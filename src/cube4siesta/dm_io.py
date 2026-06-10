"""
Read/write SIESTA .DM (density matrix) binary format.

Format (see siesta-4.1.5/Src/m_iodm.F90, subroutine write_dm, which delegates
to io_write_Sp / io_write_d2D), unformatted sequential Fortran records:

    Record 1: no_u, nspin, nsc(1:3)                    int32, 5 values
              nsc = number of supercell images along each lattice vector
              (so no_s = no_u * nsc(1) * nsc(2) * nsc(3))
    Record 2: (numd(io), io=1..no_u)                   int32, no_u values
              numd(io) = number of nonzero column entries in row io
    Then, for io = 1..no_u:
        Record:  (listd(j), j=1..numd(io))             int32
                 column indices into the *supercell* orbital set (1-based)
    Then, for ispin = 1..nspin, for io = 1..no_u:
        Record:  (dm(j,ispin), j=1..numd(io))          real*8

(The older siesta format, m_iodm_old.F, had a 2-int header without nsc; this
reader/writer targets the 4.1.5+ format actually emitted by SIESTA, which our
example .DM files use.)

The column indices in ``listd`` run over supercell orbitals (1..no_s); use
the .ORB_INDX file (see orb_indx.py) to decode each into
(atom, unit-cell orbital, cell offset).

This is the sparse, compressed-row representation SIESTA uses internally
(numd / listdptr / listd / dm). We keep the same layout so a converted DM can
be written back with an existing .DM's sparsity pattern as a template.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.io import FortranFile


@dataclass
class DMFile:
    no_u: int                 # number of orbitals in the unit cell
    nspin: int
    nsc: tuple[int, int, int]  # supercell image counts along each lattice vector
    numd: np.ndarray          # (no_u,) int, nonzeros per row
    listdptr: np.ndarray      # (no_u,) int, 0-based start of each row in flat arrays
    listd: np.ndarray         # (nnz,) int, 1-based supercell column indices
    dm: np.ndarray            # (nnz, nspin) float64, values aligned with listd

    @property
    def nnz(self) -> int:
        return int(self.numd.sum())

    @property
    def no_s(self) -> int:
        return self.no_u * self.nsc[0] * self.nsc[1] * self.nsc[2]

    def row(self, io: int) -> tuple[np.ndarray, np.ndarray]:
        """Return (listd, dm) slices for unit-cell row ``io`` (0-based)."""
        p = int(self.listdptr[io])
        n = int(self.numd[io])
        return self.listd[p : p + n], self.dm[p : p + n, :]


def _build_ptr(numd: np.ndarray) -> np.ndarray:
    listdptr = np.zeros_like(numd)
    if numd.size > 1:
        np.cumsum(numd[:-1], out=listdptr[1:])
    return listdptr


def read_dm(path: str | Path) -> DMFile:
    """Read a SIESTA .DM file (unformatted). Mirror of :func:`write_dm`."""
    with FortranFile(str(path), "r") as f:
        head = f.read_ints(dtype=np.int32)
        if head.size != 5:
            raise ValueError(
                f"header record has {head.size} ints, expected 5 (no_u, nspin, nsc(3))"
            )
        no_u, nspin = int(head[0]), int(head[1])
        nsc = (int(head[2]), int(head[3]), int(head[4]))

        numd = f.read_ints(dtype=np.int32).astype(np.int64)
        if numd.size != no_u:
            raise ValueError(f"numd record has {numd.size} ints, expected no_u={no_u}")

        listdptr = _build_ptr(numd)
        nnz = int(numd.sum())

        listd = np.empty(nnz, dtype=np.int64)
        for io in range(no_u):
            row = f.read_ints(dtype=np.int32)
            if row.size != numd[io]:
                raise ValueError(
                    f"listd row {io} has {row.size} ints, expected {numd[io]}"
                )
            p = int(listdptr[io])
            listd[p : p + row.size] = row

        dm = np.empty((nnz, nspin), dtype=np.float64)
        for ispin in range(nspin):
            for io in range(no_u):
                vals = f.read_reals(dtype=np.float64)
                if vals.size != numd[io]:
                    raise ValueError(
                        f"dm row {io} (spin {ispin}) has {vals.size} reals, "
                        f"expected {numd[io]}"
                    )
                p = int(listdptr[io])
                dm[p : p + vals.size, ispin] = vals

    return DMFile(
        no_u=no_u,
        nspin=nspin,
        nsc=nsc,
        numd=numd,
        listdptr=listdptr,
        listd=listd,
        dm=dm,
    )


def write_dm(path: str | Path, dmf: DMFile) -> None:
    """Write a SIESTA .DM file from a :class:`DMFile`."""
    no_u, nspin = dmf.no_u, dmf.nspin
    numd = np.asarray(dmf.numd, dtype=np.int64)
    if numd.size != no_u:
        raise ValueError(f"numd has {numd.size} entries, expected no_u={no_u}")
    if dmf.dm.shape != (numd.sum(), nspin):
        raise ValueError(
            f"dm shape {dmf.dm.shape} inconsistent with nnz={int(numd.sum())}, nspin={nspin}"
        )

    listdptr = _build_ptr(numd)

    nsc = np.asarray(dmf.nsc, dtype=np.int32)
    with FortranFile(str(path), "w") as f:
        f.write_record(np.array([no_u, nspin, nsc[0], nsc[1], nsc[2]], dtype=np.int32))
        f.write_record(numd.astype(np.int32))
        for io in range(no_u):
            p = int(listdptr[io])
            n = int(numd[io])
            f.write_record(np.ascontiguousarray(dmf.listd[p : p + n], dtype=np.int32))
        for ispin in range(nspin):
            for io in range(no_u):
                p = int(listdptr[io])
                n = int(numd[io])
                f.write_record(
                    np.ascontiguousarray(dmf.dm[p : p + n, ispin], dtype=np.float64)
                )


def copy_pattern(template: DMFile, nspin: int | None = None) -> DMFile:
    """
    Return a DMFile with the same sparsity pattern (numd/listdptr/listd) as
    ``template`` but zero-filled values. Use as a target for a converted DM.
    """
    ns = template.nspin if nspin is None else nspin
    return DMFile(
        no_u=template.no_u,
        nspin=ns,
        nsc=template.nsc,
        numd=template.numd.copy(),
        listdptr=template.listdptr.copy(),
        listd=template.listd.copy(),
        dm=np.zeros((template.nnz, ns), dtype=np.float64),
    )
