"""
Read OpenMX ``.scfout`` files (SCFOUT_VERSION 3, OpenMX 3.9).

This is a straight port of the binary layout in
``openmx3.9/source/read_scfout.c``. Note the file is a *raw* C binary stream
(plain ``fread`` of ints/doubles, NO Fortran record-length markers), so we read
sequentially with ``numpy.fromfile`` using an endianness detected from the
header (mirroring the ``conversionSwitch`` logic in the C reader).

Layout (1-based atom / neighbour indexing, matching the C source):

    i_vec[0..5]  : atomnum, (SpinP_switch + 4*version), Catomnum, Latomnum,
                   Ratomnum, TCpyCell
    order_max    : 1 int
    atv          : (TCpyCell+1) x 4 doubles   (cell translation vectors, Bohr;
                                               components in cols 1..3)
    atv_ijk      : (TCpyCell+1) x 4 ints       (integer cell indices, cols 1..3)
    Total_NumOrbs: atomnum ints
    FNAN         : atomnum ints                (# first-nearest neighbours)
    natn[ct][0..FNAN]  : global neighbour atom indices (1-based)
    ncn[ct][0..FNAN]   : cell index Rn into atv / atv_ijk
    tv  : 3 x 4 doubles    (lattice vectors, Bohr; rows 1..3, cols 1..3)
    rtv : 3 x 4 doubles    (reciprocal lattice)
    Gxyz[ct][0..3] : atomic coordinates (Bohr; cols 1..3), for ct=1..atomnum
    Hks  : per (spin, ct, h, i) -> TNO2 doubles   [we read & discard]
    iHks : only if SpinP_switch==3                [discard]
    OLP  : per (ct, h, i) -> TNO2 doubles         [kept]
    OLPpo: 3 x order_max x (ct,h,i)               [discard]
    OLPmo: 3 x (ct,h,i)                           [discard]
    DM   : per (spin, ct, h, i) -> TNO2 doubles   [kept]
    iDM  : spin 0..1 (ct,h,i)                     [discard]
    Solver (1 int), d_vec[10] (ChemP, E_Temp, dipoles, Valence_Electrons,
    Total_SpinS), then the embedded input file.

For a neighbour pair (ct, h): the neighbour's unit-cell atom is
``Gh = natn[ct][h]`` and its cell offset is ``atv_ijk[ncn[ct][h]]`` (integer)
/ ``atv[ncn[ct][h]]`` (Bohr). The on-site pair is h == 0.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np


@dataclass
class Scfout:
    atomnum: int
    spin_switch: int          # 0 = non-pol, 1 = collinear, 3 = non-collinear
    version: int
    tcpycell: int
    order_max: int
    total_numorbs: np.ndarray  # (atomnum+1,) int, [0] unused
    fnan: np.ndarray           # (atomnum+1,) int, [0] unused
    natn: list[np.ndarray]     # natn[ct] -> (FNAN[ct]+1,) global atom indices (1-based)
    ncn: list[np.ndarray]      # ncn[ct]  -> (FNAN[ct]+1,) cell index Rn
    atv: np.ndarray            # (TCpyCell+1, 4) doubles, cols 1..3 = xyz (Bohr)
    atv_ijk: np.ndarray        # (TCpyCell+1, 4) ints,    cols 1..3 = ijk
    tv: np.ndarray             # (4, 4) doubles, rows/cols 1..3 = lattice (Bohr)
    gxyz: np.ndarray           # (atomnum+1, 4) doubles, cols 1..3 = xyz (Bohr)
    OLP: dict                  # {(ct, h): (TNO1, TNO2) ndarray}
    DM: list                   # DM[spin] = {(ct, h): (TNO1, TNO2) ndarray}
    chemp: float
    valence_electrons: float
    total_spin: float

    @property
    def nspin(self) -> int:
        # number of independent DM spin blocks stored
        return {0: 1, 1: 2, 3: 4}[self.spin_switch]

    def cell_bohr(self) -> np.ndarray:
        """(3,3) lattice vectors as rows, Bohr."""
        return self.tv[1:4, 1:4].copy()

    def atom_xyz(self) -> np.ndarray:
        """(atomnum, 3) atomic coordinates, Bohr (0-based atom order)."""
        return self.gxyz[1:, 1:4].copy()

    def neighbor_cell_ijk(self, ct: int, h: int) -> tuple[int, int, int]:
        """Integer cell offset (i,j,k) of neighbour h of central atom ct."""
        rn = int(self.ncn[ct][h])
        v = self.atv_ijk[rn]
        return int(v[1]), int(v[2]), int(v[3])


class _Reader:
    """Sequential raw-binary reader with a fixed endianness."""

    def __init__(self, path: str | Path):
        self.fp = open(path, "rb")
        self._int = "<i4"
        self._dbl = "<f8"

    def set_endian(self, big: bool) -> None:
        self._int = ">i4" if big else "<i4"
        self._dbl = ">f8" if big else "<f8"

    def ints(self, n: int) -> np.ndarray:
        return np.fromfile(self.fp, dtype=self._int, count=n)

    def dbls(self, n: int) -> np.ndarray:
        return np.fromfile(self.fp, dtype=self._dbl, count=n)

    def close(self) -> None:
        self.fp.close()


def _read_pair_matrices(
    r: _Reader,
    spins: int,
    atomnum: int,
    fnan: np.ndarray,
    natn: list[np.ndarray],
    total_numorbs: np.ndarray,
    keep: bool,
) -> list[dict] | None:
    """
    Read the jagged per-(ct,h,i) row blocks for ``spins`` spin channels.

    For each spin, ct=1..atomnum, h=0..FNAN[ct], i=0..TNO1-1: TNO2 doubles,
    where TNO1 = Total_NumOrbs[ct], TNO2 = Total_NumOrbs[natn[ct][h]].
    Returns a list (per spin) of dicts {(ct,h): (TNO1,TNO2) array} if ``keep``,
    else consumes the bytes and returns None.
    """
    out: list[dict] | None = [dict() for _ in range(spins)] if keep else None
    for s in range(spins):
        for ct in range(1, atomnum + 1):
            tno1 = int(total_numorbs[ct])
            for h in range(int(fnan[ct]) + 1):
                gh = int(natn[ct][h])
                tno2 = int(total_numorbs[gh])
                block = r.dbls(tno1 * tno2)
                if keep:
                    out[s][(ct, h)] = block.reshape(tno1, tno2)
    return out


def read_scfout(path: str | Path) -> Scfout:
    """Parse an OpenMX ``.scfout`` (version 3) file."""
    r = _Reader(path)
    try:
        i_vec = r.ints(6)
        # Endianness check, mirroring read_scfout.c: i_vec[1] must be in a
        # plausible range; otherwise the file is the other endianness.
        if not (0 <= int(i_vec[1]) <= 3 * 4 + 3):
            r.set_endian(big=True)
            i_vec = i_vec.byteswap()
            if not (0 <= int(i_vec[1]) <= 3 * 4 + 3):
                raise ValueError("scfout: cannot determine endianness / unsupported version")

        atomnum = int(i_vec[0])
        spin_switch = int(i_vec[1]) % 4
        version = int(i_vec[1]) // 4
        tcpycell = int(i_vec[5])
        if version != 3:
            raise ValueError(f"scfout version {version} unsupported (expected 3)")

        order_max = int(r.ints(1)[0])

        atv = r.dbls((tcpycell + 1) * 4).reshape(tcpycell + 1, 4)
        atv_ijk = r.ints((tcpycell + 1) * 4).reshape(tcpycell + 1, 4)

        total_numorbs = np.empty(atomnum + 1, dtype=np.int64)
        total_numorbs[0] = 1
        total_numorbs[1:] = r.ints(atomnum)

        fnan = np.empty(atomnum + 1, dtype=np.int64)
        fnan[0] = 0
        fnan[1:] = r.ints(atomnum)

        natn: list[np.ndarray] = [np.array([0], dtype=np.int64)]
        for ct in range(1, atomnum + 1):
            natn.append(r.ints(int(fnan[ct]) + 1).astype(np.int64))
        ncn: list[np.ndarray] = [np.array([0], dtype=np.int64)]
        for ct in range(1, atomnum + 1):
            ncn.append(r.ints(int(fnan[ct]) + 1).astype(np.int64))

        tv = np.zeros((4, 4), dtype=np.float64)
        for j in (1, 2, 3):
            tv[j] = r.dbls(4)
        rtv = np.zeros((4, 4), dtype=np.float64)
        for j in (1, 2, 3):
            rtv[j] = r.dbls(4)  # read & discard (kept local for completeness)

        gxyz = np.zeros((atomnum + 1, 4), dtype=np.float64)
        for ct in range(1, atomnum + 1):
            gxyz[ct] = r.dbls(4)

        nspin_h = spin_switch + 1  # Hks/DM stored for spin = 0..SpinP_switch

        # Hks  (discard)
        _read_pair_matrices(r, nspin_h, atomnum, fnan, natn, total_numorbs, keep=False)
        # iHks (only non-collinear; discard)
        if spin_switch == 3:
            _read_pair_matrices(r, 3, atomnum, fnan, natn, total_numorbs, keep=False)
        # OLP (keep)
        olp_list = _read_pair_matrices(r, 1, atomnum, fnan, natn, total_numorbs, keep=True)
        OLP = olp_list[0]
        # OLPpo : 3 directions x order_max  (discard)
        if order_max > 0:
            _read_pair_matrices(
                r, 3 * order_max, atomnum, fnan, natn, total_numorbs, keep=False
            )
        # OLPmo : 3 directions (discard)
        _read_pair_matrices(r, 3, atomnum, fnan, natn, total_numorbs, keep=False)
        # DM (keep)
        DM = _read_pair_matrices(r, nspin_h, atomnum, fnan, natn, total_numorbs, keep=True)
        # iDM : spin 0..1 (discard)
        _read_pair_matrices(r, 2, atomnum, fnan, natn, total_numorbs, keep=False)

        solver = int(r.ints(1)[0])  # noqa: F841 (read to advance; not stored)
        d_vec = r.dbls(10)
        chemp = float(d_vec[0])
        valence_electrons = float(d_vec[8])
        total_spin = float(d_vec[9])
    finally:
        r.close()

    return Scfout(
        atomnum=atomnum,
        spin_switch=spin_switch,
        version=version,
        tcpycell=tcpycell,
        order_max=order_max,
        total_numorbs=total_numorbs,
        fnan=fnan,
        natn=natn,
        ncn=ncn,
        atv=atv,
        atv_ijk=atv_ijk,
        tv=tv,
        gxyz=gxyz,
        OLP=OLP,
        DM=DM,
        chemp=chemp,
        valence_electrons=valence_electrons,
        total_spin=total_spin,
    )


def trace_PS(sc: Scfout, spin: int = 0) -> float:
    """
    Electron count in the OpenMX basis: sum over neighbour pairs of
    DM[spin][ct,h] : OLP[ct,h] (elementwise). For non-spin-polarised data
    multiply by 2 externally if needed; here returns the raw per-channel trace.
    """
    total = 0.0
    for (ct, h), dm in sc.DM[spin].items():
        olp = sc.OLP[(ct, h)]
        total += float(np.sum(dm * olp))
    return total
