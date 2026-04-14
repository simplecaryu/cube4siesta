"""
Read/write SIESTA .RHO binary format.

Format (see siesta-4.1.5/Src/m_iorho.F):
    unformatted sequential Fortran records
    Record 1: cell(3,3)         real*8, 9 values, Fortran column-major
                                cell(j,i) = component j of lattice vector i
                                i.e. flat order [a1x,a1y,a1z, a2x,a2y,a2z, a3x,a3y,a3z]
    Record 2: mesh(1:3), nspin  int32, 4 values
    Record 3..:                 real*4, mesh(1) values per record
                                loop order: spin (outer) -> z -> y -> (x as record)

Units: cell in Bohr, rho in electrons/Bohr^3.
Single-process layout only (cube4siesta targets writer; SIESTA redistributes on read).
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.io import FortranFile


@dataclass
class RhoFile:
    cell: np.ndarray  # shape (3,3), float64, rows are lattice vectors (Bohr)
    mesh: tuple[int, int, int]
    nspin: int
    rho: np.ndarray  # shape (mesh[0], mesh[1], mesh[2], nspin), float32

    @property
    def n_electrons(self) -> float:
        volume = float(abs(np.linalg.det(self.cell)))
        dV = volume / (self.mesh[0] * self.mesh[1] * self.mesh[2])
        return float(self.rho.sum() * dV)


def write_rho(
    path: str | Path,
    cell: np.ndarray,
    rho: np.ndarray,
) -> None:
    """
    Write a SIESTA .RHO file.

    Parameters
    ----------
    path : output file path
    cell : (3,3) array, lattice vectors as rows, in Bohr
    rho  : (nx, ny, nz, nspin) array, density on mesh, electrons/Bohr^3
    """
    cell = np.asarray(cell, dtype=np.float64)
    if cell.shape != (3, 3):
        raise ValueError(f"cell must be (3,3), got {cell.shape}")

    rho = np.asarray(rho)
    if rho.ndim == 3:
        rho = rho[..., np.newaxis]
    if rho.ndim != 4:
        raise ValueError(f"rho must be 3D or 4D (nx,ny,nz[,nspin]), got {rho.shape}")

    nx, ny, nz, nspin = rho.shape
    mesh = np.array([nx, ny, nz], dtype=np.int32)

    # Flatten cell in the order Fortran writes cell(j,i):
    #   axis1 then axis2 then axis3, components x,y,z each.
    # Python rows = axes, C-order flatten -> [a1x,a1y,a1z, a2x,a2y,a2z, a3x,a3y,a3z]. Correct.
    cell_flat = np.ascontiguousarray(cell).reshape(9).astype(np.float64)

    rho_sp = np.ascontiguousarray(rho, dtype=np.float32)

    with FortranFile(str(path), "w") as f:
        f.write_record(cell_flat)
        f.write_record(np.concatenate([mesh, np.array([nspin], dtype=np.int32)]))
        for ispin in range(nspin):
            for iz in range(nz):
                for iy in range(ny):
                    # one record of mesh(1) real*4 values (the x-row)
                    f.write_record(rho_sp[:, iy, iz, ispin])


def read_rho(path: str | Path) -> RhoFile:
    """Read a SIESTA .RHO file. Mirror of write_rho."""
    with FortranFile(str(path), "r") as f:
        cell_flat = f.read_reals(dtype=np.float64)
        if cell_flat.size != 9:
            raise ValueError(f"cell record has {cell_flat.size} reals, expected 9")
        cell = cell_flat.reshape(3, 3)

        header = f.read_ints(dtype=np.int32)
        if header.size != 4:
            raise ValueError(f"header record has {header.size} ints, expected 4")
        nx, ny, nz, nspin = int(header[0]), int(header[1]), int(header[2]), int(header[3])

        rho = np.empty((nx, ny, nz, nspin), dtype=np.float32)
        for ispin in range(nspin):
            for iz in range(nz):
                for iy in range(ny):
                    row = f.read_reals(dtype=np.float32)
                    if row.size != nx:
                        raise ValueError(
                            f"row (spin={ispin},z={iz},y={iy}) has {row.size} reals, "
                            f"expected {nx}"
                        )
                    rho[:, iy, iz, ispin] = row

    return RhoFile(cell=cell, mesh=(nx, ny, nz), nspin=nspin, rho=rho)
