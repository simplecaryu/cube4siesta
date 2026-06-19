"""
Read SIESTA ``.STRUCT_OUT`` files.

The .STRUCT_OUT is a small ASCII file SIESTA writes with the final geometry:

    a1x a1y a1z            cell vectors, one per line, in Angstrom
    a2x a2y a2z
    a3x a3y a3z
    natoms
    ispec  Z  fx fy fz     one line per atom: species index, atomic number,
    ...                     and *fractional* coordinates

We only need the lattice (Angstrom) and the fractional atomic coordinates; the
Cartesian positions used for the real-space quadrature are reconstructed from
the .RHO cell (Bohr) so units stay consistent with the density grid.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np


@dataclass
class StructOut:
    cell_ang: np.ndarray   # (3, 3), lattice vectors as rows, Angstrom
    frac: np.ndarray       # (natoms, 3), fractional coordinates


def read_struct_out(path: str | Path) -> StructOut:
    """Parse a SIESTA .STRUCT_OUT file."""
    rows = [ln.split() for ln in Path(path).read_text().splitlines() if ln.strip()]
    cell_ang = np.array([[float(x) for x in rows[i]] for i in range(3)], dtype=np.float64)
    natoms = int(rows[3][0])
    frac = np.array(
        [[float(x) for x in rows[4 + a][2:5]] for a in range(natoms)],
        dtype=np.float64,
    )
    return StructOut(cell_ang=cell_ang, frac=frac)
