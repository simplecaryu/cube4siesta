"""
Gaussian cube file I/O (reader + writer).

Format reference: http://paulbourke.net/dataformats/cube/
  line 1,2  : comments
  line 3    : natoms, origin_x, origin_y, origin_z       (Bohr)
  line 4    : N1, v1x, v1y, v1z     (voxel vector along axis 1)
  line 5    : N2, v2x, v2y, v2z
  line 6    : N3, v3x, v3y, v3z
  lines 7.. : Z, zval, x, y, z       (per atom, |natoms| lines)
  data      : rho(i1, i2, i3), i3 fastest, i2 middle, i1 slowest
                (typically 6 values per line, free format)

If natoms < 0, the cube stores a set of MO volumes preceded by an index
line — not used for density, so we error out in that case.

Units: lengths in Bohr; density in e/Bohr^3. No unit conversion here.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np


@dataclass
class CubeFile:
    comment1: str
    comment2: str
    origin: np.ndarray        # (3,) Bohr
    voxel: np.ndarray         # (3,3), voxel[i] is voxel vector along axis i (Bohr)
    mesh: tuple[int, int, int]
    atoms: np.ndarray         # (natoms, 5): Z, zval, x, y, z
    data: np.ndarray          # (N1, N2, N3) float64

    @property
    def cell(self) -> np.ndarray:
        """Full cell matrix: cell[i] = mesh[i] * voxel[i]."""
        n = np.asarray(self.mesh).reshape(3, 1)
        return self.voxel * n

    @property
    def voxel_volume(self) -> float:
        return float(abs(np.linalg.det(self.voxel)))

    @property
    def n_electrons(self) -> float:
        return float(self.data.sum() * self.voxel_volume)


def read_cube(path: str | Path) -> CubeFile:
    path = Path(path)
    with path.open() as f:
        lines = f.readlines()

    comment1 = lines[0].rstrip("\n")
    comment2 = lines[1].rstrip("\n")

    tok3 = lines[2].split()
    natoms = int(tok3[0])
    if natoms < 0:
        raise NotImplementedError(
            "cube with natoms<0 (molecular orbital volumes) not supported"
        )
    origin = np.array([float(x) for x in tok3[1:4]])

    voxel = np.zeros((3, 3))
    mesh = [0, 0, 0]
    for i in range(3):
        toks = lines[3 + i].split()
        mesh[i] = int(toks[0])
        voxel[i] = [float(x) for x in toks[1:4]]
        if mesh[i] <= 0:
            raise ValueError(
                f"cube specifies Angstrom units (negative N{i+1}); "
                "convert to Bohr first"
            )

    atoms = np.zeros((natoms, 5))
    for i in range(natoms):
        toks = lines[6 + i].split()
        atoms[i, 0] = int(toks[0])
        atoms[i, 1:] = [float(x) for x in toks[1:5]]

    # Gather remaining tokens as data (free format, any whitespace / linebreaks).
    data_tokens: list[str] = []
    for line in lines[6 + natoms:]:
        data_tokens.extend(line.split())
    total = mesh[0] * mesh[1] * mesh[2]
    if len(data_tokens) < total:
        raise ValueError(
            f"cube data truncated: got {len(data_tokens)} values, expected {total}"
        )
    raw = np.asarray(data_tokens[:total], dtype=np.float64)
    # cube order: i3 fastest, then i2, then i1. reshape accordingly.
    data = raw.reshape(mesh[0], mesh[1], mesh[2])

    return CubeFile(
        comment1=comment1,
        comment2=comment2,
        origin=origin,
        voxel=voxel,
        mesh=(mesh[0], mesh[1], mesh[2]),
        atoms=atoms,
        data=data,
    )


def write_cube(
    path: str | Path,
    cube: CubeFile,
    *,
    values_per_line: int = 6,
) -> None:
    path = Path(path)
    with path.open("w") as f:
        f.write(cube.comment1 + "\n")
        f.write(cube.comment2 + "\n")
        natoms = cube.atoms.shape[0]
        f.write(f"{natoms:5d} {cube.origin[0]:12.6f} "
                f"{cube.origin[1]:12.6f} {cube.origin[2]:12.6f}\n")
        for i in range(3):
            f.write(f"{cube.mesh[i]:5d} {cube.voxel[i,0]:12.6f} "
                    f"{cube.voxel[i,1]:12.6f} {cube.voxel[i,2]:12.6f}\n")
        for i in range(natoms):
            Z = int(cube.atoms[i, 0])
            zval, x, y, z = cube.atoms[i, 1:]
            f.write(f"{Z:5d} {zval:12.6f} {x:12.6f} {y:12.6f} {z:12.6f}\n")

        # data — i3 fastest per cube convention. flatten in that order.
        flat = cube.data.reshape(-1)
        # group values_per_line per physical line, but spec says newline after
        # each i2 inner-loop row so files read correctly in all parsers.
        n1, n2, n3 = cube.mesh
        idx = 0
        for i1 in range(n1):
            for i2 in range(n2):
                # write n3 values for this (i1,i2), wrapping every values_per_line
                for start in range(0, n3, values_per_line):
                    end = min(start + values_per_line, n3)
                    chunk = flat[idx + start:idx + end]
                    f.write(" ".join(f"{v:.5e}" for v in chunk) + "\n")
                idx += n3
