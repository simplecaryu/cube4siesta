"""
VASP CHGCAR -> Gaussian cube conversion.

VASP's CHGCAR stores rho(r) * V_cell (volume-weighted) in Angstrom units.
Gaussian cube convention is rho(r) in e/Bohr^3 with Bohr lengths.
We use pymatgen only to parse the binary-ish CHGCAR (fiddly to parse by
hand) and then handle unit conversion + cube emission via our own
cube_io to keep the pipeline symmetric.

Note: pymatgen.Chgcar.to_cube() does NOT divide by volume and so produces
cubes that integrate to N_e * V (not N_e). Do not use it directly.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from .cube_io import CubeFile, write_cube

# Conversion constant consistent with pymatgen's ang_to_bohr (CODATA 2018)
_ANG_TO_BOHR = 1.8897261258369282


def chgcar_to_cube(chgcar_path: str | Path, cube_path: str | Path,
                   *, channel: str = "total") -> CubeFile:
    """
    Parse a VASP CHGCAR and write a properly normalized Gaussian cube.

    Parameters
    ----------
    chgcar_path : input CHGCAR path
    cube_path   : output .cube path
    channel     : "total" (default) or "diff" for the magnetization density of
                  a spin-polarized CHGCAR.

    Returns
    -------
    CubeFile : the in-memory cube that was written.
    """
    from pymatgen.io.vasp.outputs import Chgcar

    chg = Chgcar.from_file(str(chgcar_path))
    if channel not in chg.data:
        raise KeyError(
            f"CHGCAR has no '{channel}' channel; available: {list(chg.data.keys())}"
        )

    # CHGCAR raw values are rho(r) * V. Divide out to get e/Bohr^3.
    # Volume from the structure's lattice is in Ang^3 -> convert to Bohr^3.
    volume_ang3 = float(chg.structure.lattice.volume)
    volume_bohr3 = volume_ang3 * _ANG_TO_BOHR ** 3
    rho = np.asarray(chg.data[channel], dtype=np.float64) / volume_bohr3

    # Cell and mesh
    cell_ang = np.asarray(chg.structure.lattice.matrix)
    cell_bohr = cell_ang * _ANG_TO_BOHR
    mesh = rho.shape
    voxel = cell_bohr / np.array(mesh).reshape(3, 1)

    # Atoms (Z, zval=0 since we don't have valence info from CHGCAR alone,
    # xyz in Bohr)
    coords_bohr = np.asarray(chg.structure.cart_coords) * _ANG_TO_BOHR
    atoms = np.zeros((len(chg.structure), 5))
    for i, site in enumerate(chg.structure):
        atoms[i, 0] = site.specie.Z
        atoms[i, 1] = 0.0
        atoms[i, 2:] = coords_bohr[i]

    cube = CubeFile(
        comment1=f"VASP CHGCAR -> cube ({chg.structure.formula})",
        comment2=f"channel={channel}",
        origin=np.zeros(3),
        voxel=voxel,
        mesh=mesh,
        atoms=atoms,
        data=rho,
    )
    write_cube(cube_path, cube)
    return cube
