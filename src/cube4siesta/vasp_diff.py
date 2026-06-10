"""
Compute a difference density (rho_SCF - rho_atomic) cube from two VASP
CHGCAR files, for use with SIESTA's Rho.Restart.Diff mode.

Workflow: run VASP twice
  1) full SCF with LCHARG=.TRUE.  -> CHGCAR       (ρ_scf × V)
  2) ICHARG=1, NELM=0, LCHARG=.TRUE. same structure
                                  -> CHGCAR.init  (ρ_atomic_sum × V)

Then
    cube4siesta-vasp-diff scf=CHGCAR atomic=CHGCAR.init out=dden.cube

produces a cube whose integral is close to zero and whose values are
rho_scf - rho_atomic in e/Bohr^3. That cube can be fed into
'cube4siesta convert --diff ...' to produce a SIESTA .RHOIN.diff file.

Note that the atomic-superposition CHGCAR must be produced from the
*same* pseudopotential set and the *same* structure as the SCF run,
otherwise the subtraction doesn't cancel the atomic cores.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from .cube_io import CubeFile, write_cube
from .vasp_io import _ANG_TO_BOHR


def vasp_diff_chgcar(scf_path: str | Path,
                    atomic_path: str | Path,
                    out_path: str | Path,
                    *, channel: str = "total") -> CubeFile:
    """
    Write (rho_scf - rho_atomic) as a cube in e/Bohr^3.

    The two CHGCARs must share the same structure and mesh.
    """
    from pymatgen.io.vasp.outputs import Chgcar

    scf = Chgcar.from_file(str(scf_path))
    atomic = Chgcar.from_file(str(atomic_path))

    if scf.data[channel].shape != atomic.data[channel].shape:
        raise ValueError(
            f"mesh mismatch: scf={scf.data[channel].shape} "
            f"atomic={atomic.data[channel].shape}"
        )
    if scf.structure.lattice.volume != atomic.structure.lattice.volume:
        # pymatgen returns a float; strict equality is too picky. allow 1e-6.
        if abs(scf.structure.lattice.volume - atomic.structure.lattice.volume) > 1e-6:
            raise ValueError("scf and atomic CHGCARs have different cell volumes")

    volume_bohr3 = scf.structure.lattice.volume * _ANG_TO_BOHR ** 3
    diff = (np.asarray(scf.data[channel]) - np.asarray(atomic.data[channel])) \
           / volume_bohr3  # e/Bohr^3

    cell_bohr = np.asarray(scf.structure.lattice.matrix) * _ANG_TO_BOHR
    mesh = diff.shape
    voxel = cell_bohr / np.array(mesh).reshape(3, 1)

    coords_bohr = np.asarray(scf.structure.cart_coords) * _ANG_TO_BOHR
    atoms = np.zeros((len(scf.structure), 5))
    for i, site in enumerate(scf.structure):
        atoms[i, 0] = site.specie.Z
        atoms[i, 2:] = coords_bohr[i]

    cube = CubeFile(
        comment1=f"VASP diff density ({scf.structure.formula})",
        comment2=f"rho_scf - rho_atomic, channel={channel}",
        origin=np.zeros(3),
        voxel=voxel,
        mesh=mesh,
        atoms=atoms,
        data=diff,
    )
    write_cube(out_path, cube)
    return cube


def main_vasp_diff(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        prog="cube4siesta-vasp-diff",
        description=(
            "Build rho_scf - rho_atomic cube from two VASP CHGCARs "
            "for Rho.Restart.Diff mode."
        ),
    )
    p.add_argument("--scf", required=True, help="CHGCAR from a converged VASP SCF")
    p.add_argument("--atomic", required=True,
                   help="CHGCAR from VASP ICHARG=1 NELM=0 (atomic superposition)")
    p.add_argument("--out", required=True, help="output cube path")
    p.add_argument("--channel", default="total", choices=["total", "diff"],
                   help="VASP CHGCAR channel to subtract")
    args = p.parse_args(argv)

    cube = vasp_diff_chgcar(args.scf, args.atomic, args.out, channel=args.channel)

    dV = cube.voxel_volume
    n_diff = float(cube.data.sum() * dV)
    peak = float(np.max(np.abs(cube.data)))
    print(f"wrote {args.out}")
    print(f"  integral Δρ dV = {n_diff:+.3e}  (expected ~0)")
    print(f"  peak |Δρ|      = {peak:.3e} e/Bohr^3")
    return 0


if __name__ == "__main__":
    raise SystemExit(main_vasp_diff())
