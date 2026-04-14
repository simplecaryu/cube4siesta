"""
Command-line interface: cube4siesta convert ...

Converts a Gaussian cube density into a SIESTA .RHO file, optionally
resampling onto a target mesh (derived either from a reference .RHO
or given explicitly).
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

from .cube_io import read_cube
from .resample import resample_to_mesh
from .rho_io import read_rho, write_rho


def _parse_mesh(spec: str) -> tuple[int, int, int]:
    parts = spec.split(",")
    if len(parts) != 3:
        raise argparse.ArgumentTypeError(f"mesh must be Nx,Ny,Nz, got {spec!r}")
    return tuple(int(p) for p in parts)  # type: ignore[return-value]


def cmd_convert(args: argparse.Namespace) -> int:
    cube = read_cube(args.cube)
    print(f"[cube4siesta] loaded {args.cube}")
    print(f"    mesh = {cube.mesh}")
    print(f"    cell (Bohr) =\n{cube.cell}")
    print(f"    integral rho dV = {cube.n_electrons:.6f} electrons")

    # Target mesh
    if args.from_siesta_rho is not None:
        ref = read_rho(args.from_siesta_rho)
        target_mesh = ref.mesh
        target_cell = ref.cell
        print(f"[cube4siesta] target mesh from {args.from_siesta_rho}: {target_mesh}")
        # Warn if cell disagrees
        if not np.allclose(ref.cell, cube.cell, atol=1e-3):
            print("WARNING: reference .RHO cell differs from cube cell by > 1e-3 Bohr")
            print(f"    cube cell =\n{cube.cell}")
            print(f"    ref  cell =\n{ref.cell}")
    elif args.target_mesh is not None:
        target_mesh = args.target_mesh
        target_cell = cube.cell
        print(f"[cube4siesta] target mesh (user): {target_mesh}")
    else:
        target_mesh = cube.mesh
        target_cell = cube.cell
        print(f"[cube4siesta] target mesh = cube mesh: {target_mesh}")

    # Resample if necessary
    if target_mesh == cube.mesh:
        rho = cube.data.astype(np.float32)
        print("[cube4siesta] grids match — no resampling")
    else:
        rho = resample_to_mesh(cube.data, target_mesh, order=args.order).astype(np.float32)
        print(f"[cube4siesta] resampled {cube.mesh} -> {target_mesh} (order={args.order})")

    # Charge normalization
    volume = float(abs(np.linalg.det(target_cell)))
    dV = volume / (target_mesh[0] * target_mesh[1] * target_mesh[2])
    n_after = float(rho.sum() * dV)
    print(f"[cube4siesta] integral after resample = {n_after:.6f} electrons")

    if args.rescale_to is not None:
        if n_after <= 0:
            print("ERROR: cannot rescale, integrated rho is non-positive", file=sys.stderr)
            return 2
        scale = args.rescale_to / n_after
        rho *= scale
        print(f"[cube4siesta] rescaled by {scale:.6f} -> {args.rescale_to:.6f} electrons")

    # Add a spin axis (nspin=1)
    rho_4d = rho[..., np.newaxis]
    write_rho(args.output, target_cell, rho_4d)
    print(f"[cube4siesta] wrote {args.output}")

    if args.verify:
        back = read_rho(args.output)
        assert back.mesh == target_mesh
        diff = float(np.max(np.abs(back.rho[..., 0] - rho)))
        print(f"[cube4siesta] verify: max |read - wrote| = {diff:.3g}")
        if diff > 1e-5:
            print("WARNING: roundtrip mismatch exceeds 1e-5", file=sys.stderr)
            return 3

    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="cube4siesta")
    sub = parser.add_subparsers(dest="cmd", required=True)

    c = sub.add_parser("convert", help="Convert a Gaussian cube to SIESTA .RHO")
    c.add_argument("--cube", required=True, help="input Gaussian cube file")
    c.add_argument("--output", required=True, help="output .RHO path (e.g. system.RHOIN)")
    grp = c.add_mutually_exclusive_group()
    grp.add_argument("--target-mesh", type=_parse_mesh,
                     help="target mesh as Nx,Ny,Nz (else cube's own mesh)")
    grp.add_argument("--from-siesta-rho", help="borrow target mesh+cell from an existing .RHO")
    c.add_argument("--order", type=int, default=1,
                   help="interpolation order (1=trilinear, 3=cubic)")
    c.add_argument("--rescale-to", type=float,
                   help="rescale density so integral = this many electrons")
    c.add_argument("--verify", action="store_true",
                   help="read back the written .RHO and check roundtrip")
    c.set_defaults(func=cmd_convert)

    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
