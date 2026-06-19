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
    elif args.mesh_cutoff is not None:
        from .mesh import siesta_mesh
        target_mesh = siesta_mesh(cube.cell, mesh_cutoff=args.mesh_cutoff)
        target_cell = cube.cell
        print(f"[cube4siesta] target mesh from MeshCutoff={args.mesh_cutoff} Ry: "
              f"{target_mesh}")
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

    # Subtract a reference density if requested (build a difference density
    # from two totals, e.g. source SCF minus source atomic superposition).
    if args.subtract is not None:
        ref_path = Path(args.subtract)
        if ref_path.suffix.lower() in (".rho", ".rhoin", ".rhoref"):
            ref = read_rho(ref_path)
            ref_data = ref.rho[..., 0]
            ref_mesh = ref.mesh
        else:
            ref_cube = read_cube(ref_path)
            ref_data = ref_cube.data
            ref_mesh = ref_cube.mesh
        if tuple(ref_mesh) != tuple(target_mesh):
            ref_data = resample_to_mesh(ref_data, target_mesh, order=args.order)
        rho = (rho - ref_data).astype(np.float32)
        print(f"[cube4siesta] subtracted reference density from {args.subtract}")

    # Charge normalization
    volume = float(abs(np.linalg.det(target_cell)))
    dV = volume / (target_mesh[0] * target_mesh[1] * target_mesh[2])
    n_after = float(rho.sum() * dV)

    if args.diff or args.subtract is not None:
        # Difference density: integral should be close to zero. SIESTA's
        # Rho.Restart.Diff mode adds its own rhoatm back in on its side.
        print(f"[cube4siesta] diff-mode: integral ~ {n_after:+.3e} "
              "(expected ~0 for a difference density)")
        if abs(n_after) > 0.1:
            print("WARNING: |integral| > 0.1 for a diff cube — are you sure this"
                  " is rho - rho_atomic and not a total density?", file=sys.stderr)
        if args.rescale_to is not None:
            print("ERROR: --rescale-to is incompatible with a difference density "
                  "(it must remain charge-neutral)", file=sys.stderr)
            return 2
    else:
        print(f"[cube4siesta] integral after resample = {n_after:.6f} electrons")
        if args.rescale_to is not None:
            if n_after <= 0:
                print("ERROR: cannot rescale, integrated rho is non-positive",
                      file=sys.stderr)
                return 2
            scale = args.rescale_to / n_after
            rho *= scale
            print(f"[cube4siesta] rescaled by {scale:.6f} -> "
                  f"{args.rescale_to:.6f} electrons")

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


def cmd_gen_atom_input(args: argparse.Namespace) -> int:
    from .gen_psf import parse_openmx_vps, parse_qe_upf, write_atom_input

    if args.from_vps:
        params = parse_openmx_vps(args.from_vps)
        source = args.from_vps
    else:
        params = parse_qe_upf(args.from_upf)
        source = args.from_upf

    content = write_atom_input(params, args.output)

    print(f"[cube4siesta] parsed {source}")
    print(f"    element: {params.symbol} (Z={params.z})")
    print(f"    XC: {params.xc}  relativistic: {params.relativistic}")
    print(f"    core shells: {params.n_core_shells}")
    print(f"    valence: {params.n_valence_electrons:.0f} electrons in "
          f"{len(params.valence)} channels")
    for v in params.valence:
        print(f"      {v.n}{_L_LABELS[v.l]}  occ={v.occ:.2f}  rc={v.rc:.2f}")
    print(f"[cube4siesta] wrote {args.output}")
    print(f"\nTo generate .psf, run:  sh /path/to/atom/bin/pg.sh {args.output}")
    return 0


_L_LABELS = "spdf"


def _parse_labeled(values: list[str] | None) -> dict[str, str]:
    """Parse repeated LABEL:PATH (or LABEL=SPEC) options into a dict."""
    out: dict[str, str] = {}
    for item in values or []:
        if ":" in item:
            label, val = item.split(":", 1)
        elif "=" in item:
            label, val = item.split("=", 1)
        else:
            raise argparse.ArgumentTypeError(f"expected LABEL:VALUE, got {item!r}")
        out[label] = val
    return out


def cmd_project_dm(args: argparse.Namespace) -> int:
    from .basis import (build_openmx_orbitals, build_siesta_orbitals,
                        parse_spec, read_ion_radials, read_pao_radials)
    from .dm_io import read_dm, write_dm
    from .orb_indx import read_orb_indx
    from .project import project_openmx_to_siesta
    from .scfout_io import read_scfout

    scfout = read_scfout(args.scfout)
    oi = read_orb_indx(args.orb_indx)
    template = read_dm(args.template)
    print(f"[cube4siesta] scfout: {scfout.atomnum} atoms, spin_switch={scfout.spin_switch}, "
          f"valence_electrons={scfout.valence_electrons:.1f}")
    print(f"[cube4siesta] template .DM: no_u={template.no_u}, nspin={template.nspin}, "
          f"nsc={template.nsc}, nnz={template.nnz}")

    ion_map = _parse_labeled(args.ion)
    pao_map = _parse_labeled(args.pao)
    spec_map = _parse_labeled(args.spec)

    # species per atom (1-based) from .ORB_INDX, assuming OpenMX atom order matches
    species_of_atom: dict[int, str] = {}
    for o in oi.unit:
        species_of_atom[o.ia] = o.spec
    n_atoms = len(species_of_atom)
    if n_atoms != scfout.atomnum:
        print(f"ERROR: atom count mismatch: ORB_INDX={n_atoms}, scfout={scfout.atomnum}",
              file=sys.stderr)
        return 2

    # SIESTA radial tables per species
    radials_by_species = {lab: read_ion_radials(p) for lab, p in ion_map.items()}
    missing = {o.spec for o in oi.unit} - set(radials_by_species)
    if missing:
        print(f"ERROR: missing --ion for species {missing}", file=sys.stderr)
        return 2
    siesta_orbitals = build_siesta_orbitals(oi.unit, radials_by_species)

    # OpenMX orbitals per species, then per atom
    omx_by_species = {}
    for lab in {species_of_atom[a] for a in species_of_atom}:
        if lab not in pao_map or lab not in spec_map:
            print(f"ERROR: missing --pao/--spec for species {lab}", file=sys.stderr)
            return 2
        omx_by_species[lab] = build_openmx_orbitals(
            parse_spec(spec_map[lab]), read_pao_radials(pao_map[lab])
        )
    omx_per_atom = [omx_by_species[species_of_atom[ct]] for ct in range(1, n_atoms + 1)]

    print(f"[cube4siesta] projecting (spacing={args.spacing} Bohr) ...")
    res = project_openmx_to_siesta(
        scfout, oi, template, siesta_orbitals, omx_per_atom,
        spacing=args.spacing, verbose=args.verbose,
        purify=args.purify, qtot=args.qtot,
    )
    write_dm(args.output, res.dm)
    print(f"[cube4siesta] wrote {args.output}")
    print(f"    Tr(P_SIE S) = {res.n_electrons:.4f} electrons "
          f"(SIESTA-representable; semicore projects out)")
    if res.purified:
        print(f"    purified: natural occupations were in "
              f"[{res.occ_min:.4f}, {res.occ_max:.4f}] -> idempotent filling, "
              f"{res.n_frac} fractional state(s) at the Fermi boundary")
    print(f"    Hermiticity residual = {res.max_asym:.2e}  (k-points: {res.nk})")
    return 0


def cmd_compare_dm(args: argparse.Namespace) -> int:
    from .dm_compare import compare, format_table
    from .dm_io import read_dm
    from .orb_indx import read_orb_indx

    oi = read_orb_indx(args.orb_indx)
    ref = read_dm(args.ref)
    cand_map = _parse_labeled(args.cand)
    results = {}
    for name, path in cand_map.items():
        results[name] = compare(ref, read_dm(path), oi, sig_thresh=args.sig_thresh)
    print(f"[cube4siesta] reference: {args.ref}  (no_u={ref.no_u}, nnz={ref.nnz})")
    print(f"[cube4siesta] DM-value similarity vs reference "
          f"(sig threshold |ref|>{args.sig_thresh}):\n")
    print(format_table(results))
    return 0


def cmd_dm_seed(args: argparse.Namespace) -> int:
    from .density_fit import make_filtered_dm_seed, restart_fdf_block
    from .dm_io import write_dm

    res = make_filtered_dm_seed(
        base_dir=args.baseline,
        stem=args.stem,
        dden_path=args.dden,
        damp=args.damp,
        verbose=args.verbose,
    )
    write_dm(args.output, res.dm)
    print(f"[cube4siesta] wrote filtered DM seed -> {args.output}")
    print(f"    nnz={res.nnz}  max|D|={res.max_abs:.3f}  "
          f"filtered Delta-rho integral={res.filtered_integral:+.4f} e")
    print(f"    seed vs converged DM (pre-SCF, diagnostic only): "
          f"L2={res.direct_l2:.4f}")
    print("    -> one SIESTA diagonalisation refines this seed. Add to your fdf:\n")
    print(restart_fdf_block())
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="cube4siesta")
    sub = parser.add_subparsers(dest="cmd", required=True)

    # --- dm-seed subcommand (primary cross-code restart route) ---
    ds = sub.add_parser(
        "dm-seed",
        help="MAIN ROUTE: build a SIESTA DM restart seed by filtering a "
        "cross-code difference density onto the SIESTA basis (semicore-safe)",
    )
    ds.add_argument("--baseline", required=True,
                    help="directory holding the SIESTA baseline run "
                    "(<stem>.DM/.ORB_INDX/.STRUCT_OUT/.RHO and <species>.ion)")
    ds.add_argument("--stem", required=True,
                    help="SIESTA SystemLabel of the baseline run (e.g. vsse)")
    ds.add_argument("--dden", required=True,
                    help="cross-code difference density rho_scf - rho_atoms "
                    "(.cube, e.g. OpenMX *.dden.cube, or a .RHO)")
    ds.add_argument("--output", required=True, help="output seed .DM path")
    ds.add_argument("--damp", type=float, default=0.1,
                    help="least-squares damping for the basis filter (default 0.1)")
    ds.add_argument("--verbose", action="store_true")
    ds.set_defaults(func=cmd_dm_seed)

    c = sub.add_parser("convert", help="Convert a Gaussian cube to SIESTA .RHO")
    c.add_argument("--cube", required=True, help="input Gaussian cube file")
    c.add_argument("--output", required=True, help="output .RHO path (e.g. system.RHOIN)")
    grp = c.add_mutually_exclusive_group()
    grp.add_argument("--target-mesh", type=_parse_mesh,
                     help="target mesh as Nx,Ny,Nz (else cube's own mesh)")
    grp.add_argument("--from-siesta-rho", help="borrow target mesh+cell from an existing .RHO")
    grp.add_argument("--mesh-cutoff", type=float,
                     help="predict SIESTA mesh from MeshCutoff in Ry (no SIESTA run needed)."
                     " Uses the cube's cell vectors.")
    c.add_argument("--order", type=int, default=1,
                   help="interpolation order (1=trilinear, 3=cubic)")
    c.add_argument("--rescale-to", type=float,
                   help="rescale density so integral = this many electrons")
    c.add_argument("--diff", action="store_true",
                   help="the cube already holds a difference density "
                   "(rho_scf - rho_atomic, e.g. OpenMX *.dden.cube); convert "
                   "as-is for SIESTA's Rho.Restart.Diff mode")
    c.add_argument("--subtract", metavar="REF",
                   help="build a difference density by subtracting a reference "
                   "cube or .RHO (the source code's atomic-superposition "
                   "density) before writing; implies diff-mode checks")
    c.add_argument("--verify", action="store_true",
                   help="read back the written .RHO and check roundtrip")
    c.set_defaults(func=cmd_convert)

    # --- gen-atom-input subcommand ---
    g = sub.add_parser("gen-atom-input",
                       help="Generate ATOM .inp from a source pseudo file")
    src = g.add_mutually_exclusive_group(required=True)
    src.add_argument("--from-vps", help="OpenMX .vps file")
    src.add_argument("--from-upf", help="Quantum ESPRESSO .UPF file")
    g.add_argument("--output", required=True, help="output ATOM .inp path")
    g.set_defaults(func=cmd_gen_atom_input)

    # --- project-dm subcommand ---
    p = sub.add_parser(
        "project-dm",
        help="Project an OpenMX density matrix onto SIESTA's NAO basis (.scfout -> .DM)",
    )
    p.add_argument("--scfout", required=True, help="OpenMX .scfout file")
    p.add_argument("--orb-indx", required=True, help="SIESTA .ORB_INDX file")
    p.add_argument("--template", required=True,
                   help="SIESTA .DM whose sparsity pattern is reused for output")
    p.add_argument("--ion", action="append", metavar="LABEL:PATH",
                   help="SIESTA .ion file for a species (repeatable), e.g. Si:Si.ion")
    p.add_argument("--pao", action="append", metavar="LABEL:PATH",
                   help="OpenMX .pao file for a species (repeatable), e.g. Si:Si7.0.pao")
    p.add_argument("--spec", action="append", metavar="LABEL:SPEC",
                   help="OpenMX basis spec for a species (repeatable), e.g. Si:s2p2d1")
    p.add_argument("--spacing", type=float, default=0.15,
                   help="real-space grid spacing in Bohr for overlap quadrature")
    p.add_argument("--purify", action="store_true",
                   help="canonical (McWeeny-limit) purification: make the DM "
                   "idempotent and fix Tr(PS) to --qtot")
    p.add_argument("--qtot", type=float,
                   help="target electron count for --purify (default: raw "
                   "Tr(PS) rounded to the nearest integer)")
    p.add_argument("--output", required=True, help="output SIESTA .DM path")
    p.add_argument("--verbose", action="store_true")
    p.set_defaults(func=cmd_project_dm)

    # --- compare-dm subcommand ---
    cm = sub.add_parser("compare-dm", help="Compare SIESTA .DM files element-wise")
    cm.add_argument("--orb-indx", required=True, help="SIESTA .ORB_INDX file")
    cm.add_argument("--ref", required=True, help="reference .DM (e.g. SIESTA baseline)")
    cm.add_argument("--cand", action="append", required=True, metavar="NAME:PATH",
                    help="candidate .DM (repeatable), e.g. projection:si_omx_proj.DM")
    cm.add_argument("--sig-thresh", type=float, default=0.02,
                    help="|ref| threshold for the 'significant element' correlation")
    cm.set_defaults(func=cmd_compare_dm)

    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
