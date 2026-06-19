"""High-level automation for OpenMX ``.scfout`` -> SIESTA ``.DM`` projection.

This module provides a pragmatic wrapper around the explicit ``project-dm``
implementation: parse the embedded OpenMX input in ``.scfout``, make a minimal
SIESTA template run to obtain ``.ion/.ORB_INDX/.DM``, then call the same overlap
projection code used by the explicit CLI.
"""

from __future__ import annotations

import argparse
import os
import re
import shlex
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

from .basis import (build_openmx_orbitals, build_siesta_orbitals, parse_spec,
                    read_ion_radials, read_pao_radials)
from .dm_io import read_dm, write_dm
from .orb_indx import read_orb_indx
from .project import ProjectionResult, project_openmx_to_siesta
from .scfout_io import Scfout, read_scfout

# Element symbols, 1-based atomic numbers.
_SYMBOLS = [
    "", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg",
    "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr",
    "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr",
    "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
    "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
    "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf",
    "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po",
    "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm",
    "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr",
]
Z_OF_SYMBOL = {s: i for i, s in enumerate(_SYMBOLS) if s}


@dataclass
class OpenMXSpecies:
    label: str
    pao_tag: str          # e.g. Mo7.0-s3p2d2f1
    vps_tag: str

    @property
    def pao_file(self) -> str:
        return self.pao_tag.split("-", 1)[0] + ".pao"

    @property
    def spec(self) -> str:
        if "-" not in self.pao_tag:
            raise ValueError(f"cannot infer basis spec from OpenMX PAO tag {self.pao_tag!r}")
        return self.pao_tag.split("-", 1)[1]


@dataclass
class OpenMXInputInfo:
    system_name: str
    data_path: Path | None
    species: list[OpenMXSpecies]
    atom_labels: list[str]
    kgrid: tuple[int, int, int] = (1, 1, 1)


def _clean_openmx_embedded_lines(text: str) -> list[str]:
    """Extract original input-like lines from OpenMX's noisy embedded strings."""
    out: list[str] = []
    for raw in text.splitlines():
        line = raw.strip()
        if not line.startswith("*"):
            continue
        # OpenMX prints embedded input lines prefixed by many stars.
        line = line.lstrip("*").strip()
        if not line or set(line) == {"*"}:
            continue
        out.append(line)
    return out


def _block(lines: list[str], start: str, end: str) -> list[str]:
    inside = False
    out: list[str] = []
    for line in lines:
        if line.startswith(start):
            inside = True
            continue
        if line.startswith(end):
            break
        if inside:
            out.append(line)
    return out


def parse_openmx_input_info(text: str) -> OpenMXInputInfo:
    """Parse the small subset of OpenMX input needed for automation."""
    lines = _clean_openmx_embedded_lines(text)
    if not lines:
        raise ValueError("could not find embedded OpenMX input in .scfout")

    system_name = "system"
    data_path: Path | None = None
    kgrid = (1, 1, 1)
    for line in lines:
        toks = line.split()
        if len(toks) >= 2 and toks[0] == "System.Name":
            system_name = toks[1]
        elif len(toks) >= 2 and toks[0] == "DATA.PATH":
            data_path = Path(toks[1])
        elif len(toks) >= 4 and toks[0] == "scf.Kgrid":
            kgrid = (int(toks[1]), int(toks[2]), int(toks[3]))

    species: list[OpenMXSpecies] = []
    for line in _block(lines, "<Definition.of.Atomic.Species", "Definition.of.Atomic.Species>"):
        toks = line.split()
        if len(toks) >= 3:
            species.append(OpenMXSpecies(toks[0], toks[1], toks[2]))
    if not species:
        raise ValueError("could not parse Definition.of.Atomic.Species from .scfout")

    atom_labels: list[str] = []
    for line in _block(lines, "<Atoms.SpeciesAndCoordinates", "Atoms.SpeciesAndCoordinates>"):
        toks = line.split()
        if len(toks) >= 2 and toks[0].lstrip("+-").isdigit():
            atom_labels.append(toks[1])
    if not atom_labels:
        raise ValueError("could not parse Atoms.SpeciesAndCoordinates from .scfout")

    return OpenMXInputInfo(system_name, data_path, species, atom_labels, kgrid)


def _link_or_copy(src: Path, dst: Path) -> None:
    if dst.exists() or dst.is_symlink():
        return
    try:
        os.symlink(src, dst)
    except OSError:
        shutil.copy2(src, dst)


def write_siesta_template_fdf(
    path: Path,
    scfout: Scfout,
    info: OpenMXInputInfo,
    label: str,
    basis_size: str,
    mesh_cutoff: float,
    kgrid: tuple[int, int, int] | None = None,
    xc_functional: str = "GGA",
    xc_authors: str = "PBE",
) -> None:
    """Write a minimal SIESTA input whose job is to emit .ion/.ORB_INDX/.DM."""
    if len(info.atom_labels) != scfout.atomnum:
        raise ValueError(
            f"atom count mismatch: input has {len(info.atom_labels)}, scfout has {scfout.atomnum}"
        )
    labels = [s.label for s in info.species]
    species_index = {lab: i + 1 for i, lab in enumerate(labels)}

    tv = scfout.tv[1:4, 1:4]
    xyz = scfout.gxyz[1:, 1:4]

    with open(path, "w") as f:
        f.write(f"SystemName          {label} template\n")
        f.write(f"SystemLabel         {label}\n\n")
        f.write(f"NumberOfAtoms       {scfout.atomnum}\n")
        f.write(f"NumberOfSpecies     {len(labels)}\n\n")
        f.write("%block ChemicalSpeciesLabel\n")
        for i, lab in enumerate(labels, 1):
            z = Z_OF_SYMBOL.get(lab)
            if z is None:
                raise ValueError(f"unknown element symbol {lab!r}")
            f.write(f"  {i:3d}  {z:3d}  {lab}\n")
        f.write("%endblock ChemicalSpeciesLabel\n\n")
        f.write("LatticeConstant 1.0 Bohr\n")
        f.write("%block LatticeVectors\n")
        for row in tv:
            f.write(f"  {row[0]:20.12f} {row[1]:20.12f} {row[2]:20.12f}\n")
        f.write("%endblock LatticeVectors\n\n")
        f.write("AtomicCoordinatesFormat Bohr\n")
        f.write("%block AtomicCoordinatesAndAtomicSpecies\n")
        for lab, pos in zip(info.atom_labels, xyz):
            f.write(
                f"  {pos[0]:20.12f} {pos[1]:20.12f} {pos[2]:20.12f}"
                f"  {species_index[lab]:3d}\n"
            )
        f.write("%endblock AtomicCoordinatesAndAtomicSpecies\n\n")
        f.write(f"XC.functional        {xc_functional}\n")
        f.write(f"XC.authors           {xc_authors}\n")
        f.write(f"MeshCutoff           {mesh_cutoff:g} Ry\n")
        f.write(f"PAO.BasisSize        {basis_size}\n")
        kg = kgrid or info.kgrid
        f.write("%block kgrid_Monkhorst_Pack\n")
        f.write(f"  {kg[0]:d} 0 0  0.0\n")
        f.write(f"  0 {kg[1]:d} 0  0.0\n")
        f.write(f"  0 0 {kg[2]:d}  0.0\n")
        f.write("%endblock kgrid_Monkhorst_Pack\n")
        f.write("SpinPolarized        false\n")
        f.write("MaxSCFIterations     1\n")
        f.write("SCFMustConverge      false\n")
        f.write("WriteKbands          false\n")


def run_auto_projection(
    scfout_path: str | Path,
    output: str | Path | None = None,
    workdir: str | Path | None = None,
    siesta_cmd: str = "siesta",
    siesta_pseudo_dir: str | Path | None = None,
    pao_dir: str | Path | None = None,
    basis_size: str = "SZ",
    mesh_cutoff: float = 150.0,
    spacing: float = 0.15,
    keep_template: bool = True,
    verbose: bool = False,
) -> tuple[ProjectionResult, Path, Path]:
    """
    End-to-end automated projection.

    Returns (projection_result, output_dm_path, template_dir).
    """
    scfout_path = Path(scfout_path).resolve()
    sc = read_scfout(scfout_path)
    info = parse_openmx_input_info(sc.input_text)
    label = re.sub(r"[^A-Za-z0-9_.-]+", "_", info.system_name or scfout_path.stem)

    if sc.spin_switch not in (0, 1):
        raise NotImplementedError("auto projection supports non-pol/collinear OpenMX only")
    if sc.spin_switch == 1:
        # The underlying projection sums spin channels and writes nspin=1.
        print("WARNING: OpenMX input is collinear spin-polarized; writing spin-summed nspin=1 DM")

    workdir = Path(workdir or f"{label}_c4s_auto").resolve()
    template_dir = workdir / "siesta_template"
    template_dir.mkdir(parents=True, exist_ok=True)

    # SIESTA 4.1.5 looks for <label>.vps or <label>.psf in the run directory.
    pseudo_src_dir = Path(siesta_pseudo_dir).resolve() if siesta_pseudo_dir else scfout_path.parent
    for sp in info.species:
        found = None
        for ext in (".psf", ".vps"):
            p = pseudo_src_dir / f"{sp.label}{ext}"
            if p.exists():
                found = p
                break
        if found is None:
            raise FileNotFoundError(
                f"SIESTA pseudo for {sp.label} not found in {pseudo_src_dir}; "
                f"put {sp.label}.psf there or pass --siesta-pseudo-dir"
            )
        _link_or_copy(found, template_dir / found.name)

    fdf_path = template_dir / f"{label}_template.fdf"
    # Use the OpenMX SCF k-grid by default. This is encoded in the embedded
    # input text in .scfout (scf.Kgrid) and avoids another user-facing knob.
    write_siesta_template_fdf(fdf_path, sc, info, label, basis_size, mesh_cutoff)

    out_path = template_dir / f"{label}_template.out"
    if not (template_dir / f"{label}.DM").exists():
        with open(fdf_path) as fdf_in, open(out_path, "w") as out:
            subprocess.run(
                shlex.split(siesta_cmd), cwd=template_dir, stdin=fdf_in, stdout=out,
                stderr=subprocess.STDOUT, check=True,
            )

    template_dm = template_dir / f"{label}.DM"
    orb_indx = template_dir / f"{label}.ORB_INDX"
    if not template_dm.exists() or not orb_indx.exists():
        raise FileNotFoundError(
            f"SIESTA template run did not create {template_dm.name}/{orb_indx.name}; "
            f"see {out_path}"
        )

    qtot = _read_siesta_total_electrons(out_path)

    ion_map = {sp.label: template_dir / f"{sp.label}.ion" for sp in info.species}
    missing_ions = [str(p) for p in ion_map.values() if not p.exists()]
    if missing_ions:
        raise FileNotFoundError("SIESTA template run did not create ion files: " + ", ".join(missing_ions))

    if pao_dir is None:
        if info.data_path is None:
            raise ValueError("OpenMX DATA.PATH not found in .scfout; pass --pao-dir")
        pao_root = info.data_path / "PAO"
    else:
        pao_root = Path(pao_dir)

    pao_map = {sp.label: pao_root / sp.pao_file for sp in info.species}
    missing_paos = [str(p) for p in pao_map.values() if not p.exists()]
    if missing_paos:
        raise FileNotFoundError("OpenMX PAO files not found: " + ", ".join(missing_paos))

    oi = read_orb_indx(orb_indx)
    template = read_dm(template_dm)
    radials_by_species = {lab: read_ion_radials(path) for lab, path in ion_map.items()}
    siesta_orbitals = build_siesta_orbitals(oi.unit, radials_by_species)
    omx_by_species = {
        sp.label: build_openmx_orbitals(parse_spec(sp.spec), read_pao_radials(pao_map[sp.label]))
        for sp in info.species
    }
    omx_per_atom = [omx_by_species[lab] for lab in info.atom_labels]

    res = project_openmx_to_siesta(
        sc, oi, template, siesta_orbitals, omx_per_atom,
        spacing=spacing, verbose=verbose, purify=True, qtot=qtot,
    )

    out_dm = Path(output or (workdir / f"{label}_omx_proj_pure.DM")).resolve()
    out_dm.parent.mkdir(parents=True, exist_ok=True)
    write_dm(out_dm, res.dm)

    return res, out_dm, template_dir


def _read_siesta_total_electrons(out_path: Path) -> float:
    """Read the target electron count from the SIESTA template output."""
    pat = re.compile(r"Total number of electrons:\s*([+-]?[0-9.]+(?:[Ee][+-]?\d+)?)")
    text = out_path.read_text(errors="ignore")
    matches = pat.findall(text)
    if not matches:
        raise ValueError(
            f"could not read target electron count from {out_path}; "
            "SIESTA did not print 'Total number of electrons'"
        )
    return float(matches[-1])


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="openmx2siesta-dm",
        description=(
            "Create a SIESTA template run and project an OpenMX .scfout density "
            "matrix into a SIESTA .DM file."
        ),
    )
    parser.add_argument("--scfout", required=True, help="OpenMX .scfout file")
    parser.add_argument("--siesta-pseudo-dir",
                        help="directory containing SIESTA Label.psf/Label.vps files; "
                        "default: the .scfout directory")
    parser.add_argument("--siesta-cmd", default="siesta",
                        help="SIESTA executable command (default: siesta)")
    parser.add_argument("--pao-dir",
                        help="OpenMX PAO directory override; default: DATA.PATH/PAO from .scfout")
    parser.add_argument("--workdir", help="working directory for generated template/output")
    parser.add_argument("--basis-size", default="SZ",
                        help="SIESTA PAO.BasisSize for the template run (default: SZ)")
    parser.add_argument("--mesh-cutoff", type=float, default=150.0,
                        help="SIESTA MeshCutoff in Ry for the template run (default: 150)")
    parser.add_argument("--spacing", type=float, default=0.15,
                        help="real-space grid spacing in Bohr for overlap quadrature")
    parser.add_argument("--output", help="output SIESTA .DM path")
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args(argv)

    res, out_dm, template_dir = run_auto_projection(
        scfout_path=args.scfout,
        output=args.output,
        workdir=args.workdir,
        siesta_cmd=args.siesta_cmd,
        siesta_pseudo_dir=args.siesta_pseudo_dir,
        pao_dir=args.pao_dir,
        basis_size=args.basis_size,
        mesh_cutoff=args.mesh_cutoff,
        spacing=args.spacing,
        verbose=args.verbose,
    )

    print(f"[openmx2siesta-dm] template dir: {template_dir}")
    print(f"[openmx2siesta-dm] wrote {out_dm}")
    print(f"    Tr(P_SIE S) = {res.n_electrons:.4f} electrons")
    if res.purified:
        print(f"    purified: natural occupations were in "
              f"[{res.occ_min:.4f}, {res.occ_max:.4f}], "
              f"{res.n_frac} fractional state(s)")
    print(f"    Hermiticity residual = {res.max_asym:.2e}  (k-points: {res.nk})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
