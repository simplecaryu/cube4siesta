"""
Parse source-code pseudopotential files (OpenMX .vps, QE .UPF) and
generate ATOM input (.inp) files for producing a matching SIESTA .psf.

The goal: when restarting SIESTA from another code's density, the
SIESTA pseudopotential must have the same valence partition. This module
automates extracting the key parameters (element, XC, valence config,
cutoff radii) and writing a ready-to-run ATOM input.

Usage:
    cube4siesta gen-atom-input --from-vps V_PBE19.vps --output V.pg.inp
    cube4siesta gen-atom-input --from-upf Si.pbe-n-van.UPF --output Si.pg.inp
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

_ELEMENTS = [
    "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
    "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar",
    "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Ga", "Ge", "As", "Se", "Br", "Kr",
    "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
    "In", "Sn", "Sb", "Te", "I", "Xe",
    "Cs", "Ba",
    "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu",
    "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "Tl", "Pb", "Bi", "Po", "At", "Rn",
    "Fr", "Ra",
    "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf",
]

_L_LABELS = "spdf"


def _element_symbol(z: int) -> str:
    return _ELEMENTS[z - 1]


@dataclass
class ValenceShell:
    n: int
    l: int
    occ: float
    rc: float


@dataclass
class PseudoParams:
    z: int
    symbol: str
    xc: str          # ATOM code: "ca", "pb", "wi", etc.
    relativistic: bool
    n_core_shells: int
    valence: list[ValenceShell] = field(default_factory=list)

    @property
    def n_valence_electrons(self) -> float:
        return sum(v.occ for v in self.valence)


# ---------------------------------------------------------------------------
# OpenMX .vps parser
# ---------------------------------------------------------------------------

def parse_openmx_vps(path: str | Path) -> PseudoParams:
    text = Path(path).read_text()

    z = int(re.search(r"AtomSpecies\s+(\d+)", text).group(1))
    symbol = _element_symbol(z)

    xc_raw = re.search(r"xc\.type\s+(\w+)", text).group(1).upper()
    xc = {"GGA": "pb", "LDA": "ca"}.get(xc_raw, "pb")

    eq_raw = re.search(r"eq\.type\s+(\w+)", text).group(1).lower()
    relativistic = eq_raw in ("dirac", "sdirac")

    # Parse occupied electrons to determine core vs valence
    occ_block = re.search(
        r"<occupied\.electrons\s*(.*?)occupied\.electrons>",
        text, re.DOTALL,
    ).group(1).strip().splitlines()

    full_config: list[tuple[int, int, float]] = []
    for line in occ_block:
        parts = line.split()
        n = int(parts[0])
        for l_idx, occ_str in enumerate(parts[1:]):
            occ = float(occ_str)
            full_config.append((n, l_idx, occ))

    # Parse pseudo.NandL for valence channels and cutoff radii
    pnl_block = re.search(
        r"<pseudo\.NandL\s*(.*?)pseudo\.NandL>",
        text, re.DOTALL,
    ).group(1).strip().splitlines()

    val_channels: list[tuple[int, int, float]] = []
    for line in pnl_block:
        parts = line.split()
        # format: index  n  l  rc  flag
        n, l, rc = int(parts[1]), int(parts[2]), float(parts[3])
        val_channels.append((n, l, rc))

    # Determine which (n,l) are valence (appear in pseudo.NandL)
    val_nl = {(n, l) for n, l, _ in val_channels}
    rc_map = {(n, l): rc for n, l, rc in val_channels}

    valence_shells: list[ValenceShell] = []
    core_count = 0
    for n, l, occ in full_config:
        if (n, l) in val_nl:
            rc = rc_map.get((n, l), 2.0)
            valence_shells.append(ValenceShell(n=n, l=l, occ=occ, rc=rc))
        elif occ > 0:
            core_count += 1

    return PseudoParams(
        z=z, symbol=symbol, xc=xc,
        relativistic=relativistic,
        n_core_shells=core_count,
        valence=valence_shells,
    )


# ---------------------------------------------------------------------------
# QE .UPF parser (v1 and v2)
# ---------------------------------------------------------------------------

def parse_qe_upf(path: str | Path) -> PseudoParams:
    text = Path(path).read_text()

    # Element
    m = re.search(r"Element:\s*(\w+)", text) or re.search(
        r"^\s*(\w{1,2})\s+Element", text, re.MULTILINE,
    )
    symbol = m.group(1).strip().capitalize()
    z = _ELEMENTS.index(symbol) + 1

    # XC functional
    xc = "pb"  # default PBE
    if re.search(r"PBE|GGA", text):
        xc = "pb"
    elif re.search(r"PZ|LDA|SLA\s+PZ", text):
        xc = "ca"

    # Relativistic
    relativistic = bool(re.search(r"[Ss]calar.?[Rr]elativistic|[Ff]ull.?[Rr]elativistic", text))

    # Z valence
    m_zv = re.search(r"([\d.]+)\s+Z valence", text)
    z_val = float(m_zv.group(1)) if m_zv else 0.0

    # Extract valence shells from PP_INFO.
    # Two UPF variants:
    #   v2+: "Valence configuration:\n nl pn l occ Rcut ..."  (prefer this block)
    #   v1:  "nl pn l occ Rcut Rcut_US E_pseu"  directly in PP_INFO
    # Either way, lines look like "3S  3  0  2.00  <Rcut>  <Rcut_US>  <E>"

    val_section = re.search(
        r"Valence configuration:\s*\n\s*nl.*?\n(.*?)(?:Generation|<|$)",
        text, re.DOTALL,
    )
    if val_section:
        val_text = val_section.group(1)
    else:
        pp_info = re.search(r"<PP_INFO>(.*?)</PP_INFO>", text, re.DOTALL)
        val_text = pp_info.group(1) if pp_info else text

    # Pattern: "nL  pn  l  occ  Rcut  [Rcut_US  [E_pseu]]"
    shell_pattern = re.compile(
        r"(\d)([SPDF])\s+\d+\s+(\d+)\s+([\d.]+)\s+([\d.E+\-]+)\s+([\d.E+\-]+)"
    )
    valence_shells: list[ValenceShell] = []
    seen: set[tuple[int, int]] = set()
    for m in shell_pattern.finditer(val_text):
        n = int(m.group(1))
        l = int(m.group(3))
        occ = float(m.group(4))
        rcut1 = float(m.group(5))
        rcut2 = float(m.group(6))
        # Use the smaller of the two Rcut columns as the pseudization radius.
        # In v1 UPFs, Rcut can be huge (grid extent) while Rcut_US is the real rc.
        rc = min(rcut1, rcut2)
        key = (n, l)
        if key in seen:
            continue
        seen.add(key)
        valence_shells.append(ValenceShell(n=n, l=l, occ=occ, rc=rc))

    # Determine core shells: total Z minus valence electrons
    core_electrons = z - z_val
    core_count = 0
    e_left = core_electrons
    for n in range(1, 8):
        for l in range(n):
            if e_left <= 0:
                break
            max_occ = 2 * (2 * l + 1)
            if e_left >= max_occ:
                core_count += 1
                e_left -= max_occ

    return PseudoParams(
        z=z, symbol=symbol, xc=xc,
        relativistic=relativistic,
        n_core_shells=core_count,
        valence=valence_shells,
    )


# ---------------------------------------------------------------------------
# ATOM .inp writer
# ---------------------------------------------------------------------------

def write_atom_input(params: PseudoParams, path: str | Path) -> str:
    """Write an ATOM pg-mode input file. Returns the content as a string."""
    xc_code = params.xc
    if params.relativistic:
        xc_code += "r"

    n_val = len(params.valence)

    lines = [
        "pg",
        f"tm2    3.0",
        f"{params.symbol:<7s}{xc_code}",
        "0.0    0.0    0.0    0.0    0.0    0.0",
        f"   {params.n_core_shells}    {n_val}",
    ]
    for v in params.valence:
        lines.append(f"   {v.n}    {v.l}    {v.occ:5.2f}   0.00")

    # Cutoff radii line: gather by l (ATOM wants s, p, d, f order)
    rc_by_l: dict[int, float] = {}
    for v in params.valence:
        if v.l not in rc_by_l:
            rc_by_l[v.l] = v.rc

    rc_list = [rc_by_l.get(l, 0.0) for l in range(4)]
    # Columns 5,6 = core correction flag and radius (0 = off)
    lines.append(
        f"  {rc_list[0]:.2f}   {rc_list[1]:.2f}   {rc_list[2]:.2f}   "
        f"{rc_list[3]:.2f}   0.00   0.00"
    )

    content = "\n".join(lines) + "\n"
    Path(path).write_text(content)
    return content
