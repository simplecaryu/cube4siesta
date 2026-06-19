"""
Radial-function loaders and real spherical harmonics for SIESTA NAOs and
OpenMX PAOs, plus an :class:`Orbital` that can be evaluated on a real-space grid.

Both codes write a localised orbital as  phi(r) = R(r) * Y(dir)  with R(r)
tabulated and Y a *real* spherical harmonic. We define ONE canonical set of
normalised real harmonics (matching OpenMX's AngularF.c and SIESTA's labels)
and attach to each orbital the harmonic implied by its code's convention, so a
cross-overlap integral naturally connects matching harmonics.

Conventions
-----------
Canonical real harmonics (normalised, ∫|Y|²dΩ = 1), keyed by name:
    s; px, py, pz; dz2, dx2y2, dxy, dxz, dyz; (f left unimplemented until needed)

SIESTA: harmonic taken from the .ORB_INDX ``sym`` column
    s | py(m=-1) pz(m=0) px(m=+1) | dxy dyz dz2 dxz dx2-y2
OpenMX (openmx3.9/source/AngularF.c), orbital index within an l-shell:
    l=1: 0->px 1->py 2->pz
    l=2: 0->dz2 1->dx2y2 2->dxy 3->dxz 4->dyz
Within an atom both codes order orbitals as l (outer) -> multiplicity/zeta ->
m (inner); e.g. OpenMX ``s2p2`` => s s' px py pz px' py' pz'.

Radial tables
-------------
SIESTA .ion: per orbital a header ``l n z is_pol pop``, then ``npts delta rc``,
then ``npts`` rows of ``r value`` on a *linear* grid; value is R(r).
OpenMX .pao: ``<pseudo.atomic.orbitals.L=k>`` blocks with rows
``x  r  R_0 R_1 ... R_{Mul-1}`` (r = exp(x), logarithmic grid); column j+2 is
the j-th radial function for that l. ``s2p2d1`` selects the first 2/2/1 of them
for l = 0/1/2.
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from scipy.interpolate import CubicSpline


# --------------------------------------------------------------------------
# Canonical normalised real spherical harmonics, as functions of a unit vector
# (nx, ny, nz). Coefficients are the standard ones (verified against AngularF.c).
# --------------------------------------------------------------------------
_C_S = 0.5 / math.sqrt(math.pi)
_C_P = math.sqrt(3.0 / (4.0 * math.pi))
_C_DZ2 = math.sqrt(5.0 / (16.0 * math.pi))
_C_DX2 = math.sqrt(15.0 / (16.0 * math.pi))
_C_DXY = math.sqrt(15.0 / (4.0 * math.pi))
# f-harmonic coefficients (from OpenMX AngularF.c, l=3)
_C_F0 = 0.373176332590116    # 5z^3-3z
_C_F1 = 0.457045799464466    # x(5z^2-1), y(5z^2-1)
_C_F3 = 1.44530572132028     # z(x^2-y^2)
_C_F4 = 2.89061144264055     # xyz
_C_F5 = 0.590043589926644    # x(x^2-3y^2), y(3x^2-y^2)


def _angular(name: str, nx: np.ndarray, ny: np.ndarray, nz: np.ndarray) -> np.ndarray:
    if name == "s":
        return np.full_like(nx, _C_S)
    if name == "px":
        return _C_P * nx
    if name == "py":
        return _C_P * ny
    if name == "pz":
        return _C_P * nz
    if name == "dz2":
        return _C_DZ2 * (3.0 * nz * nz - 1.0)
    if name == "dx2y2":
        return _C_DX2 * (nx * nx - ny * ny)
    if name == "dxy":
        return _C_DXY * nx * ny
    if name == "dxz":
        return _C_DXY * nx * nz
    if name == "dyz":
        return _C_DXY * ny * nz
    # f harmonics (OpenMX index order names)
    if name == "fz3":
        return _C_F0 * (5.0 * nz * nz * nz - 3.0 * nz)
    if name == "fxz2":
        return _C_F1 * nx * (5.0 * nz * nz - 1.0)
    if name == "fyz2":
        return _C_F1 * ny * (5.0 * nz * nz - 1.0)
    if name == "fzx2y2":
        return _C_F3 * nz * (nx * nx - ny * ny)
    if name == "fxyz":
        return _C_F4 * nx * ny * nz
    if name == "fx3":
        return _C_F5 * (nx * nx * nx - 3.0 * nx * ny * ny)
    if name == "fy3":
        return _C_F5 * (3.0 * nx * nx * ny - ny * ny * ny)
    raise NotImplementedError(f"real harmonic '{name}' not implemented")


_L_OF_NAME = {
    "s": 0,
    "px": 1, "py": 1, "pz": 1,
    "dz2": 2, "dx2y2": 2, "dxy": 2, "dxz": 2, "dyz": 2,
    "fz3": 3, "fxz2": 3, "fyz2": 3, "fzx2y2": 3, "fxyz": 3, "fx3": 3, "fy3": 3,
}

# OpenMX index-within-l -> canonical harmonic name (from AngularF.c)
OPENMX_NAMES = {
    0: ["s"],
    1: ["px", "py", "pz"],
    2: ["dz2", "dx2y2", "dxy", "dxz", "dyz"],
    3: ["fz3", "fxz2", "fyz2", "fzx2y2", "fxyz", "fx3", "fy3"],
}


def siesta_sym_to_name(sym: str) -> str:
    """Map a SIESTA .ORB_INDX ``sym`` label to a canonical harmonic name."""
    s = sym.lstrip("P")              # strip polarisation-orbital prefix
    s = s.replace("d", "d", 1)       # keep as-is
    table = {
        "s": "s",
        "px": "px", "py": "py", "pz": "pz",
        "dz2": "dz2", "dx2-y2": "dx2y2", "dxy": "dxy", "dxz": "dxz", "dyz": "dyz",
    }
    if s not in table:
        raise NotImplementedError(f"SIESTA sym '{sym}' (->'{s}') not mapped")
    return table[s]


@dataclass
class Orbital:
    """
    A localised orbital phi(r) = R(r) * Y_name(dir), where R(r) is the *full*
    radial part (so for l>0 it already vanishes ~r^l at the origin) and Y_name
    is a canonical normalised real harmonic. The stored functions come from the
    code's own (normalised) basis files, so no extra normalisation is applied.
    """
    name: str                 # canonical harmonic name
    l: int
    zeta: int                 # multiplicity index (1-based)
    rcut: float
    r_min: float              # smallest tabulated radius (spline lower bound)
    _spline: CubicSpline

    def radial(self, r: np.ndarray) -> np.ndarray:
        out = np.zeros_like(r, dtype=np.float64)
        inside = r < self.rcut
        if np.any(inside):
            rc = np.clip(r[inside], self.r_min, self.rcut)
            out[inside] = self._spline(rc)
        return out

    def evaluate(self, dx: np.ndarray, dy: np.ndarray, dz: np.ndarray) -> np.ndarray:
        """Evaluate the orbital at displacements (dx,dy,dz) from its centre."""
        r = np.sqrt(dx * dx + dy * dy + dz * dz)
        val = np.zeros_like(r, dtype=np.float64)
        inside = r < self.rcut
        if not np.any(inside):
            return val
        ri = r[inside]
        rc = np.clip(ri, self.r_min, self.rcut)
        rad = self._spline(rc)
        if self.l == 0:
            ang = np.full_like(ri, _C_S)
        else:
            safe = ri > 1e-12
            nx = np.zeros_like(ri); ny = np.zeros_like(ri); nz = np.zeros_like(ri)
            nx[safe] = dx[inside][safe] / ri[safe]
            ny[safe] = dy[inside][safe] / ri[safe]
            nz[safe] = dz[inside][safe] / ri[safe]
            ang = _angular(self.name, nx, ny, nz)
            ang[~safe] = 0.0     # l>0 orbitals vanish at the centre
        val[inside] = rad * ang
        return val


def _spline_from(r: np.ndarray, R: np.ndarray) -> CubicSpline:
    # CubicSpline requires strictly increasing x; r grids here already are.
    return CubicSpline(r, R, extrapolate=False)


# --------------------------------------------------------------------------
# SIESTA .ion radial tables
# --------------------------------------------------------------------------
def read_ion_radials(path: str | Path) -> dict[tuple[int, ...], tuple[np.ndarray, np.ndarray, float]]:
    """
    Parse the ``# PAOs`` section of a SIESTA .ion file.

    Returns radial tables keyed primarily by ``(n, l, zeta)``. For ordinary
    bases where ``(l, zeta)`` is unique, a backward-compatible ``(l, zeta)``
    alias is also stored. The principal quantum number is needed for semicore
    PAO blocks such as Mo 4s and 5s, which can share the same ``l`` and ``z``.
    """
    radials: dict[tuple[int, ...], tuple[np.ndarray, np.ndarray, float]] = {}
    seen_lz: set[tuple[int, int]] = set()
    duplicate_lz: set[tuple[int, int]] = set()
    with open(path) as fh:
        lines = fh.readlines()
    i = 0
    while i < len(lines):
        if lines[i].lstrip().startswith("# PAOs"):
            i += 1
            break
        i += 1
    # multiplicity counter per l
    zeta_count: dict[int, int] = {}
    while i < len(lines):
        line = lines[i].strip()
        if not line or line.startswith("#"):
            i += 1
            continue
        head = line.split()
        # orbital header: "l n z is_pol pop"
        if len(head) >= 5 and all(_is_int(head[k]) for k in range(4)):
            l = int(head[0]); n = int(head[1]); z = int(head[2])
            i += 1
            npts, _delta, rc = (lambda t: (int(t[0]), float(t[1]), float(t[2])))(
                lines[i].split()
            )
            i += 1
            rr = np.empty(npts); RR = np.empty(npts)
            for k in range(npts):
                a, b = lines[i].split()[:2]
                rr[k] = float(a); RR[k] = float(b)
                i += 1
            value = (rr, RR, rc)
            radials[(n, l, z)] = value
            lz = (l, z)
            if lz in seen_lz:
                duplicate_lz.add(lz)
                radials.pop(lz, None)
            elif lz not in duplicate_lz:
                radials[lz] = value
            seen_lz.add(lz)
        else:
            # reached the end of the PAO section (e.g. "# KBs")
            break
    return radials


def read_ion_populations(path: str | Path) -> dict[tuple[int, int, int], float]:
    """
    Return the neutral-atom shell populations from a SIESTA .ion ``# PAOs``
    section, keyed by ``(n, l, zeta)``.

    The population printed on each PAO header line is the occupancy of the whole
    ``(n, l, z)`` shell. Dividing by ``2l + 1`` gives the per-``m`` occupancy used
    to seed a diagonal neutral-atom density matrix.
    """
    pops: dict[tuple[int, int, int], float] = {}
    with open(path) as fh:
        lines = fh.readlines()
    i = 0
    while i < len(lines) and not lines[i].lstrip().startswith("# PAOs"):
        i += 1
    for line in lines[i + 1:]:
        t = line.split()
        if len(t) >= 5 and all(_is_int(t[k]) for k in range(4)):
            l, n, z = int(t[0]), int(t[1]), int(t[2])
            try:
                pops[(n, l, z)] = float(t[4])
            except ValueError:
                pass
        elif t and t[0].startswith("#"):
            # reached the end of the PAO section (e.g. "# KBs")
            break
    return pops


def build_siesta_orbitals(unit_orbitals, radials_by_species) -> list[Orbital]:
    """
    Build the unit-cell Orbital list (aligned to .ORB_INDX unit orbital order).

    unit_orbitals : list of orb_indx.Orbital (the no_u unit-cell entries)
    radials_by_species : {species_label: read_ion_radials(...) dict}
    """
    out: list[Orbital] = []
    for o in unit_orbitals:
        table = radials_by_species[o.spec]
        key_nlz = (o.n, o.l, o.z)
        key_lz = (o.l, o.z)
        if key_nlz in table:
            rr, RR, rc = table[key_nlz]
        else:
            rr, RR, rc = table[key_lz]
        # SIESTA .ion stores f(r) with orbital = f(r) * r^l * Y_lm (see
        # atmfuncs.f: phi = phir * r**l * rly). Bake r^l into the radial so
        # R(r) is the full radial part, matching OpenMX's convention.
        #
        # SIESTA's rlylm uses real harmonics C*P_l^m*{cos,sin}(m phi) with the
        # Condon-Shortley (-1)^m phase folded into P_l^m (verified empirically:
        # px, py flip sign vs s; confirmed for d). Our canonical harmonics omit
        # CS (they match OpenMX's AngularF.c), so multiply SIESTA orbitals by
        # (-1)^m to reproduce SIESTA's actual orbital signs.
        cs_sign = (-1.0) ** o.m
        R_full = cs_sign * RR * (rr ** o.l)
        orb = Orbital(
            name=siesta_sym_to_name(o.sym),
            l=o.l,
            zeta=o.z,
            rcut=rc,
            r_min=float(rr[0]),
            _spline=_spline_from(rr, R_full),
        )
        out.append(orb)
    return out


# --------------------------------------------------------------------------
# OpenMX .pao radial tables
# --------------------------------------------------------------------------
def parse_spec(spec: str) -> dict[int, int]:
    """'s2p2d1' -> {0:2, 1:2, 2:1}."""
    lmap = {"s": 0, "p": 1, "d": 2, "f": 3}
    out: dict[int, int] = {}
    for letter, num in re.findall(r"([spdf])(\d+)", spec):
        out[lmap[letter]] = int(num)
    if not out:
        raise ValueError(f"could not parse basis spec '{spec}'")
    return out


def read_pao_radials(path: str | Path) -> dict[int, tuple[np.ndarray, np.ndarray]]:
    """
    Parse ``<pseudo.atomic.orbitals.L=k>`` blocks of an OpenMX .pao file.

    Returns {l: (r, R_all)} where R_all has shape (npts, Mul); column j is the
    j-th radial function for angular momentum l. Grid r = column 2 (= exp(x)).
    """
    radials: dict[int, tuple[np.ndarray, np.ndarray]] = {}
    with open(path) as fh:
        lines = fh.readlines()
    i = 0
    while i < len(lines):
        m = re.match(r"<pseudo\.atomic\.orbitals\.L=(\d+)", lines[i].strip())
        if not m:
            i += 1
            continue
        l = int(m.group(1))
        i += 1
        rows = []
        while i < len(lines) and not lines[i].strip().startswith("pseudo.atomic.orbitals.L="):
            tok = lines[i].split()
            if tok:
                rows.append([float(x) for x in tok])
            i += 1
        arr = np.array(rows)          # (npts, 2 + Mul)
        r = arr[:, 1].copy()          # column 2 = r
        R_all = arr[:, 2:].copy()     # the Mul radial functions
        radials[l] = (r, R_all)
    return radials


def build_openmx_orbitals(
    spec: dict[int, int],
    pao_radials: dict[int, tuple[np.ndarray, np.ndarray]],
) -> list[Orbital]:
    """
    Build one atom's OpenMX Orbital list in DM index order: l -> zeta -> m.

    spec : {l: multiplicity}, e.g. {0:2, 1:2, 2:1}
    pao_radials : read_pao_radials(...) output
    """
    out: list[Orbital] = []
    for l in sorted(spec):
        r, R_all = pao_radials[l]
        rcut = float(r[-1])
        r_min = float(r[0])
        names = OPENMX_NAMES[l]
        for mul in range(spec[l]):
            spline = _spline_from(r, R_all[:, mul])
            for name in names:
                out.append(
                    Orbital(
                        name=name, l=l, zeta=mul + 1, rcut=rcut,
                        r_min=r_min, _spline=spline,
                    )
                )
    return out


def _is_int(tok: str) -> bool:
    try:
        int(tok)
        return True
    except ValueError:
        return False
