"""
Two-centre overlap integrals between localised orbitals, by real-space grid
quadrature.

    <a|b> = ∫ phi_a(r - Ra) phi_b(r - Rb) d^3r

For orbitals on the *same* centre we use an exact 1-D radial quadrature
(the angular parts are orthonormal, so the integral is delta_{name} ∫R_aR_b r²dr).
For different centres we integrate on a regular 3-D grid covering the lens-shaped
overlap region of the two cutoff spheres. The grid spacing controls accuracy;
0.10–0.15 Bohr reproduces code-computed overlaps to ~1e-3.

These primitives are convention-agnostic: each orbital is evaluated from its own
radial table and canonical real harmonic, so signs / m-orderings of the two
codes are handled correctly as long as each Orbital was built correctly. The
recommended validation is to rebuild a code's own overlap matrix and compare to
the overlap stored by that code (see tests / scripts).
"""

from __future__ import annotations

import numpy as np

from .basis import Orbital

# numpy>=2.0 renamed trapz -> trapezoid
_trapz = getattr(np, "trapezoid", None) or np.trapz


def radial_overlap(orb_a: Orbital, orb_b: Orbital, npts: int = 4000) -> float:
    """Exact same-centre overlap: 0 unless harmonics match, else ∫R_aR_b r² dr."""
    if orb_a.name != orb_b.name:
        return 0.0
    rmax = min(orb_a.rcut, orb_b.rcut)
    rmin = max(orb_a.r_min, orb_b.r_min)
    if rmax <= rmin:
        return 0.0
    r = np.linspace(rmin, rmax, npts)
    Ra = orb_a.radial(r)
    Rb = orb_b.radial(r)
    return float(_trapz(Ra * Rb * r * r, r))


def two_center_overlap(
    orb_a: Orbital,
    ctr_a: np.ndarray,
    orb_b: Orbital,
    ctr_b: np.ndarray,
    spacing: float = 0.12,
    same_center_tol: float = 1e-6,
) -> float:
    """
    Overlap <phi_a(.-ctr_a) | phi_b(.-ctr_b)>.

    ctr_a, ctr_b : (3,) Cartesian centres (Bohr).
    spacing      : grid spacing (Bohr) for the different-centre case.
    """
    ctr_a = np.asarray(ctr_a, dtype=np.float64)
    ctr_b = np.asarray(ctr_b, dtype=np.float64)
    d = float(np.linalg.norm(ctr_a - ctr_b))

    if d < same_center_tol:
        return radial_overlap(orb_a, orb_b)

    if d >= orb_a.rcut + orb_b.rcut:
        return 0.0

    # Axis-aligned bounding box of the intersection of the two cutoff spheres.
    lo = np.maximum(ctr_a - orb_a.rcut, ctr_b - orb_b.rcut)
    hi = np.minimum(ctr_a + orb_a.rcut, ctr_b + orb_b.rcut)
    if np.any(hi <= lo):
        return 0.0

    # Grid points at cell centres (midpoint rule), spacing ~ `spacing`.
    ns = np.maximum(np.ceil((hi - lo) / spacing).astype(int), 1)
    xs = lo[0] + (np.arange(ns[0]) + 0.5) * (hi[0] - lo[0]) / ns[0]
    ys = lo[1] + (np.arange(ns[1]) + 0.5) * (hi[1] - lo[1]) / ns[1]
    zs = lo[2] + (np.arange(ns[2]) + 0.5) * (hi[2] - lo[2]) / ns[2]
    dV = np.prod((hi - lo) / ns)

    X, Y, Z = np.meshgrid(xs, ys, zs, indexing="ij")
    Xf, Yf, Zf = X.ravel(), Y.ravel(), Z.ravel()

    va = orb_a.evaluate(Xf - ctr_a[0], Yf - ctr_a[1], Zf - ctr_a[2])
    nz = va != 0.0
    if not np.any(nz):
        return 0.0
    vb = orb_b.evaluate(Xf[nz] - ctr_b[0], Yf[nz] - ctr_b[1], Zf[nz] - ctr_b[2])
    return float(np.sum(va[nz] * vb) * dV)


def onsite_overlap_matrix(orbitals: list[Orbital]) -> np.ndarray:
    """Same-centre overlap matrix for a list of orbitals on one atom (exact)."""
    n = len(orbitals)
    S = np.zeros((n, n))
    for i in range(n):
        for j in range(i, n):
            S[i, j] = S[j, i] = radial_overlap(orbitals[i], orbitals[j])
    return S


def block_overlap_matrix(
    orbs_a: list[Orbital],
    ctr_a: np.ndarray,
    orbs_b: list[Orbital],
    ctr_b: np.ndarray,
    spacing: float = 0.12,
    same_center_tol: float = 1e-6,
) -> np.ndarray:
    """
    Full overlap block <a_i | b_j> for all orbitals of two atoms at ``ctr_a``,
    ``ctr_b``. The grid is built ONCE for the atom pair and every orbital is
    evaluated once on it; the block is then a single matrix product. This is far
    cheaper than calling :func:`two_center_overlap` per (i, j).

    Returns an (len(orbs_a), len(orbs_b)) array.
    """
    ctr_a = np.asarray(ctr_a, dtype=np.float64)
    ctr_b = np.asarray(ctr_b, dtype=np.float64)
    na, nb = len(orbs_a), len(orbs_b)
    out = np.zeros((na, nb))
    d = float(np.linalg.norm(ctr_a - ctr_b))

    if d < same_center_tol:
        for i in range(na):
            for j in range(nb):
                out[i, j] = radial_overlap(orbs_a[i], orbs_b[j])
        return out

    rca = max(o.rcut for o in orbs_a)
    rcb = max(o.rcut for o in orbs_b)
    if d >= rca + rcb:
        return out

    lo = np.maximum(ctr_a - rca, ctr_b - rcb)
    hi = np.minimum(ctr_a + rca, ctr_b + rcb)
    if np.any(hi <= lo):
        return out

    ns = np.maximum(np.ceil((hi - lo) / spacing).astype(int), 1)
    xs = lo[0] + (np.arange(ns[0]) + 0.5) * (hi[0] - lo[0]) / ns[0]
    ys = lo[1] + (np.arange(ns[1]) + 0.5) * (hi[1] - lo[1]) / ns[1]
    zs = lo[2] + (np.arange(ns[2]) + 0.5) * (hi[2] - lo[2]) / ns[2]
    dV = float(np.prod((hi - lo) / ns))
    X, Y, Z = np.meshgrid(xs, ys, zs, indexing="ij")
    Xf, Yf, Zf = X.ravel(), Y.ravel(), Z.ravel()

    # Evaluate every orbital of each atom once on the shared grid.
    A = np.empty((na, Xf.size))
    for i, o in enumerate(orbs_a):
        A[i] = o.evaluate(Xf - ctr_a[0], Yf - ctr_a[1], Zf - ctr_a[2])
    # Restrict to points where some A orbital is nonzero (inside A's sphere).
    mask = np.any(A != 0.0, axis=0)
    if not np.any(mask):
        return out
    A = A[:, mask]
    Xm, Ym, Zm = Xf[mask], Yf[mask], Zf[mask]
    B = np.empty((nb, Xm.size))
    for j, o in enumerate(orbs_b):
        B[j] = o.evaluate(Xm - ctr_b[0], Ym - ctr_b[1], Zm - ctr_b[2])

    return (A @ B.T) * dV
