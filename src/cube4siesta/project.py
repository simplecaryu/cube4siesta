"""
Project an OpenMX density matrix onto a SIESTA NAO basis and write it as a
SIESTA .DM, using the overlap-projection method.

The density operator carried by OpenMX,

    rho_hat = sum_{mu,nu} sum_{cells} |phi_mu C> P_{mu,nu}(D) <phi_nu C+D| ,

is projected onto SIESTA's NAO basis. In a periodic crystal this is done
exactly on the Born–von-Karman k-grid implied by the .DM supercell size ``nsc``:

    S(k)_{ab}  = sum_R  <chi_a 0|chi_b R>      e^{i k.R}            (SIESTA overlap)
    U(k)_{a,mu}= sum_C  <chi_a 0|phi_mu C>     e^{i k.C}            (cross overlap)
    P(k)_{mu,nu}= sum_D P_{mu,nu}(D)           e^{i k.D}            (OpenMX DM)
    D(k)       = U(k) P(k) U(k)^dagger                              (rho in SIESTA basis)
    P_SIE(k)   = S(k)^{-1} D(k) S(k)^{-1}                           (contravariant DM)
    P_SIE(R)   = (1/Nk) sum_k P_SIE(k) e^{-i k.R}

This is exact in the complete-basis limit and independent of valence count: any
OpenMX density that SIESTA's basis cannot represent (e.g. semicore) simply
projects out, so Tr(P_SIE S) returns the representable electron count.

All overlap integrals are real-space grid quadratures (see overlap.py), with
each orbital evaluated from its own radial table + canonical real harmonic, so
the two codes' conventions are handled automatically.

Assumptions (validated by the electron-count / Hermiticity checks):
  * SIESTA and OpenMX describe the same crystal with the same lattice-vector
    ordering and atom ordering, so we use OpenMX geometry (gxyz, tv) and integer
    cell indices for both sides.
  * Collinear non-polarised input (spin_switch == 0) for now: the OpenMX DM[0]
    is one spin channel and the written SIESTA nspin=1 DM is 2x that.
"""

from __future__ import annotations

import itertools
from dataclasses import dataclass

import numpy as np

from .basis import Orbital
from .dm_io import DMFile, copy_pattern
from .orb_indx import OrbIndx
from .overlap import block_overlap_matrix
from .scfout_io import Scfout


def _cell_range(nsc: tuple[int, int, int]) -> list[tuple[int, int, int]]:
    """Centred integer cell triples spanning the BvK supercell of size nsc."""
    rngs = [range(-(n // 2), n - (n // 2)) for n in nsc]
    return list(itertools.product(*rngs))


def _atom_ranges(sizes: list[int]) -> list[tuple[int, int]]:
    """Contiguous orbital index ranges for atoms with given orbital counts."""
    out = []
    start = 0
    for s in sizes:
        out.append((start, start + s))
        start += s
    return out


@dataclass
class ProjectionResult:
    dm: DMFile
    n_electrons: float
    max_asym: float          # max |P_SIE(R)_ab - P_SIE(-R)_ba| (Hermiticity check)
    nk: int
    purified: bool = False
    occ_min: float = 0.0     # natural-occupation range BEFORE purification
    occ_max: float = 0.0
    n_frac: int = 0          # states left fractional at the Fermi boundary


def _sqrt_and_invsqrt(S: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Hermitian S^(1/2) and S^(-1/2) via eigendecomposition (eigenvalues clamped)."""
    s, V = np.linalg.eigh(S)
    s = np.clip(s, 1e-12, None)
    rs = np.sqrt(s)
    return (V * rs) @ V.conj().T, (V / rs) @ V.conj().T


def purify_canonical(
    Psie_k: np.ndarray,
    Sk: np.ndarray,
    n_channel: float,
    occ_cap: float,
    deg_tol: float = 1e-6,
) -> tuple[float, float, int]:
    """
    Spectral purification / occupation refill of a k-resolved density matrix,
    in place.

    Each P(k) is diagonalised in the S(k) metric to get natural occupations
    f_i(k) and S-orthonormal natural orbitals. The natural orbitals are kept
    fixed, while the occupations are refilled aufbau-style with a single global
    Fermi level across all k-points: the highest-occupation states are set to
    ``occ_cap`` until the average occupation per cell equals ``n_channel``.
    States degenerate with the Fermi boundary (within ``deg_tol``) share the
    remaining charge equally, so metals and odd electron counts are handled
    without breaking k-point symmetry.

    This is not an untargeted iterative McWeeny procedure. The target trace is
    explicit and determines the rank/filling of the final idempotent projector.
    Note: for spin-polarised sources the projection carries the spin-summed DM,
    so purification acts on total occupations capped at 2.

    Returns (occ_min, occ_max, n_frac): the natural-occupation range before
    purification and the number of fractional boundary states after.
    """
    nk, n, _ = Psie_k.shape
    occs = np.empty((nk, n))
    nat = np.empty_like(Psie_k)
    for ik in range(nk):
        S12, Sm12 = _sqrt_and_invsqrt(Sk[ik])
        Q = S12 @ Psie_k[ik] @ S12
        Q = 0.5 * (Q + Q.conj().T)
        f, W = np.linalg.eigh(Q)
        occs[ik] = f
        nat[ik] = Sm12 @ W
    occ_min, occ_max = float(occs.min()), float(occs.max())

    flat = occs.reshape(-1)
    order = np.argsort(flat)[::-1]
    fs = flat[order]
    total = n_channel * nk                      # required sum of occupations
    if total > occ_cap * len(fs) + 1e-9:
        raise ValueError(
            f"target {n_channel} electrons/channel exceeds basis capacity "
            f"{occ_cap * n}"
        )
    nfull = int(np.floor(total / occ_cap + 1e-12))
    rem = total - nfull * occ_cap
    new = np.zeros_like(flat)
    if total > 1e-12:
        f_b = fs[nfull] if rem > 1e-9 else fs[nfull - 1]
        group = np.abs(flat - f_b) < deg_tol
        new[(flat > f_b) & ~group] = occ_cap
        share = (total - new.sum()) / group.sum()
        if share < -1e-9 or share > occ_cap + 1e-9:
            raise ValueError("degenerate Fermi group cannot hold the remaining charge; "
                             "retry with a smaller deg_tol")
        new[group] = min(max(share, 0.0), occ_cap)
    n_frac = int(np.sum((new > 1e-9) & (new < occ_cap - 1e-9)))

    new = new.reshape(nk, n)
    for ik in range(nk):
        Psie_k[ik] = (nat[ik] * new[ik]) @ nat[ik].conj().T
    return occ_min, occ_max, n_frac


def project_openmx_to_siesta(
    scfout: Scfout,
    orb_indx: OrbIndx,
    template: DMFile,
    siesta_orbitals: list[Orbital],
    omx_orbitals_per_atom: list[list[Orbital]],
    spacing: float = 0.15,
    spin: int = 0,
    verbose: bool = False,
    purify: bool = False,
    qtot: float | None = None,
) -> ProjectionResult:
    """
    Build a SIESTA .DM (same sparsity as ``template``) from an OpenMX scfout.

    siesta_orbitals       : length no_u, aligned to .ORB_INDX unit orbital order
    omx_orbitals_per_atom : length atomnum (1-based atoms -> index ct-1), each the
                            OpenMX orbital list for that atom in DM order
    purify                : apply canonical purification (see purify_canonical),
                            fixing Tr(P S) to ``qtot`` (default: raw count rounded
                            to the nearest integer)
    """
    if scfout.spin_switch not in (0, 1):
        raise NotImplementedError("only spin_switch 0 (non-pol) or 1 (collinear) supported")

    no_u = template.no_u
    tv = scfout.tv[1:4, 1:4]                      # (3,3) lattice rows, Bohr
    gxyz = scfout.gxyz                            # 1-based
    atomnum = scfout.atomnum

    # ---- atom <-> orbital maps -------------------------------------------------
    s_atom_of_orb = orb_indx.atom_of_uorb()       # (no_u,) 1-based atom per SIESTA orb
    s_ranges: dict[int, tuple[int, int]] = {}
    for a in range(no_u):
        at = int(s_atom_of_orb[a])
        lo, hi = s_ranges.get(at, (a, a))
        s_ranges[at] = (min(lo, a), a + 1)

    omx_sizes = [len(omx_orbitals_per_atom[ct - 1]) for ct in range(1, atomnum + 1)]
    omx_ranges = _atom_ranges(omx_sizes)          # 0-based atom -> (lo,hi)
    no_omx = sum(omx_sizes)

    cells = _cell_range(template.nsc)
    cellset = set(cells)

    def cart(cell):
        c = np.asarray(cell, dtype=np.float64)
        return c[0] * tv[0] + c[1] * tv[1] + c[2] * tv[2]

    # ---- real-space SIESTA overlap blocks  S_{A,B}(R) -------------------------
    if verbose:
        print("building SIESTA overlap blocks ...", flush=True)
    S_blocks: dict[tuple[int, int, tuple[int, int, int]], np.ndarray] = {}
    for A in range(1, atomnum + 1):
        a0, a1 = s_ranges[A]
        orbsA = siesta_orbitals[a0:a1]
        posA = gxyz[A][1:4]
        for B in range(1, atomnum + 1):
            b0, b1 = s_ranges[B]
            orbsB = siesta_orbitals[b0:b1]
            for R in cells:
                posB = gxyz[B][1:4] + cart(R)
                blk = block_overlap_matrix(orbsA, posA, orbsB, posB, spacing=spacing)
                if np.any(blk != 0.0):
                    S_blocks[(A, B, R)] = blk

    # ---- real-space cross overlap blocks  U_{A,B}(C) -------------------------
    if verbose:
        print("building cross overlap blocks ...", flush=True)
    U_blocks: dict[tuple[int, int, tuple[int, int, int]], np.ndarray] = {}
    for A in range(1, atomnum + 1):
        a0, a1 = s_ranges[A]
        orbsA = siesta_orbitals[a0:a1]
        posA = gxyz[A][1:4]
        for B in range(1, atomnum + 1):
            orbsB = omx_orbitals_per_atom[B - 1]
            for C in cells:
                posB = gxyz[B][1:4] + cart(C)
                blk = block_overlap_matrix(orbsA, posA, orbsB, posB, spacing=spacing)
                if np.any(blk != 0.0):
                    U_blocks[(A, B, C)] = blk

    # ---- OpenMX real-space DM blocks  P_{ctA,ctB}(D) -------------------------
    # Sum spin channels to get the TOTAL density matrix. For spin_switch==0 the
    # single channel is one spin (the factor-2 below restores both); for
    # spin_switch==1 DM[0]+DM[1] is already the total (factor 1).
    spin_channels = [0] if scfout.spin_switch == 0 else [0, 1]
    P_blocks: dict[tuple[int, int, tuple[int, int, int]], np.ndarray] = {}
    for ct in range(1, atomnum + 1):
        for h in range(int(scfout.fnan[ct]) + 1):
            Gh = int(scfout.natn[ct][h])
            D = scfout.neighbor_cell_ijk(ct, h)
            if D not in cellset:
                continue                          # outside BvK range (negligible)
            blk = sum(scfout.DM[s][(ct, h)] for s in spin_channels)
            key = (ct, Gh, D)
            if key in P_blocks:
                P_blocks[key] = P_blocks[key] + blk
            else:
                P_blocks[key] = blk.copy()

    # ---- k-grid ---------------------------------------------------------------
    n0, n1, n2 = template.nsc
    kfracs = [
        (m0 / n0, m1 / n1, m2 / n2)
        for m0 in range(n0) for m1 in range(n1) for m2 in range(n2)
    ]
    nk = len(kfracs)
    two_pi = 2.0 * np.pi

    def phase(kf, cell):
        return np.exp(1j * two_pi * (kf[0] * cell[0] + kf[1] * cell[1] + kf[2] * cell[2]))

    # ---- per-k projection, accumulate P_SIE(k) -------------------------------
    if verbose:
        print(f"projecting over {nk} k-points ...", flush=True)
    Psie_k = np.empty((nk, no_u, no_u), dtype=np.complex128)
    Sk_all = np.zeros((nk, no_u, no_u), dtype=np.complex128)
    for ik, kf in enumerate(kfracs):
        Sk = Sk_all[ik]
        for (A, B, R), blk in S_blocks.items():
            a0, a1 = s_ranges[A]; b0, b1 = s_ranges[B]
            Sk[a0:a1, b0:b1] += blk * phase(kf, R)
        Uk = np.zeros((no_u, no_omx), dtype=np.complex128)
        for (A, B, C), blk in U_blocks.items():
            a0, a1 = s_ranges[A]; b0, b1 = omx_ranges[B - 1]
            Uk[a0:a1, b0:b1] += blk * phase(kf, C)
        Pk = np.zeros((no_omx, no_omx), dtype=np.complex128)
        for (ctA, ctB, D), blk in P_blocks.items():
            a0, a1 = omx_ranges[ctA - 1]; b0, b1 = omx_ranges[ctB - 1]
            Pk[a0:a1, b0:b1] += blk * phase(kf, D)

        Dk = Uk @ Pk @ Uk.conj().T
        Sinv = np.linalg.inv(Sk)
        Psie_k[ik] = Sinv @ Dk @ Sinv

    spin_factor = 2.0 if scfout.spin_switch == 0 else 1.0

    # ---- optional canonical purification --------------------------------------
    occ_min = occ_max = 0.0
    n_frac = 0
    if purify:
        ne_raw = spin_factor * np.real(
            np.einsum("kab,kba->", Psie_k, Sk_all)) / nk
        n_target = float(round(ne_raw)) if qtot is None else float(qtot)
        if verbose:
            print(f"purifying: raw Tr(PS) = {ne_raw:.4f} -> target {n_target:.4f}",
                  flush=True)
        occ_min, occ_max, n_frac = purify_canonical(
            Psie_k, Sk_all, n_target / spin_factor, 2.0 / spin_factor
        )

    # ---- inverse transform onto the template sparsity ------------------------
    out = copy_pattern(template, nspin=1)
    kf_arr = np.array(kfracs)                      # (nk,3)
    for io in range(no_u):
        listd, _ = template.row(io)
        p = int(template.listdptr[io])
        for jj, jo in enumerate(listd):
            juo, isc = orb_indx.decode_column(int(jo))
            b = juo - 1
            ph = np.exp(-1j * two_pi * (kf_arr @ np.array(isc)))   # (nk,)
            val = np.sum(Psie_k[:, io, b] * ph) / nk
            out.dm[p + jj, 0] = spin_factor * val.real

    # ---- diagnostics: electron count and Hermiticity -------------------------
    ne = spin_factor * np.real(np.einsum("kab,kba->", Psie_k, Sk_all)) / nk

    max_asym = _hermiticity_residual(out, orb_indx)

    return ProjectionResult(
        dm=out, n_electrons=float(ne), max_asym=max_asym, nk=nk,
        purified=purify, occ_min=occ_min, occ_max=occ_max, n_frac=n_frac,
    )


def _hermiticity_residual(dm: DMFile, orb_indx: OrbIndx) -> float:
    """max |P(R)_ab - P(-R)_ba| over stored entries (should be ~0)."""
    # index stored entries by (a, b, isc)
    store: dict[tuple[int, int, tuple[int, int, int]], float] = {}
    for io in range(dm.no_u):
        listd, vals = dm.row(io)
        for jo, v in zip(listd, vals[:, 0]):
            juo, isc = orb_indx.decode_column(int(jo))
            store[(io, juo - 1, isc)] = float(v)
    worst = 0.0
    for (a, b, isc), v in store.items():
        mirror = store.get((b, a, (-isc[0], -isc[1], -isc[2])))
        if mirror is not None:
            worst = max(worst, abs(v - mirror))
    return worst
