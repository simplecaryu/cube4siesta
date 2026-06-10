"""
Tests for the OpenMX -> SIESTA density-matrix projection pipeline.

The convention-critical checks (overlap engine reproduces each code's own OLP,
projection conserves electrons) run against the bulk-Si fixtures in
testdata/si_projection (see its README for provenance).
"""

from __future__ import annotations

import tempfile
from pathlib import Path

import numpy as np
import pytest

from cube4siesta.basis import (OPENMX_NAMES, _angular, parse_spec,
                               siesta_sym_to_name)
from cube4siesta.dm_io import copy_pattern, read_dm, write_dm

TD = Path(__file__).resolve().parent.parent / "testdata/si_projection"


# ----------------------------- unconditional --------------------------------
def test_parse_spec():
    assert parse_spec("s2p2d1") == {0: 2, 1: 2, 2: 1}
    assert parse_spec("s3p2d2f1") == {0: 3, 1: 2, 2: 2, 3: 1}


def test_siesta_sym_mapping():
    assert siesta_sym_to_name("s") == "s"
    assert siesta_sym_to_name("px") == "px"
    assert siesta_sym_to_name("Pdx2-y2") == "dx2y2"   # strips polarisation prefix
    assert siesta_sym_to_name("Pdz2") == "dz2"


def test_real_harmonics_orthonormal():
    """Canonical real harmonics integrate to delta_ij over the unit sphere."""
    rng = np.linspace(0, 1, 200)
    # Fibonacci sphere for a quick angular quadrature
    n = 20000
    i = np.arange(n) + 0.5
    phi = np.arccos(1 - 2 * i / n)
    theta = np.pi * (1 + 5 ** 0.5) * i
    nx = np.sin(phi) * np.cos(theta)
    ny = np.sin(phi) * np.sin(theta)
    nz = np.cos(phi)
    names = ["s", "px", "py", "pz", "dz2", "dx2y2", "dxy", "dxz", "dyz",
             "fz3", "fxz2", "fyz2", "fzx2y2", "fxyz", "fx3", "fy3"]
    Y = np.array([_angular(nm, nx, ny, nz) for nm in names])
    G = (Y @ Y.T) * (4 * np.pi / n)        # Monte-Carlo-ish surface integral
    assert np.allclose(np.diag(G), 1.0, atol=0.05)
    off = G - np.diag(np.diag(G))
    assert np.max(np.abs(off)) < 0.05


def test_dm_roundtrip():
    dm = read_dm(TD / "si.DM")
    assert dm.no_u == 26 and dm.nspin == 1 and dm.nsc == (5, 5, 5)
    with tempfile.TemporaryDirectory() as d:
        p = Path(d) / "out.DM"
        write_dm(p, dm)
        # byte-identical roundtrip
        assert p.read_bytes() == (TD / "si.DM").read_bytes()


def test_copy_pattern_zeroed():
    dm = read_dm(TD / "si.DM")
    z = copy_pattern(dm)
    assert np.array_equal(z.listd, dm.listd)
    assert np.allclose(z.dm, 0.0)


# ----------------------- bulk-Si fixture based ------------------------------
_scfout = TD / "Si.scfout"


def test_overlap_reproduces_openmx_olp():
    """Overlap engine must reproduce OpenMX's own on-site + NN overlap matrix."""
    from cube4siesta.basis import build_openmx_orbitals, read_pao_radials
    from cube4siesta.overlap import block_overlap_matrix, onsite_overlap_matrix
    from cube4siesta.scfout_io import read_scfout

    sc = read_scfout(_scfout)
    orbs = build_openmx_orbitals(parse_spec("s2p2d1"), read_pao_radials(TD / "Si7.0.pao"))
    # on-site
    assert np.max(np.abs(onsite_overlap_matrix(orbs) - sc.OLP[(1, 0)])) < 1e-4
    # nearest neighbour
    gx = sc.gxyz
    best = None
    for h in range(1, int(sc.fnan[1]) + 1):
        Gh = int(sc.natn[1][h]); rn = int(sc.ncn[1][h])
        cB = gx[Gh][1:4] + sc.atv[rn][1:4]
        dd = float(np.linalg.norm(gx[1][1:4] - cB))
        if dd > 0.1 and (best is None or dd < best[0]):
            best = (dd, h, cB)
    _, h, cB = best
    blk = block_overlap_matrix(orbs, gx[1][1:4], orbs, cB, spacing=0.15)
    assert np.max(np.abs(blk - sc.OLP[(1, h)])) < 1e-3


def _random_spd(rng, n):
    A = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    return A @ A.conj().T / n + np.eye(n)


def _dm_with_occupations(rng, S, f):
    """Hermitian P whose natural occupations in the S metric are exactly f."""
    s, V = np.linalg.eigh(S)
    Sm12 = (V / np.sqrt(s)) @ V.conj().T
    H = rng.normal(size=S.shape) + 1j * rng.normal(size=S.shape)
    _, W = np.linalg.eigh(H + H.conj().T)
    X = Sm12 @ W
    return (X * f) @ X.conj().T


def test_purify_canonical_insulating_filling():
    """Fuzzy occupations snap to {0, cap}, trace lands exactly on the target."""
    from cube4siesta.project import purify_canonical

    rng = np.random.default_rng(42)
    nk, n, cap = 4, 6, 2.0
    Sk = np.array([_random_spd(rng, n) for _ in range(nk)])
    P = np.empty((nk, n, n), dtype=complex)
    for ik in range(nk):
        # 4 nearly-full + 2 nearly-empty states, away from idempotency
        f = np.concatenate([1.7 + 0.2 * rng.random(4), 0.15 * rng.random(2)])
        P[ik] = _dm_with_occupations(rng, Sk[ik], f)

    occ_min, occ_max, n_frac = purify_canonical(P, Sk, n_channel=8.0, occ_cap=cap)

    ne = np.real(np.einsum("kab,kba->", P, Sk)) / nk
    assert abs(ne - 8.0) < 1e-9                  # trace constraint exact
    assert n_frac == 0                           # insulating: no boundary sharing
    assert 0.0 <= occ_min < 0.2 and 1.7 < occ_max < 2.0
    for ik in range(nk):                         # idempotent: P S P == cap * P
        assert np.allclose(P[ik] @ Sk[ik] @ P[ik], cap * P[ik], atol=1e-9)
        assert np.allclose(P[ik], P[ik].conj().T, atol=1e-12)


def test_purify_canonical_fractional_boundary():
    """Odd electron count -> fractional state(s) at the Fermi boundary."""
    from cube4siesta.project import purify_canonical

    rng = np.random.default_rng(7)
    nk, n, cap = 3, 6, 2.0
    Sk = np.array([_random_spd(rng, n) for _ in range(nk)])
    P = np.empty((nk, n, n), dtype=complex)
    for ik in range(nk):
        P[ik] = _dm_with_occupations(rng, Sk[ik], 2.0 * rng.random(n))

    # 7 e/cell with cap 2 -> 10 full pooled states + 1 e left over
    _, _, n_frac = purify_canonical(P, Sk, n_channel=7.0, occ_cap=cap)

    ne = np.real(np.einsum("kab,kba->", P, Sk)) / nk
    assert abs(ne - 7.0) < 1e-9
    assert n_frac >= 1


def test_purify_canonical_capacity_error():
    from cube4siesta.project import purify_canonical

    rng = np.random.default_rng(0)
    Sk = np.array([_random_spd(rng, 4)])
    P = np.array([_dm_with_occupations(rng, Sk[0], np.full(4, 1.0))])
    with pytest.raises(ValueError):
        purify_canonical(P, Sk, n_channel=9.0, occ_cap=2.0)   # > 4*2


def test_projection_conserves_electrons_and_hermitian():
    from cube4siesta.basis import (build_openmx_orbitals, build_siesta_orbitals,
                                   read_ion_radials, read_pao_radials)
    from cube4siesta.orb_indx import read_orb_indx
    from cube4siesta.project import project_openmx_to_siesta
    from cube4siesta.scfout_io import read_scfout

    sc = read_scfout(_scfout)
    oi = read_orb_indx(TD / "si.ORB_INDX")
    tmpl = read_dm(TD / "si.DM")
    sorb = build_siesta_orbitals(oi.unit, {"Si": read_ion_radials(TD / "Si.ion")})
    oorb = build_openmx_orbitals(parse_spec("s2p2d1"), read_pao_radials(TD / "Si7.0.pao"))
    res = project_openmx_to_siesta(sc, oi, tmpl, sorb, [oorb, oorb], spacing=0.2)
    assert abs(res.n_electrons - 8.0) < 0.1     # matched valence -> ~8 e
    assert res.max_asym < 1e-8                   # Hermitian P(R)_ab == P(-R)_ba
