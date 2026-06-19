"""
Build a SIESTA density-matrix *restart seed* by fitting a cross-code difference
density onto the SIESTA LCAO basis.

Motivation
----------
Transferring a converged density matrix between codes that use different
pseudopotentials is hard precisely because of semicore states: OpenMX may treat
e.g. V 3s/3p as valence (25 electrons) while the SIESTA pseudopotential freezes
them into the core (17 electrons). A direct density-matrix projection transplants
that semicore information and corrupts the valence block.

This module takes the *difference* density Δρ = ρ_crystal − ρ_atoms (which any
DFT package can emit, and in which the frozen semicore largely cancels) and
*fits* it onto the SIESTA basis:

    Δρ(r) ≈ Σ_i Σ_{(j,R)} x_ij(R) · col_{i,j,R}(r),
    col_{i,j,R}(r) = Σ_n φ_i(r−c_i−n) φ_j(r−c_j−R−n)      (n = lattice images)

The unknowns ``x`` live on the sparsity pattern of a SIESTA ``.DM`` template. The
least-squares fit is a *basis filter*: any density the SIESTA basis cannot
represent — semicore cusps near the nucleus among them — simply drops out. Adding
the SIESTA neutral-atom density matrix back gives a seed

    D_seed = D_atomic + fit(Δρ),

which SIESTA reads with ``DM.UseSaveDM`` and refines with a single
diagonalisation (``MaxSCFIterations 1``), no source density matrix and no
SIESTA source-code patch required. See :func:`restart_fdf_block`.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import scipy.sparse as sp
from scipy.ndimage import map_coordinates
from scipy.sparse.linalg import lsqr

from .basis import build_siesta_orbitals, read_ion_populations, read_ion_radials
from .cube_io import read_cube
from .dm_io import DMFile, copy_pattern, read_dm
from .orb_indx import OrbIndx, read_orb_indx
from .rho_io import read_rho
from .struct_out import read_struct_out


# --------------------------------------------------------------------------
# Geometry / basis bundle
# --------------------------------------------------------------------------
@dataclass
class Baseline:
    """A converged SIESTA baseline run: geometry, basis and .DM sparsity pattern.

    Only the *pattern* of ``dm`` (numd/listd) is used to define the unknowns; its
    values double as the converged reference for the optional direct-error
    diagnostic.
    """
    dm: DMFile
    oi: OrbIndx
    orbitals: list             # built SIESTA orbitals, one per unit-cell orbital
    frac: np.ndarray           # (natoms, 3) fractional coordinates
    atom_of_uorb: np.ndarray   # (no_u,) 1-based atom index per unit-cell orbital
    pops_by_species: dict      # {species: {(n,l,z): shell population}}


def load_baseline(base_dir: str | Path, stem: str) -> Baseline:
    """Load a SIESTA baseline from ``<base_dir>/<stem>.{DM,ORB_INDX,STRUCT_OUT}``
    plus the per-species ``<species>.ion`` files in the same directory."""
    base = Path(base_dir)
    dm = read_dm(base / f"{stem}.DM")
    oi = read_orb_indx(base / f"{stem}.ORB_INDX")
    frac = read_struct_out(base / f"{stem}.STRUCT_OUT").frac
    species = sorted({o.spec for o in oi.unit})
    radials = {s: read_ion_radials(base / f"{s}.ion") for s in species}
    pops = {s: read_ion_populations(base / f"{s}.ion") for s in species}
    orbitals = build_siesta_orbitals(oi.unit, radials)
    return Baseline(
        dm=dm,
        oi=oi,
        orbitals=orbitals,
        frac=frac,
        atom_of_uorb=np.array([o.ia for o in oi.unit], dtype=np.int64),
        pops_by_species=pops,
    )


# --------------------------------------------------------------------------
# Real-space resampling (corner-aligned)
# --------------------------------------------------------------------------
def resample_corner(src: np.ndarray, target_mesh: tuple[int, int, int]) -> np.ndarray:
    """
    Resample a periodic density onto ``target_mesh``, corner-aligned: grid point
    ``k`` sits at fractional coordinate ``k / N`` on both source and target. This
    matches the SIESTA .RHO and OpenMX .cube grid convention and avoids the
    half-voxel shift a centre-aligned resample would introduce.
    """
    ns = src.shape
    coords = np.stack(
        np.meshgrid(
            *[np.arange(target_mesh[d]) * ns[d] / target_mesh[d] for d in range(3)],
            indexing="ij",
        ),
        axis=0,
    )
    return map_coordinates(np.asarray(src, dtype=np.float64), coords,
                           order=1, mode="grid-wrap")


# --------------------------------------------------------------------------
# Design matrix: one column per .DM entry = its real-space density
# --------------------------------------------------------------------------
def build_design_matrix(
    base: Baseline,
    cell: np.ndarray,
    mesh: tuple[int, int, int],
    verbose: bool = False,
):
    """
    Build the sparse design matrix ``A`` (ngrid × nnz). Column ``c`` is the
    real-space density of density-matrix entry ``D_ij(R)`` summed over lattice
    images, evaluated on the grid defined by ``cell`` (Bohr) and ``mesh``.

    Returns ``(A_csc, dV)`` where ``dV`` is the voxel volume (Bohr^3).
    """
    nx, ny, nz = mesh
    ngrid = nx * ny * nz
    cell = np.asarray(cell, dtype=np.float64)
    atom_cart = base.frac @ cell
    no_u = base.dm.no_u

    fi, fj, fk = np.arange(nx) / nx, np.arange(ny) / ny, np.arange(nz) / nz
    FI, FJ, FK = np.meshgrid(fi, fj, fk, indexing="ij")
    gcart = np.stack([FI.ravel(), FJ.ravel(), FK.ravel()], 1) @ cell

    rcut = np.array([o.rcut for o in base.orbitals])
    a_len = np.linalg.norm(cell, axis=1)
    nr = [int(np.ceil(rcut.max() / a_len[d])) + 1 for d in range(3)]

    # localized orbital-image fields g[(u, n1, n2, n3)] = (grid idx, value)
    g: dict = {}
    for u in range(no_u):
        c0 = atom_cart[base.atom_of_uorb[u] - 1]
        orb = base.orbitals[u]
        for n1 in range(-nr[0], nr[0] + 1):
            for n2 in range(-nr[1], nr[1] + 1):
                for n3 in range(-nr[2], nr[2] + 1):
                    c = c0 + n1 * cell[0] + n2 * cell[1] + n3 * cell[2]
                    d = gcart - c
                    val = orb.evaluate(d[:, 0], d[:, 1], d[:, 2])
                    idx = np.nonzero(val)[0]
                    if idx.size:
                        g[(u, n1, n2, n3)] = (idx.astype(np.int32), val[idx])
    g_by_orb: dict = {io: [] for io in range(no_u)}
    for k, v in g.items():
        g_by_orb[k[0]].append((k, v))

    rows, cols, data = [], [], []
    col = 0
    for io in range(no_u):
        listd, _ = base.dm.row(io)
        for jj in range(listd.size):
            juo, isc = base.oi.decode_column(int(listd[jj]))
            j = juo - 1
            tmp = np.zeros(ngrid)
            touched = []
            for (u, m1, m2, m3), (idx_i, val_i) in g_by_orb[io]:
                gj = g.get((j, m1 + isc[0], m2 + isc[1], m3 + isc[2]))
                if gj is None:
                    continue
                idx_j, val_j = gj
                common, ia_, ib_ = np.intersect1d(
                    idx_i, idx_j, assume_unique=True, return_indices=True)
                if common.size:
                    tmp[common] += val_i[ia_] * val_j[ib_]
                    touched.append(common)
            if touched:
                gi = np.unique(np.concatenate(touched))
                rows.append(gi)
                cols.append(np.full(gi.size, col, np.int32))
                data.append(tmp[gi])
            col += 1
        if verbose and (io + 1) % max(1, no_u // 10) == 0:
            print(f"    design matrix: row {io + 1}/{no_u}", flush=True)

    A = sp.csr_matrix(
        (np.concatenate(data), (np.concatenate(rows), np.concatenate(cols))),
        shape=(ngrid, base.dm.nnz),
    ).tocsc()
    dV = abs(np.linalg.det(cell)) / ngrid
    return A, dV


# --------------------------------------------------------------------------
# Neutral-atom density matrix and the least-squares fit
# --------------------------------------------------------------------------
def build_atomic_dm(base: Baseline) -> np.ndarray:
    """
    SIESTA neutral-atom initial density matrix as a vector on the .DM pattern:
    diagonal, each ``m``-orbital of a shell gets ``pop(n,l,z) / (2l+1)``;
    off-diagonal and off-cell entries are zero.
    """
    x = np.zeros(base.dm.nnz)
    for io in range(base.dm.no_u):
        o = base.oi.unit[io]
        occ = base.pops_by_species[o.spec].get((o.n, o.l, o.z), 0.0) / (2 * o.l + 1)
        if occ == 0.0:
            continue
        listd, _ = base.dm.row(io)
        p = int(base.dm.listdptr[io])
        for jj in range(listd.size):
            juo, isc = base.oi.decode_column(int(listd[jj]))
            if juo - 1 == io and isc == (0, 0, 0):     # onsite diagonal entry
                x[p + jj] = occ
                break
    return x


def fit_density(
    A: sp.spmatrix,
    rho_target: np.ndarray,
    damp: float = 0.1,
    atol: float = 1e-7,
    btol: float = 1e-7,
    iter_lim: int = 2000,
) -> np.ndarray:
    """Tikhonov-damped least-squares fit of ``rho_target`` onto the columns of
    ``A`` (the basis filter). ``damp`` regularises the basis null-space."""
    return lsqr(A, rho_target, damp=damp, atol=atol, btol=btol, iter_lim=iter_lim)[0]


# --------------------------------------------------------------------------
# High-level driver
# --------------------------------------------------------------------------
@dataclass
class SeedResult:
    dm: DMFile
    direct_l2: float          # ||D_seed - D_converged|| / ||D_converged|| (pre-SCF)
    filtered_integral: float  # ∫ filtered Δρ dV (electrons); ~0 for a difference density
    max_abs: float            # max |D_seed| entry
    nnz: int


def _load_dden_on_grid(dden_path: str | Path, mesh: tuple[int, int, int]) -> np.ndarray:
    """Load a difference density (.cube or .RHO) and corner-resample to ``mesh``."""
    path = Path(dden_path)
    if path.suffix.lower() == ".cube":
        data = read_cube(path).data
    else:
        data = read_rho(path).rho[..., 0]
    if tuple(data.shape) == tuple(mesh):
        return data.ravel()
    return resample_corner(data, mesh).ravel()


def make_filtered_dm_seed(
    base_dir: str | Path,
    stem: str,
    dden_path: str | Path,
    damp: float = 0.1,
    verbose: bool = False,
) -> SeedResult:
    """
    Build a filtered density-matrix restart seed.

    Reads the SIESTA baseline ``<base_dir>/<stem>.*`` (for geometry, basis and the
    .DM sparsity pattern), filters the cross-code difference density ``dden_path``
    through that basis, and returns ``D_seed = D_atomic + fit(Δρ)`` as a
    :class:`~cube4siesta.dm_io.DMFile` ready to be written with
    :func:`~cube4siesta.dm_io.write_dm`.
    """
    base = load_baseline(base_dir, stem)
    rho = read_rho(Path(base_dir) / f"{stem}.RHO")
    if verbose:
        print(f"[dm-seed] no_u={base.dm.no_u} nnz={base.dm.nnz} "
              f"mesh={rho.mesh} nsc={base.dm.nsc}", flush=True)

    A, dV = build_design_matrix(base, rho.cell, rho.mesh, verbose=verbose)
    if verbose:
        print(f"[dm-seed] built design matrix nnz={A.nnz}", flush=True)

    rho_dd = _load_dden_on_grid(dden_path, rho.mesh)
    x_delta = fit_density(A, rho_dd, damp=damp)
    d_atom = build_atomic_dm(base)
    d_seed = d_atom + x_delta

    dmf = copy_pattern(base.dm, nspin=1)
    dmf.dm[:, 0] = d_seed

    x_true = base.dm.dm[:, 0]
    direct_l2 = float(np.linalg.norm(d_seed - x_true) / np.linalg.norm(x_true))
    filtered_integral = float((A @ x_delta).sum() * dV)

    return SeedResult(
        dm=dmf,
        direct_l2=direct_l2,
        filtered_integral=filtered_integral,
        max_abs=float(np.abs(d_seed).max()),
        nnz=base.dm.nnz,
    )


def restart_fdf_block(system_label: str = "<SystemLabel>") -> str:
    """Return the SIESTA settings that consume a seed .DM in a single,
    patch-free diagonalisation step."""
    return (
        "# --- filtered DM restart (patch-free) ---\n"
        "# Rename the seed to <SystemLabel>.DM in the run directory, then:\n"
        "DM.UseSaveDM         true\n"
        "MaxSCFIterations     1\n"
        "DM.MixingWeight      1.0\n"
        "SCFMustConverge      false\n"
        "# SIESTA rebuilds rho[D_seed] and diagonalises once; the output\n"
        f"# {system_label}.DM and bands are the one-step refined result.\n"
    )
