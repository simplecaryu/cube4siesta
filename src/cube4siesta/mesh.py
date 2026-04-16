"""
Predict the SIESTA FFT mesh dimensions from cell vectors and MeshCutoff.

Replicates the logic in siesta-4.1.5/Src/meshsubs.F:get_mesh_sizes +
fft1d.F:nfft + chkgmx.f:chkgmx so that cube4siesta can determine the
target mesh without running SIESTA.

Units: cell in Bohr, MeshCutoff in Ry (same as SIESTA fdf).
"""

from __future__ import annotations

import numpy as np


def _next_fft_size(n: int, nsm: int = 2) -> int:
    """Round n up to the next integer that is a product of 2,3,5 AND a
    multiple of nsm (SIESTA's sub-mesh size, default 2)."""
    while True:
        if _is_fft_friendly(n) and n % nsm == 0:
            return n
        n += 1


def _is_fft_friendly(n: int) -> bool:
    if n <= 0:
        return False
    for p in (2, 3, 5):
        while n % p == 0:
            n //= p
    return n == 1


def _chkgmx(rcell: np.ndarray, mesh: tuple[int, int, int]) -> float:
    """Compute the actual planewave cutoff G2max for a given mesh on a
    (possibly non-orthogonal) cell. Mirrors chkgmx.f.

    Returns G2max in Ry (the square of the cutoff wavevector).
    """
    bm = np.zeros((3, 3))
    for i in range(3):
        bm[:, i] = rcell[:, i] * mesh[i]

    gmax = np.inf
    for i1 in range(-1, 2):
        for i2 in range(-1, 2):
            for i3 in range(-1, 2):
                if i1 == 0 and i2 == 0 and i3 == 0:
                    continue
                g = bm[:, 0] * i1 + bm[:, 1] * i2 + bm[:, 2] * i3
                gmod = np.linalg.norm(g)
                r = 0.5 * gmod  # k=0
                gmax = min(gmax, r)

    return gmax * gmax - 1e-8


def siesta_mesh(
    cell: np.ndarray,
    mesh_cutoff: float = 200.0,
    nsm: int = 2,
) -> tuple[int, int, int]:
    """Predict the SIESTA FFT mesh for a given cell and MeshCutoff.

    Parameters
    ----------
    cell : (3,3) array, lattice vectors as rows, in Bohr
    mesh_cutoff : MeshCutoff in Ry (default 200, SIESTA's common default)
    nsm : sub-mesh sampling (default 2, SIESTA default)

    Returns
    -------
    (nx, ny, nz) : mesh dimensions (including sub-points, i.e. ntm)
    """
    cell = np.asarray(cell, dtype=np.float64)
    g2max = mesh_cutoff

    # Reciprocal lattice vectors (with 2pi), stored as columns.
    # SIESTA convention: rcell(:,i) · cell(:,j) = 2π δ_ij
    # With Python row-vectors (cell[i,:] = i-th vector):
    #   rcell_fortran = 2π · inv(cell_python)
    rcell = 2.0 * np.pi * np.linalg.inv(cell)

    # Initial estimate: ntm[i] = 2*sqrt(G2max) / |b_i| + 1
    ntm = [0, 0, 0]
    for i in range(3):
        vecmod = np.linalg.norm(rcell[:, i])
        ntm[i] = int(2.0 * np.sqrt(g2max) / vecmod) + 1

    # Iteratively round up to FFT-friendly + nsm-multiple,
    # then verify effective cutoff (non-orthogonal correction)
    for _ in range(50):
        for i in range(3):
            ntm[i] = _next_fft_size(ntm[i], nsm)

        real_cutoff = _chkgmx(rcell, tuple(ntm))
        if real_cutoff >= g2max:
            break
        # Cutoff too small for this mesh — increase and retry
        for i in range(3):
            ntm[i] += 1
    else:
        raise RuntimeError("siesta_mesh: could not find suitable mesh in 50 iterations")

    return (ntm[0], ntm[1], ntm[2])
