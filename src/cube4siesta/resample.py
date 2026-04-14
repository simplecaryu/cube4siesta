"""
Resample a density from the cube's native grid onto a target mesh.

MVP assumption: cube voxel vectors are parallel to the target cell's
lattice vectors, so interpolation reduces to a trilinear resample along
each axis independently. Non-orthogonal alignment is left to a follow-up.
"""

from __future__ import annotations

import numpy as np
from scipy.ndimage import map_coordinates


def resample_to_mesh(
    data: np.ndarray,
    target_mesh: tuple[int, int, int],
    *,
    order: int = 1,
) -> np.ndarray:
    """
    Resample a 3D density from its native grid shape to target_mesh.

    Parameters
    ----------
    data : (Ns1, Ns2, Ns3) array on the source grid
    target_mesh : (N1, N2, N3) target grid
    order : interpolation order (1 = trilinear, 3 = cubic spline)

    Returns
    -------
    out : (N1, N2, N3) array on the target grid

    Notes
    -----
    Uses periodic wrap (mode='grid-wrap') — appropriate for periodic DFT
    densities produced by plane-wave codes. If source and target shapes
    match exactly, returns a copy unchanged.
    """
    src = np.asarray(data, dtype=np.float64)
    n1s, n2s, n3s = src.shape
    n1, n2, n3 = target_mesh

    if (n1s, n2s, n3s) == (n1, n2, n3):
        return src.copy()

    # Coordinates in source-index space for each target voxel. Target voxel
    # centers at fractional positions (i+0.5)/N map to source fractional
    # (i+0.5)/N * Ns - 0.5 so that voxel centers align, not corners.
    i1 = (np.arange(n1) + 0.5) / n1 * n1s - 0.5
    i2 = (np.arange(n2) + 0.5) / n2 * n2s - 0.5
    i3 = (np.arange(n3) + 0.5) / n3 * n3s - 0.5
    coords = np.stack(np.meshgrid(i1, i2, i3, indexing="ij"), axis=0)

    out = map_coordinates(src, coords, order=order, mode="grid-wrap")
    return out
