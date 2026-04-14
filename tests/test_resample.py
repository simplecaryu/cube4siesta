"""Resampling sanity tests."""

from __future__ import annotations

import numpy as np
import pytest

from cube4siesta.resample import resample_to_mesh


def test_identity_when_meshes_match():
    rng = np.random.default_rng(0)
    a = rng.random((5, 7, 9))
    b = resample_to_mesh(a, (5, 7, 9))
    np.testing.assert_array_equal(a, b)


def test_uniform_preserved():
    a = np.full((6, 8, 10), 0.37)
    b = resample_to_mesh(a, (12, 16, 20))
    np.testing.assert_allclose(b, 0.37, rtol=1e-10)


def test_charge_integral_preserved_for_smooth_field():
    """A smooth field's integral should be conserved across upsampling."""
    n = 16
    x = (np.arange(n) + 0.5) / n * 2 * np.pi
    field = (np.sin(x)[:, None, None] * np.sin(x)[None, :, None]
             * np.sin(x)[None, None, :]) ** 2
    voxel_volume = 1.0 / n ** 3
    integ_src = field.sum() * voxel_volume

    up = resample_to_mesh(field, (32, 32, 32))
    voxel_volume_up = 1.0 / 32 ** 3
    integ_up = up.sum() * voxel_volume_up

    # trilinear interpolation is not exactly conservative, but for a smooth
    # oscillatory field on a fine-enough grid it should agree to a few %
    assert integ_up == pytest.approx(integ_src, rel=2e-2)
