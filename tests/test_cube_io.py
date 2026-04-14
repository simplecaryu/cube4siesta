"""Cube roundtrip + field integrity tests."""

from __future__ import annotations

import numpy as np
import pytest

from cube4siesta.cube_io import CubeFile, read_cube, write_cube


def _make_cube(n1=4, n2=5, n3=6, nat=2):
    rng = np.random.default_rng(42)
    origin = np.array([0.0, 0.0, 0.0])
    voxel = np.diag([0.2, 0.25, 0.3])
    atoms = np.array([
        [1, 1.0, 0.5, 0.5, 0.5],
        [8, 6.0, 1.0, 1.2, 1.4],
    ])[:nat]
    data = rng.random((n1, n2, n3))
    return CubeFile(
        comment1="test cube",
        comment2="generated",
        origin=origin,
        voxel=voxel,
        mesh=(n1, n2, n3),
        atoms=atoms,
        data=data,
    )


def test_cube_roundtrip(tmp_path):
    cube = _make_cube()
    path = tmp_path / "t.cube"
    write_cube(path, cube)
    out = read_cube(path)

    assert out.mesh == cube.mesh
    np.testing.assert_allclose(out.origin, cube.origin)
    np.testing.assert_allclose(out.voxel, cube.voxel)
    np.testing.assert_allclose(out.atoms, cube.atoms)
    # text formatting is %.5e -> ~5 decimal digits
    np.testing.assert_allclose(out.data, cube.data, rtol=1e-4)


def test_cube_cell_and_norm(tmp_path):
    cube = _make_cube(n1=3, n2=3, n3=3)
    cube.data[:] = 1.0  # constant field
    path = tmp_path / "uniform.cube"
    write_cube(path, cube)
    out = read_cube(path)

    expected_cell = np.diag([3 * 0.2, 3 * 0.25, 3 * 0.3])
    np.testing.assert_allclose(out.cell, expected_cell)

    # integral of constant rho=1 over volume = volume
    volume = 3 * 0.2 * 3 * 0.25 * 3 * 0.3
    assert out.n_electrons == pytest.approx(volume, rel=1e-5)


def test_cube_ordering_is_i3_fastest(tmp_path):
    """cube spec: i3 (axis-3) is the fastest-varying index in the data block."""
    n1, n2, n3 = 2, 2, 3
    data = np.arange(n1 * n2 * n3, dtype=np.float64).reshape(n1, n2, n3)
    cube = CubeFile(
        comment1="ord", comment2="test",
        origin=np.zeros(3),
        voxel=np.eye(3),
        mesh=(n1, n2, n3),
        atoms=np.zeros((0, 5)),
        data=data,
    )
    path = tmp_path / "ord.cube"
    write_cube(path, cube)

    # parse raw file: after the header (6 lines, no atoms) the first data line
    # should begin with data[0,0,0], data[0,0,1], data[0,0,2], ... (i3 fastest)
    text = path.read_text().splitlines()
    header = 6
    first_data = text[header].split()
    assert float(first_data[0]) == data[0, 0, 0]
    assert float(first_data[1]) == data[0, 0, 1]
    assert float(first_data[2]) == data[0, 0, 2]

    out = read_cube(path)
    np.testing.assert_array_equal(out.data, data)
