"""Tests for the difference-density conversion path (--diff / --subtract)."""

from __future__ import annotations

import numpy as np

from cube4siesta.cli import main
from cube4siesta.cube_io import CubeFile, write_cube
from cube4siesta.rho_io import read_rho


def _write_cube(path, data):
    n1, n2, n3 = data.shape
    cube = CubeFile(
        comment1="test cube",
        comment2="generated",
        origin=np.zeros(3),
        voxel=np.diag([0.2, 0.25, 0.3]),
        mesh=(n1, n2, n3),
        atoms=np.array([[1, 1.0, 0.5, 0.5, 0.5]]),
        data=data,
    )
    write_cube(path, cube)
    return cube


def test_convert_diff_passthrough(tmp_path):
    """--diff converts a charge-neutral field verbatim, no rescaling."""
    rng = np.random.default_rng(1)
    data = rng.normal(size=(4, 5, 6))
    data -= data.mean()                      # integral exactly 0
    cube_path = tmp_path / "dden.cube"
    _write_cube(cube_path, data)

    out = tmp_path / "t.RHOIN.diff"
    rc = main(["convert", "--cube", str(cube_path),
               "--output", str(out), "--diff", "--verify"])
    assert rc == 0
    back = read_rho(out)
    np.testing.assert_allclose(back.rho[..., 0], data, rtol=1e-4, atol=1e-7)


def test_convert_diff_rejects_rescale(tmp_path):
    cube_path = tmp_path / "dden.cube"
    _write_cube(cube_path, np.zeros((3, 3, 3)))
    rc = main(["convert", "--cube", str(cube_path),
               "--output", str(tmp_path / "t.RHOIN"),
               "--diff", "--rescale-to", "8.0"])
    assert rc == 2


def test_convert_subtract_builds_difference(tmp_path):
    """--subtract REF writes (cube - REF), e.g. scf minus atomic totals."""
    rng = np.random.default_rng(2)
    atomic = rng.random((4, 4, 4)) + 1.0
    bonding = rng.normal(size=(4, 4, 4)) * 0.01
    bonding -= bonding.mean()
    scf_path = tmp_path / "scf.cube"
    ref_path = tmp_path / "atomic.cube"
    _write_cube(scf_path, atomic + bonding)
    _write_cube(ref_path, atomic)

    out = tmp_path / "t.RHOIN.diff"
    rc = main(["convert", "--cube", str(scf_path), "--subtract", str(ref_path),
               "--output", str(out)])
    assert rc == 0
    back = read_rho(out)
    # cube text format keeps ~5 significant digits of the ~1.0 totals,
    # so the difference is only good to ~1e-5 absolute
    np.testing.assert_allclose(back.rho[..., 0], bonding, atol=5e-5)
