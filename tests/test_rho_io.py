"""
Tests for SIESTA .RHO reader/writer.

Strategy:
  1. roundtrip: write_rho -> read_rho yields identical data
  2. hand-crafted layout: manually parse the raw bytes and verify each record
  3. real-file: read the h2o.RHO produced by a live SIESTA run (if present)
"""

from __future__ import annotations

import struct
from pathlib import Path

import numpy as np
import pytest

from cube4siesta.rho_io import read_rho, write_rho


# ---------------------------------------------------------------------------
# 1. Roundtrip

def test_roundtrip_nspin1(tmp_path):
    rng = np.random.default_rng(0)
    cell = np.array([[10.0, 0.0, 0.0],
                     [0.0, 12.0, 0.0],
                     [0.0, 0.0, 15.0]])
    rho = rng.random((6, 8, 10), dtype=np.float32)

    path = tmp_path / "test.RHO"
    write_rho(path, cell, rho)
    out = read_rho(path)

    assert out.mesh == (6, 8, 10)
    assert out.nspin == 1
    np.testing.assert_allclose(out.cell, cell)
    np.testing.assert_array_equal(out.rho[..., 0], rho)


def test_roundtrip_nspin2(tmp_path):
    rng = np.random.default_rng(1)
    cell = np.eye(3) * 7.5
    rho = rng.random((4, 5, 6, 2), dtype=np.float32)

    path = tmp_path / "test_spin.RHO"
    write_rho(path, cell, rho)
    out = read_rho(path)

    assert out.mesh == (4, 5, 6)
    assert out.nspin == 2
    np.testing.assert_array_equal(out.rho, rho)


def test_non_orthogonal_cell(tmp_path):
    # triclinic-ish — verify cell element ordering survives
    cell = np.array([[5.1, 0.0, 0.0],
                     [2.5, 4.4, 0.0],
                     [0.3, 0.7, 6.2]])
    rho = np.arange(3 * 4 * 5, dtype=np.float32).reshape(3, 4, 5)

    path = tmp_path / "tri.RHO"
    write_rho(path, cell, rho)
    out = read_rho(path)

    np.testing.assert_allclose(out.cell, cell)
    assert out.cell[1, 0] == pytest.approx(2.5)  # a2_x component preserved
    assert out.cell[2, 1] == pytest.approx(0.7)  # a3_y component preserved


# ---------------------------------------------------------------------------
# 2. Raw byte-layout verification (independent of scipy.io.FortranFile)
#
# Fortran unformatted sequential records use a 4-byte length marker before
# AND after each record's payload (default GCC/ifort). We parse manually.

def _read_fortran_record(f):
    head = f.read(4)
    if len(head) < 4:
        return None
    n = struct.unpack("<i", head)[0]
    data = f.read(n)
    tail = struct.unpack("<i", f.read(4))[0]
    assert tail == n, f"record length mismatch: head={n} tail={tail}"
    return data


def test_raw_byte_layout(tmp_path):
    # small known array
    cell = np.array([[1.0, 2.0, 3.0],
                     [4.0, 5.0, 6.0],
                     [7.0, 8.0, 9.0]])
    rho = np.array([[[0.1, 0.2], [0.3, 0.4]],
                    [[0.5, 0.6], [0.7, 0.8]]], dtype=np.float32)  # shape (2,2,2)

    path = tmp_path / "raw.RHO"
    write_rho(path, cell, rho)

    with open(path, "rb") as f:
        # Record 1: 9 doubles = 72 bytes, Fortran column-major order
        r1 = _read_fortran_record(f)
        assert len(r1) == 72
        vals = struct.unpack("<9d", r1)
        # Expected flat order: a1x,a1y,a1z, a2x,a2y,a2z, a3x,a3y,a3z
        assert vals == (1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0)

        # Record 2: 4 int32 = 16 bytes
        r2 = _read_fortran_record(f)
        assert len(r2) == 16
        mesh_nspin = struct.unpack("<4i", r2)
        assert mesh_nspin == (2, 2, 2, 1)

        # Records 3..: nspin(1) * nz(2) * ny(2) = 4 records, each 2 float32 = 8 bytes
        # loop order spin -> z -> y; per record = x-row
        records = []
        while True:
            rec = _read_fortran_record(f)
            if rec is None:
                break
            assert len(rec) == 8  # nx * 4 bytes
            records.append(struct.unpack("<2f", rec))

        assert len(records) == 4
        # ispin=0, iz=0, iy=0 -> rho[:,0,0,0] = [0.1, 0.5]
        # ispin=0, iz=0, iy=1 -> rho[:,1,0,0] = [0.3, 0.7]
        # ispin=0, iz=1, iy=0 -> rho[:,0,1,0] = [0.2, 0.6]
        # ispin=0, iz=1, iy=1 -> rho[:,1,1,0] = [0.4, 0.8]
        np.testing.assert_allclose(records[0], [0.1, 0.5], rtol=1e-6)
        np.testing.assert_allclose(records[1], [0.3, 0.7], rtol=1e-6)
        np.testing.assert_allclose(records[2], [0.2, 0.6], rtol=1e-6)
        np.testing.assert_allclose(records[3], [0.4, 0.8], rtol=1e-6)


# ---------------------------------------------------------------------------
# 3. Real SIESTA output

H2O_RHO = Path("/tmp/siesta_test/h2o.RHO")


@pytest.mark.skipif(not H2O_RHO.exists(), reason="h2o.RHO not present")
def test_read_real_h2o_rho():
    out = read_rho(H2O_RHO)

    # h2o.fdf uses a cubic cell ~10 Ang and MeshCutoff; nspin=1
    assert out.nspin == 1
    assert out.cell.shape == (3, 3)
    assert all(m > 0 for m in out.mesh)

    # H2O has 8 valence electrons (O: 6, H: 1+1). integral should be close.
    n_elec = out.n_electrons
    assert 7.0 < n_elec < 9.0, f"electron count out of range: {n_elec}"

    # density should be non-negative (up to tiny FP noise)
    assert out.rho.min() > -1e-6
