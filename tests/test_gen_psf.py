"""Tests for gen_psf: parse source pseudo → ATOM input."""

from __future__ import annotations

import os
from pathlib import Path

import pytest

from cube4siesta.gen_psf import (
    parse_openmx_vps,
    parse_qe_upf,
    write_atom_input,
)

# ---------------------------------------------------------------------------
# OpenMX .vps
# ---------------------------------------------------------------------------

_VPS_DEFAULT = "/home/users2/cha/programs/openmx3.9/DFT_DATA19/VPS/V_PBE19.vps"
_VPS = Path(os.environ.get("CUBE4SIESTA_VPS", _VPS_DEFAULT))


@pytest.mark.skipif(not _VPS.exists(), reason="OpenMX V_PBE19.vps not available")
def test_parse_openmx_vps():
    p = parse_openmx_vps(_VPS)
    assert p.z == 23
    assert p.symbol == "V"
    assert p.xc == "pb"
    assert p.relativistic is True
    assert p.n_valence_electrons == pytest.approx(13.0)
    # Should have channels for 3s, 3p, 3d, 4s, 4p, 4d, 4f
    assert len(p.valence) >= 4
    # 3s channel
    s3 = [v for v in p.valence if v.n == 3 and v.l == 0]
    assert len(s3) == 1
    assert s3[0].occ == pytest.approx(2.0)
    assert s3[0].rc == pytest.approx(1.10)


@pytest.mark.skipif(not _VPS.exists(), reason="OpenMX V_PBE19.vps not available")
def test_vps_to_atom_input(tmp_path):
    p = parse_openmx_vps(_VPS)
    out = tmp_path / "V.pg.inp"
    content = write_atom_input(p, out)

    lines = content.strip().splitlines()
    assert lines[0] == "pg"
    assert "tm2" in lines[1]
    assert lines[2].startswith("V")
    assert "pbr" in lines[2]  # PBE + relativistic

    # Last line should have 6 floats (4 rc + 2 core-correction)
    rc_line = lines[-1].split()
    assert len(rc_line) == 6
    assert float(rc_line[0]) == pytest.approx(1.10)  # s cutoff


# ---------------------------------------------------------------------------
# QE .UPF
# ---------------------------------------------------------------------------

_UPF_SI = Path("/home/users2/cha/programs/wannier90-3.1.0/pseudo/Si.pbe-n-van.UPF")
_UPF_TE = Path("/home/users2/cha/programs/wannier90-3.1.0/examples/example24/input/Te.pbe-n-nc.UPF")


@pytest.mark.skipif(not _UPF_SI.exists(), reason="QE Si UPF not available")
def test_parse_qe_upf_si():
    p = parse_qe_upf(_UPF_SI)
    assert p.symbol == "Si"
    assert p.z == 14
    assert p.xc == "pb"
    assert p.n_valence_electrons == pytest.approx(4.0)
    assert any(v.l == 0 for v in p.valence)  # has s channel
    assert any(v.l == 1 for v in p.valence)  # has p channel


@pytest.mark.skipif(not _UPF_TE.exists(), reason="QE Te UPF not available")
def test_parse_qe_upf_te():
    p = parse_qe_upf(_UPF_TE)
    assert p.symbol == "Te"
    assert p.z == 52
    assert p.n_valence_electrons == pytest.approx(6.0)


@pytest.mark.skipif(not _UPF_SI.exists(), reason="QE Si UPF not available")
def test_upf_to_atom_input(tmp_path):
    p = parse_qe_upf(_UPF_SI)
    out = tmp_path / "Si.pg.inp"
    content = write_atom_input(p, out)

    lines = content.strip().splitlines()
    assert lines[0] == "pg"
    assert "Si" in lines[2]

    # Verify round-trip: the generated file should have correct n_val shells
    shell_count = int(lines[4].split()[1])
    assert shell_count == len(p.valence)
