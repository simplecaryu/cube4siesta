"""
Integration tests against real DFT output from other codes.

These tests don't run the source codes — they exercise the cube4siesta
pipeline (parse -> convert -> write .RHO) on live OpenMX cube files and
VASP CHGCAR files.

Paths come from env vars so CI / other machines can point to their own
data; tests are skipped silently when the data is absent.

Env vars:
  CUBE4SIESTA_OPENMX_CUBE  -> path to an OpenMX *.tden.cube file
  CUBE4SIESTA_VASP_CHGCAR  -> path to a VASP CHGCAR file

Defaults for the current machine:
  OpenMX  -> /home/users2/cha/work/jx_spirit/tutorial/VSSe_1T/VSSe.tden.cube
  VASP    -> /home/users2/cha/programs/vasp.6.5.0/testsuite/tests/
             SiC_HSE06_ALGO=D_RPR/CHGCAR
"""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pytest

from cube4siesta.cube_io import read_cube
from cube4siesta.rho_io import read_rho, write_rho

_OPENMX_DEFAULT = "/home/users2/cha/work/jx_spirit/tutorial/VSSe_1T/VSSe.tden.cube"
_VASP_DEFAULT = (
    "/home/users2/cha/programs/vasp.6.5.0/testsuite/tests/"
    "SiC_HSE06_ALGO=D_RPR/CHGCAR"
)


def _path_from_env(var: str, default: str) -> Path | None:
    p = Path(os.environ.get(var, default))
    return p if p.exists() else None


# ---------------------------------------------------------------------------
# OpenMX: VSSe (Janus TMD monolayer, 1T)
# ---------------------------------------------------------------------------

_OPENMX = _path_from_env("CUBE4SIESTA_OPENMX_CUBE", _OPENMX_DEFAULT)


@pytest.mark.skipif(_OPENMX is None, reason="OpenMX cube not available")
def test_openmx_vsse_cube_parse():
    """VSSe total density cube produced by OpenMX parses correctly."""
    cube = read_cube(_OPENMX)

    # 1T-VSSe is a 3-atom Janus monolayer: V (Z=23), S (Z=16), Se (Z=34)
    Zs = sorted(int(z) for z in cube.atoms[:, 0])
    assert Zs == [16, 23, 34], f"unexpected species: {Zs}"

    # hexagonal in-plane, long c-axis for a 2D slab
    assert cube.mesh[0] == cube.mesh[1], "expect square in-plane mesh for hex"
    assert cube.mesh[2] > cube.mesh[0], "expect large z-mesh for 2D slab"

    # VSSe total-valence estimate: depends on pseudo choice, but should be
    # in a physically reasonable range. OpenMX PBE pseudos for V-S-Se
    # typically give around 13-21 valence electrons.
    n_e = cube.n_electrons
    assert 10.0 < n_e < 40.0, f"N_e out of plausible range: {n_e}"


@pytest.mark.skipif(_OPENMX is None, reason="OpenMX cube not available")
def test_openmx_vsse_convert_and_roundtrip(tmp_path):
    """Convert the OpenMX cube into .RHO (same mesh) and verify charge preserves."""
    cube = read_cube(_OPENMX)
    rho_4d = cube.data.astype(np.float32)[..., np.newaxis]

    out = tmp_path / "vsse.RHOIN"
    write_rho(out, cube.cell, rho_4d)
    back = read_rho(out)

    assert back.mesh == cube.mesh
    np.testing.assert_allclose(back.cell, cube.cell, atol=1e-6)
    # real*4 truncation ~7 sig figs; charge integral should match to 5 sig figs
    assert back.n_electrons == pytest.approx(cube.n_electrons, rel=1e-4)


# ---------------------------------------------------------------------------
# VASP: SiC (FCC 2-atom via HSE testsuite)
# ---------------------------------------------------------------------------

_VASP = _path_from_env("CUBE4SIESTA_VASP_CHGCAR", _VASP_DEFAULT)


@pytest.mark.skipif(_VASP is None, reason="VASP CHGCAR not available")
def test_vasp_chgcar_to_cube_preserves_charge(tmp_path):
    """CHGCAR -> cube should give an integrable density with correct N_e."""
    pytest.importorskip("pymatgen")

    from cube4siesta.vasp_io import chgcar_to_cube

    cube_path = tmp_path / "sic.cube"
    in_mem = chgcar_to_cube(_VASP, cube_path)
    # SiC with PBE VASP pseudopotentials: Si gives 4 val. electrons, C gives 4
    assert in_mem.n_electrons == pytest.approx(8.0, abs=0.05)

    # also check reread from disk matches in-memory
    reread = read_cube(cube_path)
    assert reread.mesh == in_mem.mesh
    np.testing.assert_allclose(reread.cell, in_mem.cell, atol=1e-5)
    assert reread.n_electrons == pytest.approx(in_mem.n_electrons, rel=1e-3)


@pytest.mark.skipif(_VASP is None, reason="VASP CHGCAR not available")
def test_vasp_chgcar_full_pipeline(tmp_path):
    """VASP CHGCAR -> cube -> cube4siesta -> .RHOIN, with N_e preserved."""
    pytest.importorskip("pymatgen")

    from cube4siesta.vasp_io import chgcar_to_cube

    cube_path = tmp_path / "sic.cube"
    chgcar_to_cube(_VASP, cube_path)
    cube = read_cube(cube_path)

    rho = cube.data.astype(np.float32)[..., np.newaxis]
    rhoin = tmp_path / "sic.RHOIN"
    write_rho(rhoin, cube.cell, rho)

    back = read_rho(rhoin)
    assert back.n_electrons == pytest.approx(cube.n_electrons, rel=1e-3)
