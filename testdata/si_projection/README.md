# Si fixtures for the DM-projection tests

Bulk Si (diamond, a = 5.43 Å, 2-atom primitive cell), PBE, computed with
matched settings in OpenMX 3.9 and SIESTA 4.1.5. These five files are the
minimal set needed by `tests/test_dm_projection.py` to exercise the full
OpenMX → SIESTA density-matrix projection pipeline on a real system.

| File | Producer | Role |
|------|----------|------|
| `Si.scfout`   | OpenMX 3.9 (`HS.fileout on`) | source DM + overlap (binary scfout v3) |
| `Si7.0.pao`   | OpenMX 3.9 `DFT_DATA19/PAO` (GPL, redistributed) | OpenMX radial basis tables |
| `si.DM`       | SIESTA 4.1.5 converged run | reference DM + sparsity template |
| `si.ORB_INDX` | SIESTA 4.1.5 | orbital/supercell index map |
| `Si.ion`      | SIESTA 4.1.5 (DZP, Cornell NNIN `Si.psf`) | SIESTA radial orbital tables |

The SIESTA run inputs are in `examples/si_xcode_compare` style: DZP basis,
MeshCutoff 200 Ry, 8×8×8 Monkhorst-Pack (see `examples/DM_PROJECTION_RESULTS.md`).
