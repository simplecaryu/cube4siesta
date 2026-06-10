#!/usr/bin/env bash
# Reproduce the OpenMX -> SIESTA density-matrix projection comparison
# (method 3: overlap projection) against the SIESTA baseline and the
# density-restart method, for Si (matched valence) and VSSe (valence mismatch).
#
# Prereqs already generated in the example dirs:
#   - OpenMX .scfout  (restart with HS.fileout on)
#   - SIESTA baseline .DM / .ORB_INDX / .ion
#   - SIESTA rho-restart .DM (patched siesta, Rho.Restart from the OpenMX cube)
#   - SIESTA drho-restart .DM (Choice 3: Rho.Restart.Diff from the OpenMX
#     .dden.cube; convert WITHOUT --rescale-to, see siesta_rr_openmx_diff/)
set -euo pipefail
cd "$(dirname "$0")/.."          # repo root

# Point OPENMX_PAO at your OpenMX installation's DFT_DATA19/PAO directory.
PAO=${OPENMX_PAO:?set OPENMX_PAO to OpenMX's DFT_DATA19/PAO directory}
PY="python -m cube4siesta.cli"   # or: python -m src.cube4siesta.cli

echo "############ Si (matched valence: OMX 4e == SIESTA 4e) ############"
S=examples/si_xcode_compare
$PY project-dm \
  --scfout   $S/openmx/scfout_gen/Si.scfout \
  --orb-indx $S/siesta_baseline/si.ORB_INDX \
  --template $S/siesta_baseline/si.DM \
  --ion Si:$S/siesta_baseline/Si.ion \
  --pao Si:$PAO/Si7.0.pao --spec Si:s2p2d1 \
  --output $S/si_omx_proj.DM

$PY project-dm \
  --scfout   $S/openmx/scfout_gen/Si.scfout \
  --orb-indx $S/siesta_baseline/si.ORB_INDX \
  --template $S/siesta_baseline/si.DM \
  --ion Si:$S/siesta_baseline/Si.ion \
  --pao Si:$PAO/Si7.0.pao --spec Si:s2p2d1 \
  --purify \
  --output $S/si_omx_proj_pure.DM

$PY compare-dm \
  --orb-indx $S/siesta_baseline/si.ORB_INDX \
  --ref  $S/siesta_baseline/si.DM \
  --cand DMprojection_OMX:$S/si_omx_proj.DM \
  --cand DMprojection_purified:$S/si_omx_proj_pure.DM \
  --cand rho-restart_OMX:$S/siesta_rr_openmx/si_rr_openmx.DM \
  --cand drho-restart_OMX:$S/siesta_rr_openmx_diff/si_rrd.DM \
  --cand rho-restart_VASP:$S/siesta_rr_vasp/si_rr_vasp.DM \
  --cand rho-restart_QE:$S/siesta_rr_qe/si_rr_qe.DM

echo
echo "############ VSSe (valence MISMATCH: OMX 25e vs SIESTA 17e) ############"
V=examples/vsse_openmx
$PY project-dm \
  --scfout   $V/vsse_scfout_gen/VSSe.scfout \
  --orb-indx $V/vsse_baseline/vsse.ORB_INDX \
  --template $V/vsse_baseline/vsse.DM \
  --ion V:$V/vsse_baseline/V.ion --ion S:$V/vsse_baseline/S.ion --ion Se:$V/vsse_baseline/Se.ion \
  --pao V:$PAO/V6.0.pao --pao S:$PAO/S7.0.pao --pao Se:$PAO/Se7.0.pao \
  --spec V:s3p2d2f1 --spec S:s3p2d2f1 --spec Se:s3p2d2f1 \
  --output $V/vsse_omx_proj.DM

$PY project-dm \
  --scfout   $V/vsse_scfout_gen/VSSe.scfout \
  --orb-indx $V/vsse_baseline/vsse.ORB_INDX \
  --template $V/vsse_baseline/vsse.DM \
  --ion V:$V/vsse_baseline/V.ion --ion S:$V/vsse_baseline/S.ion --ion Se:$V/vsse_baseline/Se.ion \
  --pao V:$PAO/V6.0.pao --pao S:$PAO/S7.0.pao --pao Se:$PAO/Se7.0.pao \
  --spec V:s3p2d2f1 --spec S:s3p2d2f1 --spec Se:s3p2d2f1 \
  --purify --qtot 17 \
  --output $V/vsse_omx_proj_pure.DM

$PY compare-dm \
  --orb-indx $V/vsse_baseline/vsse.ORB_INDX \
  --ref  $V/vsse_baseline/vsse.DM \
  --cand DMprojection_OMX:$V/vsse_omx_proj.DM \
  --cand DMprojection_purified:$V/vsse_omx_proj_pure.DM \
  --cand rho-restart_OMX:$V/vsse_rr/vsse_rr.DM \
  --cand drho-restart_OMX:$V/vsse_rrd/vsse_rrd.DM
