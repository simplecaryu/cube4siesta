"""
MoS2 (Mo 4s4p semicore -> strong pp mismatch). Two questions:

(1) SEMICORE-FILTER test: is the part of OpenMX dden that the SIESTA basis
    cannot represent (the fit residual) the unwanted cross-pp semicore-change,
    or lost valence bonding? Reference = SIESTA's own valence rearrangement
    dden_SIESTA = A @ (x_true - D_atom). If basis-filtered dden is CLOSER to
    dden_SIESTA than raw OpenMX dden, projection is a beneficial filter.

(2) DM numbers for the comparison table: D_atom alone, D_atom + fit(dden).
"""
import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from rho_fit import Geometry, build_A, build_D_atomic, fit
sys.path.insert(0, "src")
from cube4siesta.rho_io import read_rho

B = Path("/home/users2/cha/work/edison/c4s-2x2/talk_dm/mos2_nnin")
SB = B / "siesta_base"
geo = Geometry(SB, stem="mos2")
rho = read_rho(SB / "mos2.RHO")
cell, mesh = rho.cell, rho.mesh
nx, ny, nz = mesh
print(f"MoS2: no_u={geo.dm.no_u} nnz={geo.dm.nnz} nsc={geo.dm.nsc} mesh={mesh}", flush=True)

A, dV = build_A(geo, cell, mesh)
D_atom = build_D_atomic(geo, SB)
xt = geo.x_true
print(f"built A nnz={A.nnz};  Ne(atomic from D_atom)={(A@D_atom).sum()*dV:.3f}", flush=True)

# OpenMX dden, already on the SIESTA grid (the density rrd actually injected)
rho_dd = read_rho(B / "siesta_rrd" / "mos2.RHOIN.diff").rho[..., 0].ravel()

# SIESTA's own valence rearrangement (reference truth)
dden_sie = A @ (xt - D_atom)

# grid coords + min-image distance to Mo (atom 1 -> row 0)
fi, fj, fk = np.arange(nx)/nx, np.arange(ny)/ny, np.arange(nz)/nz
FI, FJ, FK = np.meshgrid(fi, fj, fk, indexing="ij")
gcart = np.stack([FI.ravel(), FJ.ravel(), FK.ravel()], 1) @ cell
Mo = (geo.frac @ cell)[0]
dmin = np.full(gcart.shape[0], 1e9)
for n1 in (-1, 0, 1):
    for n2 in (-1, 0, 1):
        dmin = np.minimum(dmin, np.linalg.norm(gcart - (Mo + n1*cell[0] + n2*cell[1]), axis=1))

relerr = lambda a, b: np.linalg.norm(a - b) / np.linalg.norm(b)
def near(field, r):
    return np.linalg.norm(field[dmin < r]) / np.linalg.norm(field)

ntru = np.linalg.norm(xt)
print(f"\n[DM] D_atom alone vs converged: L2err={np.linalg.norm(D_atom-xt)/ntru:.4f}"
      f"  corr={np.corrcoef(D_atom,xt)[0,1]:+.4f}", flush=True)

for damp in (0.03, 0.1):
    x_delta = fit(A, rho_dd, damp)
    filt = A @ x_delta
    resid = rho_dd - filt
    Dn = D_atom + x_delta
    print(f"\n--- damp={damp} ---", flush=True)
    print(f"  [DM] D_atom + fit(dden): L2err={np.linalg.norm(Dn-xt)/ntru:.4f}"
          f"  corr={np.corrcoef(Dn,xt)[0,1]:+.4f}  max|D|={np.abs(x_delta).max():.2f}")
    print(f"  [filter] E_raw (raw OMX dden vs SIESTA dden) = {relerr(rho_dd, dden_sie):.4f}")
    print(f"  [filter] E_filt(filtered    vs SIESTA dden) = {relerr(filt,   dden_sie):.4f}"
          f"   -> {'HELPS' if relerr(filt,dden_sie)<relerr(rho_dd,dden_sie) else 'HURTS'}")
    for r in (0.8, 1.2):
        print(f"  [spatial r<{r}Bohr of Mo] raw_dden={near(rho_dd,r):.3f}  "
              f"SIESTA_dden={near(dden_sie,r):.3f}  discarded_resid={near(resid,r):.3f}")
