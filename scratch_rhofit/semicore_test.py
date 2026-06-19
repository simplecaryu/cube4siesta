"""
Test the hypothesis: when fitting OpenMX dden onto the SIESTA basis, the part the
fit CANNOT represent (the residual) is the unwanted cross-pp semicore-change,
NOT lost valence bonding. If true, basis projection is a beneficial filter and
the filtered dden should be CLOSER to SIESTA's own valence rearrangement than the
raw OpenMX dden.

Reference 'truth' for the valence rearrangement in SIESTA's world:
    dden_SIESTA = A @ (x_true - D_atom)
    (density of converged crystal DM minus neutral-atom DM, in the SIESTA basis)
"""
import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from rho_fit import Geometry, build_A, build_D_atomic, resample_corner, fit
sys.path.insert(0, "src")
from cube4siesta.rho_io import read_rho
from cube4siesta.cube_io import read_cube

base = Path("examples/vsse_openmx/vsse_baseline")
geo = Geometry(base)
rho = read_rho(base / "vsse.RHO")
cell, mesh = rho.cell, rho.mesh
nx, ny, nz = mesh
A, dV = build_A(geo, cell, mesh)
D_atom = build_D_atomic(geo, base)
xt = geo.x_true

# OpenMX difference density on the SIESTA grid
dden = read_cube("examples/vsse_openmx/vsse_scfout_gen/VSSe.dden.cube")
rho_dd = resample_corner(dden.data, mesh).ravel()

# SIESTA's own valence rearrangement (the reference truth)
dden_sie = A @ (xt - D_atom)

# grid Cartesian coords (same construction as build_A)
fi, fj, fk = np.arange(nx)/nx, np.arange(ny)/ny, np.arange(nz)/nz
FI, FJ, FK = np.meshgrid(fi, fj, fk, indexing="ij")
gcart = np.stack([FI.ravel(), FJ.ravel(), FK.ravel()], 1) @ cell

# min-image distance of every grid point to the V atom (atom index 1 -> row 0)
atom_cart = geo.frac @ cell
V = atom_cart[0]
dmin = np.full(gcart.shape[0], 1e9)
for n1 in (-1, 0, 1):
    for n2 in (-1, 0, 1):
        c = V + n1*cell[0] + n2*cell[1]
        d = np.linalg.norm(gcart - c, axis=1)
        dmin = np.minimum(dmin, d)

def relerr(a, b):
    return np.linalg.norm(a - b) / np.linalg.norm(b)

def near_frac(field, r):
    """fraction of the field's L2 weight within radius r of V."""
    m = dmin < r
    return np.linalg.norm(field[m]) / np.linalg.norm(field)

print(f"V at {np.round(V,2)} Bohr;  grid {mesh}")
for damp in (0.03, 0.1):
    x_delta = fit(A, rho_dd, damp)
    filt = A @ x_delta          # basis-filtered OpenMX dden
    resid = rho_dd - filt       # what the fit discarded
    print(f"\n--- damp={damp} ---")
    print(f"  E_raw  (raw OMX dden   vs SIESTA dden) = {relerr(rho_dd, dden_sie):.4f}")
    print(f"  E_filt (filtered dden  vs SIESTA dden) = {relerr(filt,   dden_sie):.4f}")
    print(f"  -> filtering {'HELPS' if relerr(filt,dden_sie)<relerr(rho_dd,dden_sie) else 'HURTS'}")
    for r in (0.8, 1.2):
        print(f"  L2 weight within r<{r} Bohr of V:  "
              f"raw_dden={near_frac(rho_dd,r):.3f}  "
              f"SIESTA_dden={near_frac(dden_sie,r):.3f}  "
              f"discarded_residual={near_frac(resid,r):.3f}")
