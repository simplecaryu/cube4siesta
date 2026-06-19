"""
Generate a basis-FILTERED dden as a SIESTA RHOIN.diff.
  argv: base_dir  stem  dden_path  out_path  [damp]
filtered dden = A @ fit(raw dden), on the SIESTA grid of base_dir/<stem>.RHO.
"""
import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from rho_fit import Geometry, build_A, build_D_atomic, fit
sys.path.insert(0, "src")
from cube4siesta.rho_io import read_rho, write_rho

base = Path(sys.argv[1]); stem = sys.argv[2]
dden_path = Path(sys.argv[3]); out = Path(sys.argv[4])
damp = float(sys.argv[5]) if len(sys.argv) > 5 else 0.1

geo = Geometry(base, stem=stem)
rho = read_rho(base / f"{stem}.RHO")
cell, mesh = rho.cell, rho.mesh
nx, ny, nz = mesh
print(f"{stem}: no_u={geo.dm.no_u} nnz={geo.dm.nnz} mesh={mesh}", flush=True)
A, dV = build_A(geo, cell, mesh)
print(f"built A nnz={A.nnz}", flush=True)

raw = read_rho(dden_path).rho[..., 0].ravel()
x_delta = fit(A, raw, damp)
filt = A @ x_delta

D_atom = build_D_atomic(geo, base)
dden_sie = A @ (geo.x_true - D_atom)
rel = lambda a, b: np.linalg.norm(a - b) / np.linalg.norm(b)
print(f"damp={damp}  integral(filt)={filt.sum()*dV:+.4f}  "
      f"E_raw={rel(raw,dden_sie):.4f}  E_filt={rel(filt,dden_sie):.4f}", flush=True)

write_rho(out, cell, filt.reshape(nx, ny, nz))
print(f"wrote {out} ({out.stat().st_size} bytes)", flush=True)
