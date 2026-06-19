"""
Write D_start = D_atom + fit(dden) as a SIESTA .DM (baseline sparsity pattern),
to test the PATCH-FREE route: standard DM.UseSaveDM restart instead of rho
injection. SIESTA's first SCF step rebuilds rho[D_start] and diagonalizes ==
the injection+diagonalization, no rho-restart patch needed.
  argv: base_dir  stem  dden_path  out.DM  [damp]
"""
import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from rho_fit import Geometry, build_A, build_D_atomic, fit
sys.path.insert(0, "src")
from cube4siesta.rho_io import read_rho
from cube4siesta.dm_io import read_dm, write_dm, copy_pattern

base = Path(sys.argv[1]); stem = sys.argv[2]
dden_path = Path(sys.argv[3]); out = Path(sys.argv[4])
damp = float(sys.argv[5]) if len(sys.argv) > 5 else 0.1

geo = Geometry(base, stem=stem)
rho = read_rho(base / f"{stem}.RHO")
A, dV = build_A(geo, rho.cell, rho.mesh)
D_atom = build_D_atomic(geo, base)
raw = read_rho(dden_path).rho[..., 0].ravel()
x_delta = fit(A, raw, damp)
D_start = D_atom + x_delta

tmpl = read_dm(base / f"{stem}.DM")
dmf = copy_pattern(tmpl, nspin=1)
dmf.dm[:, 0] = D_start
write_dm(out, dmf)
xt = geo.x_true
print(f"D_start vs converged (direct): "
      f"L2={np.linalg.norm(D_start-xt)/np.linalg.norm(xt):.4f} "
      f"(this is the BAD direct number; SIESTA's 1 diag should fix it)")
print(f"wrote seed .DM -> {out}", flush=True)
