"""Resample an OpenMX .cube density onto a reference SIESTA .RHO grid and write
a .RHO file (corner-aligned). argv: cube_path  ref_rho_path  out_rho_path"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from rho_fit import resample_corner
sys.path.insert(0, "src")
from cube4siesta.cube_io import read_cube
from cube4siesta.rho_io import read_rho, write_rho
import numpy as np

cube = read_cube(sys.argv[1])
ref = read_rho(sys.argv[2])
out = sys.argv[3]
r = resample_corner(cube.data, ref.mesh)
dV = abs(np.linalg.det(ref.cell)) / (ref.mesh[0]*ref.mesh[1]*ref.mesh[2])
print(f"cube {cube.mesh} -> grid {ref.mesh}; integral={r.sum()*dV:+.4f} (dden~0)")
write_rho(out, ref.cell, r)
print(f"wrote {out}", flush=True)
