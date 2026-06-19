"""
Prototype: fit a target real-space density rho(r) onto the SIESTA LCAO basis to
recover a density matrix D_ij(R), using the .DM sparsity pattern as the unknown
set.

    rho(r) = sum_i sum_{(j,R)} D_ij(R) * col_{i,j,R}(r)
    col_{i,j,R}(r) = sum_n phi_i(r-c_i-n) * phi_j(r-c_j-R-n)   (n = lattice images)

Two experiments, both compared to the converged SIESTA DM (vsse.DM):
  (A) round-trip : target = SIESTA's own RHO          -> tests the fit machinery
  (B) mismatch   : target = OpenMX tden cube (25 e)   -> tests pp-mismatch handling
The geometry, basis and DM sparsity pattern always come from the SIESTA baseline.
"""
from __future__ import annotations
import sys, time
from pathlib import Path
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import lsqr

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from cube4siesta.rho_io import read_rho
from cube4siesta.cube_io import read_cube
from cube4siesta.dm_io import read_dm
from cube4siesta.orb_indx import read_orb_indx
from cube4siesta.basis import read_ion_radials, build_siesta_orbitals


def read_struct_out(path):
    lines = [l.split() for l in Path(path).read_text().splitlines() if l.strip()]
    cell_ang = np.array([[float(x) for x in lines[i]] for i in range(3)])
    nat = int(lines[3][0])
    frac = np.array([[float(x) for x in lines[4 + a][2:5]] for a in range(nat)])
    return cell_ang, frac


class Geometry:
    """SIESTA baseline geometry + basis + DM pattern (the authoritative side)."""
    def __init__(self, base, stem="vsse"):
        base = Path(base)
        self.dm = read_dm(base / f"{stem}.DM")
        self.oi = read_orb_indx(base / f"{stem}.ORB_INDX")
        _, self.frac = read_struct_out(base / f"{stem}.STRUCT_OUT")
        species = sorted({o.spec for o in self.oi.unit})
        radials = {s: read_ion_radials(base / f"{s}.ion") for s in species}
        self.orbs = build_siesta_orbitals(self.oi.unit, radials)
        self.atom_of_uorb = np.array([o.ia for o in self.oi.unit])
        # true converged DM as the unknown vector (pattern order)
        self.x_true = np.concatenate(
            [self.dm.row(io)[1][:, 0] for io in range(self.dm.no_u)])


def build_A(geo, cell, mesh):
    """Sparse design matrix A (ngrid x nnz): column = density of one DM entry."""
    nx, ny, nz = mesh
    ngrid = nx * ny * nz
    atom_cart = geo.frac @ cell
    fi, fj, fk = np.arange(nx)/nx, np.arange(ny)/ny, np.arange(nz)/nz
    FI, FJ, FK = np.meshgrid(fi, fj, fk, indexing="ij")
    gcart = np.stack([FI.ravel(), FJ.ravel(), FK.ravel()], 1) @ cell

    rcut = np.array([o.rcut for o in geo.orbs])
    a_len = np.linalg.norm(cell, axis=1)
    nr = [int(np.ceil(rcut.max()/a_len[d])) + 1 for d in range(3)]

    # localized orbital-image fields g[(u,n1,n2,n3)] = (idx, val)
    g = {}
    for u in range(geo.dm.no_u):
        c0 = atom_cart[geo.atom_of_uorb[u] - 1]
        orb = geo.orbs[u]
        for n1 in range(-nr[0], nr[0]+1):
            for n2 in range(-nr[1], nr[1]+1):
                for n3 in range(-nr[2], nr[2]+1):
                    c = c0 + n1*cell[0] + n2*cell[1] + n3*cell[2]
                    d = gcart - c
                    val = orb.evaluate(d[:, 0], d[:, 1], d[:, 2])
                    idx = np.nonzero(val)[0]
                    if idx.size:
                        g[(u, n1, n2, n3)] = (idx.astype(np.int32), val[idx])
    g_by_orb = {io: [] for io in range(geo.dm.no_u)}
    for k, v in g.items():
        g_by_orb[k[0]].append((k, v))

    rows, cols, data = [], [], []
    col = 0
    for io in range(geo.dm.no_u):
        listd, _ = geo.dm.row(io)
        for jj in range(listd.size):
            juo, isc = geo.oi.decode_column(int(listd[jj]))
            j = juo - 1
            tmp = np.zeros(ngrid); touched = []
            for (u, m1, m2, m3), (idx_i, val_i) in g_by_orb[io]:
                gj = g.get((j, m1+isc[0], m2+isc[1], m3+isc[2]))
                if gj is None:
                    continue
                idx_j, val_j = gj
                common, ia_, ib_ = np.intersect1d(idx_i, idx_j, assume_unique=True,
                                                  return_indices=True)
                if common.size:
                    tmp[common] += val_i[ia_] * val_j[ib_]
                    touched.append(common)
            if touched:
                gi = np.unique(np.concatenate(touched))
                rows.append(gi); cols.append(np.full(gi.size, col, np.int32))
                data.append(tmp[gi])
            col += 1
    A = sp.csr_matrix((np.concatenate(data),
                       (np.concatenate(rows), np.concatenate(cols))),
                      shape=(ngrid, geo.dm.nnz)).tocsc()
    dV = abs(np.linalg.det(cell)) / ngrid
    return A, dV


def parse_ion_pops(path):
    """{(n,l,z): shell_population} from a SIESTA .ion PAO section."""
    out = {}; lines = Path(path).read_text().splitlines(); i = 0
    while i < len(lines) and not lines[i].lstrip().startswith("# PAOs"):
        i += 1
    for j in range(i + 1, len(lines)):
        t = lines[j].split()
        if len(t) >= 5 and all(t[k].lstrip("-").isdigit() for k in range(4)):
            l, n, z = int(t[0]), int(t[1]), int(t[2])
            try:
                out[(n, l, z)] = float(t[4])
            except ValueError:
                pass
    return out


def build_D_atomic(geo, base):
    """SIESTA neutral-atom initial DM as a vector on the .DM pattern: diagonal,
    each m-orbital of a shell gets pop(n,l,z)/(2l+1); off-diagonal/off-cell = 0."""
    species = sorted({o.spec for o in geo.oi.unit})
    pops = {s: parse_ion_pops(base / f"{s}.ion") for s in species}
    x = np.zeros(geo.dm.nnz)
    for io in range(geo.dm.no_u):
        o = geo.oi.unit[io]
        occ = pops[o.spec].get((o.n, o.l, o.z), 0.0) / (2 * o.l + 1)
        if occ == 0.0:
            continue
        listd, _ = geo.dm.row(io)
        p = int(geo.dm.listdptr[io])
        for jj in range(listd.size):
            juo, isc = geo.oi.decode_column(int(listd[jj]))
            if juo - 1 == io and isc == (0, 0, 0):     # onsite diagonal entry
                x[p + jj] = occ
                break
    return x


def resample_corner(src, target_mesh):
    """Resample a periodic density (corner-aligned: point k at frac k/N) to a
    coarser/finer mesh. i_src = i_target * N_src/N_target."""
    from scipy.ndimage import map_coordinates
    ns = src.shape
    coords = np.stack(np.meshgrid(
        *[np.arange(target_mesh[d]) * ns[d] / target_mesh[d] for d in range(3)],
        indexing="ij"), axis=0)
    return map_coordinates(src, coords, order=1, mode="grid-wrap")


def fit(A, rho_t, damp):
    return lsqr(A, rho_t, damp=damp, atol=1e-7, btol=1e-7, iter_lim=2000)[0]


def main():
    t0 = time.time()
    base = Path("examples/vsse_openmx/vsse_baseline")
    geo = Geometry(base)
    print(f"no_u={geo.dm.no_u} nnz={geo.dm.nnz} nsc={geo.dm.nsc}")

    # Build A ONCE on the SIESTA grid; reuse for both targets.
    rho = read_rho(base / "vsse.RHO")
    A, dV = build_A(geo, rho.cell, rho.mesh)
    print(f"built A (SIESTA grid {rho.mesh}) nnz={A.nnz} ({time.time()-t0:.0f}s)",
          flush=True)

    xt = geo.x_true
    nt = np.linalg.norm(xt)
    def rep(tag, x):
        err = np.linalg.norm(x - xt) / nt
        corr = np.corrcoef(x, xt)[0, 1]
        print(f"  {tag:<34s} L2err={err:.4f}  corr={corr:+.4f}  "
              f"max|D|={np.abs(x).max():.2f}", flush=True)

    # End-to-end: D_start = D_atomic + fit(Delta-rho_OpenMX), vs converged DM.
    D_atom = build_D_atomic(geo, base)
    print(f"\n[end-to-end vs converged crystal DM]  (true max|D|={np.abs(xt).max():.2f})")
    rep("D_atomic alone (neutral-atom guess)", D_atom)

    dden = read_cube("examples/vsse_openmx/vsse_scfout_gen/VSSe.dden.cube")
    rho_dd = resample_corner(dden.data, rho.mesh).ravel()
    for damp in (1e-3, 1e-2, 3e-2, 1e-1):
        x_delta = fit(A, rho_dd, damp)
        rep(f"D_atomic + fit(dDelta-rho) damp={damp:g}", D_atom + x_delta)

    # Reference: directly fitting the TOTAL OpenMX density (the broken route).
    tden = read_cube("examples/vsse_openmx/vsse_scfout_gen/VSSe.tden.cube")
    x_tot = fit(A, resample_corner(tden.data, rho.mesh).ravel(), 1e-2)
    rep("fit(total rho_OpenMX) damp=1e-2", x_tot)
    print(f"\ndone ({time.time()-t0:.0f}s)", flush=True)


if __name__ == "__main__":
    main()
