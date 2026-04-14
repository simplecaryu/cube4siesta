# SiC: VASP CHGCAR -> SIESTA Rho.Restart end-to-end

Demonstrates restarting SIESTA from a density produced by a completely
different code (VASP 6.5 HSE06 testsuite, zincblende SiC, 2 atoms,
24^3 grid).

## Prerequisites

- VASP CHGCAR at some path, e.g.
  `/home/users2/cha/programs/vasp.6.5.0/testsuite/tests/SiC_HSE06_ALGO=D_RPR/CHGCAR`
- SIESTA 4.1.5 built with the `Rho.Restart` patch (see the
  sibling `siesta-4.1.5/` repo).
- `pip install -e ..` from the cube4siesta root.

`Si.psf` and `C.psf` in this directory are copied from SIESTA's own
`Tests/Pseudos/` and included for reproducibility. Attribution: SIESTA
project (GPL).

## Run

```bash
# 1. Baseline SIESTA SCF (generates sic.RHO so cube4siesta can match grids)
mpirun -np 4 siesta < sic.fdf > sic_base.out

# 2. VASP CHGCAR -> Gaussian cube
python -c "from cube4siesta.vasp_io import chgcar_to_cube; \
    chgcar_to_cube('/path/to/CHGCAR', 'sic_vasp.cube')"

# 3. cube -> SIESTA .RHOIN (matches baseline mesh via --from-siesta-rho)
cube4siesta convert --cube sic_vasp.cube --output sic.RHOIN \
    --from-siesta-rho sic.RHO --verify

# 4. SIESTA Rho.Restart: one-shot rho -> H -> DM
mpirun -np 4 siesta < sic_rr.fdf > sic_rr.out
```

## Observations

- The SIESTA mesh (24,24,24 for `MeshCutoff=200 Ry` on this FCC cell)
  happens to match VASP's grid, so resampling is skipped.
- VASP's CHGCAR integrates to ~8 electrons (Si 4 + C 4 valence).
- After Rho.Restart one-step, Fermi level sits ~0.07 eV off SIESTA's
  converged baseline; total energy ~0.5 eV off. Differences this size
  are expected because VASP (PAW + plane waves) and SIESTA (NC psf + DZP
  LCAO) differ in their core–valence partitioning and basis completeness.
  The value is **not** a check of physical agreement — it confirms the
  restart mechanism delivers a sane starting point for a subsequent SCF.

