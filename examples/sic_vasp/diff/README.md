# SiC — difference-density restart (VASP PAW → SIESTA NC)

Companion to the main `sic_vasp/` example. Here we go one step further
and build a proper Δρ cube from two VASP runs, then use SIESTA's
`Rho.Restart.Diff` mode.

## Build the atomic-superposition CHGCAR

```bash
# In a fresh dir, same POSCAR/POTCAR/KPOINTS as the SCF run
cp ../POSCAR ../POTCAR ../KPOINTS .     # or point at the VASP testsuite dir
cp INCAR .                               # this file: ICHARG=1 NELM=0 LCHARG=.T.
mpirun -np 4 /path/to/vasp_std
cp CHGCAR CHGCAR.atomic
```

## Build the Δρ cube and convert

```bash
cube4siesta-vasp-diff \
    --scf    /path/to/VASP/scf/CHGCAR \
    --atomic CHGCAR.atomic \
    --out    sic_dden.cube

cube4siesta convert --cube sic_dden.cube --output sic.RHOIN.diff \
    --from-siesta-rho sic.RHO --diff --verify
```

## Run SIESTA

`sic_rrd.fdf` is identical to `../sic_rr.fdf` except it adds
`Rho.Restart.Diff true` and points at `sic.RHOIN.diff`.

```bash
mpirun -np 4 siesta < sic_rrd.fdf > sic_rrd.out
```

## What we observed

| Mode | E_KS (eV) | Fermi (eV) |
|---|---|---|
| Baseline SCF | −263.98 | −4.56 |
| Total-ρ restart (parent dir) | −264.54 | −4.49 |
| **Diff-ρ restart** (this dir) | **−284.06** | **−1.71** |

Contrary to VSSe, for SiC the diff-mode result is **worse** than the
straight total-ρ result. VASP PAW and the Cornell NNIN .psf give Si
and C the same valence configurations (4 electrons each), so there is
no semicore mismatch for diff mode to fix — and the residual
difference between VASP's PAW atomic density and SIESTA's NC atomic
density shows up as a ~20 eV offset. See the top-level README for the
mode-selection rule of thumb.
