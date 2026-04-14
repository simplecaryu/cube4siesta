# VSSe: OpenMX .tden.cube -> SIESTA Rho.Restart (pseudo mismatch study)

Demonstrates the pipeline on a real 2D Janus TMD density (VSSe 1T, 3
atoms, 54x54x343 cube produced by OpenMX 3.9 with PBE + spin polarization).

## Source data

OpenMX run in `/home/users2/cha/work/jx_spirit/tutorial/VSSe_1T/`,
total density `VSSe.tden.cube` (13 MB).

## Key observation: this is a **pseudopotential-mismatch** demo

The mechanism runs to completion, but the output is physically
degraded. Comparing:

| Quantity | Baseline SIESTA | VSSe rho-restart |
|----------|-----------------|-------------------|
| E_KS (eV) | -696.6 | -540.9 |
| Fermi (eV) | -4.73 | +30.4 |
| N_e       | 17 (V 5 + S 6 + Se 6) | 25 raw / 17 after rescale |

The 150 eV energy offset is **not** a bug in Rho.Restart. It is the
direct consequence of OpenMX and the Cornell NNIN V.psf using different
valence configurations:

- OpenMX's V_PBE19 treats V 3s / 3p as *valence* (13 valence electrons
  in total), so the cube density has a pronounced peak in the semicore
  region around V.
- Cornell NNIN `V-gga.psf` freezes V 3s / 3p into the core (5 valence
  electrons), so SIESTA's Hamiltonian has no states to occupy at those
  energies.

Rescaling the integral from 25 to 17 electrons with cube4siesta's
`--rescale-to 17.0` is a uniform multiplicative fix that does **not**
remove the spatially localized semicore peak — it just squashes it.
SIESTA's diagonalization then tries to fit 17 valence states to a
density that is structurally a 25-electron distribution, and the
one-shot Fermi level lands somewhere absurd.

## Conclusion

Cross-code rho restart requires pseudo-level compatibility. For a
physically useful cross-code pipeline one needs **matching core/valence
partitioning** on both sides — either use a V pseudo with the same
semicore treatment on the SIESTA side, or regenerate the OpenMX
calculation with a frozen-core V pseudopotential.

For SiC (`examples/sic_vasp/`) the VASP PAW and Cornell NNIN
`C.psf`/`Si.psf` happen to agree closely on core/valence, so that case
lands within 0.5 eV of a fully converged SIESTA SCF.

## Files

- `vsse.fdf`     - SIESTA baseline SCF input (converges in ~32 steps
                    with MixingWeight=0.03, Pulay=6, MaxSCF=300).
- `vsse_rr.fdf`  - same plus `Rho.Restart true` / `Rho.RestartFile
                    vsse.RHOIN`.
- `V.psf`, `S.psf`, `Se.psf` from Cornell NNIN Virtual Vault,
  GGA-PBE (redistributed under their public terms).

## Reproducing

```bash
# baseline
mpirun -np 4 siesta < vsse.fdf > vsse_base.out

# cube -> RHOIN (resamples 54^2 -> 24^2 and rescales N_e 25 -> 17)
cube4siesta convert \
  --cube /path/to/VSSe.tden.cube \
  --output vsse.RHOIN \
  --from-siesta-rho vsse.RHO \
  --rescale-to 17.0 --verify

# rho-restart
mpirun -np 4 siesta < vsse_rr.fdf > vsse_rr.out
```
