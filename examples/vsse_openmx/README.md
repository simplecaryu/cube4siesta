# VSSe: OpenMX .tden.cube -> SIESTA Rho.Restart (pseudo mismatch study)

Demonstrates the pipeline on a real 2D Janus TMD density (VSSe 1T, 3
atoms, 54x54x343 cube produced by OpenMX 3.9 with PBE + spin polarization).

## Source data

OpenMX run in `/home/users2/cha/work/jx_spirit/tutorial/VSSe_1T/`,
total density `VSSe.tden.cube` (13 MB).

## Key observation: the two codes are transferring *different* ρ(r)

The mechanism runs to completion, but the output is physically
degraded. Comparing:

| Quantity | Baseline SIESTA | VSSe rho-restart |
|----------|-----------------|-------------------|
| E_KS (eV) | -696.6 | -540.9 |
| Fermi (eV) | -4.73 | +30.4 |
| ∫ρ dV (N_e) | 17 (V 5 + S 6 + Se 6) | 25 raw / 17 after rescale |

A natural first reaction: "we're transferring real-space ρ(r), what
does pseudo valence *count* have to do with it?" — and that is right
in the abstract. The subtlety is that *ρ(r) itself is defined by
what each code considers valence*. ρ is a scalar field, but a
scalar field of what? "Valence electron density". Change the
partition and you change the field.

- **Not a format issue.** The Cornell NNIN `V.psf` parses fine; the
  SIESTA baseline run converges cleanly (E_KS = -696 eV, 32 SCF
  steps). The psf header format matches SIESTA's own pseudos in
  `Tests/Pseudos/` exactly (`ATM ... Troullier-Martins`).
- **A physics issue.** OpenMX's V_PBE19 pseudo treats V 3s/3p as
  *valence*, so the OpenMX cube contains a pronounced 8-electron
  peak in the semicore region around V (the same 3s/3p charge
  distribution that a full-core all-electron calculation would also
  contain). Cornell NNIN `V-gga.psf` freezes 3s/3p into the ionic
  core, so SIESTA's ρ_valence(r) for the same physical atom does
  NOT contain that peak. These are literally different scalar fields.

Rescaling the integral from 25 to 17 electrons with cube4siesta's
`--rescale-to 17.0` is a uniform multiplicative fix that does **not**
remove the spatially localized semicore peak — it just squashes it.
SIESTA's diagonalization then tries to fit 17 valence states to a
density that is structurally a 25-electron distribution, and the
one-shot Fermi level lands somewhere absurd.

## Conclusion

Cross-code ρ restart is physically meaningful only when the source
and target codes call the *same set of electrons* valence. For VSSe
with these two pseudos they don't, so the scalar field ρ_openmx(r)
and the scalar field ρ_siesta(r) that the SCF expects are
inherently different quantities — no grid resampling or
renormalization turns one into the other.

Paths to make this physical:
1. Regenerate the SIESTA-side V pseudopotential (e.g. via ATOM)
   with the same configuration as OpenMX's V_PBE19 (3s, 3p as
   valence); or
2. Rerun the OpenMX calculation with a frozen-core V pseudo that
   matches Cornell NNIN's 5-valence choice.

For SiC (`examples/sic_vasp/`) the VASP PAW and Cornell NNIN
`C.psf`/`Si.psf` happen to agree on core/valence, so that case
lands within ~0.5 eV of a fully converged SIESTA SCF — close
enough to serve as a seed for continued SCF.

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

## Update: difference-density mode resolves the mismatch

Pipeline on the diff-density branch (`feat/diff-density-restart` of
both repos):

```bash
# Convert OpenMX's already-computed ρ - ρ_atomic  (VSSe.dden.cube)
cube4siesta convert \
    --cube /path/to/VSSe.dden.cube \
    --output vsse.RHOIN.diff \
    --from-siesta-rho vsse.RHO \
    --diff --verify

# Inject it and ask SIESTA to add its own rhoatm back in
# (add  Rho.Restart.Diff true  to vsse.fdf; file provided as vsse_rr_diff.fdf)
mpirun -np 4 siesta < vsse_rr_diff.fdf > vsse_rr_diff.out
```

Results summary:

| Quantity | Baseline SCF | Total-ρ restart | **Diff-ρ restart** |
|----------|-------------|-----------------|----------------------|
| E_KS (eV) | -696.6 | -540.9 | **-693.1** |
| Fermi (eV) | -4.73 | +30.4 | **-3.37** |

The 150-eV gap closes to ~3 eV, and Fermi returns to the physically
sensible range. See `docs/issues/001-cross-pseudo-diff-density.md`
for the full analysis.
