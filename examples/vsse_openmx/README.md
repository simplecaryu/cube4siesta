# VSSe: OpenMX .tden.cube -> SIESTA Rho.Restart (pseudo mismatch study)

Demonstrates the pipeline on a real 2D Janus TMD density (VSSe 1T, 3
atoms, 54x54x343 cube produced by OpenMX 3.9 with PBE + spin polarization).

## Source data

OpenMX 3.9 run of the same structure (PBE, V_PBE19/S_PBE19/Se_PBE19
pseudopotentials), total density `VSSe.tden.cube` (13 MB) and difference
density `VSSe.dden.cube` — OpenMX writes both per default.

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

Instead of the total density, transfer only the bonding-induced
redistribution Δρ = ρ_scf − ρ_atomic (OpenMX already writes it as
`*.dden.cube`) and let SIESTA add its own atomic superposition back:

```bash
# Convert OpenMX's difference density (no rescaling — Δρ integrates to ~0)
cube4siesta convert \
    --cube /path/to/VSSe.dden.cube \
    --output vsse.RHOIN.diff \
    --from-siesta-rho vsse.RHO \
    --diff --verify

# vsse_rr_diff.fdf = vsse.fdf + Rho.Restart + Rho.Restart.Diff
mpirun -np 4 siesta < vsse_rr_diff.fdf > vsse_rr_diff.out
```

Results vs the converged SIESTA baseline:

| Quantity | Baseline SCF | Total-ρ restart | **Diff-ρ restart** |
|----------|-------------|-----------------|----------------------|
| E_KS (eV) | −696.6 | −540.9 | **−693.1** |
| Fermi (eV) | −4.73 | +30.4 | **−3.37** |
| DM rel. Frobenius vs baseline | — | 1.49 | **0.103** |
| mean \|Δband\| occupied (eV) | — | 26.9 | **0.21** |

The 150-eV energy gap closes to ~3 eV, the bands match to 0.21 eV, and
the one-shot density matrix lands within 10 % of the converged one —
all with the standard 5-valence `V.psf`, because the V semicore lives
in ρ_atomic on both sides and cancels out of Δρ.

Caveat: this only helps when there *is* a valence-partition mismatch
and the source provides a real norm-conserving difference density.
When the partitions already match, plain total-ρ restart is more
accurate (SiC: see `examples/sic_vasp/diff/`). And the diff route does
not work from PAW sources (DM-verified: VSSe-PAW Δρ gives relFrob 1.41,
no better than total-ρ's 1.49, where the NC source reaches 0.103) —
for a PAW source with mismatched valence, regenerate the SIESTA
pseudopotential instead.
