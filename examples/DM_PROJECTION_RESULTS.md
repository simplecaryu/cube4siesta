# Cross-code density-matrix transfer: overlap projection vs density restart

A third cross-code method for SIESTA, alongside (1) pseudopotential regeneration
and (2) electron-count-matched density restart: **take the density matrix that
OpenMX already stores (`.scfout`) and project it directly onto SIESTA's NAO
basis**, writing a SIESTA `.DM`.

The density-restart family has two variants, both compared below:
- **total-ρ (Choice 1)**: inject the source ρ(r) verbatim, one-shot SCF.
- **Δρ (Choice 3, diff-density)**: inject ρ_scf − ρ_atomic from the source
  (OpenMX `.dden.cube`) and let SIESTA re-add its own rhoatm
  (`Rho.Restart.Diff`, patched SIESTA one-shot).

Pseudopotential regeneration (gen-atom-input → matched 13e psf) is NOT in these
tables: it changes SIESTA's basis and electron count, so its DM lives in a
different space and has no common baseline to compare against (it is
band-validated separately: mean |Δband| 0.61 eV on VSSe).

## Method (overlap projection)

A density matrix is basis-dependent, so OpenMX's `P` cannot be reused verbatim.
We project the OpenMX density operator onto SIESTA's basis, exactly on the
Born–von-Kármán k-grid implied by the `.DM` supercell size `nsc`:

```
S(k)_ab   = Σ_R <χ_a 0|χ_b R> e^{ik·R}          (SIESTA overlap)
U(k)_aμ   = Σ_C <χ_a 0|φ_μ C> e^{ik·C}          (SIESTA↔OpenMX cross overlap)
P(k)_μν   = Σ_D P_μν(D) e^{ik·D}                (OpenMX DM, spin-summed)
P_SIE(k)  = S(k)⁻¹ [U(k) P(k) U(k)†] S(k)⁻¹
P_SIE(R)  = (1/Nk) Σ_k P_SIE(k) e^{-ik·R}
```

All overlaps are real-space grid quadratures with each orbital evaluated from its
own radial table (`.ion` / `.pao`) and real harmonic. Code conventions are handled
automatically: OpenMX uses Cartesian-positive real harmonics (AngularF.c); SIESTA
uses the Condon–Shortley `(-1)^m` convention (rlylm) — verified by reproducing each
code's own overlap matrix to ~1e-5.

The method is **valence-count agnostic**: density OpenMX carries in shells SIESTA's
basis cannot represent (e.g. semicore) simply projects out, so `Tr(P_SIE·S)` returns
the SIESTA-representable electron count.

Implemented in `src/cube4siesta/{scfout_io,orb_indx,basis,overlap,dm_io,project,dm_compare}.py`,
CLI `cube4siesta project-dm` and `cube4siesta compare-dm`. Reproduce with
`examples/dm_projection_compare.sh`.

## Results — DM-value similarity vs the SIESTA self-consistent baseline

Metric: relative Frobenius error of the DM (lower = closer to SIESTA), and the
correlation of the significant (|ref|>0.02) elements.

### Si — matched valence (OpenMX 4 e == SIESTA 4 e), DZP ↔ s2p2d1

| method                | relFrob | corr(sig) |
|-----------------------|--------:|----------:|
| total-ρ restart (OpenMX)   |  0.021  |  0.9999   |
| total-ρ restart (VASP)     |  0.096  |  0.999    |
| total-ρ restart (QE)       |  0.100  |  0.998    |
| Δρ restart (OpenMX)        |  0.215  |  0.984    |
| **DM projection (OpenMX)** | **0.363** | **0.979** |
| DM projection, purified    |  0.382  |  0.969    |

When valence partitions match, total-ρ restart wins: SIESTA re-diagonalises from a
correct density, so its DM is nearly self-consistent. The Δρ correction actively
*hurts* here (0.021 → 0.215): with no mismatch to fix, swapping the source atomic
density for SIESTA's rhoatm only injects atomic-shape noise (one-shot E_KS
−198.8 eV vs baseline −214.1 eV; run in `si_xcode_compare/siesta_rr_openmx_diff/`).
This is the DM-level confirmation of the SiC-NC total-energy result. Direct DM
projection is crudest here — it measures the pure PAO↔NAO basis-representation gap
(qualitatively faithful, corr 0.98, but ~36 % in magnitude).

### VSSe — valence MISMATCH (OpenMX 25 e vs SIESTA 17 e), SZ ↔ s3p2d2f1

OpenMX's V_PBE19 treats V 3s/3p as valence (13 e); the Cornell NNIN V pseudo freezes
them into the core (5 e). Projection: `Tr(P_SIE·S) = 16.80 e` — the V semicore
(3s²3p⁶ = 8 e) projects out automatically (25 → 16.8 ≈ 17), no pre-matching needed.

| method                     | relFrob | corr(sig) |
|----------------------------|--------:|----------:|
| Δρ restart (OpenMX)        |  0.103  |  0.996    |
| **DM projection, purified**| **0.163** | **0.991** |
| DM projection (raw)        |  0.243  |  0.974    |
| total-ρ restart (OpenMX)   |  1.487  |  0.562    |

The ranking vs the matched case **inverts**. Total-ρ restart fails (relFrob 1.49,
E_KS −540.9 eV vs the baseline −696.6 eV) because it transfers a physically
different 25-e density. Both mismatch-aware methods survive: Δρ restart is closest
(0.103 — the semicore lives in ρ_atomic on both sides, so the difference density
is partition-neutral), and DM projection follows (0.243) because it is
valence-agnostic by construction.

## Canonical purification (`--purify`)

The raw projection is not idempotent: basis incompleteness and quadrature noise
leave the natural occupations of `P_SIE·S` fuzzy (Si: [−0.017, 1.18] per
channel; VSSe: [−0.10, 2.28]) and the trace short of SIESTA's qtot (16.80 vs
17). `project-dm --purify [--qtot N]` fixes both: each P(k) is diagonalised in
the S(k) metric and the natural occupations are refilled aufbau-style with one
global Fermi level (degenerate boundary states share charge equally — metals
are safe). This is the fixed point that iterative McWeeny purification
converges to, with the trace constraint built in. Effect:

- **Mismatched valence (VSSe): purification helps a lot** — relFrob
  0.243 → 0.163, and the one-shot E_KS goes from −683.6 to −694.8 eV
  (baseline −696.6). The missing 0.2 e and occupation noise were a real error.
- **Matched valence (Si): purification slightly hurts** (0.363 → 0.382) —
  the raw projection is already the least-squares optimum; forcing
  idempotency moves it away. Use raw for ML-similarity work, purified for
  seeding/physics.

## Band verification (one-shot from the projected DM)

The projected .DM can seed a stock SIESTA one-shot — `DM.UseSaveDM true`,
`MaxSCFIterations 1` + BandLines — **no Rho.Restart patch needed at all**.
Mean |Δε| vs the baseline band structure over occupied bands:

| one-shot seed | Si (occupied) | VSSe (occupied, E_F-aligned) |
|---|---:|---:|
| DM projection (raw)      | **0.015 eV** | 0.70 eV |
| DM projection, purified  | 0.022 eV | **0.22 eV** |
| Δρ restart (reference)   | — | 0.21 eV |

Si is essentially exact either way (the relFrob-0.36 DM still produces an H
indistinguishable from converged: E_KS matches the baseline to 4 decimals).
For VSSe the purified projection reproduces the baseline bands as well as the
best density-restart route (0.22 vs 0.21 eV) while remaining pure
post-processing plus one stock SIESTA run.

## Takeaway

DM-similarity ranking (vs the same SIESTA self-consistent baseline):

| | matched valence (Si) | mismatched valence (VSSe) |
|---|---|---|
| total-ρ restart | **0.021 (best)** | 1.487 (fails) |
| Δρ restart      | 0.215 | **0.103 (best)** |
| DM projection (raw)      | 0.363 | 0.243 |
| DM projection (purified) | 0.382 | 0.163 |

- Matched valence → **total-ρ restart** is the best DM-transfer route; the Δρ
  correction only adds atomic-shape noise when there is no mismatch to fix.
- Mismatched valence → **Δρ restart** gives the closest DM, but needs the patched
  SIESTA (`Rho.Restart.Diff`) plus a one-shot run, and the source must provide
  ρ − ρ_atomic (NC source; PAW atomic shapes don't transfer).
- **DM projection** never wins on raw DM similarity but is the only route that is
  (a) pure post-processing — no SIESTA run, no patch — and (b) automatically
  valence-agnostic (semicore projects out, exact in the complete-basis limit).
  Its residual is the genuine PAO↔NAO basis-representation difference, itself the
  quantity of interest for ML-style DM-similarity studies. With `--purify` it
  additionally matches the best restart route at the band level under valence
  mismatch.

## Validation performed

- SIESTA `.DM` read/write byte-identical roundtrip.
- OpenMX `.scfout` parse: `Tr(P·S)` = 8.0 (Si), 25.0 (VSSe) exact.
- Overlap engine reproduces OpenMX's own OLP to 1e-7 (Si NN), 3e-5 (VSSe V on-site,
  validates f harmonics).
- Projection: electron count conserved, Hermiticity residual 1e-11…1e-15,
  symmetry-equivalent bonds equal.
- All numbers in the tables measured with the same `cube4siesta compare-dm` tool
  against the same baselines. Δρ DMs generated with the `Rho.Restart.Diff`
  one-shot (Choice 3, OpenMX `.dden.cube` converted without `--rescale-to`):
  Si in `si_xcode_compare/siesta_rr_openmx_diff/`, VSSe in
  `vsse_openmx/vsse_rrd/`.
- Purification: synthetic unit tests in `tests/test_dm_projection.py`
  (occupations snap to {0, cap}, trace exact, fractional Fermi boundary);
  on the real systems Tr(P·S) lands on the target to 1e-9 and the one-shot
  numbers above confirm the physics.
- Band checks: one-shot runs in `si_xcode_compare/siesta_oneshot_proj{,_pure}/`
  and `vsse_openmx/vsse_oneshot_proj{,_pure}/`.
