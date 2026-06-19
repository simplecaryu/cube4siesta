# `dm-seed` demo: VSSe (OpenMX → SIESTA, semicore mismatch)

A small, self-contained example for the main `dm-seed` route. VSSe is the
hard case: OpenMX treats V 3s/3p as valence (25 e), the SIESTA `V.psf`
freezes them into the core (17 e total for the cell).

Run it:

```bash
cube4siesta dm-seed \
    --baseline examples/vsse_openmx/dm_seed_demo \
    --stem vsse \
    --dden examples/vsse_openmx/dm_seed_demo/vsse_dden.RHO \
    --output vsse_seed.DM
```

Then refine `vsse_seed.DM` with a single SIESTA diagonalization using the
`DM.UseSaveDM` settings that `dm-seed` prints.

## Files

| File | Role |
|---|---|
| `vsse.DM` | SIESTA baseline density matrix — only its **sparsity pattern** is used (and it doubles as the converged reference for the pre-SCF diagnostic) |
| `vsse.ORB_INDX` | orbital ↔ supercell-column map |
| `vsse.STRUCT_OUT` | geometry (fractional coordinates) |
| `vsse.RHO` | defines the real-space grid for the basis fit |
| `{V,S,Se}.ion` | SIESTA basis radials + neutral-atom populations |
| `vsse_dden.RHO` | OpenMX difference density `Δρ = ρ_scf − ρ_atoms`, corner-resampled from `VSSe.dden.cube` onto the SIESTA grid (compact stand-in for the 13 MB cube) |

The baseline files come from a plain SIESTA run of `examples/vsse_openmx/vsse.fdf`;
`vsse_dden.RHO` was produced from the OpenMX `*.dden.cube` with
`cube4siesta.density_fit.resample_corner`.
