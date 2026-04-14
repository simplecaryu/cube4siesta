# cube4siesta

Restart a SIESTA calculation from a real-space charge density produced
by another DFT code (VASP, OpenMX, Quantum ESPRESSO, …).

Takes a Gaussian cube file of ρ(r) in, emits a SIESTA `.RHO`-format
file out. Pairs with a small SIESTA 4.1.5 patch that adds the
`Rho.Restart` fdf flag, so SIESTA runs one `ρ → V_KS → H → DM` step
from your external density instead of the usual atomic-superposition
initial guess.

Key features:

- **Gaussian cube** reader/writer and **SIESTA `.RHO`** binary I/O.
- **VASP `CHGCAR` → cube** converter built on pymatgen.
- **3D resampling** onto SIESTA's FFT mesh (trilinear or cubic).
- **Difference-density mode** (`--diff` + `Rho.Restart.Diff`) that
  makes the pipeline robust against cross-code pseudopotential
  mismatches — transfer only the bonding-induced ρ − ρ_atomic and let
  SIESTA regenerate the pseudo-specific atomic part itself.

---

## Why this exists

SIESTA's native SCF restart file is `.DM` (the density matrix in the
NAO basis). Other codes don't speak that format — VASP gives you
`CHGCAR`, OpenMX gives you `*.tden.cube`, etc. Converting ρ(r) back
into a valid `.DM` requires one self-consistency step's worth of
diagonalization.

This project does three things:

1. Patches SIESTA 4.1.5 to accept a real-space density at the top of
   its SCF loop: `ρ → V_H[ρ] + V_xc[ρ] → H_μν → diag → DM`, one shot.
2. Provides a Python package to convert external cube / CHGCAR files
   into the grid layout SIESTA expects.
3. Handles the subtle point that different codes' pseudopotentials
   define *different* "valence" densities, via a difference-density
   transfer path.

---

## Installation

### 1. SIESTA with the Rho.Restart patch

Clone or unpack SIESTA 4.1.5 (sibling directory `siesta-4.1.5/`) and
apply the patch by checking out the feature branch of the SIESTA side
of this project:

```bash
cd siesta-4.1.5
git checkout feat/diff-density-restart   # or main for the first version
cd Obj
sh ../Src/obj_setup.sh
cp intel.make arch.make           # or gfortran.make, edit as needed
make -j 8 siesta
```

Build tested with Intel oneAPI 2022 (`ifort`, Intel MPI, MKL). See
the SIESTA 4.1.5 manual for arch.make templates on other systems.

### 2. cube4siesta Python package

```bash
pip install -e .
```

Dependencies: `numpy`, `scipy`, and `pymatgen` (only if you use the
VASP CHGCAR converter).

Installs two console scripts:

- `cube4siesta` — main converter
- `cube4siesta-vasp-diff` — builds Δρ cubes from VASP CHGCAR pairs

---

## How injection actually works

Every run of `Rho.Restart` in SIESTA ends up calling

$$\rho_\text{inject}(r) \;=\; \rho_\text{file}(r) \;+\; [\text{optional}]\;\rho_\text{atom}^\text{SIESTA}(r),$$

where the optional atomic term is added whenever you set
`Rho.Restart.Diff true`. Everything else — the numerical quality of the
seed DM that comes out — depends only on what you put into
`ρ_file(r)`. There are three useful choices, and they are strictly
ordered by robustness:

### Choice 1 — total density (simplest)

```
ρ_file := ρ_scf^src
```

Plain `cube4siesta convert ...` with no `--diff`, and no
`Rho.Restart.Diff` on the SIESTA side. The SCF density from the source
code is injected directly; SIESTA treats it as its own valence density.

Failure mode: if the source and SIESTA disagree on *which electrons are
valence* (e.g. V with semicore vs V with frozen 3s/3p), the scalar
field doesn't match SIESTA's Hamiltonian at all.

### Choice 2 — diff against SIESTA's own atomic density (**recommended**)

```
ρ_file := ρ_scf^src − ρ_atom^SIESTA       # subtracted on the cube4siesta side
ρ_inject := ρ_file + ρ_atom^SIESTA         # added back by SIESTA
         := ρ_scf^src                      # identical to Choice 1 by construction
```

`cube4siesta convert --diff --subtract <SIESTA rhoatm cube>` with
`Rho.Restart.Diff true` on the SIESTA side. When the pseudos agree on
valence this is numerically identical to Choice 1 — the patch is a
true **no-op relative to total-ρ**. It is strictly no worse than
Choice 1 in any case.

### Choice 3 — diff against the *source*'s own atomic density (specialized)

```
ρ_file := ρ_scf^src − ρ_atom^src           # source-side subtraction (e.g. OpenMX *.dden.cube)
ρ_inject := ρ_file + ρ_atom^SIESTA
         := ρ_scf^src + (ρ_atom^SIESTA − ρ_atom^src)
```

`cube4siesta convert --diff` with `Rho.Restart.Diff true`. The residual
`(ρ_atom^SIESTA − ρ_atom^src)` lives near the nuclei and is nonzero
whenever the pseudos disagree — either on valence partition (the win
case) or on pseudized-atomic shape (the loss case, e.g. PAW vs NC).

Useful **only** when the source code has semicore valence that
SIESTA's pseudo freezes, and both sides are in the same pseudo family
(e.g. OpenMX NC with V 3s/3p valence → SIESTA NC with V 3s/3p frozen).
In that specific regime it selectively removes the semicore peak that
SIESTA can't describe; outside that regime it just injects a local
shape mismatch.

### Measured

| | Baseline SCF | Choice 1 (total) | Choice 2 (diff vs SIESTA) | Choice 3 (diff vs source) |
|---|---|---|---|---|
| SiC (VASP PAW, 4+4 val, matches SIESTA) | −263.98 | **−264.54** | **−264.54** | −284.06 |
| VSSe (OpenMX 13-val V → SIESTA 5-val V) | −696.61 | −540.92 | −540.92 | **−693.13** |

Choice 2 numerically tracks Choice 1 exactly, which is the right answer
in the matching case (SiC) and still terrible in the mismatch case
(VSSe) — so Choice 2 is safe but can't work around a semicore
disagreement. Only Choice 3 fixes VSSe, but only when the rest of the
pseudo setup is similar enough.

### Rule of thumb

| Situation | Recommendation |
|-----------|----------------|
| You trust the two pseudos agree on valence | **Choice 1 or 2** (pick whichever is easier — same numerical answer) |
| Source has semicore valence that SIESTA freezes, both sides NC | **Choice 3** |
| Source is PAW (VASP) | **Choice 1 or 2** — never Choice 3 |

## Quickstart — the universal recipe (works for any source code)

This is Choice 1 or 2 from the previous section; the two are
numerically identical when the pseudos agree on valence (the common
case), so pick whichever is easier to script.

```bash
# Step 1. Source code produces a cube of rho_scf(r) in e/Bohr^3.
#   OpenMX:  *.tden.cube is written by default.
#   VASP:    python -c "from cube4siesta.vasp_io import chgcar_to_cube; \
#                        chgcar_to_cube('CHGCAR','scf.cube')"
#   QE / ...: use their own cube output option or ase.io.cube.

# Step 2. Run SIESTA once to fix the target mesh.
#   (Any modest SCF with your.fdf will do; we just need the .RHO
#    for its cell+mesh metadata.)
mpirun -np 4 siesta < your.fdf > base.out       # writes your.RHO

# Step 3a. Choice 1 — feed the total density directly
cube4siesta convert --cube scf.cube \
    --output your.RHOIN \
    --from-siesta-rho your.RHO --verify
# fdf: Rho.Restart true, Rho.RestartFile your.RHOIN  (no Rho.Restart.Diff)

# Step 3b. Choice 2 — SIESTA-consistent diff (a no-op vs Choice 1 when
#          pseudos match, cheap insurance otherwise)
#   First extract SIESTA's own rhoatm on the same mesh:
cat your.fdf analyze.fdf > run_rhoatm.fdf   # with  AnalyzeChargeDensityOnly  true
mpirun -np 4 siesta < run_rhoatm.fdf > rhoatm.out   # writes your.RHO (= rhoatm)

cube4siesta convert --cube scf.cube \
    --output your.RHOIN.diff \
    --from-siesta-rho your.RHO \
    --diff --subtract your.RHO --verify
# fdf: Rho.Restart true, Rho.Restart.Diff true, Rho.RestartFile your.RHOIN.diff
```

See `examples/sic_vasp/` for a worked end-to-end run on zincblende
SiC (Choice 1 and Choice 2 give identical numbers there).

---

## When you need Choice 3: semicore mismatch

Use this path only when you know that your source code treats more
electrons as valence than your SIESTA `.psf` does (e.g. OpenMX
`V_PBE19` with 13 valence electrons going into a run whose Cornell NNIN
`V.psf` has only 5). In that specific case the source's own atomic
density carries the semicore peaks as part of the "atomic" reference;
subtracting it cancels those peaks out of the cube, and SIESTA's own
valence-only `rhoatm` is the correct field to re-add.

For OpenMX this requires only a single extra file — `*.dden.cube` is
written automatically:

```bash
cube4siesta convert --cube /path/to/System.dden.cube \
    --output your.RHOIN.diff \
    --from-siesta-rho your.RHO --diff --verify
# fdf: Rho.Restart true, Rho.Restart.Diff true, Rho.RestartFile your.RHOIN.diff
```

For VASP (Choice 3 almost never helps because VASP is PAW — prefer
Choice 1 or 2) you would build the analogue with two VASP runs and
`cube4siesta-vasp-diff`; see `examples/sic_vasp/diff/README.md` for an
instrumented demo that also illustrates why this is a bad idea when
the frameworks differ.

---

## The SIESTA side: fdf flags

| Flag | Default | Meaning |
|------|---------|---------|
| `Rho.Restart`       | `false` | Enable the ρ → one-step DM path. Forces `MaxSCFIterations=1`, disables DM.UseSaveDM and SCF convergence checks. |
| `Rho.RestartFile`   | `<SystemLabel>.RHOIN` | Path to the input SIESTA-format density file. |
| `Rho.Restart.Diff`  | `false` | Treat the loaded field as Δρ = ρ − ρ_atomic. SIESTA adds its own rhoatm on top. Requires `Rho.Restart=true`. |

The output after a `Rho.Restart` run is a standard `.DM` / `.EIG` /
`.FA` / `.XV` set — a clean starting point for a normal continuation
SCF (just remove the `Rho.Restart*` lines and add `DM.UseSaveDM true`).

---

## Pseudopotentials

SIESTA 4.1.5 reads `.psf` (ATOM / Froyen format). PSML is only
supported from SIESTA 5.x onward.

Where to get .psf files:

- **Cornell NNIN Virtual Vault** —
  <https://nninc.cnf.cornell.edu/psplist.php?element=SYMBOL>.
  Element pages have direct `.psf` download links (LDA and GGA PBE
  variants). This is what `testdata/pseudos/` bundles.
- **SIESTA `Tests/Pseudos/`** — a handful of well-tested ones
  shipped with SIESTA itself.
- **SIMUNE** — commercial subscription database.
- **ATOM** — regenerate from scratch when community files don't cover
  your element or valence choice.

### Cross-code caveat

**When cross-code restart matters:** both codes must agree on which
electrons are "valence", otherwise the scalar field ρ(r) they write
to disk is literally a different quantity. E.g. OpenMX's `V_PBE19` is
a 13-valence pseudo (includes 3s/3p) while Cornell NNIN's `V-gga.psf`
is 5-valence (freezes 3s/3p). Use the difference-density path
(`--diff` + `Rho.Restart.Diff`) whenever you are not sure the two
pseudos agree.

See `docs/issues/001-cross-pseudo-diff-density.md` for the full
discussion and a numerical example.

---

## CLI reference

### `cube4siesta convert`

```
cube4siesta convert --cube IN.cube --output OUT.RHOIN [options]

  --cube              Gaussian cube file (required)
  --output            SIESTA .RHO / .RHOIN path (required)

  --target-mesh Nx,Ny,Nz     target grid (default: cube's own grid)
  --from-siesta-rho REF.RHO  borrow mesh + cell from an existing SIESTA .RHO
                             (mutually exclusive with --target-mesh)
  --order 1|3                interpolation order (1 = trilinear, default)

  --diff                     treat input as a difference density
                             (skip --rescale-to, keep sign)
  --rescale-to N             rescale so ∫ρ dV = N electrons
                             (useful in total-ρ mode when charge is
                             close but not exact)
  --verify                   read back and confirm roundtrip
```

### `cube4siesta-vasp-diff`

```
cube4siesta-vasp-diff --scf CHGCAR.scf --atomic CHGCAR.atomic --out dden.cube

  --scf             CHGCAR from a converged VASP SCF
  --atomic          CHGCAR from  ICHARG=1  NELM=0  (atomic superposition)
  --out             output cube file
  --channel         "total" (default) or "diff" for spin-polarized CHGCAR
```

### Python API

```python
from cube4siesta import read_rho, write_rho
from cube4siesta.cube_io import read_cube, write_cube
from cube4siesta.resample import resample_to_mesh
from cube4siesta.vasp_io import chgcar_to_cube
from cube4siesta.vasp_diff import vasp_diff_chgcar
```

---

## Repository layout

```
cube4siesta/
├── pyproject.toml
├── README.md
├── src/cube4siesta/
│   ├── rho_io.py              SIESTA .RHO binary reader/writer
│   ├── cube_io.py             Gaussian cube reader/writer
│   ├── resample.py            3D trilinear / cubic resampling
│   ├── vasp_io.py             CHGCAR → cube (unit-corrected)
│   ├── vasp_diff.py           CHGCAR pair → Δρ cube
│   └── cli.py                 `cube4siesta convert` entry point
├── examples/
│   ├── sic_vasp/              VASP CHGCAR → SIESTA Rho.Restart (SiC)
│   └── vsse_openmx/           OpenMX dden.cube → Rho.Restart.Diff (VSSe)
├── testdata/pseudos/          GGA-PBE .psf for C, Si, V, S, Se, Mn, Te
├── docs/issues/
│   └── 001-cross-pseudo-diff-density.md
└── tests/                     pytest suite
```

---

## Status

- [x] .RHO I/O with real-SIESTA-file roundtrip tests
- [x] Gaussian cube I/O
- [x] Trilinear resampling
- [x] VASP CHGCAR → cube (unit-corrected vs pymatgen's `to_cube`)
- [x] `Rho.Restart` one-shot SIESTA patch
- [x] `Rho.Restart.Diff` difference-density mode
- [x] OpenMX VSSe end-to-end validated
- [x] VASP SiC end-to-end validated
- [ ] Spin-polarized (nspin=2) path
- [ ] Non-collinear (nspin=4) path
- [ ] Graceful mesh mismatch with in-SIESTA interpolation

---

## License

The Python code in `src/cube4siesta/` is released under the MIT
License. Redistributed pseudopotentials in `testdata/pseudos/` are
from the Cornell NNIN Virtual Vault and carry their own terms. The
SIESTA patch is derivative of SIESTA 4.1.5 and is distributed under
SIESTA's GPL.
