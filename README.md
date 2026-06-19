# cube4siesta

**Restart a SIESTA calculation from a charge density produced by
another DFT code** — VASP, OpenMX, Quantum ESPRESSO, or anything that
can write a Gaussian cube file.

If you've ever run a DFT calculation in one code and wished you could
continue or analyze that result in SIESTA without starting from
scratch, this tool is for you.

---

## What does it do?

It turns the result of your other DFT calculation into a standard
SIESTA `.DM` (density matrix) file, which you can then use for any
further SIESTA calculation — band structures, geometry optimization,
transport, etc.

### The main route: filtered DM restart (`dm-seed`)

The recommended, general-purpose route is **`dm-seed`**. You give it a
*difference density* `Δρ = ρ_scf − ρ_atoms` (which essentially any DFT
code can write) and it fits that density onto SIESTA's own basis,
producing a seed `.DM`. SIESTA then refines the seed in a **single
diagonalization** — no SIESTA source patch, and no density matrix from
the source code.

Why this is the default:

- **Semicore-safe / mismatch-capable.** Pseudopotentials that disagree
  on how many electrons are valence (e.g. OpenMX treats V 3s/3p as
  valence, a SIESTA `V.psf` freezes them into the core) are the hard
  case. Fitting onto SIESTA's basis acts as a *filter*: any density the
  basis cannot represent — the semicore cusp among it — simply drops
  out, so it never contaminates the valence block.
- **Source-DM-free.** It needs only a density, so it works from codes
  that do not expose a density matrix at all (VASP, Quantum ESPRESSO),
  not just OpenMX.
- **Patch-free.** The seed is consumed by ordinary `DM.UseSaveDM`
  settings; you do **not** need the patched SIESTA.
- **Most consistent across systems.** Across matched and mismatched
  test systems it is the route whose accuracy varies least.

> **Caveat — simple/matched cases.** When the source and SIESTA
> pseudopotentials count the *same* valence electrons (no semicore
> mismatch), the plain **total-density restart** is simpler and can be
> a touch more accurate. Reach for `dm-seed` whenever the partitions
> differ, or when you want one route that works everywhere.

### Alternative routes (kept as additional subcommands)

- **Total / difference-density restart** (`convert` + a small SIESTA
  patch): convert the charge density to SIESTA's grid format (`.RHO`),
  then a patched SIESTA reads it. Best for the *matched* case; the
  `--diff` variant also handles mismatch but needs the patch.
- **Direct density-matrix projection** (`project-dm` / `openmx2siesta-dm`,
  OpenMX only): read the density matrix OpenMX stores and project it onto
  SIESTA's basis. No patch; most accurate when semicore is weak/absent,
  but it requires a source density matrix and degrades on strong-semicore
  atoms.

Further experimental transfer methods are developed on the `dev`
branch; the routes above are the supported ones on `main`.

---

## Which route should I use?

**If in doubt, use `dm-seed`** — it is the one route that works for
every source code and every pseudopotential pairing. The table below is
for when you want to optimize for a specific situation.

The key question is whether your source code's pseudopotential counts
the **same valence electrons** as your SIESTA `.psf` (e.g. OpenMX's
default vanadium has 13 valence electrons, a typical SIESTA `V.psf`
has 5 — that is a *mismatch*). If you're not sure, run
`cube4siesta convert` on the source cube and check whether the printed
electron count matches what SIESTA expects for your structure.

| Your situation | Use this | Patch needed? | Source DM needed? |
|---|---|---|---|
| **Anything / not sure** (default) | **filtered DM restart (`dm-seed`)** | **no** | **no** (density only) |
| Valence counts **match** (simple case) | total-density restart (Quickstart below) | yes | no |
| Valence counts **differ**, source is OpenMX / QE-with-NC-pseudos | difference-density restart (`--diff`), or `dm-seed` | yes / **no** | no |
| Valence counts **differ**, source is OpenMX, semicore weak/absent | DM projection (`project-dm --purify`) | **no** | yes (OpenMX `.scfout`) |
| Valence counts **differ**, source is VASP (PAW) | `dm-seed` (density only), or regenerate the pseudo (`gen-atom-input`) | **no** / yes | no |

Measured accuracy (relative Frobenius error of the one-step density
matrix vs a fully converged SIESTA reference): `dm-seed` lands at ~1 %
on matched Si, ~8 % on strong-semicore V-series (VSSe, VSe₂) where it
is the most accurate patch-free route, and is the most *consistent*
across the set. DM projection wins when semicore is weak (Si ~1 %,
MoS₂ ~8 %) but degrades on V-series; total-ρ is best only in the
matched case. Details and band structures in
`examples/DM_PROJECTION_RESULTS.md`.

These measurements all use **OpenMX** source densities. `dm-seed` needs
only a density, so VASP/QE sources work by construction — but the
shipped, validated numbers are OpenMX-sourced; treat other-code sources
as supported-but-unbenchmarked here.

---

## Installation

> **The main `dm-seed` route does not need the SIESTA patch** (Step 1
> below) — it only uses ordinary `DM.UseSaveDM` settings. Patch SIESTA
> only if you also want the density-grid restart routes.

### Step 1: Patch SIESTA (optional; only for the density-grid routes)

cube4siesta includes a patch file for SIESTA 4.1.5. Apply it and
rebuild:

```bash
cd siesta-4.1.5
patch -p1 < /path/to/cube4siesta/siesta-4.1.5-rho-restart.patch
cd Obj
sh ../Src/obj_setup.sh
cp intel.make arch.make       # pick the template for your compiler
make -j 8 siesta
```

If you're not sure how to build SIESTA, see the
[YHKLab Wiki build guide](https://yhklab.github.io/YHKimLabWiki/site/installation/siesta-cluster/).

### Step 2: Install cube4siesta

```bash
pip install -e /path/to/cube4siesta
```

This gives you the `cube4siesta` command. You'll also need Python 3.9+,
NumPy, and SciPy (installed automatically). If you're working with VASP
files, you'll also need `pymatgen` (`pip install pymatgen`).

---

## Main route: `dm-seed` (filtered DM restart)

This is the recommended workflow. It needs **no SIESTA patch** and **no
density matrix from the source code** — just a difference density.

### 1. Get a difference density from your other code

`dm-seed` consumes `Δρ = ρ_scf − ρ_atoms`:

- **OpenMX** writes it directly as `*.dden.cube`.
- **Quantum ESPRESSO / VASP**: subtract the atomic-superposition density
  from the SCF density (e.g. `cube4siesta-vasp-diff` for VASP, see below),
  or pass the total density and let the basis filter do the work.

### 2. Run a cheap SIESTA baseline once

`dm-seed` reuses a SIESTA run's geometry, basis, `.RHO` grid, and `.DM`
sparsity pattern. Any quick SIESTA run on the *same structure* works
(even unconverged — only the layout is borrowed). It produces
`<stem>.DM`, `<stem>.ORB_INDX`, `<stem>.STRUCT_OUT`, `<stem>.RHO`, and
`<species>.ion` in one directory.

### 3. Build the seed `.DM`

```bash
cube4siesta dm-seed \
    --baseline path/to/siesta_baseline \
    --stem vsse \
    --dden VSSe.dden.cube \
    --output vsse_seed.DM
```

This fits `Δρ` onto the SIESTA basis and writes
`D_seed = D_atomic + fit(Δρ)`. A small, fully runnable example ships in
the repo under `examples/vsse_openmx/dm_seed_demo/` (the OpenMX `Δρ` is
provided as a compact `.RHO` so no 13 MB cube is needed):

```bash
cube4siesta dm-seed \
    --baseline examples/vsse_openmx/dm_seed_demo \
    --stem vsse \
    --dden examples/vsse_openmx/dm_seed_demo/vsse_dden.RHO \
    --output vsse_seed.DM
```

(VSSe is the hard case — OpenMX treats V 3s/3p as valence, the SIESTA
`V.psf` does not. The basis filter removes that semicore so it never
reaches the valence block.)

### 4. Refine the seed with one SIESTA diagonalization

Rename the seed to `<SystemLabel>.DM` in your SIESTA run directory and
add the settings `dm-seed` prints:

```fdf
DM.UseSaveDM         true
MaxSCFIterations     1
DM.MixingWeight      1.0
SCFMustConverge      false
```

SIESTA rebuilds `ρ[D_seed]`, diagonalizes once, and the resulting
`<SystemLabel>.DM`, bands, and density are your transferred result.

---

## Quickstart: your first cross-code restart

Here's the typical workflow, step by step.

### 1. Get a cube file from your other code

**If you're coming from OpenMX:**
OpenMX already writes `*.tden.cube` (total density) by default.
You're done — skip to step 2.

**If you're coming from VASP:**
VASP stores the density in a special format called CHGCAR. We need
to convert it to a cube file first:

```bash
python -c "
from cube4siesta.vasp_io import chgcar_to_cube
chgcar_to_cube('CHGCAR', 'my_density.cube')
"
```

**If you're coming from Quantum ESPRESSO:**
Use QE's own `pp.x` post-processing tool with `plot_num=0` to write
a cube file of the charge density.

### 2. Run SIESTA once (to determine the grid)

SIESTA needs to know what grid to use. The easiest way is to run a
quick SIESTA calculation on the same structure:

```bash
mpirun -np 4 siesta < my_system.fdf > baseline.out
```

This creates `my_system.RHO`, which tells cube4siesta what grid size
and cell shape SIESTA is using.

### 3. Convert the cube file

```bash
cube4siesta convert \
    --cube my_density.cube \
    --output my_system.RHOIN \
    --from-siesta-rho my_system.RHO \
    --verify
```

The `--verify` flag reads back the written file to make sure
everything looks correct. The `--from-siesta-rho` flag ensures
the grid matches SIESTA's expectations — if the grids are different,
cube4siesta will automatically resample.

### 4. Run SIESTA with the imported density

Add these two lines to your `.fdf` input file:

```
Rho.Restart       true
Rho.RestartFile   my_system.RHOIN
```

Then run SIESTA as usual:

```bash
mpirun -np 4 siesta < my_system.fdf > restart.out
```

SIESTA will read the density, build the Hamiltonian from it, and
produce a `.DM` file. You can now continue with a normal SIESTA SCF
by removing the `Rho.Restart` lines and adding `DM.UseSaveDM true`.

---

## Pseudopotentials: an important note

For the restart to be physically meaningful, **SIESTA's pseudopotential
must describe the same electrons as the source code's pseudopotential.**

For example, if OpenMX treats vanadium as having 13 valence electrons
(including 3s and 3p) but your SIESTA `.psf` file only has 5 (3d and
4s), then the two codes are literally computing different physical
quantities — and importing the density won't give sensible results.

### How to get matching pseudopotentials

cube4siesta can help! If you have the source code's pseudopotential
file, it can automatically generate the input needed to create a
matching SIESTA pseudopotential using the ATOM code:

**From an OpenMX calculation:**

```bash
cube4siesta gen-atom-input \
    --from-vps /path/to/V_PBE19.vps \
    --output V.pg.inp
```

**From a Quantum ESPRESSO calculation:**

```bash
cube4siesta gen-atom-input \
    --from-upf /path/to/Si.pbe-n-van.UPF \
    --output Si.pg.inp
```

Both commands will:
- Read the source pseudopotential file
- Extract the element, exchange-correlation functional, valence
  configuration, and cutoff radii
- Write an ATOM input file (`.inp`) ready to use

Then generate the `.psf` with the ATOM code:

```bash
sh /path/to/atom/bin/pg.sh V.pg.inp
# This creates V.pg.psf — rename and use in your SIESTA calculation
```

For a tutorial on the ATOM code and pseudopotential generation, see the
[YHKLab Wiki](https://yhklab.github.io/YHKimLabWiki/site/atom/pseudopotential/)
or the [SIESTA ATOM documentation](https://docs.siesta-project.org/projects/atom/en/latest/tutorial/ps-generation/).

**For VASP users:** VASP uses PAW pseudopotentials, which are
fundamentally different from SIESTA's norm-conserving format. Automatic
conversion is not possible (see [issue #3](https://github.com/simplecaryu/cube4siesta/issues/3)).
In practice, for common elements like Si, C, O, the standard Cornell
NNIN `.psf` files work well with VASP densities. You can download them
from the [Cornell NNIN Virtual Vault](https://nninc.cnf.cornell.edu/).

### Alternative: transfer the difference density instead

If regenerating the pseudopotential is not an option, there is a second
route that keeps your standard SIESTA `.psf`: transfer only the
**difference density** Δρ = ρ_scf − ρ_atomic and let SIESTA add its own
atomic-superposition density back (`Rho.Restart.Diff`). The
pseudo-dependent part of ρ lives almost entirely in ρ_atomic, so Δρ
transfers across different core/valence partitions where the total
density cannot.

```bash
# OpenMX writes the difference density natively (*.dden.cube):
cube4siesta convert --cube SYSTEM.dden.cube --output system.RHOIN.diff \
    --from-siesta-rho system.RHO --diff --verify

# VASP: build it from two CHGCARs (SCF run, and ICHARG=1 NELM=0 run):
cube4siesta-vasp-diff --scf CHGCAR --atomic CHGCAR.atomic --out dden.cube

# then add to your fdf:  Rho.Restart.Diff true
```

Measured on VSSe (OpenMX V has 13 valence electrons, SIESTA's V.psf
has 5): total-ρ restart misses the baseline by 156 eV, while the
diff-density restart lands within 3.5 eV, reproduces the occupied bands
to 0.21 eV, and the resulting density matrix is within 10 % of the
converged one. **Only use this when the partitions actually differ** —
with matching pseudopotentials the plain total-ρ restart is more
accurate. Note this is a norm-conserving-source technique: from PAW
sources (VASP) the diff route gives no improvement, so for a PAW
source with mismatched valence regenerate the pseudopotential instead
(see `docs/issues/001-cross-pseudo-diff-density.md` for the measured
numbers).

### Pre-packaged pseudopotentials

We include tested GGA-PBE `.psf` files for C, Si, V, S, Se, Mn, and Te
in `testdata/pseudos/`, downloaded from the Cornell NNIN Virtual Vault.

---

## The no-patch route: direct density-matrix projection (OpenMX)

If your source code is OpenMX, you can skip the density transfer (and
the SIESTA patch) entirely: OpenMX already stores its density matrix in
the `.scfout` file, and `cube4siesta project-dm` projects it directly
onto SIESTA's basis, writing a SIESTA `.DM`. This also handles valence
mismatches automatically — any density SIESTA's basis cannot represent
(e.g. semicore electrons) simply drops out of the projection.

### Convenience shortcut: `openmx2siesta-dm`

For the OpenMX projection route, this command automates the SIESTA
template run and file wrangling. It is a separate entry point from the
`cube4siesta` subcommands:

```bash
openmx2siesta-dm \
    --scfout omx/MoS2.scfout \
    --siesta-pseudo-dir siesta_pseudos/ \
    --siesta-cmd siesta \
    --output mos2_omx_proj_pure.DM
```

What it automates:

1. Reads OpenMX `DATA.PATH`, species PAO tags, geometry, and `scf.Kgrid` from
   the embedded input text in `*.scfout`.
2. Finds the corresponding OpenMX `*.pao` files under `DATA.PATH/PAO`.
3. Builds a minimal SIESTA template input using the OpenMX `scf.Kgrid`, runs
   SIESTA once, and obtains `<SystemLabel>.DM`, `<SystemLabel>.ORB_INDX`, and
   `<Element>.ion`.
4. Reads the SIESTA target electron count from the template output, purifies the
   projected matrix to that count, and writes the output `.DM`.

SIESTA 4.1.5 looks for pseudopotentials as `Label.vps` or `Label.psf` in the
run directory; it does not have an OpenMX-like `DATA.PATH`. Therefore the
shortcut takes `--siesta-pseudo-dir` and symlinks/copies `Mo.psf`, `S.psf`, etc.
into the generated template directory.

Use the explicit `cube4siesta project-dm` command below when you want full
control over every input file.

### 1. Collect three things from your OpenMX calculation

- **`SYSTEM.scfout`** — add `HS.fileout  on` to your OpenMX input and
  run (or re-run with `scf.restart on`; it finishes in one step).
- **The `.pao` files** for each element — these are in your OpenMX
  installation under `DFT_DATA19/PAO/` (e.g. `Si7.0.pao`). Use the same
  ones your OpenMX input named.
- **The basis spec** for each element — read it off the
  `Definition.of.Atomic.Species` block of your OpenMX input: for
  `Si7.0-s2p2d1` the spec is `s2p2d1`.

### 2. Run SIESTA once, cheaply, to get the basis files

The projection needs to know SIESTA's orbitals and sparsity pattern.
Any SIESTA run of the same structure and basis produces them — **it
does not need to be converged**, so you can use
`MaxSCFIterations 1` + `SCFMustConverge false`:

```bash
mpirun -np 4 siesta < my_system.fdf > template.out
# this leaves behind: my_system.DM, my_system.ORB_INDX, <Element>.ion
```

### 3. Project

```bash
cube4siesta project-dm \
    --scfout   SYSTEM.scfout \
    --orb-indx my_system.ORB_INDX \
    --template my_system.DM \
    --ion  Si:Si.ion \
    --pao  Si:/path/to/openmx/DFT_DATA19/PAO/Si7.0.pao \
    --spec Si:s2p2d1 \
    --purify \
    --output my_system_projected.DM
```

(Repeat `--ion/--pao/--spec` once per element.) The tool prints the
electron count `Tr(P·S)` — check it is close to what SIESTA expects.
`--purify` cleans the projected matrix up to an idempotent density
matrix with exactly the right electron count (add `--qtot N` to set N
explicitly); use it when you plan to run SIESTA from the result. Skip
it if you want the mathematically closest representation of the OpenMX
density matrix itself.

### 4. Use it

```bash
cp my_system_projected.DM my_system.DM
# in your fdf:  DM.UseSaveDM  true
mpirun -np 4 siesta < my_system.fdf > restart.out
```

This works with a completely standard SIESTA — `DM.UseSaveDM` is a
stock feature. To check the quality of any `.DM` against a reference:

```bash
cube4siesta compare-dm --orb-indx my_system.ORB_INDX \
    --ref converged.DM --cand projected:my_system_projected.DM
```

---

## SIESTA fdf flags

| Flag | Default | What it does |
|------|---------|--------------|
| `Rho.Restart` | `false` | Tells SIESTA to read an external density and do one diagonalization step to produce a DM. |
| `Rho.RestartFile` | `<SystemLabel>.RHOIN` | Path to the density file that cube4siesta wrote. |
| `Rho.Restart.Diff` | `false` | Treat the file as a difference density (ρ − ρ_atomic) and add SIESTA's own atomic density on top. Requires `Rho.Restart`. |

When `Rho.Restart` is on, SIESTA automatically sets
`MaxSCFIterations=1` and turns off convergence checks — this is
expected behavior, not an error.

---

## Full CLI reference

### `cube4siesta convert`

Converts a cube file into SIESTA's `.RHO` format.

```
cube4siesta convert --cube INPUT.cube --output OUTPUT.RHOIN [options]

Options:
  --from-siesta-rho FILE   Match the grid to an existing SIESTA .RHO file
                           (recommended — handles grid differences automatically)
  --target-mesh Nx,Ny,Nz   Specify the target grid manually
  --order 1|3              Interpolation quality: 1 = fast (default), 3 = smooth
  --rescale-to N           Adjust the total electron count to N
  --diff                   Input is a difference density (e.g. OpenMX
                           *.dden.cube); keep it charge-neutral, for
                           use with Rho.Restart.Diff
  --subtract REF           Build a difference density: subtract a
                           reference cube/.RHO (the source code's
                           atomic-superposition density) before writing
  --verify                 Double-check the output by reading it back
```

### `cube4siesta-vasp-diff`

Builds a difference-density cube from two VASP CHGCARs, for
`Rho.Restart.Diff` mode (requires `pymatgen`).

```
cube4siesta-vasp-diff --scf CHGCAR --atomic CHGCAR.atomic --out dden.cube

  --scf      CHGCAR from the converged SCF run
  --atomic   CHGCAR from an ICHARG=1 NELM=0 run (atomic superposition,
             same structure and POTCARs)
  --out      output cube (integral ≈ 0)
```

### `cube4siesta gen-atom-input`

Reads a source pseudopotential and writes an ATOM input file for
generating a matching SIESTA `.psf`.

```
cube4siesta gen-atom-input --from-vps FILE.vps --output NAME.pg.inp
cube4siesta gen-atom-input --from-upf FILE.UPF --output NAME.pg.inp
```

Supported formats:
- OpenMX `.vps` files (fully automatic)
- Quantum ESPRESSO `.UPF` files (v1 and v2, fully automatic)
- VASP `POTCAR` — not supported; see [issue #3](https://github.com/simplecaryu/cube4siesta/issues/3)

### `cube4siesta project-dm`

Projects an OpenMX density matrix (`.scfout`) onto SIESTA's basis and
writes a SIESTA `.DM`. No SIESTA patch required. See the walkthrough
section above for where each input comes from.

```
cube4siesta project-dm --scfout FILE.scfout --orb-indx FILE.ORB_INDX \
    --template FILE.DM --ion EL:FILE.ion --pao EL:FILE.pao \
    --spec EL:s2p2d1 [--purify [--qtot N]] --output OUT.DM

  --ion/--pao/--spec   repeat once per element (EL:value)
  --purify             make the result idempotent with exactly --qtot
                       electrons (default: nearest integer)
  --spacing            quadrature grid spacing in Bohr (default 0.15)
```

### `cube4siesta compare-dm`

Element-wise comparison of SIESTA `.DM` files (relative Frobenius
error, correlations, largest deviation) — useful to judge any
transferred DM against a converged reference.

```
cube4siesta compare-dm --orb-indx FILE.ORB_INDX --ref REF.DM \
    --cand NAME:CAND.DM [--cand NAME2:CAND2.DM ...]
```

### Python API

If you prefer to script things in Python:

```python
from cube4siesta import read_rho, write_rho        # SIESTA .RHO files
from cube4siesta.cube_io import read_cube, write_cube  # Gaussian cubes
from cube4siesta.resample import resample_to_mesh   # 3D interpolation
from cube4siesta.vasp_io import chgcar_to_cube      # VASP CHGCAR → cube
from cube4siesta.gen_psf import parse_openmx_vps, parse_qe_upf, write_atom_input
from cube4siesta.dm_io import read_dm, write_dm     # SIESTA .DM files
from cube4siesta.scfout_io import read_scfout       # OpenMX .scfout
from cube4siesta.project import project_openmx_to_siesta, purify_canonical
```

---

## Examples

The `examples/` directory contains fully worked demonstrations:

- **`examples/sic_vasp/`** — Importing a VASP charge density into
  SIESTA for zincblende SiC. Includes pseudopotentials and ready-to-run
  `.fdf` files.

- **`examples/vsse_openmx/`** — Importing an OpenMX charge density for
  a Janus VSSe monolayer. Also demonstrates what goes wrong when the
  pseudopotentials don't match (and how to fix it).

- **`examples/DM_PROJECTION_RESULTS.md`** — the full measured
  comparison of all transfer routes (density matrices and band
  structures, for matched and mismatched pseudopotentials), with
  `examples/dm_projection_compare.sh` to reproduce it (set the
  `OPENMX_PAO` environment variable to your OpenMX `DFT_DATA19/PAO`).

---

## Project status

Done:
- [x] SIESTA `.RHO` binary reader/writer (with roundtrip tests)
- [x] Gaussian cube reader/writer
- [x] Automatic grid resampling (trilinear interpolation)
- [x] VASP CHGCAR → cube converter
- [x] SIESTA `Rho.Restart` patch
- [x] Pseudopotential matching tool (`gen-atom-input`)
- [x] Cross-pseudo difference-density mode (`--diff`, `Rho.Restart.Diff`;
      [issue #1](https://github.com/simplecaryu/cube4siesta/issues/1))
- [x] Direct OpenMX→SIESTA density-matrix projection (`project-dm`)
      with canonical purification, and `compare-dm`
- [x] End-to-end tested on H₂O, Si, SiC (VASP), VSSe (OpenMX, QE)

Future:
- [ ] VASP POTCAR auto-conversion ([issue #3](https://github.com/simplecaryu/cube4siesta/issues/3))
- [ ] Spin-polarized (nspin=2) support
- [ ] Non-collinear spin (nspin=4)
- [ ] `project-dm` for sources other than OpenMX (needs the source DM
      + basis tables; QE/VASP do not expose these as directly)

---

## How to get help

If something doesn't work or you're unsure how to set up your
calculation, please [open an issue](https://github.com/simplecaryu/cube4siesta/issues).

---

## License

The Python code is released under the MIT License. Pseudopotentials in
`testdata/pseudos/` are from the Cornell NNIN Virtual Vault.
The SIESTA patch is derivative of SIESTA 4.1.5 (GPL).
