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
transport, etc. There are two ways to get there:

- **Density route** (`convert` + a small SIESTA patch): convert the
  charge density to SIESTA's grid format (`.RHO`), then a patched
  SIESTA reads it and produces the `.DM` in a single step.
- **Density-matrix route** (`project-dm`, OpenMX only): read the
  density matrix OpenMX already stores and project it directly onto
  SIESTA's basis. Pure post-processing — **no SIESTA patch needed**.

---

## Which route should I use?

The key question is whether your source code's pseudopotential counts
the **same valence electrons** as your SIESTA `.psf` (e.g. OpenMX's
default vanadium has 13 valence electrons, a typical SIESTA `V.psf`
has 5 — that is a *mismatch*). If you're not sure, run
`cube4siesta convert` on the source cube and check whether the printed
electron count matches what SIESTA expects for your structure.

| Your situation | Use this | Patch needed? | Accuracy (measured) |
|---|---|---|---|
| Valence counts **match** | total-density restart (Quickstart below) | yes | best (DM within ~2 %) |
| Valence counts **differ**, source is OpenMX / QE-with-NC-pseudos | difference-density restart (`--diff`) | yes | best under mismatch (DM ~10 %) |
| Valence counts **differ**, source is OpenMX, you can't patch SIESTA | DM projection (`project-dm --purify`) | **no** | close second (DM ~16 %, bands match the Δρ route) |
| Valence counts **differ**, source is VASP (PAW) | regenerate the SIESTA pseudo (`gen-atom-input`) | yes | no density route works from PAW under mismatch |

(The accuracy numbers are relative Frobenius errors of the resulting
density matrix vs a fully converged SIESTA reference; details and band
structures in `examples/DM_PROJECTION_RESULTS.md`.)

---

## Installation

### Step 1: Patch SIESTA

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
