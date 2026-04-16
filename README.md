# cube4siesta

**Restart a SIESTA calculation from a charge density produced by
another DFT code** — VASP, OpenMX, Quantum ESPRESSO, or anything that
can write a Gaussian cube file.

If you've ever run a DFT calculation in one code and wished you could
continue or analyze that result in SIESTA without starting from
scratch, this tool is for you.

---

## What does it do?

1. **Reads** the charge density from your other code (Gaussian cube
   format, or VASP CHGCAR).
2. **Converts** it onto SIESTA's internal grid format (`.RHO`).
3. A small **SIESTA patch** (`Rho.Restart`) tells SIESTA to use that
   density as a starting point instead of guessing from atomic
   orbitals.

The result is a standard SIESTA `.DM` file that you can use for any
further SIESTA calculation — band structures, geometry optimization,
transport, etc.

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

### Pre-packaged pseudopotentials

We include tested GGA-PBE `.psf` files for C, Si, V, S, Se, Mn, and Te
in `testdata/pseudos/`, downloaded from the Cornell NNIN Virtual Vault.

---

## SIESTA fdf flags

| Flag | Default | What it does |
|------|---------|--------------|
| `Rho.Restart` | `false` | Tells SIESTA to read an external density and do one diagonalization step to produce a DM. |
| `Rho.RestartFile` | `<SystemLabel>.RHOIN` | Path to the density file that cube4siesta wrote. |

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
  --verify                 Double-check the output by reading it back
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

### Python API

If you prefer to script things in Python:

```python
from cube4siesta import read_rho, write_rho        # SIESTA .RHO files
from cube4siesta.cube_io import read_cube, write_cube  # Gaussian cubes
from cube4siesta.resample import resample_to_mesh   # 3D interpolation
from cube4siesta.vasp_io import chgcar_to_cube      # VASP CHGCAR → cube
from cube4siesta.gen_psf import parse_openmx_vps, parse_qe_upf, write_atom_input
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

---

## Project status

Done:
- [x] SIESTA `.RHO` binary reader/writer (with roundtrip tests)
- [x] Gaussian cube reader/writer
- [x] Automatic grid resampling (trilinear interpolation)
- [x] VASP CHGCAR → cube converter
- [x] SIESTA `Rho.Restart` patch
- [x] Pseudopotential matching tool (`gen-atom-input`)
- [x] End-to-end tested on H₂O, SiC (VASP), VSSe (OpenMX)

In progress:
- [ ] Cross-pseudo difference-density mode ([PR #2](https://github.com/simplecaryu/cube4siesta/pull/2), [issue #1](https://github.com/simplecaryu/cube4siesta/issues/1))
- [ ] VASP POTCAR auto-conversion ([issue #3](https://github.com/simplecaryu/cube4siesta/issues/3))

Future:
- [ ] Spin-polarized (nspin=2) support
- [ ] Non-collinear spin (nspin=4)

---

## How to get help

If something doesn't work or you're unsure how to set up your
calculation, please [open an issue](https://github.com/simplecaryu/cube4siesta/issues).

---

## License

The Python code is released under the MIT License. Pseudopotentials in
`testdata/pseudos/` are from the Cornell NNIN Virtual Vault.
The SIESTA patch is derivative of SIESTA 4.1.5 (GPL).
