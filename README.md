# cube4siesta

Restart a SIESTA calculation from a real-space charge density produced
by another DFT code (VASP, OpenMX, Quantum ESPRESSO, …).

Takes a Gaussian cube file of ρ(r) in, emits a SIESTA `.RHO`-format
file out. Pairs with a small SIESTA 4.1.5 patch that adds the
`Rho.Restart` fdf flag, so SIESTA runs one `ρ → V_KS → H → DM` step
from your external density instead of the usual atomic-superposition
initial guess.

## Installation

### SIESTA with the Rho.Restart patch

Apply the patch from the sibling `siesta-4.1.5/` repo and rebuild:

```bash
cd siesta-4.1.5/Obj
sh ../Src/obj_setup.sh
cp intel.make arch.make       # edit for your system
make -j 8 siesta
```

### cube4siesta Python package

```bash
pip install -e .
```

Dependencies: `numpy`, `scipy`, and optionally `pymatgen` (for VASP
CHGCAR conversion).

## Quickstart

```bash
# 1. Get a cube file from the source code
#   OpenMX:  *.tden.cube is written by default
#   VASP:    python -c "from cube4siesta.vasp_io import chgcar_to_cube; \
#                        chgcar_to_cube('CHGCAR','scf.cube')"

# 2. Run SIESTA once to fix the target mesh
mpirun -np 4 siesta < your.fdf > base.out       # writes your.RHO

# 3. Convert the cube onto SIESTA's mesh
cube4siesta convert --cube scf.cube \
    --output your.RHOIN \
    --from-siesta-rho your.RHO --verify

# 4. Run SIESTA with Rho.Restart
#   Add to your.fdf:
#     Rho.Restart       true
#     Rho.RestartFile   your.RHOIN
mpirun -np 4 siesta < your.fdf > rr.out
```

The output is a standard `.DM` / `.EIG` / `.FA` set — continue with
a normal SCF by removing the `Rho.Restart` lines and adding
`DM.UseSaveDM true`.

## SIESTA fdf flags

| Flag | Default | Meaning |
|------|---------|---------|
| `Rho.Restart` | `false` | Enable the ρ → one-step DM path. Forces `MaxSCFIterations=1`. |
| `Rho.RestartFile` | `<SystemLabel>.RHOIN` | Path to the input density file. |

## CLI reference

```
cube4siesta convert --cube IN.cube --output OUT.RHOIN [options]

  --from-siesta-rho REF.RHO   borrow mesh + cell from an existing .RHO
  --target-mesh Nx,Ny,Nz      target grid (alternative to --from-siesta-rho)
  --order 1|3                  interpolation order (default: 1 = trilinear)
  --rescale-to N               rescale so ∫ρ dV = N electrons
  --verify                     read back and confirm roundtrip
```

## Examples

- `examples/sic_vasp/` — VASP CHGCAR → SIESTA (zincblende SiC)
- `examples/vsse_openmx/` — OpenMX cube → SIESTA (Janus VSSe monolayer)

## Pseudopotentials

SIESTA 4.1.5 reads `.psf` format. Pre-tested GGA-PBE pseudos for
C, Si, V, S, Se, Mn, Te are in `testdata/pseudos/` (from the
[Cornell NNIN Virtual Vault](https://nninc.cnf.cornell.edu/)).

## Status

- [x] SIESTA `.RHO` I/O with roundtrip tests
- [x] Gaussian cube I/O
- [x] Trilinear resampling
- [x] VASP CHGCAR → cube converter
- [x] `Rho.Restart` SIESTA patch
- [x] End-to-end validation (H₂O, SiC, VSSe)
- [ ] Spin-polarized (nspin=2) support
- [ ] Cross-pseudo difference-density mode (see [#1](https://github.com/simplecaryu/cube4siesta/issues/1))
