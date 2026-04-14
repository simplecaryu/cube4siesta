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

## Which mode should I use?

Two modes are available, and they are **not** ordered by quality —
each fails in a different way.

| | Total-ρ restart | Diff-ρ restart |
|---|---|---|
| CLI | `cube4siesta convert ...` | `cube4siesta convert --diff ...` |
| fdf | `Rho.Restart true` | `Rho.Restart true`<br>`Rho.Restart.Diff true` |
| What's transferred | ρ(r) itself | ρ(r) − ρ_atomic(r) |
| Atomic part provided by | the source code (implicit in ρ) | SIESTA's own `rhoatm` |

Measured on our test cases:

|  | Baseline SCF | Total-ρ | Diff-ρ |
|---|---|---|---|
| **SiC** (VASP PAW → Cornell NC, same valence) | −263.98 | **−264.54** | −284.06 |
| **VSSe** (OpenMX 13-val V → Cornell 5-val V)  | −696.61 | −540.92 | **−693.13** |

### Why diff-ρ can make things worse

Diff-ρ performs

$$\rho_\text{inject} \;=\; \underbrace{(\rho_\text{scf}^\text{src} - \rho_\text{atom}^\text{src})}_{\text{from the cube}} \;+\; \underbrace{\rho_\text{atom}^\text{SIESTA}}_{\text{SIESTA's rhoatm}}
\;=\; \rho_\text{scf}^\text{src} \;+\; \bigl(\rho_\text{atom}^\text{SIESTA} - \rho_\text{atom}^\text{src}\bigr).$$

That last bracket is zero only if the source and target atomic
densities are identical. They are identical exactly when both codes use
the *same pseudopotential family* generated the same way; otherwise
diff-ρ injects a spurious atomic-shape residual.

Measured for SiC (both codes use 4+4 valence):

- SIESTA `rhoatm` peak (Troullier-Martins NC): **0.260 e/Bohr³**
- VASP `ρ_atomic` peak (PAW augmented): **0.588 e/Bohr³**
- Pointwise peak |Δ|: **0.57 e/Bohr³** ⇒ ~20 eV total-energy error.

So "same valence" is necessary but not sufficient to make diff mode a
no-op. The frameworks have to agree on the *shape* of ρ_atomic too,
and PAW ↔ NC do not.

### Rule of thumb

| Source code | Framework | Recommended |
|---|---|---|
| VASP | PAW | **total-ρ** (diff-ρ injects PAW↔NC shape error) |
| Quantum ESPRESSO with NC pseudo | NC | total-ρ if valence matches SIESTA's; else diff-ρ |
| OpenMX | NC | total-ρ if valence matches; **diff-ρ if source has semicore that SIESTA's pseudo freezes** |
| OpenMX with same pseudo choice | NC | total-ρ (diff-ρ reduces to the same thing) |

In short: diff-ρ only earns its keep when the valence partitions
disagree *and* both sides are in the same pseudo family. In every
other case prefer total-ρ.

## Quickstart for OpenMX users

OpenMX writes `*.tden.cube` (total density) and `*.dden.cube`
(ρ − ρ_atomic) by default.

**Recommended when you suspect valence mismatch — difference density**

```bash
# 1. Run SIESTA once to fix the target mesh and obtain baseline values
mpirun -np 4 siesta < your.fdf > base.out           # writes your.RHO

# 2. Convert OpenMX's difference-density cube onto SIESTA's mesh
cube4siesta convert \
    --cube path/to/System.dden.cube \
    --output your.RHOIN.diff \
    --from-siesta-rho your.RHO \
    --diff --verify

# 3. Run SIESTA with Rho.Restart + Rho.Restart.Diff
#
#    Add to your.fdf:
#      Rho.Restart       true
#      Rho.Restart.Diff  true
#      Rho.RestartFile   your.RHOIN.diff
#
mpirun -np 4 siesta < your.fdf > rr.out
```

See `examples/vsse_openmx/` for a worked example on a Janus VSSe
monolayer.

**Simpler path — total density (only OK when pseudo choices match)**

```bash
cube4siesta convert \
    --cube path/to/System.tden.cube \
    --output your.RHOIN \
    --from-siesta-rho your.RHO --verify

# fdf: Rho.Restart true, Rho.RestartFile your.RHOIN  (leave Rho.Restart.Diff off)
```

---

## Quickstart for VASP users

VASP doesn't output a difference density directly, but you can get one
cheaply with a second no-SCF run.

```bash
# 1. Run a normal VASP SCF and keep CHGCAR
cp INCAR.my_scf INCAR
vasp_std
cp CHGCAR CHGCAR.scf

# 2. Run VASP again with ICHARG=1 and NELM=0 — atomic superposition, no SCF
cat > INCAR <<'EOF'
ICHARG = 1
NELM   = 0
LCHARG = .TRUE.
EOF
vasp_std
cp CHGCAR CHGCAR.atomic

# 3. Build the Δρ cube
cube4siesta-vasp-diff \
    --scf    CHGCAR.scf \
    --atomic CHGCAR.atomic \
    --out    dden.cube

# 4. Run SIESTA once for the baseline mesh, then convert + restart
mpirun -np 4 siesta < your.fdf > base.out

cube4siesta convert --cube dden.cube --output your.RHOIN.diff \
    --from-siesta-rho your.RHO --diff --verify

# add:  Rho.Restart true  /  Rho.Restart.Diff true  /  Rho.RestartFile your.RHOIN.diff
mpirun -np 4 siesta < your.fdf > rr.out
```

If you *know* your VASP pseudopotential and your SIESTA `.psf` agree
on core/valence (e.g. a PBE PAW with no semicore vs a plain NC .psf
on the same element), you can skip the atomic-superposition run and
inject the raw CHGCAR as a total density:

```bash
python -c "from cube4siesta.vasp_io import chgcar_to_cube; \
           chgcar_to_cube('CHGCAR', 'out.cube')"
cube4siesta convert --cube out.cube --output your.RHOIN \
    --from-siesta-rho your.RHO --verify
```

See `examples/sic_vasp/` for a worked example on zincblende SiC
(from the VASP 6.5 testsuite).

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
