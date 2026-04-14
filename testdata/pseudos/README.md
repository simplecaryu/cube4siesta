# Pseudopotentials for cube4siesta integration tests

PSML 1.1 pseudopotentials from [PseudoDojo](http://www.pseudo-dojo.org/),
set **NC-SR-04 PBE standard**.

These are used when testing the full SIESTA `Rho.Restart` pipeline against
densities produced by other DFT codes (OpenMX VSSe, VASP SiC, etc.).

## Files

| Element | Valence | Notes |
|---------|---------|-------|
| C       | 4       | VASP SiC integration test partner |
| Si      | 4       | VASP SiC integration test partner |
| V       | 13      | OpenMX VSSe (transition metal) |
| S       | 6       | OpenMX VSSe |
| Se      | 6       | OpenMX VSSe |
| Mn      | 15      | (future VASP MnTe test) |
| Te      | 6       | (future VASP MnTe test) |

## Source

Downloaded from:
```
https://www.pseudo-dojo.org/pseudos/nc-sr-04_pbe_standard_psml.tgz
```

Redistribution terms: PseudoDojo pseudopotentials are freely available under
the PseudoDojo terms of use (see http://www.pseudo-dojo.org/about.html).

## SIESTA usage

```
%block PS_File_List
V  V.psml
S  S.psml
Se Se.psml
%endblock PS_File_List
```

Note: SIESTA 4.1.5 needs to be built with PSML support for these to work.
The tarball build we made for this project links libPSML automatically
when the `libxmlf90`/`libpsml` stubs are present — check
`libSiestaXC`/`libfdf` build output for PSML messages.
