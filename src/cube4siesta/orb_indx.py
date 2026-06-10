"""
Parse SIESTA .ORB_INDX files.

The .ORB_INDX file enumerates every orbital of the *supercell* used to fold the
sparse density matrix. Its body has one line per supercell orbital:

    io  ia is  spec iao  n  l  m  z  p   sym    rc   isc1 isc2 isc3  iuo

  io    : supercell orbital index (1..no_s)
  ia    : atom index in the supercell
  is    : species index
  spec  : species label
  iao   : orbital index within the atom (1..norb_per_atom)
  n,l,m,z : principal / angular / magnetic / zeta quantum numbers
  p     : polarisation flag (T/F)
  sym   : human-readable symbol (e.g. "s", "py", "Pdxy")
  rc    : cutoff radius (Bohr)
  isc   : supercell cell offset (integer multiples of the lattice vectors)
  iuo   : index of the corresponding *unit-cell* orbital (1..no_u)

The first no_u entries (isc == 0,0,0) ARE the unit cell. This file lets us
decode the column indices stored in a .DM (which run over supercell orbitals)
into (unit-cell orbital, lattice cell offset), and gives the (l, m, z) labels we
need to attach radial functions to each orbital.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np


@dataclass
class Orbital:
    """A single supercell orbital entry from .ORB_INDX."""
    io: int          # supercell orbital index (1-based)
    ia: int          # supercell atom index (1-based)
    ispec: int       # species index (1-based)
    spec: str        # species label
    iao: int         # orbital index within atom (1-based)
    n: int
    l: int
    m: int
    z: int           # zeta
    pol: bool        # polarisation orbital?
    sym: str
    rc: float
    isc: tuple[int, int, int]   # cell offset in lattice-vector units
    iuo: int         # corresponding unit-cell orbital index (1-based)


@dataclass
class OrbIndx:
    no_u: int
    no_s: int
    orbitals: list[Orbital]   # length no_s, ordered by io

    @property
    def unit(self) -> list[Orbital]:
        """The no_u unit-cell orbitals (io = 1..no_u, isc = 0)."""
        return self.orbitals[: self.no_u]

    def atom_of_uorb(self) -> np.ndarray:
        """(no_u,) array: unit-cell atom index (1-based) for each unit-cell orbital."""
        return np.array([o.ia for o in self.unit], dtype=np.int64)

    def decode_column(self, jo: int) -> tuple[int, tuple[int, int, int]]:
        """
        Decode a 1-based supercell column index ``jo`` (as stored in .DM listd)
        into (unit-cell orbital index [1-based], cell offset isc).
        """
        o = self.orbitals[jo - 1]
        return o.iuo, o.isc


def read_orb_indx(path: str | Path) -> OrbIndx:
    """Parse a SIESTA .ORB_INDX file."""
    orbitals: list[Orbital] = []
    no_u = no_s = 0
    with open(path) as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            if "orbitals in unit cell and supercell" in line:
                parts = line.split()
                no_u, no_s = int(parts[0]), int(parts[1])
                continue
            if line.startswith("io") or line.startswith("#"):
                continue
            tok = line.split()
            # data rows start with an integer io; skip anything else (trailer text)
            if len(tok) < 16 or not tok[0].lstrip("-").isdigit():
                continue
            io = int(tok[0])
            if no_s and io > no_s:
                # past the orbital table into the trailing supercell-offset legend
                continue
            orbitals.append(
                Orbital(
                    io=io,
                    ia=int(tok[1]),
                    ispec=int(tok[2]),
                    spec=tok[3],
                    iao=int(tok[4]),
                    n=int(tok[5]),
                    l=int(tok[6]),
                    m=int(tok[7]),
                    z=int(tok[8]),
                    pol=(tok[9].upper() == "T"),
                    sym=tok[10],
                    rc=float(tok[11]),
                    isc=(int(tok[12]), int(tok[13]), int(tok[14])),
                    iuo=int(tok[15]),
                )
            )

    if no_s == 0:
        no_s = len(orbitals)
    if no_u == 0:
        no_u = sum(1 for o in orbitals if o.isc == (0, 0, 0))
    orbitals.sort(key=lambda o: o.io)
    return OrbIndx(no_u=no_u, no_s=no_s, orbitals=orbitals)
