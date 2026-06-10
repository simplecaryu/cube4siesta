"""
Compare SIESTA density matrices element-by-element.

The deliverable metric of the cross-code study is *DM-value similarity*: how close
is a converted / restarted density matrix to SIESTA's self-consistent one, in the
same NAO basis. We compare on the shared set of (orbital_a, orbital_b, cell)
entries, decoding supercell column indices via .ORB_INDX so DMs with slightly
different sparsity patterns are still comparable.

Raw correlation over all entries is dominated by the tens of thousands of
near-zero long-range tail elements, so we also report magnitude-thresholded
correlation, the relative Frobenius error, and an on-site (same-atom) breakdown.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .dm_io import DMFile
from .orb_indx import OrbIndx


def _entries(dm: DMFile, oi: OrbIndx) -> dict[tuple[int, int, tuple[int, int, int]], float]:
    """Map {(row_uorb, col_uorb, cell): value} for spin 0."""
    out: dict[tuple[int, int, tuple[int, int, int]], float] = {}
    for io in range(dm.no_u):
        listd, vals = dm.row(io)
        for jo, v in zip(listd, vals[:, 0]):
            juo, isc = oi.decode_column(int(jo))
            out[(io, juo - 1, isc)] = float(v)
    return out


@dataclass
class Comparison:
    n_shared: int
    n_ref: int
    n_cand: int
    rel_frobenius: float
    max_abs_diff: float
    corr_all: float
    corr_sig: float            # correlation on |ref| > sig_thresh
    sig_thresh: float
    n_sig: int
    rel_frob_onsite: float     # same-atom block only
    corr_onsite: float
    n_onsite: int


def compare(
    ref: DMFile,
    cand: DMFile,
    oi: OrbIndx,
    atom_of_uorb: np.ndarray | None = None,
    sig_thresh: float = 0.02,
) -> Comparison:
    """Compare candidate DM ``cand`` against reference ``ref`` (same basis)."""
    R = _entries(ref, oi)
    C = _entries(cand, oi)
    keys = [k for k in R if k in C]
    rv = np.array([R[k] for k in keys])
    cv = np.array([C[k] for k in keys])
    diff = cv - rv

    rel_frob = float(np.linalg.norm(diff) / np.linalg.norm(rv)) if rv.size else float("nan")
    corr_all = float(np.corrcoef(rv, cv)[0, 1]) if rv.size > 1 else float("nan")

    sig = np.abs(rv) > sig_thresh
    corr_sig = float(np.corrcoef(rv[sig], cv[sig])[0, 1]) if sig.sum() > 1 else float("nan")

    if atom_of_uorb is None:
        atom_of_uorb = oi.atom_of_uorb()
    onsite = np.array(
        [k[2] == (0, 0, 0) and atom_of_uorb[k[0]] == atom_of_uorb[k[1]] for k in keys]
    )
    if onsite.sum() > 1:
        rel_frob_on = float(np.linalg.norm(diff[onsite]) / np.linalg.norm(rv[onsite]))
        corr_on = float(np.corrcoef(rv[onsite], cv[onsite])[0, 1])
    else:
        rel_frob_on = corr_on = float("nan")

    return Comparison(
        n_shared=len(keys),
        n_ref=len(R),
        n_cand=len(C),
        rel_frobenius=rel_frob,
        max_abs_diff=float(np.max(np.abs(diff))) if diff.size else float("nan"),
        corr_all=corr_all,
        corr_sig=corr_sig,
        sig_thresh=sig_thresh,
        n_sig=int(sig.sum()),
        rel_frob_onsite=rel_frob_on,
        corr_onsite=corr_on,
        n_onsite=int(onsite.sum()),
    )


def format_table(name_to_cmp: dict[str, Comparison]) -> str:
    """Pretty-print a ranking table of several candidates vs one reference."""
    hdr = (
        f"{'candidate':<22}{'relFrob':>9}{'corr_sig':>10}{'corr_all':>10}"
        f"{'maxdiff':>9}{'relFrob_on':>12}{'shared':>8}"
    )
    lines = [hdr, "-" * len(hdr)]
    # rank by significant-element relative Frobenius (lower = closer to SIESTA)
    for name, c in sorted(name_to_cmp.items(), key=lambda kv: kv[1].rel_frobenius):
        lines.append(
            f"{name:<22}{c.rel_frobenius:>9.4f}{c.corr_sig:>10.4f}{c.corr_all:>10.4f}"
            f"{c.max_abs_diff:>9.4f}{c.rel_frob_onsite:>12.4f}{c.n_shared:>8d}"
        )
    return "\n".join(lines)
