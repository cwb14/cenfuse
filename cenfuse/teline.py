"""The LTR-retrotransposon line: cenrep's clade test for every clade (not only the marker), so the loci table can
place each concentrated clade's core, and the marker clade's density for the overview figure."""
from __future__ import annotations

import numpy as np

from .sat import te
from .sat.params import TE_BINS, TE_MIN_BIN


def te_line(gff, chroms: dict) -> dict:
    """{'clades': [(clade, {chrom: core})] most concentrated first, 'marker': (clade, cores, decisive) | None,
    'density': {chrom: (bin_bp, counts)} of the marker clade, 'mids': {clade: {chrom: sorted element midpoints}}}.
    Cores are cenrep's (0-based half-open)."""
    if gff is None:
        return dict(clades=[], marker=None, density={}, mids={})
    tes = te.load_te_gff(gff, chroms)
    clades = te.candidate_clades(tes, chroms)
    marker = te.marker(clades)
    dens = {}
    if marker:
        mids = tes.get(marker[0], {})
        for c, L in chroms.items():
            b = max(TE_MIN_BIN, L // TE_BINS)
            m = mids.get(c, np.zeros(0, np.int64))
            dens[c] = (b, np.bincount(np.minimum(m // b, L // b), minlength=L // b + 1))
    return dict(clades=clades, marker=marker, density=dens, mids=tes)


def composition(mids: dict, chrom: str, s: int, e: int, top: int = 3) -> tuple:
    """(element count, 'clade:n,...' of the most common clades) of LTR retrotransposons with midpoints in [s, e)."""
    n = {c: int(np.searchsorted(d[chrom], e) - np.searchsorted(d[chrom], s)) for c, d in mids.items() if chrom in d}
    n = {c: k for c, k in n.items() if k}
    best = sorted(n.items(), key=lambda x: (-x[1], x[0]))[:top]
    return sum(n.values()), ",".join(f"{c}:{k}" for c, k in best) or "."
