"""The LTR-retrotransposon line: the satellite engine's clade test for every clade (not only the marker), so the loci table can
place each concentrated clade's core, and the marker clade's density for the overview figure."""
from __future__ import annotations

from pathlib import Path

import numpy as np

from .sat import te
from .sat.params import TE_BINS, TE_MIN_BIN
from .util import log


def te_line(gff, chroms: dict) -> dict:
    """{'clades': [(clade, {chrom: core})] most concentrated first, 'marker': (clade, cores, decisive) | None,
    'density': {chrom: (bin_bp, counts)} of the marker clade, 'mids': {clade: {chrom: sorted element midpoints}}}.
    Cores are the satellite engine's (0-based half-open). A file that gives no classified element on the chromosomes
    is warned about with the reason (te.te_census): the line is then off."""
    if gff is None:
        return dict(clades=[], marker=None, density={}, mids={})
    tes = te.load_te_gff(gff, chroms)
    if not tes:
        n = te.te_census(gff, chroms)
        why = (f"{n['ltr']:,} LTR retrotransposons on the chromosomes, none classified to a lineage (clade=, or a "
               f"DANTE_LTR Final_Classification); the LTR line needs lineages such as CRM or Athila (e.g. TEsorter "
               f"with REXdb, or LTRquest)" if n["ltr"] else
               f"no LTR retrotransposon on the chromosomes ({n['foreign']:,} on other sequence names: another "
               f"assembly or naming?)" if n["foreign"] else
               "no LTR_retrotransposon feature (or DANTE_LTR transposable_element)")
        log.warning(f"tes/{Path(gff).name}: {why}; LTR line off")
    clades = te.candidate_clades(tes, chroms)
    marker = te.choose_marker(tes, chroms, clades)
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
