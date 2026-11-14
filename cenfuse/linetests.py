"""Do the lines of evidence agree with each other? One row per pair of lines, in the genome at hand.

Each line calls on its own, so two lines that place the centromere at one locus on most chromosomes, beyond chance, validate
each other without any truth; two that do not mean at least one of them is wrong in this genome. The rules test specific
pairs for themselves (the Hi-C gate, the arm fold, atline.validate). This table reports every pair the same way, with
agreement strictly defined (the intervals overlap, or lie within one Hi-C bin when one of them is a Hi-C readout) and chance
from each chromosome's own geometry: the share of placements of the two widths that would agree. It changes no call."""
from __future__ import annotations

import pandas as pd

from .params import GATE_MIN_REF, GATE_P
from .util import gap, poisson_binomial_sf

COLS = ["line_a", "line_b", "chromosomes", "agree", "expected", "p", "verdict"]
HIC_LINES = {"hic_anchored", "hic_prior_free", "fold", "alignment", "cis"}


def _iv(a, b):
    try:
        return int(a) - 1, int(b)
    except (TypeError, ValueError):
        return None


def line_calls(r) -> dict:
    """Where each line puts the centromere on one chromosome (a row of the centromeres table): {line: (start, end)}."""
    over = _iv(getattr(r, "seq_overruled_start", "."), getattr(r, "seq_overruled_end", "."))
    sat = over or (_iv(r.seq_start, r.seq_end) if r.seq_evidence in ("satellite", "weak_satellite") else None)
    out = dict(satellite=sat, ltr=_iv(r.te_start, r.te_end),
               hic_anchored=_iv(r.hic_anch_peak_start, r.hic_anch_peak_end),
               hic_prior_free=_iv(r.hic_denovo_peak_start, r.hic_denovo_peak_end),
               fold=_iv(r.fold_core_start, r.fold_core_end), alignment=_iv(r.align_core_start, r.align_core_end),
               at_trough=_iv(getattr(r, "at_start", "."), getattr(r, "at_end", ".")),
               repeats=_iv(getattr(r, "rep_start", "."), getattr(r, "rep_end", ".")),
               families=_iv(getattr(r, "fam_start", "."), getattr(r, "fam_end", ".")),
               cis=_iv(getattr(r, "cis_peak_start", "."), getattr(r, "cis_peak_end", ".")))
    return {k: v for k, v in out.items() if v is not None}


def table(t: pd.DataFrame, hic_bin: int | None) -> pd.DataFrame:
    """The pairs with at least one chromosome where both lines have a locus."""
    rows = list(t.itertuples(index=False))
    calls = [(int(r.chrom_len), line_calls(r)) for r in rows]
    lines = sorted({k for _, d in calls for k in d})
    out = []
    for i, a in enumerate(lines):
        for b in lines[i + 1:]:
            if {a, b} == {"cis", "ltr"}:     # the LTR cores chose the cis signal's sign, and hicfuse.cis_line tested the two
                continue
            tol = (hic_bin or 0) if {a, b} & HIC_LINES else 0
            k, ch = 0, []
            for L, d in calls:
                if a in d and b in d:
                    (s0, e0), (s1, e1) = d[a], d[b]
                    k += gap(s0, e0, s1, e1) <= tol
                    ch.append(min(1.0, ((e0 - s0) + (e1 - s1) + 2 * tol) / L))
            if not ch:
                continue
            n, p = len(ch), poisson_binomial_sf(k, ch)
            verdict = ("too few chromosomes" if n < GATE_MIN_REF else
                       "agree on most, beyond chance" if p <= GATE_P and 2 * k > n else
                       "agree beyond chance, on a minority" if p <= GATE_P and k >= GATE_MIN_REF else
                       "no agreement beyond chance")
            out.append(dict(line_a=a, line_b=b, chromosomes=n, agree=k, expected=round(sum(ch), 2), p=float(f"{p:.3g}"),
                            verdict=verdict))
    return pd.DataFrame(out, columns=COLS)
