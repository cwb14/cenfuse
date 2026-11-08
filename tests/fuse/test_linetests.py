"""linetests.py (future.md #65): every pair of lines, agreement against each chromosome's own chance; it changes no call."""
import pandas as pd

from cenfuse import linetests

L = 100_000_000


def _row(sat, te, anch=None, evidence="satellite"):
    f = lambda iv: (str(iv[0] + 1), str(iv[1])) if iv else (".", ".")
    (ss, se), (ts, te_), (hs, he) = f(sat), f(te), f(anch)
    return dict(chrom_len=L, seq_start=ss, seq_end=se, seq_evidence=evidence, te_start=ts, te_end=te_,
                hic_anch_peak_start=hs, hic_anch_peak_end=he, hic_denovo_peak_start=".", hic_denovo_peak_end=".",
                fold_core_start=".", fold_core_end=".", align_core_start=".", align_core_end=".")


def test_lines_that_agree_on_most_chromosomes_beyond_chance():
    rows = [_row((k * 1_000_000, k * 1_000_000 + 500_000), (k * 1_000_000 + 100_000, k * 1_000_000 + 300_000),
                 (k * 1_000_000 + 600_000, k * 1_000_000 + 850_000)) for k in range(10, 20)]
    t = linetests.table(pd.DataFrame(rows), hic_bin=250_000)
    sat_ltr = t[(t.line_a == "ltr") & (t.line_b == "satellite")].iloc[0]
    assert sat_ltr.agree == 10 and sat_ltr.verdict == "agree on most, beyond chance"
    hic = t[(t.line_a == "hic_anchored") & (t.line_b == "satellite")].iloc[0]
    assert hic.agree == 10                       # 100 kb apart: within the one Hi-C bin a Hi-C line is allowed


def test_disagreement_and_too_few():
    rows = [_row((10_000_000, 11_000_000), (60_000_000, 61_000_000)) for _ in range(8)]
    t = linetests.table(pd.DataFrame(rows), None)
    assert t.iloc[0].agree == 0 and t.iloc[0].verdict == "no agreement beyond chance"
    t2 = linetests.table(pd.DataFrame(rows[:2]), None)
    assert t2.iloc[0].verdict == "too few chromosomes"


def test_an_overruled_array_is_the_satellite_line():
    r = pd.DataFrame([dict(_row((1, 10), (500, 600), evidence="te_core"), seq_overruled_start="101", seq_overruled_end="200")])
    d = linetests.line_calls(next(r.itertuples(index=False)))
    assert d["satellite"] == (100, 200) and d["ltr"] == (500, 600)
