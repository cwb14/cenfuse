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


def test_the_repeat_lines_take_chance_over_their_span_and_drop_empty_family_windows():
    """0.19: the repeat-density and repeat-family windows skip the terminal eighths, so a pair with either takes chance over
    the chromosome outside them, as their own tests do (repline.with_hic, famline.with_repeats); a family window holding fewer
    than famline.MIN_COPIES family K-mers is no readout. Other pairs keep the whole length."""
    from cenfuse.famline import MIN_COPIES
    from cenfuse.hicfuse import END_ZONE
    rows = [dict(_row((k * 1_000_000, k * 1_000_000 + 500_000), (k * 1_000_000 + 100_000, k * 1_000_000 + 300_000)),
                 rep_start=str(k * 1_000_000 + 1), rep_end=str(k * 1_000_000 + 400_000),
                 fam_start=str(k * 1_000_000 + 1), fam_end=str(k * 1_000_000 + 400_000),
                 fam_votes=str(MIN_COPIES if k < 16 else MIN_COPIES - 1)) for k in range(10, 20)]
    t = linetests.table(pd.DataFrame(rows), None).set_index(["line_a", "line_b"])
    whole = (500_000 + 200_000) / L
    assert t.loc[("ltr", "satellite"), "expected"] == round(10 * whole, 2)
    assert t.loc[("repeats", "satellite"), "expected"] == round(10 * (500_000 + 400_000) / ((1 - 2 * END_ZONE) * L), 2)
    assert t.loc[("families", "satellite"), "chromosomes"] == 6 and t.loc[("families", "repeats"), "chromosomes"] == 6
    d = linetests.line_calls(next(pd.DataFrame([dict(rows[0], fam_votes=".")]).itertuples(index=False)))
    assert "families" not in d and "repeats" in d
    d = linetests.line_calls(next(pd.DataFrame([{k: v for k, v in rows[0].items() if k != "fam_votes"}]).itertuples(index=False)))
    assert "families" in d                                     # a table without the column: every window, as before
