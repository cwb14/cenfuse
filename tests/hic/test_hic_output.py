from __future__ import annotations

import numpy as np
import pandas as pd

from cenfuse.hic.bins import Bins
from cenfuse.hic.output import CALL_COLS, write_bed, write_bedgraph, write_tsv
from cenfuse.hic.plot import overview


def rows():
    base = {c: "" for c in CALL_COLS}
    a = dict(base, chrom="a", length=1000, callable=True, bin_size=100, peak_start=300, peak_end=400,
             core_start=200, core_end=500, peak_over_median=4.2, confidence="high", flags="")
    b = dict(base, chrom="b", length=900, callable=True, bin_size=100, peak_over_median=1.1, confidence="none",
             flags="no_rabl")
    return [a, b]


def test_bed9_only_called_rows(tmp_path):
    p = tmp_path / "x.bed"
    write_bed(rows(), p, 100)
    lines = p.read_text().splitlines()
    assert lines[0].startswith("#cenfuse-hic") and "bin_size=100" in lines[0]
    f = lines[1].split("\t")
    assert len(lines) == 2 and f[:6] == ["a", "200", "500", "a_cen", "420", "."] and f[6:8] == ["300", "400"]
    assert len(f) == 9


def test_tsv_keeps_every_row(tmp_path):
    p = tmp_path / "x.tsv"
    write_tsv(rows(), p)
    t = pd.read_csv(p, sep="\t")
    assert list(t.columns) == CALL_COLS and len(t) == 2 and t.loc[1, "confidence"] == "none"


def test_bedgraph(tmp_path):
    b = Bins({"a": 300, "b": 200}, 100)
    p = tmp_path / "s.bedgraph"
    write_bedgraph(np.array([1.0, np.nan, 2.0, 3.0, 4.0]), b, np.array([True, False]), p, "score")
    lines = [l for l in p.read_text().splitlines() if not l.startswith("track")]
    assert lines == ["a\t0\t100\t1", "a\t200\t300\t2"]


def test_overview_writes_pdf(tmp_path):
    b = Bins({"a": 1000, "b": 900}, 100)
    S = np.linspace(0.5, 3, b.N)
    overview(S, b, np.array([True, True]), rows(), {0: 350.0}, dict(n=2, negative=2, p=0.25, passed=False), tmp_path / "o.pdf")
    assert (tmp_path / "o.pdf").stat().st_size > 1000


def _bed_lines(path):
    return [l for l in path.read_text().splitlines() if not l.startswith("#")]


def test_bed_nan_ratio_scores_zero(tmp_path):
    r = rows()
    r[0] = dict(r[0], peak_over_median=float("nan"))
    p = tmp_path / "nan.bed"
    write_bed(r, p, 100)
    lines = _bed_lines(p)
    assert len(lines) == 1 and lines[0].split("\t")[4] == "0"


def test_bed_rejects_bad_labels_and_missing_core(tmp_path):
    r = rows()
    bad = [dict(r[0], confidence=float("nan")),
           dict(r[0], confidence="weird"),
           dict(r[0], confidence="medium", core_start=float("nan"))]
    p = tmp_path / "bad.bed"
    write_bed(bad, p, 100)
    assert _bed_lines(p) == []


def test_bed_medium_colour_and_score_cap(tmp_path):
    r = rows()
    r[0] = dict(r[0], confidence="medium", peak_over_median=12.3)
    p = tmp_path / "med.bed"
    write_bed(r, p, 100)
    f = _bed_lines(p)[0].split("\t")
    assert f[4] == "1000" and f[8] == "86,180,233"


def test_bed_integer_coordinates(tmp_path):
    r = rows()
    r[0] = dict(r[0], core_start=200.0, core_end=500.0, peak_start=300.0, peak_end=400.0)
    p = tmp_path / "float.bed"
    write_bed(r, p, 100)
    f = _bed_lines(p)[0].split("\t")
    assert f[1] == "200" and f[2] == "500" and f[6] == "300" and f[7] == "400"


def test_tsv_float_coordinate_kept_exact(tmp_path):
    r = rows()
    r[0] = dict(r[0], length=12345678.0, peak_start=12345678.0)
    p = tmp_path / "exact.tsv"
    write_tsv(r, p)
    text = p.read_text()
    cols = text.splitlines()[0].split("\t")
    vals = text.splitlines()[1].split("\t")
    assert vals[cols.index("length")] == "12345678"
    assert vals[cols.index("peak_start")] == "12345678"
    assert vals[cols.index("peak_over_median")] == "4.2"
    assert "e+" not in text


def test_overview_no_members_writes_nothing(tmp_path):
    b = Bins({"a": 1000, "b": 900}, 100)
    out = tmp_path / "none.pdf"
    overview(np.linspace(0.5, 3, b.N), b, np.array([False, False]), rows(), {}, dict(n=0, negative=0, p=1.0, passed=False), out)
    assert not out.exists()


def test_overview_all_nan_members_writes_nothing(tmp_path):
    b = Bins({"a": 1000, "b": 900}, 100)
    out = tmp_path / "nan.pdf"
    overview(np.full(b.N, np.nan), b, np.array([True, True]), rows(), {}, dict(n=2, negative=2, p=0.25, passed=False), out)
    assert not out.exists()
