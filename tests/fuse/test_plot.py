"""Overview figure (cenfuse.plot): the genome pages (A every chromosome to one scale, B each signal's role in each call,
C each signal's genome test, D each signal around every centromere against random places), the plants page (F, G) and
the two legends in plain words; the pages of a genome are drawn, saved and closed one at a time. Synthetic genomes:
chromosomes c1..cn of 10 Mb, each called on its satellite array at 4-5 Mb."""
import math
import re
import warnings
from types import SimpleNamespace

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.figure import Figure
from matplotlib.text import Text

from cenfuse import closeup, plot
from cenfuse import signals as sg

MM = 1 / 25.4
# a tomato-like summary row: every Hi-C readout ran, the satellite test failed, the repeat-family test is a lead
SUMMARY = dict(lines="satellite,te,hic", sat_one_per_chromosome="2/3 chromosomes hold the 3 largest clusters (null 2, "
               "P 0.87): fail", te_marker="CRM", te_decisive="True", te_verdict="CRM seed x, 32.7-fold in its cores, "
               "significant signature cores on 3/3 chromosomes; decisive", rep_gate="fail", rep_p="nan",
               fam_gate="pass", fam_p="0.0004", hic_gate_pass="0", hic_gate_support="0", hic_gate_untested="0",
               hic_chromosomes="3", cis_line="peak at the core on 2/3 chromosomes (P 2.4e-05): pass", fold_gate="fail",
               fold_on_denovo="0/3", fold_p="1", align_gate="fail", align_on_denovo="0/3", align_p="1")
FASTA = dict(lines="satellite", sat_one_per_chromosome="3/3 chromosomes hold the 3 largest clusters (null 1, P 0.001): "
             "pass")


def _fig(rdna=None, n=1, summary=None, te=False, hic=False, L=10_000_000, **by_chrom):
    """The figure dict `overview` receives, for n chromosomes (c1..cn) of length L, each called on one satellite array
    at 4-5 Mb (1-based closed in the table); by_chrom['c2'] = {...} overrides that chromosome's columns. With te, 40 CRM
    elements around each call and 20 spread over the chromosome; with hic, a 100-kb Hi-C score peaking at each call."""
    names = [f"c{i}" for i in range(1, n + 1)]
    recs = []
    for c in names:
        rec = dict(chrom=c, chrom_len=str(L), start="4000001", end="5000000", seq_start="4000001", seq_end="5000000",
                   seq_evidence="satellite", relation="sequence_only", confidence="high", call_source="sequence",
                   support="satellite", against=".", rule="sequence", flags=".", te_start=".", te_end=".",
                   hic_anch_core_start=".", hic_anch_core_end=".", hic_anch_peak_start=".", hic_anch_peak_end=".")
        rec.update(by_chrom.get(c, {}))
        recs.append(rec)
    cols = set().union(*recs)
    calls = pd.DataFrame([{**{k: "." for k in cols}, **r} for r in recs])
    arrays = pd.DataFrame(dict(chrom=names, start=[4_100_000] * n, end=[4_900_000] * n, role=["centromeric"] * n))
    f = dict(genome=SimpleNamespace(name="g", lines=["satellite"]), seq=SimpleNamespace(arrays=arrays), tel={},
             hic=None, meth=None, calls=calls, telo={}, plants=None, out="x", rdna=rdna or dict(nors=[], five_s=[]))
    if summary is not None:
        f["summary"] = summary
    lens = {c: L for c in names}
    if te:
        mids = {}
        for r in recs:
            c = (int(r["start"]) + int(r["end"])) // 2 if r["start"] != "." else L // 2
            m = np.concatenate([np.linspace(c - 400_000, c + 400_000, 40), np.linspace(0, L - 1, 20)])
            mids[r["chrom"]] = np.sort(np.clip(m, 0, L - 1).astype(np.int64))
        f["tel"] = dict(marker=("CRM", {}, True), mids={"CRM": mids}, signature=None,
                        density={c: (100_000, np.bincount(m // 100_000, minlength=L // 100_000 + 1))
                                 for c, m in mids.items()})
    if hic:
        from cenfuse.hic.bins import Bins
        size = 100_000
        b = Bins(lens, size)
        S = np.ones(b.N)
        for k, r in enumerate(recs):
            if r["start"] != ".":
                S[b.off[k] + int(r["start"]) // size] = 4.0
        f["hic"] = dict(bins=b, members=np.ones(n, bool), S=S, denovo_S=S.copy(), bin_size=size)
    return f


def _rows(f):
    hic_bin = f["hic"]["bins"].size if f["hic"] is not None else None
    return sg.rows(f["calls"], f.get("summary"), hic_bin)


def _texts(fig):
    return [t for t in fig.findobj(Text) if t.get_text()]


def _gid(fig, prefix):
    return [a for a in fig.findobj() if str(a.get_gid()).startswith(prefix)]


def _ax(fig, gid):
    return next(a for a in fig.axes if a.get_gid() == gid)


def _plants(f):
    """Five plants scored on f's chromosomes (docs/plants.md): every verdict of the plants page, p1 the reference plant
    and the others other plants, a CHG track per plant and chromosome (methylation bins) with a stretch without
    data, and the sites scored."""
    from cenfuse.meth import params as mp
    chroms, L = list(f["calls"].chrom), int(f["calls"].chrom_len.iat[0])
    plants = ["p1", "p2", "p3", "p4", "p5"]
    kinds = ["shared", "shifted", "possible_shift", "call_absent", "conflict", "unclear", "no_data"]
    recs = []
    for i, p in enumerate(plants):
        for j, c in enumerate(chroms):
            v = kinds[(i + 2 * j) % len(kinds)]
            far = v in ("shifted", "possible_shift", "conflict")
            recs.append(dict(plant=p, chrom=c, verdict=v, shift_bp="2500000" if far else ".",
                             alt_start=7_000_000 if far else -1, alt_end=7_200_000 if far else -1))
    f["plants"] = pd.DataFrame(recs)
    f["verdicts"] = {("meth", p): {"verdict": "reference" if p == "p1" else "other"} for p in plants}
    rng = np.random.default_rng(0)
    adj = {}
    for p in plants:
        adj[p] = {}
        for c in chroms:
            y = rng.normal(0, 0.02, L // mp.BIN)
            y[10:20] = np.nan
            adj[p][c] = y
    f["meth"] = dict(adj=adj, samples=plants)
    f["psites"] = pd.DataFrame([dict(chrom=c, sample=p, mid=m, test="loss" if m == 4_500_000 else "gain")
                                for c in chroms for p in plants for m in (4_500_000, 7_100_000)])
    f["libraries"] = None
    return f


def _rich():
    """A genome with every mark of the genome pages and a plants page: LTR elements, Hi-C (without the Hi-C table a
    pipeline run adds), a Hi-C call where the LTR cluster points elsewhere, a chromosome without a call, telomeres,
    rDNA and five plants."""
    f = _fig(dict(nors=[("c1", 100_000, 400_000, 30, 9_100)], five_s=[("c2", 8_000_000, 8_200_000, 80, 500)]),
             n=3, te=True, hic=True, summary=SUMMARY,
             c2=dict(support="hic", against="te", call_source="hic", rule="cis_only", flags="cis_only",
                     confidence="low", start="8000001", end="8250000", cis_peak_start="8000001",
                     cis_peak_end="8250000", te_start="4000001", te_end="5000000", seq_evidence="retrotransposon"),
             c3=dict(start=".", end=".", relation="none", confidence="none", call_source="none", rule="none",
                     support="."))
    f["telo"] = dict(unit="AAACCCT", ends=[("c1", "start", 0, 5_000), ("c1", "end", 9_995_000, 10_000_000)],
                     interstitial=[("c2", 6_000_000, 6_002_000, 30)])
    for key in ("S", "denovo_S"):                             # no Hi-C reads map inside c1's satellite array
        f["hic"][key][42:48] = np.nan
    return _plants(f)


def _paragraphs(text):
    return text.split("\n\n")


def _drawn(f):
    """f as `overview` hands it to `legends`: the genome pages drawn, then the plants page (it sets the chromosome
    panel G shows and G's methylation scale); every page closed."""
    for fig in plot._genome_pages(f, _rows(f)):
        plt.close(fig)
    pp = plot._plants_page(f)
    if pp is not None:
        plt.close(pp)
    return f


def test_short_legend_has_no_banned_words():
    """The short legend is for a biologist: the result first (how many chromosomes have a centromere, at which
    confidence, which signals agree, which point elsewhere), then one paragraph per panel letter in order (A-E, the
    plants page F-G), saying why Hi-C is blank inside satellite arrays; no word of cenfuse's jargon, no code names, no
    tool names."""
    short = plot.legends(_drawn(_rich()))[0]
    pars = _paragraphs(short)
    assert pars[0].startswith("In g, 2 of 3 chromosomes have a centromere (1 high and 1 low confidence). ")
    assert "Hi-C alone places 1" in pars[0] and "The CRM cluster points elsewhere on 1" in pars[0]
    letters = [re.match(r"\*\*([A-Z])\*\*, ", p)[1] for p in pars[1:]]
    assert letters == list("ABCDEFG")
    assert "too few reads map" in pars[1]
    hits = sg.BANNED.findall(short)
    assert not hits, hits
    assert "_" not in short and "cenfuse" not in short.lower() and "tsv" not in short


def test_short_legend_opens_with_the_result():
    """The opening sentence counts the calls by confidence and groups them by the signals that agree; a genome
    without a call is told so, and the legend names no centromere mark."""
    f = _fig(n=3, summary=FASTA)
    assert plot.legends(f)[0].split("\n\n")[0] == ("In g, 3 of 3 chromosomes have a centromere (all high "
                                                   "confidence). The satellite array alone places all 3.")
    f = _fig(n=1, summary=FASTA, c1=dict(start=".", end=".", relation="none", confidence="none", call_source="none",
                                          rule="none", support="."))
    short = plot.legends(f)[0]
    assert short.split("\n\n")[0] == "In g, no centromere was called on its 1 chromosome."
    assert "triangle" not in short and "band through every track" not in short      # neither is drawn


def test_the_long_legend_names_every_rule_that_places_a_call():
    """The long legend explains how cenfuse places a centromere in plain words: the roles (used, points elsewhere, not
    used, not used for methylation from other plants), the confidence styles, and every kind of verdict, naming each
    value of the `rule` column that docs/outputs.md lists, in parentheses: the only place code names appear. No word
    of the jargon outside those names."""
    from pathlib import Path
    doc = (Path(__file__).parents[2] / "docs" / "outputs.md").read_text()
    para = doc.split("**`rule`.**", 1)[1].split("\n\n", 1)[0]
    rules = re.findall(r"`([a-z_]+)`(?=\s\(|,)", para)
    assert len(rules) >= 31
    long = plot.legends(_drawn(_rich()))[1]
    body = "\n\n".join(_paragraphs(long)[:-1])               # the last paragraph: where the analysis ran
    inside = " ".join(re.findall(r"\(([^()]*)\)", body))
    for r in rules:
        assert re.search(rf"(?<![\w]){r}(?![\w])", inside), r
    outside = re.sub(r"\([^()]*\)", "", body)
    assert not re.search(r"\b[a-z]+_[a-z_]+\b", outside), re.findall(r"\b[a-z]+_[a-z_]+\b", outside)
    for x in ("*used*", "*points elsewhere*", "*not used*", "*not used (other plants)*", "black", "gray", "open",
              "dashed"):
        assert x in long.lower(), x
    hits = sg.BANNED.findall(long)
    assert not hits, hits
    assert _paragraphs(long)[-1].startswith("Analysis: ")


def test_legends_without_the_hic_table():
    """A figure dict whose Hi-C lacks the per-chromosome table (synthetic fixtures; the pipeline adds it) still gets
    both legends, the long one naming the Hi-C bin."""
    f = _fig(n=2, te=True, hic=True, summary=SUMMARY)
    assert "table" not in f["hic"]
    short, long = plot.legends(f)
    assert short and "100-kb" in long


def test_legend_pages_fit_the_page():
    """Legend pages are 180 mm wide and at most 240 mm tall at 10 pt: a text too long for one page goes on to the
    next, between paragraphs, every word once and in order; **x** is set bold and *x* italic, the marks dropped."""
    par = ("**A**, The satellite arrays of *Poa annua* are drawn to one scale, " + "word " * 160).strip()
    text = "\n\n".join([par] * 8)
    figs = list(plot._text_pages(text))
    try:
        assert len(figs) >= 2
        words = []
        for fig in figs:
            w, h = fig.get_size_inches() / MM
            assert w == pytest.approx(180) and h <= 240 + 1e-6
            texts = _texts(fig)
            assert all(t.get_fontsize() == 10 for t in texts)
            assert not any("*" in t.get_text() for t in texts)
            lines = {}                                        # a line's runs carry their own spaces
            for t in sorted(texts, key=lambda t: t.get_position()[0]):
                lines.setdefault(round(t.get_position()[1], 6), []).append(t.get_text())
            words += " ".join("".join(lines[y]) for y in sorted(lines, reverse=True)).split()
            bold = {t.get_text().strip() for t in texts if t.get_fontweight() == "bold"}
            italic = {t.get_text().strip() for t in texts if t.get_fontstyle() == "italic"}
            assert bold <= {"A"} and italic <= {"Poa annua", "Poa", "annua"}
        assert words == text.replace("*", "").split()
        assert {t.get_text() for t in _texts(figs[0]) if t.get_fontweight() == "bold"} == {"A"}
    finally:
        for fig in figs:
            plt.close(fig)


def test_a_word_longer_than_a_line_stays_on_the_page():
    """A path longer than a line (the provenance line of a deep run directory) is cut after its slashes, filling the
    line it starts on, and every piece of text ends inside the page's margins."""
    path = "/" + "/".join(f"v{i}" for i in range(1, 60)) + "/runs/tomato"
    text = f"Analysis: host:{path} · Script: cenfuse/plot.py · 2026-11-16"
    figs = list(plot._text_pages(text))
    try:
        assert len(figs) == 1
        fig = figs[0]
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        right = fig.bbox.x1 - plot.LEGEND_SIDE * MM * fig.dpi
        texts = _texts(fig)
        assert all(t.get_window_extent(renderer).x1 <= right + 0.5 for t in texts)
        lines = {}
        for t in sorted(texts, key=lambda t: t.get_position()[0]):
            lines.setdefault(round(t.get_position()[1], 6), []).append(t.get_text())
        lines = ["".join(lines[y]) for y in sorted(lines, reverse=True)]
        assert len(lines) >= 2 and lines[0].startswith("Analysis: host:/v1/")      # the path starts on line one
        assert "".join(lines).replace(" ", "") == text.replace(" ", "")
        assert all(x.endswith("/") or x.endswith(" ") or x == lines[-1] for x in lines)
    finally:
        for fig in figs:
            plt.close(fig)


def test_plants_page_has_no_banned_words():
    """The plants page keeps its layout in plain words: letters F and G after the genome pages' A-E, no word of the
    jargon in its key, labels or verdicts, every text on the page."""
    f = _rich()
    fig = plot._plants_page(f)                                # drawn on its own: it needs no genome page first
    try:
        texts = [t.get_text() for t in _texts(fig)]
        letters = sorted(t.get_text() for t in fig.texts if t.get_fontweight() == "bold")
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        out = [t.get_text() for t in _texts(fig) if t.get_window_extent(renderer).x1 > fig.bbox.x1
               or t.get_window_extent(renderer).x0 < 0]
    finally:
        plt.close(fig)
    assert not out, out                                       # the longer verdict words stay on the page
    assert letters == ["F", "G"]
    assert "centromere" in texts and "CRM elements" in texts
    hits = [x for x in texts if sg.BANNED.search(x)]
    assert not hits, hits


def test_overview_marks_45s_and_5s_rdna_arrays():
    """45S arrays (NORs) as filled and 5S arrays as open gray diamonds at their midpoints, in the key and the legend;
    a genome without either says so."""
    from matplotlib.colors import to_hex
    f = _fig(dict(nors=[("c1", 100_000, 400_000, 30, 9_100)], five_s=[("c1", 8_000_000, 8_200_000, 80, 500)]))
    fig = plot._pages(f)[0]
    marks = sorted((float(x.get_xdata()[0]), to_hex(x.get_markerfacecolor())) for x in _ax(fig, "A").lines
                   if x.get_marker() == "D")
    assert marks == [(250_000.0, to_hex(plot.RDNA_MARK)), (8_100_000.0, "#ffffff")]
    labels = [t.get_text() for t in fig.legends[0].get_texts()]
    assert "45S rDNA array (NOR)" in labels and "5S rDNA array" in labels
    plt.close(fig)
    assert "1 array of 45S rDNA and 1 of 5S rDNA" in plot.legends(f)[1]
    f = _fig(dict(nors=[], five_s=[]))
    fig = plot._pages(f)[0]
    assert not [x for x in _ax(fig, "A").lines if x.get_marker() == "D"]
    assert not {"45S rDNA array (NOR)", "5S rDNA array"} & {t.get_text() for t in fig.legends[0].get_texts()}
    plt.close(fig)
    assert "no 45S or 5S rDNA array" in plot.legends(f)[1]


def test_long_chromosome_names_stay_on_the_page():
    """Panel A's names are drawn whole, however long (a GWH accession has 17 characters; the left margin held about
    15), and end left of the chromosomes."""
    long = {f"c{i}": f"GWHERKE0000000{i}.1" for i in range(1, 5)}
    f = _fig(n=4)
    f["calls"]["chrom"] = f["calls"]["chrom"].map(long)
    f["seq"].arrays["chrom"] = f["seq"].arrays["chrom"].map(long)
    fig = plot._pages(f)[0]
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    a = _ax(fig, "A").get_window_extent(renderer)
    names = [t for t in _texts(fig) if t.get_text().startswith("GWH")]
    assert len(names) == 4
    for t in names:
        box = t.get_window_extent(renderer)
        assert box.x0 >= 0 and box.y0 >= 0 and box.x1 <= fig.bbox.x1 and box.y1 <= fig.bbox.y1, t.get_text()
        assert box.x1 <= a.x0 + 1, t.get_text()                   # left of the chromosome, which starts at 0
    plt.close("all")


def test_genome_pages_paginate_at_24_rows():
    """60 chromosomes: three pages of panels A and B (24, 24 and 12 rows), every name once, every page 180 mm wide and
    at most 240 mm tall; C and D once, under the last rows or on a page of their own."""
    f = _fig(n=60)
    rows = _rows(f)
    figs = list(plot._genome_pages(f, rows))
    try:
        a_pages = [fig for fig in figs if any(ax.get_gid() == "A" for ax in fig.axes)]
        assert len(a_pages) == 3
        per_page = [len(_gid(fig, "name:")) for fig in a_pages]
        assert per_page == [24, 24, 12]
        names = [t.get_text() for fig in a_pages for t in _gid(fig, "name:")]
        assert sorted(names) == sorted(r.chrom for r in rows)
        for fig in figs:
            w, h = fig.get_size_inches() / MM
            assert w == pytest.approx(180) and h <= 240 + 1e-6
        assert sum(len(_gid(fig, "letter:C")) for fig in figs) == 1
        assert sum(len(_gid(fig, "letter:D")) for fig in figs) == 1
        assert len(figs) in (3, 4)
    finally:
        for fig in figs:
            plt.close(fig)


def test_close_up_boxes_mark_the_close_up_windows():
    """Panel A's gray boxes are the close-up windows (closeup.window), cut at the chromosome ends; a centromere is one
    triangle, and a trusted signal that points elsewhere is an x at its own place. The repeat-family cluster, a lead
    only, points elsewhere in B but gets no x in A."""
    f = _fig(n=3, te=True, hic=True, summary=SUMMARY,
             c2=dict(support="hic", against="te", call_source="hic", rule="cis_only", flags="cis_only",
                     confidence="low", start="8000001", end="8250000", cis_peak_start="8000001",
                     cis_peak_end="8250000", te_start="4000001", te_end="5000000", seq_evidence="retrotransposon"),
             c3=dict(start="200001", end="400000", seq_start="200001", seq_end="400000", against="families",
                     flags="fam_elsewhere", fam_start="7000001", fam_end="7500000"))
    rows = _rows(f)
    assert rows[2].roles["fam"] == "elsewhere"
    keys = plot._keys(f, rows)
    assert keys == sg.present(rows, SUMMARY, dict(te=True, hic=True, meth=False))
    _, wins = closeup.window(rows, keys, [r.L for r in rows])
    fig = next(iter(plot._genome_pages(f, rows)))
    boxes = {b.get_gid().split(":")[1]: (b.get_x(), b.get_x() + b.get_width()) for b in _gid(fig, "window:")}
    want = {r.chrom: (max(lo, 0), min(hi, r.L)) for r, (lo, hi) in zip(rows, wins)}
    assert boxes.keys() == want.keys()
    for c in want:
        assert boxes[c] == pytest.approx(want[c])
    assert boxes["c3"][0] == 0                                    # c3's window starts before the chromosome
    tri = {m.get_gid().split(":")[1]: (float(m.get_xdata()[0]), m.get_marker()) for m in _gid(fig, "centromere:")}
    assert tri == {"c1": (4_500_000.0, "^"), "c2": (8_125_000.0, "^"), "c3": (300_000.0, "^")}
    assert [x.get_gid() for x in _gid(fig, "elsewhere:")] == ["elsewhere:c2:te"]
    assert {m.get_gid(): m.get_marker() for m in _gid(fig, "role:c3:fam")} == {"role:c3:fam": "x"}
    assert float(_gid(fig, "elsewhere:")[0].get_xdata()[0]) == 4_500_000.0
    plt.close(fig)


def test_profiles_null_is_reproducible():
    """Panel D: two calls with seed 0 give identical random-place bands, another seed another band; every random window
    lies inside its chromosome (no missing value in the band), and the centromeres stand above it."""
    f = _fig(n=4, te=True, hic=True)
    rows = _rows(f)

    def run(seed):
        fig = Figure()
        axes = {k: fig.add_subplot(1, 3, i + 1) for i, k in enumerate(("hic", "te", "sat"))}
        return plot._profiles(f, rows, axes, seed=seed)

    a, b, c = run(0), run(0), run(1)
    assert set(a) == {"hic", "te", "sat"}
    for k in a:
        for part in ("lo", "med", "hi"):
            np.testing.assert_array_equal(a[k][part], b[k][part])
        assert np.isfinite(a[k]["null"]).all(), k             # no random window runs past a chromosome end
    assert any(not np.array_equal(a[k]["med"], c[k]["med"]) for k in a)
    mid = len(a["te"]["x"]) // 2
    assert a["te"]["mean"][mid] > a["te"]["hi"][mid] and a["sat"]["mean"][mid] > a["sat"]["hi"][mid]


def test_table_rows_match_signal_tests():
    """Panel C: one row per genome test (signals.signal_tests), in the order of every view (signals.TRACK_ORDER), its
    label, question and verdict as written there, P beside the question where the summary gives one; a P the summary
    leaves as 'nan' is left out, never printed."""
    f = _fig(n=3, te=True, hic=True, summary=SUMMARY)
    rows = _rows(f)
    figs = list(plot._genome_pages(f, rows))
    want = sorted(sg.signal_tests(SUMMARY, 3, "CRM"), key=lambda t: sg.TRACK_ORDER.index(t.key))
    texts = {t.get_gid(): t.get_text().replace("\n", " ") for fig in figs for t in _texts(fig) if t.get_gid()}
    labels = [t for g, t in texts.items() if g.startswith("test-label:")]
    assert labels == [t.label for t in want]
    for t in want:
        assert texts[f"test-question:{t.key}"].startswith(t.question.rstrip(")"))      # P inside '(chance: ...)'
        assert texts[f"test-verdict:{t.key}"] == t.verdict
    assert texts["test-question:flank"].endswith("(P = 2e-5)")
    assert texts["test-question:sat"].endswith("(chance: 2; P = 0.87)")
    assert texts["test-question:rep"] == next(t.question for t in want if t.key == "rep")      # P 'nan': none
    assert not [x for fig in figs for x in _texts(fig) if re.search(r"\bnan\b", x.get_text(), re.I)]
    for fig in figs:
        plt.close(fig)


def test_b_columns_are_the_close_up_tracks():
    """Panel B: one column per close-up track (signals.present, in TRACK_ORDER), the keys plot._keys gives the
    close-ups; a signal that ran but has no track (here the arm fold and alignment, which failed, and the repeat
    density and repeat-family cluster, never consulted) is in panel C only. Tomato-like: Hi-C contact with the other
    centromeres, the flanks, the LTR family and the satellite arrays."""
    f = _fig(n=3, te=True, hic=True, summary=SUMMARY,
             c1=dict(hic_anch_core_start="6000001", hic_anch_core_end="6100000"),
             c2=dict(support="hic", against="te", call_source="hic", rule="cis_only", flags="cis_only",
                     confidence="low", start="8000001", end="8250000", cis_peak_start="8000001",
                     cis_peak_end="8250000", te_start="4000001", te_end="5000000", seq_evidence="retrotransposon"))
    rows = _rows(f)
    keys = sg.present(rows, SUMMARY, dict(te=True, hic=True, meth=False))
    assert keys == plot._keys(f, rows) == ["hic", "flank", "te", "sat"]
    fig = next(iter(plot._genome_pages(f, rows)))
    heads = {t.get_gid().split(":")[1]: t.get_text() for t in _gid(fig, "column:")}
    assert list(heads) == keys
    assert list(heads.values()) == ["Hi-C", "Hi-C flanks", "CRM", "satellite"]
    assert {t.get_gid().split(":")[1] for t in _gid(fig, "test-label:")} > set(keys)      # C: every signal that ran
    plt.close(fig)


def test_fasta_only_genome_page():
    """A FASTA-only run, the most common input: one column in B (the satellite arrays), one row in C, one profile in D
    (centromeric satellite), no Hi-C trace and no LTR elements; one page."""
    f = _fig(n=3, summary=FASTA)
    rows = _rows(f)
    figs = list(plot._genome_pages(f, rows))
    try:
        assert len(figs) == 1
        fig = figs[0]
        assert [t.get_gid() for t in _gid(fig, "column:")] == ["column:sat"]
        assert [t.get_gid() for t in _gid(fig, "test-label:")] == ["test-label:sat"]
        assert [ax.get_gid() for ax in fig.axes if str(ax.get_gid()).startswith("D:")] == ["D:sat"]
        roles = {g.get_gid(): g.get_marker() for g in _gid(fig, "role:")}
        assert roles == {f"role:c{i}:sat": "o" for i in (1, 2, 3)}
        assert not _gid(fig, "hic:") and not _gid(fig, "te:")
        assert sorted(t.get_text() for t in _gid(fig, "letter:")) == ["A", "B", "C", "D"]
    finally:
        for fig in figs:
            plt.close(fig)


def test_letters_run_in_reading_order():
    """Panel letters A, B, C, D in reading order: A left of B on the top row, then C, then D below."""
    f = _fig(n=3, te=True, hic=True, summary=SUMMARY)
    fig = next(iter(plot._genome_pages(f, _rows(f))))
    pos = {t.get_text(): t.get_position() for t in _gid(fig, "letter:")}
    assert pos["A"][1] == pytest.approx(pos["B"][1]) and pos["A"][0] < pos["B"][0]
    assert pos["B"][1] > pos["C"][1] > pos["D"][1] and pos["A"][0] == pytest.approx(pos["C"][0])
    plt.close(fig)


def test_no_banned_words_on_genome_pages():
    """No word of cenfuse's jargon on a genome page: names, column heads, the test table, axis labels and the key,
    with every Hi-C readout, a lead, telomeres and rDNA drawn."""
    f = _fig(dict(nors=[("c1", 100_000, 400_000, 30, 9_100)], five_s=[("c2", 8_000_000, 8_200_000, 80, 500)]),
             n=3, te=True, hic=True, summary=SUMMARY,
             c2=dict(support="hic", against="te", call_source="hic", rule="cis_only", flags="cis_only",
                     confidence="low", start="8000001", end="8250000", cis_peak_start="8000001",
                     cis_peak_end="8250000", te_start="4000001", te_end="5000000", seq_evidence="retrotransposon"),
             c3=dict(start=".", end=".", relation="none", confidence="none", call_source="none", rule="none",
                     support="."))
    f["telo"] = dict(unit="AAACCCT", ends=[("c1", "start", 0, 5_000), ("c1", "end", 9_995_000, 10_000_000)],
                     interstitial=[("c2", 6_000_000, 6_002_000, 30)])
    figs = list(plot._genome_pages(f, _rows(f)))
    texts = [t.get_text() for fig in figs for t in _texts(fig)]
    for fig in figs:
        plt.close(fig)
    assert "telomere array" in texts and "no call" in texts
    hits = [x for x in texts if sg.BANNED.search(x)]
    assert not hits, hits


@pytest.fixture
def quick_pages(monkeypatch):
    """Pages of 2 rows: the page lifecycle does not depend on the number of rows a page holds."""
    monkeypatch.setattr(plot, "ROWS_PER_PAGE", 2)


def test_overview_keeps_one_page_open_at_a_time(tmp_path, monkeypatch, quick_pages):
    """A genome of more pages than matplotlib's open-figure limit (20; the 670-sequence Cissus run had 28) is drawn,
    saved and closed one page at a time: no "More than 20 figures" RuntimeWarning (no warning at all), exactly one
    figure open at every save (a page is closed before the next is drawn), none left, and every page in the PDF."""
    n_genome, n_legend_pages = 25, 2                          # the genome pages; the short and the long legend
    f = _fig(n=n_genome * plot.ROWS_PER_PAGE)
    rows = _rows(f)
    n_closeup = len(list(closeup.pages(f, rows, plot._keys(f, rows))))
    for fig in closeup.pages(f, rows, plot._keys(f, rows)):
        plt.close(fig)
    n_pages = n_genome + n_closeup
    open_at_save = []

    def spy(save):
        def wrapped(self, *args, **kwargs):
            open_at_save.append(len(plt.get_fignums()))
            return save(self, *args, **kwargs)
        return wrapped

    monkeypatch.setattr(Figure, "savefig", spy(Figure.savefig))
    monkeypatch.setattr(PdfPages, "savefig", spy(PdfPages.savefig))
    plt.close("all")                                          # figures other tests left open are not this test's
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            plot.overview([f], tmp_path)
        left = plt.get_fignums()
    finally:
        plt.close("all")                                      # a leak fails this test, not the next one
    assert len(open_at_save) >= n_pages + n_legend_pages      # the spy saw every save
    # a PDF page is counted twice (PdfPages.savefig calls Figure.savefig), both times with that page alone open
    assert max(open_at_save) == 1 and left == []
    pdf = (tmp_path / "cenfuse_overview.pdf").read_bytes()
    assert len(re.findall(rb"/Type\s*/Page\b", pdf)) == n_pages + n_legend_pages     # page objects, not /Pages
    assert (tmp_path / "g.overview.png").is_file()


@pytest.mark.parametrize("where", ["draw", "save"])
def test_overview_closes_the_page_that_fails(tmp_path, monkeypatch, quick_pages, where):
    """A page whose drawing or saving fails is closed before the error reaches the caller: the pipeline logs it and
    goes on, and a figure left open would stay so for the rest of the process."""
    f = _fig(n=2 * plot.ROWS_PER_PAGE + 1)                    # three pages
    n_calls = []

    def fail_second(real):
        def wrapped(*args, **kwargs):
            n_calls.append(1)
            if len(n_calls) == 2:
                raise OSError("disk full")
            return real(*args, **kwargs)
        return wrapped

    if where == "draw":
        monkeypatch.setattr(plot, "_panel_a", fail_second(plot._panel_a))
    else:
        monkeypatch.setattr(PdfPages, "savefig", fail_second(PdfPages.savefig))
    plt.close("all")
    try:
        with pytest.raises(OSError, match="disk full"):
            plot.overview([f], tmp_path)
        left = plt.get_fignums()
    finally:
        plt.close("all")
    assert left == []


def test_cpg_profile_where_the_figure_has_cpg():
    """Panel D draws CpG methylation around the centromeres where the figure dict carries the CpG profiles
    (cdr.prepare's chroms: a 5-kb methylation profile per chromosome): a dip at every call stands below random
    places."""
    from cenfuse import cdr
    f = _fig(n=3)
    chroms = {}
    for c in ("c1", "c2", "c3"):
        y = np.full(10_000_000 // cdr.BIN, 0.8)
        y[4_000_000 // cdr.BIN:5_000_000 // cdr.BIN] = 0.2
        chroms[c] = SimpleNamespace(y=y)
    f["cpg"] = dict(chroms=chroms)
    fig = Figure()
    axes = {"cpg": fig.add_subplot(1, 2, 1), "sat": fig.add_subplot(1, 2, 2)}
    prof = plot._profiles(f, _rows(f), axes, seed=0)
    mid = len(prof["cpg"]["x"]) // 2
    assert prof["cpg"]["mean"][mid] < prof["cpg"]["lo"][mid]
    assert math.isfinite(prof["cpg"]["mean"][mid])


def test_overview_writes_both_pngs_and_page_order(tmp_path):
    """Two chromosomes: the PDF holds the genome page, the close-up page and the two legend pages, in that order; both
    PNGs exist."""
    f = _fig(n=2, summary=FASTA)
    plt.close("all")
    plot.overview([f], tmp_path)
    assert (tmp_path / "g.overview.png").is_file() and (tmp_path / "g.closeups.png").is_file()
    pdf = (tmp_path / "cenfuse_overview.pdf").read_bytes()
    assert len(re.findall(rb"/Type\s*/Page\b", pdf)) == 4
    kinds = [k for k, fig in plot._page_iter(_fig(n=2, summary=FASTA)) if plt.close(fig) is None]
    assert kinds == ["genome", "closeup"]


def test_overview_without_summary(tmp_path):
    """A caller without `summary` (older callers) still draws: the roles come from support and against alone, and
    table C says "not reported"."""
    f = _fig(n=2)
    assert "summary" not in f
    plt.close("all")
    plot.overview([f], tmp_path)
    assert (tmp_path / "g.overview.png").is_file() and (tmp_path / "g.closeups.png").is_file()
    fig = plot._pages(_fig(n=2))[0]
    try:
        assert any("not reported" in t.get_text() for t in _texts(fig))
    finally:
        plt.close(fig)


def test_c_and_d_on_their_own_page_for_24_rows_with_hic():
    """24 rows with Hi-C: C and D do not fit under A, so they take a page of their own, every page within 240 mm, C
    and D drawn once."""
    f = _fig(n=24, summary=SUMMARY, hic=True)
    figs = list(plot._genome_pages(f, _rows(f)))
    try:
        assert len(figs) == 2
        assert sum(len(_gid(fig, "letter:C")) for fig in figs) == 1
        assert sum(len(_gid(fig, "letter:D")) for fig in figs) == 1
        assert not _gid(figs[0], "letter:C") and _gid(figs[1], "letter:C")
        for fig in figs:
            w, h = fig.get_size_inches() / MM
            assert w == pytest.approx(180) and h <= 240 + 1e-6
    finally:
        for fig in figs:
            plt.close(fig)


@pytest.mark.parametrize("kw", [dict(summary=FASTA), dict(), dict(summary=SUMMARY, te=True)])
def test_one_chromosome_genome_draws_every_page(tmp_path, kw):
    f = _fig(n=1, **kw)
    plt.close("all")
    plot.overview([f], tmp_path)
    assert (tmp_path / "g.overview.png").is_file() and (tmp_path / "g.closeups.png").is_file()


def test_one_chromosome_genome_with_hic_draws_its_figure_pages():
    f = _fig(n=1, summary=SUMMARY, te=True, hic=True)
    kinds = [k for k, fig in plot._page_iter(f) if plt.close(fig) is None]
    assert kinds == ["genome", "closeup"]


def test_plants_page_hic_scale_is_that_of_the_trace_it_draws():
    """Panel G draws the reference's contact with the other centromeres on that score's own genome-wide top (the
    99.5th percentile, as the genome pages scale their trace), not on the top of the hub the genome pages may draw; and
    G draws no box for where the DNA sequence alone puts the centromere."""
    f = _rich()
    f["hic"]["denovo_S"][::10] = 40.0                         # the hub peaks far higher than the other score
    hsc = closeup._hic_scale(f, list(f["calls"].chrom), "hic")
    assert hsc != pytest.approx(closeup._hic_scale(f, list(f["calls"].chrom), "hub"))
    fig = plot._plants_page(f)
    try:
        line = _gid(fig, "G:hic")[0]
        ax = line.axes
        lo, hi = ax.get_xlim()
        x, y = plot._track(f["hic"], f["focus"], "S")
        k = (x >= lo) & (x <= hi)
        n_plants = f["plants"].plant.nunique()
        np.testing.assert_allclose(line.get_ydata(), n_plants + 0.15 + 0.9 + np.clip(y[k] / hsc, 0, 1) * 0.55)
        boxes = [r for r in ax.patches if not r.get_fill() and matplotlib.colors.to_hex(r.get_edgecolor()) == "#000000"]
        assert not boxes
    finally:
        plt.close(fig)


def test_overview_of_a_genome_without_any_call(tmp_path):
    """No call on any chromosome (a holocentric genome), with LTR elements, Hi-C and a summary: every page is drawn
    (genome page without panel D, close-ups, both legends), both PNGs written, and the legend says no centromere was
    called."""
    none = dict(start=".", end=".", relation="none", confidence="none", call_source="none", rule="none", support=".")
    f = _fig(n=3, te=True, hic=True, summary=SUMMARY, c1=none, c2=none, c3=none)
    plt.close("all")
    plot.overview([f], tmp_path)
    assert plt.get_fignums() == []
    assert (tmp_path / "g.overview.png").is_file() and (tmp_path / "g.closeups.png").is_file()
    pdf = (tmp_path / "cenfuse_overview.pdf").read_bytes()
    assert len(re.findall(rb"/Type\s*/Page\b", pdf)) >= 4
    figs = [fig for _, fig in plot._page_iter(_fig(n=3, te=True, hic=True, summary=SUMMARY, c1=none, c2=none, c3=none))]
    try:
        assert not [a for fig in figs for a in fig.axes if str(a.get_gid()).startswith("D:")]
        assert not _gid(figs[0], "centromere:")
    finally:
        for fig in figs:
            plt.close(fig)
    short = plot.legends(f)[0]
    assert short.startswith("In g, no centromere was called on any of its 3 chromosomes.")
    assert "**D**" not in short and not sg.BANNED.search(short)


def test_panel_d_counts_the_centromeres_each_signal_placed():
    """Over each panel D plot, how many of its centromeres that signal placed (signals.Row.placed): the satellite
    placed the two calls kept on their arrays, Hi-C the one it placed (c2, by the flank contact), the LTR family none.
    A signal peaks at the centromeres it placed by construction, and the legend says so."""
    f = _fig(n=3, te=True, hic=True, summary=SUMMARY,
             c2=dict(support="hic", call_source="hic", rule="cis_only", flags="cis_only", confidence="low",
                     start="8000001", end="8250000", cis_peak_start="8000001", cis_peak_end="8250000",
                     seq_evidence="retrotransposon"))
    rows = _rows(f)
    assert [r.placed for r in rows] == [{"sat"}, {"flank"}, {"sat"}]
    fig = Figure()
    axes = {k: fig.add_subplot(1, 3, i + 1) for i, k in enumerate(("hic", "te", "sat"))}
    prof = plot._profiles(f, rows, axes, seed=0)
    assert {k: (p["placed"], p["n"]) for k, p in prof.items()} == {"hic": (1, 3), "te": (0, 3), "sat": (2, 3)}
    assert {t.get_gid(): t.get_text() for t in _gid(fig, "placed:")} == {
        "placed:hic": "placed 1 of 3", "placed:te": "placed 0 of 3", "placed:sat": "placed 2 of 3"}
    d = next(p for p in _paragraphs(plot.legends(f)[0]) if p.startswith("**D**"))
    assert "by construction" in d and "placed by other signals" in d


def test_panel_d_uses_the_same_chromosomes_for_centromeres_and_random_places():
    """A chromosome shorter than panel D's span (here 1 Mb against a span of at least 5 Mb) is left out of both the
    centromeres' traces and the random places, so both cover the same chromosomes; n counts the rest."""
    f = _fig(n=3, c3=dict(chrom_len="1000000", start="400001", end="500000", seq_start="400001", seq_end="500000"))
    f["seq"].arrays.loc[2, ["start", "end"]] = [410_000, 490_000]
    rows = _rows(f)
    W = closeup.window(rows, plot._keys(f, rows), [r.L for r in rows])[0]
    assert W > 1_000_000
    fig = Figure()
    prof = plot._profiles(f, rows, {"sat": fig.add_subplot(1, 1, 1)}, seed=0, W=W)["sat"]
    assert prof["n"] == 2 and prof["obs"].shape[0] == 2 and prof["null"].shape[1] == 2
    assert np.isfinite(prof["null"]).all()
    assert plot._profile_keys(f, rows[2:], None, W) == []         # the short chromosome alone: no panel D


def test_short_legend_explains_a_dot_in_a_grayed_column():
    """A signal not trusted genome-wide (here the satellite arrays, whose test failed) can still be the evidence a
    centromere rests on: B shows a dot in its grayed column, and the legend says why; where the satellite is trusted,
    no such sentence."""
    f = _fig(n=3, summary=SUMMARY)
    b = next(p for p in _paragraphs(plot.legends(f)[0]) if p.startswith("**B**"))
    assert "A dot in a grayed column" in b and "lower confidence" in b
    b = next(p for p in _paragraphs(plot.legends(_fig(n=3, summary=FASTA))[0]) if p.startswith("**B**"))
    assert "grayed" not in b


def test_short_legend_says_what_the_hic_bar_counts():
    """Panel C's Hi-C row counts the chromosomes whose own test passed, not those Hi-C marks: the legend says so where
    that row has a bar, and not in a FASTA-only genome."""
    c = next(p for p in _paragraphs(plot.legends(_fig(n=3, summary=SUMMARY, hic=True))[0]) if p.startswith("**C**"))
    assert "For Hi-C contact with other centromeres, the bar counts the chromosomes whose own test passes." in c
    c = next(p for p in _paragraphs(plot.legends(_fig(n=3, summary=FASTA))[0]) if p.startswith("**C**"))
    assert "own test passes" not in c


def test_legends_say_a_far_signal_points_elsewhere_only_where_counted():
    """The roles paragraph: a signal points elsewhere where cenfuse counted it against the centromere or it is the
    satellite array or LTR cluster a moved or overruled centromere left behind; E: a colored bar well off the band is
    a disagreement, a gray bar a signal not used."""
    short, long = plot.legends(_drawn(_rich()))
    roles = next(p for p in _paragraphs(long) if p.startswith("**Roles.**"))
    assert "counted it against this centromere" in roles and "left behind" in roles
    assert "puts it over" not in roles
    e = next(p for p in _paragraphs(short) if p.startswith("**E**"))
    assert "well off the band" in e and "a gray bar a signal not used" in e
