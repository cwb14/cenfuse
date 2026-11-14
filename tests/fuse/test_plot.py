"""Overview plot: panel C (each line's call relative to the sequence call) scales to the calls drawn; the pages of a
genome are drawn, saved and closed one at a time."""
import re
import warnings

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pytest
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.figure import Figure

from cenfuse import plot


def test_offsets_axis_covers_sequence_only_calls():
    """FASTA-only run with every chromosome called: only rectangle patches are drawn, which do not rescale the axis
    by themselves (mouse and maize kept 0-1 Mb, showing the upper half of each call)."""
    rows = [dict(chrom=f"c{i}", seq=(1_000_000, 9_000_000), call=(1_000_000, 9_000_000), te=None, hic=None)
            for i in range(3)]
    fig, ax = plt.subplots()
    plot._offsets(ax, rows)
    fig.canvas.draw()
    lo, hi = ax.get_ylim()
    plt.close(fig)
    assert lo <= -4 and hi >= 4


def _fig(rdna, n=1):
    """The least `overview` needs, for n chromosomes (c1..cn), each with one satellite array and a sequence call."""
    import pandas as pd
    from types import SimpleNamespace
    names = [f"c{i}" for i in range(1, n + 1)]
    calls = pd.DataFrame([dict(chrom=c, chrom_len=10_000_000, start="4000001", end="5000000", seq_start="4000001",
                               seq_end="5000000", relation="sequence_only", confidence="high", call_source="sequence",
                               support="satellite", against=".", te_start=".", te_end=".", hic_anch_core_start=".",
                               hic_anch_core_end=".", hic_anch_peak_start=".", hic_anch_peak_end=".") for c in names])
    arrays = pd.DataFrame(dict(chrom=names, start=[4_100_000] * n, end=[4_900_000] * n, role=["centromeric"] * n))
    return dict(genome=SimpleNamespace(name="g", lines=["satellite"]), seq=SimpleNamespace(arrays=arrays), tel={},
                hic=None, meth=None, calls=calls, telo={}, plants=None, out="x", rdna=rdna)


def test_overview_marks_45s_and_5s_rdna_arrays():
    """45S arrays (NORs) as filled and 5S arrays as open gray diamonds at their midpoints, in the key and the legend;
    a genome without either says so."""
    from matplotlib.colors import to_hex
    f = _fig(dict(nors=[("c1", 100_000, 400_000, 30, 9_100)], five_s=[("c1", 8_000_000, 8_200_000, 80, 500)]))
    fig = plot._pages(f)[0]
    marks = sorted((float(x.get_xdata()[0]), to_hex(x.get_markerfacecolor())) for x in fig.axes[0].lines
                   if x.get_marker() == "D")
    assert marks == [(250_000.0, to_hex(plot.RDNA_MARK)), (8_100_000.0, "#ffffff")]
    labels = [t.get_text() for t in fig.legends[0].get_texts()]
    assert "45S rDNA array (NOR)" in labels and "5S rDNA array" in labels
    plt.close(fig)
    assert "1 of 45S" in plot.legends(f)[1] and "1 of 5S" in plot.legends(f)[1]
    f = _fig(dict(nors=[], five_s=[]))
    fig = plot._pages(f)[0]
    assert not [x for x in fig.axes[0].lines if x.get_marker() == "D"]
    assert not {"45S rDNA array (NOR)", "5S rDNA array"} & {t.get_text() for t in fig.legends[0].get_texts()}
    plt.close(fig)
    assert "No 45S or 5S rDNA array was found" in plot.legends(f)[1]


def _composition(f, i, flag, line):
    """Chromosome i's call placed by a FASTA line alone (call_source composition): an AT trough, a repeat-density peak or a
    repeat-family cluster, with no sequence call."""
    for k, v in dict(call_source="composition", seq_start=".", seq_end=".", confidence="low", flags=flag, support=line).items():
        f["calls"].loc[i, k] = v


def test_repeat_line_calls_have_their_own_band_key_and_columns():
    """0.16's repeat-density and repeat-family calls are drawn in the repeat lines' color and named so in the key, not as AT
    troughs, and panel B shows their lines (rep, fam) with the supporting dot on each call's row. C. neoformans from the FASTA
    (0.16.0): 14 of 14 calls rep_only, drawn as 'integrated call (AT trough)' with an empty panel B."""
    from matplotlib.colors import to_hex
    f = _fig(dict(nors=[], five_s=[]), n=3)
    _composition(f, 0, "rep_only", "repeats")
    _composition(f, 1, "fam_only", "families")
    fig = plot._pages(f)[0]
    keys = [t.get_text() for t in fig.legends[0].get_texts()]
    assert "integrated call (repeat peak)" in keys and "integrated call (AT trough)" not in keys
    bands = {round(p.get_y()): to_hex(p.get_facecolor()) for p in fig.axes[0].patches if p.get_zorder() == 0}
    band = {k: to_hex(v) for k, v in plot.BAND.items()}
    assert bands == {2: band["repeats"], 1: band["repeats"], 0: band["sequence"]}                    # rows top-down
    assert band["repeats"] not in (band["composition"], band["hic"], band["sequence"])
    cols = [t.get_text() for t in fig.axes[1].get_xticklabels()]
    assert cols == ["sat", "TE", "Hi-C", "rep", "fam", "meth"]
    dots = {(round(float(x.get_xdata()[0])), float(x.get_ydata()[0])) for x in fig.axes[1].lines if x.get_marker() == "o"}
    assert dots == {(cols.index("rep"), 2.5), (cols.index("fam"), 1.5), (cols.index("sat"), 0.5)}
    guides = sorted(float(x.get_xdata()[0]) for x in fig.axes[1].lines if to_hex(x.get_color()) == "#ebebeb")
    assert guides == list(range(len(cols)))                     # a guide under every column, not only the first four
    plt.close(fig)
    _composition(f, 2, "at_only", "at")                         # with an AT call too: both keys, each band its own color
    fig = plot._pages(f)[0]
    keys = [t.get_text() for t in fig.legends[0].get_texts()]
    assert "integrated call (repeat peak)" in keys and "integrated call (AT trough)" in keys
    bands = {round(p.get_y()): to_hex(p.get_facecolor()) for p in fig.axes[0].patches if p.get_zorder() == 0}
    assert bands[0] == band["composition"] and bands[2] == band["repeats"]
    plt.close(fig)


def test_the_long_legend_names_every_rule_that_places_a_call():
    """The 'Terms' paragraph names the rules that can place a call besides the Hi-C gate and the Rabl axis (stale since 0.12):
    the LTR cores, the cis signal, the AT troughs and (0.16) the two repeat lines; panel A's legend says what each band color
    means, and panel B's what its AT and repeat columns mean."""
    f = _fig(dict(nors=[], five_s=[]), n=2)
    _composition(f, 0, "rep_only", "repeats")
    long = plot.legends(f)[1]
    terms = long.split("\n\n")[0]
    for x in ("LTR core", "cis signal", "AT trough", "densest stretch of repeated sequence", "repeat-family cluster"):
        assert x in terms, x
    a = long.split("\n\n")[1]
    assert "wine: a repeat peak" in a and "green" not in a       # only the bands drawn are named
    assert "AT, rep and fam" in long.split("\n\n")[2]


def test_relation_labels_end_inside_the_page():
    """Panel B's longest label ("sequence only, medium") ends inside the page: the PNG has no margin to spare (10 of
    the 23 v8 gallery pages cut it at the right edge)."""
    f = _fig(dict(nors=[], five_s=[]))
    f["calls"].loc[0, "confidence"] = "medium"
    fig = plot._pages(f)[0]
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    labels = [t for t in fig.axes[1].texts if t.get_text() == "sequence only, medium"]
    assert labels and all(t.get_window_extent(r).x1 <= fig.bbox.x1 for t in labels)
    a, b = fig.axes[0].get_window_extent(r), fig.axes[1].get_window_extent(r)
    assert a.x1 < b.x0                                       # panel A stays left of panel B
    plt.close(fig)


def test_long_chromosome_names_stay_on_the_page():
    """Panel A's names and panel C's rotated names are drawn whole, however long (a GWH accession has 17 characters; the left margin
    held about 15). The two left edges stay aligned with each other."""
    f = _fig(dict(nors=[], five_s=[]), n=4)
    long = {f"c{i}": f"GWHERKE0000000{i}.1" for i in range(1, 5)}
    f["calls"]["chrom"] = f["calls"]["chrom"].map(long)
    f["seq"].arrays["chrom"] = f["seq"].arrays["chrom"].map(long)
    pages = plot._pages(f)
    fig = pages[-1]
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    width, height = fig.get_size_inches() * fig.dpi
    seen = 0
    for ax in fig.axes:
        for t in ax.texts + ax.get_xticklabels():
            if t.get_text().startswith("GWH"):
                box = t.get_window_extent(renderer)
                assert box.x0 >= 0 and box.y0 >= 0 and box.x1 <= width and box.y1 <= height, t.get_text()
                seen += 1
    assert seen == 8                                       # every name in panel A and in panel C
    assert fig.axes[0].get_position().x0 == pytest.approx(fig.axes[2].get_position().x0)   # panels A and C share their left edge
    plt.close("all")


@pytest.fixture
def quick_pages(monkeypatch):
    """Pages of 2 rows and a fixed label width. The page lifecycle depends on neither, and measuring each chromosome's
    relation label (a TextPath, 0.1 s each) would take over a minute for the 600 rows of 25 full pages."""
    monkeypatch.setattr(plot, "ROWS_PER_PAGE", 2)
    monkeypatch.setattr(plot, "_text_mm", lambda s, size: 20.0)


def test_overview_keeps_one_page_open_at_a_time(tmp_path, monkeypatch, quick_pages):
    """A genome of more pages than matplotlib's open-figure limit (20; the 670-sequence Cissus run had 28) is drawn,
    saved and closed one page at a time: no "More than 20 figures" RuntimeWarning, exactly one figure open at every
    save (a page is closed before the next is drawn), none left, and every page in the PDF."""
    n_pages, n_legend_pages = 25, 2                           # the figure pages; the short and the long legend
    f = _fig(dict(nors=[], five_s=[]), n=n_pages * plot.ROWS_PER_PAGE)
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
    f = _fig(dict(nors=[], five_s=[]), n=2 * plot.ROWS_PER_PAGE + 1)        # three pages
    n_calls = []

    def fail_second(real):
        def wrapped(*args, **kwargs):
            n_calls.append(1)
            if len(n_calls) == 2:
                raise OSError("disk full")
            return real(*args, **kwargs)
        return wrapped

    if where == "draw":
        monkeypatch.setattr(plot, "_gallery", fail_second(plot._gallery))
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


def _hic(size, anchored):
    """A Hi-C result with two 60-bin chromosomes: anchored score all NaN unless `anchored` (no sequence call to anchor
    partner windows), prior-free score peaking at bin 30."""
    import numpy as np
    from cenfuse.hic.bins import Bins
    b = Bins({"c1": 60 * size, "c2": 60 * size}, size)
    peak = np.ones(b.N)
    peak[b.off + 30] = 5.0
    return dict(bins=b, members=np.array([True, True]), S=peak.copy() if anchored else np.full(b.N, np.nan),
                denovo_S=peak, bin_size=size)


@pytest.mark.filterwarnings("ignore:All-NaN slice")     # grid columns no trace reaches (as before 0.10)
def test_score_panel_scales_with_the_bin_and_falls_back_to_the_prior_free_score():
    rows = [dict(chrom=c, call=(30 * 10_000, 31 * 10_000), relation="functional_only") for c in ("c1", "c2")]
    fig = Figure()
    ax = fig.add_subplot()
    plot._overlay(ax, rows, _hic(10_000, anchored=False))
    assert ax.get_xlim() == (-0.4, 0.4)                       # +/- 40 bins of 10 kb
    assert ax.get_ylabel().startswith("prior-free Hi-C score")
    fig = Figure()
    ax = fig.add_subplot()
    rows = [dict(chrom=c, call=(30 * 250_000, 31 * 250_000), relation="sequence_only") for c in ("c1", "c2")]
    plot._overlay(ax, rows, _hic(250_000, anchored=True))
    assert ax.get_xlim() == (-10.0, 10.0) and ax.get_ylabel().startswith("anchored Hi-C score")   # as before 0.10
