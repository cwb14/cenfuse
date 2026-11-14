#!/usr/bin/env python3
"""Draw cenfuse's performance over every benchmark genome from the frozen results table (benchmark/results.tsv).

  python -m benchmark.performance [--png docs/img/performance.png] [--pdf docs/img/performance.pdf] [--numbers N.tsv]

A: one genome to scale, its mapped centromeres and its call with FASTA, +tes, +Hi-C and +tes+Hi-C. The genome is
chosen by rule: of the genomes with truth coordinates run with all four, the one whose gain from the FASTA alone to its
all-inputs run is nearest the median gain (ties: genome id). B: every chromosome with truth pooled per input set and
centromere class (one run per genome and input set, the 'agreement' mammals left out, as benchmark/summary.py (a)), each
genome's share right as a dot, beside the chance of a right call by luck (benchmark/README.md; rows whose truth is
withheld are left out of it). C: every genome and input set, controls included. D: confidence calibration over each
genome's all-inputs run. E: calls changed by adding real or randomized Hi-C. The PDF is the figure and its legend
page; --numbers writes every number it shows."""
import argparse
import re
import sys
from datetime import date
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt                                          # noqa: E402
import numpy as np                                                       # noqa: E402
import pandas as pd                                                      # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages                     # noqa: E402
from matplotlib.collections import PatchCollection                       # noqa: E402
from matplotlib.lines import Line2D                                      # noqa: E402
from matplotlib.offsetbox import AnnotationBbox, HPacker, TextArea       # noqa: E402
from matplotlib.patches import Patch, Rectangle                          # noqa: E402

from benchmark.summary import (CALL, HERE, RIGHT, SAME_ASSEMBLY, all_inputs, by_inputs, calibration,  # noqa: E402
                               controls, one_run, pooled_rows, read_table)

REPO = HERE.parent
OUTCOMES = ["right_high", "right_medium", "right_low", "conf_wrong", "wrong_low", "none"]   # bottom to top
COLOR = {"right_high": "#00694C", "right_medium": "#009E73", "right_low": "#9DD9C5", "conf_wrong": "#D55E00",
         "wrong_low": "#F6C3A2", "none": "#DADADA"}
CHANCE, INK, GRAY, CHROM, HIC = "#8A8A8A", "#222222", "#5A5A5A", "#EDEDED", "#0072B2"
CLASSES = {"satellite": "satellite arrays", "TE_or_none": "retrotransposons or none",
           "point_regional": "point or regional"}
TYPES = [("FASTA", "FASTA"), ("tes", "TEs (tes/)"), ("Hi-C", "Hi-C (hic/)"), ("genes", "genes (genes/)"),
         ("meth", "methylation (meth/)")]
WORD = {"tes": "TEs", "Hi-C": "Hi-C", "genes": "genes", "meth": "methylation"}
LADDER = ["FASTA", "FASTA+tes", "FASTA+Hi-C", "FASTA+tes+Hi-C"]
W, MM = 180.0, 1 / 25.4                                  # figure width (mm); inches per mm
FS, FR, FL, FP = 6.5, 6.0, 7.5, 10                       # base, row label, title and panel-letter sizes (pt)


# ---------------------------------------------------------------- data
def has(inputs, t):
    """True if the input set (e.g. 'FASTA+tes+meth(CpG)') holds data type t (FASTA, tes, Hi-C, genes, meth)."""
    return any(k == t or k.startswith(t + "(") for k in inputs.split("+"))


def column_key(inputs):
    """Input sets in order of how many inputs they hold, then TEs, Hi-C, genes, methylation first."""
    return (inputs.count("+"), [not has(inputs, t) for t, _ in TYPES[1:]], inputs)


def chance(d):
    """Per row, the chance of a right call by luck: (call + truth + 2 tol) / chromosome length, at most 1 (lengths of
    1-based inclusive intervals); 0 without a call; NaN where the truth is withheld ('') or absent with a call."""
    n = {c: pd.to_numeric(d[c], errors="coerce") for c in ("call_start", "call_end", "truth_start", "truth_end",
                                                             "tol_bp", "chrom_len")}
    v = ((n["call_end"] - n["call_start"] + 1) + (n["truth_end"] - n["truth_start"] + 1) + 2 * n["tol_bp"]) \
        / n["chrom_len"]
    v = v.clip(upper=1.0)
    v[d.call_start == "."] = 0.0
    v[d.truth_start == ""] = np.nan
    return v


def labels(canonical, truth):
    """genome -> (label, lineage): canonical.tsv, else the start of truth.tsv's source (before the accession)."""
    lab = {}
    if truth and Path(truth).exists():
        t = pd.read_csv(truth, sep="\t", dtype=str, keep_default_na=False)
        lab = {g: (re.split(r" GC[AF]_|:|;", s)[0].strip() or g, "") for g, s in zip(t.genome, t.source)}
    if canonical and Path(canonical).exists():
        c = pd.read_csv(canonical, sep="\t", dtype=str, keep_default_na=False)
        lab.update({g: (lb, li) for g, lb, li in zip(c.genome, c.label, c.lineage)})
    return lab


BINOMIAL = re.compile(r"\b([A-Z][a-z]+|[A-Z]\.) ([a-z]{3,})\b")


def runs_of(label, lineage, short=True):
    """A label as (text, style) runs: Latin binomials italic, the genus abbreviated when short; common names of
    mammals and birds left alone (heuristic: 'Sumatran orangutan' is no binomial)."""
    label = re.sub(r" var\. [a-z]+", "", label)
    if lineage.startswith(("mammals", "birds")):
        return [(label, "normal")]
    out, i = [], 0
    for m in BINOMIAL.finditer(label):
        genus = m.group(1)[0] + "." if short else m.group(1)
        out += [(label[i:m.start()], "normal"), (f"{genus} {m.group(2)}", "italic")]
        i = m.end()
    return [r for r in out + [(label[i:], "normal")] if r[0]]


def tables(d):
    """perf: the rows pooled in A, B and D (summary.pooled_rows: one run per genome and input set, one id per assembly);
    every genome's rows for C (both ids of an assembly), agreement and holocentric rows; alli: perf's all-inputs runs
    (summary.all_inputs: the run with the most input types)."""
    real = d[d.control == "none"]
    perf = pooled_rows(d)
    every = one_run(real[real.truth_kind != "agreement"])
    agree = one_run(real[real.truth_kind == "agreement"])
    holo = one_run(d[d.control == "holocentric"])
    return perf, every, agree, holo, all_inputs(perf)


def pooled(perf, d):
    """B: per input set and centromere class, genomes, chromosomes, each outcome and the chance (%); checked against
    summary.by_inputs."""
    rows = []
    for (cls, inp), x in perf.groupby(["cen_class", "inputs"]):
        r = dict(cen_class=cls, inputs=inp, genomes=x.genome.nunique(), chromosomes=len(x),
                 chance_rows=int(x.chance.notna().sum()), chance_pct=100 * x.chance.mean())
        r.update({o: int((x.outcome == o).sum()) for o in OUTCOMES})
        r["right"] = r["right_high"] + r["right_medium"] + r["right_low"]
        rows.append(r)
    p = pd.DataFrame(rows).set_index(["inputs", "cen_class"])
    s = by_inputs(d)
    s = s[s.index.get_level_values(1) != "all"]
    chk = p.loc[s.index].assign(wrong=lambda y: y.conf_wrong + y.wrong_low, no_call=lambda y: y.none)
    cols = ["genomes", "chromosomes", "right", "wrong", "no_call", "conf_wrong"]
    if not (chk[cols].values == s[cols].values).all():
        sys.exit("performance: pooled counts differ from benchmark.summary.by_inputs")
    return p


def hic_changes(d):
    """E: per run with Hi-C (real, control 'none', one per genome and input set; or randomized, every run, as
    summary.controls; one id per assembly), its chromosomes and those whose call (start, end or confidence) differs
    from the same genome's run with the same inputs but no Hi-C (the first arm). The randomized total is checked
    against summary.controls."""
    d = d[~d.genome.isin(list(SAME_ASSEMBLY))]                       # one Hi-C library per assembly
    base = one_run(d[~d.control.str.contains("random_hic")])
    runs = pd.concat([one_run(d[d.inputs.str.contains("Hi-C") & (d.control == "none")]),
                      d[d.control.str.contains("random_hic")]])
    rows = []
    for (arm, g), x in runs.groupby(["arm", "genome"]):
        want = re.sub(r"\+Hi-C(\(randomized\))?", "", x.inputs.iloc[0])
        p = base[(base.genome == g) & (base.inputs == want)]
        if len(p):
            m = x.merge(p, on="chrom", suffixes=("", "_p"))
            changed = (m[CALL].values != m[[c + "_p" for c in CALL]].values).any(axis=1)
            rows.append(dict(hic="randomized" if "random" in x.control.iloc[0] else "real", arm=arm, genome=g,
                             chromosomes=len(m), changed=int(changed.sum())))
    r = pd.DataFrame(rows, columns=["hic", "arm", "genome", "chromosomes", "changed"])
    c, _ = controls(d)
    y = r[r.hic == "randomized"]
    if len(c) and (int(pd.to_numeric(c.changed, errors="coerce").fillna(0).sum()) != y.changed.sum()
                   or int(c.chromosomes.sum()) != y.chromosomes.sum()):
        sys.exit("performance: randomized Hi-C runs differ from benchmark.summary.controls")
    return r


def pick_example(perf):
    """Of the genomes with truth coordinates run with all four LADDER input sets: the one whose gain (share right in
    its all-inputs run, the one with the most input types, minus with the FASTA alone) is nearest the median gain;
    ties by genome id."""
    x = perf[~perf.genome.isin(perf.genome[perf.truth_start == ""])]
    x = x.assign(r=x.outcome.isin(RIGHT))
    sets = x.groupby("genome").inputs.agg(set)
    alls = all_inputs(x).groupby("genome").r.mean()
    fasta = x[x.inputs == "FASTA"].groupby("genome").r.mean()
    cand = sorted(g for g, s in sets.items() if set(LADDER) <= s and g in alls.index)
    if not cand:
        return None, pd.Series(dtype=float), np.nan
    gain = pd.Series({g: alls[g] - fasta[g] for g in cand})
    med = float(gain.median())
    return min(cand, key=lambda g: (round(abs(gain[g] - med), 9), g)), gain, med      # exact ties: genome id


def check_example(ex, perf):
    """Structural checks of the example's intervals: inside the chromosome; truth width against its class median."""
    num = {c: pd.to_numeric(ex[c], errors="coerce") for c in ("call_start", "call_end", "truth_start", "truth_end",
                                                              "chrom_len")}
    called = ex.call_start != "."
    bad = int((called & ~((num["call_start"] >= 1) & (num["call_start"] <= num["call_end"])
                          & (num["call_end"] <= num["chrom_len"]))).sum())
    bad += int((~((num["truth_start"] >= 1) & (num["truth_start"] <= num["truth_end"])
                  & (num["truth_end"] <= num["chrom_len"]))).sum())
    t = perf[(perf.cen_class == ex.cen_class.iloc[0]) & ~perf.truth_start.isin(["", "."])]
    cls_w = (pd.to_numeric(t.truth_end) - pd.to_numeric(t.truth_start) + 1).median()
    ex_w = (num["truth_end"] - num["truth_start"] + 1).median()
    return bad, ex_w, cls_w


def chrom_names(chroms):
    """Chromosome numbers: 'Chr7' -> '7'; otherwise the rank in table order."""
    m = [re.fullmatch(r"(?:chr|Chr|CHR|chromosome)?_?0*(\d+|[XYZWU])", c) for c in chroms]
    return [x.group(1) for x in m] if all(m) else [str(i + 1) for i in range(len(chroms))]


def short_name(inputs):
    return "FASTA" if inputs == "FASTA" else "+" + "+".join(WORD[t] if t != "meth" else "meth"
                                                            for t, _ in TYPES[1:] if has(inputs, t))


def set_name(inputs):
    return "FASTA only" if inputs == "FASTA" else "+ " + " + ".join(WORD[t] for t, _ in TYPES[1:] if has(inputs, t))


# ---------------------------------------------------------------- drawing helpers
def rich(ax, x, y, runs, size, align=(1, 0.5), color=INK):
    """Text of mixed styles (runs of (text, style)) anchored at data (x, y)."""
    box = HPacker(children=[TextArea(t, textprops=dict(style=s, size=size, color=color)) for t, s in runs],
                  align="baseline", pad=0, sep=0)
    ax.add_artist(AnnotationBbox(box, (x, y), xycoords="data", box_alignment=align, frameon=False, pad=0,
                                 annotation_clip=False))


def text_mm(fig, runs, size):
    """Width in mm of runs of (text, style) at size pt."""
    r = fig.canvas.get_renderer()
    w = 0.0
    for t, s in runs:
        a = fig.text(0, 0, t.replace(" ", " "), size=size, style=s)
        w += a.get_window_extent(r).width / fig.dpi / MM
        a.remove()
    return w


def stacked(ax, x, width, counts, total, bottom=0.0, scale=100.0, horizontal=False, edge=None):
    """A bar split by outcome (counts: outcome -> n), drawn from bottom; each n / total of scale."""
    patches, colors = [], []
    for o in OUTCOMES:
        n = counts.get(o, 0)
        if n:
            h = scale * n / total
            patches.append(Rectangle((bottom, x - width / 2), h, width) if horizontal
                           else Rectangle((x - width / 2, bottom), width, h))
            colors.append(COLOR[o])
            bottom += h
    ax.add_collection(PatchCollection(patches, facecolors=colors, edgecolors=edge or "none", linewidths=0.15))


def swarm(ax, x0, ys, half, dx, **kw):
    """Dots at ys around x0: values within 3 units of each other side by side, spaced dx, within +-half."""
    ys = np.sort(np.asarray(ys, dtype=float))
    xs, b = np.zeros(len(ys)), np.round(ys / 3)
    for v in np.unique(b):
        i = np.where(b == v)[0]
        xs[i] = (np.arange(len(i)) - (len(i) - 1) / 2) * min(dx, 2 * half / max(len(i) - 1, 1))
    ax.scatter(x0 + xs, ys, zorder=3, clip_on=False, **kw)


def bare(ax):
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_xticks([])
    ax.set_yticks([])


def yaxis(ax, label=None, ticks=(0, 50, 100), labels_on=True):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.set_xticks([])
    ax.set_yticks(list(ticks))
    ax.tick_params(axis="y", labelsize=FS, length=2, width=0.5, pad=1, labelleft=labels_on)
    if label:
        ax.set_ylabel(label, size=FS, labelpad=1.5)


def glyphs(ax, cols, size=7.0, lw=0.7):
    """Which inputs each column holds: one row per data type, filled where given."""
    for j, inp in enumerate(cols):
        on = [i for i, (t, _) in enumerate(TYPES) if has(inp, t)]
        ax.plot([j, j], [0, max(on)], color=INK, lw=lw, zorder=1, solid_capstyle="round")
        ax.scatter([j] * len(TYPES), range(len(TYPES)), s=size, lw=0, zorder=2,
                   color=[INK if i in on else "#D3D3D3" for i in range(len(TYPES))])


# ---------------------------------------------------------------- the figure
def layout(blocks, lab):
    """Positions (mm): B and its glyphs on wide columns, C on its own columns after fitted label areas."""
    fig = plt.figure(figsize=(W * MM, 4))
    for b in blocks:
        ws = [text_mm(fig, runs_of(*lab.get(g, (g, ""))), FR) for k, g, _ in b["rows"] if k == "g"]
        b["lw"] = max(ws + [0]) + 4.8                         # + the chromosome count
    gw = max(text_mm(fig, [(t, "normal")], FR) for t in [n for _, n in TYPES] + ["chromosomes"]) + 2.2
    hh = max(text_mm(fig, [(short_name(c), "normal")], FS) for b in blocks for c in b["cols"]) + 1.6
    plt.close(fig)
    L = dict(left=1.0, right=W - 1.0, gw=gw, rh=2.15, gh=2.2, hh=hh)
    ncol = sum(len(b["cols"]) for b in blocks)
    gb, gc = 4.5, 2.5
    L["cwb"] = (L["right"] - L["left"] - gw - gb * (len(blocks) - 1)) / ncol
    L["cwc"] = (L["right"] - L["left"] - sum(b["lw"] for b in blocks) - gc * (len(blocks) - 1)) / ncol
    xb, xc = L["left"] + gw, L["left"]
    for b in blocks:
        k = len(b["cols"])
        b["xb"], b["wb"] = xb, k * L["cwb"]
        b["xc"], b["wc"] = xc + b["lw"], k * L["cwc"]
        xb += b["wb"] + gb
        xc += b["lw"] + b["wc"] + gc
    L["key"], L["a"], L["a_h"] = 1.0, 7.0, 31.0
    L["b"] = L["a"] + 6.5 + L["a_h"] + 5.0
    L["bar"], L["bar_h"] = L["b"] + 4.5, 29.0
    L["g"] = L["bar"] + L["bar_h"] + 1.2
    L["c"] = L["g"] + L["gh"] * (len(TYPES) + 2) + 3.0
    L["rows"] = L["c"] + L["hh"]
    n = [len(b["rows"]) for b in blocks]
    L["d"] = L["rows"] + max(n[1:] + [0]) * L["rh"] + 7.0
    L["d_h"] = max(L["rows"] + n[0] * L["rh"] - 4.0 - L["d"] - 4.5, 22.0)
    L["H"] = max(L["rows"] + n[0] * L["rh"], L["d"] + 4.5 + L["d_h"] + 4.0) + 1.0
    return L


def draw(d, lab, order, verbose=False):
    """The figure and the numbers behind it."""
    perf, every, agree, holo, alli = tables(d)
    pool = pooled(perf, d)
    framed = all_inputs(pd.concat([every, agree, holo]))                # each genome's run with the most inputs
    framed = set(zip(framed.genome, framed.arm))
    g_ex, gain, med = pick_example(perf)
    nums = []

    def note(panel, item, group, value, source="results.tsv"):
        nums.append(dict(panel=panel, item=item, group=group, value=value, source=source))

    blocks = []                                          # per class: its input sets (columns) and C's rows
    for cls in CLASSES:
        rows_t = sorted(every[every.cen_class == cls].genome.unique(), key=order)
        rows_a = sorted(agree[agree.cen_class == cls].genome.unique(), key=order)
        if rows_t or rows_a:
            rows = [("g", g, every) for g in rows_t]
            rows += [("h", "vs another caller", None)] + [("g", g, agree) for g in rows_a] if rows_a else []
            used = pd.concat([every[every.cen_class == cls], agree[agree.cen_class == cls]]).inputs.unique()
            blocks.append(dict(cls=cls, cols=sorted(used, key=column_key), rows=rows))
    if len(holo) and blocks:                             # the controls go under the last class
        b = blocks[-1]
        b["cols"] = sorted(set(b["cols"]) | set(holo.inputs), key=column_key)
        b["rows"] += [("h", "holocentric controls", None)] + [("g", g, holo)
                                                              for g in sorted(holo.genome.unique(), key=order)]
    L = layout(blocks, lab)
    H = L["H"]
    fig = plt.figure(figsize=(W * MM, H * MM))

    def ax_mm(x, y, w, h):
        return fig.add_axes([x / W, 1 - (y + h) / H, w / W, h / H])

    def letter(x, y, s):
        fig.text(x / W, 1 - y / H, s, size=FP, weight="bold", va="top", ha="left")

    # key
    hs = [Patch(fc=COLOR[o], ec="none") for o in OUTCOMES] + [
        Line2D([], [], ls="", marker="o", ms=2.6, mfc="white", mec=INK, mew=0.5), Patch(fc=CHANCE, ec="none"),
        Patch(fc=INK, ec="none"), Patch(fc="white", ec=INK, lw=0.6)]
    ls = ["right, high", "right, medium", "right, low", "wrong, confident", "wrong, low", "no call", "one genome",
          "chance", "mapped centromere", "all-inputs run"]
    fig.legend(hs, ls, loc="upper left", bbox_to_anchor=(L["left"] / W, 1 - L["key"] / H), ncol=len(ls),
               frameon=False, fontsize=FS, handlelength=1.0, handleheight=0.8, handletextpad=0.35,
               columnspacing=0.95, borderaxespad=0, borderpad=0)

    # ---- A: the example genome, to scale, with each of the four input sets
    letter(0.3, L["a"], "A")
    x0 = L["left"] + L["gw"]
    if g_ex is not None:
        ex = perf[(perf.genome == g_ex) & perf.inputs.isin(LADDER)]
        chroms = list(dict.fromkeys(ex[ex.inputs == "FASTA"].chrom))
        names = chrom_names(chroms)
        size = pd.to_numeric(ex.drop_duplicates("chrom").set_index("chrom").chrom_len).reindex(chroms) / 1e6
        n, aw = len(chroms), L["right"] - x0
        span = 4 * n + 3 * 1.2                            # n slots per input set, 1.2 between sets
        ax = ax_mm(x0, L["a"] + 6.5, aw, L["a_h"])
        hdr = ax_mm(x0, L["a"], aw, 6.0)
        for a in (ax, hdr):
            a.set_xlim(-0.2, span)
        top = float(size.max()) * 1.02
        ax.set_ylim(top, 0)
        hdr.set_ylim(0, 1)
        bare(hdr)
        mb = top / L["a_h"]                               # Mb per mm
        cw = min(0.3, 1.6 / (aw / (span + 0.2)))          # chromosome and call width (slots), at most 1.6 mm
        lb, li = lab.get(g_ex, (g_ex, ""))
        rich(hdr, -0.2, 1.0, runs_of(lb, li, short=False) + [(f" ({n} chromosomes)", "normal")], FL, align=(0, 1))
        bad, ex_w, cls_w = check_example(ex, perf)
        for k, inp in enumerate(LADDER):
            xs = k * (n + 1.2)
            hdr.text(xs + n / 2, 0.02, set_name(inp), ha="center", va="bottom", size=FL)
            e = ex[ex.inputs == inp].set_index("chrom")
            for i, c in enumerate(chroms):
                xc = xs + i + 0.5
                ax.add_patch(Rectangle((xc - cw - 0.04, 0), cw, size[c], fc=CHROM, ec="#9A9A9A", lw=0.4))
                r = e.loc[c]
                if r.truth_start not in ("", "."):
                    t0, t1 = int(r.truth_start) / 1e6, int(r.truth_end) / 1e6
                    h = max(t1 - t0, 0.9 * mb)
                    ax.add_patch(Rectangle((xc - cw - 0.1, (t0 + t1 - h) / 2), cw + 0.12, h, fc=INK, ec="none"))
                if r.call_start != ".":
                    c0, c1 = int(r.call_start) / 1e6, int(r.call_end) / 1e6
                    h = max(c1 - c0, 1.3 * mb)
                    ax.add_patch(Rectangle((xc + 0.04, (c0 + c1 - h) / 2), cw, h, fc=COLOR[r.outcome], ec="none"))
                ax.text(xc - 0.04, top + 1.0 * mb, names[i], ha="center", va="top", size=FR, color=GRAY)
            note("A", "right", inp, int(e.outcome.isin(RIGHT).sum()))
            note("A", "chromosomes", inp, len(e))
            note("A", "confidently wrong", inp, int((e.outcome == "conf_wrong").sum()))
        yaxis(ax, "Mb", ticks=[t for t in (0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000) if t <= top]
              if top > 150 else ax.get_yticks())
        ax.spines["bottom"].set_visible(False)
        ax.set_ylim(top, 0)
        note("A", "example genome", g_ex, lb, "rule: gain nearest the median of the candidates")
        note("A", "candidates", "genomes with truth run with FASTA, +tes, +Hi-C, +tes+Hi-C", len(gain))
        note("A", "median gain", "share right: all inputs minus FASTA", round(med, 4))
        note("A", "example gain", g_ex, round(float(gain[g_ex]), 4))
        note("A", "tolerance (bp)", g_ex, int(ex.tol_bp.iloc[0]), "results.tsv tol_bp")
        note("A", "intervals outside their chromosome", g_ex, bad, "check")
        note("A", "median truth width (bp)", g_ex, int(ex_w))
        note("A", "median truth width of its class (bp)", ex.cen_class.iloc[0], int(cls_w))
        if verbose:
            print(f"example: {g_ex} (gain {gain[g_ex]:.3f}, median {med:.3f} over {len(gain)} genomes); "
                  f"{bad} intervals outside the chromosome; truth width {ex_w:.0f} vs class median {cls_w:.0f}",
                  file=sys.stderr)

    # ---- B: pooled per input set, with every genome as a dot; the glyphs beneath
    letter(0.3, L["b"], "B")
    for bi, b in enumerate(blocks):
        cols, k, cls = b["cols"], len(b["cols"]), b["cls"]
        xb = ax_mm(b["xb"], L["bar"], b["wb"], L["bar_h"])
        xb.set_xlim(-0.5, k - 0.5)
        xb.set_ylim(0, 100)
        xb.text((k - 1) / 2, 103, CLASSES[cls], ha="center", va="bottom", size=FL)
        for j, inp in enumerate(cols):
            if (inp, cls) not in pool.index:
                continue
            p = pool.loc[(inp, cls)]
            stacked(xb, j - 0.1, 0.58, p.to_dict(), p.chromosomes)
            if p.chance_rows:
                xb.add_patch(Rectangle((j + 0.31 - 0.075, 0), 0.15, p.chance_pct, fc=CHANCE, ec="none"))
            x = perf[(perf.cen_class == cls) & (perf.inputs == inp)]
            s = x.assign(r=x.outcome.isin(RIGHT)).groupby("genome").r.mean()
            swarm(xb, j - 0.1, 100 * s.values, 0.22, 0.8 / L["cwb"], s=4.5, facecolor="white", edgecolor=INK, lw=0.4)
            for o in ["genomes", "chromosomes", "right", "right_high", "right_medium", "right_low", "conf_wrong",
                      "wrong_low", "none", "chance_rows"]:
                note("B", o, f"{cls} | {inp}", int(p[o]), "results.tsv (counts as summary.by_inputs)")
            note("B", "right (%)", f"{cls} | {inp}", round(100 * p.right / p.chromosomes, 1))
            note("B", "chance (%)", f"{cls} | {inp}", round(float(p.chance_pct), 1))
            note("B", "genome share right, min-median-max (%)", f"{cls} | {inp}",
                 "-".join(f"{100 * v:.0f}" for v in (s.min(), s.median(), s.max())))
        yaxis(xb, "chromosomes (%)" if bi == 0 else None, ticks=(0, 25, 50, 75, 100), labels_on=bi == 0)
        xb.spines["bottom"].set_linewidth(0.5)
        xg = ax_mm(b["xb"], L["g"], b["wb"], L["gh"] * (len(TYPES) + 2))
        bare(xg)
        xg.set_xlim(-0.5, k - 0.5)
        xg.set_ylim(len(TYPES) + 1.5, -0.5)
        glyphs(xg, cols)
        for j, inp in enumerate(cols):
            if (inp, cls) in pool.index:
                p = pool.loc[(inp, cls)]
                for i, v in ((len(TYPES), p.genomes), (len(TYPES) + 1, p.chromosomes)):
                    xg.text(j, i, f"{int(v)}", ha="center", va="center", size=FS - 0.5, color=GRAY)
        if bi == 0:
            for i, name in enumerate([nm for _, nm in TYPES] + ["genomes", "chromosomes"]):
                xg.text(-0.5 - 0.6 / L["cwb"], i, name, ha="right", va="center", size=FR,
                        color=INK if i < len(TYPES) else GRAY)

    # ---- C: every genome and input set
    letter(0.3, L["c"], "C")
    note("C", "genomes compared with another caller", "agreement", agree.genome.nunique())
    note("C", "holocentric genomes", "holocentric", holo.genome.nunique())
    note("C", "holocentric chromosome-runs", "holocentric", len(holo))
    note("C", "holocentric confident calls", "holocentric", int(holo.confidence.isin(["high", "medium"]).sum()))
    for b in blocks:
        cols, k = b["cols"], len(b["cols"])
        xc = ax_mm(b["xc"], L["rows"], b["wc"], len(b["rows"]) * L["rh"])
        bare(xc)
        xc.set_xlim(-0.5, k - 0.5)
        xc.set_ylim(len(b["rows"]) - 0.5, -0.5)
        for j, inp in enumerate(cols):                    # column heads: the inputs besides the FASTA
            xc.text(j, -0.5 - 0.6 / L["rh"], short_name(inp), rotation=90, ha="center", va="bottom", size=FS,
                    clip_on=False)
        nw = 3.9 / L["cwc"]                               # the chromosome count's column (data units)
        for i, (kind, g, src) in enumerate(b["rows"]):
            if kind == "h":
                xc.text(-0.5 - (b["lw"] - 0.3) / L["cwc"], i + 0.12, g, ha="left", va="center", size=FR,
                        style="italic", color=GRAY)
                continue
            x = src[src.genome == g]
            rich(xc, -0.5 - nw, i, runs_of(*lab.get(g, (g, ""))), FR)
            xc.text(-0.62, i, f"{x.chrom.nunique()}", ha="right", va="center", size=FR, color=GRAY)
            note("C", "chromosomes", g, x.chrom.nunique())
            for j, inp in enumerate(cols):
                y = x[x.inputs == inp]
                if not len(y):
                    continue
                stacked(xc, i, 0.66, y.outcome.value_counts().to_dict(), len(y), bottom=j - 0.42, scale=0.84,
                        horizontal=True, edge="white" if 0.84 * L["cwc"] / len(y) > 0.3 else None)
                if (g, y.arm.iloc[0]) in framed:
                    xc.add_patch(Rectangle((j - 0.47, i - 0.45), 0.94, 0.9, fill=False, ec=INK, lw=0.6))
                note("C", "right/chromosomes", f"{g} | {inp}", f"{int(y.outcome.isin(RIGHT).sum())}/{len(y)}")

    # ---- D: calibration over each genome's all-inputs run
    bd = blocks[1] if len(blocks) > 1 else blocks[0]
    dx = bd["xc"] if len(blocks) > 1 else bd["xc"] + bd["wc"] + 12
    letter(dx - 11.0, L["d"], "D")
    cal = calibration(alli)
    xd = ax_mm(dx, L["d"] + 4.5, 34.0, L["d_h"])
    xd.set_xlim(-0.6, 2.6)
    xd.set_ylim(0, 100)
    xd.text(1, 104, "calls, all inputs", ha="center", va="bottom", size=FL)
    for j, c in enumerate(["high", "medium", "low"]):
        y = alli[alli.confidence == c]
        if not len(y):
            continue
        stacked(xd, j - 0.1, 0.58, y.outcome.value_counts().to_dict(), len(y))
        ch = 100 * y.chance.mean()
        if y.chance.notna().any():
            xd.add_patch(Rectangle((j + 0.31 - 0.075, 0), 0.15, ch, fc=CHANCE, ec="none"))
        s = y.assign(r=y.outcome.isin(RIGHT)).groupby("genome").r.mean()
        swarm(xd, j - 0.1, 100 * s.values, 0.27, 0.8 * 3.2 / 34.0, s=4.5, facecolor="white", edgecolor=INK, lw=0.4)
        xd.text(j - 0.1, 2, f"{len(y)}", ha="center", va="bottom", size=FS - 0.5,
                color="white" if c != "low" else INK)
        note("D", "calls", c, int(cal.loc[c, "calls"]), "results.tsv via summary.calibration")
        note("D", "right", c, int(cal.loc[c, "right"]), "results.tsv via summary.calibration")
        note("D", "right (%)", c, round(100 * cal.loc[c, "right"] / cal.loc[c, "calls"], 1))
        note("D", "chance (%)", c, round(float(ch), 1))
    note("D", "no call", "all inputs", int(cal.loc["no call", "calls"]), "results.tsv via summary.calibration")
    note("D", "genomes", "all inputs", alli.genome.nunique())
    yaxis(xd, "calls (%)")
    xd.set_xticks([-0.1, 0.9, 1.9])
    xd.set_xticklabels(["high", "medium", "low"], size=FS)
    xd.tick_params(axis="x", length=0, pad=1.5)
    xd.spines["bottom"].set_linewidth(0.5)

    # ---- E: calls changed by adding Hi-C, real against randomized
    hic = hic_changes(d)
    ex_ = dx + 34.0 + 17.0
    letter(ex_ - 11.0, L["d"], "E")
    we = min(24.0, L["right"] - ex_)
    xe = ax_mm(ex_, L["d"] + 4.5, we, L["d_h"])
    xe.set_xlim(-0.6, 1.6)
    xe.set_ylim(0, 100)
    xe.text(0.5, 104, "adding Hi-C", ha="center", va="bottom", size=FL)
    for j, kind in enumerate(["real", "randomized"]):
        y = hic[hic.hic == kind]
        if not len(y):
            continue
        share = 100 * y.changed.sum() / y.chromosomes.sum()
        xe.add_patch(Rectangle((j - 0.29, 0), 0.58, share, fc=HIC, ec="none"))
        per = y.groupby("genome")[["changed", "chromosomes"]].sum()          # one dot per genome, its runs pooled
        swarm(xe, j, 100 * per.changed / per.chromosomes, 0.27, 0.8 * 2.2 / we, s=4.5, facecolor="white",
              edgecolor=INK, lw=0.4)
        note("E", "runs", kind, len(y))
        note("E", "genomes", kind, len(per))
        note("E", "chromosome-runs", kind, int(y.chromosomes.sum()))
        note("E", "calls changed", kind, int(y.changed.sum()))
        note("E", "calls changed (%)", kind, round(float(share), 1))
    yaxis(xe, "calls changed (%)")
    xe.set_xticks([0, 1])
    xe.set_xticklabels(["real", "randomized"], size=FS)
    xe.tick_params(axis="x", length=0, pad=1.5)
    xe.spines["bottom"].set_linewidth(0.5)
    return fig, nums, dict(perf=perf, every=every, agree=agree, holo=holo, alli=alli, pool=pool, example=g_ex,
                           gain=gain, med=med, cal=cal, hic=hic)


# ---------------------------------------------------------------- legend page
def legend_runs(d, lab, r):
    """The legend as paragraphs of (text, style) runs; every number from the table."""
    pool, cal, alli, perf, holo, agree = r["pool"], r["cal"], r["alli"], r["perf"], r["holo"], r["agree"]

    def pct(cls, inp):
        p = pool.loc[(inp, cls)] if (inp, cls) in pool.index else None
        return "." if p is None else f"{100 * p.right / p.chromosomes:.0f}%"

    def pct_of(y):
        return f"{100 * y.changed.sum() / max(y.chromosomes.sum(), 1):.0f}% of {int(y.chromosomes.sum()):,}"

    def size(bp):
        return f"{bp / 1e6:g} Mb" if bp >= 1e6 else f"{bp / 1e3:g} kb"

    def span(t):
        t = sorted(set(pd.to_numeric(t).astype(int)))
        if len(t) < 2:
            return size(t[0]) if t else "."
        u, name = (1e6, "Mb") if t[0] >= 1e6 else (1e3, "kb")
        return f"{t[0] / u:g}-{t[-1] / u:g} {name}"

    g = r["example"]
    ex = perf[perf.genome == g]
    tol_s = size(int(ex.tol_bp.iloc[0])) if len(ex) else "."
    pr = perf[perf.cen_class == "point_regional"]
    hic_pr = pr.inputs.str.contains("Hi-C")
    twice = []                                       # an assembly under two genome ids (one Hi-C library each)
    for drop, keep in SAME_ASSEMBLY.items():
        lb = lab.get(keep, (keep, ""))
        lib = re.search(r"\((.+)\)", lb[0])
        if drop in set(r["every"].genome):
            twice += [(runs_of(re.sub(r" \(.*\)", "", lb[0]), lb[1])[0][0], "italic"),
                      (" appears twice, one assembly with two Hi-C libraries; B, D and E pool only the one with more "
                       f"read pairs{' (' + lib.group(1) + ')' if lib else ''}. ", "normal")]
    hic = r["hic"]
    real, rnd = hic[hic.hic == "real"], hic[hic.hic == "randomized"]
    hi, me, lo = (100 * cal.loc[c, "right"] / max(cal.loc[c, "calls"], 1) for c in ("high", "medium", "low"))
    N = "normal"
    paras = [
        [("What a centromere is made of decides which data find it.", "bold"),
         (f" From the FASTA alone cenfuse was right on {pct('satellite', 'FASTA')} of satellite centromeres. "
          f"Centromeres of retrotransposons or of no repeat went from {pct('TE_or_none', 'FASTA')} to "
          f"{pct('TE_or_none', 'FASTA+tes')} with the TE annotation, and the point and regional centromeres of fungi "
          f"and ", N), ("Plasmodium", "italic"),
         (f" from {pct('point_regional', 'FASTA')} to {pct('point_regional', 'FASTA+Hi-C')} with Hi-C.", N)],
        [("A", "bold"), (" One genome, ", N)] + (runs_of(*lab.get(g, (g, "")), short=False) if g else []) +
        [(f": its chromosomes to scale, its centromeres mapped by {ex.truth_kind.iloc[0] if len(ex) else ''} "
          "(black) and the call made from each set of inputs, colored by outcome (small ones drawn as a visible "
          "mark). Its gain from the FASTA alone to all its inputs is nearest the median of the "
          f"{len(r['gain'])} genomes run with all four sets. A call is right within {tol_s} of the mapped "
          "centromere, so one beside it can still miss.", N)],
        [("B", "bold"), (" Every chromosome with a mapped centromere, pooled by inputs and centromere type: right "
                         "(green, darker when more confident), wrong (vermillion, strong when confident) or no call "
                         "(gray). Dots: each genome's share right. Thin gray bars: the share a call of the same size "
                         "placed at random would get right. Beneath: the inputs of each column, and the genomes and "
                         "chromosomes it pools. In fungi and ", N), ("Plasmodium", "italic"),
         (f" a call is right within {span(pr[~hic_pr].tol_bp)} without Hi-C and within one Hi-C bin "
          f"({span(pr[hic_pr].tol_bp)}) with it.", N)],
        [("C", "bold"), (" Every genome with every set of inputs it was run with; each cell splits its chromosomes "
                         "by outcome, and a frame marks the genome's run with the most inputs, the one D pools. ", N)]
        + twice + [(f"{agree.genome.nunique()} mammals were compared with another tool's satellite locus, not "
                         f"with a mapped centromere, and are left out of B. The {holo.genome.nunique()} holocentric "
                         "genomes have no single centromere, so every call counts as wrong; of their "
                         f"{len(holo)} chromosomes over all sets of inputs, "
                         f"{int(holo.confidence.isin(['high', 'medium']).sum())} got a confident call.", N)],
        [("D", "bold"), (f" Over the {alli.genome.nunique()} genomes with a mapped centromere, each run with its "
                         f"most inputs, {hi:.0f}% of high, {me:.0f}% of medium and {lo:.0f}% of low calls were right "
                         f"(calls per bar at its foot; dots: genomes); {int(cal.loc['no call', 'calls'])} "
                         "chromosomes got no call.", N)],
        [("E", "bold"), (" Calls that adding Hi-C changed (moved, or their confidence), against the same genome "
                         f"without Hi-C: {pct_of(real)} with real Hi-C ({len(real)} runs), {pct_of(rnd)} with Hi-C "
                         f"whose contacts were randomized ({len(rnd)} runs). Dots: genomes."
                         + (" Hi-C without centromere signal changes nothing." if len(rnd) and
                            rnd.changed.sum() == 0 else "") + " Whether real Hi-C helped: B.", N)],
        [(f"Data: benchmark/results.tsv (cenfuse {', '.join(sorted(d.version.unique()))}) · Script: "
          f"benchmark/performance.py · {date.today().isoformat()}", N)],
    ]
    return paras


def legend_page(paras, size=10.0):
    """One page as wide as the figure, words laid out ragged right at size pt, its height fitted to the text."""
    fig = plt.figure(figsize=(W * MM, 4))
    r = fig.canvas.get_renderer()
    cache = {}

    def width(word, style):
        if (word, style) not in cache:
            t = fig.text(0, 0, word, size=size, style="italic" if style == "italic" else "normal",
                         weight="bold" if style == "bold" else "normal")
            cache[(word, style)] = t.get_window_extent(r).width / fig.dpi
            t.remove()
        return cache[(word, style)]

    space = width("a a", "normal") - width("aa", "normal")
    left, right, line_h = 6 * MM, (W - 6) * MM, size / 72 * 1.32
    lines = []
    for p in paras:
        toks, pending = [], False                       # (word, style, space before)
        for text, style in p:
            i = 0
            for m in re.finditer(r"\S+", text):
                toks.append((m.group(), style, bool(toks) and (pending or m.start() > i)))
                pending, i = False, m.end()
            pending = pending or i < len(text)
        line, x = [], left
        for w, style, sp in toks:
            gap = space if (sp and line) else 0.0
            if x + gap + width(w, style) > right and line:
                lines.append(line)
                line, x, gap = [], left, 0.0
            line.append((x + gap, w, style))
            x += gap + width(w, style)
        lines += [line, None]                           # None: the gap between paragraphs
    plt.close(fig)
    h = (sum(1.0 if ln else 0.45 for ln in lines) + 1.0) * line_h + 8 * MM
    fig = plt.figure(figsize=(W * MM, h))
    y = h - 6 * MM - line_h * 0.8
    for ln in lines:
        if ln is None:
            y -= 0.45 * line_h
            continue
        for x, w, style in ln:
            fig.text(x / (W * MM), y / h, w, size=size, va="baseline", ha="left",
                     style="italic" if style == "italic" else "normal", weight="bold" if style == "bold" else "normal")
        y -= line_h
    return fig


# ---------------------------------------------------------------- main
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--table", type=Path, default=HERE / "results.tsv",
                    help="the results table [benchmark/results.tsv]")
    ap.add_argument("--canonical", type=Path, default=HERE / "canonical.tsv",
                    help="genome labels, lineages and row order [benchmark/canonical.tsv]")
    ap.add_argument("--truth", type=Path, default=HERE / "truth.tsv",
                    help="the truth manifest: labels of the genomes canonical.tsv lacks [benchmark/truth.tsv]")
    ap.add_argument("--png", type=Path, default=REPO / "docs/img/performance.png",
                    help="the figure, without legend [docs/img/performance.png]")
    ap.add_argument("--pdf", type=Path, default=REPO / "docs/img/performance.pdf",
                    help="the figure and its legend page [docs/img/performance.pdf]")
    ap.add_argument("--numbers", type=Path, help="write every number on the figure and in the legend here (TSV)")
    ap.add_argument("--dpi", type=int, default=300, help="PNG resolution [300]")
    ap.add_argument("-v", "--verbose", action="store_true", help="print the example's selection and the pooled table")
    a = ap.parse_args(argv)
    plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
                         "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none", "font.size": FS,
                         "axes.linewidth": 0.5})
    if not a.table.exists():
        sys.exit(f"performance: no table {a.table}")
    d = read_table(a.table)
    d["chance"] = chance(d)
    lab = labels(a.canonical, a.truth)
    canon = list(pd.read_csv(a.canonical, sep="\t", dtype=str).genome) if a.canonical.exists() else []
    rank = {g: i for i, g in enumerate(canon)}

    def order(g):
        return (rank.get(g, len(rank)), lab.get(g, (g, ""))[0].lower(), g)

    print(f"performance: {a.table.name}: {len(d)} chromosome-runs, {d.genome.nunique()} genomes", file=sys.stderr)
    fig, nums, r = draw(d, lab, order, a.verbose)
    if a.verbose:
        print(r["pool"].drop(columns=["chance_rows"]).round(1).to_string(), file=sys.stderr)
    for p in (a.png, a.pdf):
        p.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(a.png, dpi=a.dpi)
    leg = legend_page(legend_runs(d, lab, r))
    with PdfPages(a.pdf) as pdf:
        pdf.savefig(fig)
        pdf.savefig(leg)
    plt.close("all")
    if a.numbers:
        a.numbers.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(nums).to_csv(a.numbers, sep="\t", index=False)
    print(f"performance: wrote {a.png}, {a.pdf}" + (f", {a.numbers}" if a.numbers else ""), file=sys.stderr)


if __name__ == "__main__":
    main()
