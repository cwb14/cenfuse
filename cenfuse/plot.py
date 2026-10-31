"""Overview figure, one per genome (cenfuse_overview.pdf: figure pages, a short and a long legend page per genome;
<genome>.overview.png: the first figure page at 300 dpi).

A, every chromosome to one bp scale, one row each:
- the anchored Hi-C score (blue) and the prior-free Hi-C score (gray);
- the methylation track (pooled composition-adjusted CHG);
- every satellite array, colored by its family's role;
- the TE marker clade's density;
- the sequence call (black outline) and the integrated call (band).
B, aligned to A's rows: which lines support or oppose the integrated call.
C: each line's call relative to the sequence call.
D: the anchored Hi-C score around every integrated call, one trace per chromosome.
With scored plants (docs/plants.md), a second page:
E: one row per plant, one column per chromosome: does the plant's centromere sit at the reference call?
F: the chromosome where the plants differ most, plant by plant: each plant's methylation track under the reference
   plant's evidence, at one bp scale."""
from __future__ import annotations

import logging
import textwrap
from datetime import date
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages  # noqa: E402
from matplotlib.collections import PolyCollection  # noqa: E402
from matplotlib.font_manager import FontProperties  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch, Rectangle  # noqa: E402
from matplotlib.textpath import TextPath  # noqa: E402

log = logging.getLogger("cenfuse")
MM = 1 / 25.4
WIDTH = 180 * MM
B_COL_MM = 4.6                                   # panel B: spacing of the line columns
ROWS_PER_PAGE = 24
ROLE_COLOR = {"centromeric": "#D55E00", "cen_associated": "#F4A582", "interstitial": "#CC79A7",
              "subtelomeric": "#56B4E9", "rDNA": "#009E73"}
OTHER = "#B0B0B0"
HIC, DENOVO, METH, TE = "#0072B2", "#A0A0A0", "#303030", "#E69F00"
BAND = {"hic": "#C6DBEF", "sequence": "#E8E8E8"}
LINES = ["satellite", "te", "hic", "meth"]
LINE_LABEL = {"satellite": "sat", "te": "TE", "hic": "Hi-C", "fold": "fold", "meth": "meth"}
REL_SHORT = {"concordant": "concordant", "offset": "moved", "contested": "contested",
             "functional_only": "functional only", "sequence_only": "sequence only", "none": "none"}
VERDICT_COLOR = {"shared": "#BDBDBD", "shifted": "#D55E00", "possible_shift": "#E69F00", "call_absent": "#000000",
                 "conflict": "#CC79A7", "unclear": "#F2F2F2", "no_data": "white", "no_call": "white"}
VERDICT_LABEL = {"shared": "consistent with the call", "shifted": "shifted", "possible_shift": "possible shift",
                 "call_absent": "call region absent", "conflict": "lines disagree", "unclear": "unclear",
                 "no_data": "no data"}
PROV_SHORT = {"reference": "reference plant", "other": "other plant", "unknown": "origin unknown"}
plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Liberation Sans", "Arial", "DejaVu Sans"],
                     "font.size": 6.5, "pdf.fonttype": 42, "ps.fonttype": 42, "axes.linewidth": 0.5,
                     "xtick.major.width": 0.5, "ytick.major.width": 0.5, "xtick.major.size": 2,
                     "ytick.major.size": 2})


def _iv(a, b):
    return (int(a) - 1, int(b)) if a != "." and b != "." else None


def _rows(calls):
    out = []
    for r in calls.to_dict("records"):
        out.append(dict(chrom=r["chrom"], L=int(r["chrom_len"]), call=_iv(r["start"], r["end"]),
                        seq=_iv(r["seq_start"], r["seq_end"]), relation=r["relation"], conf=r["confidence"],
                        source=r["call_source"], support=set(str(r["support"]).split(",")) - {"."},
                        against=set(str(r["against"]).split(",")) - {"."},
                        te=_iv(r["te_start"], r["te_end"]),
                        hic=_iv(r["hic_anch_core_start"], r["hic_anch_core_end"]),
                        hic_peak=_iv(r["hic_anch_peak_start"], r["hic_anch_peak_end"])))
    return out


def _track(hic, chrom, key):
    """(x bp centers, y over the chromosome median) of a Hi-C score track, or None."""
    if hic is None:
        return None
    b = hic["bins"]
    if chrom not in b.names:
        return None
    k = b.names.index(chrom)
    if not hic["members"][k]:
        return None
    y = np.asarray(hic[key][b.sl(k)], float)
    med = np.nanmedian(y) if np.isfinite(y).any() else np.nan
    if not med > 0:
        return None
    x = (b.start[b.sl(k)] + b.end[b.sl(k)]) / 2
    return x, y / med


def _meth_track(meth, chrom):
    if meth is None or chrom not in meth.get("pooled_adj", {}):
        return None
    from .meth import params as mp
    y = np.asarray(meth["pooled_adj"][chrom], float)
    return np.arange(len(y)) * mp.BIN + mp.BIN / 2, y


def _nan_runs(x, y, min_run=2):
    bad = ~np.isfinite(y)
    d = np.diff(np.concatenate([[0], bad.astype(int), [0]]))
    step = x[1] - x[0] if len(x) > 1 else 0
    return [(x[s] - step / 2, x[e - 1] + step / 2) for s, e in zip(np.flatnonzero(d == 1), np.flatnonzero(d == -1))
            if e - s >= min_run]


def _rel_label(r) -> str:
    return f"{REL_SHORT.get(r['relation'], r['relation'])}{', ' + r['conf'] if r['conf'] != 'none' else ''}"


def _text_mm(s: str, size: float) -> float:
    return TextPath((0, 0), s, size=size, prop=FontProperties(size=size)).get_extents().width / 72 * 25.4


def _gallery(ax, axb, rows, f, scale, has_hic, has_meth, Lmax, b_xmax):
    seq, tel, hic, meth = f["seq"], f["tel"], f["hic"], f["meth"]
    rdna = f.get("rdna") or {}
    arrays = seq.arrays
    min_w = Lmax * 0.0012                                   # sub-pixel arrays drawn at this width (legend)
    n = len(rows)
    lanes = {}
    top = 0.94
    if has_hic:
        lanes["hic"] = (top - 0.36, top)
        top -= 0.40
    if has_meth:
        lanes["meth"] = (top - 0.22, top)
        top -= 0.26
    has_te = bool(tel.get("density"))
    if has_te:
        lanes["te"] = (top - 0.12, top)
        top -= 0.14
    lanes["feat"] = (0.05, max(0.05 + 0.14, top))
    for i, r in enumerate(rows):
        y0 = n - 1 - i
        c, L = r["chrom"], r["L"]
        if r["call"]:
            s, e = r["call"]
            ax.add_patch(Rectangle((s, y0 + 0.02), max(e - s, min_w), 0.96, color=BAND[r["source"]], lw=0,
                                   zorder=0))
        if r["call"] and r["seq"] and r["call"] != r["seq"]:
            a = (r["seq"][0] + r["seq"][1]) / 2
            b = (r["call"][0] + r["call"][1]) / 2
            ax.annotate("", xy=(b, y0 + 0.03), xytext=(a, y0 + 0.03),
                        arrowprops=dict(arrowstyle="-|>", lw=0.6, color="black", mutation_scale=5), zorder=5)
        fl, fh = lanes["feat"]
        ax.plot([0, L], [y0 + fl + 0.02] * 2, color="#808080", lw=0.4, zorder=1)
        a = arrays[arrays.chrom == c]
        if len(a):
            verts, cols = [], []
            for s, e, role in zip(a.start.to_numpy(), a.end.to_numpy(), a.role.to_numpy()):
                w = max(e - s, min_w)
                verts.append([(s, y0 + fl + 0.02), (s + w, y0 + fl + 0.02), (s + w, y0 + fh - 0.04),
                              (s, y0 + fh - 0.04)])
                cols.append(ROLE_COLOR.get(role, OTHER))
            order = np.argsort([c_ != OTHER for c_ in cols])          # gray context first, colored on top
            ax.add_collection(PolyCollection([verts[j] for j in order], facecolors=[cols[j] for j in order],
                                             edgecolors="none", zorder=2))
        dens = tel.get("density", {}).get(c)
        if dens is not None and dens[1].max() > 0 and "te" in lanes:
            tl, th = lanes["te"]
            bw, cnt = dens
            xs = np.arange(len(cnt)) * bw
            h = np.minimum(cnt / scale["te"], 1.0) * (th - tl)
            ax.fill_between(xs, y0 + tl, y0 + tl + h, step="post", color=TE, alpha=0.8, lw=0, zorder=1.5)
        if r["seq"]:
            s, e = r["seq"]
            ax.add_patch(Rectangle((s, y0 + fl - 0.02), max(e - s, min_w), fh - fl + 0.02, fill=False,
                                   ec="black", lw=0.7, zorder=4))
        for tc, side, ts, te_, *_ in (f.get("telo") or {}).get("ends", []):
            if tc == c:
                ax.plot(0 if side == "start" else L, y0 + fh - 0.02, marker=">" if side == "start" else "<",
                        ms=3.2, color="black", mew=0, zorder=5, clip_on=False)
        for tc, ts, te_, _ in (f.get("telo") or {}).get("interstitial", []):
            if tc == c:
                ax.plot((ts + te_) / 2, y0 + fh + 0.06, marker="v", ms=3.0, color="black", mew=0, zorder=5,
                        clip_on=False)
        for arr, face in ((rdna.get("nors", ()), ROLE_COLOR["rDNA"]), (rdna.get("five_s", ()), "white")):
            for rc, rs, re_, *_ in arr:                          # 45S filled, 5S open
                if rc == c:
                    ax.plot((rs + re_) / 2, y0 + fh + 0.06, marker="D", ms=2.6, mew=0.6, mec=ROLE_COLOR["rDNA"],
                            mfc=face, ls="", zorder=5, clip_on=False)
        if has_hic:
            hl, hh = lanes["hic"]
            ax.plot([0, L], [y0 + hl] * 2, color="#D0D0D0", lw=0.3, zorder=1)
            for key, col, lw in (("denovo_S", DENOVO, 0.4), ("S", HIC, 0.7)):
                t = _track(hic, c, key)
                if t is not None:
                    yv = np.clip(t[1] / scale["hic"], 0, 1) * (hh - hl)
                    ax.plot(t[0], y0 + hl + yv, color=col, lw=lw, zorder=3, solid_joinstyle="round")
                    if key == "S":
                        for s, e in _nan_runs(t[0], t[1]):
                            ax.add_patch(Rectangle((s, y0 + hl), e - s, hh - hl, fill=False, hatch="////",
                                                   ec="#C8C8C8", lw=0, zorder=1))
        if has_meth:
            ml, mh = lanes["meth"]
            mid = (ml + mh) / 2
            ax.plot([0, L], [y0 + mid] * 2, color="#D0D0D0", lw=0.3, zorder=1)
            t = _meth_track(meth, c)
            if t is not None:
                yv = np.clip(t[1] / scale["meth"], -1, 1) * (mh - ml) / 2
                ax.plot(t[0], y0 + mid + yv, color=METH, lw=0.5, zorder=3)
                for s, e in _nan_runs(t[0], t[1], 4):
                    ax.add_patch(Rectangle((s, y0 + ml), e - s, mh - ml, fill=False, hatch="////", ec="#C8C8C8",
                                           lw=0, zorder=1))
        ax.text(-Lmax * 0.01, y0 + 0.5, c, ha="right", va="center", fontsize=6.5)
        for j, line in enumerate(f["lines"]):
            if line in r["support"]:
                axb.plot(j, y0 + 0.5, "o", ms=3.6, color="black", mew=0)
            elif line in r["against"]:
                axb.plot(j, y0 + 0.5, "x", ms=3.6, color="#D55E00", mew=0.9)
        axb.text(len(f["lines"]) - 0.35, y0 + 0.5, _rel_label(r), va="center", fontsize=6)
    ax.set_xlim(-Lmax * 0.005, Lmax * 1.005)
    ax.set_ylim(-0.1, n + 0.1)
    ax.set_yticks([])
    for sp in ("left", "right", "top"):
        ax.spines[sp].set_visible(False)
    ticks = _nice_ticks(Lmax)
    ax.set_xticks(ticks)
    ax.set_xticklabels([f"{t / 1e6:g}" for t in ticks])
    ax.set_xlabel("position (Mb)", labelpad=1)
    axb.set_xlim(-0.6, b_xmax)
    axb.set_ylim(-0.1, n + 0.1)
    axb.set_xticks(range(len(f["lines"])))
    axb.set_xticklabels([LINE_LABEL[x] for x in f["lines"]], rotation=90)
    axb.tick_params(axis="x", length=0, labeltop=True, labelbottom=False, pad=1)
    axb.set_yticks([])
    for sp in axb.spines.values():
        sp.set_visible(False)
    for j in range(len(LINES)):
        axb.axvline(j, color="#EBEBEB", lw=0.4, zorder=0)
    return lanes


def _nice_ticks(Lmax):
    for step in (1e6, 2e6, 5e6, 10e6, 20e6, 25e6, 50e6, 100e6, 200e6):
        if Lmax / step <= 8:
            return np.arange(0, Lmax + 1, step)
    return np.arange(0, Lmax + 1, 500e6)


def _offsets(ax, rows):
    """Each line's call (center and extent) relative to the sequence call center, one column per chromosome."""
    xs = range(len(rows))
    for x, r in zip(xs, rows):
        ref = r["seq"]
        if not ref:
            ax.plot(x, 0, marker="|", color="#BBBBBB", ms=4)
            continue
        m = (ref[0] + ref[1]) / 2
        rel = lambda iv: ((iv[0] - m) / 1e6, (iv[1] - m) / 1e6)
        a, b = rel(ref)
        ax.add_patch(Rectangle((x - 0.18, a), 0.36, max(b - a, 0.02), color="black", lw=0, zorder=3))
        for iv, col, dx in ((r["te"], TE, -0.32), (r["hic"], HIC, 0.32)):
            if iv:
                a, b = rel(iv)
                ax.plot([x + dx, x + dx], [a, b], color=col, lw=1.0, zorder=3)
                ax.plot(x + dx, (a + b) / 2, "o", color=col, ms=2.4, mew=0, zorder=4)
        if r["call"] and r["call"] != r["seq"]:
            a, b = rel(r["call"])
            ax.plot(x, (a + b) / 2, marker="D", ms=3.4, mfc="white", mec="black", mew=0.7, zorder=5)
    ax.autoscale_view()          # patches alone never rescale: a FASTA-only run calling every chromosome kept 0-1 Mb
    ax.axhline(0, color="#D0D0D0", lw=0.4, zorder=0)
    ax.set_xticks(list(xs))
    ax.set_xticklabels([r["chrom"] for r in rows], rotation=90, fontsize=5.5)
    ax.set_xlim(-0.8, len(rows) - 0.2)
    ax.set_ylabel("position relative to the\nsequence call center (Mb)")
    for sp in ("right", "top"):
        ax.spines[sp].set_visible(False)


def _overlay(ax, rows, hic, flank=10e6):
    """Anchored Hi-C score (x chromosome median) against distance from each integrated call's center."""
    traces = []
    for r in rows:
        t = _track(hic, r["chrom"], "S")
        if t is None or not r["call"]:
            continue
        m = (r["call"][0] + r["call"][1]) / 2
        d = (t[0] - m)
        keep = np.abs(d) <= flank
        moved = r["relation"] in ("offset", "functional_only")
        ax.plot(d[keep] / 1e6, t[1][keep], color="#D55E00" if moved else "#9E9E9E", lw=0.5 if moved else 0.35,
                alpha=0.95 if moved else 0.6, zorder=3 if moved else 2)
        traces.append((d, t[1]))
    if traces:
        b = hic["bins"].size
        grid = np.arange(-flank, flank + 1, b)
        vals = np.full((len(traces), len(grid)), np.nan)
        for i, (d, y) in enumerate(traces):
            idx = np.round((d - d[np.argmin(np.abs(d))]) / b).astype(int) + len(grid) // 2
            ok = (idx >= 0) & (idx < len(grid))
            vals[i, idx[ok]] = y[ok]
        with np.errstate(all="ignore"):
            med = np.nanmedian(vals, axis=0)
        ax.plot(grid / 1e6, med, color=HIC, lw=1.2, zorder=4)
    ax.axhline(1, color="#D0D0D0", lw=0.4, zorder=0)
    ax.axvline(0, color="#D0D0D0", lw=0.4, zorder=0)
    ax.set_xlim(-flank / 1e6, flank / 1e6)
    ax.set_xlabel("distance from the integrated call center (Mb)")
    ax.set_ylabel("anchored Hi-C score\n(x chromosome median)")
    for sp in ("right", "top"):
        ax.spines[sp].set_visible(False)


def _key(fig, y, has_hic, has_meth, roles, rdna=(False, False)):
    items = [Patch(fc=ROLE_COLOR[r], ec="none", label=f"{r.replace('_', '-')} satellite") for r in ROLE_COLOR
             if r in roles]
    items.append(Patch(fc=OTHER, ec="none", label="other satellite"))
    items += [Patch(fc="none", ec="black", lw=0.7, label="sequence call"),
              Patch(fc=BAND["sequence"], ec="none", label="integrated call (sequence)"),
              Patch(fc=BAND["hic"], ec="none", label="integrated call (Hi-C placed)"),
              Patch(fc=TE, alpha=0.8, ec="none", label="TE marker clade density")]
    if has_hic:
        items += [Line2D([], [], color=HIC, lw=0.9, label="anchored Hi-C score"),
                  Line2D([], [], color=DENOVO, lw=0.6, label="prior-free Hi-C score")]
    if has_meth:
        items.append(Line2D([], [], color=METH, lw=0.6, label="adjusted CHG (down = loss)"))
    items.append(Line2D([], [], ls="", marker=">", ms=3.4, color="black", label="telomere array"))
    items.append(Line2D([], [], ls="", marker="v", ms=3.4, color="black", label="interstitial telomeric array"))
    for has, face, label in zip(rdna, (ROLE_COLOR["rDNA"], "white"), ("45S rDNA array (NOR)", "5S rDNA array")):
        if has:
            items.append(Line2D([], [], ls="", marker="D", ms=3.0, mew=0.6, mec=ROLE_COLOR["rDNA"], mfc=face,
                                label=label))
    items += [Line2D([], [], ls="", marker="o", ms=3.4, color="black", label="line supports the call"),
              Line2D([], [], ls="", marker="x", ms=3.4, mew=0.9, color="#D55E00", label="line opposes the call"),
              Patch(fc="none", ec="#C8C8C8", hatch="////", lw=0, label="no data")]
    fig.legend(handles=items, loc="upper left", bbox_to_anchor=(0.01, y), ncol=4, frameon=False, fontsize=6,
               handlelength=1.6, columnspacing=1.2, handletextpad=0.5, borderaxespad=0)


def _pages(f):
    g = f["genome"]
    rows = _rows(f["calls"])
    has_hic, has_meth = f["hic"] is not None, f["meth"] is not None
    Lmax = max(r["L"] for r in rows)
    scale = dict(te=1.0, hic=4.0, meth=0.1)
    dens = [d[1].max() for d in f["tel"].get("density", {}).values() if len(d[1])]
    scale["te"] = max(np.percentile(dens, 90), 1) if dens else 1.0
    if has_hic:
        v = np.concatenate([t[1] for t in (_track(f["hic"], r["chrom"], "S") for r in rows) if t is not None]
                           or [np.ones(1)])
        scale["hic"] = max(float(np.nanpercentile(v, 99.5)), 2.0)
    if has_meth:
        v = np.concatenate([np.abs(t[1]) for t in (_meth_track(f["meth"], r["chrom"]) for r in rows) if t is not None]
                           or [np.full(1, 0.1)])
        scale["meth"] = max(float(np.nanpercentile(v, 99.5)), 0.02) if np.isfinite(v).any() else 0.1
    roles = set(f["seq"].arrays.role)
    rd = f.get("rdna") or {}
    rdna = (bool(rd.get("nors")), bool(rd.get("five_s")))
    f["scale"] = scale
    fold = any("fold" in (r["support"] | r["against"]) for r in rows)   # a fold column only where the fold acted
    f["lines"] = LINES[:3] + (["fold"] if fold else []) + LINES[3:]
    label_mm = 3 + 1.1 * max(len(r["chrom"]) for r in rows)            # rotated chromosome names under panel C
    # panel B ends 1 mm inside the page with its longest relation label; panel A takes the width left of it
    b_xmax = len(f["lines"]) - 0.35 + (max(_text_mm(_rel_label(r), 6) for r in rows) + 1.0) / B_COL_MM
    b_mm = (b_xmax + 0.6) * B_COL_MM
    b_left = 1 - (b_mm + 1) / (WIDTH / MM)
    pages = []
    chunks = [rows[i:i + ROWS_PER_PAGE] for i in range(0, len(rows), ROWS_PER_PAGE)]
    for pi, chunk in enumerate(chunks):
        last = pi == len(chunks) - 1
        row_mm = (6.6 if (has_hic and has_meth) else 5.4 if (has_hic or has_meth) else 4.2) + \
            (0.8 if f["tel"].get("density") else 0)
        gal_mm = row_mm * len(chunk)
        yb_mm = max(9.0, label_mm)
        bottom_mm = yb_mm + 53 if last else 0
        gtop = 27                                              # mm from the top: key, then letters, then rows
        height = (gtop + gal_mm + 12 + bottom_mm) * MM
        fig = plt.figure(figsize=(WIDTH, height))
        H = height / MM
        ax = fig.add_axes([0.11, 1 - (gtop + gal_mm) / H, b_left - 0.015 - 0.11, gal_mm / H])
        axb = fig.add_axes([b_left, 1 - (gtop + gal_mm) / H, b_mm / (WIDTH / MM), gal_mm / H])
        _gallery(ax, axb, chunk, f, scale, has_hic, has_meth, Lmax, b_xmax)
        _key(fig, 1 - 2 / H, has_hic, has_meth, roles, rdna)
        fig.text(0.01, 1 - 21 / H, "A", fontsize=9, fontweight="bold", va="bottom")
        fig.text(b_left - 0.02, 1 - 21 / H, "B", fontsize=9, fontweight="bold", va="bottom")
        if last:
            yb = yb_mm / H
            hb = 44 / H
            axc = fig.add_axes([0.11, yb, 0.5 if has_hic else 0.85, hb])
            _offsets(axc, rows)
            fig.text(0.01, yb + hb + 3 / H, "C", fontsize=9, fontweight="bold", va="bottom")
            if has_hic:
                axd = fig.add_axes([0.7, yb, 0.28, hb])
                _overlay(axd, rows, f["hic"])
                fig.text(0.63, yb + hb + 3 / H, "D", fontsize=9, fontweight="bold", va="bottom")
        pages.append(fig)
    return pages


def _focus(ptab):
    """The chromosome where the plants differ most from the reference call (shifted first), or None."""
    w = {"shifted": 3, "call_absent": 2, "possible_shift": 1, "conflict": 1}
    s = ptab.assign(w=ptab.verdict.map(w).fillna(0)).groupby("chrom").w.sum()
    s = s[s > 0]
    return None if s.empty else s.sort_values(ascending=False, kind="stable").index[0]


def _verdict_grid(ax, ptab, chroms, prov):
    plants_ = sorted(ptab.plant.unique(), key=lambda p: ({"reference": 0, "unknown": 1}.get(prov.get(p), 2), p))
    by = {(r.plant, r.chrom): r for r in ptab.itertuples(index=False)}
    for i, p in enumerate(plants_):
        y = len(plants_) - 1 - i
        for j, c in enumerate(chroms):
            r = by.get((p, c))
            v = r.verdict if r is not None else "no_data"
            ax.add_patch(Rectangle((j + 0.06, y + 0.08), 0.88, 0.84, fc=VERDICT_COLOR.get(v, "white"),
                                   ec="#BDBDBD" if v in ("unclear", "no_data", "no_call") else "none", lw=0.4,
                                   hatch="////" if v == "no_data" else None, zorder=2))
            if r is not None and v in ("shifted", "possible_shift", "conflict") and r.shift_bp not in (".", None) \
                    and np.isfinite(float(r.shift_bp)):
                ax.text(j + 0.5, y + 0.5, f"{float(r.shift_bp) / 1e6:.1f}", ha="center", va="center", fontsize=5,
                        color="white" if v == "shifted" else "black", zorder=3)
        ax.text(-0.15, y + 0.5, f"{p} ({PROV_SHORT.get(prov.get(p), 'origin unknown')})", ha="right", va="center",
                fontsize=6)
    ax.set_xlim(0, len(chroms))
    ax.set_ylim(0, len(plants_))
    ax.set_xticks(np.arange(len(chroms)) + 0.5)
    ax.set_xticklabels(chroms, rotation=90, fontsize=5.5)
    ax.tick_params(axis="x", length=0, labeltop=True, labelbottom=False, pad=1)
    ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)
    return plants_


def _plant_rows(ax, f, chrom, plants_, prov):
    """F: the reference plant's evidence on top, then one row per plant: its adjusted-CHG track (down = loss), the
    sites it was scored at, its call site and its best centromere-like site away from the call."""
    from .meth import params as mp
    ptab, sites, meth, seq, tel, hic = f["plants"], f["psites"], f["meth"], f["seq"], f["tel"], f["hic"]
    r0 = next(r for r in _rows(f["calls"]) if r["chrom"] == chrom)
    L = r0["L"]
    d = ptab[ptab.chrom == chrom]
    pts = [x for x in (r0["call"] or ()) if x is not None]
    pts += [float(v) for v in d.alt_start.tolist() + d.alt_end.tolist() if float(v) >= 0]
    lo, hi = (min(pts), max(pts)) if pts else (0, L)
    pad = max(5e6, 0.15 * (hi - lo))
    lo, hi = max(0, lo - pad), min(L, hi + pad)
    if hi - lo < 20e6:
        mid = (lo + hi) / 2
        lo, hi = max(0, mid - 10e6), min(L, mid + 10e6)
    n = len(plants_)
    ref_h = 1.6
    ytop = n + ref_h
    if r0["call"]:
        s, e = r0["call"]
        ax.add_patch(Rectangle((s, 0), e - s, ytop, color=BAND[r0["source"]], lw=0, zorder=0))
    # reference evidence
    y0 = n + 0.15
    a = seq.arrays[(seq.arrays.chrom == chrom) & (seq.arrays.end > lo) & (seq.arrays.start < hi)]
    min_w = (hi - lo) * 0.0015
    for s, e, role in zip(a.start.to_numpy(), a.end.to_numpy(), a.role.to_numpy()):
        ax.add_patch(Rectangle((s, y0), max(e - s, min_w), 0.32, fc=ROLE_COLOR.get(role, OTHER), ec="none", zorder=3))
    dens = tel.get("density", {}).get(chrom)
    if dens is not None and dens[1].max() > 0:
        bw, cnt = dens
        xs = np.arange(len(cnt)) * bw
        k = (xs >= lo - bw) & (xs <= hi)
        top = max(float(cnt[k].max()), 1.0) if k.any() else 1.0
        ax.fill_between(xs[k], y0 + 0.4, y0 + 0.4 + 0.45 * cnt[k] / top, step="post", color=TE, alpha=0.8, lw=0,
                        zorder=2)
    tr = _track(hic, chrom, "S")
    if tr is not None:
        k = (tr[0] >= lo) & (tr[0] <= hi)
        yv = np.clip(tr[1][k] / f["scale"]["hic"], 0, 1) * 0.55
        ax.plot(tr[0][k], y0 + 0.9 + yv, color=HIC, lw=0.7, zorder=3)
        for s, e in _nan_runs(tr[0][k], tr[1][k], 1):
            ax.add_patch(Rectangle((s, y0 + 0.88), e - s, 0.6, fill=False, hatch="////", ec="#C8C8C8", lw=0,
                                   zorder=1))
    if r0["seq"]:
        s, e = r0["seq"]
        ax.add_patch(Rectangle((s, y0 - 0.04), e - s, 0.4, fill=False, ec="black", lw=0.7, zorder=4))
    x_lab = lo - (hi - lo) * 0.01
    for yy, name in ((y0 + 0.16, "satellite arrays"), (y0 + 0.58, "TE marker"), (y0 + 1.12, "Hi-C contact")):
        ax.text(x_lab, yy, name, ha="right", va="center", fontsize=5.5, color="#505050")
    # plants
    vals = []
    for p in plants_:
        y = np.asarray(meth["adj"].get(p, {}).get(chrom, []), float)
        if len(y):
            x = np.arange(len(y)) * mp.BIN + mp.BIN / 2
            vals.append(np.abs(y[(x >= lo) & (x <= hi)]))
    v = np.concatenate(vals) if vals else np.zeros(1)
    sc = max(float(np.nanpercentile(v, 99)), 0.02) if np.isfinite(v).any() else 0.1
    srow = sites[sites.chrom == chrom] if sites is not None else None
    for i, p in enumerate(plants_):
        yb = n - 1 - i
        mid = yb + 0.5
        ax.plot([lo, hi], [mid, mid], color="#D0D0D0", lw=0.3, zorder=1)
        lib = (f.get("libraries") or {}).get(p)
        lt = _track(lib, chrom, "S") if lib is not None else None
        if lt is not None:
            k = (lt[0] >= lo) & (lt[0] <= hi)
            ax.plot(lt[0][k], yb + 0.5 + np.clip(lt[1][k] / f["scale"]["hic"], 0, 1) * 0.45, color=HIC, lw=0.6,
                    alpha=0.85, zorder=2.5)
        y = np.asarray(meth["adj"].get(p, {}).get(chrom, []), float)
        if len(y):
            x = np.arange(len(y)) * mp.BIN + mp.BIN / 2
            k = (x >= lo) & (x <= hi)
            ax.plot(x[k], mid + np.clip(y[k] / sc, -1, 1) * 0.42, color=METH, lw=0.5, zorder=3)
            for s, e in _nan_runs(x[k], y[k], 4):
                ax.add_patch(Rectangle((s, yb + 0.06), e - s, 0.88, fill=False, hatch="////", ec="#C8C8C8", lw=0,
                                       zorder=1))
        if srow is not None:
            mine = srow[srow["sample"] == p]
            for s in mine.itertuples(index=False):
                if lo <= s.mid <= hi:
                    ax.plot([s.mid, s.mid], [yb + 0.86, yb + 0.96], color="#808080", lw=0.5, zorder=2)
        rr = d[d.plant == p]
        if len(rr):
            rr = rr.iloc[0]
            if rr.alt_start >= 0:
                ax.plot((rr.alt_start + rr.alt_end) / 2, yb + 0.12, marker="^", ms=4.2, mew=0,
                        color=VERDICT_COLOR.get(rr.verdict if rr.verdict != "shared" else "possible_shift"), zorder=5)
            lab = VERDICT_LABEL.get(rr.verdict, rr.verdict)
            if rr.verdict in ("shifted", "possible_shift", "conflict") and rr.shift_bp not in (".", None) \
                    and np.isfinite(float(rr.shift_bp)):
                lab += f", {float(rr.shift_bp) / 1e6:.1f} Mb"
            ax.text(hi + (hi - lo) * 0.01, mid, lab, ha="left", va="center", fontsize=6,
                    color=VERDICT_COLOR.get(rr.verdict) if rr.verdict in ("shifted", "possible_shift", "conflict")
                    else "black")
        ax.text(lo - (hi - lo) * 0.01, mid, p, ha="right", va="center", fontsize=6.5)
    ax.set_xlim(lo, hi)
    ax.set_ylim(0, ytop + 0.1)
    ax.set_yticks([])
    for sp in ("left", "right", "top"):
        ax.spines[sp].set_visible(False)
    ticks = [x for x in np.arange(0, L + 1, _nice_step(hi - lo)) if lo <= x <= hi]
    ax.set_xticks(ticks)
    ax.set_xticklabels([f"{x / 1e6:g}" for x in ticks])
    ax.set_xlabel(f"{chrom} position (Mb)", labelpad=1)
    return sc


def _nice_step(span):
    for step in (1e6, 2e6, 5e6, 10e6, 20e6, 50e6):
        if span / step <= 8:
            return step
    return 100e6


def _plants_page(f):
    ptab = f.get("plants")
    if ptab is None or not len(ptab):
        return None
    prov = {n: v["verdict"] for (k, n), v in f["verdicts"].items() if k == "meth"}
    prov.update({n: v["verdict"] for (k, n), v in f["verdicts"].items() if k == "hic" and n not in prov})
    chroms = list(f["calls"].chrom)
    focus = _focus(ptab)
    n_pl = ptab.plant.nunique()
    has_f = focus is not None and f["meth"] is not None
    grid_mm = 4.6 * n_pl
    top = 17                                               # mm: verdict key, letter E, chromosome labels
    f_mm = (12 + 8.5 * n_pl + 9) if has_f else 0
    H = top + grid_mm + (13 + f_mm if has_f else 0) + 4
    fig = plt.figure(figsize=(WIDTH, H * MM))
    axe = fig.add_axes([0.22, 1 - (top + grid_mm) / H, 0.76, grid_mm / H])
    plants_ = _verdict_grid(axe, ptab, chroms, prov)
    fig.text(0.01, 1 - 10 / H, "E", fontsize=9, fontweight="bold", va="bottom")
    items = [Patch(fc=VERDICT_COLOR[v], ec="#BDBDBD" if v in ("unclear", "no_data") else "none", lw=0.4,
                   hatch="////" if v == "no_data" else None, label=VERDICT_LABEL[v])
             for v in ("shared", "shifted", "possible_shift", "call_absent", "conflict", "unclear", "no_data")
             if v in set(ptab.verdict) or v in ("shared", "shifted")]
    fig.legend(handles=items, loc="upper left", bbox_to_anchor=(0.04, 1 - 1.5 / H), ncol=len(items), frameon=False,
               fontsize=6, handlelength=1.2, columnspacing=1.0, handletextpad=0.4, borderaxespad=0)
    f["focus"] = focus
    if has_f:
        yb = 4 / H
        axf = fig.add_axes([0.13, yb + 9 / H, 0.69, (f_mm - 9) / H])
        f["focus_scale"] = _plant_rows(axf, f, focus, plants_, prov)
        fig.text(0.01, yb + f_mm / H + 5 / H, "F", fontsize=9, fontweight="bold", va="bottom")
        keys = [Line2D([], [], color=METH, lw=0.6, label="plant's adjusted CHG (down = loss)"),
                Line2D([], [], ls="", marker="^", ms=4, mew=0, color=VERDICT_COLOR["shifted"],
                       label="centromere-like away from the call"),
                Line2D([], [], color="#808080", lw=0.6, label="site scored"),
                Patch(fc=BAND["sequence"], ec="none", label="integrated call"),
                Patch(fc="none", ec="#C8C8C8", hatch="////", lw=0, label="no data")]
        if f.get("libraries"):
            keys.insert(1, Line2D([], [], color=HIC, lw=0.6, label="Hi-C (reference or the plant's library)"))
        fig.legend(handles=keys, loc="lower left", bbox_to_anchor=(0.04, yb + f_mm / H + 0.5 / H), ncol=len(keys),
                   frameon=False, fontsize=6, handlelength=1.4, columnspacing=1.0, handletextpad=0.4,
                   borderaxespad=0)
    return fig


def _species(name: str) -> str:
    return name.replace("_", " ")


def legends(f) -> tuple[str, str]:
    g, t = f["genome"], f["calls"]
    called = t[t.start != "."]
    rel = t.relation.value_counts().to_dict()
    n = len(t)
    moved = rel.get("offset", 0) + rel.get("functional_only", 0)
    hic, meth, tel = f["hic"], f["meth"], f["tel"]
    head = (f"In {_species(g.name)}, {len(called)} of {n} chromosomes have a centromere call; {moved} sit away from "
            f"(or without) a satellite array where Hi-C contact with the other centromeres"
            f"{' and methylation' if meth is not None else ''} place them, and {rel.get('contested', 0)} "
            f"stay contested. Satellite position and centromere function can part ways, so each line is shown "
            f"on its own.")
    short = [head,
             "A, Every chromosome at one scale: satellite arrays by family role, retrotransposon marker density, "
             "Hi-C contact with the other centromeres and the methylation track, with the sequence call (box) and "
             "the integrated call (band).",
             "B, Which evidence lines support or oppose each integrated call.",
             "C, Where each line places the centromere relative to the satellite call."]
    if hic is not None:
        short.append("D, Hi-C contact with the other chromosomes' centromeres peaks at the integrated calls "
                     "(orange: calls placed away from the satellite).")
    lines_used = ", ".join(g.lines)
    long = [f"Terms. Sequence call: the satellite (or retrotransposon) centromere. Anchored Hi-C score: mean "
            f"observed/expected contact of each bin with windows at the other chromosomes' sequence calls (the "
            f"chromosome's own call is never used), over the chromosome median. Integrated call: the sequence call, "
            f"unless Hi-C that agrees with the sequence calls of the other chromosomes (gate: majority, "
            f"Poisson-binomial P <= 0.01) places a strong peak elsewhere; a high sequence call moves only with a "
            f"second line (retrotransposon core or methylation) at the Hi-C peak. Where that gate fails, the arm fold "
            f"(the anti-diagonal of each chromosome's own contact map in a Rabl nucleus) can validate the prior-free "
            f"Hi-C peaks instead (majority, P <= 0.01); the two then move calls they both contradict. Lines used "
            f"here: {lines_used}."]
    a = (f"A, {n} chromosomes in assembly order, to one base-pair scale. Each satellite array "
         f"is a box colored by its family's role in this genome ({len(f['seq'].arrays):,} arrays; arrays too small "
         f"to see are drawn at the minimum visible width). The orange area is the density of the "
         f"{tel['marker'][0] + ' ' if tel.get('marker') else ''}retrotransposon marker clade.")
    if hic is not None:
        h = hic["table"]
        a += (f" Blue: anchored Hi-C score at {hic['bin_size'] / 1e3:g}-kb bins (0 to {f['scale']['hic']:.1f} x "
              f"the chromosome median over the track height, clipped above); gray: the prior-free Hi-C score. The gate passed on "
              f"{int(h.gate_pass.sum())} of {len(h)} chromosomes; the anchored peak lies on the sequence call on "
              f"{int(h.anch_agrees_seq.sum())} of {int(h.seq_start.ge(0).sum())}, the prior-free peak on "
              f"{int(h.denovo_agrees_seq.sum())}.")
    if meth is not None:
        sig = ",".join(f"{r.context} {r.cls} {r.direction}" for r in meth["signature"].itertuples()) or "none"
        a += (f" Black line: composition-adjusted CHG methylation pooled over {len(meth['samples'])} plants "
              f"(0 = like the same sequence class nearby; down = loss). The centromere signature learned from the "
              f"calls of the other chromosomes: {sig}.")
    telo = f.get("telo") or {}
    if telo.get("unit"):
        k = len(telo.get("interstitial", []))
        its = (f"downward triangles: {k} array{'s' if k != 1 else ''} of the unit inside chromosomes (interstitial "
               f"telomeric sequences)" if k else "none inside chromosomes")
        a += (f" Triangles: telomere arrays of the unit {telo['unit']} (found de novo; {len(telo['ends'])} of "
              f"{2 * n} chromosome ends); {its}.")
    else:
        a += " No telomere unit was enriched in tandem arrays at the chromosome ends."
    rd = f.get("rdna") or {}
    n45, n5 = len(rd.get("nors", ())), len(rd.get("five_s", ()))
    a += (f" Green diamonds: rDNA arrays, {n45} of 45S (filled; the 18S motif in tandem) and {n5} of 5S (open; whole "
          f"5S genes in tandem)." if n45 or n5 else " No 45S or 5S rDNA array was found.")
    a += " Hatching: no data (masked or uncovered bins). Bands: integrated calls; arrows: sequence call to moved call."
    long.append(a)
    long.append("B, One column per evidence line: filled circle, the line supports the integrated call (satellite: "
                ">= 5 kb of centromeric-family arrays inside it; TE: the marker core agrees; Hi-C: gated and the "
                "anchored peak within one bin; methylation: centromere-like under the leave-one-chromosome-out "
                "signature; fold, shown only where it acted: the call lies in the arm-fold core); cross, a gated "
                "strong Hi-C peak, a CRM core, the arm fold or a validated prior-free Hi-C peak points elsewhere. "
                "Text: relation and confidence.")
    long.append("C, Black bar: the sequence call, centered at 0. Orange: the TE marker core; blue: the anchored Hi-C "
                "core with its peak; open diamond: the integrated call when it differs from the sequence call.")
    if hic is not None:
        long.append("D, Anchored Hi-C score in the 10 Mb on each side of every integrated call, one trace per "
                    "chromosome (gray: kept sequence calls; orange: calls placed by Hi-C); blue: the median. A peak "
                    "at 0 means the calls sit where the chromosomes' centromeres meet in the nucleus.")
    ptab = f.get("plants")
    if ptab is not None and len(ptab):
        _plant_legends(f, short, long)
    long.append(f"Exceptions are listed in {g.name}.centromeres.tsv (relation, against, flags, seq_flags) and every "
                f"candidate feature in {g.name}.loci.tsv.")
    long.append(f"Analysis: {__import__('socket').gethostname()}:{f['out']} · Script: cenfuse/plot.py · "
                f"{date.today()}")
    return "\n\n".join(short), "\n\n".join(long)


def _plant_legends(f, short: list, long: list) -> None:
    ptab, sites, ver = f["plants"], f.get("psites"), f["verdicts"]
    g = f["genome"]
    cnt = ptab.verdict.value_counts().to_dict()
    n_pl = ptab.plant.nunique()
    prov = {}
    for (k, n), v in ver.items():
        prov.setdefault(n, v["verdict"])
    others = sorted(p for p in ptab.plant.unique() if prov.get(p) == "other")
    focus = f.get("focus")
    short.append(f"E, Each of {n_pl} plants, chromosome by chromosome: is its own evidence consistent with the "
                 f"reference call (gray), or does its centromere seem to sit elsewhere (orange; numbers: Mb from the "
                 f"call)? Black: the plant lacks the DNA of the call.")
    if focus is not None and f.get("meth") is not None:
        short.append(f"F, {focus}, plant by plant: each plant's methylation track (centromeric loss points down) "
                     f"beneath the reference plant's evidence; triangles mark where a plant looks centromeric away "
                     f"from the call.")
    n_tests = int(sites.test.isin(["loss", "gain"]).sum()) if sites is not None and "test" in sites else 0
    long.append(
        f"Plants and provenance. Each input was tested for whether it comes from the reference plant (the assembled "
        f"genotype): methylation samples by fixed differences (symmetric CHG sites whose rows are missing on both "
        f"strands where reads span them, so that read errors, which strike one strand, do not count) that another "
        f"sample lacks; every input by segments it lacks while another input covers them. "
        f"{len(others)} of {n_pl} plants come from another plant ({', '.join(others) or 'none'}); only "
        f"reference-plant inputs and inputs of unknown origin decide the reference calls "
        f"(methylation samples deciding here: {f['calls'].meth_samples.iat[0] if 'meth_samples' in f['calls'] else '.'}).")
    long.append(
        f"E, One row per plant (its provenance in brackets), one column per chromosome. Gray: consistent with the "
        f"reference call (centromere-like methylation there, or the plant's own Hi-C peak there or where the "
        f"reference library peaks); methylation that looks centromeric is not proof, since some non-centromeric "
        f"places look so in every plant. Orange: a dip at the call and centromere-like methylation at a site away from it, with between-plant "
        f"support (the plant less centromeric than the others at the call, or more at that site; signature-score "
        f"differences against their spread at the other chromosomes' calls, Student's t prediction interval, "
        f"Benjamini-Hochberg q <= 0.05 over {n_tests} tests), or a strong Hi-C peak of the plant's own library away "
        f"from the call where the reference library peaks elsewhere. Light orange: the methylation pattern without "
        f"between-plant support. Black: the call lies in a segment the plant lacks. Purple: the plant's lines "
        f"disagree. Pale: unclear; hatched: no data at the call. Of {len(ptab)} cells: {cnt.get('shared', 0)} at the "
        f"call, {cnt.get('shifted', 0)} shifted, {cnt.get('possible_shift', 0)} possible shifts, "
        f"{cnt.get('call_absent', 0)} call regions absent, {cnt.get('conflict', 0)} with disagreeing lines, "
        f"{cnt.get('unclear', 0)} unclear, {cnt.get('no_data', 0)} without data, {cnt.get('no_call', 0)} on "
        f"chromosomes without a reference call.")
    if focus is not None and f.get("meth") is not None:
        long.append(
            f"F, {focus}, chosen as the chromosome with the most shift verdicts. Top: the reference plant's evidence "
            f"(satellite arrays by family role, the retrotransposon marker density in orange, the anchored Hi-C score "
            f"in blue, the sequence call boxed); the band is the integrated call. Below, one row per plant: its "
            f"composition-adjusted CHG track (one scale for all plants, +-{f.get('focus_scale', 0.1):.2f}; down = "
            f"loss), gray ticks at the sites every plant was scored at, a triangle at its best centromere-like site "
            f"away from the call, and its verdict at the right. Hatching: no data (masked bins or DNA the plant lacks).")


def _text_page(text: str):
    paras = [textwrap.fill(p, 100) for p in text.split("\n\n")]
    body = "\n\n".join(paras)
    nlines = body.count("\n") + 1
    h = nlines * 4.7 + 12
    fig = plt.figure(figsize=(WIDTH, h * MM))
    fig.text(0.04, 1 - 6 / h, body, va="top", ha="left", fontsize=10, linespacing=1.35)
    return fig


def overview(figs: list, out: Path) -> None:
    out = Path(out)
    with PdfPages(out / "cenfuse_overview.pdf") as pdf:
        for f in figs:
            f.setdefault("out", out)
            pages = _pages(f)
            for i, fig in enumerate(pages):
                pdf.savefig(fig)
                if i == 0:
                    fig.savefig(out / f"{f['genome'].name}.overview.png", dpi=300)
                plt.close(fig)
            pp = _plants_page(f)
            if pp is not None:
                pdf.savefig(pp)
                pp.savefig(out / f"{f['genome'].name}.plants.png", dpi=300)
                plt.close(pp)
            for txt in legends(f):
                fig = _text_page(txt)
                pdf.savefig(fig)
                plt.close(fig)
    log.info(f"overview: {out / 'cenfuse_overview.pdf'}")
