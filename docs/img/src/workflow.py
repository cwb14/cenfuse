#!/usr/bin/env python3
"""How cenfuse works: the workflow figure of the README (docs/img/workflow.png and .pdf).

A and C are a schematic of the method (cenfuse 0.16.0); B is real data: the one-per-chromosome test of the satellite
line on rice and tomato, read from workflow_panelB.tsv beside this script. The script reads no other file, so the repo
rebuilds the figure on its own:

  python docs/img/src/workflow.py [--numbers N.tsv]     # writes docs/img/workflow.png (figure) and .pdf (+ legend)

workflow_panelB.tsv was frozen from cenfuse 0.16.0's FASTA-only benchmark runs of rice (T2T Nipponbare) and tomato
(Heinz 1706) by an extraction step kept with the benchmark runs, outside the repo: for each genome, every array of the
satellite families the satellite engine called centromeric, merged into clusters of arrays within 500 kb
(cenfuse.sat.classify.make_clusters), the chromosome lengths, and the test's verdict and 10,000 random draws
(cenfuse.sat.karyotype.satellite_test, null_covered), each checked against the runs' own test results. Its header
says what each row holds. With the same fonts and matplotlib, a rebuild is byte-identical (fixed date, no PDF
creation time).
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.patheffects as pe  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages  # noqa: E402
from matplotlib.font_manager import FontProperties, fontManager  # noqa: E402
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch, PathPatch, Polygon, Rectangle  # noqa: E402
from matplotlib.path import Path as MPath  # noqa: E402
from matplotlib.textpath import TextToPath  # noqa: E402

HERE = Path(__file__).resolve().parent             # docs/img/src
IMG = HERE.parent                                  # docs/img
PANEL_B = HERE / "workflow_panelB.tsv"
FIG_DATE = "2026-11-14"                            # when panel B's data were frozen (cenfuse 0.16.0); fixed, so a
                                                   # rebuild is byte-identical

# ---------------------------------------------------------------- style
FONT = "Arial" if any(f.name == "Arial" for f in fontManager.ttflist) else "Liberation Sans"
plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": [FONT, "DejaVu Sans"],
                     "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none"})
A_GROW = 11.0                                      # 0.16: panel A's extra row (repeats), mm
BC_GROW = 7.0                                      # 0.16: panel C's extra row (flag), mm
W, H = 180.0, 129.0 + A_GROW + BC_GROW             # page, mm
A_SHIFT = H - 126.0                                # panel A's frame (rows in mm) sits this much higher on the page
BC_DY = BC_GROW                                    # panels B and C: their band is this much taller than in 0.15
FS = 7.5                                           # smallest words on the figure (pt)
FS_NAME = 8.0                                      # line names, stage cues
FS_LET = 10.0                                      # panel letters

# palette: cenfuse/plot.py (ROLE_COLOR, HIC, TE, METH, AT_MARK, CIS_MARK, OTHER, RDNA_MARK)
SAT, KNOB, SUBTEL = "#D55E00", "#CC79A7", "#56B4E9"
HIC, TE, METH, AT, CIS = "#0072B2", "#E69F00", "#303030", "#009E73", "#6A3D9A"
REP = "#882255"                                    # 0.16 repeat lines (density, families): Tol's wine, apart from every
REP_CTX = "#BDBDBD"                                # Okabe-Ito hue already used; windows that are not the peak: gray
RDNA_MARK, OTHER = "#8C8C8C", "#B0B0B0"
TE_TXT = "#B87A00"                                 # the LTR amber, darkened for text only (contrast on white)
TE_LTR, TE_INT = "#9A6200", "#E69F00"              # young element: LTRs, internal region
OLD_LTR, OLD_INT = "#D9B46A", "#F6E2B5"            # old element
CONF = {"high": "#1E1E1E", "medium": "#7A7A7A", "low": "#BEBEBE"}   # confidence ramp (luminance only)
CHR_FILL, CHR_EDGE = "#ECECEC", "#9A9A9A"
BAND = "#F0F0F0"
CHANCE = "#CDCDCD"                                 # B: the chance distribution
INK = "#1E1E1E"
CUE = "#6B6B6B"
LINK = "#8C8C8C"


# ---------------------------------------------------------------- primitives (mm coordinates)
def bar_h(ax, x0, x1, y, h=1.5, fill=CHR_FILL, edge=CHR_EDGE, z=2, lw_=0.6, ls="-"):
    ax.add_patch(FancyBboxPatch((x0, y - h / 2), x1 - x0, h, boxstyle=f"round,pad=0,rounding_size={h / 2}",
                                fc=fill, ec=edge, lw=lw_, ls=ls, zorder=z))


def bar_v(ax, x, y0, y1, w=1.5, fill=CHR_FILL, edge=CHR_EDGE, z=2, lw_=0.6):
    ax.add_patch(FancyBboxPatch((x - w / 2, y0), w, y1 - y0, boxstyle=f"round,pad=0,rounding_size={w / 2}",
                                fc=fill, ec=edge, lw=lw_, zorder=z))


def rect(ax, x0, y0, w, h, fc, ec="none", lw_=0.6, z=3, ls="-"):
    ax.add_patch(Rectangle((x0, y0), w, h, fc=fc, ec=ec, lw=lw_, zorder=z, ls=ls))


def call_box(ax, xc, y, conf, w=2.6, h=2.3, z=6, ls="-", ec=INK, lw_=0.7):
    rect(ax, xc - w / 2, y - h / 2, w, h, CONF[conf] if conf else "white", ec=ec, lw_=lw_, z=z, ls=ls)


def tri_down(ax, x, y, s, fc, z=7):
    """downward triangle, tip at (x, y): where a tested line points"""
    ax.add_patch(Polygon([(x, y), (x - s / 2, y + s * 0.95), (x + s / 2, y + s * 0.95)], closed=True, fc=fc,
                         ec="none", zorder=z))


def unit_arrow(ax, x, y0, L, h, fc, z=5, rev=False):
    tip = min(h * 0.55, L * 0.45)
    if not rev:
        pts = [(x, y0), (x + L - tip, y0), (x + L, y0 + h / 2), (x + L - tip, y0 + h), (x, y0 + h)]
    else:
        pts = [(x + L, y0), (x + tip, y0), (x, y0 + h / 2), (x + tip, y0 + h), (x + L, y0 + h)]
    ax.add_patch(Polygon(pts, closed=True, fc=fc, ec="none", zorder=z))


def ltr_element(ax, x, y, L, h, ltr_c, int_c, z=5, rev=False):
    """LTR retrotransposon: internal region with an LTR at each end, both pointing the element's way"""
    ltr = max(L * 0.2, 0.45)
    rect(ax, x + ltr * 0.9, y - h * 0.32, L - 1.8 * ltr, h * 0.64, int_c, z=z)
    unit_arrow(ax, x, y - h / 2, ltr, h, ltr_c, z=z + 0.1, rev=rev)
    unit_arrow(ax, x + L - ltr, y - h / 2, ltr, h, ltr_c, z=z + 0.1, rev=rev)


def arrow(ax, p0, p1, color, lw_=1.0, z=4, rad=0.0, head=4.5, ls="-"):
    ax.add_patch(FancyArrowPatch(p0, p1, arrowstyle=f"-|>,head_length={head * 0.11},head_width={head * 0.07}",
                                 connectionstyle=f"arc3,rad={rad}", color=color, lw=lw_, zorder=z, ls=ls,
                                 shrinkA=0, shrinkB=0, mutation_scale=10, capstyle="butt", joinstyle="miter"))


def curve(ax, verts, color, lw_=1.0, z=4, ls="-", halo=False):
    """polyline (2 points) or cubic Bezier (4 points)"""
    codes = [MPath.MOVETO] + ([MPath.LINETO] if len(verts) == 2 else [MPath.CURVE4] * 3)
    p = PathPatch(MPath(verts, codes), fc="none", ec=color, lw=lw_, zorder=z, ls=ls, capstyle="butt")
    if halo:
        p.set_path_effects([pe.Stroke(linewidth=lw_ + 2.4, foreground=BAND if halo == "band" else "white"),
                            pe.Normal()])
    ax.add_patch(p)


def head(ax, x, y, color, s=1.3, direction="right", z=6):
    d = {"right": [(x, y), (x - s, y + s * 0.42), (x - s, y - s * 0.42)],
         "left": [(x, y), (x + s, y + s * 0.42), (x + s, y - s * 0.42)],
         "down": [(x, y), (x - s * 0.42, y + s), (x + s * 0.42, y + s)],
         "up": [(x, y), (x - s * 0.42, y - s), (x + s * 0.42, y - s)]}[direction]
    ax.add_patch(Polygon(d, closed=True, fc=color, ec="none", zorder=z))


def check(ax, x, y, s, color=INK, lw_=1.6, z=8):
    ax.plot([x - s * 0.5, x - s * 0.12, x + s * 0.55], [y + s * 0.02, y - s * 0.4, y + s * 0.5], color=color, lw=lw_,
            solid_capstyle="round", solid_joinstyle="round", zorder=z)


def cross(ax, x, y, s, color=INK, lw_=1.6, z=8):
    for sgn in (1, -1):
        ax.plot([x - s * 0.42, x + s * 0.42], [y - sgn * s * 0.42, y + sgn * s * 0.42], color=color, lw=lw_,
                solid_capstyle="round", zorder=z)


def txt(ax, x, y, s, size=FS, color=INK, ha="left", va="center", weight="normal", style="normal", z=10, **kw):
    return ax.text(x, y, s, fontsize=size, color=color, ha=ha, va=va, fontweight=weight, fontstyle=style, zorder=z,
                   **kw)


_T2P = TextToPath()


def text_mm(s, size, style="normal", weight="normal"):
    """width of a string in mm"""
    w, _, _ = _T2P.get_text_width_height_descent(s, FontProperties(family=FONT, size=size, style=style,
                                                                   weight=weight), ismath=False)
    return w * 25.4 / 72


def page_icon(ax, x, y, w, h, ec=INK, fill="white", lw_=0.9, z=5):
    f = min(w, h) * 0.32
    ax.add_patch(Polygon([(x, y), (x + w, y), (x + w, y + h - f), (x + w - f, y + h), (x, y + h)], closed=True,
                         fc=fill, ec=ec, lw=lw_, zorder=z, joinstyle="miter"))
    ax.plot([x + w - f, x + w - f, x + w], [y + h, y + h - f, y + h - f], color=ec, lw=lw_ * 0.8, zorder=z + 0.1)


def folder_icon(ax, x, y, w, h, ec, ls=(0, (2.4, 1.3)), lw_=1.0, z=5):
    tw, th = w * 0.38, h * 0.17
    ax.add_patch(Polygon([(x, y), (x + w, y), (x + w, y + h), (x + tw + th, y + h), (x + tw, y + h + th),
                          (x, y + h + th)], closed=True, fc="white", ec=ec, lw=lw_, ls=ls, zorder=z,
                         joinstyle="round"))


# ---------------------------------------------------------------- panel A: the method
X_IN = 2.0                 # inputs column, left edge
X_G0, X_G1 = 24.0, 79.0    # evidence: names from X_G0, glyphs end at X_G1
XB0 = 27.5                 # glyph chromosomes start
X_BAND0, X_BAND1 = 84.0, 106.0
X_HEAD1 = 125.0            # the test region's header (emblem) reaches over the integration column
Y_HEAD0 = 105.6
X_NODE = 98.5
X_LAD0, X_LAD1 = 109.8, 122.2       # the confidence rungs, inside the integration frame
X_FR0, X_FR1 = 107.9, 124.1
Y_LAD = 71.0
Y_TOP = 124.0
_D = A_GROW                         # rows below the new repeats row move down by this much (0.15 positions + _D)
ROW = {"sat": 101.0, "at": 89.6, "rep": 78.2, "te": 77.6 - _D, "trans": 66.2 - _D, "cis": 59.6 - _D,
       "meth": 45.0 - _D}
Y_BAND_LO = 39.2 - _D               # the test column's lower end
ROW_COLOR = {"sat": SAT, "at": AT, "rep": REP, "te": TE, "trans": HIC, "cis": CIS, "meth": METH}
ORDER = ["sat", "at", "rep", "te", "trans", "cis", "meth"]


def panel_a(ax):
    for x, s in [(X_IN + 6.0, "inputs"), ((X_G0 + X_G1) / 2 + 3, "evidence"), ((X_BAND0 + X_HEAD1) / 2, "test"),
                 (151.0, "calls")]:
        txt(ax, x, Y_TOP - 1.5, s, size=FS_NAME, color=CUE, ha="center")
    txt(ax, (X_FR0 + X_FR1) / 2, Y_LAD + 9.4, "integration", size=FS_NAME, color=CUE, ha="center")
    txt(ax, 0.4, Y_TOP - 1.5, "A", size=FS_LET, weight="bold")
    # the test region: a header holding the rule, a column every line passes through
    ax.add_patch(FancyBboxPatch((X_BAND0, Y_HEAD0), X_HEAD1 - X_BAND0, 15.0, boxstyle="round,pad=0,rounding_size=2.2",
                                fc=BAND, ec="none", zorder=0.5))
    ax.add_patch(FancyBboxPatch((X_BAND0, Y_BAND_LO), X_BAND1 - X_BAND0, Y_HEAD0 - Y_BAND_LO + 3.0,
                                boxstyle="round,pad=0,rounding_size=2.2", fc=BAND, ec="none", zorder=0.5))
    emblem(ax)
    inputs(ax)
    row_satellites(ax)
    row_at(ax)
    row_repeats(ax)
    row_te(ax)
    row_hic(ax)
    row_meth(ax)
    gates(ax)
    ladder(ax)
    calls(ax)


def emblem(ax):
    """the shared test: one mark per chromosome (pass) vs marks in proportion to length (fail);
    each chromosome is judged by what the other chromosomes show"""
    heights = [8.5, 7.4, 6.3, 5.2, 4.1]
    top = Y_HEAD0 + 12.1
    gap = 2.75
    xk = {"pass": X_BAND0 + 6.0, "fail": X_BAND0 + 25.5}
    for kind, x0 in xk.items():
        xs = [x0 + i * gap for i in range(5)]
        for i, (x, h) in enumerate(zip(xs, heights)):
            this = kind == "pass" and i == 4
            bar_v(ax, x, top - h, top, w=1.25, fill="white" if this else CHR_FILL, edge=INK if this else CHR_EDGE,
                  lw_=0.8 if this else 0.55, z=3)
        if kind == "pass":
            for i, (x, h) in enumerate(zip(xs, heights)):
                y = top - h * 0.45
                rect(ax, x - 0.95, y - 0.38, 1.9, 0.76, INK if i < 4 else "white", ec=INK, lw_=0.6, z=5)
            # learned on the other chromosomes, applied to this one
            yb, yt = top + 1.0, top - heights[4] * 0.45
            ax.plot([xs[0] - 0.7, xs[0] - 0.7, xs[3] + 0.7, xs[3] + 0.7], [yb - 0.55, yb, yb, yb - 0.55], color=INK,
                    lw=0.75, zorder=5, solid_joinstyle="miter")
            curve(ax, [((xs[1] + xs[2]) / 2, yb), ((xs[1] + xs[2]) / 2, yb + 1.6), (xs[4] + 2.9, yb + 1.0),
                       (xs[4] + 2.0, yt + 0.4)], INK, lw_=0.75, z=5)
            head(ax, xs[4] + 1.15, yt, INK, s=1.15, direction="left")
            check(ax, x0 - 3.6, top - 4.0, 2.2, lw_=1.8)
            txt(ax, x0 + 2 * gap, Y_HEAD0 + 1.6, "one each", size=FS, color=CUE, ha="center")
        else:
            marks = {0: [0.2, 0.5, 0.8], 1: [0.3, 0.72], 2: [0.6], 3: [], 4: []}
            for i, (x, h) in enumerate(zip(xs, heights)):
                for f in marks[i]:
                    rect(ax, x - 0.95, top - h * f - 0.38, 1.9, 0.76, INK, z=5)
            cross(ax, x0 - 3.6, top - 4.0, 1.9, lw_=1.8)
            txt(ax, x0 + 2 * gap, Y_HEAD0 + 1.6, "by length", size=FS, color=CUE, ha="center")
    # fail: reported only (once, dotted, to the report)
    yf = top - 6.0
    curve(ax, [(xk["fail"] + 4 * gap + 1.6, yf), (164.5, yf)], LINK, lw_=1.0, ls=(0, (1.2, 1.4)), z=2)
    head(ax, 165.6, yf, LINK, s=1.4)


def inputs(ax):
    y_fa = ROW["at"] + 0.2                       # the FASTA feeds three lines: satellites, AT troughs, repeats
    page_icon(ax, X_IN, y_fa - 4.3, 12.4, 8.6, ec=INK, lw_=1.1)
    txt(ax, X_IN + 5.6, y_fa - 0.3, "FASTA", size=FS, weight="bold", ha="center")
    for key, name, col in [("te", "tes/", TE), ("hic", "hic/", HIC), ("meth", "meth/", METH)]:
        y = ROW[key] if key != "hic" else (ROW["trans"] + ROW["cis"]) / 2
        folder_icon(ax, X_IN, y - 3.3, 12.4, 6.0, ec=col)
        txt(ax, X_IN + 6.2, y - 0.25, name, size=FS, weight="bold", ha="center", color=TE_TXT if col == TE else col)
    folder_icon(ax, X_IN, 109.6, 12.4, 6.0, ec="#8A8A8A")
    txt(ax, X_IN + 6.2, 112.35, "genes/", size=FS, ha="center", color=CUE)
    ex = X_G0 - 0.4
    arrow(ax, (X_IN + 12.9, y_fa + 2.4), (ex, ROW["sat"] - 0.1), "#555555", lw_=1.0, rad=-0.2)
    arrow(ax, (X_IN + 12.9, y_fa - 0.2), (ex, ROW["at"] - 0.2), "#555555", lw_=1.0)
    arrow(ax, (X_IN + 12.9, y_fa - 2.8), (ex, ROW["rep"] - 0.6), "#555555", lw_=1.0, rad=0.2)
    p0, p1 = np.array([X_IN + 12.9, 111.4]), np.array([XB0 + 3.6, ROW["sat"] + 3.0])
    u = (p1 - p0) / np.linalg.norm(p1 - p0)
    curve(ax, [tuple(p0), tuple(p1 - 1.2 * u)], "#8A8A8A", lw_=0.9, ls=(0, (2, 1.3)))
    arrow(ax, tuple(p1 - 1.4 * u), tuple(p1), "#8A8A8A", lw_=0.9)
    for y in (ROW["te"], (ROW["trans"] + ROW["cis"]) / 2, ROW["meth"]):
        col = {ROW["te"]: TE, ROW["meth"]: METH}.get(y, HIC)
        arrow(ax, (X_IN + 12.9, y - 0.3), (ex, y - 0.3), col, lw_=1.0)


def row_satellites(ax):
    y = ROW["sat"]
    txt(ax, X_G0, y + 7.9, "satellites", size=FS_NAME, weight="bold", color=SAT)
    bar_h(ax, XB0, X_G1, y, h=1.5)
    for xa, wa, c in [(XB0 + 0.5, 2.2, SUBTEL), (X_G1 - 2.7, 2.2, SUBTEL), (41.5, 2.4, KNOB), (53.4, 4.8, SAT),
                      (61.4, 1.0, SAT), (69.5, 1.6, KNOB)]:
        rect(ax, xa, y - 0.75, wa, 1.5, c, z=4)
    cb = (53.0, 58.6)
    rect(ax, cb[0], y - 1.3, cb[1] - cb[0], 2.6, "none", ec=INK, lw_=1.0, z=6)
    # zoom: the main array as separate tandem copies
    zx0, zx1, zy = 43.0, 70.0, y + 3.3
    n = 11
    L = (zx1 - zx0) / n
    for i in range(n):
        unit_arrow(ax, zx0 + i * L + 0.14, zy, L - 0.34, 1.8, SAT, z=5)
    ax.plot([cb[0], zx0], [y + 1.3, zy - 0.15], color=LINK, lw=0.5, zorder=3)
    ax.plot([cb[1], zx1], [y + 1.3, zy - 0.15], color=LINK, lw=0.5, zorder=3)
    # genes/ veto: a gene-rich repeat is struck out
    gx = 33.8
    rect(ax, gx - 1.3, y - 0.75, 2.6, 1.5, "#9A9A9A", z=4)
    for k in range(3):
        unit_arrow(ax, gx - 2.4 + k * 1.65, y + 1.2, 1.4, 0.95, "#6B6B6B", z=5)
    ax.plot([gx - 2.0, gx + 2.0], [y - 1.5, y + 1.5], color=INK, lw=0.9, zorder=7, solid_capstyle="round")


def gc_trace(x):
    """flat schematic GC curve with one deep trough (no data)"""
    y = 0.42 * np.sin(x * 1.25) + 0.28 * np.sin(x * 2.7 + 1.0)
    y -= 3.0 * np.exp(-0.5 * ((x - 57.0) / 0.95) ** 2)
    y -= 0.7 * np.exp(-0.5 * ((x - 38.5) / 0.9) ** 2)
    return y


def row_at(ax):
    y = ROW["at"]
    txt(ax, X_G0, y + 5.6, "AT troughs", size=FS_NAME, weight="bold", color=AT)
    yb = y - 1.2
    bar_h(ax, XB0, X_G1, yb, h=1.5)
    xs = np.linspace(XB0 + 0.5, X_G1 - 0.5, 800)
    base = y + 2.4
    g = base + gc_trace(xs)
    thr = base + gc_trace(np.array([57.0]))[0] / 2
    m = (g < thr) & (np.abs(xs - 57) < 3)
    ax.fill_between(xs, g, thr, where=m, color=AT, lw=0, zorder=3)
    ax.plot(xs, g, color="#4D4D4D", lw=0.8, zorder=4)
    rect(ax, xs[m].min(), yb - 0.75, xs[m].max() - xs[m].min(), 1.5, AT, z=4)


# repeated sequence per window along the chromosome (schematic, no data): mostly low, one window per chromosome high
REP_UNITS = [1, 0, 1, 2, 1, 0, 1, 1, 2, 0, 1, 1, 0, 1, 2, 2, 3, 5, 3, 2, 1, 1, 0, 1, 2, 1, 0, 1, 1, 0]


def row_repeats(ax):
    """0.16's two FASTA-only repeat lines in one row: the windows where repeated sequence piles up (in small genomes the
    density of repeated 21-mers; in large ones the votes of 31-mer families that cluster at one place per chromosome),
    drawn as a unit histogram, one unit per repeated k-mer or family vote; the peak window marks the chromosome"""
    y = ROW["rep"]
    txt(ax, X_G0, y + 6.1, "repeats", size=FS_NAME, weight="bold", color=REP)
    yb = y - 1.2
    bar_h(ax, XB0, X_G1, yb, h=1.5)
    n = len(REP_UNITS)
    x_lo, x_hi = XB0 + 1.2, X_G1 - 1.2
    pitch = (x_hi - x_lo) / n
    uw, uh, ug = pitch * 0.74, 0.62, 0.24
    ipk = int(np.argmax(REP_UNITS))
    for i, k in enumerate(REP_UNITS):
        xc = x_lo + (i + 0.5) * pitch
        for j in range(k):
            rect(ax, xc - uw / 2, yb + 1.1 + j * (uh + ug), uw, uh, REP if i == ipk else REP_CTX, z=4)
    xpk = x_lo + (ipk + 0.5) * pitch
    rect(ax, xpk - pitch * 0.62, yb - 0.75, pitch * 1.24, 1.5, REP, z=4)


def row_te(ax):
    y = ROW["te"]
    txt(ax, X_G0, y + 6.3, "LTR retrotransposons", size=FS_NAME, weight="bold", color=TE_TXT)
    yb = y - 1.6
    bar_h(ax, XB0, X_G1, yb, h=1.5)
    L, h, dy = 3.2, 1.0, 1.25
    core = [(52.0, 0, False), (55.5, 0, True), (59.0, 0, False), (53.7, 1, True), (57.3, 1, False), (55.4, 2, False)]
    old = [(38.4, 0, True), (41.9, 0, False), (40.1, 1, False)]
    scattered = [(30.0, 0, False), (66.4, 0, True), (72.6, 0, False)]
    for xs, r, rev in core:
        ltr_element(ax, xs, yb + 1.6 + r * dy, L, h, TE_LTR, TE_INT, rev=rev)
    for xs, r, rev in old + scattered:
        ltr_element(ax, xs, yb + 1.6 + r * dy, L, h, OLD_LTR, OLD_INT, rev=rev)
    rect(ax, 51.8, yb - 0.75, 10.6, 1.5, TE, z=4)                       # the core
    rect(ax, 51.8, yb - 1.3, 10.6, 2.6, "none", ec=INK, lw_=1.0, z=6)
    ex0, ey = 65.0, y + 3.7                                              # zoom: one element, two LTRs
    ltr_element(ax, ex0, ey, X_G1 - ex0, 1.9, TE_LTR, TE_INT)
    ax.plot([55.4 + L, ex0], [yb + 1.6 + 2 * dy + h / 2, ey - 0.95], color=LINK, lw=0.5, zorder=3)


def row_hic(ax):
    yt, yc = ROW["trans"], ROW["cis"]
    ymid = (yt + yc) / 2
    txt(ax, X_G0, ymid + 9.0, "Hi-C", size=FS_NAME, weight="bold", color=HIC)
    # the nucleus, Rabl configuration: centromeres cluster at one pole, arms trail behind
    nx, ny, nr = XB0 + 6.2, ymid + 0.4, 5.3
    ax.add_patch(Circle((nx, ny), nr, fc="white", ec="#8C8C8C", lw=0.9, zorder=2))
    apexes = [(nx - 1.4, ny + 3.5), (nx, ny + 4.0), (nx + 1.4, ny + 3.5)]
    arms = [[(-2.6, -6.0), (-1.1, -6.9)], [(-0.55, -7.3), (0.65, -6.7)], [(1.2, -6.4), (2.7, -5.6)]]
    for (cx, cy), aa in zip(apexes, arms):
        for dx, dy in aa:
            t = np.linspace(0, 1, 30)
            ax.plot(cx + dx * t, cy + dy * t, color="#8C8C8C", lw=1.5, solid_capstyle="round", zorder=3)
    for cx, cy in apexes:
        ax.add_patch(Circle((cx, cy), 0.7, fc=HIC, ec="white", lw=0.5, zorder=4))
    # a schematic contact map of two chromosomes
    S, a = 13.0, 7.4                     # map side; chromosome 1 takes a, chromosome 2 the rest
    mx0, mtop = 44.0, ymid + 6.0
    c1, c2 = 3.4, a + 2.4                # centromeres

    def M(u, v):
        return (mx0 + u, mtop - v)
    for (u0, u1) in [(0, a), (a, S)]:
        for (v0, v1) in [(0, a), (a, S)]:
            ax.add_patch(Polygon([M(u0, v0), M(u1, v0), M(u1, v1), M(u0, v1)], closed=True,
                                 fc="#E1EDF7" if u0 == v0 else "#F4F8FC", ec="white", lw=0.9, zorder=2))
    for (u0, u1) in [(0, a), (a, S)]:   # contacts fall with distance: the diagonal
        d = 0.6
        ax.add_patch(Polygon([M(u0, u0), M(u0 + d, u0), M(u1, u1 - d), M(u1, u1), M(u1 - d, u1), M(u0, u0 + d)],
                             closed=True, fc=HIC, ec="none", zorder=3))
    for c, (u0, u1) in [(c1, (0, a)), (c2, (a, S))]:   # Rabl axis: the arms fold at the centromere
        k = min(c - u0, u1 - c)
        ax.plot(*zip(M(c - k, c + k), M(c + k, c - k)), color="#78ADD8", lw=1.2, zorder=3.5, solid_capstyle="butt")
    for (u, v) in [(c1, c2), (c2, c1)]:                # trans: the centromeres meet
        ax.add_patch(Circle(M(u, v), 0.9, fc=HIC, ec="none", zorder=4))
    for c in (c1, c2):                                 # cis: the two flanks of the centromere touch
        x, yy = M(c, c)
        ax.add_patch(Polygon([(x, yy + 0.85), (x + 0.85, yy), (x, yy - 0.85), (x - 0.85, yy)], closed=True, fc=CIS,
                             ec="white", lw=0.4, zorder=5))
    for (u0, u1), c in [((0, a), c1), ((a, S), c2)]:
        bar_h(ax, mx0 + u0 + 0.2, mx0 + u1 - 0.2, mtop + 1.2, h=1.0, z=3)
        rect(ax, mx0 + c - 0.4, mtop + 0.55, 0.8, 1.3, INK, z=4)
    # readouts leave the map
    sx, sy = M(c2, c1)
    curve(ax, [(sx + 0.9, sy), (X_G1 - 4.0, sy), (X_G1 - 3.0, yt), (X_G1, yt)], HIC, lw_=1.1, z=3)
    dx_, dy_ = M(c2, c2)
    curve(ax, [(dx_ + 0.85, dy_), (X_G1 - 4.0, dy_), (X_G1 - 3.0, yc), (X_G1, yc)], CIS, lw_=1.1, z=3)
    txt(ax, X_G1 - 2.0, sy + 1.6, "trans", size=FS, color=HIC, ha="right")
    txt(ax, X_G1 - 2.0, dy_ - 1.6, "cis", size=FS, color=CIS, ha="right")


def row_meth(ax):
    y = ROW["meth"]
    txt(ax, X_G0, y + 7.1, "methylation", size=FS_NAME, weight="bold", color=METH)
    x0 = XB0 + 7.5
    base = y + 0.2
    ax.plot([x0, X_G1], [base, base], color="#8C8C8C", lw=0.6, zorder=2)
    xs = np.linspace(x0 + 1.0, X_G1 - 1.0, 19)
    dip = (xs > 53.0) & (xs < 60.5)
    for x, d in zip(xs, dip):
        top = base + (1.4 if d else 3.5 + 0.25 * np.sin(x))
        ax.plot([x, x], [base, top], color=METH, lw=0.75, zorder=3)
        ax.add_patch(Circle((x, top), 0.58, fc="white" if d else METH, ec=METH, lw=0.75, zorder=4))
    txt(ax, x0 - 0.8, base + 2.0, "CpG", size=FS, color=METH, ha="right")
    # CHG: a second, broader signal (plants)
    xc = np.linspace(x0, X_G1, 400)
    yc = base - 1.75 + 0.14 * np.sin(xc * 1.7) - 1.0 * np.exp(-0.5 * ((xc - 56.6) / 4.5) ** 2)
    ax.plot(xc, yc, color=METH, lw=0.8, zorder=3)
    txt(ax, x0 - 0.8, base - 1.75, "CHG", size=FS, color=METH, ha="right")
    # the call, with the kinetochore (the CpG dip) inside it
    rect(ax, 46.5, base - 4.75, 20.6, 1.6, CONF["high"], ec=INK, lw_=0.6, z=4)
    ax.add_patch(Circle((56.75, base - 3.95), 0.55, fc="white", ec="none", zorder=5))


def gates(ax):
    """each line meets the shared test at a node; lines test one another (links); passing lines feed integration"""
    for k in ORDER:
        x0 = X_G1 + (1.2 if k in ("sat", "at", "rep", "te", "meth") else 0.0)
        curve(ax, [(x0, ROW[k]), (X_BAND0, ROW[k])], ROW_COLOR[k], lw_=1.2, z=5)
        curve(ax, [(X_BAND0, ROW[k]), (X_NODE - 1.3, ROW[k])], ROW_COLOR[k], lw_=1.2, halo="band", z=5)
    # links (background, under the row lines): Hi-C tested against the sequence calls and the LTR cores, and its cis
    # and trans readouts against each other. Nested, so they never cross one another.
    yt, yc, ye = ROW["trans"], ROW["cis"], ROW["te"]
    link(ax, ("sat", 205), ("trans", 150), SAT, (88.0, 96.0), (88.0, yt + 5.3))
    link(ax, ("te", 262), ("trans", 103), TE, (96.6, ye - 4.0), (96.4, yt + 4.2))
    link(ax, ("trans", 215), ("cis", 145), HIC, (94.6, yt - 1.8), (94.6, yc + 1.8), both=True)
    for k in ORDER:
        x, y = X_NODE, ROW[k]
        ax.add_patch(Circle((x, y), 1.35, fc="white", ec=ROW_COLOR[k], lw=1.5, zorder=8))
        ax.add_patch(Circle((x, y), 0.5, fc=ROW_COLOR[k], ec="none", zorder=8.1))


LINK_ALPHA = 0.6


def link(ax, src, dst, color, c1, c2, both=False):
    """a thin background link between two nodes: (row, angle in degrees on the node ring) -> (row, angle);
    c1, c2 are the Bezier control points"""
    def on_ring(row, ang, r=1.75):
        return np.array([X_NODE + r * np.cos(np.radians(ang)), ROW[row] + r * np.sin(np.radians(ang))])
    p0, p3 = on_ring(*src), on_ring(*dst)
    c1, c2 = np.array(c1), np.array(c2)
    hs = 1.0
    d1 = (p3 - c2) / np.linalg.norm(p3 - c2)              # arrival direction
    d0 = (p0 - c1) / np.linalg.norm(p0 - c1)              # departure direction, reversed
    q3 = p3 - d1 * hs * 0.8
    q0 = p0 - d0 * hs * 0.8 if both else p0
    p = PathPatch(MPath([tuple(q0), tuple(c1), tuple(c2), tuple(q3)], [MPath.MOVETO] + [MPath.CURVE4] * 3),
                  fc="none", ec=color, lw=0.6, alpha=LINK_ALPHA, zorder=4.6, capstyle="butt")
    ax.add_patch(p)
    for tip, d in [(p3, d1)] + ([(p0, d0)] if both else []):
        n = np.array([-d[1], d[0]])
        base = tip - d * hs
        ax.add_patch(Polygon([tuple(tip), tuple(base + n * hs * 0.42), tuple(base - n * hs * 0.42)], closed=True,
                             fc=color, ec="none", alpha=LINK_ALPHA, zorder=4.6))


def ladder(ax):
    hgt, gap = 3.6, 0.5
    ys = {"high": Y_LAD + hgt + gap, "medium": Y_LAD, "low": Y_LAD - hgt - gap}
    for k, y in ys.items():
        rect(ax, X_LAD0, y - hgt / 2, X_LAD1 - X_LAD0, hgt, CONF[k], ec=INK if k == "low" else "none", lw_=0.5, z=5)
        txt(ax, (X_LAD0 + X_LAD1) / 2, y - 0.1, k, size=FS, ha="center", color="white" if k != "low" else INK)
    # the frame: every passing line enters one integration step, which grades the call
    ax.add_patch(FancyBboxPatch((X_FR0, Y_LAD - 7.5), X_FR1 - X_FR0, 15.0, boxstyle="round,pad=0,rounding_size=1.4",
                                fc="white", ec="#8C8C8C", lw=0.9, zorder=4.5))
    ends = np.linspace(Y_LAD + 2.1, Y_LAD - 2.1, len(ORDER))
    for k, ye in zip(ORDER, ends):
        x0, y0 = X_NODE + 1.35, ROW[k]
        curve(ax, [(x0, y0), (x0 + 5.0, y0), (X_FR0 - 4.4, ye), (X_FR0, ye)], ROW_COLOR[k], lw_=1.1, z=4)
    arrow(ax, (X_FR1 + 0.3, Y_LAD), (X_FR1 + 4.0, Y_LAD), INK, lw_=1.3)


def chromosome_v(ax, x, top, L, cen, conf, w=1.7, feats=(), kin=False):
    bar_v(ax, x, top - L, top, w=w, z=3)
    yc = top - L * cen
    call_box(ax, x, yc, conf, w=3.1, h=2.4)
    s = 1.0
    ax.add_patch(Polygon([(x, top + 1.3), (x - s * 0.5, top + 0.3), (x + s * 0.5, top + 0.3)], fc=INK, ec="none",
                         zorder=5))
    ax.add_patch(Polygon([(x, top - L - 1.3), (x - s * 0.5, top - L - 0.3), (x + s * 0.5, top - L - 0.3)], fc=INK,
                         ec="none", zorder=5))
    for kind, f in feats:
        yy = top - L * f
        if kind == "knob":
            rect(ax, x - w / 2, yy - 0.75, w, 1.5, KNOB, z=4)
        else:
            ax.add_patch(Polygon([(x + 1.6, yy + 0.85), (x + 2.45, yy), (x + 1.6, yy - 0.85), (x + 0.75, yy)],
                                 fc=RDNA_MARK if kind == "45S" else "white", ec="none" if kind == "45S" else "#6B6B6B",
                                 lw=0.75, zorder=5))
    if kin:
        ax.add_patch(Circle((x, yc), 0.6, fc="white", ec="none", zorder=7))
    return yc


KARYO_SCALE = 1.40          # 0.16: the karyotype fills the taller panel A (chromosome lengths x this)


def calls(ax):
    top, x0, dx = 103.0, 130.4, 4.7
    spec = [(40.0, 0.42, "high", [("45S", 0.12)], True), (36.0, 0.47, "high", [("knob", 0.8)], False),
            (33.0, 0.31, "medium", [], False), (30.5, 0.5, "high", [("5S", 0.74)], False),
            (28.0, 0.16, "medium", [], False), (25.0, 0.4, "high", [("knob", 0.15)], False),
            (21.0, 0.55, "low", [], False)]
    ycs = [chromosome_v(ax, x0 + i * dx, top, L * KARYO_SCALE, cen, conf, feats=f, kin=k)
           for i, (L, cen, conf, f, k) in enumerate(spec)]
    # the call in words (basis, caveat)
    bubble(ax, x0 + 6 * dx + 4.2, ycs[6] - 2.4, 10.0, 5.0, (x0 + 6 * dx + 1.75, ycs[6]))
    # outputs: report and overview
    rx, ry = 166.0, 108.4
    page_icon(ax, rx, ry, 10.0, 12.6, lw_=0.9)
    for k in range(5):
        ax.plot([rx + 1.5, rx + 7.4 - (2.4 if k == 4 else 0)], [ry + 9.2 - k * 1.9] * 2, color="#7A7A7A", lw=0.8,
                zorder=6)
    txt(ax, rx + 5.0, ry - 2.0, "report", size=FS, ha="center")
    oy = 56.0
    page_icon(ax, rx, oy, 10.0, 12.6, lw_=0.9)
    for k in range(5):
        yy = oy + 9.4 - k * 1.95
        ax.plot([rx + 1.3, rx + 8.3 - k * 1.0], [yy, yy], color=OTHER, lw=1.4, zorder=6, solid_capstyle="butt")
        ax.plot([rx + 1.3, rx + 8.3 - k * 1.0], [yy + 0.6, yy + 0.6], color=HIC, lw=0.5, zorder=6)
        rect(ax, rx + 2.4 + 0.45 * k, yy - 0.5, 1.0, 1.0, "white", ec=INK, lw_=0.5, z=7)
    txt(ax, rx + 5.0, oy - 2.0, "overview", size=FS, ha="center")
    # direct key for the other features
    k1, k2, kx = 40.2, 36.0, 128.4
    ax.add_patch(Polygon([(kx + 0.6, k1 + 0.75), (kx + 0.05, k1 - 0.35), (kx + 1.15, k1 - 0.35)], fc=INK, ec="none",
                         zorder=5))
    txt(ax, kx + 2.0, k1, "telomere", size=FS)
    dx2 = kx + 16.0
    ax.add_patch(Polygon([(dx2, k1 + 0.85), (dx2 + 0.85, k1), (dx2, k1 - 0.85), (dx2 - 0.85, k1)], fc=RDNA_MARK,
                         ec="none", zorder=5))
    ax.add_patch(Polygon([(dx2 + 2.4, k1 + 0.85), (dx2 + 3.25, k1), (dx2 + 2.4, k1 - 0.85), (dx2 + 1.55, k1)],
                         fc="white", ec="#6B6B6B", lw=0.75, zorder=5))
    txt(ax, dx2 + 4.2, k1, "45S, 5S rDNA", size=FS)
    rect(ax, kx - 0.2, k2 - 0.75, 1.7, 1.5, KNOB, z=5)
    txt(ax, kx + 2.0, k2, "knob", size=FS)
    call_box(ax, dx2 + 0.0, k2, "high", w=3.1, h=2.2)
    ax.add_patch(Circle((dx2, k2), 0.6, fc="white", ec="none", zorder=7))
    txt(ax, dx2 + 2.4, k2, "kinetochore", size=FS)


def bubble(ax, x, y, w, h, tip, small=False):
    """the call in words (basis, caveat): a speech bubble whose tail points at the call"""
    lw_ = 0.7 if small else 0.8
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle=f"round,pad=0,rounding_size={0.8 if small else 1.1}",
                                fc="white", ec=INK, lw=lw_, zorder=8))
    a, b = (x + 0.0, y + h * (0.56 if small else 0.66)), (x + 0.0, y + h * (0.22 if small else 0.30))   # the tail
    ax.add_patch(Polygon([(a[0] + 0.45, a[1]), tip, (b[0] + 0.45, b[1])], fc="white", ec="none", zorder=8.2))
    ax.plot([a[0], tip[0], b[0]], [a[1], tip[1], b[1]], color=INK, lw=lw_, zorder=8.3, solid_joinstyle="miter")
    m, dy, short = (1.0, 1.2, 1.6) if small else (1.4, 1.6, 2.2)
    for k in range(2):
        ax.plot([x + m, x + w - m - short * k], [y + h * (0.68 if small else 0.64) - k * dy] * 2, color="#7A7A7A",
                lw=0.75 if small else 0.85, zorder=9)


# ---------------------------------------------------------------- panel B: the test on real data
def read_panel_b(path=PANEL_B):
    """panel B's data (workflow_panelB.tsv): per genome, the chromosomes in drawing order with their lengths, every
    cluster of the centromeric satellite families (bp, chromosome, start, end; 0-based half-open), the N largest of
    them (the ones the test counts), the chance distribution (chromosomes covered: random draws) and the verdict"""
    t = pd.read_csv(path, sep="\t", comment="#", dtype=str)
    data = {}
    for g, d in t.groupby("genome", sort=False):
        ch = d[d.kind == "chromosome"]
        cl = d[d.kind.isin(["cluster", "largest"])]
        test = dict(zip(d[d.kind == "test"].name, d[d.kind == "test"].value))
        nu = d[d.kind == "chance"]
        data[g] = dict(chroms=dict(zip(ch.name, ch.end.astype(int))),
                       clusters=[(int(v), c, int(s), int(e)) for c, s, e, v in zip(cl.name, cl.start, cl.end, cl.value)],
                       top=[(int(v), c, int(s), int(e)) for c, s, e, v, k in zip(cl.name, cl.start, cl.end, cl.value,
                                                                                  cl.kind) if k == "largest"],
                       chance=dict(zip(nu.start.astype(int), nu.value.astype(int))),
                       test=dict(n=int(test["n"]), covered=int(test["covered"]),
                                 null_median=float(test["chance_median"]), p=float(test["p"]), status=test["status"]),
                       cluster_gap=int(test["merge_gap_bp"]))
        data[g]["covered"] = {c for _, c, _, _ in data[g]["top"]}
        assert len(data[g]["covered"]) == data[g]["test"]["covered"] and len(data[g]["top"]) == data[g]["test"]["n"]
    return data


SCALE_B = 0.30         # mm per Mb, both genomes (0.15: 0.27; 0.16 gives B a taller band)
MINH_B = 0.45          # smallest drawn cluster height (mm)


def panel_b(ax, data):
    txt(ax, 0.4, 36.4 + BC_DY, "B", size=FS_LET, weight="bold")
    top = 29.0 + BC_DY
    xs0 = {"rice": 6.0, "tomato": 43.0}
    dx = 2.9
    for g, sp in (("rice", "Oryza sativa"), ("tomato", "Solanum lycopersicum")):
        d = data[g]
        x0 = xs0[g]
        txt(ax, x0 - 1.0, 35.2 + BC_DY, sp, size=FS, style="italic")
        txt(ax, x0 - 1.0 + text_mm(sp + " ", FS, style="italic"), 35.2 + BC_DY, f"({g})", size=FS, color=CUE)
        topset = {(c, s, e) for _, c, s, e in d["top"]}
        for i, (c, L) in enumerate(d["chroms"].items()):
            x = x0 + i * dx
            bar_v(ax, x, top - L / 1e6 * SCALE_B, top, w=1.35, z=2, lw_=0.5)
            mine = [t for t in d["clusters"] if t[1] == c]
            for accent in (False, True):
                for bp, cc, s, e in mine:
                    if ((cc, s, e) in topset) != accent:
                        continue
                    y1, y2 = top - s / 1e6 * SCALE_B, top - e / 1e6 * SCALE_B
                    hh = max(y1 - y2, MINH_B)
                    rect(ax, x - 1.05, (y1 + y2) / 2 - hh / 2, 2.1, hh, SAT if accent else "#A6A6A6",
                         z=4 if accent else 3)
            cov = c in d["covered"]
            ax.add_patch(Circle((x, top + 2.0), 0.65, fc=SAT if cov else "white", ec=SAT if cov else "#8C8C8C",
                                lw=0.75, zorder=5))
    sb = 10 * SCALE_B                               # scale bar under the (shorter) rice karyotype
    xb = xs0["rice"] - 0.7
    ax.plot([xb, xb + sb], [9.0, 9.0], color=INK, lw=1.0, solid_capstyle="butt")
    txt(ax, xb + sb + 1.0, 9.0, "10 Mb", size=FS)
    # the count against chance: chromosomes covered by the 12 largest clusters
    ax0, ay0, ay1 = 86.0, 5.0, 29.0 + BC_DY
    kmin, kmax = 4, 12

    def Y(k):
        return ay0 + (k - kmin) / (kmax - kmin) * (ay1 - ay0)
    ax.plot([ax0, ax0], [Y(kmin) - 0.8, Y(kmax) + 0.8], color=INK, lw=0.6)
    for k in (4, 8, 12):
        ax.plot([ax0 - 0.8, ax0], [Y(k), Y(k)], color=INK, lw=0.6)
        txt(ax, ax0 - 1.3, Y(k), str(k), size=FS, ha="right")
    ax.text(ax0 - 5.6, (ay0 + ay1) / 2, "chromosomes covered", fontsize=FS, rotation=90, ha="center", va="center")
    for j, g in enumerate(("rice", "tomato")):
        d = data[g]
        cx = ax0 + 4.6 + j * 8.4
        tot = sum(d["chance"].values())
        for v, c in d["chance"].items():
            if v >= kmin:
                hw = 3.4 * c / tot / 0.45
                rect(ax, cx - hw, Y(v) - 1.1, 2 * hw, 2.2, CHANCE, z=3)
        obs = d["test"]["covered"]
        ax.add_patch(Circle((cx, Y(obs)), 1.05, fc=SAT, ec="white", lw=0.6, zorder=6))
        (check if d["test"]["status"] == "pass" else cross)(ax, cx + 3.0, Y(obs), 1.9 if g == "rice" else 1.6,
                                                             lw_=1.5)
        txt(ax, cx, 2.0, g, size=FS, ha="center", color=CUE)
        if g == "tomato":                          # the chance distribution, named where it is widest
            tot = sum(d["chance"].values())
            edge = max(3.4 * d["chance"].get(k, 0) / tot / 0.45 for k in (7, 8))
            txt(ax, cx + edge + 1.2, Y(7.2), "chance", size=FS, color="#8C8C8C")


# ---------------------------------------------------------------- panel C: the integration rules
def panel_c(ax):
    X0 = 115.0
    txt(ax, X0, 36.4 + BC_DY, "C", size=FS_LET, weight="bold")
    b0, b1, a0, a1 = 128.5, 150.5, 157.0, 179.0
    ys = np.linspace(29.8 + BC_DY, 4.6, 5)
    rows = list(zip(["agree", "oppose", "move", "place", "flag"], ys))
    for name, y in rows:
        txt(ax, X0 + 1.0, y, name, size=FS)
        bar_h(ax, b0, b1, y, h=1.4)
        bar_h(ax, a0, a1, y, h=1.4)
        arrow(ax, (b1 + 1.4, y), (a0 - 1.4, y), "#555555", lw_=1.1, head=4.2)
    s = 1.8
    y = rows[0][1]                       # agree: each tested line at the call raises it one level (two at most)
    call_box(ax, 139.5, y, "low")
    tri_down(ax, 138.4, y + 1.7, s, HIC)
    tri_down(ax, 140.6, y + 1.7, s, METH)
    call_box(ax, 168.0, y, "high")
    y = rows[1][1]                       # oppose: a tested line elsewhere keeps the call, one level lower
    call_box(ax, 136.0, y, "high")
    tri_down(ax, 146.0, y + 1.7, s, HIC)
    call_box(ax, 164.5, y, "medium")
    y = rows[2][1]                       # move: a high call moves only with two lines at the new place
    call_box(ax, 134.0, y, "high")
    tri_down(ax, 144.9, y + 1.7, s, HIC)
    tri_down(ax, 147.1, y + 1.7, s, TE)
    arrow(ax, (134.8, y - 1.3), (145.6, y - 1.25), INK, lw_=0.9, rad=0.32, head=3.4)
    call_box(ax, 162.5, y, None, ls=(0, (1.6, 1.0)), ec="#6B6B6B")
    call_box(ax, 173.5, y, "high")
    y = rows[3][1]                       # place: no sequence call; a tested line places one
    tri_down(ax, 142.5, y + 1.7, s, HIC)
    call_box(ax, 171.0, y, "medium")
    y = rows[4][1]                       # flag (0.16): a tested repeat line elsewhere neither moves nor lowers a low call;
    call_box(ax, 136.0, y, "low")        # it is said in the call's words (the caveat)
    tri_down(ax, 146.0, y + 1.7, s, REP)
    call_box(ax, 164.5, y, "low")
    bubble(ax, 168.4, y + 0.9, 7.6, 3.6, (165.6, y + 1.0), small=True)


# ---------------------------------------------------------------- legend page and numbers
def numbers(data):
    """every number on the figure and in the legend, with its source (written only with --numbers)"""
    r, t = data["rice"], data["tomato"]
    tsv = "docs/img/src/workflow_panelB.tsv"
    me = "docs/img/src/workflow.py"
    tot = {g: sum(d["chance"].values()) for g, d in data.items()}
    rows = [
        ("rice_chromosomes", r["test"]["n"], f"{tsv} test n; B: 'the 12 largest'"),
        ("rice_covered", r["test"]["covered"], f"{tsv} test covered (= chromosomes of the 'largest' rows)"),
        ("rice_chance_median", f"{r['test']['null_median']:g}", f"{tsv} test chance_median"),
        ("rice_P", f"{r['test']['p']:.4f}", f"{tsv} test p (10,000 draws: the smallest possible); B check = pass"),
        ("rice_clusters_drawn", len(r["clusters"]), f"{tsv} cluster + largest rows"),
        ("tomato_chromosomes", t["test"]["n"], f"{tsv} test n"),
        ("tomato_covered", t["test"]["covered"], f"{tsv} test covered"),
        ("tomato_chance_median", f"{t['test']['null_median']:g}", f"{tsv} test chance_median"),
        ("tomato_P", f"{t['test']['p']:.2f}", f"{tsv} test p; B cross = fail"),
        ("tomato_clusters_drawn", len(t["clusters"]), f"{tsv} cluster + largest rows"),
        ("cluster_merge_gap_kb", r["cluster_gap"] // 1000, f"{tsv} test merge_gap_bp"),
        ("random_draws", tot["rice"], f"{tsv} chance rows (sum per genome; tomato {tot['tomato']})"),
        ("drawn_min_cluster_Mb", round(MINH_B / SCALE_B, 1), f"{me} MINH_B / SCALE_B"),
        ("B_scale_bar_Mb", 10, f"{me} (SCALE_B mm per Mb)"),
        ("B_axis_ticks", "4, 8, 12", f"{me} (chromosomes covered axis)"),
        ("max_raise_levels", 2, "cenfuse/params.py MAX_RAISE; C: 'two at most'"),
        ("lines_of_evidence", 6, "panel A rows: 3 from the assembly (satellites, AT troughs, repeats), 3 from optional "
                                 "inputs (LTR retrotransposons, Hi-C read two ways = 2 test circles, methylation)"),
    ]
    return pd.DataFrame(rows, columns=["number", "value", "source"])


def legend_page_text(data):
    """the one legend page (making-figures skill): headline and meaning, then a few sentences per panel; no methods,
    thresholds or statistics (those are in docs/design.md); worded so that 0.16.1's leave-one-out statuses and a
    recalibrated repeat-family test leave it true"""
    r, t = data["rice"]["test"], data["tomato"]["test"]
    minmb = MINH_B / SCALE_B
    paras = [
        ("", "A centromere occurs once per chromosome, so whatever marks centromeres falls once on each chromosome, "
             "while other repeats fall roughly in proportion to chromosome length. cenfuse lets each kind of evidence act only "
             "after it shows this pattern in the genome at hand, and grades every call."),
        ("A", "Six lines of evidence find centromere-like places on their own: three from the assembly alone "
              "(satellite arrays, AT-rich troughs, repeats) and three from optional inputs (dashed: retrotransposons, "
              "Hi-C read two ways, methylation). Each is tested in this genome, mostly by whether it marks one place "
              "per chromosome rather than places spread by length (gray column; each chromosome weighed against the "
              "others); a line that fails changes no call. Lines that pass, some tested against one another (thin "
              "links), give one call per chromosome, graded high, medium or low and explained in words (bubble)."),
        ("B", f"The test on real satellites: every cluster of the centromeric satellite families, to scale (drawn at "
              f"least {minmb:.1f} Mb tall); orange, the {r['n']} largest; filled dots, the chromosomes they reach. In "
              f"rice they reach all {r['covered']} chromosomes, far more than the {r['null_median']:.0f} expected by "
              f"chance; in tomato they pile onto some chromosomes and miss others ({t['covered']} of {t['n']}), as "
              f"chance would, so its satellite calls stay low."),
        ("C", "Tested lines adjust the call from the sequence (left: before; right: after; triangles: where a tested "
              "line points, in its color). Each agreeing line raises the call one level, two at most. Hi-C, an LTR "
              "core or a CpG dip pointing elsewhere lowers it one level; strong evidence elsewhere moves it, a high "
              "call only with two lines (dashed: the old place, still reported). Where the sequence gives no call, "
              "Hi-C or an LTR core places one at medium, troughs or repeats alone at low, and a trough or repeat peak "
              "at a Hi-C call lifts it one level. A repeat line pointing away from a low call leaves it in place and "
              "flags it in words (bubble)."),
    ]
    stamp = (f"Data: docs/img/src/workflow_panelB.tsv (cenfuse 0.16.0 benchmark runs) · Script: "
             f"docs/img/src/workflow.py · {FIG_DATE}")
    return paras, stamp


def rich_lines(text, width_mm, size=10.0):
    """greedy word wrap with *italic* spans; returns lines of (word, style) tokens"""
    toks = []
    for i, part in enumerate(re.split(r"(\*[^*]+\*)", text)):
        if not part:
            continue
        style = "italic" if part.startswith("*") and part.endswith("*") else "normal"
        for w in part.strip("*").split(" "):
            if not w:
                continue
            if "/" in w and text_mm(w, size) > width_mm * 0.5:     # a long path: break after any slash
                bits = re.findall(r"[^/]*/|[^/]+$", w)
                toks.append((bits[0], style, False))
                toks += [(b, style, True) for b in bits[1:]]
            else:
                toks.append((w, style, False))
    space = text_mm(" ", size)
    lines, cur, wcur = [], [], 0.0
    for w, st, glue in toks:
        ww = text_mm(w, size, style=st)
        gap = 0.0 if (glue or not cur) else space
        if cur and wcur + gap + ww > width_mm:
            lines.append(cur)
            cur, wcur, gap = [], 0.0, 0.0
        cur.append((w, st, ww, gap))
        wcur += gap + ww
    if cur:
        lines.append(cur)
    return lines


def legend_page(pdf, paras, stamp=None, size=10.0):
    left, right, indent, top_m, bot_m = 10.0, 10.0, 6.0, 10.0, 9.0
    lh = size * 1.32 * 25.4 / 72            # line height, mm
    blocks = []
    for letter, body in paras:
        x0 = left + (indent if letter else 0)
        blocks.append((letter, x0, rich_lines(body, W - x0 - right, size)))
    if stamp:
        blocks.append(("", left, rich_lines(stamp, W - left - right, size)))
    n = sum(len(b[2]) for b in blocks)
    hh = top_m + bot_m + n * lh + (len(blocks) - 1) * lh * 0.55
    fig = plt.figure(figsize=(W / 25.4, hh / 25.4))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, W)
    ax.set_ylim(0, hh)
    ax.axis("off")
    y = hh - top_m
    for letter, x0, lines in blocks:
        if letter:
            ax.text(left, y, letter, fontsize=size, fontweight="bold", va="baseline", ha="left")
        for ln in lines:
            x = x0
            for w, st, ww, gap in ln:
                x += gap
                ax.text(x, y, w, fontsize=size, fontstyle=st, va="baseline", ha="left")
                x += ww
            y -= lh
        y -= lh * 0.55
    pdf.savefig(fig, facecolor="white")
    plt.close(fig)


# ---------------------------------------------------------------- main
def build(verbose=False):
    data = read_panel_b()
    if verbose:
        for g, d in data.items():
            print(f"{g}: {d['test']} clusters={len(d['clusters'])}", file=sys.stderr)
    fig = plt.figure(figsize=(W / 25.4, H / 25.4))
    fig.patch.set_facecolor("white")
    axes = []
    for shift in (A_SHIFT, 0.0):                  # panel A's frame, then the page frame for B and C (both in mm)
        ax = fig.add_axes([0, 0, 1, 1])
        ax.set_xlim(0, W)
        ax.set_ylim(-shift, H - shift)
        ax.set_aspect("equal")
        ax.axis("off")
        axes.append(ax)
    panel_a(axes[0])
    panel_b(axes[1], data)
    panel_c(axes[1])
    return fig, data


def main():
    ap = argparse.ArgumentParser(description="Draw cenfuse's workflow figure: PNG (figure only) and PDF (figure, then "
                                             "its legend page), from this folder's workflow_panelB.tsv.")
    ap.add_argument("--png", default=str(IMG / "workflow.png"), help="PNG to write [docs/img/workflow.png]")
    ap.add_argument("--pdf", default=str(IMG / "workflow.pdf"), help="PDF to write [docs/img/workflow.pdf]")
    ap.add_argument("--numbers", help="also write every number on the figure and in the legend, with its source")
    ap.add_argument("--dpi", type=int, default=400, help="PNG resolution [400]")
    ap.add_argument("-v", "--verbose", action="store_true", help="print panel B's test numbers")
    a = ap.parse_args()
    fig, data = build(a.verbose)
    fig.savefig(a.png, dpi=a.dpi, facecolor="white")
    from PIL import Image                          # opaque RGB: GitHub shows it on dark backgrounds too
    Image.open(a.png).convert("RGB").save(a.png, dpi=(a.dpi, a.dpi))
    paras, stamp = legend_page_text(data)
    with PdfPages(a.pdf, metadata={"CreationDate": None}) as pdf:     # no creation time: rebuilds are identical
        pdf.savefig(fig, facecolor="white")
        legend_page(pdf, paras, stamp)
    plt.close(fig)
    if a.numbers:
        numbers(data).to_csv(a.numbers, sep="\t", index=False)
    print(f"wrote {a.png} and {a.pdf}" + (f" and {a.numbers}" if a.numbers else ""))


if __name__ == "__main__":
    main()
