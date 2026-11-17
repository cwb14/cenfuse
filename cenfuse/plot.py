"""Overview figure, one per genome (cenfuse_overview.pdf: figure pages, then a short and a long legend per genome, each
on pages of at most closeup.PAGE_H; <genome>.overview.png and <genome>.closeups.png: the first genome page and the first
close-up page at 300 dpi).

Genome pages, ROWS_PER_PAGE chromosomes a page:
A, every chromosome to one bp scale, one row each: the satellite arrays painted by their family's role, the centromeric
   LTR family's elements as ticks, the Hi-C trace (the contact with the other centromeres, else the hub) as a thin line
   on one genome-wide scale, the centromere as one triangle in the confidence style (high solid black, medium gray,
   low open), an x in a signal's color at the place of each signal that points elsewhere (not one that failed its
   genome test or is a lead only), gray boxes at the close-up windows (closeup.window), telomere and rDNA arrays.
B, beside A, row by row: each signal's role in each call (signals.rows: used, points elsewhere, not used), one column
   per close-up track (_keys: signals.present, in signals.TRACK_ORDER), grayed where the signal failed its genome test;
   then the call's confidence.
C: each signal's genome test (signals.signal_tests, in TRACK_ORDER) in plain words, with P where the summary gives it,
   k of n chromosomes as a bar with chance as a tick, and the verdict.
D: each signal with a trace (Hi-C, CHG and CpG methylation, the LTR family's density, the centromeric satellite's share
   of the DNA) within half a close-up width of every centromere (faint) and their mean (bold), against N_NULL random
   places on each of the same chromosomes, windows kept inside the chromosome (the mean over chromosomes at each draw:
   median and 5-95% band); both on the chromosomes at least a close-up width long. A window that runs past a chromosome
   end has no value. Over each plot, how many of its centromeres that signal placed (signals.Row.placed): a signal
   rises at the centromeres it placed by construction.
C and D sit under the last rows of A, or on a page of their own where that page would pass closeup.PAGE_H.
With scored plants (docs/plants.md), a plants page after the close-ups (E):
F: one row per plant, one column per chromosome: does the plant's centromere sit at the reference call?
G: the chromosome where the plants differ most, plant by plant: each plant's methylation track under the reference
   plant's evidence, at one bp scale.
Legends (legends): the short one gives the result, then each panel's marks in plain words; the long one says how cenfuse
places a centromere (the roles, the confidence, each kind of verdict with its `rule` values) and what this genome holds.
Neither uses a word of signals.BANNED; `rule` values appear in the long one only, in parentheses."""
from __future__ import annotations

import logging
import math
import re
import warnings
from collections import Counter
from collections.abc import Iterator
from contextlib import contextmanager
from dataclasses import dataclass
from datetime import date
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages  # noqa: E402
from matplotlib.collections import LineCollection, PolyCollection  # noqa: E402
from matplotlib.figure import Figure  # noqa: E402
from matplotlib.font_manager import FontProperties  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch, Rectangle  # noqa: E402
from matplotlib.textpath import TextPath  # noqa: E402

from . import closeup as cu  # noqa: E402  (closeup reads this module's attributes at call time only)
from . import signals as sg  # noqa: E402
from .runtime import atomic_path  # noqa: E402

log = logging.getLogger("cenfuse")
MM = 1 / 25.4                                    # layout: inches per millimeter
WIDTH = 180 * MM
ROWS_PER_PAGE = 24                               # layout: chromosome rows per page
ROLE_COLOR = {"centromeric": "#D55E00", "cen_associated": "#F4A582", "interstitial": "#CC79A7",
              "subtelomeric": "#56B4E9", "rDNA": "#009E73"}
OTHER = "#B0B0B0"
RDNA_MARK = "#8C8C8C"                            # rDNA arrays: context, so they do not outshine the centromere calls
HIC, METH, TE = "#0072B2", "#303030", "#E69F00"
CALL_BAND = "#E8E8E8"                            # the plants page: the reference call
VERDICT_COLOR = {"shared": "#BDBDBD", "shifted": "#D55E00", "possible_shift": "#E69F00", "call_absent": "#000000",
                 "conflict": "#CC79A7", "unclear": "#F2F2F2", "no_data": "white", "no_call": "white"}
VERDICT_LABEL = {"shared": "consistent with the centromere", "shifted": "shifted", "possible_shift": "possible shift",
                 "call_absent": "centromere DNA absent", "conflict": "signals disagree", "unclear": "unclear",
                 "no_data": "no data"}
PROV_SHORT = {"reference": "reference plant", "other": "other plant", "unknown": "origin unknown"}
plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Liberation Sans", "Arial", "DejaVu Sans"],
                     "font.size": 6.5, "pdf.fonttype": 42, "ps.fonttype": 42, "axes.linewidth": 0.5,
                     "xtick.major.width": 0.5, "ytick.major.width": 0.5, "xtick.major.size": 2,
                     "ytick.major.size": 2})
# genome pages
N_NULL = 200                                     # prior: random places per chromosome in panel D, ten draws in each
                                                 # tail of the 5-95% band
ROW_HIC = 6.0                                    # layout: a chromosome row with a Hi-C trace, mm
ROW_BASE = 4.2                                   # layout: a chromosome row without one, mm
TRI_Y = 0.75                                     # layout: the centromere triangle's center above the row's foot, mm
BAR_Y = 1.5                                      # layout: the chromosome bar's lower edge above the row's foot, mm
BAR_H = 1.3                                      # layout: the chromosome bar's height, mm
TICK_H = 0.8                                     # layout: an LTR element's tick above the bar, mm
HIC_H = 2.6                                      # layout: the Hi-C trace's height above the bar (0 to the top), mm
MIN_FEAT = 0.25                                  # layout: the narrowest array drawn, mm (smaller ones widened)
COL_MM = 4.2                                     # layout: panel B's column spacing, mm (a rotated 6-pt head: 2.2 mm)
NAME_GAP = 2.0                                   # layout: chromosome names to the chromosomes' start, mm
AB_GAP = 5.0                                     # layout: panel A's right end to panel B, mm
HEAD_TOP = 2.0                                   # layout: page top to panel B's column heads, mm
AXIS_H = 8.0                                     # layout: panel A's axis under its rows (ticks, labels, title), mm
BLOCK_GAP = 5.0                                  # layout: space above panels C and D, holding their letters, mm
TEST_ROW = 4.2                                   # layout: one row of panel C, mm
TEST_BAR = 18.0                                  # layout: panel C's k-of-n bar, mm
D_H = 22.0                                       # layout: a panel D plot's height, mm
D_TOP = 3.0                                      # layout: above a panel D plot, holding n, mm
D_FOOT = 8.5                                     # layout: under a panel D plot, its two-line axis title, mm
D_GUT = 14.0                                     # layout: left of each panel D plot, its two-line axis title, mm
D_MAX_W = 36.0                                   # layout: the widest panel D plot, mm
KEY_GAP = 2.0                                    # layout: the last panel to the key, mm
WIN_BAND = "#EFEFEF"                             # a close-up's window
UNTRUSTED = ("not trusted here", "not used (other plants)", "a lead only")
CEN_ROLES = ("centromeric", "cen_associated")    # the satellite arrays panel D measures
B_HEAD = {"hic": "Hi-C", "hub": "Hi-C hub", "flank": "Hi-C flanks", "fold": "Hi-C arm fold",
          "align": "Hi-C arms align", "meth": "CHG", "cpg": "CpG", "te": "{clade}", "sat": "satellite",
          "at": "AT-rich", "rep": "repeats", "fam": "repeat families"}
D_LABEL = {"hic": "Hi-C contact\n(× chromosome median)", "hub": "Hi-C contact\n(× chromosome median)",
           "meth": "CHG methylation\nchange", "cpg": "CpG methylation\n(% methylated)",
           "te": "{clade} elements\nper Mb", "sat": "centromeric satellite\n(% of DNA)"}
D_COLOR = {"hic": HIC, "hub": HIC, "meth": METH, "cpg": METH, "te": TE, "sat": ROLE_COLOR["centromeric"]}
# per panel D signal, the placing signals (signals.Row.placed) that put a centromere where it peaks by construction;
# methylation places none
D_PLACED = {"hic": sg.HIC_KEYS, "hub": sg.HIC_KEYS, "meth": frozenset(), "cpg": frozenset(), "te": frozenset({"te"}),
            "sat": frozenset({"sat"})}


def _iv(a, b):
    return (int(a) - 1, int(b)) if a != "." and b != "." else None


def _rows(calls):
    """Each chromosome's length and call (0-based half-open), for the plants page."""
    return [dict(chrom=r["chrom"], L=int(r["chrom_len"]), call=_iv(r["start"], r["end"]))
            for r in calls.to_dict("records")]


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


def _text_mm(s: str, size: float) -> float:
    return TextPath((0, 0), s, size=size, prop=FontProperties(size=size)).get_extents().width / 72 * 25.4


def _nice_ticks(Lmax):
    for step in (1e4, 2e4, 5e4, 1e5, 2e5, 5e5, 1e6, 2e6, 5e6, 10e6, 20e6, 25e6, 50e6, 100e6, 200e6):
        if Lmax / step <= 8:
            return np.arange(0, Lmax + 1, step)
    return np.arange(0, Lmax + 1, 500e6)


# ------------------------------------------------------------------------------------------------ genome pages: data
@dataclass
class _Genome:
    """What every genome page of one genome shares: its rows, signals, tests, scales and widths (mm)."""
    f: dict
    rows: list
    summary: dict | None
    clade: str
    keys: list                           # the close-up tracks (_keys): panel B's columns, the windows' signals
    tests: list                          # panel C's rows (signals.signal_tests, in TRACK_ORDER)
    verdict: dict                        # signal -> its genome test's verdict
    W: float                             # the close-up width, bp; panel D spans it
    wins: dict                           # chromosome -> its close-up window (lo, hi)
    trace: str | None                    # the Hi-C trace: "hic", else "hub", else None
    hsc: float                           # the Hi-C trace's top (x chromosome median)
    te: bool                             # LTR elements drawn
    pkeys: list                          # panel D's signals
    Lmax: int
    row_mm: float
    a_x0: float
    a_w: float
    b_x0: float
    b_w: float
    top: float                           # page top to the first row, mm


def _profile_keys(f: dict, rows, trace, W: float) -> list:
    """Panel D's signals, in TRACK_ORDER: those with a trace on this genome's figure (none without a call on a
    chromosome at least W, the plot's span, long)."""
    if not any(r.call and r.L >= W for r in rows):
        return []
    out = set()
    if trace and f.get("hic") is not None:
        out.add(trace)
    if f.get("meth") is not None and f["meth"].get("pooled_adj"):
        out.add("meth")
    if (f.get("cpg") or {}).get("chroms"):
        out.add("cpg")
    tel = f.get("tel") or {}
    if (tel.get("signature") is not None or tel.get("marker")) and any(len(cu._mids(tel, r.chrom)) for r in rows):
        out.add("te")
    a = f["seq"].arrays
    if len(a) and a.role.isin(CEN_ROLES).any():
        out.add("sat")
    return [k for k in sg.TRACK_ORDER if k in out]


def _has(f: dict) -> dict:
    """The inputs signals.present asks about: LTR elements, Hi-C, CHG methylation."""
    return dict(te=bool((f.get("tel") or {}).get("density")), hic=f.get("hic") is not None,
                meth=f.get("meth") is not None)


def _keys(f: dict, rows) -> list:
    """A genome's close-up tracks (signals.present, in TRACK_ORDER): the close-ups' keys, panel B's columns and the
    signals whose places the close-up windows (A's gray boxes) hold. Every signal that ran is in panel C."""
    return sg.present(rows, f.get("summary"), _has(f))


def _setup(f: dict, rows) -> _Genome:
    summary = f.get("summary")
    clade = sg.clade_of(f["calls"], summary)
    keys = _keys(f, rows)
    W, wins = cu.window(rows, keys, [r.L for r in rows])
    trace = cu._trace_key(keys)
    tests = sorted(sg.signal_tests(summary, len(rows), clade), key=lambda t: sg.TRACK_ORDER.index(t.key))
    conf_w = max([cu._w("confidence")] + [cu._w(_conf_word(r)) + 3.2 for r in rows]) + 1.0
    b_w = COL_MM * len(keys) + 1.0 + conf_w
    b_x0 = cu.PAGE_W - 1.0 - b_w
    a_x0 = max(cu._w(r.chrom, cu.FS) for r in rows) + 1.0 + NAME_GAP
    heads = [_head(k, clade) for k in keys]
    return _Genome(f=f, rows=rows, summary=summary, clade=clade, keys=keys, tests=tests,
                   verdict={t.key: t.verdict for t in tests}, W=W, wins=dict(zip((r.chrom for r in rows), wins)),
                   trace=trace, hsc=cu._hic_scale(f, [r.chrom for r in rows], trace), te="te" in keys,
                   pkeys=_profile_keys(f, rows, trace, W), Lmax=max(r.L for r in rows),
                   row_mm=ROW_HIC if trace else ROW_BASE, a_x0=a_x0, a_w=b_x0 - AB_GAP - a_x0, b_x0=b_x0, b_w=b_w,
                   top=HEAD_TOP + max(cu._w(h) for h in heads) + 1.5)


def _head(key: str, clade: str) -> str:
    return B_HEAD[key].format(clade=clade)


def _conf_word(r) -> str:
    return r.conf if r.call else "no call"


def _p_text(p) -> str:
    """P as the page prints it ('0.87', '0.0037', '2e-5'); '' where the summary gives none (None or NaN)."""
    if p is None or not math.isfinite(p):
        return ""
    return f"{p:.2g}" if p >= 0.00095 else f"{p:.0e}".replace("e-0", "e-")


def _question(t: sg.TestRow) -> str:
    """The test's question with its P: '(chance: 8)' becomes '(chance: 8; P = 0.87)'."""
    p = _p_text(t.p)
    if not p:
        return t.question
    if t.question.endswith(")") and "(chance:" in t.question:
        return f"{t.question[:-1]}; P = {p})"
    return f"{t.question} (P = {p})"


# ------------------------------------------------------------------------------------------------ panel D data
def _bin_mean(y: np.ndarray, bs: float, a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Mean of the finite bins (bin k: [k bs, (k+1) bs)) whose centers lie in [a, b), per window; NaN without one."""
    ok = np.isfinite(y)
    s = np.concatenate([[0.0], np.cumsum(np.where(ok, y, 0.0))])
    c = np.concatenate([[0], np.cumsum(ok)])
    k0 = np.clip(np.ceil(a / bs - 0.5).astype(np.int64), 0, len(y))
    k1 = np.clip(np.ceil(b / bs - 0.5).astype(np.int64), 0, len(y))
    n = c[k1] - c[k0]
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(n > 0, (s[k1] - s[k0]) / np.maximum(n, 1), np.nan)


def _covered(starts: np.ndarray, ends: np.ndarray, x: np.ndarray) -> np.ndarray:
    """bp of the merged intervals [starts, ends) (sorted, disjoint) left of x."""
    cum = np.concatenate([[0], np.cumsum(ends - starts)])
    i = np.searchsorted(starts, x, side="right") - 1
    j = np.clip(i, 0, max(len(starts) - 1, 0))
    inside = np.clip(x - starts[j], 0, ends[j] - starts[j]) if len(starts) else np.zeros_like(x)
    return np.where(i >= 0, cum[np.clip(i, 0, None)] + inside, 0)


def _merged(s, e):
    o = np.argsort(s, kind="stable")
    s, e = np.asarray(s, np.int64)[o], np.asarray(e, np.int64)[o]
    out = []
    for a, b in zip(s, e):
        if out and a <= out[-1][1]:
            out[-1][1] = max(out[-1][1], b)
        else:
            out.append([a, b])
    arr = np.asarray(out, np.int64).reshape(-1, 2)
    return arr[:, 0], arr[:, 1]


def _signal_fn(f: dict, key: str, W: float):
    """(step, value) for one panel D signal: the window width (bp) and value(chrom, a, b) -> its value over each window
    [a, b) (None for a chromosome without data)."""
    if key in ("hic", "hub"):
        bs = f["hic"]["bins"].size

        def value(c, a, b):
            t = _track(f["hic"], c, cu.TRACE[key])
            return None if t is None else _bin_mean(t[1], bs, a, b)
    elif key == "meth":
        from .meth import params as mp
        bs = mp.BIN

        def value(c, a, b):
            t = _meth_track(f["meth"], c)
            return None if t is None else _bin_mean(t[1], bs, a, b)
    elif key == "cpg":
        from .cdr import BIN as bs

        def value(c, a, b):
            ch = f["cpg"]["chroms"].get(c)
            return None if ch is None else 100 * _bin_mean(np.asarray(ch.y, float), bs, a, b)
    elif key == "te":
        bs = 0

        def value(c, a, b):
            m = cu._mids(f.get("tel") or {}, c)
            return (np.searchsorted(m, b) - np.searchsorted(m, a)) / (b - a) * 1e6
    else:
        bs = 0
        arr = f["seq"].arrays
        arr = arr[arr.role.isin(CEN_ROLES)]

        def value(c, a, b):
            x = arr[arr.chrom == c]
            s, e = _merged(x.start.to_numpy(), x.end.to_numpy()) if len(x) else (np.zeros(0, np.int64),) * 2
            return (_covered(s, e, b) - _covered(s, e, a)) / (b - a) * 100
    step = max(1, round(W / 48 / bs)) * bs if bs else W / 48
    return step, value


def _profile(f: dict, called, key: str, W: float, rand) -> dict:
    """One panel D signal: x (bp from the centromere), every centromere's trace (obs), their mean, each random place's
    trace (null: draw x chromosome x position), the median and 5-95% band of the mean over chromosomes at each draw,
    the centromeres with data (n) and how many of them the signal placed (placed: D_PLACED)."""
    step, value = _signal_fn(f, key, W)
    half = int((W - step) / 2 // step)                       # every window [x - step/2, x + step/2) inside +-W/2
    x = np.arange(-half, half + 1) * step
    obs = np.full((len(called), len(x)), np.nan)
    nul = np.full((N_NULL, len(called), len(x)), np.nan)
    for i, r in enumerate(called):
        for out, pos in ((obs[i], (r.call[0] + r.call[1]) / 2 + x), (nul[:, i], rand[i][:, None] + x[None, :])):
            a, b = pos - step / 2, pos + step / 2
            v = value(r.chrom, a, b)
            if v is not None:
                out[...] = np.where((a >= 0) & (b <= r.L), v, np.nan)
    with warnings.catch_warnings(), np.errstate(all="ignore"):
        warnings.simplefilter("ignore", RuntimeWarning)
        mean = np.nanmean(obs, axis=0)
        lo, med, hi = np.nanpercentile(np.nanmean(nul, axis=1), [5, 50, 95], axis=0)
    has = np.isfinite(obs).any(axis=1)
    placed = np.array([bool(r.placed & D_PLACED[key]) for r in called], bool)
    return dict(x=x, obs=obs, mean=mean, null=nul, lo=lo, med=med, hi=hi, n=int(has.sum()),
                placed=int((has & placed).sum()))


def _top(p: dict, both: bool) -> float:
    """A panel D plot's top: the 99th percentile of the traces, at least the mean and the band (absolute values where
    the axis runs both ways); 0 without data."""
    obs, band = p["obs"][np.isfinite(p["obs"])], np.concatenate([p["mean"], p["hi"], p["lo"]])
    band = band[np.isfinite(band)]
    if both:
        obs, band = np.abs(obs), np.abs(band)
    return max([float(np.percentile(obs, 99))] * bool(len(obs)) + [float(band.max())] * bool(len(band)) + [0.0])


def _profiles(f: dict, rows, ax_by_key: dict, seed: int = 0, W: float | None = None) -> dict:
    """Panel D: each signal of ax_by_key around every centromere against N_NULL random places per chromosome (numpy
    default_rng(seed), centers drawn so each window lies inside its chromosome), drawn into its axes; W, the span (the
    close-up width when None). Centromeres and random places alike on the chromosomes at least W long, so both cover
    the same chromosomes. Returns each signal's profile (_profile)."""
    if not ax_by_key:
        return {}
    if W is None:
        W = cu.window(rows, _keys(f, rows), [r.L for r in rows])[0]
    called = [r for r in rows if r.call and r.L >= W]
    if not called:
        return {}
    rng = np.random.default_rng(seed)
    rand = [rng.uniform(W / 2, r.L - W / 2, N_NULL) for r in called]
    clade = sg.clade_of(f["calls"], f.get("summary"))
    kb = W < 2e6
    out = {}
    for key, ax in ax_by_key.items():
        p = out[key] = _profile(f, called, key, W, rand)
        xs = p["x"] / (1e3 if kb else 1e6)
        col = D_COLOR[key]
        ax.fill_between(xs, p["lo"], p["hi"], color="#DADADA", lw=0, zorder=1)
        ax.plot(xs, p["med"], color="#8C8C8C", lw=0.7, ls=(0, (2, 1)), zorder=2)
        for o in p["obs"]:
            ax.plot(xs, o, color=col, lw=0.35, alpha=0.3, zorder=2)
        ax.plot(xs, p["mean"], color=col, lw=1.3, zorder=4)
        ax.axvline(0, color="#C8C8C8", lw=0.4, zorder=0)
        top = _top(p, key == "meth") * 1.08 or 1.0
        if key == "meth":
            ax.set_ylim(-top, top)
            ax.axhline(0, color="#C8C8C8", lw=0.4, zorder=0)
        else:
            ax.set_ylim(0, min(top, 105.0) if key in ("sat", "cpg") else top)
            if key in ("hic", "hub"):
                ax.axhline(1, color="#C8C8C8", lw=0.4, zorder=0)
        ax.set_xlim(-W / 2 / (1e3 if kb else 1e6), W / 2 / (1e3 if kb else 1e6))
        ax.set_xlabel(f"distance from the\ncentromere ({'kb' if kb else 'Mb'})", labelpad=1)
        ax.set_ylabel(D_LABEL[key].format(clade=clade), labelpad=1.5)
        for sp in ("right", "top"):
            ax.spines[sp].set_visible(False)
        ax.text(1.0, 1.03, f"placed {p['placed']} of {p['n']}", transform=ax.transAxes, ha="right", va="bottom",
                fontsize=cu.FSS, gid=f"placed:{key}")
    return out


# ------------------------------------------------------------------------------------------------ genome pages: drawing
VERDICT_STYLE = {"trusted": ("black", "bold"), "not trusted here": (sg.ELSE, "normal"),
                 "a lead only": ("#555555", "bold"), "not used (other plants)": (sg.GRAY, "normal"),
                 "not reported": (sg.GRAY, "normal")}


def _letter(pg, x: float, y: float, s: str) -> None:
    pg.text(x, y, s, fontsize=9, fontweight="bold", va="top", ha="left", gid=f"letter:{s}")


def _gaps(f: dict, trace: str | None, r) -> list:
    """The stretches of r's Hi-C trace without data (the whole chromosome where it has no trace)."""
    if not trace:
        return []
    t = _track(f.get("hic"), r.chrom, cu.TRACE[trace])
    return [(0, r.L)] if t is None else [(max(s, 0), min(e, r.L)) for s, e in _nan_runs(*t)]


def _no_data(g: _Genome, r) -> list:
    """The stretches of r's Hi-C lane without data (the whole lane where the chromosome has no trace)."""
    return _gaps(g.f, g.trace, r)


def _panel_a(pg, g: _Genome, chunk, y0: float) -> None:
    """Panels A and B for one page's chromosomes, their first row y0 mm under the page top."""
    f, rm = g.f, g.row_mm
    h = rm * len(chunk)
    ax = pg.ax(g.a_x0, y0, g.a_w, h + 0.3)
    ax.set_gid("A")
    axb = pg.ax(g.b_x0, y0, g.b_w, h + 0.3)
    axb.set_gid("B")
    x0, x1 = -g.Lmax * 0.003, g.Lmax * 1.003
    ax.set_xlim(x0, x1)
    ax.set_ylim(h + 0.3, 0)                                  # mm from the first row's top
    bp_mm = (x1 - x0) / g.a_w
    arrays, tel = f["seq"].arrays, f.get("tel") or {}
    telo, rdna = f.get("telo") or {}, f.get("rdna") or {}
    te_alpha = 0.4 if g.verdict.get("te") in UNTRUSTED else 1.0
    hic_alpha = 0.45 if g.verdict.get(g.trace) in UNTRUSTED else 1.0
    for i, r in enumerate(chunk):
        c, L = r.chrom, r.L
        foot = (i + 1) * rm
        top, bot = foot - BAR_Y - BAR_H, foot - BAR_Y        # the bar's top and bottom edges
        lo, hi = g.wins[c]
        s, e = max(lo, 0), min(hi, L)
        if e > s:
            ax.add_patch(Rectangle((s, i * rm + 0.15), e - s, rm - 0.3, color=WIN_BAND, lw=0, zorder=0,
                                   gid=f"window:{c}"))
        ax.add_patch(Rectangle((0, top), L, BAR_H, fc="#F4F4F4", ec="#9A9A9A", lw=0.3, zorder=1))
        a = arrays[arrays.chrom == c]
        if len(a):
            cols = [ROLE_COLOR.get(x, OTHER) for x in a.role]
            o = np.argsort([cc != OTHER for cc in cols], kind="stable")          # gray context under the colors
            st, en = a.start.to_numpy()[o], a.end.to_numpy()[o]
            verts = [[(p, top), (q, top), (q, bot), (p, bot)]
                     for p, q in (cu._widen(p, q, MIN_FEAT * bp_mm, L) for p, q in zip(st, en))]
            ax.add_collection(PolyCollection(verts, facecolors=[cols[j] for j in o], edgecolors="none", zorder=2))
        if g.te:
            m = cu._mids(tel, c)
            if len(m):
                seg = np.stack([np.column_stack([m, np.full(len(m), top - 0.1)]),
                                np.column_stack([m, np.full(len(m), top - 0.1 - TICK_H)])], axis=1)
                ax.add_collection(LineCollection(seg, colors=TE, lw=0.45, alpha=te_alpha, zorder=2, gid=f"te:{c}"))
        if g.trace:
            base = top - 0.15
            ax.plot([0, L], [base, base], color="#D8D8D8", lw=0.3, zorder=1)
            t = _track(f.get("hic"), c, cu.TRACE[g.trace])
            if t is not None:
                ax.plot(t[0], base - np.clip(t[1] / g.hsc, 0, 1) * HIC_H, color=HIC, lw=0.5, alpha=hic_alpha,
                        zorder=3, solid_joinstyle="round", gid=f"hic:{c}")
            for p, q in _no_data(g, r):
                ax.add_patch(Rectangle((p, base - HIC_H), q - p, HIC_H, fill=False, hatch=cu.HATCH, ec=cu.NODATA,
                                       lw=0, zorder=1))
        if r.call:
            ax.plot((r.call[0] + r.call[1]) / 2, foot - TRI_Y, marker="^", ms=4.2, mew=0.6, ls="", zorder=6,
                    mfc=cu.CONF_FACE.get(r.conf, "white"), mec=cu.CONF_EDGE.get(r.conf, "black"), clip_on=False,
                    gid=f"centromere:{c}")
        for k in g.keys:                                     # where a trusted signal that points elsewhere puts it
            p = r.places.get(k)
            if r.roles.get(k) == "elsewhere" and p and g.verdict.get(k) not in UNTRUSTED:
                ax.plot((p[0] + p[1]) / 2, foot - TRI_Y, marker="x", ms=3.4, mew=1.0, color=sg.SIGNALS[k].color,
                        ls="", zorder=5.5, clip_on=False, gid=f"elsewhere:{c}:{k}")
        for tc, side, *_ in telo.get("ends", []):
            if tc == c:
                ax.plot(0 if side == "start" else L, top, marker=">" if side == "start" else "<", ms=3.2,
                        color="black", mew=0, zorder=5, clip_on=False)
        for tc, ts, te_, *_ in telo.get("interstitial", []):
            if tc == c:
                ax.plot((ts + te_) / 2, top - 0.55, marker="v", ms=3.0, color="black", mew=0, zorder=5, clip_on=False)
        for arr, face in ((rdna.get("nors", ()), RDNA_MARK), (rdna.get("five_s", ()), "white")):
            for rc, rs, re_, *_ in arr:                      # 45S filled, 5S open
                if rc == c:
                    ax.plot((rs + re_) / 2, top - 0.55, marker="D", ms=2.0, mew=0.5, mec=RDNA_MARK, mfc=face, ls="",
                            zorder=5, clip_on=False)
        ax.text(x0 - 1.0 * bp_mm, (top + bot) / 2, c, ha="right", va="center", fontsize=cu.FS, gid=f"name:{c}")
        for j, k in enumerate(g.keys):                       # B: each signal's role, then the confidence
            axb.plot(COL_MM * (j + 0.5), (top + bot) / 2, ls="", gid=f"role:{c}:{k}",
                     **cu.ROLE_GLYPH[r.roles.get(k, "not")])
        xc = COL_MM * len(g.keys) + 1.0 + 1.4
        if r.call:
            axb.plot(xc, (top + bot) / 2, marker="s", ms=4.0, mfc=cu.CONF_FACE.get(r.conf, "white"),
                     mec=cu.CONF_EDGE.get(r.conf, "black"), mew=0.7, ls="")
        axb.text(xc + 1.8, (top + bot) / 2, _conf_word(r), va="center", fontsize=cu.FSS,
                 color="black" if r.call else sg.GRAY)
    ax.set_yticks([])
    for sp in ("left", "right", "top"):
        ax.spines[sp].set_visible(False)
    kb = g.Lmax < 2e6
    ticks = _nice_ticks(g.Lmax)
    ax.set_xticks(ticks)
    ax.set_xticklabels([f"{t / (1e3 if kb else 1e6):g}" for t in ticks])
    ax.set_xlabel(f"position on the chromosome ({'kb' if kb else 'Mb'})", labelpad=1)
    axb.set_xlim(0, g.b_w)
    axb.set_ylim(h + 0.3, 0)
    axb.axis("off")
    for j, k in enumerate(g.keys):
        x, off = COL_MM * (j + 0.5), g.verdict.get(k) in UNTRUSTED
        axb.plot([x, x], [0, h], color="#E4E4E4", lw=0.4, zorder=0)
        if off:                                              # failed its genome test: the column grayed
            axb.add_patch(Rectangle((x - COL_MM * 0.42, 0), COL_MM * 0.84, h, color="#F1F1F1", lw=0, zorder=-1))
        pg.text(g.b_x0 + x, y0 - 1.0, _head(k, g.clade), rotation=90, ha="center", va="bottom", fontsize=cu.FSS,
                color="#8C8C8C" if off else "black", gid=f"column:{k}")
    pg.text(g.b_x0 + COL_MM * len(g.keys) + 1.0, y0 - 1.0, "confidence", va="bottom", fontsize=cu.FSS)


def _c_layout(g: _Genome) -> tuple[dict, list, float]:
    """Panel C's columns (x, mm), its rows (test, question lines, height) and its height under its letter."""
    xl = 4.6                                                 # a color chip, then the label
    xv = cu.PAGE_W - 1.0 - max(cu._w(t.verdict, cu.FS, "bold") for t in g.tests)
    xb = xv - 3.0 - TEST_BAR
    xq = xl + max(cu._w(t.label, cu.FS) for t in g.tests) + 3.0
    rows = []
    for t in g.tests:
        q = cu._wrap(_question(t), xb - 3.0 - xq)
        rows.append((t, q, TEST_ROW + (len(q) - 1) * cu.LINE_H))
    return dict(label=xl, question=xq, bar=xb, verdict=xv), rows, 1.0 + TEST_ROW + sum(h for *_, h in rows)


def _panel_c(pg, g: _Genome, y: float) -> float:
    """Panel C, its letter y mm under the page top; returns its foot."""
    x, rows, height = _c_layout(g)
    _letter(pg, 1.0, y, "C")
    yh = y + 1.0 + TEST_ROW / 2
    for k, lab in (("label", "signal"), ("question", "does it mark the centromeres in this genome?"),
                   ("bar", "chromosomes"), ("verdict", "verdict")):
        pg.text(x[k], yh, lab, fontsize=cu.FSS, va="center", color="#555555")
    y0 = y + 1.0 + TEST_ROW
    axc = pg.ax(x["bar"], y0, TEST_BAR, height - 1.0 - TEST_ROW)
    axc.set_gid("C")
    axc.set_xlim(0, 1)
    axc.set_ylim(height - 1.0 - TEST_ROW, 0)
    axc.axis("off")
    yy = 0.0
    for t, q, h in rows:
        ym = yy + TEST_ROW / 2                               # the first line's center
        hue = sg.SIGNALS[t.key].color
        pg.fig.add_artist(Rectangle((1.8 / cu.PAGE_W, 1 - (y0 + ym + 1.0) / pg.h), 1.8 / cu.PAGE_W, 2.0 / pg.h,
                                    color=hue, lw=0, transform=pg.fig.transFigure))
        pg.text(x["label"], y0 + ym, t.label, fontsize=cu.FS, va="center", gid=f"test-label:{t.key}")
        pg.text(x["question"], y0 + ym - cu.LINE_H / 2, "\n".join(q), fontsize=cu.FSS, va="top", linespacing=0.95,
                gid=f"test-question:{t.key}")
        color, weight = VERDICT_STYLE.get(t.verdict, ("black", "normal"))
        pg.text(x["verdict"], y0 + ym, t.verdict, fontsize=cu.FS, va="center", color=color, fontweight=weight,
                gid=f"test-verdict:{t.key}")
        if t.k is not None and t.n:
            axc.add_patch(Rectangle((0, ym - 1.0), 1, 2.0, fc="white", ec="#BBBBBB", lw=0.4))
            axc.add_patch(Rectangle((0, ym - 1.0), min(t.k / t.n, 1), 2.0, fc=hue, ec="none", alpha=0.85))
            if t.chance is not None:
                cx = min(t.chance / t.n, 1)
                axc.plot([cx, cx], [ym - 1.4, ym + 1.4], color="black", lw=0.9, solid_capstyle="butt",
                         gid=f"chance:{t.key}")
        yy += h
    return y + height


def _panel_d(pg, g: _Genome, y: float) -> float:
    """Panel D, its letter y mm under the page top; returns its foot."""
    _letter(pg, 1.0, y, "D")
    n = len(g.pkeys)
    ew = min(D_MAX_W, (cu.PAGE_W - D_GUT - 2.0 - D_GUT * (n - 1)) / n)
    axes = {}
    for j, k in enumerate(g.pkeys):
        axes[k] = pg.ax(D_GUT + j * (ew + D_GUT), y + 1.0 + D_TOP, ew, D_H)
        axes[k].set_gid(f"D:{k}")
    _profiles(g.f, g.rows, axes, seed=0, W=g.W)
    return y + 1.0 + D_TOP + D_H + D_FOOT


def _lower_h(g: _Genome) -> float:
    """Panels C and D's height, with the space above each."""
    return BLOCK_GAP + _c_layout(g)[2] + (BLOCK_GAP + 1.0 + D_TOP + D_H + D_FOOT if g.pkeys else 0.0)


def _lower(pg, g: _Genome, y: float) -> float:
    y = _panel_c(pg, g, y + BLOCK_GAP)
    return _panel_d(pg, g, y + BLOCK_GAP) if g.pkeys else y


def _key_a(g: _Genome) -> list:
    """The key of panels A and B: what the genome's rows draw."""
    f, arr = g.f, g.f["seq"].arrays
    roles = set(arr.role)
    items = [Patch(fc=ROLE_COLOR[x], ec="none", label=cu.SAT_NAME[x]) for x in ROLE_COLOR if x in roles]
    if len(arr) and (~arr.role.isin(list(ROLE_COLOR))).any():
        items.append(Patch(fc=OTHER, ec="none", label="other tandem repeat"))
    if g.te:
        items.append(Line2D([], [], ls="", marker="|", ms=4, mew=0.8, color=TE, label=f"{g.clade} element"))
    if g.trace:
        items.append(Line2D([], [], color=HIC, lw=0.8, label=f"Hi-C contact, 0–{g.hsc:.1f}× median"))
        if any(_no_data(g, r) for r in g.rows):
            items.append(Patch(fc="none", ec=cu.NODATA, hatch=cu.HATCH, lw=0, label="no data"))
    items += [Line2D([], [], ls="", marker="^", ms=4.2, mfc=cu.CONF_FACE[c], mec=cu.CONF_EDGE[c], mew=0.6,
                     label=f"centromere, {c} confidence")
              for c in ("high", "medium", "low") if any(r.call and r.conf == c for r in g.rows)]
    if any(r.roles.get(k) == "elsewhere" and r.places.get(k) and g.verdict.get(k) not in UNTRUSTED
           for r in g.rows for k in g.keys):
        items.append(Line2D([], [], ls="", marker="x", ms=3.4, mew=1.0, color="#555555",
                            label="where another signal points"))
    items.append(Patch(fc=WIN_BAND, ec="none", label="close-up in panel E"))
    telo, rdna = f.get("telo") or {}, f.get("rdna") or {}
    if telo.get("ends"):
        items.append(Line2D([], [], ls="", marker=">", ms=3.2, mew=0, color="black", label="telomere array"))
    if telo.get("interstitial"):
        items.append(Line2D([], [], ls="", marker="v", ms=3.0, mew=0, color="black",
                            label="interstitial telomeric array"))
    for has, face, label in ((rdna.get("nors"), RDNA_MARK, "45S rDNA array (NOR)"),
                             (rdna.get("five_s"), "white", "5S rDNA array")):
        if has:
            items.append(Line2D([], [], ls="", marker="D", ms=2.6, mew=0.5, mec=RDNA_MARK, mfc=face, label=label))
    return items + cu.role_key({r.roles.get(k, "not") for r in g.rows for k in g.keys})


def _key_cd(g: _Genome) -> list:
    """The key of panels C and D."""
    items = []
    if any(t.chance is not None and t.k is not None and t.n for t in g.tests):
        items.append(Line2D([], [], ls="", marker="|", ms=6, mew=1.0, color="black", label="expected by chance"))
    if g.pkeys:
        items += [Line2D([], [], color="#555555", lw=1.3, label="average of all centromeres"),
                  Line2D([], [], color="#8C8C8C", lw=0.7, ls=(0, (2, 1)), label="random places, median"),
                  Patch(fc="#DADADA", ec="none", label="random places, 5–95%")]
    return items


def _key_h(items: list) -> float:
    return cu.KEY_ROW * math.ceil(len(items) / cu._key_columns(items)) + cu.KEY_PAD if items else 0.0


def _key(pg, items: list, y: float) -> None:
    if items:
        pg.fig.legend(handles=items, loc="upper left", bbox_to_anchor=(1.5 / cu.PAGE_W, 1 - y / pg.h),
                      ncol=cu._key_columns(items), frameon=False, fontsize=cu.FSS, handlelength=1.4,
                      columnspacing=1.0, handletextpad=0.4, borderaxespad=0)


def _genome_pages(f: dict, rows: list) -> Iterator[Figure]:
    """The genome pages of one genome, each drawn when asked for (the caller saves and closes a page before asking for
    the next): panels A and B, ROWS_PER_PAGE chromosomes a page; C and D under the last rows, or on a page of their own
    where that page would pass closeup.PAGE_H; the key on every page. `rows` come from signals.rows."""
    if not rows:
        return
    g = _setup(f, rows)
    items_a, items_cd = _key_a(g), _key_cd(g)
    lower = _lower_h(g)
    chunks = [rows[i:i + ROWS_PER_PAGE] for i in range(0, len(rows), ROWS_PER_PAGE)]
    alone = False
    for pi, chunk in enumerate(chunks):
        last = pi == len(chunks) - 1
        ya = g.top + g.row_mm * len(chunk) + AXIS_H
        under = last and ya + lower + KEY_GAP + _key_h(items_a + items_cd) <= cu.PAGE_H
        items = items_a + items_cd if under else items_a
        pg = cu._Page(ya + (lower if under else 0.0) + KEY_GAP + _key_h(items))
        _letter(pg, 1.0, 1.0, "A")
        _letter(pg, g.b_x0 - 3.5, 1.0, "B")
        _panel_a(pg, g, chunk, g.top)
        _key(pg, items, (_lower(pg, g, ya) if under else ya) + KEY_GAP)
        alone = last and not under
        yield pg.fig
    if alone:
        pg = cu._Page(1.0 + lower - BLOCK_GAP + KEY_GAP + _key_h(items_cd))
        _key(pg, items_cd, _lower(pg, g, 1.0 - BLOCK_GAP) + KEY_GAP)
        yield pg.fig


def _page_iter(f: dict) -> Iterator[tuple[str, Figure]]:
    """The figure pages of one genome as (kind, page), kind "genome" (A-D) then "closeup" (E), each drawn when asked
    for. The caller saves and closes a page before asking for the next, so one is open at a time: matplotlib warns at
    20 open figures, and the 670-sequence Cissus run (28 pages) held every page, with its memory, until the last was
    drawn."""
    hic = f.get("hic")
    rows = sg.rows(f["calls"], f.get("summary"), hic["bins"].size if hic is not None else None)
    for fig in _genome_pages(f, rows):
        yield "genome", fig
    for fig in cu.pages(f, rows, _keys(f, rows)):
        yield "closeup", fig


def _pages(f):
    """Every page of one genome at once, all open (the caller closes them): to inspect pages; `overview` uses
    `_page_iter`, which holds one."""
    return [fig for _, fig in _page_iter(f)]



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


def _plant_rows(ax, f, chrom, plants_, prov, hsc):
    """G: the reference plant's evidence on top, then one row per plant: its adjusted-CHG track (down = loss), the
    sites it was scored at, its call site and its best centromere-like site away from the call. hsc: the top of the
    Hi-C traces (the contact with the other centromeres, the reference's and each plant's own)."""
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
        ax.add_patch(Rectangle((s, 0), e - s, ytop, color=CALL_BAND, lw=0, zorder=0))
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
        yv = np.clip(tr[1][k] / hsc, 0, 1) * 0.55
        ax.plot(tr[0][k], y0 + 0.9 + yv, color=HIC, lw=0.7, zorder=3, gid="G:hic")
        for s, e in _nan_runs(tr[0][k], tr[1][k], 1):
            ax.add_patch(Rectangle((s, y0 + 0.88), e - s, 0.6, fill=False, hatch="////", ec="#C8C8C8", lw=0,
                                   zorder=1))
    x_lab = lo - (hi - lo) * 0.01
    clade = sg.clade_of(f["calls"], f.get("summary"))
    for yy, name in ((y0 + 0.16, "satellite arrays"), (y0 + 0.58, f"{clade} elements"), (y0 + 1.12, "Hi-C contact")):
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
            ax.plot(lt[0][k], yb + 0.5 + np.clip(lt[1][k] / hsc, 0, 1) * 0.45, color=HIC, lw=0.6,
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
    top = max(17.0, 6 + max(_text_mm(c, 5.5) for c in chroms))  # mm: verdict key, letter E, chromosome labels
    f_mm = (12 + 8.5 * n_pl + 9) if has_f else 0
    H = top + grid_mm + (13 + f_mm if has_f else 0) + 4
    fig = plt.figure(figsize=(WIDTH, H * MM))
    axe = fig.add_axes([0.22, 1 - (top + grid_mm) / H, 0.76, grid_mm / H])
    plants_ = _verdict_grid(axe, ptab, chroms, prov)
    fig.text(0.01, 1 - 10 / H, "F", fontsize=9, fontweight="bold", va="bottom")
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
        hsc = cu._hic_scale(f, chroms, "hic")                # the top of the trace G draws
        f["focus_scale"] = _plant_rows(axf, f, focus, plants_, prov, hsc)
        fig.text(0.01, yb + f_mm / H + 5 / H, "G", fontsize=9, fontweight="bold", va="bottom")
        keys = [Line2D([], [], color=METH, lw=0.6, label="plant's CHG methylation (down = loss)"),
                Line2D([], [], ls="", marker="^", ms=4, mew=0, color=VERDICT_COLOR["shifted"],
                       label="centromere-like elsewhere"),
                Line2D([], [], color="#808080", lw=0.6, label="site scored"),
                Patch(fc=CALL_BAND, ec="none", label="centromere"),
                Patch(fc="none", ec="#C8C8C8", hatch="////", lw=0, label="no data")]
        if f.get("libraries"):
            keys.insert(1, Line2D([], [], color=HIC, lw=0.6, label="Hi-C (reference or the plant's own)"))
        fig.legend(handles=keys, loc="lower left", bbox_to_anchor=(0.04, yb + f_mm / H + 0.5 / H), ncol=len(keys),
                   frameon=False, fontsize=6, handlelength=1.4, columnspacing=1.0, handletextpad=0.4,
                   borderaxespad=0)
    return fig


# ------------------------------------------------------------------------------------------------ legends
LEGEND_PT = 10.0                                 # layout: legend text, pt (a long legend takes a page more, never less)
LEGEND_LH = 1.3                                  # layout: a legend line, in font sizes
LEGEND_SIDE = 8.0                                # layout: a legend page's left and right margins, mm
LEGEND_TB = 7.0                                  # layout: a legend page's top and bottom margins, mm
LEGEND_PAR = 0.5                                 # layout: space between two paragraphs, in lines
_MARKUP = re.compile(r"(\*\*[^*]+\*\*|\*[^*]+\*)")
_STYLE = {"normal": ("normal", "normal"), "bold": ("bold", "normal"), "italic": ("normal", "italic")}
CONF_LEVELS = ("high", "medium", "low")
# a signal that alone placed some centromeres, in the opening sentence (any Hi-C readout: "Hi-C")
ALONE = {"sat": "the satellite array", "te": "the {clade} cluster", "meth": "CHG methylation",
         "cpg": "the CpG methylation dip", "at": "the AT-rich trough", "rep": "repeat density",
         "fam": "the repeat-family cluster"}


def _species(name: str) -> str:
    return name.replace("_", " ")


def _and(xs: list, comma: bool = False) -> str:
    """'a', 'a and b', 'a, b and c' (with comma: 'a, and b', for clauses)."""
    xs = list(xs)
    if len(xs) < 2:
        return "".join(xs)
    return ", ".join(xs[:-1]) + (", and " if comma else " and ") + xs[-1]


@dataclass
class _Legend:
    """What both legends read: the genome's rows (signals.rows), its close-up tracks, panel C's verdicts, the Hi-C
    trace and panel D's signals, as the figure pages draw them."""
    f: dict
    name: str
    rows: list
    keys: list
    verdict: dict
    clade: str
    trace: str | None
    pkeys: list
    hic_bin: int | None
    tests: list                          # panel C's rows


def _legend_setup(f: dict) -> _Legend:
    hic = f.get("hic")
    hic_bin = hic["bins"].size if hic is not None else None
    rows = sg.rows(f["calls"], f.get("summary"), hic_bin)
    keys = _keys(f, rows)
    clade = sg.clade_of(f["calls"], f.get("summary"))
    trace = cu._trace_key(keys)
    tests = sg.signal_tests(f.get("summary"), len(rows), clade)
    W = cu.window(rows, keys, [r.L for r in rows])[0] if rows else 0.0
    return _Legend(f=f, name=_species(f["genome"].name), rows=rows, keys=keys,
                   verdict={t.key: t.verdict for t in tests}, clade=clade, trace=trace,
                   pkeys=_profile_keys(f, rows, trace, W) if rows else [], hic_bin=hic_bin, tests=tests)


def _result(lg: _Legend) -> str:
    """The opening: how many chromosomes have a centromere and at which confidence, the signals that agree on them
    (the commonest sets of used signals, any Hi-C readout named once), and each signal that points elsewhere."""
    n, called = len(lg.rows), [r for r in lg.rows if r.call]
    if not called:
        return f"In {lg.name}, no centromere was called on " + \
            ("its 1 chromosome." if n == 1 else f"any of its {n} chromosomes.")
    k = len(called)
    conf = Counter(r.conf for r in called)
    levels = [c for c in CONF_LEVELS if conf.get(c)]
    how = ("" if not levels else f" (all {levels[0]} confidence)" if len(levels) == 1 and conf[levels[0]] == k else
           f" ({_and([f'{conf[c]} {c}' for c in levels])} confidence)")
    out = f"In {lg.name}, {k} of {n} chromosome{'s' if n != 1 else ''} {'has' if k == 1 else 'have'} a centromere{how}."
    groups, first = Counter(), {}
    for r in called:
        used = [x for x in sg.SIGNALS if r.roles.get(x) == "used"]
        names = tuple(dict.fromkeys(sg.SIGNALS[x].short.format(clade=lg.clade) for x in used))
        groups[names] += 1
        first.setdefault(names, used[0] if used else None)
    top = groups.most_common()
    said = []
    for names, c in top[:3]:
        cnt = f"all {c}" if c == k and k > 1 else str(c)
        if len(names) > 1:
            said.append(f"{_and(names)} agree on {cnt}")
        else:
            x = first[names]
            who = ("the DNA sequence" if x is None else "Hi-C" if x in sg.HIC_KEYS else
                   ALONE.get(x, sg.SIGNALS[x].short).format(clade=lg.clade))
            said.append(f"{who} alone places {cnt}")
    rest = sum(c for _, c in top[3:])
    if rest:
        said.append(f"other sets of signals place {rest}")
    out += " " + sg._cap(_and(said, comma=True)) + "."
    far = [(sg.SIGNALS[x].where.format(clade=lg.clade), m) for x in lg.keys if sg.SIGNALS[x].where
           and (m := sum(r.roles.get(x) == "elsewhere" for r in called))]
    if far:
        out += " " + sg._cap(_and([f"the {w} points elsewhere on {m}" if i == 0 else f"the {w} on {m}"
                                   for i, (w, m) in enumerate(far)])) + "."
    return out


def _short_panels(lg: _Legend) -> list:
    """The short legend's panels A-E: what each draws, every mark without a label named."""
    f, rows = lg.f, lg.rows
    drawn = ["satellite arrays colored by where their family sits in this genome"]
    if "te" in lg.keys:
        drawn.append(f"{lg.clade} elements as orange ticks")
    if lg.trace:
        drawn.append(("Hi-C contact with the other centromeres" if lg.trace == "hic" else
                      "the Hi-C contact found from Hi-C alone") + " as a blue trace")
    a = f"**A**, Every chromosome to one scale, with {_and(drawn)}."
    if any(_gaps(f, lg.trace, r) for r in rows):
        a += (" Hatching marks Hi-C without data, often inside satellite arrays, where the repeat copies are so alike "
              "that too few reads map there.")
    called = any(r.call for r in rows)
    crosses = any(r.roles.get(x) == "elsewhere" and r.places.get(x) and lg.verdict.get(x) not in UNTRUSTED
                  for r in rows for x in lg.keys)
    a += (" The triangle under each chromosome is its centromere" +
          (", a cross in a signal's color marks where a trusted signal that disagrees puts it," if crosses else "") +
          " and gray boxes are the close-ups in E." if called else " Gray boxes are the close-ups in E.")
    b = ("**B**, Each signal's role in each centromere and the centromere's confidence (the roles are defined in the "
         "second legend).")
    if any(lg.verdict.get(x) in UNTRUSTED for x in lg.keys):
        b += " A grayed column is a signal not trusted in this genome."
    if any(r.roles.get(x) == "used" and lg.verdict.get(x) in UNTRUSTED for r in rows for x in lg.keys):
        b += (" A dot in a grayed column is a signal that, though not trusted genome-wide, is still the evidence for "
              "that centromere (for example the satellite array it sits on), at lower confidence.")
    c = ("**C**, How well each signal marks centromeres across the genome, with the share of chromosomes where it "
         "does (bar) and the share expected by chance (tick, where known). A trusted signal can move a centromere or "
         "raise or lower its confidence. An untrusted one cannot, but it may still offer a low-confidence best guess "
         "where nothing trusted places one.")
    per = [t.label for t in lg.tests if t.per_chrom and t.k is not None and t.n]
    if per:
        c += f" For {_and(per)}, the bar counts the chromosomes whose own test passes."
    if "a lead only" in lg.verdict.values():
        c += " A lead only places a centromere just where nothing else does."
    out = [a, b, c]
    if lg.pkeys:
        out.append("**D**, Each signal around every centromere (thin traces) and their average (bold), against the "
                   "average at random places on the same chromosomes (gray band, 5–95% of draws, with its median "
                   "dashed), on the chromosomes at least as long as the plot is wide. Over each plot, how many of the "
                   "centromeres drawn that signal placed. A signal rises at the centromeres it placed by construction, "
                   "so independent support is a rise or dip clear of the band at centromeres placed by other "
                   "signals.")
    out.append("**E**, One close-up per chromosome, all at one scale. " +
               ("The band through every track is the centromere, and the bar at the foot of each track is where that "
                "signal alone would put it: a colored bar well off the band is a disagreement, a gray bar a signal "
                "not used, wherever it lies. Each track is labeled with its "
                "role (gray when not used), and the italic sentence says what placed the centromere and which signals "
                "disagree, with their distance." if called else
                "The bar at the foot of each track is where that signal alone would put the centromere."))
    return out


def _how(lg: _Legend) -> list:
    """The long legend's account of how cenfuse places a centromere: the roles, the confidence and each kind of
    verdict, with the `rule` values (docs/outputs.md) in parentheses."""
    f, cl = lg.f, lg.clade
    tol = f"{sg.distance_text(lg.hic_bin)} (one Hi-C bin)" if lg.hic_bin else sg.distance_text(sg.NO_HIC_TOL)
    out = [f"**How cenfuse places a centromere.** Each signal is tested genome-wide (C). The centromere stays "
           f"where the DNA sequence puts it, on a satellite array or a cluster of {cl} elements, unless trusted Hi-C "
           f"places it elsewhere or the sequence offers none (rules in parentheses, from "
           f"{f['genome'].name}.centromeres.tsv)."]
    roles = (f"**Roles.** A signal is *used* when it placed or agrees with this centromere. It *points elsewhere* "
             f"when cenfuse counted it against this centromere, or when it is the satellite array or {cl} cluster "
             f"that a moved or overruled centromere left behind, over {tol} away, having passed its test. It is "
             f"*not used* when it failed its test, gives no place, or was counted neither for nor against this "
             f"centromere, however far its place lies. Hi-C contact with the other centromeres has a place only "
             f"where it peaks.")
    if f.get("meth") is not None:
        roles += " *Not used (other plants)* marks methylation from other plants, whose centromeres may differ."
    out.append(roles)
    out.append(
        f"**Confidence.** A centromere kept on the DNA sequence starts at the level its array or cluster earns alone, "
        f"lower where the satellite family fails its test or the array lies near a chromosome end, holds little "
        f"satellite or has a rival, and higher where the {cl} cluster agrees. It rises one level for agreeing Hi-C "
        f"and one for centromere-like CHG methylation, two at most. Where trusted signals disagree enough to contest "
        f"it but not to move it, it instead drops one level from its start, two where the arms' fold or alignment and "
        f"the place Hi-C alone singles out both disagree, one more where the CpG dip disagrees too, and agreeing "
        f"methylation gives one back. A CpG dip alone puts it one level below what the agreeing signals give. A "
        f"centromere that Hi-C moved or "
        f"placed is medium, high with a second signal at its place, and low where only the flank contact, a weaker "
        f"Hi-C peak or untested Hi-C placed it (medium with an agreeing trough or repeat peak) or one of two rules "
        f"shaped on a single genome did (cis_trans_pair, cis_end_pair). One placed from the DNA alone, by an AT-rich "
        f"trough, repeat density or a repeat-family cluster, is low, a trough medium where troughs and Hi-C agree "
        f"across the genome. High is drawn black, medium gray and low open (dashed in the close-ups).")
    out.append(
        f"**How each centromere was placed.** The italic sentence over each close-up names what placed the centromere "
        f"and what agrees. *Satellite, {cl} and Hi-C agree* keeps it on its satellite array, *Satellite array alone* "
        f"with nothing else bearing on it, *{cl} and Hi-C agree, on a long tandem repeat* on an array with too little "
        f"centromeric satellite to count, and *{cl} cluster, backed by Hi-C* on the densest {cl} cluster, where no "
        f"satellite array marks it or the cluster overruled one (sequence, hic_agrees, hic_support, hic_end, "
        f"cis_agrees, te_hic_pair, meth_support). *Second satellite cluster, chosen by Hi-C* or *by the AT-rich "
        f"trough* takes the second of two similar clusters (hic_chose_alt, at_chose_alt). *Placed by Hi-C* moves it "
        f"to a Hi-C peak or places it where the sequence offers none (hic_moved, hic_strong_peak, hic_medium_peak, "
        f"hic_denovo_only, te_hic_pair, cis_only), an agreeing AT-rich trough or repeat peak narrowing it (at_pair, "
        f"rep_pair). *Placed by two Hi-C measures* needs flank contact and contact with other chromosomes to agree "
        f"(cis_trans_pair, cis_end_pair), or arms that fold back or align at the place Hi-C alone singles out "
        f"(fold_pair, align_pair). An AT-rich trough, alone or with that place, repeat density or a repeat-family "
        f"cluster places a centromere only where nothing else does, the last two as leads (at_only, at_denovo_pair, "
        f"rep_only, fam_only), and *No centromere called* leaves none (none). After the semicolon come the trusted "
        f"signals that point elsewhere, with the gap to their place, and those that contest the centromere "
        f"without moving it lower its confidence (hic_offset, fold_offset, align_offset, hic_denovo_offset, "
        f"cis_elsewhere, cpg_elsewhere). *Overruled* marks a satellite array that a {cl} cluster replaced.")
    return out


def _genome_facts(lg: _Legend) -> str:
    """The long legend's account of what this genome holds: arrays, Hi-C bins, methylation, telomeres, rDNA."""
    f, n = lg.f, len(lg.rows)
    has = [f"{n} chromosome{'s' if n != 1 else ''} in assembly order", f"{len(f['seq'].arrays):,} satellite arrays "
           f"(the smallest widened to stay visible)"]
    if lg.hic_bin:
        has.append(f"Hi-C in {lg.hic_bin / 1e3:g}-kb bins")
    meth = f.get("meth")
    if meth is not None:
        k = len(meth.get("samples") or ())
        has.append(f"CHG methylation from {k} plant{'s' if k != 1 else ''}, each bin against the same kind of DNA "
                   f"nearby (0, no change)")
    out = f"**This genome.** It has {_and(has)}."
    telo = f.get("telo") or {}
    if telo.get("unit"):
        i = len(telo.get("interstitial", []))
        out += (f" Telomere arrays of {telo['unit']} cap {len(telo.get('ends', []))} of {2 * n} chromosome ends" +
                ("" if not i else ", and 1 lies inside a chromosome" if i == 1 else f", and {i} lie inside "
                 f"chromosomes") + ".")
    else:
        out += " No telomere repeat was found in arrays at the chromosome ends."
    rd = f.get("rdna") or {}
    n45, n5 = len(rd.get("nors", ())), len(rd.get("five_s", ()))
    out += (f" It holds {n45} array{'s' if n45 != 1 else ''} of 45S rDNA and {n5} of 5S rDNA." if n45 or n5 else
            " It holds no 45S or 5S rDNA array.")
    return out


def legends(f) -> tuple[str, str]:
    """The text of the two legends, paragraphs apart by a blank line, **x** set bold and *x* italic (_text_pages).
    Short: the result, then each panel's marks in plain words (A-E; F-G with scored plants). Long: how cenfuse places
    a centromere, with the `rule` values in parentheses, what this genome holds, the plants, and where the analysis
    ran. Neither uses a word of signals.BANNED."""
    lg = _legend_setup(f)
    short = [_result(lg)] + _short_panels(lg)
    long = _how(lg) + [_genome_facts(lg)]
    ptab = f.get("plants")
    if ptab is not None and len(ptab):
        _plant_legends(f, lg, short, long)
    long.append(f"Analysis: {__import__('socket').gethostname()}:{f['out']} · Script: cenfuse/plot.py · "
                f"{date.today()}")
    return "\n\n".join(short), "\n\n".join(long)


def _plant_legends(f, lg: _Legend, short: list, long: list) -> None:
    """Panels F and G of the plants page, in both legends."""
    ptab, ver = f["plants"], f["verdicts"]
    cnt = ptab.verdict.value_counts().to_dict()
    n_pl = ptab.plant.nunique()
    prov = {}
    for (_, n), v in ver.items():
        prov.setdefault(n, v["verdict"])
    other = sorted(p for p in ptab.plant.unique() if prov.get(p) == "other")
    focus = f.get("focus")
    hues = ["orange that its centromere seems to sit elsewhere, by its methylation with support from the other "
            "plants or by a strong peak of its own Hi-C (the number gives the distance in Mb)",
            "light orange that its methylation suggests so without that support",
            "black that the plant lacks the centromere's DNA"]
    if cnt.get("conflict"):
        hues.append("purple that its methylation and Hi-C disagree")
    short.append(f"**F**, Each of {n_pl} plants (rows) on each chromosome (columns). Gray means the plant's own "
                 f"methylation or Hi-C is consistent with this genome's centromere, {_and(hues)}.")
    if focus is not None and f.get("meth") is not None:
        short.append(f"**G**, {focus}, where the plants differ most. The reference plant's evidence is on top, then "
                     f"each plant's CHG methylation on one scale (±{f.get('focus_scale', 0.1):.2f}, loss down) with a "
                     f"triangle at its most centromere-like place elsewhere.")
    names = f" ({', '.join(other)})" if other else ""
    long.append(f"**Plants.** Each input was tested for whether it comes from the reference (assembled) plant, by "
                f"genetic differences in its reads and stretches of the assembly it lacks. {len(other)} of {n_pl} "
                f"plants come from other plants{names}, and only the reference plant's inputs and those of unknown "
                f"origin decide its centromeres. In F, gray is consistent with the centromere but no proof, since some "
                f"places look centromeric in every plant. Orange needs either methylation that dips at the centromere "
                f"and looks centromeric elsewhere, supported by the other plants at a 5% false discovery rate, or a "
                f"strong peak of the plant's own Hi-C away from the centromere, where the reference Hi-C does not "
                f"peak.")


def _tokens(par: str) -> list:
    """A paragraph's words, each a list of (text, style) pieces (a word can change style inside, as '**A**,'):
    **x** bold, *x* italic, else normal."""
    words, glued = [], False
    for part in _MARKUP.split(par):
        if not part:
            continue
        style = "bold" if part.startswith("**") else "italic" if part.startswith("*") else "normal"
        text = part[2:-2] if style == "bold" else part[1:-1] if style == "italic" else part
        for i, w in enumerate(text.split(" ")):
            if w and i == 0 and glued and words:             # no space before it: part of the previous word
                words[-1].append((w, style))
            elif w:
                words.append([(w, style)])
        glued = not text.endswith(" ")
    return words


def _width(pieces) -> float:
    return sum(cu._w(t, LEGEND_PT, *_STYLE[s]) for t, s in pieces)


def _cut(text: str, style: str, first: float, width: float) -> list:
    """A word wider than a line (a path), cut after its slashes (anywhere in a stretch wider than a line): the first
    piece at most `first` mm ('' where not even its first stretch fits there), the others at most `width`."""
    parts = []
    for x in re.findall(r"[^/]*/|[^/]+", text):
        parts += list(x) if _width([(x, style)]) > width else [x]
    out, cur, room = [], "", first
    for x in parts:
        if _width([(cur + x, style)]) > room:
            out.append(cur)
            cur, room = "", width
        cur += x
    return out + [cur]


def _lines(words: list, width: float) -> list:
    """Words into lines of at most `width` mm, as many words a line as fit (ragged right); a word wider than a line
    is cut (_cut), its first piece filling the line it starts on."""
    space = cu._w(" ", LEGEND_PT)
    out, cur, used = [], [], 0.0
    for word in words:
        w = _width(word)
        if w > width and len(word) == 1:
            (text, style), = word
            pieces = _cut(text, style, width - used - space if cur else width, width)
        else:
            pieces = None
        for i, piece in enumerate(pieces or [word]):
            if pieces is not None:
                if i and cur:                                # a cut: the rest of the word on the next line
                    out.append(cur)
                    cur, used = [], 0.0
                if not piece:
                    continue
                piece = [(piece, style)]
                w = _width(piece)
            if cur and used + space + w > width:
                out.append(cur)
                cur, used = [], 0.0
            used += (space if cur else 0.0) + w
            cur.append(piece)
    if cur:
        out.append(cur)
    return out


def _runs(line: list) -> list:
    """A line as runs of one style each, the space between two words closing the first."""
    runs = []
    for i, word in enumerate(line):
        if i:
            runs[-1][0] += " "
        for t, s in word:
            if runs and runs[-1][1] == s:
                runs[-1][0] += t
            else:
                runs.append([t, s])
    return runs


def _text_pages(text: str) -> Iterator[Figure]:
    """A legend's pages, each drawn when asked for (the caller saves and closes it before asking for the next): 10-pt
    text, ragged right, paragraphs a blank line apart in `text`; a paragraph goes whole to the next page where it does
    not fit on this one (split only where it is longer than a page), each page as tall as its text and at most
    closeup.PAGE_H."""
    lh = LEGEND_PT * LEGEND_LH * 25.4 / 72
    room = cu.PAGE_H - 2 * LEGEND_TB
    pages, cur, used = [], [], 0.0
    for par in text.split("\n\n"):
        lines = _lines(_tokens(par), cu.PAGE_W - 2 * LEGEND_SIDE)
        gap = LEGEND_PAR * lh if cur else 0.0
        if cur and used + gap + len(lines) * lh > room and len(lines) * lh <= room:
            pages.append(cur)
            cur, used, gap = [], 0.0, 0.0
        for i, ln in enumerate(lines):
            g = gap if i == 0 else 0.0
            if cur and used + g + lh > room:
                pages.append(cur)
                cur, used, g = [], 0.0, 0.0
            cur.append((g, ln))
            used += g + lh
    if cur:
        pages.append(cur)
    for page in pages:
        pg = cu._Page(2 * LEGEND_TB + sum(g + lh for g, _ in page))
        y = LEGEND_TB
        for g, ln in page:
            y += g
            x = LEGEND_SIDE
            for t, s in _runs(ln):
                weight, style = _STYLE[s]
                pg.text(x, y + 0.8 * lh, t, fontsize=LEGEND_PT, fontweight=weight, fontstyle=style, va="baseline")
                x += cu._w(t, LEGEND_PT, weight, style)
            y += lh
        yield pg.fig


@contextmanager
def _closing_new_figures():
    """Close the figures opened inside the block and still open when it ends, however it ends. After an error that is
    the page being drawn or saved (a half-drawn page never reaches the code that would close it); the pipeline logs
    the error and goes on, so it would stay open for the rest of the process."""
    before = set(plt.get_fignums())
    try:
        yield
    finally:
        for num in sorted(set(plt.get_fignums()) - before):
            plt.close(num)


def overview(figs: list, out: Path) -> None:
    out = Path(out)
    with _closing_new_figures(), atomic_path(out / "cenfuse_overview.pdf") as tmp, PdfPages(tmp) as pdf:
        for f in figs:
            f.setdefault("out", out)
            first = set()
            for kind, fig in _page_iter(f):
                pdf.savefig(fig)
                if kind not in first:                          # the first page of each kind gets a PNG
                    first.add(kind)
                    name = "overview" if kind == "genome" else "closeups"
                    with atomic_path(out / f"{f['genome'].name}.{name}.png") as png:
                        fig.savefig(png, dpi=300)
                plt.close(fig)                                 # before the next page is drawn
            pp = _plants_page(f)
            if pp is not None:
                pdf.savefig(pp)
                with atomic_path(out / f"{f['genome'].name}.plants.png") as png:
                    pp.savefig(png, dpi=300)
                plt.close(pp)
            for txt in legends(f):
                for fig in _text_pages(txt):
                    pdf.savefig(fig)
                    plt.close(fig)                             # before the next page is drawn
    log.info(f"overview: {out / 'cenfuse_overview.pdf'}")
