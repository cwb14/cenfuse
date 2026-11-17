"""Close-up pages of the overview figure (panel E): one tile per chromosome, two tiles a row, every tile the same width
in bp.

A tile: the chromosome, its window (Mb; kb in windows under 2 Mb), the confidence chip and the verdict
(signals.verdict); under it one track per signal (signals.present), labeled at the left and, in every tile, with the
signal's role in this call (used, points elsewhere, not used: grayed). The centromere is one band through every track,
styled by confidence (high solid black, medium gray, low open with dashed edges). Each track draws its own data and,
as a bar at its foot, the place where that signal alone puts the centromere (signals.Row.places): a colored bar well
off the band is a disagreement cenfuse counted, a gray one a signal not used; a place beyond the window is an arrow at
the tile edge with its position. Satellite arrays are boxes
colored by family role, packed into at most SAT_LANES rows; a long tandem repeat call (seq_evidence long_unit_tandem,
whose repeat is no satellite array) a dark-gray outlined box; a second satellite cluster a dashed box. The centromeric
LTR family's elements are ticks over their density; Hi-C is one trace (the contact with the other centromeres, else
the hub), the score over the chromosome median (1x in gray) on one genome-wide scale, the 99.5th percentile of the
trace as on the genome pages, with a caret where it runs above; CHG methylation is the pooled composition-adjusted
track (down = loss). Other signals have their bar alone. No data is hatched; nothing is drawn past a chromosome end,
which is marked. Tracks run in the order of `keys` (signals.present: signals.TRACK_ORDER).

Window (window): wmin = min(6 Mb, median chromosome length / 2), far = 2 x wmin. One width W for every tile: at least
wmin, and 1.1 x the widest span of a call with the places of the drawn signals used or pointing elsewhere within `far`
of it (and the second satellite cluster, which the verdict names); rounded up to 0.5 Mb (50 kb below 2 Mb), at most the
longest chromosome. Each window is centered on the call and shifted only to hold those places."""
from __future__ import annotations

import math
from collections.abc import Iterator
from dataclasses import dataclass
from functools import lru_cache

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection, PolyCollection
from matplotlib.figure import Figure
from matplotlib.font_manager import FontProperties
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle
from matplotlib.textpath import TextToPath

from . import plot as P                 # rcParams, _track, _meth_track, _nan_runs, ROLE_COLOR, OTHER, HIC, METH, TE
from . import signals as sg

MM = 1 / 25.4                           # layout: inches per millimeter
PAGE_W = 180.0                          # layout: page width, mm (two columns of a journal page)
PAGE_H = 240.0                          # layout: the most a page may be tall, mm (within A4 with margins)
FS = 6.5                                # layout: body text, pt
FSS = 6.0                               # layout: smallest text, pt
WIN_MAX = 6_000_000                     # layout: the narrowest window where chromosomes are over 12 Mb: a plant
                                        # centromere (0.5-3 Mb) with its flanks; smaller genomes take half the median
                                        # chromosome length
SAT_LANES = 6                           # layout: most rows of satellite arrays in a track (overlaps go to the emptiest)
CONF_FACE = {"high": "black", "medium": "#8C8C8C", "low": "white"}
CONF_EDGE = {"high": "black", "medium": "#8C8C8C", "low": "black"}
BAND_FILL = {"high": (0, 0, 0, 0.13), "medium": (0.55, 0.55, 0.55, 0.13), "low": (1, 1, 1, 0)}
BAND_EDGE = {"high": ("black", "-"), "medium": ("#8C8C8C", "-"), "low": ("black", (0, (2, 1.5)))}
# a role's glyph, in the close-ups, panel B and both keys (its ink: signals.ROLE_INK)
ROLE_GLYPH = {"used": dict(marker="o", ms=3.0, color=sg.ROLE_INK["used"], mew=0),
              "elsewhere": dict(marker="x", ms=3.0, color=sg.ROLE_INK["elsewhere"], mew=0.9),
              "not": dict(marker="_", ms=3.6, color=sg.ROLE_INK["not"], mew=1.0)}
FADE = {"trace": 0.35, "sat": 0.45}     # alpha of a signal's data where it is not used
LONG_TANDEM, END, NODATA, HATCH, TICK = "#505050", "#606060", "#C8C8C8", "//////", "#808080"
SAT_NAME = {"centromeric": "centromeric satellite", "cen_associated": "near-centromere satellite",
            "interstitial": "interstitial satellite", "subtelomeric": "subtelomeric satellite",
            "rDNA": "rDNA satellite"}
TRACE = {"hic": "S", "hub": "denovo_S"}  # the Hi-C score the Hi-C trace draws: hic's, else the hub's
NAN_RUN = {"meth": 4}                    # missing bins hatched from this many in a row (Hi-C: one bin is a visible gap
                                         # here; CHG: 50-kb bins, as on the genome pages)
H_TRACK = {"cap": 1.5, "trace": 4.8, "te": 4.0, "meth": 4.0}               # track heights, mm (else H_BAR)
H_BAR = 2.3                              # layout: a track with its bar alone, mm
LINE_H = 2.1                             # layout: one line of 6-pt text at linespacing 0.95, mm
LABEL_GAP = 0.5                          # layout: least space between two tracks' labels, mm
LANE = 0.7                               # layout: one row of satellite arrays, mm
ALT_ROW = 0.8                            # layout: the row of the second satellite cluster, mm
FOOT_BAR = 1.0                           # layout: the zone at a track's foot that holds its bar, mm
CARET = 0.3                              # layout: the Hi-C trace's top (its clip level) below the track's top, mm
TRACK_GAP = 0.5                          # layout: space between tracks, mm
LABEL_MAX = 22.0                         # layout: track labels wrap to this width, mm
COL_GAP = 3.0                            # layout: space between the two columns of tiles, mm
TOP = 2.5                                # layout: page top to the first tile, mm
HEAD_LINE = 2.4                          # layout: tile top to the chromosome name's baseline, mm
VERDICT_LINE = 2.6                       # layout: one line of the verdict, mm
FOOT = 2.0                               # layout: under the tracks, the ticks and the space to the next tile row, mm
SCALE_H = 4.0                            # layout: last tile row to the key, holding the scale bar, mm
KEY_ROW = 3.35                           # layout: one key row at 6 pt, mm (measured: 3.34 a row + 0.64)
KEY_PAD = 1.3                            # layout: margin under the key, mm
CHIP_HALF = 0.75                         # layout: half the confidence chip's width (4 pt), mm
MIN_NAME = 12                            # layout: characters a shortened name keeps before the window's span is dropped


# ------------------------------------------------------------------------------------------------ helpers
_T2P = TextToPath()


@lru_cache(maxsize=None)
def _prop(size: float, weight: str, style: str) -> FontProperties:
    return FontProperties(size=size, weight=weight, style=style)


def _w(s: str, size: float = FSS, weight: str = "normal", style: str = "normal") -> float:
    """Advance width (mm) of the widest line of a text: TextToPath's metrics, a thousandth of plot._text_mm's time."""
    prop = _prop(size, weight, style)
    return max(_T2P.get_text_width_height_descent(x, prop, ismath=False)[0] for x in s.split("\n")) / 72 * 25.4


def _wrap(s: str, width: float, size: float = FSS, **style) -> list[str]:
    """Word wrap to `width` mm: as few lines as a greedy wrap needs, two lines split where the longer is shortest."""
    out, cur = [], ""
    for word in s.split(" "):
        nxt = f"{cur} {word}" if cur else word
        if cur and _w(nxt, size, **style) > width:
            out.append(cur)
            cur = word
        else:
            cur = nxt
    out.append(cur)
    if len(out) == 2:
        words = s.split(" ")
        out = min(([" ".join(words[:i]), " ".join(words[i:])] for i in range(1, len(words))),
                  key=lambda x: max(_w(y, size, **style) for y in x))
    return out


def _shorten(s: str, width: float, size: float, weight: str = "normal") -> str:
    """s cut in the middle with an ellipsis to fit `width` mm (both ends kept: accessions differ at either end)."""
    if _w(s, size, weight) <= width:
        return s
    for keep in range(len(s) - 1, 1, -1):
        x = s[:(keep + 1) // 2] + "\u2026" + s[len(s) - keep // 2:]
        if _w(x, size, weight) <= width:
            return x
    return s[:1] + "\u2026"


def _header(name: str, span: str, room: float) -> tuple[str, str]:
    """A tile's chromosome name and window span within `room` mm: the name shortened (_shorten) to leave room for the
    span; where that would keep fewer than MIN_NAME characters, the span dropped and the name shortened to the room."""
    sw = _w(span, FSS) + 1.5
    if _w(name, FS, "bold") + sw <= room:
        return name, span
    short = _shorten(name, room - sw, FS, "bold")
    if len(short) - 1 >= min(MIN_NAME, len(name)):
        return short, span
    return _shorten(name, room, FS, "bold"), ""


def _mid(iv) -> float:
    return (iv[0] + iv[1]) / 2


def _overlaps(iv, lo, hi) -> bool:
    return iv[1] >= lo and iv[0] <= hi


def _widen(s, e, min_w, L):
    """An interval at least min_w wide (sub-pixel features stay visible), kept inside the chromosome."""
    e = min(max(e, s + min_w), L)
    return max(0, min(s, e - min_w)), e


def _step(W: float) -> float:
    """Tick spacing and scale bar: 1 Mb (100 kb in windows under 2 Mb), coarser or finer for 2 to 20 ticks a tile."""
    steps = [1e4, 2e4, 5e4, 1e5, 2e5, 5e5, 1e6, 2e6, 5e6, 1e7, 2e7, 5e7, 1e8]
    i = steps.index(1e6 if W >= 2e6 else 1e5)
    while W / steps[i] > 20 and i < len(steps) - 1:
        i += 1
    while W / steps[i] < 2 and i > 0:
        i -= 1
    return steps[i]


def _length(bp: float) -> str:
    return f"{bp / 1e6:g} Mb" if bp >= 1e6 else f"{bp / 1e3:g} kb"


def _pos(x: float, kb: bool) -> str:
    return f"{x / 1e3:,.0f} kb" if kb else f"{x / 1e6:.1f} Mb"


def _span(lo: float, hi: float, kb: bool) -> str:
    return f"{lo / 1e3:,.0f}\u2013{hi / 1e3:,.0f} kb" if kb else f"{lo / 1e6:.1f}\u2013{hi / 1e6:.1f} Mb"


def _mids(tel: dict, chrom: str) -> np.ndarray:
    """Midpoints of the centromeric LTR family's elements on one chromosome, those teline counts in its density: the
    signature's elements where there is a signature, else the marker clade's."""
    sig, marker = tel.get("signature"), tel.get("marker")
    if sig is not None:
        m = (sig.cores.get(chrom) or {}).get("mids")
    elif marker:
        m = (tel.get("mids") or {}).get(marker[0], {}).get(chrom)
    else:
        m = None
    return np.sort(np.asarray(m, np.int64)) if m is not None else np.zeros(0, np.int64)


def _pack(s, e, min_w, gap, max_lanes):
    """Rows for intervals (IGV expanded): each into the first row it fits with `gap` to spare, at most max_lanes rows
    (then the row that ends first). Returns the rows and their number."""
    lanes = np.zeros(len(s), int)
    ends = []
    for i in np.argsort(s, kind="stable"):
        de = max(e[i], s[i] + min_w)
        for k, le in enumerate(ends):
            if s[i] >= le + gap:
                lanes[i], ends[k] = k, de
                break
        else:
            if len(ends) < max_lanes:
                ends.append(de)
                lanes[i] = len(ends) - 1
            else:
                k = int(np.argmin(ends))
                lanes[i], ends[k] = k, max(ends[k], de)
    return lanes, max(len(ends), 1)


# ------------------------------------------------------------------------------------------------ window
def _held(r: sg.Row, keys, far) -> list:
    """What r's window must hold: the call, the places of the drawn signals used or pointing elsewhere within `far` of
    it, and the second satellite cluster within `far` (none without a call)."""
    if r.call is None:
        return []
    out = [r.call] + [p for k in keys if (p := r.places.get(k)) and r.roles.get(k) in ("used", "elsewhere")
                      and sg._gap(p, r.call) <= far]
    if "sat" in keys and r.alt and sg._gap(r.alt, r.call) <= far:
        out.append(r.alt)
    return out


def window(rows: list[sg.Row], keys: list[str], chrom_lens: list[int]) -> tuple[float, list[tuple[float, float]]]:
    """The shared window width W and each chromosome's window (lo, hi), by the rule in the module docstring. A
    chromosome without a call is centered on its first drawn signal's place, else on its middle."""
    wmin = min(WIN_MAX, float(np.median(chrom_lens)) / 2)
    far = 2 * wmin
    held = [_held(r, keys, far) for r in rows]
    need = max((max(e for _, e in h) - min(s for s, _ in h) for h in held if h), default=0)
    W = max(wmin, 1.1 * need)
    unit = 50_000 if W < 2e6 else 500_000
    W = float(min(math.ceil(W / unit) * unit, max(chrom_lens)))
    wins = []
    for r, h in zip(rows, held):
        anchor = r.call or next((r.places[k] for k in keys if r.places.get(k)), None)
        lo = (_mid(anchor) if anchor else r.L / 2) - W / 2
        if h:
            ulo, uhi, m = min(s for s, _ in h), max(e for _, e in h), 0.05 * W
            if ulo < lo + m:
                lo = ulo - m
            if uhi > lo + W - m:
                lo = uhi + m - W
        wins.append((lo, lo + W))
    return W, wins


# ------------------------------------------------------------------------------------------------ data
@dataclass
class _Tile:
    row: sg.Row
    lo: float
    hi: float
    ev: str                              # seq_evidence of the call
    arrays: tuple                        # satellite arrays in the window: starts, ends, roles, rows
    traces: dict                         # Hi-C and CHG tracks in the window: key -> (x, y) or None
    mids: np.ndarray                     # LTR elements on the chromosome


def _trace_key(keys) -> str | None:
    """The one Hi-C readout drawn as a trace: the contact with the other centromeres, else the hub."""
    return "hic" if "hic" in keys else "hub" if "hub" in keys else None


def _hic_scale(f: dict, chroms, trace) -> float:
    """The top of a Hi-C trace (signals key `trace`, TRACE), one for the genome: the 99.5th percentile of its values
    over the chromosomes `chroms`, at least 2x (as the genome pages scale it); higher values are clipped."""
    v = [t[1] for c in chroms if (t := P._track(f.get("hic"), c, TRACE[trace])) is not None] if trace else []
    v = np.concatenate(v) if v else np.zeros(0)
    v = v[np.isfinite(v)]
    return max(float(np.percentile(v, 99.5)), 2.0) if len(v) else 2.0


def _clipped(x, v, top):
    """Midpoints of the runs of finite values above `top`."""
    over = np.isfinite(v) & (v > top)
    d = np.diff(np.concatenate([[0], over.astype(int), [0]]))
    return [float(x[s:e].mean()) for s, e in zip(np.flatnonzero(d == 1), np.flatnonzero(d == -1))]


def _prepare(f: dict, rows, keys, W, wins, min_w):
    """Each tile's data in its window, the shared scales and what the key must name."""
    trace = _trace_key(keys)
    hsc = _hic_scale(f, [r.chrom for r in rows], trace)
    ev = {}
    if "seq_evidence" in f["calls"]:
        ev = dict(zip(f["calls"].chrom.astype(str), f["calls"].seq_evidence.astype(str)))
    A = f["seq"].arrays
    tel = f.get("tel") or {}
    tiles, nl = [], 1
    for r, (lo, hi) in zip(rows, wins):
        a = A[(A.chrom == r.chrom) & (A.end > lo) & (A.start < hi)] if len(A) else A
        s, e = a.start.to_numpy(np.int64), a.end.to_numpy(np.int64)
        lanes, k = _pack(s, e, min_w, 0.6 * min_w, SAT_LANES)
        nl = max(nl, k if len(a) else 1)
        traces = {}
        for key in keys:
            if key == trace:
                traces[key] = P._track(f.get("hic"), r.chrom, TRACE[key])
            elif key == "meth":
                traces[key] = P._meth_track(f.get("meth"), r.chrom)
        tiles.append(_Tile(r, lo, hi, ev.get(r.chrom, "."), (s, e, a.role.to_numpy(), lanes), traces,
                           _mids(tel, r.chrom) if "te" in keys else np.zeros(0, np.int64)))
    te_max = 1e-9
    te_bin = W / 32
    seen = dict(arrow=False, end=False, hatch=False, alt=False, long=False, roles=set(), other=False, clip=False)
    for t in tiles:
        r, lo, hi = t.row, t.lo, t.hi
        c0, c1 = max(lo, 0), min(hi, r.L)
        seen["end"] |= lo < 0 or hi > r.L
        seen["long"] |= (t.ev == "long_unit_tandem" and "sat" in keys and bool(p := r.places.get("sat"))
                         and _overlaps(p, lo, hi))
        if "sat" in keys and r.alt:
            seen["alt" if _overlaps(r.alt, lo, hi) else "arrow"] = True
        seen["arrow"] |= any(p and not _overlaps(p, lo, hi) for k in keys if (p := r.places.get(k)))
        roles = t.arrays[2]
        seen["roles"] |= {x for x in roles if x in P.ROLE_COLOR}
        seen["other"] |= any(x not in P.ROLE_COLOR for x in roles)
        for key, tr in t.traces.items():
            if tr is None:
                seen["hatch"] = True
                continue
            x, y = tr
            k = (x > lo) & (x < hi)
            if key == trace:
                seen["clip"] |= bool(_clipped(x[k], y[k], hsc))
            seen["hatch"] |= any(e_ > c0 and s_ < c1 for s_, e_ in P._nan_runs(x, y, NAN_RUN.get(key, 1)))
        if len(t.mids) and c1 > c0:
            cnt = _density(t.mids, c0, c1, te_bin)[1]
            te_max = max(te_max, float(cnt.max()) if len(cnt) else 0)
    msc = 0.1
    if "meth" in keys and f.get("meth") is not None:
        v = [np.abs(t.traces["meth"][1]) for t in tiles if t.traces.get("meth") is not None]
        v = np.concatenate(v) if v else np.zeros(0)
        msc = max(float(np.nanpercentile(v, 99.5)), 0.02) if np.isfinite(v).any() else 0.1
    return tiles, nl, dict(hic=hsc, meth=msc, te=te_max, te_bin=te_bin), seen


def _density(mids, c0, c1, b):
    """LTR elements per Mb in bins of b bp from c0 to c1 (a last bin under half a bin joins the one before)."""
    edges = np.arange(c0, c1, b, dtype=float)
    if len(edges) > 1 and c1 - edges[-1] < b / 2:
        edges = edges[:-1]
    edges = np.append(edges, c1)
    return edges, np.histogram(mids, bins=edges)[0] / (np.diff(edges) / 1e6)


# ------------------------------------------------------------------------------------------------ layout
class _Page:
    """A page in mm from its top-left corner."""

    def __init__(self, h: float):
        self.h = h
        self.fig = plt.figure(figsize=(PAGE_W * MM, h * MM))

    def ax(self, x, y, w, h):
        return self.fig.add_axes([x / PAGE_W, 1 - (y + h) / self.h, w / PAGE_W, h / self.h])

    def text(self, x, y, s, **kw):
        return self.fig.text(x / PAGE_W, 1 - y / self.h, s, **kw)

    def mark(self, x, y, **kw):
        self.fig.add_artist(Line2D([x / PAGE_W], [1 - y / self.h], transform=self.fig.transFigure, ls="", **kw))


def role_key(seen) -> list:
    """Key entries for the roles in `seen`, in the order used, points elsewhere, not used."""
    return [Line2D([], [], ls="", label=sg.ROLE_WORDS[k], **ROLE_GLYPH[k]) for k in ROLE_GLYPH if k in seen]


def _role_lines(text: str) -> str:
    return text.replace(" (", "\n(")                # 'not used (other plants)' on two lines


@dataclass
class _Geo:
    gut: float                           # track labels, left of the first column
    rolew: float                         # role glyph and word, every tile
    tile_w: float
    dw: float                            # data width
    head: float                          # header height
    heights: dict                        # track -> height
    tracks_h: float
    tile_h: float
    labels: dict                         # track -> wrapped label
    verdicts: dict                       # chrom -> wrapped verdict


def _geometry(rows, keys, clade, summary, nl, alt) -> _Geo:
    trace = _trace_key(keys)
    labels = {k: "\n".join(_wrap(sg.SIGNALS[k].label(clade), LABEL_MAX)) for k in keys}
    gut = max([_w(x) for x in labels.values()] + [_w("centromere")]) + 2.0
    words = {_role_lines(sg.role_label(k, r.roles.get(k, "not"), summary)) for r in rows for k in keys}
    rolew = 2.8 + max(_w(x) for x in words) + 1.5
    tile_w = (PAGE_W - 1.0 - gut - COL_GAP) / 2
    dw = tile_w - rolew
    verdicts = {r.chrom: "\n".join(_wrap(r.verdict, tile_w, style="italic")) for r in rows}
    nv = max(v.count("\n") + 1 for v in verdicts.values())
    head = HEAD_LINE + 0.5 + VERDICT_LINE * nv + 0.4
    heights, lines = {"cap": H_TRACK["cap"]}, {"cap": 1}
    for k in keys:
        heights[k] = (LANE * nl + (ALT_ROW if alt else 0) + FOOT_BAR + 0.2 if k == "sat" else
                      H_TRACK.get("trace" if k == trace else k, H_BAR))
        lines[k] = max(labels[k].count("\n") + 1, max(_role_lines(sg.role_label(k, r.roles.get(k, "not"), summary))
                                                      .count("\n") + 1 for r in rows))
    names = list(heights)
    for a, b in zip(names, names[1:]):                       # labels and roles, centered on their tracks, never touch
        need = (lines[a] + lines[b]) / 2 * LINE_H + LABEL_GAP - TRACK_GAP - heights[a] / 2
        heights[b] = max(heights[b], 2 * need)
    tracks_h = sum(heights.values()) + TRACK_GAP * (len(heights) - 1)
    return _Geo(gut, rolew, tile_w, dw, head, heights, tracks_h, head + tracks_h + FOOT, labels, verdicts)


def _key_items(rows, keys, clade, scales, seen) -> list:
    """The key, in groups that a column-by-column fill keeps together where it can: roles, the centromere, the places,
    then each track's marks."""
    items = role_key({r.roles.get(k, "not") for r in rows for k in keys})
    if any(r.call for r in rows):
        items.append(Patch(fc=BAND_FILL["high"], ec="black", lw=0.6, label="centromere, through every track"))
    items += [Patch(fc=CONF_FACE[c], ec=CONF_EDGE[c], lw=0.7, label=f"centromere, {c} confidence")
              for c in ("high", "medium", "low") if any(r.call and r.conf == c for r in rows)]
    items.append(Line2D([], [], color="#555555", lw=2.2, solid_capstyle="butt",
                        label="where this signal alone puts it"))
    if seen["arrow"]:
        items.append(Line2D([], [], ls="", marker=r"$\rightarrow$", ms=6, mew=0, color="#555555",
                            label="its place, beyond the close-up"))
    if _trace_key(keys):
        items.append(Line2D([], [], color=P.HIC, lw=0.8,
                            label=f"Hi-C, \u00d7 chromosome median (gray 1\u00d7, top {scales['hic']:.1f}\u00d7)"))
    if seen["clip"]:
        items.append(Line2D([], [], ls="", marker="^", ms=2.6, mew=0, color=P.HIC, label="Hi-C above the top"))
    if "meth" in keys:
        items.append(Line2D([], [], color=P.METH, lw=0.7, label="CHG methylation, down = loss (gray 0)"))
    if "te" in keys:
        items += [Line2D([], [], ls="", marker="|", ms=4, mew=0.8, color=P.TE, label=f"{clade} element"),
                  Patch(fc=P.TE, alpha=0.35, ec="none", label=f"{clade} elements per Mb (dashed: chromosome)")]
    items += [Patch(fc=P.ROLE_COLOR[x], ec="none", label=SAT_NAME[x]) for x in P.ROLE_COLOR if x in seen["roles"]]
    if seen["other"]:
        items.append(Patch(fc=P.OTHER, ec="none", label="other tandem repeat"))
    if seen["long"]:
        items.append(Patch(fc="none", ec=LONG_TANDEM, lw=0.8, label="long tandem repeat"))
    if seen["alt"]:
        items.append(Patch(fc="none", ec=P.ROLE_COLOR["centromeric"], lw=0.8, ls=(0, (2, 1.2)),
                           label="second satellite cluster"))
    if seen["end"]:
        items.append(Line2D([], [], color=END, lw=0.6, label="chromosome end"))
    if seen["hatch"]:
        items.append(Patch(fc="none", ec=NODATA, hatch=HATCH, lw=0, label="no data"))
    return items


def _key_columns(items) -> int:
    """The most key columns whose widths fit the page, filled as matplotlib fills them: column by column, the first
    len % ncol columns one entry longer."""
    widths = [_w(x.get_label()) + 4.0 for x in items]        # handle and its pad: 1.8 em at 6 pt
    for ncol in range(min(6, len(items)), 0, -1):
        rows, large = divmod(len(items), ncol)
        sizes = [rows + 1] * large + [rows] * (ncol - large)
        starts = np.cumsum([0] + sizes[:-1])
        used = sum(max(widths[a:a + n]) for a, n in zip(starts, sizes)) + 2.2 * (ncol - 1)   # column spacing 1 em
        if used <= PAGE_W - 1.5 - 1.0 - 1.8:                 # anchor, right margin, border padding
            return ncol
    return 1


# ------------------------------------------------------------------------------------------------ drawing
def _arrow(ax, t: _Tile, y, p, color, gid, kb):
    """A place beyond the window: an arrow at the tile edge on its side, with its position."""
    W = t.hi - t.lo
    left = p[1] < t.lo
    xa, xt = (t.lo, t.lo + 0.06 * W) if left else (t.hi, t.hi - 0.06 * W)
    ax.annotate("", xy=(xa, y), xytext=(xt, y), annotation_clip=False, zorder=7,
                arrowprops=dict(arrowstyle="-|>", lw=0.8, color=color, mutation_scale=5, shrinkA=0, shrinkB=0))
    ax.text(xt + (0.012 if left else -0.012) * W, y, _pos(_mid(p), kb), ha="left" if left else "right", va="center",
            fontsize=FSS, color=color, zorder=7, gid=gid, bbox=dict(fc="white", ec="none", pad=0.3, alpha=0.85))


def _tile(pg: _Page, t: _Tile, x0, y0, geo: _Geo, keys, scales, nl, kb, summary, label: bool):
    """One chromosome's tile, its top-left corner at (x0, y0) mm."""
    r, lo, hi, L = t.row, t.lo, t.hi, t.row.L
    trace = _trace_key(keys)
    c0, c1 = max(lo, 0), min(hi, L)
    xd = x0 + geo.rolew                                      # data area
    bp_mm = (hi - lo) / geo.dw
    min_w = 0.3 * bp_mm
    # header: name, window, confidence; the verdict
    room = geo.tile_w - (_w(r.conf, FS) + 2.0 + CHIP_HALF + 1.0 if r.call else 0.0)     # left of the chip
    name, span = _header(r.chrom, _span(c0, c1, kb), room)
    pg.text(x0, y0 + HEAD_LINE, name, fontsize=FS, fontweight="bold", va="baseline", gid=f"name:{r.chrom}")
    if span:
        pg.text(x0 + _w(name, FS, "bold") + 1.5, y0 + HEAD_LINE, span, fontsize=FSS, va="baseline", color="#666666",
                gid=f"window:{r.chrom}")
    if r.call:
        pg.text(x0 + geo.tile_w, y0 + HEAD_LINE, r.conf, fontsize=FS, va="baseline", ha="right")
        pg.mark(x0 + geo.tile_w - _w(r.conf, FS) - 2.0, y0 + HEAD_LINE - 0.85, marker="s", ms=4.0,
                mfc=CONF_FACE.get(r.conf, "white"), mec=CONF_EDGE.get(r.conf, "black"), mew=0.7,
                gid=f"chip:{r.chrom}")
    pg.text(x0, y0 + HEAD_LINE + 0.6, geo.verdicts[r.chrom], fontsize=FSS, va="top", style="italic",
            linespacing=1.05, gid=f"verdict:{r.chrom}")
    ytr = y0 + geo.head
    foot_y = geo.tracks_h + 0.3
    ax = pg.ax(xd, ytr, geo.dw, foot_y + 1.0)
    ax.set_xlim(lo, hi)
    ax.set_ylim(foot_y + 1.0, 0)                             # mm from the top of the tracks
    ax.axis("off")
    ax.set_gid(f"tile:{r.chrom}")
    # track offsets, labels and roles
    offs, y = {}, 0.0
    for k, h in geo.heights.items():
        offs[k] = (y, y + h)
        mid = ytr + y + h / 2
        if label:
            pg.text(geo.gut - 1.0, mid, "centromere" if k == "cap" else geo.labels[k], ha="right", va="center",
                    fontsize=FSS, linespacing=0.95, gid="cap" if k == "cap" else f"track:{k}")
        if k != "cap":
            role = r.roles.get(k, "not")
            pg.mark(x0 + 1.2, mid, **ROLE_GLYPH[role])
            pg.text(x0 + 2.8, mid, _role_lines(sg.role_label(k, role, summary)), va="center", fontsize=FSS,
                    linespacing=0.95, color=sg.ROLE_INK[role], gid=f"role:{r.chrom}:{k}")
        y += h + TRACK_GAP
    bottom = geo.tracks_h
    # the centromere: a mark in the confidence style, and one band through every track
    if r.call:
        s, e = _widen(*r.call, min_w, L)
        cap = offs["cap"]
        ax.add_patch(Rectangle((s, cap[0] + 0.15), e - s, cap[1] - cap[0] - 0.3, fc=CONF_FACE.get(r.conf, "white"),
                               ec=CONF_EDGE.get(r.conf, "black"), lw=0.7, zorder=5, gid=f"band:{r.chrom}"))
        ax.add_patch(Rectangle((s, cap[1]), e - s, bottom - cap[1], fc=BAND_FILL.get(r.conf, BAND_FILL["low"]),
                               ec="none", lw=0, zorder=0, gid=f"band:{r.chrom}"))
        col, ls = BAND_EDGE.get(r.conf, BAND_EDGE["low"])
        for x in (s, e):
            ax.plot([x, x], [cap[1], bottom], color=col, lw=0.45, ls=ls, zorder=0.5)
    for k in keys:
        yt, yb = offs[k]
        role = r.roles.get(k, "not")
        hue = sg.SIGNALS[k].color
        base, avail = yb - FOOT_BAR, yb - FOOT_BAR - yt
        if k == "sat":
            s, e, roles, lanes = t.arrays
            if len(s):
                verts, cols = [], []
                for a, b, lane, rl in zip(s, e, lanes, roles):
                    a, b = _widen(a, b, min_w, L)
                    verts.append([(a, yt + lane * LANE + 0.08), (b, yt + lane * LANE + 0.08),
                                  (b, yt + lane * LANE + 0.62), (a, yt + lane * LANE + 0.62)])
                    cols.append(P.ROLE_COLOR.get(rl, P.OTHER))
                o = np.argsort([c != P.OTHER for c in cols], kind="stable")       # gray context under the colors
                ax.add_collection(PolyCollection([verts[j] for j in o], facecolors=[cols[j] for j in o],
                                                 edgecolors="none", alpha=1.0 if role != "not" else FADE["sat"],
                                                 zorder=3))
            p = r.places.get("sat")
            if t.ev == "long_unit_tandem" and p and _overlaps(p, lo, hi):
                a, b = _widen(max(p[0], c0), min(p[1], c1), min_w, L)
                ax.add_patch(Rectangle((a, yt + 0.04), b - a, LANE * nl - 0.08, fill=False, ec=LONG_TANDEM, lw=0.8,
                                       zorder=4, gid=f"long:{r.chrom}"))
            if r.alt:
                ya = yt + LANE * nl + 0.1
                if _overlaps(r.alt, lo, hi):
                    a, b = _widen(max(r.alt[0], c0), min(r.alt[1], c1), min_w, L)
                    ax.add_patch(Rectangle((a, ya), b - a, ALT_ROW - 0.2, fill=False, ec=P.ROLE_COLOR["centromeric"],
                                           lw=0.8, ls=(0, (2, 1.2)), zorder=4, gid=f"alt:{r.chrom}"))
                else:
                    _arrow(ax, t, ya + (ALT_ROW - 0.2) / 2, r.alt, P.ROLE_COLOR["centromeric"], f"arrow:{r.chrom}:alt",
                           kb)
        elif k == "te":
            alpha = 1.0 if role != "not" else FADE["trace"]
            if len(t.mids) and c1 > c0:
                edges, cnt = _density(t.mids, c0, c1, scales["te_bin"])
                top = base - np.minimum(cnt / scales["te"], 1) * 0.55 * avail
                ax.fill_between(edges, base, np.append(top, top[-1]), step="post", color=hue, alpha=0.35 * alpha,
                                lw=0, zorder=2)
                avg = len(t.mids) / (L / 1e6)
                ya = base - min(avg / scales["te"], 1) * 0.55 * avail
                ax.plot([c0, c1], [ya, ya], color=hue, lw=0.5, ls=(0, (2, 1.2)), alpha=alpha, zorder=2.5)
                m = t.mids[(t.mids >= lo) & (t.mids <= hi)]
                ax.add_collection(LineCollection([[(x, yt), (x, yt + 0.38 * avail)] for x in m], colors=hue, lw=0.6,
                                                 alpha=alpha, zorder=3))
        elif k == trace or k == "meth":
            alpha = 1.0 if role != "not" else FADE["trace"]
            tr = t.traces.get(k)
            if k == "meth":
                msc, mid = scales["meth"], yt + avail / 2
                to_y = lambda v: mid - np.clip(v, -msc, msc) / msc * avail / 2      # noqa: E731
                ref, color, lw, run = 0.0, P.METH, 0.55, NAN_RUN["meth"]
            else:
                hsc, span = scales["hic"], avail - CARET                          # the top: the clip level
                to_y = lambda v: base - np.clip(v, 0, hsc) / hsc * span           # noqa: E731
                ref, color, lw, run = 1.0, hue, 0.7, NAN_RUN.get(k, 1)
            if c1 > c0:
                ax.plot([c0, c1], [float(to_y(ref))] * 2, color=NODATA, lw=0.4, zorder=1)
            if tr is None:
                if c1 > c0:
                    ax.add_patch(Rectangle((c0, yt), c1 - c0, avail, fill=False, hatch=HATCH, ec=NODATA, lw=0,
                                           zorder=1))
            else:
                x, v = tr
                step = x[1] - x[0] if len(x) > 1 else 0
                keep = (x > lo - step) & (x < hi + step)
                ax.plot(x[keep], to_y(v[keep]), color=color, lw=lw, alpha=alpha, zorder=3, solid_joinstyle="round")
                if k == trace:                                           # where the trace runs above the top
                    win = (x > lo) & (x < hi)                            # as _prepare decides the key
                    xc = _clipped(x[win], v[win], hsc)
                    if xc:
                        ax.plot(xc, [yt + CARET] * len(xc), ls="", marker="^", ms=2.6, mew=0, color=color,
                                alpha=alpha, zorder=4, gid=f"clip:{r.chrom}")
                for s, e in P._nan_runs(x, v, run):
                    s, e = max(s, c0), min(e, c1)
                    if e > s:
                        ax.add_patch(Rectangle((s, yt), e - s, avail, fill=False, hatch=HATCH, ec=NODATA, lw=0,
                                               zorder=1))
        # where this signal alone puts the centromere: a bar at the track's foot, or an arrow at the tile edge
        p = r.places.get(k)
        if p:
            color = sg.GRAY if role == "not" else hue
            yb_bar = yb - FOOT_BAR / 2
            if _overlaps(p, lo, hi):
                a, b = _widen(max(p[0], c0), min(p[1], c1), min_w, L)
                ax.plot([a, b], [yb_bar, yb_bar], color=color, lw=2.2, solid_capstyle="butt", zorder=6,
                        gid=f"place:{r.chrom}:{k}")
            else:
                _arrow(ax, t, yb - LINE_H / 2, p, color, f"arrow:{r.chrom}:{k}", kb)      # its label inside the track
    # chromosome ends inside the window; ticks along the foot
    for edge in (0, L):
        if lo < edge < hi:
            ax.plot([edge, edge], [0, bottom], color=END, lw=0.6, zorder=2)
    if c1 > c0:
        step = _step(hi - lo)
        ax.plot([c0, c1], [foot_y, foot_y], color=TICK, lw=0.4)
        ticks = np.arange(math.ceil(c0 / step) * step, c1 + 1, step)
        ax.add_collection(LineCollection([[(x, foot_y), (x, foot_y + 0.7)] for x in ticks], colors=TICK, lw=0.4))


def pages(f: dict, rows: list[sg.Row], keys: list[str], letter: str = "E") -> Iterator[Figure]:
    """The close-up pages of one genome, each drawn when asked for (the caller saves and closes a page before asking
    for the next): two tiles a row, as many rows as keep a page within 240 mm, the key and a scale bar on every page.
    `f` is the dict plot.overview receives (its 'summary', the genome's cenfuse_summary.tsv row, may be missing);
    `rows` come from signals.rows and `keys` from signals.present."""
    if not rows:
        return
    keys = list(keys)
    summary = f.get("summary")
    clade = sg.clade_of(f["calls"], summary)
    W, wins = window(rows, keys, [r.L for r in rows])
    kb = W < 2e6
    # the data width depends only on the labels, so it is known before the tiles are packed
    probe = _geometry(rows, keys, clade, summary, 1, False)
    min_w = 0.3 * W / probe.dw
    tiles, nl, scales, seen = _prepare(f, rows, keys, W, wins, min_w)
    geo = _geometry(rows, keys, clade, summary, nl, seen["alt"])
    items = _key_items(rows, keys, clade, scales, seen)
    ncol = _key_columns(items)
    key_h = KEY_ROW * math.ceil(len(items) / ncol) + KEY_PAD
    per_page = 2 * max(1, int((PAGE_H - TOP - SCALE_H - key_h) // geo.tile_h))
    step = _step(W)
    for i0 in range(0, len(tiles), per_page):
        chunk = tiles[i0:i0 + per_page]
        nrow = math.ceil(len(chunk) / 2)
        pg = _Page(TOP + nrow * geo.tile_h + SCALE_H + key_h)
        pg.text(1.0, 0.5, letter, fontsize=9, fontweight="bold", va="top", ha="left")
        for j, t in enumerate(chunk):
            ri, ci = divmod(j, 2)
            _tile(pg, t, geo.gut + ci * (geo.tile_w + COL_GAP), TOP + ri * geo.tile_h, geo, keys, scales, nl, kb,
                  summary, label=ci == 0)
        yb = TOP + nrow * geo.tile_h + 0.6
        x1 = geo.gut + 2 * geo.tile_w + COL_GAP
        x0 = x1 - step / W * geo.dw
        pg.fig.add_artist(Line2D([x0 / PAGE_W, x1 / PAGE_W], [1 - yb / pg.h] * 2, color="black", lw=1.0,
                                 transform=pg.fig.transFigure))
        pg.text(x0 - 1.0, yb, f"{_length(step)} (all close-ups)", ha="right", va="center", fontsize=FSS, gid="scale")
        pg.fig.legend(handles=items, loc="upper left", bbox_to_anchor=(1.5 / PAGE_W, 1 - (yb + SCALE_H - 0.6) / pg.h),
                      ncol=ncol, frameon=False, fontsize=FSS, handlelength=1.4, columnspacing=1.0, handletextpad=0.4,
                      borderaxespad=0)
        yield pg.fig
