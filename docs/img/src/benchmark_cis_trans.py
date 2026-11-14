#!/usr/bin/env python3
"""benchmark_cis_trans: draws docs/img/benchmark_cis_trans.pdf and .png (docs/benchmark.md, cenfuse 0.16.0).

The contact across each bin (cis) and the contact between chromosomes (trans) are two readouts of one Hi-C library; where each
chromosome's cis extreme, read outside its terminal eighths, is met by its trans peak made from the OTHER chromosomes'
extremes, the cis signal marks centromeres that carry no repeat. This script records how the figure was made, legend text
included. It CANNOT rebuild the figure from a clone of the repository alone: it reads the cenfuse 0.16.0 benchmark runs and
tables, which are not in the repository (too large; some truth unpublished).

Three steps, each writing to the work folder ($CENFUSE_FIG_WORK, default: the current directory):
1. benchmark_cis_trans_data.py: contact maps, profiles and scores read with cenfuse 0.16.0's own functions from the stored
   matrices of the release runs -> benchmark_cis_trans_cache_<entry>.pkl
2. benchmark_cis_trans_libs.py: every 0.16.0 run whose Hi-C was tested by trans contacts, with its chance
   -> benchmark_cis_trans_runs16.tsv, and the libraries re-read from their matrices -> benchmark_cis_trans_libs16.tsv
3. this script -> benchmark_cis_trans.pdf (figure + one legend page) and .png (figure only, 300 dpi) in $CENFUSE_FIG_OUT
   (default: docs/img; docs/img's PNG is that PNG, losslessly recompressed), and benchmark_cis_trans_numbers.tsv (every
   number with its source) and benchmark_cis_trans_libraries.tsv in the work folder

Data read, under the benchmark data root $CENFUSE_BENCH_DATA (default: the folder they were made in, as in
benchmark/paths.py), with V12 = v2/v3/v4/v5/v6/v7/v8/v9/v10/v10.1/v10.2/v10.3/v10.4/v11/v12 and V16 = V12/v13/v14/v15/v16:
- V16/runs/rel16/SH/{tomato,aetAL878,cneoG1} and V16/runs/rel16/TNEG/tomato: hic/<g>.cache/*.npz (contact matrix) and
  work/hic/<g>.pkl (bins, members, hub)
- every 0.16.0 release run with Hi-C (V16/bench/jobs_rel16*.tsv): cenfuse_summary.tsv (cis_line) and its stored matrix
- V16/results/docs16_rel16/release16.{runs,chrom}.tsv: calls and outcomes of arms S (FASTA only) and SH (with Hi-C)
- truth: V12/../bench/rel/score_rel.py (the benchmark's frozen truth rules); tomato's CENH3 ChIP and input reads
  (v2/.../v8/data/tomato_scau/slH1706)
- 0.13-era tables (V12/results/cis_trans_anchor*.tsv; the 0.13 figure's per-library chances, V12/figures/fig_cis_trans):
  checked against the 0.16.0 re-read; a library's exact chance is reused only where its matrix, k, n, verdict and sign
  are unchanged. Code: cenfuse 0.16.0 (V16/code/rel16).

usage: CENFUSE_FIG_WORK=<work folder> python docs/img/src/benchmark_cis_trans.py [-v]
"""
import argparse
import datetime
import os
import pickle
import socket
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(os.environ.get("CENFUSE_FIG_WORK") or Path.cwd()).resolve()   # work folder: caches and tables
SRC = Path(__file__).resolve().parent                                      # docs/img/src
OUT = Path(os.environ.get("CENFUSE_FIG_OUT") or SRC.parent)                 # the PDF and the PNG (default: docs/img)
STEM = "benchmark_cis_trans"
DATA = Path(os.environ.get("CENFUSE_BENCH_DATA") or "/data2/chris/centromere_benchmark")   # as benchmark/paths.py
V = DATA / "v2/v3/v4/v5/v6/v7/v8/v9/v10/v10.1/v10.2/v10.3/v10.4/v11/v12"
V11 = V.parent
V16 = V / "v13/v14/v15/v16"
RES = V16 / "results/docs16_rel16"          # the 0.16.0 release's score tables
RES13 = V / "results"                       # 0.13-era data-level tables: compared with the 0.16.0 re-read, never drawn
RUNS = V16 / "runs/rel16"
TOL = 1_000_000                             # 'on CENH3': score_rel's tolerance for tomato, Ae. tauschii, potato (4 x 250-kb bins)
VERBOSE = False

C_CIS, C_TRANS, C_SAT = "#6A3D9A", "#0072B2", "#D55E00"
C_CEN, C_CEN_EDGE = "#C6DFF0", "#7FA7C4"     # CENH3 / CENP-A: light blue band
C_RAND, C_END, C_MEET, C_HATCH = "#A0A0A0", "#ECECEC", "#222222", "#C2C2C2"
HALF = 4                    # hp.HALF_WIN: the crossing box half-width and the zone around an anchor (bins)
END_ZONE = 0.125            # hicfuse.END_ZONE: terminal eighths, not read
FW, FH = 180.0, 236.0       # figure size, mm



def cis_sign(line: str) -> str:
    """'peak' or 'dip' from a summary cis_line, in either wording: 0.13's final 'trans peaks at the cis peak on k/n ...' or the
    pre-release 'peak at the trans peak on k/n ...'."""
    line = str(line)
    if " at the cis " in line:
        return line.split(" at the cis ", 1)[1].split()[0]
    return line.split(" at ")[0]

def log(*a):
    if VERBOSE:
        print(*a, file=sys.stderr, flush=True)


def load(name):
    return pickle.load(open(HERE / f"benchmark_cis_trans_cache_{name}.pkl", "rb"))


def chr_label(gid, c, names):
    return f"chr{names.index(c) + 1}"


# ============================================================================== numbers
class Numbers:
    """Every number in the legends, with its source; the legend text is formatted from here."""
    def __init__(self):
        self.rows = []
        self.d = {}

    def add(self, key, value, source, note=""):
        self.d[key] = value
        self.rows.append(dict(key=key, value=value, source=source, note=note))
        return value

    def __getitem__(self, k):
        return self.d[k]

    def write(self, path):
        pd.DataFrame(self.rows).to_csv(path, sep="\t", index=False)


def fmt_p(p):
    if p >= 0.01:
        return f"{p:.2g}"
    m, e = f"{p:.1e}".split("e")
    return f"{m} × 10$^{{{int(e)}}}$"


def runs16():
    """Every 0.16.0 release run whose cis_line was tested by trans contacts (benchmark_cis_trans_libs.py table)."""
    return pd.read_csv(HERE / "benchmark_cis_trans_runs16.tsv", sep="\t")


def on_truth(G, tr, mode, sign):
    """Chromosomes whose cis extreme (mode 'int' or 'all', sign +1/-1) lies within TOL of the truth: (k, n with truth)."""
    b = G["bin"]
    ext = G["ext"][mode][sign]
    have = [c for c in ext if c in tr]

    def dist(s, e, c):
        return min(max(0, max(s, a) - min(e, z)) for a, z in tr[c])
    k = sum(dist(ext[c] * b, (ext[c] + 1) * b, c) <= TOL for c in have)
    return k, len(have)


def collect_numbers(T, R, A, C, Pt, meta):
    N = Numbers()
    ct = runs16()
    r16 = pd.read_csv(RES / "release16.runs.tsv", sep="\t")
    c16 = pd.read_csv(RES / "release16.chrom.tsv", sep="\t")
    src_ct = "benchmark_cis_trans_runs16.tsv (cis_line of runs/rel16/<arm>/<genome>/cenfuse_summary.tsv)"

    def ctrow(arm, g):
        return ct[(ct.arm == arm) & (ct.genome == g)].iloc[0]
    for arm, g, key in (("SH", "tomato", "tom"), ("SH", "cneoG1", "cneo"), ("SH", "brachy", "brachy"), ("SH", "ncra", "ncra"),
                        ("SH", "vdah", "vdah"), ("SH", "barley", "barley"), ("SH", "cneoA", "cneoA"), ("FULL", "MMU8", "mmu"),
                        ("SH", "spom", "spom"), ("NEG2", "gmX26", "rgm"), ("NEG2", "MMU8", "rmmu"), ("NEG2", "chm13", "rchm"),
                        ("SH", "aetAL878", "aet16"), ("SH", "potato", "pot16"), ("SH", "ts60", "ts60")):
        r = ctrow(arm, g)
        N.add(f"{key}_kn", f"{r.k}/{r.n}", src_ct, f"{arm} {g}: {r.cis_line}")
        N.add(f"{key}_p", float(r.p), src_ct, f"{arm} {g}, Bonferroni x2 signs")
        N.add(f"{key}_sign", r.sign, src_ct, f"{arm} {g}")
        N.add(f"{key}_verdict", r.verdict, src_ct, f"{arm} {g}")
        N.add(f"{key}_acts", int(r.acts), src_ct, f"{arm} {g}: chromosomes the line acts on")
    rnd = ct[ct.arm.isin(["NEG", "NEG2"])]
    N.add("n_rand", len(rnd), src_ct, "randomized libraries (arms NEG, NEG2) where the test ran")
    N.add("n_rand_pass", int((rnd.verdict == "pass").sum()), src_ct, "randomized libraries that pass")
    N.add("n_rand_support", int((rnd.verdict == "support").sum()), src_ct, "randomized libraries at support")
    tomx = ct[(ct.arm == "SH") & ct.genome.str.match(r"^(la\d|ts\d)")]
    N.add("tomx_n", len(tomx), src_ct, "the 13 further tomato accessions (SH arm), new in 0.16")
    for v in ("pass", "support", "fail"):
        N.add(f"tomx_{v}", int((tomx.verdict == v).sum()), src_ct, f"further tomato accessions at {v}")
    sup = tomx[tomx.verdict == "support"]
    N.add("tomx_support_k", f"{sup.k.min()}–{sup.k.max()}", src_ct, "further tomato accessions at support: chromosomes meeting")
    N.add("sver_kn", f"{ctrow('SH', 'sver').k}/{ctrow('SH', 'sver').n}", src_ct, "S. verrucosum (SH)")
    N.add("sver_verdict", ctrow("SH", "sver").verdict, src_ct, "S. verrucosum (SH)")
    for g, key in (("tomato", "tom"), ("cneoG1", "cneo"), ("brachy", "brachy"), ("ncra", "ncra"), ("vdah", "vdah")):
        for arm in ("S", "SH"):
            b = r16[(r16.arm == arm) & (r16.genome == g)].iloc[0]
            N.add(f"{key}_right{arm}", f"{int(b.right)}/{int(b.chromosomes)}", "results/docs16_rel16/release16.runs.tsv", f"{arm} {g}: right calls")
            N.add(f"{key}_confwrong{arm}", int(b.conf_wrong), "results/docs16_rel16/release16.runs.tsv", f"{arm} {g}: confidently wrong")
            N.add(f"{key}_rightmed{arm}", int(b.right_high) + int(b.right_medium), "results/docs16_rel16/release16.runs.tsv",
                  f"{arm} {g}: right at high or medium confidence")
            N.add(f"{key}_P{arm}", float(b.P), "results/docs16_rel16/release16.runs.tsv", f"{arm} {g}: P of the right count against chance")
    cw = c16[(c16.arm == "SH") & (c16.genome == "cneoG1") & (c16.ok == "W")]
    N.add("cneo_wrong_n", len(cw), "results/docs16_rel16/release16.chrom.tsv", "SH cneoG1: wrong calls")
    N.add("cneo_wrong_gap_kb", f"{cw.gap.min() / 1e3:.0f}", "results/docs16_rel16/release16.chrom.tsv", "SH cneoG1: the wrong call's distance to CENP-A (kb)")
    N.add("cneo_wrong_chrom", ",".join(chr_label("cneoG1", c, C["names"]) for c in cw.chrom), "results/docs16_rel16/release16.chrom.tsv",
          "SH cneoG1: the wrong call's chromosome")
    N.add("cneo_wrong_conf", ",".join(cw.confidence), "results/docs16_rel16/release16.chrom.tsv", "SH cneoG1: the wrong call's confidence")
    cS = c16[(c16.arm == "S") & (c16.genome == "cneoG1")]
    N.add("cneo_S_conf", ",".join(sorted(set(cS.confidence))), "results/docs16_rel16/release16.chrom.tsv", "S cneoG1: confidence of the calls")
    # the readouts with the ends kept or out, from the 0.16.0 caches; each checked against the 0.13-era anchor tables
    anc = pd.read_csv(RES13 / "cis_trans_anchor.tsv", sep="\t")
    anci = pd.read_csv(RES13 / "cis_trans_anchor_int.tsv", sep="\t")
    tr = meta["truth"]
    checks = []

    def anchor(G, gid, arm, mode, sign, key, what):
        r = G["res"][(mode, sign)]
        kr, nr = on_truth(G, tr[gid], mode, sign)
        N.add(f"{key}_kn", f"{r['k']}/{r['n']}", f"benchmark_cis_trans_cache (runs/rel16/{arm}/{gid})", what)
        N.add(f"{key}_p", float(r["p"]), f"benchmark_cis_trans_cache (runs/rel16/{arm}/{gid})", what + " (one sign)")
        N.add(f"{key}_right", f"{kr}/{nr}", f"benchmark_cis_trans_cache + truth (runs/rel16/{arm}/{gid})",
              what + ": cis extremes within 1 Mb of CENH3")
        t = anc if mode == "all" else anci
        o = t[(t.arm == arm) & (t.id == gid) & (t.sign == ("peak" if sign > 0 else "dip"))].iloc[0]
        checks.append((key, (int(o.k), int(o.n), o.ext_right) == (r["k"], r["n"], f"{kr}/{nr}") and
                       bool(np.isclose(float(o.p), float(f"{r['p']:.3g}"), rtol=1e-6, atol=0))))
    anchor(A, "aetAL878", "SH", "all", +1, "aet_whole", "SH aetAL878, ends kept, cis peaks")
    anchor(A, "aetAL878", "SH", "int", +1, "aet_int", "SH aetAL878, terminal eighths out, cis peaks")
    anchor(Pt, "potato", "SH", "all", -1, "pot_dip_whole", "SH potato, ends kept, cis dips")
    anchor(Pt, "potato", "SH", "int", -1, "pot_dip", "SH potato, terminal eighths out (as the release reads), cis dips")
    anchor(T, "tomato", "SH", "int", +1, "tom_int", "SH tomato, terminal eighths out, cis peaks")
    anchor(T, "tomato", "SH", "all", +1, "tom_whole", "SH tomato, ends kept, cis peaks")
    for sg, s in (("peak", +1), ("dip", -1)):
        anchor(R, "tomato", "TNEG", "int", s, f"tneg_{sg}", f"TNEG tomato (randomized), interior {sg}s")
    N.add("anchor_tables_reproduced", f"{sum(ok for _, ok in checks)}/{len(checks)}",
          "results/cis_trans_anchor*.tsv (0.13 era) vs the 0.16.0 caches", "rows equal in k, n, P (3 digits) and extremes on CENH3: "
          + ", ".join(f"{k} {'=' if ok else '!='}" for k, ok in checks))
    st = meta["same_state"]
    N.add("hic_state_same", all(all(v for k_, v in d.items() if k_ != "npz") for d in st.values()),
          "benchmark_cis_trans_cache_meta.pkl", "0.13 and 0.16.0 runs: same matrix file, bins, members, hub: "
          + ", ".join(f"{a} {g} {d['npz']}" for (a, g), d in st.items()))
    r = T["res"][("int", 1)]
    N.add("tom_k_fig", f"{r['k']}/{r['n']}", "benchmark_cis_trans_cache_tomato.pkl", "anchored_on_cores re-read; equals the release")
    N.add("tom_chance_sum", round(sum(r["chance"].values()), 2), "benchmark_cis_trans_cache_tomato.pkl",
          "tomato: chromosomes expected to meet by chance (sum of zone shares)")
    N.add("tneg_line", f"{R['line']['k']}/{R['line']['n']} ({'peak' if R['line']['sign'] > 0 else 'dip'}, "
          f"{R['line']['status']})", "benchmark_cis_trans_cache_tneg.pkl", "cis_trans_line on the randomized tomato library")
    rc = C["res"][("int", 1)]
    N.add("cneo_k_fig", f"{rc['k']}/{rc['n']}", "benchmark_cis_trans_cache_cneo.pkl", "re-read; equals the release")
    N.add("cneo_chance_sum", round(sum(rc["chance"].values()), 2), "benchmark_cis_trans_cache_cneo.pkl",
          "C. neoformans G1: chromosomes expected to meet by chance")
    assert f"{r['k']}/{r['n']}" == N["tom_kn"] and f"{rc['k']}/{rc['n']}" == N["cneo_kn"], "cache and release disagree"
    return N


# ============================================================================== drawing helpers
def setup():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Liberation Sans", "Arial"], "font.size": 6.5,
                         "pdf.fonttype": 42, "ps.fonttype": 42, "axes.linewidth": 0.5, "xtick.major.width": 0.5,
                         "ytick.major.width": 0.5, "xtick.major.size": 2, "ytick.major.size": 2, "xtick.labelsize": 6,
                         "ytick.labelsize": 6, "axes.labelsize": 6.5, "legend.fontsize": 6, "hatch.linewidth": 0.4,
                         "hatch.color": "#B5B5B5", "mathtext.fontset": "custom", "mathtext.it": "Liberation Sans:italic",
                         "mathtext.rm": "Liberation Sans", "mathtext.bf": "Liberation Sans:bold", "axes.unicode_minus": False})
    return plt


class Page:
    def __init__(self, plt):
        self.plt = plt
        self.fig = plt.figure(figsize=(FW / 25.4, FH / 25.4))

    def ax(self, x, y, w, h, **kw):
        """Axes at x, y (mm from the top-left corner), w x h mm."""
        return self.fig.add_axes([x / FW, 1 - (y + h) / FH, w / FW, h / FH], **kw)

    def text(self, x, y, s, **kw):
        return self.fig.text(x / FW, 1 - y / FH, s, **kw)

    def letter(self, x, y, s):
        self.fig.text(x / FW, 1 - y / FH, s, fontsize=9, fontweight="bold", va="top", ha="left")


def clean(a, keep=()):
    for s in ("top", "right", "left", "bottom"):
        a.spines[s].set_visible(s in keep)
    if "bottom" not in keep:
        a.set_xticks([])
    if "left" not in keep:
        a.set_yticks([])


def bins_x(nb, binsize):
    return (np.arange(nb) + 0.5) * binsize / 1e6


def runs(mask):
    """(start, end) index runs where mask is True."""
    out, j, n = [], 0, len(mask)
    while j < n:
        if mask[j]:
            i = j
            while j < n and mask[j]:
                j += 1
            out.append((i, j))
        else:
            j += 1
    return out


def shade_ends(a, nb, binsize, ylim):
    from matplotlib.patches import Rectangle
    e = int(np.ceil(END_ZONE * nb))
    for i, j in ((0, e), (nb - e, nb)):
        a.add_patch(Rectangle((i * binsize / 1e6, ylim[0]), (j - i) * binsize / 1e6, ylim[1] - ylim[0], fc=C_END, ec="none",
                              zorder=0))


def hatch_masked(a, good, binsize, ylim):
    from matplotlib.patches import Rectangle
    for i, j in runs(~np.asarray(good, bool)):
        a.add_patch(Rectangle((i * binsize / 1e6, ylim[0]), (j - i) * binsize / 1e6, ylim[1] - ylim[0], fill=False,
                              hatch="//////", ec=C_HATCH, lw=0, zorder=0.5))


def band(a, iv, color=C_CEN, z=0.2, min_mb=0.0):
    s, e = iv[0] / 1e6, iv[1] / 1e6
    if e - s < min_mb:                                  # too narrow to see at this scale: widened about its center
        m = (s + e) / 2
        s, e = m - min_mb / 2, m + min_mb / 2
    a.axvspan(s, e, color=color, lw=0, zorder=z)


def logy(a, lim, ticks=(1, 2), labels=True):
    a.set_yscale("log")
    a.set_ylim(*lim)
    a.set_yticks(list(ticks))
    a.set_yticks([], minor=True)
    a.set_yticklabels([f"{t:g}×" for t in ticks] if labels else [])
    a.axhline(1, color="#BBBBBB", lw=0.35, ls=(0, (1, 1.5)), zorder=0.6)


def truth_iv(tr, c):
    return min(x[0] for x in tr[c]), max(x[1] for x in tr[c])


# ============================================================================== example
def pick_example(T, meta, N):
    """Among tomato chromosomes whose trans peak meets the cis extreme and that the 0.16.0 call with Hi-C (arm SH) gets right: the
    one whose extreme height and CENH3 width sit nearest the median of those chromosomes (each scaled by its median absolute
    deviation). The 0.13 figure's rule, with the 0.16.0 call in place of the 0.13 call."""
    c16 = pd.read_csv(RES / "release16.chrom.tsv", sep="\t")
    ok = c16[(c16.arm == "SH") & (c16.genome == "tomato")].set_index("chrom").ok
    r = T["res"][("int", 1)]
    tr = meta["truth"]["tomato"]
    rows = []
    for c in T["names"]:
        j = T["ext"]["int"][1][c]
        s, e = truth_iv(tr, c)
        rows.append(dict(chrom=c, label=chr_label("tomato", c, T["names"]), meets=c in r["peaks"], right=ok[c] == "R",
                         height=float(T["prof"][c][j]), cenh3_kb=(e - s) / 1e3, ext_bin=j))
    d = pd.DataFrame(rows)
    cand = d[d.meets & d.right].copy()
    for col in ("height", "cenh3_kb"):
        med = cand[col].median()
        mad = (cand[col] - med).abs().median() or 1.0
        cand[col + "_z"] = (cand[col] - med) / mad
    cand["dist"] = np.hypot(cand.height_z, cand.cenh3_kb_z)
    cand = cand.sort_values("dist")
    # the skill: verify structurally, then choose among the clean candidates (the 0.13 rule asserted the check on its pick)
    cand["clean"] = [struct_checks(T, meta, c)["ok"] for c in cand.chrom]
    fails = cand[~cand.clean]
    def why(c):
        r = struct_checks(T, meta, c)
        return "; ".join({"box_unmasked": "unmappable bin within 4 bins of the extreme", "outside_ends": "in a terminal eighth",
                          "dist_bins": "far from CENH3"}[k] for k in ("box_unmasked", "outside_ends") if r[k] is False)
    N.add("ex_screened_out", ", ".join(f"{r.label} ({why(r.chrom)})" for r in fails.itertuples()) or "none",
          "benchmark_cis_trans.py struct_checks", "candidates failing the structural checks (not drawn as the example)")
    first = cand.iloc[0]
    clean = cand[cand.clean]
    ex = clean.iloc[0]
    N.add("ex_rule_first", first.label, "benchmark_cis_trans.py pick_example", "nearest the medians before the structural screen")
    N.add("ex_n_candidates", len(cand), "benchmark_cis_trans_cache_tomato.pkl + results/docs16_rel16/release16.chrom.tsv",
          "tomato chromosomes where the trans peak meets the cis extreme and the 0.16.0 SH call is right")
    N.add("ex_candidates", ",".join(f"{r.label} {r.dist:.2f}{'' if r.clean else ' (screened out)'}" for r in cand.itertuples()),
          "benchmark_cis_trans.py pick_example", "candidates and their distance from the medians (MAD units), nearest first")
    N.add("ex_n_clean", len(clean), "benchmark_cis_trans.py struct_checks", "candidates passing the structural checks")
    N.add("ex_label", ex.label, "benchmark_cis_trans.py pick_example", "the example chromosome")
    N.add("ex_dist", round(float(ex.dist), 2), "benchmark_cis_trans.py pick_example", "example: distance from the medians (MAD units)")
    N.add("ex_runner_up_dist", round(float(clean.iloc[1].dist), 2), "benchmark_cis_trans.py pick_example",
          "runner-up's distance (clean candidates)")
    N.add("ex_height", round(ex.height, 2), "benchmark_cis_trans_cache_tomato.pkl", "example: cis extreme height (obs./exp.)")
    N.add("ex_cenh3_kb", int(ex.cenh3_kb), "V8/data/tomato_scau/slH1706/truth.chip_truth.bed", "example: CENH3 domain width")
    N.add("cand_height_med", round(cand.height.median(), 2), "benchmark_cis_trans_cache_tomato.pkl", "median over the candidates")
    N.add("cand_cenh3_med", int(cand.cenh3_kb.median()), "truth.chip_truth.bed", "median CENH3 width over the candidates (kb)")
    N.add("ex_runner_up", clean.iloc[1].label, "benchmark_cis_trans.py pick_example", "next nearest the median (clean)")
    log(cand[["label", "height", "cenh3_kb", "dist"]].to_string())
    return ex.chrom, cand


def struct_checks(T, meta, c):
    """Structural checks before drawing: extreme within 4 bins of CENH3, outside the terminal eighths, no masked bin in its
    crossing box (bins j-4 .. j+3) nor in the zone j-4 .. j+4."""
    b = T["bin"]
    j = T["ext"]["int"][1][c]
    nb = T["nb"][c]
    s, e = truth_iv(meta["truth"]["tomato"], c)
    lo, hi = s // b, (e - 1) // b
    dist_bins = 0 if lo <= j <= hi else min(abs(j - lo), abs(j - hi))
    e8 = int(np.ceil(END_ZONE * nb))
    good = T["good"][c]
    box_ok = bool(good[j - HALF:j + HALF + 1].all())
    res = dict(dist_bins=dist_bins, outside_ends=bool(e8 <= j < nb - e8), box_unmasked=box_ok)
    res["ok"] = bool(dist_bins <= 4 and res["outside_ends"] and box_ok)
    return res


def check_example(T, meta, c, N):
    res = struct_checks(T, meta, c)
    dist_bins, box_ok = res["dist_bins"], res["box_unmasked"]
    N.add("ex_check_dist_bins", dist_bins, "benchmark_cis_trans_cache_tomato.pkl", "example: bins between cis extreme and CENH3")
    N.add("ex_check_outside_ends", res["outside_ends"], "benchmark_cis_trans_cache_tomato.pkl", "example: extreme outside 1/8 ends")
    N.add("ex_check_box_unmasked", box_ok, "benchmark_cis_trans_cache_tomato.pkl", "example: no masked bin within 4 bins")
    assert dist_bins <= 4 and res["outside_ends"] and box_ok, f"example {c} fails the structural checks: {res}"
    return res


# ============================================================================== panel A: tomato, every chromosome
YL_CIS = (0.3, 2.6)         # cis: contact across the bin, obs./exp. (log); covers every tomato value, real and randomized
YL_TR = (0.3, 3.0)          # trans: contact with the other chromosomes' extremes, obs./exp. (log)


def meet_line(P, a_top, a_bot, x_top, x_bot, color=C_MEET):
    """A connector from the top of a_top at x_top to the bottom of a_bot at x_bot (the two readouts meet)."""
    from matplotlib.patches import ConnectionPatch
    cp = ConnectionPatch((x_top, a_top.get_ylim()[1]), (x_bot, a_bot.get_ylim()[0]), "data", "data", axesA=a_top, axesB=a_bot,
                         color=color, lw=0.6, zorder=0.7, alpha=0.9)
    P.fig.add_artist(cp)


def carets(a, x, v, ylim, color=C_RAND):
    """Values beyond the shared axis are drawn at its edge and marked (not hidden)."""
    v = np.asarray(v, float)
    hi, lo = v > ylim[1], v < ylim[0]
    if hi.any():
        a.plot(x[hi], np.full(hi.sum(), ylim[1]), ls="none", marker="^", ms=1.6, color=color, mec="none", clip_on=False,
               zorder=1.5)
    if lo.any():
        a.plot(x[lo], np.full(lo.sum(), ylim[0]), ls="none", marker="v", ms=1.6, color=color, mec="none", clip_on=False,
               zorder=1.5)
    return int(hi.sum()), int(lo.sum())


def draw_pair(P, G, R, c, x0, y0, w, mbmm, hc, ht, gap, tr_iv, labels, key=("int", 1), ext_mode="int", show_rand=True,
              ylc=YL_CIS, ylt=YL_TR, lw=0.55, marker=3.2, ticklabels=True, min_band_mm=0.0, clipped=None):
    """One chromosome: the cis track over the trans track (same x), randomized partner in gray, end zones, masked bins, truth
    band, cis extreme (down triangle), trans peak (up triangle) and the connector where the trans peak meets the extreme."""
    b = G["bin"]
    nb = G["nb"][c]
    x = bins_x(nb, b)
    L = G["lengths"][c] / 1e6
    ac = P.ax(x0, y0, w, hc)
    at = P.ax(x0, y0 + hc + gap, w, ht)
    for a, yl in ((ac, ylc), (at, ylt)):
        a.set_xlim(0, L)
        shade_ends(a, nb, b, yl)
        hatch_masked(a, G["good"][c], b, yl)
        if tr_iv is not None:
            band(a, tr_iv, min_mb=min_band_mm / mbmm)
    v = G["prof"][c]
    S = G["S"][key][c]
    if show_rand and R is not None and c in R["prof"]:
        xr = bins_x(R["nb"][c], R["bin"])
        for a, vr, yl in ((ac, R["prof"][c], ylc), (at, R["S"][key][c], ylt)):
            a.plot(xr, np.clip(vr, yl[0], yl[1]), color=C_RAND, lw=lw * 0.75, zorder=1)
            n_hi, n_lo = carets(a, xr, vr, yl)
            if clipped is not None:
                clipped.append((c, a is ac, n_hi, n_lo, float(np.nanmax(vr)), float(np.nanmin(vr))))
    for a, vv, yl, col in ((ac, v, ylc, C_CIS), (at, S, ylt, C_TRANS)):
        a.plot(x, np.clip(vv, yl[0], yl[1]), color=col, lw=lw, zorder=2)
        n_hi, n_lo = carets(a, x, vv, yl, color=col)
        if clipped is not None and (n_hi or n_lo):
            clipped.append((c, "real", n_hi, n_lo, float(np.nanmax(vv)), float(np.nanmin(vv))))
    logy(ac, ylc, labels=ticklabels)
    logy(at, ylt, labels=ticklabels)
    for a in (ac, at):
        clean(a, keep=("left",))
        a.tick_params(axis="y", length=1.5, pad=1)
    j = G["ext"][ext_mode][key[1]][c]
    tp = G["trans_peak"][key][c]
    xj, xt = (j + 0.5) * b / 1e6, (tp + 0.5) * b / 1e6
    ac.plot([xj], [ylc[1]], marker="v", ms=marker, color=C_CIS, mec="none", clip_on=False, zorder=5, transform=ac.transData)
    at.plot([xt], [ylt[0]], marker="^", ms=marker, color=C_TRANS, mec="none", clip_on=False, zorder=5)
    meets = c in G["res"][key]["peaks"]
    if meets:
        meet_line(P, ac, at, xj, xt)
    if labels:
        yl_ = ac.transAxes.inverted().transform(ac.transAxes.transform((0, 0)))  # noqa: F841
        P.text(x0 - 1.2, y0 + hc + gap / 2, labels, ha="right", va="center", fontsize=6.5)
    return ac, at, meets


def calls_track(P, x0, y0, w, h, L, items, tr_iv=None):
    """Calls as boxes at their size: items [(interval or None, color, row)] with row 0 on top."""
    from matplotlib.patches import Rectangle
    a = P.ax(x0, y0, w, h)
    a.set_xlim(0, L)
    a.set_ylim(0, 2)
    if tr_iv is not None:
        band(a, tr_iv)
    for iv, col, row in items:
        y = 1.05 if row == 0 else 0.05
        a.axhline(y + 0.45, color="#D8D8D8", lw=0.3, zorder=0.3)
        if iv is not None:
            a.add_patch(Rectangle((iv[0] / 1e6, y), (iv[1] - iv[0]) / 1e6, 0.9, fc=col, ec="none", zorder=2))
    clean(a)
    return a


def panel_A(P, T, R, meta, N, y0, ex):
    tr = meta["truth"]["tomato"]
    cS, cSH = meta["calls"][("S", "tomato")], meta["calls"][("SH", "tomato")]      # 0.16.0: FASTA only; with Hi-C
    names = T["names"]
    mbmm = 80.0 / (max(T["lengths"].values()) / 1e6)
    hc, ht, gap, hk, rowh = 4.4, 4.4, 0.5, 2.0, 14.4
    cols = (10.0, 100.0)
    axes, clipped = {}, []
    for i, c in enumerate(names):
        col, row = divmod(i, 6)
        x0, yy = cols[col], y0 + row * rowh
        L = T["lengths"][c] / 1e6
        w = L * mbmm
        iv = truth_iv(tr, c)
        ac, at, meets = draw_pair(P, T, R, c, x0, yy, w, mbmm, hc, ht, gap, iv, chr_label("tomato", c, names),
                                  ticklabels=False, clipped=clipped)
        ak = calls_track(P, x0, yy + hc + gap + ht + 0.5, w, hk, L,
                         [(cS[c][:2] if cS[c] else None, C_SAT, 0), (cSH[c][:2] if cSH[c] else None, "black", 1)], iv)
        axes[c] = (ac, at, ak, x0, yy, w)
    for col in (0, 1):                                  # one x axis per column, at the bottom
        c = names[col * 6 + 5]
        ac, at, ak, x0, yy, w = axes[c]
        Lmax = max(T["lengths"][names[col * 6 + r]] for r in range(6)) / 1e6
        ax = P.ax(x0, yy + hc + gap + ht + 0.5 + hk + 1.0, Lmax * mbmm, 0.1)
        ax.set_xlim(0, Lmax)
        clean(ax, keep=("bottom",))
        ax.set_xticks(np.arange(0, Lmax + 0.1, 20))
        ax.tick_params(axis="x", length=2, pad=1)
        ya = yy + hc + gap + ht + 0.5 + hk + 1.0        # the axis line, mm from the top
        P.text(x0 + 0.5 * Lmax * mbmm, ya + 3.4, "tomato position (Mb)", ha="center", va="top", fontsize=6.5)
    rc = [x for x in clipped if x[1] is True and x[2]]
    N.add("tneg_clipped_chroms", ",".join(chr_label("tomato", x[0], names) for x in rc), "benchmark_cis_trans_cache_tneg.pkl",
          "randomized cis profile above the shared axis (carets)")
    N.add("tneg_cis_max", round(max(x[4] for x in clipped if x[1] is True), 2), "benchmark_cis_trans_cache_tneg.pkl",
          "randomized cis profile maximum (beside unmappable bins)")
    real = [x for x in clipped if x[1] == "real"]
    N.add("tom_clipped", ", ".join(f"{chr_label('tomato', x[0], names)} ({x[2] + x[3]} bin)" for x in real) or "none",
          "benchmark_cis_trans_cache_tomato.pkl", "tomato tracks with values beyond the shared axis (carets)")
    return axes, mbmm, rowh


# ============================================================================== panel B: one chromosome, both readouts
WIN = 48                    # zoom: the extreme +- 48 bins (12 Mb at 250 kb)
MAXD = 24                   # triangle: pairs up to 24 bins (6 Mb) apart


def offsets(G):
    off, o = {}, 0
    for c in G["names"]:
        off[c] = o
        o += G["nb"][c]
    return off


def panel_B(P, T, R, meta, chip, N, c, x0, y0, w):
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    from matplotlib.patches import Polygon, Rectangle
    b = T["bin"]
    binmb = b / 1e6
    nb = T["nb"][c]
    j = T["ext"]["int"][1][c]
    i0, i1 = max(0, j - WIN), min(nb, j + WIN + 1)
    xlim = (i0 * binmb, i1 * binmb)
    mbmm = w / (xlim[1] - xlim[0])
    tr_iv = truth_iv(meta["truth"]["tomato"], c)
    lab = chr_label("tomato", c, T["names"])
    # --- triangle: observed/expected contact within the chromosome, each pair drawn above its midpoint
    tri_h = MAXD / 2 * binmb * mbmm
    at = P.ax(x0, y0, w, tri_h)
    X = T["oe"][c]
    n_ = i1 - i0
    I, J = np.meshgrid(np.arange(i0, i1 + 1), np.arange(i0, i1 + 1), indexing="ij")
    Xr, Yr = (I + J) / 2 * binmb, (J - I) / 2 * binmb
    Cm = np.full((n_, n_), np.nan)
    for a_ in range(n_):
        for b_ in range(a_ + 1, min(n_, a_ + MAXD + 1)):
            Cm[a_, b_] = X[i0 + a_, i0 + b_]
    at.add_patch(Polygon([(i0 * binmb, 0), (i1 * binmb, 0), ((i1 - MAXD / 2) * binmb, MAXD / 2 * binmb),
                          ((i0 + MAXD / 2) * binmb, MAXD / 2 * binmb)], closed=True, fill=False, hatch="//////", lw=0,
                         ec=C_HATCH, zorder=0))
    norm_c = LogNorm(vmin=0.5, vmax=2.0)
    cmap_c = plt.get_cmap("PRGn_r").copy()
    cmap_c.set_bad((1, 1, 1, 0))
    at.pcolormesh(Xr, Yr, np.ma.masked_invalid(Cm), cmap=cmap_c, norm=norm_c, shading="flat", rasterized=True,
                  edgecolors="none", zorder=1)
    dia = np.array([(j - HALF / 2, HALF / 2), (j, 0), (j + HALF / 2, HALF / 2), (j, HALF)]) * binmb
    at.add_patch(Polygon(dia, closed=True, fill=False, ec="black", lw=0.7, zorder=5))
    at.set_xlim(*xlim)
    at.set_ylim(0, MAXD / 2 * binmb)
    clean(at, keep=("left",))
    at.set_yticks([0, 1, 2, 3])
    at.set_yticklabels(["0", "2", "4", "6"])
    at.tick_params(axis="y", length=1.5, pad=1)
    at.set_ylabel("distance (Mb)", labelpad=1, fontsize=6)
    # --- trans rows: the chromosome's contact with each other chromosome's anchor window
    off = offsets(T)
    others = [y for y in T["names"] if y != c and y in T["partner_oe"]]
    rh = 1.15
    yr = y0 + tri_h + 1.6
    ap = P.ax(x0, yr, w, rh * len(others))
    Mrow = np.vstack([T["partner_oe"][y][off[c] + i0: off[c] + i1] for y in others])
    norm_t = LogNorm(vmin=1.0, vmax=4.0)
    cmap_t = plt.get_cmap("Blues").copy()
    cmap_t.set_bad((1, 1, 1, 0))
    cmap_t.set_under("white")
    ap.add_patch(Rectangle((xlim[0], 0), xlim[1] - xlim[0], len(others), fill=False, hatch="//////", lw=0, ec=C_HATCH,
                           zorder=0))
    ap.pcolormesh(np.arange(i0, i1 + 1) * binmb, np.arange(len(others) + 1), np.ma.masked_invalid(Mrow), cmap=cmap_t,
                  norm=norm_t, shading="flat", rasterized=True, edgecolors="none", zorder=1)
    for k_ in range(1, len(others)):
        ap.axhline(k_, color="white", lw=0.25, zorder=2)
    ap.set_xlim(*xlim)
    ap.set_ylim(len(others), 0)
    clean(ap, keep=("left",))
    ap.set_yticks([0.5, len(others) - 0.5])
    ap.set_yticklabels([chr_label("tomato", others[0], T["names"]), chr_label("tomato", others[-1], T["names"])], fontsize=6)
    ap.tick_params(axis="y", length=0, pad=1)
    P.text(x0 - 6.5, yr + rh * len(others) / 2, "other\nchromosomes'\nextremes", ha="right", va="center", fontsize=6,
           color=C_TRANS, linespacing=0.95)
    # --- the two readouts as lines, then CENH3 ChIP and the calls
    hc = 6.6
    yl = yr + rh * len(others) + 1.6
    G1 = dict(T)
    ac = P.ax(x0, yl, w, hc)
    at2 = P.ax(x0, yl + hc + 0.6, w, hc)
    x = bins_x(nb, b)
    for a, v, vr, col, ylim in ((ac, T["prof"][c], R["prof"].get(c), C_CIS, YL_CIS),
                                (at2, T["S"][("int", 1)][c], R["S"][("int", 1)].get(c), C_TRANS, YL_TR)):
        a.set_xlim(*xlim)
        hatch_masked(a, T["good"][c], b, ylim)
        band(a, tr_iv)
        if vr is not None:
            a.step(bins_x(R["nb"][c], R["bin"]), vr, where="mid", color=C_RAND, lw=0.45, zorder=1)
        a.step(x, v, where="mid", color=col, lw=0.7, zorder=2)
        logy(a, ylim)
        clean(a, keep=("left",))
        a.tick_params(axis="y", length=1.5, pad=1)
    P.text(x0 - 6.5, yl + hc / 2, "cis", ha="right", va="center", fontsize=6.5, color=C_CIS)
    P.text(x0 - 6.5, yl + hc + 0.6 + hc / 2, "trans", ha="right", va="center", fontsize=6.5, color=C_TRANS)
    tp = T["trans_peak"][("int", 1)][c]
    xj, xt = (j + 0.5) * binmb, (tp + 0.5) * binmb
    ac.plot([xj], [YL_CIS[1]], marker="v", ms=3.6, color=C_CIS, mec="none", clip_on=False, zorder=5)
    at2.plot([xt], [YL_TR[0]], marker="^", ms=3.6, color=C_TRANS, mec="none", clip_on=False, zorder=5)
    meet_line(P, ac, at2, xj, xt)
    # CENH3 ChIP / input, 100-kb bins (the truth rule's own track)
    yc = yl + 2 * hc + 0.6 + 1.6
    hchip = 4.2
    ah = P.ax(x0, yc, w, hchip)
    e = chip[c] if chip is not None else np.ones(int(np.ceil(T["lengths"][c] / 1e5)))
    xc = (np.arange(len(e)) + 0.5) * 0.1
    sel = (xc >= xlim[0] - 0.1) & (xc <= xlim[1] + 0.1)
    band(ah, tr_iv)
    ah.fill_between(xc[sel], 1, e[sel], step="mid", color=C_CEN_EDGE, lw=0, zorder=2)
    ah.set_xlim(*xlim)
    ah.set_yscale("log")
    top = max(2.0, float(np.nanmax(e[sel])))
    ah.set_ylim(0.8, top * 1.15)
    ah.set_yticks([1, 10])
    ah.set_yticks([], minor=True)
    ah.set_yticklabels(["1×", "10×"])
    clean(ah, keep=("left",))
    ah.tick_params(axis="y", length=1.5, pad=1)
    P.text(x0 - 6.5, yc + hchip / 2, "CENH3\nChIP / input", ha="right", va="center", fontsize=6, linespacing=0.95)
    # calls
    cS, cSH = meta["calls"][("S", "tomato")][c], meta["calls"][("SH", "tomato")][c]
    ak = P.ax(x0, yc + hchip + 0.8, w, 3.6)
    ak.set_xlim(*xlim)
    ak.set_ylim(0, 2)
    band(ak, tr_iv)
    for iv, col, y in ((cS, C_SAT, 1.05), (cSH, "black", 0.05)):
        ak.axhline(y + 0.45, color="#D8D8D8", lw=0.3, zorder=0.3)
        if iv:
            ak.add_patch(Rectangle((iv[0] / 1e6, y), (iv[1] - iv[0]) / 1e6, 0.9, fc=col, ec="none", zorder=2))
    clean(ak)
    P.text(x0 - 1.0, yc + hchip + 0.8 + 0.9, "FASTA only", ha="right", va="center", fontsize=6, color="#444444")
    P.text(x0 - 1.0, yc + hchip + 0.8 + 2.7, "with Hi-C", ha="right", va="center", fontsize=6, color="#444444")
    axx = P.ax(x0, yc + hchip + 0.8 + 3.6 + 0.6, w, 0.1)
    axx.set_xlim(*xlim)
    clean(axx, keep=("bottom",))
    t0 = np.ceil(xlim[0] / 2) * 2
    axx.set_xticks(np.arange(t0, xlim[1] + 0.01, 2))
    axx.tick_params(axis="x", length=2, pad=1)
    axx.set_xlabel(f"tomato {lab} position (Mb)", labelpad=1)
    # colorbars, in the empty corners above the triangle's flanks
    cbw = 13.0
    cb1 = P.ax(x0 + 6.0, y0 - 5.0, cbw, 1.2)
    cb2 = P.ax(x0 + 32.0, y0 - 5.0, cbw, 1.2)
    for a, cm, nm, ticks, lab_, ha in ((cb1, cmap_c, norm_c, (0.5, 1, 2), "cis contact, obs./exp.", "center"),
                                       (cb2, cmap_t, norm_t, (1, 2, 4), "trans contact, obs./exp.", "center")):
        g = np.geomspace(nm.vmin, nm.vmax, 128)
        a.pcolormesh(g, [0, 1], g[None, :-1], cmap=cm, norm=nm, shading="flat", rasterized=True)
        a.set_xscale("log")
        a.set_xlim(nm.vmin, nm.vmax)
        a.set_yticks([])
        a.set_xticks(list(ticks))
        a.set_xticks([], minor=True)
        a.set_xticklabels([("½×" if t == 0.5 else f"{t:g}×") for t in ticks], fontsize=6)
        a.tick_params(axis="x", length=1.2, pad=0.8)
        for s in a.spines.values():
            s.set_linewidth(0.3)
        a.text(0.5, 1.5, lab_, transform=a.transAxes, ha="center", va="bottom", fontsize=6)
    bottom = yc + hchip + 0.8 + 3.6 + 0.6
    N.add("ex_trans_peak_offset_bins", int(tp - j), "benchmark_cis_trans_cache_tomato.pkl", "example: trans peak minus cis extreme")
    N.add("ex_chip_max", round(top, 1), "benchmark_cis_trans_cache_chip.pkl", "example: max CENH3 ChIP/input in the window")
    if cS:
        satm = (xc >= cS[0] / 1e6) & (xc <= cS[1] / 1e6)
        N.add("ex_chip_sat", round(float(np.nanmax(e[satm])), 1), "benchmark_cis_trans_cache_chip.pkl",
              "example: max CENH3 ChIP/input over the FASTA-only call's satellite array")
    return dict(top=y0, bottom=bottom, xlim=xlim, x0=x0, w=w)


# ============================================================================== panel C: schematic (flat cartoon)
def panel_C(P, x0, y0, w, h):
    from matplotlib.patches import Arc, FancyBboxPatch, Rectangle
    a = P.ax(x0, y0, w, h)
    a.set_xlim(0, w)
    a.set_ylim(h, 0)
    a.set_aspect("equal")
    clean(a)
    rows = [(5.0, 0.42, 1.00), (15.5, 0.62, 0.86), (26.0, 0.33, 0.93)]       # y, centromere position, relative length
    L0, xl = w - 9.0, 6.5
    cen_xy = []
    for y, cp, rl in rows:
        L = L0 * rl
        e8 = 0.125 * L
        for xs in (xl, xl + L - e8):
            a.add_patch(Rectangle((xs, y - 2.3), e8, 4.6, fc=C_END, ec="none", zorder=0))
        a.add_patch(FancyBboxPatch((xl, y - 0.7), L, 1.4, boxstyle="round,pad=0,rounding_size=0.7", fc="#CFCFCF", ec="none",
                                   zorder=1))
        xc = xl + cp * L
        a.add_patch(Rectangle((xc - 0.5, y - 0.75), 1.0, 1.5, fc=C_CIS, ec="none", zorder=2))
        for r in (1.3, 2.3, 3.3):                       # flanks meeting across the centromere: contacts that cross it
            a.add_patch(Arc((xc, y - 0.7), 2 * r, 2 * r * 0.9, theta1=0, theta2=180, color=C_CIS, lw=0.6, zorder=3))
        for xe, sgn in ((xl, 1), (xl + L, -1)):           # the ends fold too
            a.add_patch(Arc((xe + sgn * 1.6, y - 0.7), 3.2, 2.6, theta1=0, theta2=180, color="#8C8C8C", lw=0.5, zorder=3))
        cen_xy.append((xc, y))
    for (x1, y1), (x2, y2) in ((cen_xy[0], cen_xy[1]), (cen_xy[1], cen_xy[2]), (cen_xy[0], cen_xy[2])):
        a.plot([x1, x2], [y1 + 0.8, y2 - 0.8] if y2 > y1 else [y1, y2], color=C_TRANS, lw=0.8, ls=(0, (2, 1.2)), zorder=2.5)
    ends = [(xl, y) for y, _, _ in rows]
    for (x1, y1), (x2, y2) in zip(ends[:-1], ends[1:]):
        a.plot([x1 - 1.2, x2 - 1.2], [y1, y2], color="#8C8C8C", lw=0.7, ls=(0, (2, 1.2)), zorder=2.5)
    xc0, yc0 = cen_xy[0]
    a.text(xc0 + 4.2, rows[0][0] - 2.9, "cis", color=C_CIS, fontsize=6.5, ha="left", va="center")
    a.text((cen_xy[1][0] + cen_xy[2][0]) / 2 + 3.0, (rows[1][0] + rows[2][0]) / 2 + 0.6, "trans", color=C_TRANS, fontsize=6.5,
           ha="left", va="center")
    a.text(xl - 1.8, (rows[0][0] + rows[1][0]) / 2, "ends", color="#6E6E6E", fontsize=6, ha="right", va="center")
    a.text(xl + L0 * rows[2][2], rows[2][0] + 3.0, "ends: not read", fontsize=6, ha="right", va="top", color="#6E6E6E")
    xc2, yc2 = cen_xy[2]
    a.text(xc2, yc2 + 3.0, "centromere", fontsize=6, ha="center", va="top", color=C_CIS)
    return a


# ============================================================================== panel D: Ae. tauschii, ends kept
def panel_D(P, A, meta, N, x0, y0, w, rowh):
    tr = meta["truth"]["aetAL878"]
    names = A["names"]
    mbmm = w / (max(A["lengths"].values()) / 1e6)
    key = ("all", 1)
    ylc = (0.12, 2.0)               # covers every value but exact zeros beside unmappable bins (carets)
    ylt = (0.12, 3.6)
    hc = (rowh - 2.6) / 2
    last = None
    for i, c in enumerate(names):
        L = A["lengths"][c] / 1e6
        ac, at, meets = draw_pair(P, A, None, c, x0, y0 + i * rowh, L * mbmm, mbmm, hc, hc, 0.4, truth_iv(tr, c),
                                  chr_label("aetAL878", c, names), key=key, ext_mode="all", show_rand=False, ylc=ylc, ylt=ylt,
                                  lw=0.45, marker=2.8, ticklabels=False, min_band_mm=0.6)
        last = (ac, at)
    Lmax = max(A["lengths"].values()) / 1e6
    ax = P.ax(x0, y0 + len(names) * rowh - 1.2, Lmax * mbmm, 0.1)
    ax.set_xlim(0, Lmax)
    clean(ax, keep=("bottom",))
    ax.set_xticks(np.arange(0, Lmax + 1, 200))
    ax.tick_params(axis="x", length=2, pad=1)
    ax.set_xlabel(r"$\it{Ae.\ tauschii}$ position (Mb)", labelpad=1)
    ext = A["ext"]["all"][1]
    e8 = {c: int(np.ceil(END_ZONE * A["nb"][c])) for c in names}
    in_end = sum(1 for c in names if ext[c] < e8[c] or ext[c] >= A["nb"][c] - e8[c])
    met = A["res"][key]["peaks"]
    met_end = sum(1 for c in met if ext[c] < e8[c] or ext[c] >= A["nb"][c] - e8[c])
    N.add("aet_ext_in_ends", f"{in_end}/{len(names)}", "benchmark_cis_trans_cache_aet.pkl", "whole-chromosome cis peaks in a terminal 1/8")
    N.add("aet_met_in_ends", f"{met_end}/{len(met)}", "benchmark_cis_trans_cache_aet.pkl", "of the meetings, in a terminal 1/8")
    N.add("aet_bin_mb", A["bin"] / 1e6, "benchmark_cis_trans_cache_aet.pkl", "Ae. tauschii Hi-C bin (Mb)")
    return last


# ============================================================================== panel F: C. neoformans G1, 14 chromosomes
def panel_F(P, Cn, meta, N, x0s, y0, ws, rowh):
    tr = meta["truth"]["cneoG1"]
    cSH = meta["calls"][("SH", "cneoG1")]                 # 0.16.0, with Hi-C
    names = Cn["names"]
    Lmax = max(Cn["lengths"].values()) / 1e6
    mbmm = ws / Lmax
    key = ("int", 1)
    ylc = (0.2, 2.3)
    ylt = (0.2, 16.0)
    hc = (rowh - 3.4) / 2
    for i, c in enumerate(names):
        col, row = divmod(i, 7)
        L = Cn["lengths"][c] / 1e6
        x0, yy = x0s[col], y0 + row * rowh
        ac, at, meets = draw_pair(P, Cn, None, c, x0, yy, L * mbmm, mbmm, hc, hc, 0.3, truth_iv(tr, c),
                                  chr_label("cneoG1", c, names), key=key, show_rand=False, ylc=ylc, ylt=ylt, lw=0.45,
                                  marker=2.6, ticklabels=False, min_band_mm=0.6)
        for a in (ac, at):
            a.set_yticks([])
            a.spines["left"].set_visible(False)
        ak = calls_track(P, x0, yy + 2 * hc + 0.3 + 0.3, L * mbmm, 0.9, L, [(cSH[c][:2] if cSH[c] else None, "black", 1)],
                         None)
        band(ak, truth_iv(tr, c), min_mb=0.6 / mbmm)
        ak.set_ylim(0, 1.05)
    for col in (0, 1):
        Lc = max(Cn["lengths"][names[col * 7 + r]] for r in range(7)) / 1e6
        ax = P.ax(x0s[col], y0 + 7 * rowh - 1.4, Lc * mbmm, 0.1)
        ax.set_xlim(0, Lc)
        clean(ax, keep=("bottom",))
        ax.set_xticks(np.arange(0, Lc + 0.01, 1.0))
        ax.tick_params(axis="x", length=2, pad=1)
    P.text(x0s[0] + (x0s[1] + Lc * mbmm - x0s[0]) / 2, y0 + 7 * rowh + 3.6,
           r"$\it{C.\ neoformans}$ position (Mb)", ha="center", va="top", fontsize=6.5)
    N.add("cneo_bin_kb", Cn["bin"] / 1e3, "benchmark_cis_trans_cache_cneo.pkl", "C. neoformans Hi-C bin (kb)")


# ============================================================================== panel E: every Hi-C library against chance
HOLO = {"celegans", "bmori", "rbrev", "cnig"}
PLANTS = {"tomato", "potato", "barley", "brachy", "aetAL878", "colcc", "sbBTx623", "ppat", "csat9930", "rice", "gmX26",
          "sver"}                                    # + the 13 further tomato accessions (la*, ts*), matched by name
FUNGI = {"cneoG1", "cneoA", "ncra", "vdah", "spom", "scer", "ztri", "msym", "ctro", "pfal"}
NICE = {"tomato": "tomato", "cneoG1": "C. neoformans", "brachy": "Brachypodium", "ncra": "N. crassa", "vdah": "V. dahliae",
        "gmX26": "soybean", "MMU8": "rhesus", "chm13": "human", "barley": "barley", "cneoA": "C. neoformans A", "spom": "S. pombe",
        "ts60": "tomato TS60"}
ITAL = {"cneoG1", "ncra", "vdah", "cneoA", "spom", "brachy"}
SUFFIX = {"cneoG1": " G1"}            # upright: a cell-cycle stage, not part of the species name


def libraries(N):
    """One row per contact map tested by trans contacts in the 0.16.0 release, with its chance: reused from the 0.13 figure's
    re-read where the same matrix gave the same k, n, verdict and sign (benchmark_cis_trans_libs.py table), re-read here for the
    new ones (benchmark_cis_trans_libs16.tsv), else the coverage-based estimate (chance_approx.py; light bars)."""
    sys.path.insert(0, str(SRC))
    import benchmark_cis_trans_chance as CA
    m = runs16()
    f = HERE / "benchmark_cis_trans_libs16.tsv"
    d = pd.read_csv(f, sep="\t") if f.exists() else pd.DataFrame(columns=["arm", "genome", "k", "n", "status", "sign", "npz",
                                                                          "expected"])
    if "error" in d:
        d = d[d.error.isna()]
    d = d.drop_duplicates(["arm", "genome"])[["arm", "genome", "k", "n", "status", "sign", "npz", "expected"]]
    m = m.merge(d, on=["arm", "genome"], how="left", suffixes=("", "_re"))
    rr = m.dropna(subset=["k_re"])
    agree = ((rr.k == rr.k_re) & (rr.n == rr.n_re) & (rr.verdict == rr.status) & (rr.sign == rr.sign_re) &
             (rr.npz == rr.npz_re))
    N.add("libs_reread", f"{int(agree.sum())}/{len(rr)}", "benchmark_cis_trans_libs16.tsv vs benchmark_cis_trans_runs16.tsv",
          "runs re-read from their 0.16.0 matrices whose k, n, verdict, sign and matrix equal the release's")
    assert bool(agree.all()), "a re-read library disagrees with the release"
    m["chance_src"] = np.where(m.expected.notna(), "0.13 re-read, same matrix and result",
                      np.where(m.expected_re.notna(), "0.16.0 re-read", "estimate"))
    m["expected"] = m.expected.fillna(m.expected_re)
    m["k_rel"], m["n_rel"], m["sign_rel"] = m.k, m.n, m.sign
    m["kind"] = np.where(m.arm.isin(["NEG", "NEG2"]), "randomized",
                np.where(m.genome.isin(HOLO), "holocentric",
                np.where(m.genome.isin(PLANTS) | m.genome.str.match(r"^(la\d|ts\d)"), "plants",
                np.where(m.genome.isin(FUNGI), "fungi", "animals"))))
    # one row per library: arms that reuse one contact map (held-out reruns HTH, rhesus FULL and CPG) give identical rows
    m["arm_rank"] = m.arm.map({"SH": 0, "FULL": 1, "CPG": 2, "HSH": 3, "HTH": 4}).fillna(5)
    n_runs = len(m)
    m = m.sort_values(["arm_rank", "genome"], kind="stable").drop_duplicates(["genome", "kind", "cis_line"], keep="first")
    m["chance_exact"] = m.expected.notna()
    for ix in m.index[~m.chance_exact]:
        e, n = CA.approx_expected(m.at[ix, "arm"], m.at[ix, "genome"])
        m.at[ix, "expected"] = e
    N.add("libs_runs", n_runs, "benchmark_cis_trans_runs16.tsv", "0.16.0 runs tested by trans contacts")
    N.add("chance_exact_n", int(m.chance_exact.sum()), "benchmark_cis_trans_runs16.tsv + _libs16.tsv",
          "libraries with the exact chance")
    N.add("chance_reused13_n", int((m.chance_src == "0.13 re-read, same matrix and result").sum()),
          "benchmark_cis_trans_runs16.tsv (from V12/figures/fig_cis_trans/fig_cis_trans_libs*.tsv)",
          "libraries whose exact chance is reused from the 0.13 re-read (same matrix file, same k/n/verdict/sign)")
    N.add("chance_reread16_n", int((m.chance_src == "0.16.0 re-read").sum()), "benchmark_cis_trans_libs16.tsv",
          "libraries whose chance was re-read from the 0.16.0 matrices (new in 0.16)")
    N.add("chance_approx_n", int((~m.chance_exact).sum()), "chance_approx.py", "libraries with the coverage-based estimate")
    m["share"] = m.k_rel / m.n_rel
    m["chance_share"] = m.expected / m.n_rel
    m[["arm", "genome", "kind", "sign", "k", "n", "p", "verdict", "acts", "expected", "chance_src", "share", "chance_share"]].to_csv(
        HERE / f"{STEM}_libraries.tsv", sep="\t", index=False)
    return m


def panel_E(P, Lb, N, x0, y0, w, h):
    from matplotlib.lines import Line2D
    hs = [Line2D([], [], marker="^", ls="none", ms=3, mfc="black", mec="black"),
          Line2D([], [], marker="^", ls="none", ms=3, mfc="#9C9C9C", mec="#9C9C9C"),
          Line2D([], [], marker="^", ls="none", ms=3, mfc="white", mec="#9C9C9C"),
          Line2D([], [], marker="^", ls="none", ms=3, mfc="white", mec="black"),
          Line2D([], [], marker="v", ls="none", ms=3, mfc="#9C9C9C", mec="#9C9C9C"),
          Line2D([], [], color="#404040", lw=0.9), Line2D([], [], color="#404040", lw=0.9, alpha=0.4),
          Line2D([], [], color="#9A9A9A", lw=0.6, ls=(0, (3, 2)))]
    kl = P.ax(x0 - 9.0, y0, w + 9.0, 6.0)
    kl.set_axis_off()
    names = ["pass (acts)", "support", "fail", "both signs pass", "dip sign", "chance", "chance, estimated", "majority"]
    if bool(Lb.chance_exact.all()):                  # every chance re-read exactly: no estimate to key
        hs, names = hs[:6] + hs[7:], names[:6] + names[7:]
    kl.legend(hs, names, loc="upper left", ncol=3, frameon=False, fontsize=6,
              handlelength=1.2, handletextpad=0.3, columnspacing=0.6, labelspacing=0.2, borderaxespad=0, borderpad=0)
    a = P.ax(x0, y0 + 9.5, w, h - 9.5)
    order = ["plants", "fungi", "animals", "holocentric", "randomized"]
    lab = {"plants": "plants", "fungi": "fungi,\nprotist", "animals": "animals", "holocentric": "holo-\ncentric",
           "randomized": "randomized Hi-C"}
    xs, pos = [], 0.0
    gaps = {"plants": 6.5, "fungi": 13.0, "animals": 13.0, "holocentric": 6.5, "randomized": 0.0}   # room for the group labels
    centers = {}
    for g in order:
        sub = Lb[Lb.kind == g].sort_values(["share", "n_rel"], ascending=[False, False])
        x = pos + np.arange(len(sub))
        centers[g] = (x[0], x[-1]) if len(sub) else (pos, pos)
        for xi, r in zip(x, sub.itertuples()):
            xs.append((xi, r))
        pos += len(sub) + gaps[g]
    a.set_xlim(-1.5, pos + 0.5)
    a.set_ylim(-3, 104)
    a.axhline(50, color="#9A9A9A", lw=0.5, ls=(0, (3, 2)), zorder=0.5)
    sty = {"pass": dict(fc="black", ec="black"), "support": dict(fc="#9C9C9C", ec="#9C9C9C"),
           "ambiguous": dict(fc="white", ec="black"), "fail": dict(fc="white", ec="#9C9C9C")}
    for xi, r in xs:
        if np.isfinite(r.chance_share):
            a.plot([xi - 0.45, xi + 0.45], [100 * r.chance_share] * 2, color="#404040", lw=0.8, solid_capstyle="butt",
                   alpha=1.0 if r.chance_exact else 0.4, zorder=1)
        st = sty[r.verdict]
        a.scatter([xi], [100 * r.share], marker="^" if r.sign_rel == "peak" else "v", s=6, linewidths=0.4,
                  facecolors=st["fc"], edgecolors=st["ec"], zorder=3)
    # the passing libraries, labeled at about their own height to the right of the fungi, joined to their points
    lab_x = centers["fungi"][1] + 2.2
    passing = sorted([(xi, r) for xi, r in xs if r.verdict in ("pass", "ambiguous")], key=lambda t: -t[1].share)
    ys, last = [], 112.0
    for xi, r in passing:
        yv = min(100 * r.share, last - 7.0)
        ys.append(yv)
        last = yv
    for (xi, r), yl in zip(passing, ys):
        nm = NICE.get(r.genome, r.genome)
        txt = (r"$\it{" + nm.replace(" ", r"\ ") + "}$" if r.genome in ITAL else nm) + SUFFIX.get(r.genome, "") + \
            (" (both signs)" if r.verdict == "ambiguous" else "")
        a.text(lab_x, yl, txt, fontsize=6, ha="left", va="center", zorder=4)
        a.plot([xi + 0.5, lab_x - 0.4], [100 * r.share, yl], color="#9A9A9A", lw=0.35, zorder=0.8)
    for g in order:
        lo, hi = centers[g]
        a.plot([lo - 0.4, hi + 0.4], [-6.0, -6.0], color="#777777", lw=0.5, clip_on=False)
        a.text((lo + hi) / 2, -8.5, lab[g], ha="center", va="top", fontsize=6, linespacing=0.95)
    for sp in ("top", "right", "bottom"):
        a.spines[sp].set_visible(False)
    a.set_xticks([])
    a.set_yticks([0, 25, 50, 75, 100])
    a.tick_params(axis="y", length=2, pad=1)
    a.set_ylabel("chromosomes where the trans\npeak meets the cis extreme (%)", labelpad=1, linespacing=0.95)
    n_by = Lb.groupby("kind").size().to_dict()
    for g in order:
        N.add(f"libs_n_{g}", int(n_by.get(g, 0)), "benchmark_cis_trans_runs16.tsv (one row per contact map)",
              f"Hi-C libraries drawn in panel E: {g}")
    N.add("libs_n", int(len(Lb)), "benchmark_cis_trans_runs16.tsv (one row per contact map)", "libraries drawn in panel E")
    rnd = Lb[Lb.kind == "randomized"]
    N.add("rand_median_share", round(100 * rnd.share.median(), 1), "benchmark_cis_trans_runs16.tsv", "randomized: median % meeting")
    N.add("rand_median_chance", round(100 * rnd.chance_share.median(), 1), "benchmark_cis_trans_runs16.tsv (chance)",
          "randomized: median chance %")
    N.add("rand_max_share", round(100 * rnd.share.max(), 1), "benchmark_cis_trans_runs16.tsv", "randomized: highest % meeting")
    for v in ("pass", "support", "fail", "ambiguous"):
        N.add(f"libs_{v}", int((Lb.verdict == v).sum()), "benchmark_cis_trans_runs16.tsv (one row per contact map)",
              f"libraries drawn in panel E at {v}")
    return a


# ============================================================================== the key over panel A
def key_row(P, y):
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    rows = [[(Line2D([], [], color=C_CIS, lw=1.0), "cis contact, obs./exp."),
             (Line2D([], [], color=C_TRANS, lw=1.0), "trans contact, obs./exp."),
             (Line2D([], [], color=C_RAND, lw=0.9), "randomized Hi-C, same axes"),
             (Line2D([], [], color="#9A9A9A", lw=0.6, ls=(0, (1, 1.5))), "1×")],
            [(Line2D([], [], marker="v", ls="none", ms=3.5, mfc=C_CIS, mec="none"), "cis extreme"),
             (Line2D([], [], marker="^", ls="none", ms=3.5, mfc=C_TRANS, mec="none"), "trans peak"),
             (Line2D([], [], color=C_MEET, lw=0.8), "the two meet"),
             (Patch(fc=C_CEN, ec="none"), "CENH3 (ChIP)"),
             (Patch(fc=C_SAT, ec="none"), "call, FASTA only"),
             (Patch(fc="black", ec="none"), "call, with Hi-C"),
             (Patch(fc=C_END, ec="none"), "chromosome ends, not read"),
             (Patch(fc="white", ec=C_HATCH, hatch="//////", lw=0), "no data")]]
    for k, items in enumerate(rows):
        a = P.ax(4, y + k * 3.6, FW - 6, 3.2)
        a.set_axis_off()
        a.legend([h for h, _ in items], [t for _, t in items], loc="upper left", ncol=len(items), frameon=False,
                 fontsize=6, handlelength=1.5, handletextpad=0.4, columnspacing=1.1, borderaxespad=0, borderpad=0)


# ============================================================================== legend pages
def IT(s):
    """Italic species name in legend text (mathtext). Its spaces become mathtext spaces with no literal space character
    ('\\;'), so the word wrapper of text_page, which splits on spaces, never breaks a name across lines."""
    return "$\\it{" + s.replace(" ", "\\;") + "}$"


def text_page(pdf, paras, fs, w_in, stamp=None):
    """A legend page as wide as the figure, height fitted to the text (wrapped by measuring rendered widths); 10 pt throughout."""
    import matplotlib.pyplot as plt
    margin = 0.3
    fig = plt.figure(figsize=(w_in, 11))
    r = fig.canvas.get_renderer()
    maxw = (w_in - 2 * margin) * 72

    def width(s):
        t = fig.text(0, 0, s, fontsize=fs)
        bb = t.get_window_extent(renderer=r)
        t.remove()
        return bb.width * 72 / fig.dpi
    lines = []
    for lab, txt in paras:
        words = txt.split(" ")
        cur, first = "", True
        prefix = f"$\\bf{{{lab}}}$  " if lab else ""
        for wd in words:
            trial = (cur + " " + wd) if cur else wd
            if width((prefix if first else "") + trial) > maxw and cur:
                lines.append((prefix if first else "") + cur)
                cur, first = wd, False
            else:
                cur = trial
        lines.append((prefix if first else "") + cur)
        lines.append("")
    if stamp:
        cur = ""
        for tok in stamp.replace("/", "/\x00").replace(" · ", " · \x00").split("\x00"):
            if cur and width(cur + tok) > maxw:
                lines.append(cur)
                cur = "    " + tok.lstrip()
            else:
                cur += tok
        lines.append(cur)
    lh = fs * 1.34 / 72
    H = margin * 2 + lh * len(lines)
    plt.close(fig)
    fig = plt.figure(figsize=(w_in, H))
    y = 1 - margin / H
    for ln in lines:
        fig.text(margin / w_in, y, ln, fontsize=fs, va="top", ha="left")
        y -= lh / H
    pdf.savefig(fig)
    plt.close(fig)


def extra_numbers(T, Cn, meta, N):
    """Numbers read from the 0.16.0 per-chromosome table and run outputs for the legend, and the check that the truth and the
    calls drawn are the ones the release scored (release16.chrom.tsv: start/end 1-based closed, as centromeres.tsv; the calls
    drawn are 0-based half-open, BED-like; truth BED 0-based half-open)."""
    import csv
    c16 = pd.read_csv(RES / "release16.chrom.tsv", sep="\t")
    n_checked = 0
    for arm, g in (("S", "tomato"), ("SH", "tomato"), ("S", "cneoG1"), ("SH", "cneoG1")):
        tr = meta["truth"][g]
        calls = meta["calls"][(arm, g)]
        t = c16[(c16.arm == arm) & (c16.genome == g)].set_index("chrom")
        assert set(t.index) == set(calls), (arm, g)
        for c, r in t.iterrows():
            s, e = truth_iv(tr, c)
            assert r.truth == f"{s}-{e}", (arm, g, c, r.truth, s, e)
            iv = calls[c]
            assert iv is not None and iv[0] + 1 == int(r.start) and iv[1] == int(r.end), (arm, g, c, iv, r.start, r.end)
            gap = min(max(0, max(iv[0], a) - min(iv[1], z)) for a, z in tr[c])
            assert abs(gap - float(r.gap)) <= 1, (arm, g, c, gap, r.gap)
            assert (r.ok == "R") == (gap <= float(r.tol)), (arm, g, c, gap, r.tol, r.ok)
            n_checked += 1
    N.add("truth_checked", n_checked, "results/docs16_rel16/release16.chrom.tsv vs V11/bench/rel/score_rel.py truth + runs/rel16 centromeres.tsv",
          "chromosome-runs whose truth interval, call coordinates, distance (within 1 bp) and right/wrong equal the release's")
    t = c16[(c16.arm == "SH") & (c16.genome == "tomato")]
    fl = t["flags"].astype(str)
    N.add("tom_moved", int(((t.relation == "offset") & fl.str.contains("cis_only")).sum()), "results/docs16_rel16/release16.chrom.tsv",
          "SH tomato: calls moved to the cis extreme (relation offset, flag cis_only)")
    N.add("tom_confirmed", int(fl.str.contains("cis_agrees").sum()), "results/docs16_rel16/release16.chrom.tsv",
          "SH tomato: FASTA calls confirmed by the cis extreme (cis_agrees)")
    ep = t[fl.str.contains("cis_end_pair")]
    N.add("tom_endpair", " and ".join(chr_label("tomato", c, T["names"]) for c in ep.chrom), "results/docs16_rel16/release16.chrom.tsv",
          "SH tomato: low calls in a terminal eighth moved to their own trans pair (cis_end_pair), new in 0.16")
    N.add("tom_endpair_right", f"{int((ep.ok == 'R').sum())}/{len(ep)}", "results/docs16_rel16/release16.chrom.tsv", "SH tomato: end-pair moves right")
    s = c16[(c16.arm == "S") & (c16.genome == "tomato") & (c16.ok == "W")]
    N.add("tom_S_gap_mb", f"{s.gap.min() / 1e6:.1f}–{s.gap.max() / 1e6:.0f}", "results/docs16_rel16/release16.chrom.tsv",
          "S tomato: wrong FASTA-only calls' distance to CENH3 (Mb)")
    N.add("tom_conf_SH", ", ".join(f"{v} {k}" for k, v in t.confidence.value_counts().items()), "results/docs16_rel16/release16.chrom.tsv",
          "SH tomato: confidence of the calls")
    ex = N["ex_chrom"]
    gS = c16[(c16.arm == "S") & (c16.genome == "tomato") & (c16.chrom == ex)].iloc[0]
    N.add("ex_old_gap_mb", round(gS.gap / 1e6, 1), "results/docs16_rel16/release16.chrom.tsv", "example: the FASTA-only call to CENH3 (Mb)")
    N.add("ex_S_ok", gS.ok, "results/docs16_rel16/release16.chrom.tsv", "example: the FASTA-only call right (R) or wrong (W)")
    gH = c16[(c16.arm == "SH") & (c16.genome == "tomato") & (c16.chrom == ex)].iloc[0]
    N.add("ex_SH", f"{gH['ok']} {gH['confidence']} {gH['relation']} {gH['flags']}", "results/docs16_rel16/release16.chrom.tsv",
          "example: the call with Hi-C")
    hubs = [chr_label("tomato", c, T["names"]) for c, f in T["peak_flags"][("int", 1)].items() if "hub" in str(f)
            and c not in T["res"][("int", 1)]["peaks"]]
    N.add("tom_hub_misses", " and ".join(hubs), "benchmark_cis_trans_cache_tomato.pkl", "tomato: trans maximum on a contact hub")
    f = RUNS / "SH/cneoG1/cneoG1.centromeres.tsv"
    rows = list(csv.DictReader(open(f), delimiter="\t"))
    ends = sum(1 for r in rows if r["hic_denovo_peak_start"] not in (".", "") and
               min(int(r["hic_denovo_peak_start"]), int(r["chrom_len"]) - int(r["hic_denovo_peak_start"])) <
               END_ZONE * int(r["chrom_len"]))
    N.add("cneo_denovo_ends", f"{ends}/{len(rows)}", "runs/rel16/SH/cneoG1/cneoG1.centromeres.tsv",
          "C. neoformans G1: prior-free Hi-C peaks in a terminal eighth")


def legends(N, ex_label):
    """One legend page: the result, then one or two sentences per panel (biology only; statistics and rules in memo.md)."""
    tom, cn, aet = "tomato", IT("C. neoformans"), IT("Ae. tauschii")
    k_tom, n_tom = N["tom_kn"].split("/")
    k_cn, n_cn = N["cneo_kn"].split("/")
    rS, rH = N["tom_rightS"].split("/")[0], N["tom_rightSH"].split("/")[0]
    cS, cH = N["cneo_rightS"].split("/")[0], N["cneo_rightSH"].split("/")[0]
    pct = N.add("ex_chip_sat_pct", round(100 * N["ex_chip_sat"] / N["ex_chip_max"]), "benchmark_cis_trans_cache_chip.pkl",
                "example: CENH3 ChIP/input over the FASTA-only call, % of the peak in the window")
    low_tom = int(N["tom_conf_SH"].split(" low")[0].split(", ")[-1]) if " low" in N["tom_conf_SH"] else 0
    num = {0: "none", 1: "one", 2: "two", 3: "three", 4: "four", 5: "five", 6: "six", 7: "seven", 8: "eight", 9: "nine", 10: "ten",
           11: "eleven", 12: "twelve", 13: "thirteen", 14: "fourteen"}
    paras = [
        ("", "Two readouts of one Hi-C library find centromeres that no repeat marks: where the contact across a stretch of "
             "chromosome peaks and the contact with the other chromosomes' peaks rises at the same place, that place is the "
             "centromere. No annotation is needed."),
        ("A", f"In {tom}, the contact across each chromosome's interior peaks once (purple), and the contact with the other "
              f"chromosomes' peaks (blue) rises at the same place, on CENH3, on {k_tom} of {n_tom} chromosomes "
              f"({N['tom_chance_sum']:.2f} expected by chance); on {N['tom_hub_misses']} the blue maximum sits on a contact "
              f"hub. A randomized library (gray) shows neither. With Hi-C, all {rH} calls land on CENH3 ({low_tom} at low "
              f"confidence), against {rS} from the FASTA alone; on {N['tom_endpair']}, whose FASTA calls lay at a chromosome "
              f"end, the call moves to where that chromosome's own readouts meet (a rule fitted on this genome, so here it is not an independent test). "
              f"Carets, values beyond the axis."),
        ("B", f"On tomato {ex_label} (boxed in A), a typical chromosome, the two flanks of the centromere contact each other more than "
              f"expected (diamond), and the same bins contact the centromeric peaks of all {int(n_tom) - 1} other "
              f"chromosomes (blue rows): both readouts mark one CENH3 domain. The FASTA-only call sits on a satellite array "
              f"{N['ex_old_gap_mb']:.1f} Mb away, whose CENH3 signal is {pct:.0f}% of the peak's."),
        ("C", "Schematic: centromeres fold locally (cis) and cluster in the nucleus (trans). Chromosome ends do both too, so "
              "only the interior of each chromosome is read."),
        ("D", f"In {aet}, read with the ends kept, both readouts peak at chromosome ends and agree there on "
              f"{N['aet_whole_kn'].split('/')[0]} of {N['aet_whole_kn'].split('/')[1]} chromosomes, away from CENH3; "
              f"with the ends excluded, {num.get(int(N['aet_int_kn'].split('/')[0]), N['aet_int_kn'].split('/')[0])} agree."),
        ("E", f"Each of {N['libs_n']} Hi-C libraries is one point (bar: the share expected by chance): black where the two "
              f"readouts meet on more than half of its chromosomes, beyond chance, so the line acts (pass); gray where they "
              f"meet beyond chance on half or fewer (support); open where they do not beat chance (fail). They agree on "
              f"most chromosomes, far beyond chance, in {tom} (Heinz 1706 and TS60), {IT('Brachypodium')}, {cn}, "
              f"{IT('N. crassa')} and {IT('V. dahliae')}; {IT('S. pombe')} agrees at both extremes, so no call is made. "
              f"They agree beyond chance on a minority in {N['tomx_support']} of {N['tomx_n']} further tomato accessions, "
              f"barley, {cn} in anaphase and rhesus. Randomized libraries never agree on a majority (on a minority beyond "
              f"chance in {N['n_rand_support']} of {N['n_rand']}), and holocentric genomes fail. In potato the dips sit on "
              f"CENH3 on {N['pot_dip_right'].split('/')[0]} of {N['pot_dip_right'].split('/')[1]} chromosomes, but trans "
              f"contacts do not confirm them."),
        ("F", f"In {cn} (G1 Hi-C, 10-kb bins), the readouts meet on {k_cn} of {n_cn} chromosomes "
              f"({N['cneo_chance_sum']:.2f} expected by chance). With Hi-C, {cH} of {n_cn} calls lie on CENP-A at medium "
              f"confidence and one ({N['cneo_wrong_chrom']}) lies {N['cneo_wrong_gap_kb']} kb beyond it; from the FASTA "
              f"alone all {cS} lie on it, at low confidence."),
    ]
    return paras


# ============================================================================== build
def build():
    plt = setup()
    from matplotlib.backends.backend_pdf import PdfPages
    from matplotlib.patches import ConnectionPatch, Rectangle
    T, R, A, Cn, Pt = load("tomato"), load("tneg"), load("aet"), load("cneo"), load("pot")
    meta = load("meta")
    assert meta["cenfuse"][0] == "0.16.0" and all(x["run"].startswith(str(RUNS)) for x in (T, R, A, Cn, Pt)), "not 0.16.0"
    chip = load("chip")
    N = collect_numbers(T, R, A, Cn, Pt, meta)
    ex, cand = pick_example(T, meta, N)
    N.add("ex_chrom", ex, "benchmark_cis_trans.py pick_example", "the example chromosome's sequence name")
    check_example(T, meta, ex, N)
    extra_numbers(T, Cn, meta, N)
    Lb = libraries(N)
    P = Page(plt)
    # ---- A
    key_row(P, 0.8)
    yA = 13.5
    axesA, mbmmA, rowh = panel_A(P, T, R, meta, N, yA, ex)
    P.letter(1.0, yA - 1.0, "A")
    # ---- B (example) and C (schematic)
    yB = 114.5
    xB, wB = 27.0, 91.0
    bB = panel_B(P, T, R, meta, chip, N, ex, xB, yB, wB)
    P.letter(1.0, yB - 8.6, "B")
    # zoom box on the example row in A, joined to B
    ac, at, ak, x0, yy, w = axesA[ex]
    xl = bB["xlim"]
    zx0, zx1 = x0 + xl[0] * mbmmA, x0 + xl[1] * mbmmA
    bh = 4.4 + 0.5 + 4.4 + 0.5 + 2.0 + 1.6
    box = P.ax(zx0, yy - 1.0, zx1 - zx0, bh)
    box.set_xlim(0, 1)
    box.set_ylim(0, 1)
    box.patch.set_alpha(0)
    for s in box.spines.values():
        s.set_linewidth(0.6)
    box.set_xticks([])
    box.set_yticks([])
    if ex == T["names"][5]:                               # the example on A's bottom row: connectors run straight down to B
        for (fx, bx) in ((zx0, xB), (zx1, xB + wB)):
            P.fig.add_artist(ConnectionPatch((fx / FW, 1 - (yy - 1.0 + bh) / FH),
                                             (bx / FW, 1 - (yB - 0.4) / FH), "figure fraction", "figure fraction",
                                             color="#8A8A8A", lw=0.5, zorder=0))
    else:                                                 # higher rows: connectors would cross the rows below; tag the box
        P.text(zx1 + 0.6, yy - 1.0, "B", fontsize=7, fontweight="bold", ha="left", va="top")
    xC = 127.0
    panel_C(P, xC, yB - 4.0, FW - xC - 1.0, 40.0)
    P.letter(xC - 4.0, yB - 8.6, "C")
    # ---- D, E, F
    y3 = 177.5
    rowh3 = 7.3
    panel_D(P, A, meta, N, 10.0, y3, 40.0, rowh3)
    P.letter(1.0, y3 - 1.0, "D")
    panel_E(P, Lb, N, 68.0, y3, 52.0, 47.0)
    P.letter(55.0, y3 - 1.0, "E")
    panel_F(P, Cn, meta, N, (134.0, 161.5), y3, 19.0, rowh3)
    P.letter(124.0, y3 - 1.0, "F")
    paras = legends(N, N["ex_label"])
    N.write(HERE / f"{STEM}_numbers.tsv")
    stamp = (f"Script: docs/img/src/{STEM}.py · data: the cenfuse 0.16.0 benchmark runs (not in the repository) · "
             f"{datetime.date.today().isoformat()}")
    meta_pdf = {"Title": "benchmark_cis_trans: two readouts of one Hi-C library test each other", "Author": None,
                "Subject": None, "Keywords": None, "Creator": None, "Producer": None}
    OUT.mkdir(parents=True, exist_ok=True)
    P.fig.savefig(OUT / f"{STEM}.png", dpi=300)
    with PdfPages(OUT / f"{STEM}.pdf", metadata=meta_pdf) as pdf:
        pdf.savefig(P.fig, dpi=300)
        text_page(pdf, paras, 10, FW / 25.4, stamp)
    plt.close(P.fig)
    print(f"wrote {OUT / STEM}.pdf, .png and {HERE / STEM}_numbers.tsv; example {N['ex_label']} ({ex})", flush=True)


def main():
    global VERBOSE
    ap = argparse.ArgumentParser(description="Build benchmark_cis_trans (PDF: figure + one legend page; PNG; numbers TSV).")
    ap.add_argument("-v", "--verbose", action="store_true", help="print per-step progress and the example ranking")
    VERBOSE = ap.parse_args().verbose
    build()


if __name__ == "__main__":
    main()
