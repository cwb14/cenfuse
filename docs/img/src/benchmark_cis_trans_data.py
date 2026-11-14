#!/usr/bin/env python3
"""Matrix-derived values for benchmark_cis_trans, cached in the work folder
($CENFUSE_FIG_WORK, default: the current directory) as benchmark_cis_trans_cache_<entry>.pkl.

usage: benchmark_cis_trans_data.py [--only tomato tneg aet cneo pot chip meta]
For each library drawn (tomato and its randomized partner, Ae. tauschii, C. neoformans G1) the stored contact matrix of the
cenfuse 0.16.0 release run (V16/runs/rel16/<arm>/<genome>) is loaded exactly as bench/cis_trans_anchor.py does (pickle of <run>/work/hic/<g>.pkl, the cached npz, Bins/coarsen, ice, balanced,
block_means) and read with cenfuse 0.16.0's own functions (V16/code/rel16): hicfuse.cis_profile, interior, anchored_on_cores, cis_trans_line, and
hic.score.partner_oe / cluster_score for the trans score itself (anchored_on_cores computes it but does not return it). No Hi-C
engine is rerun. meta: the truth (score_rel G + read_truth), the 0.16.0 calls from the FASTA alone (arm S) and with Hi-C
(arm SH), and a check that the anchor tables of 0.13 (results/cis_trans_anchor*.tsv) read the same matrix and Hi-C state as
the 0.16.0 runs. Also bins the tomato CENH3 ChIP and input fragments per 100 kb as the held-out truth rule does
(v2/scripts/chip_truth.py: CPM + 0.5 ratio, running mean over 3 bins)."""
import argparse
import os
import time
import pickle
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(os.environ.get("CENFUSE_FIG_WORK") or Path.cwd()).resolve()   # work folder: caches and tables
SRC = Path(__file__).resolve().parent                                      # docs/img/src
DATA = Path(os.environ.get("CENFUSE_BENCH_DATA") or "/data2/chris/centromere_benchmark")   # as benchmark/paths.py
V = DATA / "v2/v3/v4/v5/v6/v7/v8/v9/v10/v10.1/v10.2/v10.3/v10.4/v11/v12"
V11 = V.parent
V16 = V / "v13/v14/v15/v16"
CODE = V16 / "code/rel16"
V8 = DATA / "v2/v3/v4/v5/v6/v7/v8"
sys.path.insert(0, str(V11 / "bench/rel"))
import score_rel as SR  # noqa: E402  (frozen truth rules; imports an old cenfuse for one util)
for _k in [k for k in sys.modules if k.startswith("cenfuse")]:
    del sys.modules[_k]                                   # then use the 0.16.0 code for everything else
sys.path.insert(0, str(CODE))
import scipy.sparse as sp  # noqa: E402
sys.path.insert(0, str(SRC))
import benchmark_cis_trans_npz as npzfast  # noqa: E402
import cenfuse  # noqa: E402
from cenfuse import hicfuse  # noqa: E402
assert Path(cenfuse.__file__).resolve().is_relative_to(CODE) and cenfuse.__version__ == "0.16.0", (cenfuse.__file__,
                                                                                                 cenfuse.__version__)
from cenfuse.hic import params as hp  # noqa: E402
from cenfuse.hic.balance import balanced, ice  # noqa: E402
from cenfuse.hic.bins import Bins, coarsen  # noqa: E402
from cenfuse.hic.score import block_means, cluster_score, partner_oe  # noqa: E402

MAXD = 24                    # triangle: pairs up to 24 bins apart (6 Mb at 250 kb)
CHIP_BIN = 100_000           # the held-out truth rule's bin


def run_dir(arm, g):
    return V16 / "runs/rel16" / arm / g


def run_dir13(arm, g):
    d = V / "runs/rel13f" / arm / g
    return d if (d / "DONE").exists() else V / "runs/rel13" / arm / g


def same_state(arm, g):
    """The 0.13 run read by results/cis_trans_anchor*.tsv and the 0.16.0 run hold the same Hi-C state: the same matrix file
    (inode) and the same bins, base size, members and hub in the pickle (the anchor tables read nothing else from a run)."""
    out = {}
    hs = [pickle.load(open(d / "work/hic" / f"{g}.pkl", "rb")) for d in (run_dir13(arm, g), run_dir(arm, g))]
    npz = [d / "hic" / f"{g}.cache" / Path(h["matrix_cache"]).name for d, h in zip((run_dir13(arm, g), run_dir(arm, g)), hs)]
    out["npz_same_file"] = os.stat(npz[0]).st_ino == os.stat(npz[1]).st_ino and npz[0].name == npz[1].name
    a, b = (h["bins"] for h in hs)
    out["bins_same"] = list(a.names) == list(b.names) and np.array_equal(a.lengths, b.lengths) and a.size == b.size
    out["base_same"] = int(hs[0]["base_size"]) == int(hs[1]["base_size"])
    out["members_same"] = np.array_equal(np.asarray(hs[0]["members"]), np.asarray(hs[1]["members"]))
    out["hub_same"] = np.array_equal(np.asarray(hs[0]["hub"]), np.asarray(hs[1]["hub"]))
    out["npz"] = npz[1].name
    return out


def load_state(p, gid):
    h = pickle.load(open(p / "work/hic" / f"{gid}.pkl", "rb"))
    b = h["bins"]
    seqs = {c: int(L) for c, L in zip(b.names, b.lengths)}
    npz = p / "hic" / f"{gid}.cache" / Path(h["matrix_cache"]).name
    M = npzfast.load_npz(npz).tocsr()
    base = Bins(seqs, int(h["base_size"]))
    if b.size != base.size:
        M = coarsen(M, base, b)
    bias, good, _, _ = ice(M)
    W = balanced(M, bias, good).tocsr()
    members = h["members"]
    mean = block_means(W, b, good, members)
    st = dict(W=W, bins=b, good=good, mean=mean, members=members, hub=h["hub"], rdna={}, bias=bias)
    return st, h, npz


def trans_score(st, anchors):
    """The trans score of anchored_on_cores for anchors {chrom: (start, end)}: per bin, the mean over the OTHER chromosomes y
    of the O/E of its contact with y's partner window (visible_window at y's anchor). Returns (S over all bins, oe per y)."""
    b, good, members = st["bins"], st["good"], st["members"]
    idx = {n: i for i, n in enumerate(b.names)}
    win = {idx[c]: hicfuse.visible_window(b, good, idx[c], *iv) for c, iv in anchors.items() if c in idx and members[idx[c]]}
    oe = partner_oe(st["W"], b, good, st["mean"], {k: w for k, w in win.items() if len(w)})
    return cluster_score(oe, b, good, members), oe, win


def flags_at(st, k, j):
    hub = st.get("hub")
    hub = np.zeros(st["bins"].N, dtype=bool) if hub is None else hub
    fl = hicfuse.call_flags(k, dict(peak=j, core_lo=j, core_hi=j), st["bins"], hub, {}, True, 0.0, True, good=st["good"])
    return ",".join(fl)


def oe_dense(W, b, good, k, maxd=MAXD):
    """Observed/expected contact of chromosome k up to maxd bins apart (expected: its mean balanced contact at each distance
    over mappable pairs, as cis_profile's expectation). NaN where either bin is masked."""
    sl = b.sl(k)
    A = W[sl, :][:, sl].toarray()
    g = good[sl]
    nb = len(g)
    X = np.full((nb, nb), np.nan, dtype=np.float32)
    for d in range(1, min(nb, maxd + 1)):
        x = np.diag(A, d)
        m = g[:nb - d] & g[d:]
        e = x[m].mean() if m.any() else np.nan
        if not e > 0:
            continue
        i = np.arange(nb - d)
        X[i, i + d] = np.where(m, x / e, np.nan)
    return X


def read_genome(arm, gid, whole=False, maps=False):
    p = run_dir(arm, gid)
    t0 = time.time()
    st, h, npz = load_state(p, gid)
    print(gid, arm, "loaded", round(time.time() - t0, 1), "s", flush=True)
    b, W, good, members = st["bins"], st["W"], st["good"], st["members"]
    names = list(b.names)
    out = dict(arm=arm, gid=gid, run=str(p), npz=str(npz), bin=int(b.size), names=names,
               lengths={c: int(L) for c, L in zip(b.names, b.lengths)}, nb={c: int(n) for c, n in zip(b.names, b.nb)},
               good={c: good[b.sl(k)].copy() for k, c in enumerate(names)},
               hub={c: np.asarray(st["hub"])[b.sl(k)].copy() for k, c in enumerate(names)},
               members={c: bool(members[k]) for k, c in enumerate(names)}, stored_cis=h.get("cis"))
    prof, ext = {}, {"int": {+1: {}, -1: {}}, "all": {+1: {}, -1: {}}}
    for k, c in enumerate(names):
        if not members[k]:
            continue
        v = hicfuse.cis_profile(W, b, good, k)
        prof[c] = v
        if not np.isfinite(v).any():
            continue
        vi = hicfuse.interior(v)
        for g in (+1, -1):
            if np.isfinite(vi).any():
                ext["int"][g][c] = int(np.nanargmax(g * vi))
            ext["all"][g][c] = int(np.nanargmax(g * v))
    out["prof"], out["ext"] = prof, ext
    modes = ("int", "all") if whole else ("int",)
    out["S"], out["res"], out["trans_peak"], out["peak_flags"] = {}, {}, {}, {}
    for mode in modes:
        for g in (+1, -1):
            anchors = {}
            for c, j in ext[mode][g].items():
                k = names.index(c)
                o = int(b.off[k]) + j
                anchors[c] = (int(b.start[o]), int(b.end[o]))
            S, oe, win = trans_score(st, anchors)
            r = hicfuse.anchored_on_cores(st, anchors, label=f"{mode} {g}")
            key = (mode, g)
            out["S"][key] = {c: S[b.sl(k)].copy() for k, c in enumerate(names)}
            out["res"][key] = r
            tp, fl = {}, {}
            for c in anchors:
                k = names.index(c)
                v = S[b.sl(k)]
                if np.isfinite(v).any():
                    tp[c] = int(np.nanargmax(v))
                    fl[c] = flags_at(st, k, tp[c])
            out["trans_peak"][key], out["peak_flags"][key] = tp, fl
            if maps and mode == "int" and g == +1:      # every partner's O/E row, whole genome (any example can be drawn)
                out["partner_oe"] = {names[y]: oe[y].astype(np.float32) for y in oe}
                out["partner_windows"] = {names[y]: (int(w.min() - b.off[y]), int(w.max() - b.off[y])) for y, w in win.items()
                                          if len(w)}
    out["line"] = hicfuse.cis_trans_line(st)
    if maps:
        out["oe"] = {c: oe_dense(W, b, good, k) for k, c in enumerate(names) if members[k]}
    return out


def chip_track(fai, chip, inp, chroms):
    """CENH3 ChIP / input per 100-kb bin, as v2/scripts/chip_truth.py (fragments, any MAPQ; CPM + 0.5 ratio; 3-bin mean)."""
    lengths = {}
    for line in open(fai):
        p = line.split("\t")
        if p[0] in chroms:
            lengths[p[0]] = int(p[1])

    def counts(path):
        """Fragments per 100-kb bin (positions are 1-based fragment positions), streamed with awk (pandas is ~30x slower)."""
        import subprocess
        c = {k: np.zeros(int(np.ceil(L / CHIP_BIN)), dtype=np.int64) for k, L in lengths.items()}
        cmd = f"zcat {path} | awk -F'\\t' '{{n[$1\"\\t\"int(($2-1)/{CHIP_BIN})]++}} END{{for(k in n) print k\"\\t\"n[k]}}'"
        out = subprocess.run(["bash", "-c", cmd], capture_output=True, text=True, check=True).stdout
        for line in out.splitlines():
            ch, b, n = line.split("\t")
            if ch in c:
                c[ch][min(int(b), len(c[ch]) - 1)] += int(n)
        return c
    C, I = counts(chip), counts(inp)
    tot, itot = sum(v.sum() for v in C.values()), sum(v.sum() for v in I.values())
    out = {}
    for k in lengths:
        e = (C[k] / tot * 1e6 + 0.5) / (I[k] / itot * 1e6 + 0.5)
        out[k] = pd.Series(e).rolling(3, center=True, min_periods=1).mean().to_numpy()
    return out


def calls(path):
    c = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    return {r["chrom"]: (int(r["start"]) - 1, int(r["end"]), r["confidence"], r["relation"], r["flags"])
            if r["start"] not in (".", "") else None for r in c.to_dict("records")}


def truth(gid):
    kind, path, rule, tol, nb = SR.G[gid]
    return SR.read_truth(path)


def main():
    ap = argparse.ArgumentParser(description="Cache the matrix-derived values of benchmark_cis_trans (one pickle per entry).")
    ap.add_argument("--only", nargs="*", default=["tomato", "tneg", "aet", "cneo", "pot", "chip", "meta"],
                    help="entries to compute: tomato tneg aet cneo pot chip meta (default: all)")
    a = ap.parse_args()

    def save(name, obj):
        f = HERE / f"benchmark_cis_trans_cache_{name}.pkl"
        tmp = f.with_suffix(".tmp")
        pickle.dump(obj, open(tmp, "wb"))
        os.replace(tmp, f)
        print(name, "written", flush=True)
    jobs = dict(tomato=("SH", "tomato", True), tneg=("TNEG", "tomato", True), aet=("SH", "aetAL878", False),
                cneo=("SH", "cneoG1", False), pot=("SH", "potato", False))
    for name in a.only:
        if name in jobs:
            arm, gid, maps = jobs[name]
            r = read_genome(arm, gid, whole=True, maps=maps)
            print(name, "cis_trans_line:", r["line"]["k"], r["line"]["n"], r["line"]["status"],
                  "stored:", (r["stored_cis"] or {}).get("k"), (r["stored_cis"] or {}).get("status"), flush=True)
            save(name, r)
        elif name == "chip":
            d = V8 / "data/tomato_scau/slH1706"
            save("chip", chip_track(d / "slH1706.fa.fai", d / "chip.pos.gz", d / "input.pos.gz",
                                    [f"CM1750{i}.1" for i in range(67, 79)]))
        elif name == "meta":
            save("meta", dict(truth={g: truth(g) for g in ("tomato", "aetAL878", "cneoG1", "potato")},
                              calls={(arm, g): calls(run_dir(arm, g) / f"{g}.centromeres.tsv")
                                     for arm in ("S", "SH") for g in ("tomato", "cneoG1", "aetAL878")},
                              same_state={(arm, g): same_state(arm, g) for arm, g in
                                          (("SH", "tomato"), ("TNEG", "tomato"), ("SH", "aetAL878"), ("SH", "potato"))},
                              cenfuse=(cenfuse.__version__, str(Path(cenfuse.__file__).resolve()))))


if __name__ == "__main__":
    main()
