#!/usr/bin/env python3
"""Every Hi-C library's cis-trans test in the cenfuse 0.16.0 release, with its chance (the expected number of chromosomes whose
trans peak meets the cis extreme), for panel E of benchmark_cis_trans.

usage:
  benchmark_cis_trans_libs.py table              build benchmark_cis_trans_runs16.tsv (in the work folder): one row per 0.16.0
                                                 release run whose cenfuse_summary.tsv cis_line says "tested by trans contacts"
                                                 (runs from bench/jobs_rel16.tsv + jobs_rel16_extra*.tsv), with its k, n, P,
                                                 verdict, sign and contact-matrix hash, and the exact chance reused from the 0.13
                                                 figure's re-read where the same matrix gave the same k, n, verdict and sign
  benchmark_cis_trans_libs.py read OUT [--reverse OTHER]
                                                 re-read the test from the stored matrix (cenfuse 0.16.0's own functions:
                                                 cis_profile, interior, anchored_on_cores, cis_trans_line) for every run of the
                                                 table without a reused chance; resumable (rows already in OUT or OTHER are
                                                 skipped); --reverse takes the jobs from the end and stops where OTHER has begun
Loads each run as bench/cis_trans_anchor.py does (pickle of <run>/work/hic/<g>.pkl, the cached npz, Bins/coarsen, ice, balanced,
block_means). One process, no pool (at most 2 python processes run for this figure). No Hi-C engine is rerun."""
import argparse
import os
import pickle
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(os.environ.get("CENFUSE_FIG_WORK") or Path.cwd()).resolve()   # work folder: caches and tables
SRC = Path(__file__).resolve().parent                                      # docs/img/src
DATA = Path(os.environ.get("CENFUSE_BENCH_DATA") or "/data2/chris/centromere_benchmark")   # as benchmark/paths.py
V12 = DATA / "v2/v3/v4/v5/v6/v7/v8/v9/v10/v10.1/v10.2/v10.3/v10.4/v11/v12"
V16 = V12 / "v13/v14/v15/v16"
CODE = V16 / "code/rel16"
RUNS = V16 / "runs/rel16"
OLD = V12 / "figures/fig_cis_trans"          # the 0.13 figure: its libraries re-read with their chance (read only)
TABLE = HERE / "benchmark_cis_trans_runs16.tsv"
sys.path.insert(0, str(CODE))
sys.path.insert(0, str(SRC))


def run_dir(arm, g):
    return RUNS / arm / g


def npz_name(arm, g):
    """The run's contact-matrix cache (hash-named npz), read from its pickle."""
    h = pickle.load(open(run_dir(arm, g) / "work/hic" / f"{g}.pkl", "rb"))
    return Path(h["matrix_cache"]).name


def parse_line(s):
    """k, n, P, verdict, sign, acts from a cis_line 'trans peaks at the cis <sign> on k/n chromosomes (P p): <verdict>, tested
    by trans contacts; acts on a'."""
    m = re.match(r"trans peaks at the cis (peak|dip) on (\d+)/(\d+) chromosomes \(P ([^)]+)\): (\w+), tested by trans contacts; "
                 r"acts on (\d+)", s)
    if not m:
        raise ValueError(f"cis_line not understood: {s}")
    sign, k, n, p, verdict, acts = m.groups()
    return dict(sign=sign, k=int(k), n=int(n), p=float(p), verdict=verdict, acts=int(acts))


def build_table():
    jobs = pd.concat([pd.read_csv(V16 / f"bench/{f}", sep="\t", header=None, usecols=[0, 1], names=["arm", "genome"])
                      for f in ("jobs_rel16.tsv", "jobs_rel16_extra1.tsv", "jobs_rel16_extra2.tsv")]).drop_duplicates()
    rows = []
    for r in jobs.itertuples(index=False):
        f = run_dir(r.arm, r.genome) / "cenfuse_summary.tsv"
        if not f.exists() or not (run_dir(r.arm, r.genome) / "DONE").exists():
            print("no finished run:", r.arm, r.genome, file=sys.stderr)
            continue
        s = pd.read_csv(f, sep="\t", dtype=str, keep_default_na=False).iloc[0]
        line = s.get("cis_line", ".")
        if "tested by trans contacts" not in str(line):
            continue
        rows.append(dict(arm=r.arm, genome=r.genome, cis_line=line, **parse_line(line), hic_bin=s.get("hic_bin", "."),
                         npz=npz_name(r.arm, r.genome)))
    t = pd.DataFrame(rows)
    old = pd.concat([pd.read_csv(OLD / f, sep="\t") for f in ("fig_cis_trans_libs.tsv", "fig_cis_trans_libs_b.tsv")],
                    ignore_index=True)
    old = old[old.get("error").isna()] if "error" in old else old
    old = old.drop_duplicates(["arm", "genome"]).set_index(["arm", "genome"])
    t["expected"], t["chance_from"] = np.nan, ""
    for ix, r in t.iterrows():
        key = (r.arm, r.genome)
        if key not in old.index:
            t.at[ix, "chance_from"] = "new: re-read"
            continue
        o = old.loc[key]
        same = (o.npz == r.npz and int(o.k) == r.k and int(o.n) == r.n and o.status == r.verdict and o.sign == r.sign)
        if same:
            t.at[ix, "expected"] = float(o.expected)
            t.at[ix, "chance_from"] = "0.13 re-read (same matrix, same k/n/verdict/sign)"
        else:
            t.at[ix, "chance_from"] = (f"changed: re-read (0.13 {o.npz} {int(o.k)}/{int(o.n)} {o.status} {o.sign}; "
                                       f"0.16 {r.npz} {r.k}/{r.n} {r.verdict} {r.sign})")
    t.to_csv(TABLE, sep="\t", index=False)
    print(f"{len(t)} runs tested by trans contacts; chance reused for {int(t.expected.notna().sum())}; "
          f"to re-read: {int(t.expected.isna().sum())}", flush=True)
    print(t[t.expected.isna()][["arm", "genome", "chance_from"]].to_string(index=False), flush=True)


def load_state(p, gid):
    import benchmark_cis_trans_npz as npzfast
    from cenfuse.hic.balance import balanced, ice
    from cenfuse.hic.bins import Bins, coarsen
    from cenfuse.hic.score import block_means
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
    state = dict(W=W, bins=b, good=good, mean=mean, members=members, hub=h["hub"], rdna={}, bias=bias)
    return state, h, npz


def one(arm, gid):
    import cenfuse
    from cenfuse import hicfuse
    from cenfuse.util import poisson_binomial_sf
    assert Path(cenfuse.__file__).resolve().is_relative_to(CODE), cenfuse.__file__
    p = run_dir(arm, gid)
    try:
        st, h, npz = load_state(p, gid)
    except Exception as e:                               # noqa: BLE001
        return dict(arm=arm, genome=gid, error=str(e)[:120])
    b, W, good, members = st["bins"], st["W"], st["good"], st["members"]
    ext = {+1: {}, -1: {}}
    for k, c in enumerate(b.names):                      # as cis_trans_line: interior extremes per sign
        if not members[k]:
            continue
        v = hicfuse.interior(hicfuse.cis_profile(W, b, good, k))
        if not np.isfinite(v).any():
            continue
        off = int(b.off[k])
        for g in (+1, -1):
            j = int(np.nanargmax(g * v))
            ext[g][c] = (int(b.start[off + j]), int(b.end[off + j]))
    row = dict(arm=arm, genome=gid, run=str(p), npz=npz.name, bin=int(b.size), chromosomes=int(members.sum()),
               cenfuse=cenfuse.__version__)
    for g, nm in ((+1, "peak"), (-1, "dip")):
        r = hicfuse.anchored_on_cores(st, ext[g], label=f"cis {nm}")
        ch = r["chance"]
        m, n = len(r["peaks"]), len(ch)
        pb = min(1.0, 2 * poisson_binomial_sf(m, list(ch.values()))) if n else 1.0
        row.update({f"{nm}_k": m, f"{nm}_n": n, f"{nm}_expected": float(sum(ch.values())), f"{nm}_p2": pb})
    line = hicfuse.cis_trans_line(st)                    # the release's own verdict, re-read
    row.update(sign="peak" if line["sign"] > 0 else "dip", k=line["k"], n=line["n"], p=line["p"], status=line["status"],
               acts=len(line["peaks"]))
    row["expected"] = row[f"{row['sign']}_expected"]
    return row


def done_in(*fs):
    d = set()
    for f in fs:
        if f and os.path.exists(f):
            x = pd.read_csv(f, sep="\t")
            d |= set(zip(x.arm, x.genome))
    return d


def read(out, other=None):
    t = pd.read_csv(TABLE, sep="\t")
    jobs = list(dict.fromkeys(zip(t[t.expected.isna()].arm, t[t.expected.isna()].genome)))
    if other:
        jobs = jobs[::-1]
    print(f"{len(jobs)} libraries to re-read; {len(done_in(out, other))} done", flush=True)
    for j in jobs:
        if other and j in done_in(other):                # the forward reader got here first: the two have met
            break
        if j in done_in(out, other):
            continue
        row = one(*j)
        pd.DataFrame([row]).to_csv(out, sep="\t", index=False, mode="a", header=not os.path.exists(out))
        print(row.get("arm"), row.get("genome"), row.get("sign"), row.get("k"), row.get("n"), row.get("status"),
              row.get("error", ""), flush=True)
    print("done", flush=True)


def main():
    ap = argparse.ArgumentParser(description="Panel E of benchmark_cis_trans: every 0.16.0 Hi-C library's cis-trans test and chance.")
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("table", help="build benchmark_cis_trans_runs16.tsv from the release summaries")
    rd = sub.add_parser("read", help="re-read the libraries without a reused chance from their matrices (resumable)")
    rd.add_argument("out", help="output TSV (appended)")
    rd.add_argument("--reverse", metavar="OTHER", help="take jobs from the end; stop where OTHER (the forward reader's TSV) has begun")
    a = ap.parse_args()
    if a.cmd == "table":
        build_table()
    else:
        read(a.out, a.reverse)


if __name__ == "__main__":
    main()
