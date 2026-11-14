"""Approximate chance of the cis-trans test without the contact matrix: per chromosome, the 9-bin zone (cis extreme +- 4 bins)
over the chromosome's mappable bins (ICE's first rule on the stored coverage: >= 0.2 x the median of nonzero bins). Used only for
libraries whose matrix was not re-read (none in the 0.16.0 figure: every chance is exact; kept so a rerun on a partial re-read
still draws every library)."""
import os
import pickle
from pathlib import Path

import numpy as np

DATA = Path(os.environ.get("CENFUSE_BENCH_DATA") or "/data2/chris/centromere_benchmark")   # as benchmark/paths.py
V16 = DATA / "v2/v3/v4/v5/v6/v7/v8/v9/v10/v10.1/v10.2/v10.3/v10.4/v11/v12/v13/v14/v15/v16"


def run_dir(arm, g):
    return V16 / "runs/rel16" / arm / g


def approx_expected(arm, g, min_frac=0.2, zone=9):
    h = pickle.load(open(run_dir(arm, g) / "work/hic" / f"{g}.pkl", "rb"))
    b, cov, members = h["bins"], np.asarray(h["coverage"], float), h["members"]
    if len(cov) != b.N:
        return np.nan, np.nan
    pos = cov[cov > 0]
    good = cov >= min_frac * np.median(pos)
    exp, n = 0.0, 0
    for k in range(len(b.names)):
        if not members[k]:
            continue
        gk = int(good[b.sl(k)].sum())
        if gk <= 2 * 4:
            continue
        exp += min(1.0, zone / gk)
        n += 1
    return exp, n
