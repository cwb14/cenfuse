#!/usr/bin/env python3
"""A release's runs against a reference release, arm by arm, with the frozen scorer: the arms whose counts changed and
every call that moved or changed confidence.

  python -m benchmark.score_dev RUNS JOBS OUT_PREFIX [--ref ROOT] [--redo ROOT] [--manifest truth.tsv]
RUNS: the release root (RUNS/<arm>/<id>). --ref: the reference root [paths.PREV_RUNS, 0.12's runs]; a run of --redo
replaces the reference run of its arm (default paths.PREV_REDO, the 20 runs 0.12 redid, only when --ref is not given).
Writes OUT_PREFIX.changes.tsv (if a call changed) and OUT_PREFIX.chrom.tsv (the release's scored chromosomes).
Port of v12 bench/score_dev.py (its TAG is RUNS here; it wrote results/TAG.*)."""
import argparse
from pathlib import Path

import pandas as pd

from benchmark import paths
from benchmark import score as S


def table(f):
    return pd.read_csv(f, sep="\t", dtype=str, keep_default_na=False).set_index("chrom")


def changed(col):
    return col.str.split(">").map(lambda z: z[0] != z[1])


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("runs", type=Path, help="release root: RUNS/<arm>/<id>")
    ap.add_argument("jobs", help="job table (arm, id, ...)")
    ap.add_argument("prefix", help="output prefix: PREFIX.changes.tsv, PREFIX.chrom.tsv")
    ap.add_argument("--ref", type=Path, help="reference root [paths.PREV_RUNS]")
    ap.add_argument("--redo", type=Path, help="root of reference runs redone [paths.PREV_REDO without --ref]")
    ap.add_argument("--manifest", type=Path, default=paths.MANIFEST, help="truth manifest [benchmark/truth.tsv]")
    a = ap.parse_args(argv)
    ref, redo = a.ref or paths.PREV_RUNS, a.redo if a.ref else (a.redo or paths.PREV_REDO)
    for k, v in (("display.width", 250), ("display.max_rows", 2000), ("display.max_columns", 60),
                 ("display.max_colwidth", 60)):
        pd.set_option(k, v)
    T, warned = S.load_manifest(a.manifest), set()
    rows, chg, allrows = [], [], []
    for p in S.jobs_of(a.jobs):
        arm, gid = p[0], p[1]
        new = a.runs / arm / gid
        if not (new / "DONE").exists():
            continue
        if ((new / "DONE").read_text().split() or ["?"])[0] != "0":
            print("FAILED", arm, gid)
            continue
        for cen in sorted(new.glob("*.centromeres.tsv")):
            g = cen.name[: -len(".centromeres.tsv")]
            o = ref / arm / gid
            if redo and (redo / arm / gid / "DONE").exists():
                o = redo / arm / gid
            x = table(cen)
            y = table(o / cen.name) if (o / cen.name).exists() else None
            scored = g in T and S.usable(T[g], warned)
            if scored:
                d, h = S.score_genome(arm, g, new, cen, T[g])
                s = S.summarize(d, arm, g, T[g].kind, h)
                d["arm"] = arm                   # as v12 (gives an empty table the column)
                allrows.append(d)
            if y is None:
                continue
            if scored:
                od, oh = S.score_genome(arm, g, o, o / cen.name, T[g])
                os_ = S.summarize(od, arm, g, T[g].kind, oh)
                rows.append(dict(arm=arm, genome=g, n=s["chromosomes"], right=f"{os_['right']}>{s['right']}",
                                 high=f"{os_['right_high']}>{s['right_high']}",
                                 med=f"{os_['right_medium']}>{s['right_medium']}",
                                 low=f"{os_['right_low']}>{s['right_low']}",
                                 conf_wrong=f"{os_['conf_wrong']}>{s['conf_wrong']}",
                                 wrong_low=f"{os_['wrong_low']}>{s['wrong_low']}", none=f"{os_['none']}>{s['none']}",
                                 dR=s["right"] - os_["right"], dW=s["conf_wrong"] - os_["conf_wrong"]))
                out = {r.chrom: r for r in d.itertuples()}
                oo = {r.chrom: r for r in od.itertuples()}
            for c in x.index:
                if c not in y.index:
                    continue
                b0, b1 = y.loc[c], x.loc[c]
                if (b0.start, b0.end, b0.confidence) != (b1.start, b1.end, b1.confidence):
                    chg.append(dict(arm=arm, g=g, chrom=c, old=f"{b0.start}-{b0.end} {b0.confidence}",
                                    new=f"{b1.start}-{b1.end} {b1.confidence}",
                                    old_outcome=oo[c].outcome if scored and c in oo else "-",
                                    new_outcome=out[c].outcome if scored and c in out else "-", flags=b1["flags"]))
    R = pd.DataFrame(rows)
    if len(R):
        ch = R[(R.dR != 0) | (R.dW != 0) | changed(R.right) | changed(R.high) | changed(R.med) | changed(R.low)
               | changed(R.none)]
        print(ch.drop(columns=["dR", "dW"]).to_string(index=False))
        print(f"arms scored {len(R)}; TOTAL dRight {R.dR.sum()} dConfWrong {R.dW.sum()}")
    C = pd.DataFrame(chg)
    if len(C):
        C.to_csv(f"{a.prefix}.changes.tsv", sep="\t", index=False)
        print(f"\n{len(C)} calls changed (position or confidence):")
        print(C.to_string(index=False))
    if allrows:
        pd.concat(allrows).to_csv(f"{a.prefix}.chrom.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
