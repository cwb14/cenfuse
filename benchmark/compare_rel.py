#!/usr/bin/env python3
"""Each release run against its reference run (the run the docs' numbers came from), per chromosome.

  python -m benchmark.compare_rel JOBS REL_ROOT OUT_PREFIX
JOBS: a job table whose column 6 is the reference run (paths.JOBS_CMP: 0.12's runs). For every <g>.centromeres.tsv of
a finished run REL_ROOT/<arm>/<id>, compares the call columns that every version since 0.8 writes. Writes
OUT_PREFIX.diffs.tsv (one row per differing chromosome and column, old and new values) and OUT_PREFIX.summary.tsv (per
run: chromosomes, changed chromosomes per column group), and prints the runs with a difference.
Port of v12 bench/rel/compare_rel.py."""
import argparse
from pathlib import Path

import pandas as pd

from benchmark import paths
from benchmark.score import jobs_of

CALL = ["start", "end", "confidence", "relation", "call_source"]
EVID = ["support", "against", "flags"]
SEQ = ["seq_start", "seq_end", "seq_confidence", "seq_flags", "seq_alt_start", "seq_alt_end", "main_family"]
CPG = ["cpg_kt_start", "cpg_kt_end"]
MORPH = ["morphology", "telomeres", "gaps", "call_gaps", "t2t"]
GROUPS = {"call": CALL, "evidence": EVID, "seq": SEQ, "cpg": CPG, "morph": MORPH}


def table(f):
    return pd.read_csv(f, sep="\t", dtype=str, keep_default_na=False).set_index("chrom")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("jobs", help="job table; column 6 = the reference run")
    ap.add_argument("root", type=Path, help="release root: ROOT/<arm>/<id>")
    ap.add_argument("prefix", help="output prefix: PREFIX.diffs.tsv, PREFIX.summary.tsv")
    a = ap.parse_args(argv)
    diffs, summ = [], []
    for p in jobs_of(a.jobs):
        arm, gid, ref = p[0], p[1], paths.local(p[5])
        new = a.root / arm / gid
        if not (new / "DONE").exists() or ((new / "DONE").read_text().split() or ["?"])[0] != "0":
            continue
        for f in sorted(new.glob("*.centromeres.tsv")):
            g = f.name[: -len(".centromeres.tsv")]
            fo = ref / f.name
            if not fo.exists():
                summ.append(dict(arm=arm, id=gid, genome=g, note="no reference table"))
                continue
            old, cur = table(fo), table(f)
            row = dict(arm=arm, id=gid, genome=g, chromosomes=len(cur),
                       missing=len(old.index.symmetric_difference(cur.index)))
            for grp, cols in GROUPS.items():
                changed = set()
                for c in cols:
                    if c not in old.columns or c not in cur.columns:
                        continue
                    for ch in old.index.intersection(cur.index):
                        if old.at[ch, c] != cur.at[ch, c]:
                            changed.add(ch)
                            diffs.append(dict(arm=arm, id=gid, genome=g, chrom=ch, column=c, old=old.at[ch, c],
                                              new=cur.at[ch, c]))
                row[grp] = len(changed)
            summ.append(row)
    s = pd.DataFrame(summ)
    s.to_csv(f"{a.prefix}.summary.tsv", sep="\t", index=False)
    pd.DataFrame(diffs).to_csv(f"{a.prefix}.diffs.tsv", sep="\t", index=False)
    cols = [c for c in GROUPS if c in s]
    ch = s[(s[cols].fillna(0).sum(axis=1) > 0) | s.get("note", pd.Series(dtype=str)).notna()] if len(s) else s
    print(f"{len(s)} genome runs compared; {len(ch)} with differences")
    with pd.option_context("display.width", 200, "display.max_rows", 300):
        print(ch.to_string(index=False))


if __name__ == "__main__":
    main()
