#!/usr/bin/env python3
"""Summarize the frozen benchmark table (benchmark/results.tsv): one row per scored chromosome-run.

  python -m benchmark.summary [--table results.tsv] [--canonical canonical.tsv]

Pooled counts take each assembly once per input set: real inputs, genomes with truth (the 5 mammals scored against
another caller's satellite locus, truth_kind 'agreement', left out), one run per genome and input set (maize T before
Tt, HG002's ONT before its PacBio methylation: the first arm), and one genome id per assembly (SAME_ASSEMBLY). Prints
(a) per input set and centromere class: genomes, chromosome-runs, right, wrong, no call, confidently wrong;
(b) confidence calibration over each genome's all-inputs run (its run with the most input types; ties: the table's
all_inputs flag, then the first arm): the README's genomes (those in canonical.tsv) and every genome with truth; (c) the
controls: randomized Hi-C (calls changed against the same genome's run without Hi-C) and holocentric genomes (any
confident call is wrong)."""
import argparse
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
RIGHT = ["right_high", "right_medium", "right_low"]
WRONG = ["conf_wrong", "wrong_low"]
CALL = ["call_start", "call_end", "confidence"]
# One assembly benchmarked under two genome ids: C. neoformans H99 (one FASTA, one truth BED) with an anaphase (cneoA)
# and a G1 (cneoG1) Hi-C library; their runs without Hi-C are identical. Pooled counts keep the library with more read
# pairs, a choice blind to outcomes (hic/<id>.qc.tsv of the 0.16.0 SH runs: cneoA 20,293,703 pairs kept of 20,366,561
# read, cneoG1 20,033,192 of 20,076,461). Key: the id left out of pooled counts; value: the id kept for its assembly.
SAME_ASSEMBLY = {"cneoG1": "cneoA"}


def read_table(path):
    """results.tsv as strings; the leading '#' lines are its header notes."""
    with open(path) as fh:
        skip = 0
        for line in fh:
            if not line.startswith("#"):
                break
            skip += 1
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False, skiprows=skip)


def counts(x, keys):
    """genomes, chromosome-runs, right, wrong, no call and confidently wrong per group of keys."""
    o = x.outcome
    x = x.assign(right=o.isin(RIGHT), wrong=o.isin(WRONG), no_call=o == "none", conf_wrong=o == "conf_wrong")
    return x.groupby(keys).agg(genomes=("genome", "nunique"), chromosomes=("chrom", "size"), right=("right", "sum"),
                               wrong=("wrong", "sum"), no_call=("no_call", "sum"), conf_wrong=("conf_wrong", "sum"))


def one_run(x):
    """One run per genome and input set: the first arm."""
    return x[x.arm == x.groupby(["genome", "inputs"]).arm.transform("min")]


def pooled_rows(d):
    """The rows pooled counts use: real inputs, genomes with truth ('agreement' left out), one run per genome and input
    set, one genome id per assembly."""
    x = d[(d.control == "none") & (d.truth_kind != "agreement") & ~d.genome.isin(list(SAME_ASSEMBLY))]
    return one_run(x)


def all_inputs(x):
    """The rows of each genome's all-inputs run in x: its run with the most input types; ties: the table's all_inputs
    flag, then the first arm."""
    r = x.drop_duplicates(["genome", "arm"])
    r = r.assign(k=r.inputs.str.count(r"\+"), f=r.all_inputs == "1")
    r = r.sort_values(["genome", "k", "f", "arm"], ascending=[True, False, False, True]).drop_duplicates("genome")
    return x[pd.MultiIndex.from_frame(x[["genome", "arm"]]).isin(pd.MultiIndex.from_frame(r[["genome", "arm"]]))]


def by_inputs(d):
    """(a) per input set and centromere class, then all classes, inputs ordered by how many they hold."""
    x = pooled_rows(d)
    t = counts(x, ["inputs", "cen_class"])
    a = counts(x, ["inputs"])
    a.index = pd.MultiIndex.from_tuples([(i, "all") for i in a.index], names=t.index.names)
    t = pd.concat([t, a])
    n = {i: (i.count("+"), i) for i in t.index.get_level_values(0)}
    return t.loc[sorted(t.index, key=lambda k: (n[k[0]], k[1] == "all", k[1]))].astype(int)


def calibration(x):
    """(b) calls, right, wrong and the share right per confidence level, and the chromosome-runs with no call."""
    rows = []
    for c in ("high", "medium", "low"):
        y = x[x.confidence == c]
        r = int(y.outcome.isin(RIGHT).sum())
        rows.append(dict(confidence=c, calls=len(y), right=r, wrong=len(y) - r,
                         share_right=f"{100 * r / len(y):.1f}%" if len(y) else "."))
    rows.append(dict(confidence="no call", calls=int((x.outcome == "none").sum()), right="", wrong="",
                     share_right=""))
    return pd.DataFrame(rows).set_index("confidence")


def controls(d):
    """(c) randomized Hi-C: each run against the same genome's run with the same inputs but no Hi-C; holocentric."""
    rows = []
    d = d[~d.genome.isin(list(SAME_ASSEMBLY))]                                  # one Hi-C library per assembly
    real = d[~d.control.str.contains("random_hic")]
    for (arm, g), x in d[d.control.str.contains("random_hic")].groupby(["arm", "genome"], sort=True):
        want = x.inputs.iloc[0].replace("+Hi-C(randomized)", "")
        p = real[(real.genome == g) & (real.inputs == want)]
        arms = sorted(p.arm.unique())                                           # several: the first, as in (a)
        changed = "."
        if arms:
            m = x.merge(p[p.arm == arms[0]], on="chrom", suffixes=("", "_p"))
            changed = int((m[CALL].values != m[[c + "_p" for c in CALL]].values).any(axis=1).sum())
        rows.append(dict(control=x.control.iloc[0], arm=arm, genome=g, chromosomes=len(x),
                         confident=int(x.confidence.isin(["high", "medium"]).sum()),
                         partner=arms[0] if arms else "none", changed=changed))
    r = pd.DataFrame(rows)
    h = d[d.control == "holocentric"]
    holo = h.assign(confident=h.confidence.isin(["high", "medium"])).groupby("inputs").agg(
        genomes=("genome", "nunique"), chromosomes=("chrom", "size"), confident=("confident", "sum"))
    return r, holo


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--table", type=Path, default=HERE / "results.tsv",
                    help="the results table [benchmark/results.tsv]")
    ap.add_argument("--canonical", type=Path, default=HERE / "canonical.tsv",
                    help="genome labels; its genomes are the README's calibration set [benchmark/canonical.tsv]")
    a = ap.parse_args(argv)
    d = read_table(a.table)
    canon = set(pd.read_csv(a.canonical, sep="\t", dtype=str).genome)
    w = dict(line_width=200)
    print(f"{a.table.name}: {len(d)} chromosome-runs, {d.groupby(['arm', 'genome']).ngroups} runs, "
          f"{d.genome.nunique()} genomes\n")
    print("(a) per input set and centromere class (real inputs; genomes with truth, 'agreement' left out)")
    print(by_inputs(d).to_string(**w), "\n")
    x = all_inputs(pooled_rows(d))
    for name, y in (("the README's genomes (canonical.tsv)", x[x.genome.isin(canon)]), ("every genome", x)):
        right = int(y.outcome.isin(RIGHT).sum())
        print(f"(b) calibration, all inputs, {name}: {y.genome.nunique()} genomes, {len(y)} chromosomes with truth, "
              f"{right} right ({100 * right / max(len(y), 1):.1f}%), "
              f"{int((y.outcome == 'conf_wrong').sum())} confidently wrong")
        print(calibration(y).to_string(**w), "\n")
    r, holo = controls(d)
    if len(r):
        print(f"(c) randomized Hi-C: {len(r)} runs, {r.genome.nunique()} genomes, {int(r.chromosomes.sum())} "
              f"chromosome-runs; calls changed against the run without Hi-C: "
              f"{int(pd.to_numeric(r.changed, errors='coerce').fillna(0).sum())} "
              f"({int((r.partner == 'none').sum())} runs without a partner); confident calls "
              f"{int(r.confident.sum())}")
        print(r.groupby(["control", "arm"]).agg(runs=("genome", "size"), chromosomes=("chromosomes", "sum"),
                                                confident=("confident", "sum")).to_string(**w), "\n")
    if len(holo):
        print("(c) holocentric genomes (any confident call is wrong)")
        print(holo.astype(int).to_string(**w))


if __name__ == "__main__":
    main()
