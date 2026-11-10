#!/usr/bin/env python3
"""Tables of docs/benchmark.md from one release's scores and runs; nothing is typed in.

  python -m benchmark.md_013 OUTDIR [--scores PREFIX] [--runs ROOT] [--jobs JOBS] [--cis-survey TSV]
Writes OUTDIR/<NAME>.md (one locus per chromosome, LTR signature, AT troughs, the cis signal, the cis signal tested by
trans contacts, which end; the long ones as a table of the informative rows and a sentence listing the others) and
OUTDIR/one_per_chromosome.tsv, OUTDIR/cis_trans.tsv. Reads the release's runs table, the
runs' summaries, each FASTA-only arm's <g>_satellite_test.json (satellite dir = job column 5) and the cis survey. The
tables show this release only (--prev and --versions are accepted and ignored, so older chains still run).
Port of v12 bench/md_013.py (the module keeps the name of the release that introduced it)."""
import argparse
import json
from pathlib import Path

import pandas as pd

from benchmark import paths
from benchmark.md_tables import canonical, label

OTHER_LABEL = {"csat9930": "cucumber 9930", "vvPN": "grape PN40024"}    # genomes outside canonical.tsv
REAL = ("SH", "HSH", "FULL")     # CPG repeats FULL's Hi-C, HTH the SH library with a tes/ that gave no signature; NEG,
                                 # NEG2, TNEG are the randomized libraries
INP = {"SH": "FASTA, Hi-C", "HSH": "FASTA, Hi-C", "HTH": "FASTA, `tes/`, Hi-C", "FULL": "FASTA, genes, Hi-C"}


def cell(r):
    if r is None:
        return ""
    if r.kind == "holocentric":
        return f"{r.conf_wrong} confident / {r.chromosomes}"
    return f"{r.right}/{r.chromosomes}" + (f" ({r.conf_wrong})" if r.conf_wrong else "")


def bp(x):
    """250000 -> '250 kb', 2000000 -> '2 Mb'; anything else unchanged."""
    try:
        v = int(float(x))
    except (TypeError, ValueError):
        return x
    return f"{v / 1e6:g} Mb" if v >= 1_000_000 else f"{v / 1e3:g} kb" if v >= 1000 else f"{v} bp"


def runs_table(prefix):
    t = pd.read_csv(f"{prefix}.runs.tsv", sep="\t", dtype={"hic_gate_pass": str})
    return {(r.genome, r.arm): r for r in t.itertuples(index=False)}


def write_tables(out, scores, rel, jobs_file, can, cis_survey):
    out.mkdir(parents=True, exist_ok=True)
    N = runs_table(scores)
    jobs = {tuple(x.split("\t")[:2]): x.rstrip("\n").split("\t") for x in open(jobs_file)}

    def lab(g):
        return label(g, can, OTHER_LABEL)

    def w(name, lines):
        (out / f"{name}.md").write_text("\n".join(lines).rstrip("\n") + "\n")

    def summary(arm, g):
        f = rel / arm / g / "cenfuse_summary.tsv"
        if not f.exists():
            return None
        return pd.read_csv(f, sep="\t", dtype=str, keep_default_na=False).iloc[0]

    # ---- one locus per chromosome: the FASTA-only arm of every genome with a centromeric family; the genomes that
    # fail as a table, those that pass as a sentence (no family called centromeric: the test does not run)
    L = ["| genome | chromosomes | covered by the N largest clusters | by chance | P | FASTA only |",
         "|---|---|---|---|---|---|"]
    rows, passed = [], []
    for g in can.index:
        if (g, "S") not in N:
            continue
        f = paths.local(jobs.get(("S", g))[4]) / f"{g}_satellite_test.json"
        if not f.exists():
            continue
        t = json.loads(f.read_text())
        if t["status"] == "pass":
            passed.append((lab(g), t))
        elif t["status"] != "none":
            L.append(f"| {lab(g)} | {t['n']} | {t['covered']} | {t['null_median']:.0f} | {t['p']:.2g} "
                     f"| {cell(N.get((g, 'S')))} |")
        rows.append(dict(genome=g, **{k: t[k] for k in ("n", "covered", "null_median", "p", "status")},
                         conf_wrong=getattr(N.get((g, "S")), "conf_wrong", None)))
    w("TABLE_ONE_PER_CHROMOSOME", L)
    ps = [t["p"] for _, t in passed]
    w("ONE_PER_CHROMOSOME_PASS", [f"It passes in {len(passed)} genomes (P {min(ps):.2g}-{max(ps):.2g}; chromosomes "
                                  "covered by the N largest clusters / N): " + ", ".join(
                                      f"{x} {t['covered']}/{t['n']}" for x, t in passed) + "." if passed else
                                  "It passes in no genome."])
    pd.DataFrame(rows).to_csv(out / "one_per_chromosome.tsv", sep="\t", index=False)

    # ---- the LTR signature: genomes with tes/
    L = ["| genome | signature | fold in its cores | decides? | FASTA only | + `tes/` | + `tes/` + Hi-C |",
         "|---|---|---|---|---|---|---|"]
    for g in can.index:
        if (g, "T") not in N:
            continue
        s = summary("T", g)
        hic = N.get((g, "TH")) or N.get((g, "THn"))
        if s.get("te_decisive") == "True":
            dec = "yes"
        elif s.get("te_marker", ".") != ".":
            dec = "no: " + s.get("te_verdict", "").rpartition("; ")[2].replace("flags only", "").strip(" ()")
        else:
            dec = "none"
        L.append(f"| {lab(g)} | {s.get('te_marker', '.')} | {s.get('te_specificity', '.')} | {dec} | "
                 f"{cell(N.get((g, 'S')))} | {cell(N.get((g, 'T')))} | {cell(hic)} |")
    w("TABLE_SIGNATURE", L)

    # ---- AT troughs: FASTA-only arm and the arm with Hi-C
    L = ["| genome | window | P | test | FASTA only | + Hi-C | troughs vs Hi-C |", "|---|---|---|---|---|---|---|"]
    for g in can.index:
        s = summary("S", g)
        if s is None or s.get("at_gate", ".") in (".", ""):
            continue
        h = summary("SH", g)
        L.append(f"| {lab(g)} | {bp(s.get('at_window', '.'))} | {s.get('at_p', '.')} | {s.get('at_gate')} | "
                 f"{cell(N.get((g, 'S')))} | {cell(N.get((g, 'SH')))} | {'' if h is None else h.get('at_hic', '.')} |")
    w("TABLE_AT", L)

    # ---- the cis signal: every genome with Hi-C and truth (the survey), and the test against the LTR cores; the
    # genomes with neither extreme on truth beyond chance (P > 0.05) and no LTR-core test run, as a sentence
    L = ["| genome | chromosomes | peak on truth | chance | P | dip on truth | chance | P | against the LTR cores |",
         "|---|---|---|---|---|---|---|---|---|"]
    flat = []
    for r in pd.read_csv(cis_survey, sep="\t").itertuples(index=False):
        arm = "TH" if (r.id, "TH") in N else "THn" if (r.id, "THn") in N else None
        s = summary(arm, r.id) if arm else None
        x = s.get("cis_line", ".") if s is not None else ""
        if min(r.p_peak, r.p_dip) > 0.05 and x in (".", ""):
            flat.append(lab(r.id))
            continue
        L.append(f"| {lab(r.id)} | {r.n} | {r.peak} | {r.chance_peak:.2f} | {r.p_peak:.2g} | {r.dip} | {r.chance_dip:.2f} | "
                 f"{r.p_dip:.2g} | {x if x != '.' else ('not run: Hi-C test passed' if s is not None else '')} |")
    w("TABLE_CIS", L)
    w("CIS_FLAT", [f"Neither extreme lies on truth beyond chance (P > 0.05) in the other {len(flat)} genomes: "
                   + ", ".join(flat) + "."])

    # ---- the cis signal tested by trans contacts: real libraries as rows, randomized ones counted
    L = ["| genome (inputs) | Hi-C bin | trans peaks at the cis extreme | verdict | acts on | right |",
         "|---|---|---|---|---|---|"]
    rows, rnd, failed = [], [], []
    for (g, arm), r in sorted(N.items(), key=lambda kv: (kv[0][1] != "SH", kv[0][1], lab(kv[0][0]))):
        s = summary(arm, g)
        if s is None or "trans contacts" not in str(s.get("cis_line", "")):
            continue
        line = str(s.get("cis_line"))
        head, _, tail = line.partition(" chromosomes (P ")
        what = head.replace("trans peaks at the cis ", "")            # 'peak on 10/12'
        p_, _, rest = tail.partition("): ")
        verdict = rest.partition(",")[0]
        acts = line.rpartition("acts on ")[2]
        rows.append(dict(genome=g, arm=arm, cis_line=line, verdict=verdict, acts=acts, right=cell(r)))
        if arm in REAL and verdict == "fail":
            failed.append(f"{lab(g)} ({INP[arm]})")
        elif arm in REAL:
            L.append(f"| {lab(g)} ({INP[arm]}) | {bp(s.get('hic_bin', '.'))} | {what.replace(' on ', ': ')} (P {p_}) | "
                     f"{verdict} | {acts} | {cell(r)} |")
        elif arm in ("NEG", "NEG2", "TNEG"):
            rnd.append((lab(g), arm, verdict, what))
    w("TABLE_CIS_TRANS", L)
    w("CIS_TRANS_FAIL", [f"It fails, and changes nothing, in {len(failed)} more runs: " + ", ".join(failed) + "."])
    sup = [f"{x[0]} ({x[3]})" for x in rnd if x[2] == "support"]
    npass = sum(x[2] == "pass" for x in rnd)              # v12 typed "never passes"; counted here
    w("CIS_TRANS_RANDOM", [f"On {len(rnd)} randomized libraries (partners re-paired, coverage kept) the test "
                           f"{'never passes' if not npass else f'passes on {npass}'}; it reaches support on "
                           f"{len(sup)} ({'; '.join(sup) if sup else 'none'}), where it does not act."])
    pd.DataFrame(rows).to_csv(out / "cis_trans.tsv", sep="\t", index=False)

    # ---- which end: calls in a terminal eighth, per genome with Hi-C
    L = ["| genome | calls at a chromosome end | at the higher end | test | raised by it |", "|---|---|---|---|---|"]
    for g in can.index:
        arm = next((a for a in ("FULL", "SH") if (g, a) in N), None)
        if arm is None:
            continue
        f = rel / arm / g / f"{g}.centromeres.tsv"
        if not f.exists():
            continue
        t = pd.read_csv(f, sep="\t", dtype=str, keep_default_na=False)
        if "hic_end_k" not in t or not (t.hic_end_n.ne(".") & t.hic_end_n.ne("")).any():
            continue
        x = t[t.hic_end_n.ne(".") & t.hic_end_n.ne("")]
        k, n = pd.to_numeric(x.hic_end_k, errors="coerce").max(), pd.to_numeric(x.hic_end_n, errors="coerce").max()
        if not n or n != n or n < 3:
            continue
        raised = int(t["flags"].str.contains("hic_end", regex=False).sum())
        L.append(f"| {lab(g)} ({'genes, ' if arm == 'FULL' else ''}Hi-C) | {int(n) + 1} | {int(k)} of {int(n)} others | "
                 f"{'pass' if raised else 'fail or not needed'} | {raised} |")
    w("TABLE_END", L)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("out", type=Path, help="output folder")
    ap.add_argument("--scores", default=paths.SCORES, help="this release's score.py prefix [paths.SCORES]")
    ap.add_argument("--runs", type=Path, default=paths.RUNS, help="release root of the scored runs [paths.RUNS]")
    ap.add_argument("--jobs", default=paths.JOBS, help="job table (column 5: satellite dir) [paths.JOBS]")
    ap.add_argument("--canonical", type=Path, default=paths.CANONICAL, help="genome table [benchmark/canonical.tsv]")
    ap.add_argument("--cis-survey", default=paths.CIS_SURVEY, help="cis-signal survey table [paths.CIS_SURVEY]")
    ap.add_argument("--prev", help="ignored: the tables show one release (accepted so older chains still run)")
    ap.add_argument("--versions", nargs=2, metavar=("OLD", "NEW"), help="ignored, as --prev")
    a = ap.parse_args(argv)
    write_tables(a.out, a.scores, a.runs, a.jobs, canonical(a.canonical), a.cis_survey)
    print("written to", a.out)


if __name__ == "__main__":
    main()
