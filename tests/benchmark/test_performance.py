"""benchmark/performance.py: the performance figure from a synthetic results table."""
import re

import numpy as np
import pandas as pd

from benchmark import performance, summary

COLS = ["genome", "arm", "version", "inputs", "control", "cen_class", "truth_kind", "all_inputs", "chrom", "chrom_len",
        "call_start", "call_end", "confidence", "relation", "rule", "outcome", "truth_start", "truth_end",
        "score_rule", "tol_bp", "fam_gate", "fam_start", "fam_end", "kt_start", "kt_end", "kt_on_truth"]
F, T, H, TH = "FASTA", "FASTA+tes", "FASTA+Hi-C", "FASTA+tes+Hi-C"
RIGHT = ("101", "200", "high", "right_high")                 # truth 151-250 on every chromosome
LOW = ("101", "200", "low", "right_low")
WRONG = ("701", "800", "low", "wrong_low")
NONE = (".", ".", "none", "none")


def run(g, arm, inputs, cls, calls, alli="0", control="none", kind="CENH3", truth=("151", "250")):
    """One run: a row per chromosome (calls: one (start, end, confidence, outcome) each)."""
    rows = []
    for i, (a, z, conf, out) in enumerate(calls):
        t0, t1 = truth if out != "none" else (".", ".")
        if truth[0] == "":
            t0 = t1 = ""
        rows.append(dict(genome=g, arm=arm, version="0.0.0", inputs=inputs, control=control, cen_class=cls,
                         truth_kind=kind, all_inputs=alli, chrom=f"chr{i + 1}", chrom_len="1000", call_start=a,
                         call_end=z, confidence=conf, outcome=out, truth_start=t0, truth_end=t1, tol_bp="10"))
    return rows


def table(tmp_path):
    rows = []
    for g, gains in (("sat", [RIGHT, RIGHT]), ("mid", [WRONG, LOW]), ("te", [RIGHT, RIGHT])):
        cls = "satellite" if g == "sat" else "TE_or_none"
        fasta = [RIGHT, RIGHT] if g == "sat" else [WRONG, WRONG]          # gains: sat 0, mid 0.5, te 1
        rows += run(g, "S", F, cls, fasta) + run(g, "T", T, cls, gains) + run(g, "SH", H, cls, fasta)
        rows += run(g, "TH", TH, cls, gains, alli="1")
    rows += run("ape", "S", F, "satellite", [RIGHT, NONE], alli="1", truth=("", ""))     # truth withheld
    for g in ("cneoA", "cneoG1"):          # one assembly, two Hi-C libraries: pooled once (summary.SAME_ASSEMBLY)
        rows += run(g, "S", F, "point_regional", [NONE, LOW])
        rows += run(g, "SH", H, "point_regional", [RIGHT, RIGHT], alli="1")
    rows += run("agr", "S", F, "satellite", [RIGHT, WRONG], alli="1", kind="agreement")
    rows += run("holo", "S", F, "holocentric", [WRONG, NONE], control="holocentric", kind="holocentric")
    rows += run("sat", "NEG2", "FASTA+Hi-C(randomized)", "satellite", [RIGHT, RIGHT], control="random_hic")
    t = tmp_path / "results.tsv"
    d = pd.DataFrame(rows).reindex(columns=COLS).fillna(".")
    t.write_text("# a header note\n" + d.to_csv(sep="\t", index=False))
    c = tmp_path / "canonical.tsv"
    c.write_text("genome\tlabel\tlineage\nsat\tBrassica rapa X1\tplants\nape\tSumatran orangutan\tmammals\n")
    return t, c


def test_chance_counts_one_based_lengths():
    d = pd.DataFrame(dict(call_start=["101", ".", "101", "1"], call_end=["200", ".", "200", "1000"],
                          truth_start=["151", ".", "", "1"], truth_end=["250", ".", "", "1000"],
                          tol_bp=["10"] * 4, chrom_len=["1000"] * 4))
    c = performance.chance(d)
    assert np.isclose(c[0], (100 + 100 + 20) / 1000) and c[1] == 0.0 and np.isnan(c[2]) and c[3] == 1.0


def test_example_rule_and_figure(tmp_path):
    t, c = table(tmp_path)
    d = summary.read_table(t)
    d["chance"] = performance.chance(d)
    perf, every, agree, holo, alli = performance.tables(d)
    assert "cneoG1" in set(every.genome) and "cneoG1" not in set(perf.genome) | set(alli.genome)
    g, gain, med = performance.pick_example(perf)
    assert sorted(gain.index) == ["mid", "sat", "te"] and med == 0.5 and g == "mid"      # the median gain
    assert performance.runs_of("Brassica rapa X1", "plants") == [("B. rapa", "italic"), (" X1", "normal")]
    assert performance.runs_of("Sumatran orangutan", "mammals") == [("Sumatran orangutan", "normal")]
    png, pdf, num = tmp_path / "p.png", tmp_path / "p.pdf", tmp_path / "n.tsv"
    performance.main(["--table", str(t), "--canonical", str(c), "--truth", str(tmp_path / "none.tsv"),
                      "--png", str(png), "--pdf", str(pdf), "--numbers", str(num), "--dpi", "60"])
    assert png.stat().st_size > 0
    assert len(re.findall(rb"/Type\s*/Page\b", pdf.read_bytes())) == 2               # the figure, the legend
    n = pd.read_csv(num, sep="\t", dtype=str)
    v = n.set_index(["panel", "item", "group"]).value
    assert v[("A", "example genome", "mid")] == "mid"
    assert v[("B", "chromosomes", "satellite | FASTA")] == "4"            # sat and the ape; 'agr' left out
    assert v[("B", "chance_rows", "satellite | FASTA")] == "2"            # the ape's truth is withheld
    assert v[("E", "calls changed", "randomized")] == "0" and v[("E", "runs", "real")] == "7"   # SH, TH x 3; cneoA SH
    assert v[("B", "chromosomes", "point_regional | FASTA")] == "2" and v[("C", "chromosomes", "cneoG1")] == "2"
    assert v[("D", "genomes", "all inputs")] == "5"                       # sat, mid, te, the ape, cneoA
