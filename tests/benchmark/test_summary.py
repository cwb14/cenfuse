"""benchmark/summary.py: the frozen results table, per input set, calibration and controls (a synthetic table)."""
from benchmark import summary

COLS = ["genome", "arm", "version", "inputs", "control", "cen_class", "truth_kind", "all_inputs", "chrom", "chrom_len",
        "call_start", "call_end", "confidence", "relation", "rule", "outcome", "truth_start", "truth_end",
        "score_rule", "tol_bp", "fam_gate", "fam_start", "fam_end", "kt_start", "kt_end", "kt_on_truth", "truth_all"]
ROWS = [  # genome arm inputs control cen_class truth_kind all_inputs call confidence outcome
    ("g1", "S", "FASTA", "none", "satellite", "CENP-A", "0", "101", "200", "high", "right_high"),
    ("g1", "SH", "FASTA+Hi-C", "none", "satellite", "CENP-A", "1", "501", "600", "medium", "conf_wrong"),
    ("g1", "NEG2", "FASTA+Hi-C(randomized)", "random_hic", "satellite", "CENP-A", "0", "101", "300", "high",
     "right_high"),
    ("g2", "S", "FASTA", "holocentric", "holocentric", "holocentric", "0", "5", "50", "low", "wrong_low"),
    ("g3", "S", "FASTA", "none", "satellite", "agreement", "1", "101", "200", "high", "right_high"),
    ("g1", "Sb", "FASTA", "none", "satellite", "CENP-A", "0", "101", "200", "high", "right_high"),   # a second run
]


def write(tmp_path):
    t = tmp_path / "results.tsv"
    lines = ["# header note one", "# header note two", "\t".join(COLS)]
    for g, arm, inp, ctl, cls, kind, alli, a, z, conf, out in ROWS:
        r = dict(zip(COLS, ["."] * len(COLS)))
        r.update(genome=g, arm=arm, version="0.0.0", inputs=inp, control=ctl, cen_class=cls, truth_kind=kind,
                 all_inputs=alli, chrom="chr1", chrom_len="1000", call_start=a, call_end=z, confidence=conf,
                 outcome=out, truth_start="" if g == "g1" else ".", truth_end="" if g == "g1" else ".")
        lines.append("\t".join(r[c] for c in COLS))
    t.write_text("\n".join(lines) + "\n")
    c = tmp_path / "canonical.tsv"
    c.write_text("genome\tlabel\ng1\tgenome one\n")
    return t, c


def test_counts_calibration_and_controls(tmp_path, capsys):
    t, c = write(tmp_path)
    d = summary.read_table(t)
    assert len(d) == 6 and d.truth_start.iloc[0] == ""          # header notes skipped, withheld truth kept empty
    a = summary.by_inputs(d)
    assert a.loc[("FASTA", "all")].tolist() == [1, 1, 1, 0, 0, 0]     # holocentric, agreement, second run out
    assert a.loc[("FASTA+Hi-C", "satellite")].tolist() == [1, 1, 0, 1, 0, 1]
    assert ("FASTA+Hi-C(randomized)", "all") not in a.index
    cal = summary.calibration(d[d.all_inputs == "1"])
    assert cal.loc["medium", "calls"] == 1 and cal.loc["high", "right"] == 1
    r, holo = summary.controls(d)
    assert r.partner.tolist() == ["S"] and r.changed.tolist() == [1]   # the randomized run moved its call
    assert holo.loc["FASTA"].tolist() == [1, 1, 0]
    summary.main(["--table", str(t), "--canonical", str(c)])
    out = capsys.readouterr().out
    assert "the README's genomes (canonical.tsv): 1 genomes, 1 chromosomes with truth, 0 right" in out
    assert "every genome: 1 genomes, 1 chromosomes" in out            # g3 is 'agreement': not truth


def test_one_assembly_once_and_the_all_inputs_run():
    """Pooled rows keep one genome id per assembly; a genome's all-inputs run is its run with the most input types, ties
    broken by the table's all_inputs flag."""
    import pandas as pd
    rows = [("cneoA", "S", "FASTA", "0"), ("cneoG1", "S", "FASTA", "0"),        # one assembly under two ids
            ("g", "S", "FASTA", "1"), ("g", "HT", "FASTA+tes", "0"),            # the flag on fewer inputs
            ("h", "M", "FASTA+tes+meth(CpG,CHG,CHH)", "1"), ("h", "THn", "FASTA+tes+Hi-C", "0")]   # a tie: the flag
    d = pd.DataFrame(rows, columns=["genome", "arm", "inputs", "all_inputs"]).assign(
        control="none", truth_kind="CENH3", chrom="chr1", outcome="right_low")
    assert set(summary.pooled_rows(d).genome) == {"cneoA", "g", "h"} and summary.SAME_ASSEMBLY["cneoG1"] == "cneoA"
    a = summary.all_inputs(summary.pooled_rows(d))
    assert dict(zip(a.genome, a.arm)) == {"cneoA": "S", "g": "HT", "h": "M"}
