"""benchmark/md_tables.py and md_013.py: the calibration table, the telomere table, bin sizes."""
import pandas as pd

from benchmark import md_013, md_tables


def test_calibration_counts_each_level_and_the_chromosomes_without_a_call():
    tr = pd.DataFrame(dict(genome=["a", "a", "a", "b", "b"], confidence=["high", "high", "medium", "low", "none"],
                           ok=["R", "W", "R", "W", "none"],
                           outcome=["right_high", "conf_wrong", "right_medium", "wrong_low", "none"]))
    t = md_tables.calibration(tr)
    assert t.startswith("2 genomes, each with all its inputs: 5 chromosomes with truth, 2 right (40.0%), "
                        "1 confidently wrong.")
    assert "| `high` | 2 | 1 | 1 | 50.0% |" in t
    assert "| `medium` | 1 | 1 | 0 | 100.0% |" in t
    assert "| `low` | 1 | 0 | 1 | 0.0% |" in t
    assert "| no call | 1 | | | |" in t


def test_telomere_rows_one_per_assembly_and_complete_ones_in_a_sentence(tmp_path):
    for g, unit, ends in (("cneoA", "AACCCCCT", "13"), ("cneoG1", "AACCCCCT", "13"), ("pot", ".", "0"),
                          ("t2t", "AACCCT", "28")):
        d = tmp_path / "S" / g
        d.mkdir(parents=True)
        pd.DataFrame([dict(genome=g, chromosomes="14", telomere_unit=unit, telomere_ends=ends, t2t_chromosomes="3/14",
                           telomere_interstitial="0", rdna_arrays="2", rdna_chromosomes="c1,c2")]
                     ).to_csv(d / "cenfuse_summary.tsv", sep="\t", index=False)
    labels = {"cneoA": "*C. neoformans* (anaphase Hi-C)", "cneoG1": "*C. neoformans* (G1 Hi-C)", "pot": "potato",
              "t2t": "goat"}
    tab, full = md_tables.telomeres(tmp_path, ["cneoA", "cneoG1", "pot", "no_run", "t2t"], labels.get)
    t = tab.splitlines()
    assert len(t) == 4                                   # header, rule, one row per incomplete assembly with a run
    assert t[2] == "| *C. neoformans* | 14 | AACCCCCT | 13/28 | 3/14 | 0 | 2 (2) |"
    assert t[3].startswith("| potato | 14 | none found | 0/28 |")
    assert full == "An array at every chromosome end in 1 assemblies: goat (3/14 complete)."


def test_bin_sizes_read_as_kb_and_mb():
    assert md_013.bp("250000") == "250 kb" and md_013.bp(2_000_000) == "2 Mb" and md_013.bp(".") == "."


def test_rule_calibration_splits_each_rule_by_confidence(tmp_path):
    d = tmp_path / "S" / "g"
    d.mkdir(parents=True)
    pd.DataFrame(dict(chrom=["c1", "c2", "c3", "c4"], rule=["sequence", "sequence", "hic_denovo_only", "at_only"],
                      seq_evidence=["satellite", "satellite", "none", "none"])).to_csv(d / "g.centromeres.tsv", sep="\t",
                                                                                       index=False)
    tr = pd.DataFrame(dict(genome=["g"] * 4, arm=["S"] * 4, chrom=["c1", "c2", "c3", "c4"],
                           confidence=["low", "low", "low", "none"], ok=["R", "W", "R", "none"]))
    t = md_tables.rule_calibration(tr, tmp_path).splitlines()
    assert t[0].startswith("| rule that set the call |") and len(t) == 4          # the chromosome without a call is left out
    assert "| the sequence call as it is (satellite) (`sequence`) | `low` | 2 | 1 | 50.0% |" in t
    assert "| Hi-C's prior-free peak alone (`hic_denovo_only`) | `low` | 1 | 1 | 100.0% |" in t
    assert {"rep_only", "rep_pair", "fam_only", "cis_end_pair", "fold_offset"} <= set(md_tables.RULE_TEXT)
