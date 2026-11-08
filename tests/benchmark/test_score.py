"""The frozen benchmark scorer (benchmark/score.py) on tiny synthetic runs and truths: every rule, the tolerance with
and without Hi-C, the chance formula, the Poisson-binomial P, missing truth, and the repository's manifest."""
import math

import pandas as pd
import pytest

from benchmark import paths
from benchmark import score as S

CEN_COLS = ["chrom", "chrom_len", "start", "end", "confidence", "relation", "flags", "support"]
CHROM_COLS = ["arm", "genome", "chrom", "chrom_len", "kind", "rule", "tol", "hic_bin", "start", "end", "confidence",
              "relation", "flags", "support", "against", "seq_start", "seq_end", "seq_confidence", "cpg_kt_start",
              "cpg_kt_end", "morphology", "ci", "telomeres", "t2t", "ok", "gap", "truth", "chance", "outcome"]
RUNS_COLS = ["arm", "genome", "kind", "chromosomes", "right", "right_high", "right_medium", "right_low", "conf_wrong",
             "wrong_low", "none", "chance", "P", "hic_bin", "hic_gate_pass", "hic_gate_support", "hic_gate_untested",
             "fold_gate", "align_gate", "hic_axis", "hic_architecture", "cpg_gate_pass", "hic_sequence_conflict"]
L = 10_000_000


def bed(path, rows):
    path.write_text("".join("\t".join(map(str, r)) + "\n" for r in rows))
    return path


def manifest(tmp_path, rows):
    """rows: (genome, kind, rule, tol_bp, hic_bins, truth) -> a truth.tsv"""
    p = tmp_path / "truth.tsv"
    p.write_text("genome\tkind\trule\ttol_bp\thic_bins\ttruth\tsource\tnotes\n"
                 + "".join("\t".join(map(str, r)) + "\tsynthetic\t.\n" for r in rows))
    return p


def run_dir(root, arm, gid, calls, hic_bin=".", done="0"):
    """A finished cenfuse run: DONE, cenfuse_summary.tsv, <gid>.centromeres.tsv; calls: (chrom, start, end, confidence),
    1-based inclusive like cenfuse's table, '.' for no call."""
    d = root / arm / gid
    d.mkdir(parents=True)
    if done is not None:
        (d / "DONE").write_text(done + "\n")
    (d / "cenfuse_summary.tsv").write_text(f"genome\thic_bin\thic_gate_pass\n{gid}\t{hic_bin}\t0\n")
    rows = ["\t".join(CEN_COLS)] + ["\t".join(map(str, (ch, L, s, e, c, "sequence_only", ".", "satellite")))
                                    for ch, s, e, c in calls]
    (d / f"{gid}.centromeres.tsv").write_text("\n".join(rows) + "\n")
    return d


def score(tmp_path, rule, tol, bins, truth, calls, hic_bin=".", kind="CENH3", files=None):
    """Score one synthetic run of genome 'g'; truth: BED rows (written to t.bed), or files: ';'-joined names."""
    tmp_path.mkdir(parents=True, exist_ok=True)
    if files is None:
        files = "." if truth is None else bed(tmp_path / "t.bed", truth).name
    T = S.load_manifest(manifest(tmp_path, [("g", kind, rule, tol, bins if bins else ".", files)]), data=tmp_path)
    run = run_dir(tmp_path / "runs", "S", "g", calls, hic_bin)
    d, hic = S.score_genome("S", "g", run, run / "g.centromeres.tsv", T["g"])
    return d.set_index("chrom"), hic


def test_overlap_needs_one_base(tmp_path):
    d, _ = score(tmp_path, "overlap", 0, None, [("c1", 100, 200), ("c2", 100, 200), ("c9", 0, 5)],
                 [("c1", 200, 300, "high"), ("c2", 201, 300, "high"), ("c3", 1, 50, "high")])
    assert d.at["c1", "outcome"] == "right_high"          # 1-based 200 = 0-based 199: one base inside [100, 200)
    assert (d.at["c2", "ok"], d.at["c2", "gap"], d.at["c2", "outcome"]) == ("W", 0, "conf_wrong")   # touching only
    assert "c3" not in d.index                            # no truth on c3: not scored


def test_gap_within_tolerance(tmp_path):
    d, _ = score(tmp_path, "gap", 500, None, [(c, 1000, 2000) for c in ("c1", "c2", "c3")] + [("c3", 5000, 6000)],
                 [("c1", 2501, 3000, "high"), ("c2", 2502, 3000, "medium"), ("c3", 5101, 5200, "low")])
    assert (d.at["c1", "gap"], d.at["c1", "ok"]) == (500, "R")
    assert (d.at["c2", "gap"], d.at["c2", "outcome"]) == (501, "conf_wrong")
    assert (d.at["c3", "truth"], d.at["c3", "outcome"]) == ("5000-6000", "right_low")   # the nearest interval
    assert set(d.tol) == {500}


def test_gap_zero_counts_touching(tmp_path):
    d, _ = score(tmp_path, "gap", 0, None, [("c1", 100, 200)], [("c1", 201, 300, "low")])
    assert (d.at["c1", "gap"], d.at["c1", "outcome"]) == (0, "right_low")


def test_four_hic_bins_with_hic_tol_without(tmp_path):
    truth = [("c1", 1_000_000, 2_000_000), ("c2", 1_000_000, 2_000_000)]
    calls = [("c1", 2_400_001, 2_500_000, "high"), ("c2", 2_400_002, 2_500_000, "high")]
    d, hic = score(tmp_path / "a", "gap", 1_000_000, 4, truth, calls, hic_bin=100_000)
    assert set(d.tol) == {400_000} and set(d.hic_bin) == {100_000} and hic["hic_bin"] == "100000"
    assert list(d.ok) == ["R", "W"]                       # gaps 400,000 and 400,001
    d, _ = score(tmp_path / "b", "gap", 1_000_000, 4, truth, calls)
    assert set(d.tol) == {1_000_000} and set(d.hic_bin) == {"."} and list(d.ok) == ["R", "R"]


def test_bin1_one_hic_bin_else_250kb(tmp_path):
    truth = [("c1", 100_000, 110_000), ("c2", 100_000, 110_000)]
    calls = [("c1", 120_001, 125_000, "high"), ("c2", 120_002, 125_000, "medium")]
    d, _ = score(tmp_path / "a", "bin1", 250_000, 1, truth, calls, hic_bin=10_000)
    assert set(d.tol) == {10_000} and list(d.outcome) == ["right_high", "conf_wrong"]
    d, _ = score(tmp_path / "b", "bin1", 250_000, 1, truth, calls)
    assert set(d.tol) == {250_000} and list(d.ok) == ["R", "R"]


def test_mammal_w_then_f_last_line_wins(tmp_path):
    bed(tmp_path / "g.W.bed", [("c1", 0, 100), ("c2", 1000, 2000)])
    bed(tmp_path / "g.F.bed", [("c1", 500, 600), ("c1", 700, 800)])
    d, _ = score(tmp_path, "mammal", 0, None, None, [("c1", 501, 600, "high"), ("c2", 1500, 1600, "high")],
                 kind="CENP-A", files="g.W.bed;g.F.bed")
    assert S.load_truth(S.load_manifest(tmp_path / "truth.tsv", data=tmp_path)["g"]) == {
        "c1": [(700, 800)], "c2": [(1000, 2000)]}
    assert (d.at["c1", "truth"], d.at["c1", "ok"]) == ("700-800", "W")    # F's last line replaced the others
    assert (d.at["c2", "truth"], d.at["c2", "ok"]) == ("1000-2000", "R")  # W kept where F has no line


def test_holocentric_every_confident_call_wrong(tmp_path):
    d, hic = score(tmp_path, "holo", 0, None, None, [("c1", 1, 100, "high"), ("c2", 1, 100, "medium"),
                                                      ("c3", 1, 100, "low"), ("c4", ".", ".", ".")], kind="holocentric")
    assert list(d.outcome) == ["conf_wrong", "conf_wrong", "wrong_low", "none"]
    assert list(d.ok) == ["W", "W", "W", "none"] and set(d.chance) == {0.0} and set(d.truth) == {"."}
    s = S.summarize(d.reset_index(), "S", "g", "holocentric", hic)
    assert (s["conf_wrong"], s["wrong_low"], s["none"], s["right"]) == (2, 1, 1, 0) and math.isnan(s["P"])


def test_censat_active_hor_only(tmp_path):
    rows = [("track name=censat",), ("#chrom", "start", "end", "name"), ("c1", 100, 200, "hor(S1C1H1L)"),
            ("c1", 1000, 2000, "active_hor(S1C1H1L)"), ("c2", 100, 200, "active_hor(S2C2H1L)"),
            ("c3", 100, 200, "hor(S3C3H1L)")]
    d, _ = score(tmp_path, "active", 0, None, rows, [("c1", 101, 200, "high"), ("c2", 101, 200, "high"),
                                                     ("c3", 101, 200, "high")], kind="CenSat active")
    assert d.at["c1", "outcome"] == "conf_wrong"          # on the inactive array only
    assert d.at["c2", "outcome"] == "right_high"
    assert "c3" not in d.index                            # no active array: no truth on c3


def test_no_call(tmp_path):
    d, _ = score(tmp_path, "gap", 0, None, [("c1", 100, 200)], [("c1", ".", ".", ".")])
    assert (d.at["c1", "ok"], d.at["c1", "outcome"], d.at["c1", "chance"], d.at["c1", "truth"]) == ("none", "none", 0.0, ".")
    assert math.isnan(d.at["c1", "gap"])


def test_chance_and_poisson_binomial_p(tmp_path):
    d, hic = score(tmp_path, "gap", 1000, None, [("c1", 1000, 2000), ("c2", 1000, 2000), ("c3", 0, 10)],
                   [("c1", 1001, 3000, "high"), ("c2", 1001, 3000, "high"), ("c3", 1, L, "low")])
    assert d.at["c1", "chance"] == pytest.approx((2000 + 1000 + 2 * 1000) / L)   # (call + truth + 2 tol) / length
    assert d.at["c3", "chance"] == 1.0                                          # at most 1
    s = S.summarize(d.reset_index(), "S", "g", "CENH3", hic)
    p = 5000 / L
    assert s["right"] == 3 and s["chance"] == round(2 * p + 1, 2)
    assert s["P"] == pytest.approx(p * p)                 # P(K >= 3): both small chances and the certain one


def test_missing_truth_skipped_with_warning(tmp_path, capsys):
    bed(tmp_path / "ok.bed", [("c1", 100, 200)])
    m = manifest(tmp_path, [("ok", "CENH3", "gap", 0, ".", tmp_path / "ok.bed"),
                            ("gone", "CENH3", "gap", 0, ".", tmp_path / "nowhere.bed")])
    for g in ("ok", "gone"):
        run_dir(tmp_path / "rel", "S", g, [("c1", 101, 200, "high")])
        run_dir(tmp_path / "rel", "T", g, [("c1", 101, 200, "high")])
    jobs = tmp_path / "jobs.tsv"
    jobs.write_text("".join(f"{a}\t{g}\t.\t1\t.\t.\n" for a in ("S", "T") for g in ("ok", "gone")))
    chrom, runs = S.score_jobs(jobs, tmp_path / "rel", S.load_manifest(m))
    assert [(r["arm"], r["genome"]) for r in runs] == [("S", "ok"), ("T", "ok")]
    err = capsys.readouterr().err
    assert err.count("warning: truth of gone missing") == 1 and "nowhere.bed" in err


def test_pending_and_failed_runs_not_scored(tmp_path, capsys):
    bed(tmp_path / "t.bed", [("c1", 100, 200)])
    m = manifest(tmp_path, [("g", "CENH3", "gap", 0, ".", tmp_path / "t.bed")])
    run_dir(tmp_path / "rel", "S", "g", [("c1", 101, 200, "high")], done=None)
    run_dir(tmp_path / "rel", "T", "g", [("c1", 101, 200, "high")], done="1")
    jobs = tmp_path / "jobs.tsv"
    jobs.write_text("S\tg\t.\t1\t.\t.\nT\tg\t.\t1\t.\t.\n")
    assert S.score_jobs(jobs, tmp_path / "rel", S.load_manifest(m)) == ([], [])
    err = capsys.readouterr().err
    assert "pending: S/g" in err and "failed: T/g exit 1" in err


def test_cli_writes_the_frozen_columns(tmp_path, capsys):
    bed(tmp_path / "t.bed", [("c1", 100, 200), ("c2", 100, 200)])
    m = manifest(tmp_path, [("g", "CENH3", "gap", 0, ".", tmp_path / "t.bed")])
    run_dir(tmp_path / "rel", "S", "g", [("c1", 101, 200, "high"), ("c2", 5001, 6000, "low")])
    (tmp_path / "jobs.tsv").write_text("S\tg\t.\t1\t.\t.\n")
    S.main([str(tmp_path / "jobs.tsv"), str(tmp_path / "rel"), str(tmp_path / "out"), "--manifest", str(m)])
    c = pd.read_csv(tmp_path / "out.chrom.tsv", sep="\t")
    r = pd.read_csv(tmp_path / "out.runs.tsv", sep="\t")
    assert list(c.columns) == CHROM_COLS and list(r.columns) == RUNS_COLS
    assert list(c.outcome) == ["right_high", "wrong_low"] and r.right.tolist() == [1]
    assert "right_high" in capsys.readouterr().out


def test_repository_manifest():
    T = S.load_manifest(paths.MANIFEST)
    assert {"agen", "sver", "chm13", "hg002", "celegans", "cneoA", "cneoG1"} <= set(T)
    for g, t in T.items():
        assert t.source not in ("", "."), g
        assert (t.rule == "mammal") == (t.kind == "CENP-A" and t.hic_bins is None and t.tol == 0
                                        and t.files[0].parent == paths.V2 / "data/truth"), g
        assert t.hic_bins is None or t.rule in ("gap", "bin1"), g
        assert t.rule != "bin1" or (t.tol, t.hic_bins) == (250_000, 1), g
        assert all(f.is_absolute() for f in t.files), g
    assert T["barley"].tol == 8_000_000 and T["barley"].hic_bins == 4
    assert (T["agen"].rule, T["agen"].tol, T["agen"].hic_bins) == ("gap", 1_000_000, None)


def test_local_moves_job_paths_under_the_data_root(monkeypatch, tmp_path):
    monkeypatch.setattr(paths, "DATA", tmp_path)
    assert paths.local(paths.DEFAULT_DATA / "v2/runs/x") == tmp_path / "v2/runs/x"
    assert paths.local("/elsewhere/x") == paths.Path("/elsewhere/x")
