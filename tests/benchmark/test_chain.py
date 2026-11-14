"""The chain's release checks on a tiny synthetic release (benchmark/paths.py, compare_rel.py, score_dev.py, rows.py,
headroom.py): a reference run right on two chromosomes, a release run that moves one call off its truth."""
import pandas as pd

from benchmark import compare_rel, headroom, paths, rows, score_dev
from benchmark import score as S

L = 10_000_000
COLS = ["chrom", "chrom_len", "start", "end", "confidence", "relation", "flags", "support", "seq_start", "seq_end"]


def run(root, arm, gid, calls):
    """A finished run; calls: (chrom, start, end, confidence, seq_start, seq_end), 1-based inclusive as cenfuse writes."""
    d = root / arm / gid
    d.mkdir(parents=True)
    (d / "DONE").write_text("0\n")
    (d / "cenfuse_summary.tsv").write_text(f"genome\thic_bin\thic_gate_pass\n{gid}\t.\t0\n")
    lines = ["\t".join(COLS)] + ["\t".join(map(str, (ch, L, s, e, c, "sequence_only", ".", "satellite", ss, se)))
                                 for ch, s, e, c, ss, se in calls]
    (d / f"{gid}.centromeres.tsv").write_text("\n".join(lines) + "\n")
    return d


def release(tmp_path):
    """Truth (BED, 0-based half-open) on c1 and c2 of 'rice', rule gap with tolerance 0; the release moves c2's call off
    its truth while its satellite line stays on it. Returns the manifest and the job table (column 6 = reference run)."""
    (tmp_path / "t.bed").write_text("c1\t1000\t2000\nc2\t1000\t2000\n")
    m = tmp_path / "truth.tsv"
    m.write_text("genome\tkind\trule\ttol_bp\thic_bins\ttruth\tsource\tnotes\n"
                 f"rice\tCENH3\tgap\t0\t.\t{tmp_path / 't.bed'}\tsynthetic\t.\n")
    right = [("c1", 1001, 2000, "high", 1001, 2000), ("c2", 1001, 2000, "high", 1001, 2000)]
    ref = run(tmp_path / "ref", "S", "rice", right)
    run(tmp_path / "rel", "S", "rice", [right[0], ("c2", 5001, 6000, "medium", 1001, 2000)])
    jobs = tmp_path / "jobs.tsv"
    jobs.write_text(f"S\trice\t.\t1\t.\t{ref}\n")
    return m, jobs


def test_shell_defaults_leave_the_optional_inputs_empty(capsys):
    paths.main(["--shell"])
    out = capsys.readouterr().out.splitlines()
    assert 'PREV_REDO="${PREV_REDO-}"' in out and 'EXTRA_REF="${EXTRA_REF-}"' in out
    assert f'RUNS="${{RUNS-{paths.RUNS}}}"' in out and paths.RUNS.is_relative_to(paths.DATA)
    paths.main([])
    assert any(x.startswith("TABLES\t") for x in capsys.readouterr().out.splitlines())


def test_compare_rel_and_score_dev_list_the_moved_call(tmp_path, capsys):
    m, jobs = release(tmp_path)
    compare_rel.main([str(jobs), str(tmp_path / "rel"), str(tmp_path / "cmp")])
    d = pd.read_csv(tmp_path / "cmp.diffs.tsv", sep="\t", dtype=str)
    assert set(d.chrom) == {"c2"} and set(d.column) == {"start", "end", "confidence"}
    s = pd.read_csv(tmp_path / "cmp.summary.tsv", sep="\t")
    assert (s.loc[0, "chromosomes"], s.loc[0, "call"], s.loc[0, "seq"]) == (2, 1, 0)
    score_dev.main([str(tmp_path / "rel"), str(jobs), str(tmp_path / "dev"), "--ref", str(tmp_path / "ref"),
                    "--manifest", str(m)])
    c = pd.read_csv(tmp_path / "dev.changes.tsv", sep="\t", dtype=str)
    assert c[["chrom", "old_outcome", "new_outcome"]].values.tolist() == [["c2", "right_high", "conf_wrong"]]
    assert "TOTAL dRight -1 dConfWrong 1" in capsys.readouterr().out


def test_rows_side_by_side_and_headroom(tmp_path, capsys):
    m, jobs = release(tmp_path)
    for root, prefix in (("rel", "new"), ("ref", "old")):
        S.main([str(jobs), str(tmp_path / root), str(tmp_path / prefix), "--manifest", str(m)])
    capsys.readouterr()
    rows.main([str(tmp_path / "new"), str(tmp_path / "old")])
    out = capsys.readouterr().out
    assert " **rice" in out and "2 | 1 (1/0/0) | cw 1" in out and "2 | 2 (2/0/0) | cw 0" in out
    headroom.main([str(tmp_path / "rel"), str(jobs), str(tmp_path / "head.tsv"), "--manifest", str(m)])
    h = pd.read_csv(tmp_path / "head.tsv", sep="\t", dtype=str).set_index("chrom")
    assert (h.at["c1", "final"], h.at["c2", "final"], h.at["c2", "seq"]) == ("R", "W", "R")   # a line right, the call not
