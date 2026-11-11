"""The CpG line (cdr.py): dips, the permutation test, both input formats, the gate, and the line on a small genome with
planted dips (and its shifted negative control)."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pysam
import pytest

sys.path.insert(0, str(Path(__file__).parent))

from cenfuse import cdr  # noqa: E402
from cenfuse.meth import io as mio  # noqa: E402

B = cdr.BIN


def plateau(n=200, level=0.8, noise=0.03, dip=None, seed=1):
    rng = np.random.default_rng(seed)
    y = level + noise * rng.standard_normal(n)
    if dip:
        y[dip[0]:dip[1]] = 0.15 + noise * rng.standard_normal(dip[1] - dip[0])
    return y


def flat(n=200):
    return np.full(n, 0.8), np.ones(n, bool)


def test_dip_is_a_run_below_half_its_background():
    y = plateau(dip=(90, 110))
    bg, m = flat()
    d = cdr.dips(y, bg, m)
    assert len(d) == 1 and d[0][1:] == (90, 110) and d[0][0] < 0.3        # depth: methylation over background
    assert cdr.dips(plateau(dip=(90, 92)), bg, m) == []                  # 10 kb: narrower than MIN_RUN
    m2 = m.copy()
    m2[95:] = False                                                       # the run ends where the arrays end
    assert cdr.dips(y, bg, m2)[0][1:] == (90, 95)
    two = plateau(dip=(90, 110))
    two[40:50] = 0.3                                                      # a shallower dip
    assert [x[1:] for x in cdr.dips(two, bg, m)] == [(90, 110), (40, 50)]  # deepest first


def test_one_missing_bin_is_bridged_two_break_the_run():
    y = plateau(dip=(90, 110))
    bg, m = flat()
    y[100] = np.nan
    assert [x[1:] for x in cdr.dips(y, bg, m)] == [(90, 110)]
    y[101] = np.nan
    assert sorted(x[1:] for x in cdr.dips(y, bg, m)) == [(90, 100), (102, 110)]


def test_profile_pools_calls_over_three_bins(monkeypatch):
    monkeypatch.setattr(cdr, "MIN_SITES", 2)
    b = cdr.bins(np.array([100, 5100, 10100, 10200, 99_999]), np.array([1.0, 0.0, 2.0, 0.0, 1.0]),
                 np.array([1.0, 1.0, 2.0, 2.0, 1.0]), 3)                  # the last row lies beyond the sequence
    assert b.shape == (3, 3) and list(b[2]) == [1, 1, 2]
    y = cdr.profile(b)
    assert np.allclose(y, [1 / 2, 3 / 6, 2 / 5])                          # (methylated, valid) calls of bins i-1..i+1


def test_array_bins_merge_overlaps():
    a = pd.DataFrame(dict(chrom=["c", "c", "c", "d"], start=[0, 5000, 12000, 0], end=[7500, 9000, 13000, 15000],
                          role="centromeric"))
    assert list(cdr.array_bins(a, "c", 3)) == [True, True, False]        # 5000/5000, 4000/5000, 1000/5000


def _tabix(path: Path, rows):
    raw = path.with_suffix("")
    raw.write_text("".join("\t".join(map(str, r)) + "\n" for r in rows))
    pysam.tabix_compress(str(raw), str(path), force=True)
    pysam.tabix_index(str(path), preset="bed", force=True)
    return path


def test_both_formats_read(tmp_path):
    bg = _tabix(tmp_path / "a.bedgraph.gz", [("c1", 100, 102, 80.0), ("c1", 200, 202, 20.0)])
    t = cdr.Track(bg)
    assert t.kind == "bedgraph" and t.scale == 100.0
    p, m, n = t.sites("c1", 0, 1000)
    assert np.allclose(m, [0.8, 0.2]) and np.allclose(n, [1, 1])                    # a bedGraph site = one call
    bm = _tabix(tmp_path / "b.bed.gz", [
        ("c1", 100, 101, "m", 10, "+", 100, 101, "0,0,0", 10, 60.0, 5, 4, 1, 0, 0, 0, 0),
        ("c1", 100, 101, "h", 10, "+", 100, 101, "0,0,0", 10, 10.0, 1, 4, 5, 0, 0, 0, 0),   # 5hmC row: skipped
        ("c1", 300, 301, "m", 2, "+", 300, 301, "0,0,0", 2, 50.0, 1, 1, 0, 0, 0, 0, 0),     # 2 calls: kept
        ("c1", 400, 401, "m", 0, "+", 400, 401, "0,0,0", 0, 0.0, 0, 0, 0, 0, 0, 0, 0)])     # no call: dropped
    t = cdr.Track(bm)
    p, m, n = t.sites("c1", 0, 1000)
    assert t.kind == "bedmethyl" and list(p) == [100, 300] and np.allclose(m, [6, 1]) and np.allclose(n, [10, 2])
    assert len(t.sites("nochrom", 0, 10)[0]) == 0


def test_encode_bedmethyl_reads_like_modkit(tmp_path):
    """ENCODE bedMethyl (11 columns: coverage, percent) gives the calls of the modkit rows above."""
    en = _tabix(tmp_path / "e.bed.gz", [("c1", 100, 101, ".", 10, "+", 100, 101, "0,0,0", 10, 60),
                                        ("c1", 300, 301, ".", 2, "+", 300, 301, "0,0,0", 2, 50),
                                        ("c1", 400, 401, ".", 0, "+", 400, 401, "0,0,0", 0, 0)])
    t = cdr.Track(en)
    p, m, n = t.sites("c1", 0, 1000)
    assert t.kind == "encode" and list(p) == [100, 300] and np.allclose(m, [6, 1]) and np.allclose(n, [10, 2])


def test_unreadable_table_is_named(tmp_path, caplog):
    cx = _tabix(tmp_path / "cx.bed.gz", [("c1", 100 + i, 101 + i, "+", 3, 2, "CG") for i in range(20)])
    arrays = pd.DataFrame(dict(chrom=["c1"], start=[0], end=[1000], role=["centromeric"]))
    with caplog.at_level("WARNING"):
        cdr.run([cx], {"c1": (0, 1000, "high")}, arrays)
    assert any("no row could be read" in r.message and str(cx) in r.message for r in caplog.records)


def test_cg_only_samples_table(tmp_path):
    f = _tabix(tmp_path / "s.bedgraph.gz", [("c1", 1, 3, 0.5)])
    (tmp_path / "m.tsv").write_text(f"sample\tCG\ns1\t{f.name}\n")
    with pytest.raises(mio.InputError):
        mio.read_samples(tmp_path / "m.tsv")                     # the methylation engine still needs CHG
    s = mio.read_samples(tmp_path / "m.tsv", require_chg=False)
    assert s[0].files == {"CG": f}
    (tmp_path / "n.tsv").write_text("sample\tCHH\ns1\tx\n")
    with pytest.raises(mio.InputError):
        mio.read_samples(tmp_path / "n.tsv", require_chg=False)


def test_gate():
    t = pd.DataFrame(dict(chrom=[f"c{i}" for i in range(6)], confidence=["high"] * 6, pos_p=[0.01] * 5 + [0.5]))
    g = cdr.gate(t)
    assert g.gate_pass.all() and (g.gate_n == 5).all()
    assert not cdr.gate(t.assign(pos_p=0.5)).gate_pass.any()
    g = cdr.gate(t.iloc[:3])                                     # two references each: no test
    assert not g.gate_pass.any() and (g.gate_n == 2).all()
    one = t.assign(pos_p=[0.01] + [0.5] * 5)                     # one positional dip in five: chance
    assert not cdr.gate(one).gate_pass.any()


L, CALL, NULLS = 4_000_000, (1_000_000, 2_000_000), (3_000_000, 3_400_000)


def small_genome(tmp_path, shift=False, n=6, seed=3, everywhere=False, moved=(), keep_array=1.0):
    """n chromosomes of L bp: a centromeric array over CALL with a 100-kb dip at 1.4 Mb, an interstitial array over
    NULLS without one; a CpG every 100 bp (bedGraph fractions). shift: positions moved by L/2 (the negative control).
    everywhere: methylated 0.8 along the whole chromosome with a 100-kb dip every 400 kb (dips are not positional).
    moved: chromosomes whose dip lies in the interstitial array (3.1-3.2 Mb) instead. keep_array: the share of the
    arrays' 50-kb blocks that keep their rows (a coverage hole elsewhere)."""
    rng = np.random.default_rng(seed)
    rows, arrays, calls = [], [], {}
    for i in range(n):
        c = f"c{i}"
        pos = np.arange(50, L, 100)
        y = 0.3 + 0.05 * rng.standard_normal(len(pos))           # unmethylated-ish background
        sat = ((pos >= CALL[0]) & (pos < CALL[1])) | ((pos >= NULLS[0]) & (pos < NULLS[1]))
        y[sat] = 0.8 + 0.05 * rng.standard_normal(sat.sum())
        if c in moved:
            y[(pos >= 3_100_000) & (pos < 3_200_000)] = 0.15
        else:
            y[(pos >= 1_400_000) & (pos < 1_500_000)] = 0.15
        if everywhere:
            y = 0.8 + 0.05 * rng.standard_normal(len(pos))
            y[(pos % 400_000) >= 300_000] = 0.15
        p = (pos + L // 2) % L if shift else pos
        keep = ~sat | (rng.random(len(pos) // 500 + 1)[pos // 50_000] < keep_array)
        o = np.argsort(p[keep])
        rows += [(c, int(a), int(a) + 2, round(float(v), 3)) for a, v in zip(p[keep][o], y[keep][o])]
        arrays += [(c, CALL[0], CALL[1], "centromeric"), (c, NULLS[0], NULLS[1], "interstitial")]
        calls[c] = (CALL[0], CALL[1], "high")
    f = _tabix(tmp_path / f"cg{int(shift)}.bedgraph.gz", rows)
    return f, calls, pd.DataFrame(arrays, columns=["chrom", "start", "end", "role"])


def test_run_finds_planted_kinetochores(tmp_path):
    f, calls, arrays = small_genome(tmp_path)
    res = cdr.run([f], calls, arrays)
    t = res["table"].set_index("chrom")
    assert t.gate_pass.all() and (t.dip_depth < cdr.FRAC).all() and (t.pos_p <= cdr.DIP_P).all()
    assert ((t.kt_start - 1_400_000).abs() <= B).all() and ((t.kt_end - 1_500_000).abs() <= B).all()
    assert (t.alt_start < 0).all() and (t.array_data == 1).all()
    assert (t.kt_meth < 0.3).all() and (t.array_meth > 0.7).all()
    out = cdr.columns(pd.DataFrame(dict(chrom=list(calls))), res)
    assert (out.cpg_gate == "pass").all() and (out.cpg_kt_start.astype(int) == t.kt_start.to_numpy() + 1).all()


def test_shifted_track_fails_the_gate(tmp_path):
    f, calls, arrays = small_genome(tmp_path, shift=True)
    t = cdr.run([f], calls, arrays)["table"]
    assert not t.gate_pass.any() and (t.kt_start < 0).all()


def test_dips_everywhere_are_not_positional(tmp_path):
    """Dips in every array, also far from the calls (methylation from elsewhere makes them too): the gate fails."""
    f, calls, arrays = small_genome(tmp_path, everywhere=True)
    t = cdr.run([f], calls, arrays)["table"]
    assert (t.dip_depth < cdr.FRAC).all() and not t.gate_pass.any()


def test_malformed_rows_are_skipped(tmp_path):
    f = _tabix(tmp_path / "m.bedgraph.gz", [("c1", 100, 102, 0.5), ("c1", 200, 202, "x"), ("c1", 300, 302, 0.7)])
    t = cdr.Track(f)
    p, m, n = t.sites("c1")
    assert list(p) == [100, 300] and t.bad == 1


def test_columns_without_the_line():
    out = cdr.columns(pd.DataFrame(dict(chrom=["a", "b"])), None)
    assert (out[cdr.COLS] == ".").all().all()


def test_pipeline_cg_only_changes_no_call(tmp_path):
    """A CG-only samples table (no tes/ needed): the CpG line runs, every call stays as without methylation."""
    import fuse_synth
    from cenfuse.pipeline import Run
    truth = fuse_synth.make(tmp_path / "a", meth=False, plants=False)
    rows, fa = [], pysam.FastaFile(str(truth["dir"] / "Syn.fa"))
    seqs = {c: bytearray(fa.fetch(c).encode()) for c in fuse_synth.CHROMS}
    fa.close()
    for c in fuse_synth.CHROMS:
        for p in range(50, fuse_synth.L, 400):                   # a CpG at every row (the assembly check)
            seqs[c][p:p + 2] = b"CG"
            rows.append((c, p, p + 2, 0.8 if abs(p - fuse_synth.bp(fuse_synth.F[c])) > 50_000 else 0.1))
    with open(truth["dir"] / "Syn.fa", "w") as f:
        for c in fuse_synth.CHROMS:
            f.write(f">{c}\n{seqs[c].decode()}\n")
    pysam.faidx(str(truth["dir"] / "Syn.fa"))
    Run(truth["dir"], tmp_path / "o0", 4, satellite_dir=truth["sat"]).run()
    (truth["dir"] / "meth").mkdir(exist_ok=True)
    f = _tabix(truth["dir"] / "meth" / "cg.bedgraph.gz", rows)
    (truth["dir"] / "meth" / "Syn.tsv").write_text(f"sample\tCG\ns1\t{f.name}\n")
    Run(truth["dir"], tmp_path / "o1", 4, satellite_dir=truth["sat"]).run()
    a = pd.read_csv(tmp_path / "o0" / "Syn.centromeres.tsv", sep="\t", dtype=str)
    b = pd.read_csv(tmp_path / "o1" / "Syn.centromeres.tsv", sep="\t", dtype=str)
    keep = ["chrom", "start", "end", "relation", "confidence", "support", "against", "flags"]
    assert a[keep].equals(b[keep]) and (a[cdr.COLS] == ".").all().all()
    assert b.cpg_gate.isin(["pass", "fail", "untested"]).all()
    s = pd.read_csv(tmp_path / "o1" / "cenfuse_summary.tsv", sep="\t")
    assert "cpg_gate_pass" in s and s.lines.iat[0] == "satellite,te,hic,meth"


def test_bedgraph_assembly_check(tmp_path):
    """Rows on CpGs pass (also one base early, as some public tracks); rows elsewhere do not."""
    seq = ("AT" * 20 + "CG") * 500
    (tmp_path / "g.fa").write_text(">c1\n" + seq + "\n")
    pysam.faidx(str(tmp_path / "g.fa"))
    cpg = [i for i in range(len(seq) - 1) if seq[i:i + 2] == "CG"]
    win = [("c1", 0, len(seq))]
    on = _tabix(tmp_path / "on.bedgraph.gz", [("c1", i, i + 2, 0.5) for i in cpg])
    early = _tabix(tmp_path / "early.bedgraph.gz", [("c1", i - 1, i + 1, 0.5) for i in cpg])
    off = _tabix(tmp_path / "off.bedgraph.gz", [("c1", i + 7, i + 9, 0.5) for i in cpg])
    assert cdr.bedgraph_on_cpg(on, tmp_path / "g.fa", win)[0] == 1.0
    assert cdr.bedgraph_on_cpg(early, tmp_path / "g.fa", win)[0] == 1.0
    assert cdr.bedgraph_on_cpg(off, tmp_path / "g.fa", win)[0] < 0.1
    assert cdr.track_kind(on) == "bedgraph"


def test_a_malformed_first_row_loses_only_itself(tmp_path):
    f = _tabix(tmp_path / "f.bedgraph.gz", [("c1", 100, 102, ""), ("c1", 200, 202, 0.5), ("c1", 300, 302, 0.7)])
    t = cdr.Track(f)
    assert list(t.sites("c1")[0]) == [200, 300] and t.bad == 1


def test_low_coverage_dip_is_found_by_pooling_calls(tmp_path):
    """About one ONT call per CpG (MMU8 thinned to 1/16). Per-site fractions with a 3-call floor keep almost no site;
    pooling each bin's calls (methylated over valid) still shows the 100-kb dip at 1.4 Mb, on every chromosome."""
    rng = np.random.default_rng(5)
    rows, arrays, calls = [], [], {}
    for i in range(8):
        c = f"c{i}"
        pos = np.arange(50, L, 100)
        frac = np.full(len(pos), 0.3)
        sat = ((pos >= CALL[0]) & (pos < CALL[1])) | ((pos >= NULLS[0]) & (pos < NULLS[1]))
        frac[sat] = 0.8
        frac[(pos >= 1_400_000) & (pos < 1_500_000)] = 0.15
        n = rng.poisson(1.0, len(pos))
        m = rng.binomial(n, frac)
        keep = n > 0
        rows += [(c, int(a), int(a) + 1, "m", int(k), "+", int(a), int(a) + 1, "0,0,0", int(k), round(100 * x / k, 2),
                  int(x), int(k - x), 0, 0, 0, 0, 0) for a, k, x in zip(pos[keep], n[keep], m[keep])]
        arrays += [(c, CALL[0], CALL[1], "centromeric"), (c, NULLS[0], NULLS[1], "interstitial")]
        calls[c] = (CALL[0], CALL[1], "high")
    f = _tabix(tmp_path / "low.bed.gz", rows)
    t = cdr.run([f], calls, pd.DataFrame(arrays, columns=["chrom", "start", "end", "role"]))["table"].set_index("chrom")
    assert t.gate_pass.all() and ((t.kt_start - 1_400_000).abs() <= 2 * B).all() and ((t.kt_end - 1_500_000).abs() <= 2 * B).all()


def test_contest_points_to_the_other_array(tmp_path):
    """c0's call holds no dip while its interstitial array holds one; the other chromosomes' calls carry theirs (the
    gate passes, leave-one-out): the line contests c0 and names the other array's dip; no other call is contested."""
    f, calls, arrays = small_genome(tmp_path, moved=("c0",))
    res = cdr.run([f], calls, arrays)
    t = res["table"].set_index("chrom")
    assert t.gate_pass.all() and t.at["c0", "kt_start"] < 0 and t.at["c0", "pos_p"] == 1.0
    assert abs(t.at["c0", "alt_start"] - 3_100_000) <= B and abs(t.at["c0", "alt_end"] - 3_200_000) <= B
    assert (t.drop(index="c0").alt_start < 0).all()
    ch = res["chroms"]["c0"]
    assert ch.elsewhere(*CALL) == (t.at["c0", "alt_start"], t.at["c0", "alt_end"])
    assert res["chroms"]["c1"].elsewhere(*CALL) is None
    out = cdr.columns(pd.DataFrame(dict(chrom=list(calls))), res).set_index("chrom")
    assert int(out.at["c0", "cpg_alt_start"]) == t.at["c0", "alt_start"] + 1 and out.at["c1", "cpg_alt_start"] == "."


def test_pbcpg_rows_are_read_with_their_coverage(tmp_path):
    """pb-CpG-tools (9 columns: score, type, coverage, estimated modified and unmodified counts): calls, not a fraction."""
    f = _tabix(tmp_path / "p.bed.gz", [("c1", 100, 101, 60.0, "Total", 10, 6, 4, 60.0),
                                       ("c1", 300, 301, 25.0, "hap1", 4, 1, 3, 25.0)])
    t = cdr.Track(f)
    p, m, n = t.sites("c1", 0, 1000)
    assert t.kind == "pbcpg" and list(p) == [100, 300] and np.allclose(m, [6, 1]) and np.allclose(n, [10, 4])


def test_coverage_hole_in_the_arrays_is_reported(tmp_path, caplog):
    """Most array rows missing (as in tracks piled up at MAPQ 60: reads in arrays map at MAPQ 20-59): a warning."""
    f, calls, arrays = small_genome(tmp_path, keep_array=0.3)
    with caplog.at_level("WARNING"):
        res = cdr.run([f], calls, arrays)
    assert res["array_data"] < cdr.HOLE and any("MAPQ" in r.message for r in caplog.records)
    caplog.clear()
    f, calls, arrays = small_genome(tmp_path)
    with caplog.at_level("WARNING"):
        res = cdr.run([f], calls, arrays)
    assert res["array_data"] >= cdr.HOLE and not any("MAPQ" in r.message for r in caplog.records)
