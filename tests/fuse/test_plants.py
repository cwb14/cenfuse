"""Provenance and per-plant scoring on tiny synthetic cases (docs/plants.md)."""
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pysam
import pytest

from cenfuse import plants, provenance
from cenfuse.inputs import hic_libraries
from cenfuse.methfuse import merge_sites
from cenfuse.util import CenfuseError

STEP = 100


def _genome(tmp: Path, L: int = 300_000, seed: int = 1) -> tuple[Path, np.ndarray]:
    """Random A/T with a CAG every STEP bp (a CHG cytosine on each strand)."""
    rng = np.random.default_rng(seed)
    seq = np.where(rng.random(L) < 0.5, ord("A"), ord("T")).astype(np.uint8)
    pos = np.arange(50, L - 3, STEP)
    seq[pos], seq[pos + 1], seq[pos + 2] = ord("C"), ord("A"), ord("G")
    fa = tmp / "g.fa"
    fa.write_bytes(b">c1\n" + seq.tobytes() + b"\n")
    pysam.faidx(str(fa))
    return fa, pos


def _bedmethyl(tmp: Path, name: str, rows: list) -> Path:
    path = tmp / f"{name}.CHG.bedmethyl"
    with open(path, "w") as f:
        for x, strand in sorted(rows):
            f.write(f"c1\t{x}\t{x + 1}\tm\t10\t{strand}\t{x}\t{x + 1}\t255,0,0\t10\t50.00\t5\t5\t0\t0\t0\t0\t0\n")
    pysam.tabix_index(str(path), preset="bed", force=True)
    return Path(f"{path}.gz")


def test_motifs():
    seq = "TTCAGTTCTGTTCGTT"
    assert list(provenance.symmetric_chg(seq, 100)) == [102, 107]
    pos = np.array([102, 104, 107, 109, 112])
    strand = np.array(["+", "-", "+", "-", "+"])
    assert list(provenance.on_motif(seq, 100, pos, strand, "CHG")) == [True, True, True, True, False]
    assert list(provenance.on_motif(seq, 100, np.array([112]), np.array(["+"]), "CG")) == [True]


def test_fixed_differences_need_both_strands(tmp_path):
    fa, pos = _genome(tmp_path)
    rng = np.random.default_rng(2)
    inner = rng.choice(pos[10:-10], 110, replace=False)
    fixed, plus_lost, minus_lost = inner[:30], inner[30:70], inner[70:]
    rows = [(int(p), "+") for p in pos if p not in set(fixed) | set(plus_lost)]
    rows += [(int(p) + 2, "-") for p in pos if p not in set(fixed) | set(minus_lost)]
    path = _bedmethyl(tmp_path, "S", rows)
    g = provenance.sample_genotype(("S", path, fa, [("c1", 0, 300_000)], {"c1": 0}))
    assert g.on_motif == 1.0
    assert set(g.fixed) == set(fixed.tolist()) and (g.plus_only, g.minus_only) == (40, 40)
    assert len(g.tested) >= len(pos) - 4
    noise = g.noise * len(g.tested)
    assert 0.4 < noise < 1.0                               # 40 x 40 / n ~ 0.5 both-missing by chance


def test_another_assembly_stops_the_run(tmp_path):
    fa, pos = _genome(tmp_path)
    good = _bedmethyl(tmp_path, "good", [(int(p), "+") for p in pos])
    bad = _bedmethyl(tmp_path, "bad", [(int(p) + 1, "+") for p in pos])
    smp = lambda n, p: SimpleNamespace(name=n, files={"CHG": p})
    assert provenance.check_assembly([smp("good", good)], fa, [("c1", 0, 300_000)]) == []
    with pytest.raises(CenfuseError, match="another assembly"):
        provenance.check_assembly([smp("bad", bad)], fa, [("c1", 0, 300_000)])


def _encode(tmp: Path, name: str, rows: list) -> Path:
    """ENCODE bedMethyl (11 columns: coverage, percent)."""
    path = tmp / f"{name}.encode.bed"
    with open(path, "w") as f:
        for x, strand in sorted(rows):
            f.write(f"c1\t{x}\t{x + 1}\t.\t10\t{strand}\t{x}\t{x + 1}\t0,0,0\t10\t50\n")
    pysam.tabix_index(str(path), preset="bed", force=True)
    return Path(f"{path}.gz")


def test_encode_cg_files_are_checked_against_the_assembly(tmp_path):
    rng = np.random.default_rng(3)
    seq = np.where(rng.random(300_000) < 0.5, ord("A"), ord("T")).astype(np.uint8)
    pos = np.arange(50, 299_990, STEP)
    seq[pos], seq[pos + 1] = ord("C"), ord("G")
    fa = tmp_path / "cg.fa"
    fa.write_bytes(b">c1\n" + seq.tobytes() + b"\n")
    pysam.faidx(str(fa))
    smp = lambda n, p: SimpleNamespace(name=n, files={"CG": p})
    good = _encode(tmp_path, "good", [(int(p), "+") for p in pos])
    bad = _encode(tmp_path, "bad", [(int(p) + 2, "+") for p in pos])
    assert provenance.check_assembly([smp("good", good)], fa, [("c1", 0, 300_000)]) == []
    with pytest.raises(CenfuseError, match="another assembly"):          # before: no row read, the check passed over
        provenance.check_assembly([smp("bad", bad)], fa, [("c1", 0, 300_000)])


def test_chg_files_must_be_modkit(tmp_path):
    fa, pos = _genome(tmp_path)
    mk = _bedmethyl(tmp_path, "mk", [(int(p), "+") for p in pos[:20]])
    en = _encode(tmp_path, "en", [(int(p), "+") for p in pos[:20]])
    s1, s2, s3 = (SimpleNamespace(name="s1", files={"CHG": mk}), SimpleNamespace(name="s2", files={"CHG": en, "CG": en}),
                  SimpleNamespace(name="s3", files={"CHG": en}))
    kept, dropped = provenance.drop_non_modkit_chg([s1, s2, s3])
    assert [s.name for s in kept] == ["s1", "s2"] and s2.files == {"CG": en} and s1.files == {"CHG": mk}
    assert [(n, w) for n, _, w in dropped] == [("s2", "ENCODE bedMethyl (coverage, percent)"),
                                               ("s3", "ENCODE bedMethyl (coverage, percent)")]


def _gt(name, tested, fixed, loss=450):          # 450 single-strand losses in 100,000 sites: ~2 both-missing by chance
    return provenance.Genotype(name, 0, 1.0, np.array(sorted(tested)), np.array(sorted(fixed)), loss, loss)


def test_genotype_verdicts():
    rng = np.random.default_rng(3)
    sites = np.arange(100_000)
    noise = lambda: set(rng.choice(sites, 2, replace=False).tolist())
    shared = set(rng.choice(sites, 300, replace=False).tolist())
    ref = _gt("ref", sites, noise())
    o1 = _gt("o1", sites, shared | set(rng.choice(sites, 900, replace=False).tolist()) | noise() | set(range(10)))
    o2 = _gt("o2", sites, shared | set(rng.choice(sites, 900, replace=False).tolist()) | noise())
    v = provenance.genotype_verdicts([ref, o1, o2], provenance.private_tests([ref, o1, o2]))
    assert v["ref"]["verdict"] == "reference" and v["o1"]["verdict"] == v["o2"]["verdict"] == "other"
    one = provenance.genotype_verdicts([o1], provenance.private_tests([o1]))
    assert one["o1"]["verdict"] == "unknown" and "one methylation sample" in one["o1"]["reason"]
    twin = _gt("twin", sites, set(o1.fixed.tolist()))      # same genotype as o1: neither can be placed
    v2 = provenance.genotype_verdicts([o1, twin], provenance.private_tests([o1, twin]))
    assert v2["o1"]["verdict"] == v2["twin"]["verdict"] == "unknown"
    small = _gt("small", sites[:300], set(), loss=1)       # matches the reference, too few sites to tell
    v3 = provenance.genotype_verdicts([small, o1, o2], provenance.private_tests([small, o1, o2]))
    assert v3["small"]["verdict"] == "unknown" and "power" in v3["small"]["reason"]


def test_missing_segments_and_call_absent():
    n = 30
    flat = lambda: {"c1": np.full(n, 100.0)}
    prof = {("meth", "A"): flat(), ("meth", "B"): flat(), ("hic", "H"): flat()}
    prof[("meth", "A")]["c1"][10:16] = 1.0                 # 6 windows absent
    prof[("meth", "A")]["c1"][20:22] = 1.0                 # too short
    prof[("meth", "B")]["c1"][25:28] = prof[("meth", "A")]["c1"][25:28] = 1.0
    prof[("hic", "H")]["c1"][25:28] = 1.0                  # absent in every input: not judged
    skip = {"c1": np.zeros(n, bool)}
    skip["c1"][[0, n - 1, 12]] = True                      # a skipped window inside the run is transparent
    seg = provenance.missing_segments(prof, skip, 1_000_000)
    a = seg[("meth", "A")]
    assert len(a) == 1 and (a[0]["start"], a[0]["end"], a[0]["windows"]) == (10_000_000, 16_000_000, 5)
    assert a[0]["covered_by"] == "meth:B" and seg[("meth", "B")] == [] and seg[("hic", "H")] == []    # same type
    calls = {"c1": (12_200_000, 12_800_000)}
    assert provenance.call_absent(prof, skip, seg, calls, 1_000_000) == {("A", "c1"): True}
    assert provenance.call_absent(prof, skip, seg, {"c1": (17_100_000, 17_600_000)}, 1_000_000) == {}
    for k in prof:                                         # a window no input covers cannot be judged: no break
        prof[k]["c1"][13] = 1.0
    a2 = provenance.missing_segments(prof, skip, 1_000_000)[("meth", "A")]
    assert len(a2) == 1 and (a2[0]["start"], a2[0]["end"], a2[0]["windows"]) == (10_000_000, 16_000_000, 4)
    prof[("meth", "A")]["c1"][14] = 300.0                 # a pile-up (> 2x the median) cannot be judged either
    a4 = provenance.missing_segments(prof, skip, 1_000_000)[("meth", "A")]
    assert [(x["start"], x["end"], x["windows"]) for x in a4] == [(10_000_000, 16_000_000, 3)]
    prof[("meth", "A")]["c1"][14] = 60.0                  # the plant has coverage: the run ends there
    a3 = provenance.missing_segments(prof, skip, 1_000_000)[("meth", "A")]
    assert [(x["start"], x["end"]) for x in a3] == []      # 10-12 and 15 are each shorter than three windows
    prof[("meth", "B")]["c1"][10:16] = 1.0                 # only the Hi-C input covers A's segment now: another
    assert provenance.missing_segments(prof, skip, 1_000_000)[("meth", "A")] == []      # data type is no evidence


def test_satellite_dominated_windows_are_not_judged(tmp_path):
    fa, _ = _genome(tmp_path, L=6_000_000)
    arrays = pd.DataFrame(dict(chrom=["c1", "c1", "c1"], start=[1_100_000, 1_500_000, 3_200_000],
                               end=[1_400_000, 1_900_000, 3_220_000]))
    share = provenance.satellite_share(arrays, "c1", 6_000_000, 1_000_000)
    assert np.allclose(share[:4], [0, 0.7, 0, 0.02])
    sk = provenance.skip_windows({"c1": 6_000_000}, 1_000_000, {"c1": [(4_500_000, 4_600_000)]}, fa, arrays)["c1"]
    assert list(sk) == [True, True, False, False, True, True]      # ends, satellite > half, the call


def test_combine_verdicts():
    geno = {"A": dict(verdict="reference", reason="r"), "B": dict(verdict="other", reason="o"),
            "C": dict(verdict="unknown", reason="u")}
    seg = {("meth", "A"): [dict(chrom="c1", start=0, end=3_000_000)]}
    v = provenance.combine(["A", "B", "C"], ["B", "hic"], geno, seg)
    assert v[("meth", "A")]["verdict"] == "other" and "conflict" in v[("meth", "A")]["reason"]
    assert v[("meth", "B")]["verdict"] == "other" and v[("hic", "B")]["verdict"] == "other"
    assert "by name" in v[("hic", "B")]["reason"]
    assert v[("meth", "C")]["verdict"] == v[("hic", "hic")]["verdict"] == "unknown"


def _sites(n_chrom=10, plants_=("P1", "P2", "P3", "P4"), rng=None):
    """Two sites per chromosome: the call (S1, inside) and an alternative satellite site (S2); shared centromeres."""
    rng = rng or np.random.default_rng(5)
    rows = []
    for k in range(1, n_chrom + 1):
        c = f"c{k}"
        for site, (s, e), origin, base in (("S1", (9_500_000, 10_500_000), "sequence_call", 1.4),
                                           ("S2", (29_500_000, 30_500_000), "satellite:sat", 0.05)):
            for p in plants_:
                sc = base + rng.normal(0, 0.08)
                rows.append(dict(chrom=c, site=f"{c}:{site}", start=s, end=e, mid=(s + e) // 2, origins=origin,
                                 sample=p, label="centromere_like" if sc >= 0.5 else "dip", score=sc, value=-0.1,
                                 with_data=1))
    return pd.DataFrame(rows)


def _set(st, chrom, site, plant, score):
    i = (st.chrom == chrom) & (st.site == f"{chrom}:{site}") & (st["sample"] == plant)
    st.loc[i, "score"] = score
    st.loc[i, "label"] = "centromere_like" if score >= 0.5 else ("dip" if np.isfinite(score) else "no_signature_data")


def test_meth_verdicts():
    st = _sites()
    calls = {f"c{k}": (9_800_000, 10_200_000) for k in range(1, 11)}
    _set(st, "c3", "S1", "P2", 0.1)                         # P2 moved its c3 centromere to the satellite site
    _set(st, "c3", "S2", "P2", 1.4)
    for p in ("P1", "P2", "P3", "P4"):                      # c5: every plant at the alternative (the call is odd)
        _set(st, "c5", "S1", p, 0.1)
        _set(st, "c5", "S2", p, 1.4)
    _set(st, "c7", "S1", "P4", np.nan)                      # no data at the call
    tab, tested = plants.score([f"c{k}" for k in range(1, 11)], calls,
                               {("meth", p): dict(verdict="other") for p in ("P1", "P2", "P3", "P4")}, st,
                               {("P3", "c9"): True}, {})
    v = tab.set_index(["plant", "chrom"]).verdict
    assert v[("P2", "c3")] == "shifted" and v[("P3", "c9")] == "call_absent" and v[("P4", "c7")] == "no_data"
    assert all(v[(p, "c5")] == "possible_shift" for p in ("P1", "P2", "P3", "P4"))
    rest = v.drop([("P2", "c3"), ("P3", "c9"), ("P4", "c7")] + [(p, "c5") for p in ("P1", "P2", "P3", "P4")])
    assert set(rest) == {"shared"}
    r = tab[(tab.plant == "P2") & (tab.chrom == "c3")].iloc[0]
    assert r.loss_q <= 0.05 and r.gain_q <= 0.05 and r.shift_bp == 29_500_000 - 10_200_000
    counts = plants.call_counts(tab, ["c3"]).iloc[0]
    assert (counts.plants_shared, counts.plants_shifted) == (3, 1)


def test_own_candidates_carry_no_gain_test():
    st = _sites()
    st.loc[st.site.str.endswith("S2") & (st["sample"] == "P2"), "origins"] = "meth:P2"
    s = plants.site_tests(st, {f"c{k}": (9_800_000, 10_200_000) for k in range(1, 11)})
    assert not s[(s["sample"] == "P2") & s.site.str.endswith("S2")].test.eq("gain").any()
    assert s[(s["sample"] == "P1") & s.site.str.endswith("S2")].test.eq("gain").all()


def test_hic_verdicts():
    b = 250_000
    row = dict(callable=True, gate_pass=True, anch_level="high", anch_strong=True, anch_peak_start=20 * b,
               anch_peak_end=21 * b, anch_core_start=19 * b, anch_core_end=22 * b)
    call = (5 * b, 6 * b)
    assert plants.hic_verdict(row, call, b)["hic_verdict"] == "shifted"
    ref = dict(anch_peak_start=20 * b, anch_peak_end=21 * b, anch_level="high")
    assert plants.hic_verdict(row, call, b, ref, b)["hic_verdict"] == "shared"      # the reference peaks there too
    assert plants.hic_verdict({**row, "gate_pass": False}, call, b)["hic_verdict"] == "uninformative"
    assert plants.hic_verdict({**row, "anch_peak_start": 5 * b, "anch_peak_end": 6 * b}, call, b)["hic_verdict"] == "shared"
    assert plants.hic_verdict({**row, "anch_strong": False}, call, b)["hic_verdict"] == "unclear"
    assert plants.combine(["shifted", "shared"]) == "conflict" and plants.combine(["no_data", "uninformative"]) == "no_data"
    assert plants.combine(["possible_shift", None]) == "possible_shift" and plants.combine(["unclear", "shared"]) == "shared"


def test_merge_sites():
    W = 1_000_000
    g = merge_sites([("call", 10_000_000, 10_400_000), ("hic", 10_300_000, 10_550_000),
                     ("sat", 14_000_000, 14_100_000), ("meth:P1", 14_200_000, 14_300_000)], W, 50_000_000)
    assert [x["origins"] for x in g] == [["call", "hic"], ["sat", "meth:P1"]]
    assert g[0]["mid"] == 10_200_000 and (g[0]["start"], g[0]["end"]) == (9_700_000, 10_700_000)


def test_hic_library_table(tmp_path):
    h = tmp_path / "hic"
    h.mkdir()
    (h / "g.pairs.gz").write_bytes(b"")
    assert [(x.name, x.path.name) for x in hic_libraries(h, "g")] == [("hic", "g.pairs.gz")]
    (h / "P3.pairs.gz").write_bytes(b"")
    (h / "g.tsv").write_text("sample\tfile\nref\tg.pairs.gz\nP3\tP3.pairs.gz\n")
    assert [x.name for x in hic_libraries(h, "g")] == ["ref", "P3"]
    for bad, msg in (("sample\tpath\nx\tP3.pairs.gz\n", "columns"), ("sample\tfile\nx y\tP3.pairs.gz\n", "bad sample"),
                     ("sample\tfile\nx\tnone.pairs.gz\n", "no file"),
                     ("sample\tfile\nx\tP3.pairs.gz\nx\tg.pairs.gz\n", "duplicate")):
        (h / "g.tsv").write_text(bad)
        with pytest.raises(CenfuseError, match=msg):
            hic_libraries(h, "g")


def test_support_uses_a_bin_with_signature_data():
    """Review finding: the most extreme bin was chosen first, then signature data was asked for there; a sample with
    data elsewhere in the window was labeled no_signature_data."""
    from cenfuse.meth import params as mp
    from cenfuse.methfuse import support
    names, nb = ["LTR:CRM", "other"], 10
    num, vok = np.zeros((2, nb)), np.full((2, nb), 1000.0)
    vok[0, 3], num[1, 3] = 10, -300                   # bin 3: the deepest dip, but too few CRM calls
    num[0, 6] = -200                                  # bin 6: CRM loss like the anchors'
    sig = pd.DataFrame(dict(context=["CHG"], cls=["LTR:CRM"], mean_dev=[-0.2], direction=["down"]))
    prep = SimpleNamespace(geo=SimpleNamespace(w=1))
    r = support({"CHG": (num, vok)}, sig, names, 0, nb * mp.BIN, nb, prep, "down")
    assert r["label"] == "centromere_like" and r["bin_start"] == 6 * mp.BIN and r["score"] == pytest.approx(1.0)
    vok[0, :] = 10                                    # no bin has signature data: the extreme bin, labeled so
    r = support({"CHG": (num, vok)}, sig, names, 0, nb * mp.BIN, nb, prep, "down")
    assert r["label"] == "no_signature_data" and r["bin_start"] == 3 * mp.BIN


def test_call_absent_needs_the_call_inside_the_segment():
    """Review finding: a call one window beyond a segment, where no input has reads, was called absent."""
    n = 30
    prof = {("meth", "A"): {"c1": np.full(n, 100.0)}, ("meth", "B"): {"c1": np.full(n, 100.0)}}
    prof[("meth", "A")]["c1"][10:16] = 1.0
    prof[("meth", "B")]["c1"][15] = 1.0                    # window 15: no input has reads
    skip = {"c1": np.zeros(n, bool)}
    seg = provenance.missing_segments(prof, skip, 1_000_000)
    assert [(s["start"], s["end"]) for s in seg[("meth", "A")]] == [(10_000_000, 15_000_000)]
    assert provenance.call_absent(prof, skip, seg, {"c1": (15_200_000, 15_800_000)}, 1_000_000) == {}
    assert provenance.call_absent(prof, skip, seg, {"c1": (12_200_000, 12_800_000)}, 1_000_000) == {("A", "c1"): True}


def test_the_deciding_library_is_the_one_the_gate_trusts():
    """Review finding: among libraries of unknown origin the most reads won, whatever the gate said."""
    from cenfuse.pipeline import deciding_library
    lib = lambda gate, reads: dict(table=pd.DataFrame(dict(gate_pass=gate)), coverage=np.array([reads]))
    hics = {"a": lib([True] * 5 + [False] * 5, 1_000_001), "b": lib([True] * 10, 1_000_000)}
    v = {("hic", "a"): dict(verdict="unknown"), ("hic", "b"): dict(verdict="unknown")}
    assert deciding_library(hics, v) == "b"
    v[("hic", "a")] = dict(verdict="reference")
    assert deciding_library(hics, v) == "a"                 # a reference library first
    v[("hic", "a")] = v[("hic", "b")] = dict(verdict="other")
    assert deciding_library(hics, v) is None


def test_no_tests_without_a_call():
    st = _sites()
    calls = {f"c{k}": (9_800_000, 10_200_000) for k in range(1, 11)}
    calls["c4"] = None
    s = plants.site_tests(st, calls)
    assert not s[s.chrom == "c4"].test.isin(["loss", "gain"]).any() and s[s.chrom == "c5"].test.eq("loss").any()
