"""Unit tests: the decision table, confidence arithmetic, the gate statistic, loci merging, the chromosome set and
reading satellite-engine outputs."""
import itertools

import numpy as np
import pandas as pd
import pytest

from cenfuse import fuse, hicfuse, inputs, loci, seqcalls
from cenfuse.util import level_shift, poisson_binomial_sf

BIN = 250_000


def ev(seq=(1_000_000, 1_400_000), conf="high", peak=None, core=None, gate=True, strong=True, te=None, meth=None,
       sat=True, level="high"):
    s = dict(start=seq[0] if seq else -1, end=seq[1] if seq else -1, evidence="satellite" if seq else "none",
             confidence=conf if seq else "none", flags=".", main_family="satA")
    h = None
    if peak is not None:
        core = core or (peak, peak + BIN)
        h = dict(bin_size=BIN, gate_pass=gate, anch_strong=strong, anch_level=level, anch_peak_start=peak,
                 anch_peak_end=peak + BIN, anch_core_start=core[0], anch_core_end=core[1])
    t = dict(clade="CRM", start=te[0], end=te[1], decisive=True) if te else None
    sat_fn = (lambda a, b: max(0, min(b, seq[1]) - max(a, seq[0]))) if (sat and seq) else (lambda a, b: 0)
    return fuse.Evidence("c", 50_000_000, s, t, h, meth, sat_fn)


def test_no_hic_keeps_the_sequence_call():
    d = fuse.decide(ev())
    assert d["call"] == (1_000_000, 1_400_000) and d["relation"] == "sequence_only" and d["confidence"] == "high"
    assert d["support"] == ["satellite"]


def test_failed_gate_reports_only():
    d = fuse.decide(ev(peak=30_000_000, gate=False))
    assert d["call"] == (1_000_000, 1_400_000) and d["relation"] == "sequence_only" and not d["against"]


def test_agreeing_hic_raises_one_level():
    d = fuse.decide(ev(conf="low", peak=1_500_000))                  # one bin from the call: agrees
    assert d["relation"] == "concordant" and d["confidence"] == "medium" and "hic" in d["support"]


def test_strong_disagreement_moves_low_and_medium_calls():
    for conf in ("low", "medium"):
        d = fuse.decide(ev(conf=conf, peak=30_000_000, core=(29_500_000, 30_750_000)))
        assert d["call"] == (29_500_000, 30_750_000) and d["relation"] == "offset" and d["source"] == "hic"
        assert d["confidence"] == "medium"


def test_high_call_is_contested_without_a_second_line():
    d = fuse.decide(ev(peak=30_000_000))
    assert d["call"] == (1_000_000, 1_400_000) and d["relation"] == "contested" and d["confidence"] == "medium"
    assert d["against"] == ["hic"] and "hic_offset" in d["flags"]


def test_high_call_moves_with_te_or_methylation_at_the_hic_core():
    d = fuse.decide(ev(peak=30_000_000, te=(29_800_000, 30_400_000)))
    assert d["relation"] == "offset" and d["confidence"] == "high" and d["second"] == ["te"]
    m = {"hic_peak": dict(label="centromere_like", score=0.9), "sequence_call": dict(label="dip", score=0.1)}
    d = fuse.decide(ev(peak=30_000_000, meth=m))
    assert d["relation"] == "offset" and d["second"] == ["meth"] and "meth" in d["support"]


def test_methylation_at_both_places_is_not_a_second_line():
    m = {"hic_peak": dict(label="centromere_like", score=0.9), "sequence_call": dict(label="centromere_like", score=0.8)}
    assert fuse.decide(ev(peak=30_000_000, meth=m))["relation"] == "contested"


def test_weak_or_hub_peaks_never_move_a_call():
    d = fuse.decide(ev(conf="low", peak=30_000_000, strong=False))
    assert d["call"] == (1_000_000, 1_400_000) and "hic_disagrees_weak" in d["flags"] and not d["against"]


def test_functional_only_and_none():
    d = fuse.decide(ev(seq=None, peak=30_000_000, te=(29_900_000, 30_300_000)))
    assert d["relation"] == "functional_only" and d["confidence"] == "high" and d["support"] == ["te", "hic"]
    assert fuse.decide(ev(seq=None, peak=30_000_000, strong=False))["relation"] == "none"
    assert fuse.decide(ev(seq=None))["call"] is None


def test_methylation_alone_supports_but_never_moves():
    m = {"sequence_call": dict(label="centromere_like", score=0.8)}
    d = fuse.decide(ev(conf="low", meth=m))
    assert d["relation"] == "concordant" and d["confidence"] == "medium" and d["call"] == (1_000_000, 1_400_000)


def test_raises_are_capped():
    m = {"sequence_call": dict(label="centromere_like", score=0.8)}
    assert fuse.decide(ev(conf="low", peak=1_000_000, meth=m))["confidence"] == "high"
    assert level_shift("high", 2) == "high" and level_shift("low", -1) == "low" and level_shift("none", 1) == "none"


def test_poisson_binomial_matches_enumeration():
    p = [0.1, 0.3, 0.05, 0.5]
    for k in range(6):
        brute = sum(np.prod([q if b else 1 - q for q, b in zip(p, bits)])
                    for bits in itertools.product([0, 1], repeat=4) if sum(bits) >= k)
        assert poisson_binomial_sf(k, p) == pytest.approx(brute)


def _gate_table(agree, conf=None):
    n = len(agree)
    conf = conf or ["high"] * n
    return pd.DataFrame(dict(chrom=[f"c{i}" for i in range(n)], anch_peak_start=[0] * n, chance=[0.05] * n,
                             seq_start=[0] * n, seq_confidence=conf, anch_agrees_seq=agree))


def test_gate_is_leave_one_out_majority_and_significant():
    t = hicfuse.gate(_gate_table([True] * 6 + [False] * 2))
    assert t.gate_pass.all() and (t.gate_n == 7).all()
    assert list(t.gate_k) == [5] * 6 + [6] * 2              # a chromosome's own agreement never counts for it
    assert not hicfuse.gate(_gate_table([True] * 3 + [False] * 5)).gate_pass.any()     # no majority
    assert not hicfuse.gate(_gate_table([True, True, False])).gate_pass.any()          # too few references


def test_gate_falls_back_to_medium_references():
    t = hicfuse.gate(_gate_table([True] * 6, conf=["high", "medium", "medium", "medium", "medium", "medium"]))
    assert t.gate_pass.all() and (t.gate_n == 5).all()


def test_gate_takes_low_references_only_when_no_call_is_high_or_medium():
    t = hicfuse.gate(_gate_table([True] * 6 + [False] * 2, conf=["low"] * 8).assign(anch_flags=""))
    assert t.gate_pass.all() and (t.gate_refs == "low").all() and (t.gate_n == 7).all()
    t = hicfuse.gate(_gate_table([True] * 8, conf=["high"] + ["low"] * 7).assign(anch_flags=""))
    assert not t.gate_pass.any() and (t.gate_refs == "reference").all()     # one high call: low ones never count
    # Rabl nuclei: wrong calls at the chromosome ends would agree with the ends' trans cluster; with low references an
    # agreement at a telomeric (or hub, rDNA) peak does not count
    t = hicfuse.gate(_gate_table([True] * 8, conf=["low"] * 8).assign(anch_flags=["telomeric"] * 5 + [""] * 3))
    assert not t.gate_pass.any() and (t.gate_k <= 3).all()


def test_with_low_references_only_a_clean_peak_supports():
    e = ev(conf="low", peak=1_200_000)
    e.hic.update(gate_refs="low", anch_flags="")
    assert fuse.decide(e)["confidence"] == "medium"                        # two independent lines agree
    e.hic["anch_flags"] = "telomeric"
    d = fuse.decide(e)
    assert d["confidence"] == "low" and "hic" not in d["support"]
    e.hic["gate_refs"] = "reference"                                       # validated by high calls: as before
    assert fuse.decide(e)["confidence"] == "medium"


def test_loci_merge_and_kinds():
    m = [(0, 100, "satellite_cluster", "satellite", "subT", dict(role="subtelomeric", rdna=".")),
         (5_000, 9_000, "sequence_call", "satellite", "satellite", {}),
         (8_000, 12_000, "integrated_call", "call", "concordant", {}),
         (11_500, 13_000, "hic_anchored_core", "hic", "high", dict(strong=True, gate=True)),
         (50_000, 51_000, "satellite_cluster", "satellite", "R", dict(role="rDNA", rdna="45S_spacer")),
         (80_000, 81_000, "te_core", "te", "CRM", dict(marker=True))]
    groups = loci.merge(m, tol=0)
    assert [len(g) for g in groups] == [1, 3, 1, 1]
    assert [loci.kind(g) for g in groups] == ["subtelomeric_block", "centromere", "rDNA", "alt_centromere_candidate"]
    assert [len(g) for g in loci.merge(m, tol=30_000)] == [4, 1, 1]   # only the Hi-C member is padded (to 43 kb)
    dense = [(i * 50_000, i * 50_000 + 20_000, "satellite_cluster", "satellite", "f", dict(role="interstitial"))
             for i in range(20)]
    assert len(loci.merge(dense, tol=500_000)) == 20                   # dense satellites never chain


def test_chromosome_set_drops_unplaced_scaffolds():
    fai = {**{f"chr{i}": 50_000_000 + i for i in range(5)}, **{f"scaf{i}": 300_000 for i in range(40)}}
    assert set(inputs.chromosome_set(fai)) == {f"chr{i}" for i in range(5)}
    chrom_only = {f"chr{i}": 10_000_000 * (i + 1) for i in range(4)}
    assert inputs.chromosome_set(chrom_only) == chrom_only


def test_view_lists_chromosomes_only(tmp_path):
    src = tmp_path / "g"
    (src / "tes").mkdir(parents=True)
    fa = src / "A.fa"
    fa.write_text(">c1\nACGT\n>s1\nAC\n")
    (src / "A.fa.fai").write_text("c1\t4\t4\t4\t5\ns1\t2\t13\t2\t3\n")
    g = inputs.discover(src)[0]
    g.chroms = {"c1": 4}
    view = inputs.make_view([g], src, tmp_path / "view")
    assert (view / "A.fa").resolve() == fa.resolve() and (view / "A.fa.fai").read_text().startswith("c1\t4")
    assert "s1" not in (view / "A.fa.fai").read_text() and (view / "tes").is_symlink()


def test_satellite_outputs_are_read_0_based(tmp_path):
    (tmp_path / "G_centromere_summary.tsv").write_text(
        "chrom\tchrom_len\tcen_start\tcen_end\tevidence\tconfidence\tflags\tmain_family\tcen_sat_bp\tte_clade\t"
        "te_core_start\tte_core_end\nc1\t1000\t101\t200\tsatellite\thigh\t.\tf1\t90\tCRM\t151\t300\n"
        "c2\t900\t.\t.\tnone\tnone\t.\t.\t.\t.\t.\t.\n")
    (tmp_path / "G_satellite_arrays.gff").write_text(
        "##gff-version 3\nc1\tcenfuse\tsatellite_DNA\t101\t150\t.\t+\t.\tID=a;family=f1;role=centromeric;"
        "monomer_len=10;n_monomers=5;in_centromere=true\n")
    s = seqcalls.read(tmp_path, "G", {"c1": 1000, "c2": 900, "c3": 800})
    c1 = s.call("c1")
    assert (c1.start, c1.end, c1.te_start, c1.te_end) == (100, 200, 150, 300)
    assert s.call("c2") is None and s.call("c3") is None
    assert s.arrays.iloc[0][["start", "end", "role"]].tolist() == [100, 150, "centromeric"]


def test_telomere_unit_found_de_novo(tmp_path):
    import random
    from cenfuse import telomere
    rnd = random.Random(1)
    seqs = {}
    for i in range(4):
        body = "".join(rnd.choice("ACGT") for _ in range(600_000))
        left = "CCCTAAA" * 600 if i < 3 else ""                           # 4.2-kb arrays on 3 of 4 chromosomes
        seqs[f"c{i}"] = left + body + ("TTTAGGG" * 600 if i < 3 else "")
    fa = tmp_path / "t.fa"
    fa.write_text("".join(f">{c}\n{s}\n" for c, s in seqs.items()))
    import pysam
    pysam.faidx(str(fa))
    r = telomere.detect(fa, {c: len(s) for c, s in seqs.items()})
    assert r["unit"] == telomere.canonical("TTTAGGG") == "AAACCCT"
    assert sorted((c, side) for c, side, *_ in r["ends"]) == sorted((f"c{i}", s) for i in range(3) for s in ("start", "end"))
    assert all(abs(x[4] - 4200) <= 7 and x[5] == 0 for x in r["ends"])
    none = tmp_path / "n.fa"
    none.write_text("".join(f">{c}\n{s[4200:-4200] if c != 'c3' else s}\n" for c, s in seqs.items()))
    pysam.faidx(str(none))
    assert telomere.detect(none, {c: len(s) - (8400 if c != "c3" else 0) for c, s in seqs.items()})["unit"] is None


def test_unexplained_loci_are_not_called_satellite():
    assert loci.kind([(0, 10, "meth_candidate", "meth", "P1", dict(label="dip"))]) == "meth_candidate"
    assert loci.kind([(0, 10, "hic_denovo_core", "hic", "low", {})]) == "hic_peak"


def test_peaks_within_the_method_resolution_neither_support_nor_oppose():
    d = fuse.decide(ev(peak=1_400_000 + 3 * BIN))                  # 3 bins past a high call: 'near'
    assert d["relation"] == "sequence_only" and d["confidence"] == "high" and not d["against"]
    assert "hic_near" in d["flags"] and "hic" not in d["support"]
    d = fuse.decide(ev(conf="low", peak=1_400_000 + 3 * BIN))       # nor does it move a low call
    assert d["call"] == (1_000_000, 1_400_000)
    d = fuse.decide(ev(peak=1_400_000 + 5 * BIN))                  # beyond 4 bins: contested
    assert d["relation"] == "contested"



def test_review_te_second_line_is_symmetric_and_decisive_only():
    # CRM core 1.1 Mb from the call, 4 Mb from the Hi-C peak: it does not support the peak (500-kb bins)
    e = ev(seq=(10_000_000, 10_400_000), peak=16_000_000, core=(14_000_000, 18_000_000), te=(11_500_000, 12_000_000))
    e.hic["bin_size"] = 500_000
    d = fuse.decide(e)
    assert d["relation"] == "contested" and "te" not in d["second"]
    e = ev(peak=30_000_000, te=(29_800_000, 30_400_000))
    e.te["decisive"] = False                                  # a flags-only clade is never a second line
    assert fuse.decide(e)["relation"] == "contested"


def test_review_no_methylation_data_at_the_call_is_not_a_second_line():
    for lab in ("no_data", "no_signature_data", ""):
        m = {"hic_peak": dict(label="centromere_like", score=0.9), "sequence_call": dict(label=lab, score=float("nan"))}
        assert fuse.decide(ev(peak=30_000_000, meth=m))["relation"] == "contested", lab


def test_review_insignificant_peaks_never_support():
    d = fuse.decide(ev(conf="low", peak=1_500_000, level="none"))
    assert d["relation"] == "sequence_only" and d["confidence"] == "low" and "hic" not in d["support"]


def test_review_contested_with_methylation_support_keeps_its_level():
    m = {"sequence_call": dict(label="centromere_like", score=0.8), "hic_peak": dict(label="dip", score=0.1)}
    d = fuse.decide(ev(peak=30_000_000, meth=m))
    assert d["relation"] == "contested" and d["confidence"] == "high"


def test_review_a_core_reaching_the_call_does_not_disagree():
    d = fuse.decide(ev(conf="medium", peak=3_000_000, core=(900_000, 3_250_000)))   # peak 6 bins off, core covers call
    assert d["call"] == (1_000_000, 1_400_000) and d["relation"] == "sequence_only" and "hic_near" in d["flags"]


def test_review_telomere_longest_array_wins_over_a_nearer_decoy():
    from cenfuse import telomere
    unit = "TTTAGGG"
    real = unit * 600                                          # 4.2 kb
    decoy = unit * 5                                           # 35 bp, 2 kb nearer the end
    seq = "A" * 3000 + real + "C" * 2000 + decoy + "G" * 100
    a = telomere.end_array(seq, telomere.canonical(unit), "end")
    assert a is not None and abs(a[2] - 4200) <= 7 and a[0] == 3000


def test_short_telomere_counts_only_at_the_sequence_end(tmp_path):
    import random
    import pysam
    from cenfuse import telomere
    unit, cu = "TTAGGG", telomere.canonical("TTAGGG")
    short = unit * 27                                          # 162 bp (Neurospora NC12 telomeres: 157-207 bp)
    body = "".join(random.Random(2).choice("ACGT") for _ in range(20_000))
    a = telomere.end_array(body + short, cu, "end")
    assert a is not None and a[1] == len(body) + len(short) and abs(a[2] - 162) <= 6
    assert telomere.end_array(short + body, cu, "start")[0] == 0
    assert telomere.end_array(body[:15_000] + short + body[15_000:], cu, "end") is None    # 5 kb inside: interstitial
    assert telomere.end_array(body + unit * 10, cu, "end") is None                          # 60 bp: too short
    rnd = random.Random(3)                                     # a genome whose every telomere is short
    seqs = {f"c{i}": "CCCTAA" * 27 + "".join(rnd.choice("ACGT") for _ in range(450_000)) + short for i in range(3)}
    fa = tmp_path / "s.fa"
    fa.write_text("".join(f">{c}\n{s}\n" for c, s in seqs.items()))
    pysam.faidx(str(fa))
    r = telomere.detect(fa, {c: len(s) for c, s in seqs.items()})
    assert r["unit"] == cu and len(r["ends"]) == 6 and all(x[5] == 0 for x in r["ends"])


def test_interstitial_telomeric_arrays(tmp_path, monkeypatch):
    import random
    import pysam
    from cenfuse import telomere
    rnd = random.Random(4)
    body = "".join(rnd.choice("ACGT") for _ in range(600_000))
    its, short, near_end = "TTAGGG" * 150, "CCCTAA" * 25, "TTAGGG" * 100     # 900 bp, 150 bp, 600 bp
    seq = (body[:50_000] + near_end + body[50_000:299_000] + its + body[299_000:450_000] + short
           + body[450_000:])
    fa = tmp_path / "i.fa"
    fa.write_text(f">c1\n{seq}\n")
    pysam.faidx(str(fa))
    monkeypatch.setattr(telomere, "SCAN_CHUNK", 100_000)              # the 900-bp array spans a chunk edge
    start = 50_000 + 600 + 249_000                                    # 299,600-300,500: across the 300-kb chunk edge
    assert seq[start:start + 900] == its and start < 300_000 < start + 900
    out = telomere.interstitial(fa, {"c1": len(seq)}, telomere.canonical("TTAGGG"))
    assert len(out) == 1                                               # not the terminal-zone or the 150-bp array
    c, s, e, bp = out[0]
    assert (c, s, e) == ("c1", start, start + 900) and abs(bp - 900) <= 6
    m = (s, e, "telomere_interstitial", "telomere", "AACCCT", dict(bp=bp))
    assert loci.kind([m]) == "interstitial_telomeric"
    assert loci.kind([m, (0, 10, "telomere_array", "telomere", "AACCCT", dict(bp=10))]) == "telomere"


def test_review_smoothing_stays_inside_chromosome_pairs():
    import scipy.sparse as sp
    from cenfuse.hic.bins import Bins
    from cenfuse.hicfuse import oe_dense
    b = Bins({"A": 5 * BIN, "B": 5 * BIN, "C": 5 * BIN}, BIN)
    W = np.ones((15, 15))
    W[5, 10:15] = W[10:15, 5] = 50.0                          # B's first bin contacts C strongly
    mean = np.ones((3, 3))
    np.fill_diagonal(mean, np.nan)
    S = oe_dense(sp.csr_matrix(W), b, np.ones(15, bool), mean, np.ones(3, bool))
    assert S[4, 12] == pytest.approx(1.0)                      # A's last bin x C: untouched by B's first bin
    assert S[5, 12] > 10


def test_genome_level_conflict_between_hic_and_sequence():
    n = 8
    t = pd.DataFrame(dict(seq_start=[0] * n, callable=[True] * n, denovo_confidence=["high"] * n,
                          denovo_agrees_seq=[False] * n, gate_pass=[False] * n))
    assert fuse.sequence_conflict(t)["conflict"]
    t["denovo_agrees_seq"] = [True] * 6 + [False] * 2
    assert not fuse.sequence_conflict(t)["conflict"]
    t["denovo_agrees_seq"] = False
    t["denovo_confidence"] = ["low"] * n                      # an unconfident prior-free Hi-C is no conflict
    assert not fuse.sequence_conflict(t)["conflict"]
