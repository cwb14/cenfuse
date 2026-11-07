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


def test_overruled_array_is_read_0_based(tmp_path):
    """seq_overruled: the satellite array a CRM core replaced, read like the other intervals (1-based in the summary,
    0-based here), -1 where nothing was replaced, on a chromosome without a call, and on one with no row."""
    head = ("chrom\tchrom_len\tcen_start\tcen_end\tevidence\tconfidence\tflags\tmain_family\tcen_sat_bp\tte_clade\t"
            "te_core_start\tte_core_end\talt_start\talt_end\toverruled_start\toverruled_end\toverruled_family\n")
    (tmp_path / "G_centromere_summary.tsv").write_text(
        head + "c1\t1000\t151\t300\tretrotransposon\tmedium\tsatellite_away_from_core\t.\t0\tCRM\t151\t300\t.\t.\t"
        "51\t120\tG_sat01_171\n"
        "c2\t900\t101\t200\tsatellite\thigh\t.\tf1\t90\t.\t.\t.\t.\t.\t.\t.\t.\n"
        "c3\t900\t.\t.\tnone\tnone\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\n")
    (tmp_path / "G_satellite_arrays.gff").write_text("##gff-version 3\n")
    calls = seqcalls.read(tmp_path, "G", {"c1": 1000, "c2": 900, "c3": 900, "c4": 800}).calls.set_index("chrom")
    assert calls.loc["c1", ["overruled_start", "overruled_end"]].tolist() == [50, 120]
    assert calls.loc["c1", ["alt_start", "alt_end"]].tolist() == [-1, -1]            # never the alternative cluster
    for c in ("c2", "c3", "c4"):
        assert calls.loc[c, ["overruled_start", "overruled_end"]].tolist() == [-1, -1], c


def test_read_calls_without_overruled_columns(tmp_path):
    """A summary made before these columns existed (satellite engine 0.8, or an older --satellite-dir) is read, and
    nothing is overruled anywhere."""
    (tmp_path / "S.tsv").write_text(
        "chrom\tchrom_len\tcen_start\tcen_end\tevidence\tconfidence\tflags\tmain_family\tcen_sat_bp\tte_clade\t"
        "te_core_start\tte_core_end\talt_start\talt_end\nc1\t1000\t151\t300\tretrotransposon\tmedium\t"
        "satellite_away_from_core\t.\t0\tCRM\t151\t300\t.\t.\nc2\t900\t.\t.\tnone\tnone\t.\t.\t.\t.\t.\t.\t.\t.\n")
    calls = seqcalls.read_calls(tmp_path / "S.tsv", {"c1": 1000, "c2": 900})
    assert list(calls.columns[-2:]) == ["overruled_start", "overruled_end"]
    assert calls.overruled_start.tolist() == [-1, -1] and calls.overruled_end.tolist() == [-1, -1]
    assert calls.start.tolist() == [150, -1]                                          # the rest is read as before


def test_overruled_array_is_reported_and_never_moves_the_call(tmp_path):
    """The array a CRM core replaced goes to seq_overruled_*, from the summary to the fused table, and nowhere else:
    fusion lets validated Hi-C move a call to the alternative cluster (hic_chose_alt), so a replaced array kept in
    seq_alt_* would change calls. Here the strong Hi-C peak sits on the replaced array: the decision is the one made
    without it."""
    from cenfuse.pipeline import Run
    head = ("chrom\tchrom_len\tcen_start\tcen_end\tevidence\tconfidence\tflags\tmain_family\tcen_sat_bp\tte_clade\t"
            "te_core_start\tte_core_end\talt_start\talt_end\toverruled_start\toverruled_end\toverruled_family\n")
    row = ("c1\t50000000\t1000001\t1400000\tretrotransposon\thigh\tsatellite_away_from_core\t.\t0\tCRM\t1000001\t"
           "1400000\t.\t.\t{}\n")
    tables, decisions = [], []
    for k, over in enumerate(("30000001\t30250000\tG_sat01_171", ".\t.\t.")):
        d = tmp_path / str(k)
        d.mkdir()
        (d / "G_centromere_summary.tsv").write_text(head + row.format(over))
        (d / "G_satellite_arrays.gff").write_text("##gff-version 3\n")
        evs = Run.evidences(None, seqcalls.read(d, "G", {"c1": 50_000_000}), {}, None)
        evs[0].hic = dict(bin_size=BIN, gate_pass=True, anch_strong=True, anch_level="high", anch_peak_start=30_000_000,
                          anch_peak_end=30_000_000 + BIN, anch_core_start=30_000_000, anch_core_end=30_000_000 + BIN)
        decisions.append(fuse.decide(evs[0]))
        tables.append(fuse.table(evs, decisions[-1:]))
    with_over, without = decisions
    same = ("call", "relation", "confidence", "source", "support", "against", "flags", "second")
    assert {k: with_over[k] for k in same} == {k: without[k] for k in same}
    assert "hic_chose_alt" not in with_over["flags"]
    assert with_over["call"] == (1_000_000, 1_400_000) and with_over["relation"] == "contested"
    r = tables[0].iloc[0]
    assert (r.seq_overruled_start, r.seq_overruled_end) == (30_000_001, 30_250_000)   # 1-based like seq_start
    assert (r.seq_alt_start, r.seq_alt_end, r.offset_bp) == (".", ".", ".")
    r = tables[1].iloc[0]
    assert (r.seq_overruled_start, r.seq_overruled_end) == (".", ".")
    cols = list(tables[0].columns)
    assert cols[cols.index("seq_alt_end") + 1:][:2] == ["seq_overruled_start", "seq_overruled_end"]


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


# ---------------------------------------------------------------------------------- support tier (0.10, hicfuse.gate)
def test_support_tier_is_significant_agreement_without_a_majority():
    # 5 of 12 references agree at chance 0.05 each, no majority: among the other 11, an agreeing chromosome sees 4
    # (P 0.0015) and the others 5 -> support tier on every chromosome
    t = hicfuse.gate(_gate_table([True] * 5 + [False] * 7))
    assert not t.gate_pass.any() and t.gate_support.all()
    assert list(t.gate_k) == [4] * 5 + [5] * 7
    # 4 of 12: an agreeing chromosome sees 3 of 11 (P 0.015 > 0.01): its own agreement never counts for it
    t = hicfuse.gate(_gate_table([True] * 4 + [False] * 8))
    assert not t.gate_support[:4].any() and t.gate_support[4:].all()
    # GATE_MIN_REF agreements is the floor, however small the chance: 2 among the others is not a pattern
    t = hicfuse.gate(_gate_table([True] * 3 + [False] * 9).assign(chance=0.001))
    assert not t.gate_support[:3].any() and t.gate_support[3:].all()
    # a passed gate is never also the support tier; too few references: neither
    assert not hicfuse.gate(_gate_table([True] * 6 + [False] * 2)).gate_support.any()
    assert not hicfuse.gate(_gate_table([True, True, False])).gate_support.any()
    # at a high chance agreement is not significant: 4 of 12 at chance 0.3 -> fail
    t = _gate_table([True] * 4 + [False] * 8).assign(chance=0.3)
    assert not hicfuse.gate(t).gate_support.any()
    # only agreements at significant peaks count (randomized Hi-C: peaks at unmappable arrays' edges, level low)
    t = _gate_table([True] * 5 + [False] * 7).assign(anch_level=["low"] * 5 + ["high"] * 7)
    assert not hicfuse.gate(t).gate_support.any()
    t = _gate_table([True] * 5 + [False] * 7).assign(anch_level=["medium"] * 5 + ["none"] * 7)
    assert hicfuse.gate(t).gate_support.all()


def _weak(conf="medium", peak=1_500_000, strong=True, level="high", **kw):
    e = ev(conf=conf, peak=peak, gate=False, strong=strong, level=level, **kw)
    e.hic.update(gate_support=True)
    return e


def test_support_tier_needs_the_chromosomes_own_peak_significant():
    d = fuse.decide(_weak(conf="low", level="low"))
    assert d["confidence"] == "low" and "hic_support" not in d["flags"] and "hic" not in d["support"]
    assert fuse.decide(_weak(conf="low", level="medium"))["confidence"] == "medium"


def test_support_tier_raises_an_agreeing_call_one_level():
    for conf, up in (("low", "medium"), ("medium", "high"), ("high", "high")):
        d = fuse.decide(_weak(conf=conf))
        assert d["relation"] == "concordant" and d["confidence"] == up and "hic_support" in d["flags"]
        assert "hic" in d["support"] and not d["against"]
    assert fuse.gate_status(dict(gate_pass=False, gate_support=True, gate_n=12)) == "support"


def test_support_tier_never_opposes_moves_or_places():
    for conf in ("low", "medium", "high"):                 # a strong peak far away: kept, not contested
        d = fuse.decide(_weak(conf=conf, peak=30_000_000))
        assert d["call"] == (1_000_000, 1_400_000) and d["relation"] == "sequence_only"
        assert d["confidence"] == conf and not d["against"] and d["flags"] == []
    d = fuse.decide(_weak(seq=None, peak=30_000_000))      # no sequence call: Hi-C places nothing on this tier
    assert d["call"] is None and d["relation"] == "none"


def test_support_tier_with_low_references_needs_a_clean_peak():
    e = _weak(conf="low")
    e.hic.update(gate_refs="low", anch_flags="telomeric")
    d = fuse.decide(e)
    assert d["confidence"] == "low" and "hic_support" not in d["flags"]
    e.hic["anch_flags"] = ""
    assert fuse.decide(e)["confidence"] == "medium"


def test_support_tier_adds_to_methylation_within_the_cap():
    m = {"sequence_call": dict(label="centromere_like", score=0.9), "hic_peak": dict(label="", score=np.nan)}
    d = fuse.decide(_weak(conf="low", meth=m))
    assert d["confidence"] == "high" and set(d["support"]) >= {"hic", "meth"}


# ------------------------------------------------------------------- Hi-C alone where nothing can test it (0.10)
def _alone(seq=None, gate_n=0, conf="high", flags="", fold=None):
    e = ev(seq=seq, peak=-1, gate=False)                           # no anchored peak (no references to anchor)
    e.hic.update(anch_peak_start=-1, anch_peak_end=-1, anch_core_start=-1, anch_core_end=-1, anch_level="none",
                 gate_n=gate_n, gate_support=False, gate_genome_untested=gate_n < 3,
                 denovo_peak_start=150_000, denovo_peak_end=160_000,
                 denovo_core_start=140_000, denovo_core_end=170_000, denovo_confidence=conf, denovo_flags=flags)
    e.fold = fold
    return e


def test_hic_alone_places_a_low_call_where_no_test_is_possible():
    d = fuse.decide(_alone())
    assert d["call"] == (140_000, 170_000) and d["relation"] == "functional_only" and d["confidence"] == "low"
    assert d["flags"] == ["hic_denovo_only"] and d["support"] == ["hic"] and d["source"] == "hic"
    assert fuse.decide(_alone(flags="no_rabl"))["flags"] == ["hic_denovo_only"]     # Rabl test: uninformative


def test_hic_alone_never_where_a_sequence_call_or_a_test_exists():
    d = fuse.decide(_alone(seq=(1_000_000, 1_400_000), conf="high"))
    assert d["call"] == (1_000_000, 1_400_000) and d["relation"] == "sequence_only" and not d["flags"]
    assert fuse.decide(_alone(gate_n=3))["call"] is None              # the gate ran (and failed): no call
    # 0.10.1: a passed axis line (fold or alignment) that does not pair here leaves the untested low call in place
    # (0.10: no call; N. crassa LG VI); where it pairs, the two 3D lines give a medium call instead
    folded = dict(status="pass", agrees=False, core_start=5_000_000, core_end=6_000_000, point=5_500_000)
    d = fuse.decide(_alone(fold=folded))
    assert d["flags"] == ["hic_denovo_only", "fold_offset"] and d["confidence"] == "low"   # the fold's place is said
    assert d["against"] == ["fold"] and "fold" not in d["support"]
    paired = dict(folded, agrees=True, core_start=100_000, core_end=200_000)
    d = fuse.decide(_alone(fold=paired))
    assert d["flags"] == ["fold_pair"] and d["confidence"] == "medium" and d["call"] == (140_000, 170_000)


@pytest.mark.parametrize("conf, flags", [("medium", ""), ("low", ""), ("high", "telomeric"), ("high", "hub"),
                                         ("high", "rDNA"), ("high", "cycle_unstable"), ("high", "not_converged")])
def test_hic_alone_needs_a_clean_high_prior_free_call(conf, flags):
    assert fuse.decide(_alone(conf=conf, flags=flags))["call"] is None


# ------------------------------------------------------------------- long-unit arrays are features unless a line agrees
def _long(e):
    """e with its sequence call a long-unit tandem array, as the satellite engine reports it: low, not a satellite."""
    e.seq.update(evidence="long_unit_tandem", main_family="long_unit", confidence="low")
    e.sat_bp = lambda a, b: 0
    return e


def test_a_long_unit_array_no_line_supports_is_no_call():
    """24 such calls from the FASTA alone: 2 right, 16 wrong, 6 without truth."""
    d = fuse.decide(_long(ev()))
    assert d["call"] is None and d["relation"] == "none" and d["confidence"] == "none" and d["source"] == "none"
    assert d["flags"] == ["long_unit_unsupported"]


def test_the_array_of_a_dropped_call_stays_in_the_sequence_columns():
    e = _long(ev())
    r = fuse.row(e, fuse.decide(e))
    assert (r["seq_start"], r["seq_end"], r["seq_evidence"]) == (1_000_001, 1_400_000, "long_unit_tandem")
    assert (r["start"], r["end"], r["relation"]) == (".", ".", "none")


def test_a_long_unit_call_stands_where_tes_or_validated_hic_or_methylation_agree():
    near = dict(clade="Tekay", start=1_100_000, end=1_300_000, decisive=False)
    e = _long(ev())
    e.te = near
    assert fuse.decide(e)["call"] == (1_000_000, 1_400_000) and "te" in fuse.decide(e)["support"]
    d = fuse.decide(_long(ev(peak=1_500_000)))                                # one bin from the call: agrees
    assert d["call"] == (1_000_000, 1_400_000) and d["relation"] == "concordant" and "hic" in d["support"]
    m = {"sequence_call": dict(label="centromere_like", score=0.8)}
    d = fuse.decide(_long(ev(meth=m)))
    assert d["call"] == (1_000_000, 1_400_000) and "meth" in d["support"]


def test_satellite_bp_inside_is_no_support_for_a_long_unit_call():
    e = _long(ev())
    e.sat_bp = lambda a, b: 100_000
    assert fuse.decide(e)["call"] is None


def test_a_hic_peak_merely_near_a_long_unit_call_is_no_support():
    """Beyond one bin, within four: neither supports nor opposes (hic_near). The array is dropped; a strong validated
    peak then places the call, as on any chromosome without a sequence call, and a weak one does not."""
    d = fuse.decide(_long(ev(peak=2_000_000)))
    assert d["call"] == (2_000_000, 2_250_000) and d["relation"] == "functional_only"
    assert "long_unit_unsupported" in d["flags"]
    d = fuse.decide(_long(ev(peak=2_000_000, strong=False)))
    assert d["call"] is None and "long_unit_unsupported" in d["flags"]


def test_an_arm_axis_at_a_long_unit_call_is_no_support():
    e = _long(ev())
    e.hic = dict(bin_size=BIN, gate_pass=False, anch_strong=False, anch_level="none", anch_peak_start=-1, anch_peak_end=-1,
                 anch_core_start=-1, anch_core_end=-1, denovo_peak_start=-1, denovo_peak_end=-1, denovo_core_start=-1,
                 denovo_core_end=-1, denovo_flags="", gate_support=False, gate_genome_untested=False)
    e.fold = dict(status="pass", agrees=False, core_start=1_000_000, core_end=1_500_000, point=1_200_000)
    d = fuse.decide(e)
    assert d["call"] is None and "fold" not in d["support"]


def test_hic_alone_can_place_a_call_once_a_long_unit_array_is_dropped():
    e = _long(_alone(seq=(1_000_000, 1_400_000)))
    d = fuse.decide(e)
    assert d["call"] == (140_000, 170_000) and d["relation"] == "functional_only" and d["confidence"] == "low"
    assert d["flags"] == ["hic_denovo_only", "long_unit_unsupported"]


def test_a_long_unit_call_a_validated_peak_moves_is_moved_as_before():
    d = fuse.decide(_long(ev(peak=30_000_000, core=(29_500_000, 30_750_000))))
    assert d["call"] == (29_500_000, 30_750_000) and d["relation"] == "offset" and "long_unit_unsupported" not in d["flags"]


def test_a_satellite_call_with_no_support_stands():
    d = fuse.decide(ev(conf="low", sat=False))
    assert d["call"] == (1_000_000, 1_400_000) and d["flags"] == []


def test_a_dropped_array_is_a_long_unit_feature_in_the_loci():
    fused = dict(call=None, relation="none", seq_evidence="long_unit_tandem", seq=(5_000, 9_000), long_unit_dropped=True)
    mem = loci.members("c", 100_000, fused, pd.DataFrame(columns=["chrom"]), {}, None, None)
    assert [(m[2], m[3]) for m in mem] == [("long_unit_array", "satellite")]
    assert loci.kind(mem) == "long_unit_array"
    kept = dict(fused, long_unit_dropped=False)
    assert loci.kind(loci.members("c", 100_000, kept, pd.DataFrame(columns=["chrom"]), {}, None, None)) == "sequence_centromere"


def test_untested_prior_free_peak_elsewhere_only_flags_a_sequence_call():
    e = _alone(seq=(1_000_000, 1_400_000))
    e.hic.update(denovo_peak_start=30_000_000, denovo_peak_end=30_010_000, denovo_core_start=29_990_000,
                 denovo_core_end=30_030_000)
    d = fuse.decide(e)
    assert d["call"] == (1_000_000, 1_400_000) and d["confidence"] == "high" and d["flags"] == ["hic_denovo_elsewhere"]
    assert not d["against"]
    e.hic.update(gate_n=5, gate_genome_untested=False)                 # a test ran: no note
    assert "hic_denovo_elsewhere" not in fuse.decide(e)["flags"]
    near = _alone(seq=(140_000, 170_000))                              # the prior-free core reaches the call
    assert fuse.decide(near)["flags"] == []



def test_hic_alone_is_decided_for_the_genome():
    import pandas as pd
    from cenfuse.pipeline import hic_rows
    # three references: their own leave-one-out n is 2, but the gate ran on every other chromosome -> no Hi-C alone
    t = pd.DataFrame(dict(chrom=["a", "b", "c", "d"], gate_n=[2, 2, 2, 3], denovo_flags="", denovo_core_start=-1))
    assert not any(r["gate_genome_untested"] for r in hic_rows(t).values())
    e = _alone(gate_n=2)
    e.hic["gate_genome_untested"] = False
    assert fuse.decide(e)["call"] is None
    t = pd.DataFrame(dict(chrom=["a", "b"], gate_n=[0, 0], denovo_flags="", denovo_core_start=-1))
    assert all(r["gate_genome_untested"] for r in hic_rows(t).values())


def test_prior_free_core_beside_a_nor_is_not_clean():
    from cenfuse.pipeline import denovo_nor_flag
    h = dict(chrom="c", bin_size=BIN, denovo_core_start=1_000_000, denovo_core_end=1_250_000, denovo_flags="",
             masked=())
    assert denovo_nor_flag(h, [("c", 1_400_000, 1_600_000)])["denovo_flags"] == "rDNA"      # one bin away
    assert denovo_nor_flag(h, [("c", 5_000_000, 5_200_000)])["denovo_flags"] == ""
    assert denovo_nor_flag(h, [("other", 1_000_000, 1_200_000)])["denovo_flags"] == ""
    e = _alone()
    e.hic.update(denovo_flags="rDNA")
    assert fuse.decide(e)["call"] is None


def test_support_tier_yields_to_a_passed_fold():
    folded = dict(status="pass", agrees=True, core_start=1_000_000, core_end=1_500_000, point=1_250_000)
    e = _weak(conf="low")
    e.hic.update(denovo_peak_start=1_200_000, denovo_peak_end=1_450_000, denovo_core_start=1_200_000,
                 denovo_core_end=1_450_000)
    e.fold = folded
    d = fuse.decide(e)
    assert "hic_support" not in d["flags"] and d["call"] == (1_000_000, 1_400_000)


def test_hic_does_not_raise_a_call_it_cannot_tell_from_its_competing_cluster():
    """0.10.1: a low call with a competing cluster beside it (no mappable bin between the two: the mouse major- and
    minor-satellite arrays) is not raised by a peak that agrees with both; a competitor Hi-C can separate still is."""
    near = ev(conf="low", peak=1_500_000)
    near.seq.update(alt_start=1_400_000, alt_end=1_900_000)            # abuts the call: no mappable bin between
    d = fuse.decide(near)
    assert d["confidence"] == "low" and "hic_alt_ambiguous" in d["flags"] and "hic" not in d["support"]
    blind = ev(conf="low", peak=1_500_000)                             # the competitor across a blind stretch
    blind.seq.update(alt_start=5_000_000, alt_end=5_500_000)
    blind.hic["masked"] = ((1_400_000, 5_000_000),)                    # everything between is unmappable
    assert "hic_alt_ambiguous" in fuse.decide(blind)["flags"]
    far = ev(conf="low", peak=1_500_000)
    far.seq.update(alt_start=20_000_000, alt_end=20_500_000)           # far: Hi-C tells them apart
    d = fuse.decide(far)
    assert d["confidence"] == "medium" and "hic_alt_ambiguous" not in d["flags"] and "hic" in d["support"]
    w = _weak(conf="low", level="medium")                              # the support tier, likewise
    w.seq.update(alt_start=1_400_000, alt_end=1_900_000)
    d = fuse.decide(w)
    assert d["confidence"] == "low" and "hic_support" not in d["flags"] and "hic_alt_ambiguous" in d["flags"]


# ---------------------------------------------------------------- LTR cores and Hi-C (hicfuse.anchored_on_cores)
def _paired(conf, core, status="support", gate=False):
    e = ev(conf=conf, peak=30_000_000, gate=gate, strong=False, level="low")
    e.hic = dict(e.hic, gate_support=False, te_pair=dict(core=core, status=status, clade="CRM"))
    return e


def test_a_pair_at_the_call_supports_it():
    d = fuse.decide(_paired("low", (1_100_000, 1_200_000)))
    assert d["call"] == (1_000_000, 1_400_000) and d["relation"] == "concordant" and d["confidence"] == "medium"
    assert "te_hic_pair" in d["flags"] and {"te", "hic"} <= set(d["support"])


def test_a_pair_elsewhere_moves_a_low_call():
    d = fuse.decide(_paired("low", (30_000_000, 31_000_000)))
    assert d["call"] == (30_000_000, 31_000_000) and d["relation"] == "offset" and d["source"] == "hic"
    assert d["confidence"] == "medium" and "te_hic_pair" in d["flags"] and {"te", "hic"} <= set(d["support"])


def test_a_pair_on_a_minority_leaves_a_confident_call():
    for conf in ("medium", "high"):
        d = fuse.decide(_paired(conf, (30_000_000, 31_000_000)))
        assert d["call"] == (1_000_000, 1_400_000) and d["confidence"] == conf and "te_hic_pair" not in d["flags"]


def test_a_pair_on_most_chromosomes_moves_any_call():
    d = fuse.decide(_paired("medium", (30_000_000, 31_000_000), "pass"))
    assert d["call"] == (30_000_000, 31_000_000) and d["confidence"] == "medium"
    d = fuse.decide(_paired("high", (30_000_000, 31_000_000), "pass"))
    assert d["call"] == (30_000_000, 31_000_000) and d["confidence"] == "high" and "moved_with_te" in d["flags"]


def test_a_pair_cannot_raise_a_call_it_cannot_tell_from_its_competitor():
    e = _paired("medium", (1_100_000, 1_200_000), "pass")
    e.seq = dict(e.seq, alt_start=1_400_000, alt_end=1_500_000)                # no mappable bin between the call and the competitor
    d = fuse.decide(e)
    assert d["call"] == (1_000_000, 1_400_000) and d["confidence"] == "medium" and "te_hic_pair" not in d["flags"]
    assert "hic" not in d["support"] and "hic_alt_ambiguous" in d["flags"]


def test_the_ltr_lines_leave_a_chromosome_whose_own_gate_passed():
    """They run where the gate fails on most chromosomes, but a chromosome can pass it alone."""
    e = ev(conf="medium", peak=1_100_000, gate=True, strong=True, level="high")
    base = fuse.decide(e)
    e.hic = dict(e.hic, te_pair=dict(core=(30_000_000, 31_000_000), status="pass", clade="CRM"),
                 cis=dict(peak=(30_000_000, 30_250_000), status="pass"))
    d = fuse.decide(e)
    assert base["relation"] == "concordant" and (d["call"], d["relation"], d["confidence"]) == (base["call"], base["relation"], base["confidence"])
    assert not {"te_hic_pair", "cis_elsewhere", "cis_agrees"} & set(d["flags"])


def test_a_pair_places_a_call_where_there_is_none():
    e = _paired("low", (30_000_000, 31_000_000))
    e.seq = dict(e.seq, start=-1, end=-1, evidence="none", confidence="none")
    d = fuse.decide(e)
    assert d["call"] == (30_000_000, 31_000_000) and d["relation"] == "functional_only" and d["confidence"] == "medium"
