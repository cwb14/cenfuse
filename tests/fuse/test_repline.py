"""The repeat-density line (cenfuse.repline, 0.16) and its rule in fuse.decide. Synthetic FASTA only."""
import numpy as np
import pysam
import pytest

from cenfuse import fuse, repline

RNG = np.random.default_rng(5)
NT = np.array(list("ACGT"))
L = 400_000
N_CHROMS = 12


def rnd(n):
    return RNG.choice(NT, size=n)


def revcomp(s: str) -> str:
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def genome(tmp_path, planted=True, ends=False, arrays_at=None):
    """N_CHROMS chromosomes of random sequence; planted: on each, a cluster of 8 copies of one 3-kb element (shared by every
    chromosome: high copy number) at an interior place; ends: the same cluster at every chromosome's last 20 kb instead;
    arrays_at: {chrom: (start, end)} where a 50-copy tandem array of a 170-bp unit sits (the satellite line's)."""
    te = "".join(rnd(3_000))
    unit = "".join(rnd(170))
    f = tmp_path / "g.fa"
    with open(f, "w") as o:
        for i in range(N_CHROMS):
            c = f"c{i}"
            a = list(rnd(L))
            if planted or ends:
                at = L - 25_000 if ends else 100_000 + 15_000 * i
                for j in range(8):
                    cp = te if j % 2 else revcomp(te)
                    a[at + j * 3_000:at + (j + 1) * 3_000] = list(cp)
            if arrays_at and c in arrays_at:
                s = arrays_at[c][0]
                a[s:s + 50 * 170] = list(unit * 50)
            o.write(f">{c}\n{''.join(a)}\n")
    pysam.faidx(str(f))
    chroms = {f"c{i}": L for i in range(N_CHROMS)}
    return f, chroms


def test_canonical_kmers_read_both_strands_alike():
    s = "".join(rnd(200))
    code = lambda x: repline.LUT[np.frombuffer(x.encode(), dtype=np.uint8)]
    a, ok_a = repline.canonical_kmers(code(s))
    b, ok_b = repline.canonical_kmers(code(revcomp(s)))
    assert ok_a.all() and ok_b.all() and sorted(a.tolist()) == sorted(b.tolist())
    n = "".join(rnd(30)) + "N" + "".join(rnd(30))
    _, ok = repline.canonical_kmers(code(n))
    assert not ok[10:31].any() and ok[:10].all() and ok[31:].all()          # every K-mer holding the N is invalid


def test_planted_clusters_are_the_densest_bins(tmp_path):
    f, chroms = genome(tmp_path)
    b = repline.repeat_bins(f, chroms, chroms)
    x = b["c3"]
    s = (100_000 + 15_000 * 3) // 1000
    dens = x["l2"] / np.maximum(x["kpos"], 1)
    assert dens[s + 2] > 6 and np.median(dens) < 0.1                       # 96 copies: log2 ~6.6; unique sequence ~0
    assert (x["called"] == 1000).sum() == L // 1000 and (x["outside"] == x["called"]).all()


def test_one_cluster_per_chromosome_passes_and_places_it(tmp_path):
    f, chroms = genome(tmp_path)
    r = repline.scan(repline.repeat_bins(f, chroms, chroms), chroms)
    assert r["status"] == "pass" and r["p"] <= repline.REP_ALPHA
    for i in range(N_CHROMS):
        x = repline.on_chromosome(r, f"c{i}")
        at = 100_000 + 15_000 * i
        assert x["start"] < at + 24_000 and x["end"] > at and x["density"] > 3


def test_random_sequence_and_end_clusters_do_not_pass(tmp_path):
    f, chroms = genome(tmp_path, planted=False)
    assert repline.scan(repline.repeat_bins(f, chroms, chroms), chroms)["status"] != "pass"
    e = tmp_path / "e"
    e.mkdir()
    f, chroms = genome(e, planted=False, ends=True)
    r = repline.scan(repline.repeat_bins(f, chroms, chroms), chroms)
    assert r["status"] != "pass"                                             # the terminal eighths are not read


def test_satellite_arrays_are_left_out(tmp_path):
    """A tandem array is the satellite line's: its positions are not read (its copies still count genome-wide)."""
    arrays = {f"c{i}": (300_000, 300_000 + 50 * 170) for i in range(N_CHROMS)}
    f, chroms = genome(tmp_path, planted=False, arrays_at=arrays)
    with_mask = repline.repeat_bins(f, chroms, chroms, {c: [v] for c, v in arrays.items()})
    without = repline.repeat_bins(f, chroms, chroms)
    j = 300_000 // 1000 + 2
    assert without["c0"]["l2"][j] > 0 and with_mask["c0"]["kpos"][j] == 0 and with_mask["c0"]["outside"][j] == 0
    assert repline.scan(with_mask, chroms)["status"] != "pass"
    p = repline.profile(with_mask["c0"], 10, L)
    assert np.isinf(p[j - 2]) and np.isinf(p[0]) and np.isinf(p[-1])     # masked window; terminal eighths


def test_samples_large_genomes_with_exact_counts(tmp_path, monkeypatch):
    f, chroms = genome(tmp_path)
    monkeypatch.setattr(repline, "SAMPLE_ABOVE", 1_000)                     # force the hash sample
    b = repline.repeat_bins(f, chroms, chroms)
    total = sum(int(x["kpos"].sum()) for x in b.values())
    assert 0.02 * N_CHROMS * L < total < 0.05 * N_CHROMS * L                # about 1 in 32 positions kept
    r = repline.scan(b, chroms)
    assert r["status"] == "pass"


# ---------------------------------------------------------------- fuse
def _ev(seq=None, rep=None, at=None):
    s = dict(start=seq[0], end=seq[1], evidence="satellite", confidence="low", flags=".", main_family="A") if seq else \
        dict(start=-1, end=-1, evidence="none", confidence="none", flags=".", main_family=".")
    return fuse.decide(fuse.Evidence("c", 1_000_000, s, None, None, None, lambda a, b: 0, at=at, rep=rep))


REP = dict(start=400_000, end=450_000, density=6.5, status="pass", window=20_000)


def test_a_passing_repeat_line_places_only_what_nothing_else_places():
    d = _ev(rep=REP)
    assert d["call"] == (400_000, 450_000) and d["confidence"] == "low" and d["rule"] == "rep_only"
    assert "repeats" in d["support"] and d["source"] == "composition"
    assert _ev(rep=dict(REP, status="fail"))["call"] is None
    d = _ev(seq=(100_000, 200_000), rep=REP)                               # a sequence call stands
    assert d["call"] == (100_000, 200_000) and "rep_only" not in d["flags"]
    at = dict(start=600_000, end=603_000, gc=0.1, status="pass", window=1000, validated=False, support=False, agrees=False)
    d = _ev(rep=REP, at=at)                                                # the AT trough first, as before
    assert d["call"] == (600_000, 603_000) and "at_only" in d["flags"] and "rep_only" not in d["flags"]


def test_the_row_reports_the_repeat_line():
    e = fuse.Evidence("c", 1_000_000, dict(start=-1, end=-1, evidence="none", confidence="none", flags=".", main_family="."),
                      None, None, None, lambda a, b: 0, rep=REP)
    r = fuse.row(e, fuse.decide(e))
    assert (r["rep_gate"], r["rep_start"], r["rep_end"], r["rep_density"]) == ("pass", 400_001, 450_000, 6.5)
    assert set(fuse.FUSED_COLS) >= {"rep_gate", "rep_start", "rep_end", "rep_density"}


# ---------------------------------------------------------------- the repeat peak as second line for a Hi-C call (rep_pair)
HB = 10_000


def _hic_alone(peak):
    """Hi-C alone (a genome where nothing tests it) with a clean high prior-free peak at `peak`."""
    return dict(callable=True, bin_size=HB, gate_pass=False, gate_support=False, gate_genome_untested=True, anch_strong=False,
                anch_level="none", anch_peak_start=-1, anch_peak_end=-1, anch_core_start=-1, anch_core_end=-1, masked=(),
                denovo_peak_start=peak, denovo_peak_end=peak + HB, denovo_core_start=peak, denovo_core_end=peak + HB,
                denovo_confidence="high", denovo_flags="")


def _dec(h, rep):
    s = dict(start=-1, end=-1, evidence="none", confidence="none", flags=".", main_family=".")
    return fuse.decide(fuse.Evidence("c", 1_000_000, s, None, h, None, lambda a, b: 0, rep=rep))


def test_a_tested_repeat_peak_raises_the_hic_call_it_agrees_with():
    tested = dict(REP, validated=True, agrees=True)
    d = _dec(_hic_alone(420_000), tested)
    assert d["call"] == (400_000, 450_000) and d["confidence"] == "medium" and "rep_pair" in d["flags"]
    assert {"hic", "repeats"} <= set(d["support"]) and "hic_denovo_only" in d["flags"] and d["rule"] == "rep_pair"
    assert _dec(_hic_alone(420_000), dict(REP, support=True, agrees=True))["confidence"] == "medium"
    for rep in (dict(REP, agrees=True), dict(REP, validated=True, agrees=False)):            # untested; not by this call
        d = _dec(_hic_alone(420_000), rep)
        assert d["confidence"] == "low" and "rep_pair" not in d["flags"]
    d = _dec(_hic_alone(420_000), dict(tested, status="untested"))                          # the line passed no test
    assert d["confidence"] == "low" and "rep_pair" not in d["flags"] and d["call"] == (420_000, 430_000)
    d = _dec(_hic_alone(800_000), tested)                                                      # the peak lies elsewhere
    assert d["call"] == (800_000, 810_000) and d["confidence"] == "low" and "rep_pair" not in d["flags"]


def test_a_cis_call_pairs_with_a_repeat_peak_tested_against_the_cis_extremes():
    h = dict(callable=True, bin_size=HB, gate_pass=False, gate_support=False, gate_genome_untested=False, anch_strong=False,
             anch_level="none", anch_peak_start=-1, anch_peak_end=-1, anch_core_start=-1, anch_core_end=-1, masked=(),
             cis=dict(peak=(430_000, 440_000), status="pass", tested_by="trans", interior=True))
    d = _dec(h, dict(REP, cis_validated=True, cis_agrees=True))
    assert d["call"] == (400_000, 450_000) and d["confidence"] == "medium" and "rep_pair" in d["flags"] and "cis_only" in d["flags"]
    d = _dec(h, dict(REP, validated=True, agrees=False, cis_agrees=False))
    assert d["confidence"] == "low" and "rep_pair" not in d["flags"]


def test_with_hic_tests_the_peaks_at_the_repeat_test_window():
    import pandas as pd
    calls = {f"c{i}": dict(start=100_000 + i * 1000, end=120_000 + i * 1000, gc=6.0, first=True) for i in range(6)}
    res = dict(scales={20_000: dict(calls=calls, p=0.001)}, best=20_000, p=0.002, status="pass")
    lengths = {f"c{i}": 1_000_000 for i in range(6)}
    tab = pd.DataFrame([dict(chrom=c, callable=True, denovo_peak_start=v["start"] + 5000, denovo_peak_end=v["start"] + 15000,
                             bin_size=HB) for c, v in calls.items()])
    out = repline.with_hic(res, dict(table=tab), lengths)
    assert out["hic"]["status"] == "pass" and out["cis"] is None
    r = repline.on_chromosome(out, "c2")
    assert r["validated"] and r["agrees"] and r["start"] == 102_000
    assert repline.with_hic(res, None, lengths) is res
