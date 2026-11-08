"""Trans arm alignment (fold.align_scan, 0.10.1): the Rabl axis read from inter-chromosomal contacts, its validation
against the prior-free Hi-C peaks, and its place beside the arm fold in the fusion rules."""
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp

from cenfuse import fold, fuse
from cenfuse.hic.balance import balanced, ice
from cenfuse.hic.bins import Bins
from cenfuse.hic.score import block_means


SIZE = 10_000
LENGTHS = (23, 81, 32, 90, 58, 27, 85, 56, 44, 75, 67, 88, 62, 78, 70, 60)
CENS = (15, 24, 11, 45, 15, 15, 50, 10, 35, 43, 44, 15, 26, 63, 32, 55)


def rabl_counts(lengths=LENGTHS, cens=CENS, align=True, depth=60.0, seed=1):
    """Poisson counts of a Rabl-like nucleus: cis decays with distance; trans = 1 + 3 exp(-|d_i - d_j| / 4) with d
    = bins from the own centromere (loci at equal distances from their centromeres meet: the trans 'X'); align=False:
    flat trans (no axis). Returns (M, Bins)."""
    b = Bins({f"chr{k + 1}": L * SIZE for k, L in enumerate(lengths)}, SIZE)
    d = np.abs(b.local - np.array(cens)[b.chrom_of]).astype(float)
    same = b.chrom_of[:, None] == b.chrom_of[None, :]
    trans = 1 + (3 * np.exp(-np.abs(d[:, None] - d[None, :]) / 4) if align else 0)
    cis = 20 * np.exp(-np.abs(b.local[:, None] - b.local[None, :]) / 3)
    lam = depth * np.where(same, cis, trans)
    U = np.triu(np.random.default_rng(seed).poisson(lam)).astype(float)
    return sp.csr_matrix(U + U.T - np.diag(np.diag(U))), b


def _scan(align, centers=None):
    M, b = rabl_counts(align=align)
    bias, good, _, _ = ice(M)
    W = balanced(M, bias, good)
    members = b.nb >= 20
    centers = {k: c for k, c in enumerate(CENS)} if centers is None else centers
    calls = fold.align_scan(W, b, good, block_means(W, b, good, members), members, centers)
    table = pd.DataFrame([dict(chrom=n, callable=bool(members[k]), denovo_confidence="high",
                               denovo_peak_start=CENS[k] * SIZE, denovo_peak_end=(CENS[k] + 1) * SIZE)
                          for k, n in enumerate(b.names)])
    return calls, table, b


@pytest.fixture(scope="module")
def rabl():
    return _scan(True)


def test_alignment_core_holds_the_centromere_in_a_rabl_nucleus(rabl):
    calls, table, b = rabl
    hits = [calls[n]["core_start"] <= CENS[k] * SIZE < calls[n]["core_end"] for k, n in enumerate(b.names)
            if n in calls]
    assert len(hits) >= 14 and sum(hits) >= 0.8 * len(hits)
    v = fold.validate(calls, table, label="trans arm alignment")
    assert v["status"] == "pass" and v["p"] <= fold.GATE_P / 2


def test_alignment_fails_without_3d_structure():
    """Flat trans contacts (no Rabl axis): the cores hold the peaks at chance; the test fails."""
    calls, table, _ = _scan(False)
    assert fold.validate(calls, table, label="trans arm alignment")["status"] != "pass"


def test_alignment_needs_two_centered_partners(rabl):
    _, _, b = rabl
    out = fold.align_scan(sp.csr_matrix((b.N, b.N)), b, np.ones(b.N, bool), np.ones((len(b.names), len(b.names))),
                          b.nb >= 20, {0: CENS[0]})
    assert out == {}


def test_telocentric_centers_are_not_scanned():
    """A chromosome whose centromere sits at its end: the alignment center is only searched inside 12.5-87.5%."""
    cens = (2,) + CENS[1:]
    M, b = rabl_counts(cens=cens)
    bias, good, _, _ = ice(M)
    W = balanced(M, bias, good)
    members = b.nb >= 20
    calls = fold.align_scan(W, b, good, block_means(W, b, good, members), members, dict(enumerate(cens)))
    c = calls.get("chr1")
    L = LENGTHS[0] * SIZE
    assert c is None or fold.CI_MIN * L <= c["point"] <= (1 - fold.CI_MIN) * L


def test_validation_uses_the_bonferroni_level_and_the_exact_chance():
    """Each Rabl-axis readout passes at GATE_P / 2. A peak agrees when it overlaps or touches the core, so a peak bin
    agrees from one bin before the core to one after it: chance = (core + 2 peak) / scanned = 0.15 here, and 4/6
    agreements (P 0.0059) pass at 0.01 but not at 0.005."""
    folds = {f"c{i}": dict(point=50, core_start=0, core_end=130, scanned=1000) for i in range(6)}
    t = pd.DataFrame([dict(chrom=f"c{i}", callable=True, denovo_confidence="high", denovo_peak_start=p,
                           denovo_peak_end=p + 10) for i, p in enumerate([50] * 4 + [900] * 2)])
    v = fold.validate(folds, t)
    assert v["k"] == 4 and fold.GATE_P / 2 < v["p"] <= fold.GATE_P and v["status"] == "fail"
    assert abs(v["p"] - 0.00589) < 2e-4
    t.loc[0, "denovo_peak_start"], t.loc[0, "denovo_peak_end"] = 130, 140      # touching the core's end: agrees
    assert fold.validate(folds, t)["agrees"]["c0"]


# ---------------------------------------------------------------------------------------- fusion: the axis line
BIN = 250_000


def ev(fold_status="fail", align_status="pass", agrees=True, seq=(500_000, 1_000_000), conf="low", untested=False,
       fold_p=0.02):
    s = dict(start=seq[0] if seq else -1, end=seq[1] if seq else -1, evidence="satellite" if seq else "none",
             confidence=conf if seq else "none", flags=".", main_family="A", alt_start=-1, alt_end=-1)
    dn = 12_000_000
    h = dict(bin_size=BIN, gate_pass=False, gate_n=0 if untested else 5, gate_genome_untested=untested,
             anch_strong=False, anch_level="none", anch_flags="", anch_peak_start=-1, anch_peak_end=-1,
             anch_core_start=-1, anch_core_end=-1, masked=(), denovo_confidence="high", denovo_flags="",
             denovo_peak_start=dn, denovo_peak_end=dn + BIN, denovo_core_start=dn - BIN, denovo_core_end=dn + 2 * BIN)
    core = dict(point=12_500_000, core_start=11_000_000, core_end=14_000_000, scanned=18_000_000)
    f = dict(core, status=fold_status, agrees=False)
    a = dict(core, status=align_status, agrees=agrees, fold_status=fold_status, fold_p=fold_p)
    return fuse.Evidence("c", 25_000_000, s, None, h, None, lambda x, y: 0, f, a)


def test_alignment_acts_only_where_the_fold_leans_the_same_way():
    """The alignment depends on the prior-free partner centers (a holocentric's symmetric compartments mimic it); it
    stands in for a failed fold only where the fold's own test leans the same way (P <= AXIS_LEAN_P)."""
    d = fuse.decide(ev(fold_p=0.5))
    assert d["call"] == (500_000, 1_000_000) and not [f for f in d["flags"] if "align" in f]
    assert "align" not in d["support"] + d["against"]
    d = fuse.decide(ev(fold_p=0.05))                        # with or without gate references
    assert "align_pair" in d["flags"] and d["confidence"] == "medium"


def test_alignment_takes_the_folds_place_where_the_fold_fails():
    d = fuse.decide(ev(untested=True))
    assert d["call"] == (11_750_000, 12_500_000) and d["relation"] == "offset" and d["confidence"] == "medium"
    assert "align_pair" in d["flags"] and {"hic", "align"} <= set(d["support"])
    row = fuse.row(ev(untested=True), d)
    assert row["align_gate"] == "pass" and row["align_core_start"] == 11_000_001 and row["fold_gate"] == "fail"
    d = fuse.decide(ev(conf="high", untested=True))
    assert "moved_with_align" in d["flags"]


def test_a_passed_fold_decides_over_the_alignment():
    d = fuse.decide(ev(fold_status="pass", align_status="pass", untested=True))     # fold: no pair here; opposes
    assert "align_pair" not in d["flags"] and "fold_offset" in d["flags"] and d["relation"] == "contested"


def test_a_passed_fold_is_the_genomes_axis_line_also_where_it_makes_no_call():
    e = ev(untested=True)
    e.fold = None                                   # no fold call on this chromosome, but the genome's fold passed
    e.align = dict(e.align, fold_status="pass")
    d = fuse.decide(e)
    assert not [f for f in d["flags"] if "align" in f] and "align" not in d["support"]


def test_no_axis_line_without_a_passed_test():
    d = fuse.decide(ev(align_status="fail", untested=True))
    assert d["call"] == (500_000, 1_000_000) and not {f for f in d["flags"] if "align" in f or "fold" in f}


def test_hic_alone_medium_where_the_alignment_pairs_low_where_not():
    d = fuse.decide(ev(seq=None, untested=True))
    assert d["flags"] == ["align_pair"] and d["confidence"] == "medium" and d["relation"] == "functional_only"
    d = fuse.decide(ev(seq=None, untested=True, agrees=False))
    assert d["flags"] == ["hic_denovo_only"] and d["confidence"] == "low"


def test_calls_at_the_prior_free_core_take_the_provisional_methylation_item():
    """0.10.1's future.md entry 17 (closed): a call placed at the prior-free core (axis pair) reads methylation at that core."""
    e = ev(seq=None, untested=True)
    e.meth = {"provisional": dict(label="centromere_like", score=0.9), "hic_peak": dict(label="dip", score=0.1)}
    d = fuse.decide(e)
    assert d["meth_call"]["label"] == "centromere_like" and "meth" in d["support"]


def test_architecture_label():
    from cenfuse.output import architecture
    base = dict(hic_chromosomes=10, hic_gate_untested=0, hic_gate_pass=0, hic_gate_support=0, fold_gate="fail",
                align_gate="fail")
    assert architecture(base) == "not_detected"
    assert architecture(dict(base, hic_gate_support=10)) == "weak_centromere_clustering"
    assert architecture(dict(base, hic_gate_pass=10)) == "centromere_clustering"
    assert architecture(dict(base, fold_gate="pass")) == "rabl_axis"
    assert architecture(dict(base, align_gate="pass", fold_p="0.5")) == "not_detected"   # fold does not lean
    assert architecture(dict(base, align_gate="pass", fold_p="0.013")) == "rabl_axis"
    alone = dict(base, hic_gate_untested=10)
    assert architecture(alone) == "untested" and architecture(dict(alone, align_gate="pass", fold_p="0.02")) == \
        "rabl_axis"


def cluster_counts(width=5, depth=60.0, seed=2):
    """Centromeres cluster over a flat +/- width-bin domain, arms do not align (no trans 'X')."""
    b = Bins({f"chr{k + 1}": L * SIZE for k, L in enumerate(LENGTHS)}, SIZE)
    d = np.abs(b.local - np.array(CENS)[b.chrom_of])
    same = b.chrom_of[:, None] == b.chrom_of[None, :]
    inside = d <= width
    trans = 1 + 3.0 * (inside[:, None] & inside[None, :])
    cis = 20 * np.exp(-np.abs(b.local[:, None] - b.local[None, :]) / 3)
    U = np.triu(np.random.default_rng(seed).poisson(depth * np.where(same, cis, trans))).astype(float)
    return sp.csr_matrix(U + U.T - np.diag(np.diag(U))), b


@pytest.mark.parametrize("width", [4, 5])
def test_a_centromere_cluster_without_aligned_arms_does_not_pass(width):
    """Code review (v10.1): with only the cluster (no arm alignment) the first form passed 15/16; the alignment now
    leaves out the cluster on both sides (partner windows and prior-free cores, k's own +/- HALF_WIN)."""
    M, b = cluster_counts(width)
    bias, good, _, _ = ice(M)
    W = balanced(M, bias, good)
    members = b.nb >= 20
    cores = {k: (c - width, c + width) for k, c in enumerate(CENS)}
    calls = fold.align_scan(W, b, good, block_means(W, b, good, members), members, dict(enumerate(CENS)), cores)
    table = pd.DataFrame([dict(chrom=n, callable=bool(members[k]), denovo_confidence="high",
                               denovo_peak_start=CENS[k] * SIZE, denovo_peak_end=(CENS[k] + 1) * SIZE)
                          for k, n in enumerate(b.names)])
    assert fold.validate(calls, table, label="trans arm alignment")["status"] != "pass"


def test_alignment_axis_uses_only_clean_prior_free_peaks():
    """Code review (v10.1): under the alignment, a flagged prior-free peak neither pairs nor opposes (the fold keeps
    0.10's rule)."""
    for flag in ("hub", "rDNA", "cycle_unstable", "not_converged"):
        e = ev()
        e.hic["denovo_flags"] = flag
        d = fuse.decide(e)
        assert "align_pair" not in d["flags"] and "hic_denovo_offset" not in d["flags"]
    e = ev(fold_status="pass", align_status="fail")
    e.fold = dict(e.fold, agrees=True)
    e.hic["denovo_flags"] = "rDNA"
    assert "fold_pair" in fuse.decide(e)["flags"]                  # barley chr2H (rDNA-flagged peak) moved right


def test_genome_axis_decision():
    assert fuse.genome_axis(None) is None
    assert fuse.genome_axis(dict(status="pass", p=1e-5, align=dict(status="fail"))) == "fold"
    assert fuse.genome_axis(dict(status="fail", p=0.013, align=dict(status="pass"))) == "align"
    assert fuse.genome_axis(dict(status="fail", p=0.06, align=dict(status="pass"))) is None
