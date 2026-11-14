"""0.13: the cis signal tested by Hi-C's own trans contacts (hicfuse.cis_trans_line), the support tier choosing the competing cluster,
the low calls joining too few confident references (hicfuse.gate, 'mixed'), and an AT trough paired with a prior-free Hi-C peak."""
import numpy as np
import pandas as pd
import scipy.sparse as sp

from cenfuse import fuse, hicfuse
from cenfuse.hic.bins import Bins
from cenfuse.hic.score import block_means

SIZE = 10_000


def _both_state(sign=+1, n=12, nb=200, trans=30.0, cis=3.0, cen=None, partners=None, other=None, seed=4):
    """Cis maps with distance decay (the contact crossing each centromere rises, sign +1, or falls, -1) and trans background, the
    centromere bins of the chromosomes in `partners` (default: all) in strong trans contact. `other`: {chrom index: bin} a second
    locus per chromosome with the opposite cis sign and its own trans cluster."""
    seqs = {f"c{i}": nb * SIZE for i in range(n)}
    cen = cen or {k: 40 + 9 * k for k in range(n)}
    b = Bins(seqs, SIZE)
    rng = np.random.default_rng(seed)
    d = sp.lil_matrix((b.N, b.N))
    for k in range(n):
        off = int(b.off[k])
        for i in range(nb):
            for j in range(i + 1, min(nb, i + 12)):
                x = 200.0 / (j - i) * rng.gamma(20, 1 / 20)
                if i < cen[k] <= j:
                    x *= cis if sign > 0 else 1 / cis
                if other and i < other[k] <= j:
                    x *= 1 / cis if sign > 0 else cis
                d[off + i, off + j] = d[off + j, off + i] = x
    dense = d.toarray()
    bg = rng.poisson(0.3, size=(b.N, b.N)).astype(float)
    bg = np.triu(bg, 1) + np.triu(bg, 1).T
    same = b.chrom_of[:, None] == b.chrom_of[None, :]
    dense = np.where(same, dense, bg)
    ps = list(range(n)) if partners is None else partners
    for locus in [cen] + ([other] if other else []):
        for x in ps:
            for y in ps:
                if x != y:
                    dense[b.off[x] + locus[x] - 1:b.off[x] + locus[x] + 2, b.off[y] + locus[y] - 1:b.off[y] + locus[y] + 2] += trans
    W = sp.csr_matrix(dense)
    good, members = np.ones(b.N, dtype=bool), np.ones(n, dtype=bool)
    return dict(W=W, bins=b, good=good, mean=block_means(W, b, good, members), members=members, bias=np.ones(b.N)), cen


def test_cis_extremes_that_trans_contacts_confirm_pass_with_their_sign():
    for sign in (+1, -1):
        st, cen = _both_state(sign)
        r = hicfuse.cis_trans_line(st)
        assert r["status"] == "pass" and r["sign"] == sign and r["tested_by"] == "trans" and r["interior"]
        assert set(r["peaks"]) == set(st["bins"].names) and all(r["pairs"].values())
        assert all(abs(r["peaks"][f"c{k}"][0] // SIZE - c) <= 1 for k, c in cen.items())
        assert hicfuse.cis_acts(r)


def test_without_trans_clustering_the_cis_signal_alone_does_not_pass():
    st, _ = _both_state(trans=0.0)
    r = hicfuse.cis_trans_line(st)
    assert r["status"] == "fail" and not r["peaks"] and not hicfuse.cis_acts(r)


def test_extremes_at_the_chromosome_ends_never_anchor():
    """The cis signal and the trans cluster both at the ends of every chromosome (ends fold and cluster whatever the centromere: Ae.
    tauschii, potato): unrestricted extremes would pass; read outside the terminal eighths they do not."""
    nb = 200
    ends = {k: 5 for k in range(12)}
    st, _ = _both_state(cen=ends, nb=nb)
    unrestricted = {}
    for k, c in enumerate(st["bins"].names):
        v = hicfuse.cis_profile(st["W"].tocsr(), st["bins"], st["good"], k)
        j = int(np.nanargmax(v))
        unrestricted[c] = (j * SIZE, (j + 1) * SIZE)
    assert hicfuse.anchored_on_cores(st, unrestricted)["status"] == "pass"
    r = hicfuse.cis_trans_line(st)
    assert r["status"] == "fail" and not r["peaks"]


def test_support_alone_does_not_act():
    """Trans contacts confirm the cis extremes on a minority only (4 of 12): beyond chance, but agreement on a minority can come from
    the coverage the two readouts share (randomized soybean Hi-C), so nothing acts."""
    st, _ = _both_state(partners=[0, 1, 2, 3])
    r = hicfuse.cis_trans_line(st)
    assert r["status"] == "support" and r["k"] == 4 and not r["peaks"] and not hicfuse.cis_acts(r)


def test_both_signs_confirmed_leave_the_sign_unknown():
    st, _ = _both_state(other={k: 150 - 3 * k for k in range(12)})
    r = hicfuse.cis_trans_line(st)
    assert r["status"] == "ambiguous" and not r["peaks"] and not hicfuse.cis_acts(r)


def test_the_other_sign_at_support_only_does_not_silence_a_passing_sign():
    """Support can come from coverage alone (randomized soybean), so it never vetoes a majority on the other sign."""
    st, cen = _both_state(other={k: 150 - 3 * k for k in range(12)})
    other = {k: 150 - 3 * k for k in range(12)}
    b = st["bins"]
    W = st["W"].toarray()
    for x in range(12):                       # the other sign's loci cluster in trans on 4 chromosomes only
        for y in range(12):
            if x != y and not (x < 4 and y < 4):
                W[b.off[x] + other[x] - 1:b.off[x] + other[x] + 2, b.off[y] + other[y] - 1:b.off[y] + other[y] + 2] -= 30.0
    st = dict(st, W=sp.csr_matrix(W))
    st["mean"] = block_means(st["W"], b, st["good"], st["members"])
    r = hicfuse.cis_trans_line(st)
    assert r["status"] == "pass" and r["sign"] == +1 and len(r["peaks"]) == 12


def test_a_chromosome_whose_others_fail_is_left_out():
    """Leave one out, as the gate: 7 of 12 agree, so each agreeing chromosome sees 6 of 11 among the others (still a majority),
    while with 6 of 12 agreeing the others of an agreeing chromosome hold 5 of 11 (no majority) and it is left out."""
    st, _ = _both_state(partners=list(range(7)))
    r = hicfuse.cis_trans_line(st)
    assert r["status"] == "pass" and r["k"] == 7 and set(r["peaks"]) == set(st["bins"].names)
    st, _ = _both_state(partners=list(range(6)))
    r = hicfuse.cis_trans_line(st)
    assert r["k"] == 6 and r["status"] != "pass" and not r["peaks"]


def test_the_cis_entry_carries_what_tested_it():
    from cenfuse.pipeline import cis_entry
    trans = dict(status="pass", tested_by="trans", interior=True, peaks={"a": (1, 2)}, pairs={"a": True})
    assert cis_entry(trans, "a") == dict(peak=(1, 2), status="pass", tested_by="trans", interior=True, pair=True,
                                         trans_pair=False)
    assert cis_entry(trans, "b") is None
    assert cis_entry(dict(trans, status="support"), "a") is None
    old = dict(status="support", peaks={"a": (1, 2)})                         # a 0.12 result: the LTR-tested line
    assert cis_entry(old, "a") == dict(peak=(1, 2), status="support", tested_by="ltr", interior=False, pair=False,
                                       trans_pair=False)


def test_ltr_tested_cis_acts_at_support_and_trans_tested_only_at_pass():
    assert hicfuse.cis_acts(dict(status="support", tested_by="ltr"))
    assert hicfuse.cis_acts(dict(status="support"))                       # 0.12 results carry no tested_by
    assert not hicfuse.cis_acts(dict(status="support", tested_by="trans"))
    assert hicfuse.cis_acts(dict(status="pass", tested_by="trans"))
    assert not hicfuse.cis_acts(None)


def _cis_ev(conf, peak, seq, interior=True, status="pass"):
    h = dict(callable=True, bin_size=SIZE, gate_pass=False, gate_support=False, anch_strong=False, anch_level="none", anch_peak_start=-1,
             anch_peak_end=-1, anch_core_start=-1, anch_core_end=-1, masked=(),
             cis=dict(peak=peak, status=status, tested_by="trans" if interior else "ltr", interior=interior))
    s = dict(start=seq[0], end=seq[1], evidence="satellite", confidence=conf, flags=".", main_family="A") if seq else \
        dict(start=-1, end=-1, evidence="none", confidence="none", flags=".", main_family=".")
    return fuse.decide(fuse.Evidence("c", 50_000_000, s, None, h, None, lambda a, b: 0))


def test_a_line_read_outside_the_ends_has_no_say_on_a_call_at_an_end():
    end_call, inner = (1_000_000, 1_400_000), (30_000_000, 30_010_000)       # 1.2 Mb into a 50-Mb chromosome: a terminal eighth
    for conf in ("low", "medium", "high"):
        d = _cis_ev(conf, inner, end_call)
        assert d["call"] == end_call and d["confidence"] == conf and d["relation"] == "sequence_only"
        assert not {"cis_only", "cis_elsewhere", "cis_agrees"} & set(d["flags"]) and "hic" not in d["against"]
        assert "cis_blind" in d["flags"]                                    # said in the caveat
    d = _cis_ev("low", (1_100_000, 1_110_000), end_call)                    # the extreme beside the call: nothing to say
    assert "cis_blind" not in d["flags"]
    d = _cis_ev("low", inner, end_call, interior=False)                      # the LTR-tested line reads every bin
    assert d["call"] == inner and "cis_only" in d["flags"]


def test_the_trans_tested_line_moves_an_interior_low_call_and_places_one():
    d = _cis_ev("low", (30_000_000, 30_010_000), (12_000_000, 12_400_000))
    assert d["call"] == (30_000_000, 30_010_000) and d["confidence"] == "low" and "cis_only" in d["flags"]
    d = _cis_ev("low", (12_100_000, 12_110_000), (12_000_000, 12_400_000))
    assert d["relation"] == "concordant" and d["confidence"] == "medium" and "cis_agrees" in d["flags"]
    d = _cis_ev(None, (30_000_000, 30_010_000), None)
    assert d["call"] == (30_000_000, 30_010_000) and d["relation"] == "functional_only" and d["confidence"] == "low"


# ---------------------------------------------------------------- the support tier chooses the competing cluster
BIN = 250_000


def _alt_ev(gate_pass=False, support=True, strong=True, level="high", flags="", peak=20_000_000, conf="low", refs="low"):
    s = dict(start=1_000_000, end=1_400_000, evidence="satellite", confidence=conf, flags="competing_cluster", main_family="satA",
             alt_start=20_000_000, alt_end=20_300_000)
    h = dict(callable=True, bin_size=BIN, gate_pass=gate_pass, gate_support=support, gate_refs=refs, anch_strong=strong,
             anch_level=level, anch_flags=flags, anch_peak_start=peak, anch_peak_end=peak + BIN, anch_core_start=peak,
             anch_core_end=peak + BIN, masked=())
    return fuse.decide(fuse.Evidence("c", 50_000_000, s, None, h, None, lambda a, b: 300_000))


def test_the_support_tier_chooses_the_competing_cluster_with_a_strong_peak():
    d = _alt_ev()
    assert d["call"] == (20_000_000, 20_300_000) and "hic_chose_alt" in d["flags"] and "hic_support" in d["flags"]
    assert d["relation"] == "concordant" and d["confidence"] == "medium" and "hic" in d["support"]


def test_the_support_tier_needs_a_strong_clean_peak_to_choose():
    for kw in (dict(strong=False, level="medium"), dict(flags="telomeric"), dict(support=False)):
        d = _alt_ev(**kw)                     # a telomeric peak never agrees under low references (clean peaks only)
        assert d["call"] == (1_000_000, 1_400_000) and "hic_chose_alt" not in d["flags"]


def test_tested_ltr_or_cis_readings_decide_instead_of_the_support_tier_choice():
    """Hi-C read on the LTR cores at the original call (te_pair) and a strong support-tier peak at the competing cluster: the
    readings tested against the LTR cores decide; the choice would otherwise be undone by them, leaving self-contradictory flags."""
    s = dict(start=12_000_000, end=12_400_000, evidence="satellite", confidence="low", flags="competing_cluster", main_family="satA",
             alt_start=20_000_000, alt_end=20_300_000)
    h = dict(callable=True, bin_size=BIN, gate_pass=False, gate_support=True, gate_refs="low", anch_strong=True, anch_level="high",
             anch_flags="", anch_peak_start=20_000_000, anch_peak_end=20_250_000, anch_core_start=20_000_000, anch_core_end=20_250_000,
             masked=(), te_pair=dict(core=(12_100_000, 12_300_000), status="pass", clade="CRM"))
    d = fuse.decide(fuse.Evidence("c", 50_000_000, s, None, h, None, lambda a, b: 300_000))
    assert "hic_chose_alt" not in d["flags"] and "te_hic_pair" in d["flags"]
    assert d["call"] == (12_000_000, 12_400_000) and d["relation"] == "concordant"
    h2 = dict(h, te_pair=None, cis=dict(peak=(30_000_000, 30_010_000), status="pass", tested_by="trans", interior=True))
    d = fuse.decide(fuse.Evidence("c", 50_000_000, s, None, h2, None, lambda a, b: 300_000))
    assert "hic_chose_alt" not in d["flags"] and "cis_only" in d["flags"] and d["call"] == (30_000_000, 30_010_000)


def test_a_passed_gate_still_chooses_as_before():
    d = _alt_ev(gate_pass=True, support=False, refs="reference")
    assert d["call"] == (20_000_000, 20_300_000) and "hic_chose_alt" in d["flags"] and "hic_support" not in d["flags"]
    assert d["relation"] == "concordant" and d["confidence"] == "medium"


# ---------------------------------------------------------------- too few confident references (future.md #30)
def _gate_table(agree, conf, flags=None):
    n = len(agree)
    return pd.DataFrame(dict(chrom=[f"c{i}" for i in range(n)], anch_peak_start=[0] * n, chance=[0.05] * n, seq_start=[0] * n,
                             seq_confidence=conf, anch_agrees_seq=agree, anch_flags=flags or [""] * n))


def test_low_calls_join_too_few_confident_references():
    t = hicfuse.gate(_gate_table([True] * 8, ["high", "high"] + ["low"] * 6))
    assert t.gate_pass.all() and (t.gate_refs == "mixed").all() and (t.gate_n == 7).all()
    t = hicfuse.gate(_gate_table([True] * 8, ["high", "high"] + ["low"] * 6, ["", ""] + ["telomeric"] * 5 + [""]))
    assert not t.gate_pass.any()                     # a low reference agrees only at a clean peak
    t = hicfuse.gate(_gate_table([True] * 8, ["high"] * 2 + ["medium"] * 2 + ["low"] * 4))
    assert (t.gate_refs == "reference").all() and (t.gate_n.isin([3, 4])).all()   # enough confident references: as before


def test_with_mixed_references_only_a_clean_peak_supports_a_low_call():
    def one(conf, flags):
        s = dict(start=1_000_000, end=1_400_000, evidence="satellite", confidence=conf, flags=".", main_family="satA")
        h = dict(callable=True, bin_size=BIN, gate_pass=True, gate_refs="mixed", anch_strong=False, anch_level="medium",
                 anch_flags=flags, anch_peak_start=1_200_000, anch_peak_end=1_450_000, anch_core_start=1_200_000,
                 anch_core_end=1_450_000, masked=())
        return fuse.decide(fuse.Evidence("c", 50_000_000, s, None, h, None, lambda a, b: 300_000))
    assert one("low", "")["confidence"] == "medium"
    assert one("low", "telomeric")["confidence"] == "low"
    assert one("high", "telomeric")["relation"] == "concordant"


# ---------------------------------------------------------------- an AT trough and a prior-free Hi-C peak (support)
def _at_ev(support=True, agrees=True, dflags="", dconf="medium", validated=False, status="fail"):
    s = dict(start=-1, end=-1, evidence="none", confidence="none", flags=".", main_family=".")
    h = dict(callable=True, bin_size=SIZE, gate_pass=False, gate_support=False, gate_genome_untested=True, anch_strong=False,
             anch_level="none", anch_peak_start=-1, anch_peak_end=-1, anch_core_start=-1, anch_core_end=-1, masked=(),
             denovo_peak_start=500_000, denovo_peak_end=510_000, denovo_core_start=500_000, denovo_core_end=510_000,
             denovo_confidence=dconf, denovo_flags=dflags)
    at = dict(start=495_000, end=520_000, gc=0.43, status=status, window=10_000, validated=validated, support=support, agrees=agrees)
    return fuse.decide(fuse.Evidence("c", 2_000_000, s, None, h, None, lambda a, b: 0, at=at))


def test_a_trough_and_a_prior_free_peak_that_agree_place_a_low_call():
    d = _at_ev()
    assert d["call"] == (495_000, 520_000) and d["confidence"] == "low" and "at_denovo_pair" in d["flags"]
    assert d["source"] == "composition" and set(d["support"]) == {"at", "hic"}


def test_no_pair_without_agreement_support_or_a_clean_peak():
    for kw in (dict(agrees=False), dict(support=False), dict(dflags="telomeric"), dict(dconf="none")):
        assert _at_ev(**kw)["call"] is None


# ---------------------------------------------------------------- 0.16: the cis extreme tested twice; the own trans pair (#27, #32)
def _cis2(conf, peak, seq, **cis):
    c = dict(peak=peak, status="pass", tested_by="ltr", interior=False)
    c.update(cis)
    h = dict(callable=True, bin_size=SIZE, gate_pass=False, gate_support=False, anch_strong=False, anch_level="none", anch_peak_start=-1,
             anch_peak_end=-1, anch_core_start=-1, anch_core_end=-1, masked=(), cis=c)
    s = dict(start=seq[0], end=seq[1], evidence="retrotransposon", confidence=conf, flags=".", main_family=".")
    return fuse.decide(fuse.Evidence("c", 50_000_000, s, None, h, None, lambda a, b: 0))


def test_a_cis_extreme_tested_twice_and_its_own_trans_pair_moves_a_contradicted_call_low():
    """Tomato with tes/ and Hi-C: the cis signal passes against the LTR cores and by trans contacts; where the chromosome's own trans
    peak lies at its extreme, the confident call it contradicts moves there, low (a rule shaped on one genome)."""
    far, ext = (18_000_000, 19_000_000), (25_500_000, 25_510_000)
    for conf in ("medium", "high"):
        d = _cis2(conf, ext, far, trans_pair=True)
        assert d["call"] == ext and d["confidence"] == "low" and d["rule"] == "cis_trans_pair" and d["relation"] == "offset"
        assert "cis_trans_pair" in d["flags"] and "cis_elsewhere" not in d["flags"] and "hic" in d["support"]
        d = _cis2(conf, ext, far, trans_pair=False)                          # one test only: kept, contested, one level lower
        assert d["call"] == far and d["relation"] == "contested" and "cis_elsewhere" in d["flags"] and "cis_trans_pair" not in d["flags"]
    d = _cis2("medium", ext, far, trans_pair=True, status="support")         # the LTR-tested line on a minority: no move
    assert d["call"] == far and "cis_trans_pair" not in d["flags"]
    d = _cis2("medium", (18_500_000, 18_510_000), far, trans_pair=True)     # at the call: it agrees, as before
    assert d["call"] == far and "cis_agrees" in d["flags"]


def test_a_low_call_at_an_end_moves_to_the_interior_extreme_only_with_its_own_trans_pair():
    end_call, inner = (61_000_000, 61_300_000), (27_000_000, 27_010_000)    # tomato CM175073.1: 61 Mb into 69.5 Mb, a terminal eighth
    h = lambda pair: dict(callable=True, bin_size=SIZE, gate_pass=False, gate_support=False, anch_strong=False, anch_level="none",
                          anch_peak_start=-1, anch_peak_end=-1, anch_core_start=-1, anch_core_end=-1, masked=(),
                          cis=dict(peak=inner, status="pass", tested_by="trans", interior=True, pair=pair))
    s = lambda conf: dict(start=end_call[0], end=end_call[1], evidence="satellite", confidence=conf, flags=".", main_family="A")
    d = fuse.decide(fuse.Evidence("c", 69_500_000, s("low"), None, h(True), None, lambda a, b: 0))
    assert d["call"] == inner and d["confidence"] == "low" and d["rule"] == "cis_end_pair" and "cis_end_pair" in d["flags"]
    d = fuse.decide(fuse.Evidence("c", 69_500_000, s("low"), None, h(False), None, lambda a, b: 0))     # no pair: an acrocentric
    assert d["call"] == end_call and "cis_blind" in d["flags"]
    for conf in ("medium", "high"):                                         # a confident call at an end is left as it is
        d = fuse.decide(fuse.Evidence("c", 69_500_000, s(conf), None, h(True), None, lambda a, b: 0))
        assert d["call"] == end_call and d["confidence"] == conf and "cis_blind" in d["flags"]


def test_cis_entry_says_where_the_ltr_tested_extreme_is_its_own_trans_pair():
    from cenfuse.pipeline import cis_entry
    trans = dict(status="pass", peaks={"a": (10, 20), "b": (50, 60)}, pairs={"a": True, "b": True, "c": True})
    cis = dict(status="pass", tested_by="ltr", peaks={"a": (10, 20), "b": (30, 40), "c": (70, 80)}, trans=trans)
    assert cis_entry(cis, "a")["trans_pair"]                    # same extreme, acting there, a pair
    assert not cis_entry(cis, "b")["trans_pair"]                # the two tests read different extremes
    assert not cis_entry(cis, "c")["trans_pair"]                # the trans test does not act on c (left out of its peaks)
    assert not cis_entry(dict(cis, trans=dict(trans, status="support")), "a")["trans_pair"]
    assert not cis_entry(dict(cis, trans=None), "a")["trans_pair"]


# ---------------------------------------------------------------- 0.16: a degenerate Hi-C null is told apart (future.md #69)
def test_a_degenerate_null_is_flagged_apart_and_never_counts_as_tested():
    from cenfuse.hic.confidence import DEGENERATE, level
    from cenfuse.hic.qc import call_flags
    from cenfuse.hic.bins import Bins
    assert level(1.0, True, float("nan")) == "medium" and level(1.0, True, DEGENERATE) == "medium"   # same level either way
    assert level(1.0, False, DEGENERATE) == "low"
    b = Bins({"c": 1_000_000}, 10_000)
    pc = dict(peak=50, core_lo=49, core_hi=51)
    hub = np.zeros(b.N, dtype=bool)
    assert "null_untestable" in call_flags(0, pc, b, hub, {}, True, float("nan"), True)
    f = call_flags(0, pc, b, hub, {}, True, DEGENERATE, True)
    assert "null_degenerate" in f and "null_untestable" not in f
    assert not {"null_untestable", "null_degenerate"} & set(call_flags(0, pc, b, hub, {}, True, 0.01, True))
    assert fuse._clean_denovo(dict(denovo_confidence="medium", denovo_flags="null_untestable"))
    assert not fuse._clean_denovo(dict(denovo_confidence="medium", denovo_flags="null_degenerate"))


def test_a_hic_call_inside_a_decisive_core_does_not_count_the_core_against_it():
    """A Hi-C-placed call (cis_only, cis_trans_pair) that lies inside a decisive LTR core: the core is not a line against it,
    and the caveat does not say the core lies elsewhere (the same agreement rule as for sequence calls)."""
    from cenfuse import fuse, interpret
    core, ext, far = (25_400_000, 25_600_000), (25_500_000, 25_600_000), (18_000_000, 19_000_000)
    te = dict(clade="CRM", start=core[0], end=core[1], decisive=True, strong=True)
    for conf, extra in (("low", {}), ("medium", {"trans_pair": True})):
        cis = dict(peak=ext, status="pass", tested_by="ltr", interior=False, **extra)
        h = dict(callable=True, bin_size=100_000, gate_pass=False, gate_support=False, anch_strong=False, anch_level="none",
                 anch_peak_start=-1, anch_peak_end=-1, anch_core_start=-1, anch_core_end=-1, masked=(), cis=cis)
        s = dict(start=far[0], end=far[1], evidence="satellite", confidence=conf, flags=".", main_family="A")
        ev = fuse.Evidence("c", 50_000_000, s, te, h, None, lambda a, b: 0)
        d = fuse.decide(ev)
        assert d["call"] == ext and d["rule"] in ("cis_only", "cis_trans_pair")
        assert "te" not in d["against"]
        t = interpret.columns(fuse.table([ev], [d]), dict(unit={}, sat_bp=lambda c, a, b: 0, marker="CRM", decisive=True,
                                                         hic=True, tes=True, its=[], te_why=""))
        assert "lies elsewhere" not in t.iloc[0].caveat
