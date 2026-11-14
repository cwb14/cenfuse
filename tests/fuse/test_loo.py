"""0.16.1: a tested line acts on a chromosome only at the status its test reaches without that chromosome (cenfuse.loo). Where a
test passes only thanks to chromosome c's own agreement, the line no longer raises, moves or places c's call; where the other
chromosomes pass it alone, it still does."""
import numpy as np
import pandas as pd
import scipy.sparse as sp

from cenfuse import atline, fold, fuse, hicfuse, loo, repline
from cenfuse.hic.bins import Bins
from cenfuse.hic.score import block_means
from cenfuse.pipeline import axis_entries, cis_entry, te_pair_entry

SIZE = 10_000


def test_the_stricter_status_and_the_genome_status_where_the_test_read_nothing():
    assert loo.stricter("pass", "support") == "support" and loo.stricter("support", "pass") == "support"
    assert loo.stricter("pass", "fail") == "fail" and loo.stricter("pass", "ambiguous") == "ambiguous"
    assert loo.stricter("fail", "untested") == "fail" and loo.stricter("support", "support") == "support"
    agree = {"a": True, "b": True, "c": False}
    by = loo.by_chromosome("pass", agree, lambda names: "pass" if all(agree[x] for x in names) else "fail")
    assert by == {"a": "fail", "b": "fail", "c": "pass"}
    res = dict(status="support", status_by_chrom={"a": "fail"})
    assert loo.status_on(res, "a") == "fail" and loo.status_on(res, "z") == "support" and loo.status_on(None, "a") is None
    statuses = ("pass", "support", "fail", "untested", "ambiguous")          # never above the genome: no control gains an action
    assert all(loo.RANK.get(loo.stricter(g, o), 0) <= loo.RANK.get(g, 0) for g in statuses for o in statuses)


# ---------------------------------------------------------------- Hi-C read on the LTR cores (te_hic_pair)
def _cores_state(n_right, n=12, nb=200, seed=3):
    """n chromosomes of 200 bins whose centromeres contact each other in trans; the LTR cores of the first n_right lie at them,
    the others 60 bins away (a zone of 10 bins: chance 0.05 per chromosome)."""
    seqs = {f"c{i}": nb * SIZE for i in range(n)}
    cen = {k: 20 + 10 * k for k in range(n)}
    b = Bins(seqs, SIZE)
    d = np.random.default_rng(seed).poisson(0.5, size=(b.N, b.N)).astype(float)
    d = np.triu(d) + np.triu(d, 1).T
    for x, cx in cen.items():
        for y, cy in cen.items():
            if x != y:
                d[b.off[x] + cx - 1:b.off[x] + cx + 2, b.off[y] + cy - 1:b.off[y] + cy + 2] += 200.0
    W = sp.csr_matrix(d)
    good, members = np.ones(b.N, dtype=bool), np.ones(n, dtype=bool)
    st = dict(W=W, bins=b, good=good, mean=block_means(W, b, good, members), members=members)
    off = lambda k: 0 if k < n_right else 60
    return st, {f"c{k}": ((cen[k] + off(k)) * SIZE, (cen[k] + off(k) + 2) * SIZE) for k in range(n)}


def _te_ev(entry, seq=(1_500_000, 1_700_000), conf="low"):
    h = dict(callable=True, bin_size=SIZE, gate_pass=False, gate_support=False, anch_strong=False, anch_level="none",
             anch_peak_start=-1, anch_peak_end=-1, anch_core_start=-1, anch_core_end=-1, masked=())
    if entry is not None:
        h["te_pair"] = entry
    s = dict(start=seq[0], end=seq[1], evidence="satellite", confidence=conf, flags=".", main_family="A") if seq else \
        dict(start=-1, end=-1, evidence="none", confidence="none", flags=".", main_family=".")
    return fuse.decide(fuse.Evidence("c0", 2_000_000, s, None, h, None, lambda a, b: 0))


def test_hic_at_the_ltr_cores_on_a_minority_does_not_act_where_its_own_agreement_made_the_test():
    """4 of 12 peaks at their cores: beyond chance (P 0.002, support), but without an agreeing chromosome 3 of 11 are not (P 0.015):
    the line no longer moves that chromosome's low call or places a call there (TS60 with tes/ and Hi-C). A chromosome that does not
    agree keeps the genome's status, which acts nowhere without its peak at its core."""
    st, cores = _cores_state(4)
    r = hicfuse.anchored_on_cores(st, cores)
    assert r["status"] == "support" and r["k"] == 4
    assert {c for c, s in r["status_by_chrom"].items() if s == "fail"} == {"c0", "c1", "c2", "c3"}
    assert {s for c, s in r["status_by_chrom"].items() if c not in ("c0", "c1", "c2", "c3")} == {"support"}
    sig_cores = {c: dict(start=a, end=b) for c, (a, b) in cores.items()}
    assert te_pair_entry(r, "c0", sig_cores, "CRM") is None and te_pair_entry(r, "c5", sig_cores, "CRM") is None
    d = _te_ev(te_pair_entry(r, "c0", sig_cores, "CRM"))
    assert d["call"] == (1_500_000, 1_700_000) and "te_hic_pair" not in d["flags"]
    assert _te_ev(te_pair_entry(r, "c0", sig_cores, "CRM"), seq=None)["call"] is None
    old = dict(core=cores["c0"], status=r["status"], clade="CRM")                    # 0.16.0: the genome's status acted here
    assert _te_ev(old)["rule"] == "te_hic_pair"


def test_hic_at_the_ltr_cores_that_the_other_chromosomes_pass_alone_still_acts():
    st, cores = _cores_state(12)
    r = hicfuse.anchored_on_cores(st, cores)
    assert r["status"] == "pass" and set(r["status_by_chrom"].values()) == {"pass"}
    sig_cores = {c: dict(start=a, end=b) for c, (a, b) in cores.items()}
    e = te_pair_entry(r, "c0", sig_cores, "CRM")
    assert e == dict(core=cores["c0"], status="pass", clade="CRM")
    d = _te_ev(e)
    assert d["call"] == cores["c0"] and d["rule"] == "te_hic_pair" and d["confidence"] == "medium"


def test_hic_at_the_ltr_cores_passes_on_a_majority_only_with_the_chromosome_itself():
    """3 of 5 at their cores pass (a majority, P 0.001); without an agreeing chromosome 2 of 4 hold no majority and too few
    agreements: the line no longer acts there (Arabidopsis Col-CC with tes/ and Hi-C: 3 of 5)."""
    st, cores = _cores_state(3, n=5)
    r = hicfuse.anchored_on_cores(st, cores)
    assert r["status"] == "pass" and r["k"] == 3
    assert [r["status_by_chrom"][f"c{i}"] for i in range(5)] == ["fail"] * 3 + ["pass"] * 2


# ---------------------------------------------------------------- the cis signal tested against the LTR cores
def _cis_state(n=12, nb=400, seed=4):
    """Cis maps with distance decay; the contact crossing each centromere (bin 30 + 25 k) rises threefold."""
    seqs = {f"c{i}": nb * SIZE for i in range(n)}
    cen = {k: 30 + 25 * k for k in range(n)}
    b = Bins(seqs, SIZE)
    rng = np.random.default_rng(seed)
    rows, cols, vals = [], [], []
    for k in range(n):
        off = int(b.off[k])
        for i in range(nb):
            for j in range(i + 1, min(nb, i + 12)):
                x = 200.0 / (j - i) * rng.gamma(20, 1 / 20)
                if i < cen[k] <= j:
                    x *= 3.0
                rows += [off + i, off + j]
                cols += [off + j, off + i]
                vals += [x, x]
    W = sp.csr_matrix((vals, (rows, cols)), shape=(b.N, b.N))
    return dict(W=W, bins=b, good=np.ones(b.N, dtype=bool), members=np.ones(n, dtype=bool), bias=np.ones(b.N)), cen


def _cis_ev(entry, call):
    h = dict(callable=True, bin_size=SIZE, gate_pass=False, gate_support=False, anch_strong=False, anch_level="none",
             anch_peak_start=-1, anch_peak_end=-1, anch_core_start=-1, anch_core_end=-1, masked=())
    if entry is not None:
        h["cis"] = entry
    s = dict(start=call[0], end=call[1], evidence="retrotransposon", confidence="low", flags=".", main_family=".")
    return fuse.decide(fuse.Evidence("c", 4_000_000, s, None, h, None, lambda a, b: 0))


def test_the_cis_line_on_a_minority_does_not_raise_a_call_where_its_own_agreement_made_the_test(monkeypatch):
    """3 of 12 cis peaks at the LTR cores: beyond chance (support), but the other chromosomes of an agreeing one hold 2, too few:
    its call at the cis peak is no longer raised. Chances are the zones' shares (the simulation is the subject of other tests)."""
    monkeypatch.setattr(hicfuse, "cis_chance", lambda *a, **k: {+1: 0.0, -1: 0.0})
    st, cen = _cis_state()
    at = {k: c if k < 3 else c + 80 if c < 300 else c - 80 for k, c in cen.items()}          # the others' cores 80 bins away
    cores = {f"c{k}": ((c - 1) * SIZE, (c + 1) * SIZE) for k, c in at.items()}
    r = hicfuse.cis_line(st, cores)
    assert r["status"] == "support" and r["sign"] == +1 and r["k"] == 3
    assert [r["status_by_chrom"][f"c{k}"] for k in range(12)] == ["fail"] * 3 + ["support"] * 9
    assert cis_entry(r, "c0") is None
    e = cis_entry(r, "c5")
    assert e["status"] == "support" and e["peak"] == r["peaks"]["c5"]
    call = (r["peaks"]["c0"][0] - 50_000, r["peaks"]["c0"][1] + 50_000)
    d = _cis_ev(cis_entry(r, "c0"), call)
    assert d["confidence"] == "low" and "cis_agrees" not in d["flags"]
    d = _cis_ev(dict(peak=r["peaks"]["c0"], status=r["status"], tested_by="ltr", interior=False), call)   # 0.16.0
    assert d["confidence"] == "medium" and "cis_agrees" in d["flags"]


def test_the_cis_line_that_the_other_chromosomes_pass_alone_still_acts(monkeypatch):
    monkeypatch.setattr(hicfuse, "cis_chance", lambda *a, **k: {+1: 0.0, -1: 0.0})
    st, cen = _cis_state()
    cores = {f"c{k}": ((c - 1) * SIZE, (c + 1) * SIZE) for k, c in cen.items()}
    r = hicfuse.cis_line(st, cores)
    assert r["status"] == "pass" and set(r["status_by_chrom"].values()) == {"pass"}
    call = (r["peaks"]["c0"][0] - 50_000, r["peaks"]["c0"][1] + 50_000)
    d = _cis_ev(cis_entry(r, "c0"), call)
    assert d["confidence"] == "medium" and "cis_agrees" in d["flags"]


# ---------------------------------------------------------------- the cis signal tested by trans contacts
def _both_state(n, partners, nb=200, seed=4):
    """As tests/fuse/test_cis_trans.py: cis peaks at each centromere, those of `partners` in strong trans contact."""
    seqs = {f"c{i}": nb * SIZE for i in range(n)}
    cen = {k: 40 + 9 * k for k in range(n)}
    b = Bins(seqs, SIZE)
    rng = np.random.default_rng(seed)
    dense = np.zeros((b.N, b.N))
    for k in range(n):
        off = int(b.off[k])
        for i in range(nb):
            for j in range(i + 1, min(nb, i + 12)):
                x = 200.0 / (j - i) * rng.gamma(20, 1 / 20) * (3.0 if i < cen[k] <= j else 1.0)
                dense[off + i, off + j] = dense[off + j, off + i] = x
    bg = rng.poisson(0.3, size=(b.N, b.N)).astype(float)
    bg = np.triu(bg, 1) + np.triu(bg, 1).T
    dense = np.where(b.chrom_of[:, None] == b.chrom_of[None, :], dense, bg)
    for x in partners:
        for y in partners:
            if x != y:
                dense[b.off[x] + cen[x] - 1:b.off[x] + cen[x] + 2, b.off[y] + cen[y] - 1:b.off[y] + cen[y] + 2] += 30.0
    W = sp.csr_matrix(dense)
    good, members = np.ones(b.N, dtype=bool), np.ones(n, dtype=bool)
    return dict(W=W, bins=b, good=good, mean=block_means(W, b, good, members), members=members, bias=np.ones(b.N))


def test_the_trans_tested_cis_line_acts_where_the_other_chromosomes_pass_without_it():
    """7 of 13 own trans peaks at their cis extremes pass (a majority); an agreeing chromosome's others hold 6 of 12, no majority, so
    the line acts only on the six whose others hold 7 of 12. Its result says so per chromosome."""
    st = _both_state(13, range(7), seed=2)            # a seed where no other chromosome agrees by chance
    r = hicfuse.cis_trans_line(st)
    assert r["status"] == "pass" and r["k"] == 7
    agree = {f"c{i}" for i in range(7)}
    assert {c for c, s in r["status_by_chrom"].items() if s == "pass"} == set(st["bins"].names) - agree
    assert set(r["peaks"]) == set(st["bins"].names) - agree
    assert cis_entry(r, "c0") is None and cis_entry(r, "c8")["status"] == "pass"


# ---------------------------------------------------------------- the arm fold and the trans arm alignment
BIN = 250_000


def _fold_table(peaks):
    return pd.DataFrame([dict(chrom=f"c{i}", callable=True, denovo_confidence="medium", denovo_peak_start=p,
                              denovo_peak_end=p + BIN) for i, p in enumerate(peaks)])


def _folds(n):
    return {f"c{i}": dict(point=50 * BIN, core_start=48 * BIN, core_end=53 * BIN, scanned=90 * BIN) for i in range(n)}


def _fold_ev(f, a=None, seq=(500_000, 1_000_000), conf="low"):
    s = dict(start=seq[0], end=seq[1], evidence="satellite", confidence=conf, flags=".", main_family="A", alt_start=-1, alt_end=-1)
    dn = 12_500_000
    h = dict(bin_size=BIN, gate_pass=False, gate_n=5, anch_strong=False, anch_level="low", anch_flags="", anch_peak_start=-1,
             anch_peak_end=-1, anch_core_start=-1, anch_core_end=-1, masked=(), denovo_peak_start=dn, denovo_peak_end=dn + BIN,
             denovo_core_start=dn - BIN, denovo_core_end=dn + 2 * BIN, denovo_confidence="high", denovo_flags="")
    return fuse.decide(fuse.Evidence("c0", 25_000_000, s, None, h, None, lambda a, b: 0, f, a))


def test_the_fold_on_a_majority_only_with_the_chromosome_itself_moves_no_call_there():
    """3 of 5 fold cores hold the prior-free peak (P 0.004 <= GATE_P / 2, a majority): without an agreeing chromosome 2 of 4 are no
    majority, so the fold and the prior-free peak no longer move that chromosome's call (fold_pair)."""
    v = fold.validate(_folds(5), _fold_table([50 * BIN] * 3 + [10 * BIN] * 2))
    assert v["status"] == "pass" and v["k"] == 3
    assert [v["status_by_chrom"][f"c{i}"] for i in range(5)] == ["fail"] * 3 + ["pass"] * 2
    assert v["p_by_chrom"]["c0"] > v["p"] and v["p_by_chrom"]["c4"] == v["p"]
    fo = dict(v, calls=_folds(5), align=None)
    f, a = axis_entries(fo, "c0")
    assert f["status"] == "fail" and a is None
    d = _fold_ev(f)
    assert d["call"] == (500_000, 1_000_000) and "fold_pair" not in d["flags"]
    d = _fold_ev(dict(f, status=v["status"]))                                         # 0.16.0: the genome's status
    assert "fold_pair" in d["flags"]


def test_a_fold_that_the_other_chromosomes_pass_alone_still_moves_the_call():
    v = fold.validate(_folds(6), _fold_table([50 * BIN] * 6))
    assert set(v["status_by_chrom"].values()) == {"pass"}
    f, _ = axis_entries(dict(v, calls=_folds(6), align=None), "c0")
    d = _fold_ev(f)
    assert "fold_pair" in d["flags"] and d["call"] == (12_250_000, 13_000_000)


def test_the_alignment_reads_the_folds_lean_without_the_chromosome():
    """The trans arm alignment stands in for a failed fold where the fold leans the same way (P <= AXIS_LEAN_P): that lean, too, is
    read on the other chromosomes."""
    fv = fold.validate(_folds(8), _fold_table([50 * BIN] * 2 + [10 * BIN] * 6))       # 2 of 8: fails, P 0.12 genome-wide
    av = fold.validate(_folds(8), _fold_table([50 * BIN] * 8), label="trans arm alignment")
    fo = dict(fv, calls={}, align=dict(av, calls=_folds(8)))
    _, a0 = axis_entries(fo, "c0")
    _, a5 = axis_entries(fo, "c5")
    assert a0["fold_status"] == "fail" and a0["fold_p"] == fv["p_by_chrom"]["c0"] > fv["p"]
    assert a5["fold_p"] == fv["p"] and a0["status"] == "pass"


# ---------------------------------------------------------------- the AT troughs and the repeat peaks tested against Hi-C
def _res(n, at=100_000):
    calls = {f"c{i}": dict(start=at, end=at + 3_000, gc=0.1, first=True) for i in range(n)}
    return dict(scales={1000: dict(calls=calls, p=0.01)}, best=1000, p=0.01, status="fail")


def test_troughs_validated_by_hic_only_with_the_chromosome_itself_raise_no_hic_call_there():
    """3 of 3 troughs agree with the prior-free peaks: a pass, which without one of them is 2 of 2, untested. The trough no longer
    pairs with that chromosome's Hi-C call (at_pair); with 6 of 6 the other chromosomes pass alone and it still does."""
    lengths = {f"c{i}": 1_000_000 for i in range(6)}
    for n, want in ((3, False), (6, True)):
        res = _res(n)
        peaks = {f"c{i}": (101_000, 101_500) for i in range(n)}
        v = atline.validate(res, peaks, dict.fromkeys(peaks, 100), lengths)
        assert v["status"] == "pass" and v["k"] == n
        out = dict(res, window=1000, hic=v, cis=None)
        x = atline.on_chromosome(out, "c0")
        assert x["agrees"] and x["validated"] is want
        h = dict(callable=True, bin_size=1000, gate_pass=False, gate_support=False, gate_genome_untested=True, anch_strong=False,
                 anch_level="none", anch_peak_start=-1, anch_peak_end=-1, anch_core_start=-1, anch_core_end=-1, masked=(),
                 denovo_peak_start=101_000, denovo_peak_end=101_500, denovo_core_start=101_000, denovo_core_end=101_500,
                 denovo_confidence="high", denovo_flags="")
        s = dict(start=-1, end=-1, evidence="none", confidence="none", flags=".", main_family=".")
        d = fuse.decide(fuse.Evidence("c0", 1_000_000, s, None, h, None, lambda a, b: 0, at=x))
        assert ("at_pair" in d["flags"]) is want and d["confidence"] == ("medium" if want else "low")


def test_a_trough_test_on_a_minority_supports_only_where_the_others_do():
    """3 of 6 agree: support genome-wide; an agreeing chromosome's others hold 2 of 5, too few (fail); a disagreeing one's others
    hold 3 of 5, a majority, yet it keeps the genome's support (its disagreement counts against the line, as in the gate)."""
    lengths = {f"c{i}": 1_000_000 for i in range(6)}
    res = _res(6)
    peaks = {f"c{i}": (101_000, 101_500) if i < 3 else (800_000, 800_500) for i in range(6)}
    v = atline.validate(res, peaks, dict.fromkeys(peaks, 100), lengths)
    assert v["status"] == "support"
    assert [v["status_by_chrom"][f"c{i}"] for i in range(6)] == ["fail"] * 3 + ["support"] * 3
    out = dict(res, window=1000, hic=v, cis=None)
    assert not atline.on_chromosome(out, "c0")["support"] and atline.on_chromosome(out, "c4")["support"]
    rep = dict(res, best=1000, hic=v, cis=None)
    assert not repline.on_chromosome(rep, "c0")["support"] and repline.on_chromosome(rep, "c4")["support"]
    cis = dict(res, window=1000, hic=None, cis=v)
    assert not atline.on_chromosome(cis, "c1")["cis_support"] and atline.on_chromosome(cis, "c5")["cis_support"]
