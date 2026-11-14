"""Hi-C anchored on LTR cores (hicfuse.anchored_on_cores), the cis signal (hicfuse.cis_line) and the table of line agreements (linetests)."""
import numpy as np
import pandas as pd
import scipy.sparse as sp

from cenfuse import hicfuse, linetests
from cenfuse.hic.bins import Bins
from cenfuse.hic.score import block_means

SIZE = 10_000
NB = 60
SEQS = {f"c{i}": NB * SIZE for i in range(8)}
CEN = {k: 12 + 4 * k for k in range(8)}        # the centromere bin of each chromosome


def state(strength=30.0, seed=1):
    """A symmetric contact map: background 1 + noise, centromeres of different chromosomes in strong contact (+/- 1 bin)."""
    b = Bins(SEQS, SIZE)
    rng = np.random.default_rng(seed)
    d = rng.poisson(1.0, size=(b.N, b.N)).astype(float)
    d = np.triu(d) + np.triu(d, 1).T
    for x, cx in CEN.items():
        for y, cy in CEN.items():
            if x != y:
                for i in (-1, 0, 1):
                    for j in (-1, 0, 1):
                        d[b.off[x] + cx + i, b.off[y] + cy + j] += strength
    W = sp.csr_matrix(d)
    good, members = np.ones(b.N, dtype=bool), np.ones(8, dtype=bool)
    return dict(W=W, bins=b, good=good, mean=block_means(W, b, good, members), members=members)


def cores(offset=0):
    return {f"c{k}": ((c + offset) * SIZE, (c + offset + 2) * SIZE) for k, c in CEN.items()}


def test_peaks_at_the_cores_pass():
    r = hicfuse.anchored_on_cores(state(), cores(1))
    assert r["status"] == "pass" and r["k"] == 8 and r["p"] < 1e-6 and set(r["peaks"]) == set(SEQS)


def test_cores_elsewhere_fail():
    r = hicfuse.anchored_on_cores(state(), cores(20))
    assert r["status"] == "fail" and r["k"] <= 2


def test_cores_right_on_a_minority_support():
    """Beyond chance needs a zone that is small on the chromosome: 12 chromosomes of 200 bins, 4 cores right."""
    n, nb = 12, 200
    seqs = {f"c{i}": nb * SIZE for i in range(n)}
    cen = {k: 20 + 10 * k for k in range(n)}
    b = Bins(seqs, SIZE)
    d = np.random.default_rng(3).poisson(0.5, size=(b.N, b.N)).astype(float)
    d = np.triu(d) + np.triu(d, 1).T
    for x, cx in cen.items():
        for y, cy in cen.items():
            if x != y:
                d[b.off[x] + cx - 1:b.off[x] + cx + 2, b.off[y] + cy - 1:b.off[y] + cy + 2] += 200.0
    W = sp.csr_matrix(d)
    good, members = np.ones(b.N, dtype=bool), np.ones(n, dtype=bool)
    st = dict(W=W, bins=b, good=good, mean=block_means(W, b, good, members), members=members)
    c = {f"c{k}": ((cen[k] + (0 if k < 4 else 60)) * SIZE, (cen[k] + (2 if k < 4 else 62)) * SIZE) for k in range(n)}
    r = hicfuse.anchored_on_cores(st, c)
    assert r["status"] == "support" and r["k"] == 4 and set(r["peaks"]) == {"c0", "c1", "c2", "c3"}


def test_no_contact_signal_fails():
    for seed in range(5):          # anchors that chase noise would find each other: they stay at the cores
        assert hicfuse.anchored_on_cores(state(strength=0.0, seed=seed), cores(1))["status"] == "fail"


def _telomere_state(seed=2):
    """No centromere signal; the last two bins of every chromosome contact each other (telomeres clustered in the nucleus)."""
    b = Bins(SEQS, SIZE)
    d = np.random.default_rng(seed).poisson(1.0, size=(b.N, b.N)).astype(float)
    d = np.triu(d) + np.triu(d, 1).T
    for x in range(8):
        for y in range(8):
            if x != y:
                d[b.off[x] + NB - 2:b.off[x] + NB, b.off[y] + NB - 2:b.off[y] + NB] += 30.0
    W = sp.csr_matrix(d)
    good, members = np.ones(b.N, dtype=bool), np.ones(8, dtype=bool)
    return dict(W=W, bins=b, good=good, mean=block_means(W, b, good, members), members=members)


def test_cores_at_the_chromosome_ends_do_not_count():
    """Telomeres cluster whatever the centromere: cores at the ends would agree with any such map (held-out R. breviuscula)."""
    r = hicfuse.anchored_on_cores(_telomere_state(), {f"c{k}": ((NB - 2) * SIZE, NB * SIZE) for k in range(8)})
    assert r["status"] == "fail" and r["k"] == 0 and r["n"] == 8


def test_interior_cores_still_pass_beside_clustered_telomeres():
    st = state()
    W = st["W"].toarray()
    b = st["bins"]
    for x in range(8):
        for y in range(8):
            if x != y:
                W[b.off[x] + NB - 2:b.off[x] + NB, b.off[y] + NB - 2:b.off[y] + NB] += 30.0
    st["W"] = sp.csr_matrix(W)
    st["mean"] = block_means(st["W"], b, st["good"], st["members"])
    r = hicfuse.anchored_on_cores(st, cores(1))
    assert r["status"] == "pass" and r["k"] == 8


# ---------------------------------------------------------------- linetests
def row(chrom, seq, te=None, hic=None, at=None, ev="satellite", over=None, L=100_000_000):
    d = dict(chrom=chrom, chrom_len=L, seq_evidence=ev, seq_start=seq[0], seq_end=seq[1], seq_overruled_start=".",
             seq_overruled_end=".", te_start=".", te_end=".", hic_anch_peak_start=".", hic_anch_peak_end=".",
             hic_denovo_peak_start=".", hic_denovo_peak_end=".", fold_core_start=".", fold_core_end=".", align_core_start=".",
             align_core_end=".", at_start=".", at_end=".")
    if te:
        d.update(te_start=te[0], te_end=te[1])
    if hic:
        d.update(hic_anch_peak_start=hic[0], hic_anch_peak_end=hic[1])
    if at:
        d.update(at_start=at[0], at_end=at[1])
    if over:
        d.update(seq_overruled_start=over[0], seq_overruled_end=over[1])
    return d


def test_agreement_is_counted_against_chance():
    t = pd.DataFrame([row(f"c{i}", (10_000_001, 11_000_000), te=(10_500_001, 10_600_000)) for i in range(10)])
    r = linetests.table(t, None).iloc[0]
    assert (r.line_a, r.line_b, r.chromosomes, r.agree) == ("ltr", "satellite", 10, 10)
    assert r.expected < 0.5 and r.p < 1e-6 and r.verdict == "agree on most, beyond chance"


def test_disagreeing_lines_are_said():
    t = pd.DataFrame([row(f"c{i}", (10_000_001, 11_000_000), te=(60_000_001, 60_100_000)) for i in range(10)])
    r = linetests.table(t, None).iloc[0]
    assert r.agree == 0 and r.verdict == "no agreement beyond chance"


def test_a_hic_line_is_allowed_a_bin():
    t = pd.DataFrame([row(f"c{i}", (10_000_001, 11_000_000), hic=(11_100_001, 11_350_000)) for i in range(5)])
    assert linetests.table(t, 250_000).iloc[0].agree == 5
    assert linetests.table(t, 50_000).iloc[0].agree == 0


def test_the_satellite_line_is_the_overruled_array_where_a_core_replaced_it():
    r = row("c", (60_000_001, 60_100_000), te=(60_000_001, 60_100_000), ev="retrotransposon", over=(10_000_001, 11_000_000))
    x = linetests.line_calls(pd.DataFrame([r]).itertuples(index=False).__next__())
    assert x["satellite"] == (10_000_000, 11_000_000) and x["ltr"] == (60_000_000, 60_100_000)


def test_few_chromosomes_are_said():
    t = pd.DataFrame([row("c", (1, 2), te=(1, 2))])
    assert linetests.table(t, None).iloc[0].verdict == "too few chromosomes"


# ---------------------------------------------------------------- which end (hicfuse.end_test)
def _table(calls):
    rows = [dict(chrom=f"c{i}", callable=True, seq_start=s, seq_end=e) for i, (s, e) in enumerate(calls)]
    return pd.DataFrame(rows)


def _ends(n, high_at_call, tot=100, L=SIZE * 100):
    """n chromosomes of 100 bins, each with its call at the start; the score high at the start where high_at_call says."""
    seqs = {f"c{i}": L for i in range(n)}
    b = Bins(seqs, SIZE)
    S = np.ones(b.N)
    for k in range(n):
        v = S[b.sl(k)]
        v[:5 if high_at_call[k] else 0] = 3.0
        if not high_at_call[k]:
            v[-5:] = 3.0
    return b, S


def test_the_end_holding_the_call_scores_higher_on_most_chromosomes():
    b, S = _ends(8, [True] * 7 + [False])
    t = hicfuse.end_test(_table([(0, 200_000)] * 8), S, b)
    assert t.end_higher.tolist() == [True] * 7 + [False]
    assert t.end_k.iat[0] == 6 and t.end_n.iat[0] == 7 and not t.end_pass.iat[0]      # 6 of 7: P 0.06
    b, S = _ends(12, [True] * 12)
    t = hicfuse.end_test(_table([(0, 200_000)] * 12), S, b)
    assert t.end_pass.all() and t.end_k.iat[0] == 11 and t.end_near.iat[0] == 3.0 and t.end_other.iat[0] == 1.0


def test_a_call_at_the_far_end_reads_the_far_end():
    b, S = _ends(4, [False] * 4)
    t = hicfuse.end_test(_table([(SIZE * 100 - 200_000, SIZE * 100)] * 4), S, b)
    assert t.end_higher.all() and t.end_near.iat[0] == 3.0


def test_interior_calls_are_not_tested():
    b, S = _ends(4, [True] * 4)
    t = hicfuse.end_test(_table([(400_000, 500_000)] * 4), S, b)
    assert t.end_near.isna().all() and t.end_n.iat[0] == 0 and not t.end_pass.any()


def test_ends_cannot_be_told_at_random():
    from cenfuse import fuse
    h = dict(callable=True, bin_size=SIZE, gate_pass=False, gate_support=False, anch_strong=False, anch_level="none", anch_peak_start=-1,
             anch_peak_end=-1, anch_core_start=-1, anch_core_end=-1, masked=(), end_pass=True, end_higher=True)
    s = dict(start=0, end=200_000, evidence="satellite", confidence="medium", flags=".", main_family="A")
    d = fuse.decide(fuse.Evidence("c", 1_000_000, s, None, h, None, lambda a, b: 300_000))
    assert d["relation"] == "concordant" and d["confidence"] == "high" and "hic_end" in d["flags"] and "hic" in d["support"]
    for x in (dict(h, end_pass=False), dict(h, end_higher=False)):
        d = fuse.decide(fuse.Evidence("c", 1_000_000, s, None, x, None, lambda a, b: 300_000))
        assert "hic_end" not in d["flags"] and d["confidence"] == "medium"
    d = fuse.decide(fuse.Evidence("c", 1_000_000, dict(s, alt_start=200_000, alt_end=300_000), None, h, None, lambda a, b: 300_000))
    assert "hic_end" not in d["flags"]                                                   # no mappable DNA between the two


# ---------------------------------------------------------------- the cis signal (hicfuse.cis_line)
def _cis_state(sign, n=10, nb=120, base=200.0, seed=4):
    """Cis maps with distance decay; at each centromere the contact crossing the bin rises (sign +1) or falls (-1)."""
    seqs = {f"c{i}": nb * SIZE for i in range(n)}
    cen = {k: 20 + 6 * k for k in range(n)}
    b = Bins(seqs, SIZE)
    rng = np.random.default_rng(seed)
    rows, cols, vals = [], [], []
    for k in range(n):
        off = int(b.off[k])
        for i in range(nb):
            for j in range(i + 1, min(nb, i + 12)):
                x = base / (j - i) * rng.gamma(20, 1 / 20)
                if i < cen[k] <= j:
                    x *= 3.0 if sign > 0 else 0.2
                rows += [off + i, off + j]
                cols += [off + j, off + i]
                vals += [x, x]
    W = sp.csr_matrix((vals, (rows, cols)), shape=(b.N, b.N))
    st = dict(W=W, bins=b, good=np.ones(b.N, dtype=bool), members=np.ones(n, dtype=bool), bias=np.ones(b.N))
    return st, cen


def test_the_sign_comes_from_the_ltr_cores():
    for sign in (+1, -1):
        st, cen = _cis_state(sign)
        cores = {f"c{k}": ((c - 1) * SIZE, (c + 1) * SIZE) for k, c in cen.items()}
        r = hicfuse.cis_line(st, cores)
        assert r["status"] == "pass" and r["sign"] == sign and r["k"] == 10
        assert all(abs(r["peaks"][f"c{k}"][0] // SIZE - c) <= 1 for k, c in cen.items())


def test_cores_elsewhere_fail_the_cis_test():
    st, cen = _cis_state(+1)
    cores = {f"c{k}": ((c + 40) * SIZE, (c + 42) * SIZE) for k, c in cen.items()}
    assert hicfuse.cis_line(st, cores)["status"] == "fail"


def _noise_state(low=0.15, n=10, nb=120, level=300.0, seed=6):
    """Counts with distance decay and coverage, no locus effect: coverage falls in a stretch around each core, balanced by ICE."""
    from cenfuse.hic.balance import balanced, ice
    seqs = {f"c{i}": nb * SIZE for i in range(n)}
    cen = {k: 20 + 6 * k for k in range(n)}
    b = Bins(seqs, SIZE)
    rng = np.random.default_rng(seed)
    cov = np.ones(b.N)
    for k, c in cen.items():
        cov[int(b.off[k]) + c - 2:int(b.off[k]) + c + 3] = low
    rows, cols, vals = [], [], []
    for k in range(n):
        off = int(b.off[k])
        for i in range(nb):
            for j in range(i + 1, min(nb, i + 12)):
                x = rng.poisson(level * cov[off + i] * cov[off + j] / (j - i))
                rows += [off + i, off + j]
                cols += [off + j, off + i]
                vals += [x, x]
    M = sp.csr_matrix((vals, (rows, cols)), shape=(b.N, b.N), dtype=float)
    bias, good, _, _ = ice(M, min_frac=0.01)
    return dict(W=balanced(M, bias, good), bins=b, good=good, members=np.ones(n, dtype=bool), bias=bias), cen


def test_noise_at_low_coverage_does_not_pass():
    """The extremes fall in the low-coverage stretches, here the cores, more often than their share of the chromosome: without
    cis_chance this map would pass."""
    st, cen = _noise_state()
    c = {f"c{k}": ((x - 1) * SIZE, (x + 1) * SIZE) for k, x in cen.items()}
    r = hicfuse.cis_line(st, c)
    assert r["status"] == "fail" and r["k"] >= 3
    lo, hi = hicfuse._core_zone(st["bins"], 0, *c["c0"])
    sim = hicfuse.cis_chance(st["W"], st["bins"], st["good"], st["bias"], 0, lo, hi)
    assert sim[+1] > 2 * (hi - lo + 1) / 112 and sim[-1] > 2 * (hi - lo + 1) / 112


def test_a_real_signal_passes_despite_low_coverage():
    """Coverage falls at the cores, and the contact crossing them rises threefold. Balancing flattens the rise and leaves shoulders,
    so either extreme can mark the locus: the test is that the locus is found."""
    from cenfuse.hic.balance import balanced, ice
    st, cen = _cis_state(+1, base=300.0)
    b = st["bins"]
    cov = np.ones(b.N)
    for k, c in cen.items():
        cov[int(b.off[k]) + c - 2:int(b.off[k]) + c + 3] = 0.3
    D = sp.diags(cov)
    M = (D @ st["W"] @ D).tocsr()
    M.data = np.random.default_rng(2).poisson(M.data).astype(float)
    bias, good, _, _ = ice(M, min_frac=0.01)
    s2 = dict(st, W=balanced(M, bias, good), good=good, bias=bias)
    c = {f"c{k}": ((x - 1) * SIZE, (x + 1) * SIZE) for k, x in cen.items()}
    r = hicfuse.cis_line(s2, c)
    assert r["status"] == "pass" and r["k"] >= 8


def test_masked_bins_make_no_extreme():
    """Beside a masked bin the box holds other distances than usual, and contact falls steeply with distance: observed over expected
    does not move."""
    nb = 60
    b = Bins({"c0": nb * SIZE, "c1": nb * SIZE}, SIZE)
    d = np.zeros((b.N, b.N))
    for off in (0, nb):
        for i in range(nb):
            for j in range(i + 1, nb):
                d[off + i, off + j] = d[off + j, off + i] = 1000.0 / (j - i) ** 1.5
    good = np.ones(b.N, dtype=bool)
    good[25:28] = False
    d[~good, :] = 0
    d[:, ~good] = 0
    v = hicfuse.cis_profile(sp.csr_matrix(d), b, good, 0)
    assert np.nanmax(np.abs(v - 1)) < 1e-9 and np.isnan(v[:4]).all() and np.isnan(v[-4:]).all()
    assert np.isfinite(v[20:32]).all()


def _cis_ev(conf, peak, status="support", seq=(1_000_000, 1_400_000)):
    from cenfuse import fuse
    h = dict(callable=True, bin_size=SIZE, gate_pass=False, gate_support=False, anch_strong=False, anch_level="none", anch_peak_start=-1,
             anch_peak_end=-1, anch_core_start=-1, anch_core_end=-1, masked=(), cis=dict(peak=peak, status=status))
    s = dict(start=seq[0], end=seq[1], evidence="satellite", confidence=conf, flags=".", main_family="A") if seq else \
        dict(start=-1, end=-1, evidence="none", confidence="none", flags=".", main_family=".")
    return fuse.decide(fuse.Evidence("c", 50_000_000, s, None, h, None, lambda a, b: 0))


def test_cis_at_the_call_raises_it():
    d = _cis_ev("low", (1_100_000, 1_110_000))
    assert d["relation"] == "concordant" and d["confidence"] == "medium" and "cis_agrees" in d["flags"] and "hic" in d["support"]


def test_cis_elsewhere_moves_only_a_low_call_and_it_stays_low():
    d = _cis_ev("low", (30_000_000, 30_010_000), status="pass")
    assert d["call"] == (30_000_000, 30_010_000) and d["confidence"] == "low" and "cis_only" in d["flags"]
    d = _cis_ev("low", (30_000_000, 30_010_000))           # 0.16: at the support tier it raises what it agrees with, moves none
    assert d["call"] == (1_000_000, 1_400_000) and d["confidence"] == "low" and "cis_only" not in d["flags"]
    assert "cis_elsewhere" in d["flags"] and d["relation"] == "sequence_only"
    for conf in ("medium", "high"):
        d = _cis_ev(conf, (30_000_000, 30_010_000))
        assert d["call"] == (1_000_000, 1_400_000) and d["confidence"] == conf and "cis_elsewhere" in d["flags"] and "hic" in d["against"]


def test_a_passing_cis_line_elsewhere_lowers_a_call_without_moving_it():
    for conf, down in (("medium", "low"), ("high", "medium")):
        d = _cis_ev(conf, (30_000_000, 30_010_000), status="pass")
        assert d["call"] == (1_000_000, 1_400_000) and d["confidence"] == down and d["relation"] == "contested"
        assert "cis_elsewhere" in d["flags"] and "hic" in d["against"]


def test_cis_within_the_resolution_of_hi_c_neither_supports_nor_opposes():
    # SIZE-sized bins: 1,400,000 + 3 bins lies beyond the 1 bin of agreement and within the 4 bins of opposition
    for conf in ("low", "medium", "high"):
        d = _cis_ev(conf, (1_430_000, 1_440_000), status="pass")
        assert d["call"] == (1_000_000, 1_400_000) and d["confidence"] == conf and "cis_near" in d["flags"]
        assert "hic" not in d["support"] and "hic" not in d["against"] and d["relation"] != "contested"


def test_cis_cannot_raise_a_call_it_cannot_tell_from_its_competitor():
    from cenfuse import fuse
    h = dict(callable=True, bin_size=SIZE, gate_pass=False, gate_support=False, anch_strong=False, anch_level="none", anch_peak_start=-1,
             anch_peak_end=-1, anch_core_start=-1, anch_core_end=-1, masked=(), cis=dict(peak=(1_100_000, 1_110_000), status="pass"))
    s = dict(start=1_000_000, end=1_400_000, evidence="satellite", confidence="low", flags=".", main_family="A", alt_start=1_400_000, alt_end=1_500_000)
    d = fuse.decide(fuse.Evidence("c", 50_000_000, s, None, h, None, lambda a, b: 0))
    assert d["confidence"] == "low" and "cis_agrees" not in d["flags"] and "hic" not in d["support"] and "hic_alt_ambiguous" in d["flags"]


def test_cis_places_a_low_call_where_there_is_none():
    d = _cis_ev("none", (30_000_000, 30_010_000), seq=None)
    assert d["call"] == (30_000_000, 30_010_000) and d["confidence"] == "low" and d["relation"] == "functional_only"


def test_the_null_and_the_cis_simulations_do_not_depend_on_the_thread_count():
    """future.md #62: the null replicates' positions are drawn first, in one serial order, and each chromosome's cis simulations
    draw from its own stream, so the thread pool changes no P value. The expected values are those of 0.14.0's serial functions
    (code/rel14g) on these maps."""
    st, cen = _noise_state()
    c = {f"c{k}": ((x - 1) * SIZE, (x + 1) * SIZE) for k, x in cen.items()}
    for t in (1, 4):
        r = hicfuse.cis_line(st, c, t)
        assert {k: r[k] for k in ("k", "n", "p", "sign", "status")} == \
            {"k": 9, "n": 10, "p": 0.6802696115190116, "sign": -1, "status": "fail"}
    s = state(strength=3.0)
    W, b, good = s["W"], s["bins"], s["good"]
    mean = block_means(W, b, good, s["members"])
    ratios = {k: 1.45 + 0.04 * k for k in CEN}
    hits = {0: 0, 1: 60, 2: 23, 3: 0, 4: 4, 5: 10, 6: 5, 7: 24}            # replicates reaching each ratio, of 60
    for t in (1, 4):
        p = hicfuse.null_all(W, b, good, mean, s["members"], dict(CEN), dict(CEN), ratios, n=60, threads=t)
        assert p == {k: (1 + h) / 61 for k, h in hits.items()}


def test_anchored_reuses_the_engines_bootstrap_only_for_its_own_windows(monkeypatch):
    """future.md #62: the anchored test takes the bootstrap the engine ran on its windows (state['boot_extra']) when they are the
    windows it would use, and its own otherwise; both give the same table."""
    s = state(strength=10.0)
    N = s["bins"].N
    st = dict(s, M=s["W"], hub=np.zeros(N, dtype=bool), rdna={})
    seq = {f"c{k}": ((c - 1) * SIZE, (c + 2) * SIZE, "high") for k, c in CEN.items()}
    win = hicfuse._partners(st["bins"], st["members"], st["good"], seq)[1]
    own, _ = hicfuse.anchored(st, seq, 2)                                     # no engine bootstrap: its own
    calls = []
    real = hicfuse.bootstrap_all
    monkeypatch.setattr(hicfuse, "bootstrap_all", lambda *a, **k: calls.append(1) or real(*a, **k))
    pre = dict(st, boot_extra={"anchored": real(st["M"], st["bins"], st["members"], win, 2)}, extra_windows={"anchored": win})
    reused, _ = hicfuse.anchored(pre, seq, 2)
    assert not calls and reused.equals(own)
    other = {k: w[:-1] for k, w in win.items()}                               # windows that are not the ones it would use
    moved, _ = hicfuse.anchored(dict(pre, extra_windows={"anchored": other}), seq, 2)
    assert calls and moved.equals(own)
