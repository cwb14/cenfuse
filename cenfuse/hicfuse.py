"""The Hi-C line.

1. The Hi-C engine's prior-free call, unchanged (outputs under OUT/hic/).
2. The anchored test, on the same balanced matrix. Partner windows sit at the other chromosomes' sequence calls: the
   2 HALF_WIN + 1 mappable bins nearest the call (the midpoint bin +/- the Hi-C engine's HALF_WIN when the call is
   mappable; the array's mappable flanks when it is not). Each chromosome's score is the mean O/E with those windows; its own call is
   never a partner, so the test is leave-one-chromosome-out by construction. Which trans cluster dominates the
   genome then no longer matters: the score asks whether a place contacts the other chromosomes' (sequence-defined)
   centromeres. Every anchored peak gets the Hi-C engine's own tests: Poisson bootstrap with the windows fixed,
   leave-one-partner-out, and the position-aware null.
3. The gate, per chromosome x: on the other chromosomes with a reference-grade sequence call, does the anchored peak
   land on that call more often than chance (Poisson-binomial tail <= GATE_P, and on a majority)? Only then is
   Hi-C allowed to change x's call. Below it, the support tier: agreements beyond chance (the same tail, at least
   GATE_MIN_REF of them) but not on a majority, counting only agreements at significant anchored peaks (Hi-C level
   medium or high: at least two of the engine's three tests). Hi-C then carries centromere information but its peak
   is often elsewhere, so it may only support x's call where its own significant peak agrees (fuse.decide): a peak
   lands on a call by chance on ~1-2% of chromosomes, so an agreement is strong evidence, while a peak elsewhere is
   weak evidence against. Where 0.9's gate failed only for want of a majority (rhesus, goat, giraffe, horse, EquCab3,
   moss, chicken), an agreeing peak marked a right call 61 of 64 times. Without the majority, agreements must stand on
   their own: with randomized Hi-C (partners re-paired, coverage kept) peaks still land on the calls 3-6 times per
   genome (chance ~1; low-coverage bins at unmappable arrays' edges carry noisy maxima), almost all at level low or
   none. Since 0.13 a strong peak on the sequence line's own competing cluster may also choose that cluster (fuse).

Distances are counted in mappable bins (`visible_gap`). Hi-C sees nothing inside a megabase array of near-identical
repeat units: no read maps uniquely there, so its bins are masked. In T2T mammals the CENP-A domain lay inside such
an unmappable block on 139/157 chromosomes and the anchored peak sat at the block's edge, 1-5 Mb away (benchmark
v2). For agreement (support, the gate) up to HIC_FAR_BINS unmappable bins between a peak and a call (the method's
resolution) count as no distance and longer blind stretches count beyond that; for opposition unmappable bins never
count. Hi-C cannot tell which side of a long blind stretch the centromere is on, so there it neither supports nor
opposes (T2T bonobo: a call 81 unmappable bins away, at the far edge of a 20-Mb pericentric array, was the only wrong
call that full credit let Hi-C support, of 92 with CENP-A truth)."""
from __future__ import annotations

from argparse import Namespace
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

from .hic import cli as hic_cli
from .hic import params as hp
from .hic.confidence import bootstrap, leave_one_out, level
from .hic.io import InputError
from .hic.locate import peak_core
from .hic.qc import call_flags
from .hic.score import cluster_score, partner_oe, window_bins
from .params import ANCHOR_LEVELS, FC_LEVEL, GATE_MIN_REF, GATE_P, HIC_FAR_BINS, HIC_TOL_BINS, NOT_STRONG_FLAGS, \
    STRONG_HIC_LEVEL
from .util import CenfuseError, gap, log, poisson_binomial_sf

HIC_COLS = ["chrom", "callable", "bin_size", "seq_start", "seq_end", "seq_confidence", "partner",
            "anch_peak_start", "anch_peak_end", "anch_core_start", "anch_core_end", "anch_ratio", "anch_boot_in_core",
            "anch_lopo_stable", "anch_null_p", "anch_level", "anch_flags", "anch_strong", "anch_agrees_seq", "chance",
            "gate_k", "gate_n", "gate_p", "gate_pass", "gate_support", "gate_refs", "masked",
            "end_near", "end_other", "end_higher", "end_k", "end_n", "end_p", "end_pass",
            "denovo_peak_start", "denovo_peak_end", "denovo_core_start", "denovo_core_end", "denovo_confidence",
            "denovo_flags", "denovo_agrees_seq", "denovo_agrees_anch"]


def run_denovo(pairs, fai, prefix, satellites, threads: int, seq: dict | None = None) -> dict:
    """The Hi-C engine's `call`, as its own command runs it; returns its rows and in-memory state. seq: the sequence calls,
    whose partner windows (anchored) the engine's bootstrap scores on its own resamples (state['boot_extra'])."""
    files = hic_cli.pair_files(pairs)              # several = replicate libraries of one plant, merged (0.10.1)
    extra = None if seq is None else (lambda bins, members, good: {"anchored": _partners(bins, members, good, seq)[1]})
    args = Namespace(pairs=str(files[0]) if len(files) == 1 else [str(f) for f in files], fai=str(fai), o=str(prefix),
                     repeats=None,
                     satellites=[str(s) for s in satellites] or None, t=max(1, threads), verbose=False,
                     chromosome_set=False,          # the .fai is cenfuse's chromosome set already
                     extra_windows=extra)
    try:
        return hic_cli.run_call(args)
    except InputError as e:
        raise CenfuseError(f"Hi-C {pairs}: {e}") from None


def bootstrap_all(M, bins, members, win: dict, threads: int, n: int = hp.N_BOOT, seed: int = hp.SEED) -> dict:
    """The Hi-C engine's bootstrap (Poisson resamples, re-balanced, partner windows fixed: hic.confidence.bootstrap) with
    explicit partner windows (bins a resample masks drop out of them), and the peak of every callable chromosome, not only of
    those holding a window. A run computes it inside the engine's own bootstrap, on the same resamples (run_denovo: seq)."""
    return bootstrap(M, bins, members, {}, n=n, seed=seed, threads=threads, extra={"w": win})[1]["w"]


def _same_windows(a, b) -> bool:
    return a is not None and set(a) == set(b) and all(np.array_equal(a[k], b[k]) for k in a)


def lopo_all(oe, bins, good, members) -> dict:
    """Peak of every callable chromosome with each partner window left out in turn (the Hi-C engine's leave_one_out, for every
    member chromosome)."""
    return leave_one_out(oe, bins, good, members, {int(k): None for k in np.where(members)[0]})


def null_all(W, bins, good, mean, members, est, peaks, ratios, n: int = hp.N_NULL, seed: int = hp.SEED,
             threads: int = 1) -> dict:
    """The Hi-C engine's position-aware null, for every chromosome with a peak: partner windows moved to random bins at least
    NULL_EXCL from that partner's anchor and its ends; P = (1 + replicates reaching the observed peak ratio within
    NULL_REACH bins of the peak) / (1 + n). NaN when untestable. Every replicate's positions are drawn first, in the order
    of one serial loop, and the replicates are then scored on the thread pool (future.md #62): the P values do not depend
    on the thread count."""
    rng = np.random.default_rng(seed)
    elig = {}
    for y, c in est.items():
        loc = np.arange(hp.NULL_EXCL, int(bins.nb[y]) - hp.NULL_EXCL)
        loc = loc[np.abs(loc - c) >= hp.NULL_EXCL]
        if len(loc):
            elig[y] = loc
    if len(elig) < 2:
        return {x: float("nan") for x in peaks}
    cents = [{y: int(rng.choice(loc)) for y, loc in elig.items()} for _ in range(n)]

    def one(cent) -> dict:
        """{x: whether the replicate reaches x's peak ratio (None: no finite score on x)}."""
        S = cluster_score(partner_oe(W, bins, good, mean, window_bins(bins, good, cent)), bins, good, members)
        out = {}
        for x in peaks:
            v = S[bins.sl(x)]
            if not np.isfinite(v).any():
                out[x] = None
                continue
            med = np.nanmedian(v)
            lo, hi = max(peaks[x] - hp.NULL_REACH, 0), min(peaks[x] + hp.NULL_REACH + 1, len(v))
            win = v[lo:hi]
            out[x] = bool(med > 0 and np.isfinite(win).any() and np.nanmax(win) / med >= ratios[x])
        return out

    with ThreadPoolExecutor(max(threads, 1)) as ex:
        reps = list(ex.map(one, cents))
    hits = {x: sum(r[x] is True for r in reps) for x in peaks}
    seen = {x: sum(r[x] is not None for r in reps) for x in peaks}
    return {x: (1 + hits[x]) / (1 + n) if np.isfinite(ratios[x]) and seen[x] else float("nan") for x in peaks}


def _bin_range(bins, k: int, s: int, e: int) -> tuple:
    nb = int(bins.nb[k])
    return min(nb - 1, s // bins.size), min(nb - 1, max(s, e - 1) // bins.size)


def visible_dist(bins, good, k: int, s: int, e: int, cap: int | None = None) -> np.ndarray:
    """Per bin of chromosome k: the number of mappable bins strictly between it and [s, e) (0 inside), plus the
    unmappable bins between them beyond the first `cap` (None: unmappable bins never count)."""
    nb, off = int(bins.nb[k]), int(bins.off[k])
    lo, hi = _bin_range(bins, k, s, e)
    cg = np.concatenate([[0], np.cumsum(good[off:off + nb])])       # cg[j] = mappable bins in [0, j)
    j = np.arange(nb)
    d = np.zeros(nb, dtype=np.int64)
    L, Rt = j < lo, j > hi
    d[L] = cg[lo] - cg[j[L] + 1]
    d[Rt] = cg[j[Rt]] - cg[hi + 1]
    if cap is not None:
        between = np.zeros(nb, dtype=np.int64)
        between[L] = lo - j[L] - 1
        between[Rt] = j[Rt] - hi - 1
        d = d + np.maximum(0, between - d - cap)
    return d


def visible_window(bins, good, k: int, s: int, e: int, half: int = hp.HALF_WIN) -> np.ndarray:
    """Global bins of the partner window of call [s, e) on chromosome k. When most bins of the Hi-C engine's window
    (midpoint +/- half) are mappable, its mappable bins; otherwise the 2 half + 1 mappable bins nearest the call
    (mappable-bin distance, then distance to the midpoint), at most `half` mappable bins outside it: for a call inside
    an unmappable array, the array's flanks, the nearest places Hi-C can see (T2T macaques: 0/22 engine windows held
    a mappable bin). Padding every window, mostly mappable ones too, lost anchored agreements in CHM13, chimpanzee
    and S. orangutan; keeping every window that holds any mappable bin lost 11 of 97 in six imported mammals."""
    off, nb = int(bins.off[k]), int(bins.nb[k])
    mid = min(nb - 1, ((s + e) // 2) // bins.size)
    w = off + np.arange(max(mid - half, 0), min(mid + half, nb - 1) + 1)
    if 2 * int(good[w].sum()) > len(w):          # mostly mappable: the Hi-C engine's own window
        return w[good[w]]
    d = visible_dist(bins, good, k, s, e)
    j = np.flatnonzero(good[off:off + nb] & (d <= half))
    order = np.lexsort((np.abs(j - mid), d[j]))
    return np.sort(off + j[order[:2 * half + 1]])


def masked_runs(bins, good, k: int) -> tuple:
    """Unmappable stretches of chromosome k as ((start, end), ...) in bp."""
    off, nb = int(bins.off[k]), int(bins.nb[k])
    bad = ~good[off:off + nb]
    out, j = [], 0
    while j < nb:
        if bad[j]:
            i = j
            while j < nb and bad[j]:
                j += 1
            out.append((int(bins.start[off + i]), int(bins.end[off + j - 1])))
        else:
            j += 1
    return tuple(out)


def visible_gap(a0: int, a1: int, b0: int, b1: int, masked, cap_bp: int | None = None) -> int:
    """bp between two intervals that Hi-C can see: the gap minus the unmappable stretches inside it, of which at most
    cap_bp are discounted (None: all)."""
    g = gap(a0, a1, b0, b1)
    if g == 0 or not masked:
        return g
    lo, hi = min(a1, b1), max(a0, b0)
    m = sum(max(0, min(e, hi) - max(s, lo)) for s, e in masked)
    return max(0, g - (m if cap_bp is None else min(m, cap_bp)))


def _partners(bins, members, good, seq: dict) -> tuple:
    """({chromosome index: midpoint bin of its sequence call}, {chromosome index: partner window bins}) for calls at
    ANCHOR_LEVELS (all calls when fewer than GATE_MIN_REF qualify). Windows without a mappable bin are dropped."""
    idx = {n: i for i, n in enumerate(bins.names)}

    def pick(levels):
        est, win = {}, {}
        for c, (s, e, conf) in seq.items():
            k = idx.get(c)
            if k is not None and members[k] and s >= 0 and conf in levels:
                w = visible_window(bins, good, k, s, e)
                if len(w):
                    est[k] = min(int(bins.nb[k]) - 1, ((s + e) // 2) // bins.size)
                    win[k] = w
        return est, win
    est, win = pick(ANCHOR_LEVELS)
    return (est, win) if len(est) >= GATE_MIN_REF else pick(("high", "medium", "low"))


def _agrees_bin(bins, good, k: int, j: int, s: int, e: int) -> bool:
    """The peak bin j lies within HIC_TOL_BINS mappable bins of [s, e) (unmappable bins free up to HIC_FAR_BINS)."""
    return bool(visible_dist(bins, good, k, s, e, HIC_FAR_BINS)[j] <= HIC_TOL_BINS)


def _chance(S, bins, good, k: int, s: int, e: int) -> float:
    """Share of chromosome k's scored bins whose peak would agree with [s, e) (a uniformly placed peak), in the same
    mappable-bin geometry as the agreement."""
    sl = bins.sl(k)
    ok = np.isfinite(S[sl])
    if not ok.any():
        return float("nan")
    near = visible_dist(bins, good, k, s, e, HIC_FAR_BINS) <= HIC_TOL_BINS
    return float((near & ok).sum() / ok.sum())


SUPPORT_LEVELS = ("medium", "high")      # anchored peaks that count for the support tier: >= 2 of the 3 engine tests
END_ZONE = 0.125                         # end_test: Levan's telocentric bound (fold.CI_MIN)
END_SHARE = 0.05                         # ... the terminal 5% of a chromosome's scored bins stands for an end (a heuristic, about one
END_MIN_BINS = 4                         # partner half-window at the bins chosen); at least the half-window itself
CIS_SIM = 200                            # cis_chance: contact maps simulated per chromosome (resolution 0.005, as the engine's null)


def clean_peak(flags) -> bool:
    """An anchored peak without a hub, rDNA or telomeric flag (those attract trans contacts whatever the centromere)."""
    return not set(str(flags or "").split(",")) & set(NOT_STRONG_FLAGS)


def gate(table: pd.DataFrame) -> pd.DataFrame:
    """Leave-one-chromosome-out gate (see module doc). Reference: chromosomes with a high sequence call (plus medium
    when fewer than GATE_MIN_REF are high) and a scored anchored peak. Hi-C must also agree on a majority of the whole
    reference set: with a few references, leaving out one that disagrees can turn 2 of 4 into 2 of 3 (T2T sheep,
    partners re-paired at random: 2 gate passes; on real Hi-C the rule changed no gate in 14 mammals).

    A genome without any high or medium call (no satellite family dominates, so the satellite engine caps every call:
    siamang without genes, chicken, EquCab3) has no reference-grade call, and Hi-C could never act where the sequence
    line is least sure. There the low calls are the references, and an agreement counts only at a clean peak: chromosome
    ends, hubs and NORs cluster in 3D whatever the centromere, so in a Rabl nucleus wrong calls at the ends would
    otherwise agree with Hi-C. Two independent lines that agree beyond chance on most chromosomes validate each other;
    gate_refs records which references were used. Where some calls are confident but too few for a leave-one-out test (fewer
    than GATE_MIN_REF + 1), the low calls join them on the same terms ('mixed'; future.md #30: cucumber with tes/ had 2 high
    calls, and Hi-C went untested on all 7 chromosomes although it had been tested on the low calls without tes/)."""
    t = table.copy()
    scored = t.anch_peak_start.ge(0) & t.chance.notna() & t.seq_start.ge(0)
    ref = t[scored & t.seq_confidence.eq("high")]
    if len(ref) < GATE_MIN_REF + 1:
        ref = t[scored & t.seq_confidence.isin(["high", "medium"])]
    refs, agree = "reference", ref.anch_agrees_seq
    if not t.seq_confidence.isin(["high", "medium"]).any():
        ref = t[scored & t.seq_confidence.eq("low")]
        refs, agree = "low", ref.anch_agrees_seq & ref.anch_flags.map(clean_peak)
    elif len(ref) < GATE_MIN_REF + 1:
        ref = t[scored & t.seq_confidence.isin(["high", "medium", "low"])]
        clean = ref.anch_flags.map(clean_peak) if "anch_flags" in ref else True
        refs, agree = "mixed", ref.anch_agrees_seq & (ref.seq_confidence.ne("low") | clean)
    majority = 2 * int(agree.sum()) > len(ref)
    sig = agree & ref.anch_level.isin(SUPPORT_LEVELS) if "anch_level" in ref else agree
    ks, ns, ps, passed, support = [], [], [], [], []
    for c in t.chrom:
        r = ref.chrom != c
        k, n, k_sig = int(agree[r].sum()), int(r.sum()), int(sig[r].sum())
        ch = ref[r].chance.to_numpy()
        p = poisson_binomial_sf(k, ch) if n else float("nan")
        ok = bool(n >= GATE_MIN_REF and 2 * k > n and p <= GATE_P and majority)
        ks.append(k)
        ns.append(n)
        ps.append(p)
        passed.append(ok)
        support.append(bool(not ok and n >= GATE_MIN_REF and k_sig >= GATE_MIN_REF
                            and poisson_binomial_sf(k_sig, ch) <= GATE_P))
    t["gate_k"], t["gate_n"], t["gate_p"], t["gate_pass"] = ks, ns, ps, passed
    t["gate_support"] = support
    t["gate_refs"] = refs
    return t


def _end_score(v: np.ndarray) -> tuple:
    """Mean anchored score of the first and of the last END_SHARE of a chromosome's scored bins (at least END_MIN_BINS)."""
    ok = np.flatnonzero(np.isfinite(v))
    if len(ok) < 2 * END_MIN_BINS:
        return float("nan"), float("nan")
    n = max(END_MIN_BINS, int(len(ok) * END_SHARE))
    return float(np.mean(v[ok[:n]])), float(np.mean(v[ok[-n:]]))


def end_test(t: pd.DataFrame, S, bins) -> pd.DataFrame:
    """Which end? (future.md #15) A call in the terminal eighth of a chromosome (Levan's telocentric bound, fold.CI_MIN) sits where a
    peak cannot count: chromosome ends cluster in 3D whatever the centromere, so such a peak is flagged telomeric. But the two ends of a
    chromosome can be told apart: the anchored score, made from the other chromosomes' calls only, is about twice as high at the end
    holding the call in goat (28 of 30), sheep (22/23), horse (18/19) and mouse (19/19). Per chromosome x: on the OTHER chromosomes
    with such a call, is the call's end the higher one, on most of them, beyond the 1/2 of chance (binomial, P <= GATE_P)?
    Chicken does not pass (12 of 19)."""
    t = t.copy()
    names = {n: i for i, n in enumerate(bins.names)}
    near, other, higher = [], [], []
    for r in t.itertuples(index=False):
        a, b, h = np.nan, np.nan, False
        k = names[r.chrom]
        L = int(bins.lengths[k])
        if r.seq_start >= 0 and r.callable:
            mid = (r.seq_start + r.seq_end) / 2
            if min(mid, L - mid) < END_ZONE * L:
                lo, hi = _end_score(S[bins.sl(k)])
                a, b = (lo, hi) if mid < L / 2 else (hi, lo)
                h = bool(a > b)
        near.append(a), other.append(b), higher.append(h)
    t["end_near"], t["end_other"], t["end_higher"] = near, other, higher
    tested = np.isfinite(t.end_near.to_numpy(dtype=float))
    ks, ns, ps, passed = [], [], [], []
    hi = t.end_higher.to_numpy(dtype=bool)
    for i in range(len(t)):
        r_ = tested.copy()
        r_[i] = False
        n, k = int(r_.sum()), int(hi[r_].sum())
        p = binom_tail(k, n) if n else float("nan")
        ks.append(k), ns.append(n), ps.append(p), passed.append(bool(n >= GATE_MIN_REF and 2 * k > n and p <= GATE_P))
    t["end_k"], t["end_n"], t["end_p"], t["end_pass"] = ks, ns, ps, passed
    return t


def binom_tail(k: int, n: int) -> float:
    """P(X >= k), X ~ Binomial(n, 1/2)."""
    return poisson_binomial_sf(k, [0.5] * n)


def anchored(state: dict, seq: dict, threads: int) -> tuple[pd.DataFrame, dict]:
    """state: the Hi-C engine's in-memory state; seq: {chrom: (start, end, confidence)} of the sequence calls (start -1 =
    none). Returns (per-chromosome table without the de novo columns, {'S': anchored score, 'est': partners})."""
    M, W, bins, good, mean, members = (state[k] for k in ("M", "W", "bins", "good", "mean", "members"))
    est, win = _partners(bins, members, good, seq)
    if len(est) < 2:
        log.warning("Hi-C: fewer than 2 chromosomes with a sequence call to anchor partner windows; anchored test "
                    "skipped")
    oe = partner_oe(W, bins, good, mean, win) if len(est) >= 2 else {}
    S = cluster_score(oe, bins, good, members) if oe else np.full(bins.N, np.nan)
    pcs = {}
    for k in (int(x) for x in np.where(members)[0]):
        pc = peak_core(S, bins, k)
        if pc is not None and np.isfinite(pc["ratio"]):
            pcs[k] = pc
    pre = state.get("boot_extra", {}).get("anchored")          # computed by the engine on its resamples (run_denovo: seq)
    same = pre is not None and _same_windows(state.get("extra_windows", {}).get("anchored"), win)
    boot = (pre if same else bootstrap_all(M, bins, members, win, threads)) if pcs else {}
    lopo = lopo_all(oe, bins, good, members) if pcs else {}
    pnull = null_all(W, bins, good, mean, members, est, {k: v["peak"] for k, v in pcs.items()},
                     {k: v["ratio"] for k, v in pcs.items()}, threads=threads) if pcs else {}
    rows = []
    for k, name in enumerate(bins.names):
        s, e, conf = seq.get(name, (-1, -1, "none"))
        r = dict(chrom=name, callable=bool(members[k]), bin_size=bins.size, seq_start=s, seq_end=e,
                 seq_confidence=conf, partner=k in est, anch_peak_start=-1, anch_peak_end=-1, anch_core_start=-1,
                 anch_core_end=-1, anch_ratio=np.nan, anch_boot_in_core=np.nan, anch_lopo_stable=False,
                 anch_null_p=np.nan, anch_level="none", anch_flags="", anch_strong=False, anch_agrees_seq=False,
                 chance=np.nan, masked=masked_runs(bins, good, k))
        if k in pcs:
            pc, g = pcs[k], lambda j, k=k: int(bins.off[k] + j)
            b = boot.get(k, np.array([], int))
            inside = float(np.mean((b >= pc["core_lo"]) & (b <= pc["core_hi"]))) if len(b) else 0.0
            lp = lopo.get(k, np.array([], int))
            stable = bool(len(lp) and ((lp >= pc["core_lo"]) & (lp <= pc["core_hi"])).all())
            lev = level(inside, stable, pnull.get(k, np.nan))
            flags = call_flags(k, pc, bins, state["hub"], state["rdna"], True, pnull.get(k, np.nan), True, good=good)
            r.update(anch_peak_start=int(bins.start[g(pc["peak"])]), anch_peak_end=int(bins.end[g(pc["peak"])]),
                     anch_core_start=int(bins.start[g(pc["core_lo"])]), anch_core_end=int(bins.end[g(pc["core_hi"])]),
                     anch_ratio=round(pc["ratio"], 3), anch_boot_in_core=round(inside, 3), anch_lopo_stable=stable,
                     anch_null_p=pnull.get(k, np.nan), anch_level=lev, anch_flags=",".join(flags),
                     anch_strong=bool(lev == STRONG_HIC_LEVEL and not set(flags) & set(NOT_STRONG_FLAGS)))
            if s >= 0:
                r.update(anch_agrees_seq=_agrees_bin(bins, good, k, pc["peak"], s, e),
                         chance=_chance(S, bins, good, k, s, e))
        rows.append(r)
    t = end_test(gate(pd.DataFrame(rows)), S, bins)
    return t, dict(S=S, est=est, win=win)


def _core_zone(bins, k: int, s: int, e: int) -> tuple:
    """(first, last) local bin of an LTR core and HALF_WIN bins on each side: where a peak counts as on the core."""
    lo, hi = _bin_range(bins, k, s, e)
    return max(0, lo - hp.HALF_WIN), min(int(bins.nb[k]) - 1, hi + hp.HALF_WIN)


def anchored_on_cores(state: dict, cores: dict, label: str = "the LTR cores") -> dict:
    """The anchored test with LTR cores in place of the sequence calls (future.md #11: where the sequence line is wrong
    genome-wide, Hi-C anchored on it finds nothing). Partner windows at the cores. A chromosome agrees when its peak, the top of a
    score made from the other chromosomes' anchors only, lies in its own core's zone; the zone's share of the chromosome is the
    chance. The anchors stay where the LTR line put them: moved to the best score within their zones they find each other in any
    map (randomized soybean Hi-C: 1 of 20 peaks on the cores, 4 of 20 moved, P 0.002). A peak counts only where it is clean
    (clean_peak, as for the gate's low references): chromosome ends, hubs and NORs attract trans contacts whatever the centromere,
    so cores at the chromosome ends agree with any nucleus that clusters its telomeres (held-out R. breviuscula, holocentric: Athila
    cores at the ends of 4 of 5 chromosomes, 3 agreed, P 5e-4). Returns dict(k, n, p, status, peaks {chrom: (start, end)} of the
    agreeing chromosomes, chance {chrom: the zone's share}): pass (on most chromosomes, Poisson-binomial P <= GATE_P); support
    (beyond chance, on fewer: the rule of the support tier); fail. cis_trans_line reads the same test on other anchors (label)."""
    W, bins, good, mean, members = (state[k] for k in ("W", "bins", "good", "mean", "members"))
    hub = state.get("hub")
    hub = np.zeros(bins.N, dtype=bool) if hub is None else hub
    rdna = state.get("rdna") or {}
    idx = {n: i for i, n in enumerate(bins.names)}
    win = {idx[c]: visible_window(bins, good, idx[c], *iv) for c, iv in cores.items() if c in idx and members[idx[c]]}
    S = cluster_score(partner_oe(W, bins, good, mean, {k: w for k, w in win.items() if len(w)}), bins, good, members)
    peaks, chance = {}, {}
    for c, iv in cores.items():
        k = idx.get(c)
        if k is None or not members[k]:
            continue
        v = S[bins.sl(k)]
        if not np.isfinite(v).any():
            continue
        lo, hi = _core_zone(bins, k, *iv)
        j = int(np.nanargmax(v))
        chance[c] = float(np.isfinite(v[lo:hi + 1]).sum() / np.isfinite(v).sum())
        fl = call_flags(k, dict(peak=j, core_lo=j, core_hi=j), bins, hub, rdna, True, 0.0, True, good=good)
        if lo <= j <= hi and clean_peak(",".join(fl)):
            g = int(bins.off[k]) + j
            peaks[c] = (int(bins.start[g]), int(bins.end[g]))
    n, m = len(chance), len(peaks)
    p = poisson_binomial_sf(m, list(chance.values())) if n else 1.0
    status = "fail" if p > GATE_P or m < GATE_MIN_REF else "pass" if 2 * m > n else "support"
    (log.info if label == "the LTR cores" else log.debug)(        # the cis signs' raw tests: cis_trans_line logs the verdict
        f"Hi-C anchored on {label}: the peak lies at its anchor on {m}/{n} chromosomes (P {p:.2g}): {status}")
    return dict(k=m, n=n, p=p, status=status, peaks=peaks, chance=chance)


def _crossing(diagonals: dict, n: int, half: int) -> np.ndarray:
    """Per bin of a chromosome of n bins: the sum, over the half bins on its left and the half bins on its right, of a value
    per pair of bins. diagonals: {distance d: value of each pair (a, a + d)}; such a pair lies across bin a + t for t from
    max(1, d - half + 1) to min(half, d)."""
    v = np.zeros(n)
    for d, x in diagonals.items():
        for t in range(max(1, d - half + 1), min(half, d) + 1):
            v[t:n - d + t] += x
    return v


def _cis_parts(W, bins, good, k: int, half: int) -> tuple:
    """Of chromosome k, per distance d within the crossing box: (balanced contact of each pair of bins, whether both are
    mappable, the mean contact over mappable pairs), and per bin what its box would hold at those means."""
    sl, n = bins.sl(k), int(bins.nb[k])
    B, g = W[sl, :][:, sl], good[sl]
    diag, ok, mean = {}, {}, {}
    for d in range(1, 2 * half):
        diag[d] = B.diagonal(d)
        ok[d] = g[:n - d] & g[d:]
        mean[d] = float(diag[d][ok[d]].mean()) if ok[d].any() else 0.0
    expected = _crossing({d: ok[d] * mean[d] for d in ok}, n, half)
    expected[:half], expected[n - half:] = 0.0, 0.0                 # a bin needs its whole box on the chromosome
    return diag, ok, mean, expected


def cis_profile(W, bins, good, k: int, half: int = hp.HALF_WIN) -> np.ndarray:
    """Per bin of chromosome k: the balanced contact between the half bins on its left and the half bins on its right (how much
    contact crosses it within the partner half-window) over what the chromosome's mean contact at those distances predicts; NaN
    where either side has no mappable bin. Beside a masked bin the box holds an unusual mix of distances, and contact falls
    steeply with distance, so a plain mean has its extremes there whatever the locus. W: CSR."""
    n = int(bins.nb[k])
    if n <= 2 * half:
        return np.full(n, np.nan)
    diag, _, _, expected = _cis_parts(W, bins, good, k, half)
    v = np.full(n, np.nan)
    use = expected > 0
    v[use] = _crossing(diag, n, half)[use] / expected[use]
    return v


def cis_chance(W, bins, good, bias, k: int, lo: int, hi: int, half: int = hp.HALF_WIN, n: int = CIS_SIM) -> dict:
    """{+1: chance that the highest, -1: chance that the lowest value of cis_profile lies in bins lo..hi of chromosome k, when
    no locus is special}: profiles from n contact maps drawn (Poisson) given each bin's coverage and the chromosome's mean
    contact at each distance. A stretch's share of the chromosome is not that chance, since a profile is noisiest where coverage
    is low, and LTR cores lie there: randomized Brachypodium Hi-C has extremes at the cores on 3 of 5 chromosomes, zones
    7% of the chromosome and 55% by simulation. (1 + hits) / (1 + n), as the engine's null."""
    nb = int(bins.nb[k])
    hits = {+1: 0, -1: 0}
    if nb <= 2 * half:
        return {x: 1.0 for x in hits}
    b = np.asarray(bias)[bins.sl(k)]
    _, ok, mean, expected = _cis_parts(W, bins, good, k, half)
    use = expected > 0
    if not use.any():
        return {x: 1.0 for x in hits}
    pos = np.flatnonzero(use)
    weight, counts = {}, {}
    for d in ok:
        with np.errstate(divide="ignore", invalid="ignore"):
            weight[d] = np.where(ok[d], 1.0 / (b[:nb - d] * b[d:]), 0.0)    # balanced = count * weight
            counts[d] = np.where(weight[d] > 0, mean[d] / weight[d], 0.0)
    rng = np.random.default_rng([hp.SEED, k])
    for _ in range(n):
        v = _crossing({d: rng.poisson(counts[d]) * weight[d] for d in ok}, nb, half)[use] / expected[use]
        hits[+1] += bool(lo <= pos[int(np.argmax(v))] <= hi)
        hits[-1] += bool(lo <= pos[int(np.argmin(v))] <= hi)
    return {x: (1 + c) / (1 + n) for x, c in hits.items()}


def cis_line(state: dict, cores: dict, threads: int = 1) -> dict:
    """A cis readout of the centromere, tested against the LTR cores (future.md #12). The contact crossing a bin (cis_profile) has
    its extreme at the centromere in many genomes, but which extreme differs: a peak in tomato, rice and Brachypodium, a dip in
    potato. So the LTR line chooses the sign: for chromosome x, the sign whose extreme lies at the cores of the OTHER chromosomes most
    significantly. Genome-wide, the better sign's agreement with the cores (the zone of _core_zone) must beat chance, Bonferroni
    over the two signs; the chance per chromosome is the larger of the zone's share and cis_chance. Returns dict(k, n, p, sign,
    status, peaks {chrom: (start, end)}, at_core {chrom: bool}): pass (on most chromosomes), support (beyond chance, on fewer),
    fail."""
    W, bins, good, members, bias = (state[k] for k in ("W", "bins", "good", "members", "bias"))
    W = W.tocsr()
    idx = {n: i for i, n in enumerate(bins.names)}
    todo = []
    for c, iv in cores.items():
        k = idx.get(c)
        if k is None or not members[k]:
            continue
        v = cis_profile(W, bins, good, k)
        if not np.isfinite(v).any():
            continue
        lo, hi = _core_zone(bins, k, *iv)
        ok = np.isfinite(v)
        todo.append((c, k, v, lo, hi, float(ok[lo:hi + 1].sum() / ok.sum())))
    with ThreadPoolExecutor(max(threads, 1)) as ex:      # each chromosome draws from its own stream: any order gives the same
        sims = list(ex.map(lambda x: cis_chance(W, bins, good, bias, x[1], x[3], x[4]), todo))
    best, chance = {}, {}
    for (c, k, v, lo, hi, share), sim in zip(todo, sims):
        chance[c] = {x: max(share, sim[x]) for x in sim}
        best[c] = {+1: int(np.nanargmax(v)), -1: int(np.nanargmin(v))}
        best[c] = {g: (j, lo <= j <= hi) for g, j in best[c].items()}

    def test(sign, leave=None):
        names = [c for c in best if c != leave]
        m = sum(best[c][sign][1] for c in names)
        return m, len(names), poisson_binomial_sf(m, [chance[c][sign] for c in names]) if names else 1.0
    if not best:
        return dict(k=0, n=0, p=1.0, sign=0, status="fail", peaks={}, at_core={})
    sign = min((+1, -1), key=lambda g: (test(g)[2], -g))
    m, n, p = test(sign)
    p = min(1.0, 2 * p)
    status = "fail" if p > GATE_P or m < GATE_MIN_REF else "pass" if 2 * m > n else "support"
    peaks, at_core = {}, {}
    for c in best:
        g = min((+1, -1), key=lambda x: (test(x, c)[2], -x))          # the sign the other chromosomes choose
        j, on = best[c][g]
        off = int(bins.off[idx[c]]) + j
        peaks[c], at_core[c] = (int(bins.start[off]), int(bins.end[off])), bool(on)
    log.info(f"Hi-C cis signal against the LTR cores: its {'peak' if sign > 0 else 'dip'} lies at the core on {m}/{n} chromosomes "
             f"(P {p:.2g}): {status}")
    return dict(k=m, n=n, p=p, sign=sign, status=status, peaks=peaks, at_core=at_core, tested_by="ltr", interior=False)


def interior(v: np.ndarray) -> np.ndarray:
    """A per-bin profile with its terminal eighths (END_ZONE, Levan's telocentric bound) set to NaN: chromosome ends cluster in
    3D whatever the centromere (the rule of the which-end test and of the telomeric flag), so a contact readout cannot test a
    place there."""
    v = np.array(v, dtype=float)
    e = int(np.ceil(END_ZONE * len(v)))
    if e:
        v[:e], v[len(v) - e:] = np.nan, np.nan
    return v


def cis_trans_line(state: dict) -> dict:
    """The cis signal tested by Hi-C's own trans contacts (future.md #26), for a genome without LTR cores: two readouts of one
    library, the contact crossing each bin and the contact between chromosomes, test each other where the gate fails. For each sign,
    each chromosome's anchor is the extreme of its cis profile (cis_profile) outside its terminal eighths (interior), and the
    anchored test (anchored_on_cores) asks whether each chromosome's trans peak, made from the other chromosomes' anchors only,
    lies at its own anchor; the zone's share is the chance, Bonferroni over the two signs. Tomato without tes/: the peak at its
    anchor on 10 of 12 chromosomes, P 1.2e-13 for the two signs, all 10 on CENH3; C. neoformans in G1, 12 of 14. The terminal eighths are out because
    there both readouts see the chromosome ends, which fold and cluster whatever the centromere (Ae. tauschii and potato: every
    agreement of an unrestricted profile lay in a terminal eighth, none on CENH3).
    The two readouts share the library's coverage, so agreement on a minority of chromosomes can come from coverage alone
    (randomized soybean Hi-C: 6 of 20 dips, P 3e-5; no randomized library reached a majority): the line acts only at pass, on a
    majority beyond chance, and only on a chromosome whose OTHER chromosomes pass that test too (leave one out, as the gate). Both
    signs passing leaves the sign unknown: no action. Returns cis_line's fields (peaks: the chromosomes where it acts) with
    tested_by 'trans', interior True, pairs {chrom: its own trans peak lies at its anchor} and status pass / support / fail /
    ambiguous."""
    W, bins, good, members = (state[k] for k in ("W", "bins", "good", "members"))
    W = W.tocsr()
    ext = {+1: {}, -1: {}}
    for k, c in enumerate(bins.names):
        if not members[k]:
            continue
        v = interior(cis_profile(W, bins, good, k))
        if not np.isfinite(v).any():
            continue
        off = int(bins.off[k])
        for g in (+1, -1):
            j = int(np.nanargmax(g * v))
            ext[g][c] = (int(bins.start[off + j]), int(bins.end[off + j]))
    empty = dict(k=0, n=0, p=1.0, sign=0, status="fail", peaks={}, at_core={}, pairs={}, tested_by="trans", interior=True)
    if len(ext[+1]) < GATE_MIN_REF:
        return empty
    st = dict(state, W=W)
    res = {g: anchored_on_cores(st, ext[g], label=f"the cis {'peaks' if g > 0 else 'dips'}") for g in (+1, -1)}

    def test(g, leave=None):
        r = res[g]
        names = [c for c in r["chance"] if c != leave]
        m = sum(c in r["peaks"] for c in names)
        p = min(1.0, 2 * poisson_binomial_sf(m, [r["chance"][c] for c in names])) if names else 1.0
        return m, len(names), p

    def verdict(m, n, p):
        return "fail" if p > GATE_P or m < GATE_MIN_REF or n < GATE_MIN_REF else "pass" if 2 * m > n else "support"
    sign = min((+1, -1), key=lambda g: (test(g)[2], -g))
    m, n, p = test(sign)
    status = verdict(m, n, p)
    if status != "fail" and verdict(*test(-sign)) == "pass":
        status = "ambiguous"                      # both signs pass: which extreme is the centromere is unknown
    peaks, pairs = {}, {}
    for c in res[sign]["chance"]:
        pairs[c] = c in res[sign]["peaks"]
        if status == "pass" and verdict(*test(sign, c)) == "pass":
            peaks[c] = ext[sign][c]
    log.info(f"Hi-C cis signal tested by trans contacts: trans peaks lie at the cis {'peak' if sign > 0 else 'dip'} on {m}/{n} "
             f"chromosomes (P {p:.2g}): {status}; acts on {len(peaks)}")
    return dict(k=m, n=n, p=p, sign=sign, status=status, peaks=peaks, at_core=dict(pairs), pairs=pairs, tested_by="trans",
                interior=True)


def cis_acts(cis: dict | None) -> bool:
    """Whether a cis line's peaks may act on calls: the LTR-tested line at pass or support; the trans-tested one at pass only
    (cis_trans_line keeps in its peaks only the chromosomes where it acts)."""
    if not cis:
        return False
    if cis.get("tested_by") == "trans":
        return cis.get("status") == "pass"
    return cis.get("status") in ("pass", "support")


def add_denovo(t: pd.DataFrame, rows: list, bins) -> pd.DataFrame:
    """The prior-free Hi-C calls next to the anchored ones (reported; never used to decide)."""
    by = {r["chrom"]: r for r in rows}
    idx = {n: i for i, n in enumerate(bins.names)}
    cols = {c: [] for c in HIC_COLS if c.startswith("denovo_")}
    for r in t.itertuples(index=False):
        d = by.get(r.chrom, {})
        ok = d.get("confidence", "none") != "none" and d.get("peak_start", "") != ""
        ps, pe = (int(d["peak_start"]), int(d["peak_end"])) if ok else (-1, -1)
        cs, ce = (int(d["core_start"]), int(d["core_end"])) if ok else (-1, -1)
        cols["denovo_peak_start"].append(ps)
        cols["denovo_peak_end"].append(pe)
        cols["denovo_core_start"].append(cs)
        cols["denovo_core_end"].append(ce)
        cols["denovo_confidence"].append(d.get("confidence", "none"))
        cols["denovo_flags"].append(d.get("flags", ""))
        tol = HIC_TOL_BINS * bins.size
        cols["denovo_agrees_seq"].append(bool(ok and r.seq_start >= 0 and gap(ps, pe, r.seq_start, r.seq_end) <= tol))
        cols["denovo_agrees_anch"].append(bool(ok and r.anch_peak_start >= 0 and
                                               gap(ps, pe, r.anch_peak_start, r.anch_peak_end) <= tol))
    for c, v in cols.items():
        t[c] = v
    return t[HIC_COLS]


def add_nors(rdna: dict, nors) -> dict:
    """The Hi-C engine's rDNA intervals {chrom: (starts, ends)} (rDNA-tagged satellite arrays) plus the 45S arrays."""
    out = {c: (np.asarray(a), np.asarray(b)) for c, (a, b) in (rdna or {}).items()}
    for c, st, en, *_ in nors or ():
        a, b = out.get(c, (np.zeros(0, np.int64), np.zeros(0, np.int64)))
        out[c] = (np.append(a, st), np.append(b, en))
    return out


def nor_peaks(t: pd.DataFrame, nors) -> pd.DataFrame:
    """Anchored peaks that could oppose or place a call (strong, or clean medium; on a chromosome without a sequence
    call or more than HIC_FAR_BINS mappable bins from it) within HIC_TOL_BINS mappable bins of a 45S array (rdna.py)
    get the rDNA flag and are not strong. NORs contact one another across chromosomes (the nucleolus), so they mimic
    a centromere cluster; their arrays are unmappable, so the peak sits beside them (CHM13 chr15, chr22: p-arm peaks
    contested correct calls). Peaks whose core overlaps a NOR are flagged by the Hi-C engine itself (add_nors)."""
    by = {}
    for c, st, en, *_ in nors or ():
        by.setdefault(c, []).append((st, en))
    if not by or not len(t):
        return t
    t = t.copy()
    hit = []
    for i, r in t.iterrows():
        if r.anch_peak_start < 0 or r.chrom not in by or "rDNA" in str(r.anch_flags).split(","):
            continue
        fl = set(str(r.anch_flags).split(",")) - {""}
        if not (r.anch_strong or (r.anch_level == "medium" and not fl & set(NOT_STRONG_FLAGS))):
            continue
        cap = HIC_FAR_BINS * r.bin_size
        if r.seq_start >= 0 and visible_gap(r.anch_peak_start, r.anch_peak_end, r.seq_start, r.seq_end,
                                             r.masked) <= cap:                    # cannot oppose: nothing to test
            continue
        if any(visible_gap(r.anch_peak_start, r.anch_peak_end, a, b, r.masked, cap) <= HIC_TOL_BINS * r.bin_size
               for a, b in by[r.chrom]):
            t.at[i, "anch_flags"] = ",".join([f for f in str(r.anch_flags).split(",") if f] + ["rDNA"])
            t.at[i, "anch_strong"] = False
            hit.append(r.chrom)
    if hit:
        log.info(f"Hi-C: anchored peaks beside 45S arrays on {', '.join(hit)}: flagged rDNA (NOR), never oppose or "
                 f"place a call")
    return t


def run(pairs, fai, prefix: Path, satellites, seq: dict, threads: int, classes: dict | None = None,
        nors=None, te_cores: dict | None = None) -> dict:
    """Both Hi-C analyses for one genome, and the 3D association of each feature class. Returns dict(table,
    S (anchored score), denovo_S, bins, members, hub, rabl, stats, bin_size, telomere_S, classes, coverage (contacts
    per bin, for the karyotype test of provenance), base_size, matrix_cache (the engine's finest matrix, for the
    arm fold; None for .cool input), te_anchored and cis (the two tests against the LTR cores, or None)). nors: the genome's 45S
    arrays (rdna.detect), flagged rDNA. te_cores: {chrom: (start, end)} of the LTR line: where the gate fails on most chromosomes,
    the anchored test runs again on them (anchored_on_cores) and the cis signal is tested against them (cis_line); without
    LTR cores the cis signal is tested by the trans contacts (cis_trans_line)."""
    res = run_denovo(pairs, fai, prefix, satellites, threads, seq=seq)
    st = res["state"]
    st["rdna"] = add_nors(st["rdna"], nors)
    t, extra = anchored(st, seq, threads)
    te_anchored, cis = None, None
    if te_cores and 2 * int(t.gate_pass.sum()) <= len(t):
        te_anchored = anchored_on_cores(st, te_cores)
        cis = cis_line(st, te_cores, threads)
    elif 2 * int(t.gate_pass.sum()) <= len(t):          # no LTR cores: the cis signal is tested by the trans contacts instead
        cis = cis_trans_line(st)
    t = nor_peaks(t, nors)
    t = add_denovo(t, res["rows"], st["bins"])
    fc = feature_contacts(st, classes) if classes else pd.DataFrame(columns=FC_COLS)
    if len(fc):
        sig = fc[fc.p_enriched <= FC_LEVEL]
        log.info(f"Hi-C feature classes: {len(fc)} tested; self-associated (P <= {FC_LEVEL}): "
                 f"{', '.join(f'{r.feature_class} {r.fold:.2f}x' for r in sig.itertuples()) or 'none'}")
    n_pass = int(t.gate_pass.sum())
    log.info(f"Hi-C {'+'.join(Path(f).name for f in hic_cli.pair_files(pairs))}: bin {st['size']} bp; anchored "
             f"peaks on the sequence call "
             f"{int(t.anch_agrees_seq.sum())}/{int(t.seq_start.ge(0).sum())}; de novo "
             f"{int(t.denovo_agrees_seq.sum())}; gate passes on {n_pass}/{len(t)} chromosomes")
    return dict(table=t, S=extra["S"], est=extra["est"], denovo_S=st["loc"].S, bins=st["bins"],
                members=st["members"], hub=st["hub"], rabl=res["rabl"], stats=res["stats"], bin_size=st["size"],
                telomere_S=st["T"], classes=fc, coverage=np.asarray(st["M"].sum(axis=1)).ravel(),
                base_size=st["base_size"], matrix_cache=st["matrix_cache"], te_anchored=te_anchored, cis=cis)


# ---------------------------------------------------------------------------------------------------------------
# Feature classes in 3D: does a class of features (the sequence centromeres, chromosome ends, telomere arrays, each
# satellite family's blocks, each concentrated LTR clade's cores) contact itself across chromosomes more than the same
# number of places drawn at random along the same chromosomes? Centromeres are the positive control in Rabl nuclei.
FC_PERM = 1000             # permutations: P resolution 0.001
FC_MIN = 3                 # instances needed, on >= 2 chromosomes (fewer give no trans pair to average)
FC_MAX_INST = 400          # larger classes are tested on a fixed random sample of this many (cost; the mean over
                           # ~80,000 trans pairs is already precise)
FC_MAX_BINS = 15_000       # dense genome-wide O/E above this many bins is too large; the test is skipped
FC_COLS = ["feature_class", "kind", "label", "instances", "chromosomes", "trans_pairs", "observed_oe",
           "null_median_oe", "fold", "p_enriched"]


def oe_dense(W, bins, good, mean, members) -> np.ndarray | None:
    """Trans observed/expected for every bin pair (W / its chromosome-pair block mean), box-smoothed over 3x3 bins
    within each chromosome-pair block (one bin of placement error); cis pairs and masked bins NaN. float32."""
    N = bins.N
    if N > FC_MAX_BINS:
        return None
    from scipy.ndimage import uniform_filter
    D = W.toarray().astype(np.float32)
    m = mean[bins.chrom_of][:, bins.chrom_of].astype(np.float32)
    with np.errstate(invalid="ignore", divide="ignore"):
        D = np.where(m > 0, D / m, np.nan)
    ok_bin = good & members[bins.chrom_of]
    D[~ok_bin, :] = np.nan
    D[:, ~ok_bin] = np.nan
    S = np.full_like(D, np.nan)
    ks = [int(k) for k in np.where(members)[0]]
    for a in ks:                       # smooth inside each chromosome-pair block: never across a chromosome boundary
        for b in ks:
            if a == b:
                continue
            blk = D[bins.sl(a), bins.sl(b)]
            fin = np.isfinite(blk)
            num = uniform_filter(np.where(fin, blk, 0).astype(np.float32), size=3, mode="constant")
            den = uniform_filter(fin.astype(np.float32), size=3, mode="constant")
            with np.errstate(invalid="ignore", divide="ignore"):
                S[bins.sl(a), bins.sl(b)] = np.where(den > 0, num / den, np.nan)
    S[~ok_bin, :] = np.nan
    S[:, ~ok_bin] = np.nan
    return S


def _rep_bin(bins, ok_bin, k: int, s: int, e: int):
    """Global bin standing for [s, e) on chromosome k: its midpoint bin, or the nearest usable bin inside it."""
    lo, hi = s // bins.size, max(s, e - 1) // bins.size
    mid = ((s + e) // 2) // bins.size
    cand = sorted(range(lo, hi + 1), key=lambda j: abs(j - mid))
    for j in cand:
        g = int(bins.off[k] + min(j, int(bins.nb[k]) - 1))
        if ok_bin[g]:
            return g
    return None


def feature_contacts(state: dict, classes: dict, n_perm: int = FC_PERM, seed: int = hp.SEED) -> pd.DataFrame:
    """classes: {name: dict(kind, label, iv=[(chrom, start, end)])}. One row per class with enough instances."""
    bins, good, members = state["bins"], state["good"], state["members"]
    S = oe_dense(state["W"], bins, good, state["mean"], members)
    if S is None:
        log.warning(f"Hi-C feature classes: {bins.N} bins is above {FC_MAX_BINS}; 3D association not tested")
        return pd.DataFrame(columns=FC_COLS)
    idx = {n: i for i, n in enumerate(bins.names)}
    ok_bin = good & members[bins.chrom_of] & np.isfinite(np.nanmax(np.where(np.isfinite(S), S, -np.inf), axis=1))
    pool = {k: np.flatnonzero(ok_bin[bins.sl(k)]) + int(bins.off[k]) for k in range(len(bins.names))}
    rng = np.random.default_rng(seed)
    rows = []
    for name, cl in classes.items():
        reps, chs = [], []
        if cl["kind"] == "ends":                         # the outermost usable bin at each end of each chromosome
            for k in np.where(members)[0]:
                p_ = pool[int(k)]
                if len(p_) >= 2:
                    reps += [int(p_[0]), int(p_[-1])]
                    chs += [int(k), int(k)]
        for c, s, e in cl["iv"]:
            k = idx.get(c)
            if k is None or not members[k]:
                continue
            g = _rep_bin(bins, ok_bin, k, int(s), int(e))
            if g is not None:
                reps.append(g)
                chs.append(k)
        reps, chs = np.array(reps, int), np.array(chs, int)
        if len(reps) < FC_MIN or len(set(chs)) < 2:
            continue
        n_all = len(reps)
        if n_all > FC_MAX_INST:
            keep = np.sort(np.random.default_rng(seed).choice(n_all, FC_MAX_INST, replace=False))
            reps, chs = reps[keep], chs[keep]
        trans = chs[:, None] != chs[None, :]
        iu = np.triu(trans, 1)

        def stat(r):
            v = S[np.ix_(r, r)][iu]
            return float(np.nanmean(v)) if np.isfinite(v).any() else np.nan
        obs = stat(reps)
        null = np.empty(n_perm)
        groups = [(np.flatnonzero(chs == k), pool[k]) for k in np.unique(chs)]
        r = np.empty(len(chs), int)
        for i in range(n_perm):
            for at, pk in groups:
                r[at] = pk[rng.integers(len(pk), size=len(at))]
            null[i] = stat(r)
        null = null[np.isfinite(null)]
        med = float(np.median(null)) if len(null) else np.nan
        p = (1 + int((null >= obs).sum())) / (1 + len(null)) if len(null) and np.isfinite(obs) else np.nan
        rows.append(dict(feature_class=name, kind=cl["kind"], label=cl["label"], instances=n_all,
                         chromosomes=len(set(chs)), trans_pairs=int(iu.sum()), observed_oe=round(obs, 4),
                         null_median_oe=round(med, 4), fold=round(obs / med, 3) if med and med > 0 else np.nan,
                         p_enriched=p))
    return pd.DataFrame(rows, columns=FC_COLS)
