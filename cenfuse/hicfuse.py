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
   Hi-C allowed to change x's call.

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
from .hic.balance import balanced, ice
from .hic.confidence import level, poisson_resample
from .hic.io import InputError
from .hic.locate import peak_core
from .hic.qc import call_flags
from .hic.score import block_means, cluster_score, partner_oe, window_bins
from .params import ANCHOR_LEVELS, GATE_MIN_REF, GATE_P, HIC_FAR_BINS, HIC_TOL_BINS, NOT_STRONG_FLAGS, \
    STRONG_HIC_LEVEL
from .util import CenfuseError, gap, log, poisson_binomial_sf

HIC_COLS = ["chrom", "callable", "bin_size", "seq_start", "seq_end", "seq_confidence", "partner",
            "anch_peak_start", "anch_peak_end", "anch_core_start", "anch_core_end", "anch_ratio", "anch_boot_in_core",
            "anch_lopo_stable", "anch_null_p", "anch_level", "anch_flags", "anch_strong", "anch_agrees_seq", "chance",
            "gate_k", "gate_n", "gate_p", "gate_pass", "gate_refs", "masked",
            "denovo_peak_start", "denovo_peak_end", "denovo_core_start", "denovo_core_end", "denovo_confidence",
            "denovo_flags", "denovo_agrees_seq", "denovo_agrees_anch"]


def run_denovo(pairs, fai, prefix, satellites, threads: int) -> dict:
    """The Hi-C engine's `call`, as its own command runs it; returns its rows and in-memory state."""
    args = Namespace(pairs=str(pairs), fai=str(fai), o=str(prefix), repeats=None,
                     satellites=[str(s) for s in satellites] or None, t=max(1, threads), verbose=False)
    try:
        return hic_cli.run_call(args)
    except InputError as e:
        raise CenfuseError(f"Hi-C {pairs}: {e}") from None


def _score_windows(M, bins, members, win: dict):
    """The Hi-C engine's _score_fixed with explicit partner windows (each resample re-balanced; bins it masks drop
    out of the windows)."""
    bias, good, _, _ = ice(M)
    W = balanced(M, bias, good)
    mean = block_means(W, bins, good, members)
    return cluster_score(partner_oe(W, bins, good, mean, {y: w[good[w]] for y, w in win.items()}), bins, good, members)


def bootstrap_all(M, bins, members, win: dict, threads: int, n: int = hp.N_BOOT, seed: int = hp.SEED) -> dict:
    """The Hi-C engine's bootstrap (Poisson resamples, re-balanced, partner windows fixed), but the peak of every
    callable chromosome, not only of those holding a window."""
    ks = [int(k) for k in np.where(members)[0]]

    def one(ss):
        S = _score_windows(poisson_resample(M, np.random.default_rng(ss)), bins, members, win)
        return {k: int(np.nanargmax(S[bins.sl(k)])) if np.isfinite(S[bins.sl(k)]).any() else -1 for k in ks}

    with ThreadPoolExecutor(max(threads, 1)) as ex:
        reps = list(ex.map(one, np.random.SeedSequence(seed).spawn(n)))
    return {k: np.array([r[k] for r in reps]) for k in ks}


def lopo_all(oe, bins, good, members) -> dict:
    """Peak of every callable chromosome with each partner window left out in turn (the Hi-C engine's leave_one_out).
    The same arithmetic as cluster_score(exclude={x: {y}}) restricted to x's own bins: scoring every chromosome for
    every (x, y) pair is cubic in the partners, and a 25-partner genome took an hour."""
    out = {}
    for x in (int(k) for k in np.where(members)[0]):
        sl = bins.sl(x)
        g = good[sl]
        ys = [y for y in oe if y != x]
        peaks = []
        for y in ys:
            rest = [z for z in ys if z != y]
            if not rest:
                peaks.append(-1)
                continue
            A = np.vstack([oe[z][sl] for z in rest])
            ok = ~np.isnan(A)
            n = ok.sum(axis=0)
            v = np.full(A.shape[1], np.nan)
            np.divide(np.where(ok, A, 0.0).sum(axis=0), n, out=v, where=n > 0)
            v[~g] = np.nan
            peaks.append(int(np.nanargmax(v)) if np.isfinite(v).any() else -1)
        out[x] = np.array(peaks, dtype=int)
    return out


def null_all(W, bins, good, mean, members, est, peaks, ratios, n: int = hp.N_NULL, seed: int = hp.SEED) -> dict:
    """The Hi-C engine's position-aware null, for every chromosome with a peak: partner windows moved to random bins at least
    NULL_EXCL from that partner's anchor and its ends; P = (1 + replicates reaching the observed peak ratio within
    NULL_REACH bins of the peak) / (1 + n). NaN when untestable."""
    rng = np.random.default_rng(seed)
    elig = {}
    for y, c in est.items():
        loc = np.arange(hp.NULL_EXCL, int(bins.nb[y]) - hp.NULL_EXCL)
        loc = loc[np.abs(loc - c) >= hp.NULL_EXCL]
        if len(loc):
            elig[y] = loc
    if len(elig) < 2:
        return {x: float("nan") for x in peaks}
    hits, seen = dict.fromkeys(peaks, 0), dict.fromkeys(peaks, 0)
    for _ in range(n):
        cent = {y: int(rng.choice(loc)) for y, loc in elig.items()}
        S = cluster_score(partner_oe(W, bins, good, mean, window_bins(bins, good, cent)), bins, good, members)
        for x in peaks:
            v = S[bins.sl(x)]
            if not np.isfinite(v).any():
                continue
            seen[x] += 1
            med = np.nanmedian(v)
            lo, hi = max(peaks[x] - hp.NULL_REACH, 0), min(peaks[x] + hp.NULL_REACH + 1, len(v))
            win = v[lo:hi]
            if med > 0 and np.isfinite(win).any() and np.nanmax(win) / med >= ratios[x]:
                hits[x] += 1
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
    gate_refs records which references were used."""
    t = table.copy()
    scored = t.anch_peak_start.ge(0) & t.chance.notna() & t.seq_start.ge(0)
    ref = t[scored & t.seq_confidence.eq("high")]
    if len(ref) < GATE_MIN_REF + 1:
        ref = t[scored & t.seq_confidence.isin(["high", "medium"])]
    refs, agree = "reference", ref.anch_agrees_seq
    if not t.seq_confidence.isin(["high", "medium"]).any():
        ref = t[scored & t.seq_confidence.eq("low")]
        refs, agree = "low", ref.anch_agrees_seq & ref.anch_flags.map(clean_peak)
    majority = 2 * int(agree.sum()) > len(ref)
    ks, ns, ps, passed = [], [], [], []
    for c in t.chrom:
        r = ref.chrom != c
        k, n = int(agree[r].sum()), int(r.sum())
        r = ref[r]
        p = poisson_binomial_sf(k, r.chance.to_numpy()) if n else float("nan")
        ks.append(k)
        ns.append(n)
        ps.append(p)
        passed.append(bool(n >= GATE_MIN_REF and 2 * k > n and p <= GATE_P and majority))
    t["gate_k"], t["gate_n"], t["gate_p"], t["gate_pass"] = ks, ns, ps, passed
    t["gate_refs"] = refs
    return t


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
    boot = bootstrap_all(M, bins, members, win, threads) if pcs else {}
    lopo = lopo_all(oe, bins, good, members) if pcs else {}
    pnull = null_all(W, bins, good, mean, members, est, {k: v["peak"] for k, v in pcs.items()},
                     {k: v["ratio"] for k, v in pcs.items()}) if pcs else {}
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
            flags = call_flags(k, pc, bins, state["hub"], state["rdna"], True, pnull.get(k, np.nan), True)
            r.update(anch_peak_start=int(bins.start[g(pc["peak"])]), anch_peak_end=int(bins.end[g(pc["peak"])]),
                     anch_core_start=int(bins.start[g(pc["core_lo"])]), anch_core_end=int(bins.end[g(pc["core_hi"])]),
                     anch_ratio=round(pc["ratio"], 3), anch_boot_in_core=round(inside, 3), anch_lopo_stable=stable,
                     anch_null_p=pnull.get(k, np.nan), anch_level=lev, anch_flags=",".join(flags),
                     anch_strong=bool(lev == STRONG_HIC_LEVEL and not set(flags) & set(NOT_STRONG_FLAGS)))
            if s >= 0:
                r.update(anch_agrees_seq=_agrees_bin(bins, good, k, pc["peak"], s, e),
                         chance=_chance(S, bins, good, k, s, e))
        rows.append(r)
    t = gate(pd.DataFrame(rows))
    return t, dict(S=S, est=est, win=win)


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
        nors=None) -> dict:
    """Both Hi-C analyses for one genome, and the 3D association of each feature class. Returns dict(table,
    S (anchored score), denovo_S, bins, members, hub, rabl, stats, bin_size, telomere_S, classes, coverage (contacts
    per bin, for the karyotype test of provenance)). nors: the genome's 45S arrays (rdna.detect), flagged rDNA."""
    res = run_denovo(pairs, fai, prefix, satellites, threads)
    st = res["state"]
    st["rdna"] = add_nors(st["rdna"], nors)
    t, extra = anchored(st, seq, threads)
    t = nor_peaks(t, nors)
    t = add_denovo(t, res["rows"], st["bins"])
    fc = feature_contacts(st, classes) if classes else pd.DataFrame(columns=FC_COLS)
    if len(fc):
        sig = fc[fc.p_enriched <= GATE_P]
        log.info(f"Hi-C feature classes: {len(fc)} tested; self-associated (P <= {GATE_P}): "
                 f"{', '.join(f'{r.feature_class} {r.fold:.2f}x' for r in sig.itertuples()) or 'none'}")
    n_pass = int(t.gate_pass.sum())
    log.info(f"Hi-C {Path(pairs).name}: bin {st['size']} bp; anchored peaks on the sequence call "
             f"{int(t.anch_agrees_seq.sum())}/{int(t.seq_start.ge(0).sum())}; de novo "
             f"{int(t.denovo_agrees_seq.sum())}; gate passes on {n_pass}/{len(t)} chromosomes")
    return dict(table=t, S=extra["S"], est=extra["est"], denovo_S=st["loc"].S, bins=st["bins"],
                members=st["members"], hub=st["hub"], rabl=res["rabl"], stats=res["stats"], bin_size=st["size"],
                telomere_S=st["T"], classes=fc, coverage=np.asarray(st["M"].sum(axis=1)).ravel())


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
