"""The Hi-C line.

1. The Hi-C engine's prior-free call, unchanged (outputs under OUT/hic/).
2. The anchored test, on the same balanced matrix. Partner windows sit at the other chromosomes' sequence calls
   (midpoint bin +/- the Hi-C engine's HALF_WIN). Each chromosome's score is the mean O/E with those windows; its own call is
   never a partner, so the test is leave-one-chromosome-out by construction. Which trans cluster dominates the
   genome then no longer matters: the score asks whether a place contacts the other chromosomes' (sequence-defined)
   centromeres. Every anchored peak gets the Hi-C engine's own tests: Poisson bootstrap with the windows fixed,
   leave-one-partner-out, and the position-aware null.
3. The gate, per chromosome x: on the other chromosomes with a reference-grade sequence call, does the anchored peak
   land on that call more often than chance (Poisson-binomial tail <= GATE_P, and on a majority)? Only then is
   Hi-C allowed to change x's call."""
from __future__ import annotations

from argparse import Namespace
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

from .hic import cli as hic_cli
from .hic import params as hp
from .hic.confidence import _score_fixed, level, poisson_resample
from .hic.io import InputError
from .hic.locate import peak_core
from .hic.qc import call_flags
from .hic.score import cluster_score, partner_oe, window_bins
from .params import ANCHOR_LEVELS, GATE_MIN_REF, GATE_P, HIC_TOL_BINS, NOT_STRONG_FLAGS, STRONG_HIC_LEVEL
from .util import CenfuseError, gap, log, poisson_binomial_sf

HIC_COLS = ["chrom", "callable", "bin_size", "seq_start", "seq_end", "seq_confidence", "partner",
            "anch_peak_start", "anch_peak_end", "anch_core_start", "anch_core_end", "anch_ratio", "anch_boot_in_core",
            "anch_lopo_stable", "anch_null_p", "anch_level", "anch_flags", "anch_strong", "anch_agrees_seq", "chance",
            "gate_k", "gate_n", "gate_p", "gate_pass",
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


def bootstrap_all(M, bins, members, est, threads: int, n: int = hp.N_BOOT, seed: int = hp.SEED) -> dict:
    """The Hi-C engine's bootstrap (Poisson resamples, re-balanced, partner windows fixed at est), but the peak of every
    callable chromosome, not only of those holding a window."""
    ks = [int(k) for k in np.where(members)[0]]

    def one(ss):
        S = _score_fixed(poisson_resample(M, np.random.default_rng(ss)), bins, members, est)
        return {k: int(np.nanargmax(S[bins.sl(k)])) if np.isfinite(S[bins.sl(k)]).any() else -1 for k in ks}

    with ThreadPoolExecutor(max(threads, 1)) as ex:
        reps = list(ex.map(one, np.random.SeedSequence(seed).spawn(n)))
    return {k: np.array([r[k] for r in reps]) for k in ks}


def lopo_all(oe, bins, good, members) -> dict:
    """Peak of every callable chromosome with each partner window left out in turn (the Hi-C engine's leave_one_out)."""
    out = {}
    for x in (int(k) for k in np.where(members)[0]):
        peaks = []
        for y in oe:
            if y == x:
                continue
            v = cluster_score(oe, bins, good, members, exclude={x: {y}})[bins.sl(x)]
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


def _partners(bins, members, seq: dict) -> dict:
    """{chromosome index: midpoint bin of its sequence call} for calls at ANCHOR_LEVELS (all calls when fewer than
    GATE_MIN_REF qualify)."""
    idx = {n: i for i, n in enumerate(bins.names)}

    def pick(levels):
        out = {}
        for c, (s, e, conf) in seq.items():
            k = idx.get(c)
            if k is not None and members[k] and s >= 0 and conf in levels:
                out[k] = min(int(bins.nb[k]) - 1, ((s + e) // 2) // bins.size)
        return out
    est = pick(ANCHOR_LEVELS)
    return est if len(est) >= GATE_MIN_REF else pick(("high", "medium", "low"))


def _agrees_bin(bins, k: int, j: int, s: int, e: int) -> bool:
    g = int(bins.off[k] + j)
    return gap(int(bins.start[g]), int(bins.end[g]), s, e) <= HIC_TOL_BINS * bins.size


def _chance(S, bins, k: int, s: int, e: int) -> float:
    """Share of chromosome k's scored bins whose peak would agree with [s, e) (a uniformly placed peak)."""
    sl = bins.sl(k)
    ok = np.isfinite(S[sl])
    if not ok.any():
        return float("nan")
    st, en = bins.start[sl], bins.end[sl]
    tol = HIC_TOL_BINS * bins.size
    near = np.maximum(0, np.maximum(st, s) - np.minimum(en, e)) <= tol
    return float((near & ok).sum() / ok.sum())


def gate(table: pd.DataFrame) -> pd.DataFrame:
    """Leave-one-chromosome-out gate (see module doc). Reference: chromosomes with a high sequence call (plus medium
    when fewer than GATE_MIN_REF are high) and a scored anchored peak."""
    t = table.copy()
    scored = t.anch_peak_start.ge(0) & t.chance.notna() & t.seq_start.ge(0)
    ref = t[scored & t.seq_confidence.eq("high")]
    if len(ref) < GATE_MIN_REF + 1:
        ref = t[scored & t.seq_confidence.isin(["high", "medium"])]
    ks, ns, ps, passed = [], [], [], []
    for c in t.chrom:
        r = ref[ref.chrom != c]
        k, n = int(r.anch_agrees_seq.sum()), len(r)
        p = poisson_binomial_sf(k, r.chance.to_numpy()) if n else float("nan")
        ks.append(k)
        ns.append(n)
        ps.append(p)
        passed.append(bool(n >= GATE_MIN_REF and 2 * k > n and p <= GATE_P))
    t["gate_k"], t["gate_n"], t["gate_p"], t["gate_pass"] = ks, ns, ps, passed
    return t


def anchored(state: dict, seq: dict, threads: int) -> tuple[pd.DataFrame, dict]:
    """state: the Hi-C engine's in-memory state; seq: {chrom: (start, end, confidence)} of the sequence calls (start -1 =
    none). Returns (per-chromosome table without the de novo columns, {'S': anchored score, 'est': partners})."""
    M, W, bins, good, mean, members = (state[k] for k in ("M", "W", "bins", "good", "mean", "members"))
    est = _partners(bins, members, seq)
    if len(est) < 2:
        log.warning("Hi-C: fewer than 2 chromosomes with a sequence call to anchor partner windows; anchored test "
                    "skipped")
    oe = partner_oe(W, bins, good, mean, window_bins(bins, good, est)) if len(est) >= 2 else {}
    S = cluster_score(oe, bins, good, members) if oe else np.full(bins.N, np.nan)
    pcs = {}
    for k in (int(x) for x in np.where(members)[0]):
        pc = peak_core(S, bins, k)
        if pc is not None and np.isfinite(pc["ratio"]):
            pcs[k] = pc
    boot = bootstrap_all(M, bins, members, est, threads) if pcs else {}
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
                 chance=np.nan)
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
                r.update(anch_agrees_seq=_agrees_bin(bins, k, pc["peak"], s, e), chance=_chance(S, bins, k, s, e))
        rows.append(r)
    t = gate(pd.DataFrame(rows))
    return t, dict(S=S, est=est)


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


def run(pairs, fai, prefix: Path, satellites, seq: dict, threads: int, classes: dict | None = None) -> dict:
    """Both Hi-C analyses for one genome, and the 3D association of each feature class. Returns dict(table,
    S (anchored score), denovo_S, bins, members, hub, rabl, stats, bin_size, telomere_S, classes, coverage (contacts
    per bin, for the karyotype test of provenance))."""
    res = run_denovo(pairs, fai, prefix, satellites, threads)
    st = res["state"]
    t, extra = anchored(st, seq, threads)
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
