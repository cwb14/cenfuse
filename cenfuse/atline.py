"""The AT-trough line: centromeres as the deepest GC troughs of the genome, one per chromosome.

Many regional centromeres sit on the most AT-rich stretch of their chromosome (Plasmodium, Malassezia, Neurospora, blast
fungi), whatever their sequence. That leaves a mark in the FASTA alone, but AT-rich stretches of other kinds exist in every
genome, so a trough means nothing until the genome shows its troughs behave like centromeres: a centromere occurs once per
chromosome, anything else falls in proportion to length. So the test is about the karyotype, not about any one locus
(cenfuse.sat.karyotype).

For each window size (their width is not known in advance, from 3 kb in Malassezia to 300 kb in Neurospora): rank the
genome's independent troughs by depth, and take the deepest class the depths themselves set apart (the largest step between
consecutive troughs, the class holding between N, the chromosomes long enough for that window, and 1.5 N troughs: the
centromeres and a minority of others; read off the depths alone, blind to the chromosomes, so it needs no correction for having
been looked for). Under the null the deepest troughs fall as a random deal: each next trough on a chromosome in proportion to
the windows it has left, a trough there taking its room. One centromere per chromosome makes the class fall more evenly than
such a deal can; the test is the lower tail of how evenly it falls (the Freeman-Tukey statistic of the counts per chromosome
against each chromosome's share, as in a test of a too-good fit; its mid-P over the deals, the statistic taking few values on few
chromosomes). The window size with the smallest P is kept, corrected for the sizes tried. The call on a chromosome is its
deepest trough at that size, widened to the half-depth run around it.

Where Hi-C is given, the prior-free Hi-C peaks test the troughs directly (agreement), an independent line. Alone, the line
gives `low` calls on chromosomes without a sequence call, and in a passing genome chooses a `low` sequence call's competing
cluster when its deepest window lies there and not on the call (fuse._at_places). It never moves a confident satellite or TE
call: satellites are AT-rich or GC-rich whatever their centromeres do."""
from __future__ import annotations

import numpy as np
import pysam

from . import loo
from .sat.karyotype import SEED
from .util import gap, log, poisson_binomial_sf

BIN = 1_000             # GC per 1 kb: the finest window (Malassezia's 3-5 kb centromeres need it), and 1 MB of memory per
                        # 500 Mb of sequence at 2 x 2 bytes per bin
SCALES = (1_000, 2_000, 5_000, 10_000, 20_000, 50_000, 100_000, 200_000, 500_000)   # the known regional centromeres span
                        # them; each step about doubles
MIN_WINDOWS = 10        # a chromosome takes part at a size when it holds at least 10 windows (a trough needs flanks)
MIN_CALLED = 0.8        # a window counts with >= 80% of its bases called (not N)
APART = 5               # independent troughs lie >= 5 windows apart: one regional centromere gives several troughs
                        # closer than that (C. neoformans at 10 kb: 13 of the top 14 lie on centromeres, on 9 chromosomes)
CLASS_MAX = 1.5         # the deepest class of troughs holds at most 1.5 N of them: a chromosome's one centromere and a minority
                        # of others (chromosome ends, AT-rich at both ends, make two per chromosome: Chlamydomonas; heuristic)
PACKING = 0.7476        # Renyi's parking constant: independent troughs >= r windows apart found deepest first fill a stretch
                        # at 0.7476 of the closest possible packing (checked: 1.4-1.7 x n / (2r + 1) on real chromosomes)
AT_DRAWS = 10_000       # null deals per window: P resolution 1e-4, as the karyotype test's
MAX_CHROMS = 200        # more sequences than this are contigs, not chromosomes: no test (the deal of many short sequences passed
                        # noise genomes of 300 contigs: no room to tell one trough per chromosome from one per contig)
AT_ALPHA = 0.05         # conventional level, after the correction for the sizes tried. The line alone gives `low` calls
MIN_AGREE = 3           # agreements needed to speak of a test at all (cenfuse's GATE_MIN_REF)
CHUNK = 10_000_000      # bp read at a time: streaming, never a whole chromosome string above this


def gc_bins(fasta, chroms: dict) -> dict:
    """{chrom: (GC count per BIN, called bases per BIN)} as int32 arrays, read in CHUNK pieces."""
    out = {}
    with pysam.FastaFile(str(fasta)) as fa:
        for c, L in chroms.items():
            nb = L // BIN + 1
            g, m = np.zeros(nb, dtype=np.int32), np.zeros(nb, dtype=np.int32)
            for s in range(0, L, CHUNK):
                e = min(L, s + CHUNK)
                a = np.frombuffer(fa.fetch(c, s, e).upper().encode(), dtype=np.uint8)
                i = (np.arange(len(a)) + s) // BIN
                g += np.bincount(i, (a == 67) | (a == 71), minlength=nb).astype(np.int32)
                m += np.bincount(i, a != 78, minlength=nb).astype(np.int32)
            out[c] = (g, m)
    return out


def _windows(g: np.ndarray, m: np.ndarray, k: int) -> np.ndarray:
    """GC fraction of every window of k bins (step 1 bin), +inf where too few bases are called."""
    if len(g) < k:
        return np.zeros(0)
    cg, cm = np.concatenate([[0], np.cumsum(g, dtype=np.int64)]), np.concatenate([[0], np.cumsum(m, dtype=np.int64)])
    G, M = cg[k:] - cg[:-k], cm[k:] - cm[:-k]
    f = np.full(len(G), np.inf)
    ok = M >= MIN_CALLED * k * BIN
    f[ok] = G[ok] / M[ok]
    return f


def troughs(f: np.ndarray, k: int, n: int) -> list:
    """Up to n deepest independent troughs [(GC fraction, first bin)] of window profile f, at least APART windows apart."""
    f = f.copy()
    out = []
    for _ in range(n):
        if not len(f):
            break
        j = int(np.argmin(f))
        if not np.isfinite(f[j]):
            break
        out.append((float(f[j]), j))
        f[max(0, j - APART * k):j + APART * k + 1] = np.inf
    return out


def half_depth_run(f: np.ndarray, j: int, k: int) -> tuple:
    """(first, last + 1) bin of the run around trough j in which GC stays under half the depth of the trough below the
    chromosome's median window (the half-maximum core used for every other line)."""
    med = float(np.median(f[np.isfinite(f)]))
    cut = f[j] + (med - f[j]) / 2
    lo = hi = j
    while lo > 0 and f[lo - 1] < cut:
        lo -= 1
    while hi < len(f) - 1 and f[hi + 1] < cut:
        hi += 1
    return lo, hi + k


def _capacity(tr: dict, prof: dict, k: int, kept: int) -> tuple:
    """(windows, capacity) per chromosome of tr: its called windows, and the independent troughs it can hold: all of them where its
    list ended before `kept`, else from its windows at the packing density of troughs >= APART windows apart."""
    room = np.array([int(np.isfinite(prof[c]).sum()) for c in tr], dtype=np.int64)
    r = APART * k
    cap = [len(x) if len(x) < kept else max(len(x), int(PACKING * w / (r + 1)) + 1) for w, x in zip(room, tr.values())]
    return room, np.array(cap, dtype=np.int64)


def _deal(room: np.ndarray, capacity: np.ndarray, m: int) -> np.ndarray:
    """(AT_DRAWS, m) array: the chromosome of each of the m deepest troughs of the genome when they fall at random: each next trough
    on a chromosome in proportion to the windows it has left, a trough there taking the room of one of its capacity[c] places (the
    deepest trough lands anywhere; a chromosome with little room left takes fewer: the greedy choice of independent troughs)."""
    rng = np.random.default_rng(SEED)
    m = min(m, int(capacity.sum()))
    per_place = room / capacity
    left, out, rows = np.tile(capacity, (AT_DRAWS, 1)), np.zeros((AT_DRAWS, m), dtype=np.int32), np.arange(AT_DRAWS)
    for t in range(m):
        weight = (left * per_place).cumsum(axis=1)
        pick = (weight > (rng.random(AT_DRAWS) * weight[:, -1])[:, None]).argmax(axis=1)
        out[:, t] = pick
        left[rows, pick] -= 1
    return out


def deepest_class(depths: np.ndarray, n: int) -> int:
    """The size M of the deepest class of the genome's troughs (depths sorted deepest first): where the troughs leave their largest
    step in depth, for M from N to CLASS_MAX * N. Read off the depths alone, blind to the chromosomes, so the test at
    that M needs no correction for having looked for it."""
    top = min(int(CLASS_MAX * n), len(depths) - 1)
    if top <= n:
        return min(n, len(depths))
    return n + int(np.argmax(depths[n:top + 1] - depths[n - 1:top]))


def evenness(counts: np.ndarray, m: int, share: np.ndarray) -> np.ndarray:
    """How evenly m troughs fall on the chromosomes: the Freeman-Tukey deviates of the counts (..., N) from m x each chromosome's
    share, summed over chromosomes; small = more even than in proportion to length, which one centromere per chromosome makes.
    The double-root form (sqrt(c) + sqrt(c + 1) - sqrt(4e + 1), Bishop, Fienberg & Holland 1975) is stable for counts of 0, 1, 2."""
    z = np.sqrt(counts) + np.sqrt(counts + 1.0) - np.sqrt(4.0 * m * share + 1.0)
    return (z * z).sum(axis=-1)


def _counts(labels: np.ndarray, n: int) -> np.ndarray:
    """Counts per chromosome of each row of labels (rows x m, chromosome numbers < n): (rows, n)."""
    rows = len(labels)
    i = labels.ravel().astype(np.int64) + np.repeat(np.arange(rows, dtype=np.int64), labels.shape[1]) * n
    return np.bincount(i, minlength=rows * n).reshape(rows, n)


def _window_test(depths: np.ndarray, labels: np.ndarray, room: np.ndarray, capacity: np.ndarray) -> dict:
    """One window: depths (GC of the kept troughs) and their chromosome numbers (0 to n - 1, n = len(capacity)); room and capacity
    per chromosome. Returns dict(p, best_possible, m, null_median): the mid-P that a random deal gives counts at the deepest class
    as even as observed (few chromosomes give few possible values, and the observed one should not count in full);
    best_possible, the same for the evenest deal drawn; m, the class size; null_median, the chromosomes the N deepest cover in the
    deals (for the log)."""
    n = len(capacity)
    order = np.lexsort((np.random.default_rng(SEED).random(len(depths)), depths))      # ties broken at random, reproducibly
    depths, labels = depths[order], labels[order]
    m = deepest_class(depths, n)
    share = room / room.sum()
    dealt = _deal(room, capacity, m)
    sim = evenness(_counts(dealt, n), m, share)
    obs = evenness(np.bincount(labels[:m], minlength=n), m, share)
    tol = 1e-9
    mid = lambda x: (1 + int((sim < x - tol).sum()) + 0.5 * int((np.abs(sim - x) <= tol).sum())) / (1 + len(sim))
    top_n = dealt[:, :n]
    return dict(p=mid(obs), best_possible=mid(sim.min()), m=m,
                null_median=float(np.median((np.diff(np.sort(top_n, axis=1), axis=1) != 0).sum(axis=1) + 1)))


def scan(bins: dict, chroms: dict) -> dict:
    """{'scales': {window: dict(n, m (size of the deepest class), covered (chromosomes holding the N deepest), p (of the test
    of one per chromosome), best_possible, null_median, calls {chrom: dict(start, end, gc, first (held among the N deepest))})},
    'best': the window with the smallest P or None, 'p': its P corrected for the sizes tried (Sidak), 'status': pass / fail /
    untested (even the most even deal could not pass)}."""
    res = scan_profiles(lambda use, k: {c: _windows(*bins[c], k) for c in use}, chroms, SCALES)
    if res["best"] is not None:
        log.info(f"AT troughs: window {res['best']} bp: the deepest {res['scales'][res['best']]['n']} cover "
                 f"{res['scales'][res['best']]['covered']} chromosomes (null {res['scales'][res['best']]['null_median']:.0f}), "
                 f"the deepest {res['scales'][res['best']]['m']} fall one per chromosome at P {res['p']:.2g} after "
                 f"{len(res['scales'])} sizes: {res['status']}")
    return res


def scan_profiles(profiles, chroms: dict, sizes=SCALES, alpha: float = AT_ALPHA) -> dict:
    """The one-per-chromosome test of scan on any per-window profile, lower = deeper (repline.py reads repeat density this way):
    profiles(use, k) -> {chrom: profile of windows of k bins, step 1 bin, +inf where not read}; sizes: the window sizes tried;
    alpha: the level of the corrected P. Returns scan's dict."""
    if len(chroms) > MAX_CHROMS:
        return dict(scales={}, best=None, p=1.0, status="untested")
    scales = {}
    for w in sizes:
        k = w // BIN
        use = [c for c, L in chroms.items() if L >= MIN_WINDOWS * w]
        if len(use) < 3 or 2 * len(use) <= len(chroms):
            continue
        kept = int(np.ceil(CLASS_MAX * len(use))) + 2                  # the deepest M + 1 of the genome are among each chromosome's
        prof = profiles(use, k)                                        # deepest `kept`
        tr = {c: x for c in use if (x := troughs(prof[c], k, kept))}
        if len(tr) < 3:
            continue
        n = len(tr)
        depths = np.array([g for x in tr.values() for g, _ in x])
        labels = np.array([i for i, x in enumerate(tr.values()) for _ in x], dtype=np.int32)
        room, capacity = _capacity(tr, prof, k, kept)
        t = _window_test(depths, labels, room, capacity)
        first = {c for _, c in sorted((g, c) for c, x in tr.items() for g, _ in x)[:n]}
        calls = {}
        for c, x in tr.items():
            gc, j = x[0]
            lo, hi = half_depth_run(prof[c], j, k)
            calls[c] = dict(start=lo * BIN, end=min(chroms[c], hi * BIN), gc=gc, first=c in first, trough_start=j * BIN,
                            trough_end=min(chroms[c], (j + k) * BIN))
        scales[w] = dict(n=n, m=t["m"], covered=len(first), p=t["p"], best_possible=t["best_possible"],
                         null_median=t["null_median"], calls=calls)
    if not scales:
        return dict(scales={}, best=None, p=1.0, status="untested")
    best = min(scales, key=lambda w: (scales[w]["p"], w))
    m = len(scales)
    padj = 1 - (1 - scales[best]["p"]) ** m
    possible = 1 - (1 - min(s["best_possible"] for s in scales.values())) ** m
    status = "pass" if padj <= alpha else "untested" if possible > alpha else "fail"
    return dict(scales=scales, best=best, p=padj, status=status)


def validate(res: dict, peaks: dict, tol: dict, lengths: dict, alpha: float = 0.01) -> dict:
    """Do the troughs agree with an independent line's loci (Hi-C's prior-free peaks) beyond chance on most chromosomes, at
    any window size (Bonferroni over the sizes)? peaks: {chrom: (start, end)}; tol: {chrom: bp allowed between them}. Returns
    dict(window, k, n, p, status, agrees {chrom: bool}, by_window {window: agrees}, status_by_chrom): pass; support (beyond chance
    but not on most chromosomes: the two lines share information and may support each other where they agree, the rule of Hi-C's
    support tier); fail; untested (fewer than 3 chromosomes to compare). A chromosome without a peak is absent from agrees: it says
    nothing. status_by_chrom: the status at which the troughs may act on each chromosome with a peak, the same test without it
    (loo)."""
    v = _validate(res, peaks, tol, lengths, alpha)
    v["status_by_chrom"] = loo.by_chromosome(v["status"], peaks, lambda names: _validate(
        res, {c: peaks[c] for c in names}, tol, lengths, alpha)["status"])
    return v


def _validate(res: dict, peaks: dict, tol: dict, lengths: dict, alpha: float) -> dict:
    best, by_window = None, {}
    for w, s in res["scales"].items():
        agree, ch = {}, []
        for c, (a, b) in peaks.items():
            x = s["calls"].get(c)
            if x is None:
                continue
            agree[c] = gap(a, b, x["start"], x["end"]) <= tol[c]
            ch.append(min(1.0, ((b - a) + (x["end"] - x["start"]) + 2 * tol[c]) / lengths[c]))
        by_window[w] = agree
        n, k = len(agree), int(sum(agree.values()))
        if n < 3:
            continue
        p = poisson_binomial_sf(k, ch)
        if best is None or (p, w) < (best["p"], best["window"]):
            best = dict(window=w, k=k, n=n, p=p, agrees=agree)
    if best is None:
        return dict(window=None, k=0, n=0, p=1.0, status="untested", agrees={}, by_window=by_window)
    sig = min(1.0, best["p"] * len(res["scales"])) <= alpha and best["k"] >= MIN_AGREE
    best["status"] = ("pass" if 2 * best["k"] > best["n"] else "support") if sig else "fail"
    best["by_window"] = by_window
    return best


def with_hic(res: dict | None, hic: dict | None, chroms: dict) -> dict | None:
    """res with its test against the deciding library's prior-free peaks (validate, one Hi-C bin of tolerance), and the window
    size the rest reads: the window that test chose when it passed; else the FASTA test's when that passed, since that window's
    troughs are what passed; else the window at which Hi-C supports them, where only agreeing chromosomes are read.

    Where the prior-free peaks neither validate nor support the troughs and the library's cis line acts (hicfuse.cis_acts: its
    extremes tested against the LTR cores or the library's own trans contacts), the troughs are tested the same way against the
    cis extremes (0.16, future.md #33): the readout that places Hi-C's calls there (C. neoformans in G1, whose prior-free peaks
    sit at the clustered telomeres). Its result is res['cis'], read like res['hic']."""
    if res is None or not res["scales"]:
        return res
    out = dict(res, window=res["best"], hic=None, cis=None)
    if hic is not None:
        pk, tol = {}, {}
        bins = {}
        for r in hic["table"].itertuples(index=False):
            if r.callable:
                bins[r.chrom] = int(r.bin_size)
            if r.callable and r.denovo_peak_start >= 0:
                pk[r.chrom], tol[r.chrom] = (int(r.denovo_peak_start), int(r.denovo_peak_end)), int(r.bin_size)
        out["hic"] = v = validate(res, pk, tol, chroms)
        log.info(f"AT troughs vs Hi-C's prior-free peaks: window {v['window']} bp: {v['k']}/{v['n']} (P {v['p']:.2g}): "
                 f"{v['status']}")
        if v["status"] == "pass" or (v["status"] == "support" and res["status"] != "pass"):
            out["window"] = v["window"]
        from .hicfuse import cis_acts                       # here: hicfuse imports much of the Hi-C engine
        cis = hic.get("cis")
        if v["status"] not in ("pass", "support") and cis_acts(cis):
            cpk = {c: (int(a), int(b)) for c, (a, b) in (cis.get("peaks") or {}).items() if c in bins}
            out["cis"] = vc = validate(res, cpk, {c: bins[c] for c in cpk}, chroms)
            log.info(f"AT troughs vs Hi-C's cis extremes: window {vc['window']} bp: {vc['k']}/{vc['n']} (P {vc['p']:.2g}): "
                     f"{vc['status']}")
            if vc["status"] == "pass" or (vc["status"] == "support" and res["status"] != "pass"):
                out["window"] = vc["window"]
    return out


def on_chromosome(res: dict | None, chrom: str) -> dict | None:
    """The fields fuse.Evidence.at reads for one chromosome, or None without a trough there."""
    if res is None or not res.get("scales") or res.get("window") is None:
        return None
    x = res["scales"][res["window"]]["calls"].get(chrom)
    if x is None:
        return None
    v = res.get("hic") or {}
    vc = res.get("cis") or {}
    sv, svc = loo.status_on(v, chrom), loo.status_on(vc, chrom)       # the tests without this chromosome (loo)
    return dict(start=int(x["start"]), end=int(x["end"]), gc=float(x["gc"]), status=res["status"], window=res["window"],
                validated=sv == "pass", support=sv == "support",
                agrees=bool((v.get("by_window") or {}).get(res["window"], {}).get(chrom)),
                cis_validated=svc == "pass", cis_support=svc == "support",
                cis_agrees=bool((vc.get("by_window") or {}).get(res["window"], {}).get(chrom)),
                trough_start=int(x.get("trough_start", x["start"])), trough_end=int(x.get("trough_end", x["end"])))
