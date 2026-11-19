"""The repeat-density line: centromeres as the densest stretch of high-copy sequence, one per chromosome (0.16).

Where centromeres are built from transposons or their remnants (the Tcn clusters of Cryptococcus, VdLTRE9 in Verticillium,
the Osser insertions of Chlorella), each chromosome's densest stretch of sequence repeated elsewhere in the genome marks its
centromere, without any annotation; an LTR annotation without lineages (fungal Ty3, future.md #47) says nothing there. Repeats
of other kinds (rDNA, segmental duplications, subtelomeric families, pericentromeric pile-ups) exist in every genome, so the
readout means nothing until the genome shows its peaks behave like centromeres: the one-per-chromosome test of the AT line
(atline.py), run on this readout with the same machinery (troughs of the negated profile, deepest class, room-aware deals,
Freeman-Tukey lower tail, Sidak over the window sizes).

Readout (REPLOG): per 1-kb bin, the mean over its positions of log2 of the genome-wide count of the position's canonical
K-mer: single-copy sequence adds 0, and there is no copy-number cutoff. Positions inside the satellite engine's arrays are
left out (tandem arrays are the satellite line's, which has its own test; with them in, six dense blocks of holocentric
C. elegans brought the test to its level); a window needs half its length outside them. Windows of 10-500 kb (one transposon
copy fills a smaller window). Windows whose midpoint lies in a terminal eighth are left out: chromosome ends gather repeats
whatever the centromere (Malassezia: the deepest class was the last 10 kb of every chromosome; tomato 7 of 12), the bound of
the cis-trans line and the which-end test (hicfuse.END_ZONE); a centromere in a terminal eighth is out of this line's reach.
Genomes above SAMPLE_ABOVE are read from a fixed 1/SAMPLE share of K-mers (FracMinHash: every copy of a kept K-mer is
kept, so its count is exact; a 1-kb bin then holds about 1000/SAMPLE kept positions, enough over 10-kb windows).

Calls: the densest window at the chosen size, widened to its half-height run. Like the AT line, it places only chromosomes
nothing else places, `low` (fuse._at_places: rep_only). Where its peaks agree beyond chance with the Hi-C readout that placed a
call (with_hic: Hi-C's prior-free peaks, or the cis extremes), a peak at a Hi-C call on a chromosome without a sequence call
gives that call its interval and raises it one level (fuse._refine: rep_pair). It never moves or opposes a sequence call
itself; since 0.19 its peaks test the repeat-family windows (famline.with_repeats), which can let the family line move one.
Development (rules fixed before computing): C. neoformans P 0.0004, 14 of 14 peaks on CENP-A; V. dahliae, Chlorella and
soybean pass, all right; 2,100 block shuffles of 8 genomes passed 6.5% of the time at the nominal 0.05, so the line is held
at 0.01 (REP_ALPHA).
Chlamydomonas (ZeppL LINE clusters) passes too, but its chromosomes all carry low satellite calls, which the line leaves.
C. elegans, C. nigoni, B. mori, M. sympodialis and tomato are silent; soybean and rice, whose centromeres are satellite
arrays, are silent once the arrays are left out."""
from __future__ import annotations

import functools
import os
import pickle

import numpy as np
import pysam

from . import atline, loo
from .hicfuse import END_ZONE
from .runtime import process_pool
from .util import log

K = 21                    # K-mer length: odd (no K-mer is its own reverse complement); chance recurrence of a 21-mer is at
                          # most 2e-4 even in a 0.5-Gb genome, and copies at 90% identity still share 11% of their 21-mers
SAMPLE_ABOVE = 200_000_000    # bp of FASTA above which K-mers are sampled (exact counting holds ~12 bytes per position)
SAMPLE = 32               # keep 1 K-mer in SAMPLE by hash (FracMinHash): about 31 kept positions per 1-kb bin
SCALES = tuple(w for w in atline.SCALES if w >= 10_000)   # one transposon copy (1-10 kb) fills a smaller window
MIN_OUTSIDE = 0.5         # a window needs half its length outside the satellite arrays: a majority
REP_ALPHA = 0.01          # conventional level, after the correction for the sizes tried: stricter than the AT line's 0.05
                          # because block shuffles passed the test at 0.05 6.5% of the time (2,100 shuffles of 8 genomes)
PIECE = 1_000_000         # bp read at a time (a multiple of atline.BIN): streaming, never a whole chromosome string
LUT = np.full(256, 4, dtype=np.uint8)
for _i, _b in enumerate(b"ACGT"):
    LUT[_b] = _i
    LUT[_b + 32] = _i     # soft-masked (lower-case) bases are bases


def fmix64(x: np.ndarray) -> np.ndarray:
    """murmur3's 64-bit finalizer: a fixed, well-mixed hash of each uint64 (wraps modulo 2^64)."""
    x = x ^ (x >> np.uint64(33))
    x = x * np.uint64(0xFF51AFD7ED558CCD)
    x = x ^ (x >> np.uint64(33))
    x = x * np.uint64(0xC4CEB9FE1A85EC53)
    return x ^ (x >> np.uint64(33))


def canonical_kmers(code: np.ndarray, k: int = K) -> tuple:
    """(canonical 2-bit K-mer, valid) of each start 0..len(code)-k; code 0-3 = ACGT, 4 = anything else (invalid)."""
    L = len(code) - k + 1
    if L <= 0:
        return np.zeros(0, np.uint64), np.zeros(0, bool)
    bad = code > 3
    c = np.where(bad, 0, code).astype(np.uint64)
    comp = np.uint64(3) - c
    fwd, rc, two = np.zeros(L, np.uint64), np.zeros(L, np.uint64), np.uint64(2)
    for j in range(k):                       # forward: the first base most significant
        np.left_shift(fwd, two, out=fwd)
        np.bitwise_or(fwd, c[j:j + L], out=fwd)
    for j in range(k - 1, -1, -1):           # reverse complement: the complement of the last base most significant
        np.left_shift(rc, two, out=rc)
        np.bitwise_or(rc, comp[j:j + L], out=rc)
    cs = np.concatenate([[0], np.cumsum(bad, dtype=np.int64)])
    return np.minimum(fwd, rc), (cs[k:] - cs[:-k]) == 0


def _inside(ivs, s: int, e: int) -> np.ndarray:
    """Which positions s..e-1 lie inside the sorted intervals ivs [(start, end)]."""
    out = np.zeros(e - s, bool)
    for a, b in ivs or ():
        if b > s and a < e:
            out[max(a, s) - s:min(b, e) - s] = True
    return out


def pieces(fasta, fai: dict, size: int, read, args, threads: int = 1):
    """(piece, read(fa, *piece)) for every piece (name, L, s, e, *args(name)) of the whole FASTA, [s, e) of at most size bp,
    in the FASTA's order (fai: {sequence: length}, file order). With threads > 1, worker processes read the pieces, each
    with its own handle on the FASTA, in jobs of consecutive pieces of about size bp (short sequences share one), and the
    results come back in the same order: what a caller builds from them does not depend on the number of threads."""
    todo = [(name, L, s, min(L, s + size), *args(name)) for name, L in fai.items() for s in range(0, L, size)]
    if threads > 1:
        jobs, n = [], size
        for p in todo:
            if n + p[3] - p[2] > size:
                jobs.append([])
                n = 0
            jobs[-1].append(p)
            n += p[3] - p[2]
        if len(jobs) > 1:
            with process_pool(min(threads, len(jobs))) as ex:
                for job, res in zip(jobs, ex.map(_job, [(read, str(fasta), j) for j in jobs])):
                    yield from zip(job, pickle.loads(res))
            return
    with pysam.FastaFile(str(fasta)) as fa:
        for p in todo:
            yield p, read(fa, *p)


@functools.lru_cache(maxsize=1)
def _fasta(pid: int, path: str):
    """This process's handle on the FASTA, kept for its next jobs (opening one reads the whole index); pid: a forked worker
    opens its own (a handle shares its file offset across a fork)."""
    return pysam.FastaFile(path)


def _job(x) -> bytes:
    """read on each piece of one job of pieces(), in a worker process, pickled for pieces() to unpickle in the thread that
    keeps the results (unpickled by the pool's own result thread, they would leave its heap fragmented: about their size
    again in freed memory that the counting after them cannot reuse)."""
    read, fasta, job = x
    fa = _fasta(os.getpid(), fasta)
    return pickle.dumps([read(fa, *p) for p in job], protocol=pickle.HIGHEST_PROTOCOL)


def _read(fa, name: str, L: int, s: int, e: int, base: int, ivs, sample: int) -> tuple:
    """Piece [s, e) of a sequence of length L: its counted K-mers (those starting in it, a hash sample of 1 in sample), the
    bin of each (-1: counted, read nowhere: other sequences, arrays) and, on an analyzed chromosome (base: its first bin,
    else -1; ivs: its arrays), (the piece's first bin, called bases per bin, called bases outside the arrays per bin), bins
    numbered over all analyzed chromosomes as base is."""
    raw = np.frombuffer(fa.fetch(name, s, min(L, e + K - 1)).encode(), dtype=np.uint8)
    km, ok = canonical_kmers(LUT[raw])
    km, ok = km[:e - s], ok[:e - s]                       # K-mers starting in [s, e)
    pos = np.flatnonzero(ok)
    km = km[pos]
    if sample > 1:                                         # the hash (murmur3 fmix64) is fixed: same K-mers every run
        sel = (fmix64(km) % np.uint64(sample)) == 0
        km, pos = km[sel], pos[sel]
    bin_id = np.full(len(km), -1, np.int32)                # -1: counted, read nowhere (other sequences, arrays)
    if base < 0:
        return km, bin_id, None
    up = raw[:e - s] & np.uint8(0xDF)                      # upper case
    j = (np.arange(e - s) + s) // atline.BIN
    lo, hi = s // atline.BIN, (e - 1) // atline.BIN + 1
    isn = up != 78
    ins = _inside(ivs, s, e)
    called = np.bincount(j - lo, isn, minlength=hi - lo).astype(np.int64)
    outside = np.bincount(j - lo, isn & ~ins, minlength=hi - lo).astype(np.int64)
    # a K-mer that touches an array is the array's: one starting up to K - 1 bases before it holds its unit too
    touch = _inside([(max(0, a - K + 1), b) for a, b in ivs or ()], s, e)
    read = ~touch[pos]                                     # an array's K-mers count copies, are not read
    bin_id[read] = ((pos[read] + s) // atline.BIN + base).astype(np.int32)
    return km, bin_id, (base + lo, called, outside)


def repeat_bins(fasta, fai: dict, chroms: dict, arrays: dict | None = None, threads: int = 1) -> dict:
    """{chrom: dict(called, outside, kpos, l2)} per atline.BIN bin of each analyzed chromosome: called (non-N) bases, called
    bases outside the satellite arrays, the counted K-mer positions outside them, and the sum of log2 of their genome-wide
    counts. fai: {sequence: length} of the whole FASTA (every sequence counts K-mers: unplaced ones and organelles share
    repeats with the chromosomes); arrays: {chrom: [(start, end)]}. Exact counts up to SAMPLE_ABOVE bp, else a fixed hash
    sample of 1 in SAMPLE. threads: processes reading the FASTA in PIECE-bp pieces (pieces(): the same result at any
    number)."""
    arrays = {c: sorted(v) for c, v in (arrays or {}).items()}
    total = int(sum(fai.values()))
    sample = SAMPLE if total > SAMPLE_ABOVE else 1
    nb = {c: L // atline.BIN + 1 for c, L in chroms.items()}
    off, o = {}, 0
    for c in chroms:
        off[c], o = o, o + nb[c]
    called, outside = np.zeros(o, np.int64), np.zeros(o, np.int64)
    K_all, B_all = [], []
    args = lambda name: (off.get(name, -1), arrays.get(name) if name in off else None, sample)     # noqa: E731
    for _, (km, bin_id, bins) in pieces(fasta, fai, PIECE, _read, args, threads):
        if bins is not None:
            i, cl, ou = bins
            called[i:i + len(cl)] += cl
            outside[i:i + len(ou)] += ou
        K_all.append(km)
        B_all.append(bin_id)
    km, b = np.concatenate(K_all), np.concatenate(B_all)
    del K_all, B_all
    _, inv, cnt = np.unique(km, return_inverse=True, return_counts=True)     # every copy counts, in any sequence or array
    ok = b >= 0
    lg = np.log2(cnt[inv[ok]]).astype(np.float64)
    kpos = np.bincount(b[ok], minlength=o).astype(np.int64)
    l2 = np.bincount(b[ok], weights=lg, minlength=o)
    out = {}
    for c in chroms:
        sl = slice(off[c], off[c] + nb[c])
        out[c] = dict(called=called[sl], outside=outside[sl], kpos=kpos[sl], l2=l2[sl])
    log.debug(f"repeats: {len(km):,} K-mer positions counted (1/{sample} sample) over {total:,} bp")
    return out


def _win(x: np.ndarray, k: int) -> np.ndarray:
    cx = np.concatenate([[0], np.cumsum(x, dtype=np.float64 if x.dtype.kind == "f" else np.int64)])
    return cx[k:] - cx[:-k]


def profile(b: dict, k: int, L: int) -> np.ndarray:
    """The trough profile of one chromosome at k bins (step 1 bin): minus the window's mean log2 K-mer count, +inf where under
    atline.MIN_CALLED of its bases are called, under MIN_OUTSIDE lie outside the arrays, or its midpoint lies in a terminal
    eighth (END_ZONE). The densest repeats are the deepest troughs."""
    if len(b["called"]) < k:
        return np.zeros(0)
    M, O, P, S = _win(b["called"], k), _win(b["outside"], k), _win(b["kpos"], k), _win(b["l2"], k)
    f = np.full(len(M), np.inf)
    ok = (M >= atline.MIN_CALLED * k * atline.BIN) & (O >= MIN_OUTSIDE * k * atline.BIN) & (P > 0)
    f[ok] = -S[ok] / P[ok]
    j = np.arange(len(f))
    mid = (j * atline.BIN + np.minimum(L, (j + k) * atline.BIN)) / 2
    f[np.minimum(mid, L - mid) < END_ZONE * L] = np.inf
    return f


def scan(bins: dict, chroms: dict) -> dict:
    """atline.scan on this readout (atline.scan_profiles): {'scales', 'best', 'p', 'status'}; each call's `gc` field holds the
    window's mean log2 K-mer count (positive: denser repeats)."""
    res = atline.scan_profiles(lambda use, k: {c: profile(bins[c], k, chroms[c]) for c in use}, chroms, SCALES, REP_ALPHA)
    for s in res["scales"].values():
        for x in s["calls"].values():
            x["gc"] = -x["gc"]
    if res["best"] is not None:
        s = res["scales"][res["best"]]
        log.info(f"repeat density: window {res['best']} bp: the densest {s['n']} cover {s['covered']} chromosomes (null "
                 f"{s['null_median']:.0f}), the densest {s['m']} fall one per chromosome at P {res['p']:.2g} after "
                 f"{len(res['scales'])} sizes: {res['status']}")
    return res


def with_hic(res: dict | None, hic: dict | None, chroms: dict, name: str = "repeat density") -> dict | None:
    """res (this line's or the repeat-family line's, famline.py; `name` in the log) with the peaks tested against the deciding
    library's Hi-C readouts, at the window the line's test chose
    (atline.validate, one Hi-C bin): against the prior-free peaks (res['hic']) and, where the cis line acts and the prior-free
    peaks say nothing, against the cis extremes (res['cis']), as atline.with_hic tests the troughs. A Hi-C call that a peak
    tested this way agrees with holds two lines at one locus (fuse: rep_pair). 0.19: this line's windows (and the family
    line's) skip the terminal eighths (END_ZONE), so a window placed at random meets a peak with chance (width + width + 2
    bins) / ((1 - 2 END_ZONE) L), not over the whole length, where it is a quarter too small (as famline.with_repeats); the
    troughs, which may lie anywhere, keep the whole length (atline.with_hic)."""
    if res is None or not res.get("scales") or res.get("best") is None or hic is None:
        return res
    one = dict(res, scales={res["best"]: res["scales"][res["best"]]})
    span = {c: (1 - 2 * END_ZONE) * L for c, L in chroms.items()}
    pk, tol, bins = {}, {}, {}
    for r in hic["table"].itertuples(index=False):
        if r.callable:
            bins[r.chrom] = int(r.bin_size)
            if r.denovo_peak_start >= 0:
                pk[r.chrom], tol[r.chrom] = (int(r.denovo_peak_start), int(r.denovo_peak_end)), int(r.bin_size)
    out = dict(res, hic=atline.validate(one, pk, tol, span), cis=None)
    log.info(f"{name} vs Hi-C's prior-free peaks: {out['hic']['k']}/{out['hic']['n']} (P {out['hic']['p']:.2g}): "
             f"{out['hic']['status']}")
    from .hicfuse import cis_acts                       # here: hicfuse imports much of the Hi-C engine
    cis = hic.get("cis")
    if out["hic"]["status"] not in ("pass", "support") and cis_acts(cis):
        cpk = {c: (int(a), int(b)) for c, (a, b) in (cis.get("peaks") or {}).items() if c in bins}
        out["cis"] = atline.validate(one, cpk, {c: bins[c] for c in cpk}, span)
        log.info(f"{name} vs Hi-C's cis extremes: {out['cis']['k']}/{out['cis']['n']} (P {out['cis']['p']:.2g}): "
                 f"{out['cis']['status']}")
    return out


def on_chromosome(res: dict | None, chrom: str) -> dict | None:
    """The fields fuse.Evidence.rep reads for one chromosome (start, end, density, status, window; validated/support/agrees
    against the prior-free peaks and cis_validated/cis_support/cis_agrees against the cis extremes, with_hic), or None."""
    if res is None or not res.get("scales") or res.get("best") is None:
        return None
    x = res["scales"][res["best"]]["calls"].get(chrom)
    if x is None:
        return None
    v, vc = res.get("hic") or {}, res.get("cis") or {}
    sv, svc = loo.status_on(v, chrom), loo.status_on(vc, chrom)       # the tests without this chromosome (loo)
    return dict(start=int(x["start"]), end=int(x["end"]), density=float(x["gc"]), status=res["status"], window=res["best"],
                validated=sv == "pass", support=sv == "support",
                agrees=bool((v.get("by_window") or {}).get(res["best"], {}).get(chrom)),
                cis_validated=svc == "pass", cis_support=svc == "support",
                cis_agrees=bool((vc.get("by_window") or {}).get(res["best"], {}).get(chrom)))
