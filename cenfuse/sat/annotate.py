"""S3: genome-wide nhmmer monomer calls (chromosome chunks on a thread pool), overlap resolution, array chaining."""
import gc
import threading
import time
from collections import deque
from concurrent.futures import ThreadPoolExecutor
from functools import partial

import numpy as np
import pandas as pd
import pyhmmer
import pysam
from pyhmmer.easel import Alphabet, DigitalSequenceBlock, TextSequence
from pyhmmer.plan7 import Background, LongTargetsPipeline

from .params import (ANNOT_CHUNK, FULL_MONOMER_COV, MULTIMER_MIN_COV, MULTIMER_RATIO, NHMMER_REPORT_E, OVERHANG_MIN_BP,
                     OVERLAP_MAX_COV, TRACT_MAX_FRAC, TRACT_MIN_BP, TRACT_MIN_MODEL)
from .util import Progress, log, run

ABC = Alphabet.dna()
MON_COLS = ["chrom", "start", "end", "strand", "family", "bitscore", "evalue",
            "hmm_from", "hmm_to", "model_len", "identity"]
ARR_COLS = ["array_id", "chrom", "start", "end", "strand", "family", "n_monomers", "mean_identity", "monomer_len"]
# what each family's searches cost (search_cost_table). CPU times differ between runs, so this is a work table and never
# a results table: the results of fresh runs stay identical (README, Reproducible)
COST_COLS = ["family", "model_len", "searches", "cpu_h", "cpu_share", "wall_h", "hits_reported", "monomers_kept"]


_PIPE = threading.local()


def search_one(hmm, blk):
    """TopHits of one model against one sequence block, on this thread's own LongTargetsPipeline, made once and cleared
    after every search as pyhmmer's nhmmer workers do; the same hits as pyhmmer.hmmer.nhmmer([hmm], blk, cpus=1,
    E=NHMMER_REPORT_E). A pipeline made per search sits in a reference cycle, so its workspace (model length x the
    256-kb block: ~1.7 GB for an 1,853-bp model) waited for Python's full garbage collection, which a large process
    runs rarely: RSS grew with every search done (chicken GGswu annotate: 277 GB at 24 threads). Reused, a thread
    holds one workspace, its largest so far."""
    pipe = getattr(_PIPE, "pipe", None)
    if pipe is None:
        pipe = _PIPE.pipe = LongTargetsPipeline(ABC, E=NHMMER_REPORT_E)
    try:
        return pipe.search_hmm(hmm, blk)
    finally:
        pipe.clear()


def release_pipeline() -> None:
    """Drop this thread's pipeline and collect it (and those of finished worker threads): free their workspaces."""
    _PIPE.pipe = None
    gc.collect()


def _identity(a) -> float:
    """Exact matches / alignment columns (HMM consensus vs target, case-insensitive)."""
    h = np.frombuffer(a.hmm_sequence.upper().encode(), np.uint8)
    t = np.frombuffer(a.target_sequence.upper().encode(), np.uint8)
    return float((h == t).mean())


def run_significant(df: pd.DataFrame, L: int, evalue: float, gap_factor: float) -> np.ndarray:
    """Which hits are significant alone or as part of a run of whole copies. A copy of a short unit cannot be
    significant alone (a 40-bp model scores <= ~20 bits; E 1e-5 on a 500-Mb chromosome needs ~47), but a tandem run
    of copies is no chance event: after one hit, a random hit starts within the array gap W = gap_factor x unit + 1
    with probability ~W x E_i / L. A run of r whole copies (>= FULL_MONOMER_COV of the model; same family and strand;
    gaps <= gap_factor x unit, as chain_arrays joins monomers) therefore has E = prod(E_i) x (W / L)^(r - 1).
    Fragments count only on their own E: weak partial hits of long models cluster in repeats (not independent), and
    an array is made of whole copies. df: one sequence of length L, E-values for that length (evalue), hmm_from,
    hmm_to, unit in model_len. Returns a mask aligned to df: own E or run E <= evalue."""
    if df.empty:
        return np.zeros(0, bool)
    alone = (np.log10(np.maximum(df.evalue.to_numpy(), 1e-300)) <= np.log10(evalue))
    whole = (df.hmm_to - df.hmm_from + 1).to_numpy() >= FULL_MONOMER_COV * df.model_len.to_numpy()
    out = alone.copy()
    if not whole.any():
        return out
    d = df.loc[whole, ["family", "strand", "start", "end", "evalue", "model_len"]].reset_index(drop=True)
    d = d.assign(_i=np.flatnonzero(whole)).sort_values(["family", "strand", "start"], kind="stable")
    s, e = d.start.to_numpy(), d.end.to_numpy()
    unit = d.model_len.to_numpy().astype(float)
    key = (d.family.astype(str) + d.strand.astype(str)).to_numpy()
    same = np.concatenate([[False], key[1:] == key[:-1]])
    reach = pd.Series(e).groupby(np.cumsum(~same)).cummax().to_numpy()         # furthest end so far in each key
    brk = ~same | np.concatenate([[True], s[1:] - reach[:-1] > gap_factor * unit[1:]])
    run = np.cumsum(brk)
    loge = np.log10(np.maximum(d.evalue.to_numpy(), 1e-300))
    g = pd.DataFrame({"run": run, "loge": loge, "w": np.log10((gap_factor * unit + 1) / L)})
    tot = g.groupby("run").agg(s=("loge", "sum"), n=("loge", "size"), w=("w", "first"))
    run_e = (tot.s + (tot.n - 1) * tot.w).to_numpy()[run - 1]
    out[d._i.to_numpy()] |= run_e <= np.log10(evalue)
    return out


def _rows(hits, mlen: dict, chrom: str, offset: int, scale: float) -> list:
    """MON_COLS rows of the hits HMMER would report (E <= NHMMER_REPORT_E per search; pyhmmer's TopHits also holds
    weaker ones); E-values scaled to the chromosome."""
    rows = []
    for hs in hits:
        fam = hs.query.name
        for h in hs:
            if h.evalue > NHMMER_REPORT_E:
                continue
            d = h.best_domain
            a = d.alignment
            lo, hi = sorted((a.target_from, a.target_to))
            rows.append((chrom, offset + lo - 1, offset + hi, d.strand, fam, float(h.score), float(h.evalue * scale),
                         a.hmm_from, a.hmm_to, mlen[fam], _identity(a)))
    return rows


def significant(rows: list, L: int, evalue: float, gap_factor: float) -> pd.DataFrame:
    """Hits of one sequence (E-values for length L) that are significant alone or in a run; chunk-overlap duplicates
    are dropped first so that no copy counts twice."""
    df = pd.DataFrame(rows, columns=MON_COLS).drop_duplicates(["start", "end", "strand", "family"])
    return df[run_significant(df, L, evalue, gap_factor)].reset_index(drop=True)


def search_chrom(hmms, chrom: str, seq: str, threads: int, evalue: float, gap_factor: float = 2.0) -> list:
    """One nhmmer search of a whole sequence: the reference the chunked search (search_genome) is tested against."""
    if not hmms:
        return []
    blk = DigitalSequenceBlock(ABC, [TextSequence(name=chrom, sequence=seq.upper()).digitize(ABC)])
    rows = _rows(pyhmmer.hmmer.nhmmer(hmms, blk, cpus=threads, E=NHMMER_REPORT_E), {h.name: h.M for h in hmms},
                 chrom, 0, 1.0)
    return list(significant(rows, len(seq), evalue, gap_factor).itertuples(index=False, name=None))


def _chunks(chroms: dict, size: int, overlap: int):
    for c, n in chroms.items():
        s = 0
        while s < n:
            e = min(n, s + size)
            yield c, s, e
            if e == n:
                break
            s = e - overlap


def low_complexity(seqs: dict) -> dict:
    """{name: (n, 2) array of the [start, end) low-complexity stretches of the sequence} by dustmasker at BLAST's
    defaults (level 20, window 64). Names hold no whitespace."""
    out, name = {k: [] for k in seqs}, None
    text = "".join(f">{k}\n{v}\n" for k, v in seqs.items())
    for line in run(["dustmasker", "-outfmt", "interval"], input=text).stdout.splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0]
        elif line:
            a, _, b = line.partition(" - ")
            out[name].append((int(a), int(b) + 1))
    return {k: np.array(v, dtype=np.int64).reshape(-1, 2) for k, v in out.items()}


def neutralize(hmms) -> tuple:
    """(models, {name: its low-complexity tracts}) with those tracts made neutral in long models. A tract makes a model match
    every low-complexity stretch of the genome, so its filters pass nearly everything and each window costs model x window
    DP: 25-50 s per 2-Mb chunk with no hit at all (Lespedeza). At a neutral position the model emits the background, so it
    scores 0 for any base: the tract triggers nothing. Hits of the tract alone go, which mark no copy. A copy usually aligns
    through the tract on the evidence of the rest of the model; where it does not, repair_hits puts it together again.
    Applies to long models with a tract (TRACT_MIN_MODEL, TRACT_MIN_BP) that is not most of the model (TRACT_MAX_FRAC); any
    other model is returned as it is, and the given models are not changed. E-value parameters stay those of the whole
    model."""
    long = {h.name: h.consensus.upper() for h in hmms if h.M >= TRACT_MIN_MODEL}
    if not long:
        return hmms, {}
    found = low_complexity(long)
    background = list(Background(ABC).residue_frequencies)
    out, tracts = [], {}
    for h in hmms:
        iv = found.get(h.name)
        if iv is None or not TRACT_MIN_BP <= int((iv[:, 1] - iv[:, 0]).sum()) < TRACT_MAX_FRAC * h.M:
            out.append(h)
            continue
        neutral = h.copy()
        for a, b in iv.tolist():
            for k in range(a + 1, b + 1):                 # match states count from 1
                for x, freq in enumerate(background):
                    neutral.match_emissions[k, x] = freq
        out.append(neutral)
        tracts[h.name] = iv
    return out, tracts


def _tract_bp(tracts: np.ndarray, lo: int, hi: int) -> int:
    """bp of the tracts ([start, end) from 0) inside model positions lo..hi (from 1, inclusive)."""
    return int(np.maximum(0, np.minimum(tracts[:, 1], hi) - np.maximum(tracts[:, 0], lo - 1)).sum())


def _same_copy(a, b, tracts: np.ndarray, scoring: int) -> bool:
    """Whether hits a and then b along the target (one strand) are pieces of one copy: colinear with the model, a tract between
    them, as far apart on the target as on the model except for the tract (simple sequence varies in length between copies),
    and missing no more scoring positions than a whole copy may."""
    lo, hi = (a[8], b[7]) if a[3] == "+" else (b[8], a[7])              # the model positions either side of the gap
    gap_model, gap_target = hi - lo - 1, b[1] - a[2]
    if gap_model < -OVERHANG_MIN_BP or gap_target < -OVERHANG_MIN_BP:
        return False
    tract = _tract_bp(tracts, lo - OVERHANG_MIN_BP, hi + OVERHANG_MIN_BP)
    return (tract > 0 and abs(gap_target - gap_model) <= tract + OVERHANG_MIN_BP
            and gap_model - tract <= (1 - FULL_MONOMER_COV) * scoring)


def join_split_copies(rows: list, tracts: np.ndarray) -> list:
    """MON_COLS rows of one neutral model on one chunk, with the pieces of one copy (_same_copy) made one row. Where a model holds
    an internal duplicate, a jump to it outscores a walk over a tract that earns nothing, so a copy can come out in pieces
    (potato's 1,362-bp family at 99% identity: 9-514 and 608-1360). Pieces chain by their summed score; scores add up, the best
    E-value stands, identity is weighted by length. Rows not joined keep their order, ahead of the joined."""
    joined, used = [], set()
    M = rows[0][9] if rows else 0
    scoring = M - int((tracts[:, 1] - tracts[:, 0]).sum())
    for strand in ("+", "-"):
        idx = sorted((i for i, r in enumerate(rows) if r[3] == strand), key=lambda i: rows[i][1:3])
        best, prev = [rows[i][5] for i in idx], [-1] * len(idx)
        for j, i in enumerate(idx):
            for k in range(j - 1, -1, -1):
                if rows[idx[k]][1] < rows[i][1] - 2 * M:                 # no piece of this copy starts further back
                    break
                if best[k] + rows[i][5] > best[j] and _same_copy(rows[idx[k]], rows[i], tracts, scoring):
                    best[j], prev[j] = best[k] + rows[i][5], k
        done = set()
        for j in sorted(range(len(idx)), key=lambda x: (-best[x], x)):
            if j in done or prev[j] == -1 or prev[j] in done:
                continue
            chain, k = [], j
            while k != -1 and k not in done:
                chain.append(k)
                k = prev[k]
            done.update(chain)
            ps = [rows[idx[k]] for k in chain]
            used.update(idx[k] for k in chain)
            bp = [r[2] - r[1] for r in ps]
            joined.append((ps[0][0], min(r[1] for r in ps), max(r[2] for r in ps), strand, ps[0][4], sum(r[5] for r in ps),
                           min(r[6] for r in ps), min(r[7] for r in ps), max(r[8] for r in ps), M,
                           sum(r[10] * n for r, n in zip(ps, bp)) / sum(bp)))
    return [r for i, r in enumerate(rows) if i not in used] + joined


def whole_by_scoring_positions(rows: list, tracts: np.ndarray) -> list:
    """MON_COLS rows with a hit over FULL_MONOMER_COV of the positions that can score reported as covering its model, so that it
    counts as a whole copy. A neutral position earns nothing, so a local alignment stops at a tract near a model end; nhmmer's
    correction for biased composition also charges a hit for the biased residues it holds, which can stop it at a tract with
    little beyond (Chlamydomonas' 924-bp family: 59 bp, a 279-bp tract, 586 bp; turning that correction off is no answer: 45%
    more hits in potato, none of them seen by the plain model). The coordinates on the genome stay those of the alignment."""
    out = []
    for r in rows:
        M = r[9]
        scoring = M - int((tracts[:, 1] - tracts[:, 0]).sum())
        covered = r[8] - r[7] + 1 - _tract_bp(tracts, r[7], r[8])
        out.append((*r[:7], 1, M, *r[9:]) if covered >= FULL_MONOMER_COV * scoring else r)
    return out


def _search_work(job) -> int:
    """Progress unit of one (chunk x model) search: chunk bp x model length, what the cost of nhmmer's filters grows
    with (the cost still varies with how satellite-dense the chunk is, so the time left is only an estimate)."""
    _, s, e, h = job
    return (e - s) * h.M


def search_genome(hmms, genome, threads: int, evalue: float, chunk: int = ANNOT_CHUNK, overlap: int = 10_000,
                  ref_len: int = 0, gap_factor: float = 2.0, cost: dict = None) -> dict:
    """{chrom: overlap-resolved monomer hits}. pyhmmer's nhmmer parallelizes over query models only, so a few
    long models would leave most threads idle; instead every (chunk x model) job runs on a thread pool (pyhmmer
    releases the GIL). One job per model, not per chunk: a satellite-dense 2-Mb chunk searched with all ~40 models
    took > 40 min on one thread while the other 24 idled (Arachis); hits and E-values are the same, since each
    model's search is independent. Chunk E-values are rescaled by chromosome/chunk length, so the cutoff equals
    that of one whole-chromosome search. overlap must exceed the longest model; duplicate and edge-truncated
    hits in overlaps are removed by resolve_overlaps. Sequences shorter than ref_len are scaled as if ref_len long,
    so short contigs get the cutoff of a chromosome instead of a looser one. Hits are kept if significant alone or
    as part of a run of copies (run_significant), so short units are found too. The search is the longest stage
    (Panax: > 25 h), so it logs its progress at INFO: each finished search is counted in a done callback, in the order
    searches finish, while their hits are collected in submission order (as Executor.map gave them: the order decides
    which of two duplicate hits is kept). Which family costs what was not recorded, and a few families take most of the
    time (CHM13: the 172-bp alpha satellite 62% of the CPU), so with `cost` (a new dict, filled here with
    {family: dict(model_len, searches, cpu_s, wall_s, hits)}) every finished search adds to its family's totals the
    thread CPU time and the wall time its nhmmer call took on its worker thread (not the time it queued, nor the
    conversion of its hits to rows) and the hits `_rows` kept. wall_s is the sum of the searches' own wall times: the
    time the threads were busy with the family, not the elapsed time of the run. Long models with a low-complexity tract are searched with
    that tract neutral (neutralize)."""
    out = {c: [] for c in genome.chroms}
    if hmms:
        mlen = {h.name: h.M for h in hmms}
        local = threading.local()
        hmms, tracts = neutralize(hmms)
        if tracts:
            log.info(f"annotate: {genome.name}: {len(tracts)} of {len(hmms)} models searched with their low-complexity "
                     f"positions neutral")

        def work(job):
            c, s, e, h = job
            if not hasattr(local, "fa"):
                local.fa, local.key = pysam.FastaFile(str(genome.fasta)), None
            if local.key != (c, s, e):                    # the thread's last chunk is reused when it comes again
                local.blk = DigitalSequenceBlock(ABC, [TextSequence(name=c, sequence=local.fa.fetch(c, s, e).upper())
                                                       .digitize(ABC)])
                local.key = (c, s, e)
            scale = max(genome.chroms[c], ref_len) / (e - s)
            wall0, cpu0 = time.perf_counter(), time.thread_time()       # CPU read inside the wall interval: CPU <= wall
            tophits = search_one(h, local.blk)
            cpu, wall = time.thread_time() - cpu0, time.perf_counter() - wall0
            rows = _rows([tophits], mlen, c, s, scale)
            if h.name in tracts:
                rows = whole_by_scoring_positions(join_split_copies(rows, tracts[h.name]), tracts[h.name])
            return c, rows, (cpu, wall)

        overlap = max(overlap, 2 * max(mlen.values()))
        jobs = [(c, s, e, h) for c, s, e in _chunks(genome.chroms, chunk, overlap) for h in hmms]
        progress = Progress(f"annotate: {genome.name}", sum(map(_search_work, jobs)))
        if cost is not None:                                  # every model has its totals, even one never searched
            cost.update({h.name: dict(model_len=h.M, searches=0, cpu_s=0.0, wall_s=0.0, hits=0) for h in hmms})
        cost_lock = threading.Lock()

        def account(job, fut):
            """Done callback, run by the thread that finished a search: the one place a finished search is accounted,
            for the progress and for the cost of its family (under a lock: callbacks of different threads run at
            once). A search that was cancelled (the run stopped at an error) or failed did no work."""
            if not fut.cancelled() and fut.exception() is None:
                progress.add(_search_work(job))
                if cost is not None:
                    _, rows, (cpu, wall) = fut.result()
                    with cost_lock:
                        total = cost[job[3].name]
                        total["searches"] += 1
                        total["cpu_s"] += cpu
                        total["wall_s"] += wall
                        total["hits"] += len(rows)

        step, t0 = max(1, len(jobs) // 20), time.time()
        with ThreadPoolExecutor(max(1, threads)) as ex:
            futs = deque()                                    # in submission order; a future is dropped once collected
            for job in jobs:
                futs.append(ex.submit(work, job))
                futs[-1].add_done_callback(partial(account, job))
            try:
                for i in range(1, len(jobs) + 1):             # collected in submission order, as ex.map gave them
                    c, rows, _ = futs.popleft().result()
                    out[c].extend(rows)
                    if i % step == 0 or i == len(jobs):       # finer progress with -v
                        log.debug(f"annotate: {genome.name}: {i:,}/{len(jobs):,} searches ({time.time() - t0:.0f} s)")
            finally:
                for fut in futs:                              # after an error, as ex.map did: do not run the rest
                    fut.cancel()
        release_pipeline()
    return {c: resolve_overlaps(drop_multimer_hits(significant(r, max(genome.chroms[c], ref_len), evalue, gap_factor),
                                                   evalue=evalue))
            for c, r in out.items()}


def _union(s: np.ndarray, e: np.ndarray):
    """Merged union of intervals sorted by start: (starts, ends, bp before each interval)."""
    if not len(s):
        return s, e, np.zeros(1, np.int64)
    run = np.maximum.accumulate(e)
    first = np.concatenate([[True], s[1:] > run[:-1]])
    us = s[first]
    ue = np.maximum.reduceat(e, np.flatnonzero(first))
    return us, ue, np.concatenate([[0], np.cumsum(ue - us)])


def _covered_below(x: np.ndarray, us, ue, cum) -> np.ndarray:
    """bp of the union below each position x."""
    k = np.searchsorted(us, x, side="left")
    over = np.where(k > 0, np.maximum(ue[np.maximum(k - 1, 0)] - x, 0), 0)
    return cum[k] - over


def drop_multimer_hits(df: pd.DataFrame, ratio: float = MULTIMER_RATIO, min_cov: float = MULTIMER_MIN_COV,
                       evalue: float = None) -> pd.DataFrame:
    """Describe satellites by their shortest repeat unit: drop a hit when hits of families whose model is
    >= ratio x shorter cover >= min_cov of it. Composite or non-integer multimer families (CentBd x7 + spacer,
    CentO x2.3) otherwise outscore the monomer hits on the monomer arrays, because bitscore grows with length.
    With evalue, only shorter hits significant on their own count: copies significant only as a run
    (run_significant) fill sequence nothing else explains but never erase a hit that is significant alone (Avena:
    run-only 38-41-bp copies erased 87/130-bp arrays and then lost the overlap contest, leaving sequence unannotated).
    Vectorized: per distinct long-model length, the union of the shorter-model hits is built once and each hit's
    coverage read from prefix sums (Avena C genomes: tens of millions of hits)."""
    if df.empty:
        return df
    df = df.sort_values(["start", "end"]).reset_index(drop=True)
    s, e, m = df.start.to_numpy(), df.end.to_numpy(), df.model_len.to_numpy()
    strong = np.ones(len(df), bool) if evalue is None else df.evalue.to_numpy() <= evalue
    keep = np.ones(len(df), bool)
    cand = m >= ratio * m.min()                          # only hits of long models can be multimers
    for t in np.unique(m[cand]):
        short = (m * ratio <= t) & strong
        if not short.any():
            continue
        us, ue, cum = _union(s[short], e[short])
        idx = np.flatnonzero(cand & (m == t))
        cov = _covered_below(e[idx], us, ue, cum) - _covered_below(s[idx], us, ue, cum)
        keep[idx[(cov > 0) & (cov >= min_cov * (e[idx] - s[idx]))]] = False
    return df[keep].reset_index(drop=True)


def resolve_overlaps(df: pd.DataFrame, max_cov: float = OVERLAP_MAX_COV) -> pd.DataFrame:
    """Greedy by bitscore: accept a hit only if <= max_cov of it is already covered by accepted hits
    (coverage bitmap over the chromosome span). A pairwise ">50% of the shorter" rule let hits of phase-shifted
    subfamilies that overlap each of two accepted neighbors by <= 50% survive (A. thaliana CEN180: +27% monomers)."""
    if df.empty:
        return df
    df = df.sort_values(["start", "end"]).reset_index(drop=True)
    s, e = df.start.to_numpy(), df.end.to_numpy()
    lo = int(s.min())
    occ = np.zeros(int(e.max()) - lo, bool)
    keep = np.zeros(len(df), bool)
    for i in np.argsort(-df.bitscore.to_numpy(), kind="stable"):
        a, b = s[i] - lo, e[i] - lo
        if occ[a:b].sum() <= max_cov * (b - a):
            occ[a:b] = True
            keep[i] = True
    return df[keep].reset_index(drop=True)


def chain_arrays(mon: pd.DataFrame, fam_len: dict, p, groups: dict = None):
    """Monomers on one strand whose families share a chain group (close subfamilies; default: the family itself)
    with gap <= array_gap_factor x the longest member model form an array; arrays need >= min_array_monomers.
    An array's family is its majority family by monomer bp; monomers keep their own family.
    Returns (monomers in arrays with array_id, arrays)."""
    if mon.empty:
        return mon.assign(array_id=pd.Series(dtype="int64")), pd.DataFrame(columns=ARR_COLS)
    groups = groups or {}
    mon = mon.assign(_grp=mon.family.map(lambda f: groups.get(f, f)))
    mon = mon.sort_values(["chrom", "_grp", "strand", "start"]).reset_index(drop=True)
    ids = np.full(len(mon), -1, np.int64)
    rows, aid = [], 0
    ident, fam, blen = mon.identity.to_numpy(), mon.family.to_numpy(), (mon.end - mon.start).to_numpy()
    for (c, _, st), g in mon.groupby(["chrom", "_grp", "strand"], sort=False):
        idx = g.index.to_numpy()
        s, e = g.start.to_numpy(), g.end.to_numpy()
        maxgap = p.array_gap_factor * max(fam_len[f] for f in set(fam[idx]))
        brk = np.where(s[1:] - np.maximum.accumulate(e)[:-1] > maxgap)[0] + 1
        for block in np.split(np.arange(len(idx)), brk):
            if len(block) < p.min_array_monomers:
                continue
            ii = idx[block]
            ids[ii] = aid
            major = pd.Series(blen[ii]).groupby(fam[ii]).sum().idxmax()
            rows.append((aid, c, int(s[block].min()), int(e[block].max()), st, major, len(block),
                         float(ident[ii].mean()), int(fam_len[major])))
            aid += 1
    mon["array_id"] = ids
    return mon[mon.array_id >= 0].drop(columns="_grp").copy(), pd.DataFrame(rows, columns=ARR_COLS)


def search_cost_table(cost: dict, monomers: pd.DataFrame) -> pd.DataFrame:
    """COST_COLS table from the cost search_genome filled and the monomers table annotate_genome writes, costliest
    family first (equal costs keep the order of the models). cpu_share is a family's part of the CPU of all searches
    (0 for every family if none was measured); wall_h sums the searches' own wall times (the threads' busy time, not the
    elapsed time); hits_reported counts what nhmmer reported (E <= NHMMER_REPORT_E, before the significance test,
    duplicate removal in chunk overlaps, multimer drop and overlap resolution), monomers_kept the family's rows in
    the final monomers table, after chaining into arrays."""
    total = sum(c["cpu_s"] for c in cost.values())
    kept = monomers.family.value_counts()
    rows = [(f, c["model_len"], c["searches"], c["cpu_s"] / 3600, c["cpu_s"] / total if total > 0 else 0.0,
             c["wall_s"] / 3600, c["hits"], int(kept.get(f, 0))) for f, c in cost.items()]
    return pd.DataFrame(rows, columns=COST_COLS).sort_values("cpu_h", ascending=False, kind="stable", ignore_index=True)


def costliest(table: pd.DataFrame, n: int = 3) -> list:
    """[(family, cpu_share)] of the n families whose searches took the most CPU, costliest first; [] when no CPU was
    measured (nothing to rank). The log names three: enough to show where the time went without a table dump."""
    top = table.sort_values("cpu_share", ascending=False, kind="stable").head(n)
    return [(f, float(s)) for f, s in zip(top.family, top.cpu_share) if s > 0]


def annotate_genome(genome, hmms, fam_len: dict, threads: int, p, out_mon, out_arr, groups: dict = None, out_cost=None):
    """Search the genome with every model, resolve overlaps, chain arrays; write the monomers and arrays tables and,
    with out_cost, the table of what each family's searches cost (search_cost_table). The search dominates the run
    time and is hard to predict (Panax: 7.9 h predicted, > 25 h actual), so the log names the families that took most
    of it."""
    cost = {}
    per = search_genome(hmms, genome, threads, p.nhmmer_evalue, ref_len=p.evalue_ref_len, gap_factor=p.array_gap_factor,
                        cost=cost)
    parts = [per[c] for c in genome.chroms if len(per[c])]
    log.debug(f"annotate: {genome.name}: monomer hits on {len(parts)}/{len(genome.chroms)} sequences")
    mon = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame(columns=MON_COLS)
    mon, arrays = chain_arrays(mon, fam_len, p, groups)
    order = {c: i for i, c in enumerate(genome.chroms)}
    mon = mon.assign(ci=mon.chrom.map(order)).sort_values(["ci", "start", "end"]).drop(columns="ci")
    mon = mon.reset_index(drop=True)
    mon.insert(0, "mon_idx", np.arange(len(mon), dtype=np.int64))
    arrays = arrays.assign(ci=arrays.chrom.map(order)).sort_values(["ci", "start"]).drop(columns="ci")
    mon.to_csv(out_mon, sep="\t", index=False)
    arrays.to_csv(out_arr, sep="\t", index=False)
    log.info(f"annotate: {genome.name}: {len(mon)} monomers in {len(arrays)} arrays")
    table = search_cost_table(cost, mon)
    if out_cost is not None:
        table.to_csv(out_cost, sep="\t", index=False, float_format="%.9g")
    top = costliest(table)
    if top:
        log.info(f"annotate: {genome.name}: search CPU {table.cpu_h.sum():.1f} h; costliest families: "
                 + ", ".join(f"{f} ({100 * share:.0f}%)" for f, share in top))
    return mon, arrays.reset_index(drop=True)
