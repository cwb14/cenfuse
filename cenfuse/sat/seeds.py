"""S1: tandem-repeat seeds on chunked chromosomes, from FasTAN (default) or TRF, cached by genome md5. Where FasTAN
crashes on a piece, fastan.chunk_rows recovers it and the crashed intervals are listed beside the cache file."""
import gzip
import hashlib
import json
import os
import shutil
import socket
import tempfile
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import pandas as pd
import pysam

from . import fastan, trf
from .util import log, md5_file, union_bp

ENGINES = {"fastan": fastan, "trf": trf}
POLL = 15   # seconds between checks while another process holds a genome's seed lock
COLS = ["chrom", "start", "end", "period", "copies", "consensus_size",
        "pct_match", "pct_indel", "score", "entropy", "consensus"]   # start/end 0-based half-open
CRASH_COLS = ["chrom", "start", "end", "signal", "action"]   # FasTAN crash table = fastan.chunk_rows events


def chunks(chroms: dict, size: int, overlap: int):
    for c, n in chroms.items():
        s = 0
        while s < n:
            e = min(n, s + size)
            yield c, s, e
            if e == n:
                break
            s = e - overlap


def _run_chunk(job):
    """(fasta, seed rows, crash events) of one chunk. Only FasTAN crashes by signal on some arrays and recovers from
    it (fastan.chunk_rows splits the piece); TRF has no such failure to route around."""
    fasta, chrom, s, e, tmpdir, engine, opt, overlap = job
    with pysam.FastaFile(str(fasta)) as fa:
        seq = fa.fetch(chrom, s, e)
    events = []
    if engine == "fastan":
        rows = fastan.chunk_rows(chrom, seq, s, tmpdir, opt, events, overlap)
    else:
        rows = ENGINES[engine].chunk_rows(chrom, seq, s, tmpdir, opt)
    return str(fasta), rows, events


def dedupe(rows: list, tol: float = 0) -> list:
    """Drop chunk-overlap duplicates: same chrom, period within max(1 bp, tol x period), overlap >= 50% of the
    shorter -> keep the higher score."""
    rows = sorted(rows, key=lambda r: (r[0], r[1], r[2], r[3]))
    kept, active, cur = [], [], None
    for r in rows:
        if r[0] != cur:
            cur, active = r[0], []
        active = [i for i in active if kept[i][2] > r[1]]
        dup = None
        for i in active:
            k = kept[i]
            if abs(k[3] - r[3]) <= max(1, tol * min(k[3], r[3])) and \
                    min(k[2], r[2]) - max(k[1], r[1]) >= 0.5 * min(k[2] - k[1], r[2] - r[1]):
                dup = i
                break
        if dup is None:
            kept.append(r)
            active.append(len(kept) - 1)
        elif r[8] > kept[dup][8]:
            kept[dup] = r
    return kept


def genome_key(g) -> str:
    """FASTA md5; plus the chromosome set when the genome uses a subset of the FASTA's sequences (cenfuse hands
    the satellite engine a view whose .fai lists chromosomes only), so a subset and the whole file never share a
    cache entry. The real file's own .fai lists its sequences; without one (the user's FASTA is a link and its .fai
    sits beside the link) the set cannot be checked against them, so it always counts: else a changed set would
    reuse the old set's seeds."""
    md5 = md5_file(g.fasta)
    own = Path(os.path.realpath(g.fasta) + ".fai")
    if own.exists():
        with open(own) as f:
            names = {line.split("\t", 1)[0] for line in f if line.strip()}
        subset = names != set(g.chroms)
    else:
        subset = True
    if subset:
        md5 += "-" + hashlib.md5("\n".join(sorted(g.chroms)).encode()).hexdigest()[:8]
    return md5


def cache_path(cache_dir, md5: str, p) -> Path:
    """Cache file for one genome: keyed by the FASTA md5, the engine and the settings that shape the result."""
    key = hashlib.md5(json.dumps(ENGINES[p.seed_engine].KEY + [p.seed_chunk, p.seed_overlap]).encode()).hexdigest()[:8]
    return Path(cache_dir) / f"{md5}.{key}.{p.seed_engine}.tsv.gz"


def lock_path(path: Path) -> Path:
    return Path(str(path) + ".lock")


def crash_path(cache_file) -> Path:
    """The FasTAN crash table of a cached seed table. It lives beside the cache file, so a shared seed cache keeps it
    and every run that reuses the seeds finds it."""
    return Path(f"{cache_file}.fastan_crashes.tsv")


def _alive(pid: int) -> bool:
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    return True


def _owner(lock: Path) -> tuple:
    """(host, pid) of a lock; a lock written by satellite engine < 0.5 holds only the pid (single-host use)."""
    x = lock.read_text().split()
    return (socket.gethostname(), int(x[0])) if len(x) == 1 else (x[0], int(x[1]))


def _acquire(lock: Path) -> bool:
    """Create the lock file atomically ('host pid'); steal it only if its owner ran on this host and is dead. A
    process on another host (shared filesystem, SLURM nodes) cannot be checked from here, so its lock is never
    stolen: remove a stale one by hand."""
    try:
        fd = os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
    except FileExistsError:
        try:
            host, pid = _owner(lock)
        except (OSError, ValueError, IndexError):
            return False                          # being written right now; retry next poll
        if host == socket.gethostname() and pid and not _alive(pid):
            log.warning(f"seeds: removing stale lock {lock.name} (pid {pid} is gone)")
            lock.unlink(missing_ok=True)
            return _acquire(lock)
        return False
    os.write(fd, f"{socket.gethostname()} {os.getpid()}".encode())
    os.close(fd)
    return True


def _write_crashes(table: Path, genome: str, events: list) -> None:
    """Write the crash table (one row per crash) and log the genome's summary: chunks crashed, pieces split, bp in
    skipped pieces. Skipped pieces overlap their neighbours, so the bp are an upper bound of what has no seeds."""
    flat = [x for chunk in events for x in chunk]
    table.write_text("\t".join(CRASH_COLS) + "\n" + "".join("\t".join(map(str, x)) + "\n" for x in flat))
    skipped = {}
    for chrom, start, end, _, action in flat:
        if action == "skipped":
            skipped.setdefault(chrom, []).append((start, end))
    n_split = sum(x[4] == "split" for x in flat)
    log.info(f"seeds: {genome}: FasTAN crashed on {len(events)} chunk(s); {n_split} piece(s) split, "
             f"{sum(map(union_bp, skipped.values())):,} bp skipped; intervals in {table}")


def _write_cache(path: Path, g, rows: list, tol: float, events=()) -> None:
    """Cache one genome's seeds. events: the FasTAN crashes, one list per chunk it crashed on (fastan.chunk_rows). The
    crash table is written first: the cache file is what other runs wait for, so it must not appear without it. A table
    left by an earlier run of this cache file is dropped when nothing crashed this time."""
    kept = dedupe(rows, tol)
    if events:
        _write_crashes(crash_path(path), g.name, events)
    else:
        crash_path(path).unlink(missing_ok=True)
    part = Path(f"{path}.{os.getpid()}.part")
    with gzip.open(part, "wt") as f:
        f.write("\t".join(COLS) + "\n")
        for r in kept:
            f.write("\t".join(map(str, r)) + "\n")
    part.rename(path)
    log.info(f"seeds: {g.name} {len(kept)} tandem repeats")


def _compute(todo: dict, cache_dir: Path, threads: int, p) -> None:
    """Jobs run genome by genome and results arrive in job order, so each genome is cached as soon as its
    last chunk is back: a later failure loses only the unfinished genomes."""
    engine = ENGINES[p.seed_engine]
    tmp = tempfile.mkdtemp(prefix="seeds_tmp_", dir=cache_dir)
    try:
        opt = engine.setup(tmp)
        jobs, left, target = [], {}, {}
        for path, g in todo.items():
            mine = [(g.fasta, c, s, e, tmp, p.seed_engine, opt, p.seed_overlap)
                    for c, s, e in chunks(g.chroms, p.seed_chunk, p.seed_overlap)]
            jobs += mine
            left[str(g.fasta)], target[str(g.fasta)] = len(mine), (path, g)
        log.info(f"seeds: {p.seed_engine}, {len(todo)} genome(s), {len(jobs)} chunks, {threads} workers")
        rows = {fa: [] for fa in left}
        crashes = {fa: [] for fa in left}                           # the events of each chunk FasTAN crashed on
        for fa in [fa for fa, n in left.items() if not n]:          # only zero-length sequences: nothing to run
            _write_cache(*target[fa], rows.pop(fa), engine.DUP_TOL)
        with ProcessPoolExecutor(threads) as ex:
            for i, (fa, r, ev) in enumerate(ex.map(_run_chunk, jobs, chunksize=1), 1):
                rows[fa].extend(r)
                if ev:
                    crashes[fa].append(ev)
                left[fa] -= 1
                log.debug(f"seeds: chunk {i}/{len(jobs)}")
                if not left[fa]:
                    _write_cache(*target[fa], rows.pop(fa), engine.DUP_TOL, crashes.pop(fa))
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def seed_genomes(genomes, cache_dir, threads: int, p) -> dict:
    """Find the tandem repeats of every genome not cached yet (md5-keyed). Returns {genome name: cached .tsv.gz
    path}. Safe for concurrent runs sharing the cache: a per-genome lock file makes one process compute a genome
    while others wait for its result (stale locks of dead processes are stolen)."""
    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)
    out, todo = {}, {}
    for g in genomes:
        path = cache_path(cache_dir, genome_key(g), p)
        out[g.name] = path
        if path.exists():
            log.info(f"seeds: {g.name} cached ({path.name})")
        else:
            todo.setdefault(path, g)          # md5-identical genomes share one run
    waiting_logged = set()
    while True:
        todo = {path: g for path, g in todo.items() if not path.exists()}
        if not todo:
            return out
        mine = {path: g for path, g in todo.items() if _acquire(lock_path(path))}
        if mine:
            try:
                _compute(mine, cache_dir, threads, p)
            finally:
                for path in mine:
                    lock_path(path).unlink(missing_ok=True)
            continue
        for path, g in todo.items():
            if path not in waiting_logged:
                log.info(f"seeds: {g.name} is being computed by another run; waiting for {path.name} "
                         f"(lock {lock_path(path).name})")
                waiting_logged.add(path)
        time.sleep(POLL)


def load_seeds(path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype={"chrom": str, "consensus": str})
