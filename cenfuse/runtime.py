"""How a run is kept: content fingerprints, atomic writes, lock files (an output folder, a shared seed cache), what a run
costs, and worker pools whose log records reach this process.

Kept out of every stage's code key (pipeline.KEYLESS), although the engines import it: nothing here changes a result, so
editing it recomputes nothing."""
from __future__ import annotations

import contextlib
import hashlib
import logging
import logging.handlers
import multiprocessing
import os
import shutil
import socket
import threading
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

log = logging.getLogger("cenfuse")

WHOLE_FILE = 32 << 20     # files up to 32 MiB are hashed whole (an .fai, a GFF of genes, most tables: well under a second)
SAMPLE = 1 << 20          # larger ones (pairs, LTR annotations, bedMethyl): size + four 1-MiB samples at fixed offsets
HEARTBEAT = 60            # s between the owner's touches of a lock it holds
STALE_AFTER = 15 * 60     # a lock whose owner cannot be checked (another host, a garbled lock) is stale after 15 missed
                          # heartbeats (heuristic)
GUARD_STALE = 60          # s: a takeover guard older than this was left by a process that died while taking a lock over


def file_fingerprint(path, nbytes: int = SAMPLE) -> str:
    """Fingerprint of a file's content, not of its name, place or time stamp: a moved or copied folder (cp without -p) keeps
    its checkpoints, and an input rewritten with other content does not. Files up to WHOLE_FILE are hashed whole; larger
    ones by their size and nbytes at the start, at one and two thirds, and at the end (a change elsewhere that keeps the
    size goes unseen: these are tool outputs, regenerated whole, not edited in place; the FASTA, which can be edited in place,
    is keyed by its whole md5: fasta_md5)."""
    p = Path(path)
    size = p.stat().st_size
    h = hashlib.md5(f"{size}:".encode())
    with open(p, "rb") as f:
        if size <= WHOLE_FILE:
            for block in iter(lambda: f.read(1 << 20), b""):
                h.update(block)
        else:
            for off in (0, size // 3, 2 * size // 3, size - nbytes):
                f.seek(off)
                h.update(f.read(nbytes))
    return h.hexdigest()[:16]


_MD5 = {}


def fasta_md5(path) -> str:
    """md5 of a whole file (the FASTA a stage reads: a same-length edit such as hard-masking or polishing keeps its .fai and
    its size, and must still recompute what reads it), once per process and file (about 1 s per 0.5 GB)."""
    p = Path(path).resolve()
    st = p.stat()
    key = (str(p), st.st_size, st.st_mtime_ns)
    if key not in _MD5:
        h = hashlib.md5()
        with open(p, "rb") as f:
            for block in iter(lambda: f.read(1 << 22), b""):
                h.update(block)
        _MD5[key] = h.hexdigest()
    return _MD5[key]


def _tmp_name(path: Path, pid: int, tid: int) -> Path:
    return path.with_name(f".{path.stem}.{pid}.{tid}.tmp{path.suffix}")


def _drop_stale_tmp(path: Path) -> None:
    """Remove the temporary files of earlier writes of path whose process is gone (a run killed mid-write, by the OOM killer
    say, leaves one; a Hi-C matrix cache's can take gigabytes)."""
    head, tail = f".{path.stem}.", f".tmp{path.suffix}"
    try:
        names = [n for n in os.listdir(path.parent) if n.startswith(head) and n.endswith(tail)]
    except OSError:
        return
    for n in names:
        mid = n[len(head):len(n) - len(tail)].split(".")
        if len(mid) == 2 and mid[0].isdigit() and int(mid[0]) != os.getpid() and not _alive(int(mid[0])):
            try:
                os.unlink(path.parent / n)
            except OSError:
                pass


@contextlib.contextmanager
def atomic_path(path):
    """A temporary path beside path, moved over path when the block ends without an error (removed otherwise). A reader, a
    resume after a kill, or a hard-linked copy of the folder (cp -al, rsync --link-dest) never sees a half-written file, and
    the copy's own file is never rewritten: writing in place would change both. The temporary name keeps path's suffix, so a
    writer that reads the format from it (pandas' compression, matplotlib) writes the same content (a gzip header names the
    temporary file, as its time stamp already differed between runs). A writer that makes no file (an empty PdfPages) leaves
    path as it was."""
    path = Path(path)
    _drop_stale_tmp(path)
    tmp = _tmp_name(path, os.getpid(), threading.get_ident())
    try:
        yield tmp
        if tmp.exists():
            os.replace(tmp, path)
    finally:
        try:
            tmp.unlink()
        except FileNotFoundError:
            pass


@contextlib.contextmanager
def atomic_open(path, mode: str = "w", **kw):
    """open(path, mode) for writing, atomically (atomic_path)."""
    with atomic_path(path) as tmp, open(tmp, mode, **kw) as f:
        yield f


def atomic_csv(df, path, **kw) -> None:
    """df.to_csv(path, **kw), atomically (atomic_path)."""
    with atomic_path(path) as tmp:
        df.to_csv(tmp, **kw)


def atomic_copy(src, dst) -> None:
    """shutil.copyfile(src, dst), atomically (atomic_path)."""
    with atomic_path(dst) as tmp:
        shutil.copyfile(src, tmp)


def atomic_write_text(path, text: str) -> None:
    """Write text to a temporary file beside path, then move it over path (atomic_path)."""
    with atomic_path(path) as tmp:
        tmp.write_text(text)


def unshare_tree(folder) -> int:
    """Unlink every file under folder that has other hard links, before a stage recomputes into it: the stage's tools (BLAST,
    FasTAN, MAFFT) write their files in place, which would change the other copy of a folder copied with cp -al. Only a
    recomputed stage's own folder is passed, so what is unlinked is out of date anyway. Returns the number of files."""
    n = 0
    for root, _, files in os.walk(folder):
        for name in files:
            p = os.path.join(root, name)
            try:
                if os.lstat(p).st_nlink > 1:
                    os.unlink(p)
                    n += 1
            except FileNotFoundError:
                pass
    return n


# ------------------------------------------------------------------ lock files
def _alive(pid: int) -> bool:
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    return True


def _start_ticks(pid: int) -> str | None:
    """The process's start time in clock ticks since boot (Linux /proc/<pid>/stat, field 22): with the pid, what tells one
    process from a later one that reuses its pid. None where the system does not say."""
    try:
        stat = Path(f"/proc/{pid}/stat").read_text()
        return stat[stat.rindex(")") + 2:].split()[19]            # fields after the command name, which may hold spaces
    except (OSError, ValueError, IndexError):
        return None


def lock_owner(lock: Path) -> tuple:
    """(host, pid) written in a lock file ('host pid ...'); a lock of satellite engine < 0.5 holds only the pid."""
    x = Path(lock).read_text().split()
    return (socket.gethostname(), int(x[0])) if len(x) == 1 else (x[0], int(x[1]))


def _lock_ticks(lock: Path) -> str | None:
    """The owner's start ticks a 0.14 lock records (its third field, 'ticks=N'), else None (older locks hold none)."""
    try:
        x = Path(lock).read_text().split()
        t = x[2][len("ticks="):] if len(x) > 2 and x[2].startswith("ticks=") else ""
        return t if t.isdigit() else None
    except OSError:
        return None


def lock_stale(lock: Path) -> str:
    """Why a lock is stale, or '' if its owner may still run. On this host: its pid is gone, or the pid now belongs to a
    process that started after the lock (a reused pid). On another host (a shared filesystem, cluster nodes), where the pid
    cannot be checked: no heartbeat for STALE_AFTER seconds. A lock that cannot be read (being written, or left garbled by a
    process killed between creating and writing it) is stale only once it is that old."""
    lock = Path(lock)
    try:
        age = time.time() - lock.stat().st_mtime
    except OSError:
        return ""
    try:
        host, pid = lock_owner(lock)
    except (OSError, ValueError, IndexError):
        return f"unreadable for {age / 60:.0f} min" if age > STALE_AFTER else ""
    if host == socket.gethostname():
        if not _alive(pid):
            return f"pid {pid} on this host is gone"
        made, now = _lock_ticks(lock), _start_ticks(pid)
        if made is not None and now is not None and made != now:
            return f"pid {pid} on this host now belongs to another process (a reused pid)"
        return ""
    return f"no heartbeat from {host} (pid {pid}) for {age / 60:.0f} min" if age > STALE_AFTER else ""


def _create(lock: Path, what: str) -> bool:
    try:
        fd = os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
    except FileExistsError:
        return False
    ticks = _start_ticks(os.getpid()) or "."
    os.write(fd, f"{socket.gethostname()} {os.getpid()} ticks={ticks} {time.strftime('%Y-%m-%dT%H:%M:%S%z')} {what}".strip().encode())
    os.close(fd)
    return True


def try_lock(lock: Path, what: str = "") -> bool:
    """Create lock atomically, holding 'host pid start-time what'. An existing lock is taken over only when stale
    (lock_stale), with a warning, and under a takeover guard (lock.takeover, created exclusively): of several processes that
    find one stale lock, one replaces it, and none can remove a live lock another just made. Returns whether this process
    holds the lock now."""
    lock = Path(lock)
    if _create(lock, what):
        return True
    if not lock.exists():                             # released between the two looks
        return _create(lock, what)
    if not lock_stale(lock):
        return False
    guard = lock.with_name(f"{lock.name}.takeover")
    try:
        os.close(os.open(guard, os.O_CREAT | os.O_EXCL | os.O_WRONLY))
    except FileExistsError:
        try:                                          # a guard left by a process that died while taking over
            if time.time() - guard.stat().st_mtime > GUARD_STALE:
                guard.unlink(missing_ok=True)
        except OSError:
            pass
        return False
    try:
        if not lock.exists():                         # released meanwhile
            return _create(lock, what)
        why = lock_stale(lock)                        # again, under the guard: it may be another run's live lock by now
        if not why:
            return False
        log.warning(f"removing stale lock {lock} ({why})")
        lock.unlink(missing_ok=True)
        return _create(lock, what)
    finally:
        guard.unlink(missing_ok=True)


def release_lock(lock: Path) -> None:
    """Remove lock if this process holds it (never another run's)."""
    try:
        if lock_owner(lock) == (socket.gethostname(), os.getpid()):
            Path(lock).unlink(missing_ok=True)
    except (OSError, ValueError, IndexError):
        pass


# ------------------------------------------------------------------ what a run costs (future.md #5, #61)
SAMPLE_EVERY = 30         # s between samples of the process tree's memory (a peak shorter than this can be missed)


def _children(pid: int) -> list:
    out = []
    try:
        for t in os.listdir(f"/proc/{pid}/task"):
            with open(f"/proc/{pid}/task/{t}/children") as f:
                out += [int(x) for x in f.read().split()]
    except (OSError, ValueError):
        pass
    return out


def _memory_kb(pid: int) -> tuple:
    """(kB, measure): the proportional set size ('pss': pages shared with other processes counted in part, and forked workers
    share most of theirs), or where the kernel gives none (before Linux 4.14) the resident size ('rss': a shared page counted in
    every process that maps it, so a sum over forked workers overstates)."""
    for path, key, kind in ((f"/proc/{pid}/smaps_rollup", "Pss:", "pss"), (f"/proc/{pid}/status", "VmRSS:", "rss")):
        try:
            with open(path) as f:
                for line in f:
                    if line.startswith(key):
                        return int(line.split()[1]), kind
        except (OSError, ValueError, IndexError):
            continue
    return 0, None


def _can_list_children() -> bool:
    return os.path.exists(f"/proc/self/task/{threading.get_native_id()}/children")


def tree_memory(pid: int | None = None) -> tuple:
    """(kB, measures) of a process and of every process it started, still running; (None, set()) without /proc (macOS)."""
    pid = os.getpid() if pid is None else pid
    if not os.path.exists(f"/proc/{pid}"):
        return None, set()
    seen, todo, total, kinds = set(), [pid], 0, set()
    while todo:
        p = todo.pop()
        if p not in seen:
            seen.add(p)
            kb, kind = _memory_kb(p)
            total += kb
            kinds.add(kind)
            todo += _children(p)
    return total, kinds - {None}


def tree_memory_kb(pid: int | None = None) -> int | None:
    """Memory of a process and of every process it started, still running (kB), or None without /proc (macOS)."""
    return tree_memory(pid)[0]


class Usage:
    """While active, sample the memory of this process and of every process it starts (pool workers, FasTAN, BLAST, ...) every
    SAMPLE_EVERY seconds (a daemon thread): GNU time's maximum counts only the largest single process. cost() gives the wall
    time, the CPU time of this process and its finished children, and the peak of the samples (None without /proc)."""

    def __init__(self, every: float = SAMPLE_EVERY):
        self.every, self.peak, self.kinds, self._stop = every, None, set(), threading.Event()
        self._t = threading.Thread(target=self._sample, daemon=True)
        self.t0, self.c0 = time.time(), os.times()

    def _sample(self):
        while True:
            kb, kinds = tree_memory()
            if kb is not None:
                self.peak = max(self.peak or 0, kb)
                self.kinds |= kinds
            if self._stop.wait(self.every):
                return

    def memory_measure(self) -> str:
        """What peak_memory_gb measures, said as it is."""
        if self.peak is None:
            return "not measured (no /proc)"
        what = ("proportional set size" if self.kinds == {"pss"} else
                "resident size (shared pages counted in every process: an overestimate)" if "rss" in self.kinds else "memory")
        tree = "the process tree's" if _can_list_children() else "this process's (its children cannot be listed here)"
        return f"{tree} {what}, sampled every {self.every:g} s"

    def cost(self) -> dict:
        c = os.times()
        cpu = sum(getattr(c, k) - getattr(self.c0, k) for k in ("user", "system", "children_user", "children_system"))
        return dict(wall_minutes=round((time.time() - self.t0) / 60, 2), cpu_hours=round(cpu / 3600, 4),
                    peak_memory_gb=None if self.peak is None else round(self.peak / 2 ** 20, 2), memory=self.memory_measure())

    def __enter__(self):
        self._t.start()
        return self

    def __exit__(self, *exc):
        self._stop.set()
        self._t.join(timeout=5)
        return False


class Heartbeat:
    """While active, touch the given lock files every HEARTBEAT seconds (a daemon thread), so runs on other hosts can tell a
    live owner from a dead one (lock_stale). A lock this process no longer holds (taken over by another run) is not
    touched, and said once."""

    def __init__(self, *locks):
        self.locks, self._stop = [Path(x) for x in locks], threading.Event()
        self._t = threading.Thread(target=self._beat, daemon=True)
        self._lost = set()

    def _beat(self):
        me = (socket.gethostname(), os.getpid())
        while not self._stop.wait(HEARTBEAT):
            for lock in self.locks:
                try:
                    if lock_owner(lock) != me:
                        if lock not in self._lost:
                            self._lost.add(lock)
                            log.warning(f"lock {lock} is no longer held by this run (another run took it over)")
                        continue
                    os.utime(lock)
                except (OSError, ValueError, IndexError):
                    pass

    def __enter__(self):
        self._t.start()
        return self

    def __exit__(self, *exc):
        self._stop.set()
        self._t.join(timeout=5)
        return False


# ------------------------------------------------------------------ worker pools (future.md #60)
LOGGERS = ("cenfuse", "cenfuse-sat", "cenfuse-hic", "cenfuse-meth")     # util.ENGINE_LOGGERS (runtime imports no cenfuse module)


class _Forward:
    """QueueListener target: a worker's record, handled by this process's logger of the same name (its handlers, its format)."""

    def handle(self, record):
        logging.getLogger(record.name).handle(record)


def _worker_logging(q, levels: dict) -> None:
    for name, level in levels.items():
        lg = logging.getLogger(name)
        lg.handlers[:] = [logging.handlers.QueueHandler(q)]
        lg.setLevel(level)
        lg.propagate = False


LISTENER_STOP = 10        # s to wait for the log forwarding to stop: a worker killed while it held the queue's lock (SIGTERM,
                          # the OOM killer) leaves it blocked for good, and a run must not hang on that


@contextlib.contextmanager
def process_pool(workers: int, mp_context=None):
    """A ProcessPoolExecutor whose workers' log records (FasTAN's crash lines, say) reach this process's handlers. Under fork a
    worker inherits the handlers, so nothing more is needed (and nothing is added: a worker killed while sending a record on a
    queue could block it). Under spawn (macOS) or forkserver (Python 3.14's Linux default) a worker starts with no handlers and
    its records were lost: there they come over a queue, whose listener is stopped within LISTENER_STOP seconds whatever a
    worker did."""
    ctx = mp_context or multiprocessing.get_context()
    if ctx.get_start_method() == "fork":
        with ProcessPoolExecutor(workers, mp_context=ctx) as ex:
            yield ex
        return
    q = ctx.Queue()
    listener = logging.handlers.QueueListener(q, _Forward())
    listener.start()
    levels = {n: logging.getLogger(n).getEffectiveLevel() for n in LOGGERS}
    try:
        with ProcessPoolExecutor(workers, mp_context=ctx, initializer=_worker_logging, initargs=(q, levels)) as ex:
            yield ex
    finally:
        _stop_listener(listener, q)


def _stop_listener(listener, q) -> None:
    """QueueListener.stop without waiting forever: the sentinel is queued, the thread (a daemon) given LISTENER_STOP seconds."""
    t = listener._thread
    try:
        listener.enqueue_sentinel()
    except Exception:                    # a broken queue: nothing more can come
        pass
    if t is not None:
        t.join(LISTENER_STOP)
        if t.is_alive():
            log.warning(f"worker log forwarding did not stop within {LISTENER_STOP} s (a worker died while sending a "
                        f"record); some worker log lines may be missing")
    listener._thread = None
    q.cancel_join_thread()               # never wait at exit for a feeder blocked on that queue
