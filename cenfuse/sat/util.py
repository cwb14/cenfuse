"""Shared helpers: logging, subprocesses, sequences, FASTA text, intervals."""
import gzip
import hashlib
import logging
import shutil
import subprocess
import threading
import time
import zlib
from contextlib import contextmanager

from ..runtime import atomic_open

log = logging.getLogger("cenfuse-sat")


class SatelliteError(RuntimeError):
    """Fatal, user-facing error (bad input, missing tool, failed command)."""


@contextmanager
def open_text(path):
    """A plain or gzip (.gz) text input. A file that cannot be read to its end (a cut download, binary data) stops with
    an error naming it."""
    try:
        with (gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)) as f:
            yield f
    except (EOFError, OSError, UnicodeDecodeError, zlib.error) as e:
        raise SatelliteError(f"{path}: unreadable or truncated ({e})") from None


def setup_logging(verbose: bool) -> None:
    h = logging.StreamHandler()
    h.setFormatter(logging.Formatter("%(asctime)s %(levelname)s %(message)s", "%Y-%m-%d %H:%M:%S"))
    log.handlers[:] = [h]
    log.setLevel(logging.DEBUG if verbose else logging.INFO)
    log.propagate = False


def fmt_duration(seconds: float) -> str:
    """'45 s', '12 min', '3 h 05 min', '2 d 4 h': two units at most, rounded to the smaller one, so a time left is not
    stated more finely than it can be estimated."""
    s = max(0, int(seconds + 0.5))
    if s < 60:
        return f"{s} s"
    m = (s + 30) // 60
    if m < 60:
        return f"{m} min"
    if m < 24 * 60:
        return f"{m // 60} h {m % 60:02d} min"
    d, h = divmod((s + 1800) // 3600, 24)
    return f"{d} d {h} h"


_FLOAT_NOISE = 1e-12          # 3 x 0.1 != 0.3, and sums of float work are off by a few ulps: a tenth is not missed


class Progress:
    """INFO lines for a long job whose work is done in pieces on many threads: one line each time another `step` of
    the work is done, and a heartbeat line when `heartbeat` seconds (default: an hour) pass without one: due at least
    hourly, but checked only when a piece of work finishes (add), so while every thread is busy with a long piece the
    log can stay quiet for longer. Each line gives the elapsed time and the time left: the work left is assumed to
    cost per unit what the work done did, which jobs of unequal cost only approximate, so the line says it is an
    estimate.
    The caller chooses the unit of work and adds it as each piece finishes, in whatever order they finish. A total
    of 0 (nothing to do) is silent."""

    def __init__(self, label: str, total: float, step: float = 0.1, heartbeat: float = 3600.0, clock=time.monotonic):
        self.label, self.total, self.step, self.heartbeat, self.clock = label, total, step, heartbeat, clock
        self.done, self.reached, self.finished = 0, 0, False        # work done; steps already reported; 100% reported
        self.start = self.last = clock()                            # start of the job; time of the last line
        self.lock = threading.Lock()

    def add(self, work: float) -> None:
        """Count `work` (in the unit of `total`) as done; log a line if a step or the heartbeat is due."""
        if self.total <= 0 or work <= 0:
            return
        with self.lock:                                             # lines are logged inside it, so they stay in order
            if self.finished:
                return
            self.done += work
            now, frac = self.clock(), min(1.0, self.done / self.total)
            steps = int(frac / self.step + _FLOAT_NOISE)
            self.finished = frac >= 1 - _FLOAT_NOISE
            if steps <= self.reached and not self.finished and now - self.last < self.heartbeat:
                return
            self.reached, self.last = steps, now
            elapsed = now - self.start
            log.info(f"{self.label}: {100 * frac:.0f}% of the work done in {fmt_duration(elapsed)}; about "
                     f"{fmt_duration(elapsed * (1 - frac) / frac)} left (estimate)")


def require_tools(names) -> None:
    missing = [n for n in names if shutil.which(n) is None]
    if missing:
        raise SatelliteError("required tools not on PATH: " + ", ".join(missing))


def run(cmd, **kw) -> subprocess.CompletedProcess:
    cmd = [str(c) for c in cmd]
    p = subprocess.run(cmd, capture_output=True, text=True, **kw)
    if p.returncode != 0:
        raise SatelliteError(f"command failed (exit {p.returncode}): {' '.join(cmd)}\n{p.stderr[-2000:]}")
    return p


def md5_file(path, bufsize: int = 1 << 24) -> str:
    h = hashlib.md5()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(bufsize), b""):
            h.update(b)
    return h.hexdigest()


_RC = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(s: str) -> str:
    return s.translate(_RC)[::-1]


def canonical_rotation(seq: str) -> str:
    """Lexicographically smallest rotation of a circular monomer (upper-cased).
    It starts with the smallest k-mer, so it is the spec's 'smallest 16-mer' rule with ties resolved."""
    s = seq.upper()
    return min((s[i:] + s[:i] for i in range(len(s))), default=s)


def write_fasta(path, records) -> None:
    with atomic_open(path) as f:
        for name, seq in records:
            f.write(f">{name}\n{seq}\n")


def parse_fasta_text(text: str) -> list:
    out, name, buf = [], None, []
    for line in text.splitlines():
        if line.startswith(">"):
            if name is not None:
                out.append((name, "".join(buf)))
            name, buf = line[1:].split()[0], []
        elif line.strip():
            buf.append(line.strip())
    if name is not None:
        out.append((name, "".join(buf)))
    return out


def merge_intervals(iv, gap: int = 0) -> list:
    """Merge 0-based half-open intervals separated by <= gap bp. Returns sorted [start, end] lists."""
    out = []
    for s, e in sorted(iv):
        if out and s - out[-1][1] <= gap:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return out


def union_bp(iv) -> int:
    return sum(e - s for s, e in merge_intervals(iv))


def clipped_bp(iv, lo: int, hi: int) -> int:
    """Union bp of intervals clipped to [lo, hi)."""
    return union_bp([(max(s, lo), min(e, hi)) for s, e in iv if min(e, hi) > max(s, lo)])


def overlap(a0: int, a1: int, b0: int, b1: int) -> int:
    return max(0, min(a1, b1) - max(a0, b0))


