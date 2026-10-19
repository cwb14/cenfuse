"""Shared helpers: errors, logging, intervals, the Poisson-binomial tail, fingerprints."""
from __future__ import annotations

import hashlib
import json
import logging
import sys
from pathlib import Path

import numpy as np

log = logging.getLogger("cenfuse")
ENGINE_LOGGERS = ("cenfuse", "cenrep", "cenhic", "cenmeth")


class CenfuseError(Exception):
    """Bad input or a failed precondition: reported as 'cenfuse: error: ...' (exit 2)."""


def setup_logging(verbose: bool) -> None:
    """One handler for cenfuse and the three engines; messages keep their engine's name."""
    h = logging.StreamHandler(sys.stderr)
    h.setFormatter(logging.Formatter("[%(name)s %(asctime)s] %(levelname)s: %(message)s", "%H:%M:%S"))
    for name in ENGINE_LOGGERS:
        lg = logging.getLogger(name)
        lg.handlers[:] = [h]
        lg.setLevel(logging.DEBUG if verbose else logging.INFO)
        lg.propagate = False


def gap(a0: int, a1: int, b0: int, b1: int) -> int:
    """bp between two half-open intervals (0 when they overlap or touch)."""
    return max(0, max(a0, b0) - min(a1, b1))


def poisson_binomial_sf(k: int, p) -> float:
    """P(K >= k) for K = a sum of independent Bernoulli(p_i) (exact, by dynamic programming)."""
    p = np.clip(np.asarray(p, float), 0.0, 1.0)
    if k <= 0:
        return 1.0
    if k > len(p):
        return 0.0
    dist = np.zeros(len(p) + 1)
    dist[0] = 1.0
    for q in p:
        dist[1:] = dist[1:] * (1 - q) + dist[:-1] * q
        dist[0] *= 1 - q
    return float(min(1.0, dist[k:].sum()))


def mcnemar_exact(b: int, c: int) -> float:
    """Two-sided exact McNemar P for discordant pair counts b and c (binomial, p = 0.5)."""
    from scipy.stats import binomtest
    n = b + c
    return 1.0 if n == 0 else float(binomtest(min(b, c), n, 0.5).pvalue)


def file_fingerprint(path, nbytes: int = 1 << 20) -> str:
    """Fast fingerprint (name, size, mtime, first and last MiB), as the engines use for their caches."""
    p = Path(path)
    st = p.stat()
    h = hashlib.md5(f"{p.resolve()}:{st.st_size}:{st.st_mtime_ns}".encode())
    with open(p, "rb") as f:
        h.update(f.read(nbytes))
        if st.st_size > 2 * nbytes:
            f.seek(-nbytes, 2)
            h.update(f.read())
    return h.hexdigest()[:16]


def obj_fingerprint(obj) -> str:
    return hashlib.md5(json.dumps(obj, sort_keys=True, default=str).encode()).hexdigest()[:16]


def level_shift(level: str, steps: int, levels=("low", "medium", "high")) -> str:
    """Move a confidence level by steps (clipped to the scale); 'none' stays 'none'."""
    if level not in levels:
        return level
    return levels[int(np.clip(levels.index(level) + steps, 0, len(levels) - 1))]
