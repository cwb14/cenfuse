"""S1 engine: FasTAN (Myers) tandem arrays, with each array's interval, mean unit, divergence of adjacent copies and
the cuts between copies. FasTAN reports no unit sequence, so a seed's consensus is the partial-order alignment (abPOA)
of its copies."""
import os
import platform
import random
import re
import shutil
import signal
import subprocess
import sys
import tempfile

import numpy as np
import pyabpoa

from .params import CONS_UNITS, FASTAN_ARGS, LEN_TOL, RELATED_LEN_TOL, SEED_MAX_PERIOD, Params
from .util import SatelliteError, log, run

TOOLS = ["FAtoGDB", "FasTAN", "ONEview"]
# FasTAN 0.8 aborts ("malloc(): corrupted top size", SIGABRT) on some large satellite arrays (Oryza minuta
# GWHESQK00000003:29,940,000-30,582,500, a 642-kb array). These signals mean the program itself failed, so the piece
# is split and retried; a kill (SIGKILL: the OOM killer; SIGTERM, SIGINT: the scheduler or the user) is not FasTAN's
# fault and must stop the run.
CRASH_SIGNALS = (signal.SIGABRT, signal.SIGSEGV, signal.SIGBUS, signal.SIGFPE, signal.SIGILL)
_KILLED_BY = re.compile(r"command failed \(exit -(\d+)\)")     # util.run: a command that died by signal N exits -N
DUP_TOL = LEN_TOL   # a unit is a mean over the copies aligned, so the same array seen from two chunks can differ
KEY = ["fastan", 2, FASTAN_ARGS, CONS_UNITS, RELATED_LEN_TOL, SEED_MAX_PERIOD, DUP_TOL]   # cache key; bump 2 if rows change
NAN = float("nan")
_POA = None
_FIXED = None     # the prefix that runs FasTAN with address-space randomization off (fixed_layout)


def fixed_layout() -> list:
    """['setarch', <arch>, '-R'] where it works, else []. FasTAN 0.8's copy cuts depend on the process's memory layout:
    with address-space randomization, ~1% of arrays get other cuts on each run (T2T moss, one 5-Mb chunk: 17 of 982,
    13 of them inside the array; four runs, four outputs; with randomization off, four identical outputs), so unit
    consensuses, families and calls shifted between runs. With a fixed layout reruns are identical on one system."""
    global _FIXED
    if _FIXED is None:
        _FIXED = []
        cmd = ["setarch", platform.machine(), "-R"]
        if sys.platform.startswith("linux") and shutil.which("setarch"):
            try:
                if subprocess.run(cmd + ["true"], capture_output=True).returncode == 0:
                    _FIXED = cmd
            except OSError:
                pass
        if not _FIXED:
            log.warning("seeds: cannot run FasTAN with address randomization off (setarch -R); FasTAN 0.8 then cuts "
                        "~1% of arrays differently on each run, so satellite families can differ between reruns")
    return _FIXED


def parse_ano(text: str) -> list:
    """ONEview text of a FasTAN .1ano -> [(beg, end, unit, pct_diff, cuts)], 0-based half-open; cuts are the copy
    boundaries relative to beg (P line), none for a near tandem (two copies with a short gap)."""
    rows = []
    for line in text.splitlines():
        t = line.split()
        if not t:
            continue
        try:
            if t[0] == "M":
                rows.append([int(t[2]), int(t[3]), 0, 0, []])
            elif rows and t[0] == "L":
                rows[-1][2] = int(t[2])
            elif rows and t[0] == "X":
                rows[-1][3] = int(t[1])
            elif rows and t[0] == "P":
                rows[-1][4] = [int(x) for x in t[2:2 + int(t[1])]]
        except (ValueError, IndexError):                 # a format change, not a bad record: stop
            raise SatelliteError(f"FasTAN output: unparsable .1ano line {line[:80]!r} (FasTAN/ONEview version?)") from None
    return [tuple(r) for r in rows]


def run_fastan(seq: str, tmpdir) -> str:
    """FasTAN on one sequence -> the ONEview text of its .1ano (FasTAN masks only from a GDB, hence FAtoGDB)."""
    work = tempfile.mkdtemp(dir=tmpdir)
    try:
        with open(os.path.join(work, "s.fa"), "w") as f:
            f.write(f">s\n{seq}\n")
        run(["FAtoGDB", "s.fa"], cwd=work)
        run([*fixed_layout(), "FasTAN", *FASTAN_ARGS, "s.1gdb"], cwd=work)
        return run(["ONEview", "-h", "s.1ano"], cwd=work).stdout
    finally:
        shutil.rmtree(work, ignore_errors=True)


def crash_signal(err: SatelliteError) -> int | None:
    """N if err says a command (FasTAN, FAtoGDB or ONEview) died by crash signal N (CRASH_SIGNALS), else None. Only the
    status that opens the message counts, not the command's stderr after it."""
    m = _KILLED_BY.match(str(err))
    n = int(m.group(1)) if m else None
    return n if n in CRASH_SIGNALS else None


def trims_last_unit(tmpdir) -> bool:
    """True if the installed FasTAN ends each interval where the lagging copy of its first-wave alignment ends, one
    unit before the array does: the 0.8 release does (fixed upstream on 2026-08-09). Probe: ten exact copies of a
    150-bp unit at [2000, 3500) between random flanks."""
    r = random.Random(1)
    rs = lambda n: "".join(r.choice("ACGT") for _ in range(n))
    unit = rs(150)
    ends = [e for _, e, u, _, _ in parse_ano(run_fastan(rs(2000) + unit * 10 + rs(2000), tmpdir)) if u == 150]
    if not ends:
        raise SatelliteError("FasTAN self-test found no array (10 copies of a 150-bp unit); check FasTAN, FAtoGDB and "
                          "ONEview (FASTGA)")
    return max(ends) < 3500 - 75


def setup(tmpdir) -> bool:
    trims = trims_last_unit(tmpdir)
    if trims:
        log.info("seeds: this FasTAN build ends each array one copy early (0.8 release); ends are extended by one unit")
    return trims


def unit_consensus(units: list, period: int) -> str:
    """Majority unit of an array: abPOA consensus of up to CONS_UNITS copies spread along it. Copies off the period by
    more than RELATED_LEN_TOL are misparsed (a missed or extra cut) and left out. One usable copy is returned as is;
    none gives ''."""
    global _POA
    ok = [u for u in units if abs(len(u) - period) <= RELATED_LEN_TOL * period]
    if len(ok) > CONS_UNITS:
        ok = [ok[i] for i in np.linspace(0, len(ok) - 1, CONS_UNITS).round().astype(int)]
    if len(ok) < 2:
        return ok[0] if ok else ""
    if _POA is None:                                 # most frequent base per column (gap-majority columns dropped):
        _POA = pyabpoa.msa_aligner(cons_algrm="MF")  # the majority rule of families.msa_consensus
    return _POA.msa(ok, out_cons=True, out_msa=False).cons_seq[0]


class ChunkFailed(SatelliteError):
    """FasTAN crashed on every piece of a chunk longer than twice the split overlap: no part of it was searched."""


def chunk_rows(chrom: str, seq: str, offset: int, tmpdir, trims: bool, events: list | None = None,
               overlap: int = Params.seed_overlap) -> list:
    """Seed-table rows (seeds.COLS) of the FasTAN arrays of one chunk at offset, units up to SEED_MAX_PERIOD (longer
    ones are the long-unit detector's). pct_match = 100 - FasTAN's % differences between adjacent copies; score = the
    matched bp; FasTAN measures no indel share or entropy (NaN). trims: see trims_last_unit.
    A piece FasTAN dies on by a crash signal is cut in two halves that overlap by `overlap` bp and each half is
    retried, recursively (seeds.dedupe merges what both halves report, as for chunk overlaps); with the overlap
    >= 2 x SEED_MAX_PERIOD no array of a seed-length unit is lost at the cut. A piece of at most 2 x overlap that still
    crashes is skipped with a warning. Each crash is appended to events as (chrom, start, end, signal name, 'split' |
    'skipped'), 0-based absolute. Any other failure (an exit status, a kill) stops the run, and so does a chunk longer
    than 2 x overlap of which no piece was searched (ChunkFailed): every piece with sequence crashed down to the skip
    floor. A piece of N only has nothing to search and counts as searched."""
    events = [] if events is None else events
    first, searched = len(events), []
    rows = _piece_rows(chrom, seq.upper(), offset, tmpdir, trims, events, overlap, searched)
    if not searched and len(events) > first and len(seq) > 2 * overlap:
        # FasTAN's known crash is local to one array (Oryza minuta lost 0 bp after 4 splits), so a chunk where every
        # piece with sequence crashed means FasTAN is broken (one crashing array beside N is skipped, not stopped)
        sigs = ", ".join(sorted({e[3] for e in events[first:]}))
        raise ChunkFailed(f"FasTAN crashed ({sigs}) on every piece of {chrom}:{offset:,}-{offset + len(seq):,}, down "
                          f"to {2 * overlap:,} bp, so none of it was searched: its known crash is local to one array, "
                          f"so FasTAN itself is failing; check it, or rerun with --seed-engine trf")
    return rows


def _piece_rows(chrom: str, seq: str, offset: int, tmpdir, trims: bool, events: list, overlap: int,
                searched: list) -> list:
    """chunk_rows on one piece of a chunk (seq upper-cased); the offset of each piece FasTAN finishes, or that has no
    bases to search, goes to searched."""
    if not any(b in seq for b in "ACGT"):
        searched.append(offset)                      # nothing to search, so nothing failed
        return []                                    # FAtoGDB rejects a sequence without bases
    try:
        ano = run_fastan(seq, tmpdir)
    except SatelliteError as e:
        sig = crash_signal(e)
        if sig is None:
            raise                                    # an exit status or a kill: not FasTAN's bug to route around
        name, stop = signal.Signals(sig).name, offset + len(seq)
        where = f"{chrom}:{offset:,}-{stop:,}"
        if len(seq) <= 2 * overlap:
            log.warning(f"seeds: FasTAN crashed ({name}) on {where} ({len(seq):,} bp), too short to split further: "
                        f"skipped, no seeds from this piece")
            events.append((chrom, offset, stop, name, "skipped"))
            return []
        log.info(f"seeds: FasTAN crashed ({name}) on {where}; retrying its two halves, overlapping by {overlap:,} bp")
        events.append((chrom, offset, stop, name, "split"))
        h, half = len(seq) // 2, overlap // 2
        return _piece_rows(chrom, seq[:h + half], offset, tmpdir, trims, events, overlap, searched) + \
            _piece_rows(chrom, seq[h - half:], offset + h - half, tmpdir, trims, events, overlap, searched)
    searched.append(offset)
    out = []
    for beg, end, unit, diff, cuts in parse_ano(ano):
        if unit > SEED_MAX_PERIOD:
            continue
        if trims:                                    # add the last copy; an array never crosses an assembly gap
            gap = seq.find("N", end, end + unit)
            end = min(len(seq), end + unit) if gap < 0 else gap
        cons = unit_consensus([seq[beg + a:beg + b] for a, b in zip(cuts, cuts[1:])], unit) or seq[beg:beg + unit]
        span = end - beg
        out.append((chrom, offset + beg, offset + end, unit, round(span / unit, 1), len(cons), 100 - diff, NAN,
                    round(span * (100 - diff) / 100), NAN, cons))
    return out
