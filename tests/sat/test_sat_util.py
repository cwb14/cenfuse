import logging
import re
import threading
import time

import pytest

from cenfuse.sat.util import (Progress, SatelliteError, canonical_rotation, clipped_bp, fmt_duration, log,
                              merge_intervals, parse_fasta_text, revcomp, run, union_bp)


def test_canonical_rotation_is_rotation_invariant():
    s = "GATTACAGATTTCC"
    rots = {canonical_rotation(s[i:] + s[:i]) for i in range(len(s))}
    assert len(rots) == 1 and len(next(iter(rots))) == len(s)


def test_merge_intervals_gap():
    assert merge_intervals([(10, 20), (0, 5), (22, 30)], gap=2) == [[0, 5], [10, 30]]
    assert merge_intervals([(0, 5), (5, 8)]) == [[0, 8]]


def test_union_and_clipped_bp():
    iv = [(0, 10), (5, 15), (20, 30)]
    assert union_bp(iv) == 25
    assert clipped_bp(iv, 8, 25) == 12


def test_revcomp_and_fasta_parse():
    assert revcomp("AACGTn") == "nACGTT"
    assert parse_fasta_text(">a x\nAC\nGT\n>b\nTT\n") == [("a", "ACGT"), ("b", "TT")]


def test_run_raises_on_failure():
    with pytest.raises(SatelliteError):
        run(["false"])


# ---- Progress: INFO lines with an estimate for long jobs (the genome-wide search showed none for up to 25 h)

class _Clock:
    """A clock the test sets by hand (Progress takes any callable that returns seconds)."""

    def __init__(self):
        self.t = 0.0

    def __call__(self):
        return self.t


@pytest.fixture
def lines(caplog, monkeypatch):
    """The INFO messages of the satellite engine's logger (it does not propagate to caplog unless told to)."""
    monkeypatch.setattr(log, "propagate", True)
    caplog.set_level(logging.INFO, logger="cenfuse-sat")
    return lambda: [r.getMessage() for r in caplog.records if r.name == "cenfuse-sat" and r.levelno == logging.INFO]


def _pcts(messages):
    """Percentages of the progress lines; each line must have exactly the documented form."""
    form = re.compile(r"search: (\d+)% of the work done in [^;]+; about [^;]+ left \(estimate\)")
    assert all(form.fullmatch(m) for m in messages), messages
    return [int(form.fullmatch(m).group(1)) for m in messages]


@pytest.mark.parametrize("total, work", [(100, 5), (1.0, 0.05), (7.0, 0.35)])      # float sums: 20 x 0.05 != 1.0
def test_progress_logs_each_tenth_once(lines, total, work):
    """One line per tenth of the work, whatever the unit and however the float sums round, and none after 100%."""
    clock = _Clock()
    prog = Progress("search", total, clock=clock)
    for _ in range(20):
        clock.t += 60
        prog.add(work)
    assert _pcts(lines()) == list(range(10, 101, 10))
    clock.t += 7200                                  # finished: not even a heartbeat after the last line
    prog.add(work)
    assert len(lines()) == 10


def test_progress_reports_a_jump_over_several_tenths_once(lines):
    """A piece of work larger than a tenth is not reported tenth by tenth."""
    clock = _Clock()
    prog = Progress("search", 100, clock=clock)
    prog.add(35)                                     # three tenths at once: one line, at the true percentage
    prog.add(4)                                      # 39%: the fourth tenth is not reached yet
    prog.add(1)
    assert _pcts(lines()) == [35, 40]


def test_progress_heartbeat_without_a_milestone(lines):
    """When no tenth is reached for an hour, one line says where the job is; the next one comes an hour later."""
    clock = _Clock()
    prog = Progress("search", 100, clock=clock)
    prog.add(1)                                      # t = 0: 1% done, no tenth crossed, nothing to say
    assert lines() == []
    clock.t = 3601                                   # an hour with no line: say where the search is
    prog.add(1)
    assert _pcts(lines()) == [2]
    clock.t = 7200                                   # the hour counts from the last line, not from the start
    prog.add(1)
    assert len(lines()) == 1
    clock.t = 7201
    prog.add(1)
    assert _pcts(lines()) == [2, 4]


def test_progress_estimate_of_time_left(lines):
    """The work left costs per unit what the work done cost: 25% in 100 s leaves 300 s. The line says it is an
    estimate."""
    clock = _Clock()
    prog = Progress("search", 100, clock=clock)
    clock.t = 100
    prog.add(25)
    assert "about 5 min left (estimate)" in lines()[0]
    clock.t = 3600                                   # exact form of a line, on round figures
    prog.add(25)                                     # 50%: 1 h done, 1 h left
    assert lines()[1] == "search: 50% of the work done in 1 h 00 min; about 1 h 00 min left (estimate)"


def test_progress_with_zero_total_is_silent(lines):
    """A run with no work (no families, no sequences) must not divide by zero or say anything."""
    clock = _Clock()
    prog = Progress("search", 0, clock=clock)
    prog.add(0)
    clock.t = 10_000
    prog.add(5)
    assert lines() == []


class _Yielding(int):
    """One unit of work whose addition gives up the GIL half-way, as an unlucky thread switch would."""

    def __radd__(self, other):
        time.sleep(0)
        return int(other) + int(self)


def test_progress_is_thread_safe(lines):
    """Searches finish on many worker threads: no work is lost and no tenth is reported twice or out of order. The
    work unit yields inside the addition, so a Progress without its lock loses updates here (and misses 100%)."""
    prog = Progress("search", 8000)
    start = threading.Barrier(8)

    def worker():
        start.wait()
        for _ in range(1000):
            prog.add(_Yielding(1))

    threads = [threading.Thread(target=worker) for _ in range(8)]
    for t in threads:
        t.start()
    for t in threads:
        t.join()
    assert _pcts(lines()) == list(range(10, 101, 10))


@pytest.mark.parametrize("seconds, text", [
    (0, "0 s"), (45, "45 s"), (59.4, "59 s"), (59.6, "1 min"), (100, "2 min"), (720, "12 min"),
    (3599, "1 h 00 min"), (3600, "1 h 00 min"), (11100, "3 h 05 min"), (86399, "1 d 0 h"), (187200, "2 d 4 h"),
    (-5, "0 s")])
def test_fmt_duration(seconds, text):
    """At most two units, rounded to the smaller one: an estimate is not stated more finely than it is known."""
    assert fmt_duration(seconds) == text
