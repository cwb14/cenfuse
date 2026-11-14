import logging
import signal

import numpy as np
import pytest

from cenfuse.sat import fastan
from cenfuse.sat.util import SatelliteError, log
from sat_common import crashing_fastan, dna, have

need_fastan = pytest.mark.skipif(not have("FasTAN", "FAtoGDB", "ONEview"), reason="FasTAN tools missing")

ANO = """g
S 1 c
C 5500

M 0 1999 3352
L 3 150
X 0
P 10 70 220 370 520 670 820 970 1120 1270 1420
M 0 4000 5700
L 3 848
X 2
"""


def rand(rng, n):
    return "".join(rng.choice(list("ACGT"), n))


def indel_copy(rng, s, rate):
    """A copy of s with substitutions, deletions and insertions (70/15/15) at total rate."""
    out = []
    for ch in s:
        x = rng.random()
        if x < 0.7 * rate:
            out.append(rng.choice([b for b in "ACGT" if b != ch]))
        elif x < 0.85 * rate:
            continue
        elif x < rate:
            out.append(ch + rng.choice(list("ACGT")))
        else:
            out.append(ch)
    return "".join(out)


def edit_in(a, t):
    """Edit distance of all of a against the best-matching substring of t (test-only reference)."""
    prev = [0] * (len(t) + 1)
    for i, ca in enumerate(a, 1):
        cur = [i] + [0] * len(t)
        for j, ct in enumerate(t, 1):
            cur[j] = min(prev[j - 1] + (ca != ct), prev[j] + 1, cur[j - 1] + 1)
        prev = cur
    return min(prev)


def test_parse_ano_reads_intervals_units_divergence_and_cuts():
    assert fastan.parse_ano(ANO) == [
        (1999, 3352, 150, 0, [70, 220, 370, 520, 670, 820, 970, 1120, 1270, 1420]),
        (4000, 5700, 848, 2, []),                 # a near tandem: no partition line, no cuts of its own
    ]


def test_unparsable_fastan_output_fails_clearly():
    from cenfuse.sat.util import SatelliteError
    with pytest.raises(SatelliteError, match="M 0 1999"):
        fastan.parse_ano("M 0 1999\nL 3 150\n")


def test_consensus_recovers_the_unit_from_divergent_copies():
    rng = np.random.default_rng(3)
    unit = rand(rng, 171)
    copies = [indel_copy(rng, unit, 0.10) for _ in range(40)]
    assert min(edit_in(c, unit) for c in copies[:5]) >= 5          # no single copy is the unit
    assert edit_in(fastan.unit_consensus(copies, 171), unit) <= 1


def test_consensus_leaves_out_misparsed_copies():
    """A missed cut gives a double copy; it must not stretch the consensus to two units."""
    rng = np.random.default_rng(4)
    unit = rand(rng, 120)
    good = [indel_copy(rng, unit, 0.05) for _ in range(10)]
    double = [indel_copy(rng, unit, 0.05) + indel_copy(rng, unit, 0.05) for _ in range(25)]
    c = fastan.unit_consensus(good + double, 120)
    assert abs(len(c) - 120) <= 3 and edit_in(c, unit) <= 2


def test_consensus_of_one_copy_is_that_copy_and_of_none_is_empty():
    assert fastan.unit_consensus(["ACGTTGCA", "A" * 30], 8) == "ACGTTGCA"
    assert fastan.unit_consensus(["A" * 30], 8) == ""


@need_fastan
def test_planted_array_has_its_coordinates_period_and_unit(tmp_path):
    """Whatever the FasTAN build (the 0.8 release ends intervals one unit early), the row spans the planted array in
    chromosome coordinates and its consensus is the planted unit."""
    rng = np.random.default_rng(5)
    unit = rand(rng, 171)
    seq = rand(rng, 3000) + "".join(indel_copy(rng, unit, 0.05) for _ in range(30)) + rand(rng, 3000)
    end = len(seq) - 3000
    rows = fastan.chunk_rows("c1", seq, 10_000, tmp_path, fastan.trims_last_unit(tmp_path))
    hit = [r for r in rows if abs(r[3] - 171) <= 5]
    assert len(hit) == 1
    chrom, s, e, period, copies, csize, pct, _, _, _, cons = hit[0]
    assert chrom == "c1" and abs(s - 13_000) <= 20 and abs(e - (10_000 + end)) <= 20
    assert abs(copies - (e - s) / period) < 0.1 and csize == len(cons) and pct >= 85
    assert edit_in(cons, unit + unit) <= 3


@need_fastan
def test_added_last_copy_never_crosses_an_assembly_gap(tmp_path):
    """The last copies (90 bp) are shorter than the mean unit (~98 bp), so adding one mean unit to the reported end
    overshoots the array; the array ends where the gap starts."""
    rng = np.random.default_rng(6)
    unit = rand(rng, 90)
    left = rand(rng, 2000) + "".join(unit[:45] + rand(rng, 10) + unit[45:] for _ in range(15)) + unit * 5
    seq = left + "N" * 100 + rand(rng, 2000)
    hit = [r for r in fastan.chunk_rows("c1", seq, 0, tmp_path, True) if 85 <= r[3] <= 105]
    assert hit and max(r[2] for r in hit) == len(left)


@need_fastan
def test_units_beyond_the_seed_range_are_left_to_the_long_unit_detector(tmp_path):
    rng = np.random.default_rng(7)
    unit = rand(rng, 2500)
    seq = rand(rng, 3000) + "".join(indel_copy(rng, unit, 0.02) for _ in range(4)) + rand(rng, 3000)
    raw = fastan.parse_ano(fastan.run_fastan(seq, tmp_path))
    assert any(abs(u - 2500) <= 50 for _, _, u, _, _ in raw)                  # FasTAN reports it ...
    assert all(r[3] <= 2000 for r in fastan.chunk_rows("c1", seq, 0, tmp_path, False))   # ... the engine drops it


@need_fastan
def test_sequence_without_bases_gives_no_rows(tmp_path):
    assert fastan.chunk_rows("c1", "N" * 5000, 0, tmp_path, False) == []
    assert fastan.chunk_rows("c1", "", 0, tmp_path, False) == []


def test_fasTAN_runs_with_a_fixed_memory_layout_where_possible():
    """FasTAN 0.8's copy cuts depend on the memory layout; on Linux with setarch it runs with randomization off."""
    import platform
    import shutil
    import sys
    pre = fastan.fixed_layout()
    if sys.platform.startswith("linux") and shutil.which("setarch"):
        assert pre == ["setarch", platform.machine(), "-R"]
    else:
        assert pre == []


def test_without_setarch_fastan_runs_as_is_and_warns(monkeypatch, caplog):
    """Without setarch (macOS, or no util-linux) FasTAN runs with the default memory layout, and cenfuse warns that
    satellite families can then differ between reruns."""
    import logging
    import shutil
    monkeypatch.setattr(fastan, "_FIXED", None)
    monkeypatch.setattr(shutil, "which", lambda name, *a, **k: None)
    with caplog.at_level(logging.WARNING, logger="cenfuse-sat"):
        assert fastan.fixed_layout() == []
    assert "setarch -R" in caplog.text


@need_fastan
def test_reruns_give_the_same_cuts(tmp_path):
    """Imperfect arrays of several units between random DNA: three runs, one output (layout-dependent cuts were seen
    on ~1% of real arrays)."""
    rng = np.random.default_rng(5)
    parts = [rand(rng, 3000)]
    for unit_len, n, rate in ((75, 4, 0.12), (170, 12, 0.08), (1200, 3, 0.10), (60, 30, 0.15)):
        unit = rand(rng, unit_len)
        parts += ["".join(indel_copy(rng, unit, rate) for _ in range(n)), rand(rng, 3000)]
    seq = "".join(parts)
    outs = {fastan.run_fastan(seq, tmp_path) for _ in range(3)}
    assert len(outs) == 1


# ---- crash recovery: FasTAN 0.8 aborts (SIGABRT) on some large satellite arrays; chunk_rows splits the piece and
# retries (these tests replace run_fastan by sat_common.crashing_fastan: no FasTAN needed)

def test_crash_on_a_long_piece_is_split_and_recovered(monkeypatch, tmp_path, caplog):
    """A piece FasTAN aborts on is cut in two halves that overlap by `overlap` bp, each retried. Here FasTAN dies on
    any piece over 1 Mb: the pieces it finishes tile the chunk, rows keep absolute coordinates, and every crash is
    logged and recorded."""
    monkeypatch.setattr(log, "propagate", True)
    seq = dna(np.random.default_rng(1), 4_000_000)
    done = crashing_fastan(monkeypatch, {"c1": seq}, lambda name, a, b: b - a > 1_000_000)
    events = []
    with caplog.at_level(logging.INFO, logger="cenfuse-sat"):
        rows = fastan.chunk_rows("c1", seq, 7_000_000, tmp_path, False, events, 10_000)
    pieces = sorted((a, b) for _, a, b in done)
    assert pieces[0][0] == 0 and pieces[-1][1] == 4_000_000                          # no gap at either end
    assert all(nxt[0] == cur[1] - 10_000 for cur, nxt in zip(pieces, pieces[1:]))    # neighbors overlap by `overlap`
    assert sorted((r[1], r[2]) for r in rows) == [(7_000_000 + a + 100, 7_000_000 + a + 300) for a, _ in pieces]
    assert len(events) == len(pieces) - 1                                            # every crash made two pieces
    assert events[0] == ("c1", 7_000_000, 11_000_000, "SIGABRT", "split")
    assert events[1] == ("c1", 7_000_000, 9_005_000, "SIGABRT", "split")             # the cut is at 9,000,000: both
    assert ("c1", 8_995_000, 11_000_000, "SIGABRT", "split") in events              # halves reach 5,000 bp past it
    assert all(e[3:] == ("SIGABRT", "split") and 7_000_000 <= e[1] < e[2] <= 11_000_000 for e in events)
    info = [r.getMessage() for r in caplog.records if r.levelno == logging.INFO]
    assert len(info) == len(events) and "SIGABRT" in info[0] and "c1:7,000,000-11,000,000" in info[0]


def test_a_piece_that_always_crashes_is_skipped(monkeypatch, tmp_path, caplog):
    """FasTAN dies on every piece that holds one position (7,500,000, 2 kb before the chunk's end). A piece down to
    2 x overlap is skipped, with a warning that names its interval; the rest of the chunk is still searched."""
    monkeypatch.setattr(log, "propagate", True)
    seq = dna(np.random.default_rng(2), 502_000)
    done = crashing_fastan(monkeypatch, {"c1": seq}, lambda name, a, b: a <= 500_000 < b)
    events = []
    with caplog.at_level(logging.INFO, logger="cenfuse-sat"):
        rows = fastan.chunk_rows("c1", seq, 7_000_000, tmp_path, False, events, 10_000)
    skipped = [e for e in events if e[4] == "skipped"]
    assert len(skipped) == 1
    chrom, start, end, signame, _ = skipped[0]
    assert signame == "SIGABRT" and end - start <= 2 * 10_000 and start <= 7_500_000 < end
    warnings = [r.getMessage() for r in caplog.records if r.levelno == logging.WARNING]
    assert len(warnings) == 1 and f"c1:{start:,}-{end:,}" in warnings[0]
    searched = np.zeros(len(seq), bool)
    for _, a, b in done:
        searched[a:b] = True
    lost = np.flatnonzero(~searched) + 7_000_000
    assert lost.min() >= start and lost.max() < end                                  # only the skipped piece is lost
    assert len(rows) == len(done) > 0


@pytest.mark.parametrize("n, actions", [(20_000, ["skipped"]), (20_001, ["split", "skipped", "skipped"])])
def test_a_piece_is_split_only_while_longer_than_twice_the_overlap(monkeypatch, tmp_path, n, actions):
    """FasTAN dies on every piece. At exactly 2 x overlap the chunk is skipped as it is; one base more and it is cut
    in two halves (15,000 and 15,001 bp), which are then skipped, and since no piece of a chunk longer than
    2 x overlap was searched, the run stops (ChunkFailed, below)."""
    seq = dna(np.random.default_rng(4), n)
    crashing_fastan(monkeypatch, {"c1": seq}, lambda name, a, b: True)
    events = []
    if n > 2 * 10_000:
        with pytest.raises(fastan.ChunkFailed):
            fastan.chunk_rows("c1", seq, 0, tmp_path, False, events, 10_000)
    else:
        assert fastan.chunk_rows("c1", seq, 0, tmp_path, False, events, 10_000) == []
    assert [e[4] for e in events] == actions


# ---- a chunk FasTAN searched no piece of stops the run: its known crash is local to one array, so a whole chunk
# failing means FasTAN itself is failing, and a run without the chunk's seeds could read as a satellite-free genome

def _searched(done: list, n: int) -> np.ndarray:
    covered = np.zeros(n, bool)
    for _, a, b in done:
        covered[a:b] = True
    return covered


def test_a_chunk_where_every_piece_crashes_stops_the_run(monkeypatch, tmp_path, caplog):
    """FasTAN dies on every piece of a 200-kb chunk at 7 Mb: every piece is split down to the skip floor and skipped,
    so none of the chunk was searched. Instead of going on without its seeds, the run stops with an error that names
    the chunk and suggests TRF; the crashes up to then are recorded as usual (the skipped pieces cover the chunk)."""
    monkeypatch.setattr(log, "propagate", True)
    seq = dna(np.random.default_rng(5), 200_000)
    done = crashing_fastan(monkeypatch, {"c1": seq}, lambda name, a, b: True)
    events = []
    with caplog.at_level(logging.INFO, logger="cenfuse-sat"):
        with pytest.raises(SatelliteError, match=r"every piece of c1:7,000,000-7,200,000") as err:
            fastan.chunk_rows("c1", seq, 7_000_000, tmp_path, False, events, 10_000)
    assert isinstance(err.value, fastan.ChunkFailed) and "SIGABRT" in str(err.value)
    assert "--seed-engine trf" in str(err.value)
    assert done == [] and events[0] == ("c1", 7_000_000, 7_200_000, "SIGABRT", "split")
    skipped = [(a - 7_000_000, b - 7_000_000) for _, a, b, _, action in events if action == "skipped"]
    assert all(b - a <= 2 * 10_000 for a, b in skipped)
    assert _searched([("c1", a, b) for a, b in skipped], 200_000).all()
    assert len([r for r in caplog.records if r.levelno == logging.WARNING]) == len(skipped)


def test_one_crashing_piece_among_searched_ones_gives_rows_as_before(monkeypatch, tmp_path):
    """FasTAN dies only on the pieces that hold one position, the middle of a 200-kb chunk (so both halves of the
    first cut hold it): the pieces at the skip floor that hold it are skipped and the rest of the chunk is searched,
    so the chunk gives its rows and the run goes on."""
    seq = dna(np.random.default_rng(6), 200_000)
    done = crashing_fastan(monkeypatch, {"c1": seq}, lambda name, a, b: a <= 100_000 < b)
    events = []
    rows = fastan.chunk_rows("c1", seq, 7_000_000, tmp_path, False, events, 10_000)
    assert len(rows) == len(done) > 0
    skipped = [("c1", a - 7_000_000, b - 7_000_000) for _, a, b, _, action in events if action == "skipped"]
    assert skipped and all(a <= 100_000 < b and b - a <= 2 * 10_000 for _, a, b in skipped)
    assert not (~_searched(done, 200_000) & ~_searched(skipped, 200_000)).any()     # only skipped pieces are lost


def test_a_short_crashing_sequence_is_skipped_with_a_warning(monkeypatch, tmp_path, caplog):
    """A lone sequence (or chunk) of at most 2 x overlap that crashes is skipped with a warning, as before: it is too
    short to split, and a crash on one short piece says nothing about the rest of the genome."""
    monkeypatch.setattr(log, "propagate", True)
    seq = dna(np.random.default_rng(7), 15_000)
    crashing_fastan(monkeypatch, {"c1": seq}, lambda name, a, b: True)
    events = []
    with caplog.at_level(logging.WARNING, logger="cenfuse-sat"):
        assert fastan.chunk_rows("c1", seq, 0, tmp_path, False, events, 10_000) == []
    assert events == [("c1", 0, 15_000, "SIGABRT", "skipped")]
    assert [r.getMessage() for r in caplog.records if "c1:0-15,000" in r.getMessage()]


def test_a_mostly_n_chunk_whose_only_array_crashes_is_skipped_not_stopped(monkeypatch, tmp_path, caplog):
    """A 300-kb chunk that is N apart from its first 15 kb, where FasTAN dies: a crash local to one array, the case
    the stop must not catch. A piece of N only has nothing to search and counts as searched, so the pieces with
    sequence are skipped with a warning each and the run goes on, without seeds from them."""
    monkeypatch.setattr(log, "propagate", True)
    seq = dna(np.random.default_rng(8), 15_000) + "N" * 285_000
    done = crashing_fastan(monkeypatch, {"c1": seq}, lambda name, a, b: True)
    events = []
    with caplog.at_level(logging.WARNING, logger="cenfuse-sat"):
        assert fastan.chunk_rows("c1", seq, 0, tmp_path, False, events, 10_000) == []
    skipped = [(a, b) for _, a, b, _, action in events if action == "skipped"]
    assert done == [] and skipped and all(b - a <= 2 * 10_000 and a < 15_000 for a, b in skipped)
    assert len([r for r in caplog.records if r.levelno == logging.WARNING]) == len(skipped)


@pytest.mark.parametrize("code", [1, -signal.SIGKILL, -signal.SIGTERM, -signal.SIGINT])
def test_non_crash_failures_still_stop_the_run(monkeypatch, tmp_path, code):
    """Only a crash is routed around. An exit status, or a kill (the OOM killer, the scheduler, Ctrl-C), stops the
    run at the first piece: splitting would hide a real problem and keep a killed job going."""
    calls = []

    def run(piece, tmpdir):
        calls.append(len(piece))
        raise SatelliteError(f"command failed (exit {code}): FasTAN -m -p -T1 s.1gdb")

    monkeypatch.setattr(fastan, "run_fastan", run)
    events = []
    with pytest.raises(SatelliteError, match=f"exit {code}"):
        fastan.chunk_rows("c1", "ACGT" * 50_000, 0, tmp_path, False, events, 10_000)
    assert calls == [200_000] and events == []


@pytest.mark.parametrize("msg, want", [
    *[pytest.param(f"command failed (exit -{s.value}): FasTAN -m -p -T1 s.1gdb\nmalloc(): corrupted top size", s.value,
                   id=s.name) for s in (signal.SIGABRT, signal.SIGSEGV, signal.SIGBUS, signal.SIGFPE, signal.SIGILL)],
    pytest.param("command failed (exit -9): FasTAN -m -p -T1 s.1gdb", None, id="SIGKILL-the-OOM-killer"),
    pytest.param("command failed (exit -15): FasTAN -m -p -T1 s.1gdb", None, id="SIGTERM"),
    pytest.param("command failed (exit -2): FasTAN -m -p -T1 s.1gdb", None, id="SIGINT"),
    pytest.param("command failed (exit 1): FasTAN -m -p -T1 s.1gdb", None, id="exit-status-1"),
    pytest.param("command failed (exit 6): FasTAN -m -p -T1 s.1gdb", None, id="exit-status-6-is-not-signal-6"),
    pytest.param("command failed (exit 1): FAtoGDB s.fa\ncommand failed (exit -6): x", None, id="opening-status-only"),
    pytest.param("FasTAN output: unparsable .1ano line 'M 0 1999'", None, id="other-error"),
    pytest.param("command failed (exit x): FasTAN", None, id="unparsable-status"),
    pytest.param("", None, id="empty")])
def test_crash_signal_parses_only_crash_signals(msg, want):
    """SIGABRT, SIGSEGV, SIGBUS, SIGFPE and SIGILL mean the program itself failed (signal numbers differ by platform,
    hence the names); every other failure, kills included, gives None."""
    assert fastan.crash_signal(SatelliteError(msg)) == want


def test_crash_at_sequence_edges_and_in_mostly_n_pieces(monkeypatch, tmp_path):
    """FasTAN dies on whatever piece holds the first or the last base of a 400-kb sequence whose middle is N (a
    chromosome end beside a gap). The recursion ends at both edges: exactly those two pieces are skipped, no piece
    without bases ever reaches FasTAN (FAtoGDB would reject it), and every base outside them is searched."""
    rng = np.random.default_rng(3)
    seq = dna(rng, 20_000) + "N" * 360_000 + dna(rng, 20_000)
    n = len(seq)
    done = crashing_fastan(monkeypatch, {"c1": seq}, lambda name, a, b: a <= 0 < b or a <= n - 1 < b)
    events = []
    rows = fastan.chunk_rows("c1", seq, 0, tmp_path, False, events, 10_000)
    skipped = sorted(e for e in events if e[4] == "skipped")
    assert len(skipped) == 2 and skipped[0][1] == 0 and skipped[1][2] == n          # the first and the last piece only
    assert all(e[2] - e[1] <= 2 * 10_000 for e in skipped)
    bases = np.ones(n, bool)
    bases[20_000:380_000] = False
    searched, skip_mask = np.zeros(n, bool), np.zeros(n, bool)
    for _, a, b in done:
        searched[a:b] = True
    for _, a, b, _, _ in skipped:
        skip_mask[a:b] = True
    assert not (bases & ~searched & ~skip_mask).any()                                # nothing else is lost
    assert len(rows) == len(done) > 0
