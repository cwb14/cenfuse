"""The checkpoint layer (future.md #59, #60): content fingerprints, code keys, --force, the output-folder lock, the satellite
engine's input keys, the run manifest."""
import json
import os
import shutil
import socket
import time

import pysam
import pytest

from cenfuse import manifest, pipeline
from cenfuse import runtime as util
from cenfuse.pipeline import Run, _normalized, code_key
from cenfuse.runtime import file_fingerprint
from cenfuse.util import CenfuseError


def _fasta(d, name="g", n=2, L=20_000):
    d.mkdir(parents=True, exist_ok=True)
    f = d / f"{name}.fa"
    with open(f, "w") as h:
        for i in range(n):
            h.write(f">chr{i + 1}\n" + ("ACGTTGCA" * (L // 8)) + "\n")
    pysam.faidx(str(f))
    return f


# ------------------------------------------------------------------ fingerprints
def test_fingerprint_is_content_not_place_or_time(tmp_path):
    a = tmp_path / "a.txt"
    a.write_text("x" * 1000)
    b = tmp_path / "sub" / "b.txt"
    b.parent.mkdir()
    shutil.copy(a, b)                                     # cp without -p: new path, new mtime
    os.utime(b, (1, 1))
    assert file_fingerprint(a) == file_fingerprint(b)
    b.write_text("x" * 999 + "y")
    assert file_fingerprint(a) != file_fingerprint(b)


def test_fingerprint_of_a_large_file_samples_it(tmp_path, monkeypatch):
    monkeypatch.setattr(util, "WHOLE_FILE", 1000)
    a = tmp_path / "a.bin"
    a.write_bytes(bytes(5000))
    f0 = file_fingerprint(a, nbytes=100)
    a.write_bytes(bytes(4999) + b"\x01")                  # the last byte: inside the end sample
    assert file_fingerprint(a, nbytes=100) != f0


# ------------------------------------------------------------------ code keys
def test_normalized_source_ignores_comments_docstrings_and_layout(tmp_path):
    a, b, c = tmp_path / "a.py", tmp_path / "b.py", tmp_path / "c.py"
    a.write_text('"""Module doc."""\nX = 1  # a comment\n\n\ndef f(y):\n    """Doc."""\n    return y + X\n')
    b.write_text('"""Another doc."""\nX = 1\ndef f(y):\n    # a new comment\n    return y + X\n')
    c.write_text('X = 2\ndef f(y):\n    return y + X\n')
    assert _normalized(a) == _normalized(b) != _normalized(c)


def test_code_key_covers_the_import_closure():
    seen, todo = set(), [pipeline._PKG / "fuse.py"]
    while todo:
        f = todo.pop()
        if f not in seen:
            seen.add(f)
            todo += list(pipeline._imports(f))
    names = {str(f.relative_to(pipeline._PKG)) for f in seen}
    assert {"fuse.py", "hicfuse.py", "ageline.py", "sat/te.py", "params.py", "util.py"} <= names
    assert code_key("fuse.py") == code_key("fuse.py") and code_key("fuse.py") != code_key("telomere.py")


# ------------------------------------------------------------------ cenfuse's own checkpoints and --force
@pytest.mark.parametrize("force, hit", [(None, True), ("telomere", False), ("all", False), ("rdna", True)])
def test_force_recomputes_that_stage_or_all(tmp_path, force, hit):
    _fasta(tmp_path / "in")
    r = Run(tmp_path / "in", tmp_path / "out", force=force)
    r._store("telomere", "g", dict(k=1), {"v": 42})
    assert (r._cached("telomere", "g", dict(k=1)) == {"v": 42}) is hit
    assert r._cached("telomere", "g", dict(k=2)) is None             # another key: recomputed


def test_unreadable_checkpoint_is_recomputed_with_a_warning(tmp_path, caplog):
    _fasta(tmp_path / "in")
    r = Run(tmp_path / "in", tmp_path / "out")
    r._store("telomere", "g", dict(k=1), {"v": 42})
    (r.work / "telomere" / "g.pkl").write_bytes(b"not a pickle")
    with caplog.at_level("WARNING", logger="cenfuse"):
        assert r._cached("telomere", "g", dict(k=1)) is None
    assert "unreadable" in caplog.text
    assert json.loads((r.work / "telomere" / "g.done").read_text()) == dict(k=1)
    assert not list((r.work / "telomere").glob(".*tmp"))             # stamps are written atomically


# ------------------------------------------------------------------ the output-folder lock
def test_a_second_run_on_one_folder_stops(tmp_path):
    _fasta(tmp_path / "in")
    out = tmp_path / "out"
    out.mkdir()
    now = time.strftime("%Y-%m-%dT%H:%M:%S")                    # a live run: this process, which started before its lock
    (out / pipeline.LOCK).write_text(f"{socket.gethostname()} {os.getpid()} {now} cenfuse run x")
    with pytest.raises(CenfuseError, match="in use by another cenfuse run"):
        Run(tmp_path / "in", out).run()
    assert (out / pipeline.LOCK).exists()                            # another run's lock is never removed


def test_a_stale_lock_is_taken_over(tmp_path, monkeypatch, caplog):
    _fasta(tmp_path / "in")
    out = tmp_path / "out"
    out.mkdir()
    (out / pipeline.LOCK).write_text(f"{socket.gethostname()} 999999999 2026-01-01T00:00:00 cenfuse run x")   # dead pid here
    seen = {}

    def fake(self, t0):
        seen["lock"] = (out / pipeline.LOCK).read_text().split()[1]
        return []
    monkeypatch.setattr(Run, "_run", fake)
    with caplog.at_level("WARNING", logger="cenfuse"):
        Run(tmp_path / "in", out).run()
    assert "stale lock" in caplog.text and seen["lock"] == str(os.getpid())
    assert not (out / pipeline.LOCK).exists()                        # released at the end
    info = json.loads((out / manifest.NAME).read_text())
    assert info["status"] == "done" and info["cenfuse"] and info["inputs"]["g"]["chromosomes"] == 2


def test_a_failed_run_records_its_error_and_releases_the_lock(tmp_path, monkeypatch):
    _fasta(tmp_path / "in")
    out = tmp_path / "out"

    def boom(self, t0):
        raise CenfuseError("bad input")
    monkeypatch.setattr(Run, "_run", boom)
    with pytest.raises(CenfuseError):
        Run(tmp_path / "in", out).run()
    info = json.loads((out / manifest.NAME).read_text())
    assert info["status"] == "failed" and "bad input" in info["error"] and not (out / pipeline.LOCK).exists()


def test_manifest_names_tools_and_packages(tmp_path):
    info = dict(tools={t: manifest.tool_version(t, a) for t, a in manifest.TOOLS.items()}, packages=manifest.packages())
    assert info["packages"]["numpy"] and set(info["tools"]) == set(manifest.TOOLS)
    for v in info["tools"].values():
        assert set(v) == {"path", "version", "conda"}


# ------------------------------------------------------------------ the satellite engine's input keys
def _sat_pipeline(cdir, tmp_path, force=None):
    from cenfuse.sat.params import Params
    from cenfuse.sat.pipeline import Pipeline
    return Pipeline(cdir, 1, None, Params(), force=force, out_dir=tmp_path / "sat")


def _stamp_all(pipe):
    from cenfuse.sat.pipeline import STAGES
    for st in STAGES:
        (pipe.work / f"{st}.done").write_text(json.dumps(pipe.fingerprint(st)))


def test_adding_tes_recomputes_classify_onward_and_keeps_the_search(tmp_path, caplog):
    from cenfuse.sat.pipeline import STAGES
    cdir = tmp_path / "cx"
    _fasta(cdir)
    _stamp_all(_sat_pipeline(cdir, tmp_path))
    (cdir / "tes").mkdir()
    (cdir / "tes" / "g.gff3").write_text("##gff-version 3\nchr1\tx\tLTR_retrotransposon\t1\t100\t.\t+\t.\tID=a;clade=CRM\n")
    with caplog.at_level("WARNING", logger="cenfuse-sat"):
        pipe = _sat_pipeline(cdir, tmp_path)
    assert "recomputing from classify" in caplog.text and "inputs g tes" in caplog.text
    left = {st for st in STAGES if (pipe.work / f"{st}.done").exists()}
    assert left == {"seeds", "families", "annotate"}


def test_a_changed_search_parameter_still_stops(tmp_path):
    from dataclasses import replace
    from cenfuse.sat.params import Params
    from cenfuse.sat.pipeline import Pipeline
    from cenfuse.sat.util import SatelliteError
    cdir = tmp_path / "cx"
    _fasta(cdir)
    _stamp_all(_sat_pipeline(cdir, tmp_path))
    p = replace(Params(), min_fam_span=Params().min_fam_span + 1)
    with pytest.raises(SatelliteError, match="families: parameters min_fam_span"):
        Pipeline(cdir, 1, None, p, out_dir=tmp_path / "sat")


def test_unchanged_inputs_resume(tmp_path):
    from cenfuse.sat.pipeline import STAGES
    cdir = tmp_path / "cx"
    _fasta(cdir)
    _stamp_all(_sat_pipeline(cdir, tmp_path))
    time.sleep(0.01)
    pipe = _sat_pipeline(cdir, tmp_path)
    assert all((pipe.work / f"{st}.done").exists() for st in STAGES)


def test_a_reused_pid_or_a_garbled_lock_is_stale(tmp_path):
    """A live pid whose process started after the lock is another process (a reused pid); a lock that cannot be read is
    stale once it is STALE_AFTER old (a run killed between creating and writing it)."""
    lock = tmp_path / "x.lock"
    lock.write_text(f"{socket.gethostname()} {os.getpid()} ticks=1 run")                        # this pid, another start time
    assert "reused pid" in util.lock_stale(lock)
    lock.write_text(f"{socket.gethostname()} {os.getpid()} ticks={util._start_ticks(os.getpid())} run")
    assert util.lock_stale(lock) == ""                                                       # this very process: live
    lock.write_text(f"{socket.gethostname()} {os.getpid()} 2001-01-01T00:00:00 run")           # no ticks (older format): live
    assert util.lock_stale(lock) == ""
    lock.write_text("")
    assert util.lock_stale(lock) == ""                                                       # fresh: being written
    old = time.time() - util.STALE_AFTER - 60
    os.utime(lock, (old, old))
    assert "unreadable" in util.lock_stale(lock)
    assert util.try_lock(lock) and lock.read_text().split()[1] == str(os.getpid())


def test_a_stale_lock_is_taken_over_by_one_process_only(tmp_path, monkeypatch):
    """Two processes that both find one stale lock: the takeover guard lets one replace it; the other backs off and never
    removes the new owner's lock (code review, 0.14)."""
    lock = tmp_path / "x.lock"
    lock.write_text(f"{socket.gethostname()} 999999999 2001-01-01T00:00:00 run")              # a dead pid: stale
    guard = lock.with_name(lock.name + ".takeover")
    guard.write_text("")                                                                     # another process is taking over
    assert not util.try_lock(lock) and "999999999" in lock.read_text()
    guard.unlink()
    assert util.try_lock(lock) and not guard.exists()
    mine = lock.read_text()
    assert not util.try_lock(lock) and lock.read_text() == mine                             # live now: kept


def test_a_module_edit_recomputes_exactly_the_stages_whose_code_holds_it(tmp_path, monkeypatch):
    """The code half of future.md #59's check: an edit to what a module does changes the key of every stage whose import
    closure holds it, and of no other; a comment or docstring edit changes none."""
    pkg = tmp_path / "cenfuse"
    shutil.copytree(pipeline._PKG, pkg, ignore=shutil.ignore_patterns("__pycache__", "data"))
    monkeypatch.setattr(pipeline, "_PKG", pkg)

    def keys():
        for f in (pipeline.code_key, pipeline._stage_method):
            f.cache_clear()
        return dict(telomere=code_key("telomere.py", stage="telomere"), rdna=code_key("rdna.py", stage="rdna"),
                    hic=code_key("hicfuse.py", packages=("hic",), stage="hic"))
    k0 = keys()
    (pkg / "telomere.py").write_text((pkg / "telomere.py").read_text() + "\n# a comment\n")
    assert keys() == k0
    (pkg / "telomere.py").write_text((pkg / "telomere.py").read_text() + "\n_X_ = 1\n")
    k1 = keys()
    assert k1["telomere"] != k0["telomere"] and k1["rdna"] == k0["rdna"] and k1["hic"] == k0["hic"]
    (pkg / "util.py").write_text((pkg / "util.py").read_text() + "\n_Y_ = 2\n")       # in every stage's key
    k2 = keys()
    assert all(k2[s] != k1[s] for s in k1)
    p = pkg / "pipeline.py"                                                          # one stage's own method
    p.write_text(p.read_text().replace("hit = telomere.detect(g.fasta, g.chroms)", "hit = telomere.detect(g.fasta, dict(g.chroms))"))
    k3 = keys()
    assert k3["telomere"] != k2["telomere"] and k3["rdna"] == k2["rdna"]
    for f in pipeline.KEYLESS:                                    # infrastructure: imported by every engine, in no key
        (pkg / f).write_text((pkg / f).read_text() + "\n_Z_ = 3\n")
    assert keys() == k3


def test_a_same_length_fasta_edit_changes_the_assembly_key(tmp_path):
    """Hard-masking keeps the .fai byte for byte; the stages that read the sequence must still recompute (code review, 0.14)."""
    from types import SimpleNamespace
    f = _fasta(tmp_path / "a")
    g = SimpleNamespace(fasta=f)
    k0, fai0 = Run._assembly(g), (tmp_path / "a" / "g.fa.fai").read_text()
    s = f.read_text()
    f.write_text(s[:100] + "N" * 1000 + s[1100:])
    pysam.faidx(str(f))
    util._MD5.clear()                         # a resume is a new process: the per-process md5 cache starts empty
    assert (tmp_path / "a" / "g.fa.fai").read_text() == fai0 and Run._assembly(g) != k0


def test_a_live_forked_owner_is_never_stale(tmp_path):
    """A lock made by a forked process that is still running is that process's, whatever its clock or /proc entry say."""
    lock = tmp_path / "x.lock"
    r, w = os.pipe()
    pid = os.fork()
    if pid == 0:                                              # the child takes the lock, says so, and waits
        ok = util.try_lock(lock)
        os.write(w, b"1" if ok else b"0")
        time.sleep(3)
        os._exit(0)
    assert os.read(r, 1) == b"1"
    try:
        assert util.lock_stale(lock) == "" and not util.try_lock(lock)
    finally:
        os.waitpid(pid, 0)
    assert "is gone" in util.lock_stale(lock)                 # once it ended, the lock is stale


def test_the_manifest_never_reads_its_own_run(tmp_path):
    (tmp_path / "satellite").mkdir()
    manifest.write(tmp_path, dict(cenfuse="9", host=socket.gethostname(), pid=os.getpid(), status="running"))
    assert manifest.read(tmp_path / "satellite") is None
    manifest.write(tmp_path, dict(cenfuse="9", host=socket.gethostname(), pid=1, status="done"))
    assert manifest.read(tmp_path / "satellite")["cenfuse"] == "9"


def test_a_stamp_of_before_0_14_is_trusted_only_for_an_older_fasta(tmp_path):
    """Pre-0.14 satellite stamps hold no FASTA md5. Accepted (and upgraded to the full key) while the FASTA is older than them;
    when the FASTA was written after them, its content cannot be vouched for and the search stops for --force (code review)."""
    from cenfuse.sat.pipeline import STAGES
    from cenfuse.sat.util import SatelliteError
    cdir = tmp_path / "cx"
    f = _fasta(cdir)
    pipe = _sat_pipeline(cdir, tmp_path)
    for st in STAGES:                                          # a 0.13 run's stamps: no 'fasta' key
        fp = pipe.fingerprint(st)
        fp.pop("fasta")
        (pipe.work / f"{st}.done").write_text(json.dumps(fp))
    old = time.time() - 3600
    os.utime(f, (old, old))
    pipe = _sat_pipeline(cdir, tmp_path)                       # FASTA older: accepted, upgraded
    assert all("fasta" in json.loads((pipe.work / f"{st}.done").read_text()) for st in STAGES)
    for st in STAGES:                                          # back to 0.13 stamps, then the FASTA rewritten
        fp = json.loads((pipe.work / f"{st}.done").read_text())
        fp.pop("fasta")
        (pipe.work / f"{st}.done").write_text(json.dumps(fp))
    time.sleep(0.05)
    f.write_text(f.read_text())
    with pytest.raises(SatelliteError, match="FASTA is newer than this checkpoint"):
        _sat_pipeline(cdir, tmp_path)


def test_threads_must_be_a_positive_whole_number(capsys):
    from cenfuse import cli
    for bad in ("0", "abc"):
        with pytest.raises(SystemExit):
            cli.main(["run", "/nonexistent", "-t", bad])
        err = capsys.readouterr().err
        assert "must be" in err and "_threads" not in err


def test_atomic_writes_leave_a_hard_linked_copy_alone(tmp_path):
    """future.md #59 (0.14): a write goes to a temporary file beside its target and is moved over it, so the other name of a
    hard-linked file (a folder copied with cp -al) keeps its content; a failed write leaves the target as it was; the
    temporary name keeps the suffix, so pandas still compresses a .gz."""
    import pandas as pd
    a = tmp_path / "a.tsv"
    a.write_text("old\n")
    b = tmp_path / "b.tsv"
    os.link(a, b)
    util.atomic_write_text(b, "new\n")
    assert a.read_text() == "old\n" and b.read_text() == "new\n" and os.stat(a).st_nlink == 1
    with pytest.raises(RuntimeError):
        with util.atomic_open(b) as f:
            f.write("half")
            raise RuntimeError
    assert b.read_text() == "new\n" and sorted(p.name for p in tmp_path.iterdir()) == ["a.tsv", "b.tsv"]
    util.atomic_csv(pd.DataFrame({"x": [1, 2]}), tmp_path / "t.tsv.gz", sep="\t", index=False)
    assert (tmp_path / "t.tsv.gz").read_bytes()[:2] == b"\x1f\x8b"
    assert pd.read_csv(tmp_path / "t.tsv.gz", sep="\t").x.tolist() == [1, 2]
    src = tmp_path / "src.txt"
    src.write_text("copied\n")
    util.atomic_copy(src, a)
    assert a.read_text() == "copied\n"


def test_unshare_tree_unlinks_only_shared_files(tmp_path):
    d = tmp_path / "stage"
    (d / "sub").mkdir(parents=True)
    (d / "own.txt").write_text("x")
    (d / "sub" / "shared.txt").write_text("y")
    os.link(d / "sub" / "shared.txt", tmp_path / "other.txt")
    (d / "link").symlink_to(tmp_path / "other.txt")
    assert util.unshare_tree(d) == 1
    assert (d / "own.txt").exists() and not (d / "sub" / "shared.txt").exists() and (d / "link").is_symlink()
    assert (tmp_path / "other.txt").read_text() == "y"


@pytest.mark.skipif(not os.path.exists("/proc/self/status"), reason="needs /proc (Linux)")
def test_usage_counts_the_memory_of_child_processes():
    """run_info.json's cost (future.md #5, #61): the peak memory of the whole process tree, which GNU time's maximum (the
    largest single process) misses, and the CPU time of the finished children."""
    import subprocess
    import sys
    with util.Usage(every=0.1) as use:
        subprocess.run([sys.executable, "-c", "x = bytearray(400 * 2 ** 20)\nimport time\nt = time.process_time()\n"
                        "while time.process_time() - t < 2.0: pass"], check=True)
    c = use.cost()
    assert c["peak_memory_gb"] >= 0.35 and c["wall_minutes"] >= 0.03
    assert c["cpu_hours"] >= 1.5 / 3600                                       # the child's CPU, once it has ended
    assert util.tree_memory_kb() > 0
    if os.path.exists("/proc/self/smaps_rollup"):
        assert c["memory"].startswith("the process tree's proportional set size")


def _log_in_worker(x):
    import logging
    logging.getLogger("cenfuse-sat").warning(f"from worker {x}")
    return x


def _die_holding_the_log_queue(x):
    import logging
    import signal
    h = logging.getLogger("cenfuse-sat").handlers[0]                       # the QueueHandler process_pool gave this worker
    h.queue._wlock.acquire()                                                # as if killed while sending a record
    os.kill(os.getpid(), signal.SIGKILL)


def test_a_worker_killed_while_sending_a_log_record_does_not_hang_the_run(monkeypatch):
    """Code review of 0.14: a worker killed while it held the log queue's lock left the listener blocked for good, so the run
    hung (and its heartbeat kept the folder's lock fresh). Now the run fails as the pool fails, within LISTENER_STOP seconds."""
    import multiprocessing
    from concurrent.futures.process import BrokenProcessPool
    monkeypatch.setattr(util, "LISTENER_STOP", 1)
    t0 = time.time()
    with pytest.raises(BrokenProcessPool):
        with util.process_pool(1, multiprocessing.get_context("spawn")) as ex:
            list(ex.map(_die_holding_the_log_queue, [1]))
    assert time.time() - t0 < 30


@pytest.mark.parametrize("method", ["fork", "spawn", "forkserver"])
def test_worker_log_records_reach_the_parent_under_every_start_method(method, tmp_path):
    """future.md #60: a worker started by spawn (macOS) or forkserver (Python 3.14's Linux default) has no handlers of its own;
    runtime.process_pool sends its records to the parent, which writes them. Under fork the worker inherits the handlers and
    writes its lines itself, as before. Either way they reach the run's log (a file here, stderr in a run)."""
    import logging
    import multiprocessing
    if method not in multiprocessing.get_all_start_methods():
        pytest.skip(f"no {method} here")
    lg = logging.getLogger("cenfuse-sat")
    keep = lg.handlers[:], lg.propagate, lg.level
    h = logging.FileHandler(tmp_path / "run.log")
    lg.handlers[:], lg.propagate = [h], False                 # as util.setup_logging configures the engine loggers
    lg.setLevel(logging.INFO)
    try:
        with util.process_pool(2, multiprocessing.get_context(method)) as ex:
            assert list(ex.map(_log_in_worker, [1, 2])) == [1, 2]
    finally:
        h.close()
        lg.handlers[:], lg.propagate = keep[0], keep[1]
        lg.setLevel(keep[2])
    assert sorted((tmp_path / "run.log").read_text().split("\n")[:-1]) == ["from worker 1", "from worker 2"]


# ------------------------------------------------------------------ stage gates
@pytest.mark.parametrize("n, runs", [(pipeline.atline.MAX_CHROMS, True), (pipeline.atline.MAX_CHROMS + 1, False)])
def test_fasta_lines_do_not_read_a_genome_their_test_cannot_judge(tmp_path, monkeypatch, n, runs):
    """The AT-trough, repeat-density and repeat-family stages read the whole genome; their one-per-chromosome test refuses
    more than atline.MAX_CHROMS sequences (contigs, not chromosomes), so beyond that the stages do not run, and the
    report's notes say why, once (0.16.0 read Ae. geniculata-sized genomes three times for a certain `untested`)."""
    import pandas as pd
    from types import SimpleNamespace
    d = tmp_path / "in"
    d.mkdir()
    names = [f"s{i}" for i in range(n)]
    (d / "g.fa").write_text("".join(f">{c}\nA\n" for c in names))
    (d / "g.fa.fai").write_text("".join(f"{c}\t300000\t{i}\t300000\t300001\n" for i, c in enumerate(names)))
    (d / "g.chromosomes.txt").write_text("\n".join(names) + "\n")
    run = Run(d, tmp_path / "out")
    g = run.genomes[0]
    assert len(g.chroms) == n
    seq = SimpleNamespace(calls=pd.DataFrame(dict(chrom=names, start=-1, evidence="none", confidence="none", alt_start=-1)),
                          arrays=pd.DataFrame(columns=["chrom", "start", "end"]))

    class Read(Exception):
        pass

    def read(*a, **k):
        raise Read
    for mod, fn in ((pipeline.atline, "gc_bins"), (pipeline.repline, "repeat_bins"), (pipeline.famline, "sample")):
        monkeypatch.setattr(mod, fn, read)
    for stage in (run.stage_attrough, run.stage_repeats, run.stage_repfam):
        if runs:
            with pytest.raises(Read):
                stage(g, seq)
        else:
            assert stage(g, seq) is None
    notes = [x for x in g.notes if "not run" in x]
    assert notes == ([] if runs else [f"the AT-trough, repeat-density and repeat-family lines are not run: {n} sequences, "
                                      f"more than their one-per-chromosome test takes ({pipeline.atline.MAX_CHROMS})"])


@pytest.mark.parametrize("name", ["NA", "null", "None"])
def test_the_methylation_key_covers_every_file_the_engine_reads(tmp_path, name):
    """The methylation stage's key and run_info.json fingerprint the files of meth/<g>.tsv; they read the table as the engine
    does (meth.io.read_samples: names as written, only an empty cell is missing), so a file named NA, null or None is not
    left out of them while the engine reads it."""
    from types import SimpleNamespace
    (tmp_path / name).write_text("x")
    (tmp_path / "g_satellite_arrays.gff").write_text("")
    table = tmp_path / "g.tsv"
    table.write_text(f"sample\tCHG\tCG\np1\t{name}\t\n")
    g = SimpleNamespace(name="g", meth=table, tes=None, genes=None, exclude=None)
    assert Run._meth_inputs_fp(None, g, tmp_path)["files"] == [file_fingerprint(tmp_path / name)]
    assert manifest._meth_files(table) == [str(tmp_path / name)]


def test_control_tables_keep_names_that_pandas_would_read_as_missing(tmp_path):
    """controls._chunks reads fields as written (na_filter=False): a sequence named NA stays NA."""
    from cenfuse import controls
    f = tmp_path / "t.tsv"
    f.write_text("NA\t1\nnull\t2\n")
    rows = [r for ch in controls._chunks(f, ["chrom", "pos"], 0, ints=("pos",)) for r in ch.itertuples(index=False)]
    assert [(r.chrom, r.pos) for r in rows] == [("NA", 1), ("null", 2)]
