from pathlib import Path

import numpy as np
import pandas as pd
import pysam
import pytest

from cenfuse.sat.genomes import Genome, discover_genomes
from cenfuse.sat.params import Params
from cenfuse.sat.seeds import COLS, cache_path, chunks, dedupe, genome_key, seed_genomes
from cenfuse.sat.util import union_bp
from sat_common import have


def test_chunks_cover_with_overlap():
    assert list(chunks({"a": 25, "b": 7}, 10, 3)) == [
        ("a", 0, 10), ("a", 7, 17), ("a", 14, 24), ("a", 21, 25), ("b", 0, 7)]


def test_dedupe_keeps_best_of_overlapping_duplicates():
    r1 = ("c", 100, 1000, 150, 6.0, 150, 90, 1, 500, 1.9, "A")
    r2 = ("c", 120, 1000, 151, 5.8, 151, 90, 1, 600, 1.9, "B")   # same array seen from the next chunk
    r3 = ("c", 100, 1000, 300, 3.0, 300, 90, 1, 400, 1.9, "C")   # different period: kept
    assert sorted(r[-1] for r in dedupe([r1, r2, r3])) == ["B", "C"]


def test_dedupe_period_tolerance_scales_with_the_unit():
    """A FasTAN unit is a mean over the copies aligned, so the same 1-kb array seen from two chunks can differ by a
    few bp; at tol 0 (TRF) only 1 bp is allowed."""
    r1 = ("c", 0, 10_000, 1001, 10.0, 1001, 90, 0, 9000, 0, "A")
    r2 = ("c", 0, 9_000, 1004, 9.0, 1004, 90, 0, 8100, 0, "B")
    assert [r[-1] for r in dedupe([r1, r2], 0.05)] == ["A"]
    assert sorted(r[-1] for r in dedupe([r1, r2])) == ["A", "B"]


def test_trf_cache_name_is_unchanged(tmp_path):
    """TRF caches written by earlier satellite-engine versions (hours per large genome) must still be found."""
    assert cache_path(tmp_path, "abc", Params(seed_engine="trf")).name == "abc.25429390.trf.tsv.gz"


def test_cache_key_depends_on_engine_and_chunking(tmp_path):
    a = cache_path(tmp_path, "abc", Params())
    assert a == cache_path(tmp_path, "abc", Params()) and a.name.endswith(".fastan.tsv.gz")
    assert a != cache_path(tmp_path, "abc", Params(seed_chunk=1_000_000))
    assert a != cache_path(tmp_path, "abc", Params(seed_engine="trf"))


# ---- genome_key: the FASTA md5, plus a suffix for the chromosome set where it can differ from the file's sequences.
# Literals: THREE_MD5 = md5 of THREE_SEQS; SET_* = md5("\n".join(sorted(names)))[:8] for the sets a, a+b and a+b+c. The
# keys of the indexed layouts are what genome_key returned before the link-indexed fix (computed with that code).

THREE_SEQS = ">a\nACGT\n>b\nACGT\n>c\nACGT\n"
THREE_MD5 = "9b8d2a54af5e98d5b9a4509765fbd8d6"
SET_A, SET_AB, SET_ABC = "0cc175b9", "8cdeb444", "c32b2057"


def _view_fasta(tmp_path, layout: str):
    """The FASTA of a cenfuse view of one genome: inputs.make_view links it to the user's FASTA, resolved. The user's
    FASTA holds a, b and c and has its .fai beside the path cenfuse was given. layout: 'file' (a file, indexed beside
    it), 'link_to_indexed' (a link to an indexed file; the index is linked beside it too) or 'link_indexed' (a link to
    a file without .fai: samtools faidx ran on the link, so the .fai sits beside the link)."""
    real = tmp_path / "shared" / "ref.fa"
    real.parent.mkdir()
    real.write_text(THREE_SEQS)
    user = real if layout == "file" else tmp_path / "genomes" / "S.fa"
    if user != real:
        user.parent.mkdir()
        user.symlink_to(real)
    pysam.faidx(str(user if layout == "link_indexed" else real))
    if layout == "link_to_indexed":
        Path(f"{user}.fai").symlink_to(f"{real}.fai")
    assert Path(f"{user}.fai").exists()                              # cenfuse needs the .fai beside the path it gets
    assert Path(f"{real}.fai").exists() == (layout != "link_indexed")    # and the layout is what its name says
    (tmp_path / "view").mkdir()
    (tmp_path / "view" / "S.fa").symlink_to(user.resolve())
    return tmp_path / "view" / "S.fa"


def _key(fasta, names: str) -> str:
    return genome_key(Genome("S", fasta, {c: 4 for c in names}))


def test_link_indexed_fasta_keys_separate_chromosome_sets(tmp_path):
    """The user's FASTA is a link and its .fai sits beside the link, so the real file has none to compare the
    chromosome set with. The set must still decide the key: otherwise `--force satellite` after a chromosome-set change
    finds the old set's seeds in the cache and reuses them."""
    fa = _view_fasta(tmp_path, "link_indexed")
    keys = [_key(fa, names) for names in ("abc", "a", "ab")]
    assert len(set(keys)) == 3
    assert keys[1:] == [f"{THREE_MD5}-{SET_A}", f"{THREE_MD5}-{SET_AB}"]      # the keys the indexed layouts give
    assert _key(fa, "ba") == keys[2]                                          # the set decides, not its order


def test_link_indexed_fasta_with_every_sequence_gets_the_set_suffix(tmp_path):
    """Without the real file's .fai nothing says the chromosomes are all of the file's sequences, so the whole file
    gets the suffix too (and never shares an entry with a subset made by a cenfuse that left the suffix out)."""
    key = _key(_view_fasta(tmp_path, "link_indexed"), "abc")
    assert key != THREE_MD5 and key == f"{THREE_MD5}-{SET_ABC}"


@pytest.mark.parametrize("layout", ["file", "link_to_indexed"])
def test_keys_are_unchanged_when_the_real_file_is_indexed(tmp_path, layout):
    """Seed caches made before the fix stay valid (hours per large genome): the whole file keeps the plain md5, a
    subset gets the suffix."""
    fa = _view_fasta(tmp_path, layout)
    assert _key(fa, "abc") == THREE_MD5
    assert _key(fa, "a") == f"{THREE_MD5}-{SET_A}"
    assert _key(fa, "ab") == f"{THREE_MD5}-{SET_AB}"


@pytest.mark.skipif(not have("FasTAN", "FAtoGDB", "ONEview"), reason="FasTAN tools missing")
def test_fastan_finds_planted_satellite_across_chunk_boundary(synth_seeds, synth_complex):
    cdir, truth = synth_complex
    df = synth_seeds["Pa_A"]
    s, e = truth["genomes"]["Pa_A"]["Pa_A01"]["cen"]      # 130,000-169,000 crosses the 145-150 kb chunk overlap
    k = (df.period / 150).round()
    mult = (k >= 1) & ((df.period / 150 - k).abs() <= 0.03 * k)
    hit = df[(df.chrom == "Pa_A01") & mult & (df.start < e) & (df.end > s)]
    assert union_bp(zip(hit.start, hit.end)) > 0.9 * (e - s)
    assert (hit.consensus.str.len() - hit.period).abs().max() <= 3
    for p in hit.itertuples():                            # no chunk duplicates left
        dup = hit[(hit.index != p.Index) & ((hit.period - p.period).abs() <= 0.05 * p.period)
                  & (hit.start < p.end) & (hit.end > p.start)]
        for q in dup.itertuples():
            ov = min(p.end, q.end) - max(p.start, q.start)
            assert ov < 0.5 * min(p.end - p.start, q.end - q.start)


@pytest.mark.skipif(not have("FasTAN", "FAtoGDB", "ONEview"), reason="FasTAN tools missing")
def test_seed_cache_is_reused(synth_complex, synth_seed_cache, test_params):
    cdir, _ = synth_complex
    genomes = discover_genomes(cdir)
    before = {p: p.stat().st_mtime for p in synth_seed_cache.glob("*.fastan.tsv.gz")}
    seed_genomes(genomes, synth_seed_cache, 2, test_params)
    after = {p: p.stat().st_mtime for p in synth_seed_cache.glob("*.fastan.tsv.gz")}
    assert before == after and len(after) == 4


def _one_genome(tmp_path):
    from cenfuse.sat.genomes import Genome
    fa = tmp_path / "g.fa"
    fa.write_text(">c\n" + "ACGTTGCA" * 500 + "\n")
    pysam.faidx(str(fa))
    return Genome("g", fa, {"c": 4000})


def test_seeds_wait_for_lock_held_by_live_process(tmp_path, test_params, monkeypatch):
    """A second process must not recompute (or corrupt) a genome another live process is running."""
    import gzip
    import os
    import threading
    import time
    import cenfuse.sat.seeds as seeds
    from cenfuse.sat.util import md5_file
    monkeypatch.setattr(seeds, "POLL", 0.1)
    g = _one_genome(tmp_path)
    cache = tmp_path / "cache"
    cache.mkdir()
    path = seeds.cache_path(cache, md5_file(g.fasta), test_params)
    lock = seeds.lock_path(path)
    lock.write_text(str(os.getpid()))                       # held by a live process (this one)
    res = {}
    t = threading.Thread(target=lambda: res.update(seeds.seed_genomes([g], cache, 1, test_params)))
    t.start()
    time.sleep(0.5)
    assert t.is_alive()                                      # still waiting on the lock
    with gzip.open(path, "wt") as f:                          # the other process finishes
        f.write("\t".join(seeds.COLS) + "\nFROM_OTHER_PROCESS\n")
    lock.unlink()
    t.join(5)
    assert not t.is_alive() and "FROM_OTHER_PROCESS" in gzip.open(res["g"], "rt").read()


@pytest.mark.skipif(not have("FasTAN", "FAtoGDB", "ONEview"), reason="FasTAN tools missing")
def test_seeds_steal_stale_lock(tmp_path, test_params):
    import socket
    import cenfuse.sat.seeds as seeds
    from cenfuse.sat.util import md5_file
    g = _one_genome(tmp_path)
    cache = tmp_path / "cache"
    cache.mkdir()
    path = seeds.cache_path(cache, md5_file(g.fasta), test_params)
    seeds.lock_path(path).write_text(f"{socket.gethostname()} 999999999")   # dead pid on this host
    out = seeds.seed_genomes([g], cache, 1, test_params)
    assert out["g"].exists() and not seeds.lock_path(path).exists()


def test_seeds_never_steal_a_lock_held_on_another_host(tmp_path, test_params, monkeypatch):
    """On a shared filesystem (SLURM nodes) a PID from another host cannot be checked here: its lock is never stolen."""
    import threading
    import time
    import cenfuse.sat.seeds as seeds
    from cenfuse.sat.util import md5_file
    monkeypatch.setattr(seeds, "POLL", 0.1)
    g = _one_genome(tmp_path)
    cache = tmp_path / "cache"
    cache.mkdir()
    path = seeds.cache_path(cache, md5_file(g.fasta), test_params)
    lock = seeds.lock_path(path)
    lock.write_text("some-other-node 999999999")             # dead here, but on another host
    t = threading.Thread(target=lambda: seeds.seed_genomes([g], cache, 1, test_params), daemon=True)
    t.start()
    time.sleep(0.5)
    assert t.is_alive() and lock.read_text() == "some-other-node 999999999"


def test_seed_lock_records_host_and_pid(tmp_path):
    import os
    import socket
    import cenfuse.sat.seeds as seeds
    lock = tmp_path / "x.lock"
    assert seeds._acquire(lock) and lock.read_text() == f"{socket.gethostname()} {os.getpid()}"


@pytest.mark.skipif(not have("FasTAN", "FAtoGDB", "ONEview"), reason="FasTAN tools missing")
def test_seeds_cache_each_genome_as_soon_as_it_is_done(tmp_path, test_params):
    """A failure in a later genome must not lose the finished genomes of the same run (hours on real data)."""
    from cenfuse.sat.genomes import Genome
    from cenfuse.sat.seeds import _compute
    rng = np.random.default_rng(3)
    fa = tmp_path / "a.fa"
    fa.write_text(">c1\n" + "".join(rng.choice(list("ACGT"), 30_000)) + "\n")
    pysam.faidx(str(fa))
    good = Genome("a", fa, {"c1": 30_000})
    bad = Genome("b", tmp_path / "missing.fa", {"c1": 30_000})
    pa, pb = tmp_path / "a.fastan.tsv.gz", tmp_path / "b.fastan.tsv.gz"
    with pytest.raises(Exception):
        _compute({pa: good, pb: bad}, tmp_path, 2, test_params)
    assert pa.exists() and not pb.exists()
    assert list(pd.read_csv(pa, sep="\t").columns) == COLS


def test_genome_without_chunks_is_cached(tmp_path, test_params):
    """Zero-length sequences give no chunks; the cache must still be written, or seed_genomes loops. (The trf
    engine needs no tool when there is nothing to scan.)"""
    from dataclasses import replace
    from cenfuse.sat.genomes import Genome
    from cenfuse.sat.seeds import _compute
    path = tmp_path / "z.trf.tsv.gz"
    _compute({path: Genome("z", tmp_path / "z.fa", {"c1": 0})}, tmp_path, 1, replace(test_params, seed_engine="trf"))
    assert path.exists() and pd.read_csv(path, sep="\t").empty


# ---- FasTAN crashes: recovered by splitting the piece (fastan.chunk_rows) and recorded beside the cached seed table

@pytest.mark.parametrize("overlap, chr1_starts, chr2_starts", [
    (10_000, [100, 95_100, 190_100], [100, 190_100]),
    (20_000, [100, 90_100, 180_100], [100, 180_100])])      # the overlap is seed_overlap: the halves' and the chunks'
def test_crash_table_beside_the_cache(tmp_path, monkeypatch, caplog, overlap, chr1_starts, chr2_starts):
    """FasTAN aborts on the first 200-kb chunk of chr1 (here: on any chr1 piece over 150 kb). The run still finishes,
    with the seeds of that chunk's two halves, and the crashed interval is listed beside the cached seed table, where
    every run that reuses the cache finds it."""
    import logging
    import cenfuse.sat.seeds as seeds
    from cenfuse.sat import fastan
    from cenfuse.sat.genomes import Genome
    from cenfuse.sat.util import log
    from sat_common import crashing_fastan, dna
    monkeypatch.setattr(log, "propagate", True)
    rng = np.random.default_rng(4)
    chroms = {c: dna(rng, 300_000) for c in ("chr1", "chr2")}
    fa = tmp_path / "g.fa"
    fa.write_text("".join(f">{c}\n{s}\n" for c, s in chroms.items()))
    pysam.faidx(str(fa))
    g = Genome("g", fa, {c: len(s) for c, s in chroms.items()})
    monkeypatch.setattr(fastan, "setup", lambda tmpdir: False)         # the real setup runs FasTAN on a test sequence
    crashing_fastan(monkeypatch, chroms, lambda name, a, b: name == "chr1" and b - a > 150_000)   # forked workers keep it
    with caplog.at_level(logging.INFO, logger="cenfuse-sat"):
        cache = seeds.seed_genomes([g], tmp_path / "cache", 2, Params(seed_chunk=200_000, seed_overlap=overlap))["g"]
    assert cache.exists()
    lines = [x.split("\t") for x in seeds.crash_path(cache).read_text().splitlines()]
    assert lines == [["chrom", "start", "end", "signal", "action"], ["chr1", "0", "200000", "SIGABRT", "split"]]
    assert [r for r in caplog.records if r.levelno == logging.INFO
            and r.getMessage().startswith("seeds: g: FasTAN crashed on 1 chunk(s)")]
    df = seeds.load_seeds(cache)      # the stand-in finds one array at [100, 300) of each piece it finishes
    assert df[df.chrom == "chr1"].start.tolist() == chr1_starts      # the two halves of chunk 1, then chunk 2
    assert df[df.chrom == "chr2"].start.tolist() == chr2_starts


def test_a_chunk_fastan_searched_no_piece_of_stops_the_run(tmp_path, monkeypatch):
    """FasTAN aborts on every piece of chr1's first 200-kb chunk, down to the skip floor. The run stops with an error
    that names the chunk, the crash table and --seed-engine trf. The table is written (the chunk's crashes, for a
    report) and no seed table is cached, so a rerun recomputes the genome."""
    import cenfuse.sat.seeds as seeds
    from cenfuse.sat import fastan
    from cenfuse.sat.genomes import Genome
    from cenfuse.sat.util import SatelliteError
    from sat_common import crashing_fastan, dna
    rng = np.random.default_rng(9)
    chroms = {c: dna(rng, 300_000) for c in ("chr1", "chr2")}
    fa = tmp_path / "g.fa"
    fa.write_text("".join(f">{c}\n{s}\n" for c, s in chroms.items()))
    pysam.faidx(str(fa))
    g = Genome("g", fa, {c: len(s) for c, s in chroms.items()})
    monkeypatch.setattr(fastan, "setup", lambda tmpdir: False)
    crashing_fastan(monkeypatch, chroms, lambda name, a, b: name == "chr1" and b <= 200_000)   # chunk 1 of chr1 only
    p = Params(seed_chunk=200_000, seed_overlap=10_000)
    cache = seeds.cache_path(tmp_path / "cache", seeds.genome_key(g), p)
    with pytest.raises(SatelliteError) as err:
        seeds.seed_genomes([g], tmp_path / "cache", 2, p)
    msg, table = str(err.value), seeds.crash_path(cache)
    assert "every piece of chr1:0-200,000" in msg and str(table) in msg and "--seed-engine trf" in msg
    assert not cache.exists() and not seeds.lock_path(cache).exists()
    rows = [x.split("\t") for x in table.read_text().splitlines()]
    assert rows[0] == seeds.CRASH_COLS and rows[1] == ["chr1", "0", "200000", "SIGABRT", "split"]
    assert {r[4] for r in rows[1:]} == {"split", "skipped"}
    assert all(r[0] == "chr1" and int(r[2]) <= 200_000 for r in rows[1:])


def test_write_cache_lists_the_crashes_and_summarizes_them(tmp_path, monkeypatch, caplog):
    """events: one list per chunk FasTAN crashed on. The table has a row per event; the log line counts the chunks, the
    pieces split and the bp of skipped pieces (skipped pieces that overlap are counted once)."""
    import logging
    import cenfuse.sat.seeds as seeds
    from cenfuse.sat.genomes import Genome
    from cenfuse.sat.util import log
    monkeypatch.setattr(log, "propagate", True)
    path = tmp_path / "x.fastan.tsv.gz"
    chunk1 = [("c1", 0, 200_000, "SIGABRT", "split"), ("c1", 100_000, 115_000, "SIGABRT", "skipped"),
              ("c1", 110_000, 125_000, "SIGABRT", "skipped")]
    chunk2 = [("c2", 0, 15_000, "SIGSEGV", "skipped")]
    with caplog.at_level(logging.INFO, logger="cenfuse-sat"):
        seeds._write_cache(path, Genome("g", tmp_path / "g.fa", {}), [], 0.05, [chunk1, chunk2])
    assert seeds.crash_path(path).read_text() == (
        "chrom\tstart\tend\tsignal\taction\n"
        "c1\t0\t200000\tSIGABRT\tsplit\nc1\t100000\t115000\tSIGABRT\tskipped\nc1\t110000\t125000\tSIGABRT\tskipped\n"
        "c2\t0\t15000\tSIGSEGV\tskipped\n")
    line = [r.getMessage() for r in caplog.records if "FasTAN crashed" in r.getMessage()]
    assert len(line) == 1 and line[0].startswith("seeds: g: FasTAN crashed on 2 chunk(s)")
    assert "1 piece(s) split" in line[0] and "40,000 bp skipped" in line[0]       # 25,000 (union on c1) + 15,000 on c2


def test_the_crash_table_is_written_before_the_cache_file(tmp_path, monkeypatch):
    """Other runs wait for the cache file and then look for the table beside it: if the cache file cannot be written,
    the table is already there, so the file never shows up without it."""
    import gzip
    import cenfuse.sat.seeds as seeds
    from cenfuse.sat.genomes import Genome

    def refuse(*a, **k):
        raise OSError("disk full")

    monkeypatch.setattr(gzip, "open", refuse)
    path = tmp_path / "x.fastan.tsv.gz"
    with pytest.raises(OSError):
        seeds._write_cache(path, Genome("g", tmp_path / "g.fa", {}), [], 0.05,
                           [[("c1", 0, 200_000, "SIGABRT", "split")]])
    assert seeds.crash_path(path).exists() and not path.exists()


def test_a_crash_table_left_by_an_earlier_run_is_dropped_when_nothing_crashed(tmp_path):
    """The table describes the cache file beside it: recomputing the cache without a crash must not keep the old one."""
    import cenfuse.sat.seeds as seeds
    from cenfuse.sat.genomes import Genome
    path = tmp_path / "x.fastan.tsv.gz"
    seeds.crash_path(path).write_text("chrom\tstart\tend\tsignal\taction\nc1\t0\t200000\tSIGABRT\tsplit\n")
    seeds._write_cache(path, Genome("g", tmp_path / "g.fa", {}), [], 0.05)
    assert path.exists() and not seeds.crash_path(path).exists()


def test_stage_seeds_links_the_crash_table_only_when_there_is_one(tmp_path, monkeypatch):
    """work/seeds/<genome>.fastan_crashes.tsv is a link to the table beside the cached seed table; a genome FasTAN never
    crashed on has none, and a link left by an earlier run does not outlive its table."""
    import cenfuse.sat.seeds as seeds
    from cenfuse.sat.pipeline import Pipeline
    for n in "ab":
        (tmp_path / f"{n}.fa").write_text(">c\nACGT\n")
        (tmp_path / f"{n}.fa.fai").write_text("c\t4\t3\t4\t5\n")
    cache = tmp_path / "cache"
    cache.mkdir()
    paths = {n: cache / f"{n}.fastan.tsv.gz" for n in "ab"}
    for p in paths.values():
        p.write_text("")
    seeds.crash_path(paths["a"]).write_text("chrom\tstart\tend\tsignal\taction\n")
    monkeypatch.setattr(seeds, "seed_genomes", lambda *a, **k: paths)
    pl = Pipeline(tmp_path, 1, cache, Params())
    (pl.work / "seeds").mkdir()
    (pl.work / "seeds" / "b.fastan_crashes.tsv").symlink_to(tmp_path / "gone.tsv")    # left by an earlier run
    pl.stage_seeds()
    assert (pl.work / "seeds" / "a.fastan_crashes.tsv").resolve() == seeds.crash_path(paths["a"]).resolve()
    assert not (pl.work / "seeds" / "b.fastan_crashes.tsv").is_symlink()
    assert (pl.work / "seeds" / "b.seeds.tsv.gz").resolve() == paths["b"].resolve()   # the seed table link is unchanged
