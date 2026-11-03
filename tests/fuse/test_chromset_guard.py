"""Resume guard for a changed chromosome set (pipeline.Run._check_chromosome_set).

The satellite engine's checkpoints (seeds, families, calls) are made on one set of chromosomes, but its resume check
covers only the code version and the parameters. cenfuse 0.9 changes the set on assemblies with a gradual scaffold tail
and takes a <genome>.chromosomes.txt, so a resumed folder would reuse seeds, families and calls made on other
sequences. The run stops instead and says to redo the satellite stages. The engine is replaced by a stand-in that
writes its seeds checkpoint and nothing else (or fails right after it)."""
import json
import logging
import os

import numpy as np
import pytest

from cenfuse.pipeline import Run
from cenfuse.sat.pipeline import Pipeline
from cenfuse.sat.util import SatelliteError
from cenfuse.util import CenfuseError, log

# 24 chromosomes and a 7.7-Mb scaffold, then scaffolds falling gradually to 1 kb: 0.8 takes all 669 sequences, 0.9 25
GRADUAL = {f"s{i}": int(n) for i, n in enumerate(np.concatenate([np.linspace(60e6, 15e6, 24), [7.68e6],
                                                                   np.geomspace(1.6e6, 1e3, 644)]))}
SET_09 = [f"s{i}" for i in range(25)]
# three chromosomes, then a tenfold drop to two scaffolds: both rules take the three
CLEAN = {"c1": 300_000, "c2": 300_000, "c3": 300_000, "s1": 5_000, "s2": 5_000}
RECORD = "satellite/work/chromosome_set.json"
HINT = "rerun with --force satellite"


def _assembly(d, name, lengths):
    """d/name.fa (a stub: the run reads only the .fai) and d/name.fa.fai with these sequence lengths."""
    (d / f"{name}.fa").write_text("".join(f">{c}\nA\n" for c in lengths))
    (d / f"{name}.fa.fai").write_text("".join(f"{c}\t{n}\t{i}\t{n}\t{n + 1}\n"
                                              for i, (c, n) in enumerate(lengths.items())))


def _earlier_summary(d, genome, chroms):
    """An earlier run's satellite outputs for --satellite-dir: a centromere summary with a row per sequence."""
    d.mkdir()
    (d / f"{genome}_centromere_summary.tsv").write_text(
        "chrom\tchrom_len\tcen_start\tcen_end\tevidence\tconfidence\n" + "".join(f"{c}\t1000\t.\t.\tnone\tnone\n"
                                                                                for c in chroms))
    return d


class Engine(list):
    """One entry per run of the stand-in engine: the chromosome-set record on disk while it ran (None: there was none).
    With fail_after_seeds set, the runs write the seeds checkpoint and then fail, as an engine does that dies in the
    stage after the seeds."""
    fail_after_seeds = False


@pytest.fixture
def engine(monkeypatch):
    """The satellite engine's run replaced by a stand-in that writes the seeds checkpoint, with the real fingerprint so
    the engine's own resume check accepts it."""
    runs = Engine()

    def run(self):
        record = self.work / "chromosome_set.json"
        runs.append(record.read_text() if record.exists() else None)
        (self.work / "seeds.done").write_text(json.dumps(self.fingerprint("seeds")))
        if runs.fail_after_seeds:
            raise SatelliteError("families: tool crashed")

    monkeypatch.setattr(Pipeline, "run", run)
    return runs


def _record(out) -> dict:
    return json.loads((out / RECORD).read_text())


def _snapshot(out) -> dict:
    """Every file and link under out with its content (a link: its target)."""
    return {str(p.relative_to(out)): os.readlink(p) if p.is_symlink() else p.read_bytes()
            for p in sorted(out.rglob("*")) if p.is_symlink() or p.is_file()}


def test_same_set_resumes(tmp_path, engine):
    _assembly(tmp_path, "S", GRADUAL)
    out = tmp_path / "out"
    Run(tmp_path, out).stage_satellite()
    assert _record(out) == {"S": sorted(SET_09)}                     # {genome: sorted names}
    assert (out / "satellite/work/seeds.done").exists()
    before = (out / RECORD).read_text()
    Run(tmp_path, out).stage_satellite()                             # the checkpoints stand on the same 25: resumes
    assert len(engine) == 2 and (out / RECORD).read_text() == before


def test_changed_set_stops_with_the_force_hint(tmp_path, engine):
    _assembly(tmp_path, "S", GRADUAL)
    out = tmp_path / "out"
    Run(tmp_path, out).stage_satellite()
    before = (out / RECORD).read_text()
    (tmp_path / "S.chromosomes.txt").write_text("".join(f"s{i}\n" for i in range(20)))      # 5 chromosomes fewer
    with pytest.raises(CenfuseError) as e:
        Run(tmp_path, out).stage_satellite()
    msg = str(e.value)
    assert "S: 25 chromosomes when they were made, 20 now (S.chromosomes.txt): 5 dropped, 0 added" in msg
    assert str(out / "satellite/work") in msg and msg.endswith(HINT)
    assert len(engine) == 1 and (out / RECORD).read_text() == before         # the engine did not run; nothing recorded


def test_same_count_but_other_chromosomes_is_a_change(tmp_path, engine):
    """A list that swaps one chromosome for another keeps the count: the message must still say what differs."""
    _assembly(tmp_path, "S", GRADUAL)
    out = tmp_path / "out"
    Run(tmp_path, out).stage_satellite()
    (tmp_path / "S.chromosomes.txt").write_text("".join(f"s{i}\n" for i in [*range(24), 30]))
    with pytest.raises(CenfuseError, match=r"25 chromosomes when they were made, 25 now .*: 1 dropped, 1 added"):
        Run(tmp_path, out).stage_satellite()


@pytest.mark.parametrize("force", ["families", "plot", "telomere"])
def test_a_force_that_keeps_the_seeds_does_not_allow_the_change(tmp_path, engine, force):
    """--force families (or a later stage, or a stage of cenfuse itself) resumes from the old seeds. The stop comes
    before the engine is set up: its force would delete the checkpoints of the stages it redoes."""
    _assembly(tmp_path, "S", GRADUAL)
    out = tmp_path / "out"
    Run(tmp_path, out).stage_satellite()
    (out / "satellite/work/families.done").write_text("{}")
    (tmp_path / "S.chromosomes.txt").write_text("s0\ns1\ns2\n")
    before = _snapshot(out)
    with pytest.raises(CenfuseError, match=HINT):
        Run(tmp_path, out, force=force).stage_satellite()
    assert len(engine) == 1 and _snapshot(out) == before             # a stopped run leaves the folder as it was


def test_a_genome_added_or_removed_is_a_change(tmp_path, engine):
    _assembly(tmp_path, "S", CLEAN)
    out = tmp_path / "out"
    Run(tmp_path, out).stage_satellite()
    _assembly(tmp_path, "T", CLEAN)                                  # added to the genome directory
    with pytest.raises(CenfuseError) as e:
        Run(tmp_path, out).stage_satellite()
    assert str(e.value).endswith("chromosome set: T: not among the genomes of the checkpoints, 3 chromosomes now; "
                                 + HINT)                             # S is unchanged and not named
    Run(tmp_path, out, force="satellite").stage_satellite()          # redone with both
    assert set(_record(out)) == {"S", "T"}
    for f in ("T.fa", "T.fa.fai"):                                   # removed again
        (tmp_path / f).unlink()
    with pytest.raises(CenfuseError, match=r"T: 3 chromosomes when they were made, no longer in .*; " + HINT):
        Run(tmp_path, out).stage_satellite()


def test_pre09_folder_same_set_resumes_and_records(tmp_path, engine):
    """A folder made before 0.9 has seeds.done and no record. Its set was the 0.8 rule's, here the same as 0.9's:
    resume, and record the set (before the engine runs, see the failure tests below) so the next run reads it."""
    _assembly(tmp_path, "S", CLEAN)
    out = tmp_path / "out"
    Run(tmp_path, out).stage_satellite()
    (out / RECORD).unlink()
    Run(tmp_path, out).stage_satellite()
    assert _record(out) == {"S": ["c1", "c2", "c3"]}
    assert [json.loads(r) for r in engine] == [{"S": ["c1", "c2", "c3"]}] * 2        # on disk when the engine ran


def test_pre09_folder_with_a_gradual_tail_stops(tmp_path, engine):
    """0.8 took all 669 sequences of this assembly, 0.9 takes 25: the old seeds and families stand on the 669."""
    _assembly(tmp_path, "S", GRADUAL)
    out = tmp_path / "out"
    Run(tmp_path, out).stage_satellite()
    (out / RECORD).unlink()
    with pytest.raises(CenfuseError) as e:
        Run(tmp_path, out).stage_satellite()
    assert ("S: 669 chromosomes when they were made (cenfuse 0.8's rule), 25 now (length rule): 644 dropped, "
            "0 added") in str(e.value) and str(e.value).endswith(HINT)
    assert len(engine) == 1 and not (out / RECORD).exists()


def test_pre09_folder_with_a_chromosome_list_that_differs_from_the_08_set_stops(tmp_path, engine):
    _assembly(tmp_path, "S", CLEAN)
    out = tmp_path / "out"
    Run(tmp_path, out).stage_satellite()
    (out / RECORD).unlink()
    (tmp_path / "S.chromosomes.txt").write_text("c1\nc2\n")
    with pytest.raises(CenfuseError) as e:
        Run(tmp_path, out).stage_satellite()
    assert "S: 3 chromosomes when they were made (cenfuse 0.8's rule), 2 now (S.chromosomes.txt)" in str(e.value)


def test_pre09_folder_of_a_classify_from_target_resumes(tmp_path, engine):
    """A --classify-from target keeps the 0.8 rule (inputs.discover), so its old set is its current one."""
    _assembly(tmp_path, "S", CLEAN)
    _assembly(tmp_path, "T", GRADUAL)
    out = tmp_path / "out"
    Run(tmp_path, out, classify_from="T=S").stage_satellite()
    (out / RECORD).unlink()
    Run(tmp_path, out, classify_from="T=S").stage_satellite()
    assert len(_record(out)["T"]) == 669 and len(_record(out)["S"]) == 3


@pytest.mark.parametrize("force", ["satellite", "all", "seeds"])
def test_force_satellite_proceeds_and_rewrites_the_record(tmp_path, engine, force):
    _assembly(tmp_path, "S", GRADUAL)
    out = tmp_path / "out"
    Run(tmp_path, out).stage_satellite()
    (tmp_path / "S.chromosomes.txt").write_text("s0\ns1\ns2\n")
    Run(tmp_path, out, force=force).stage_satellite()                # the seeds are redone on the new set
    assert len(engine) == 2 and _record(out) == {"S": ["s0", "s1", "s2"]}
    Run(tmp_path, out).stage_satellite()                             # and the folder now resumes on it
    assert len(engine) == 3


def test_an_engine_that_cannot_be_set_up_leaves_the_record(tmp_path, engine):
    """The record follows the engine's construction, which applies the force wipe and the stale-fingerprint check: a
    run that stops before it has wiped nothing, so the record must still describe the checkpoints on disk (here a
    seed engine that does not exist stops the forced run)."""
    _assembly(tmp_path, "S", GRADUAL)
    out = tmp_path / "out"
    Run(tmp_path, out).stage_satellite()
    before = (out / RECORD).read_text()
    (tmp_path / "S.chromosomes.txt").write_text("s0\ns1\ns2\n")
    with pytest.raises(CenfuseError, match="satellite engine: seed engine 'nope'"):
        Run(tmp_path, out, seed_engine="nope", force="satellite").stage_satellite()
    assert (out / RECORD).read_text() == before and (out / "satellite/work/seeds.done").exists()
    with pytest.raises(CenfuseError, match=HINT):                    # the old seeds are still guarded
        Run(tmp_path, out).stage_satellite()


def test_a_first_run_that_fails_after_the_seeds_resumes_on_the_same_set(tmp_path, engine):
    """The engine dies in the stage after the seeds: seeds.done stands, made on the 25. The rerun on the same set must
    resume, not fall back to the 0.8 comparison (669 vs 25 here) as for a folder made before 0.9."""
    _assembly(tmp_path, "S", GRADUAL)
    out = tmp_path / "out"
    engine.fail_after_seeds = True
    with pytest.raises(CenfuseError, match="satellite engine: families: tool crashed"):
        Run(tmp_path, out).stage_satellite()
    engine.fail_after_seeds = False
    Run(tmp_path, out).stage_satellite()                             # resumes: no stop
    assert len(engine) == 2 and _record(out) == {"S": sorted(SET_09)}


def test_a_forced_run_that_fails_after_the_seeds_still_stops_a_revert(tmp_path, engine):
    """Set A finished. A forced run on set B redid the seeds and died in the next stage; the user goes back to A. The
    seeds on disk were made on B: the guard must stop, not resume on them (the old record said A)."""
    _assembly(tmp_path, "S", GRADUAL)
    out = tmp_path / "out"
    Run(tmp_path, out).stage_satellite()                             # A: the rule's 25, finished
    (tmp_path / "S.chromosomes.txt").write_text("s0\ns1\ns2\n")      # B
    engine.fail_after_seeds = True
    with pytest.raises(CenfuseError, match="satellite engine: families: tool crashed"):
        Run(tmp_path, out, force="satellite").stage_satellite()
    (tmp_path / "S.chromosomes.txt").unlink()                        # back to A
    engine.fail_after_seeds = False
    with pytest.raises(CenfuseError) as e:
        Run(tmp_path, out).stage_satellite()
    assert "S: 3 chromosomes when they were made, 25 now (length rule): 0 dropped, 22 added" in str(e.value)
    assert str(e.value).endswith(HINT) and len(engine) == 2          # the engine did not run again
    Run(tmp_path, out, force="satellite").stage_satellite()          # the way out
    assert len(engine) == 3 and _record(out) == {"S": sorted(SET_09)}


def test_a_record_without_a_seeds_checkpoint_is_not_checked(tmp_path, engine):
    """A run stopped inside the seeds stage leaves the old record and no seeds.done: nothing was checkpointed on the
    old set, so the new set goes ahead."""
    _assembly(tmp_path, "S", GRADUAL)
    out = tmp_path / "out"
    Run(tmp_path, out).stage_satellite()
    (out / "satellite/work/seeds.done").unlink()
    (tmp_path / "S.chromosomes.txt").write_text("s0\ns1\ns2\n")
    Run(tmp_path, out).stage_satellite()
    assert _record(out) == {"S": ["s0", "s1", "s2"]}


def test_an_unreadable_record_stops_with_the_force_hint(tmp_path, engine):
    _assembly(tmp_path, "S", CLEAN)
    out = tmp_path / "out"
    Run(tmp_path, out).stage_satellite()
    for text in ('{"S": ["c1", ', "[]", '{"S": 3}'):                 # cut off, not a mapping, not a list of names
        (out / RECORD).write_text(text)
        with pytest.raises(CenfuseError, match=r"chromosome_set\.json: cannot be read.*" + HINT):
            Run(tmp_path, out).stage_satellite()
    assert len(engine) == 1


def test_satellite_dir_import_on_another_set_warns(tmp_path, caplog, monkeypatch):
    """--satellite-dir takes the earlier calls as they are, with a warning when they were made on other sequences
    (here all 669, the 0.8 set); nothing is recorded."""
    _assembly(tmp_path, "S", GRADUAL)
    earlier = _earlier_summary(tmp_path / "earlier", "S", list(GRADUAL))
    monkeypatch.setattr(log, "propagate", True)
    with caplog.at_level(logging.WARNING, logger="cenfuse"):
        got = Run(tmp_path, tmp_path / "out", satellite_dir=earlier).stage_satellite()
    assert got == earlier.resolve()
    warned = [r.getMessage() for r in caplog.records if r.levelno == logging.WARNING]
    assert len(warned) == 1
    assert "S: " in warned[0] and "S_centromere_summary.tsv" in warned[0] and "669" in warned[0]
    assert "25 chromosomes" in warned[0] and "644 rows on sequences outside the set" in warned[0]
    assert "another chromosome set" in warned[0] and "--satellite-dir" in warned[0]
    assert not (tmp_path / "out" / RECORD).exists()


def test_satellite_dir_import_with_chromosomes_missing_warns(tmp_path, caplog, monkeypatch):
    _assembly(tmp_path, "S", GRADUAL)
    earlier = _earlier_summary(tmp_path / "earlier", "S", SET_09[:20])
    monkeypatch.setattr(log, "propagate", True)
    with caplog.at_level(logging.WARNING, logger="cenfuse"):
        Run(tmp_path, tmp_path / "out", satellite_dir=earlier).stage_satellite()
    warned = [r.getMessage() for r in caplog.records if r.levelno == logging.WARNING]
    assert len(warned) == 1 and "5 chromosomes without a row" in warned[0] and "0 rows" in warned[0]


def test_satellite_dir_import_on_the_same_set_is_silent(tmp_path, caplog, monkeypatch):
    _assembly(tmp_path, "S", GRADUAL)
    earlier = _earlier_summary(tmp_path / "earlier", "S", reversed(SET_09))      # row order does not matter
    monkeypatch.setattr(log, "propagate", True)
    with caplog.at_level(logging.WARNING, logger="cenfuse"):
        Run(tmp_path, tmp_path / "out", satellite_dir=earlier).stage_satellite()
    assert not [r for r in caplog.records if r.levelno >= logging.WARNING]


def test_satellite_dir_import_without_a_usable_summary_warns_nothing(tmp_path, caplog, monkeypatch):
    """A missing or foreign summary is not this check's to report: seqcalls.read stops the run with the reason when it
    reads the file."""
    _assembly(tmp_path, "S", CLEAN)
    earlier = tmp_path / "earlier"
    earlier.mkdir()
    (earlier / "S_centromere_summary.tsv").write_text("not\ta\tsummary\n")
    monkeypatch.setattr(log, "propagate", True)
    with caplog.at_level(logging.WARNING, logger="cenfuse"):
        Run(tmp_path, tmp_path / "out", satellite_dir=earlier).stage_satellite()
        (earlier / "S_centromere_summary.tsv").unlink()
        Run(tmp_path, tmp_path / "out", satellite_dir=earlier).stage_satellite()
    assert not [r for r in caplog.records if r.levelno >= logging.WARNING]
