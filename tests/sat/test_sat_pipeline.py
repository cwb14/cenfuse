from dataclasses import fields

import pytest

from cenfuse.sat.params import Params
from cenfuse.sat.util import SatelliteError
from cenfuse.sat.pipeline import STAGE_PARAMS, STAGES


def test_every_parameter_belongs_to_exactly_one_stage():
    mapped = [x for s in STAGES for x in STAGE_PARAMS[s]]
    assert sorted(mapped) == sorted(f.name for f in fields(Params)) and len(mapped) == len(set(mapped))


def test_every_parameter_has_a_reason():
    from cenfuse.sat.params import REASON
    assert set(REASON) == {f.name for f in fields(Params)} and all(REASON.values())


def _complex(d, names, nseq=1):
    for n in names:
        chroms = [f"{n}_{i}" for i in range(nseq)]
        (d / f"{n}.fa").write_text("".join(f">{c}\nACGT\n" for c in chroms))
        (d / f"{n}.fa.fai").write_text("".join(f"{c}\t4\t{i}\t4\t5\n" for i, c in enumerate(chroms)))
    return d


@pytest.mark.parametrize("spec, msg", [
    ("X=a", "unknown genome"), ("a=X", "unknown genome"), ("a=a", "borrows from itself"),
    ("a=b,b=c", "is itself borrowed"), ("a=b,a=c", "more than once"), ("a", "TARGET=SOURCE"),
    ("a=", "TARGET=SOURCE"), ("=b", "TARGET=SOURCE"), ("a=b+", "TARGET=SOURCE")])
def test_classify_from_is_validated(tmp_path, spec, msg):
    from cenfuse.sat.pipeline import Pipeline
    cdir = _complex(tmp_path, "abc")
    with pytest.raises(SatelliteError, match=msg):
        Pipeline(cdir, 1, tmp_path / "cache", Params(classify_from=spec))


def test_classify_from_is_normalised_and_sources_come_first(tmp_path):
    from cenfuse.sat.pipeline import Pipeline
    cdir = _complex(tmp_path, ["Ctg", "a", "b"])
    pl = Pipeline(cdir, 1, tmp_path / "cache", Params(classify_from=" Ctg = b+a+b "))
    assert pl.p.classify_from == "Ctg=a+b" and pl.borrow == {"Ctg": ["a", "b"]}
    assert [g.name for g in pl.genomes] == ["a", "b", "Ctg"]


def test_seed_engine_changes_every_fingerprint(tmp_path):
    """Seeds feed every stage: a checkpoint made with one engine is never reused with the other."""
    from cenfuse.sat.pipeline import Pipeline
    cdir = _complex(tmp_path, ["a"])
    a = Pipeline(cdir, 1, tmp_path / "cache", Params())
    b = Pipeline(cdir, 1, tmp_path / "cache", Params(seed_engine="trf"))
    assert all(a.fingerprint(s) != b.fingerprint(s) for s in STAGES)


def test_default_cache_keeps_using_an_older_runs_trf_cache(tmp_path):
    """Runs of satellite engine < 0.9 cached TRF in OUTDIR/work/trf_cache (hours per large genome); a rerun must find it."""
    from cenfuse.sat.pipeline import Pipeline
    cdir = _complex(tmp_path, ["a"])
    assert Pipeline(cdir, 1, None, Params()).cache == cdir / "centromere" / "work" / "seed_cache"
    (cdir / "centromere" / "work" / "trf_cache").mkdir()
    assert Pipeline(cdir, 1, None, Params(seed_engine="trf")).cache == cdir / "centromere" / "work" / "trf_cache"


def test_unknown_seed_engine_fails_clearly(tmp_path):
    from cenfuse.sat.pipeline import Pipeline
    with pytest.raises(SatelliteError, match="seed engine 'nope'"):
        Pipeline(_complex(tmp_path, ["a"]), 1, tmp_path / "cache", Params(seed_engine="nope"))


def test_classify_from_changes_classify_fingerprint_only(tmp_path):
    from cenfuse.sat.pipeline import Pipeline
    cdir = _complex(tmp_path, ["Ctg", "a"])
    a = Pipeline(cdir, 1, tmp_path / "cache", Params())
    b = Pipeline(cdir, 1, tmp_path / "cache", Params(classify_from="Ctg=a"))
    assert a.fingerprint("annotate") == b.fingerprint("annotate")
    assert a.fingerprint("classify") != b.fingerprint("classify")


def test_many_sequences_warn_unless_borrowed(tmp_path, caplog, monkeypatch):
    from cenfuse.sat.pipeline import Pipeline
    from cenfuse.sat.util import log
    monkeypatch.setattr(log, "propagate", True)
    cdir = _complex(tmp_path, ["Ctg", "a"], nseq=60)
    Pipeline(cdir, 1, tmp_path / "cache", Params(classify_from="Ctg=a"))
    assert [r for r in caplog.records if "60 sequences" in r.message and r.message.startswith("a:")]
    assert not [r for r in caplog.records if r.message.startswith("Ctg:")]


def test_standalone_run_names_sequences_too_short_to_be_chromosomes(tmp_path, monkeypatch, capsys):
    """The one-per-chromosome test and family prevalence count every sequence as a chromosome, so `python -m cenfuse.sat
    run` names the sequences cenfuse's chromosome rule leaves out (here a 20-kb organelle beside three yeast-sized
    chromosomes, which turns the satellite test from untested to fail); a contig assembly (--classify-from) is exempt."""
    from cenfuse.sat import cli
    from cenfuse.sat.pipeline import Pipeline
    lengths = {"I": 5_600_000, "II": 4_500_000, "III": 2_500_000, "mt": 20_000}
    for g in ("Y", "Ctg"):
        (tmp_path / f"{g}.fa").write_text("".join(f">{c}\nACGT\n" for c in lengths))
        (tmp_path / f"{g}.fa.fai").write_text("".join(f"{c}\t{L}\t0\t4\t5\n" for c, L in lengths.items()))
    from cenfuse.sat.util import log
    monkeypatch.setattr(Pipeline, "run", lambda self: None)
    monkeypatch.setattr(log, "handlers", list(log.handlers))  # main() sets up the engine's logger: restored after
    monkeypatch.setattr(log, "propagate", log.propagate)
    monkeypatch.setattr(log, "level", log.level)
    assert cli.main(["run", str(tmp_path), "-t", "1", "--classify-from", "Ctg=Y"]) == 0
    err = capsys.readouterr().err
    assert "Y: 1 of 4 sequences" in err and "'mt'" in err and "Ctg:" not in err


def test_output_stage_names_the_costliest_searches_with_their_final_role(tmp_path, caplog, monkeypatch):
    """The three searches that cost most CPU, by the name and role the family ends with: a dispersed family that cost
    much decides nothing, and one without arrays in this genome is 'absent' (as in satellite_families.tsv). A folder
    made before the cost table existed, and a run with no families (a table of its header only), log nothing."""
    import logging
    import pandas as pd
    from cenfuse.sat.annotate import MON_COLS, search_cost_table
    from cenfuse.sat.pipeline import Pipeline
    from cenfuse.sat.util import log
    monkeypatch.setattr(log, "propagate", True)
    pl = Pipeline(_complex(tmp_path, ["a"]), 1, tmp_path / "cache", Params())
    g = pl.genomes[0]
    fams = pd.DataFrame({"family": ["fam001", "fam002"], "role": ["centromeric", "dispersed"]})
    names = {f"fam00{i}": n for i, n in enumerate(["S_sat01_150", "S_sat02_90", "S_sat03_3000", "S_sat04_40"], 1)}
    with caplog.at_level(logging.INFO, logger="cenfuse-sat"):
        pl._log_costliest(g, fams, names)
        assert not caplog.records
        (pl.work / "annotate").mkdir()
        table = pl._p("annotate", "a.search_cost.tsv")
        search_cost_table({}, pd.DataFrame(columns=MON_COLS)).to_csv(table, sep="\t", index=False)
        pl._log_costliest(g, fams, names)
        assert not caplog.records
        hours = {"fam001": 1, "fam002": 6, "fam003": 3, "fam004": 0}
        cost = {f: dict(model_len=100, searches=1, cpu_s=3600.0 * h, wall_s=3600.0 * h, hits=0)
                for f, h in hours.items()}
        search_cost_table(cost, pd.DataFrame(columns=MON_COLS)).to_csv(table, sep="\t", index=False)
        pl._log_costliest(g, fams, names)
    assert [r.getMessage() for r in caplog.records] == [
        "output: a: costliest searches: S_sat02_90 (dispersed, 60%), S_sat03_3000 (absent, 30%), "
        "S_sat01_150 (centromeric, 10%)"]


@pytest.mark.parametrize("table, with_role", [
    ("", True),                                                  # an empty file (a run killed while writing it)
    ("family\tmodel_len\nfam001\t100\n", True),                  # no cpu_share column
    ("family\tcpu_share\nfam001\tmost\n", True),                 # a share that is not a number
    ("family\tcpu_share\nfam001\t0.5\n", False),                 # a families table without roles
    (None, True)],                                               # a folder where the table should be
    ids=["empty", "no-cpu-share", "share-not-a-number", "no-roles", "a-folder"])
def test_a_cost_table_or_role_that_cannot_be_read_only_warns(tmp_path, caplog, monkeypatch, table, with_role):
    """The costliest-search line is a diagnostic: a malformed cost table, or families without roles, give one warning
    that names the table, and the output stage goes on."""
    import logging
    import pandas as pd
    from cenfuse.sat.pipeline import Pipeline
    from cenfuse.sat.util import log
    monkeypatch.setattr(log, "propagate", True)
    pl = Pipeline(_complex(tmp_path, ["a"]), 1, tmp_path / "cache", Params())
    (pl.work / "annotate").mkdir()
    path = pl._p("annotate", "a.search_cost.tsv")
    path.mkdir() if table is None else path.write_text(table)
    fams = pd.DataFrame({"family": ["fam001"], **({"role": ["centromeric"]} if with_role else {})})
    with caplog.at_level(logging.INFO, logger="cenfuse-sat"):
        pl._log_costliest(pl.genomes[0], fams, {"fam001": "a_sat01_100"})
    assert [r.levelname for r in caplog.records] == ["WARNING"]
    assert "a.search_cost.tsv" in caplog.records[0].getMessage()


def test_the_output_stage_completes_with_a_malformed_cost_table(tmp_path, caplog, monkeypatch):
    """A run whose cost table is malformed (no cpu_share column) still finishes the output stage: output.done is
    written and the plot stage runs. The stages up to hor are checkpointed, and the output writers and the plot are
    stand-ins, so only the stage's own steps run."""
    import json
    import logging
    import pandas as pd
    from cenfuse.sat import output, pipeline
    from cenfuse.sat.pipeline import Pipeline
    from cenfuse.sat.util import log
    monkeypatch.setattr(log, "propagate", True)
    pl = Pipeline(_complex(tmp_path, ["a"]), 1, tmp_path / "cache", Params())
    for st in STAGES[:STAGES.index("output")]:
        (pl.work / st).mkdir()
        (pl.work / f"{st}.done").write_text(json.dumps(pl.fingerprint(st)))
    pl._p("annotate", "a.search_cost.tsv").write_text("family\tmodel_len\nfam001\t100\n")
    plotted = []
    monkeypatch.setattr(pipeline, "require_tools", lambda names: None)
    monkeypatch.setattr(pl, "fam_table", lambda: pd.DataFrame({"fam_id": ["fam001"]}))
    monkeypatch.setattr(pl, "genome_data", lambda g: dict(arrays=None, fams=pd.DataFrame(
        {"family": ["fam001"], "role": ["centromeric"]})))
    monkeypatch.setattr(output, "final_names", lambda fam, arrays, name: {"fam001": "a_sat01_100"})
    monkeypatch.setattr(output, "write_genome", lambda *a: None)
    monkeypatch.setattr(output, "write_families", lambda *a: None)
    monkeypatch.setattr(Pipeline, "stage_plot", lambda self: plotted.append(1))
    with caplog.at_level(logging.INFO, logger="cenfuse-sat"):
        pl.run()
    assert (pl.work / "output.done").exists() and (pl.work / "plot.done").exists() and plotted
    warned = [r.getMessage() for r in caplog.records if r.levelno == logging.WARNING]
    assert len(warned) == 1 and "a.search_cost.tsv" in warned[0]


def test_tables_read_back_keep_sequence_names_verbatim(tmp_path):
    """Sequence names pandas reads as missing (NA, None, nan, null) or as numbers (01) come back verbatim from the tables
    the engine writes and reads again, and so do the chromosomes a family defines."""
    import pandas as pd
    from cenfuse.sat.pipeline import Pipeline
    pl = Pipeline(_complex(tmp_path, ["a"]), 1, tmp_path / "cache", Params())
    g = pl.genomes[0]
    names = ["NA", "None", "nan", "null", "N/A", "01"]
    t = pd.DataFrame(dict(chrom=names, start=range(6), consensus="ACGT"))
    for stage, f in (("annotate", "a.monomers.tsv.gz"), ("annotate", "a.arrays.tsv"), ("seeds", "a.seeds.tsv.gz")):
        (pl.work / stage).mkdir(exist_ok=True)
        t.to_csv(pl._p(stage, f), sep="\t", index=False)
    (pl.work / "centromeres").mkdir()
    pd.DataFrame(dict(family=["fam001", "fam002"], defines=["01", "NA"], bp=[1, 2])).to_csv(
        pl._p("centromeres", "a.families.tsv"), sep="\t", index=False)
    for read in (pl.monomers, pl.arrays, pl.seed_table):
        assert list(read(g).chrom) == names
    assert list(pl.fams(g).defines) == ["01", "NA"]


def test_a_moved_output_folder_finds_its_seed_tables(tmp_path):
    """The seed table is a link into the seed cache. With the cache inside the output folder (the default) the link is
    relative, so a moved folder keeps it; an absolute link that a move broke is re-linked to the cache entry of the same
    key; with that entry gone too, the error names the way out."""
    import gzip
    import os
    import shutil
    from cenfuse.sat import seeds
    from cenfuse.sat.pipeline import Pipeline
    cdir = _complex(tmp_path, ["a"])
    pl = Pipeline(cdir, 1, None, Params(), out_dir=tmp_path / "out1")
    g = pl.genomes[0]
    key = seeds.cache_path(pl.cache, seeds.genome_key(g), pl.p)
    key.parent.mkdir(parents=True)
    with gzip.open(key, "wt") as f:
        f.write("chrom\tstart\tend\tconsensus\na_0\t0\t4\tAC\n")
    (pl.work / "seeds").mkdir()
    pl.stage_seeds()                                                    # cached: no seed engine runs
    link = pl._p("seeds", "a.seeds.tsv.gz")
    assert not os.path.isabs(os.readlink(link))
    shutil.move(tmp_path / "out1", tmp_path / "out2")
    pl = Pipeline(cdir, 1, None, Params(), out_dir=tmp_path / "out2")
    assert list(pl.seed_table(g).chrom) == ["a_0"]
    link = pl._p("seeds", "a.seeds.tsv.gz")
    link.unlink()
    link.symlink_to(tmp_path / "out1" / "work" / "seed_cache" / key.name)  # an absolute link into the old place
    assert list(pl.seed_table(g).chrom) == ["a_0"] and link.exists()
    (pl.cache / key.name).unlink()
    link.unlink()
    link.symlink_to(tmp_path / "out1" / "work" / "seed_cache" / key.name)
    with pytest.raises(SatelliteError, match="seed table of a not found.*--seed-cache DIR or --force seeds"):
        pl.seed_table(g)
