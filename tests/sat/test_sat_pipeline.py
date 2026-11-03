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
