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
