import logging
import re
import shutil
import subprocess
import sys
from dataclasses import replace

import pandas as pd
import pytest

from cenfuse.sat.params import Params
from cenfuse.sat.util import SatelliteError
from sat_common import TEST_PARAMS, have
from sat_gffcheck import check_gff
from sat_synth import build_complex

TOOLS = ("FasTAN", "FAtoGDB", "ONEview", "makeblastdb", "blastn", "mafft", "cd-hit-est", "minimap2")
pytestmark = pytest.mark.skipif(not have(*TOOLS), reason="external tools missing")


@pytest.fixture(scope="module", params=["fastan", "trf"])
def synth_run(request, tmp_path_factory):
    """The whole pipeline on the synthetic complex, once per seed engine: every truth below holds for both."""
    from cenfuse.sat.pipeline import Pipeline
    from cenfuse.sat.seeds import ENGINES
    if not have(*ENGINES[request.param].TOOLS):
        pytest.skip(f"{request.param} tools missing")
    root = tmp_path_factory.mktemp(f"e2e_{request.param}")
    cdir = root / "Synth"
    truth = build_complex(cdir)
    p = Params(**dict(TEST_PARAMS, seed_engine=request.param))
    Pipeline(cdir, 4, root / "cache", p).run()
    return cdir, truth, root, p


def _summary(cdir, g):
    return pd.read_csv(cdir / "centromere" / f"{g}_centromere_summary.tsv", sep="\t", dtype=str).set_index("chrom")


def _gff(cdir, g):
    from cenfuse.sat.genomes import load_fai
    return check_gff(cdir / "centromere" / f"{g}_centromeric_repeat.gff", load_fai(cdir / f"{g}.fa"))


# satellite_families.tsv columns, written out: comparing with FAMTAB_COLS itself could not fail if cost columns were
# added to it
FAMTAB_COLS_WRITTEN = ["family", "fam_id", "superfamily", "consensus_len", "hmm_mode", "rdna", "genome", "class", "role",
                       "strength", "defines", "bp", "n_arrays", "prevalence", "concentration", "subtel_frac", "gene_ratio",
                       "median_main_bp", "eligible", "dominant", "in_ltr", "ltr_clade", "tandem_frac", "whole_frac",
                       "locus_fold", "dominant_satellite", "te_rejected", "reason"]


def test_outputs_exist_and_gffs_valid(synth_run):
    from cenfuse.sat.annotate import COST_COLS
    cdir, truth, _, _ = synth_run
    out = cdir / "centromere"
    for f in ("satellite_families.tsv", "satellite_families.fa", "satellite_families.hmm", "centromere_overview.pdf"):
        assert (out / f).stat().st_size > 0
    for g in truth["genomes"]:
        assert _gff(cdir, g)
        assert (out / "work" / "annotate" / f"{g}.search_cost.tsv").stat().st_size > 0       # a fresh run writes it
    # CPU times differ between runs, so they stay out of the results tables (README, Reproducible)
    cols = pd.read_csv(out / "satellite_families.tsv", sep="\t", nrows=0).columns.tolist()
    assert cols == FAMTAB_COLS_WRITTEN and not set(cols) & (set(COST_COLS) - {"family"})


def test_search_cost_table_of_a_run(synth_run):
    """work/annotate/<genome>.search_cost.tsv: one row per family, its monomers add up to the monomers table."""
    from cenfuse.sat.annotate import COST_COLS
    cdir, truth, _, _ = synth_run
    work = cdir / "centromere" / "work"
    fam_ids = set(pd.read_csv(work / "families" / "families.tsv", sep="\t").fam_id)
    for g in truth["genomes"]:
        cost = pd.read_csv(work / "annotate" / f"{g}.search_cost.tsv", sep="\t")
        assert cost.columns.tolist() == COST_COLS and len(cost) == len(fam_ids) and set(cost.family) == fam_ids
        assert cost.monomers_kept.sum() == len(pd.read_csv(work / "annotate" / f"{g}.monomers.tsv.gz", sep="\t"))
        assert abs(cost.cpu_share.sum() - 1) < 1e-6


def test_output_stage_logs_the_costliest_searches(synth_run, caplog, monkeypatch):
    """The output stage names each genome's costliest searches by final family name and role, those of
    satellite_families.tsv (only the output and plot stages run again here)."""
    from cenfuse.sat.pipeline import Pipeline
    from cenfuse.sat.util import log
    cdir, truth, root, p = synth_run
    monkeypatch.setattr(log, "propagate", True)
    with caplog.at_level(logging.INFO, logger="cenfuse-sat"):
        Pipeline(cdir, 2, root / "cache", p, force="output").run()
    msgs = [r.getMessage() for r in caplog.records]
    fam = pd.read_csv(cdir / "centromere" / "satellite_families.tsv", sep="\t")
    for g in truth["genomes"]:
        line = [m for m in msgs if m.startswith(f"output: {g}: costliest searches: ")]
        assert len(line) == 1, (g, msgs)
        top = re.findall(r"(\S+) \((\w+), (\d+)%\)", line[0])
        role = dict(zip(fam[fam.genome == g].family, fam[fam.genome == g].role))
        assert 1 <= len(top) <= 3 and all(role[name] == r for name, r, _ in top), line
        shares = [int(s) for _, _, s in top]
        assert shares == sorted(shares, reverse=True)


def test_centromere_boundaries_and_none(synth_run):
    cdir, truth, _, _ = synth_run
    for g, chroms in truth["genomes"].items():
        s = _summary(cdir, g)
        for c, t in chroms.items():
            if t["cen"] is None:
                assert s.loc[c, "evidence"] == "none"
                continue
            assert abs(int(s.loc[c, "cen_start"]) - 1 - t["cen"][0]) <= 150, (g, c)
            assert abs(int(s.loc[c, "cen_end"]) - t["cen"][1]) <= 150, (g, c)


def test_decoy_family_is_other(synth_run):
    cdir, _, _, _ = synth_run
    f = pd.read_csv(cdir / "centromere" / "satellite_families.tsv", sep="\t")
    dec = f[(f.consensus_len - 90).abs() <= 3]
    assert len(dec) and set(dec["class"]) == {"other"}
    cen = f[(f.consensus_len - 150).abs() <= 3]
    assert set(cen["class"]) == {"centromeric"}


def test_tandem_hor_chr1_and_duplication_chr2(synth_run):
    cdir, _, _, _ = synth_run
    u = pd.read_csv(cdir / "centromere" / "Pa_A_hor_units.tsv", sep="\t")
    c1 = u[u.chrom == "Pa_A01"]
    assert ((c1.n_monomers == 4) & (c1.max_tandem_copies >= 14)).any()
    c2 = u[u.chrom == "Pa_A02"]
    assert ((c2.n_monomers >= 2) & (c2.occurrences >= 2) & (c2.max_tandem_copies == 1)).any()


def test_chr3_gap_microsat_and_long_unit(synth_run):
    cdir, _, _, p = synth_run
    assert "gaps:1" in _summary(cdir, "Pa_A").loc["Pa_A03", "flags"]
    tr = [a for c, s, e, t, a in _gff(cdir, "Pa_A") if t == "tandem_repeat" and c == "Pa_A03"]
    assert any(a["period"] == "2" and a["unit_source"] == p.seed_engine for a in tr)
    assert any(a["unit_source"] == "selfaln" and abs(int(a["period"]) - 3000) <= 150 for a in tr)


def test_resume_and_param_change(synth_run, monkeypatch):
    """Review I1: a checkpoint is reused only if the code version and the parameters of that stage and all earlier stages
    are unchanged. A stale search stage (seeds, families, annotate) stops the run and --force must name it (or an earlier
    one); a stale stage from classify on (cheap: a parameter, its code, or the genes/ and tes/ inputs) is recomputed from
    there, with a warning that says what changed (0.14, future.md #59)."""
    import cenfuse.sat.pipeline as pl
    from cenfuse.sat.pipeline import STAGES, Pipeline
    cdir, _, root, same = synth_run
    work = cdir / "centromere" / "work"
    Pipeline(cdir, 2, root / "cache", same).run()                          # all checkpoints: no-op
    changed = replace(same, cluster_gap=60_000)                            # a classify parameter: recomputed from classify
    pipe = Pipeline(cdir, 2, root / "cache", changed)
    assert pipe.force == "classify"
    assert {st for st in STAGES if (work / f"{st}.done").exists()} == {"seeds", "families", "annotate"}
    Pipeline(cdir, 2, root / "cache", same).run()                          # back to the fixture's outputs
    monkeypatch.setitem(pl.STAGE_VERSION, "classify", "9.9.9")             # classify code changed: recomputed too
    assert Pipeline(cdir, 2, root / "cache", same, force="plot").force == "classify"   # a later force cannot skip it
    Pipeline(cdir, 2, root / "cache", same).run()
    monkeypatch.undo()
    Pipeline(cdir, 2, root / "cache", same).run()
    monkeypatch.setitem(pl.STAGE_VERSION, "seeds", "9.9.9")                # seed code changed: the search is not redone unasked
    with pytest.raises(SatelliteError, match="--force seeds"):
        Pipeline(cdir, 2, root / "cache", same)
    monkeypatch.undo()
    other = replace(same, seed_engine={"fastan": "trf", "trf": "fastan"}[same.seed_engine])
    with pytest.raises(SatelliteError, match="--force seeds"):               # checkpoints of the other engine are stale
        Pipeline(cdir, 2, root / "cache", other)
    assert all((work / f"{st}.done").exists() for st in STAGES)             # a stop deletes nothing


def test_a_hard_linked_copy_is_never_rewritten(synth_run, tmp_path):
    """A folder copied with hard links (cp -al) and recomputed from the search on must leave the original's files as they
    were: each recomputed stage unlinks its folder's shared files first (its tools write in place), and every output is
    written to a temporary file and moved over its target (future.md #59, found by the 0.14 verification pass)."""
    import hashlib
    import os
    from cenfuse.sat.pipeline import Pipeline
    cdir, _, root, same = synth_run
    copy = tmp_path / "Synth"
    shutil.copytree(cdir, copy, copy_function=os.link, symlinks=True)
    files = [p for p in cdir.rglob("*") if p.is_file() and not p.is_symlink()]
    before = {p: hashlib.md5(p.read_bytes()).hexdigest() for p in files}
    Pipeline(copy, 2, root / "cache", same, force="families").run()
    changed = [str(p.relative_to(cdir)) for p in files if hashlib.md5(p.read_bytes()).hexdigest() != before[p]]
    assert not changed, changed
    assert (copy / "centromere" / "satellite_families.tsv").stat().st_nlink == 1     # the copy has its own outputs now


def test_every_array_is_reported_with_family_and_role(synth_run):
    """<genome>_satellite_arrays.gff holds every array of every kept family, with the satellite_families.tsv name and
    the family's role in that genome."""
    from cenfuse.sat.classify import ROLES
    from cenfuse.sat.genomes import load_fai
    cdir, truth, _, _ = synth_run
    out = cdir / "centromere"
    fam = pd.read_csv(out / "satellite_families.tsv", sep="\t")
    assert "role" in fam.columns and set(fam.role) <= set(ROLES) | {"absent", "unassessed"}
    for g in truth["genomes"]:
        gff = out / f"{g}_satellite_arrays.gff"
        assert check_gff(gff, load_fai(cdir / f"{g}.fa"))
        lines = [l.split("\t") for l in gff.read_text().splitlines() if not l.startswith("#")]
        n = len(pd.read_csv(out / "work" / "annotate" / f"{g}.arrays.tsv", sep="\t"))
        assert len(lines) == n > 0 and all(l[2] == "satellite_DNA" for l in lines)
        att = [dict(kv.split("=", 1) for kv in l[8].split(";")) for l in lines]
        assert {a["family"] for a in att} <= set(fam.family)
        role = dict(zip(fam[fam.genome == g].family, fam[fam.genome == g].role))
        assert all(a["role"] == role[a["family"]] for a in att)


def test_rdna_column_in_family_table(synth_run):
    cdir, _, _, _ = synth_run
    f = pd.read_csv(cdir / "centromere" / "satellite_families.tsv", sep="\t", dtype=str)
    assert "rdna" in f.columns and set(f.rdna) == {"."}


def test_bgzipped_inputs_give_identical_results(tmp_path):
    """The cleaned polyploid genomes are bgzip FASTA (.fa.gz + .fai + .gzi): every stage reads them through
    pysam, so results must equal the plain-FASTA run."""
    import pysam
    from cenfuse.sat.pipeline import Pipeline
    src = tmp_path / "src"
    truth = build_complex(src)
    cdir, gdir = tmp_path / "fa" / "Synth", tmp_path / "gz" / "Synth"     # same name: it prefixes family names
    cdir.mkdir(parents=True)
    gdir.mkdir(parents=True)
    for g in truth["genomes"]:
        for ext in (".fa", ".fa.fai"):
            shutil.copy(src / f"{g}{ext}", cdir / f"{g}{ext}")
        gz = gdir / f"{g}.fa.gz"
        pysam.tabix_compress(str(src / f"{g}.fa"), str(gz))
        pysam.faidx(str(gz))
    for d in (cdir, gdir):
        Pipeline(d, 4, tmp_path / "cache", Params(**TEST_PARAMS)).run()
    for g in truth["genomes"]:
        assert _summary(gdir, g).equals(_summary(cdir, g)), g
        assert (gdir / "centromere" / f"{g}_centromeric_repeat.gff").read_text() == \
            (cdir / "centromere" / f"{g}_centromeric_repeat.gff").read_text(), g
    fam = lambda d: pd.read_csv(d / "centromere" / "satellite_families.tsv", sep="\t")
    assert fam(gdir).equals(fam(cdir))


def test_cli_help_is_short():
    """Thresholds are internal constants (params.py), not options."""
    r = subprocess.run([sys.executable, "-m", "cenfuse.sat", "run", "--help"], capture_output=True, text=True)
    assert r.returncode == 0 and "GENOME_DIR" in r.stdout and "--outdir" in r.stdout
    assert "--seed-engine {fastan,trf}" in r.stdout and "--trf-cache" not in r.stdout     # old name: hidden alias
    assert "--cluster-gap" not in r.stdout and r.stdout.count("\n  -") <= 8


def test_cli_needs_a_genome_dir_and_fails_clearly(tmp_path):
    run = lambda *a: subprocess.run([sys.executable, "-m", "cenfuse.sat", "run", *a], capture_output=True, text=True)
    assert run().returncode == 2
    r = run(str(tmp_path / "empty"))
    assert r.returncode == 1 and "no *.fa genomes" in r.stderr
    r = run("--complex-dir", str(tmp_path / "empty"))                      # hidden alias of GENOME_DIR
    assert r.returncode == 1 and "no *.fa genomes" in r.stderr
    r = run(str(tmp_path / "empty"), "--trf-cache", str(tmp_path / "c"))   # hidden alias of --seed-cache
    assert r.returncode == 1 and "no *.fa genomes" in r.stderr
    assert run(str(tmp_path / "empty"), "--seed-engine", "nope").returncode == 2
