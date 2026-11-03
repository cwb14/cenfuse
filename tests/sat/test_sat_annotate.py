import itertools
import logging
import re
import time
from types import SimpleNamespace

import pandas as pd
import pytest

from cenfuse.sat.annotate import MON_COLS, chain_arrays, drop_multimer_hits, resolve_overlaps
from sat_common import have


def _mon(rows):
    return pd.DataFrame(rows, columns=MON_COLS)


def test_resolve_overlaps_keeps_best_and_adjacent():
    rows = [("c", 0, 150, "+", "f1", 100.0, 1e-20, 1, 150, 150, 0.9),
            ("c", 2, 150, "+", "f2", 80.0, 1e-15, 1, 148, 150, 0.85),   # same monomer, other family: dropped
            ("c", 145, 295, "+", "f1", 90.0, 1e-18, 1, 150, 150, 0.9)]   # 5-bp overlap only: kept
    out = resolve_overlaps(_mon(rows))
    assert out.family.tolist() == ["f1", "f1"] and out.start.tolist() == [0, 145]


def test_chain_arrays_splits_on_gap_and_strand(test_params):
    rows = [("c", i * 150, (i + 1) * 150, "+", "f", 90.0, 1e-20, 1, 150, 150, 0.9) for i in range(5)]
    rows += [("c", 1150, 1300, "+", "f", 90.0, 1e-20, 1, 150, 150, 0.9)]          # gap 400 > 300: singleton
    rows += [("c", 2000 + i * 150, 2150 + i * 150, "-", "f", 90.0, 1e-20, 1, 150, 150, 0.8) for i in range(3)]
    mon, arr = chain_arrays(_mon(rows), {"f": 150}, test_params)
    assert sorted(arr.n_monomers.tolist()) == [3, 5]
    assert set(arr.strand) == {"+", "-"} and len(mon) == 8


@pytest.mark.skipif(not have("FasTAN", "FAtoGDB", "ONEview", "blastn", "mafft"), reason="tools missing")
def test_annotate_synthetic_centromere(synth_annot, synth_complex, synth_families):
    mon, arr = synth_annot
    fam, _ = synth_families
    cen_fam = fam[(fam.consensus_len - 150).abs() <= 3].fam_id.iloc[0]
    s, e = synth_complex[1]["genomes"]["Pa_A"]["Pa_A01"]["cen"]
    a = arr[(arr.chrom == "Pa_A01") & (arr.family == cen_fam)]
    assert abs(a.start.min() - s) <= 150 and abs(a.end.max() - e) <= 150
    assert (a.strand == "-").any()                                    # the inverted block
    assert len(mon[(mon.chrom == "Pa_A01") & (mon.family == cen_fam)]) >= 250
    assert mon.mon_idx.tolist() == list(range(len(mon)))


@pytest.mark.skipif(not have("FasTAN", "FAtoGDB", "ONEview", "blastn", "mafft"), reason="tools missing")
def test_chunked_search_matches_whole_chromosome(synth_complex, synth_families):
    """Speed fix: nhmmer parallelises over query models only, so long models left 31 of 32 threads idle.
    Chunked search (chunks x all models on a thread pool, E-values rescaled to the chromosome) must give
    the same monomer calls as one whole-chromosome search."""
    import pyhmmer
    import pysam
    from cenfuse.sat.annotate import search_chrom, search_genome
    from cenfuse.sat.genomes import discover_genomes
    fam, wd = synth_families
    with pyhmmer.plan7.HMMFile(wd / "families.hmm") as f:
        hmms = list(f)
    g = [x for x in discover_genomes(synth_complex[0]) if x.name == "Pa_A"][0]
    with pysam.FastaFile(str(g.fasta)) as fa:
        whole = resolve_overlaps(_mon(search_chrom(hmms, "Pa_A01", fa.fetch("Pa_A01"), 4, 1e-5)))
    chunked = search_genome(hmms, g, 4, 1e-5, chunk=40_000, overlap=2_000)["Pa_A01"]
    key = lambda d: set(zip(d.start, d.end, d.strand, d.family))
    a, b = key(whole), key(chunked)
    assert len(a ^ b) <= 0.01 * len(a), (len(a), len(b), len(a ^ b))


def test_drop_multimer_hits_prefers_shortest_unit():
    """Real data: composite/multimer families (CentBd x7 + spacer = 1549 bp; CentO 155 x2.3 = 353 bp) whose long
    HMM hits outscore the monomer hits on the monomer arrays. A hit is dropped when hits of families with a
    >= 1.9x shorter model cover >= 50% of it; a long hit without such cover is kept."""
    rows = [("c", 0, 1549, "+", "long", 900.0, 1e-200, 1, 1549, 1549, 0.9)]
    rows += [("c", 170 + 157 * i, 327 + 157 * i, "+", "mono", 120.0, 1e-30, 1, 157, 157, 0.92) for i in range(7)]
    rows += [("c", 5000, 6549, "+", "long", 900.0, 1e-200, 1, 1549, 1549, 0.9)]        # no monomer cover
    rows += [("c", 5100, 5257, "+", "mono", 120.0, 1e-30, 1, 157, 157, 0.92)]           # 10% cover only
    out = drop_multimer_hits(_mon(rows))
    assert (out.family == "mono").sum() == 8
    assert out[out.family == "long"].start.tolist() == [5000]


def test_phase_shifted_hits_are_not_double_counted():
    """Real data (A. thaliana, 4 phase-shifted CEN180 subfamilies): hits overlapping each of two accepted
    neighbours by <= 50% survived the pairwise clash rule, inflating monomer counts by ~27% and putting
    overlapping monomers into HOR strings. A hit is accepted only if <= 20% of it is covered by already
    accepted (higher-scoring) hits."""
    rows = [("c", 0, 178, "+", "A", 100.0, 1e-30, 1, 178, 178, 0.95),
            ("c", 178, 356, "+", "A", 100.0, 1e-30, 1, 178, 178, 0.95),
            ("c", 89, 267, "+", "B", 95.0, 1e-28, 1, 178, 178, 0.93),        # half-phase shifted: 100% covered
            ("c", 150, 285, "+", "C", 60.0, 1e-15, 20, 154, 178, 0.9),       # bridging partial
            ("c", 350, 528, "+", "B", 90.0, 1e-27, 1, 178, 178, 0.92)]       # 6 bp overlap: kept
    out = resolve_overlaps(_mon(rows))
    assert list(zip(out.start, out.family)) == [(0, "A"), (178, "A"), (350, "B")]


def test_chain_arrays_across_related_families(test_params):
    """Interleaved monomers of related subfamilies form one array (majority family) instead of the minority ones
    being dropped as singletons."""
    fams = ["A", "A", "B", "A", "A", "B", "A"]
    rows = [("c", 155 * i, 155 * i + (165 if f == "B" else 155), "+", f, 90.0, 1e-20, 1, 155, 155, 0.9)
            for i, f in enumerate(fams)]
    mon, arr = chain_arrays(_mon(rows), {"A": 155, "B": 165}, test_params, groups={"A": "A", "B": "A"})
    assert len(arr) == 1 and arr.n_monomers.iloc[0] == 7 and arr.family.iloc[0] == "A" and len(mon) == 7
    mon, arr = chain_arrays(_mon(rows), {"A": 155, "B": 165}, test_params)
    assert (mon.family == "B").sum() == 0                                  # per-family chaining drops them


@pytest.mark.skipif(not have("FasTAN", "FAtoGDB", "ONEview", "blastn", "mafft"), reason="tools missing")
def test_evalue_ref_len_scales_short_sequences(synth_complex, synth_families):
    """Contigs: E-values scale with sequence length, so a 1-Mb contig got a cutoff ~8-13 bits looser than a
    100-Mb chromosome. evalue_ref_len scales every shorter sequence as if it had that length."""
    import pyhmmer
    from cenfuse.sat.annotate import search_genome
    from cenfuse.sat.genomes import discover_genomes
    fam, wd = synth_families
    with pyhmmer.plan7.HMMFile(wd / "families.hmm") as f:
        hmms = list(f)
    g = [x for x in discover_genomes(synth_complex[0]) if x.name == "Pa_A"][0]
    run = lambda n: search_genome(hmms, g, 4, 1e-5, chunk=100_000, overlap=2_000, ref_len=n)["Pa_A01"]
    base, same, big = run(0), run(1_000), run(3_000_000_000)
    key = lambda d: dict(zip(zip(d.start, d.end, d.family), d.evalue))
    assert key(same) == key(base)                                    # shorter than the chromosome: no effect
    b, s = key(base), key(big)
    shared = set(b) & set(s)
    assert shared and all(abs(s[k] / b[k] / 10_000 - 1) < 1e-6 for k in shared if b[k] > 0)
    from cenfuse.sat.annotate import run_significant
    assert run_significant(big, 3_000_000_000, 1e-5, 2.0).all()      # every hit passes alone or in its run
    assert len(big) <= len(base)                                     # a stricter cutoff never adds hits


def _drop_multimer_hits_reference(df, ratio=1.9, min_cov=0.5):
    """The per-hit loop of satellite engine <= 0.6.1, kept as the reference for the vectorised version."""
    import numpy as np
    from cenfuse.sat.util import union_bp
    if df.empty:
        return df
    df = df.sort_values(["start", "end"]).reset_index(drop=True)
    s, e, m = df.start.to_numpy(), df.end.to_numpy(), df.model_len.to_numpy()
    maxlen = int((e - s).max())
    keep = np.ones(len(df), bool)
    for i in np.where(m >= ratio * m.min())[0]:
        lo, hi = np.searchsorted(s, s[i] - maxlen), np.searchsorted(s, e[i])
        iv = [(max(s[j], s[i]), min(e[j], e[i])) for j in range(lo, hi)
              if j != i and m[j] * ratio <= m[i] and e[j] > s[i]]
        if iv and union_bp(iv) >= min_cov * (e[i] - s[i]):
            keep[i] = False
    return df[keep].reset_index(drop=True)


def test_drop_multimer_hits_matches_the_per_hit_reference():
    """Avena C genomes: tens of millions of hits, all of models >= 1.9x the 38-bp shortest, made the per-hit loop
    take hours. The vectorised version must keep exactly the same hits."""
    import numpy as np
    rng = np.random.default_rng(9)
    for trial in range(30):
        n = int(rng.integers(1, 400))
        lens = rng.choice([38, 48, 90, 178, 356, 553, 1655], size=n)
        start = rng.integers(0, 20_000, n)
        span = (lens * rng.uniform(0.6, 1.0, n)).astype(int) + 1
        rows = [("c", int(a), int(a + b), rng.choice(["+", "-"]), f"f{L}", 10.0, 1e-9, 1, L, int(L), 0.9)
                for a, b, L in zip(start, span, lens)]
        df = pd.DataFrame(rows, columns=MON_COLS)
        a, b = drop_multimer_hits(df), _drop_multimer_hits_reference(df)
        assert a.equals(b), trial


def test_run_only_copies_never_drop_a_hit_significant_alone():
    """Copies significant only as a run fill sequence nothing else explains; they never erase a hit significant alone."""
    rows = [("c", 0, 87, "+", "f87", 60.0, 1e-20, 1, 87, 87, 0.9),
            ("c", 0, 40, "+", "f39", 18.0, 1e-2, 1, 39, 39, 0.9),         # two weak copies (significant only as a run)
            ("c", 42, 82, "+", "f39", 18.0, 1e-2, 1, 39, 39, 0.9)]
    df = pd.DataFrame(rows, columns=MON_COLS)
    assert "f87" not in set(drop_multimer_hits(df).family)               # without evalue: all copies count
    assert "f87" in set(drop_multimer_hits(df, evalue=1e-5).family)      # weak copies cannot override it
    strong = df.assign(evalue=[1e-20, 1e-7, 1e-7])
    assert "f87" not in set(drop_multimer_hits(strong, evalue=1e-5).family)


# ---- progress of the genome-wide search (annotate dominates the run time and was silent without -v)

PROGRESS = re.compile(r"annotate: Pa_A: (\d+)% of the work done in [^;]+; about [^;]+ left \(estimate\)")


def _pa_a(synth_complex, synth_families):
    """The models of the synthetic families and the genome Pa_A."""
    import pyhmmer
    from cenfuse.sat.genomes import discover_genomes
    with pyhmmer.plan7.HMMFile(synth_families[1] / "families.hmm") as f:
        hmms = list(f)
    return hmms, [x for x in discover_genomes(synth_complex[0]) if x.name == "Pa_A"][0]


@pytest.mark.skipif(not have("FasTAN", "FAtoGDB", "ONEview", "blastn", "mafft"), reason="tools missing")
def test_search_logs_progress_at_info(caplog, monkeypatch, synth_complex, synth_families):
    """Without -v the log showed 'annotate: start' and then nothing for up to 25 h (Panax). The search now says at INFO
    how much of its work is done, how long that took and, as an estimate, how long is left; it ends with 100%."""
    from cenfuse.sat.annotate import search_genome
    from cenfuse.sat.util import log
    monkeypatch.setattr(log, "propagate", True)
    hmms, g = _pa_a(synth_complex, synth_families)
    with caplog.at_level(logging.INFO, logger="cenfuse-sat"):
        search_genome(hmms, g, 2, 1e-5, chunk=40_000, overlap=2_000)
    msgs = [r.getMessage() for r in caplog.records if r.name == "cenfuse-sat"]
    assert msgs and all(PROGRESS.fullmatch(m) for m in msgs), msgs
    pcts = [int(PROGRESS.fullmatch(m).group(1)) for m in msgs]
    assert pcts == sorted(set(pcts)) and pcts[-1] == 100, pcts


@pytest.mark.skipif(not have("FasTAN", "FAtoGDB", "ONEview", "blastn", "mafft"), reason="tools missing")
def test_search_keeps_its_verbose_lines(caplog, monkeypatch, synth_complex, synth_families):
    """With -v the finer 'i/n searches' lines stay beside the INFO progress; the last one says all n are done."""
    from cenfuse.sat.annotate import search_genome
    from cenfuse.sat.util import log
    monkeypatch.setattr(log, "propagate", True)
    hmms, g = _pa_a(synth_complex, synth_families)
    with caplog.at_level(logging.DEBUG, logger="cenfuse-sat"):
        search_genome(hmms, g, 2, 1e-5, chunk=40_000, overlap=2_000)
    fine = [r.getMessage() for r in caplog.records if r.levelno == logging.DEBUG]
    assert fine and re.fullmatch(r"annotate: Pa_A: ([\d,]+)/\1 searches \(\d+ s\)", fine[-1]), fine
    assert any(PROGRESS.fullmatch(r.getMessage()) for r in caplog.records if r.levelno == logging.INFO)


def _stand_ins(monkeypatch, tmp_path, search_one, chunks=4):
    """Stand-ins for a search of `chunks` chunks x 3 models (no nhmmer): `search_one(model)` replaces the search, and
    what it returns becomes one hit, at the same place in every job and scored by the model's rank. The chromosome is
    random sequence cut into 1-kb chunks. Returns (hmms, genome, keyword arguments of search_genome)."""
    import numpy as np
    import pysam
    import cenfuse.sat.annotate as annotate
    from sat_common import dna
    fa = tmp_path / "g.fa"
    fa.write_text(">c\n" + dna(np.random.default_rng(1), 800 * chunks + 200) + "\n")
    pysam.faidx(str(fa))
    hmms = [SimpleNamespace(name=f"m{i}", M=100, rank=i) for i in range(3)]
    monkeypatch.setattr(annotate, "search_one", lambda h, blk: search_one(h))
    monkeypatch.setattr(annotate, "_rows", lambda hits, mlen, c, off, scale: [
        (c, 100, 150, "+", "fam", float(hits[0].rank), 1e-9, 1, 100, 100, 1.0)])
    return hmms, SimpleNamespace(name="g", fasta=fa, chroms={"c": 800 * chunks + 200}), dict(chunk=1_000, overlap=200)


def test_search_collects_hits_in_submission_order(monkeypatch, tmp_path):
    """Searches finish in any order, but their hits are collected in submission order (chunk by chunk, model by
    model) as ex.map gave them: when two chunks find the same hit, duplicate removal keeps the first one. Here model
    0 is the slowest, so it finishes last; the surviving duplicate must still be the one of (chunk 0, model 0)."""
    from cenfuse.sat.annotate import search_genome

    def search_one(h):
        time.sleep((3 - h.rank) * 0.02)
        return h

    hmms, g, kw = _stand_ins(monkeypatch, tmp_path, search_one)
    hits = search_genome(hmms, g, 4, 1e-5, **kw)["c"]
    assert hits.bitscore.tolist() == [0.0]


def test_a_failed_search_stops_the_run_and_counts_no_work(caplog, monkeypatch, tmp_path):
    """An error in one search must reach the caller without running the searches still queued (ex.map cancelled them;
    otherwise hours of search go on after the error), and no progress line claims work that was never done."""
    from cenfuse.sat.annotate import search_genome
    from cenfuse.sat.util import log
    monkeypatch.setattr(log, "propagate", True)
    count, calls = itertools.count(), []

    def search_one(h):
        calls.append(h)
        if next(count) == 0:                         # the first search to start fails (next() is atomic)
            raise RuntimeError("nhmmer died")
        time.sleep(0.02)
        return h

    hmms, g, kw = _stand_ins(monkeypatch, tmp_path, search_one, chunks=40)             # 120 searches
    with caplog.at_level(logging.INFO, logger="cenfuse-sat"), pytest.raises(RuntimeError, match="nhmmer died"):
        search_genome(hmms, g, 2, 1e-5, **kw)
    assert len(calls) < 60
    assert not [r for r in caplog.records if "100% of the work done" in r.getMessage()]


def test_search_with_nothing_to_search_logs_no_progress(caplog, monkeypatch, tmp_path):
    """A run with zero families (no models) or no sequences has no work: the search says nothing, divides by nothing."""
    from cenfuse.sat.annotate import search_genome
    from cenfuse.sat.util import log
    monkeypatch.setattr(log, "propagate", True)
    hmms, g, kw = _stand_ins(monkeypatch, tmp_path, lambda h: h)
    with caplog.at_level(logging.INFO, logger="cenfuse-sat"):
        no_models = search_genome([], g, 2, 1e-5, **kw)
        no_sequences = search_genome(hmms, SimpleNamespace(name="g", fasta=g.fasta, chroms={}), 2, 1e-5, **kw)
    assert no_models["c"].empty and no_sequences == {}
    assert not [r for r in caplog.records if r.name == "cenfuse-sat"]
