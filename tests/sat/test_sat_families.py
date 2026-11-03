import numpy as np
import pytest

from cenfuse.sat.families import (RRNA_CM, absorb_multimers, base_unit, cluster_seeds, dedup_families, frame_rotation,
                             msa_consensus, chain_groups, pad_hit, rdna_tags, seed_match, trim_overhang)
from sat_common import have
from sat_synth import mutate, rseq


def test_seed_match_rotation_same_period():
    assert seed_match(150, 150, [(97.0, 1, 150, 40, 189)], 0.9, 0.8)


def test_seed_match_rejects_low_identity_same_period():
    assert not seed_match(150, 150, [(88.0, 1, 150, 40, 189)], 0.9, 0.8)


def test_seed_match_multimer_needs_tiling_and_integer_ratio():
    tiled = [(85.0, 1, 150, 1 + 150 * i, 150 * (i + 1)) for i in range(6)]   # 6 copies on a 900-bp dimer
    assert seed_match(150, 450, tiled, 0.9, 0.8)
    assert not seed_match(150, 450, tiled[:1], 0.9, 0.8)      # covers 1/3 of the circle
    assert not seed_match(150, 400, tiled, 0.9, 0.8)          # 400/150 is not an integer ratio


def test_pad_hit_extends_to_full_model():
    # + strand: hit covers model 3..148 of 150 at target 101..246 -> pad 2 left, 2 right
    assert pad_hit(101, 246, "+", 3, 148, 150, 1000) == (99, 248)
    # - strand: model 1 sits at the high coordinate; model 3..148 at 101..246 -> 2 right, 2 left
    assert pad_hit(101, 246, "-", 3, 148, 150, 1000) == (99, 248)
    # clipped to the sequence
    assert pad_hit(1, 146, "+", 5, 150, 150, 1000) == (1, 146)


def test_msa_consensus_drops_gappy_columns():
    assert msa_consensus(["AC-GT", "ACTGT", "AC-GA"]) == "ACGT"


@pytest.mark.skipif(not have("blastn", "makeblastdb"), reason="blast missing")
def test_cluster_seeds_rotation_multimer_unrelated(tmp_path):
    rng = np.random.default_rng(5)
    a, b = rseq(rng, 150), rseq(rng, 150)
    seqs = {"A": a,
            "Arot": mutate(rng, a[40:] + a[:40], 0.02),
            "A3": "".join(mutate(rng, a, 0.05) for _ in range(3)),
            "B": b,
            "Asub": mutate(rng, a, 0.14)}
    assign = cluster_seeds(seqs, ["A", "Arot", "A3", "B", "Asub"], tmp_path, 2, 0.90, 0.80)
    assert assign["Arot"] == "A" and assign["A3"] == "A"
    assert assign["B"] == "B" and assign["Asub"] == "Asub"


def test_cluster_seeds_logs_its_rounds_at_info(tmp_path, caplog, monkeypatch):
    """Acanthus: 76 min of seed clustering with nothing in the default log. Each round says how many seeds are left
    (BLAST stubbed out: no seed is related to another, so every one founds a cluster; batches of 2 give 3 rounds)."""
    import logging
    from cenfuse.sat import families
    from cenfuse.sat.util import log
    monkeypatch.setattr(log, "propagate", True)
    monkeypatch.setattr(families, "BLAST_BATCH", 2)
    monkeypatch.setattr(families, "blast_pairs", lambda *a, **k: {})
    seqs = {f"s{i}": "ACGT" * 40 for i in range(5)}
    with caplog.at_level(logging.INFO, logger="cenfuse-sat"):
        assign = cluster_seeds(seqs, list(seqs), tmp_path, 1, 0.90, 0.80)
    assert assign == {s: s for s in seqs}
    assert [r.getMessage() for r in caplog.records if r.name == "cenfuse-sat"] == [
        "families: round 1: 2 new representatives, 3 seeds left",
        "families: round 2: 2 new representatives, 1 seeds left",
        "families: round 3: 1 new representatives, 0 seeds left"]


@pytest.mark.skipif(not have("FasTAN", "FAtoGDB", "ONEview", "blastn", "mafft"), reason="tools missing")
def test_build_families_on_synthetic(synth_families):
    fam, workdir = synth_families
    near = lambda n: fam[(fam.consensus_len - n).abs() <= 3]
    assert len(near(150)) == 1 and near(150).hmm_mode.iloc[0] == "profile"   # CEN + its HOR multimer
    assert len(near(90)) == 1                                                  # subtelomeric decoy
    assert len(near(600)) == 0                                                 # multimer merged into CEN
    assert (workdir / "families.hmm").stat().st_size > 0


@pytest.mark.skipif(not have("blastn", "makeblastdb"), reason="blast missing")
def test_absorb_multimers_folds_divergent_hor_multimers(tmp_path):
    """TRF often reports HOR-rich arrays at the HOR period; when some monomers of the HOR copy are too
    divergent for the 80% multimer rule, the multimer becomes its own family (real data: CEN180 x5/x7).
    Absorption at superfamily identity folds it back into its monomer family."""
    rng = np.random.default_rng(8)
    a, b = rseq(rng, 180), rseq(rng, 180)
    cons = {"A": a,
            "A5": "".join(mutate(rng, a, 0.12 if i % 2 else 0.25) for i in range(5)),   # some copies ~75%
            "B": b,
            "B2": "".join(mutate(rng, b, 0.22) for _ in range(2)),
            "C": rseq(rng, 900)}                                                     # unrelated 900 bp
    parent = absorb_multimers(cons, tmp_path, 2, 0.70)
    assert parent == {"A5": "A", "B2": "B"}


@pytest.mark.skipif(not have("blastn", "makeblastdb"), reason="blast missing")
def test_monomer_must_not_join_multimer_rep_below_family_identity(tmp_path):
    """Review C1: TRF often calls the largest array at a multimer period, so the first representative is a
    dimer; a monomer of a related but distinct subfamily (87%, CentBd/CentBs-like) must not join it at the
    80% multimer threshold. Result must not depend on seed order."""
    rng = np.random.default_rng(12)
    a = rseq(rng, 160)
    seqs = {"A2": a + mutate(rng, a, 0.03), "B1": mutate(rng, a, 0.13), "A1": mutate(rng, a, 0.02)}
    for order in (["A2", "B1", "A1"], ["A1", "B1", "A2"]):
        assign = cluster_seeds(seqs, order, tmp_path / "_".join(order), 2, 0.90, 0.80)
        assert assign["B1"] == "B1", order
        assert assign["A1"] == assign["A2"]


@pytest.mark.skipif(not have("cmsearch", "cmemit"), reason="infernal missing")
def test_rdna_tags_flags_5s_rdna_units(tmp_path):
    """Review C2: pericentromeric 5S rDNA arrays look like satellites; tag them from the Rfam 5S model."""
    import subprocess
    out = subprocess.run(["cmemit", "-c", "--seed", "1", str(RRNA_CM)], capture_output=True, text=True).stdout
    five_s = "".join(out.split(">5_8S_rRNA")[0].splitlines()[1:]).upper().replace("U", "T")
    rng = np.random.default_rng(3)
    unit = rseq(rng, 200) + five_s + rseq(rng, 180)                   # ~500-bp 5S rDNA repeat unit
    tags = rdna_tags({"f5s": unit, "sat": rseq(rng, 178)}, tmp_path, 2)
    assert tags == {"f5s": "5S_rRNA", "sat": "."}


@pytest.mark.skipif(not have("blastn", "makeblastdb"), reason="blast missing")
def test_frame_rotation_aligns_related_family_to_reference(tmp_path):
    """Real data: related subfamilies (98-100% identical) got monomer frames 110/178 and 113/157 bp apart from
    independent lexicographic rotations, so their hits half-overlapped and 10-18% of array bp stayed unannotated.
    A related family is rotated (and oriented) into the frame of the most similar earlier family."""
    from cenfuse.sat.util import canonical_rotation, revcomp
    rng = np.random.default_rng(21)
    a = rseq(rng, 180)
    b = mutate(rng, a[70:] + a[:70], 0.05)
    for i, q in enumerate((b, revcomp(b))):
        c = frame_rotation(q, {"A": a}, tmp_path / f"q{i}", 1, 0.70)
        assert sum(x == y for x, y in zip(c[:40], a[:40])) >= 32, i
    u = rseq(rng, 180)
    assert frame_rotation(u, {"A": a}, tmp_path / "u", 1, 0.70) == canonical_rotation(u)
    assert frame_rotation(u, {}, tmp_path / "e", 1, 0.70) == canonical_rotation(u)


@pytest.mark.skipif(not have("blastn", "makeblastdb"), reason="blast missing")
def test_dedup_families_merges_near_identical(tmp_path):
    rng = np.random.default_rng(22)
    x = rseq(rng, 178)
    cons = {"f1": x, "f2": mutate(rng, x[50:] + x[:50], 0.02), "f3": mutate(rng, x, 0.13), "f4": rseq(rng, 178)}
    assert dedup_families(cons, ["f1", "f2", "f3", "f4"], tmp_path, 1, 0.95) == {"f2": "f1"}
    # real data: CentBr1 vs CentBr2 consensus ~95% identical must stay separate at the 0.97 default
    from cenfuse.sat.params import Params
    pair = {"b1": x, "b2": mutate(np.random.default_rng(23), x, 0.05)}
    assert dedup_families(pair, ["b1", "b2"], tmp_path / "p", 1, Params().dedup_ident) == {}


@pytest.mark.skipif(not have("blastn", "makeblastdb"), reason="blast missing")
def test_trim_overhang_cuts_terminal_redundancy(tmp_path):
    """Real data (Msookensis): a 796-bp family = Cent728 unit + the first 68 bp of the next unit (terminally
    redundant) outscored the real 730-bp model on every Cent728 array once both shared a frame. Such a consensus
    is trimmed to its true period; a normal monomer, an internal motif repeat and a short chance match are not."""
    rng = np.random.default_rng(31)
    u = rseq(rng, 730)
    assert trim_overhang(u + mutate(rng, u[:68], 0.02), tmp_path / "a", 1) == u
    assert trim_overhang(u, tmp_path / "b", 1) == u
    w = u[:500] + u[:40] + u[540:]                                     # internal 40-bp duplication at 500
    assert trim_overhang(w, tmp_path / "d", 1) == w
    v = rseq(rng, 300)
    assert trim_overhang(v + v[:12], tmp_path / "c", 1) == v + v[:12]


@pytest.mark.skipif(not have("blastn", "makeblastdb"), reason="blast missing")
def test_chain_groups_join_close_subfamilies_only(tmp_path):
    """Real data (Oryza): CentO-165 monomers interleaved among CentO-155 were dropped as per-family singletons.
    Arrays are chained per chain group: subfamilies >= 90% identical (length within 15%); CentBd/CentBs-like
    87% pairs and unrelated families stay apart."""
    rng = np.random.default_rng(41)
    a = rseq(rng, 155)
    b = mutate(rng, a[:80] + "ATTTATAGGC" + a[80:], 0.01)          # 165-bp insertion variant
    c = mutate(rng, a, 0.13)
    d = rseq(rng, 155)
    g = chain_groups({"a": a, "b": b, "c": c, "d": d}, ["a", "b", "c", "d"], tmp_path, 1, 0.90)
    assert g["a"] == g["b"] == "a" and g["c"] == "c" and g["d"] == "d"


@pytest.mark.skipif(not have("cmsearch", "cmemit"), reason="infernal missing")
def test_rdna_neighbourhood_flags_45s_spacer_repeats(tmp_path):
    """Real data (P. infirma chr1A NOR): the 45S rDNA intergenic-spacer subrepeat (Poa_sat18_125) forms satellite-like
    arrays without any rRNA gene in its consensus, so rdna_tags missed it and it was called centromeric. Arrays whose
    flanks carry 45S rRNA genes (18S/5.8S/25S) are rDNA-associated; a satellite with plain flanks is not."""
    import subprocess
    import pandas as pd
    import pysam
    from cenfuse.sat.families import rdna_neighbourhood
    out = subprocess.run(["cmemit", "-c", "--seed", "1", str(RRNA_CM)], capture_output=True, text=True).stdout
    recs = {r.splitlines()[0].split()[0]: "".join(r.splitlines()[1:]) for r in out.split(">")[1:]}
    ssu = next(s for n, s in recs.items() if n.startswith("SSU")).upper().replace("U", "T")
    rng = np.random.default_rng(5)
    igs, sat = rseq(rng, 125), rseq(rng, 146)
    seq, rows = [], []
    pos = 0
    for i in range(6):                                    # 6 rDNA units: 18S gene + spacer with an IGS array
        for part, fam in ((rseq(rng, 2000), None), (ssu, None), (rseq(rng, 3000), None), (igs * 20, "igs")):
            if fam:
                rows.append((len(rows), "c1", pos, pos + len(part), "+", fam, 20, 0.9, 125))
            seq.append(part)
            pos += len(part)
    for i in range(6):                                    # satellite arrays with random flanks
        seq.append(rseq(rng, 40_000))
        pos += 40_000
        rows.append((len(rows), "c1", pos, pos + 146 * 30, "+", "sat", 30, 0.9, 146))
        seq.append(sat * 30)
        pos += 146 * 30
    seq.append(rseq(rng, 40_000))
    fa = tmp_path / "g.fa"
    fa.write_text(">c1\n" + "".join(seq) + "\n")
    pysam.faidx(str(fa))
    from cenfuse.sat.annotate import ARR_COLS
    arrays = pd.DataFrame(rows, columns=ARR_COLS)
    res = rdna_neighbourhood(fa, {"c1": len("".join(seq))}, arrays, ["igs", "sat"], tmp_path / "w", 2,
                             flank=15_000, n_arrays=5)
    assert res["igs"] == (5, 5) and res["sat"] == (5, 0)


def _seeds(rows):
    import pandas as pd
    return pd.DataFrame(rows, columns=["seed_id", "genome", "chrom", "start", "end", "period", "consensus"])


def test_base_unit_skips_consensus_shorter_than_its_period():
    """TRF can drop alignment columns at indels: A. arenosa pAa seed with period 179 but a 170-bp consensus. As the
    'shortest member' it founded a family whose 170-bp model missed 9 bp of every 179-bp monomer."""
    rng = np.random.default_rng(3)
    unit = rseq(rng, 179)
    mem = _seeds([("s1", "G", "c1", 0, 50_000, 179, unit[:170]),           # largest span, truncated consensus
                  ("s2", "G", "c1", 60_000, 90_000, 179, unit),            # full period
                  ("s3", "G", "c2", 0, 40_000, 537, unit * 3),             # trimer
                  ("s4", "G", "c2", 50_000, 60_000, 168, unit[:168])])     # small deletion-variant array
    base, seq = base_unit(mem, {})
    assert base.seed_id == "s2" and seq == unit


def test_base_unit_takes_a_genome_period_when_no_consensus_spans_it(tmp_path):
    import pysam
    rng = np.random.default_rng(4)
    unit = rseq(rng, 179)
    chrom = rseq(rng, 1000) + unit * 20 + rseq(rng, 1000)
    fa = tmp_path / "G.fa"
    fa.write_text(f">c1\n{chrom}\n")
    pysam.faidx(str(fa))
    mem = _seeds([("s1", "G", "c1", 1000, 1000 + 179 * 20, 179, unit[:170])])
    base, seq = base_unit(mem, {"G": fa})
    assert base.seed_id == "s1" and seq == unit                       # one period of the array, from the genome


def test_base_unit_keeps_the_monomer_when_only_multimers_span_their_period(tmp_path):
    """A. arenosa fam014: the only monomer-period seed had a truncated consensus; its full-period cluster mates were
    1220-bp multimers. The model must still be founded on one monomer period, taken from the genome."""
    import pysam
    rng = np.random.default_rng(5)
    unit = rseq(rng, 179)
    chrom = rseq(rng, 1000) + unit * 30 + rseq(rng, 1000)
    fa = tmp_path / "G.fa"
    fa.write_text(f">c1\n{chrom}\n")
    pysam.faidx(str(fa))
    mem = _seeds([("s1", "G", "c1", 1000, 1000 + 179 * 30, 1253, unit * 7),            # full-period 7-mer
                  ("s2", "G", "c1", 1000 + 179 * 5, 1000 + 179 * 25, 179, unit[:170])])  # truncated monomer
    base, seq = base_unit(mem, {"G": fa})
    assert base.seed_id == "s2" and seq == unit


@pytest.mark.skipif(not have("blastn", "mafft"), reason="tools missing")
def test_failed_refinement_keeps_the_pass1_profile(tmp_path, test_params, monkeypatch):
    """Arachis fam018: pass 1 sampled enough monomers, its consensus collapsed (20 bp) and pass 2 found none; the
    family then got a useless single-sequence model. A failed refinement now keeps the pass-1 profile."""
    import random
    import cenfuse.sat.families as F
    rng = np.random.default_rng(3)
    unit = rseq(rng, 120)
    copies = [mutate(rng, unit, 0.05) for _ in range(30)]
    calls = iter([copies, []])
    monkeypatch.setattr(F, "sample_monomers", lambda *a, **k: next(calls))
    hmm, cons, mode, n = F.build_hmm("f", unit, [], {}, tmp_path, 2, test_params, random.Random(1))
    assert mode == "profile" and n == 30 and abs(hmm.M - 120) <= 3


def test_microsatellite_seeds_do_not_found_families(test_params):
    """Arachis: TRF reported the telomere-like 10-bp repeat AGAAAACCCT at period 40 and 30; it founded 'satellite'
    families that were called centromeric. A consensus repeating itself below min_period is a microsatellite."""
    import pandas as pd
    from cenfuse.sat.families import internal_period, select_seeds
    rng = np.random.default_rng(5)
    micro = "AGAAAACCCT" * 4
    real = rseq(rng, 40)
    assert internal_period(micro, 29) == 10 and internal_period(real, 29) is None
    assert internal_period(rseq(rng, 178), 29) is None
    trf = pd.DataFrame({"chrom": ["c", "c"], "start": [0, 10_000], "end": [5_000, 15_000], "period": [40, 40],
                        "copies": [125.0, 125.0], "consensus_size": [40, 40], "pct_match": [95, 95], "pct_indel": [0, 0],
                        "score": [5000, 5000], "entropy": [1.9, 1.9], "consensus": [micro, real]})
    s = select_seeds({"g": trf}, test_params)
    assert s.consensus.tolist() == [real]


@pytest.mark.skipif(not have("mafft"), reason="mafft missing")
def test_mafft_alignment_does_not_depend_on_threads(tmp_path):
    """>= 200 monomers: --auto refines iteratively (FFT-NS-i), whose multithreaded refinement varied between runs; the
    alignment must not change with the run or the thread count."""
    from cenfuse.sat.families import mafft
    rng = np.random.default_rng(3)
    unit = rseq(rng, 300)
    seqs = []
    for _ in range(220):
        s = mutate(rng, unit, 0.2)
        cut = int(rng.integers(0, 280))
        seqs.append(s[:cut] + s[cut + int(rng.integers(0, 20)):])          # a deletion per copy
    a = mafft(seqs, tmp_path / "t1", 1)
    assert all(mafft(seqs, tmp_path / f"t16_{i}", 16) == a for i in range(2))


def test_build_families_founds_families_on_arrayed_clusters_only(synth_complex, synth_seeds, test_params, tmp_path,
                                                                 monkeypatch):
    """build_families carries forward only the clusters arrayed() returns: with none, no family."""
    from cenfuse.sat import families
    from cenfuse.sat.genomes import discover_genomes
    monkeypatch.setattr(families, "arrayed", lambda seeds, p: set())
    fam = families.build_families(synth_seeds, discover_genomes(synth_complex[0]), tmp_path, 2, test_params)
    assert len(fam) == 0


def test_only_clusters_with_an_array_seed_found_families():
    """>= min_array_monomers copies (LEN_TOL short allowed: 2.85 of 3) in one seed at least; two-copy duplications
    alone do not found a satellite family (T2T sheep: four such families, all dispersed, took 93% of the search)."""
    import pandas as pd
    from cenfuse.sat.families import arrayed
    from cenfuse.sat.params import Params
    seeds = pd.DataFrame({"rep": ["r1", "r1", "r2", "r2", "r3"], "copies": [2.0, 2.7, 2.0, 2.86, 40.0]})
    assert arrayed(seeds, Params()) == {"r2", "r3"}
    assert arrayed(seeds, Params(min_array_monomers=2)) == {"r1", "r2", "r3"}


def test_needs_rebuild():
    """A model built before an earlier family was final is the serial one unless that family's consensus length is
    within RELATED_LEN_TOL of one of the build's frame_rotation queries (the only references it would have consulted);
    the boundary is inclusive, as in frame_rotation."""
    from cenfuse.sat.families import needs_rebuild
    assert needs_rebuild([170, 340], [160])
    assert not needs_rebuild([170, 340], [500])
    assert not needs_rebuild([170], [])
    assert needs_rebuild([100], [115])
    assert not needs_rebuild([100], [116])


def _fake_build_hmm(calls: list, slow: str = "fam001", fail: str = None):
    """Stand-in for build_hmm with its signature: records (family, references it was given), sleeps 0.3 s for `slow`,
    raises for `fail`, queries its unit length as frame_rotation queries a consensus, and returns a result that
    depends on the references within RELATED_LEN_TOL of that query (what a frame rotation can depend on) and on its
    random generator; one note per build. It records no frame rotation, so a flagged build is rebuilt."""
    import time
    from cenfuse.sat.params import RELATED_LEN_TOL
    from cenfuse.sat.util import SatelliteError

    def fake(fam_id, rep_seq, intervals, fastas, workdir, threads, p, rng, refs=None, queries=None, notes=None,
             frames=None):
        calls.append((fam_id, list(refs)))
        if fam_id == slow:
            time.sleep(0.3)
        if fam_id == fail:
            raise SatelliteError(f"{fam_id}: build failed")
        queries.append(len(rep_seq))
        frame = [f for f, c in refs.items() if abs(len(c) / len(rep_seq) - 1) <= RELATED_LEN_TOL]
        notes.append(f"families: {fam_id}: frame of {frame}")
        return f"hmm {fam_id} {frame}", rep_seq, "profile", rng.randint(0, 10 ** 9)
    return fake


def test_build_models_rebuilds_when_a_late_reference_matters(monkeypatch, tmp_path, caplog):
    """At 3 threads all three families start before any is final. fam002 (175 bp) started without fam001 (170 bp,
    within 15%), so it is rebuilt once with fam001 as its reference; fam003 (900 bp) started without both, but no
    reference that short could have changed it, so it is kept. Results, and the notes logged in family order at
    acceptance, equal those of one thread (the serial order)."""
    import logging
    import random
    from cenfuse.sat import families
    from cenfuse.sat.util import log
    monkeypatch.setattr(log, "propagate", True)
    monkeypatch.setattr(families, "THREADS_PER_BUILD", 1)          # 3 threads = 3 concurrent builds
    calls = []
    monkeypatch.setattr(families, "build_hmm", _fake_build_hmm(calls))
    jobs = [("fam001", "A" * 170, []), ("fam002", "C" * 175, []), ("fam003", "G" * 900, [])]
    with caplog.at_level(logging.INFO, logger="cenfuse-sat"):
        par = families.build_models(jobs, {}, tmp_path / "t3", 3, None)
    notes = [r.getMessage() for r in caplog.records if r.name == "cenfuse-sat" and r.levelno == logging.WARNING]
    assert [c for c in calls if c[0] == "fam002"] == [("fam002", []), ("fam002", ["fam001"])]
    assert [c for c in calls if c[0] == "fam003"] == [("fam003", [])]
    calls.clear()
    caplog.clear()
    with caplog.at_level(logging.INFO, logger="cenfuse-sat"):
        ser = families.build_models(jobs, {}, tmp_path / "t1", 1, None)
    assert calls == [("fam001", []), ("fam002", ["fam001"]), ("fam003", ["fam001", "fam002"])]
    assert par == ser
    assert [r["fid"] for r in par] == ["fam001", "fam002", "fam003"]
    assert [r["n"] for r in par] == [random.Random(i).randint(0, 10 ** 9) for i in (1, 2, 3)]   # as the serial loop
    assert notes == [r.getMessage() for r in caplog.records if r.name == "cenfuse-sat"] == [
        "families: fam001: frame of []", "families: fam002: frame of ['fam001']", "families: fam003: frame of []"]


def test_build_models_never_references_an_unaccepted_family(monkeypatch, tmp_path):
    """At 2 threads fam003 starts when fam002's first build ends, while fam001 (slow) is not final: fam002 is finished
    but not accepted (it must be rebuilt with fam001), so it is not a reference of fam003. Only accepted families are
    final; a finished one may still be rebuilt."""
    from cenfuse.sat import families
    monkeypatch.setattr(families, "THREADS_PER_BUILD", 1)          # 2 threads = 2 concurrent builds
    calls = []
    monkeypatch.setattr(families, "build_hmm", _fake_build_hmm(calls))
    jobs = [("fam001", "A" * 170, []), ("fam002", "C" * 175, []), ("fam003", "G" * 900, [])]
    families.build_models(jobs, {}, tmp_path, 2, None)
    fam003 = [refs for f, refs in calls if f == "fam003"]
    assert len(fam003) == 1 and "fam002" not in fam003[0]
    assert [refs for f, refs in calls if f == "fam002"] == [[], ["fam001"]]


def test_build_models_surfaces_errors(monkeypatch, tmp_path):
    """A family whose build raises stops the stage with that error within 10 s, at one thread or many, with jobs still
    waiting for a worker, instead of waiting forever for that family to be accepted. build_models runs in a thread so
    that a hang fails the test instead of stalling the suite."""
    import threading
    from cenfuse.sat import families
    from cenfuse.sat.util import SatelliteError
    monkeypatch.setattr(families, "THREADS_PER_BUILD", 1)          # 2 and 8 threads: builds run concurrently
    monkeypatch.setattr(families, "build_hmm", _fake_build_hmm([], fail="fam002"))
    jobs = [(f"fam{i:03d}", "A" * (170 + 50 * i), []) for i in range(1, 6)]
    for threads in (1, 2, 8):
        raised = []

        def call():
            try:
                families.build_models(jobs, {}, tmp_path / f"t{threads}", threads, None)
            except Exception as e:                           # noqa: BLE001 -- the test inspects it
                raised.append(e)

        t = threading.Thread(target=call, daemon=True)
        t.start()
        t.join(10)
        assert not t.is_alive(), f"no error after 10 s at {threads} threads"
        assert len(raised) == 1 and isinstance(raised[0], SatelliteError), (threads, raised)
        assert "fam002: build failed" in str(raised[0])


@pytest.mark.skipif(not have("FasTAN", "FAtoGDB", "ONEview", "blastn", "mafft"), reason="tools missing")
def test_families_identical_at_1_and_4_threads(synth_seeds, synth_complex, test_params, tmp_path, caplog,
                                               monkeypatch):
    """Concurrent family models are the serial ones: families.tsv, families.fa and the seed table byte for byte,
    families.hmm apart from its DATE and COM lines (build time, command), and the stage's log lines in the same order.
    The synthetic fam003 (98 bp) is within 15% of fam002 (90 bp); at 4 threads (one per build) it starts before fam002
    is final, is flagged and its frame rotation replayed with fam002: unrelated, the rotation is unchanged and the
    build kept. 16 threads on 2 families (THREADS_PER_BUILD as shipped): more threads than families, two builds of 8
    threads each."""
    import gzip
    import logging
    from dataclasses import replace
    from cenfuse.sat import families
    from cenfuse.sat.genomes import discover_genomes
    from cenfuse.sat.util import log
    monkeypatch.setattr(log, "propagate", True)
    shipped = families.THREADS_PER_BUILD
    monkeypatch.setattr(families, "THREADS_PER_BUILD", 1)          # 4 threads = 3 concurrent builds
    genomes = discover_genomes(synth_complex[0])
    calls, replays = [], []
    real, real_replay = families.build_hmm, families.rotations_unchanged

    def spy(fam_id, *a, **k):
        calls.append(fam_id)
        return real(fam_id, *a, **k)

    def replay_spy(frames, refs, workdir, ident):
        replays.append((workdir.name, list(refs), real_replay(frames, refs, workdir, ident)))
        return replays[-1][2]

    monkeypatch.setattr(families, "build_hmm", spy)
    monkeypatch.setattr(families, "rotations_unchanged", replay_spy)

    def run(name, threads, p):
        calls.clear()
        replays.clear()
        caplog.clear()
        d = tmp_path / name
        with caplog.at_level(logging.INFO, logger="cenfuse-sat"):
            families.build_families(synth_seeds, genomes, d, threads, p)
        hmm = [x for x in (d / "families.hmm").read_text().splitlines(True) if not x.startswith(("DATE", "COM"))]
        lines = [(r.levelname, r.getMessage()) for r in caplog.records if r.name == "cenfuse-sat"]
        return ((d / "families.tsv").read_bytes(), (d / "families.fa").read_bytes(), "".join(hmm),
                gzip.decompress((d / "seeds.tsv.gz").read_bytes()), lines), list(calls), list(replays)

    one, built1, replayed1 = run("t1", 1, test_params)
    four, built4, replayed4 = run("t4", 4, test_params)
    assert one == four
    assert built1 == ["fam001", "fam002", "fam003"] and replayed1 == []
    assert sorted(built4) == ["fam001", "fam002", "fam003"]                   # flagged, replayed, kept
    assert replayed4 == [("fam003", ["fam001", "fam002"], True)]
    assert ("WARNING", "families: fam003: only 0 monomers sampled; single-sequence HMM") in one[4]
    monkeypatch.setattr(families, "THREADS_PER_BUILD", shipped)
    two = replace(test_params, max_families=2)
    assert run("two1", 1, two)[0] == run("two16", 16, two)[0]


def test_model_workers_split_threads_into_builds_of_threads_per_build():
    """Builds of at least THREADS_PER_BUILD threads (fewer when fewer are given), as many at once as the threads allow
    and there are families; once that number is fixed the threads are shared out among them, so with few families no
    thread idles (MAFFT is still faster at 16 threads than at 4)."""
    from cenfuse.sat.families import THREADS_PER_BUILD, model_workers
    assert THREADS_PER_BUILD == 4
    assert [model_workers(t, 40) for t in (1, 3, 4, 16)] == [(1, 1), (1, 3), (1, 4), (4, 4)]
    assert model_workers(16, 2) == (2, 8) and model_workers(16, 1) == (1, 16)    # more threads than families
    assert model_workers(64, 40) == (16, 4) and model_workers(18, 40) == (4, 4)  # 2 left: fewer than a build's 4
    assert model_workers(7, 40) == (1, 7) and model_workers(10, 2) == (2, 5)


@pytest.mark.skipif(not have("blastn", "makeblastdb"), reason="blast missing")
def test_build_models_replays_the_frame_rotations_of_a_flagged_family(monkeypatch, tmp_path, test_params):
    """fam002 starts before fam001 (within 15% of its length) is final, so it is flagged; its two recorded frame
    rotations are replayed in order with fam001 (real BLAST). The first (a 120-bp piece: no family near that length)
    is unchanged either way. An unrelated sequence of fam001's length keeps its canonical rotation: the build is kept,
    not rebuilt. A relative of fam001 in another frame is rotated into fam001's frame by the second call: the build is
    rebuilt with fam001. Either way the result is the serial one, and the replay writes only to its own folders."""
    from cenfuse.sat import families
    from cenfuse.sat.util import canonical_rotation
    rng = np.random.default_rng(7)
    a = canonical_rotation(rseq(rng, 170))           # fam001's model: its own canonical rotation (no reference)
    r = "TTTTTT" + a[6:]                              # without a's smallest k-mer the relative's canonical rotation
    relative = mutate(rng, r[60:] + r[:60], 0.02)     # starts elsewhere: out of phase with fam001 (as TRF seeds are)
    unrelated = rseq(rng, 175)
    ident = test_params.superfam_ident
    rotated = {q: families.frame_rotation(q, {"fam001": a}, tmp_path / f"pre{i}", 1, ident)
               for i, q in enumerate((unrelated, relative))}
    assert rotated[unrelated] == canonical_rotation(unrelated) and rotated[relative] != canonical_rotation(relative)
    calls = []

    def build_hmm(fam_id, rep_seq, intervals, fastas, workdir, threads, p, rng, refs=None, queries=None, notes=None,
                  frames=None):
        """build_hmm reduced to its two reference-dependent steps, recorded: real frame rotations of a 120-bp piece of
        the unit (as of a trimmed seed) and of the whole unit (as of the pass-1 consensus)."""
        calls.append((fam_id, list(refs)))
        out = []
        for i, q in enumerate((rep_seq[:120], rep_seq)):
            out.append(families.frame_rotation(q, refs, workdir / f"frame{i}", threads, p.superfam_ident,
                                               queries=queries))
            frames.append((q, out[-1]))
        return f"hmm {out}", out[-1], "profile", rng.randint(0, 10 ** 9)

    monkeypatch.setattr(families, "build_hmm", build_hmm)
    monkeypatch.setattr(families, "THREADS_PER_BUILD", 1)          # 2 threads = 2 concurrent builds
    for name, second, rebuilt in (("unrelated", unrelated, False), ("relative", relative, True)):
        jobs = [("fam001", a, []), ("fam002", second, [])]
        calls.clear()
        par = families.build_models(jobs, {}, tmp_path / name, 2, test_params)
        assert [c for c in calls if c[0] == "fam002"] == [("fam002", [])] + [("fam002", ["fam001"])] * rebuilt, name
        wd = tmp_path / name / "fam002"         # BLAST runs only with a candidate: for the whole unit, after fam001
        assert (wd / "replay1").is_dir() and (wd / "frame1").exists() == rebuilt, name     # frame1: the rebuild's
        calls.clear()
        assert par == families.build_models(jobs, {}, tmp_path / f"{name}_serial", 1, test_params), name
