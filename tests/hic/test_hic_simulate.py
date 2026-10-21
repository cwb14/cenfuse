"""Tests for the mis-mapping null."""
from __future__ import annotations

import argparse
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
import pysam
import pytest

import cenfuse.hic.simulate as sim
from cenfuse.hic import params
from cenfuse.hic.io import InputError
from cenfuse.hic.simulate import best_segments, pair_mates, simulate_pairs

HAVE_MM2 = shutil.which("minimap2") is not None


def tiny_genome(d: Path, seed=3):
    """3 random chromosomes of 2 Mb; a 20-kb segment of chr1 copied into chr2 with 1% substitutions (mis-mapping
    bait: reads with sequencing errors can then map uniquely, MAPQ >= 30, to the wrong copy)."""
    rng = np.random.default_rng(seed)
    seqs = {f"chr{i}": "".join(rng.choice(list("ACGT"), 2_000_000)) for i in (1, 2, 3)}
    seg = np.array(list(seqs["chr1"][500_000:520_000]))
    mut = rng.random(len(seg)) < 0.01
    seg[mut] = [{"A": "C", "C": "G", "G": "T", "T": "A"}[b] for b in seg[mut]]
    seqs["chr2"] = seqs["chr2"][:1_000_000] + "".join(seg) + seqs["chr2"][1_020_000:]
    fa = d / "g.fa"
    with open(fa, "w") as f:
        for n, s in seqs.items():
            f.write(f">{n}\n")
            for i in range(0, len(s), 60):
                f.write(s[i:i + 60] + "\n")
    pysam.faidx(str(fa))
    return fa


def _genome_read(f, chrom, pos, strand, n):
    comp = str.maketrans("ACGT", "TGCA")
    if strand == "+":
        return f.fetch(chrom, pos, pos + n)
    return f.fetch(chrom, pos - n + 1, pos + 1).translate(comp)[::-1]


def test_simulated_reads_match_the_genome(tmp_path):
    fa = tiny_genome(tmp_path)
    fq1, fq2, truth = simulate_pairs(fa, {"chr1": 2_000_000, "chr3": 2_000_000}, 500, np.array([[100, 90]]),
                                     tmp_path / "s", err=0.0)
    t = pd.read_csv(truth, sep="\t")
    f = pysam.FastaFile(str(fa))
    for fq, pc, sc, n in ((fq1, "pos1", "strand1", 100), (fq2, "pos2", "strand2", 90)):
        reads = Path(fq).read_text().split("\n")[1::4]
        assert len(reads) == len(t)
        for (_, r), s in zip(t.iterrows(), reads):
            assert len(s) == n and _genome_read(f, r.chrom, r[pc], r[sc], n) == s


def test_errors_substitute_at_the_requested_rate(tmp_path):
    fa = tiny_genome(tmp_path)
    fq1, _, truth = simulate_pairs(fa, {"chr1": 2_000_000}, 2000, np.array([[100, 100]]), tmp_path / "s", err=0.05)
    t = pd.read_csv(truth, sep="\t")
    reads = Path(fq1).read_text().split("\n")[1::4]
    f = pysam.FastaFile(str(fa))
    mm = tot = 0
    for (_, r), s in zip(t.iterrows(), reads):
        g = _genome_read(f, r.chrom, r.pos1, r.strand1, 100)
        diff = [(x, y) for x, y in zip(g, s) if x != y]
        mm += len(diff)
        tot += 100
    assert 0.03 < mm / tot < 0.07


def test_chunking_keeps_keys_increasing(tmp_path, monkeypatch):
    fa = tiny_genome(tmp_path)
    monkeypatch.setattr(params, "SIM_CHUNK", 50)
    _, _, truth = simulate_pairs(fa, {"chr1": 2_000_000, "chr3": 2_000_000}, 400, np.array([[100, 100]]),
                                 tmp_path / "s")
    k = pd.read_csv(truth, sep="\t").k.to_numpy()
    assert len(k) > 300 and (np.diff(k) == 1).all()


def test_best_segments(tmp_path):
    rows = ["sim.1\t0\t+\tchr1\t1000\t1100\t60", "sim.2\t0\t-\tchr1\t2000\t2100\t60",
            "sim.3\t50\t+\tchr2\t500\t550\t60", "sim.3\t0\t-\tchr1\t900\t950\t40"]
    f = tmp_path / "m.tsv"
    f.write_text("\n".join(rows) + "\n")
    b = best_segments(f).set_index("k")
    assert b.loc[1, "pos"] == 1000 and b.loc[2, "pos"] == 2099          # '-' strand: 5' end = tend - 1
    assert (b.loc[3, "chrom"], b.loc[3, "pos"], b.loc[3, "mq"]) == ("chr1", 949, 40)   # smallest qstart wins
    (tmp_path / "e.tsv").write_text("")
    e = best_segments(tmp_path / "e.tsv")
    assert list(e.columns) == ["k", "chrom", "pos", "st", "mq"] and e.empty


def test_pair_mates_rules():
    m1 = pd.DataFrame(dict(k=[1, 2, 3, 4], chrom=["a", "a", "a", "a"], pos=[100, 100, 100, 100], st=["+"] * 4,
                           mq=[60, 60, 10, 60]))
    m2 = pd.DataFrame(dict(k=[1, 2, 3, 4], chrom=["b", "a", "b", "b"], pos=[5000, 500, 5000, 5000], st=["-"] * 4,
                           mq=[60, 60, 60, 60]))
    p = pair_mates(m1, m2)
    # k=2: cis < 1 kb dropped; k=3: MAPQ 10 dropped; k=4 duplicates k=1 -> dropped
    assert list(p.k) == [1]


@pytest.mark.skipif(not HAVE_MM2, reason="minimap2 not on PATH")
def test_artefacts_come_from_the_duplicated_segment(tmp_path):
    fa = tiny_genome(tmp_path)
    fq1, fq2, truth = simulate_pairs(fa, {"chr1": 2_000_000, "chr2": 2_000_000, "chr3": 2_000_000}, 20_000,
                                     np.array([[120, 120]]), tmp_path / "s")
    idx = sim._build_index(fa, 6_000_000, 4, tmp_path / "g.mmi")         # one index, reused for both mates
    assert idx.stat().st_size > 0 and "-I 1000000000" in Path(f"{idx}.log").read_text()
    tabs = [best_segments(sim._map(idx, fq, 4, Path(str(fq) + ".tsv"))) for fq in (fq1, fq2)]
    truth_t = pd.read_csv(truth, sep="\t").rename(columns={"pos1": "t1", "pos2": "t2"})
    p = pair_mates(*tabs).merge(truth_t, on="k")
    art = p[p.chrom1 != p.chrom2]
    near = lambda c, x: ((c == "chr1") & (abs(x - 510_000) < 11_000)) | ((c == "chr2") & (abs(x - 1_010_000) < 11_000))
    # smoke test: the bait may yield zero artefacts, in which case the 'near the segment' check below is vacuous
    assert len(art) / len(p) < 0.01                       # random sequence: artefacts are rare
    assert len(art) == 0 or (near(art.chrom1, art.pos1) | near(art.chrom2, art.pos2)).mean() >= 0.9


def test_map_failure_raises_input_error(tmp_path):
    fq = tmp_path / "r with space.fq"
    fq.write_text("@sim.0\nACGT\n+\nIIII\n")
    with pytest.raises(InputError, match="minimap2 failed"):
        sim._map(tmp_path / "missing.mmi", fq, 1, tmp_path / "o ut.tsv")
    with pytest.raises(InputError, match="minimap2 failed"):
        sim._build_index(tmp_path / "missing.fa", 1000, 1, tmp_path / "i dx.mmi")


def test_index_holds_the_whole_genome_in_one_part():
    """minimap2 splits its index every -I bases (default 8 Gb) and computes MAPQ per part, so -I must exceed the
    genome length: whole Gb, strictly above it."""
    assert sim.index_bases(1_805_882_605) == 2_000_000_000
    assert sim.index_bases(3_000_000_000) == 4_000_000_000
    assert sim.index_bases(14_500_000_000) == 15_000_000_000 > 8_000_000_000
    assert sim.index_bases(1) == 1_000_000_000


def _centromeres(path, rows):
    cols = ["chrom", "length", "callable", "bin_size", "peak_start", "peak_end", "core_start", "core_end"]
    pd.DataFrame(rows, columns=cols).to_csv(path, sep="\t", index=False)


def test_run_simulate_input_errors(tmp_path, monkeypatch):
    monkeypatch.setattr(sim.shutil, "which", lambda x: "/bin/minimap2")
    fa = tiny_genome(tmp_path)
    o = tmp_path / "x"
    pairs = tmp_path / "p"
    args = argparse.Namespace(fasta=fa, pairs=pairs, o=str(o), fastq=None, n=10, t=1)
    with pytest.raises(InputError, match="p: no such file"):            # every input is checked first
        sim.run_simulate(args)
    pairs.write_text("chr1\t10\t+\tchr1\t5000\t-\n")
    with pytest.raises(InputError, match="r1.fq: no such file"):
        sim.run_simulate(argparse.Namespace(**{**vars(args), "fastq": tmp_path / "r1.fq"}))
    cool = tmp_path / "x.mcool"
    cool.write_bytes(b"")
    with pytest.raises(InputError, match="not .cool"):
        sim.run_simulate(argparse.Namespace(**{**vars(args), "pairs": cool}))
    with pytest.raises(InputError, match="centromeres.tsv"):
        sim.run_simulate(args)
    pd.DataFrame({"chrom": ["chr1"], "bin_size": [1000]}).to_csv(f"{o}.centromeres.tsv", sep="\t", index=False)
    with pytest.raises(InputError, match="missing columns"):
        sim.run_simulate(args)
    _centromeres(f"{o}.centromeres.tsv", [["chr1", 2_000_000, 1, 10_000, None, None, None, None]])
    with pytest.raises(InputError, match="no chromosome"):
        sim.run_simulate(args)


def test_run_simulate_counts_artefacts_per_core(tmp_path, monkeypatch):
    rng = np.random.default_rng(5)
    fa = tmp_path / "g.fa"
    L = 200_000
    with open(fa, "w") as f:
        for c in ("chr1", "chr2", "chr3"):
            f.write(f">{c}\n" + "".join(rng.choice(list("ACGT"), L)) + "\n")
    pysam.faidx(str(fa))
    o = tmp_path / "x"
    # cores (bin 10 kb, window = peak bin +/- 4): chr1 [50,70) kb, win [10,100); chr2 [100,120), win [60,150);
    # chr3 [20,40), win [-20,70)
    _centromeres(f"{o}.centromeres.tsv", [["chr1", L, 20, 10_000, 50_000, 60_000, 50_000, 70_000],
                                         ["chr2", L, 20, 10_000, 100_000, 110_000, 100_000, 120_000],
                                         ["chr3", L, 20, 10_000, 20_000, 30_000, 20_000, 40_000]])
    pairs = [("chr1", 1000 + 700 * i, "chr1", 60_000 + 700 * i) for i in range(100)]         # 100 cis
    pairs += [("chr1", 55_000 + i, "chr2", 80_000 + i) for i in range(40)]                    # chr1 core x chr2 win
    pairs += [("chr2", 105_000 + i, "chr3", 30_000 + i) for i in range(10)]                   # chr2 x chr3 cores
    (tmp_path / "real.pairs").write_text("".join(f"{a}\t{b}\t+\t{c}\t{d}\t-\n" for a, b, c, d in pairs))
    n = 200
    art1 = [(1, 55_000 + 10 * i, 2, 80_000 + 10 * i) for i in range(5)]       # -> chr1 core only
    art2 = [(2, 105_000 + 10 * i, 3, 30_000 + 10 * i) for i in range(3)]      # -> chr2 and chr3 cores

    calls = dict(index=[], map=[], n_truth=None)

    def fake_index(fasta, total, threads, out):
        calls["index"].append((Path(fasta), total))
        Path(out).write_bytes(b"index")
        return Path(out)

    def fake_map(index, fq, threads, out):
        calls["map"].append(Path(index))
        t = pd.read_csv(Path(out).parent / "sim.truth.tsv", sep="\t")
        calls["n_truth"] = len(t)
        mate = 1 if "_1." in str(fq) else 2
        rows = []
        for r in t.itertuples():
            chrom, pos, st = r.chrom, (r.pos1 if mate == 1 else r.pos2), (r.strand1 if mate == 1 else r.strand2)
            rows.append([f"sim.{r.k}", chrom, pos, st])
        key = {}
        # first chr1-truth pairs become art1, first chr2-truth pairs become art2
        k1 = t[t.chrom == "chr1"].k.to_numpy()[:5]
        k2 = t[t.chrom == "chr2"].k.to_numpy()[:3]
        for k, a in zip(k1, art1):
            key[k] = (f"chr{a[0]}", a[1]) if mate == 1 else (f"chr{a[2]}", a[3])
        for k, a in zip(k2, art2):
            key[k] = (f"chr{a[0]}", a[1]) if mate == 1 else (f"chr{a[2]}", a[3])
        with open(out, "w") as f:
            for q, chrom, pos, st in rows:
                k = int(q.split(".")[1])
                if k in key:
                    chrom, pos = key[k]
                ts, te = (pos, pos + 100) if st == "+" else (pos - 99, pos + 1)
                f.write(f"{q}\t0\t{st}\t{chrom}\t{ts}\t{te}\t60\n")
        return out

    monkeypatch.setattr(sim, "_build_index", fake_index)
    monkeypatch.setattr(sim, "_map", fake_map)
    monkeypatch.setattr(sim.shutil, "which", lambda x: "/bin/minimap2")
    sim.run_simulate(argparse.Namespace(fasta=fa, pairs=str(tmp_path / "real.pairs"), o=str(o), fastq=None,
                                        n=n, t=1))
    # one index of the whole FASTA, built once and used for both mates
    assert calls["index"] == [(fa, 3 * L)] and len(calls["map"]) == 2 and len(set(calls["map"])) == 1
    # large intermediates are deleted once the table is written; the minimap2 logs stay
    work = Path(f"{o}.simulate")
    assert not [p.name for p in work.iterdir() if not p.name.endswith(".log")]
    txt = Path(f"{o}.simulate.tsv").read_text().splitlines()
    assert "minimap2 -x sr or bwa mem" in "\n".join(txt[:3]) and "MAPQ >= 30" in "\n".join(txt[:3])
    out = pd.read_csv(f"{o}.simulate.tsv", sep="\t", comment="#").set_index("chrom")
    scale = 100 / (calls["n_truth"] - 8)                 # real cis / correctly mapped simulated cis
    assert f"scale to real {scale:.3g}" in txt[0]
    assert out.loc["chr1", "real_cen_contacts"] == 40 and out.loc["chr2", "real_cen_contacts"] == 10
    assert out.loc["chr3", "real_cen_contacts"] == 10
    for c, k in (("chr1", 5), ("chr2", 3), ("chr3", 3)):
        assert out.loc[c, "expected_artefact_cen_contacts"] == pytest.approx(k * scale, rel=1e-3)


@pytest.mark.parametrize("n", ["0", "-5", "x"])
def test_simulate_needs_at_least_one_pair(tmp_path, n, capsys):
    from cenfuse.hic.cli import main
    with pytest.raises(SystemExit) as e:
        main(["simulate", "--fasta", str(tmp_path / "g.fa"), "--pairs", str(tmp_path / "p"), "-o", str(tmp_path / "x"),
              "-n", n])
    assert e.value.code == 2 and "-n" in capsys.readouterr().err


def test_short_chromosomes_are_skipped(tmp_path, caplog):
    fa = tiny_genome(tmp_path)
    with caplog.at_level("WARNING", logger="cenfuse-hic"):
        _, _, truth = simulate_pairs(fa, {"chr1": 2_000_000, "tiny": 500}, 50, np.array([[100, 100]]), tmp_path / "s")
    assert set(pd.read_csv(truth, sep="\t").chrom) == {"chr1"} and "skipping tiny" in caplog.text
