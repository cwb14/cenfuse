"""Bins for small chromosomes (0.10): candidate sizes from the chromosome lengths, and a yeast-like Rabl genome
(0.2-0.6-Mb chromosomes) called at fine bins. Genomes whose chromosomes all have 20 bins at 250 kb keep 0.9's ladder."""
import numpy as np
import pandas as pd
import pytest

from cenfuse.hic import params
from cenfuse.hic.bins import candidate_sizes
from cenfuse.hic.cli import main

from hic_synth import make_genome

STD = [250_000, 500_000, 1_000_000, 2_000_000]


def test_large_chromosomes_keep_the_09_ladder():
    assert candidate_sizes([248e6, 45e6, 5_000_000]) == STD             # shortest exactly 20 x 250 kb
    assert candidate_sizes([]) == STD


def test_yeast_like_genome_goes_down_to_the_shortest_chromosome():
    yeast = [230_218, 813_184, 316_620, 1_531_933, 576_874, 270_161, 1_090_940, 562_643, 439_888, 745_751, 666_816,
             1_078_177, 924_431, 784_333, 1_091_291, 948_066]
    assert candidate_sizes(yeast) == [10_000, 50_000] + STD             # 230 kb / 20 = 11.5 kb -> 10 kb
    assert candidate_sizes(yeast + [85_779]) == candidate_sizes(yeast)   # < 20 x 5 kb: never callable, ignored


def test_intermediate_genomes():
    pfal = [640_851, 947_102, 1_067_971, 3_291_936]
    assert candidate_sizes(pfal)[0] == 10_000                          # 50 kb x 20 = 1 Mb > 640 kb
    pombe = [5_579_133, 4_539_804, 2_452_883]
    assert candidate_sizes(pombe)[0] == 50_000
    bird = [197e6, 150e6, 111e6, 2_500_000] + [5e6] * 30                # microchromosomes, ~1.1 Gb
    assert candidate_sizes(bird)[0] == 50_000                           # 197 Mb / 50 kb = 3,940 bins <= 4,000
    chicken = [200_044_509, 150e6, 111e6, 2_501_204] + [5e6] * 30       # GGswu chr1 needs 4,001 bins at 50 kb
    assert candidate_sizes(chicken) == STD


def test_memory_caps_keep_large_genomes_at_250_kb():
    # a 3-Gb genome with one 1-Mb chromosome: 50 kb would put 5,000 bins on chr1 (> MAX_CHROM_BINS)
    assert candidate_sizes([250e6] * 12 + [1_000_000]) == STD
    # many medium chromosomes: the genome cap (MAX_GENOME_BINS) binds
    many = [4_000_000] * 700                                            # 2.8 Gb in 4-Mb chromosomes
    assert candidate_sizes(many)[0] == 250_000
    assert sum(many) / 50_000 > params.MAX_GENOME_BINS


@pytest.fixture(scope="module")
def yeast(tmp_path_factory):
    """16 chromosomes of 23-90 bins at 10 kb (0.23-0.9 Mb), centromeres inside, clustered (Rabl); no hub, no
    unattached chromosome."""
    d = tmp_path_factory.mktemp("yeast")
    lengths = (23, 81, 32, 90, 58, 27, 85, 56, 44, 75, 67, 88, 62, 78, 70, 60)
    cens = (15, 24, 11, 45, 15, 15, 50, 10, 35, 43, 44, 15, 26, 63, 32, 55)
    truth = make_genome(d, lengths=lengths, cens=cens, n_trans=600_000, n_cis=300_000, unattached=(), hub=None,
                        sat_shift=None, rdna_on=None, size=10_000)
    main(["call", "--pairs", str(truth["pairs"]), "--fai", str(truth["fai"]), "-o", str(d / "o" / "y"), "-t", "4"])
    t = pd.read_csv(d / "o" / "y.centromeres.tsv", sep="\t", keep_default_na=False).set_index("chrom")
    return truth, t


def test_yeast_like_genome_is_called_at_10_kb(yeast):
    truth, t = yeast
    assert set(t.bin_size.astype(int)) == {10_000}
    assert t.callable.astype(str).eq("True").all()                       # the 230-kb chromosome too (23 bins)


def test_yeast_like_centromeres_found(yeast):
    truth, t = yeast
    hit = 0
    for nm, c in zip(truth["names"], truth["cens"]):
        r = t.loc[nm]
        if r.confidence != "none" and r.peak_start != "":
            hit += int(abs(int(r.peak_start) // 10_000 - c) <= 1)
    assert hit >= 14, hit


def test_cool_resolutions_and_the_sequence_floor(monkeypatch):
    from cenfuse.hic import cli
    yeast = {f"c{i}": L for i, L in enumerate([230_218, 813_184, 1_531_933])}
    monkeypatch.setattr(cli, "cool_resolutions", lambda p: [5_000])          # one fine resolution only
    assert cli.bin_sizes("y.cool", yeast) == [5_000]
    assert cli.sequence_floor([5_000], yeast.values()) == 100_000
    big = {f"c{i}": int(6.2e6) if i == 0 else int(150e6) for i in range(20)}  # a 6.2-Mb chromosome in a large genome
    monkeypatch.setattr(cli, "cool_resolutions", lambda p: [10_000 * 2 ** k for k in range(8)])  # 10 kb .. 1.28 Mb
    sizes = cli.bin_sizes("g.mcool", big)
    assert sizes == [320_000, 640_000, 1_280_000]
    assert cli.sequence_floor(sizes, big.values()) == 5_000_000             # 0.9's floor: the 6.2-Mb one stays


def test_standalone_engine_counts_chromosomes_only():
    from cenfuse.hic import cli
    tair = {"Chr1": 30_427_671, "Chr2": 19_698_289, "Chr3": 23_459_830, "Chr4": 18_585_056, "Chr5": 26_975_502,
            "ChrM": 366_924, "ChrC": 154_478}
    assert cli.bin_sizes("p.pairs", tair, chromosome_set=True) == STD       # organelles out: 0.9's ladder
    assert cli.bin_sizes("p.pairs", tair, chromosome_set=False)[0] == 10_000   # cenfuse's set is taken as given
