"""Negative controls (cenfuse.controls, future.md #34): Hi-C re-paired (NEG2) and relocated (NEG), the dispersion check
that catches clumped contacts, GC blocks and LTR positions shuffled, and the command line. Tiny synthetic fixtures."""
import gzip
import logging
from collections import Counter

import numpy as np
import pandas as pd
import pysam
import pytest

from cenfuse import controls as C
from cenfuse import inputs
from cenfuse.sat import signature, te
from cenfuse.util import ENGINE_LOGGERS, CenfuseError

CHROMS = {"chrA": 50_000, "chrB": 30_000, "chrC": 20_000}
COLS = "readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type mapq1 mapq2".split()


def header(chroms=CHROMS, cols=COLS):
    return (["## pairs format v1.0.0", "#sorted: chr1-chr2-pos1-pos2", "#shape: upper triangle",
             "#genome_assembly: test"] + [f"#chromsize: {c} {L}" for c, L in chroms.items()]
            + ["#samheader: @PG\tID:bwa\tPN:bwa", "#columns: " + " ".join(cols)])


def hic_rows(n, seed=3, cis=0.8):
    """Hi-C-like pairs sorted as pairtools sorts them: mostly cis, with contacts decaying with distance."""
    rng = np.random.default_rng(seed)
    names, rows = list(CHROMS), []
    for k in range(n):
        c1 = names[rng.choice(3, p=[0.5, 0.3, 0.2])]
        p1 = int(rng.integers(1, CHROMS[c1] + 1))
        if rng.random() < cis:
            c2, p2 = c1, int(np.clip(p1 + rng.choice([-1, 1]) * rng.exponential(2000), 1, CHROMS[c1]))
        else:
            c2 = names[rng.integers(3)]
            p2 = int(rng.integers(1, CHROMS[c2] + 1))
        rows.append([f"SRR1.{k}", c1, p1, c2, p2, "+-"[rng.integers(2)], "+-"[rng.integers(2)], "UU",
                     str(rng.choice([60, 30, 7])), str(rng.choice([60, 42]))])
    rows.sort(key=lambda r: (r[1], r[3], r[2], r[4]))
    return rows


def write_pairs(path, rows, head=None):
    with (gzip.open(path, "wt") if str(path).endswith(".gz") else open(path, "w")) as f:
        f.write("\n".join(head or header()) + "\n" + "".join("\t".join(map(str, r)) + "\n" for r in rows))
    return path


def read_pairs(path, cols=COLS):
    """(header lines, data as a DataFrame of strings)."""
    with (gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)) as f:
        lines = f.read().splitlines()
    return ([x for x in lines if x.startswith("#")],
            pd.DataFrame([x.split("\t") for x in lines if not x.startswith("#")], columns=cols))


def ends(d):
    """Every read end: (chromosome, position, strand, MAPQ)."""
    return Counter(zip(d.chrom1, d.pos1, d.strand1, d.mapq1)) + Counter(zip(d.chrom2, d.pos2, d.strand2, d.mapq2))


def coverage(d, size=1000):
    """Read ends per bin of `size` bp on each chromosome."""
    out = {c: np.zeros(L // size + 1, int) for c, L in CHROMS.items()}
    for c, p in list(zip(d.chrom1, d.pos1)) + list(zip(d.chrom2, d.pos2)):
        out[c][int(p) // size] += 1
    return out


@pytest.fixture(autouse=True)
def keep_logging():
    """cenfuse.controls.main configures the engine loggers (handlers, level, no propagation): put them back."""
    saved = [(lg, lg.handlers[:], lg.level, lg.propagate) for lg in map(logging.getLogger, ENGINE_LOGGERS)]
    yield
    for lg, h, lv, pr in saved:
        lg.handlers[:], lg.level, lg.propagate = h, lv, pr


@pytest.fixture(scope="module")
def pairs(tmp_path_factory):
    return write_pairs(tmp_path_factory.mktemp("pairs") / "in.pairs.gz", hic_rows(600))


@pytest.fixture(scope="module")
def big(tmp_path_factory):
    """3,000 pairs: enough for the dispersion check to choose 32 genome bins."""
    return write_pairs(tmp_path_factory.mktemp("big") / "big.pairs.gz", hic_rows(3000))


# ---------------------------------------------------------------------------------------- re-paired Hi-C (NEG2)
def test_repair_keeps_every_read_end_and_so_every_bins_coverage(pairs, tmp_path):
    r = C.repair_pairs(pairs, tmp_path / "out.pairs.gz", seed=1)
    _, d0 = read_pairs(pairs)
    _, d1 = read_pairs(tmp_path / "out.pairs.gz")
    assert r["pairs"] == len(d1) == len(d0) == 600 and r["status"] == "pass"
    assert ends(d1) == ends(d0)
    assert all((coverage(d1)[c] == coverage(d0)[c]).all() for c in CHROMS)
    assert list(d1.readID) == [f"r{k}" for k in range(600)] and set(d1.pair_type) == {"UU"}
    assert Counter(zip(d1.chrom1, d1.chrom2)) != Counter(zip(d0.chrom1, d0.chrom2))    # partners did change


def test_repair_and_relocate_are_deterministic_given_a_seed(pairs, tmp_path):
    for f in (C.repair_pairs, C.relocate_pairs):
        a, b, c = (tmp_path / f"{f.__name__}{k}.pairs.gz" for k in "abc")
        f(pairs, a, seed=5), f(pairs, b, seed=5), f(pairs, c, seed=6)
        assert a.read_bytes() == b.read_bytes()                  # gzip without name or time stamp: the same bytes
        assert read_pairs(a)[1].values.tolist() != read_pairs(c)[1].values.tolist()


# ---------------------------------------------------------------------------------------- relocated Hi-C (NEG)
def test_relocate_keeps_chromosome_pairs_and_every_other_column(pairs, tmp_path):
    r = C.relocate_pairs(pairs, tmp_path / "out.pairs.gz", seed=1)
    _, d0 = read_pairs(pairs)
    _, d1 = read_pairs(tmp_path / "out.pairs.gz")
    assert r["status"] == "pass" and r["null"] == "relocate"
    assert list(zip(d1.chrom1, d1.chrom2)) == list(zip(d0.chrom1, d0.chrom2))        # each contact keeps its pair
    keep = ["readID", "strand1", "strand2", "pair_type", "mapq1", "mapq2"]
    assert d1[keep].equals(d0[keep])
    for c, p in (("chrom1", "pos1"), ("chrom2", "pos2")):
        pos, L = d1[p].astype(int).to_numpy(), d1[c].map(CHROMS).to_numpy()
        assert ((pos >= 1) & (pos <= L)).all()
        assert (pos != d0[p].astype(int).to_numpy()).mean() > 0.99


# ---------------------------------------------------------------------------------------- shared: header, formats
@pytest.mark.parametrize("make", [C.repair_pairs, C.relocate_pairs])
def test_header_lines_are_kept_with_sorted_made_none(pairs, tmp_path, make):
    make(pairs, tmp_path / "out.pairs.gz")
    h0, h1 = read_pairs(pairs)[0], read_pairs(tmp_path / "out.pairs.gz")[0]
    assert h1 == [("#sorted: none" if x.startswith("#sorted:") else x) for x in h0]
    assert h1[1] == "#sorted: none" and "#samheader: @PG\tID:bwa\tPN:bwa" in h1


def test_plain_pairs_and_a_minimal_column_set(tmp_path):
    cols = ["readID", "chrom1", "pos1", "chrom2", "pos2"]
    rows = [r[:5] for r in hic_rows(300)]
    src = write_pairs(tmp_path / "in.pairs", rows, header(cols=cols))
    C.repair_pairs(src, tmp_path / "rep.pairs")
    C.relocate_pairs(src, tmp_path / "rel.pairs")
    d0, rep, rel = read_pairs(src, cols)[1], read_pairs(tmp_path / "rep.pairs", cols)[1], \
        read_pairs(tmp_path / "rel.pairs", cols)[1]
    assert Counter(zip(rep.chrom1, rep.pos1)) + Counter(zip(rep.chrom2, rep.pos2)) == \
        Counter(zip(d0.chrom1, d0.pos1)) + Counter(zip(d0.chrom2, d0.pos2))
    assert list(rel.readID) == list(d0.readID) and list(rel.chrom2) == list(d0.chrom2)


def permute_pairs(src, dst, mode, seed, chunk):
    """v2/scripts/permute_pairs.py (2026-10-21), the generator of the benchmark's NEG and NEG2 controls, as it reads and
    draws, with its chunk size (2,000,000) as a parameter."""
    sizes, head, cols = {}, [], None
    with gzip.open(src, "rt") as f:
        for line in f:
            if not line.startswith("#"):
                break
            if line.startswith("#chromsize:"):
                _, c, n = line.split()
                sizes[c] = int(n)
            if line.startswith("#columns:"):
                cols = line.split()[1:]
            head.append("#sorted: none\n" if line.startswith("#sorted:") else line)
    rng = np.random.default_rng(seed)

    def chunks():
        return pd.read_csv(src, sep="\t", header=None, names=cols, skiprows=len(head), dtype=str, chunksize=chunk)

    with gzip.open(dst, "wt") as out:
        out.write("".join(head))
        if mode == "relocate":
            for ch in chunks():
                for c, p in (("chrom1", "pos1"), ("chrom2", "pos2")):
                    lens = ch[c].map(sizes)
                    ch[p] = (np.floor(rng.random(len(ch)) * lens.to_numpy(float)).astype(np.int64) + 1).astype(str)
                ch.to_csv(out, sep="\t", header=False, index=False)
            return
        names = list(sizes)
        code = {c: i for i, c in enumerate(names)}
        have = [c for c in ("strand1", "strand2", "mapq1", "mapq2") if c in cols]
        e = {k: [] for k in ("chrom", "pos", "strand", "mapq")}
        for ch in chunks():
            for cc, pc, sc, qc in C.END_COLS:
                e["chrom"].append(ch[cc].map(code).to_numpy(np.int32))
                e["pos"].append(ch[pc].to_numpy(np.int64))
                e["strand"].append(ch[sc].to_numpy(str) if sc in have else np.full(len(ch), "."))
                e["mapq"].append(ch[qc].to_numpy(str) if qc in have else np.full(len(ch), "."))
        e = {k: np.concatenate(v) for k, v in e.items()}
        n2 = len(e["pos"]) - len(e["pos"]) % 2
        order = rng.permutation(len(e["pos"]))[:n2]
        i, j = order[0::2], order[1::2]
        for s in range(0, len(i), chunk):
            ii, jj = i[s:s + chunk], j[s:s + chunk]
            d = {"readID": [f"r{k}" for k in range(s, s + len(ii))],
                 "chrom1": np.array(names)[e["chrom"][ii]], "pos1": e["pos"][ii],
                 "chrom2": np.array(names)[e["chrom"][jj]], "pos2": e["pos"][jj],
                 "strand1": e["strand"][ii], "strand2": e["strand"][jj], "pair_type": "UU",
                 "mapq1": e["mapq"][ii], "mapq2": e["mapq"][jj]}
            pd.DataFrame({c: d.get(c, ".") for c in cols}).to_csv(out, sep="\t", header=False, index=False)


@pytest.mark.parametrize("mode", ["repair", "relocate"])
def test_a_seed_gives_permute_pairs_bytes_across_chunks(pairs, tmp_path, monkeypatch, mode):
    """The stored NEG/NEG2 controls stay reproducible: the same bytes as permute_pairs.py, here over 7 chunks."""
    monkeypatch.setattr(C, "CHUNK", 97)
    permute_pairs(pairs, tmp_path / "ref.pairs.gz", mode, 1, 97)
    (C.repair_pairs if mode == "repair" else C.relocate_pairs)(pairs, tmp_path / "new.pairs.gz", seed=1)
    with gzip.open(tmp_path / "ref.pairs.gz", "rb") as a, gzip.open(tmp_path / "new.pairs.gz", "rb") as b:
        assert a.read() == b.read()


def bad_input(tmp_path, kind):
    rows, head = hic_rows(50), header()
    if kind == "no chromsize":
        head = [x for x in head if not x.startswith("#chromsize")]
    elif kind == "no columns":
        head = [x for x in head if not x.startswith("#columns")]
    elif kind == "no pos2":
        head, rows = header(cols=COLS[:4]), [r[:4] for r in rows]
    elif kind == "bad chromsize":
        head[4] = "#chromsize: chrA"
    elif kind == "unknown chromosome":
        rows[10][3] = "chrZ"
    elif kind == "position beyond the chromosome":
        rows[10][2] = CHROMS[rows[10][1]] + 1
    elif kind == "short line":
        rows[10] = rows[10][:9]
    elif kind == "long first line":
        rows[0] = rows[0] + ["extra"]
    elif kind == "not an integer":
        rows[10][4] = "12x"
    return write_pairs(tmp_path / "bad.pairs.gz", rows, head)


@pytest.mark.parametrize("kind", ["no chromsize", "no columns", "no pos2", "bad chromsize", "unknown chromosome",
                                  "position beyond the chromosome", "short line", "long first line", "not an integer"])
def test_bad_pairs_stop_with_a_clear_error(tmp_path, kind):
    src = bad_input(tmp_path, kind)
    with pytest.raises(CenfuseError, match="bad.pairs.gz"):
        C.repair_pairs(src, tmp_path / "out.pairs.gz")
    assert not (tmp_path / "out.pairs.gz").exists() and not (tmp_path / "out.pairs.gz.part").exists()
    if kind not in ("position beyond the chromosome", "not an integer"):    # relocate replaces positions unread
        with pytest.raises(CenfuseError):
            C.relocate_pairs(src, tmp_path / "out.pairs.gz")


def test_outputs_that_cannot_be_written_stop(pairs, tmp_path):
    with pytest.raises(CenfuseError, match="no such directory"):
        C.repair_pairs(pairs, tmp_path / "missing" / "out.pairs.gz")
    with pytest.raises(CenfuseError, match="overwrite the input"):
        C.relocate_pairs(pairs, pairs)
    with pytest.raises(CenfuseError, match="no such file"):
        C.repair_pairs(tmp_path / "none.pairs.gz", tmp_path / "out.pairs.gz")


# ---------------------------------------------------------------------------------------- the dispersion check
@pytest.mark.parametrize("make", [C.repair_pairs, C.relocate_pairs])
def test_dispersion_of_random_controls_sits_at_one(big, tmp_path, make):
    """The null measured: 30 seeds all pass, their mean is 1 and their spread no wider than the bound's sd says."""
    r = [make(big, tmp_path / f"s{s}.pairs.gz", seed=s) for s in range(1, 31)]
    x = np.array([v["dispersion"] for v in r])
    assert {v["status"] for v in r} == {"pass"} and {v["bins"] for v in r} == {32}
    assert abs(x.mean() - 1) < 0.08
    assert x.std() < 1.5 * r[0]["sd"] and x.max() < r[0]["bound"] - 5 * r[0]["sd"]


def test_null_moments_match_simulated_controls():
    """The expected count and the variance of every bin pair, as dispersion() computes them, against 4,000 simulated
    controls of a small genome in 8 bins of 400 bp: s1 fills bins 0-1, s2 lies inside bin 2, s3 spans bins 2-7."""
    L = np.array([800, 300, 2100])
    off, G, B = np.concatenate([[0], np.cumsum(L)]), int(L.sum()), 8
    T = {0: 60, 1: 20, 4: 40, 8: 80, 2: 15, 7: 10}                    # oriented sequence pair i * 3 + j -> contacts
    c1 = np.repeat([k // 3 for k in T], list(T.values()))
    c2 = np.repeat([k % 3 for k in T], list(T.values()))
    n, R, rng = len(c1), C._relocate_null(T, L, B, G), np.random.default_rng(4)

    def counts(b1, b2):
        return np.bincount(b1 * B + b2, minlength=B * B).reshape(B, B)

    rel, rep = [], []
    ends = np.concatenate([(off[c1] + rng.integers(0, L[c1])) * B // G, (off[c2] + rng.integers(0, L[c2])) * B // G])
    for _ in range(4000):
        b1 = (off[c1] + np.floor(rng.random(n) * L[c1]).astype(np.int64)) * B // G
        b2 = (off[c2] + np.floor(rng.random(n) * L[c2]).astype(np.int64)) * B // G
        rel.append(C._moments(counts(b1, b2), n, R)[0])
        o = rng.permutation(2 * n)
        rep.append(C._moments(counts(ends[o[0::2]], ends[o[1::2]]), n)[0])
    nulls = ((rel, C._moments(counts(b1, b2), n, R)), (rep, C._moments(counts(ends[:n], ends[n:]), n)))
    for sims, (_, exp, var) in nulls:
        sims = np.array(sims)
        assert (np.abs(sims.mean(axis=0) - exp) <= 4 * np.sqrt(var / 4000) + 1e-6).all()
        v = var > 1                                     # (bin 2 holds s2's contacts with itself whatever the draw)
        assert np.allclose(sims.var(axis=0)[v], var[v], rtol=0.1)
        assert (np.abs(sims.var(axis=0)[~v] - var[~v]) < 0.5).all()


def test_relocated_contacts_of_sequences_inside_one_bin_do_not_dilute_the_index(tmp_path):
    """Cis contacts of a scaffold shorter than a bin land in one bin pair whatever the draw (no variance): they are
    left out, so the index stays at 1 (the expected count as its denominator gave 0.21 on a synthetic genome of
    100,000 scaffolds)."""
    scaffolds = {f"scf{k}": 300 for k in range(200)}
    rng = np.random.default_rng(8)
    rows = hic_rows(3000) + [[f"s{k}", c, int(rng.integers(1, 301)), c, int(rng.integers(1, 301)), "+", "-", "UU",
                              "60", "60"] for k, c in enumerate(rng.choice(list(scaffolds), 2000))]
    src = write_pairs(tmp_path / "in.pairs.gz", rows, header(CHROMS | scaffolds))
    x = [C.relocate_pairs(src, tmp_path / f"r{s}.pairs.gz", seed=s)["dispersion"] for s in range(1, 11)]
    assert abs(np.mean(x) - 1) < 0.1


def clumped(src, dst):
    """A broken re-pairing: the ends sorted along the genome and paired with their neighbors (each bin pairs with
    itself), coverage kept."""
    _, d = read_pairs(src)
    order = {c: i for i, c in enumerate(CHROMS)}
    e = sorted([(order[c], int(p), c, s, q) for c, p, s, q in ends(d).elements()])
    rows = [[f"r{k}", a[2], a[1], b[2], b[1], a[3], b[3], "UU", a[4], b[4]]
            for k, (a, b) in enumerate(zip(e[0::2], e[1::2]))]
    return write_pairs(dst, rows)


def test_dispersion_flags_clumped_and_structured_contacts(big, tmp_path):
    bad = C.dispersion(clumped(big, tmp_path / "clumped.pairs.gz"), "repair")
    assert bad["status"] == "fail" and bad["dispersion"] > 3 * bad["bound"]   # the bin pairs tested are empty
    real = C.dispersion(big, "repair")                   # unrandomized Hi-C: mostly cis, near the diagonal
    assert real["status"] == "fail" and real["dispersion"] > 5
    C.relocate_pairs(big, tmp_path / "rel.pairs.gz")
    assert C.dispersion(tmp_path / "rel.pairs.gz", "relocate")["status"] == "pass"
    assert C.dispersion(tmp_path / "rel.pairs.gz", "repair")["status"] == "fail"    # relocation keeps cis contacts


def test_repair_stops_when_its_pairing_is_not_random(big, tmp_path, monkeypatch):
    """A pairing that keeps the order of the read ends (as a regular shuffle would) is caught, and no output is
    written under the control's name."""
    monkeypatch.setattr(C, "_pairing", lambda m, rng: (np.arange(0, m - 1, 2), np.arange(1, m, 2)))
    with pytest.raises(C.ControlError, match="clumped"):
        C.repair_pairs(big, tmp_path / "out.pairs.gz")
    assert not (tmp_path / "out.pairs.gz").exists() and (tmp_path / "out.pairs.gz.part").exists()


def test_relocate_stops_when_positions_are_not_random(big, tmp_path, monkeypatch):
    class Middle:                                        # every end at its chromosome's middle
        def __init__(self, seed):
            pass

        def random(self, n):
            return np.full(n, 0.5)

    monkeypatch.setattr(C.np.random, "default_rng", Middle)
    with pytest.raises(C.ControlError, match="relocated contacts are clumped"):
        C.relocate_pairs(big, tmp_path / "out.pairs.gz")


def test_too_few_pairs_are_written_untested(tmp_path, caplog, monkeypatch):
    src = write_pairs(tmp_path / "tiny.pairs.gz", hic_rows(20))
    monkeypatch.setattr(logging.getLogger("cenfuse"), "propagate", True)
    with caplog.at_level(logging.WARNING, logger="cenfuse"):
        r = C.repair_pairs(src, tmp_path / "out.pairs.gz")
    assert r["status"] == "untested" and (tmp_path / "out.pairs.gz").exists()
    assert any("too few contacts" in m for m in caplog.messages)


# ---------------------------------------------------------------------------------------- GC blocks
def write_fasta(path, seqs, width=70):
    with open(path, "w") as f:
        for name, s in seqs.items():
            f.write(f">{name}\n" + "".join(s[i:i + width] + "\n" for i in range(0, len(s), width)))
    pysam.faidx(str(path))
    return path


def genome(lengths, seed=2):
    """Random sequence with soft-masked and N stretches, so that composition and case can be checked."""
    rng = np.random.default_rng(seed)
    out = {}
    for name, L in lengths.items():
        s = rng.choice(list("ACGT"), size=L)
        s[L // 3:L // 3 + 40] = np.char.lower(s[L // 3:L // 3 + 40])
        s[L // 2:L // 2 + 25] = "N"
        out[name] = "".join(s)
    return out


def fetch_all(path):
    with pysam.FastaFile(str(path)) as fa:
        return {c: fa.fetch(c) for c in fa.references}


GENOME = {"chr1": 2500, "chr2": 1700, "chr3": 800, "scf1": 50}


def test_gc_shuffle_keeps_lengths_and_composition_and_moves_whole_blocks(tmp_path):
    seqs = genome(GENOME)
    src = write_fasta(tmp_path / "g.fa", seqs)
    chroms = ["chr1", "chr2", "chr3"]
    r = C.shuffle_gc_blocks(src, tmp_path / "s.fa", seed=1, block=100, chroms=chroms)
    out = fetch_all(tmp_path / "s.fa")
    assert r["blocks"] == 50 and {c: len(s) for c, s in out.items()} == GENOME
    assert out["scf1"] == seqs["scf1"]                               # outside the set: copied
    a, b = "".join(seqs[c] for c in chroms), "".join(out[c] for c in chroms)
    assert Counter(a) == Counter(b)                                  # composition, case and N included
    blocks = lambda s: sorted(s[i:i + 100] for i in range(0, len(s), 100))
    assert blocks(a) == blocks(b) and a != b                         # the same blocks in another order


def test_gc_shuffle_with_a_short_last_block_keeps_lengths_and_composition(tmp_path):
    lengths = {"chr1": 2531, "chr2": 1702, "chr3": 811}
    seqs = genome(lengths, seed=4)
    src = write_fasta(tmp_path / "g.fa", seqs)
    C.shuffle_gc_blocks(src, tmp_path / "s.fa", seed=3, block=100, chroms=list(lengths))
    out = fetch_all(tmp_path / "s.fa")
    assert {c: len(s) for c, s in out.items()} == lengths
    assert Counter("".join(seqs.values())) == Counter("".join(out.values()))


def test_gc_shuffle_is_deterministic_and_writes_indexed_bgzip(tmp_path):
    src = write_fasta(tmp_path / "g.fa", genome(GENOME))
    for name in ("a.fa", "b.fa", "c.fa", "d.fa.gz"):
        C.shuffle_gc_blocks(src, tmp_path / name, seed=7 if name != "c.fa" else 8, block=100,
                            chroms=["chr1", "chr2", "chr3"])
    assert (tmp_path / "a.fa").read_bytes() == (tmp_path / "b.fa").read_bytes()
    assert (tmp_path / "a.fa").read_bytes() != (tmp_path / "c.fa").read_bytes()
    assert (tmp_path / "d.fa.gz.fai").exists() and (tmp_path / "d.fa.gz.gzi").exists()
    assert fetch_all(tmp_path / "d.fa.gz") == fetch_all(tmp_path / "a.fa")
    assert all(len(x) <= C.FASTA_WIDTH for x in (tmp_path / "a.fa").read_text().splitlines())


def test_gc_shuffle_takes_the_chromosome_set_a_run_would_use(tmp_path):
    seqs = genome(GENOME)
    src = write_fasta(tmp_path / "g.fa", seqs)
    (tmp_path / "g.chromosomes.txt").write_text("chr1\nchr3\n")
    r = C.shuffle_gc_blocks(src, tmp_path / "s.fa", seed=1, block=100)
    out = fetch_all(tmp_path / "s.fa")
    assert r["source"] == "g.chromosomes.txt" and r["chromosomes"] == 2
    assert out["chr2"] == seqs["chr2"] and out["scf1"] == seqs["scf1"]
    assert Counter(out["chr1"] + out["chr3"]) == Counter(seqs["chr1"] + seqs["chr3"])
    (tmp_path / "g.chromosomes.txt").unlink()                        # no list: the length rule
    r = C.shuffle_gc_blocks(src, tmp_path / "t.fa", seed=1, block=100)
    rule = inputs.chromosome_set({c: len(s) for c, s in seqs.items()})
    out = fetch_all(tmp_path / "t.fa")
    assert r["chromosomes"] == len(rule) and all(out[c] == seqs[c] for c in seqs if c not in rule)


def test_gc_shuffle_needs_an_indexed_fasta(tmp_path):
    (tmp_path / "g.fa").write_text(">chr1\nACGT\n")
    with pytest.raises(CenfuseError, match="samtools faidx"):
        C.shuffle_gc_blocks(tmp_path / "g.fa", tmp_path / "s.fa")


# ---------------------------------------------------------------------------------------- LTR positions
SEQS = {"chr1": 200_000, "chr2": 150_000, "chr3": 100_000}


def ltrquest_gff(path, n=40, seed=5, extra=()):
    """An LTRquest-like annotation: elements with clade, family and age; their LTRs and a domain under them; a gene
    that is no LTR retrotransposon; '###' between elements."""
    rng = np.random.default_rng(seed)
    lines = ["##gff-version 3"] + [f"##sequence-region\t{c}\t1\t{L}" for c, L in SEQS.items()]
    fams = [("CRM", "fam1"), ("Tekay", "fam2"), ("Athila", "fam3")]
    for k in range(n):
        c = list(SEQS)[k % 3]
        w = int(rng.integers(3000, 12000))
        s = int(rng.integers(1, SEQS[c] - w))
        clade, fam = fams[k % 3]
        e, i = f"LTRRT_{k}", s + w // 3
        lines += [f"{c}\tLTRquest\tLTR_retrotransposon\t{s}\t{s + w - 1}\t.\t+\t.\tID={e};Name={c}:{s}-{s + w - 1};"
                  f"clade={clade};family={fam};K2P_T={k * 1000}",
                  f"{c}\tLTRquest\tlong_terminal_repeat\t{s}\t{s + 499}\t.\t+\t.\tID={e}.lLTR;Parent={e}",
                  f"{c}\tLTRquest\tprotein_match\t{i}\t{i + 800}\t.\t+\t.\tID={e}.RT.1;Parent={e};Name=RT;"
                  f"clade={clade}",
                  f"{c}\tLTRquest\tlong_terminal_repeat\t{s + w - 500}\t{s + w - 1}\t.\t+\t.\tID={e}.rLTR;Parent={e}",
                  "###"]
    lines += ["chr2\tmaker\tgene\t5000\t9000\t.\t-\t.\tID=g1", "chr2\tmaker\texon\t5000\t6000\t.\t-\t.\tParent=g1",
              *extra]
    path.write_text("\n".join(lines) + "\n")
    return path


def write_fai(path, seqs=SEQS):
    path.write_text("".join(f"{c}\t{L}\t0\t60\t61\n" for c, L in seqs.items()))
    return path


def features(path):
    rows = [x.split("\t") for x in path.read_text().splitlines() if x and not x.startswith("#")]
    return pd.DataFrame([(r[0], r[2], int(r[3]), int(r[4]), r[8]) for r in rows],
                        columns=["seq", "type", "start", "end", "attributes"])


def test_ltr_shuffle_moves_elements_with_their_parts(tmp_path):
    gff = ltrquest_gff(tmp_path / "t.gff3")
    r = C.shuffle_ltr(gff, tmp_path / "s.gff3", write_fai(tmp_path / "g.fa.fai"), seed=1)
    a, b = features(gff), features(tmp_path / "s.gff3")
    assert r["moved"] == r["trees"] == 40 and r["left"] == 0
    assert (a[["seq", "type", "attributes"]] == b[["seq", "type", "attributes"]]).all().all()
    assert ((b.end - b.start) == (a.end - a.start)).all()           # every feature keeps its length
    shift = b.start - a.start
    tree = a.attributes.str.extract(r"(?:ID|Parent)=(LTRRT_\d+)")[0]
    moved = tree.notna()
    assert (shift[moved].groupby(tree[moved]).nunique() == 1).all()  # an element and its parts move as one
    assert (shift[moved] != 0).mean() > 0.9 and (shift[~moved] == 0).all()    # the gene stays
    assert ((b.start >= 1) & (b.end <= b.seq.map(SEQS))).all()
    keep = lambda p: [x for x in p.read_text().splitlines() if x.startswith("#")]
    assert keep(gff) == keep(tmp_path / "s.gff3")


def test_ltr_shuffle_draws_starts_as_control_te_did(tmp_path):
    gff = ltrquest_gff(tmp_path / "t.gff3")
    C.shuffle_ltr(gff, tmp_path / "s.gff3", write_fai(tmp_path / "g.fa.fai"), seed=3)
    a, b = features(gff), features(tmp_path / "s.gff3")
    el = a.type == "LTR_retrotransposon"
    width = (a.end - a.start + 1)[el].to_numpy()                     # control_te.py: end - start, 0-based half-open
    L = a.seq[el].map(SEQS).to_numpy()
    st = (np.random.default_rng(3).random(int(el.sum())) * np.maximum(1, L - width)).astype(np.int64)
    assert (b.start[el].to_numpy() == st + 1).all()


def test_cenfuse_reads_the_shuffled_annotation_as_the_original(tmp_path):
    gff = ltrquest_gff(tmp_path / "t.gff3")
    C.shuffle_ltr(gff, tmp_path / "s.gff3", write_fai(tmp_path / "g.fa.fai"), seed=1)
    t0, t1 = signature.load_table(gff, SEQS), signature.load_table(tmp_path / "s.gff3", SEQS)
    assert (t0.chrom == t1.chrom).all() and (t0.family == t1.family).all() and (t0.lineage == t1.lineage).all()
    assert ((t0.end - t0.start) == (t1.end - t1.start)).all() and (t0.start != t1.start).any()
    m0, m1 = te.load_te_gff(gff, SEQS), te.load_te_gff(tmp_path / "s.gff3", SEQS)
    assert {k: {c: len(x) for c, x in v.items()} for k, v in m0.items()} == \
        {k: {c: len(x) for c, x in v.items()} for k, v in m1.items()}


def test_ltr_shuffle_moves_a_whole_feature_tree(tmp_path):
    """EDTA's layout: a repeat_region with its TSDs, the LTR retrotransposon under it, and LTRs under that."""
    gff = tmp_path / "edta.gff3"
    gff.write_text("##gff-version 3\n"
                   "chr1\tEDTA\trepeat_region\t1000\t9010\t.\t+\t.\tID=repeat_region_1\n"
                   "chr1\tEDTA\ttarget_site_duplication\t1000\t1004\t.\t+\t.\tID=lTSD_1;Parent=repeat_region_1\n"
                   "chr1\tEDTA\tGypsy_LTR_retrotransposon\t1005\t9005\t.\t+\t.\tID=LTRRT_1;Parent=repeat_region_1;"
                   "clade=CRM\n"
                   "chr1\tEDTA\tlong_terminal_repeat\t1005\t1500\t.\t+\t.\tID=lLTR_1;Parent=LTRRT_1\n"
                   "chr1\tEDTA\tlong_terminal_repeat\t8510\t9005\t.\t+\t.\tID=rLTR_1;Parent=LTRRT_1\n"
                   "chr1\tEDTA\ttarget_site_duplication\t9006\t9010\t.\t+\t.\tID=rTSD_1;Parent=repeat_region_1\n")
    C.shuffle_ltr(gff, tmp_path / "s.gff3", write_fai(tmp_path / "g.fa.fai"), seed=2)
    a, b = features(gff), features(tmp_path / "s.gff3")
    assert (b.start - a.start).nunique() == 1 and (b.start != a.start).all()
    assert b.start.min() >= 1 and b.end.max() <= SEQS["chr1"]


def test_ltr_shuffle_is_deterministic_and_reads_gzip(tmp_path):
    gff = ltrquest_gff(tmp_path / "t.gff3")
    gz = tmp_path / "t.gff3.gz"
    gz.write_bytes(gzip.compress(gff.read_bytes()))
    fai = write_fai(tmp_path / "g.fa.fai")
    for name, src, seed in (("a.gff3", gff, 4), ("b.gff3", gz, 4), ("c.gff3", gff, 5)):
        C.shuffle_ltr(src, tmp_path / name, fai, seed=seed)
    assert (tmp_path / "a.gff3").read_text() == (tmp_path / "b.gff3").read_text()
    assert (tmp_path / "a.gff3").read_text() != (tmp_path / "c.gff3").read_text()


def test_ltr_shuffle_leaves_unknown_sequences_and_stops_on_mismatches(tmp_path, caplog, monkeypatch):
    other = ("scaffold_9\tLTRquest\tLTR_retrotransposon\t100\t5000\t.\t+\t.\tID=X1;clade=CRM;family=fam1",)
    gff = ltrquest_gff(tmp_path / "t.gff3", extra=other)
    fai = write_fai(tmp_path / "g.fa.fai")
    monkeypatch.setattr(logging.getLogger("cenfuse"), "propagate", True)
    with caplog.at_level(logging.WARNING, logger="cenfuse"):
        r = C.shuffle_ltr(gff, tmp_path / "s.gff3", fai)
    assert r["left"] == 1 and r["moved"] == 40
    assert any("scaffold_9" in m and "stay in place" in m for m in caplog.messages)
    assert "scaffold_9\tLTRquest\tLTR_retrotransposon\t100\t5000" in (tmp_path / "s.gff3").read_text()
    with pytest.raises(CenfuseError, match="another assembly or naming"):
        C.shuffle_ltr(gff, tmp_path / "s.gff3", write_fai(tmp_path / "x.fai", {"1": 10, "2": 20}))
    with pytest.raises(CenfuseError, match="outside the sequence"):
        C.shuffle_ltr(gff, tmp_path / "s.gff3", write_fai(tmp_path / "y.fai", {**SEQS, "chr1": 1000}))
    (tmp_path / "genes.gff3").write_text("chr1\tmaker\tgene\t10\t90\t.\t+\t.\tID=g1\n")
    with pytest.raises(CenfuseError, match="no LTR retrotransposon"):
        C.shuffle_ltr(tmp_path / "genes.gff3", tmp_path / "s.gff3", fai)


# ---------------------------------------------------------------------------------------- command line
def test_command_line_makes_every_control(pairs, tmp_path):
    fa = write_fasta(tmp_path / "g.fa", genome(GENOME))
    gff = ltrquest_gff(tmp_path / "t.gff3")
    fai = write_fai(tmp_path / "t.fai")
    assert C.main(["hic-repair", str(pairs), str(tmp_path / "rep.pairs.gz"), "--seed", "2"]) == 0
    assert C.main(["hic-relocate", str(pairs), str(tmp_path / "rel.pairs.gz"), "-v"]) == 0
    assert C.main(["gc-shuffle", str(fa), str(tmp_path / "s.fa")]) == 0
    assert C.main(["ltr-shuffle", str(gff), str(tmp_path / "s.gff3"), "--fai", str(fai)]) == 0
    assert all((tmp_path / f).exists() for f in ("rep.pairs.gz", "rel.pairs.gz", "s.fa", "s.fa.fai", "s.gff3"))
    with gzip.open(tmp_path / "rep.pairs.gz", "rb") as x:
        C.repair_pairs(pairs, tmp_path / "again.pairs.gz", seed=2)
        with gzip.open(tmp_path / "again.pairs.gz", "rb") as y:
            assert x.read() == y.read()


@pytest.mark.parametrize("argv", [["ltr-shuffle", "a.gff3", "b.gff3"], ["hic-repair", "a", "b", "--fai", "x.fai"],
                                  ["hic-repair", "a", "b", "--seed", "-1"], ["hic-shuffle", "a", "b"], ["hic-repair"]])
def test_command_line_usage_errors_exit_2(argv):
    with pytest.raises(SystemExit) as e:
        C.main(argv)
    assert e.value.code == 2


def test_command_line_reports_bad_input(tmp_path, capsys):
    assert C.main(["hic-repair", str(tmp_path / "none.pairs.gz"), str(tmp_path / "out.pairs.gz")]) == 2
    assert "cenfuse: error:" in capsys.readouterr().err
