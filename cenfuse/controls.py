"""Negative controls: an input with the signal a line reads destroyed and the rest kept (future.md #34). Run cenfuse
on a copy of GENOME_DIR with the control in place of the input: no test should pass and no call should change. Every
control is deterministic given its seed (numpy's default_rng).

  hic-repair    the benchmark's NEG2: read ends stay where they are (chromosome, position, strand, MAPQ) and are
                re-paired at random across all ends. Every bin keeps its coverage, so mappability, copy-number and
                mis-mapping biases survive; every contact, cis or trans, is random.
  hic-relocate  NEG: every read end moves to a uniform random position on its own chromosome. Each contact keeps its
                chromosome pair (the cis share and the contacts between chromosomes are kept); positions, coverage
                included, carry no structure.
  gc-shuffle    the AT-trough line's control: the chromosomes, end to end, cut into blocks of 100 kb that are put back
                in random order and cut again at the chromosome lengths. Local base composition stays; any pattern of
                one locus per chromosome goes.
  ltr-shuffle   the LTR line's control: every LTR retrotransposon, with the features under it, moves to a uniform
                random position on its own sequence, keeping its length and its attributes (clade, family, age).

The Hi-C controls are v2/scripts/permute_pairs.py (--mode repair | relocate) as a module: a seed writes the same bytes
as that script (checked on the stored NEG2 control of C. tropicalis), because the same numbers are drawn in the same
order. Both test their own output with dispersion() and stop when its contacts are clumped: the re-pairing of
v11/bench/repair_pairs.sh (`shuf --random-source=<(yes 7)`, withdrawn) gave clumped contacts that nothing caught until a
figure measured them. The GC and LTR shuffles follow v11/bench/at_shuffle.py and control_te.py, which drew theirs in
memory; these write them out as inputs for a run."""
from __future__ import annotations

import argparse
import gzip
import io
import os
import sys
import zlib
from array import array
from contextlib import contextmanager
from pathlib import Path

import numpy as np
import pandas as pd

from .util import CenfuseError, log, setup_logging

CHUNK = 2_000_000       # pairs read and written at a time. hic-repair numbers the read ends by these blocks (each
                        # block's first ends, then its second ends), as permute_pairs.py did: a seed gives its controls
                        # only with this value
END_COLS = (("chrom1", "pos1", "strand1", "mapq1"), ("chrom2", "pos2", "strand2", "mapq2"))
DISP_BINS = 1024        # the finest genome-wide bins of the dispersion check: a power of 2, so that halving keeps the
                        # bins nested; 1024 bins make 524,800 bin pairs
DISP_MIN_E = 5          # a bin pair is tested when its count's null variance is >= 5 (Cochran's rule for Pearson's
                        # statistic: >= 5 expected for Poisson counts)
DISP_MIN_CELLS = 10     # fewer testable bin pairs: untested, with a warning
DISP_Z = 10             # clumped above 1 + 10 null sd (the sd of the index for independent counts; the spread measured
                        # is 0.85-1.1x it). Random controls sit at 1: 40 seeds of C. tropicalis Hi-C (11 M pairs, 1024
                        # bins) gave 0.9999 +- 0.0017 re-paired and 1.0002 +- 0.0047 relocated, at most +2.2 sd; 200
                        # seeds of a 3,000-pair fixture at most +3.0 sd; the 14 stored NEG2 controls -1.4 to +0.5 sd,
                        # the stored NEG of mouse MMU8 +1.7 sd. The clumped `shuf` partners (v11/data/neg) scored
                        # 176-1,712, over 88,000 sd
GC_BLOCK = 100_000      # at_shuffle.py's blocks (100 bins of 1 kb): the benchmark's GC control
FASTA_WIDTH = 60        # bases per FASTA line written (samtools faidx's)
COPY = 10_000_000       # bp read at a time when a sequence is copied
KINDS = ("hic-repair", "hic-relocate", "gc-shuffle", "ltr-shuffle")


class ControlError(CenfuseError):
    """A control failed its own check (clumped contacts): the randomization is broken."""


# ------------------------------------------------------------------------------------------------------- 4DN pairs
def _is_gzip(path) -> bool:
    with open(path, "rb") as f:
        return f.read(2) == b"\x1f\x8b"


def read_header(path) -> tuple[dict, list, list]:
    """(chromosome lengths of the '#chromsize:' lines, in file order; the header lines, '#sorted:' made 'none'; the
    '#columns:' names) of a 4DN .pairs[.gz] file, as permute_pairs.py reads them. Other header lines are kept as they
    are ('#shape: upper triangle' too, which a randomized file no longer is; cenfuse does not read it)."""
    sizes, header, cols, first = {}, [], None, ""
    try:
        with (gzip.open(path, "rt") if _is_gzip(path) else open(path)) as f:
            for line in f:
                if not line.startswith("#"):
                    first = line
                    break
                if line.startswith("#chromsize:"):
                    p = line.split()
                    if len(p) != 3 or not p[2].isdigit() or int(p[2]) < 1:
                        raise CenfuseError(f"{path}: malformed header line {line.strip()!r} (want #chromsize: NAME "
                                           f"LENGTH)")
                    sizes[p[1]] = int(p[2])
                elif line.startswith("#columns:"):
                    cols = line.split()[1:]
                header.append("#sorted: none\n" if line.startswith("#sorted:") else line)
    except (OSError, EOFError, zlib.error, UnicodeDecodeError) as e:
        raise CenfuseError(f"{path}: cannot read ({e})") from None
    if not sizes or not cols or not {"chrom1", "pos1", "chrom2", "pos2"} <= set(cols):
        raise CenfuseError(f"{path}: not a 4DN .pairs file: needs '#chromsize:' header lines and a '#columns:' line "
                           f"naming chrom1 pos1 chrom2 pos2 (pairtools writes them)")
    if len(set(cols)) < len(cols):
        raise CenfuseError(f"{path}: a column is named twice in {' '.join(cols)!r}")
    nf = len(first.rstrip("\r\n").split("\t"))
    if first.strip() and nf != len(cols):        # pandas would read a longer first line as an index, shifting columns
        raise CenfuseError(f"{path}: the first data line has {nf} fields, the '#columns:' line names {len(cols)}")
    return sizes, header, cols


def _chunks(path, cols: list, n_header: int, usecols=None, ints=()):
    """The data lines in blocks of CHUNK rows, fields as written (no NA parsing); columns `ints` as int64. A line
    with too few fields (an empty last field), or too many, or an unparsable integer stops the run."""
    use = usecols or cols
    last = use[-1] if use[-1] not in ints else None     # an empty integer field already fails to parse
    try:
        reader = pd.read_csv(path, sep="\t", header=None, names=cols, usecols=usecols, skiprows=n_header,
                             index_col=False, dtype={c: np.int64 if c in ints else str for c in use},
                             na_filter=False, chunksize=CHUNK, compression="gzip" if _is_gzip(path) else None)
        for ch in reader:
            if last is not None and len(ch):
                empty = (ch[last] == "").to_numpy()
                if empty.any():
                    raise CenfuseError(f"{path}: a data line near line {n_header + int(np.flatnonzero(empty)[0]) + 1} "
                                       f"has too few fields")
            yield ch
    except (ValueError, OSError, EOFError, zlib.error) as e:
        raise CenfuseError(f"{path}: unreadable pairs data ({str(e).strip()})") from None


def _chrom_index(s: pd.Series, names: pd.Index, path) -> np.ndarray:
    i = names.get_indexer(s)
    if (i < 0).any():
        c = s.iat[int(np.flatnonzero(i < 0)[0])]
        raise CenfuseError(f"{path}: chromosome {c!r} has no '#chromsize:' line"
                           + (" (unmapped ends: keep mapped pairs only, e.g. pairtools select)" if c == "!" else ""))
    return i


class _Codes:
    """Strings as 2-byte codes and back: the strands and MAPQs of all read ends, held in memory."""

    def __init__(self):
        self.index = {}

    def encode(self, values) -> np.ndarray:
        codes, uniq = pd.factorize(values)
        lut = np.array([self.index.setdefault(u, len(self.index)) for u in uniq], dtype=np.int64)
        if len(self.index) > 65_535:
            raise CenfuseError("more than 65,535 distinct strand or MAPQ values: not a pairs file")
        return lut[codes].astype(np.uint16) if len(codes) else np.zeros(0, np.uint16)

    def decode(self, codes: np.ndarray) -> np.ndarray:
        return np.array(list(self.index), dtype=object)[codes]


@contextmanager
def _text_out(path, gz: bool):
    """A text file to write: gzip at level 3 (permute_pairs.py's) without a name or time stamp in its header, so that a
    seed gives the same bytes; or plain."""
    with open(path, "wb") as raw:
        if gz:
            with gzip.GzipFile(filename="", mode="wb", fileobj=raw, compresslevel=3, mtime=0) as z, \
                    io.TextIOWrapper(z, encoding="utf-8", newline="") as f:
                yield f
        else:
            with io.TextIOWrapper(raw, encoding="utf-8", newline="") as f:
                yield f


def _paths(src, dst) -> tuple[Path, Path, Path]:
    """(input, output, the temporary output beside it), after checks."""
    src, dst = Path(src), Path(dst)
    if not src.is_file():
        raise CenfuseError(f"{src}: no such file")
    if not dst.parent.is_dir():
        raise CenfuseError(f"{dst.parent}: no such directory")
    if dst.exists() and dst.resolve() == src.resolve():
        raise CenfuseError(f"{dst}: the output would overwrite the input")
    return src, dst, dst.with_name(dst.name + ".part")


def _pairing(m: int, rng) -> tuple[np.ndarray, np.ndarray]:
    """Partners for m read ends: a uniform random permutation taken two by two (one end left out when m is odd). The
    draws of permute_pairs.py's rng.permutation(m), in half its memory below 2^31 ends."""
    order = np.arange(m, dtype=np.int32 if m < 2 ** 31 else np.int64)
    rng.shuffle(order)
    n = m // 2
    return order[0:2 * n:2], order[1:2 * n:2]


def repair_pairs(in_pairs, out_pairs, seed: int = 1) -> dict:
    """NEG2: every read end keeps its chromosome, position, strand and MAPQ, and partners are re-paired at random
    across all ends (readIDs become r<k>, pair_type UU, other columns '.'). Holds all ends in memory: about 40 bytes
    per pair at the peak (measured 2.8 GB and 5 min for 53 M pairs; 1.2 GB for 11 M). Writes OUT only after
    dispersion() passed; returns its result with the pair count."""
    src, dst, tmp = _paths(in_pairs, out_pairs)
    sizes, header, cols = read_header(src)
    names = pd.Index(list(sizes))
    lens = np.array(list(sizes.values()), dtype=np.int64)
    have = {c for c in ("strand1", "strand2", "mapq1", "mapq2") if c in cols}
    use = [c for c in cols if c in {"chrom1", "pos1", "chrom2", "pos2"} | have]
    ctype, ptype = (np.int16 if len(names) < 2 ** 15 else np.int32), (np.int32 if lens.max() < 2 ** 31 else np.int64)
    strand, mapq = _Codes(), _Codes()
    parts = {"chrom": [], "pos": [], "strand": [], "mapq": []}
    done = 0
    for ch in _chunks(src, cols, len(header), usecols=use, ints=("pos1", "pos2")):
        for cc, pc, sc, qc in END_COLS:      # the chunk's first ends, then its second ends: permute_pairs.py's order
            ci = _chrom_index(ch[cc], names, src)
            p = ch[pc].to_numpy(np.int64)
            bad = (p < 1) | (p > lens[ci])
            if bad.any():
                k = int(np.flatnonzero(bad)[0])
                raise CenfuseError(f"{src}: pair {done + k + 1:,}: {pc} {p[k]} lies outside {names[ci[k]]} (1-"
                                   f"{lens[ci[k]]:,} by its '#chromsize:' line)")
            parts["chrom"].append(ci.astype(ctype))
            parts["pos"].append(p.astype(ptype))
            parts["strand"].append(strand.encode(ch[sc] if sc in have else np.full(len(ch), ".")))
            parts["mapq"].append(mapq.encode(ch[qc] if qc in have else np.full(len(ch), ".")))
        done += len(ch)
        log.debug(f"hic-repair: {done:,} pairs read")
    if not done:
        raise CenfuseError(f"{src}: no pairs to re-pair")
    e = {k: np.concatenate(parts.pop(k)) for k in list(parts)}
    i, j = _pairing(len(e["pos"]), np.random.default_rng(seed))
    n = len(i)
    labels = np.array(list(names), dtype=object)
    try:
        with _text_out(tmp, dst.name.endswith(".gz")) as out:
            out.write("".join(header))
            for s in range(0, n, CHUNK):
                ii, jj = i[s:s + CHUNK], j[s:s + CHUNK]
                d = {"readID": [f"r{k}" for k in range(s, s + len(ii))],
                     "chrom1": labels[e["chrom"][ii]], "pos1": e["pos"][ii],
                     "chrom2": labels[e["chrom"][jj]], "pos2": e["pos"][jj],
                     "strand1": strand.decode(e["strand"][ii]), "strand2": strand.decode(e["strand"][jj]),
                     "pair_type": "UU", "mapq1": mapq.decode(e["mapq"][ii]), "mapq2": mapq.decode(e["mapq"][jj])}
                pd.DataFrame({c: d.get(c, ".") for c in cols}).to_csv(out, sep="\t", header=False, index=False,
                                                                       lineterminator="\n")
                log.debug(f"hic-repair: {s + len(ii):,} pairs written")
    except BaseException:
        tmp.unlink(missing_ok=True)
        raise
    return _finish(tmp, dst, "repair", n, seed)


def relocate_pairs(in_pairs, out_pairs, seed: int = 1) -> dict:
    """NEG: every read end moves to a uniform random position on its own chromosome (1-based, 1 to its length); the
    other columns are kept as they are. Streams. Writes OUT only after dispersion() passed; returns its result with the
    pair count."""
    src, dst, tmp = _paths(in_pairs, out_pairs)
    sizes, header, cols = read_header(src)
    names = pd.Index(list(sizes))
    lens = np.array(list(sizes.values()), dtype=np.float64)
    rng = np.random.default_rng(seed)
    n = 0
    try:
        with _text_out(tmp, dst.name.endswith(".gz")) as out:
            out.write("".join(header))
            for ch in _chunks(src, cols, len(header)):
                for c, p in (("chrom1", "pos1"), ("chrom2", "pos2")):
                    L = lens[_chrom_index(ch[c], names, src)]
                    ch[p] = np.floor(rng.random(len(ch)) * L).astype(np.int64) + 1
                ch.to_csv(out, sep="\t", header=False, index=False, lineterminator="\n")
                n += len(ch)
                log.debug(f"hic-relocate: {n:,} pairs written")
    except BaseException:
        tmp.unlink(missing_ok=True)
        raise
    if n == 0:
        tmp.unlink(missing_ok=True)
        raise CenfuseError(f"{src}: no pairs to relocate")
    return _finish(tmp, dst, "relocate", n, seed)


def _finish(tmp: Path, dst: Path, null: str, n: int, seed: int) -> dict:
    """Check the control written at tmp; move it to dst when its contacts are not clumped."""
    r = dispersion(tmp, null)
    what = (f"dispersion {r['dispersion']:.3f}, 1 for random contacts (bound {r['bound']:.3f}; {r['cells']:,} bin "
            f"pairs of {r['bins']} genome bins)")
    if r["status"] == "fail":
        raise ControlError(f"{dst}: the {'re-paired' if null == 'repair' else 'relocated'} contacts are clumped, "
                           f"{what}: the randomization is broken (a bug: please report it). The file is left at {tmp}")
    if r["status"] == "untested":
        log.warning(f"{dst}: too few contacts to test for clumps ({r['cells']} bin pairs expect >= {DISP_MIN_E})")
    os.replace(tmp, dst)
    log.info(f"hic-{null}: {n:,} pairs written to {dst} (seed {seed}); "
             + (f"contacts not clumped: {what}" if r["status"] == "pass" else "not tested for clumps"))
    return dict(r, pairs=n)


# ------------------------------------------------------------------------------------------------------ dispersion
def dispersion(pairs, null: str = "repair") -> dict:
    """Are the contacts of a randomized .pairs[.gz] file clumped? A dispersion index of its contacts per pair of
    genome-wide bins: the mean over bin pairs of (observed - expected)^2 / variance, both under the null:
      repair    partners at random given every bin's coverage (c: read ends per bin, M: all ends): expected
                c_a c_b / (M - 1) contacts between bins a and b and c_a (c_a - 1) / 2(M - 1) within bin a; variance
                the expected count times (1 - c_a/M)(1 - c_b/M), the correction for drawing partners without return;
      relocate  ends uniform on their sequences given each sequence pair's contacts: binomial counts, each contact
                landing in a bin pair with the product of its two sequences' length shares there.
    About 1 for a correct control, at any bin size; clumped contacts give far more (about the mean clump size where
    clumps fall in tested bin pairs). Bins: the '#chromsize:' sequences end to end, cut into 2^k equal bins (k <= 10),
    the k with the most bin pairs whose variance is >= DISP_MIN_E; only those are tested. Streams.
    Returns dict(status pass | fail | untested, dispersion, bound (1 + DISP_Z sd), sd (of the index for independent
    counts, sqrt(sum(2 + 1/variance)) / cells), bins, cells, null, pairs)."""
    if null not in ("repair", "relocate"):
        raise ValueError(f"null must be 'repair' or 'relocate', not {null!r}")
    sizes, header, cols = read_header(pairs)
    names = pd.Index(list(sizes))
    L = np.array(list(sizes.values()), dtype=np.int64)
    off = np.concatenate([[0], np.cumsum(L)])
    G, B, nc = int(off[-1]), DISP_BINS, len(L)
    O = np.zeros(B * B, dtype=np.int64)                 # oriented counts: (bin of end 1, bin of end 2)
    T = {}                                              # contacts per oriented sequence pair (relocate)
    n = 0
    for ch in _chunks(pairs, cols, len(header), usecols=["chrom1", "pos1", "chrom2", "pos2"], ints=("pos1", "pos2")):
        c1, c2 = _chrom_index(ch.chrom1, names, pairs), _chrom_index(ch.chrom2, names, pairs)
        b1 = (off[c1] + np.clip(ch.pos1.to_numpy(), 1, L[c1]) - 1) * B // G
        b2 = (off[c2] + np.clip(ch.pos2.to_numpy(), 1, L[c2]) - 1) * B // G
        O += np.bincount(b1 * B + b2, minlength=B * B)
        if null == "relocate":
            k, m = np.unique(c1.astype(np.int64) * nc + c2, return_counts=True)
            for a, b in zip(k.tolist(), m.tolist()):
                T[a] = T.get(a, 0) + b
        n += len(ch)
    O = O.reshape(B, B)
    R = _relocate_null(T, L, B, G) if null == "relocate" else None
    best, k = None, B
    while k >= 2 and n:
        r = _index(O, k, n, R)
        log.debug(f"dispersion ({null}): {k} bins, {r['cells']:,} testable bin pairs, index {r['dispersion']:.3f}")
        if best is None or r["cells"] > best["cells"]:
            best = r
        k //= 2
    if best is None or best["cells"] < DISP_MIN_CELLS:
        return dict(status="untested", dispersion=float("nan"), bound=float("nan"), sd=float("nan"),
                    bins=best["bins"] if best else 0, cells=best["cells"] if best else 0, null=null, pairs=n)
    bound = 1 + DISP_Z * best["sd"]
    return dict(best, status="fail" if best["dispersion"] > bound else "pass", bound=bound, null=null, pairs=n)


def _relocate_null(T: dict, L: np.ndarray, B: int, G: int) -> tuple:
    """(P, T, D) of the relocate null: the share of each sequence's length in each of B genome bins (sparse, sequences
    x bins; a bin may hold several sequences), the contacts per oriented sequence pair, and its diagonal (cis)."""
    from scipy import sparse
    nc = len(L)
    off = np.concatenate([[0], np.cumsum(L)])
    first = -((-np.arange(B + 1, dtype=np.int64) * G) // B)       # first genome position of each bin
    r, c, v = [], [], []
    for i in range(nc):
        s, e = off[i], off[i + 1]
        a = np.arange(s * B // G, (e - 1) * B // G + 1)
        r.append(np.full(len(a), i))
        c.append(a)
        v.append((np.minimum(e, first[a + 1]) - np.maximum(s, first[a])) / L[i])
    P = sparse.csr_matrix((np.concatenate(v), (np.concatenate(r), np.concatenate(c))), shape=(nc, B))
    k = np.fromiter(T.keys(), dtype=np.int64, count=len(T))
    m = np.fromiter(T.values(), dtype=np.float64, count=len(T))
    Tm = sparse.csr_matrix((m, (k // nc, k % nc)), shape=(nc, nc))
    return P, Tm, sparse.diags(Tm.diagonal()).tocsr()


def _moments(Ok: np.ndarray, n: int, R=None) -> tuple:
    """(observed, expected, variance) of the contacts in every unordered pair of the bins of Ok (oriented counts,
    k x k), in upper-triangle order, under the null. R: the relocate null with k bins (P with k columns), else the
    repair null (n pairs)."""
    k = len(Ok)
    iu = np.triu_indices(k)
    diag = iu[0] == iu[1]
    obs = (Ok + Ok.T)[iu].astype(np.float64)
    obs[diag] /= 2
    if R is None:                                       # repair: a uniform random matching of the M ends
        c = (Ok.sum(axis=0) + Ok.sum(axis=1)).astype(np.float64)
        M = 2 * n
        p, w = 1 / (M - 1), 1 / ((M - 1) * (M - 3))   # P(two ends matched); P(two given disjoint pairs) - p^2 = 2pw
        ca, cb = c[iu[0]], c[iu[1]]
        exp = ca * cb * p
        exp[diag] = (ca * (ca - 1) * p / 2)[diag]
        f = 1 - p * (ca + cb - 1) + 2 * (ca - 1) * (cb - 1) * w     # variance / expected, between bins
        f[diag] = (1 - p * (2 * ca - 3) + (ca - 2) * (ca - 3) * w)[diag]     # within a bin
        return obs, exp, exp * f
    P, Tm, D = R                                        # relocate: a contact of sequences i -> j lands in bins a, b
    P2 = P.multiply(P).tocsr()                          # with q = p_i(a) p_j(b) + p_i(b) p_j(a); variance sum q(1 - q)
    E = (P.T @ Tm @ P).toarray()                        # oriented: end 1 in a, end 2 in b
    S = (P2.T @ Tm @ P2).toarray()
    X = (P2.T @ D @ P2).toarray()                       # the cross term of q^2: one sequence spans both bins (cis)
    exp = (E + E.T)[iu]
    exp[diag] /= 2
    q2 = (S + S.T + 2 * X)[iu]
    q2[diag] = np.diag(S)
    return obs, exp, exp - q2


def _index(O: np.ndarray, k: int, n: int, R=None) -> dict:
    """The dispersion index at k bins (merged from O's B), over the unordered bin pairs whose null variance is
    >= DISP_MIN_E. R: the relocate null (_relocate_null, B bins), else the repair null."""
    from scipy import sparse
    B, f = O.shape[0], O.shape[0] // k
    Ok = O.reshape(k, f, k, f).sum(axis=(1, 3))
    if R is not None:
        R = (R[0] @ sparse.csr_matrix((np.ones(B), (np.arange(B), np.arange(B) // f)), shape=(B, k)), R[1], R[2])
    obs, exp, var = _moments(Ok, n, R)
    t = var >= DISP_MIN_E
    K = int(t.sum())
    if not K:
        return dict(bins=k, cells=0, dispersion=float("nan"), sd=float("nan"))
    x2 = float(((obs[t] - exp[t]) ** 2 / var[t]).sum())
    return dict(bins=k, cells=K, dispersion=x2 / K, sd=float(np.sqrt((2 + 1 / var[t]).sum()) / K))


# ------------------------------------------------------------------------------------------------------- GC blocks
def _fasta_index(fasta: Path) -> dict:
    """name -> length from the FASTA's .fai (as `cenfuse run` requires; bgzip also needs .gzi)."""
    from .sat.genomes import load_fai
    from .sat.util import SatelliteError
    if not fasta.is_file():
        raise CenfuseError(f"{fasta}: no such file")
    if fasta.name.endswith(".gz") and not Path(f"{fasta}.gzi").exists():
        raise CenfuseError(f"{fasta}: a compressed FASTA must be bgzip with .fai + .gzi (run: samtools faidx {fasta})")
    try:
        return load_fai(fasta)
    except SatelliteError as e:
        raise CenfuseError(str(e)) from None


def _chromosome_set(fasta: Path, fai: dict) -> tuple[dict, str]:
    """The chromosomes `cenfuse run` would analyze in the FASTA's folder: <genome>.chromosomes.txt, else the length
    rule; and what set them."""
    from .inputs import CHROM_LIST, RULE_SOURCE, chromosome_set, read_chromosome_list
    from .sat.genomes import _fasta_name
    listed = fasta.parent / f"{_fasta_name(fasta) or fasta.name}{CHROM_LIST}"
    if listed.is_file():
        return read_chromosome_list(listed, fai), listed.name
    return chromosome_set(fai), RULE_SOURCE


class _FastaWriter:
    """FASTA lines of FASTA_WIDTH bases; bgzip when the name ends in .gz (indexed after, as cenfuse needs)."""

    def __init__(self, path: Path):
        import pysam
        self.f = pysam.BGZFile(str(path), "wb") if path.name.endswith(".gz") else open(path, "wb")
        self.col = 0

    def header(self, name: str) -> None:
        self.f.write((b"\n" if self.col else b"") + b">" + name.encode() + b"\n")
        self.col = 0

    def write(self, s: str) -> None:
        b, W, out = s.encode("latin-1"), FASTA_WIDTH, []
        if self.col:                                    # finish the open line
            t = b[:W - self.col]
            out.append(t)
            self.col += len(t)
            b = b[len(t):]
            if self.col == W:
                out.append(b"\n")
                self.col = 0
        full = len(b) - len(b) % W
        if full:
            a = np.frombuffer(b[:full], dtype=np.uint8).reshape(-1, W)
            out.append(np.hstack([a, np.full((len(a), 1), 10, np.uint8)]).tobytes())
        if full < len(b):
            out.append(b[full:])
            self.col = len(b) - full
        self.f.write(b"".join(out))

    def close(self) -> None:
        if self.col:
            self.f.write(b"\n")
        self.f.close()


def shuffle_gc_blocks(fasta, out_fasta, seed: int = 1, block: int = GC_BLOCK, chroms=None) -> dict:
    """The control of the AT-trough line (at_shuffle.py on the sequence itself): the chromosomes, end to end in .fai
    order, are cut into blocks of `block` bp from the start (the last one shorter), the blocks put back in a random
    order and the result cut again at the chromosome lengths. Each chromosome keeps its length and the set its base
    composition (case too); a block may straddle two chromosomes, as at_shuffle.py's did. chroms: the names to shuffle
    among (default: the chromosome set `cenfuse run` would use); other sequences are copied unchanged, so a run on the
    output picks the same set. Writes OUT (bgzip if .gz) and its index. Streams, one block at a time.
    at_shuffle.py drew its 20 shuffles per genome from one stream (seed 7) on 1-kb GC bins, so no seed here repeats
    them."""
    import pysam
    if block < 1:
        raise CenfuseError(f"block must be >= 1 bp, not {block}")
    fasta, out, _ = _paths(fasta, out_fasta)
    fai = _fasta_index(fasta)
    if chroms is None:
        chroms, source = _chromosome_set(fasta, fai)
    else:
        source = "the names given"
        unknown = [c for c in chroms if c not in fai]
        if unknown:
            raise CenfuseError(f"{fasta}: not in its .fai: {', '.join(map(str, unknown[:5]))}")
    names = [c for c in fai if c in chroms]
    pos = {c: i for i, c in enumerate(names)}
    soff = np.concatenate([[0], np.cumsum([fai[c] for c in names])]).astype(np.int64)
    G = int(soff[-1])
    if not G:
        raise CenfuseError(f"{fasta}: no chromosome sequence to shuffle")
    nb = -(-G // block)
    order = np.random.default_rng(seed).permutation(nb)
    land = np.concatenate([[0], np.cumsum(np.minimum(block, G - order * block))])   # where each block lands
    made = [out, Path(f"{out}.fai"), Path(f"{out}.gzi")]
    for p in made:                                      # no stale index beside a new FASTA
        p.unlink(missing_ok=True)
    try:
        with pysam.FastaFile(str(fasta)) as fa:
            def source_seq(u: int, v: int):             # the chromosomes end to end, [u, v)
                c = int(np.searchsorted(soff, u, side="right")) - 1
                while u < v:
                    e = min(v, int(soff[c + 1]))
                    yield fa.fetch(names[c], u - int(soff[c]), e - int(soff[c]))
                    u, c = e, c + 1

            w = _FastaWriter(out)
            try:
                for name, Lc in fai.items():
                    w.header(name)
                    if name not in pos:
                        for s in range(0, Lc, COPY):
                            w.write(fa.fetch(name, s, min(Lc, s + COPY)))
                        continue
                    x, y = int(soff[pos[name]]), int(soff[pos[name] + 1])
                    k = int(np.searchsorted(land, x, side="right")) - 1
                    while x < y:
                        e = min(y, int(land[k + 1]))
                        u = int(order[k]) * block + x - int(land[k])
                        for piece in source_seq(u, u + e - x):
                            w.write(piece)
                        x, k = e, k + 1
            finally:
                w.close()
        pysam.faidx(str(out))
    except BaseException:
        for p in made:
            p.unlink(missing_ok=True)
        raise
    with pysam.FastaFile(str(out)) as fo:
        if dict(zip(fo.references, fo.lengths)) != fai:
            raise CenfuseError(f"{out}: lengths differ from {fasta}.fai after the shuffle (is the .fai stale?)")
    log.info(f"gc-shuffle: {len(names)} chromosomes ({G:,} bp; set from {source}) cut into {nb:,} blocks of {block:,} "
             f"bp and shuffled (seed {seed}): {out}")
    if source.endswith(".chromosomes.txt"):
        log.info(f"gc-shuffle: copy {source} beside {out.name}, named for its genome, so that a run uses the same set")
    return dict(chromosomes=len(names), bp=G, blocks=nb, block=block, source=source)


# --------------------------------------------------------------------------------------------------- LTR positions
def _gff_lines(path):
    """(feature fields or None, line, malformed) for every line of a GFF3[.gz]: a feature has 9 tab-separated fields;
    a line with fewer is malformed; comments, blank lines and the sequences after '##FASTA' are neither."""
    fasta = False
    with (gzip.open(path, "rt") if _is_gzip(path) else open(path)) as f:
        for line in f:
            if fasta or line.startswith("#") or not line.strip():
                fasta = fasta or line.startswith("##FASTA")
                yield None, line, False
                continue
            p = line.rstrip("\r\n").split("\t")
            yield (p, line, False) if len(p) >= 9 else (None, line, True)


def _attrs(col: str) -> dict:
    return {k.strip().lower(): v for k, _, v in (kv.partition("=") for kv in col.split(";")) if v}


def _is_ltr(ftype: str, a: dict) -> bool:
    """An LTR retrotransposon as cenfuse's LTR line reads one (sat.te): a type holding 'LTR_retrotransposon' (LTRquest,
    EDTA), or a DANTE_LTR transposable_element of an LTR lineage."""
    return "LTR_retrotransposon" in ftype or (ftype == "transposable_element"
                                              and "|LTR|" in a.get("final_classification", ""))


def shuffle_ltr(gff, out_gff, fai, seed: int = 1) -> dict:
    """The control of the LTR line (control_te.py, written out): every feature tree (a top-level feature and all
    under it by ID/Parent) holding an LTR retrotransposon moves as one to a uniform random position on its sequence,
    drawn as control_te.py drew an element's start: start0 = floor(u * max(1, L - width)), u ~ U(0, 1), one draw per
    tree in the order of their first LTR line. LTRquest's elements are trees of their own (their LTRs and domains under
    them), so with a .fai of the chromosome set this repeats control_te.py's shuffle of the same seed. Columns other
    than start and end, other features and comment lines are kept (Name= may still hold the old coordinates). Trees on
    sequences absent from the .fai stay in place (warned). Reads the GFF twice; holds one entry per feature."""
    from .hic.io import InputError, read_fai
    gff, out, _ = _paths(gff, out_gff)
    try:
        lengths = read_fai(fai)
    except (InputError, OSError) as e:
        raise CenfuseError(str(e)) from None
    ids, up = {}, []                                    # union-find over features: ID -> node, node -> parent node

    def node(key=None):
        if key is not None and key in ids:
            return ids[key]
        up.append(len(up))
        if key is not None:
            ids[key] = len(up) - 1
        return len(up) - 1

    def find(x):
        while up[x] != x:
            up[x] = up[up[x]]
            x = up[x]
        return x

    seq_code = {}
    line_node, line_seq, start, end, ltr = array("q"), array("q"), array("q"), array("q"), []
    short = 0
    try:
        for k, (p, line, bad) in enumerate(_gff_lines(gff)):
            short += bad
            if p is None:
                continue
            a = _attrs(p[8])
            x = node(a.get("id"))
            for par in a.get("parent", "").split(","):
                if par:
                    ra, rb = find(x), find(node(par))
                    up[max(ra, rb)] = min(ra, rb)
            try:
                s0, e0 = int(p[3]), int(p[4])
            except ValueError:
                raise CenfuseError(f"{gff}: line {k + 1}: start and end are not integers") from None
            line_node.append(x)
            line_seq.append(seq_code.setdefault(p[0], len(seq_code)))
            start.append(s0)
            end.append(e0)
            ltr.append(_is_ltr(p[2], a))
    except (OSError, EOFError, zlib.error, UnicodeDecodeError) as e:
        raise CenfuseError(f"{gff}: cannot read ({e})") from None
    seqs = list(seq_code)
    ltr = np.array(ltr, dtype=bool)
    if not ltr.any():
        raise CenfuseError(f"{gff}: no LTR retrotransposon (a type holding LTR_retrotransposon, or a DANTE_LTR "
                           f"transposable_element)")
    root = np.array([find(x) for x in line_node], dtype=np.int64)
    trees, tree = np.unique(root, return_inverse=True)
    sq, st, en = np.frombuffer(line_seq, np.int64), np.frombuffer(start, np.int64), np.frombuffer(end, np.int64)
    nt = len(trees)
    lo, hi = np.full(nt, np.iinfo(np.int64).max), np.full(nt, np.iinfo(np.int64).min)
    np.minimum.at(lo, tree, st)
    np.maximum.at(hi, tree, en)
    smin, smax = np.full(nt, len(seqs)), np.full(nt, -1)
    np.minimum.at(smin, tree, sq)
    np.maximum.at(smax, tree, sq)
    first = np.full(nt, len(ltr))                       # each tree's first LTR line: the order of the draws
    np.minimum.at(first, tree[ltr], np.flatnonzero(ltr))
    holds = first < len(ltr)
    L = np.array([lengths.get(s, -1) for s in seqs], dtype=np.int64)
    known = L[smin] > 0
    move = holds & known
    if not move.any():
        raise CenfuseError(f"{gff}: none of its LTR retrotransposons lies on a sequence of {fai} (another assembly or "
                           f"naming?)")
    bad = move & (smin != smax)
    if bad.any():
        raise CenfuseError(f"{gff}: a feature tree spans several sequences ({seqs[smin[bad][0]]}, "
                           f"{seqs[smax[bad][0]]}): IDs used twice?")
    bad = move & ((lo < 1) | (hi > L[smin]))
    if bad.any():
        t = int(np.flatnonzero(bad)[0])
        raise CenfuseError(f"{gff}: a feature at {seqs[smin[t]]}:{lo[t]}-{hi[t]} lies outside the sequence "
                           f"({L[smin[t]]:,} bp in {fai}): another assembly?")
    mv = np.flatnonzero(move)
    mv = mv[np.argsort(first[mv], kind="stable")]
    width = hi[mv] - lo[mv] + 1
    u = np.random.default_rng(seed).random(len(mv))
    shift = np.zeros(nt, dtype=np.int64)
    shift[mv] = (u * np.maximum(1, L[smin[mv]] - width)).astype(np.int64) + 1 - lo[mv]
    moved, left = int((ltr & move[tree]).sum()), int((ltr & ~move[tree]).sum())     # LTR retrotransposons
    if left:
        miss = sorted({seqs[s] for s in smin[holds & ~known]})
        log.warning(f"ltr-shuffle: {left:,} LTR retrotransposons on {len(miss)} sequences absent from {fai} "
                    f"({', '.join(miss[:3])}{', ...' if len(miss) > 3 else ''}) stay in place")
    if short:
        log.warning(f"ltr-shuffle: {short:,} lines with fewer than 9 tab-separated fields copied unchanged")
    i = 0
    try:
        with _text_out(out, out.name.endswith(".gz")) as f:
            for p, line, _ in _gff_lines(gff):
                if p is not None:
                    d = int(shift[tree[i]])
                    i += 1
                    if d:
                        p[3], p[4] = str(int(p[3]) + d), str(int(p[4]) + d)
                        line = "\t".join(p) + line[len(line.rstrip("\r\n")):]
                f.write(line)
    except BaseException:
        out.unlink(missing_ok=True)
        raise
    log.info(f"ltr-shuffle: {moved:,} LTR retrotransposons ({len(mv):,} feature trees) moved within "
             f"{len(set(smin[mv].tolist()))} sequences (seed {seed}): {out}")
    return dict(moved=moved, left=left, trees=len(mv))


# ---------------------------------------------------------------------------------------------------- command line
def main(argv=None) -> int:
    """cenfuse control KIND IN OUT [--seed N] [--fai FAI] [-v] (argv: the arguments after 'control')."""
    ap = argparse.ArgumentParser(
        prog="cenfuse control", formatter_class=argparse.RawDescriptionHelpFormatter,
        description="Negative controls: an input with the signal one line reads destroyed and the rest kept. Run\n"
                    "cenfuse on a copy of GENOME_DIR with the control in place of the input: no test should pass\n"
                    "and no call should change. Deterministic given the seed.",
        epilog="kinds:\n"
               "  hic-repair    re-pair read ends at random; each end keeps chromosome, position, strand and MAPQ, so\n"
               "                every bin keeps its coverage (the benchmark's NEG2). Holds all ends in memory\n"
               "                (about 40 bytes per pair: 2.8 GB and 5 min for 53 M pairs)\n"
               "  hic-relocate  move each read end to a uniform random position on its chromosome; chromosome pairs\n"
               "                are kept (NEG)\n"
               "  gc-shuffle    shuffle 100-kb blocks of the chromosomes among them, lengths kept; writes OUT's .fai\n"
               "  ltr-shuffle   move each LTR retrotransposon (with its LTRs and domains) to a uniform random\n"
               "                position on its sequence; needs --fai\n"
               "hic-repair and hic-relocate test their output for clumped contacts and stop if they find them.")
    ap.add_argument("kind", choices=KINDS, metavar="{" + ",".join(KINDS) + "}", help="the control to make")
    ap.add_argument("inp", metavar="IN", help="4DN .pairs[.gz] (hic-*: #chromsize and #columns header, e.g. from "
                                              "pairtools); FASTA with .fai (gc-shuffle); GFF3[.gz] (ltr-shuffle)")
    ap.add_argument("out", metavar="OUT", help="the control, same format (.gz: gzip; a FASTA: bgzip and indexed)")
    ap.add_argument("--seed", type=int, default=1, metavar="N", help="random seed, >= 0 [1, the benchmark's]")
    ap.add_argument("--fai", help="ltr-shuffle only, required: the genome's .fai (sequence lengths)")
    ap.add_argument("-v", "--verbose", action="store_true", help="progress per chunk and the check at every bin size")
    a = ap.parse_args(argv)
    if a.seed < 0:
        ap.error("--seed must be >= 0")
    if a.kind == "ltr-shuffle" and not a.fai:
        ap.error("ltr-shuffle needs --fai (the genome's .fai)")
    if a.kind != "ltr-shuffle" and a.fai:
        ap.error("--fai is used by ltr-shuffle only")
    setup_logging(a.verbose)
    try:
        if a.kind == "hic-repair":
            repair_pairs(a.inp, a.out, a.seed)
        elif a.kind == "hic-relocate":
            relocate_pairs(a.inp, a.out, a.seed)
        elif a.kind == "gc-shuffle":
            shuffle_gc_blocks(a.inp, a.out, a.seed)
        else:
            shuffle_ltr(a.inp, a.out, a.fai, a.seed)
    except CenfuseError as e:
        print(f"cenfuse: error: {e}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
