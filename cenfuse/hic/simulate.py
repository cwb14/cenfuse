"""Mis-mapping null: cis-only read pairs simulated from the assembly, mapped and paired like the real library.

Every inter-chromosomal pair that survives is a mapping artefact. Counts are scaled to the real library by its
cis-pair count and compared with the real centromere contacts of each `python -m cenfuse.hic call` core."""
from __future__ import annotations

import gzip
import logging
import shlex
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import pysam

from . import params
from .io import InputError, PairStats, check_readable, is_cool, iter_pairs, read_fai

log = logging.getLogger("cenfuse-hic")
COMP = np.frombuffer(bytes.maketrans(b"ACGTN", b"TGCAN"), dtype=np.uint8)      # lookup tables, not parameters
BASES = np.frombuffer(b"ACGT", dtype=np.uint8)


def read_lengths(fastq, n=params.SIM_FASTQ_READS):
    op = gzip.open if str(fastq).endswith(".gz") else open
    out = []
    with op(fastq, "rt") as f:
        for i, line in enumerate(f):
            if i % 4 == 1:
                out.append(min(len(line.strip()), params.SIM_MAXLEN))
                if len(out) >= n:
                    break
    L = np.array(out)
    return np.column_stack([L, np.roll(L, 1)])


def _reads(seq, pos, strand, ln, rng, err):
    off = np.arange(params.SIM_MAXLEN)
    start = np.where(strand == 0, pos, pos - ln + 1)
    m = seq[np.clip(start[:, None] + off[None, :], 0, len(seq) - 1)]
    for i in np.where(strand == 1)[0]:
        m[i, :ln[i]] = COMP[m[i, :ln[i]][::-1]]
    e = rng.random(m.shape) < err
    if e.any():
        m[e] = BASES[(np.searchsorted(BASES, m[e]) + rng.integers(1, 4, int(e.sum()))) % 4]
    return [m[i, :ln[i]].tobytes().decode() if ln[i] > 0 else "N" for i in range(len(pos))]


def simulate_pairs(fasta, chroms, n, lengths, out_prefix, err=params.SIM_ERR, seed=params.SEED):
    """n cis-only pairs over `chroms` (length-weighted): mate 1 uniform, mate 2 a log-uniform distance (1 kb ..
    chromosome length) away, random strands, mate lengths drawn jointly from `lengths` (rows len1, len2)."""
    rng = np.random.default_rng(seed)
    fa = pysam.FastaFile(str(fasta))
    names = []
    ml, chunk = params.SIM_MAXLEN, params.SIM_CHUNK
    for c, L in chroms.items():
        if L < params.SIM_MIN_CHROM:
            log.warning(f"simulate: skipping {c} ({L} bp < {params.SIM_MIN_CHROM} bp)")
        else:
            names.append(c)
    if not names:
        raise InputError("no chromosome long enough to simulate from")
    w = np.array([chroms[c] for c in names], float)
    per = rng.multinomial(n, w / w.sum())
    fq1, fq2, tr = Path(f"{out_prefix}_1.fq"), Path(f"{out_prefix}_2.fq"), Path(f"{out_prefix}.truth.tsv")
    k0 = 0
    with open(fq1, "w") as f1, open(fq2, "w") as f2, open(tr, "w") as ft:
        ft.write("k\tchrom\tpos1\tstrand1\tpos2\tstrand2\n")
        for c, m in zip(names, per):
            if m == 0:
                continue
            seq = np.frombuffer(fa.fetch(c).upper().encode(), dtype=np.uint8)
            L = len(seq)
            for done in range(0, int(m), chunk):
                mm = min(chunk, int(m) - done)
                p1 = rng.integers(ml, L - ml, mm)
                d = (10 ** rng.uniform(np.log10(params.MIN_CIS_DIST), np.log10(L), mm)).astype(np.int64)
                s = np.where(rng.random(mm) < 0.5, -1, 1)
                p2 = p1 + s * d
                flip = (p2 < ml) | (p2 >= L - ml)
                p2[flip] = p1[flip] - s[flip] * d[flip]
                ok = (p2 >= ml) & (p2 < L - ml)
                p1, p2 = p1[ok], p2[ok]
                ln = lengths[rng.integers(0, len(lengths), len(p1))]
                s1, s2 = rng.integers(0, 2, len(p1)), rng.integers(0, 2, len(p1))
                r1, r2 = _reads(seq, p1, s1, ln[:, 0], rng, err), _reads(seq, p2, s2, ln[:, 1], rng, err)
                ks = k0 + np.arange(len(p1))
                f1.write("".join(f"@sim.{k}\n{a}\n+\n{'I' * len(a)}\n" for k, a in zip(ks, r1)))
                f2.write("".join(f"@sim.{k}\n{b}\n+\n{'I' * len(b)}\n" for k, b in zip(ks, r2)))
                ft.write("".join(f"{k}\t{c}\t{a}\t{'+-'[x]}\t{b}\t{'+-'[y]}\n"
                                 for k, a, x, b, y in zip(ks, p1, s1, p2, s2)))
                k0 += len(p1)
    fa.close()
    return fq1, fq2, tr


def best_segments(path):
    """5'-most segment per read (as `bwa mem -5` does); position = 5' end. Input: minimap2 PAF columns
    1,3,5,6,8,9,12 (qname, qstart, strand, tname, tstart, tend, mapq)."""
    cols = ["k", "chrom", "pos", "st", "mq"]
    try:
        t = pd.read_csv(path, sep="\t", header=None, names=["q", "qs", "st", "chrom", "ts", "te", "mq"],
                        dtype={"q": str, "chrom": str})
    except pd.errors.EmptyDataError:
        return pd.DataFrame({"k": pd.Series(dtype=np.int64), "chrom": pd.Series(dtype=str),
                             "pos": pd.Series(dtype=np.int64), "st": pd.Series(dtype=str),
                             "mq": pd.Series(dtype=np.int64)})[cols]
    t["k"] = t.q.str.rsplit(".", n=1).str[1].astype(np.int64)
    t = t.sort_values(["k", "qs"], kind="stable").drop_duplicates("k")
    t["pos"] = np.where(t.st == "+", t.ts, t.te - 1)
    return t[["k", "chrom", "pos", "st", "mq"]]


def pair_mates(m1, m2, q=params.MAPQ):
    """Both mates MAPQ >= q; cis pairs < MIN_CIS_DIST (1 kb) apart dropped; exact duplicates dropped (the
    usual Hi-C pair filters)."""
    p = m1.merge(m2, on="k", suffixes=("1", "2"))
    p = p[(p.mq1 >= q) & (p.mq2 >= q)]
    p = p[~((p.chrom1 == p.chrom2) & ((p.pos1 - p.pos2).abs() < params.MIN_CIS_DIST))]
    a = p[["chrom1", "pos1", "st1"]].astype(str).agg(":".join, axis=1)
    b = p[["chrom2", "pos2", "st2"]].astype(str).agg(":".join, axis=1)
    key = np.where(a <= b, a + "|" + b, b + "|" + a)
    return p.loc[~pd.Series(key, index=p.index).duplicated()].reset_index(drop=True)


def index_bases(total: int) -> int:
    """minimap2 -I value that keeps the whole genome in one index part: whole Gb, strictly above `total` bp.
    minimap2 splits its index every -I bases (default 8 Gb) and computes MAPQ within each part, so on a larger
    genome a read with copies in two parts would get a falsely high MAPQ."""
    return -(-(int(total) + 1) // 10**9) * 10**9


def _run(cmd, logf):
    r = subprocess.run(["bash", "-c", "set -o pipefail; " + cmd])
    if r.returncode != 0:
        tail = "\n".join(logf.read_text().splitlines()[-5:]) if logf.exists() else ""
        raise InputError(f"minimap2 failed (exit {r.returncode}); log {logf}:\n{tail}")


def _build_index(fasta, total, threads, out):
    """One minimap2 -x sr index of the whole FASTA (`total` bp) in a single part; built once, used for both mates."""
    out = Path(out)
    logf = Path(f"{out}.log")
    q = shlex.quote
    _run(f"minimap2 -x sr -I {index_bases(total)} -t {int(threads)} -d {q(str(out))} {q(str(fasta))} "
         f"2> {q(str(logf))}", logf)
    return out


def _map(index, fq, threads, out):
    out = Path(out)
    logf = Path(f"{out}.log")
    q = shlex.quote
    _run(f"minimap2 -cx sr --secondary=no -t {int(threads)} {q(str(index))} {q(str(fq))} 2> {q(str(logf))} "
         f"| cut -f1,3,5,6,8,9,12 > {q(str(out))}", logf)
    return out


def run_simulate(args):
    check_readable([args.fasta, args.pairs, *([args.fastq] if args.fastq else [])])
    if is_cool(args.pairs):
        raise InputError(f"{args.pairs}: simulate needs the real pairs file (its cis-pair count scales the null), "
                         "not .cool/.mcool")
    if shutil.which("minimap2") is None:
        raise InputError("minimap2 not found on PATH (needed to map the simulated reads)")
    cen = Path(f"{args.o}.centromeres.tsv")
    if not cen.exists():
        raise InputError(f"{cen} not found: run `python -m cenfuse.hic call -o {args.o}` first")
    fai_path = Path(str(args.fasta) + ".fai")
    if not fai_path.exists():
        raise InputError(f"{fai_path} not found: index the FASTA (samtools faidx)")
    calls = pd.read_csv(cen, sep="\t")
    need = {"chrom", "bin_size", "peak_start", "core_start", "core_end"}
    if need - set(calls.columns):
        raise InputError(f"{cen}: missing columns {sorted(need - set(calls.columns))}")
    called = calls[calls.core_start.notna() & calls.peak_start.notna()]
    if called.empty:
        raise InputError(f"{cen}: no chromosome has a called centromere core")
    size = int(calls.bin_size.dropna().iloc[0])
    fai = read_fai(fai_path)
    seqs = {c: fai[c] for c in calls.chrom if c in fai}
    lengths = (read_lengths(args.fastq) if args.fastq
               else np.array([[params.SIM_DEFAULT_LEN, params.SIM_DEFAULT_LEN]]))
    work = Path(f"{args.o}.simulate")
    work.mkdir(parents=True, exist_ok=True)
    log.info(f"simulate: {args.n} cis-only pairs")
    fq1, fq2, truth = simulate_pairs(args.fasta, seqs, args.n, lengths, work / "sim")
    total = sum(fai.values())
    log.info(f"simulate: minimap2 -x sr index of the whole FASTA ({total / 1e9:.2f} Gb, one part)")
    idx = _build_index(args.fasta, total, args.t, work / "genome.mmi")
    m1 = best_segments(_map(idx, fq1, args.t, work / "m1.tsv"))
    m2 = best_segments(_map(idx, fq2, args.t, work / "m2.tsv"))
    p = pair_mates(m1, m2).merge(pd.read_csv(truth, sep="\t").rename(columns={"pos1": "t1", "pos2": "t2"}), on="k")
    tol = params.SIM_TRUTH_TOL
    sim_cis_ok = int(((p.chrom1 == p.chrom2) & (p.chrom1 == p.chrom) & ((p.pos1 - p.t1).abs() <= tol)
                      & ((p.pos2 - p.t2).abs() <= tol)).sum())
    art = p[p.chrom1 != p.chrom2]
    st = PairStats()
    real = []
    win = {}
    for r in called.itertuples():
        c = int(r.peak_start) // size
        win[r.chrom] = ((c - params.HALF_WIN) * size, (c + params.HALF_WIN + 1) * size)
    core = {r.chrom: (int(r.core_start), int(r.core_end)) for r in called.itertuples()}

    def cen_counts(c1, p1, c2, p2):
        out = {}
        for ca, pa, cb, pb in ((c1, p1, c2, p2), (c2, p2, c1, p1)):
            for chrom, (s, e) in core.items():
                m = (ca == chrom) & (pa >= s) & (pa < e) & (cb != chrom)
                hit = np.zeros(int(m.sum()), bool)
                for oc, (ws, we) in win.items():
                    if oc != chrom:
                        hit |= (cb[m] == oc) & (pb[m] >= ws) & (pb[m] < we)
                out[chrom] = out.get(chrom, 0) + int(hit.sum())
        return out

    names = np.array(list(seqs))
    for ci1, q1, ci2, q2 in iter_pairs(args.pairs, seqs, st):
        real.append(cen_counts(names[ci1], q1, names[ci2], q2))
    real_tot = {c: sum(d.get(c, 0) for d in real) for c in core}
    scale = st.cis / max(sim_cis_ok, 1)
    a = cen_counts(art.chrom1.to_numpy(), art.pos1.to_numpy(), art.chrom2.to_numpy(), art.pos2.to_numpy())
    rows = [dict(chrom=c, real_cen_contacts=real_tot[c], expected_artefact_cen_contacts=a.get(c, 0) * scale,
                 artefact_share=a.get(c, 0) * scale / max(real_tot[c], 1)) for c in core]
    t = pd.DataFrame(rows)
    with open(f"{args.o}.simulate.tsv", "w") as f:
        f.write(f"# simulated {args.n} pairs; kept {len(p)}; artefact trans {len(art)} "
                f"({len(art) / max(len(p), 1):.4%}); scale to real {scale:.3g}; read lengths "
                f"{'from ' + str(args.fastq) if args.fastq else f'{params.SIM_DEFAULT_LEN} bp (assumed)'}, "
                f"capped at {params.SIM_MAXLEN} bp\n"
                "# valid only if the real library was mapped with a similar short-read mapper "
                "(minimap2 -x sr or bwa mem) and filtered to both-ends MAPQ >= 30\n")
        t.to_csv(f, sep="\t", index=False, float_format="%.4g")
    for f in (fq1, fq2, truth, work / "m1.tsv", work / "m2.tsv", idx):      # large intermediates; the logs stay
        Path(f).unlink(missing_ok=True)
    log.info(f"simulate done: artefacts {len(art)} of {len(p)} kept pairs; per-core shares in {args.o}.simulate.tsv")
