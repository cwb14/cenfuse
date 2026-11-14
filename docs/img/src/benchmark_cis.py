#!/usr/bin/env python3
"""benchmark_cis: draws docs/img/benchmark_cis.pdf and .png (docs/benchmark.md, cenfuse 0.16.0).

The two flanks of a centromere contact each other, so a cis Hi-C signal marks centromeres that no repeat marks. This script
records how the figure was made, legend text included. It CANNOT rebuild the figure from a clone of the repository alone:
it reads the cenfuse 0.16.0 benchmark runs and tables, which are not in the repository (too large; some truth unpublished).

Data read, under the benchmark data root $CENFUSE_BENCH_DATA (default: the folder they were made in, as in
benchmark/paths.py), with V11 = v2/v3/v4/v5/v6/v7/v8/v9/v10/v10.1/v10.2/v10.3/v10.4/v11 and V16 = V11/v12/v13/v14/v15/v16:
- V16/runs/rel16/TH/{tomato,potato}: hic/<g>.cache/*.npz (contact matrix) and work/hic/<g>.pkl (bins, the engine's cis result)
- V16/runs/rel16/TNEG/<g>, 10 genomes: the same for the randomized-Hi-C partners, and cenfuse_summary.tsv (cis_line)
- V16/runs/rel16/T/{tomato,potato}/<g>.te_signature.tsv (LTR cores, signature families); each run's input folder from
  V16/bench/jobs_rel16.tsv (tes/<g>.gff3 elements, the FASTA index) and tomato's satellite folder (satellite arrays, families)
- V16/results/docs16_rel16/release16.{runs,chrom}.tsv: calls and outcomes of arms T (no Hi-C) and TH (with Hi-C)
- panels E, F: V11/results/cis_survey4.tsv and the stored 0.11 runs it was computed on (V11/../runs/REL11/{SH,FULL}/<g>; the
  Hi-C engine is unchanged since 0.10.1); labels and lineages from V11/results/canonical.tsv
- truth: V11/bench/explib.py (the benchmark's frozen truth rules and files); code: cenfuse 0.16.0 (V16/code/rel16)

usage: CENFUSE_FIG_WORK=<work folder> python docs/img/src/benchmark_cis.py [--refresh] [-v]
Writes benchmark_cis.pdf (figure + one legend page) and benchmark_cis.png (figure only, 300 dpi) to $CENFUSE_FIG_OUT
(default: docs/img), and benchmark_cis_numbers.tsv (every number with its source) and a cache of the matrix-derived values
to the work folder ($CENFUSE_FIG_WORK, default: the current directory). docs/img's PNG is that PNG, losslessly recompressed.
"""
import argparse
import datetime
import os
import pickle
import socket
import sys
from pathlib import Path

os.environ.setdefault('OMP_NUM_THREADS', '1')
import numpy as np
import pandas as pd

SRC = Path(__file__).resolve().parent                                  # docs/img/src
OUT = Path(os.environ.get('CENFUSE_FIG_OUT') or SRC.parent)            # the PDF and the PNG (default: docs/img)
HERE = Path(os.environ.get('CENFUSE_FIG_WORK') or Path.cwd()).resolve()  # work folder: cache, numbers table
STEM = 'benchmark_cis'
DATA = Path(os.environ.get('CENFUSE_BENCH_DATA') or '/data2/chris/centromere_benchmark')  # as benchmark/paths.py
VV = DATA / 'v2/v3/v4/v5/v6/v7/v8/v9/v10/v10.1/v10.2/v10.3/v10.4/v11'
V16 = VV / 'v12/v13/v14/v15/v16'                               # the 0.16.0 release: runs, tables, code
REL = V16 / 'runs/rel16'
TAB = V16 / 'results/docs16_rel16'
CODE = V16 / 'code/rel16'                                       # identical to the repo's 0.16.0 commit
sys.path.insert(0, str(VV / 'bench'))
import explib as E                                            # frozen truth rules (score_rel), elements()
for _k in [k for k in sys.modules if k.startswith('cenfuse')]:
    del sys.modules[_k]                                         # use the frozen 0.16.0 release code, not explib's copy
sys.path.insert(0, str(CODE))
import scipy.sparse as sp
from cenfuse import hicfuse
from cenfuse.hic import params as hp
from cenfuse.hic.balance import balanced, ice
from cenfuse.hic.bins import Bins, coarsen
from cenfuse.util import poisson_binomial_sf
import cenfuse as _cf
assert Path(_cf.__file__).resolve().parent == (CODE / 'cenfuse').resolve(), _cf.__file__
assert _cf.__version__ == '0.16.0', _cf.__version__

HALF = hp.HALF_WIN          # 4 bins each side: the crossing box of cis_profile
TOLB = 4                    # bins allowed between an extreme and truth (as results/cis_survey4.tsv)
NSIM = 100                  # simulated maps per chromosome (as results/cis_survey4.tsv)
PILE = 12                   # pile-up half window, bins (3 Mb at 250 kb)
FHALF = 12                  # small multiples: bins each side of the truth center
TNEG = REL / 'TNEG'
SURVEY = VV / 'results/cis_survey4.tsv'           # survey of the stored 0.11 runs (engine unchanged since 0.10.1)
VERBOSE = False
CACHE = HERE / f'{STEM}_cache.pkl'


class Cache(dict):
    """A dict that remembers which keys this run asked for, so stale entries (old partner runs) can be dropped."""
    def __init__(self, *a):
        super().__init__(*a)
        self.used = set()

    def __contains__(self, k):
        self.used.add(k)
        return super().__contains__(k)

    def __getitem__(self, k):
        self.used.add(k)
        return super().__getitem__(k)

    def pruned(self):
        return Cache({k: v for k, v in self.items() if k in self.used})


def stash(cache):
    """Write the cache after every new entry (a partner run moved away mid-run then costs nothing already done)."""
    tmp = CACHE.with_suffix('.tmp')
    with open(tmp, 'wb') as f:
        pickle.dump(cache, f)
    os.replace(tmp, CACHE)


def log(*a, force=False):
    if VERBOSE or force:
        print(*a, file=sys.stderr, flush=True)


# --------------------------------------------------------------------------------------------- loading and statistics
def sig(p):
    p = Path(p)
    s = p.stat()
    return (str(p.resolve()), s.st_size, int(s.st_mtime))


def matrix_file(run, gid):
    h = pickle.load(open(run / 'work/hic' / f'{gid}.pkl', 'rb'))
    return h, run / 'hic' / f'{gid}.cache' / Path(h['matrix_cache']).name


def load(run, gid):
    """The engine's bins, balanced matrix and mask for one run (as bench/cis_plot.py)."""
    h, f = matrix_file(run, gid)
    b = h['bins']
    seqs = {c: int(L) for c, L in zip(b.names, b.lengths)}
    M0 = sp.load_npz(f).tocsr()
    base = Bins(seqs, int(h['base_size']))
    M = coarsen(M0, base, b) if b.size != base.size else M0
    bias, good, _, _ = ice(M)
    W = balanced(M, bias, good).tocsr()
    return dict(h=h, b=b, M=M, M0=M0, base=base, bias=bias, good=good, W=W, file=f)


def sim_profiles(W, b, good, bias, k, n=NSIM, half=HALF):
    """Profiles of n contact maps drawn as in hicfuse.cis_chance (Poisson counts from each bin's coverage and the
    chromosome's mean contact at each distance; same random stream), so the band is that null drawn."""
    nb = int(b.nb[k])
    bb = np.asarray(bias)[b.sl(k)]
    _, ok, mean, expected = hicfuse._cis_parts(W, b, good, k, half)
    use = expected > 0
    weight, counts = {}, {}
    for d in ok:
        with np.errstate(divide='ignore', invalid='ignore'):
            weight[d] = np.where(ok[d], 1.0 / (bb[:nb - d] * bb[d:]), 0.0)
            counts[d] = np.where(weight[d] > 0, mean[d] / weight[d], 0.0)
    rng = np.random.default_rng([hp.SEED, k])
    S = np.full((n, nb), np.nan)
    for i in range(n):
        S[i, use] = hicfuse._crossing({d: rng.poisson(counts[d]) * weight[d] for d in ok}, nb, half)[use] / expected[use]
    return S


def oe_dense(W, b, good, k, maxd):
    """Observed/expected contact of chromosome k (expected: the chromosome's mean balanced contact at each distance)."""
    sl = b.sl(k)
    A = W[sl, :][:, sl].toarray()
    g = good[sl]
    nb = len(g)
    ok = np.outer(g, g)
    X = np.full((nb, nb), np.nan)
    for d in range(1, min(nb, maxd + 1)):
        x, m = np.diag(A, d), np.diag(ok, d)
        e = x[m].mean() if m.any() else np.nan
        if not e > 0:
            continue
        i = np.arange(nb - d)
        X[i, i + d] = np.where(m, x / e, np.nan)
        X[i + d, i] = X[i, i + d]
    return X


def truth_zone(tr, c, b, k, tolb=TOLB):
    s0, e0 = min(x[0] for x in tr[c]), max(x[1] for x in tr[c])
    n = int(b.nb[k])
    return s0, e0, max(0, (s0 - tolb * b.size) // b.size), min(n - 1, (e0 + tolb * b.size) // b.size)


def extremes_on_truth(L, gid, chance=False):
    """Per chromosome with truth: the cis profile's peak and dip, whether each lies on truth (TOLB bins), and (chance=True)
    the simulated chance of that, as bench/cis_survey4.py."""
    b, good, W, bias, h = L['b'], L['good'], L['W'], L['bias'], L['h']
    tr, _ = E.truth(gid, b.size)
    idx = {n: i for i, n in enumerate(b.names)}
    rows = []
    for c in tr:
        k = idx.get(c)
        if k is None or not h['members'][k]:
            continue
        v = hicfuse.cis_profile(W, b, good, k)
        if not np.isfinite(v).any():
            continue
        s0, e0, lo, hi = truth_zone(tr, c, b, k)
        jp, jm = int(np.nanargmax(v)), int(np.nanargmin(v))
        r = dict(chrom=c, peak=lo <= jp <= hi, dip=lo <= jm <= hi, peak_val=float(v[jp]), dip_val=float(v[jm]))
        if chance:
            ok = np.isfinite(v)
            share = float(ok[lo:hi + 1].sum() / ok.sum())
            s = hicfuse.cis_chance(W, b, good, bias, k, lo, hi, n=NSIM)
            r.update(sp_=max(share, s[1]), sm_=max(share, s[-1]))
        rows.append(r)
    return pd.DataFrame(rows)


def dispersion(L):
    """Pearson dispersion of the partner's raw counts at distances of 1-7 bins (the engine's bins, so sparse maps such as
    barley's still have counts to test) against what a uniform shuffle of the same pairs gives (each pair's expected count
    r_i r_j / total): 1 for a uniform shuffle (Poisson counts), far above 1 when contacts are clumped."""
    M, b = L['M'], L['b']
    r = np.asarray(M.sum(axis=1)).ravel()
    T = r.sum()
    num, den = 0.0, 0
    for k in range(len(b.names)):
        sl, n = b.sl(k), int(b.nb[k])
        if n < 10:
            continue
        B = M[sl, :][:, sl]
        for d in range(1, 8):
            x = B.diagonal(d)
            e = r[sl][:n - d] * r[sl][d:] / T
            m = e > 0.5
            num += float(((x[m] - e[m]) ** 2 / e[m]).sum())
            den += int(m.sum())
    return num / den if den >= 20 else float('nan')


def partner_run(gid):
    """The randomized-Hi-C run of a genome: V16/runs/rel16/TNEG/<g>, once DONE."""
    run = TNEG / gid
    if (run / 'DONE').exists():
        return run
    raise SystemExit(f'randomized partner {run} has no DONE file; rerun this script when it is done')


# --------------------------------------------------------------------------------------------- the cached computations
def genome_views(gid, run, partner, cache):
    """Profiles, simulated bands and the example's O/E block for tomato / potato, real and partner."""
    key = ('views', gid, sig(matrix_file(run, gid)[1]), sig(matrix_file(partner, gid)[1]), NSIM)
    if key in cache:
        return cache[key]
    log(f'{gid}: profiles, simulated bands, pile-ups', force=True)
    R, N = load(run, gid), load(partner, gid)
    b = R['b']
    tr, _ = E.truth(gid, b.size)
    out = dict(bin=b.size, names=list(b.names), lengths=[int(x) for x in b.lengths], members=list(R['h']['members']),
               cis=R['h']['cis'], prof={}, good={}, band={}, part={}, part_bin=N['b'].size, simmax={})
    for k, c in enumerate(b.names):
        if not R['h']['members'][k]:
            continue
        out['prof'][c] = hicfuse.cis_profile(R['W'], b, R['good'], k)
        out['good'][c] = R['good'][b.sl(k)].copy()
        S = sim_profiles(R['W'], b, R['good'], R['bias'], k)
        out['band'][c] = np.nanquantile(S, [0.025, 0.975], axis=0)
        out['simmax'][c] = float(np.nanmedian(np.nanmax(S, axis=1)))
        kn = list(N['b'].names).index(c)
        out['part'][c] = hicfuse.cis_profile(N['W'], N['b'], N['good'], kn)
    # pile-ups of O/E around the truth center and around random positions of the same chromosomes (bench/cis_pileup.py)
    rng = np.random.default_rng(7)
    S_, N_, Sr, Nr = (np.zeros((2 * PILE + 1,) * 2) for _ in range(4))
    nchr = 0
    oe = {}
    for k, c in enumerate(b.names):
        if not R['h']['members'][k] or c not in tr:
            continue
        X = oe_dense(R['W'], b, R['good'], k, 2 * PILE + 2)
        oe[c] = X
        nb = X.shape[0]
        s0, e0 = min(x[0] for x in tr[c]), max(x[1] for x in tr[c])
        cb = ((s0 + e0) // 2) // b.size

        def add(cc, S, Nn):
            lo, hi = cc - PILE, cc + PILE + 1
            if lo < 0 or hi > nb:
                return
            blk = X[lo:hi, lo:hi]
            f = np.isfinite(blk)
            S[f] += blk[f]
            Nn[f] += 1
        add(cb, S_, N_)
        nchr += 1
        for _ in range(5):
            r = int(rng.integers(PILE, max(PILE + 1, nb - PILE)))
            if abs(r - cb) > 2 * PILE:
                add(r, Sr, Nr)
    P = np.where(N_ > 0, S_ / np.maximum(N_, 1), np.nan)
    Q = np.where(Nr > 0, Sr / np.maximum(Nr, 1), np.nan)
    box = lambda A: float(np.nanmean(A[PILE - HALF:PILE, PILE:PILE + HALF]))
    same = lambda A: float(np.nanmean(np.concatenate([A[PILE - HALF:PILE, PILE - HALF:PILE][np.triu_indices(HALF, 1)],
                                                        A[PILE:PILE + HALF, PILE:PILE + HALF][np.triu_indices(HALF, 1)]])))
    out['pile'] = dict(P=P, Q=Q, n=nchr, nrand=int(Nr[PILE, PILE + 1]), box=box(P), box_rand=box(Q), same=same(P),
                       same_rand=same(Q))
    out['oe'] = oe
    out['partner_dispersion'] = dispersion(N)
    out['ext_real'] = extremes_on_truth(R, gid)
    out['ext_part'] = extremes_on_truth(N, gid)
    cache[key] = out
    stash(cache)
    return out


def survey_profiles(cache):
    """Every genome of results/cis_survey4.tsv (the stored 0.11 release runs it used, REL11; not rerun for 0.16): each
    chromosome's profile around its truth center, in bins."""
    sv = pd.read_csv(SURVEY, sep='\t')
    out = {}
    for r in sv.itertuples():
        run = E.stored(r.arm, r.id)
        key = ('survey', r.arm, r.id, sig(matrix_file(run, r.id)[1]), FHALF)
        if key not in cache:
            log(f'survey: {r.arm}/{r.id}', force=True)
            L = load(run, r.id)
            b = L['b']
            tr, _ = E.truth(r.id, b.size)
            idx = {n: i for i, n in enumerate(b.names)}
            rows, widths = [], []
            for c in tr:
                k = idx.get(c)
                if k is None or not L['h']['members'][k]:
                    continue
                v = hicfuse.cis_profile(L['W'], b, L['good'], k)
                if not np.isfinite(v).any():
                    continue
                s0, e0 = min(x[0] for x in tr[c]), max(x[1] for x in tr[c])
                cb = ((s0 + e0) // 2) // b.size
                w = np.full(2 * FHALF + 1, np.nan)
                for j in range(-FHALF, FHALF + 1):
                    if 0 <= cb + j < len(v):
                        w[j + FHALF] = v[cb + j]
                rows.append(w)
                widths.append((e0 - s0) / b.size)
            cache[key] = dict(bin=b.size, traces=np.array(rows), width=float(np.median(widths)) if widths else 0.0)
            stash(cache)
        out[r.id] = cache[key]
    return sv, out


def partner_survey(sv, cache):
    """The same count on every randomized partner (V16/runs/rel16/TNEG), with its own simulated chance."""
    rows = []
    for g in sorted(p.name for p in TNEG.iterdir() if p.is_dir()):
        if g not in set(sv.id):
            continue
        run = partner_run(g)
        key = ('partner', g, sig(matrix_file(run, g)[1]), NSIM, TOLB, 'disp2')
        if key not in cache:
            log(f'partner: {g}', force=True)
            L = load(run, g)
            d = extremes_on_truth(L, g, chance=True)
            cache[key] = dict(n=len(d), peak=int(d.peak.sum()), dip=int(d.dip.sum()), chance_peak=float(d.sp_.mean()),
                              chance_dip=float(d.sm_.mean()),
                              p_peak=poisson_binomial_sf(int(d.peak.sum()), list(d.sp_)),
                              p_dip=poisson_binomial_sf(int(d.dip.sum()), list(d.sm_)), dispersion=dispersion(L),
                              bin=L['b'].size)
            stash(cache)
        s = pd.read_csv(run / 'cenfuse_summary.tsv', sep='\t', dtype=str, keep_default_na=False).iloc[0]
        rows.append(dict(id=g, **cache[key], cis_line=s['cis_line'], run=str(run)))
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------------------------- annotations
JOBS = V16 / 'bench/jobs_rel16.tsv'      # arm, genome, input dir, satellite dir, 0.15 run, group


def jobs():
    J = {}
    for l in open(JOBS):
        p = l.rstrip('\n').split('\t')
        J[(p[0], p[1])] = p
    return J


def satellites(gid, J):
    d = Path(J[('T', gid)][3])                 # the satellite folder the 0.16 T run used
    fam = pd.read_csv(d / 'satellite_families.tsv', sep='\t')
    cen = set(fam.loc[fam['class'] == 'centromeric', 'family'])
    rows = []
    for l in open(d / f'{gid}_satellite_arrays.gff'):
        if l[0] == '#':
            continue
        p = l.rstrip('\n').split('\t')
        a = dict(x.split('=', 1) for x in p[8].split(';') if '=' in x)
        rows.append((p[0], int(p[3]) - 1, int(p[4]), a.get('family', '.'), a.get('family', '.') in cen))
    return pd.DataFrame(rows, columns=['chrom', 's', 'e', 'family', 'cen']), sorted(cen), d


def ltr_signature(gid, J, chroms):
    t = pd.read_csv(REL / 'T' / gid / f'{gid}.te_signature.tsv', sep='\t')
    gff = Path(J[('T', gid)][2]) / 'tes' / f'{gid}.gff3'
    el = E.elements(gff, set(chroms))
    keep = []
    for r in t.itertuples():
        fams = set(str(r.families).split(',')) if isinstance(r.families, str) else set()
        keep.append(el[(el.chrom == r.chrom) & el.family.isin(fams)])
    return t.set_index('chrom'), pd.concat(keep), gff


def calls(path, arm, gid):
    t = pd.read_csv(path, sep='\t', dtype=str, keep_default_na=False)
    t = t[(t.arm == arm) & (t.genome == gid)]
    return {r.chrom: (int(r.start), int(r.end), r.ok, r.outcome, r.truth) for r in t.itertuples() if r.start not in ('', '.')}


def check_truth(gid, hic_bin):
    """The truth drawn (explib, frozen rules) must be the truth that scored the 0.16 release, and each call's outcome in
    release16.chrom.tsv must follow from it: gap (call as 1-based closed -> 0-based half-open) within the tolerance."""
    tr, tol = E.truth(gid, hic_bin)
    t = pd.read_csv(TAB / 'release16.chrom.tsv', sep='\t', dtype=str, keep_default_na=False)
    t = t[t.arm.isin(['T', 'TH']) & (t.genome == gid)]
    assert len(t), gid
    n = 0
    for r in t.itertuples():
        if r.truth in ('', '.'):
            continue
        c = r.chrom
        assert c in tr, (gid, c)
        s0, e0 = min(x[0] for x in tr[c]), max(x[1] for x in tr[c])
        assert r.truth in (f'{s0}-{e0}', ','.join(f'{a}-{b}' for a, b in tr[c])), (gid, c, r.truth, tr[c])
        assert int(float(r.tol)) == int(tol), (gid, r.tol, tol)
        if r.start in ('', '.'):
            assert r.outcome == 'none', (gid, c, r.outcome)
            continue
        s, e = int(r.start) - 1, int(r.end)
        g_ = min(max(0, max(s, a) - min(e, b)) for a, b in tr[c])
        assert abs(g_ - float(r.gap)) < 1, (gid, c, g_, r.gap)
        ok = 'R' if g_ <= tol else 'W'
        assert ok == r.ok, (gid, c, ok, r.ok)
        assert r.outcome.startswith('right') == (ok == 'R'), (gid, c, r.outcome)
        n += 1
    return n


def fasta_gaps(fa, chrom, s, e):
    """N-runs in fa[chrom][s:e] via the .fai (no full load)."""
    fai = {l.split('\t')[0]: l.split('\t') for l in open(str(fa) + '.fai')}
    _, L, off, lb, lB = fai[chrom][:5]
    L, off, lb, lB = int(L), int(off), int(lb), int(lB)
    s, e = max(0, s), min(L, e)
    with open(fa, 'rb') as f:
        f.seek(off + (s // lb) * lB + s % lb)
        raw = f.read((e - s) + (e - s) // lb * (lB - lb) + lB)
    seq = raw.replace(b'\n', b'').replace(b'\r', b'')[:e - s].upper()
    runs, i = [], seq.find(b'N')
    while i >= 0:
        j = i
        while j < len(seq) and seq[j:j + 1] == b'N':
            j += 1
        runs.append((s + i, s + j))
        i = seq.find(b'N', j)
    return runs


# --------------------------------------------------------------------------------------------- names
CAN = pd.read_csv(VV / 'results/canonical.tsv', sep='\t').set_index('genome')
SHORT = {'aetAL878': '$\\it{Ae.\\ tauschii}$', 'barley': 'barley', 'brachy': '$\\it{Brachypodium}$',
         'CC1': '$\\it{B.\\ rapa}$', 'colcc': '$\\it{A.\\ thaliana}$', 'gmX26': 'soybean', 'potato': 'potato',
         'ppat': 'moss', 'rice': 'rice', 'sbBTx623': 'sorghum', 'tomato': 'tomato',
         'scer': '$\\it{S.\\ cerevisiae}$', 'spom': '$\\it{S.\\ pombe}$', 'ncra': '$\\it{N.\\ crassa}$',
         'cneoA': '$\\it{C.\\ neoformans}$ A', 'cneoG1': '$\\it{C.\\ neoformans}$ G1', 'ctro': '$\\it{C.\\ tropicalis}$',
         'msym': '$\\it{M.\\ sympodialis}$', 'vdah': '$\\it{V.\\ dahliae}$', 'ztri': '$\\it{Z.\\ tritici}$',
         'pfal': '$\\it{P.\\ falciparum}$', 'ggswu': 'chicken', 'chm13': 'human', 'mPanTro3': 'chimpanzee',
         'mPanPan1': 'bonobo', 'mGorGor1': 'gorilla', 'mPonAbe1': 'Sumatran orang.', 'mPonPyg2': 'Bornean orang.',
         'mSymSyn1': 'siamang', 'MFA8': 'crab-eating mac.', 'MMU8': 'rhesus macaque', 'mmT2T': 'mouse',
         'horse': 'horse EquCab3', 'horse1': 'horse T2T', 'goat1': 'goat', 'sheep1': 'sheep', 'mGirTip1': 'giraffe',
         'Sscrofa': 'pig'}
FSHORT = {'cneoA': '$\\it{C.\\ neoformans}$', 'cneoG1': '$\\it{C.\\ neoformans}$', 'mPonAbe1': 'Sum. orangutan',
          'mPonPyg2': 'Bor. orangutan', 'MFA8': 'cynomolgus', 'MMU8': 'rhesus'}
FNOTE = {'cneoA': 'anaphase', 'cneoG1': 'G1'}
PLAIN = {k: v.replace('$\\it{', '').replace('}$', '').replace('\\ ', ' ') for k, v in SHORT.items()}


def lineage(g):
    s = str(CAN.at[g, 'lineage']) if g in CAN.index else '?'
    return 'mammals' if s.startswith('mammals') else 'fungi' if s.startswith('fungi') else s


# --------------------------------------------------------------------------------------------- figure
PURPLE, GRAY, BLUE, ORANGE, VERM = '#6A3D9A', '#8C8C8C', '#A6CEE3', '#E69F00', '#D55E00'
SATGRAY = '#B0B0B0'
YL = (0.25, 2.25)                      # every profile track: one log scale
FW, FH = 180.0, 236.0


def build(args):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    from matplotlib.colors import LogNorm
    from matplotlib.patches import Polygon, Rectangle
    from matplotlib.lines import Line2D
    from matplotlib import ticker
    plt.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'Liberation Sans'], 'font.size': 6.5,
                         'pdf.fonttype': 42, 'ps.fonttype': 42, 'axes.linewidth': 0.5, 'xtick.major.width': 0.5,
                         'ytick.major.width': 0.5, 'xtick.major.size': 2, 'ytick.major.size': 2, 'xtick.labelsize': 6,
                         'ytick.labelsize': 6, 'axes.labelsize': 6.5, 'mathtext.default': 'regular',
                         'hatch.linewidth': 0.4, 'hatch.color': '#BBBBBB', 'mathtext.fontset': 'custom',
                         'mathtext.rm': 'Liberation Sans', 'mathtext.it': 'Liberation Sans:italic',
                         'mathtext.bf': 'Liberation Sans:bold', 'mathtext.sf': 'Liberation Sans',
                         'mathtext.cal': 'Liberation Sans', 'mathtext.tt': 'Liberation Mono'})
    cache = Cache() if args.refresh or not CACHE.exists() else Cache(pickle.load(open(CACHE, 'rb')))
    J = jobs()
    runs = {g: REL / 'TH' / g for g in ('tomato', 'potato')}
    part = {g: partner_run(g) for g in ('tomato', 'potato')}
    G = {g: genome_views(g, runs[g], part[g], cache) for g in ('tomato', 'potato')}
    sv, SV = survey_profiles(cache)
    PS = partner_survey(sv, cache)
    stash(dict(cache.pruned()))
    NUM = []

    def num(panel, what, value, source, release='0.16.0'):
        NUM.append(dict(panel=panel, quantity=what, value=value, source=source, release=release))
        return value

    # ---------------- tomato annotations
    T = G['tomato']
    tchroms = [c for c, m in zip(T['names'], T['members']) if m]
    tlen = dict(zip(T['names'], T['lengths']))
    trT, _ = E.truth('tomato', T['bin'])
    sigT, elT, gffT = ltr_signature('tomato', J, tchroms)
    satT, cenfam, satdir = satellites('tomato', J)
    cHi = calls(TAB / 'release16.chrom.tsv', 'TH', 'tomato')     # 0.16 call with Hi-C (lower row)
    cNo = calls(TAB / 'release16.chrom.tsv', 'T', 'tomato')      # 0.16 call without Hi-C: FASTA + tes/ (upper row)
    chrname = {c: f'chr{i + 1}' for i, c in enumerate(tchroms)}
    P = G['potato']
    pchroms = [c for c, m in zip(P['names'], P['members']) if m]
    plen = dict(zip(P['names'], P['lengths']))
    trP, _ = E.truth('potato', P['bin'])
    sigP = pd.read_csv(REL / 'T/potato/potato.te_signature.tsv', sep='\t').set_index('chrom')
    nchk = {g: check_truth(g, T['bin'] if g == 'tomato' else P['bin']) for g in ('tomato', 'potato')}
    log(f'truth and outcomes agree with release16.chrom.tsv: {nchk}', force=True)

    # ---------------- the example (panel B): rule from the spec
    at_core = T['cis']['at_core']
    cand = []
    for c in tchroms:
        v = T['prof'][c]
        e = T['ext_real'].set_index('chrom').loc[c]
        if e.peak and at_core.get(c):
            s0, e0 = min(x[0] for x in trT[c]), max(x[1] for x in trT[c])
            cand.append(dict(chrom=c, name=chrname[c], height=float(np.nanmax(v)), width=(e0 - s0) / 1e6))
    cand = pd.DataFrame(cand)
    for col in ('height', 'width'):
        med = cand[col].median()
        mad = (cand[col] - med).abs().median()
        cand[f'z_{col}'] = (cand[col] - med) / mad
    cand['dist'] = np.hypot(cand.z_height, cand.z_width)
    cand = cand.sort_values('dist')
    fa = Path(J[('T', 'tomato')][2]) / 'tomato.fa'

    def structure(c):
        """Masked bins within 2 Mb (8 bins) of the peak bin, and assembly gaps (N runs) within 2 Mb of it."""
        v_, g_ = T['prof'][c], T['good'][c]
        j_ = int(np.nanargmax(v_))
        near = [j for j in range(max(0, j_ - 2 * HALF), min(len(g_), j_ + 2 * HALF + 1)) if not g_[j]]
        gp = fasta_gaps(fa, c, j_ * T['bin'] - 2_000_000, (j_ + 1) * T['bin'] + 2_000_000)
        return len(near), len(gp)
    cand['masked'], cand['gaps'] = zip(*[structure(c) for c in cand.chrom])
    clean_ = cand[(cand.masked == 0) & (cand.gaps == 0)]
    EX = clean_.iloc[0].chrom                                    # nearest the medians among the structurally clean
    skipped = cand[cand.dist < clean_.iloc[0].dist]
    log(cand.to_string(index=False), force=True)
    s0x, e0x = min(x[0] for x in trT[EX]), max(x[1] for x in trT[EX])
    vx = T['prof'][EX]
    jpx = int(np.nanargmax(vx))
    masked_near, ngaps = structure(EX)
    rel_truth = cHi[EX][4]
    assert rel_truth == f'{s0x}-{e0x}', (rel_truth, s0x, e0x)        # the band drawn is the truth that scored the call
    assert cHi[EX][2] == 'R' and cNo[EX][2] == 'R', (cHi[EX], cNo[EX])  # the example shows no miss (skill rule)

    # =========================================================================== page 1: the figure
    fig = plt.figure(figsize=(FW / 25.4, FH / 25.4))

    def ax_mm(x, y, w, h, **kw):
        return fig.add_axes([x / FW, 1 - (y + h) / FH, w / FW, h / FH], **kw)

    def letter(x, y, s):
        fig.text(x / FW, 1 - y / FH, s, fontsize=9, fontweight='bold', va='top', ha='left')

    def clean(a, left=False, bottom=False):
        for s in ('top', 'right'):
            a.spines[s].set_visible(False)
        a.spines['left'].set_visible(left)
        a.spines['bottom'].set_visible(bottom)
        if not left:
            a.set_yticks([])
        if not bottom:
            a.set_xticks([])
        a.patch.set_alpha(0)

    def logy(a, labels='ends'):
        a.set_yscale('log')
        a.set_ylim(*YL)
        a.yaxis.set_major_locator(ticker.FixedLocator([0.5, 1, 2]))
        a.yaxis.set_minor_locator(ticker.NullLocator())
        a.set_yticklabels({'ends': ['½×', '', '2×'], 'all': ['½×', '1×', '2×']}.get(labels, ['', '', '']))

    def hatch_missing(a, good, binmb, v=None):
        bad = ~good if v is None else (~good | ~np.isfinite(v))
        j = 0
        while j < len(bad):
            if bad[j]:
                i = j
                while j < len(bad) and bad[j]:
                    j += 1
                a.add_patch(Rectangle((i * binmb, YL[0]), (j - i) * binmb, YL[1] - YL[0], hatch='//////',
                                      fill=False, lw=0, zorder=0.5))
            j += 1

    def profile_track(a, c, D, L, labels='none', partner=True, mark=+1):
        """mark +1: the highest bin of each profile marked above the track (real filled, partner open); -1: the lowest,
        below."""
        binmb = D['bin'] / 1e6
        v = D['prof'][c]
        x = (np.arange(len(v)) + 0.5) * binmb
        lo, hi = D['band'][c]
        a.fill_between(x, lo, hi, color='#D9D9D9', lw=0, zorder=1)
        a.axhline(1, color='#9E9E9E', lw=0.35, zorder=1.2)
        ym, mk = (YL[1] * 1.25, 'v') if mark > 0 else (YL[0] / 1.25, '^')
        pick = np.nanargmax if mark > 0 else np.nanargmin
        if partner:
            w = D['part'][c]
            xp = (np.arange(len(w)) + 0.5) * D['part_bin'] / 1e6
            a.plot(xp, np.where(w > 0, w, YL[0] * 0.5), color=GRAY, lw=0.45, zorder=2, solid_joinstyle='round')
            a.plot(xp[pick(w)], ym, marker=mk, ms=2.5, mfc='white', mec=GRAY, mew=0.5, clip_on=False, zorder=6)
        a.plot(x, v, color=PURPLE, lw=0.7, zorder=3, solid_joinstyle='round')
        a.plot(x[pick(v)], ym, marker=mk, ms=2.5, mfc=PURPLE, mec=PURPLE, mew=0.3, clip_on=False, zorder=7)
        hatch_missing(a, D['good'][c], binmb, v)
        a.set_xlim(0, L / 1e6)
        logy(a, labels)

    def band(a, ivs, alpha=0.75):
        for (bs, be) in ivs:
            a.axvspan(bs / 1e6, be / 1e6, color=BLUE, lw=0, zorder=0, alpha=alpha)

    def ltr_track(a, c, el, sigt, lw=0.35):
        e = el[el.chrom == c]
        a.vlines(e.mid / 1e6, 0.38, 1.0, color=ORANGE, lw=lw)
        r = sigt.loc[c]
        a.add_patch(Rectangle((r.core_start / 1e6, 0.0), (r.core_end - r.core_start) / 1e6, 0.26, color=ORANGE, lw=0))

    def sat_track(a, c, sat, lw=0.25):
        s_ = sat[sat.chrom == c]
        a.axhline(0.5, color='#CCCCCC', lw=0.3, zorder=1)
        for rr in s_.itertuples():
            colr = VERM if rr.cen else SATGRAY
            a.add_patch(Rectangle((rr.s / 1e6, 0.1), (rr.e - rr.s) / 1e6, 0.8, fc=colr, ec=colr, lw=lw,
                                  zorder=3 if rr.cen else 2))

    def call_box(a, c, cc, colr, y0, h):
        if c in cc:
            a.add_patch(Rectangle((cc[c][0] / 1e6, y0), (cc[c][1] - cc[c][0]) / 1e6, h, fc=colr, ec=colr, lw=0.4))

    # ---------------------------------------------------------------- A: tomato, 12 chromosomes, one bp scale
    MBMM = 75.0 / 94.0                      # mm per Mb in A and D
    PT, LT, ST, CT, GAP, BG = 4.6, 1.1, 1.1, 1.4, 0.25, 1.45
    BLOCK = PT + GAP + LT + GAP + ST + GAP + CT + BG
    yA = 9.6
    letter(0, 0.3, 'A')
    xcol = (13.0, 105.0)
    key = [(Line2D([], [], color=PURPLE, lw=1), 'contact across bin, obs./exp.'),
           (Line2D([], [], color=GRAY, lw=0.8), 'randomized Hi-C, same axis'),
           (Rectangle((0, 0), 1, 1, fc='#D9D9D9', lw=0), 'simulated maps, 95% per bin'),
           (Line2D([], [], color=PURPLE, lw=0, marker='v', ms=3.2), 'extreme bin (open: randomized)'),
           (Rectangle((0, 0), 1, 1, fc=BLUE, lw=0), 'CENH3 (ChIP)'),
           (Line2D([], [], color=ORANGE, lw=0, marker='|', ms=5, mew=0.8), 'LTR signature element'),
           (Line2D([], [], color=ORANGE, lw=2.2), 'LTR core'),
           (Rectangle((0, 0), 1, 1, fill=False, hatch='//////', lw=0.3, ec='#BBBBBB'), 'no data'),
           (Rectangle((0, 0), 1, 1, fc=VERM, lw=0), 'centromeric satellite family'),
           (Rectangle((0, 0), 1, 1, fc=SATGRAY, lw=0), 'other satellite arrays'),
           (Rectangle((0, 0), 1, 1, fc='#7F7F7F', lw=0), 'call, no Hi-C (upper)'),
           (Rectangle((0, 0), 1, 1, fc='black', lw=0), 'call, with Hi-C (lower)')]
    kax = ax_mm(6.0, 0.6, 174, 8.0)
    kax.axis('off')
    kax.legend([k[0] for k in key], [k[1] for k in key], ncol=4, loc='upper left', frameon=False, fontsize=6,
               handlelength=1.5, handleheight=0.7, columnspacing=1.2, labelspacing=0.22, borderaxespad=0,
               handletextpad=0.4)
    zoomA = None
    for i, c in enumerate(tchroms):
        col, row = divmod(i, 6)
        x0 = xcol[col]
        y0 = yA + row * BLOCK
        L = tlen[c]
        w = L / 1e6 * MBMM
        a = ax_mm(x0, y0, w, PT)
        profile_track(a, c, T, L, labels='ends' if row == 0 else 'none')
        clean(a, left=True)
        a.spines['left'].set_linewidth(0.4)
        a.tick_params(axis='y', length=1.5, pad=0.8, labelsize=5.5)
        fig.text((x0 - 4.6) / FW, 1 - (y0 + PT / 2) / FH, chrname[c], ha='right', va='center', fontsize=6.5)
        yy = y0 + PT + GAP
        al = ax_mm(x0, yy, w, LT)
        ltr_track(al, c, elT, sigT)
        yy += LT + GAP
        asat = ax_mm(x0, yy, w, ST)
        sat_track(asat, c, satT)
        yy += ST + GAP
        ac = ax_mm(x0, yy, w, CT)
        call_box(ac, c, cNo, '#7F7F7F', 0.55, 0.42)
        call_box(ac, c, cHi, 'black', 0.0, 0.42)
        for aa in (al, asat, ac):
            aa.set_xlim(0, L / 1e6)
            aa.set_ylim(0, 1)
            clean(aa)
        for aa in (a, al, asat, ac):
            band(aa, trT[c])
        if c == EX:
            zoomA = (x0, y0, w, L)
    for col in (0, 1):
        ax_ = ax_mm(xcol[col], yA + 6 * BLOCK - BG + 0.5, 94 * MBMM, 0.01)
        clean(ax_, bottom=True)
        ax_.set_xlim(0, 94)
        ax_.xaxis.set_major_locator(ticker.MultipleLocator(20))
        ax_.tick_params(axis='x', labelsize=6, pad=1)
        ax_.set_xlabel('tomato position (Mb)', fontsize=6.5, labelpad=1)
    yAend = yA + 6 * BLOCK + 5.5
    bmid = (s0x + e0x) / 2
    bwin = (max(0, int(bmid - 7_000_000)), min(tlen[EX], int(bmid + 7_000_000)))
    x0z, y0z, wz, Lz = zoomA
    zx0 = x0z + bwin[0] / 1e6 * MBMM
    zx1 = x0z + bwin[1] / 1e6 * MBMM
    fig.add_artist(Rectangle((zx0 / FW, 1 - (y0z + BLOCK - BG + 0.25) / FH), (zx1 - zx0) / FW,
                             (BLOCK - BG + 1.3) / FH, fill=False, lw=0.6, ec='black', transform=fig.transFigure))
    fig.text((zx1 + 0.5) / FW, 1 - (y0z - 0.2) / FH, 'B', fontsize=6.5, fontweight='bold', va='top')

    # ---------------------------------------------------------------- B: one chromosome as a contact map
    yB = yAend + 2.0
    letter(0, yB, 'B')
    bx0, bw = 13.0, 75.0
    mbmm_b = bw / ((bwin[1] - bwin[0]) / 1e6)
    binmb = T['bin'] / 1e6
    X = T['oe'][EX]
    i0, i1 = bwin[0] // T['bin'], -(-bwin[1] // T['bin'])
    maxd = 24                                     # 6 Mb
    tri_h = maxd / 2 * binmb * mbmm_b
    at = ax_mm(bx0, yB + 1.5, bw, tri_h)
    norm = LogNorm(vmin=0.5, vmax=2.0)
    cmap = plt.get_cmap('PRGn_r').copy()
    cmap.set_bad((1, 1, 1, 0))
    gX = T['good'][EX]
    n_ = i1 - i0
    I, Jg = np.meshgrid(np.arange(i0, i1 + 1), np.arange(i0, i1 + 1), indexing='ij')
    Xr, Yr = (I + Jg) / 2 * binmb, (Jg - I) / 2 * binmb
    C = np.full((n_, n_), np.nan)
    for a_ in range(n_):
        for b_ in range(a_ + 1, min(n_, a_ + maxd + 1)):
            C[a_, b_] = X[i0 + a_, i0 + b_]
    at.add_patch(Polygon([(i0 * binmb, 0), (i1 * binmb, 0), ((i1 - maxd / 2) * binmb, maxd / 2 * binmb),
                          ((i0 + maxd / 2) * binmb, maxd / 2 * binmb)], closed=True, fill=False, hatch='//////', lw=0))
    at.pcolormesh(Xr, Yr, np.ma.masked_invalid(C), cmap=cmap, norm=norm, shading='flat', rasterized=True,
                  edgecolors='none')
    jx = jpx
    dia = np.array([(jx - HALF / 2, HALF / 2), (jx, 0), (jx + HALF / 2, HALF / 2), (jx, HALF)]) * binmb
    at.add_patch(Polygon(dia, closed=True, fill=False, ec='black', lw=0.7, zorder=5))
    at.set_xlim(bwin[0] / 1e6, bwin[1] / 1e6)
    at.set_ylim(0, maxd / 2 * binmb)
    at.set_aspect('equal', adjustable='box', anchor='S')
    clean(at)
    at.spines['left'].set_visible(True)
    at.spines['left'].set_linewidth(0.4)
    at.set_yticks([0, 1, 2, 3])
    at.set_yticklabels(['0', '2', '4', '6'])
    at.tick_params(axis='y', labelsize=5.5, length=1.5, pad=1)
    at.set_ylabel('distance\n(Mb)', fontsize=6, labelpad=1)
    yb = yB + 1.5 + tri_h + 0.8
    PB = 8.0
    ap = ax_mm(bx0, yb, bw, PB)
    v = T['prof'][EX]
    lo_, hi_ = T['band'][EX]
    edges = np.arange(len(v) + 1) * binmb
    ap.fill_between(edges[:-1], lo_, hi_, step='post', color='#D9D9D9', lw=0, zorder=1)
    ap.axhline(1, color='#9E9E9E', lw=0.35, zorder=1.2)
    w_ = T['part'][EX]
    ap.step(np.arange(len(w_)) * T['part_bin'] / 1e6, np.where(w_ > 0, w_, YL[0] * 0.5), where='post', color=GRAY,
            lw=0.5, zorder=2)
    ap.step(edges[:-1], v, where='post', color=PURPLE, lw=0.8, zorder=3)
    hatch_missing(ap, gX, binmb, v)
    ap.plot([jx * binmb, jx * binmb], [v[jx], YL[1] * 3.5], color='black', lw=0.5, ls=(0, (2, 1.5)), clip_on=False,
            zorder=4)
    ap.set_xlim(bwin[0] / 1e6, bwin[1] / 1e6)
    logy(ap, 'all')
    clean(ap, left=True)
    ap.tick_params(axis='y', labelsize=5.5, length=1.5, pad=1)
    ap.set_ylabel('obs./exp.', fontsize=6, labelpad=1)
    yy = yb + PB + 0.5
    tracks = []
    for nm, h_ in (('LTR', 1.5), ('satellite', 1.5), ('no Hi-C', 1.0), ('with Hi-C', 1.0)):
        a_ = ax_mm(bx0, yy, bw, h_)
        tracks.append(a_)
        fig.text((bx0 - 1) / FW, 1 - (yy + h_ / 2) / FH, nm, ha='right', va='center', fontsize=5.5, color='#444444')
        yy += h_ + 0.55
    ltr_track(tracks[0], EX, elT, sigT, lw=0.5)
    sat_track(tracks[1], EX, satT, lw=0.3)
    call_box(tracks[2], EX, cNo, '#7F7F7F', 0.1, 0.8)
    call_box(tracks[3], EX, cHi, 'black', 0.1, 0.8)
    for a_ in tracks:
        a_.set_xlim(bwin[0] / 1e6, bwin[1] / 1e6)
        a_.set_ylim(0, 1)
        clean(a_)
    for a_ in [ap] + tracks:
        band(a_, trT[EX])
    axx = ax_mm(bx0, yy, bw, 0.01)
    clean(axx, bottom=True)
    axx.set_xlim(bwin[0] / 1e6, bwin[1] / 1e6)
    axx.xaxis.set_major_locator(ticker.MultipleLocator(2))
    axx.tick_params(axis='x', labelsize=6, pad=1)
    axx.set_xlabel(f'tomato {chrname[EX]} position (Mb)', fontsize=6.5, labelpad=1)
    yBend = yy + 5.5

    # ---------------------------------------------------------------- C: pile-ups, and what a fold or a domain would look like
    xC = 99.0
    letter(xC - 4.0, yB, 'C')
    tw = 25.0
    th = tw / 2

    def tri_map(a, A, binsz):
        n2 = A.shape[0]
        I2, J2 = np.meshgrid(np.arange(n2 + 1), np.arange(n2 + 1), indexing='ij')
        Xr2, Yr2 = (I2 + J2) / 2 - n2 / 2, (J2 - I2) / 2
        C2 = np.where(np.triu(np.ones_like(A, bool), 1), A, np.nan)
        a.add_patch(Polygon(np.array([(-n2 / 2, 0), (n2 / 2, 0), (0, n2 / 2)]) * binsz, closed=True, fill=False,
                            hatch='//////', lw=0))
        a.pcolormesh(Xr2 * binsz, Yr2 * binsz, np.ma.masked_invalid(C2), cmap=cmap, norm=norm, shading='flat',
                     rasterized=True, edgecolors='none')
        h2 = HALF                    # crossing box of the center bin: rows c-4..c-1, columns c..c+3
        d2 = np.array([(-h2 / 2 - 0.5, h2 / 2), (-0.5, 0), (h2 / 2 - 0.5, h2 / 2), (-0.5, h2)]) * binsz
        a.add_patch(Polygon(d2, closed=True, fill=False, ec='black', lw=0.5, zorder=5))
        a.set_xlim(-n2 / 2 * binsz, n2 / 2 * binsz)
        a.set_ylim(0, n2 / 2 * binsz)
        a.set_aspect('equal', adjustable='box', anchor='S')
        clean(a)
    ytri = yB + 4.2
    xt = (xC + 9.0, xC + 9.0 + tw + 2.5)
    for ri, g in enumerate(('tomato', 'potato')):
        D = G[g]
        for ci, (A_, t_) in enumerate(((D['pile']['P'], 'centered on CENH3'), (D['pile']['Q'], 'random positions'))):
            a = ax_mm(xt[ci], ytri + ri * (th + 2.5), tw, th)
            tri_map(a, A_, D['bin'] / 1e6)
            if ri == 0:
                fig.text((xt[ci] + tw / 2) / FW, 1 - (ytri - 0.5) / FH, t_, ha='center', va='bottom', fontsize=6)
            if ci == 0:
                fig.text((xt[0] - 1.0) / FW, 1 - (ytri + ri * (th + 2.5) + th * 0.55) / FH,
                         f"{g}\n{D['pile']['n']} chr.", ha='right', va='center', fontsize=6)
            if ri == 1:
                a.spines['bottom'].set_visible(True)
                a.spines['bottom'].set_linewidth(0.4)
                a.set_xticks([-3, 0, 3])
                a.set_xticklabels(['−3', '0', '3'])
                a.tick_params(axis='x', labelsize=5.5, length=1.5, pad=1)
                a.set_xlabel('Mb from the center', fontsize=6, labelpad=1)
    cby = ytri + 2 * th + 2.5 + 7.5
    cax = ax_mm(xt[0] + 6, cby, 2 * tw + 2.5 - 12, 1.4)
    cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax, orientation='horizontal')
    cb.set_ticks([0.5, 1, 2])
    cb.set_ticklabels(['½×', '1×', '2×'])
    cb.ax.tick_params(labelsize=5.5, length=1.5, pad=1, width=0.4)
    cb.outline.set_linewidth(0.4)
    cb.ax.xaxis.set_minor_locator(ticker.NullLocator())
    cb.set_label('contact, obs./exp.', fontsize=6, labelpad=1)
    # schematic column: flat cartoons, set apart
    xs0 = xt[1] + tw + 4.0
    sw = 180.0 - xs0 - 0.5
    fig.text((xs0 + sw / 2) / FW, 1 - (ytri - 0.5) / FH, 'expected shapes', ha='center', va='bottom', fontsize=6,
             color='#444444')
    for si, (nm, kind) in enumerate((('fold', 'fold'), ('insulated\ndomain', 'domain'))):
        ys_ = ytri + si * (th + 2.5)
        tws = sw * 0.62
        a = ax_mm(xs0, ys_ + th - tws / 2, tws, tws / 2)
        a.add_patch(Polygon([(-1, 0), (1, 0), (0, 1)], closed=True, fc='white', ec='#999999', lw=0.4))
        if kind == 'fold':
            a.add_patch(Polygon([(-0.08, 0.08), (0.08, 0.08), (0.08, 0.92), (-0.08, 0.92)], fc=PURPLE, lw=0))
        else:
            a.add_patch(Polygon([(-0.24, 0), (0.24, 0), (0, 0.24)], fc=PURPLE, lw=0))
            a.add_patch(Polygon([(-0.24, 0), (0, 0.24), (-0.38, 0.62), (-0.62, 0.38)], fc='#A6DBA0', lw=0))
            a.add_patch(Polygon([(0.24, 0), (0, 0.24), (0.38, 0.62), (0.62, 0.38)], fc='#A6DBA0', lw=0))
        a.set_xlim(-1, 1)
        a.set_ylim(0, 1)
        a.axis('off')
        b_ = ax_mm(xs0 + tws + 0.6, ys_ + th - tws / 2, sw - tws - 0.6, tws / 2)
        b_.axis('off')
        b_.set_xlim(0, 1)
        b_.set_ylim(0, 1)
        if kind == 'fold':
            tt = np.linspace(0, 1, 30)
            b_.plot(0.05 + 0.4 * tt, 0.95 - 0.8 * tt, color='black', lw=0.6)
            b_.plot(0.95 - 0.4 * tt, 0.95 - 0.8 * tt, color='black', lw=0.6)
            b_.plot([0.45, 0.55], [0.15, 0.15], color='black', lw=0.6)
            b_.plot([0.5], [0.15], 'o', color=BLUE, ms=2.6, mec='black', mew=0.3)
            for q in (0.3, 0.5, 0.7):
                b_.plot([0.05 + 0.4 * q, 0.95 - 0.4 * q], [0.95 - 0.8 * q] * 2, color=PURPLE, lw=0.5, ls=(0, (1, 1)))
        else:
            b_.plot([0.0, 0.3], [0.45, 0.45], color='black', lw=0.6)
            b_.plot([0.7, 1.0], [0.45, 0.45], color='black', lw=0.6)
            th_ = np.linspace(0, 2 * np.pi, 120)
            b_.plot(0.5 + 0.2 * np.cos(th_) * (1 + 0.15 * np.sin(5 * th_)), 0.45 + 0.4 * np.sin(th_) *
                    (1 + 0.15 * np.sin(5 * th_)), color='black', lw=0.6)
            b_.plot([0.5], [0.45], 'o', color=BLUE, ms=2.6, mec='black', mew=0.3)
        fig.text((xs0 + tws / 2) / FW, 1 - (ys_ + th - tws / 2 - 0.3) / FH, nm, ha='center', va='bottom',
                 fontsize=5.8)
    yCend = cby + 6.0

    # ---------------------------------------------------------------- D: potato, one row per chromosome, the same scales
    yD = max(yBend, yCend) + 2.0
    letter(0, yD, 'D')
    PD, GD = PT, 0.35
    xD = 13.0
    for i, c in enumerate(pchroms):
        y0 = yD + 1.2 + i * (PD + GD)
        L = plen[c]
        a = ax_mm(xD, y0, L / 1e6 * MBMM, PD)
        profile_track(a, c, P, L, labels='ends' if i == 0 else 'none', mark=-1)
        clean(a, left=True)
        a.spines['left'].set_linewidth(0.4)
        a.tick_params(axis='y', length=1.5, pad=0.8, labelsize=5.5)
        band(a, trP.get(c, []))
        if c in sigP.index:
            r = sigP.loc[c]
            a.add_patch(Rectangle((r.core_start / 1e6, 0), (r.core_end - r.core_start) / 1e6, 0.1, color=ORANGE, lw=0,
                                  transform=a.get_xaxis_transform(), zorder=4))
        fig.text((xD - 4.6) / FW, 1 - (y0 + PD / 2) / FH, c.replace('chr0', 'chr'), ha='right', va='center',
                 fontsize=6)
    yDend = yD + 1.2 + len(pchroms) * (PD + GD)
    ax_ = ax_mm(xD, yDend + 0.3, 94 * MBMM, 0.01)
    clean(ax_, bottom=True)
    ax_.set_xlim(0, 94)
    ax_.xaxis.set_major_locator(ticker.MultipleLocator(20))
    ax_.tick_params(axis='x', labelsize=6, pad=1)
    ax_.set_xlabel('potato position (Mb)', fontsize=6.5, labelpad=1)
    yDend += 6.0

    # ---------------------------------------------------------------- E: every genome: extreme on truth vs simulated chance
    xE = 110.0
    letter(xE - 14.0, yD, 'E')
    aE = ax_mm(xE, yD + 4.0, 68.0, yDend - yD - 11.0)
    aE.plot([0, 1], [0, 1], color='#BBBBBB', lw=0.6, zorder=1)
    for r in sv.itertuples():
        for sign, mk in (('peak', '^'), ('dip', 'v')):
            aE.plot(getattr(r, f'chance_{sign}'), getattr(r, sign) / r.n, marker=mk, ms=3.3, mfc=PURPLE, mec='white',
                    mew=0.3, ls='none', zorder=3)
    for r in PS.itertuples():
        for sign, mk in (('peak', '^'), ('dip', 'v')):
            aE.plot(getattr(r, f'chance_{sign}'), getattr(r, sign) / r.n, marker=mk, ms=3.5, mfc='none', mec=GRAY,
                    mew=0.6, ls='none', zorder=2)
    # names of the genomes beyond chance: (label x, label y, ha) in data units, a leader line when moved
    pos = {('tomato', 'peak'): (0.13, 1.0, 'right'), ('vdah', 'dip'): (0.215, 1.10, 'right'),
           ('spom', 'peak'): (0.29, 1.10, 'center'), ('spom', 'dip'): None, ('colcc', 'peak'): (0.43, 1.10, 'center'),
           ('cneoG1', 'peak'): (0.46, 0.93, 'left'), ('potato', 'dip'): (0.11, 0.727, 'left'),
           ('rice', 'peak'): (0.36, 0.667, 'left'), ('aetAL878', 'dip'): (0.085, 0.571, 'left'),
           ('msym', 'dip'): (0.16, 0.50, 'left'), ('barley', 'dip'): (0.07, 0.429, 'left'),
           ('ppat', 'peak'): (0.215, 0.30, 'left')}
    for r in sv.itertuples():
        for sign in ('peak', 'dip'):
            if getattr(r, f'p_{sign}') > 0.05:
                continue
            x_, y_ = getattr(r, f'chance_{sign}'), getattr(r, sign) / r.n
            p_ = pos.get((r.id, sign), (x_ + 0.02, y_, 'left'))
            if p_ is None:
                continue
            lx, ly, ha = p_
            aE.text(lx, ly, SHORT[r.id], fontsize=5.5, ha=ha, va='center', color='#333333', zorder=5)
            if abs(ly - y_) > 0.04:
                aE.plot([x_, lx], [y_ + 0.025, ly - 0.03], color='#AAAAAA', lw=0.35, zorder=1)
    aE.set_xlim(0, 0.7)
    aE.set_ylim(-0.03, 1.15)
    aE.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    aE.xaxis.set_major_formatter(ticker.FuncFormatter(lambda v_, p_: f'{v_ * 100:.0f}'))
    aE.yaxis.set_major_formatter(ticker.FuncFormatter(lambda v_, p_: f'{v_ * 100:.0f}'))
    aE.set_xlabel('chance from simulated maps (% of chromosomes)', fontsize=6.5, labelpad=1)
    aE.set_ylabel('chromosomes with the extreme\non the centromere (%)', fontsize=6.5, labelpad=1)
    for s_ in ('top', 'right'):
        aE.spines[s_].set_visible(False)
    for s_ in ('left', 'bottom'):
        aE.spines[s_].set_linewidth(0.5)
    aE.tick_params(labelsize=6, length=2, pad=1)
    aE.text(0.69, 0.60, 'chance', fontsize=5.5, color='#999999', ha='right', va='top')
    hE = [Line2D([], [], marker='^', ms=3.4, mfc=PURPLE, mec='white', mew=0.3, ls='none'),
          Line2D([], [], marker='v', ms=3.4, mfc=PURPLE, mec='white', mew=0.3, ls='none'),
          Line2D([], [], marker='^', ms=3.4, mfc='none', mec=GRAY, mew=0.6, ls='none')]
    aE.legend(hE, [f'peak  ({len(sv)} libraries)', 'dip', f'randomized Hi-C ({len(PS)})'], loc='lower right',
              bbox_to_anchor=(1.0, 0.03), frameon=False, fontsize=5.8, handletextpad=0.2, labelspacing=0.25,
              borderaxespad=0)

    # ---------------------------------------------------------------- F: every genome, mean profile around truth
    yF = yDend + 2.0
    letter(0, yF, 'F')
    order = ['aetAL878', 'barley', 'brachy', 'rice', 'sbBTx623', 'colcc', 'CC1', 'gmX26', 'potato', 'tomato', 'ppat',
             'scer', 'spom', 'ncra', 'cneoA', 'cneoG1', 'ctro', 'msym', 'vdah', 'ztri', 'pfal', 'ggswu',
             'chm13', 'mPanTro3', 'mPanPan1', 'mGorGor1', 'mPonAbe1', 'mPonPyg2', 'mSymSyn1', 'MFA8', 'MMU8', 'mmT2T',
             'horse', 'horse1', 'goat1', 'sheep1', 'mGirTip1', 'Sscrofa']
    order = [g for g in order if g in SV] + [g for g in sv.id if g not in order]
    groups = [('plants', [g for g in order if lineage(g) == 'plants']),
              ('fungi,\nprotist,\nbird', [g for g in order if lineage(g) in ('fungi', 'protists', 'birds')]),
              ('mammals', [g for g in order if lineage(g) == 'mammals'])]
    NC = 11
    BINS = {}                                       # Hi-C bin size -> genomes, for legend F
    cw, chh, gapx, rowh = 13.4, 5.2, 1.6, 8.2
    xF = 13.0
    rowi = 0
    FL = (0.4, 2.5)
    for gname, gl in groups:
        chunks = [gl[i:i + NC] for i in range(0, len(gl), NC)]
        for ci_, chunk in enumerate(chunks):
            y0 = yF + 0.8 + rowi * rowh
            if ci_ == 0:
                fig.text(1.0 / FW, 1 - (y0 + 2.0 + chh / 2) / FH, gname, ha='left', va='center', fontsize=6,
                         color='#333333', linespacing=1.0)
            for j, g in enumerate(chunk):
                d = SV[g]
                a = ax_mm(xF + j * (cw + gapx), y0 + 2.0, cw, chh)
                xx = np.arange(-FHALF, FHALF + 1)
                hw = max(d['width'] / 2, 0.5)
                a.axvspan(-hw, hw, color=BLUE, lw=0, alpha=0.75, zorder=0)
                for tr_ in d['traces']:
                    a.plot(xx, tr_, color=PURPLE, lw=0.25, alpha=0.22, zorder=1)
                if len(d['traces']):
                    with np.errstate(all='ignore'):
                        m = np.nanmean(d['traces'], axis=0)
                    a.plot(xx, m, color=PURPLE, lw=0.8, zorder=3)
                    for j2 in np.flatnonzero(~np.isfinite(m)):
                        a.add_patch(Rectangle((xx[j2] - 0.5, FL[0]), 1.0, FL[1] - FL[0], hatch='//////', fill=False,
                                              lw=0, zorder=0.6))
                a.axhline(1, color='#9E9E9E', lw=0.3, zorder=0.5)
                a.set_yscale('log')
                a.set_ylim(*FL)
                a.set_xlim(-FHALF, FHALF)
                a.yaxis.set_major_locator(ticker.FixedLocator([0.5, 1, 2]))
                a.yaxis.set_minor_locator(ticker.NullLocator())
                a.set_yticklabels(['½×', '', '2×'] if j == 0 else [])
                a.xaxis.set_major_locator(ticker.FixedLocator([-10, 0, 10]))
                last_row = (gname == 'mammals' and ci_ == len(chunks) - 1)
                a.set_xticklabels(['−10', '0', '10'] if last_row else [])
                for s_ in ('top', 'right'):
                    a.spines[s_].set_visible(False)
                for s_ in ('left', 'bottom'):
                    a.spines[s_].set_linewidth(0.4)
                a.tick_params(labelsize=5, length=1.2, pad=0.8, width=0.4)
                a.patch.set_alpha(0)
                bl = f"{d['bin'] / 1e3:.0f} kb" if d['bin'] < 1e6 else f"{d['bin'] / 1e6:.0f} Mb"
                a.set_title(FSHORT.get(g, SHORT[g]), fontsize=5.5, pad=0.8, loc='left')
                BINS.setdefault(int(d['bin']), []).append(g)
                if g in FNOTE:
                    a.text(0.99, 0.97, FNOTE[g], transform=a.transAxes, ha='right', va='top', fontsize=4.6,
                           color='#555555', zorder=6)
                if last_row and j == 0:
                    a.set_xlabel('bins from the centromere center', fontsize=6, labelpad=0.5, loc='left')
            rowi += 1
    yFend = yF + 0.8 + rowi * rowh + 3.5
    log(f'figure height used: {yFend:.1f} mm of {FH}', force=True)

    # =========================================================================== numbers
    R11 = 'REL11 survey (0.11 runs; Hi-C engine unchanged since 0.10.1)'     # every other row: 0.16.0
    S4 = 'V11/results/cis_survey4.tsv (stored 0.11 runs, V11/../runs/REL11)'
    TRT = f'{TAB.relative_to(V16)}/release16.runs.tsv'
    TCH = f'{TAB.relative_to(V16)}/release16.chrom.tsv'
    st = sv.set_index('id')
    num('A', 'tomato chromosomes with the cis peak on CENH3 (4-bin rule)', f"{st.at['tomato', 'peak']}/{st.at['tomato', 'n']}", S4, R11)
    num('A', 'tomato chance of that (mean per chromosome, simulated maps)', st.at['tomato', 'chance_peak'], S4, R11)
    num('A', 'tomato P (Poisson-binomial), peak on CENH3', st.at['tomato', 'p_peak'], S4, R11)
    er = T['ext_real']
    num('A', 'tomato peak on CENH3, recomputed here from runs/rel16/TH/tomato (the figure draws these)',
        f"{int(er.peak.sum())}/{len(er)}", 'this script (cenfuse 0.16.0 hicfuse.cis_profile on runs/rel16/TH/tomato)')
    hts = np.array([np.nanmax(T['prof'][c]) for c in tchroms])
    num('A', 'tomato peak height, median (range), obs./exp.', f'{np.median(hts):.2f} ({hts.min():.2f}-{hts.max():.2f})',
        'runs/rel16/TH/tomato matrix, this script')
    smx = np.array([T['simmax'][c] for c in tchroms])
    num('A', 'tomato: highest value of a simulated profile, median over chromosomes', f'{np.median(smx):.2f}',
        'runs/rel16/TH/tomato matrix, simulated maps as hicfuse.cis_chance, this script')
    ep = T['ext_part']
    num('A', 'tomato randomized partner: peak on CENH3', f"{int(ep.peak.sum())}/{len(ep)}",
        f'{part["tomato"].relative_to(V16)} matrix, this script')
    num('A', 'tomato randomized partner: dispersion of contacts (1 = uniform shuffle)', f"{T['partner_dispersion']:.1f}",
        f'{part["tomato"].relative_to(V16)} matrix, this script')
    sT = pd.read_csv(runs['tomato'] / 'cenfuse_summary.tsv', sep='\t', dtype=str, keep_default_na=False).iloc[0]
    sP = pd.read_csv(runs['potato'] / 'cenfuse_summary.tsv', sep='\t', dtype=str, keep_default_na=False).iloc[0]
    num('A', 'tomato cis line against the LTR cores', sT['cis_line'], 'runs/rel16/TH/tomato/cenfuse_summary.tsv cis_line')
    num('A', 'tomato chromosomes whose cis peak lies at the LTR core (engine at_core)',
        f"{sum(bool(v_) for v_ in at_core.values())}/{len(at_core)}: " +
        ', '.join(chrname[c] for c in tchroms if at_core.get(c)), 'runs/rel16/TH/tomato/work/hic/tomato.pkl cis.at_core')
    num('D', 'potato cis line against the LTR cores', sP['cis_line'], 'runs/rel16/TH/potato/cenfuse_summary.tsv cis_line')
    cores_on = sum(E.dist(int(sigT.at[c, 'core_start']), int(sigT.at[c, 'core_end']), c, trT) <= 1_000_000 for c in tchroms)
    num('A', 'tomato LTR cores within 1 Mb of CENH3', f'{cores_on}/12', 'runs/rel16/T/tomato/tomato.te_signature.tsv + truth')
    for c in tchroms:
        num('A', f'tomato {chrname[c]} ({c}) LTR core (bp, as drawn)', f"{int(sigT.at[c, 'core_start'])}-{int(sigT.at[c, 'core_end'])}",
            'runs/rel16/T/tomato/tomato.te_signature.tsv')
    ncen = int(satT[satT.cen].chrom.nunique())
    num('A', 'tomato chromosomes holding an array of the centromeric satellite family', f'{ncen}/12 ({", ".join(cenfam)})',
        f'{satdir}/tomato_satellite_arrays.gff, satellite_families.tsv (the satellite folder of the 0.16 T run, bench/jobs_rel16.tsv)')
    num('A', 'tomato satellite arrays drawn', len(satT), f'{satdir}/tomato_satellite_arrays.gff')
    num('A', 'tomato LTR signature elements drawn', len(elT), f'{gffT} (families of runs/rel16/T/tomato/tomato.te_signature.tsv)')
    r16 = pd.read_csv(TAB / 'release16.runs.tsv', sep='\t')
    grun = {(a_, g_): r16[(r16.arm == a_) & (r16.genome == g_)].iloc[0] for a_ in ('T', 'TH', 'TNEG')
            for g_ in ('tomato', 'potato')}
    for (a_, g_), r_ in grun.items():
        lab = {'T': 'without Hi-C (FASTA + tes/)', 'TH': 'with Hi-C', 'TNEG': 'with randomized Hi-C'}[a_]
        num('A' if g_ == 'tomato' else 'D', f'{g_} cenfuse 0.16 call {lab}: right / chromosomes (high, medium, low); '
            f'wrong with high or medium confidence; wrong low', f'{r_.right}/{r_.chromosomes} ({r_.right_high}, '
            f'{r_.right_medium}, {r_.right_low}); {r_.conf_wrong}; {r_.wrong_low}', f'{TRT} arm {a_}')
    for c in tchroms:
        num('A', f'tomato {chrname[c]} ({c}) calls drawn: no Hi-C | with Hi-C (start-end, 1-based; outcome)',
            f'{cNo[c][0]}-{cNo[c][1]} {cNo[c][3]} | {cHi[c][0]}-{cHi[c][1]} {cHi[c][3]}', f'{TCH} arms T, TH')
    num('A', 'truth rows checked against explib (truth string, gap, R/W, outcome) for tomato and potato, arms T and TH',
        f"{nchk['tomato']} tomato, {nchk['potato']} potato: all agree", f'{TCH}; V11/bench/explib.py truth (assertion)')
    num('B', 'example chromosome (rule: peak on CENH3 and at the LTR core; nearest the median of peak height and CENH3 width, MAD-scaled)',
        f'{chrname[EX]} ({EX})', 'this script')
    for rr in cand.itertuples():
        num('B', f'candidate {rr.name}: peak height, CENH3 width (Mb), distance, masked bins and gaps within 2 Mb',
            f'{rr.height:.2f}, {rr.width:.1f}, {rr.dist:.2f}, {rr.masked}, {rr.gaps}', 'this script')
    num('B', 'example: masked bins within 8 bins (2 Mb) of the peak', masked_near, 'runs/rel16/TH/tomato matrix (ICE mask)')
    num('B', 'example: assembly gaps (N runs) within 2 Mb of the peak bin', ngaps, f'{fa}')
    for rr in skipped.itertuples():
        num('B', f'skipped {rr.name} (nearer the medians): masked bins, gaps within 2 Mb of its peak', f'{rr.masked}, {rr.gaps}',
            'runs/rel16/TH/tomato matrix; FASTA')
    num('B', 'example: CENH3 domain (ChIP truth used for scoring)', rel_truth, f'{TCH} truth; explib.truth')
    num('B', 'example: peak bin', f'{jpx * T["bin"]:,}-{(jpx + 1) * T["bin"]:,} ({vx[jpx]:.2f}x)', 'this script')
    num('B', 'example: window drawn (bp)', f'{bwin[0]:,}-{bwin[1]:,}', 'this script (CENH3 center +- 7 Mb)')
    for g in ('tomato', 'potato'):
        pl = G[g]['pile']
        num('C', f'{g}: chromosomes in the pile-up / random windows', f"{pl['n']} / {pl['nrand']}", 'this script')
        num('C', f'{g}: mean obs./exp. in the 4+4 crossing box at the center vs random positions',
            f"{pl['box']:.2f} vs {pl['box_rand']:.2f}", 'this script')
        num('C', f'{g}: mean obs./exp. of same-side pairs within 4 bins of the center vs random',
            f"{pl['same']:.2f} vs {pl['same_rand']:.2f}", 'this script')
    num('D', 'potato chromosomes with the cis dip on CENH3 (4-bin rule)', f"{st.at['potato', 'dip']}/{st.at['potato', 'n']}", S4, R11)
    num('D', 'potato chance (mean per chromosome)', st.at['potato', 'chance_dip'], S4, R11)
    num('D', 'potato P, dip on CENH3', st.at['potato', 'p_dip'], S4, R11)
    num('D', 'potato peak on CENH3', f"{st.at['potato', 'peak']}/{st.at['potato', 'n']} (P {st.at['potato', 'p_peak']})", S4, R11)
    pd_ = G['potato']['ext_real']
    num('D', 'potato dip on CENH3, recomputed here from runs/rel16/TH/potato (the figure draws these)',
        f"{int(pd_.dip.sum())}/{len(pd_)}", 'this script (cenfuse 0.16.0 hicfuse.cis_profile on runs/rel16/TH/potato)')
    dips = pd_[pd_.dip].dip_val
    num('D', 'potato dip depth on CENH3, median (range)', f'{dips.median():.2f} ({dips.min():.2f}-{dips.max():.2f})',
        'runs/rel16/TH/potato matrix, this script')
    epp = G['potato']['ext_part']
    num('D', 'potato randomized partner: dip on CENH3', f"{int(epp.dip.sum())}/{len(epp)}",
        f'{part["potato"].relative_to(V16)} matrix, this script')
    sig_g = sv[(sv.p_peak <= 0.05) | (sv.p_dip <= 0.05)]
    num('E', 'genomes with the peak or the dip on truth beyond chance (P <= 0.05, uncorrected for two signs)',
        f'{len(sig_g)}/{len(sv)}: ' + ', '.join(PLAIN[g] for g in sig_g.id), S4, R11)
    for lin in ('plants', 'fungi', 'protists', 'birds', 'mammals'):
        ids = [g for g in sv.id if lineage(g) == lin]
        k_ = sum(g in set(sig_g.id) for g in ids)
        num('E', f'{lin}: genomes beyond chance', f'{k_}/{len(ids)}', S4 + ' ; V11/results/canonical.tsv lineage', R11)
    for r in sv.itertuples():
        num('E', f'{PLAIN[r.id]} ({r.arm}/{r.id}): peak, dip on truth / n; chance; P',
            f'{r.peak}, {r.dip} / {r.n}; {r.chance_peak}, {r.chance_dip}; {r.p_peak}, {r.p_dip}', S4, R11)
    for r in PS.itertuples():
        num('E', f'randomized {PLAIN[r.id]}: peak, dip on truth / n; chance; P; dispersion',
            f'{r.peak}, {r.dip} / {r.n}; {r.chance_peak:.3f}, {r.chance_dip:.3f}; {r.p_peak:.2g}, {r.p_dip:.2g}; {r.dispersion:.1f}',
            f'{Path(r.run).relative_to(V16)} matrix, this script')
        num('E', f'randomized {PLAIN[r.id]}: cis line against the LTR cores', r.cis_line or '.',
            f'{Path(r.run).relative_to(V16)}/cenfuse_summary.tsv')
    tested = PS[PS.cis_line.str.contains('P ', regex=False)]
    pmin = min(float(x.split('(P ')[1].split(')')[0]) for x in tested.cis_line)
    n_fail = int(tested.cis_line.str.endswith('fail').sum())
    num('E', 'randomized partners tested by the cis line: fail / tested; smallest P', f'{n_fail}/{len(tested)} fail; {pmin:.2g}',
        'runs/rel16/TNEG/*/cenfuse_summary.tsv cis_line')
    _pm = float(min(PS.p_peak.min(), PS.p_dip.min()))
    num('E', f'randomized partners: smallest truth-rule P of {2 * len(PS)} (two signs each); chance of one so small',
        f'{_pm:.2g}; {1 - (1 - _pm) ** (2 * len(PS)):.2f}', 'this script')
    yes7 = PS[PS.dispersion > 2]
    num('E', 'randomized partners with clumped contacts (dispersion > 2)', ', '.join(yes7.id) or 'none', 'this script')
    for g in sv.id:
        d = SV[g]
        m = np.nanmean(d['traces'], axis=0) if len(d['traces']) else np.full(2 * FHALF + 1, np.nan)
        num('F', f'{PLAIN[g]}: mean obs./exp. at the truth center bin; chromosomes', f'{m[FHALF]:.2f}; {len(d["traces"])}',
            f'V11/../runs/REL11/{st.at[g, "arm"]}/{g} matrix (as results/cis_survey4.tsv), this script', R11)

    # =========================================================================== legend (one page, biology first)
    # the profiles drawn (0.16 runs) must give the survey's counts (stored 0.11 runs of the same matrices)
    assert int(er.peak.sum()) == int(st.at['tomato', 'peak']) and len(er) == int(st.at['tomato', 'n']), er
    assert int(pd_.dip.sum()) == int(st.at['potato', 'dip']) and len(pd_) == int(st.at['potato', 'n']), pd_
    n_sig = len(sig_g)
    mamm = [g for g in sv.id if lineage(g) == 'mammals']
    plants = [g for g in sv.id if lineage(g) == 'plants']
    fungi = [g for g in sv.id if lineage(g) == 'fungi']
    sig_pl = [g for g in plants if g in set(sig_g.id)]
    sig_fu = [g for g in fungi if g in set(sig_g.id)]
    tNo, tHi = grun[('T', 'tomato')], grun[('TH', 'tomato')]
    pNo, pHi = grun[('T', 'potato')], grun[('TH', 'potato')]
    part_sig = [(r.id, 'peak' if r.p_peak <= r.p_dip else 'dip') for r in PS.itertuples() if min(r.p_peak, r.p_dip) <= 0.05]
    p_part_min = float(min(PS.p_peak.min(), PS.p_dip.min()))
    p_one = 1 - (1 - p_part_min) ** (2 * len(PS))
    nw = lambda k: {0: 'no', 1: 'one', 2: 'two', 3: 'three', 4: 'four', 5: 'five', 10: 'ten'}.get(int(k), str(int(k)))
    if not part_sig:
        E_part = f'Open marks: randomized Hi-C of {len(PS)} of these genomes, all near chance.'
    elif len(part_sig) == 1:
        E_part = (f'Open marks: randomized Hi-C of {len(PS)} of these genomes stays near chance; its one outlier, a '
                  f'{part_sig[0][1]} in {PLAIN[part_sig[0][0]]}, is what chance alone yields among {2 * len(PS)} tries '
                  f'about one time in {nw(round(1 / float(f"{p_one:.1g}")))}.')
    else:
        E_part = (f'Open marks: randomized Hi-C of {len(PS)} of these genomes; {len(part_sig)} of {2 * len(PS)} tries '
                  f'pass chance.')
    no_miss = int(tNo.chromosomes) - int(tNo.right)
    bstr = lambda b: f'{b / 1e3:.0f} kb' if b < 1e6 else f'{b / 1e6:.0f} Mb'
    nm_ = lambda gs: list(dict.fromkeys(FSHORT.get(x, SHORT[x]) for x in gs))
    and_ = lambda xs: xs[0] if len(xs) == 1 else ', '.join(xs[:-1]) + ' and ' + xs[-1]
    common = max(BINS, key=lambda b: len(BINS[b]))
    F_bins = (f'Hi-C bins: {bstr(common)}; ' + '; '.join(f'{bstr(b)} in {and_(nm_(BINS[b]))}' for b in
              sorted(BINS, reverse=True) if b != common) + '.')
    num('F', 'Hi-C bin sizes of the survey', F_bins, S4, R11)
    GEN = {'cneoA': 'cneo', 'cneoG1': 'cneo'}       # C. neoformans H99: two Hi-C libraries (anaphase, G1), one genome
    sig_set = set(sig_g.id)

    def by_genome(ids):                             # a genome counts once: beyond chance if any of its libraries is
        gs = {}
        for g_ in ids:
            gs.setdefault(GEN.get(g_, g_), []).append(g_ in sig_set)
        return len(gs), sum(any(v) for v in gs.values())
    n_gen, n_gen_sig = by_genome(sv.id)
    nf_gen, nf_gen_sig = by_genome(fungi)
    rest = [g for g in sv.id if lineage(g) in ('protists', 'birds')]
    assert sorted(rest) == ['ggswu', 'pfal'], rest
    assert len(plants) + nf_gen + len(rest) + len(mamm) == n_gen, (len(plants), nf_gen, len(rest), len(mamm), n_gen)
    assert len(sig_pl) + nf_gen_sig + sum(g in sig_set for g in rest + mamm) == n_gen_sig
    assert not any(g in sig_set for g in rest + mamm), 'legend E names no beyond-chance mammal, Plasmodium or chicken'
    assert 'cneoG1' in sig_set and 'cneoA' not in sig_set, 'legend E: G1 beats chance, anaphase does not'
    num('E', 'genomes beyond chance, C. neoformans H99 counted once (either library)', f'{n_gen_sig}/{n_gen}', S4, R11)
    num('E', 'fungal genomes beyond chance, C. neoformans counted once', f'{nf_gen_sig}/{nf_gen}', S4, R11)
    num('E', 'Hi-C libraries drawn (C. neoformans: anaphase and G1)', len(sv), S4, R11)
    nmam = sum(g in set(sig_g.id) for g in mamm)
    legend = [
        ('', 'The two sides of every tomato centromere contact each other about twice as often as equally distant DNA '
             'elsewhere, so a contact map alone marks centromeres that no repeat marks. Potato shows a deficit instead, '
             'and the signal recurs in some plants and fungi but in no mammal.'),
        ('A', f'On all {len(er)} tomato chromosomes, contact between the two sides of each point peaks at the CENH3 '
              f'domain; '
              f'randomized Hi-C and simulated maps show no such peak. The centromeric satellite and the retrotransposons '
              f'that mark centromeres often lie elsewhere, so calls from the genome and its transposons alone miss '
              f'{nw(no_miss)} centromeres, {nw(tNo.conf_wrong)} of them confidently (upper row); with Hi-C, all '
              f'{tHi.right} are placed (lower row).'),
        ('B', f'Tomato {chrname[EX]} (box in A) as a contact map: the excess is a compact patch of contact between the '
              f'two sides of the centromere, directly above it. The diamond holds the contacts that make the profile '
              f'value beneath.'),
        ('C', 'Averaged over chromosomes, tomato centromeres gain contact across the centromere and lose contact with '
              'the rest of the chromosome, more like a small insulated domain than arms folded together; potato '
              'centromeres lose contact instead, and random positions stay flat. The cartoons are expectations, and '
              'the mechanism is a hypothesis.'),
        ('D', f'In potato the same measure dips at most centromeres instead of peaking (triangles: lowest bin; orange '
              f'bar: the retrotransposon core). Here Hi-C adds no centromere: calls place {pHi.right} of '
              f'{pHi.chromosomes} with it and {pNo.right} without.'),
        ('E', f'In {n_gen_sig} of {n_gen} genomes with Hi-C and a centromere map, the peak or the dip falls on the '
              f'centromere more often than chance (line): {len(sig_pl)} of {len(plants)} plants, {nf_gen_sig} of '
              f'{nf_gen} fungi, and none of the {len(mamm)} mammals, $\\it{{Plasmodium}}$ or chicken. Each point is '
              f'one Hi-C library; $\\it{{C.\\ neoformans}}$ has two and counts once, as its G1 library beats chance '
              f'and its anaphase library does not. {E_part}'),
        ('F', 'The average profile around the centromeres of each genome (faint lines: single chromosomes; blue: '
              'centromere width; hatched: no data) shows which genomes peak, which dip and which stay flat. ' + F_bins),
    ]
    stamp = (f'Script: docs/img/src/{STEM}.py · data: the cenfuse 0.16.0 benchmark runs (not in the repository) · '
             f'{datetime.date.today().isoformat()}')
    for lab_, txt_ in legend:
        num('legend', f'words in paragraph {lab_ or "headline"}', len(txt_.split()), 'this script')

    # =========================================================================== write
    NUMS = pd.DataFrame(NUM)
    NUMS.to_csv(HERE / f'{STEM}_numbers.tsv', sep='\t', index=False)
    meta = {'Title': 'benchmark_cis: the contact crossing the centromere (cenfuse 0.16.0)', 'Author': None,
            'Subject': None, 'Keywords': None, 'Creator': None}
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / f'{STEM}.png', dpi=300)
    with PdfPages(OUT / f'{STEM}.pdf', metadata=meta) as pdf:
        pdf.savefig(fig, dpi=400)                  # rasterized contact maps at 400 dpi, everything else vector
        text_page(pdf, legend, 10, FW / 25.4, stamp)
    plt.close(fig)
    log(f'wrote {OUT / STEM}.pdf, .png and {HERE / STEM}_numbers.tsv; example {chrname[EX]} ({EX}); partner dispersion '
        f'tomato {T["partner_dispersion"]:.1f}, potato {G["potato"]["partner_dispersion"]:.1f}', force=True)


def IT(s):
    """Italic species name in legend text (mathtext; spaces escaped so a name wraps as one word)."""
    return '$\\it{' + s.replace(' ', '\\ ') + '}$'


def text_page(pdf, paras, fs, w_in, stamp=None):
    """A legend page as wide as the figure, height fitted to the text (wrapped by measuring rendered widths)."""
    import matplotlib.pyplot as plt
    margin = 0.25
    fig = plt.figure(figsize=(w_in, 11))
    r = fig.canvas.get_renderer()
    maxw = (w_in - 2 * margin) * 72
    lines = []

    def width(s):
        t = fig.text(0, 0, s, fontsize=fs)
        bb = t.get_window_extent(renderer=r)
        t.remove()
        return bb.width * 72 / fig.dpi

    for lab, txt in paras:
        words = txt.split(' ')
        cur, first = '', True
        prefix = f'$\\bf{{{lab}}}$  ' if lab else ''
        for wd in words:
            trial = (cur + ' ' + wd) if cur else wd
            if width((prefix if first else '') + trial) > maxw and cur:
                lines.append((prefix if first else '') + cur)
                cur, first = wd, False
            else:
                cur = trial
        lines.append((prefix if first else '') + cur)
        lines.append('')
    stamps = []
    if stamp:                                                    # one provenance line, wrapped at '/' or ' · '
        cur = ''
        for tok in stamp.replace('/', '/\x00').replace(' · ', ' · \x00').split('\x00'):
            if cur and width(cur + tok) * 0.8 > maxw:
                stamps.append(cur)
                cur = '    ' + tok
            else:
                cur += tok
        stamps.append(cur)
    lines += stamps
    lh = fs * 1.32 / 72
    H = margin * 2 + lh * len(lines)
    plt.close(fig)
    fig = plt.figure(figsize=(w_in, H))
    y = 1 - margin / H
    for ln in lines:
        fig.text(margin / w_in, y, ln, fontsize=fs if ln not in stamps else fs * 0.8, va='top', ha='left')
        y -= lh / H
    pdf.savefig(fig)
    plt.close(fig)


def main():
    global VERBOSE
    ap = argparse.ArgumentParser(description='Build benchmark_cis (cenfuse 0.16.0): PDF (figure + one legend page), PNG '
                                             '(figure only, 300 dpi) and the numbers table, next to this script.')
    ap.add_argument('--refresh', action='store_true', help='recompute every matrix-derived value (ignore the cache)')
    ap.add_argument('-v', '--verbose', action='store_true', help='print per-step progress')
    args = ap.parse_args()
    VERBOSE = args.verbose
    build(args)


if __name__ == '__main__':
    main()
