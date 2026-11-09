# benchmark/

The benchmark of cenfuse, inside the repository ([future.md #34](../docs/future.md#34-the-benchmark-cannot-be-rerun-from-the-repository)):
a truth manifest, the frozen scorer as a tested module, and the release chain that scores stored runs, writes every
table of the docs and fills `README.md` and `docs/benchmark.md` from templates. It is a port of the 0.13 development
folder's scripts (`v12/bench/finish_rel13f.sh` and what it calls). On the 0.13.0 release runs it reproduces every
stored output byte for byte, and the filled docs equal the repository's (see Verification).

**The data stay outside the repository.** Truth BEDs, stored runs, job tables and the earlier releases' tables live
under the data root, `/data2/chris/centromere_benchmark`, or `$CENFUSE_BENCH_DATA` if set. The manifest only points to
them. Some truth is unpublished: the ape CDRs are the Logsdon lab's curated set (v1.3), "contact the lab before
publication" (benchmark.md credits Logsdon et al. 2024; future.md #37).

## Files

| file | what |
|---|---|
| `truth.tsv` | the manifest: one row per genome (58 of the frozen scorer + `agen`, `sver` and 13 tomatoes added in 0.14/0.15): `kind`, `rule`, `tol_bp`, `hic_bins`, `truth` (paths relative to the data root, `;`-separated; `.` for holocentric), `source` (assembly, data accessions, paper, how the BED was made), `notes` (caveats) |
| `paths.py` | every root the chain reads, in one place: the data root, the development folders, the 0.13.0 release's runs, job tables and reference tables. `python -m benchmark.paths` lists them |
| `score.py` | the frozen scorer (`v12/bench/rel/score_rel.py`), reading the manifest |
| `compare_rel.py` | release runs vs their reference runs, column by column |
| `rows.py` | the docs' result rows, beside the previous release's |
| `md_tables.py`, `md_013.py` | the tables of `docs/benchmark.md` |
| `readme_table.py` | the README's results table |
| `rabl_table.py` | the Rabl-axis table |
| `lines_alone.py` | each line alone vs the integrated call |
| `score_k2_rel.py` | the karyotype check K2 (v7 `bench/score_k.py`, inlined) |
| `randomized.py` | the randomized Hi-C controls: tests passed, calls changed |
| `score_dev.py` | calls and counts changed against a reference release |
| `headroom.py` | where a line's own locus is right but the final call is not |
| `precision.py` | how close the right calls are: width and midpoint against the truth domain |
| `fill_docs.py` | fills `docs_src/*.template.md` with the computed tables and the hand-written texts |
| `release.sh` | runs the chain for a release root |
| `docs_src/` | the doc templates and their hand-written texts (copied from `v12/docs_src`; edit them here, not the generated docs) |
| `canonical.tsv` | genome labels, lineages and the arm each doc table shows (copied from `v11/results`) |

Each module runs as `python -m benchmark.<name> --help` from the repository root.

## Rerun

Use the `cenfuse` mamba env (numpy, pandas, pytest). That env holds a stale editable install of cenfuse 0.8.0 (#34),
so run from the repository root with `PYTHONPATH` set to it. `release.sh` stops when python imports cenfuse from
anywhere else.

```bash
cd /data2/chris/centromere_benchmark/v2/v3/v4/v5/v6/v7/v8/v9/v10/v10.1/v10.2/v10.3/v10.4/v11/v12/v13/cenfuse
export PYTHONPATH=$PWD OMP_NUM_THREADS=1 NUMPY_MADVISE_HUGEPAGE=0
PY=/home/chris/bin/mambaforge/envs/cenfuse/bin/python
V12=/data2/chris/centromere_benchmark/v2/v3/v4/v5/v6/v7/v8/v9/v10/v10.1/v10.2/v10.3/v10.4/v11/v12
$PY -m benchmark.paths                                                    # roots in use
mkdir -p out && $PY -m benchmark.score $V12/bench/jobs_rel13.tsv $V12/runs/rel13f out/release13   # scores (5 s)
PY=$PY bash benchmark/release.sh out/rel 13                               # the whole chain (1 min)
diff out/rel/docs/README.md README.md && diff out/rel/docs/docs/benchmark.md docs/benchmark.md
DOCS=$PWD PY=$PY bash benchmark/release.sh out/rel 13                     # rewrite the repo's docs from docs_src/
$PY -m pytest -q tests/benchmark -p no:cacheprovider                      # 17 tests
```

`release.sh OUT [TAG]` writes `OUT/release$TAG.*` (scores), `cmp$TAG.*`, `rows_release$TAG.txt`, `md$TAG/` (the doc
tables), `lines_alone$TAG.*`, `k2_release$TAG.*`, `randomized_release$TAG.*`, `dev$TAG.vs_prev.*`,
`dev$TAG.vs_extra.*`, `headroom$TAG.*` and the filled docs in `OUT/docs/`. A failing table is reported and the
chain goes on (exit 1 at the end); failed scoring stops it. For another release, set the inputs in the environment:
`RUNS` (release root, `RUNS/<arm>/<id>`), `JOBS`, `JOBS_CMP` (column 6 = the previous release's run of each arm),
`JOBS_RABL`, `PREV` (the previous release's scores, a prefix), `PREV_RUNS`, `PREV_REDO`, `EXTRA_REF` (empty: skipped),
`CIS_SURVEY`, `VERSIONS` (`"0.13 0.14"` for md_013's headers). Job tables hold absolute paths; under
`$CENFUSE_BENCH_DATA` they are moved to the new root (`paths.local`).

## Scoring rules

Frozen when each benchmark was first run; `score.py` reproduces them exactly.

| rule | right when | tolerance | genomes |
|---|---|---|---|
| `mammal` | the call overlaps the CENP-A interval by >= 1 bp | none | 9 T2T mammals with CENP-A truth |
| `overlap` | >= 1-bp overlap with a truth interval | none | 5 mammals scored against the satellome locus (kind `agreement`); *Arabidopsis* |
| `active` | >= 1-bp overlap with a CenSat row named `active_hor*` | none | HG002 |
| `gap`, `hic_bins` 4 | gap <= 4 x the run's Hi-C bin; without Hi-C, gap <= `tol_bp` | 1 Mb (rice, potato, horse EquCab3); barley 8 Mb (4 x its 2-Mb bin; 4 Mb in the arm with 1-Mb bins) | rice, potato, barley, horse EquCab3 |
| `gap` | gap <= `tol_bp`, with or without Hi-C | 1 Mb; 100 kb (moss, *Chlamydomonas*, *Chlorella*); 0 (chicken, mouse) | the other plants and algae, chicken, mouse, `agen`, `sver` |
| `bin1` | gap <= 1 Hi-C bin of the run; without Hi-C, gap <= 250 kb | 10 or 50 kb with Hi-C | fungi, *Plasmodium* |
| `holo` | never: no single centromere exists | | *C. elegans*, *C. nigoni*, *B. mori*, *R. breviuscula* |

- Calls are 1-based inclusive (cenfuse's tables), truth 0-based half-open (BED). The gap is
  `max(0, max(a, ts) - min(z, te))`: 0 when the intervals overlap or touch, so `gap` with tolerance 0 counts a
  touching call as right, which the overlap rules do not. With several truth intervals on a chromosome, the nearest
  is used, a hit first.
- `mammal` reads tier FW: the W file (ChIP), then the F file (CDRs) over it, keeping one interval per chromosome, from
  the last line read. Each mammal has only one of the two files today.
- Only chromosomes with truth are scored; a holocentric control scores every chromosome, a confident call there
  being confidently wrong.
- Outcomes: `right_high`, `right_medium`, `right_low`, `conf_wrong` (a wrong `high` or `medium` call), `wrong_low`,
  `none` (no call). Confident = `high` or `medium`.
- Chance of a right call by luck = (call + truth + 2 x tolerance) / chromosome length, at most 1 (0 for no call and
  for holocentric controls); `P` = Poisson-binomial tail of the right count over the run's chromosomes.

Where `docs/benchmark.md` ("How calls are scored") says less than the scorer does:
- Small genomes without Hi-C are scored within 250 kb (`bin1`), against 10-50 kb with Hi-C. *S. cerevisiae*'s AT
  troughs are "right" on 12 of 16 that way and overlap the centromere on 2 (future.md #37).
- "Plants: within 4 Hi-C bins (1 Mb without Hi-C)" holds for rice, potato, barley and EquCab3 only. The other plants
  are scored at a fixed 1 Mb (100 kb for moss and the algae) with or without Hi-C, so *Ae. tauschii*, whose Hi-C bin
  is 2 Mb, is held to half a bin.
- Barley's 8 Mb applies to its arms without Hi-C as well.
- Chicken and mouse are scored at gap 0 ("overlap" in the docs), which also counts touching intervals.

## Truth

Each row of `truth.tsv` gives the assembly, the data and their accessions, the paper and how the BED was made. In
short: CENP-A CDRs (CHM13: Logsdon et al. 2024, Supp. Table 7; five great apes: the lab's unpublished curated set) or
CENP-A ChIP (crab-eating macaque, pig and the T2T horse, all from other animals; EquCab3 from its own animal); CENH3
ChIP of the assembled line or of a related one (sorghum V9 on BTx623; *Arabidopsis* Col-0 on Col-CC), mostly turned
into domains by one frozen rule (100-kb bins, ChIP/input, half-maximum run, >= 4x; v2 `scripts/map_chip.sh` and
`chip_truth.py`); published domains (rice, potato, moss, *B. rapa* TU1 and CX1 whose ChIP is too weak); CUT&Tag
(*Chlamydomonas*); the authors' CENP-A windows (chicken); the CENP-A minor-satellite array (mouse); the CenSat active
arrays (HG002); published CEN, CENP-A or CenH3 coordinates (fungi, *Plasmodium*); another caller's satellome locus for
the 5 mammals without functional data (scored as agreement, not truth).

New for 0.14, held out (not in the 0.13 release): `agen`, *Ae. geniculata* GCA_053567775.1, CENH3 ChIP peaks of
Yunusbaev et al. 2026 Table S2 at Mb resolution (development data for two cenrep LTR rules, future.md #9); `sver`,
*S. verrucosum* GCA_031230405.1, CENH3 ChIP by the frozen rule with an input stand-in from another accession. Both are
scored at 1 Mb. A genome whose truth file is missing (`sver` until its truth is computed) is skipped with one
warning; its runs are not scored.

To add a genome: append a row to `truth.tsv`, with the rule fixed before its runs are scored and its source given.

## Caveats (future.md #37)

- Overlap is lenient for wide calls: mammal calls are a median 3.6 Mb wide against CENP-A domains of 0.1-0.5 Mb. No
  width-normalized or midpoint score is reported.
- The chance model draws calls uniformly along the chromosome, which is optimistic where centromeres sit near a known
  end (acrocentric goat, sheep, horse, mouse).
- Weak or indirect truth: pig and horse ChIP is from other animals (pig's 3 confident misses lie where that ChIP's
  second peak falls inside the call); siamang, sheep, goat, rhesus and giraffe are scored against another caller;
  the Bornean orangutan's CM054652.2 may be haplotype-specific (#38); *N. crassa* LG VI carries an inversion;
  *C. tropicalis* CENP-A was carried over from an older assembly; *V. dahliae*'s Hi-C was used by its authors to infer
  its centromeres; *C. neoformans* counts twice (two Hi-C states).
- The truth rule's 100-kb bins and 4x have no stated reason, and the rule decides most plant truths (cucumber: 1 of 7
  chromosomes passes).
- Shallow controls: the holocentric and held-out Hi-C libraries are prefixes of the data (*R. breviuscula* 3.3 M pairs,
  *B. mori* 5.9 M, *Ae. tauschii* 0.84 M, sorghum 3.9 M, mouse 15.7 M).
- Attribution: the ape CDRs and the ape DiMeLo-seq and CUT&RUN tracks are unpublished (see above).
- No longer reported: the head-to-head with the independent satellite caller (0.1.0's B7) and the unique-read Hi-C
  control (C3, #52).
- Typed-in results: `readme_table.py` writes a few results as text ("a kinetochore on 40", "all 20 at the right end",
  "unchanged"), and the hand-written texts in `docs_src/` quote numbers. Check them against the tables each release.

## Changes from the v12 scripts

None changes an output of the 0.13 release.
- The truth table comes from `truth.tsv`, not code; the *B. rapa* fallback to Table S15 when the ChIP BED is empty is
  written into the manifest (TU1, CX1). Every path comes from `paths.py`; every script takes its inputs as arguments.
- A missing truth file skips its genome with a warning (v12 never checked), and scoring with no run left exits with a
  message (v12 crashed).
- `md_013` writes `one_per_chromosome.tsv` and `cis_trans.tsv` into its output folder (v12 wrote them into
  `v12/results`), takes the two releases' names for its headers (`--versions`), and counts the randomized libraries
  where the cis test passes instead of typing "never passes".
- `score_dev` takes a release root and an output prefix instead of a tag under `v12/runs` and `v12/results`.
- `score_k2_rel` inlines v7's `k2()` and reads the truth files of the manifest (the same files).
- `rabl_table` reads `canonical.tsv` from here; v12's first version read a nonexistent `v12/results/canonical.tsv`,
  failed in the release chain, and was rerun by hand with the right path.
- Blank and `#` lines of job tables are skipped; `pairs_kept` reads a sorted list of Hi-C QC tables (each run has one).

## Verification (2026-11-07)

On `v12/bench/jobs_rel13.tsv` and `v12/runs/rel13f`, `release.sh` reproduced byte for byte the 45 stored outputs of
`v12/results` that it rewrites (46 comparisons): `release13f.*` and its copies `release13.{chrom,runs}.tsv`, `cmp13.*`,
`rows_release13_vs_12.txt`, all 21 files of `md13/`, `one_per_chromosome13.tsv`, `cis_trans13.tsv`,
`lines_alone13.*`, `k2_release13.*`, `randomized_release13.*`, `rel13f.vs12.txt`, `rel13f.vs13.txt`,
`rel13f.changes.tsv`, `rel13f.chrom.tsv` and `headroom13.*`. The filled `README.md` and `docs/benchmark.md` equal the
repository's. The only difference, in `md13.log`, is v12's failed first rabl_table run (above) and the output path.
Also identical: the frozen scorer and `score.py` in `ref` mode, on the 0.11 and the 0.12 runs of every arm (3,504 and
3,677 chromosome rows); and the manifest against the frozen truth table, for all 58 genomes (kind, rule, tolerance,
bins, files and the intervals read). Seven mutations of the scoring rules each fail the tests.

## Still outside the repository (the rest of #34)

- The negative-control generators: `v2/scripts/permute_pairs.py` (re-paired and relocated Hi-C, `--seed 1`),
  `v11/bench/make_tneg.py`, `at_shuffle.py`, `control_te.py`; not `v11/bench/repair_pairs.sh` (DO NOT USE).
- The truth-making scripts (v2 `scripts/prep_truth.py`, `map_chip.sh`, `chip_truth.py`; the fungal builders of v10
  and v10.1), the cis survey that `TABLE_CIS` reads (`v11/bench/cis_survey4.py`), and the analysis scripts the docs
  cite (`v11/bench/pa1a_age.py`, `young_core.py`, `poly_one.py`).
- The runs themselves: the job tables, the satellite replays (`--satellite-dir`), checkpoint seeding
  (`v12/bench/seed_run.py`), and the gallery (`v12/bench/rel/make_gallery.py`, not in the chain).

## 0.15

The 0.15.0 tables were made from `v14/runs/rel15b` by `v14/bench/docs15.sh` (outside the repository): `release.sh` against the 0.14
release (`JOBS` = the 260 arms of 0.14, `JOBS_CMP` = those arms with column 6 at the 0.14 runs), and, for the genomes outside that
job table (agen, sver and the 13 tomatoes, which `truth.tsv` now names), the row-making scripts of v13 `bench/`. A truth-free check of LTR-line
rules, homoeolog agreement on the polyploid panel (future.md #68), is in `v14/bench/polyeval.py`.
