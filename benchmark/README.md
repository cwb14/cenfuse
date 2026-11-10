# benchmark/

For developers: the truth manifest, the frozen scorer, and the scripts that score a release's stored runs, write every
table of the docs and fill `README.md` and `docs/benchmark.md` from templates. Results and how to read them are in
[docs/benchmark.md](../docs/benchmark.md).

**The data stay outside the repository.** Truth BEDs, stored runs and job tables live under the data root,
`/data2/chris/centromere_benchmark` (or `$CENFUSE_BENCH_DATA`). The manifest only points to them. Some truth is
unpublished: the ape CDRs are the Logsdon lab's curated set (v1.3), "contact the lab before publication"
(docs/benchmark.md credits Logsdon et al. 2024).

## Files

| file | what |
|---|---|
| `truth.tsv` | the manifest, one row per genome: `kind`, `rule`, `tol_bp`, `hic_bins`, `truth` (paths relative to the data root, `;`-separated; `.` for holocentric), `source` (assembly, accessions, paper, how the BED was made), `notes` |
| `score.py` | the frozen scorer, reading the manifest |
| `paths.py` | every root the scripts read; `python -m benchmark.paths` lists them. The defaults are the 0.13.0 release's; later releases set them in the environment (below) |
| `release.sh` | the whole chain for one release: scores, every table, the filled docs |
| `fill_docs.py` | fills `docs_src/*.template.md` with the computed tables and the hand-written texts |
| `docs_src/` | the doc templates and hand-written texts: edit these, not `README.md` or `docs/benchmark.md` |
| `canonical.tsv` | genome labels, lineages and the arm each doc table shows |
| `compare_rel.py`, `rows.py`, `score_dev.py` | release runs against reference runs; result rows beside the previous release; calls changed |
| `md_tables.py`, `md_013.py`, `readme_table.py`, `rabl_table.py`, `precision.py` | the tables of `docs/benchmark.md` and the README |
| `lines_alone.py`, `headroom.py`, `score_k2_rel.py`, `randomized.py` | each line alone; where a line is right but the call is not; the karyotype check; the randomized-Hi-C controls |

Each module runs as `python -m benchmark.<name> --help` from the repository root.

## Regenerate the docs

Use the `cenfuse` mamba env with `PYTHONPATH` set to the repository: that env holds a stale editable install of
cenfuse 0.8.0, and `release.sh` stops when python imports cenfuse from anywhere else.

After editing a template or a hand-written text, refill from the release's stored tables (seconds):

```bash
cd <repo>
PY=/home/chris/bin/mambaforge/envs/cenfuse/bin/python
T=/data2/chris/centromere_benchmark/v2/v3/v4/v5/v6/v7/v8/v9/v10/v10.1/v10.2/v10.3/v10.4/v11/v12/v13/v14/results/docs15_audit/md15
PYTHONPATH=$PWD $PY -m benchmark.fill_docs $T .          # rewrites README.md and docs/benchmark.md
PYTHONPATH=$PWD $PY -m pytest -q tests/benchmark -p no:cacheprovider
```

The whole chain for a release: `bash benchmark/release.sh OUT TAG`, with the release's inputs in the environment:
`RUNS` (release root, `RUNS/<arm>/<id>`), `JOBS`, `JOBS_CMP` (column 6 = the previous release's run of each arm),
`JOBS_RABL`, `PREV` (the previous release's scores, a prefix), `PREV_RUNS`, `PREV_REDO`, `EXTRA_REF` (empty: skipped),
`CIS_SURVEY`, `DOCS` (where to write the filled docs; the repository root rewrites its own).
It writes the scores, every table (`OUT/mdTAG/`) and the filled docs; a failing table is reported and the chain goes on
(exit 1 at the end); failed scoring stops it. The 0.15.0 docs were made by `v14/bench/docs15.sh rel15c --fill`
(outside the repository), which runs `release.sh` against the 0.14 release and adds the rows of the genomes outside its
job table (*Ae. geniculata*, *S. verrucosum*, 13 tomatoes) with v13 `bench/` scripts.

## Scoring rules

Frozen when each benchmark was first run; `score.py` reproduces them exactly.

| rule | right when | tolerance | genomes |
|---|---|---|---|
| `mammal` | the call overlaps the CENP-A interval by >= 1 bp | none | 9 T2T mammals with CENP-A truth |
| `overlap` | >= 1-bp overlap with a truth interval | none | 5 mammals scored against the satellome locus (kind `agreement`); *Arabidopsis* |
| `active` | >= 1-bp overlap with a CenSat row named `active_hor*` | none | HG002 |
| `gap`, `hic_bins` 4 | gap <= 4 x the run's Hi-C bin; without Hi-C, gap <= `tol_bp` | 1 Mb (rice, potato, horse EquCab3); barley 8 Mb | rice, potato, barley, horse EquCab3 |
| `gap` | gap <= `tol_bp`, with or without Hi-C | 1 Mb; 100 kb (moss, *Chlamydomonas*, *Chlorella*); 0 (chicken, mouse) | the other plants and algae, chicken, mouse |
| `bin1` | gap <= 1 Hi-C bin of the run; without Hi-C, gap <= 250 kb | 10 or 50 kb with Hi-C | fungi, *Plasmodium* |
| `holo` | never: no single centromere exists | | *C. elegans*, *C. nigoni*, *B. mori*, *R. breviuscula* |

- Calls are 1-based inclusive, truth 0-based half-open (BED). The gap is `max(0, max(a, ts) - min(z, te))`: 0 when the
  intervals overlap or touch, so `gap` with tolerance 0 counts a touching call as right. With several truth intervals
  on a chromosome, the nearest is used, a hit first.
- `mammal` reads the ChIP file, then the CDR file over it, keeping one interval per chromosome.
- Only chromosomes with truth are scored; a holocentric control scores every chromosome.
- Outcomes: `right_high`, `right_medium`, `right_low`, `conf_wrong` (a wrong `high` or `medium` call), `wrong_low`,
  `none`. Confident = `high` or `medium`.
- Chance of a right call by luck = (call + truth + 2 x tolerance) / chromosome length, at most 1; `P` = the
  Poisson-binomial tail of the right count over the run's chromosomes.
- Where docs/benchmark.md simplifies: plants other than rice, potato, barley and EquCab3 are held to a fixed 1 Mb
  whatever their Hi-C bin (*Ae. tauschii*, with 2-Mb bins, to half a bin); small genomes without Hi-C are scored
  within 250 kb, which is lenient (*S. cerevisiae*'s AT troughs are "right" on 12 of 16 that way and overlap the
  centromere on 2).

## Truth

Each row of `truth.tsv` gives the assembly, the data and their accessions, the paper and how the BED was made. In
short: CENP-A CDRs (CHM13: Logsdon et al. 2024; five great apes: the lab's unpublished curated set); CENP-A ChIP
(crab-eating macaque, pig and the T2T horse, from other animals; EquCab3 from its own animal); CENH3 ChIP of the
assembled line or a related one, turned into domains by one frozen rule (100-kb bins, ChIP/input, half-maximum run,
>= 4x); published domains (rice, potato, moss, *B. rapa* TU1 and CX1); CUT&Tag (*Chlamydomonas*); the authors'
CENP-A windows (chicken); the CENP-A minor-satellite array (mouse); the CenSat active arrays (HG002); published
CEN, CENP-A or CenH3 coordinates (fungi, *Plasmodium*); *Ae. geniculata*: CENH3 ChIP peaks of Yunusbaev et al. 2026
(Table S2, Mb resolution); *S. verrucosum*: CENH3 ChIP by the frozen rule, with input from
another accession; another caller's satellome locus for 5 mammals without functional data (agreement, not
truth). A genome whose truth file is missing is skipped with a warning.

To add a genome: append a row to `truth.tsv`, with its rule fixed before its runs are scored and its source given.

## Caveats

- Overlap is lenient for wide calls: mammal calls are a median 3.6 Mb wide against CENP-A domains of 0.1-0.5 Mb
  (`precision.py` reports width and midpoint distance).
- The chance model draws calls uniformly along the chromosome, which is optimistic where centromeres sit near one end
  (acrocentric goat, sheep, horse, mouse).
- Weak or indirect truth: pig and horse ChIP is from other animals; siamang, sheep, goat, rhesus and giraffe are scored
  against another caller; the Bornean orangutan's CM054652.2 may be haplotype-specific; *N. crassa* LG VI carries an
  inversion; *C. tropicalis* CENP-A was carried over from an older assembly; *V. dahliae*'s Hi-C was used by its
  authors to infer its centromeres; *C. neoformans* counts twice (two Hi-C states).
- The truth rule's 100-kb bins and 4x have no stated reason, though 3x, 5x or 50-kb bins change no outcome on the
  chromosomes both rules give truth.
- Shallow controls: some Hi-C libraries are prefixes of the data (*R. breviuscula* 3.3 M pairs, *B. mori* 5.9 M,
  *Ae. tauschii* 0.84 M, sorghum 3.9 M, mouse 15.7 M).
- `readme_table.py` and the hand-written texts in `docs_src/` type in a few results ("a kinetochore on 40",
  "unchanged"): check them against the tables at each release.

## Still outside the repository

- The negative-control generators used for the benchmark (`cenfuse control` reproduces the Hi-C ones byte for byte).
- The truth-making scripts (v2 `scripts/prep_truth.py`, `map_chip.sh`, `chip_truth.py`; the fungal builders of v10 and
  v10.1), the cis survey (`v11/bench/cis_survey4.py`) and the analysis scripts the docs cite.
- The runs themselves: job tables, satellite replays (`--satellite-dir`), and the gallery
  (`v14/bench/make_gallery15.py`).
- The truth-free check of LTR rules on 106 allotetraploids (homoeolog agreement): `v14/bench/polyeval.py`.
