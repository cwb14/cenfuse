# benchmark/

For developers: every scored result of the 0.16.1 release, the truth manifest, the frozen scorer, and the chain that
scores a release's runs and refills [docs/benchmark.md](../docs/benchmark.md), which presents the results.

## Files

| file | what |
|---|---|
| `results.tsv` | one row per scored chromosome-run of the 0.16.1 release (the 0.16.0 runs, 51 of them rerun with 0.16.1): call, confidence, rule, outcome, truth, the repeat-family window and the CpG kinetochore; 1-based inclusive coordinates throughout (truth start = BED start + 1); its header lines give the counts and define the columns |
| `summary.py` | summaries of `results.tsv`: right, wrong and no call per input set and centromere class, confidence calibration, and the controls (`python -m benchmark.summary`) |
| `truth.tsv` | the manifest, one row per genome: `kind`, `rule`, `tol_bp`, `hic_bins`, `truth` (files under the data root, `;`-separated; `.` for holocentric), `source` (assembly, accessions, paper, how the BED was made), `notes` |
| `score.py` | the frozen scorer, reading the manifest |
| `canonical.tsv` | genome labels, lineages and the arm each doc table shows |
| `paths.py` | every input of the chain, under the data root; defaults to the 0.16.1 release (`python -m benchmark.paths` lists them) |
| `release.sh` | the chain for one release: scores, every computed table, the filled docs/benchmark.md |
| `fill_docs.py`, `docs_src/benchmark.template.md` | the template of docs/benchmark.md and its filler: edit the template, never docs/benchmark.md |
| `md_tables.py`, `md_013.py`, `rabl_table.py`, `precision.py` | the computed tables of docs/benchmark.md |
| `lines_alone.py`, `score_k2_rel.py`, `randomized.py`, `headroom.py` | numbers the docs cite: each line alone, the karyotype check, the randomized-Hi-C controls, where a line is right but the call is not |
| `compare_rel.py`, `score_dev.py`, `rows.py` | a release against the previous one: changed calls, changed outcomes, result rows side by side |

Each module runs as `python -m benchmark.<name> --help` from the repository root.

## Score a release and refill the docs

Use a python with numpy and pandas that imports cenfuse from this repository (`PYTHONPATH=$PWD`): `release.sh` stops if
it imports cenfuse from anywhere else. Point `CENFUSE_BENCH_DATA` at the data root (below).

```bash
export CENFUSE_BENCH_DATA=/path/to/centromere_benchmark
PYTHONPATH=$PWD python -m benchmark.paths                       # every input and its default (the 0.16.1 release)
RUNS=/path/to/runs PREV=/path/to/release16 PY=python bash benchmark/release.sh OUT 17
```

`release.sh OUT TAG` writes the scores (`OUT/releaseTAG.*`), every computed table (`OUT/mdTAG/`), the checks against the
previous release and the filled `OUT/docs/docs/benchmark.md`; `DOCS=<repository root>` rewrites the repository's own. Any input of
`paths.py` can be set in the environment (`RUNS`, `JOBS`, `JOBS_CMP`, `JOBS_RABL`, `PREV`, `PREV_RUNS`, `PREV_REDO`,
`EXTRA_REF`, `CIS_SURVEY`); an empty value skips the optional ones. A failing table is reported and the chain goes on
(exit 1 at the end); failed scoring stops it. Eight placeholders stay unfilled: the rows of the genomes outside the job
table (*Ae. geniculata*, *S. verrucosum*, the 13 tomatoes and the tomatoes' own Hi-C) come from scripts kept with the runs.

After editing the template, refill from a release's stored tables (seconds), then run the tests:

```bash
T=$(PYTHONPATH=$PWD python -m benchmark.paths | awk '$1 == "TABLES" {print $2}')     # the 0.16.1 tables
PYTHONPATH=$PWD python -m benchmark.fill_docs "$T" .          # rewrites docs/benchmark.md; prints any unfilled name
PYTHONPATH=$PWD python -m pytest -q tests/benchmark -p no:cacheprovider
```

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

- The gap is the number of bases strictly between call and truth: 0 when they overlap or touch, so `gap` with
  tolerance 0 counts a touching call as right. `score.py` reads calls 1-based inclusive and truth BEDs 0-based
  half-open, so the gap is `max(0, max(a, ts) - min(z, te))` with `a` = call start - 1. With several truth intervals on
  a chromosome, the nearest is used, a hit first.
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
>= 4x); published domains (rice, potato, moss, *B. rapa* TU1 and CX1, 11 of the 13 tomatoes); CUT&Tag
(*Chlamydomonas*); the authors' CENP-A windows (chicken); the CENP-A minor-satellite array (mouse); the CenSat active
arrays (HG002); published CEN, CENP-A or CenH3 coordinates (fungi, *Plasmodium*); *Ae. geniculata*: CENH3 ChIP peaks of
Yunusbaev et al. 2026 (Table S2, Mb resolution); *S. verrucosum*: CENH3 ChIP by the frozen rule, with input from another
accession; another caller's satellome locus for 5 mammals without functional data (agreement, not truth). A genome
whose truth file is missing is skipped with a warning.

To add a genome: append a row to `truth.tsv`, with its rule fixed before its runs are scored and its source given.

## Caveats

- Overlap is lenient for wide calls: mammal calls are a median 3.6 Mb wide against CENP-A domains of 0.1-0.5 Mb
  (`precision.py` reports width and midpoint distance).
- The chance model draws calls uniformly along the chromosome, which is optimistic where centromeres sit near one end
  (acrocentric goat, sheep, horse, mouse).
- Weak or indirect truth: pig and horse ChIP is from other animals; siamang, sheep, goat, rhesus and giraffe are scored
  against another caller; the Bornean orangutan's CM054652.2 may be haplotype-specific; *N. crassa* LG VI carries an
  inversion; *C. tropicalis* CENP-A was carried over from an older assembly; *V. dahliae*'s Hi-C was used by its
  authors to infer its centromeres; *C. neoformans* counts twice (one assembly, two Hi-C states).
- The truth rule's 100-kb bins and 4x have no stated reason, though 3x, 5x or 50-kb bins change no outcome on the
  chromosomes both rules give truth.
- Shallow controls: some Hi-C libraries are prefixes of the data (*R. breviuscula* 3.3 M pairs, *B. mori* 5.9 M,
  *Ae. tauschii* 0.84 M, sorghum 3.9 M, mouse 15.7 M).
- `results.tsv` shows one truth interval per chromosome-run, the one the scorer used for the call. 19 chromosomes have
  several (18 of HG002's, *Ae. geniculata* chr05M), so another window checked against the interval shown can disagree
  with the scorer: chr05M's repeat-family window overlaps its other CENH3 interval and is right.
- The template's hand-written text types in many results: check it against the tables at each release.

## Outside the repository

The runs and their inputs stayed on the server where cenfuse was built, in one nested folder per release (`v2` to `v16`,
as `paths.py` spells out); `CENFUSE_BENCH_DATA` points the scripts at a copy. They stay out because of size and rights:

- **Too large for git:** the assemblies, Hi-C pairs, methylation tracks, LTR annotations and truth BEDs, and the runs
  themselves (318 for 0.16.0, with their checkpoints and the stored satellite searches they replay).
- **Unpublished:** the five great apes' CDRs are the Logsdon lab's curated set (v1.3), "contact the lab before
  publication" (docs/benchmark.md credits Logsdon et al. 2024). `results.tsv` keeps their outcomes and leaves their
  coordinates empty.
- **Scripts tied to that data:** the job tables, the doc chain's steps after `release.sh` (the genomes outside its job
  table), the truth-making scripts (ChIP mapping, domain calls), the cis survey, the homoeolog-agreement check of the
  106-allotetraploid panel, the gallery's maker, the analysis scripts the docs cite, and the generators of the
  negative controls (`cenfuse control` reproduces the randomized-Hi-C ones byte for byte).

`results.tsv` holds every chromosome-run of the release, so the docs' results can be checked against it without any of
these. A few checks rest on earlier runs, which docs/benchmark.md labels and lists in its opening: among them the Hi-C
depth series, the CpG line's kinetochore details (against DiMeLo-seq), its misplaced-methylation controls and its survey
of 83 vertebrates, and the shuffled-LTR and shuffled-GC controls.
