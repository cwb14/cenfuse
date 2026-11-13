# Future development

Open problems after cenfuse 0.15.0, for planning the next version (not user documentation). Each entry gives the problem as it
stands, the evidence, what would settle it, and where in the code to start. Entry numbers are stable (the code, tests and docs cite
`future.md #N`); a closed entry keeps its number ([Closed entries](#closed-entries)); the next new entry is 73.

## How to read this file

**Status.** *open*; *partly done*; *built, unchecked* (in the code, untested by truth); *tested, not adopted* (reopen only with new
evidence); *waiting for data*; *user's decision* (an outward-facing step).

**Words.** *Arm*: one benchmark run of a genome with one input set; *chromosome-run*: one chromosome in one arm. *Release runs*: every
arm rerun for a release, most replaying a stored satellite search (`--satellite-dir`). *The Hi-C test*: the gate. A line *decides*
when it may replace or move a call (design.md). *Confident*: `high` or `medium`. *Headroom*: chromosome-runs where some line's own
locus is right but the final call is not (`benchmark/headroom.py`).

**Where things live.** In the repo: `cenfuse/`, `tests/`, `docs/`, `benchmark/` (truth manifest, frozen scorer, table chain
`benchmark/release.sh`, doc templates in `benchmark/docs_src/`). README.md and docs/benchmark.md are generated from those templates
(`benchmark/fill_docs.py`) and constants.md from the code (`python -m cenfuse.consts`): edit the sources. Outside the repo, under
`/data2/chris/centromere_benchmark/` (`$CENFUSE_BENCH_DATA`): truth files, stored runs and analysis scripts, in the nested development
folders `v2/v3/.../v14` (v2-v7 = 0.1.0 development, v8 = 0.8, v9 = 0.9, v10 = 0.10.0, v10.1-v10.3 = 0.10.1, v10.4 = 0.11, v11 =
0.12, v12 = 0.13, v13 = 0.14, v14 = 0.15 with the repo at `v14/cenfuse`; cenrep and cenhic at the top). 0.15's runs:
`v14/runs/rel15c/<arm>/<id>`; its record: `v14/memo.md`. The polyploid panel: `/data2/chris/polyploid_genomes_v2/` (#68).

## Rules for any change

An entry is done only when these hold (design.md, Principles).

1. **Write the check first** (the release PLAN.md): what should change, on which genomes, and what must not.
2. **Give every constant a reason**: a convention (P ≤ 0.01), a majority, a value measured in the data, or a named heuristic. No
   threshold options; never a value chosen because it matches known answers.
3. **A rule prompted by one case** needs a reason independent of that case, a recheck on cases not used to make it, and a row in
   benchmark.md's table of rules changed after a result was seen. A rule built on one genome waits for a second, or stays a named
   heuristic that makes no confident call (#70, #72).
4. **No data are held out** (user's decision, 2026-11-08): every benchmark genome, the 13 tomatoes, *S. verrucosum*, *Ae. geniculata*
   and the polyploid panel are development data. Rule 3 does the work that holding data back did.
5. **Replay on the polyploid panel too** (#68). Homoeolog agreement has no truth, but it caught a bad 0.15 rule the benchmark missed.
6. **Keep the controls silent** (`cenfuse control`; the methylation one has no generator, #34): randomized Hi-C (57 runs) passes
   no test and changes no call; shuffled LTR
   positions (270 shuffles of 54 genomes) give no signature; holocentric genomes (196 chromosome-runs) get no confident call;
   misplaced methylation passes no test; GC-block shuffles pass the AT test at most at its 5% level (3.9% of 3,336).
7. **Do no harm**: no right call becomes wrong and no confident error is added, or each change is listed (in the release record and commit message).
   Right, then unsure (`low`), never confidently wrong.
8. **Report a missed criterion as written**; do not reword it.
9. **Test every rule** on a small synthetic fixture. Develop the LTR line on LTRquest annotations only.
10. **On the shared server**: stay within the thread budget (100), detach long jobs (`setsid nohup`), never edit a running script in
    place, and find your own jobs by full path, never with `pkill -f` on a pattern.

## Where to start

0.15.0: 2,236 of 2,624 chromosome-runs with truth right in the real arms (85.2%; 0.14: 82.9%; confident errors unchanged). Taking, on
each, a line whose own locus is right would reach 91.3%; most headroom is in AT troughs (91 of 161, nearly all *C. neoformans*) and
competing clusters (30). Confident misses: 12 in the release runs ([benchmark.md](benchmark.md#every-confident-miss)), and 16 with `tes/` in the genomes added
since 0.14 (2 in *Ae. geniculata* and *S. verrucosum*, 14 in the 13 tomatoes). In order:

1. **#27, #46, #72**: wrong LTR cores in non-CRM genomes (tomatoes: cores right on 117 of 156; 14 of 72 confident calls wrong).
2. **#35**: the public sets (Solanaceae with Hi-C or ChIP, wheat, oat, soybean lines, cattle, gibbon, a same-animal horse).
3. **#45**: calls at the edge of a wide LTR core (5 confident errors); tie-break the core's maximum first.
4. **#28, #39, #47**: *C. neoformans*, whose troughs are right but fail the test and whose LTR elements have no lineage.
5. **#37, #38, #50**: the 7 mammal misses of the 12 (6 scored against other animals' ChIP or another caller).
6. **#63**: the remaining large functions.
7. **#32**: the own-trans-pair rule, once a second genome shows the case.
8. **#60, #65**: CI's first result.

## All entries

Closed, one line each under [Closed entries](#closed-entries): 11, 12, 17, 18, 19, 20, 21, 23, 24, 25, 26, 36, 59, 64.

| # | area | problem | status | next step |
|---|---|---|---|---|
| [1](#1-calls-on-chains-of-partial-copies) | satellites | calls on chains of partial copies | partly done | CENH3 where a chain dominates without a CRM core (cotton) |
| [2](#2-the-40-family-cap) | satellites | the 40-family cap | open | a genome where it costs a call |
| [3](#3-satellites-specific-to-subgenomes-vs-the-honesty-cap) | satellites | satellites specific to subgenomes vs the honesty cap | tested, not adopted | CENH3 for an allopolyploid with such satellites |
| [4](#4-the-decisive-crm-test-with-subgenome-specific-lineages) | LTR | the decisive CRM test with subgenome-specific lineages | waiting for data | CENH3 in the named species |
| [5](#5-a-faster-search) | speed | a faster search | partly done | every arm searched again; the cost of real hits |
| [6](#6-faster-seed-clustering) | speed | faster seed clustering | open | a two-pass BLAST |
| [7](#7-the-fastan-crash-reported-upstream) | engineering | the FasTAN crash, reported upstream | user's decision | file the report |
| [8](#8-large-unplaced-scaffolds-and-dot-chromosomes) | inputs | large unplaced scaffolds and dot chromosomes | partly done | the other edge; organelles; known karyotypes |
| [9](#9-cenh3-truths-for-allopolyploids) | benchmark | CENH3 truths for allopolyploids | partly done | the public allopolyploid sets (#35) |
| [10](#10-the-ltr-retrotransposon-line) | LTR | the LTR retrotransposon line | partly done | a per-rule switch; `weak_core`'s cost |
| [13](#13-bird-microchromosomes-in-hi-c) | Hi-C | bird microchromosomes in Hi-C | waiting for data | deeper chicken Hi-C |
| [14](#14-the-rabl-axis-where-the-fold-has-no-power) | Hi-C | the Rabl axis where the fold has no power | partly done | partners on confident calls; Rabl tissue |
| [15](#15-centromeres-at-the-edge-of-a-terminal-unmappable-array) | Hi-C | centromeres at the edge of a terminal unmappable array | partly done | chicken; the sheep X call |
| [16](#16-the-chg-line-has-changed-no-call) | methylation | the CHG line has changed no call | open | same-plant CHG with CENH3 (data exist) |
| [22](#22-pa1a-in-other-p-annua-plants) | plants | Pa1A in other *P. annua* plants | waiting for data | CENH3 or deeper reads from those plants |
| [27](#27-wrong-ltr-cores-that-look-right) | LTR | wrong LTR cores that look right | open | an independent line; non-CRM genomes beyond *Solanum* |
| [28](#28-where-at-troughs-do-not-work) | AT | where AT troughs do not work | open | a second feature that *C. neoformans* passes |
| [29](#29-what-the-contact-crossing-a-centromere-is) | Hi-C | what the contact crossing a centromere is | waiting for data | the gap artifact; mutant or cell-cycle Hi-C |
| [30](#30-a-hi-c-test-with-too-few-confident-references) | Hi-C | a Hi-C test with too few confident references | built, unchecked | a genome where it acts, with truth |
| [31](#31-a-centromere-that-left-its-satellite-and-its-ltr-families) | LTR | a centromere that left its satellite and its LTR families | built, unchecked | CENH3 of *P. annua* |
| [32](#32-calls-in-a-terminal-eighth-that-an-interior-line-contradicts) | Hi-C | calls in a terminal eighth that an interior line contradicts | open | a second genome for the own-trans-pair rule |
| [33](#33-cis-and-trans-at-one-locus-as-a-confident-call) | Hi-C | cis and trans at one locus as a confident call | open | the AT trough's finer interval |
| [34](#34-the-benchmark-cannot-be-rerun-from-the-repository) | benchmark | the benchmark cannot be rerun from the repository | partly done | the scripts still outside; a fresh run of every arm |
| [35](#35-truth-on-disk-or-public-not-yet-used) | benchmark | truth on disk or public, not yet used | partly done | fetch and score the public sets |
| [37](#37-scoring-rules-and-controls-to-tighten) | benchmark | scoring rules and controls to tighten | partly done | the end-conditioned chance model; full-depth controls |
| [38](#38-which-array-holds-the-kinetochore) | satellites | which array holds the kinetochore | open | settle the orangutan truth |
| [39](#39-satellite-arrays-below-the-engines-absolute-sizes) | satellites | satellite arrays below the engine's absolute sizes | open | a family-free tandem readout, tested one per chromosome |
| [40](#40-a-centromeric-satellite-smaller-than-a-shared-familys-cluster) | satellites | a centromeric satellite smaller than a shared family's cluster | tested, not adopted | a readout that does not need the largest array |
| [41](#41-right-calls-held-at-low-by-the-honesty-cap) | satellites | right calls held at low by the honesty cap | open | a stability test of dominance |
| [42](#42-the-one-per-chromosome-test-near-its-level) | satellites | the one-per-chromosome test near its level | open | potato with genes; oat |
| [43](#43-assemblies-that-lack-or-misassemble-the-centromeric-array) | satellites | assemblies that lack or misassemble the centromeric array | open | the 56 repaired ape centromeres as positives |
| [44](#44-holocentric-genomes-one-low-call-per-chromosome) | satellites | holocentric genomes: one low call per chromosome | open | score the satellite holocentromeres with CENH3 |
| [45](#45-confident-calls-at-the-edge-of-a-wide-ltr-core) | LTR | confident calls at the edge of a wide LTR core | open | tied maxima; agreement judged by chance |
| [46](#46-crm-decides-by-its-name) | LTR | CRM decides by its name | partly done | the same test for CRM; tobacco, eggplant |
| [47](#47-centromeric-retroelements-the-ltr-line-cannot-read) | LTR | centromeric retroelements the LTR line cannot read | open | fungal lineage labels; non-LTR families |
| [48](#48-how-much-the-ltr-line-depends-on-the-annotation) | LTR | how much the LTR line depends on the annotation | open | rerun the thinning on current code |
| [49](#49-hi-c-cannot-see-inside-megabase-arrays) | Hi-C | Hi-C cannot see inside megabase arrays | open | multi-mapped reads or array-unique k-mers |
| [50](#50-satellite-free-centromeres-in-mammals) | Hi-C | satellite-free centromeres in mammals | open | same-animal horse data |
| [51](#51-hi-c-tests-at-their-threshold) | Hi-C | Hi-C tests at their threshold | open | the depth series on current code; the bin |
| [52](#52-hi-c-alone-untested) | Hi-C | Hi-C alone, untested | partly done | an effect-size guard; the read-class split |
| [53](#53-where-the-cis-trans-line-passes-every-extreme-is-trusted) | Hi-C | where the cis-trans line passes, every extreme is trusted | tested, not adopted | the diagnostic on gate-passing genomes |
| [54](#54-the-support-tiers-evidence) | Hi-C | the support tier's evidence | open | score its raises against functional truth |
| [55](#55-cpg-kinetochores-outside-primates) | methylation | CpG kinetochores outside primates | open | HGSVC3 CDRs; a mouse negative |
| [56](#56-telomeres-and-rdna-arrays-the-benchmark-does-not-count) | features | telomeres and rDNA arrays the benchmark does not count | partly done | single-end telomere arrays; a 2-copy rDNA rule |
| [57](#57-per-plant-verdicts-that-depend-on-other-chromosomes) | plants | per-plant verdicts that depend on other chromosomes | open | a leave-one-call-out check |
| [58](#58-report-wording-and-undocumented-outputs) | reporting | report wording and undocumented outputs | partly done | one flag-to-sentence table |
| [60](#60-reproducibility-across-environments) | engineering | reproducibility across environments | partly done | the macOS difference; CI's first result |
| [61](#61-memory) | engineering | memory | open | profile the rhesus peak |
| [62](#62-hi-c-engine-speed) | engineering | Hi-C engine speed | partly done | the library loop; a 2n = 80 polyploid timed |
| [63](#63-the-integration-cascade-and-other-large-functions) | code | the integration cascade and other large functions | partly done | the report's sentence functions, `Run._run`, `hicfuse.py` |
| [65](#65-tests-that-pass-without-running) | code | tests that pass without running | partly done | a slow test at production values |
| [66](#66-statements-the-013-outputs-no-longer-bear-out) | docs | statements the 0.13 outputs no longer bear out | partly done | the old figures; the audit as a release step |
| [67](#67-fragments-and-solo-ltrs-as-a-core-locator) | LTR | fragments and solo LTRs as a core locator | tested, not adopted | a larger gain than the age choice's |
| [68](#68-homoeolog-agreement-on-the-polyploid-panel) | benchmark | homoeolog agreement on the polyploid panel | built, unchecked | a panel species with CENH3 truth |
| [69](#69-why-a-null-test-could-not-be-run) | Hi-C | why a null test could not be run | open | a distinct flag for the degenerate case |
| [70](#70-the-at-trough-choosing-between-two-satellites-rests-on-one-genome) | AT | the AT trough choosing between two satellites rests on one genome | open | a second genome with two satellites per centromere |
| [71](#71-age-choice-and-family-composition) | LTR | age choice and family composition | open | the age test within one family |
| [72](#72-minority-cores-beyond-solanum) | LTR | minority cores beyond *Solanum* | open | non-CRM genomes with truth outside *Solanum* |

## Benchmark, truth and validation

### 9. CENH3 truths for allopolyploids

Partly done. One allotetraploid has CENH3 truth on its own assembly: *Ae. geniculata* (Yunusbaev 2026, Table S2, Mb resolution;
`polyploid_genomes/_cenrep_batch/review/Aegilops_geniculata.cenh3.tsv`); it shaped two cenrep LTR rules. 0.15: 0 of 14 from the FASTA
(all `low`), 13 of 14 with `tes/`, 1 confidently wrong (chr05U, #45). Open: more allopolyploids, chosen to settle #1, #3, #4 and #10
(wheat T2T, *Avena sterilis*, *A. maroccana*, *Fragaria × ananassa*: #35); provisional subgenome truth carried from a progenitor by
unique flanks (v6 `bench/flank_truth.py`: 13 of 14 per subgenome in a *Cryptococcus* hybrid), valid only where centromeres did not
move after polyploidization.

### 34. The benchmark cannot be rerun from the repository

Partly done. In the repo since 0.14: the truth manifest (`benchmark/truth.tsv`), the frozen scorer as a tested module, the table chain
and templates (`benchmark/release.sh`; it reproduced 0.13's tables byte for byte), the controls (`cenfuse control`), and a stop when
the imported cenfuse is not the code under test. Still outside, with hard-coded `/data2/chris/...` paths:
- 0.15's doc chain (`v14/bench/docs15.sh`): after the repo chain it scores the genomes added since 0.14 and writes their sections
  (v13 `bench/score_new.py`, `md_014.py`, `tom14_md.py`, `tom14_te.py`, `cmp_runs.py`, `telomere_recall.py`; v14 `tom15_detail.py`,
  `md_015.py`). `benchmark/paths.py` still defaults to 0.13.0's runs.
- The panel check (v14 `bench/polyeval.py`, `replay_T.sh`, `run_P.sh`); the truth-making scripts (`map_chip.sh`, `chip_truth.py`,
  the fetch scripts); the analysis scripts the docs cite (v11 `bench/pa1a_age.py`, `young_core.py`, `poly_one.py`); the
  misplaced-methylation control's generator (never found).

Release runs replay searches mostly made before 0.11 (#5); fresh unseeded runs equal the replays where tried (0.15: *V. dahliae*, *M.
sympodialis*, `v14/runs/fresh15c`). **Settle:** a clean checkout regenerates benchmark.md from the stored runs; a fresh run of every
arm reproduces the published counts or lists each difference.

### 35. Truth on disk or public, not yet used

Partly done (*Ae. geniculata*, *S. verrucosum* and 13 tomatoes scored since 0.14). Fetch in table order; confirm each antibody and
target first (the horse CUT&Tag target and the pig antibody are unconfirmed); score with the frozen ChIP rule (#37). Sheep and goat
CENP-A ChIP gives no truth ([Tried and set aside](#tried-and-set-aside-not-in-designmd)).

| data | what it would test | where it is recorded |
|---|---|---|
| tobacco, eggplant CENH3 (accessions to find) | CRM that is not centromeric (#46) | cenrep's literature notes |
| wheat Chinese Spring (T2T GCA_040256815.2; ChIP SRR9089557/8, DRR821424/5) | a hexaploid; CRM plateaus (#45) | `v8/data/published/manifest_ena_runs.tsv` |
| *Avena sterilis* (GWHERDU00000000; CRR1339975/6); *A. maroccana*, *A. atlantica*, *A. eriantha*, 'Starter' (CRA007177) | oat's non-centromeric satellite on every chromosome (#42); allopolyploids | `v8/data/published/manifest_gsa_runs.tsv` |
| 6 accessions of Shi et al. 2026 Table 8 published without deposited ChIP | Solanaceae LTR rules (#46, #72), if how their centromeres were found is confirmed | Shi 2026, Supp Table 8 |
| 4 more soybean lines (ChiHei23, QiHuang34, ZH13, Williams82); the pan-centromere of 27 accessions (PNAS 2023, Dataset S1, on disk, 0.1-0.5-Mb resolution; needs each accession's 2020 pan-genome assembly) | CRM plateaus (#45); centromeres that move within a species (#31) | `v8/data/published/download_manifest_cenh3.tsv` |
| cattle CENP-A CUT&RUN and input (SRR28462359/8) | an all-acrocentric ruminant; the which-end rule (#15) | v2 `data/research/research_truth.md` |
| hoolock gibbon CENP-A CUT&RUN | the siamang's failure mode | same |
| siamang same-cell-line CENP-A DiMeLo-seq (bigWig, 15.1 GB) | the siamang's 2 confident disagreements; `low` references (#30) | same |
| horse: same-animal assembly (GCA_052323815.1), Hi-C (SRR36069624), CUT&Tag/ChIP with IgG (SRR36096936-9, PRJNA1248178); CENP-A ChIP of about 15 horses (96 runs) | horse chr11 (#50); the horse Hi-C that fails; centromere sliding | same |
| pig CENP-A ChIP and input of the assembly's breed, Duroc liver (DRR1139646-9; at GSA PRJCA032804 and DDBJ only; antibody not stated) | pig's 3 confident misses (#37) | same |
| the repaired ape assemblies (56 centromeres misassembled in v2.1, with repair BEDs), both haplotypes | about twice the ape truth; misassembly as a positive control (#43) | same |
| HGSVC3 CDRs (67 genomes, 2,156 centromeres; Gao 2026) | the CpG line (#55) | figshare 10.6084/m9.figshare.29921990 |
| holocentric CENH3: *R. breviuscula* ChIP (SRR27965399/8), *Chionographis japonica* (GSM7120336), *Luzula sylvatica* (358 domains) | the satellite line on satellite holocentromeres (#44) | v4 `data/research/hard_genomes.md` |
| same-plant CHG with CENH3: *S. verrucosum* CPC54, potato DM FAST5 (about 625 GB), rice T2T-NIP modBAMs | the CHG line (#16) | same |
| Hi-C: *Chlorella* NS4-2 (CRR559578), the primary Bd21 library (CRR611509, 289 M pairs), deeper libraries for the shallow controls (#37) | Hi-C on LTR cores in an alga; the Rabl axis in a grass (#14) | v10 `data/hic/hic_runs_plants.tsv` |

### 37. Scoring rules and controls to tighten

Partly done (0.14: call widths, domain widths and midpoint distances reported, `benchmark/precision.py`; the scorer's tolerances
stated; the truth rule at 3×, 5× and 50-kb bins moves no outcome, v13 `bench/truth_sensitivity.py`). Open:
- **The chance model** draws calls uniformly along the chromosome, optimistic where centromeres sit near a known end (goat, sheep,
  horse, mouse). Condition it on end distance.
- **250 kb without Hi-C** for fungi and *Plasmodium*, whose centromeres span 1-100 kb: *S. cerevisiae*'s AT troughs are "right" on 12
  of 16 that way and overlap on 2. Score small genomes without Hi-C at a fine fixed tolerance. Barley is scored at 8 Mb without Hi-C.
- **Weak truth.** Pig and horse ChIP come from other animals; pig's 3 confident misses lie where that ChIP's second peak (3.0-4.4×)
  falls, and the karyotype agrees with the calls (chr11 bi-armed but its domain at CI 0.002; chr17 acrocentric but at CI 0.24).
  Siamang, sheep, goat, rhesus and giraffe are scored against another caller. Bornean orangutan CM054652.2 may be haplotype-specific
  (#38); *N. crassa* LG VI carries an inversion (Galazka 2016); *C. tropicalis* CENP-A comes from an older assembly; *V. dahliae*'s
  authors inferred its centromeres from this Hi-C; *C. neoformans* counts twice in the 47-genome calibration. The ChIP truth rule
  (100-kb bins, ≥ 4×) has no stated reason for either number; in cucumber 1 of 7 chromosomes passes it.
- **Shallow controls** (library prefixes): *R. breviuscula* 3.3 M pairs, *B. mori* 5.9 M, *Ae. tauschii* 0.84 M (another accession),
  sorghum 3.9 M, mouse 15.7 M; "no holocentric genome passes" rests on them (#51). Deeper: *B. mori* SRR25988387 (369 M), *R.
  breviuscula* SRR17429201 (863 M), soybean SRR37686999 (415 M), *B. rapa* CC1 SRR8633037 (368 M), Col-0 SRR2626163 (143 M).
- **Attribution** (user's decision): benchmark.md credits the ape CDRs to Logsdon et al. 2024; the v2 notes record an unpublished
  curated set (v1.3, "contact the lab before publication"), and the ape DiMeLo-seq and CUT&RUN tracks as unpublished.
- **Settle:** each scorer change reported beside the old numbers for one release; no README claim rests on overlap alone (an animal
  call is about ten times its CENP-A domain).

### 68. Homoeolog agreement on the polyploid panel

Built in 0.15 as a check, no truth. Each of the 106 allotetraploids (`/data2/chris/polyploid_genomes_v2/<species>`: `cenfuse_in/`, a
stored run, `homoeologs/<phasing>/*.anchors.tsv`) is replayed from classify (v14 `bench/replay_T.sh` with `FIX_STAMPS=1`, as the
stages carry 0.12 stamps; `bench/run_P.sh`). A call lifted through the collinear anchors agrees within 5% of the chromosome of its
homoeolog's (`bench/polyeval.py`; per rule `bench/review8.py`). 0.15's rules: 663 of 1,103 pairs agree over 105 species, 657 under
0.14's; of 26 changed pairs, 14 agree against 8. It caught a first form of a 0.15 rule that let a signature with unspecific cores
decide in *Arachis* (11.6-fold seed, 1.5-fold cores; agreement on the changed pairs 11 → 2 of 14). A disagreement can be a right call
on an evolutionarily new centromere. **Settle:** agreement calibrated against CENH3 on one panel species (#9).

## Satellites

### 1. Calls on chains of partial copies

Partly done. A family can concentrate at one place without being tandem (under half its bp in tandem, under half in whole-unit hits):
a chain of pieces of a longer repeat, usually a retrotransposon (*Glycine dolichocarpa*: 44 whole hits of 25,732, 65% inside SIRE;
10-Mb calls). Such a call is at most `medium`, `low` unless a CRM core (or a core of a lineage the family is not made of) overlaps it,
and a narrower decisive core inside it becomes the call (`chain_to_core`; design.md). The cap rests on a principle, not on truth. On
the panel 41 such families in 29 species carry 409 of 2,355 calls (benchmark.md says 459: #66). No benchmark genome has a dominant
chain; its chain calls were 29 right, 4 wrong of 33 with truth (0.12), `chain_to_core` 7 of 7. *Ae. geniculata* with `tes/` elects a
1,216-bp chain family inside CRM cores: 13 of 14 right, none `high`. Never measured: whether short (< 1 kb) single-sequence families
lose diverged copies (`min_hmm_seqs` 20).
- **Settle:** CENH3 for a genome whose chain dominates without a CRM core (6 cotton species: 23-26 of 26 chromosomes). **Code:**
  `sat/classify.py` (roles, LOCUS_FOLD), `sat/centromere.py` (`partial_chain`, `chain_to_core`).

### 2. The 40-family cap

Open. At most 40 clusters per complex are modeled; the dropped ones are listed (`satellite/dropped_clusters.tsv`). Over 71 runs that
dropped any: 2,376 clusters, 1,553 on one chromosome; only 4 reach 20 kb on a chromosome without a call (all *B. carinata*), 349 lie
on a chromosome with a `low` call. No case yet where the cap costs a centromere. Changing the family set is not neutral: overlapping
hits go to the higher-scoring family (`resolve_overlaps`; 0.8: moss's fam003 became centromeric and Chr08's call moved). If a genome
shows a cost: a budget priced from `search_cost.tsv` instead of a count.
- **Settle:** *B. carinata* without the cap (roles, main families, arrays, calls). **Code:** `sat/families.py` (`max_families`).

### 3. Satellites specific to subgenomes vs the honesty cap

Tested, not adopted: dominance judged within subgenomes lifts no more capped calls than random splits do (design.md, Choices tested
and rejected). Reopen with CENH3 for an allopolyploid with subgenome-specific satellites (*Fragaria × ananassa*, Jin 2025: #35).

### 38. Which array holds the kinetochore

Open. Where a chromosome holds several candidate arrays, the sequence line takes the largest cluster; no tested line can choose
another.
- **Bornean orangutan CM054652.2** (`high`, 14.3 Mb off; the only confident ape miss): CENP-A on 4 small alpha-satellite arrays (73
  kb), the call on a 3.4-Mb cluster. Hi-C's peak lies 0.26 Mb from the CDR but is `low` (1 of 3 engine tests) and cannot oppose
  `high`; the called array dips in CpG too. The truth may be haplotype-specific (CenSat v2.0 marks the called array active; the CDR
  region was redrawn in 2025; hap2's DiMeLo-seq shows CENP-A at both sites).
- **Bonobo CM055489.2** (`low`, wrong): the true array is the reported competing cluster, but only a strong Hi-C peak may choose it
  (of the two weaker peaks that rule kept out, this one would have been right and CHM13 chr21's wrong).
- **Mouse**: since 0.15 the AT trough chooses the minor satellite (`at_chose_alt`, 18 of 20 from the FASTA, `low`; #70). The rule
  that makes it the alternative was shaped on mouse alone. Chr2, 5 and 14 share one byte-identical end (9,259-807,536).
- **Approach.** Settle the orangutan truth (2025 CDR, CenSat v2.0, both haplotypes). A second line for a weak Hi-C peak against
  `high` (a CpG dip at the peak; CpG alone chose badly: design.md and [Tried and set
  aside](#tried-and-set-aside-not-in-designmd)). For two satellites at one end, a sequence feature of the CENP-A satellite (CENP-B
  boxes in mouse; an idea).
- **Settle:** the orangutan, bonobo and mouse calls; no right ape call lost; CHM13 chr21 stays right. **Code:** `sat/centromere.py`
  (`satellite_call`, `rival_bp`), `fuse.py` (`hic_chose_alt`, `hic_alt_ambiguous`), `cdr.py`.

### 39. Satellite arrays below the engine's absolute sizes

Open. One prior, S_MIN = 20 kb, scales the satellite engine's size rules in every genome (the "prior" fields of `sat/params.py`, 5 kb
to 5 Mb), and units under 30 bp or over 30 kb are not searched (`min_period`, `LU_MAX_UNIT`). So the satellite line places 0 of 14
*C. neoformans*, 0 of 8 *V. dahliae*, 0 of 21 *Z. tritici*. In 0.10.1's comparison (optimistic: readouts picked by the truth) TRASH2
placed 7 of 14 *C. neoformans* with 0 wrong (a 225-bp unit, one ~2.2-kb array in 7 centromeres only), and satellome's largest tandem
region per chromosome 10 of 14 *C. neoformans* and 7 of 8 *V. dahliae*.
- **Approach.** A family-free readout tested like the AT troughs (each chromosome's densest tandem region, accepted where the N
  densest fall one per chromosome beyond chance); size floors scaled to the genome from a reason in the data.
- **Settle:** *C. neoformans* placed from the FASTA; no change where satellites are large; shuffled and holocentric genomes fail;
  the TRASH2/satellome comparison rerun on current code. **Code:** `sat/params.py`, `atline.py`, `sat/karyotype.py`.

### 40. A centromeric satellite smaller than a shared family's cluster

Tested in 0.15, not adopted. A family is centromeric only where it forms a chromosome's largest centromeric cluster, so a
chromosome-specific satellite smaller than a shared family's cluster is never called: potato's centromeres sit on kb-unit families
(St24, 981 bp, 299 kb at Cen1) while a 157-bp family dominates, so potato is 2 of 11 from the FASTA. Taking the largest cluster of
families on at most 3 chromosomes places 2 of 11 potato, 0 of 12 tomato and 3 of 12 *S. verrucosum* (16 of 33 chicken, 4 of 17
*Chlamydomonas*; v14 `bench/sat40.py`): these satellites are rarely the chromosome's largest array. Reopen with a readout that does
not need the largest array (#39). **Code:** `sat/classify.py`, `sat/karyotype.py`.

### 41. Right calls held at low by the honesty cap

Open. Where no family dominates, every satellite call is `low`, even when mostly right: chicken 27 of 33 (chromosome-specific
satellites), moss 16 of 26 from the FASTA, *Chlamydomonas* 6 of 17. Lifting the cap where the one-per-chromosome test passes was
rejected (chicken +27 right, +5 wrong; pig 8 and 9). Dominance itself is unstable (0.8's fresh run split pig's satellites differently:
FASTA-only pig 13 → 8 of 17), and nothing measures how far a small change of the family set moves it. A bird lead: 22 of chicken's 50
interstitial telomeric arrays lie within 100 kb of a CENP-A window (1.4 expected, P 4e-21; Huang 2023); as a line only chicken beats
chance (6 right, 0 wrong), so it needs a second bird T2T genome.
- **Settle:** a stability test of dominance (resampled seeds, families pooled by kind); chicken's right calls rise and no wrong one
  does. **Code:** `sat/centromere.py` (`no_dominant_satellite`), `telomere.py`.

### 42. The one-per-chromosome test near its level

Open. The test passes on any family that happens to lie once per chromosome (design.md, Known risks), and genomes near its level flip:
- **Potato** fails from the FASTA (7 of 12, P 0.87) but passes with `genes/` (10 of 12, P 0.044), whose gene veto changes the family
  set, so two wrong calls become confident (chr04 3.4 Mb off, chr08 1.7 Mb off): an input that should only help adds 2 errors.
- ***S. chilense*** passes (10 of 12, P 0.045) on mostly non-centromeric arrays; only the honesty cap keeps its calls `low`.
- **Oat**: a 335-360-bp family on every chromosome is not centromeric (cenrep 0.8.9) yet would pass (#35).
- **Siamang** fails (18 of 26, P 0.2): its arrays form two clusters over 500 kb apart on 8 chromosomes.
- **Few chromosomes**: no power at 4 or fewer; the acrocentric test needs about 16; no small acrocentric karyotype tried.
- **Approach.** Report the test's P beside each call that depends on it; measure how the outcome moves with the family set (with and
  without the gene veto; families pooled by kind). **Settle:** no confident error in potato with `genes/`; oat's FASTA-only calls
  against CENH3. **Code:** `sat/karyotype.py`, `sat/genomes.py` (the gene veto).

### 43. Assemblies that lack or misassemble the centromeric array

Open. Where the array is missing, a non-centromeric satellite elsewhere can be called with confidence (cenrep 0.8.9: *B. napus*
A06/A09, `high`); gaps are reported (`call_gaps`) but change nothing. The 56 ape centromeres misassembled in NCBI v2.1, excluded from
the truth and with published repair intervals, are free positive controls. A collapsed array makes spanning pairs farther apart than
their coordinates, mimicking the cis dip (#29); mouse chr2, 5 and 14 share one templated end (#38).
- **Approach.** Report per call what an incomplete centromere looks like (a gap in or beside the call, a short array, a call away from
  its karyotype position). **Settle:** the ape repairs flagged; T2T assemblies not. **Code:** `morphology.py`, `interpret.py`.

### 44. Holocentric genomes: one low call per chromosome

Open. Holocentric genomes get no confident call (0 of 196 chromosome-runs) but still one `low` "centromere" per chromosome, and
nothing says the genome may be holocentric. A flag from `no_dominant_satellite` was rejected (also on 25 of 26 siamang, 25 of 32
EquCab3 calls); the one-per-chromosome test fails there (P 0.13-0.73) and in many monocentric genomes. Never scored: the satellite line
on satellite holocentromeres (*R. breviuscula*, *Chionographis japonica*, *Luzula sylvatica*: #35). An LTR signature can form at the
ends (*R. breviuscula*: Athila, 5.8-fold, cores within 1 Mb of an end on 4 of 5); Hi-C there is guarded, but a non-CRM signature
above 10-fold or a CRM seed would place calls at the ends from `tes/` alone.
- **Approach.** Several independent signs, each tested on the three genomes with CENH3: no dominant satellite, a failed
  one-per-chromosome test, Hi-C's holocentric signature, CENH3-domain-like arrays along the chromosome. **Settle:** the four
  holocentric controls flagged, no monocentric genome flagged. **Code:** `sat/signature.py`, `interpret.py`.

## LTR retrotransposons

### 4. The decisive CRM test with subgenome-specific lineages

Waiting for data. CRM core counts differ between subgenomes in 3 of 87 species at P ≤ 0.01 (0.9 expected), and the marker is decisive
in 2 of them anyway (*B. napus* A 8/8, B 2/8; *Euphrasia confusa* 10/11 vs 3/11); at P ≤ 0.05, 9 species (4.4 expected, P 0.03), 5
still decisive (*Cerastium alpinum* only flags: 11/18 vs 0/18). No evidence yet that a decisive marker is often lost. **Settle:** those
species against CENH3 (#35), random splits of diploids as the null; then test CRM concentration within subgenomes.

### 10. The LTR retrotransposon line

Partly done. Open questions on the branches that move, overrule or lower a call (design.md, What a decisive core does):
- **Which rule earned which gain.** Rejected families are reported (`seq_rejected_*`), but no benchmark-only switch per rule attributes
  the gains.
- **The head-to-head test** (a confident call is replaced only if the core is 10-fold enriched and holds significantly more of the
  signature's elements than the call; otherwise the call drops a level, `te_discordant`): its cited cases (*B. rapa* A09, soybean
  chr1/chr6) have not arisen since 0.13. In 0.15 it keeps two right calls with truth, grape chr11 (`high`, HT and HTH arms) and *Ae.
  geniculata* chr05M (`low`). Replacing a `low` call needs no such test (of the moves it would have blocked, 2 right, 1 wrong:
  soybean chr13, #27).
- **`weak_core`** (a thin core of a decisive signature takes a chromosome whose satellite failed, `low`): it cost *S. verrucosum* 2
  right satellite calls (cores of 2 and 4 elements, chr10 and chr12) and is right on 7 of 15 tomato chromosomes. No confident error,
  but right `low` calls lost. **`weak_satellite_beside_core`**: one confident error (*S. verrucosum* chr06, 1.4 Mb), tomatoes 7 of 7.
- **`te_concordant`** is a weak flag (#45).
- **Settle:** the per-rule switch replayed on every `tes/` arm and on the panel (#68), each branch's right and wrong counts reported.
  **Code:** `sat/centromere.py`, `sat/signature.py`.

### 27. Wrong LTR cores that look right

Open. A signature learned from the other chromosomes can mark where centromeric families accumulated beside the centromere, and no
statistic of the signature tells those cores from the right ones.
- **Tomatoes.** Under 0.15's rules the cores are right on 117 of 156 chromosomes of the 13 tomatoes (102 of 151 under 0.14's) and the
  benchmark tomato's wrong cores fell from 5 to 4. With `tes/`, 14 of 72 confident calls are wrong, 5 of the 18 that merely agree with
  the core (benchmark.md, Thirteen tomatoes). Without Hi-C two benchmark-tomato misses stay `medium` (CM175074.1, CM175078.1); with
  Hi-C the cis peak, which is right, lowers them.
- **Others.** Soybean chr13: a CRM core 1.8 Mb from CENH3, `medium`. Potato chr08 (#45). *Arabidopsis*: Athila (13.3-fold) has its
  cores 0.46-1.38 Mb beside CENH3 on all 5 chromosomes and only flags because CEN180 passes the one-per-chromosome test. In 6 of 27
  polyploids above 10-fold (three *Sphagnum*, *Lotus subbiflorus*, *Malus toringoides*, *Solanum stoloniferum*) cores and FASTA-only
  calls mostly disagree, no truth.
- **Tried in 0.15, none separates them:** seed specificity, the core's share or fold, a second comparable peak (accuracy 0.85 with
  one, 0.83 without; reported as `core_alt`), a super-majority of the core's share; fragments as a locator (#67).
- **Leads.** A core beside the dominant satellite array rather than over it (*Arabidopsis*); Hi-C anchored on the cores; the cis
  signal tested by trans contacts where LTR cores exist too (not built in 0.13 for its machinery, not for a bad result: tomato with
  `tes/` and Hi-C would go from 10 to 12 of 12).
- **Settle:** non-CRM signatures with truth beyond *Solanum* (#35, #72); *Arabidopsis* with CEN180 masked must not move onto the
  Athila cores; randomized TNEG tomato silent. **Code:** `sat/signature.py`, `hicfuse.py`.

### 31. A centromere that left its satellite and its LTR families

Built, unchecked (no CENH3). *P. annua* Pa1A: earlier analyses (no ChIP) and Hi-C (peak 168.0 Mb; the test passes on 14 of 14) put
its centromere at 163.0-169.4 Mb, 21 Mb from its largest satellite, around which the CRM families learned from the other chromosomes
concentrate (145.1-148.3 Mb). Insertion age agrees with Hi-C: the signature's copies are younger at the Hi-C site (median 0 vs 620,000
years, P 2e-13; within one CRM family 0 vs 157,000, P 2e-7, v11 `bench/pa1a_age.py`), as young CRM colonize present maize centromeres
and older ones their flanks (Wolfgruber et al. 2009).
- **State in 0.15.** The age choice between density peaks (`signature.age_choice`; design.md) makes the young peak Pa1A's core, so the
  call is 163.0-169.4 Mb, `high`, Hi-C agreeing. The 0.13 tie-breaker (`ageline.py`, `moved_with_age`: younger copies at a validated
  strong Hi-C peak let Hi-C move a `high` call) now acts on no benchmark run and has never been scored against truth.
- **Limits.** The premise is tested on the other chromosomes (young copies in the cores more often than old: Fisher P 5e-10 here). The
  young copies' core alone is a worse locator than the density core (128 vs 155 of 172 with truth; v11 `bench/young_core.py`). A
  younger family can mimic a move (#71). The README counts Pa1A among cenfuse's strengths ("no ChIP yet"), ahead of the evidence.
- **Settle:** CENH3 ChIP or CUT&Tag of *P. annua*; genomes with known evolutionarily new centromeres and Hi-C (horse chr11, #50; the
  orangutan chr10 neocentromere polymorphism, de Gennaro 2026); within-species centromere tables (#35). **Code:** `sat/signature.py`
  (`age_choice`), `ageline.py`, `fuse.py`.

### 45. Confident calls at the edge of a wide LTR core

Open; 0.15's LTR rules do not reach these calls. A call agrees with a core when the gap is at most the wider of the two
(`te_concordant`), and an agreeing core lifts both the honesty cap and the one-per-chromosome cap. Of 130 such calls with truth in the
0.13 FASTA + `tes/` arms, 125 are right and 5 confidently wrong, 5 of the 7 confident errors of those arms. In 3 of the 4 below the
core holds CENH3 and the call sits at its edge:

| chromosome | confidence | call (Mb) | LTR core (Mb) | CENH3 (Mb) |
|---|---|---|---|---|
| maize Mo17 CM039154.1 (chr5) | `high` | 107.19-107.40 | 103.54-108.83 | 103.5-105.2 |
| *Ae. tauschii* Chr5 | `high` | 199.36-201.42 | 194.02-201.03 | 195.4-198.1 |
| *Ae. tauschii* Chr7 | `medium` | 340.27-340.50 (45 kb of satellite) | 340.47-349.60 | 343.6-344.9 |
| potato chr08 | `high` | 15.54-16.92 | 16.46-17.41 (5 elements) | 18.6-18.7 |

The fifth is tomato CM175074.1 (its core misses CENH3: #27); *Ae. geniculata* chr05U has the same shape (2.1 Mb off, in a 9.3-Mb core
whose other edge reaches CENH3). Maize chr5 is the plants' only confident error from the FASTA alone (no Hi-C). Exploratory: where the
core is at least 3 times wider than the call, 11 right and 4 wrong; otherwise 114 and 1.
- **Tested in 0.14, not adopted: agreement at the core's maximum** (within one TE bin, only where agreement lifts a cap; v13 PLAN.md
  C2, `results/core_peak/`): it lowers 2 of the 5 errors and 3 right calls (potato chr05, tomato CM175070.1, CM175077.1). Overlap with
  the peak bin: 4 right and 3 wrong lowered; any tied maximum: 2 and 1.
- **Ties first.** 66 of the 130 core profiles have tied maxima (small counts, a 5-bin box filter), 32 with a tie over a bin from the
  one `np.argmax` takes, so a core's position can depend on which equal maximum comes first. A tie-break with a reason (unsmoothed
  counts, a smooth kernel, the tied run with more elements) needs its own check before any rule reads the maximum.
- **Approach.** Judge agreement by chance, not by gap (`linetests.py` computes each chromosome's chance of agreement).
- **Settle:** the five calls fall to `low` or move onto CENH3; none of the 125 right calls drops below `medium`; CRM-plateau genomes
  not used to make the rule (wheat, the soybean lines: #35). **Code:** `sat/centromere.py` (`te_concordant`), `sat/te.py` (`agrees`),
  `linetests.py`.

### 46. CRM decides by its name

Partly done. CRM decides from 3 significant cores with no fold level; other lineages need 10-fold, specific cores (median core fold ≥
10, 0.15) and, to act on a confident call, a core with most of the signature's elements on its chromosome (`TE_MAJOR_SHARE`). CRM's
standing rests on every CRM core short of a majority lying on CENH3 (24 of 24), but those come from thinned annotations of genomes
whose CRM is centromeric (rice, potato, maize, *Brachypodium*). Non-centromeric CRM exists: tobacco (CENH3 binds Tekay; CRM clusters
beside), eggplant, soybean chr13 (#27); LTRquest labels 14 CHM13 ERVs CRM. 0.12's seed, which ranks a 3-core CRM candidate first, was
not retested against the *Arabidopsis* case that sank a v6 form ([Tried and set aside](#tried-and-set-aside-not-in-designmd)). The
reverse case, non-CRM lineages held back, eased in 0.15: a minority core replaces a call `low` on its own and stays `low`
(`minority_core`: right on 41 of 58 tomato chromosomes, none confident; #72).
- **Approach.** The same test for CRM as for other lineages (a fold level, or the decisive test measured per genome), or CRM's
  standing only where its cores pass the one-per-chromosome test. **Settle:** no right call lost in the CRM-decisive genomes; no
  tobacco or eggplant CENH3 domain overruled by CRM (#35). **Code:** `sat/signature.py` (`seed_of`, `TE_SEED_MIN_CORES`).

### 47. Centromeric retroelements the LTR line cannot read

Open. The LTR line reads classified LTR elements only. In *Cryptococcus* 89-94% of LTRquest's elements lie in the CENP-A domains, but
none has a REXdb lineage (unknown, chromo-outgroup), so the FASTA and `tes/` place 0 of 14 (Hi-C does). Pooling unclassified elements
(v6) qualified in fish genomes without truth (design.md); all LTR elements pooled per chromosome (0.15) are right on 6 of 14 *C.
neoformans* and 3 of 8 *V. dahliae* but never significant (too few independent insertions; the 3-site minimum is what a pooled test
would give up). *Chlamydomonas*' centromeres are clusters of ZeppL, a LINE bound by two CENH3 paralogs; LTRquest finds 669 elements,
613 without lineage; 6 of 17, all `low`.
- **Approach.** Lineage labels covering fungal Ty3 and chromoviruses ("chromo-outgroup" may be the centromeric ones); a non-LTR
  annotation read by the family signature; each needs a second genome with truth. **Settle:** *Cryptococcus* (both nuclear states,
  no Hi-C) and *Chlamydomonas* placed; fish and *Oryza* unchanged; shuffled positions silent. **Code:** `sat/te.py` (reads
  `LTR_retrotransposon` only), `sat/signature.py`.

### 48. How much the LTR line depends on the annotation

Open. The decisive test needs significant cores on most chromosomes; potato sits at its threshold (8 of 12 in 0.13; a majority is 7).
benchmark.md's annotation table says "measured on 0.11" but holds v6's lineage-line numbers (half of potato's CRM elements, EDTA or
DANTE_LTR each lose the gain from `tes/`; rice drops from 8 cores to 3, 1, 1 at 50, 25, 10%). The family signature (0.12 on) was
never thinned. LTRquest writes no version into its GFF3, so a run records only the `tes/` file's fingerprint. The README marks
`--proteins` untested. **Settle:** the thinning rerun on current code with LTRquest only (v6 `bench/td_thin.py`), reporting where each
signature stops deciding; the README's command benchmarked with `--proteins`, or the flag dropped. **Code:** `sat/signature.py`,
`teline.py`.

### 67. Fragments and solo LTRs as a core locator

Tested in 0.15, not adopted. The signature families' LTRs mapped back with minimap2 (leave one chromosome out) peak on CENH3 on 248 of
299 chromosomes with an element core, against 237 for the element core and 248 for the age choice (251 both; v14 `bench/frag4.py`):
+3 for a new external stage per genome, gaining in tomatoes and losing in potato, soybean and sorghum. Reopen with a larger gain.

### 71. Age choice and family composition

Open (0.15 code review). `signature.age_choice` compares the ages of all recruited families' copies at two peaks; families burst at
different times, so a peak made of a younger family reads as younger with no move of the centromere. In 0.15 it moved 27 cores
over 32 genomes with truth and `tes/` (12 wrong to right, 1 right to wrong; v14 `bench/chk_moved.py`); on the panel the 62 pairs it
touches agree 32 times, as before. **Settle:** the test within each peak's commonest family, where copies suffice; compare the moves.
**Code:** `sat/signature.py`.

### 72. Minority cores beyond *Solanum*

Open (0.15 code review). `minority_core` was built on the tomatoes (#46). Outside *Solanum* one scored chromosome: cucumber chr1 (HT,
HTH), right in 0.14, 2.1 Mb off in 0.15, still `low`, the only right calls 0.15 lost on the release arms. On the panel it acts on 14
calls in 3 species (*Arabidopsis suecica*, *Malus toringoides*, *Sphagnum papillosum*); the 11 pairs they touch agree in 7 against 3
under 0.14's rules. It makes no confident call. **Settle:** non-CRM genomes with truth outside *Solanum* (#35). **Code:**
`sat/centromere.py`, `sat/signature.py`.

## Hi-C

### 13. Bird microchromosomes in Hi-C

Waiting for data. Chicken's 2.5-7-Mb microchromosomes stay out of the Hi-C line because its library (the first 40 M read pairs of
SRR15423486; 26.1 M unique) is too shallow for bins under 250 kb. The bin cap is not the limit: at `MAX_CHROM_BINS` 8,000 (0.14) the
depth rule still chose 250 kb, the gate passed on 0 of 40 and no call changed. **Settle:** deeper chicken Hi-C, scored against the
CENP-A windows. **Code:** `hic/params.py`, `hic/bins.py`, `fold.py`.

### 14. The Rabl axis where the fold has no power

Partly done. The trans arm alignment passes in *C. neoformans* (anaphase, 12/14) and *M. sympodialis* (6/8), but the arm fold, which
must lean the same way first (the holocentric guard), has no power in these small nuclei (1/14; 0/8), so the alignment raises nothing.
In holocentric *C. elegans* the alignment reached 2/5 (P 0.095), and a position-matched null that rejects it also rejects barley. G1
*C. neoformans* loses the cluster, yet anchored on CENP-A 8 of 14 chromosomes still align (P 0.0014). Other lines now place these
genomes (`at_pair`, the cis-trans line, `at_denovo_pair`, Hi-C alone: #52); with 6 chromosomes *C. tropicalis*' Rabl test cannot reach
P ≤ 0.01 (floor 0.0156).
- **Approach.** Partners placed on the confident sequence calls as a corroboration test (prototype v10.1
  `explore/anchored_align.py`: G1 8/14, anaphase 12/14, randomized 0/14), with a per-genome test like the gate's; a
  compartment-adjusted alignment (the E1 sign at centromeres is opposite in plants and apes, cenhic 0.2); the fold at the 5-kb rung.
- **Settle:** barley, *Neurospora* and the yeasts pass; *C. elegans*, *Rhynchospora*, *Bombyx* and randomized Hi-C fail. More data:
  *Cryptococcus* relatives, *Toxoplasma*, *T. brucei*, other *Candida*; Rabl tissue for grasses (*B. distachyon* root tips, not
  leaves; barley endosperm loses Rabl; Bd21 CRR611509, 289 M pairs, never used). The README gives no tissue or cell-cycle advice.
  **Code:** `fold.py` (`align_profile`), `fuse.py` (`fold_pair`, `align_pair`).

### 15. Centromeres at the edge of a terminal unmappable array

Partly done. The which-end test (`hic_end`, 0.12; design.md) raises a call at the end with the higher anchored score where that holds
on most other chromosomes (goat, sheep, horse, mouse; 54 of the 55 scored raises right). Open: chicken does not pass (12 of 19): its
microchromosomes' arrays are unmappable for 1-1.25 Mb at the very end, so their peaks are `telomeric` and the calls stay `low`. The
one raise not scored right, sheep CP162266.1 (X, `high`, 83.4 Mb from another caller's locus), agrees with the karyotype (CI 0.054,
acrocentric), so the other caller is probably wrong (#37). **Code:** `hicfuse.py` (`end_test`).

### 29. What the contact crossing a centromere is

Waiting for data. Piled up on CENH3, tomato's two flanks contact each other about twice as much as DNA at the same distance elsewhere
(within 0.25 Mb; the figure gives 1.79 against 0.97 at random positions, and no saved file holds the 2.4 once quoted); potato shows a
smaller central peak between deep dips. Both lose contact with the arms out to 6 Mb, a mean effect only (per chromosome the lowest
arm contact is not the centromere: tomato 0/12, potato 0/11, rice 0/12). *S. pombe* has heterochromatin-dependent inter-arm contacts
near its centromeres (Mizuguchi et al. 2014). A fold of the flanks or an insulated domain could explain all of it: a hypothesis.
- **An artifact first.** Pairs spanning a collapsed or missing array are farther apart than their coordinates, a short-range
  depletion that fades with distance; the dips (potato, *Verticillium*) were never compared with the gaps and collapsed arrays.
- **Chance** comes from simulated maps; on moss's re-paired control the dip lies at the LTR cores on 11 of 26 (about 5 expected, P
  0.011, just short; about one in five of 20 control tests gets this low). If the simulation runs short, rebalance each simulated map.
- **Settle:** the profile at several box sizes; Hi-C of heterochromatin or condensin mutants, or other cell-cycle stages, in a genome
  with the signal. **Code:** `hicfuse.py` (`cis_profile`, `cis_chance`).

### 30. A Hi-C test with too few confident references

Built, unchecked. When the confident references are fewer than the gate needs (4), the `low` calls join them, agreements counted at
clean peaks only (`hic_gate_refs` `mixed`, 0.13). It acts only on cucumber 9930 with `tes/` (Hi-C passes on 3 of 7, no call changes).
The older rule for a genome without confident calls (v5) acts only in the siamang (no CENP-A truth), where it raises the two calls that
confidently disagree with the satellome locus (CM054514.2 by 1.0 Mb, CM054529.2 by 15.0 Mb). **Settle:** a genome where it acts, with
truth; the siamang's DiMeLo-seq track (#35). **Code:** `hicfuse.py` (`gate`).

### 32. Calls in a terminal eighth that an interior line contradicts

Open; a candidate waits for a second genome (rule 3). The cis-trans line is read outside the terminal eighths, so it has no say on a
call inside one: tomato without `tes/` keeps 2 wrong `low` calls 8.6% and 11.8% from an end (CM175073.1, CM175075.1; `cis_blind`)
while their interior extremes lie on CENH3. The candidate (0.14) is the chromosome's own trans pair: the test takes each chromosome's
trans peak over its whole length (`hicfuse.anchored_on_cores`, clean peaks only, so an end never counts), and where that pair holds,
the strongest clustering lies at the interior extreme. Both tomato calls are such pairs (calls at 61.06-61.31 and 61.76-62.21 Mb;
extremes 27.00-27.25 in CENH3 26.1-27.5 Mb, 20.25-20.50 in 19.6-20.9 Mb); an acrocentric chromosome among metacentric ones has no pair
and keeps its call. Tomato is the only genome with `cis_blind`. **Settle:** tomato moves its 2 calls; goat, sheep and mouse keep
theirs. **Code:** `hicfuse.py` (`cis_trans_line`), `fuse.py`.

### 33. Cis and trans at one locus as a confident call

Open. Where the cis-trans line passes, a chromosome whose own trans peak lies at its cis extreme holds two readouts at one locus (as
`te_hic_pair` does), yet stays `low`: in G1 *C. neoformans* 10 of 12 pairs lie on CENP-A and 2 lie 2 bins (19-20 kb) off, which as
`medium` would be confident errors. The FASTA-only AT trough is finer (16-40 kb), overlaps CENP-A on 12 of 14 and agrees with 11 of the
12 right cis calls (on CP003829.1 it lies on CEN10, the cis bin 19 kb off), but troughs pair only with prior-free peaks. **Approach:**
pair the cis-trans line with the trough and take its interval; or a sub-bin fit of the trans spot (Centurion's 2D Gaussian: 2.5 kb
mean deviation at 40-kb bins in *S. cerevisiae*). **Settle:** *V. dahliae* and *N. crassa* at 10-kb bins; randomized controls silent.

### 49. Hi-C cannot see inside megabase arrays

Open. At about 100 M pairs and MAPQ 30, the CENP-A domain lies in bins without a uniquely mapped contact on 139 of 157 T2T mammal
chromosomes (blind blocks of median 2.5 Mb in chimpanzee to 19 Mb in gorilla); the anchored peak overlaps CENP-A on 5 of 157, a median
3.2 Mb away at the blind block's edge. A passing test certifies that Hi-C sees the array, not the kinetochore. MAPQ ≥ 1 changes nothing
(near-identical copies map at MAPQ 0). **Settle:** remap one ape library keeping multi-mapped reads, or through array-unique k-mers, and
measure what is recoverable; Pore-C if public data exist for a T2T genome with CENP-A truth. **Code:** `hic/io.py`, `hic/bins.py`,
`hicfuse.py` (`visible_window`).

### 50. Satellite-free centromeres in mammals

Open. Horse chr11 (CM128059.1), the classic evolutionarily new centromere and the benchmark's only satellite-free mammal centromere
with truth, gets no call in any arm (nor EquCab3's chr11 and chr31). The rule for its case exists (`hic_medium_peak`, a `low` call),
but horse Hi-C fails its test at every depth; the library comes from another horse, and chr11's CENP-A slides about 500 kb between
horses. The cis signal works in no mammal. A gap-rank test (does the call's gap have the top flank score among its chromosome's gaps?)
reached P ≤ 0.01 in horse (9 of 11), pig (8/11) and rhesus (13/16), where the gate fails, and 0 of 14 on randomized Hi-C; never
adopted (v3). **Settle:** the same-animal horse set and the multi-horse CENP-A union (#35), then the gap-rank test with a stated role
(a `low` call where there is no sequence call); other mammals and controls unchanged. **Code:** `fuse.py`, `hicfuse.py` (`gate`).

### 51. Hi-C tests at their threshold

Open. The tests flip rather than fade with depth: both orangutans pass at full depth and fail at half (the bin doubles to 500 kb); pig
fails at 100% and 50%, passes on 20 chromosomes at 25% and 10% (500-kb and 1-Mb bins), fails at 5%, its calls unchanged; rice passes
at 50% and fails at 25%. A verdict that improves when reads are discarded suggests the bin decides as much as the signal. The depth
series comes from older code; the share of about 8 random halves that pass was never tried (a split-half flag was rejected: design.md).
**Settle:** the depth series rerun on current code; the gate's chance model checked across bin sizes at fixed depth; the share of
passing halves reported as README depth guidance, not a threshold. **Code:** `hic/bins.py`, `hicfuse.py` (`gate`).

### 52. Hi-C alone, untested

Partly done. In a genome without sequence calls, a clean `high` prior-free peak becomes a `low` call that nothing tests; the strongest
trans cluster may be something else, and a holocentric genome without sequence calls would get such calls (design.md, Known risks).
Since 0.15 a peak held at `medium` only by `null_untestable` counts too (*C. tropicalis* 0 → 6 of 7; on randomized libraries the flag
never occurs without `not_converged`, which keeps a peak unclean; #69).
- **Sharp single-bin peaks**: potato's 7 at 12-27 times the chromosome median are `high` and wrong; one carries `hub`.
- **An effect-size readout** (median peak over median, ≥ 2) separated *C. elegans* (1.21) and randomized libraries (1.05-1.51) from
  yeast (5.96) and *Plasmodium* (5.40); never tried as a guard.
- ***Z. tritici***: 19 of 21, all `low`; 8 of its 21 chromosomes are accessory and may dilute every "most chromosomes" test.
- **Dormant controls** that `cenfuse run` never runs: the read-class split (`--repeats`: trans contacts of reads in repeats against
  unique ones; passed in 6 of 6 mammals as criterion C3 of the first benchmark) and `simulate` (cis-only pairs from the assembly).
- **Settle:** the effect-size guard, a single-bin flag and the read-class split by default (satellite arrays and LTR cores as the
  repeat class): holocentric and randomized libraries give no Hi-C-alone call, the yeasts and *Plasmodium* unchanged, *Z. tritici*
  scored on its core chromosomes. **Code:** `fuse.py` (`hic_denovo_only`), `hic/qc.py`, `hicfuse.py` (`repeats=None`),
  `hic/simulate.py`.

### 53. Where the cis-trans line passes, every extreme is trusted

Tested, not adopted; the diagnostic is open. Where it passes, the line acts on every chromosome, not only where that chromosome's own
trans peak meets its extreme. Some extremes are far off (*N. crassa* 911 and 82 kb on 2 of 7, *V. dahliae* 112 kb on 1 of 8); no call
changed only because those chromosomes had Hi-C or AT calls. In v10 anchoring at the other chromosomes' extremes agreed in rice on 8
of 12 (P 2e-10), none right; the line now skips terminal eighths (where rice's extremes lay) and runs only where the gate fails, so it
never ran on rice, soybean or the apes. Acting only at pairs (0.14) would lose 6 right calls and fix none (G1 *C. neoformans*:
CP003821.1 and CP003826.1 are right non-pairs, its 2 wrong calls are pairs; tomato: CM175070.1 and CM175077.1 right non-pairs).
**Settle:** the line as a diagnostic on gate-passing genomes with truth: rice, soybean and the apes pass nothing; the off extremes
would not act. **Code:** `hicfuse.py` (`cis_trans_line`), `fuse.py`.

### 54. The support tier's evidence

Open. The tier lets Hi-C that finds centromeres beyond chance, but not on most chromosomes, raise the calls its significant peaks agree
with (and since 0.13 choose a competing cluster). Its evidence when adopted (agreeing peaks right 61 of 64 times) counts agreement
with another caller as right in goat (12/12), rhesus (8/8) and giraffe (4/4); goat's 8 raises to `high` stand on that caller. The
pre-registered bar (≥ 90% of raises right) was missed: 14 of 16, 88% (pig's two `high` misses). Its risk, raising a wrong call on a
non-centromeric trans cluster (potato: 1 of 64), is not stated in design.md. **Settle:** its raises scored against functional truth
only (cattle and goat CENP-A first: #35), ≥ 90% right or its role narrowed; the risk stated in design.md. **Code:** `fuse.py`
(`hic_support`).

### 69. Why a null test could not be run

Open (0.15 code review). `fuse._clean_denovo` counts a prior-free peak held at `medium` by `null_untestable` as having passed every test
it could take (#52). The Hi-C engine sets that flag for two causes it does not separate: fewer than two partner chromosomes with bins
to resample (*C. tropicalis*; what the rule means), and a non-finite observed ratio or a chromosome never scored in a replicate (a
degenerate case). **To do:** a distinct flag for the degenerate case, in the engine's stage; it could add one `low` call, only in a
genome that cannot test Hi-C at all (no benchmark run seen). **Code:** `hic/confidence.py`, `hic/qc.py`, `fuse.py`.

## AT troughs

### 28. Where AT troughs do not work

Open. The 0.15 test (design.md, AT troughs) passes in *P. falciparum*, *M. sympodialis*, *N. crassa*, *V. dahliae* (deepest trough
right on 13/14, 8/8, 7/7, 8/8) and mouse, on 3.9% of 3,336 GC-shuffled genomes and on no holocentric genome. Where it fails:
- ***C. neoformans*** (P 0.36): its deepest 10-kb troughs are right on 12 of 14 but tie with other AT-rich loci (0.423-0.452 vs
  0.43-0.46), so nothing in the depths sets a class apart. It holds most of the AT headroom (91 of 161); with Hi-C the cis-trans line
  (G1) and `at_denovo_pair` (anaphase) place it. *S. pombe*: `untested` with 3 chromosomes, its troughs right.
- **Rightly failed**: *Z. tritici* (P 0.69; 0 of 21 troughs on the centromere), *C. tropicalis* (P 0.15; 0 of 7).
- ***Chlamydomonas*** no longer passes (P 0.09): its ends make two deep troughs per chromosome (21 of the 30 deepest in a terminal 5%).
- **Latent false passes**: *B. rapa* TC1 (P 0.002), tomato LA0426 (P 0.04); the line acts there only on a chromosome without a call.
- **Settle:** a second feature beside the trough (a gene desert; without `genes/`, low k-mer complexity) so that *C. neoformans* passes
  from the FASTA and *Z. tritici*, *C. tropicalis*, *Chlamydomonas*, holocentric and shuffled genomes do not. **Code:** `atline.py`;
  comparisons in `v14/wt_at/analysis/`.

### 70. The AT trough choosing between two satellites rests on one genome

Open (0.15 code review); a named heuristic (design.md, AT troughs). `at_chose_alt` lets a tested trough choose a `low` call's competing
cluster, though two satellites differ in composition whatever their centromeres do. Only T2T mouse exercises it (17 chromosomes in 3
arms, all right, all `low`); it never fired on the panel. Of the 136 scored `low` calls with a competing cluster, the 2 outside mouse
where exactly one cluster is right and the trough lies on exactly one (*S. tuberosum*, *S. verrucosum*, whose AT test fails, so the
rule does not act) have the trough on the wrong one; outside mouse the first choice is right on 27 of 46 such chromosomes. It stays
because the call stays `low`, the flag names it and the first choice stays reported (`seq_start`, `seq_end`). **Settle:** a genome with
two satellites at one centromere per chromosome, CENH3 or CENP-A truth and a passing AT test. **Code:** `fuse.py`, `atline.py`.

## Methylation

### 16. The CHG line has changed no call

Open. Only *P. annua* and *Arabidopsis* have CHG data here, and the line has never supported or moved a call on scored data. *P.
annua*'s five samples are other plants, so they do not decide, though the signal is real (the dip is the chromosome's lowest point on
12 of 14). *Arabidopsis*' expected domain (315 kb, from genome size) is far below its CEN180 arrays (2.2-3.5 Mb), so the reference
ring (1.5-4 times the domain) falls inside the arrays. In the way:
- **D, the domain size**, is a maize prior for every genome (22.2 Mb over Mo17's 10 centromeres, scaled by genome size, slope 0.97,
  Plackova et al. 2021); it sets the window, ring and candidate spacing, depends on the chromosomes listed, and is computed twice
  (`pipeline.py`, `methfuse.py`). A D measured from the anchor calls was never tested on a second same-plant case.
- **A per-site floor** (MIN_CALLS_SITE 3), which the CpG line dropped for per-bin pooling (macaque kinetochores kept at 1/16 depth).
- **Fraction-only input**: the reader takes modkit or ENCODE bedMethyl only, so the public *Solanum* bedGraphs cannot be read.
- **Settle:** same-plant CHG with CENH3 truth (#35; *S. verrucosum* CPC54 first), with D from the genome's confident calls (median
  width, else the prior, reported) and CHG pooled per bin. **Code:** `meth/params.py` (D_REF_BP, D_REF_G, D_SLOPE, MIN_CALLS_SITE),
  `meth/track.py`, `meth/io.py`, `methfuse.py`.

### 55. CpG kinetochores outside primates

Open. The CpG line was designed on the primates it is benchmarked on (rhesus, *Arabidopsis* and HG002 were not used to design it),
and its specificity rests on rotated tracks and one plant.
- **No negative**: mouse should fail, but its ONT reads carry no modification tags (a HiFi 5mC track would do). No grass, other
  vertebrate or holocentric genome run. Unused: the HGSVC3 CDRs (#35), giraffe's ONT methylation (21 modBAMs) and rhesus's.
- **A gate without a majority**: with five chromosomes the gate has little power (2 of 4 at chance 0.05: P 0.014); wide calls
  weaken it too (#73). CG-only samples decide without a provenance test.
- **Small arrays** hypomethylated as a whole (bonobo and Bornean orangutan minicentromeres) are never refined (10 bins with data
  needed); a wrong array cannot be told (#38). Since the spread depth and the width (2026-11, cpg-robust), chimpanzee, gorilla
  and siamang pass on R9 Guppy tracks (25 of 25, 25 of 25, 26 of 26).
- **Genomes without truth that 0.15.0 passed**: the duck passed 0.15.0 only because its rotations skipped offsets near the
  call; with every offset, 0.15.0's own test passes none of its chromosomes (5 of 44 calls positional instead of 24), and 21 of
  its 47 calls span at least a sixth of their (micro)chromosome, where no exact test within the chromosome reaches a small P.
  The width alone would pass it (43 chromosomes); the three measures together pass 2: minP pays for the measures that carry
  nothing there. The house wren likewise (0.15.0 with every offset: 0 of 40), which the depths pass.
- **A perturbation**: CHM13 clones whose methylation was changed by dox moved CENP-A (PRJNA1270043, fiber-seq and DiMeLo-seq).
- **Settle:** kinetochores within 50 kb of the HGSVC3 CDRs; a mouse HiFi 5mC track fails the gate; in the perturbed clones the called
  kinetochore follows CENP-A; rotated controls silent. **Code:** `cdr.py`.

### 73. The CpG test in wide calls

Partly done (2026-11): gorilla passes with both technologies since the width (0.15.0's statistic, read on each bin alone)
joined the depths in the test; its kinetochores lie on DiMeLo-seq CENP-A on 9 of 9 truth chromosomes. Open for genomes of
small chromosomes (duck, #55). The text below is from before.

The positional test compares a call's deepest dips with those that methylation from anywhere on the chromosome makes on
the same array bins, so its power falls as calls widen. Gorilla's dips lie on DiMeLo-seq CENP-A on 10 of 10 truth chromosomes,
yet its 3.5-8 Mb calls keep the test below the gate with ONT R9 and with PacBio HiFi (the depth-only test of the cpg-vgp
branch passed the HiFi track narrowly, on 15 of 25). The macaques (calls a median 9-12 Mb) fail too, though there the dips
miss CENP-A on Anvil (pig-tailed macaque), so failing may be right.
- **Tried:** Stouffer's method for the gate (rescues nothing, lowers siamang); testing only the largest array run (arrays fill
  90-99% of calls).
- **Leads:** the active array alone, where the satellite line can name it (HOR structure); a test of the dip's place within the
  array rather than its depth against the chromosome; for genomes of microchromosomes, rotations across chromosomes (which
  changes the question: chromosomes differ in methylation as a whole).
- **Settle:** gorilla passes with both technologies while every rotated control stays silent. **Code:** `cdr.py`.

## Inputs, chromosome set and other features

### 8. Large unplaced scaffolds and dot chromosomes

Partly done. On 121 polyploid assemblies with known chromosome lists, the length rule keeps 31 sequences outside the lists (1.5-117 Mb)
in 12 assemblies, and extra sequences dilute the prevalence tests (on the lists *Cissus* gains 14 `high` calls, *Mimulus* 16;
*Malus toringoides* moves 19 calls). It also drops a dot chromosome alone after a tenfold gap (*D. melanogaster* chr4 without the Y),
dots before a tapering scaffold tail (a simulated 37-chromosome bird karyotype keeps 31), and small chromosomes spanning a tenfold
range; it keeps organelles (yeast's mitochondrion, NC_001224.1, enters every *S. cerevisiae* denominator; headers carry accessions
only). MAX_CHROMOSOMES 50 (a warning) and MAX_SEQUENCES 5,000 (no one-per-chromosome test) express one concept 100-fold apart. Since
0.13 complete molecules left out (telomere arrays at both ends) are named (`left_out_telomeric`), never added; it fired in none of the
0.15 release runs and cannot act without a telomere unit (#56).
- **Settle:** large scaffolds kept without a telomere, centromere-like array or genes flagged; an organelle test from the sequence
  (the organelle rRNA context `rdna.py` reads, or size and copy number); the 121 lists as truth, *D. melanogaster*, bird karyotypes and
  the simulated karyotypes as regression tests, the rule's 383 assemblies unchanged. **Code:** `inputs.py` (`chromosome_set`,
  `left_out`), `meth/track.py` (`chromosome_count`).

### 56. Telomeres and rDNA arrays the benchmark does not count

Partly done: benchmark.md ([Other features](benchmark.md#other-features)) lists every genome's telomere unit, ends with an array, complete chromosomes and 45S arrays.
Open: a unit is learned only from arrays at two ends or more (`telomere.MIN_ENDS`), so potato (one 378-bp TTTAGGG array at chr11's
start), *S. pombe* (at most 8 copies of a GGTTAC-like k-mer per 10 kb) and *C. nigoni* (one 240-bp TTAGGC array, under
`MIN_ARRAY`) get none, and a single-end array is not reported. NOR recall was never measured (goat: arrays on 3 of its 5 NOR
chromosomes); the 45S scan needs 3 motif copies, so a collapsed 2-unit locus is missed (*Cryptococcus* JEC21). Four T2T mammal ends
hold an array the truth rule (terminal 10 kb) does not count, never inspected.
**Settle:** single-end arrays reported; a 2-copy rDNA rule from the organelle context the 45S scan reads; telomere and NOR recall in
the benchmark tables. **Code:** `telomere.py`, `rdna.py` (`MIN_UNITS`).

## Several plants

### 22. Pa1A in other P. annua plants

Waiting for data. Young CRM insertions of four of the five other plants cluster 9-15 Mb from Pa1A's centromere (163.0-169.4 Mb, from
earlier analyses), which suggests other centromeres; J10's lie inside it, and C11 and J10 show a plant-only methylation dip at their
clusters. 0.15 calls Pa1A at 163.0-169.4 Mb (`high`, #31), with plant verdicts 4 `shared` (C11, J10, P3, R1) and 1 `unclear` (I3).
Over all 70 plant-chromosomes: 59 `shared`, 5 `unclear`, 2 `possible_shift` (J10 Pa4B, Pa5A), 2 `call_absent`, 2 `no_data`, none
`shifted`, at 4-8× reads (`v14/runs/rel15c/ALL/poa/annua.plants.tsv`). Nothing to build until data can settle it. **Settle:** CENH3
ChIP or deeper long reads from those plants.

### 57. Per-plant verdicts that depend on other chromosomes

Open. A plant's verdicts change when another chromosome's call moves: when 0.13 moved Pa1A, J10 Pa4B went from `shifted` (loss q
0.015) to `possible_shift` (q 0.32) at the same site, and R1 Pa3B from `shared` to `unclear`. The cause (the multiple-testing pool or
the between-plant test) was not traced, and plants.md does not say it. Untested since 0.1 (plants.md, Untested): the noise model on two
samples of one plant, coverage below about 2.5×, the per-plant Hi-C path (pairs cannot be genotyped, so another plant's library
decides as `unknown`). **Settle:** a leave-one-call-out check of each verdict; one sample's reads split in two must come out as the
reference plant. **Code:** `plants.py`, `provenance.py`.

## Reporting and outputs

### 58. Report wording and undocumented outputs

Partly done (0.14: failed sections recorded in `incomplete`; every column and the Hi-C engine flags documented, with a test). Open:
- **One table of flag → basis phrase, caveat phrase and settling input**, read by `basis` and `caveat`, tested flag by flag (#63).
- **`tes/` advice where it cannot act**: caveats advise `tes/` in mammals and fungi too (`interpret._could`; the siamang's confident
  misses, every *C. tropicalis* no-call); cenfuse does not know the lineage.
- **`hic_sequence_conflict`** reads as blame on the sequence line but fires where it is right (horse T2T 29/32, pig, goat, sheep,
  rhesus, barley 6/7). **`hic_rabl`** is the engine's telomere-pole test (it passes on non-Rabl rice and most T2T mammals), yet
  outputs.md lists it among the Rabl-axis tests.
- **The HOR summary**: the strongest class period equals the top tandem unit on 6 of 46 HG002 calls (`sat/hor.py`, unchanged since
  0.10).
- **No per-run homolog check** (homoeologs, haplotypes), though cheap (cenrep 0.8.9: one goat's assemblies chose different main
  families); the panel check (#68) runs outside cenfuse.
- **Replaced calls**: their own confidence is not stored; `basis` calls a replaced long-unit array a "weak satellite array".
- **Code:** `interpret.py`, `output.py`, `sat/hor.py`.

## Speed, memory and reproducibility

### 5. A faster search

Partly done: since 0.11 long models search with low-complexity tracts neutral, 1.00-4.93× faster on 14 genomes (the table is in docs/benchmark.md at commit 8059fb7). Open:
- **What it does to results.** Over 14 full reruns no chromosome with truth (143) changed, but only 8 of those genomes have truth and
  the rest of the benchmark replays older searches (#34). Overlapping families trade the lead in *Cuscuta* and *Panax* (19 calls change
  main family, 2 change place): a neutral model scores less, so a close family without a tract wins the overlaps (cotton: 906 copies
  move to a relative). **Settle:** every arm rerun from the search; overlaps compared only on positions both models can score.
- **The cost of real hits.** Chain families hit nearly everywhere; dense short-unit arrays hold the rest (alpha satellite: 62% of
  CHM13's CPU; mouse: 14 dispersed models took 95% of a 2-Mb chunk). Tried: [Tried and set aside](#tried-and-set-aside-not-in-designmd).
  **Settle:** a search that treats the seed arrays apart, with monomer and array tables equal to the stored ones (CEN180 on Col-CEN:
  66,131 in Naish 2021, 67,241 by cenrep).
- **The runtime model** of the progress estimate (model length × genome size) fails on heavy genomes (*Panax*: 7.9 h predicted, over
  25 h actual; R² 0.72 over 104 genomes). Runs record their cost (`run_info.json`).
- **Code:** `sat/annotate.py` (`search_cost.tsv`, `resolve_overlaps`), `sat/util.py` (progress).

### 6. Faster seed clustering

Open. BLAST takes 98% of seed clustering, driven by multimeric seeds (*Acanthus* 76 min at 16 threads; median over 106 genomes 0.7
min), and the clustering reads little of its output (*Acanthus*' heaviest round: 7.0% of the HSPs). The thread count changes no family
(0.14: *S. verrucosum* and *Acanthus* at 4 and 16 threads; v13 `bench/exp_threads.py`), though BLAST's hit lists depend on it below
what the clustering reads (`sat/families.py`'s comment says the hits are identical). Tried: aligning only the consulted pairs (10×
faster, but hit lists differ for 0.2-0.5% of them, so not exact; v9 `bench/variant_blast.py`); a fixed partition of queries over
single-thread BLAST processes (0.13: 0.79× on the heaviest round, not identical). **Settle:** a two-pass search (the best HSP
per pair, then full tiling only for consulted pairs) with identical clusters and representatives on every stored seed set. **Code:**
`sat/families.py` (`cluster_seeds`, `blast_pairs`).

### 7. The FasTAN crash, reported upstream

User's decision; not filed yet. A 642.5-kb window of *Oryza minuta* (GWHESQK00000003:29,940,000-30,582,500), about 90% one satellite
array (15-mer spacings of 243 and 342 bp), aborts FasTAN 0.8 every time (`malloc(): corrupted top size`); neither half crashes alone.
*Sonchus oleraceus* crashed at another site. cenfuse recovers by splitting the piece (757 rows after 4 splits, nothing skipped). **To
do:** file with the reproducer (`FAtoGDB`, then `FasTAN -m -p -T1` on the window) and cenfuse's crash table; keep the splitting.

### 60. Reproducibility across environments

Partly done in 0.14 (design.md, Engineering notes): `run_info.json` records the versions of cenfuse, Python, packages and every tool,
the inputs' fingerprints and the cost; `environment.lock.txt` pins the exact environment; the seed cache keys on FasTAN's binary and
layout state; `--satellite-dir` warns on another version; workers log under every start method. Open: on macOS, without `setarch`,
about 1% of FasTAN's array cuts differ (a warning only): measure what that does to calls and state it. CI
(`.github/workflows/tests.yml`): main is now pushed; record its first result. **Code:** `manifest.py`, `sat/seeds.py`, `sat/fastan.py`.

### 61. Memory

Open. The satellite search holds about one workspace per thread (up to ~2.4 GB for 2-kb units) with no cap, so `-t 32` on a 64-GB
node can be killed. Peaks reach 74 GB (rhesus, 12 threads), well above 12 workspaces, unexplained; the 13 tomatoes (0.8-1.2 Gb, 16
threads) peaked at 7-46 GB, so genome size and threads do not set it. The CHG line holds about threads × 1-2 GB; the 3D feature-class
test builds dense N × N matrices and is skipped above 15,000 bins. Tried: [Tried and set aside](#tried-and-set-aside-not-in-designmd).
**Settle:** profile the rhesus peak; then a `--max-memory` resource option limiting concurrent long-model searches, finishing with
identical outputs; the feature-class statistic from the sparse matrix. **Code:** `sat/annotate.py`, `hicfuse.py` (`feature_contacts`).

### 62. Hi-C engine speed

Partly done in 0.14 (one bootstrap, a faster resampler, nulls and cis simulations on the thread pool, one leave-one-out; every
number unchanged; the 139 Hi-C stages 4.6 h instead of 8.5 h). Open: each Hi-C library reruns the whole stage; a 2n = 80 polyploid
was never timed on real data. **Settle:** identical tables, faster on chicken and a 2n = 80 polyploid. **Code:**
`hic/confidence.py`, `hicfuse.py`.

## Code and docs

### 63. The integration cascade and other large functions

Partly done in 0.14: `fuse._decide` is one function per rule in a declared order (`CASCADE`, `PHASES`), each call names its rule
(`rule`), and an equivalence test checks it against 0.14.0's single function. Open:
- `interpret.genome_section` (223 lines), `caveat` (124) and `basis` (123) need coordinated edits for every new flag (#58);
  `pipeline.Run._run` is 189 lines.
- `hicfuse.py` holds four concerns in 862 lines and drives the Hi-C engine through a hand-built argparse Namespace (a library call
  instead).
- Binomial tails written several ways and the Hi-C confidence tests twice (one `stats.py`); the chromosome-set rule lives in the
  methylation engine; the LTR signature is computed twice per run and `tes/` parsed about six times.
- The four CLIs differ in exit codes, logging and the meaning of `-o`; the engines' CLIs are undocumented.
- Bare numeric literals inside functions (the constants test covers module-level ones only: closed #64).
- **Settle:** byte-identical benchmark outputs before and after; `test_fuse_rules.py` as the per-rule spec. **Code:** `fuse.py`,
  `interpret.py`, `pipeline.py`, `hicfuse.py`.

### 65. Tests that pass without running

Partly done in 0.14: `CENFUSE_REQUIRE_TOOLS=1` turns a skip for a missing tool into a failure; unit tests cover the checkpoint layer,
locks, manifest, constants, `cli.py`, `linetests.py`, `seqcalls.py`, `teline.py` and parts of `loci.py`, `methfuse.py` and
`provenance.py`; 1,176 tests pass (0.15). Open: the satellite tests use test-scale values for 13 of the 51 parameters (one slow test at
production values); the parallel family build is exact only while `frame_rotation` writes nothing a later step reads, a premise
guarded by a docstring (make it a test); the report's sentences are tested only end to end (#58); no fixture holds IUPAC codes; CI's
first result (#60). **Code:** `tests/conftest.py`, `tests/sat/sat_common.py`, `tests/fuse/test_fuse_e2e.py`.

### 66. Statements the 0.13 outputs no longer bear out

Partly done (0.14 fixed the statements 0.13's audit found). Open:
- **The audit as a release step**: the docs checked against each release's outputs, docstrings citing benchmark sections, not numbers.
- **Old figures** (as of 0.15.0): three of benchmark.md's figures date from 0.10.1 (`benchmark_cpg`, `benchmark_hic`,
  `benchmark_hic2`) and four from 0.12-0.13; each caption names its release. The workflow (`docs/img/workflow.png`, `.pdf`)
  was redrawn for 0.15.0 (`v14/figures/workflow/workflow.py`).
- **Counts that disagree**: potato's dips on CENH3 on 8 of 11 in the cis survey and 9 of 11 in the cis-trans text;
  `sat/families.py` on BLAST hits at any thread count (#6). (The 459 vs 409 chain-call count and the mislabeled annotation
  table were removed or relabeled in benchmark.md.)

## Leads for biology

Hypotheses cenfuse's outputs could test; none is a finding yet.

- **Young insertions mark where a centromere is now.** CRM copies are younger at Pa1A's Hi-C site than at its satellite (#31), and
  choosing the youngest comparable LTR peak moved 27 cores over 32 genomes with truth, 12 wrong to right and 1 right to wrong (#71).
  Test: centromeres that moved within a species with CENH3 tables (soybean, tomatoes: #35) should carry young copies at the new site.
- **Which lineage CENH3 binds.** Tobacco's CENH3 binds Tekay, not the CRM beside it; a Tekay signature marks all 13 tomatoes (18- to
  53-fold). Tekay may be the CENH3-bound lineage across Solanaceae (#46).
- **The Rabl configuration and condensin II.** The Rabl readouts pass in fungi, *Plasmodium* and barley, in no territorial nucleus;
  condensin-II-deficient lineages keep Rabl (Hoencamp et al. 2021). Predictions: Diptera (no CAP-G2) pass; human cells without CAP-H2
  gain the alignment (public Hi-C); declustering mutants (*S. pombe* csi1) lose it (#14).
- **What clusters in a mammal nucleus.** CHM13's off-target trans peaks sit on pericentromeric HSat2/3 and β-satellite, not the active
  array; clustering differs between mammals at equal depth (pig, horse, sheep, giraffe fail; goat, rhesus reach support).
  `hic_architecture` labels about 50 genomes, without a compartment-matched null (#52).
- **The cis signal's sign**: a peak in tomato, rice and *C. neoformans*, a dip in potato and *Verticillium*, absent in mammals (#29).
- **Telomere repeats at bird centromeres**: chicken's interstitial telomeric arrays lie at its CENP-A windows (P 4e-21; Huang et al.
  2023); other bird T2T assemblies would show whether this is general (#41).
- **5S rDNA beside HSA1q42 genes** in every primate, horse ECA1 and pig SSC14 (v6): a 5S locus perhaps conserved across mammals.
- **Methylation and centromere position**: whether the CpG kinetochore follows CENP-A in CHM13 clones with changed methylation
  (PRJNA1270043) tests the CpG line's premise (#55).

## Tried and set aside (not in design.md)

Measured ideas that design.md's "Choices tested and rejected" omits or records without reason or numbers. Reopen only with new
evidence; most belong in that table.

| idea | result | where |
|---|---|---|
| `malloc_trim` after each search; model-major job order; `MALLOC_ARENA_MAX=2`; freeing hit rows | 229 GB; 125 GB; no change; no change (the fix was one pyhmmer pipeline per thread) | v4-v5 (#61) |
| smaller pyhmmer blocks | 5× less memory and 3× faster on a CEN180 job, but about 40% of raw hits keep their position (multimer models align at any phase) | v5 (#5) |
| tandem-copy models for short units | 4 copies: 43 s per 100 kb (1 copy: 0.1 s) | cenrep (#5) |
| a sample screen to skip dispersed families (why it failed) | 128 families skipped on 216 stored genomes with no false skip, none on real runs: under 7 chromosomes cannot reach P ≤ 0.01; array-poor families cannot be measured (moss's fam009, 41% of its CPU); dense chaining families give too few breaks; removing a family shifts the others' arrays. Code: `v8/patches/v8-screen_screen.py` | v8 (#5) |
| a BLASTN prefilter for long models; a minimap2 prefilter | lost 12% of the hit bp; never tested | v9 (#5) |
| BLAST `-mt_mode`; single-thread query jobs | identical in 12 of 13 reruns, mixed times; 1.6-1.9× per call, 3 of 873,000 hits differ | v9 (#6) |
| a CRM minority as the genome's marker | two right `high` *Arabidopsis* calls fell to `low`; 41 genomes changed | v6 (#46) |
| a share-in-core cutoff for specificity | a cut in a 0.18-0.24 gap fitted to five clades | v6 |
| sites in the core's binomial test too | cores moved in 113 of 220 replay genomes | v8 |
| element counts without a margin | failed its design case (*B. rapa* A09: call 16, core 18) | v8 |
| a core that may favor any cluster it touches | lost *Ae. geniculata* chr02U's CENH3 call | v3 (#9) |
| every candidate clade must disagree before a satellite is rejected | Triticeae chromosomes give 6-8 clades; nothing is ever rejected | cenrep |
| a two-fold concentration margin between clades | labels *Aegilops* genomes inconsistently | cenrep |
| agreement as "the midpoint of either interval inside the other" | 13 of 14 on *Ae. geniculata*, but overrides CENH3-bound satellite calls (CentO-C 4/8 in *O. minuta* C) | cenrep |
| a minimum fold for recruited families | right cores 155, 155, 148, 131 of 172 at none, 5, 10, 20-fold; barley 7 → 3 at 10-fold | v11 |
| an LTR core choosing the competing cluster | no benchmark case to build on (moss's 5 of 5 were solved another way) | v7, v11 |
| the dispersed test without the locus fold | moved 160 batch calls off retrotransposon-derived arrays (*Coffea*, *Glycine*, *Rosa*) | v3 |
| gate chance from the engine's null replicates | removed 1 of 2 spurious passes on randomized sheep; the whole-set majority removed both | v3 |
| full credit for blind bins; replacing only empty partner windows | supported a wrong bonobo call; right anchored peaks 97 → 86 of 139 | v3 |
| both flanks of an unmappable array rising | 75 of 82 centromeric gaps, but also 67 of 96 other gaps under 5 Mb | v3 |
| the fold center by the lower confidence bound | biased to the middle (barley 2/7) | v3 |
| numpy huge pages | the Hi-C bootstrap 30× slower on fragmented memory (now off at import) | v3 |
| peaks restricted to B-compartment bins | gorilla 0 → 0.60, *P. annua* 0.93 → 0.36: the E1 sign at centromeres is opposite in plants and apes | cenhic (#14) |
| arm-scale segregation across the centromere | CHM13 6 of 23, plants at chance | v10 |
| the compartment eigenvector or the trans fraction as a locator | at or near chance; the trans fraction works on a minority, in mammals only | v10 |
| an effect-size test of nucleus type | separated holocentric and randomized libraries from yeast, not Type I nuclei from Type II | v10 (#52) |
| a confidence cap on `hic_sequence_conflict` | the flag fires where the sequence calls are right | v8 (#58) |
| the arms' depletion as a locator | tomato 0 of 12, potato 0/11, rice 0/12 per chromosome (a mean effect only) | v11 |
| the engine's peak levels required for core-anchored Hi-C | cleaned the controls, but left no real hits outside genomes that pass anyway | v11 |
| the Hi-C table swapped for the core-anchored one | potato 6 of 12, short of a majority | v11 |
| CUSUM on the GC deficit; the window by the largest trough gap | *M. sympodialis* 0 of 8; wrong windows in *C. neoformans* and *S. pombe* | v11 (#28) |
| lifting the honesty cap where one-per-chromosome passes | chicken 27 right, 5 wrong; pig 8 right, 9 wrong | v11 (#41) |
| interstitial telomeric arrays as a line | only chicken beats chance (9 calls carry one, 3.9 expected, P 6e-4) | v11 (#41) |
| a genome-level holocentric flag | also on 25 of 26 siamang and 25 of 32 EquCab3 calls | v4 (#44) |
| CpG: a cluster whose valley floor is at most half the call's | 6 changes in 7 primates: 1 right (CM054652.2), 2 wrong | v4 (#38) |
| a CpG null from other satellite DNA; CpG bins as the mean of per-CpG fractions | infeasible (MFA8: 3.7 Mb outside the calls); MMU8 0 of 21 kinetochores at 1/16 depth | v4, v6 |
| a data-derived CHG domain size | rejected in v5, numbers not kept | v5 (#16) |
| 45S arrays by the Rfam SSU HMM; even spacing for 45S arrays | stalled on low-complexity sequence; only 100 of 272 arrays are regular | v3, v10.3 |
| 5S from the satellite families' Rfam tags | misses units over 2 kb (mammals, birds) | v6 |
| sheep and goat CENP-A ChIP as truth | the antibody (CST #2187) recognizes phospho-Ser7 of CENP-A, absent in those species: at most 1.1-1.7× | v2 (#35) |
| loci split into dense runs; sparse clusters dropped | chimpanzee +5,757 loci; barley lost 108 of 156 subtelomeric loci | v8 |

## Closed entries

One line each; the full text is in git history (before 2026-11-10).

- **11. Hi-C anchors where the sequence line is wrong genome-wide.** 0.12: better calls (the LTR signature) make better anchors, and
  where the test still fails Hi-C is tested against the LTR cores (design.md).
- **12. Cis signals of satellite-free centromeres.** 0.12: the cis profile, its sign and level tested against the LTR cores (design.md);
  without cores, #26.
- **17. Knobs and the size floor.** 0.11: `classify.related`'s period test works both ways, so maize's knob family no longer counts
  CentC's chromosomes; reopen with another knob-rich line with CENH3 (size floors: #39).
- **18. Organelle insertions called 45S arrays.** 0.11: a 45S motif hit counts only with nuclear context (design.md).
- **19. Long-unit calls are mostly wrong.** 0.11: an unsupported long-unit array is a feature, not the call (`long_unit_unsupported`).
- **20. The rDNA spacer tag at its threshold.** 0.14: the 11 *B. rapa* accessions run as one complex tag the family in all 11, 110 of
  110 calls right (v13 `runs/joint14`); several accessions of one species are run together.
- **21. Moved calls whose sentences disagree.** 0.11: `caveat` reads the evidence `basis` reports.
- **23. A family's unit length can overstate its copies.** 0.11: `basis` gives the array's own spacing (`seq_spacing`).
- **24. Telomeres outside the rules.** 0.11: 13-29-bp units and loose stretches at sequence ends; genomes with no unit: #56.
- **25. Non-CRM centromeric LTR retrotransposons.** 0.12: the family-level signature, in which another lineage decides at 10-fold
  (design.md); CRM's own standing: #46.
- **26. The cis signal without `tes/`.** 0.13: the library's own trans contacts test the cis extremes (`cis_trans_line`, design.md);
  what remains: #32, #33, #53.
- <a name="36-no-held-out-test-of-the-013-rules"></a>**36. No held-out test of the 0.13 rules.** 0.15, the user's decision
  (2026-11-08): no data are held out; rules made on few genomes are checked under Rules for any change, 3-5. What the LTR branches
  that move a call did on *S. verrucosum*, *Ae. geniculata* and the 13 tomatoes is in #10.
- **59. Checkpoints that miss changed inputs or code.** 0.14: one checkpoint rule (code, parameters, input content), `--force all`, an
  output-folder lock, atomic writes (design.md, Engineering notes); a changed search stage stops for `--force` by design.
- **64. Constants without a reason in design.md.** 0.14: every module-level constant has a reason, enforced by
  tests/fuse/test_constants.py, and docs/constants.md is generated from the code; bare literals inside functions: #63.
