# Future development

Open problems and planned work after cenfuse 0.13.0. The entries come from the benchmark ([benchmark.md](benchmark.md)), a replay of
106 allotetraploid plant genomes (2,355 calls), and the development notes of every version, from the predecessors (cenrep, cenhic) to
0.13, each checked against the 0.13.0 code (2026-11-07). An entry gives the evidence, an approach, the check that would settle it, and
where in the code to start. Entry numbers are stable, because the code, the tests and the other docs cite them (`future.md #31`): a
closed entry keeps its number, and a new one takes the next free number.

## How to read this file

| status | meaning |
|---|---|
| open | nothing built yet |
| partly done | part built, the rest open |
| built, unchecked | the rule is in the code, and no truth has tested it |
| tested, not adopted | tried and measured, kept out for the reason given; reopen only with new evidence |
| waiting for data | nothing to build until new data exist |
| user's decision | an outward-facing step that the user decides |
| closed | done; kept for its number and its evidence ([Closed entries](#closed-entries)) |

**Words.** An *arm* is one benchmark run of a genome with one input set (FASTA + Hi-C, say); a *chromosome-run* is one chromosome in
one arm. The *release runs* are the arms rerun for a release; most *replay* the satellite engine from stored outputs (`--satellite-dir`)
instead of searching again. *The Hi-C test* is the gate (design.md, Hi-C step 6). A line *decides* when it may replace or move a call
(design.md, What may decide). *Confident* means `high` or `medium`. *Headroom* counts the chromosome-runs where some line's own locus is
right but the final call is not.

**Where things live.** future.md, design.md, outputs.md and plants.md are edited in the repo. README.md and docs/benchmark.md are
generated: `v12/bench/fill_docs.py` fills `v12/docs_src/*.template.md` with the tables of `v12/results/md13/`, so an edit made in the
repo is lost at the next fill (#34). The scorer, the truth files, the negative-control generators and the analysis scripts cited below
also live outside the repo, in the development folders `v2/` to `v12/` of `/data2/chris/centromere_benchmark/` (v8 = 0.8, v9 = 0.9,
v10 = 0.10.0, v10.1-v10.3 = 0.10.1, v10.4 = 0.11, v11 = 0.12, v12 = 0.13; v2-v7 are 0.1.0 development; cenrep and cenhic sit at the top).

## Rules for any change

Every release so far has kept these (design.md, Principles; benchmark.md, section 12). An entry is done only when they hold.

1. **Write the check first.** Before a run, write down what the change should do, on which genomes, and what must not change (each
   release's PLAN.md).
2. **Give every constant a reason**: a convention (P ≤ 0.01), a majority, a value measured in the data, or a named heuristic. Add no
   threshold options, and choose no value because it matches known answers.
3. **Recheck a rule prompted by one case.** It needs a reason independent of that case, a recheck on cases not used to make it, and a
   row in benchmark.md's table of rules changed after a result was seen. A rule made on one genome waits for a second case.
4. **Keep the controls silent.** Randomized Hi-C (57 runs) passes no test and changes no call. Shuffled LTR positions (157 annotations)
   give no signature. Holocentric genomes get no confident call. Misplaced methylation passes no test. Shuffled GC passes the AT test
   at most at its level (5%).
5. **Do no harm.** No right call becomes wrong and no confident error is added on the earlier benchmarks, and every changed call is
   listed (benchmark.md, section 13). The goal order is: right, then unsure (`low`), never confidently wrong.
6. **Report a missed criterion as written**; do not reword it.
7. **Test every rule** on a small synthetic fixture. Develop the LTR line on LTRquest annotations only.
8. **On the shared server**, stay within the thread budget the user sets (100 in 0.13's notes), detach long jobs (`setsid nohup`), and
   match your own jobs by full path, never by a pattern.

## Where to start

The 0.13 release is right on 2,174 of the 2,624 chromosome-runs with truth in its real arms (82.9%). Taking, on each, a line whose own
locus is right would reach 91.7% (`v12/bench/headroom.py`, `v12/results/headroom13.tsv`). Most of that headroom lies in AT troughs
(117 of the 231 chromosome-runs, 113 where only the trough is right) and in competing clusters (64, mouse 38). Twelve confident misses
remain (benchmark.md, section 11), and until this rewrite 11 of them were in no entry. In order:

1. **#34 Put the benchmark in the repository**: the scorer, a truth manifest, the control generators and the doc templates. Every other
   check depends on it.
2. **#35 Score the truth already on disk**, then fetch the public sets. On disk: *Ae. geniculata*, an allotetraploid with CENH3 on all
   14 chromosomes, and *S. verrucosum*. Public: wheat, oat, 14 tomato lines, 4 soybean lines, cattle, gibbon, a same-animal horse set and a
   second pig ChIP. They bear on #1, #3, #4, #9, #10, #27, #36, #45 and #50.
3. **#36 Freeze a held-out set for the 0.13 rules** before writing the next rule. No 0.13 rule has one.
4. **#45 Confident calls at the edge of a wide LTR core**: 5 of the 7 confident errors with `tes/` (maize chr5, *Ae. tauschii* Chr5 and
   Chr7, potato chr08, tomato CM175074.1).
5. **#28 and #39 Fungal centromeres the FASTA already marks.** Each chromosome's deepest AT trough lies on the centromere in *V.
   dahliae* (8/8) and *C. neoformans* (12/14). Small tandem arrays below the satellite engine's 20-kb priors mark 7-10 of 14 *C.
   neoformans* centromeres. From the FASTA, cenfuse calls none of them.
6. **#59 Fix checkpoints that miss changed inputs.** After `tes/` or `genes/` is added to an earlier run, the run resumes on stale
   satellite calls unless `classify` is forced.
7. **#60 Write a run manifest and a lock file.** "Fresh runs are identical" holds only within one environment, and no output records
   the tools.
8. **#37, #38 and #50 The mammal misses.** Seven of the 12 confident misses are scored against truth from other animals or from another
   caller (pig 3, siamang 2, sheep 1), or against a disputed array (Bornean orangutan). Better truth exists for pig, siamang and horse
   (#35).
9. **#63 Split `fuse._decide`** (217 lines, 161 branches) into ordered rules with a `rule` column before adding more rules.
10. **#38 Find which array holds the kinetochore.** Mouse's competing cluster is the right one on 19 of 20 chromosomes.

## All entries

| # | area | problem | status | next step |
|---|---|---|---|---|
| [1](#1-calls-on-chains-of-partial-copies) | satellites | calls on chains of partial copies | partly done | CENH3 where a chain dominates (#35) |
| [2](#2-the-40-family-cap) | satellites | the 40-family cap | open | a genome where it costs a call |
| [3](#3-satellites-specific-to-subgenomes-vs-the-honesty-cap) | satellites | satellites specific to subgenomes vs the honesty cap | tested, not adopted | CENH3 for an allopolyploid with such satellites |
| [4](#4-the-decisive-crm-test-with-subgenome-specific-lineages) | LTR | the decisive CRM test with subgenome-specific lineages | waiting for data | CENH3 in the named species |
| [5](#5-a-faster-search) | speed | a faster search | partly done | the whole benchmark searched again; the cost of real hits |
| [6](#6-faster-seed-clustering) | speed | faster seed clustering | open | a two-pass BLAST; BLAST at two thread counts |
| [7](#7-the-fastan-crash-reported-upstream) | engineering | the FasTAN crash, reported upstream | user's decision | file the report |
| [8](#8-large-unplaced-scaffolds-and-dot-chromosomes) | inputs | large unplaced scaffolds and dot chromosomes | partly done | the other edge; organelles; known karyotypes |
| [9](#9-cenh3-truths-for-allopolyploids) | benchmark | CENH3 truths for allopolyploids | partly done | score *Ae. geniculata* (on disk) |
| [10](#10-the-ltr-retrotransposon-line) | LTR | the LTR retrotransposon line | partly done | out-of-sample test (#36); the head-to-head test acts on no call |
| [11](#11-hi-c-anchors-where-the-sequence-line-is-wrong-genome-wide) | Hi-C | Hi-C anchors where the sequence line is wrong genome-wide | closed | |
| [12](#12-cis-signals-of-satellite-free-centromeres) | Hi-C | cis signals of satellite-free centromeres | closed | |
| [13](#13-bird-microchromosomes-in-hi-c) | Hi-C | bird microchromosomes in Hi-C | open | the per-chromosome bin cap; deeper bird Hi-C |
| [14](#14-the-rabl-axis-where-the-fold-has-no-power) | Hi-C | the Rabl axis where the fold has no power | partly done | partners on confident calls; Rabl tissue |
| [15](#15-centromeres-at-the-edge-of-a-terminal-unmappable-array) | Hi-C | centromeres at the edge of a terminal unmappable array | partly done | chicken; the sheep `high` call |
| [16](#16-the-chg-line-has-changed-no-call) | methylation | the CHG line has changed no call | open | same-plant CHG with CENH3 (data exist) |
| [17](#17-knobs-and-the-size-floor) | satellites | knobs and the size floor | partly done | other knob-rich lines with CENH3 |
| [18](#18-organelle-insertions-called-45s-arrays) | rDNA | organelle insertions called 45S arrays | closed | |
| [19](#19-long-unit-calls-are-mostly-wrong) | satellites | long-unit calls are mostly wrong | closed | |
| [20](#20-the-rdna-spacer-tag-at-its-threshold) | satellites | the rDNA spacer tag at its threshold | open | run the 11 *B. rapa* accessions together |
| [21](#21-moved-calls-whose-sentences-disagree) | reporting | moved calls whose sentences disagree | closed | |
| [22](#22-pa1a-in-other-p-annua-plants) | plants | Pa1A in other *P. annua* plants | waiting for data | CENH3 or deeper reads from those plants |
| [23](#23-a-familys-unit-length-can-overstate-its-copies) | satellites | a family's unit length can overstate its copies | closed | |
| [24](#24-telomeres-outside-the-rules) | telomeres | telomeres outside the rules | closed | |
| [25](#25-non-crm-centromeric-ltr-retrotransposons) | LTR | non-CRM centromeric LTR retrotransposons | closed | |
| [26](#26-the-cis-signal-without-tes) | Hi-C | the cis signal without `tes/` | closed | potato and the ends: #32, #53 |
| [27](#27-wrong-ltr-cores-that-look-right) | LTR | wrong LTR cores that look right | open | non-CRM signatures with truth (wild tomatoes) |
| [28](#28-where-at-troughs-do-not-work) | AT | where AT troughs do not work | open | a test that *V. dahliae* and *C. neoformans* pass |
| [29](#29-what-the-contact-crossing-a-centromere-is) | Hi-C | what the contact crossing a centromere is | waiting for data | mutant or cell-cycle Hi-C |
| [30](#30-a-hi-c-test-with-too-few-confident-references) | Hi-C | a Hi-C test with too few confident references | built, unchecked | a genome where it acts, with truth |
| [31](#31-a-centromere-that-left-its-satellite-and-its-ltr-families) | LTR | a centromere that left its satellite and its LTR families | built, unchecked | CENH3 of *P. annua* |
| [32](#32-calls-in-a-terminal-eighth-that-an-interior-line-contradicts) | Hi-C | calls in a terminal eighth that an interior line contradicts | open | a line that can see the ends |
| [33](#33-cis-and-trans-at-one-locus-as-a-confident-call) | Hi-C | cis and trans at one locus as a confident call | open | the AT trough's finer interval |
| [34](#34-the-benchmark-cannot-be-rerun-from-the-repository) | benchmark | the benchmark cannot be rerun from the repository | open | scorer, truth, controls and templates into the repo |
| [35](#35-truth-on-disk-or-public-not-yet-used) | benchmark | truth on disk or public, not yet used | waiting for compute | *Ae. geniculata* and *S. verrucosum* first |
| [36](#36-no-held-out-test-of-the-013-rules) | benchmark | no held-out test of the 0.13 rules | open | freeze a held-out set |
| [37](#37-scoring-rules-and-controls-to-tighten) | benchmark | scoring rules and controls to tighten | open | midpoint scores; full-depth controls |
| [38](#38-which-array-holds-the-kinetochore) | satellites | which array holds the kinetochore | open | settle the orangutan truth; mouse as a test bed |
| [39](#39-satellite-arrays-below-the-engines-absolute-sizes) | satellites | satellite arrays below the engine's absolute sizes | open | a family-free tandem readout, tested one per chromosome |
| [40](#40-a-centromeric-satellite-smaller-than-a-shared-familys-cluster) | satellites | a centromeric satellite smaller than a shared family's cluster | open | chromosome-specific families tested together |
| [41](#41-right-calls-held-at-low-by-the-honesty-cap) | satellites | right calls held at low by the honesty cap | open | a stability test of dominance |
| [42](#42-the-one-per-chromosome-test-near-its-level) | satellites | the one-per-chromosome test near its level | open | potato with genes; oat |
| [43](#43-assemblies-that-lack-or-misassemble-the-centromeric-array) | satellites | assemblies that lack or misassemble the centromeric array | open | the 56 repaired ape centromeres as positives |
| [44](#44-holocentric-genomes-one-low-call-per-chromosome) | satellites | holocentric genomes: one low call per chromosome | open | score the satellite holocentromeres with CENH3 |
| [45](#45-confident-calls-at-the-edge-of-a-wide-ltr-core) | LTR | confident calls at the edge of a wide LTR core | open | agreement by chance, not by gap |
| [46](#46-crm-decides-by-its-name) | LTR | CRM decides by its name | open | Solanaceae with CENH3 |
| [47](#47-centromeric-retroelements-the-ltr-line-cannot-read) | LTR | centromeric retroelements the LTR line cannot read | open | fungal lineage labels; non-LTR families |
| [48](#48-how-much-the-ltr-line-depends-on-the-annotation) | LTR | how much the LTR line depends on the annotation | open | rerun the thinning on 0.13 |
| [49](#49-hi-c-cannot-see-inside-megabase-arrays) | Hi-C | Hi-C cannot see inside megabase arrays | open | multi-mapped reads or array-unique k-mers |
| [50](#50-satellite-free-centromeres-in-mammals) | Hi-C | satellite-free centromeres in mammals | open | same-animal horse data |
| [51](#51-hi-c-tests-at-their-threshold) | Hi-C | Hi-C tests at their threshold | open | the depth series on 0.13; the bin |
| [52](#52-hi-c-alone-untested) | Hi-C | Hi-C alone, untested | open | guards and the dormant mis-mapping controls |
| [53](#53-where-the-cis-trans-line-passes-every-extreme-is-trusted) | Hi-C | where the cis-trans line passes, every extreme is trusted | open | act only where a chromosome's own pair holds |
| [54](#54-the-support-tiers-evidence) | Hi-C | the support tier's evidence | open | score its raises against functional truth |
| [55](#55-cpg-kinetochores-outside-primates) | methylation | CpG kinetochores outside primates | open | HGSVC3 as a held-out set; a mouse negative |
| [56](#56-telomeres-and-rdna-arrays-the-benchmark-does-not-count) | features | telomeres and rDNA arrays the benchmark does not count | open | recall per genome |
| [57](#57-per-plant-verdicts-that-depend-on-other-chromosomes) | plants | per-plant verdicts that depend on other chromosomes | open | a leave-one-call-out check |
| [58](#58-report-wording-and-undocumented-outputs) | reporting | report wording and undocumented outputs | open | one flag-to-sentence table |
| [59](#59-checkpoints-that-miss-changed-inputs-or-code) | engineering | checkpoints that miss changed inputs or code | open | one checkpoint system keyed on code and inputs |
| [60](#60-reproducibility-across-environments) | engineering | reproducibility across environments | open | a run manifest; a lock file |
| [61](#61-memory) | engineering | memory | open | a memory budget |
| [62](#62-hi-c-engine-speed) | engineering | Hi-C engine speed | open | one leave-one-out and null implementation |
| [63](#63-the-integration-cascade-and-other-large-functions) | code | the integration cascade and other large functions | open | ordered rules with a `rule` column |
| [64](#64-constants-without-a-reason-in-designmd) | code | constants without a reason in design.md | open | generate the Constants table from the code |
| [65](#65-tests-that-pass-without-running) | code | tests that pass without running | open | fail on missing tools in CI |
| [66](#66-statements-the-013-outputs-no-longer-bear-out) | docs | statements the 0.13 outputs no longer bear out | open | a docs audit per release |

## Benchmark, truth and validation

### 9. CENH3 truths for allopolyploids

**Evidence.** No allotetraploid in the scored benchmark has CENH3 truth, but one exists on disk: *Ae. geniculata* (UUMM, 14
chromosomes), with CENH3 ChIP-seq peaks on the same assembly (Yunusbaev 2026, Table S2, at Mb resolution), in
`polyploid_genomes/_cenrep_batch/review/Aegilops_geniculata.cenh3.tsv`. 0.10.1 scored it 10 of 14 in its satellite-batch section,
together with *S. verrucosum* (3 of 12 overlapping, 7 of 11 within 1 Mb). The docs rewrite of 0.10.1 dropped that section, 0.12 then
wrote that no allotetraploid has truth, and no release since has scored either genome. *Ae. geniculata* shaped two rules: a core may
favor only the largest cluster (letting it favor "any cluster the core touches" lost chr02U's CENH3 call), and the beside-core test
turned chr02M right. It is therefore development data, not held out. Where a CRM core replaced a satellite array (143 of 2,212 audited
allotetraploid calls, 0.8's run), homoeolog agreement favors the core (26 vs 14 pairs, P = 0.08) and gene density favors the array
(80 vs 34, a biased proxy). Public CENH3 ChIP exists for other allopolyploids (#35): hexaploid wheat with a T2T assembly, *Avena
sterilis* and *A. maroccana*.
**Approach.** Score 0.13 (FASTA only, and with `tes/`) on *Ae. geniculata* and *S. verrucosum*. Then take the public sets of #35,
chosen to settle entries 1, 3, 4, 5 and 10 (`seq_overruled_*` lists the replaced arrays). Truth carried from a progenitor into a
subgenome by unique flanks would give provisional truth (v6 `bench/flank_truth.py` reached 13 of 14 per subgenome in a *Cryptococcus*
hybrid), but centromeres can move after polyploidization.
**Check.** Score as in the benchmark (the Mb-resolution table needs the 1-Mb tolerance); no rerun of the release needed.

### 34. The benchmark cannot be rerun from the repository

**Evidence.** Every number in benchmark.md and the README comes from scripts and data outside the repo, which holds only
`cenfuse/`, `docs/` and `tests/`:
- **The scorer.** The frozen scorer is `v11/bench/rel/score_rel.py` (it still imports code from `v10.2/code/rel`). Its truth table reads
  BEDs from eight version folders and holds the tolerance rules the README summarizes; `v12/bench/score_dev.py` patches it at run time.
- **The negative controls.** Re-paired or relocated Hi-C comes from `v2/scripts/permute_pairs.py --mode repair|relocate --seed 1`,
  TNEG arms from `v11/bench/make_tneg.py`, GC-block shuffles from `at_shuffle.py` and LTR shuffles from `control_te.py`. The generator
  of the misplaced-methylation control was not found. Beside them, `v11/bench/repair_pairs.sh` is marked DO NOT USE: its `shuf` made
  clumped contacts, caught only when a figure measured their dispersion (55-276, against 1.0).
- **The release chain.** `v12/bench/finish_rel13f.sh` runs the scoring, the comparisons, `md_013.py` and `headroom.py`. Then
  `fill_docs.py` writes README.md and docs/benchmark.md from `v12/docs_src/`.
- **Scripts the docs cite.** `bench/pa1a_age.py`, `bench/young_core.py` (entry 31) and `bench/poly_one.py` (benchmark.md) exist only
  in `v11/bench`.
- **Replayed runs.** The runs replay the satellite engine (`--satellite-dir`) and seed their checkpoints from the previous release
  (`seed_run.py`). Apart from 14 genomes in 0.11, the benchmark's satellite searches predate 0.11 (entry 5).

All the scripts hard-code `/data2/chris/...` paths. The `cenfuse` mamba env imports an editable 0.8.0 install from the v8 folder, so
running its `cenfuse` command outside the 0.13 checkout uses 0.8 code silently.
**Approach.** A `benchmark/` folder, or a companion repo, holding:
- a truth manifest (genome, assembly accession, truth source and accession, how the BED was made, tolerance rule) with the BEDs or fetch scripts;
- the scorer as a tested module;
- the control generators as `cenfuse.controls`, with seeds and a dispersion check;
- the doc templates and table scripts, with paths taken from one config file;
- the cited analysis scripts.

Release scripts should stop when the imported cenfuse is not the code under test.
**Check.** From a clean checkout, regenerate benchmark.md's tables from the stored runs and get the same text. One fresh, unseeded run
of every arm should reproduce the published counts, or list every difference.

### 35. Truth on disk or public, not yet used

**Evidence.** The development notes found these sets and never scored them. Most bear on confident misses or on rules that have no
held-out truth:

| data | what it would test | where it is recorded |
|---|---|---|
| *Ae. geniculata* CENH3 (14 chromosomes, same assembly) | an allotetraploid: #1, #3, #4, #9, #10 | `polyploid_genomes/_cenrep_batch/review/` |
| *S. verrucosum* CENH3 (ChIP SRR18548893; a `tes/` file exists) | the Solanaceae satellite line: #40, #45, #46 | `polyploid_genomes/_cenrep_batch/validation/solanum_cenh3/` |
| wheat Chinese Spring (T2T GCA_040256815.2; ChIP SRR9089557/8, DRR821424/5) | a hexaploid; CRM plateaus (#45) | `v8/data/published/manifest_ena_runs.tsv` |
| *Avena sterilis* (GWHERDU00000000; CRR1339975/6); *A. maroccana*, *A. atlantica*, *A. eriantha*, 'Starter' (CRA007177) | oat's dominant non-centromeric satellite (#42); allopolyploids | `v8/data/published/manifest_gsa_runs.tsv` |
| 14 tomato and wild-tomato lines (Shi 2026: 240 centromeres in 20 accessions; SCAU CRA035506) | a held-out test of the tomato rules (#26, #27, #32, #33) | `v8/data/published/download_manifest_cenh3.tsv` |
| 4 more soybean lines (ChiHei23, QiHuang34, ZH13, Williams82) and centromeres of about 27 accessions (PNAS 2023) | CRM plateaus (#45); centromeres that move within a species (#31) | same manifest |
| cattle CENP-A CUT&RUN and input (SRR28462359/8) | an all-acrocentric ruminant: the acrocentric and which-end rules | v2 `data/research/research_truth.md` |
| hoolock gibbon CENP-A CUT&RUN | the siamang's failure mode | same |
| siamang same-cell-line CENP-A DiMeLo-seq (bigWig, 15.1 GB) | the siamang's 2 confident disagreements; the gate on `low` references (#30) | same |
| horse: a same-animal assembly (GCA_052323815.1), Hi-C (SRR36069624) and CUT&Tag/ChIP with IgG (SRR36096936-9, PRJNA1248178); CENP-A ChIP of about 15 horses (96 runs) | horse chr11 (#50); the horse Hi-C that fails its test; centromere sliding | same |
| pig CENP-A ChIP and input from the assembly's breed, Duroc liver (DRR1139646/7; PRJCA032804) | pig's 3 confident misses (#37) | same |
| the repaired ape assemblies (56 centromeres misassembled in v2.1, with repair BEDs), both haplotypes | about twice the ape truth; misassembly as a positive control (#43) | same |
| HGSVC3 CDRs (67 genomes, 2,156 centromeres; Gao 2026) | a held-out test of the CpG line (#55) | figshare 10.6084/m9.figshare.29921990 |
| holocentric CENH3: *R. breviuscula* ChIP (SRR27965399/8), *Chionographis japonica* (GSM7120336), *Luzula sylvatica* (358 domains) | what the satellite line reports on satellite holocentromeres (#44) | v4 `data/research/hard_genomes.md` |
| same-plant CHG with CENH3: *S. verrucosum* CPC54, potato DM FAST5 (about 625 GB), rice T2T-NIP modBAMs | the CHG line (#16) | same |
| Hi-C: *Chlorella* NS4-2 (CRR559578), the primary Bd21 library (CRR611509, 289 M pairs), and deeper libraries for the shallow controls (#37) | Hi-C on LTR cores in an alga; the Rabl axis in a grass | v10 `data/hic/hic_runs_plants.tsv` |

**Approach.** Score the two genomes already on disk first (one rerun each). Fetch the rest in table order, confirming each antibody and
target first: the horse CUT&Tag target and the pig antibody are unconfirmed. Sheep and goat ChIP gives no truth: its antibody (CST
#2187) recognizes phospho-Ser7 of CENP-A, which those species lack, and enrichment never passes 1.1-1.7×.
**Check.** Score each set with the frozen ChIP rule (#37), and add its genomes to the held-out set (#36) before any rule is changed on them.

### 36. No held-out test of the 0.13 rules

**Evidence.** The genomes held out from 0.12 were in view while 0.13 was built ("none is held out from it", benchmark.md), so no 0.13
rule has a held-out test, and each acts on one or two genomes:

| 0.13 rule | where it acts |
|---|---|
| the cis signal tested by trans contacts | tomato and G1 *C. neoformans*; it passes without changing a call in *Brachypodium*, *N. crassa* and *V. dahliae* |
| the support tier's competing cluster | moss only |
| the trough paired with a prior-free peak | *C. neoformans* in anaphase only |
| `mixed` references | cucumber only, with no call changed |
| the age tie-breaker | Pa1A only, without truth |
| `left_out_telomeric` | nowhere (0 of 260 runs) |

Three of these rules were chosen from 0.12's headroom table, yet they are missing from benchmark.md's table of rules changed after a
result was seen. The LTR line is in-sample too. In 0.9 all 20 gains from `tes/` were on genomes the rules were built on (4 genomes up, 1
down, 12 equal; P = 0.375 by genome). The 0.12 held-out set moved no call (*B. rapa* ×7, grape, cucumber, fungi). So every branch that
moves or overrules a call has been checked only on development genomes: `seq_overruled_*`, `te_discordant`, `weak_core`,
`chain_to_core`, and a non-CRM signature deciding. The 10-fold level at which a non-CRM lineage decides comes from 42 polyploids
without truth. The clean-peak guard on Hi-C anchored on LTR cores (0.12) was rechecked on no untouched genome.
**Approach.** Freeze a held-out set from #35 before the next rule is written:
- genomes with satellite-free or beside-satellite centromeres and Hi-C, run without `tes/`, for the cis line tested by trans;
- CRM-plateau genomes (wheat, the soybean lines: #45);
- non-CRM signatures with truth (the wild tomatoes: #27);
- *D. melanogaster* and bird karyotypes (#8).

Write each rule's expected effect on them first, as the release PLAN.md does. Add the three headroom-chosen rules to the
rules-changed table now.
**Check.** The held-out scores, reported as written, with the release.

### 37. Scoring rules and controls to tighten

**Evidence.**
- *Overlap is lenient for wide calls.* Mammal calls are a median 3.6 Mb wide against CENP-A domains of 0.1-0.5 Mb, so a few-kb overlap
  counts as right (in 0.1.0 a 19.9-Mb bonobo call was right by 32.6 kb). No width-normalized or midpoint score is reported.
- *The tolerance without Hi-C.* Fungi and *Plasmodium* are scored within one Hi-C bin, and the frozen scorer uses 250 kb for an arm
  without Hi-C (`score_rel.py`, `bin1`), which benchmark.md's scoring section does not say. *S. cerevisiae*'s AT troughs are "right" on 12
  of 16 that way and overlap the centromere on 2. Barley is scored at 8 Mb (4 bins of its 2-Mb Hi-C) also in its arms without Hi-C.
- *The chance model* draws calls uniformly along the chromosome. That is optimistic where centromeres sit near a known end (acrocentric
  goat, sheep, horse, mouse). 0.1.0 said so; the docs no longer do.
- *Weak or indirect truth.*
  - Pig and horse truth is ChIP from other animals. Pig's three confident misses lie where that ChIP's second peak (3.0-4.4×) falls
    inside the call, and the karyotype agrees with the call: chr11 is bi-armed but its domain sits at CI 0.002; chr17 is acrocentric
    but its domain sits at CI 0.24.
  - Siamang, sheep, goat, rhesus and giraffe are scored against another caller (the satellome locus). The sheep X call it disputes
    agrees with the karyotype check (CI 0.054, acrocentric).
  - Bornean orangutan CM054652.2 may be haplotype-specific (#38).
  - *N. crassa* LG VI carries an inversion (Galazka 2016). *C. tropicalis* CENP-A was carried over from an older assembly.
    *V. dahliae*'s Hi-C was used by its authors to infer its centromeres, so it is no independent Hi-C test.
  - *C. neoformans* counts twice (two Hi-C states) among the 47 genomes of the calibration.
- *The truth rule* (100-kb bins, ChIP/input, half-maximum run, ≥ 4×) gives no reason for 100 kb or 4×, and it decides most plant
  truths: in cucumber only 1 of 7 chromosomes passes it.
- *Shallow controls.* The holocentric and held-out Hi-C libraries are prefixes of the data: *R. breviuscula* 3.3 M pairs, *B. mori*
  5.9 M, *Ae. tauschii* 0.84 M (another accession), sorghum 3.9 M, mouse 15.7 M. Deeper runs exist: *B. mori* SRR25988387 (369 M),
  *R. breviuscula* SRR17429201 (863 M), soybean SRR37686999 (415 M), *B. rapa* CC1 SRR8633037 (368 M), Col-0 SRR2626163 (143 M).
  "No holocentric genome passes" rests on these libraries, and the tests lose power with depth (#51).
- *Attribution.* benchmark.md credits the ape CDRs to Logsdon et al. 2024. The v2 notes record them as the Logsdon lab's curated set
  (v1.3), unpublished ("contact the lab before publication"), and the ape DiMeLo-seq and CUT&RUN tracks as unpublished too.
- *Gone from the report.* The head-to-head with the independent satellite caller (0.1.0's criterion B7; it would tie at 90/92 today)
  and 0.1.0's unique-read Hi-C control (criterion C3, #52) are no longer reported.

**Approach.**
- Report a midpoint distance and the call width beside overlap.
- Score small genomes without Hi-C at a fixed fine tolerance.
- Use a chance model conditioned on end distance for acrocentric chromosomes.
- Run the truth rule at 3× and 5× and with 50-kb bins, as a sensitivity check.
- Remap the shallow controls at full depth.
- Settle the CDR attribution with the lab.

**Check.** Report each scorer change beside the old numbers for one release; no README claim may depend on a lenient rule.

## Satellites

### 1. Calls on chains of partial copies

**Evidence.** A family can concentrate in one place without being a tandem satellite: under half of its bp lie in tandem and under
half in whole-unit hits, so it is a chain of pieces of a longer repeat, usually a retrotransposon. *Glycine dolichocarpa*'s 1,790-bp
family has 44 whole hits of 25,732, 65% inside SIRE elements, and its calls are 10 Mb wide. In the 106 allotetraploids, 41 such
families in 29 species carry 409 of 2,355 calls (benchmark.md's rules table says 459; one of the two is wrong), 209 of them `high`
before 0.11. They include most of the 131 calls of 0.10.1's entry 1, which stood on single-sequence models (214 such models in 68 of
105 species). Those models are nearly all at least 1 kb long, and a better profile would not make them tandem; whether short (< 1 kb)
single-sequence families lose diverged copies was never measured (`min_hmm_seqs` 20). A chain call (under 20 kb of real satellite
inside) is at most `medium`. It is `low` unless a CRM core overlaps it, or a core of another lineage that the family is not made of
(design.md, chain of partial hits). No call moved. The new search changes these counts little: over its 14 reruns, 75 chain calls
become 74 (entry 5).
That rule rests on a principle (only a decisive LTR signature may decide among retrotransposons), not on truth. No benchmark genome has
a chain as its dominant family. In 0.12 the benchmark's 35 chain calls (chicken, potato, moss, tomato) are 29 right and 4 wrong of the
33 with truth, all `low` or `medium` from the sequence line; 3 in moss rise to `high` with Hi-C and `tes/`, all right. Two proxies
cannot settle it. In gene density, calls on CRM chains lie between TE cores and real satellites (0.135, 0.054, 0.009). They match their
homoeolog's place within 5% of the chromosome about as often as the others (0.60, 0.70, 0.66). Two held-out complexes of cenrep failed
the same way, both with "no centromeric satellite" by CENH3 or FISH: *Glycine* (a 1.8-kb tandem family) and *Musa* (a 62-bp repeat from
inside SIRE). A tandem repeat cut from inside an element is no chain, so it is only flagged (`te_derived`).
**Done in 0.12.** Where a decisive core lies inside a chain call and is narrower, it becomes the call (`chain_to_core`; 7 of 7 right on
the genomes with truth, before and after). Such cores are tighter and more gene-poor than the calls they overlap (gene density 0.021 vs
0.053), with about the same homoeolog agreement. A core that straddles the chain call's edge is left as it is (the geometry of #45).
**Check still open.** CENH3 for a genome where a chain dominates (6 cotton species have one on 23-26 of 26 chromosomes); *Ae.
geniculata* and the public allopolyploid sets (#35).
**Code.** `sat/classify.py` (roles, LOCUS_FOLD), `sat/centromere.py` (`partial_chain`, `chain_to_core`).

### 2. The 40-family cap

**Evidence.** At most 40 clusters per complex are modeled, since search time grows with models × genome size; since 0.11 the dropped
ones are listed (`satellite/dropped_clusters.tsv`). Over 71 runs that dropped any: 2,376 clusters, none on half of the chromosomes, 3
on three or more, 1,553 on one; 760 hold under 20 kb per chromosome. Only 4 reach 20 kb on a chromosome without a call, all in *B.
carinata*; 349 lie on a chromosome with a `low` call. So far nothing suggests the cap costs a centromere (0.9's case, *Sonchus
oleraceus*, dropped 60 clusters of 20 kb or more, the largest 107 kb). Changing the family set is not neutral. Overlapping hits go to the
higher-scoring family (`sat/annotate.py`, `resolve_overlaps`), so adding or removing a family moves copies to its neighbors: in 0.8,
moss's fam003 gained 90 kb of arrays and became centromeric, and Chr08's call moved.
**Approach.** A cost budget instead of a count (model length × genome size, priced from `search_cost.tsv`), if a genome ever shows
the cap costing a call.
**Check.** Rerun *B. carinata* without the cap, and compare roles, main families, every family's arrays and the calls.
**Code.** `sat/families.py` (`max_families`, `dropped_table`).

### 3. Satellites specific to subgenomes vs the honesty cap

**Tested in 0.11, not implemented.** When no family is the largest centromeric cluster on most chromosomes, every call is capped
`low`, and in an allopolyploid each subgenome can have its own satellite. Judging dominance within the subgenomes of 20 capped
polyploids lifts 27 of their 292 capped calls, but random splits of the same chromosomes lift 9.3; one species of 20 beats its own
random splits (*O. schlechteri*, P 0.005). Random splits of diploids lift up to 7.9 calls (siamang), 5.7 (pig) and 4.7 (chicken), so
a lift this size is what splitting does anyway.
**Check, to reopen it.** CENH3 truth for an allopolyploid with subgenome-specific satellites: *Fragaria × ananassa* (CENH3 ChIP, Jin
2025, with chromosome-specific satellite types) or *Ae. geniculata* (#35).

### 17. Knobs and the size floor

**Done in 0.11 for maize; one genome.** Nothing recognizes a knob as such. Until 0.11, maize Mo17's 180-bp knob family missed the
centromeric role only because its median largest block (74 kb) is under the 100-kb floor (`strong_cluster_bp`). At a 70-kb floor,
FASTA only, it took chr4, chr6 and chr8. `classify.related` tested relatedness one way only, by unit length (155 vs 180 bp passed,
180 vs 155 failed), so CentC's chromosomes counted for the knob. The test now works both ways: at floors of 100, 70 and 50 kb the knob
stays `interstitial` and CentC keeps all ten, `high` (in a replay; the 0.10.1 code at 70 kb gives the three knob calls). The same change put *B. rapa* CX1
right ([20](#20-the-rdna-spacer-tag-at-its-threshold)). Relatedness by sequence alone, proposed here before, stays out: it changes the
main family on 90 of 3,459 chromosomes and removes dominance in 6 genomes with no truth behind it, and it works against satellites
specific to subgenomes. The floor is one of the engine's absolute priors (#39). Maize's remaining confident miss, chr5, is a different
problem (#45).
**Check.** Other knob-rich lines with CENH3 truth.

### 20. The rDNA spacer tag at its threshold

**Evidence.** A family is tagged as an rDNA spacer when 45S genes flank at least half of its 10 largest arrays. *B. rapa*'s 238-bp
family reaches 5-6 of 10 in 5 accessions (tagged) and 3-4 in 6 (not tagged). In CX1, at 4 of 10, it stays eligible. Since 0.11 the
two-way period test ([17](#17-knobs-and-the-size-floor)) takes it out of CX1's centromeric class, which puts A05 right and narrows A06,
but the tag still differs between accessions.
**Approach.** The code already propagates the tag when the accessions are analyzed together (several FASTAs in one folder): a family
tagged in one genome of a complex is tagged in all (`sat/classify.py`, `apply_rdna(elsewhere)`; the case it was built on is *Lemna*).
Nothing new is needed but the joint run.
**Check.** The 11 accessions run as one complex, scored against CENH3.

### 38. Which array holds the kinetochore

**Evidence.** Where a chromosome holds several candidate arrays, the sequence line takes the largest cluster, and no tested line can
choose another:
- **Bornean orangutan CM054652.2** (`high`, 14.3 Mb off; the only confident ape miss). CENP-A sits on 4 small alpha-satellite arrays
  (73 kb), the call on a 3.4-Mb cluster. Hi-C's peak lies 0.26 Mb from the CDR but is `low` (it passes 1 of the engine's 3 tests), so it
  cannot oppose a `high` call. The deepest CpG bin of the chromosome lies at the CDR, but the called array has a dip too, and the CpG
  line puts its kinetochore there. The truth may be haplotype-specific: CenSat v2.0 marks the called array active, the CDR region was
  redrawn in 2025, and hap2's DiMeLo-seq CENP-A shows signal at both sites.
- **Bonobo CM055489.2** (`low`, wrong). The true array is the reported competing cluster. Only a strong peak may choose it. The two
  weaker peaks that rule kept out were 1 right (this one) and 1 wrong (CHM13 chr21), so the rule rests on two cases.
- **Mouse** (all 20 chromosomes acrocentric). Every call sits on the major satellite at the right end, within 45 kb of the CENP-A minor
  satellite, and the minor satellite is the reported alternative on 19 of 20 (1/20 overlap; 38 of the 231 chromosome-runs of
  headroom). No mappable bin separates the two, so Hi-C cannot choose (`hic_alt_ambiguous`). The rule that makes the minor satellite
  the alternative (its terminal block outside the call) was shaped on mouse alone, after a first form put the alternative there on 17 of
  20; the replays hold no second case of an unrelated pericentric satellite. Three mouse chromosomes (chr2, 5, 14) share one
  byte-identical centromeric-end sequence (9,259-807,536).

Tried for CpG choosing between arrays (design.md, in part):
- the least-methylated window lies on the CDR on 3 of 23;
- the most prominent valley across clusters picks the wrong cluster on 8-10 of 25;
- a cluster whose valley floor is at most half the call's made 6 changes in 7 primates: 1 right (CM054652.2), 2 wrong.

**Approach.**
- Settle the orangutan truth first, from the 2025 CDR, CenSat v2.0 and both haplotypes.
- Give a weak Hi-C peak that opposes a `high` call a second line, such as a CpG dip at the peak, as `moved_with_age` is for a strong one.
  Test it on the ape loci where two references disagree; there cenhic's peaks sided with the published active array 4 of 4 times.
- For two satellites at one end, look for a sequence feature of the CENP-A satellite (in mouse the minor satellite carries CENP-B boxes
  and the major does not; an idea, not yet tried), with mouse as the test bed.

**Check.** The orangutan, bonobo and mouse calls; no right ape call lost; CHM13 chr21 stays right.
**Code.** `sat/centromere.py` (`satellite_call`, `rival_bp`), `fuse.py` (`hic_chose_alt`, `hic_alt_ambiguous`), `cdr.py`.

### 39. Satellite arrays below the engine's absolute sizes

**Evidence.** One absolute prior, S_MIN = 20 kb, scales the satellite engine's size rules whatever the genome. `sat/params.py` marks
these fields "prior": `min_fam_span` and `chrom_min_bp` 20 kb, `low_sat_bp` 50 kb, `strong_cluster_bp` 100 kb, `fallback_min_bp` 5 kb,
`cluster_gap` 500 kb and `diffuse_min_len` 5 Mb. Units under 30 bp or over 30 kb are not searched at all (`min_period`, "misses 5-12-bp
satellites, e.g. Drosophila"; `LU_MAX_UNIT`). So the fungi get no centromeric family. From the FASTA, cenfuse calls 0 of 14 *C.
neoformans*, 0 of 8 *V. dahliae* and 0 of 21 *Z. tritici* centromeres, while tandem-repeat callers at their defaults find many:
- TRASH2 places 7 of 14 *C. neoformans* with 0 wrong (a 225-bp unit forms one ~2.2-kb array in 7 centromeres and nowhere else);
- satellome's largest tandem region per chromosome places 10 of 14 *C. neoformans* and 7 of 8 *V. dahliae*.

These figures come from 0.10.1's comparison, where each tool's readout was picked by the truth, so they are optimistic.
**Approach.** A family-free readout tested like the AT troughs: take each chromosome's densest tandem-repeat region, and accept it
where the N densest fall one per chromosome beyond chance. Scale the size floors to the genome or chromosome, from a reason in the data
rather than a fit.
**Check.** *C. neoformans* and *V. dahliae* placed from the FASTA; no call changed in genomes with large satellites; shuffled and
holocentric genomes do not pass. Rerun the comparison with TRASH2 and satellome on 0.13, since cenfuse's side has changed (the AT line).
**Code.** `sat/params.py` (S_MIN and the 'prior' fields), `atline.py` (the test to reuse), `sat/karyotype.py`.

### 40. A centromeric satellite smaller than a shared family's cluster

**Evidence.** A family is centromeric only where it forms a chromosome's largest centromeric cluster, so a chromosome-specific
satellite smaller than a shared family's cluster on the same chromosome is never called. Potato's centromeres sit on
chromosome-specific kb-unit families while a genome-wide 157-bp family dominates. The 981-bp St24 holds 299 kb at Cen1, and the 1,044-
and 948-bp families do the same on their chromosomes. All of them stay `interstitial` ("never the main or co-located cluster on a
chromosome", prevalence 0.17), so potato is 2 of 11 from the FASTA (10 of 11 with `tes/`). Other *Solanum* centromeres sit on
chromosome-specific satellites with 1-2-kb units too, as in *S. verrucosum* (3, later 4, of 12 in the batch).
**Approach.** Test chromosome-specific families together: take the largest array of each family present on one or two chromosomes,
and accept them where those arrays fall one per chromosome beyond chance. That is the one-per-chromosome test over a set of families,
not one family.
**Check.** Potato from the FASTA against CENH3; *S. verrucosum* (#35); no change where the centromeric family is shared.
**Code.** `sat/classify.py` (roles), `sat/karyotype.py`.

### 41. Right calls held at low by the honesty cap

**Evidence.** Where no family dominates, every satellite call is `low`, even when the calls are mostly right: chicken 27 of 33 (its
centromeres use chromosome-specific satellites, and chr5, chr27 and Z have none), moss 16 of 26 from the FASTA, *Chlamydomonas* 6 of
17. Over the benchmark, `low` calls are right 97 of 141 times. Lifting the cap where the one-per-chromosome test passes was measured

and not adopted: chicken would gain 27 right and 5 wrong, pig 8 and 9 (about 71% right, near `low` and far below `medium`). Dominance
itself is unstable. In 0.8's fresh run, pig's 73/74-bp satellite split differently and its 336-bp satellite lost genome-wide
dominance, so FASTA-only pig fell from 13 of 17 to 8 of 17, all `low` (5 right `medium` calls lowered). The families stage is now
deterministic, but nothing measures how far a small change of the family set moves dominance or the cap.
One FASTA-only feature agrees with chicken's centromeres: 22 of its 50 interstitial telomeric arrays lie within 100 kb of a CENP-A
window (1.4 expected, P 4e-21), the telomere-repeat centromere cores of Huang 2023. Tested as a line, though, only chicken beats chance
(9 calls carry an array, 3.9 expected, P 6e-4; 6 right, 0 wrong).
**Approach.** Before trusting the cap on a genome, test how stable its dominance is (resampled seeds, or families pooled by kind). For
birds, telomeric cores could support a call, if a second bird genome shows them (other bird T2T assemblies).
**Check.** Chicken's 27 right calls rise without any wrong call rising; pig's family set stays the same under resampling.
**Code.** `sat/centromere.py` (`no_dominant_satellite`), `telomere.py` (`interstitial`).

### 42. The one-per-chromosome test near its level

**Evidence.** The test passes on any family that happens to lie once per chromosome (design.md, Known risks), and genomes near its
level flip with small changes:
- **Potato** fails from the FASTA (7 of 12, P 0.87). With `genes/` the vetoed gene-rich families change the family set and it passes
  (10 of 12, P 0.044), so two wrong calls stay `high` (chr04 3.4 Mb off, chr08 1.7 Mb off). An input that should only help adds 2
  confident errors; design.md's genes row omits them.
- **Oat.** A 335-360-bp family on every chromosome is not centromeric (cenrep 0.8.9), yet it would pass and give confident FASTA-only
  calls. *Avena* is in the polyploid set, not the benchmark.
- **Siamang** fails (18 of 26, P 0.2): on 8 chromosomes its arrays fall into two clusters more than 500 kb apart.
- **Few chromosomes.** With 4 or fewer chromosomes the test has no power (`untested`). The acrocentric test needs about 16, and no
  smaller acrocentric karyotype (fish, rodents) has been tried.

**Approach.** Report the test's P beside every call that depends on it. Measure how the outcome moves with the family set: with and
without the gene veto, and with families pooled by kind. Score *Avena* with and without `tes/` (#35).
**Check.** No confident error in potato with `genes/`; oat's FASTA-only calls scored against CENH3.
**Code.** `sat/karyotype.py`, `sat/genomes.py` (the gene veto).

### 43. Assemblies that lack or misassemble the centromeric array

**Evidence.**
- **A missing array.** Where the assembly lacks a chromosome's centromeric array, a non-centromeric satellite elsewhere on it can be
  called with confidence (cenrep 0.8.9: *B. napus* A06/A09 on an AT-rich 156-bp repeat, `high`). Gaps are reported (`call_gaps`) but
  change nothing.
- **Misassembled arrays.** The 56 ape centromeres misassembled in the NCBI v2.1 assemblies, excluded from the truth and with published
  repair intervals, are free positive controls for any signal that might flag one.
- **A collapsed array.** Pairs that span it are farther apart than their coordinates, which would mimic the cis dip at short distances
  (#29).
- **Templated arrays.** In mouse, three chromosomes share one byte-identical centromeric-end sequence (#38): unmappable for Hi-C and
  identical for the satellite line.

**Approach.** Report, per call, what an incomplete centromere looks like: a gap inside or beside the call, a short array, a call away
from its karyotype position. Test which of these flags the 56 repaired ape centromeres.
**Check.** The ape repairs flagged; T2T assemblies not flagged.
**Code.** `morphology.py` (completeness), `interpret.py`.

### 44. Holocentric genomes: one low call per chromosome

**Evidence.** On holocentric genomes cenfuse makes no confident call (0 of 196 chromosome-runs), but it still reports one `low`
"centromere" per chromosome, and nothing says the genome may be holocentric. A genome-level flag from `no_dominant_satellite` was
rejected: it is on every call of *C. elegans* and *R. breviuscula*, but also on 25 of 26 siamang and 25 of 32 EquCab3 calls. The
one-per-chromosome test fails in holocentric genomes (P 0.13-0.73), but in many monocentric ones too. What the satellite line reports on
satellite-based holocentromeres was never scored: *R. breviuscula* (Tyba arrays along every chromosome; CENH3 ChIP exists),
*Chionographis japonica* (a CENH3 track on the same assembly) and *Luzula sylvatica* (358 CENH3 domains). An LTR signature can form at
the chromosome ends of a holocentric genome (*R. breviuscula*: Athila, 5.8-fold, cores within 1 Mb of an end on 4 of 5). Hi-C on those
cores is guarded by the clean-peak rule, but a non-CRM signature above 10-fold, or a CRM seed, would place calls at the ends from `tes/`
alone.
**Approach.** A positive report needs several independent signs, each tested on the three genomes with CENH3: no dominant satellite,
a failed one-per-chromosome test, Hi-C's holocentric signature, and CENH3-domain-like arrays along the chromosome. LTR cores in a
terminal eighth should not decide, as the cis line already excludes them.
**Check.** The four holocentric controls flagged; no monocentric genome flagged.
**Code.** `sat/signature.py`, `interpret.py`.

## LTR retrotransposons

### 4. The decisive CRM test with subgenome-specific lineages

**Evidence** (CRM counts only, no truth). CRM core counts differ between subgenomes in 3 of 87 species at P ≤ 0.01 (0.9 expected), but
the marker is decisive in 2 of the 3 anyway (*B. napus*: A 8/8, B 2/8; *Euphrasia confusa*: 10/11 vs 3/11). At P ≤ 0.05 the count is
9 species (4.4 expected, P 0.03), the only hint that imbalance is more common than chance. In 5 of the 9 the marker decides anyway; in
the others it only flags (*Cerastium alpinum*: 11/18 vs 0/18). There is no evidence yet that a decisive marker is often lost, and in 0.9
the test was judged "not worth the code before truth shows" one is.
**Approach.** Test CRM concentration within subgenomes.
**Check.** Those species against CENH3 truth (#35), with random splits of diploids as the null, as in entry 3.

### 10. The LTR retrotransposon line

A check on 17 genomes with `tes/` and truth (190 chromosomes) pointed to these changes; their status in 0.12, on the same genomes,
then what 0.13 shows:

- **Show which families the CRM line rejected. Done in 0.11:** columns `seq_rejected_family`, `seq_rejected_start` and
  `seq_rejected_end`, and a caveat naming the family. 19 of the 30 calls that moved with `tes/` had such a family. Which rule earned
  which gain is still unknown: a benchmark-only switch per rule (not a user option) would attribute them.
- **One standard for replacing a `low` call. Measured, left as it is.** A core replaces a disagreeing `low` call without the
  head-to-head test a confident call must pass. Among such moves the test would have blocked, 2 are right and 1 wrong (soybean chr13,
  still confidently wrong: entry 27). That is no case for a stricter rule, and still too few cases to settle it.
- **An inconclusive head-to-head test. Left as it is.** When a core disagrees with a call and cannot overrule it, the call drops one
  level (`te_discordant`). In 0.12 the `tes/` arms held 3 such calls with truth, all in tomato, all `low` and all wrong: no case for
  keeping the level.
- **The head-to-head test acts on no benchmark call in 0.13.** Its stated cases (*B. rapa* A09, soybean chr1 and chr6, cited in
  design.md, benchmark.md and `sat/centromere.py`) no longer arise. A09 is `high` and `te_concordant` in all four accessions, soybean
  chr1 and chr6 carry no LTR flag, and the only `te_discordant` calls left are 3 `low` tomato calls. A replay with the test off would
  show whether it still protects anything.
- **`te_concordant` is a weak flag. Measured, left as it is in 0.12.** Two estimates agree if the gap between them is at most the
  wider of the two, so a wide core or call agrees with much of its chromosome. The calls that keep the flag are right 125 of 130 times,
  and the widest are no worse (35 of 35 above a tenth of the chromosome). But the 5 wrong ones are 5 of the 7 confident errors of the
  `tes/` arms (#45), so a stricter rule would gain calibration, not right calls. What the agreement is worth against chance is reported
  per genome (`asm.line_tests.tsv`).
- **In-sample.** The gains come from development genomes, and the held-out genomes of 0.12 exercise none of the branches that move a
  call (#36).

### 27. Wrong LTR cores that look right

**Evidence.**
- **Tomato.** Its signature has significant cores on all 12 chromosomes, and 5 of them miss CENH3. Without Hi-C two of those misses
  stay `medium` (CM175074.1, CM175078.1). With Hi-C they fall to `low` because the cis peak lies elsewhere, and on both the cis peak is
  the right one.
- **Soybean chr13** is similar: a CRM core 1.8 Mb from CENH3, at `medium`.
- **Potato chr08.** Its 5-element CRM core also misses CENH3 (#45).
- ***Arabidopsis*.** Its Athila signature (13.3-fold, above the level at which a non-CRM lineage may decide) has its cores beside CENH3
  on all 5 chromosomes, 0.46-1.38 Mb from the domain. It only flags because CEN180 passes the one-per-chromosome test.
- **Polyploids.** In 6 of 27 polyploids above 10-fold (three *Sphagnum*, *Lotus subbiflorus*, *Malus toringoides*, *Solanum
  stoloniferum*), the cores and the FASTA-only calls mostly disagree, with no truth to say which is wrong (*M. toringoides*' calls also
  shift with its chromosome set: 19 calls moved on its project list).
- **A lead.** In tomato the core of the youngest half of the signature's copies lies on CENH3 on 9 of 12, against 7 of 12 for the
  density core. It was chosen after seeing the result, so it is untested.

**Approach.** Nothing in the LTR line tells these cores apart so far: their fold and share look like the right cores'. Three leads:
- a core that sits beside the dominant satellite array rather than over it (*Arabidopsis*);
- the Hi-C test anchored on the cores, wherever Hi-C exists;
- with Hi-C, the cis signal tested by trans contacts also where LTR cores exist. Two tests of one signal could then move a confident
  call: tomato with `tes/` and Hi-C would go from 10 to 12 of 12. This was not built in 0.13, because of the extra machinery, not
  because of a bad result.

**Check.** More genomes with a non-CRM signature and truth (the wild tomatoes of #35); *Arabidopsis* with CEN180 masked must not move
its calls onto the Athila cores; randomized TNEG tomato stays silent.

### 31. A centromere that left its satellite and its LTR families

**Evidence.** Earlier analyses (no ChIP) place *P. annua* Pa1A's functional centromere at 163.0-169.4 Mb, where Hi-C peaks (168.0 Mb),
21 Mb from its largest satellite; 0.11 called it there. In 0.12 the CRM families learned from the other 13 chromosomes concentrated
around the satellite (145.1-148.3 Mb, 30 of their 146 elements on Pa1A). The satellite and the LTR line then agreed on the old site and
the call stayed there (`medium`, `contested`, `hic_offset`), although the Hi-C test passes on 14 of 14 chromosomes; with randomized
Hi-C the same call is `high`. A family that marks centromeres on the other chromosomes marks where centromeric families accumulated, and
a centromere that has moved leaves that behind: only Hi-C sees the new site.
Insertion age, an independent line, agrees, in two analyses:
- exploratory (`bench/pa1a_age.py`, CRM copies): median about 0 years at the Hi-C site (89% under 100,000 years; 62 copies) against
  325,000 at the satellite (39% young; 28 copies); within one family, 0 vs 157,000 years (P 2e-7); 61 of Pa1A's CRM copies under 50,000
  years old lie at 163-172 Mb, 2 at the satellite;
- the 0.13 rule's own test (signature families within 4 bins of each site): median 0 vs 620,000 years (P 2e-13).

CRM inserts into active centromeric chromatin, so young copies mark where the centromere is now. In maize, likewise, recently active
CRM subgroups colonize the present centromeres and older ones their flanks
([Wolfgruber et al. 2009](https://pmc.ncbi.nlm.nih.gov/articles/PMC2776974/)). No CENH3 data place Pa1A's centromere. In the scored
benchmark, Hi-C contests a sequence call (`hic_offset`) on one chromosome only, CHM13 CP068258.2 in two runs, and there the sequence
call is right. It has no `tes/` ages, so the tie-breaker has never been scored.
**Done in 0.13** (`ageline.py`). Where validated Hi-C's strong peak lies far from a `high` call, copies of the chromosome's signature
families (learned from the other chromosomes) that are younger at the peak than at the call become the second line that lets Hi-C move
the call (`moved_with_age`). The test is a one-sided Mann-Whitney at P ≤ 0.01, with at least 3 dated copies at each site and both sites
widened by 4 bins. It applies only in a genome where the premise holds: on the other chromosomes, the signature's younger copies lie
in the cores more often than its older ones (Fisher, P ≤ 0.01; *P. annua* leaving Pa1A out: P 5e-10). As a locator on its own, the core
of the young copies does worse than the density core (on CENH3 for 128 vs 155 of 172 chromosomes with truth; `bench/young_core.py`), so
age only breaks a tie between sites that another line raised. In the benchmark it acts only on Pa1A, which 0.13 moves 20.7 Mb to the
Hi-C site (`high`). The README counts this among cenfuse's strengths ("no ChIP yet"), which runs ahead of the evidence.
**Check.** CENH3 ChIP or CUT&Tag of *P. annua*. Genomes with known evolutionarily new centromeres and Hi-C: horse chr11 (#50), the
orangutan chr10 neocentromere polymorphism (de Gennaro 2026). Within-species centromere tables (soybean, tomato: #35). The randomized
Hi-C controls stay silent.

### 45. Confident calls at the edge of a wide LTR core

**Evidence.** A call agrees with an LTR core when the gap between them is at most the wider of the two (`te_concordant`, entry 10),
and an agreeing core lifts both the honesty cap and the one-per-chromosome cap (`sat/centromere.py`). That is how *Ae. tauschii* and
potato, whose satellites fail that test, keep confident satellite calls. Of 130 such calls with truth in the 0.13 FASTA + `tes/` arms,
125 are right and 5 are confidently wrong, and those 5 are 5 of the 7 confident errors of those arms:

| chromosome | confidence | call (Mb) | LTR core (Mb) | CENH3 (Mb) |
|---|---|---|---|---|
| maize Mo17 CM039154.1 (chr5) | `high` | 107.19-107.40 | 103.54-108.83 | 103.5-105.2 |
| *Ae. tauschii* Chr5 | `high` | 199.36-201.42 | 194.02-201.03 | 195.4-198.1 |
| *Ae. tauschii* Chr7 | `medium` | 340.27-340.50 (45 kb of satellite) | 340.47-349.60 | 343.6-344.9 |
| potato chr08 | `high` | 15.54-16.92 | 16.46-17.41 (5 elements) | 18.6-18.7 |

The fifth is tomato CM175074.1 (`medium`), whose core misses CENH3 too (#27). In three of the four in the table the core holds CENH3
and the call sits at its edge. Maize chr5's caveat is empty and its `basis` says the CRM core agrees: 118 of the chromosome's 158
signature elements lie in the core, 11 in the call. It is the plants' only confident error from the FASTA alone, and it has been one
since 0.8. Exploratory, not a rule: where the core is at least 3 times wider than the call, the counts are 11 right and 4 wrong;
otherwise 114 right and 1 wrong. 0.9 already named the pattern: "damage occurs where a CRM-rich pericentromere is broader than the
centromere". With Hi-C, potato chr08 falls to `medium` because the cis extreme lies within one bin of CENH3 (`cis_elsewhere`), but a
confident call is not moved by one tested line against another. Maize has no Hi-C (only HiChIP), so its chr5 has the shape of Pa1A
(entry 31) without a Hi-C peak to move it.
**Approach.** Judge agreement by chance, not by gap. `linetests.py` already computes each chromosome's chance that two placements
agree. A call at the edge of a core much wider than itself agrees little beyond chance, so it should not be raised; or it should meet the
head-to-head test that a call beside a core meets. Inside an agreeing wide core, the young copies' core (entry 31) is a candidate
locator, untested.
**Check.** The five calls fall to `low` or move onto CENH3; none of the 125 right calls drops below `medium`. Held-out CRM-plateau
genomes (wheat, the soybean lines, *Ae. geniculata*: #35), not these five, decide the rule.
**Code.** `sat/centromere.py` (`te_concordant`), `sat/te.py` (`agrees`), `linetests.py`.

### 46. CRM decides by its name

**Evidence.** CRM is the one lineage that decides from 3 significant cores on, with no fold level (`sat/signature.py`); other lineages
need 10-fold. The reason given is that every CRM core short of a majority lay on CENH3 (24 of 24). But those 24 come from weaker or
thinned annotations of genomes whose CRM is centromeric (rice, potato, maize, *Brachypodium*), so the sample holds no genome with a
non-centromeric CRM. Such genomes exist:
- in tobacco, CENH3 binds Tekay and avoids CRM, which clusters beside the centromeres (cenrep's literature notes);
- in eggplant, CRM lies outside the centromeres (same notes);
- soybean chr13's CRM core lies 1.8 Mb from CENH3 (entry 27).

In v6, a rule that let a CRM clade with more cores than chance take the marker ahead of a lineage with cores on most chromosomes
lowered two right `high` *Arabidopsis* calls (TAIR10), and it was rejected. 0.12's seed puts a CRM candidate with 3 cores ahead of any
lineage, with a floor that excludes that case, and was not retested. Mammalian ERVs that LTRquest labels CRM (14 in CHM13) could seed
one.
**Approach.** Give CRM the same test as other lineages (a fold level, or the decisive test measured in each genome), or keep its
standing only where its cores pass the one-per-chromosome test. Test on Solanaceae with CENH3 (*S. verrucosum*, tobacco).
**Check.** No right call lost in the CRM-decisive genomes of the benchmark; no tobacco or eggplant CENH3 domain overruled by CRM.
**Code.** `sat/signature.py` (`seed_of`, `TE_SEED_MIN_CORES`).

### 47. Centromeric retroelements the LTR line cannot read

**Evidence.** The LTR line uses classified LTR elements only. Two genomes show what that misses:
- ***Cryptococcus*.** 89-94% of LTRquest's elements lie in the CENP-A domains, but none has a REXdb lineage (unknown, chromo-outgroup),
  so from the FASTA and `tes/` cenfuse places 0 of 14; only Hi-C reaches them. v6 tried pooling unclassified elements three times. The
  guarded form placed JEC21 10 of 14 and a hybrid's D subgenome 10 of 13, all `low` and none wrong, and it never lifted the honesty
  cap. But it still qualified in 9 fish genomes without truth and in *A. suecica*, and its 20-kb bin floor cost *O. punctata* a CRM call.
  design.md records only the first form's failure.
- ***Chlamydomonas*.** Its centromeres are clusters of ZeppL, a LINE, and two CENH3 paralogs bind ZeppL. LTRquest finds 669 elements,
  613 of unknown lineage, and no line sees the LINEs: 6 of 17, all `low`. It holds 22 of the chromosome-runs where no line is right.

**Approach.** Lineage labels that cover fungal Ty3 and chromoviruses (the "chromo-outgroup" elements may be the centromeric ones). A
non-LTR annotation read by the same family-level signature (families, not classes, carry the test). Both need a second genome with truth
before they may act.
**Check.** *Cryptococcus* (both nuclear states, without Hi-C) and *Chlamydomonas* placed; fish and *Oryza* unchanged; shuffled positions
give no signature.
**Code.** `sat/te.py` (the readers take `LTR_retrotransposon` features only), `sat/signature.py`.

### 48. How much the LTR line depends on the annotation

**Evidence.** The decisive test needs significant cores on most chromosomes, and potato's sits at its threshold: 8 of 12 significant
signature cores in 0.13, where a majority is 7. benchmark.md's annotation table says it was measured on 0.11, but its numbers equal v6's,
measured with the lineage line. There, potato's 214 LTRquest CRM elements gave exactly a majority; half of them, EDTA's intact set or
DANTE_LTR each lost the gain from `tes/`; and rice fell from 8 cores to 3, 1 and 1 at 50, 25 and 10% of its elements. The 0.12 family
signature, which recruits from the other chromosomes, was never thinned. The README prescribes "exactly this command" with `--proteins`,
which no benchmark annotation used, and the LTRquest version (1.0.1 in v6) is recorded nowhere.
**Approach.** Rerun the thinning on 0.13 with LTRquest only (v6 `bench/td_thin.py` thins elements by hash). Report the depth at which
each genome's signature stops deciding. Record the LTRquest version and flags in the outputs. Benchmark the README's command with
`--proteins`, or drop the flag.
**Check.** The annotation table regenerated from 0.13 code.
**Code.** `sat/signature.py`, `teline.py`.

## Hi-C

### 13. Bird microchromosomes in Hi-C

**Evidence.** Hi-C bins finer than 250 kb are capped by the largest chromosome (4,000 bins). Chicken's chr1 (200,044,509 bp) needs
4,001 bins at 50 kb, one over the cap, so its 2.5-7-Mb microchromosomes stay out of the Hi-C line. The genome cap (50,000) would allow
about 21,000. The chicken library is a prefix of the data: the first 40 M read pairs of SRR15423486, 26.1 M unique pairs.
**Approach.** First, a one-constant experiment: a per-chromosome cap set by memory rather than by 4,000 bins. Check which readers use the
dense per-chromosome block (the fold, balancing, the cis profile), because a sparse or banded fold scan lifts the cap only if they all
can work without it.
**Check.** Chicken with deeper Hi-C; calls on the microchromosomes against the CENP-A windows.
**Code.** `hic/params.py` (`MAX_CHROM_BINS`), `hic/bins.py`, `fold.py`.

### 14. The Rabl axis where the fold has no power

**Evidence.** The trans arm alignment passes in *C. neoformans* (anaphase, 12/14) and *M. sympodialis* (6/8), and Hi-C's prior-free calls
there are right on 10 of 10 and 7 of 7. But the arm fold, which must lean the same way before the alignment may act (the holocentric
guard), has no power in these small nuclei (1/14; 0/8), so the alignment alone raises nothing. In holocentric *C. elegans* the alignment
reached 2/5 (P 0.095), and a position-matched null that rejects it also rejects barley. In G1 *C. neoformans* the prior-free peaks sit at
the chromosome ends, so the alignment's partners cannot be placed (0/14); anchored on CENP-A, 8 of 14 chromosomes still align (P
0.0014), so the alignment weakens without the cluster rather than vanishing.
**In 0.12 and 0.13.**
- AT troughs that agree with Hi-C raise *M. sympodialis*' 8 calls and 4 of *C. neoformans*' to `medium` (`at_pair`, 0.12).
- The cis signal tested by trans contacts places G1 *C. neoformans* at 12 of 14 (`low`, entry 26).
- A trough paired with a prior-free peak brings *C. neoformans* in anaphase from 10 to 13 of 14 (`at_denovo_pair`).
- *C. tropicalis* (3C-seq, 11.4 M pairs): its prior-free peaks are right on 6 of 7, but with 6 chromosomes the engine's Rabl test cannot
  reach P ≤ 0.01 (floor 0.0156), and nothing passes (0 of 7).

**Approach.** Three directions:
- Place the partners on the sequence line's confident calls, as a corroboration test. A prototype exists (v10.1
  `explore/anchored_align.py`: G1 8/14, P 0.0014; anaphase 12/14, P 2e-6; randomized 0/14); a per-genome test like the gate's comes first.
- A compartment-adjusted alignment, its expected value taken from each bin pair's cis eigenvector classes (as in saddle plots). Beware:
  restricting peaks to B-compartment bins helped apes and broke *P. annua*, because the sign of E1 at centromeres is opposite in plants
  and apes (cenhic 0.2).
- The fold read at the 5-kb rung in small nuclei.

**Check.** These must pass: barley, *Neurospora*, the yeasts. These must fail: *C. elegans*, *Rhynchospora*, *Bombyx* (small-chromosome
holocentrics), and randomized Hi-C. The stored Hi-C arms must be unchanged except where predicted. Truth and Hi-C from *Cryptococcus*
relatives, *Toxoplasma*, *T. brucei* and other *Candida*. A Rabl tissue for the grasses: *B. distachyon* root tips are Rabl and leaves
are not, and barley endosperm loses Rabl. The primary Bd21 library (CRR611509, 289 M pairs) was never used, and the README gives no
advice on tissue or cell-cycle stage.
**Code.** `fold.py` (`align_profile`), `fuse.py` (`fold_pair`, `align_pair`).

### 15. Centromeres at the edge of a terminal unmappable array

**Evidence.** Chicken's acrocentric microchromosomes carry their centromeric array at the very end, unmappable for 1-1.25 Mb. Their
Hi-C peaks lie at its inner edge, which in mappable bins is the end, so they stay `telomeric`, and chicken's calls stay `low`.
**Done in 0.12** (`hic_end`). Telomere clustering lifts both ends of a chromosome, a centromere only one. The anchored score is about
twice as high at the end holding the call in goat (28 of 30), sheep (22/23), horse (18/19) and mouse (19/19). So a call at the higher end
rises a level where that holds on most of the other chromosomes (binomial, P ≤ 0.01). No randomized library passes. Of its 56 raises,
the 55 scored are right 54 times.
**Still open.**
- Chicken does not pass (12 of 19), so its calls stay `low`.
- The one raise that is not right is sheep CP162266.1 (X), `high` and 83.4 Mb from another caller's locus. The karyotype check agrees
  with cenfuse's call (CI 0.054, acrocentric), so the other caller is probably wrong; sheep has no functional truth (#35, #37).
- design.md rejects "the other end's height" for telling a right end peak from a wrong one, but adopts this test, which also compares
  the two ends. The difference (a per-peak discriminator vs a genome-level binomial over calls that only raises) should be stated there.

**Code.** `hicfuse.py` (`end_test`).

### 29. What the contact crossing a centromere is

**Evidence.** Piled up on CENH3, tomato's two flanks contact each other more than DNA at the same distance elsewhere (2.4-fold within
0.25 Mb, back to background by 1 Mb; no saved file holds the 2.4, and the figure's numbers give 1.79 against 0.97 at random
positions). Potato shows a smaller central peak between deep dips. Both lose contact with the arms out to at least 6 Mb (a mean effect
only: per chromosome, the lowest contact with the arms is not the centromere, tomato 0/12, potato 0/11, rice 0/12). *S. pombe* has
inter-arm contacts close to its centromeres that depend on heterochromatin ([Mizuguchi et al. 2014](https://doi.org/10.1038/nature13833)).
A fold of the two flanks, or a domain insulated from the arms, could give all of these, the sign depending on the domain's width against
the 1 Mb the profile reads. That is a hypothesis: nothing here perturbs it. A plain artifact must be ruled out first. Pairs that span a
collapsed or missing array are farther apart than their coordinates, which gives an apparent depletion at short distances that fades
with distance; true insulation persists. The dips (potato, *Verticillium*) were never tested for that.
The chance the test uses comes from contact maps simulated with each bin's coverage and the chromosome's mean contact at each distance.
On moss's re-paired control map the dip lies at the LTR cores on 11 of 26 chromosomes, where the simulation expects about 5 (P 0.011, just
short of passing). Among 20 control tests a value this small turns up about one time in five. Still, if the simulation runs short in
some genomes, balancing each simulated map again (instead of reusing the real map's balancing) would be the first thing to try.
**Approach.** Rerun the profile at several box sizes and look for genomes where the sign changes with the box. Compare the dips with
the assembly's gaps and collapsed arrays.
**Check.** Hi-C of heterochromatin or condensin mutants, or of other cell-cycle stages, in a genome with this signal.
**Code.** `hicfuse.py` (`cis_profile`, `cis_chance`).

### 30. A Hi-C test with too few confident references

**Evidence** (0.12). The Hi-C test took the `high` calls as references (with the `medium` ones if too few), and fell back to the `low`
calls at clean peaks only when the genome had no `high` or `medium` call at all. Cucumber 9930 (held out) without `tes/` has 7 `low` calls,
and Hi-C was tested on them. With `tes/` the Tekay signature raises 2 calls to `high`, so 2 references remained, too few, and the test was
`untested` on all 7 chromosomes. No call differed: on the `low` references Hi-C had changed no level either.
**Done in 0.13.** When the confident references are fewer than the test needs (4), the `low` calls join them, with agreements counted at
clean peaks only for the `low` ones (`hic_gate_refs` `mixed`). In cucumber Hi-C now passes on 3 of 7 chromosomes and changes no call. No
benchmark call with truth changes, and the randomized controls stay silent.
**Check still open.** A genome where this lets Hi-C act, with truth. The older `low`-reference rule (v5) acts in one genome only, the
siamang, which has no CENP-A truth. There it raises the two calls that confidently disagree with the satellome locus (CM054514.2 by 1.0
Mb, CM054529.2 by 15.0 Mb). The siamang's DiMeLo-seq track (#35) would test it.
**Code.** `hicfuse.py` (`gate`).

### 32. Calls in a terminal eighth that an interior line contradicts

**Evidence.** The cis signal tested by trans contacts (entry 26) is read outside the terminal eighths, because there both readouts see
the chromosome ends. So it has no say on a call inside one. Tomato without `tes/` keeps 2 wrong `low` satellite calls 8.6% and 11.8% from
an end (CM175073.1, CM175075.1), while its interior extremes lie on CENH3. Moving such a call would be right in tomato and wrong for an
acrocentric chromosome among metacentric ones, which the line cannot see.
**Approach.** A second line that can see the ends: the which-end test (the anchored score at the two ends) where it passes, or an LTR
core.
**Check.** Tomato moves its 2 calls; the acrocentric mammals (goat, sheep, mouse) keep theirs.

### 33. Cis and trans at one locus as a confident call

**Evidence.** Where the cis signal passes its test by trans contacts, a chromosome whose own trans peak lies at its cis extreme holds two
readouts at one locus, as `te_hic_pair` holds an LTR core and Hi-C. Such calls stay `low` in 0.13. In G1 *C. neoformans*, 10 of the 12
pairs lie on CENP-A and the other 2 lie 2 bins (19-20 kb) from it, wrong by the 1-bin rule, so as `medium` they would be confident
errors. The FASTA-only AT trough is a finer interval (16-40 kb). It overlaps CENP-A on 12 of 14 chromosomes and agrees with 11 of the 12
right cis calls, and on CP003829.1 it lies on CEN10 where the cis bin is 19 kb off. AT troughs are paired with Hi-C's prior-free peaks
(which sit at the ends in G1: 1/14), but not with the cis-trans line.
**Approach.** Pair the cis-trans line with the AT trough: two lines that can test each other, with the trough's interval for the call.
Or `medium` for a pair whose bin is coarse relative to the centromere. A sub-bin fit of the trans centromere spot is another route: in
*S. cerevisiae*, Centurion's 2D-Gaussian fit has a mean deviation of 2.5 kb at 40-kb bins.
**Check.** Small genomes with truth at 10-kb bins (*V. dahliae* and *N. crassa*, where both lines pass); the randomized controls stay
silent.

### 49. Hi-C cannot see inside megabase arrays

**Evidence.** With about 100 M pairs at MAPQ 30, the CENP-A domain lies in bins without a uniquely mapped contact on 139 of 157 T2T mammal
chromosomes (blind blocks of median 2.5 Mb in chimpanzee to 19 Mb in gorilla). The anchored peak overlaps CENP-A on 5 of 157 and lies
within 4 bins on 32, a median 3.2 Mb away, at the edge of the blind block (20 of 22 in human). So a passing Hi-C test certifies that Hi-C
sees the array, not the kinetochore: the pre-registered "peaks within 4 bins of CENP-A wherever the test passes" failed in four apes, and
offsets inside an array are invisible (design.md, Known risks). MAPQ ≥ 1 changes nothing, because near-identical copies map at MAPQ 0.
Never tried: contacts assigned through array-unique k-mers, multi-mapping-aware assignment, or long-read contact data (Pore-C). v2
predicted those "should move the maxima onto the CENP-A domain".
**Approach.** Remap one ape library keeping multi-mapped reads (or through array-unique k-mers) to measure what is recoverable. Add
Pore-C as an input type if public data exist for a T2T genome with CENP-A truth.
**Check.** Peaks inside the CENP-A domain on most ape chromosomes; randomized libraries stay silent.
**Code.** `hic/io.py`, `hic/bins.py` (masking), `hicfuse.py` (`visible_window`).

### 50. Satellite-free centromeres in mammals

**Evidence.** Horse chr11 (CM128059.1) is the classic evolutionarily new centromere and the benchmark's only satellite-free mammal
centromere with truth, and it gets no call in any arm (neither do EquCab3's chr11 and chr31). In 0.1.0 its anchored Hi-C peak lay inside
the CENP-A domain but at engine level `medium`, so no call was made. The rule for that case exists now (`hic_medium_peak`, a `low` call),
but horse Hi-C has failed its test at every depth since v3's mappable windows ("no window rule tried restored it without losing the ape
and macaque gains"). The library comes from another horse, and CENP-A on chr11 slides about 500 kb between horses. The cis signal works
in no mammal. A gap-rank test (does the call's gap have the top flank score among its chromosome's gaps?) reached P ≤ 0.01 in horse (9 of
11), pig (8/11) and rhesus (13/16), where the gate fails, and found 0 of 14 on randomized Hi-C. It was never adopted (v3, analysis only).
**Approach.** The same-animal horse set and the multi-horse CENP-A union (#35). Then a test for Hi-C where the gate fails (the gap-rank
test) with a stated role: place a `low` call on a chromosome without a sequence call.
**Check.** Horse chr11 placed; no call changes in the other mammals; randomized and holocentric controls silent.
**Code.** `fuse.py` (`hic_medium_peak`), `hicfuse.py` (`gate`).

### 51. Hi-C tests at their threshold

**Evidence.** The tests do not lose power gradually as depth falls:
- both orangutans pass at full depth and fail at half, where the bin doubles to 500 kb;
- pig fails at 100% and 50% of its pairs, passes on 20 chromosomes at 25% and 10% (500-kb and 1-Mb bins), then fails at 5%, with its calls
  unchanged throughout;
- rice passes at 50% and fails at 25%.

A test whose verdict improves when reads are discarded suggests that the bin decides as much as the signal: each depth chooses its own
bin, and agreement is counted in bins. A split-half flag was rejected (stable genomes agreed 13 of 13, borderline ones disagreed with each
other, and Hi-C time doubles). The share of several (about 8) random halves that pass, proposed with it, was never tried. The thinned
series and the robustness claims come from older code. The libraries thinned to 2-25% "passed none in 0.11 and were not rerun".
benchmark.md says barley falls back to "3 confidently wrong" FASTA-only calls, which 0.12's one-per-chromosome test made `low`.
**Approach.** Rerun the depth series on 0.13. Check the gate's chance model across bin sizes at fixed depth. Report the share of random
halves that pass as depth guidance, not as a threshold.
**Check.** Pig's passes at lower depth explained; depth guidance in the README from 0.13 numbers.
**Code.** `hic/bins.py` (the bin choice), `hicfuse.py` (`gate`).

### 52. Hi-C alone, untested

**Evidence.** In a genome without sequence calls, a clean `high` prior-free peak becomes a `low` call that nothing tests. The strongest
trans cluster may be something other than the centromeres, and a holocentric genome without sequence calls would get such calls too
(design.md, Known risks).
- **Sharp single-bin peaks.** Potato's 7 single-bin peaks at 12-27 times the chromosome median are `high` and wrong; one carries `hub`,
  which needs both 2× coverage and 2× trans contacts.
- **An effect-size readout.** The median peak over the median (≥ 2) separated *C. elegans* (1.21) and the randomized libraries (1.05-1.51)
  from yeast (5.96) and *Plasmodium* (5.40). It was rejected as a test of nucleus type and never tried as a guard on Hi-C-alone calls.
- ***C. tropicalis*** (entry 14): right peaks, and no test that can license them.
- ***Z. tritici*** is right on 19 of 21 with Hi-C, all `low`. 8 of its 21 chromosomes are accessory, which may dilute every "most
  chromosomes" test.
- **Dormant mis-mapping controls.** The Hi-C engine has two that `cenfuse run` never runs, and neither is documented. One is the
  read-class split (`--repeats`: trans contacts of reads in repeats against those in unique sequence); as criterion C3 of the first
  benchmark it passed in 6 of 6 held-out mammals. The other is `simulate`, cis-only pairs simulated from the assembly.

**Approach.**
- Add the effect-size guard and a single-bin-hotspot flag on Hi-C-alone calls.
- Run the read-class split by default, with the satellite arrays and LTR cores as the repeat class: a free control in every genome.
- Run the majority tests on core chromosomes where the accessory ones are known.

**Check.** Holocentric and randomized libraries give no Hi-C-alone call; the yeasts and *Plasmodium* unchanged; *C. tropicalis* and
*Z. tritici* scored.
**Code.** `fuse.py` (`hic_denovo_only`), `hic/qc.py` (flags), `hicfuse.py` (`repeats=None` where the engine is called), `hic/simulate.py`.

### 53. Where the cis-trans line passes, every extreme is trusted

**Evidence.** The cis signal tested by trans contacts acts on every chromosome of a genome that passes (leave one out), not only where that
chromosome's own trans peak meets its extreme.
- **Off extremes.** Where it passes, some extremes are far off: *N. crassa* 911 kb and 82 kb on two of 7 chromosomes, *V. dahliae* 112 kb
  on one of 8. No call changed, only because those chromosomes already had Hi-C or AT calls, which the cis rules never touch.
- **An untested false positive.** In v10, anchoring trans contacts at the other chromosomes' cis extremes agreed in rice on 8 of 12 (P
  2e-10), none of them right. The line now reads outside the terminal eighths, and rice's extremes then lay in a terminal eighth on 12 of
  12. But it runs only where the gate fails, so it has never run on rice, soybean or the apes, and design.md's "none seen outside the
  terminal eighths" is untested there.
- **Where LTR cores exist** the trans test is not run (entry 27).

**Approach.** Let the line act on a chromosome only where its own pair holds (the pairs it already records). Run it as a diagnostic on
the gate-passing genomes with truth.
**Check.** The off extremes of *N. crassa* and *V. dahliae* do not act; rice, soybean and the apes pass nothing; randomized libraries
silent.
**Code.** `hicfuse.py` (`cis_trans_line`), `fuse.py` (the cis rules).

### 54. The support tier's evidence

**Evidence.** The tier lets Hi-C that finds centromeres beyond chance, but not on most chromosomes, raise the calls its significant peaks
agree with (and, since 0.13, choose a competing cluster). Its evidence, "right 61 of 64" (design.md), counts agreement with another
satellite caller as right in goat (12/12), rhesus (8/8) and giraffe (4/4), and goat's 8 raises to `high` stand on that caller. The
pre-registered held-out bar (≥ 90% of raised calls right) was missed: 14 of 16, 88% (pig's two `high` misses). Its known risk, raising a
wrong call that sits on a non-centromeric trans cluster (potato: 1 of 64 stored agreements), left design.md in the 0.10.1 rewrite.
**Approach.** Score its raises against functional truth only (cattle and goat CENP-A data would come first; #35). Restore the risk to
design.md.
**Check.** The tier's raises right on ≥ 90% with functional truth, or its role narrowed.
**Code.** `fuse.py` (`hic_support`).

## AT troughs

### 28. Where AT troughs do not work

**Evidence.** The test passes in *P. falciparum*, *M. sympodialis* and *N. crassa*. It fails elsewhere, and its failures are of two kinds:
- **Right troughs, failed test.** Each chromosome's deepest trough lies on the centromere in *V. dahliae* (8 of 8) and *C. neoformans*
  (12 of 14), but the genome-level count of the N deepest troughs fails: *V. dahliae* has three deep troughs Mb apart on one chromosome,
  and a *C. neoformans* centromere holds two sub-troughs. These 20 calls are computed and discarded; AT troughs hold 117 of the 231
  chromosome-runs of headroom (113 where only the trough is right, a count inflated by the 250-kb tolerance without Hi-C, #37).
- **Wrong troughs.** In *Z. tritici* (0 of 21 overlapping) and *C. tropicalis* (0 of 7) the troughs miss, so a test that passed there
  would make wrong calls.

Point centromeres say nothing (yeast: 125 bp), and genomes with few chromosomes cannot reach the level (*S. pombe*: `untested`,
though its 3 troughs are right). The test passes on about 4% of genomes with shuffled GC blocks, as its 5% level allows, and then gives
`low` calls. It also passes on *Chlamydomonas* (P 0.0125, 500-kb window) with its troughs right on 6 of 13. It does not run there only
because most chromosomes have a sequence call: a latent false pass. *P. falciparum* NC_037281.1 is its one wrong call (482 kb off,
`low`).
**Tested in 0.13, not changed.** Counting the troughs of one AT-rich run once (each trough's half-depth run excluded before the next)
changes no genome. *C. neoformans*' 30-65-kb centromeres hold two sub-troughs with GC-normal sequence between them (CEN11: 0.865 and
0.925 Mb), so the run breaks, and *V. dahliae*'s troughs are Mb apart. Also tried before 0.12: CUSUM on the GC deficit (*M. sympodialis*
0 of 8) and choosing the window by the largest gap between the first and second trough (wrong windows in *C. neoformans* and *S. pombe*).
In *C. neoformans* the cis signal tested by trans contacts now places them (entry 26), and an AT trough agreeing with a weak prior-free
peak places 3 more (`at_denovo_pair`).
**Approach.** A per-chromosome test: each chromosome's deepest trough against its own chance, combined over chromosomes, instead of
the genome's N deepest. A second feature with the trough would also help: a gene desert, or without `genes/` a run of low k-mer
complexity.
**Check.** *V. dahliae* and *C. neoformans* pass from the FASTA; *Z. tritici*, *C. tropicalis*, *Chlamydomonas*, the holocentric genomes
and the shuffled genomes do not.
**Code.** `atline.py`; `v12/bench/at_diag.py` and `v12/results/at_merge.tsv` hold the per-chromosome numbers.

## Methylation

### 16. The CHG line has changed no call

**Evidence.** Only *P. annua* and *Arabidopsis* have CHG data here, and across all 3,773 chromosomes of the 0.12 release runs the line
supported no call and backed no move; it has never acted on scored data. *P. annua*'s five samples come from other plants, so they do
not decide. *Arabidopsis* finds no signature: its expected domain, from genome size, is 315 kb, while its CEN180 arrays span 2.2-3.5 Mb,
so the reference ring (1.5-4 times the domain) falls inside the arrays. The signal is real in *P. annua*: the dip is the chromosome's
lowest point on 12 of 14. Four things stand in the way:
- **The domain size D** comes from a maize prior applied to every genome (22.2 Mb over Mo17's 10 centromeres, scaled by genome size with
  slope 0.97 from Plackova et al. 2021). It sets the window, the ring and the candidate spacing. A D measured in the data (from the
  anchor calls) was the lead in v4, not adopted for want of a second same-plant case; v5 rejected "a data-derived CHG domain size" with
  no numbers kept. D also takes the number of chromosomes listed, so `asm.chromosomes.txt` must name all of them, and it is computed in
  two places (`pipeline.py`, `methfuse.py`).
- **A per-site floor.** The CHG engine drops sites with fewer than 3 calls (MIN_CALLS_SITE), the floor the CpG line removed when pooling
  per bin kept rhesus and crab-eating macaque kinetochores at 1/16 of the reads.
- **Fraction-only input.** The CHG reader takes modkit or ENCODE bedMethyl only, so the public *Solanum* tracks (4-column
  fraction-only bedGraphs) cannot be read.
- **No data.** Same-plant CHG with CENH3 truth exists (#35: *S. verrucosum* CPC54, potato DM FAST5, rice T2T-NIP) and was never used.

**Approach.** D from the genome's own confident calls (the median call width) where there are enough, else the prior, reporting which
was used; CHG pooled per bin like CpG; fraction-only tracks read without genotyping.
**Check.** CHG from the assembled plant of a genome with CENH3 truth (CPC54 first).
**Code.** `meth/params.py` (D_REF_BP, D_REF_G, D_SLOPE, MIN_CALLS_SITE), `meth/track.py`, `meth/io.py`, `methfuse.py`.

### 55. CpG kinetochores outside primates

**Evidence.** The CpG line was designed on the primates it is benchmarked on (rhesus, *Arabidopsis* and HG002 held out), and its
specificity rests on rotated tracks and one plant. The gaps:
- **Negatives and held-out sets.** No mammal without CDRs serves as a negative: mouse is predicted to fail, its ONT reads carry no
  modification tags, but a HiFi 5mC track would do. No grass, other vertebrate or holocentric genome has been run. The HGSVC3 CDRs (67
  genomes) are unused as a held-out set, and giraffe's ONT methylation (21 modBAMs) and rhesus's are unused for CDR truth.
- **A gate without a majority.** Gorilla passes on positional dips at 7 of 25 references, and a kinetochore then needs only the shuffle
  test on its own call.
- **Small arrays.** Arrays hypomethylated as a whole (the bonobo and Bornean orangutan minicentromeres) are never refined, since a call
  needs 10 bins with data, and the line cannot tell a wrong array (#38).
- **Old tracks.** Chimpanzee and siamang fail the gate on R9 Guppy tracks (positional dips at 1 of 25 and 4 of 22).
- **Provenance and power.** CG-only samples decide without a provenance test. With five chromosomes the gate has little power: 2 of 4
  references at chance 0.05 gives P 0.014.
- **A perturbation.** CHM13 clones whose methylation was changed by dox moved CENP-A (PRJNA1270043: fiber-seq and DiMeLo-seq), a direct
  test of the line's premise.

**Approach.** HGSVC3 as a held-out set; a mouse HiFi 5mC negative; the perturbed clones (does the called kinetochore follow CENP-A?); a
per-call positional requirement, or a recorded reason for its absence.
**Check.** Kinetochores within 50 kb of the HGSVC3 CDRs; mouse fails the gate; rotated controls silent.
**Code.** `cdr.py`.

## Inputs, chromosome set and other features

### 8. Large unplaced scaffolds and dot chromosomes

**Evidence.** On 121 polyploid assemblies with known chromosome lists, the length rule keeps 31 sequences outside the lists (1.5-117 Mb,
median 7 Mb) in 12 assemblies. Extra sequences dilute the prevalence tests: on the project lists *Cissus* goes from 0 to 14 `high` calls
and *Mimulus* from 9 to 25, and *Malus toringoides* moved 19 calls. Five list assemblies were never scored on 0.9's sets (*Cissus
quadrangularis*, *M. toringoides*, *Oryza schlechteri*, *Pontederia crassipes*, *Armoracia rusticana*). The rule fails at three other edges:
- it drops a dot chromosome that stands alone after a tenfold gap (*D. melanogaster* chromosome 4 without the Y);
- it drops dot chromosomes followed by a tapering scaffold tail (a simulated bird-like karyotype of 37 chromosomes and 300 scaffolds
  tapering from 0.4 Mb: the rule keeps 31 and drops six dots of 0.5-1.4 Mb; with a clean gap it keeps 37);
- it cuts as a tail small chromosomes that themselves span a tenfold range.

It keeps organelles: the yeast mitochondrion (NC_001224.1, 86 kb) is one of the 17 "chromosomes" of every *S. cerevisiae* run and enters
every denominator (telomere ends, the majorities, the one-per-chromosome N). Two limits on sequence counts, MAX_CHROMOSOMES 50 (a
warning) and MAX_SEQUENCES 5,000 (no one-per-chromosome test), express one concept 100-fold apart.
**Done in 0.13, in part.** The longest sequences left out of the set (up to 100, at least 200 kb) that hold arrays of the genome's
telomere unit at both ends are named in the report and in `left_out_telomeric`. These are complete molecules, such as a dot chromosome
after a tenfold gap or a B chromosome; they are flagged, never added. The flag fired in none of the 260 release runs, and it cannot act
in a genome without a telomere unit (potato, *S. pombe*: #56); 0.9's version also asked for a tandem array, so a dot whose telomeres are
unassembled is missed.
**Approach, still open.** The other edge: large scaffolds kept in the set without a telomere, a centromere-like array or genes (Hi-C
contacts are another signal). An organelle test from the sequence (the organelle SSU context `rdna.py` already reads is one route).
**Check.** The 121 lists as truth; *D. melanogaster* and bird karyotypes, with the simulated karyotypes as regression tests; the 383
assemblies of the rule's evidence as a no-change control.
**Code.** `inputs.py` (`chromosome_set`, `left_out`), `meth/track.py` (where the length rule lives today).

### 56. Telomeres and rDNA arrays the benchmark does not count

**Evidence.** benchmark.md names telomere misses only in chicken and *Bombyx*, but the release runs find:
- no unit at all in potato (0 of 24 ends) and *S. pombe* (0/6);
- 7 of 34 ends in yeast (TG1-3 found as ACACCC; the 34 include the mitochondrion's 2);
- 3 of 14 in barley, 3 of 10 in *R. breviuscula* and 8 of 10 in *Arabidopsis* (its NOR arms).

Without a unit, `left_out_telomeric` cannot act. NOR recall was never measured: goat has arrays on 3 of its 5 NOR chromosomes. The 45S
scan needs 3 motif copies, so a collapsed 2-unit rDNA locus is missed (*Cryptococcus* JEC21). Four T2T mammal ends hold an array that
the truth rule (terminal 10 kb) does not count, and they were never inspected.
**Approach.** Report telomere recall per genome; find why potato has no unit (assembly ends without repeats, or detection); a 2-copy
rDNA rule, using the organelle context the 45S scan already reads.
**Check.** Telomere and NOR recall in the benchmark tables.
**Code.** `telomere.py`, `rdna.py` (`MIN_UNITS`).

## Several plants

### 22. Pa1A in other P. annua plants

**Evidence.** Young CRM insertions of four of the five other plants cluster 9-15 Mb from the functional centromere that earlier analyses
found on Pa1A (163.0-169.4 Mb), which suggests other centromeres; J10's lie inside it. C11 and J10 show a plant-only methylation dip at
their clusters. 0.13 calls Pa1A at the Hi-C site (167.5-170.0 Mb, `high`, `moved_with_age`: entry 31). The plant verdicts there are 3
`shared` (J10, P3, R1) and 2 `unclear` (C11, I3). Over all 70 plant-chromosomes they are 58 `shared`, 6 `unclear`, 2 `possible_shift`
(J10 Pa4B, J10 Pa5A), 2 `call_absent`, 2 `no_data`, and none `shifted`, at 4-8× reads (0.12 had J10 Pa4B `shifted`). No verdict settles it.
**Approach.** None in the code until the data can settle it.
**Check.** CENH3 ChIP or deeper long reads from those plants.

### 57. Per-plant verdicts that depend on other chromosomes

**Evidence.** Moving Pa1A's call in 0.13 changed verdicts on other chromosomes. J10 Pa4B went from `shifted` (loss q 0.015) to
`possible_shift` (q 0.32) with the same call site, and R1 Pa3B from `shared` to `unclear`. plants.md does not say that a verdict depends
on other chromosomes' calls; the cause (the multiple-testing pool, or the between-plant test) was not traced. Untested since 0.1 (plants.md,
Untested):
- the genotype noise model on two samples of one plant;
- coverage below about 2.5×;
- the per-plant Hi-C path, which has no other-plant library. Hi-C pairs cannot be genotyped, so another plant's library decides as `unknown`.

Pa2A's missing segments come out in 3 pieces each.
**Approach.** A leave-one-call-out check of each verdict's stability; a split of one sample's reads (or one plant sequenced twice) for
the noise model.
**Check.** Verdicts stable when an unrelated call moves; the split sample comes out as the reference plant.
**Code.** `plants.py`, `provenance.py`.

## Reporting and outputs

### 58. Report wording and undocumented outputs

**Evidence.**
- **Caveats.** They advise `tes/` wherever it is absent, also in mammals and fungi, where it cannot act (the siamang's confident misses;
  every *C. tropicalis* no-call; `interpret.py`, `_could`). `hic_near` has no sentence.
- **`hic_sequence_conflict`.** It warns that one line is wrong genome-wide, which reads as blame on the sequence line, but it fires in
  genomes whose sequence calls are right (horse T2T 29/32, pig, goat, sheep, rhesus, barley 6/7). A confidence cap on it was rejected for
  that reason (v8).
- **`hic_rabl`.** It is the Hi-C engine's telomere-pole test, which passed on non-Rabl rice and on nearly every T2T mammal, yet
  outputs.md lists it among the Rabl-axis tests.
- **Silent failures.** A failure in the report or the figure becomes one warning and the run ends with exit 0, the report silently
  missing that genome's section.
- **Undocumented outputs.** outputs.md documents the main tables column by column (all 105 columns of the calls table match), but not
  the loci, feature, plants, provenance and `te_signature` tables. It also omits 5 Hi-C engine flags (`no_rabl`, `null_untestable`,
  `not_converged`, `cycle_unstable`, `no_contacts`), two of which decide whether a prior-free peak is clean.
- **The HOR summary describes the repeat poorly.** The strongest class period equals the top tandem unit on 6 of 46 HG002 calls, and
  one changed class breaks a 16-monomer HOR into 5- and 4-monomer units (`sat/hor.py`, unchanged since 0.10).
- **No per-run check that homologous genomes agree** (homoeologs, haplotypes), though agreement is a benchmark proxy and cheap: in cenrep
  0.8.9, one goat's assemblies and sheep's haplotypes chose different main families.
- **Replaced calls.** The replaced call's own confidence is not stored, and `basis` calls a replaced long-unit array a "weak satellite
  array".

**Approach.** One table of flag → basis phrase, caveat phrase and settling input, read by `basis` and `caveat` and tested flag by flag.
Write failed sections into the report and `cenfuse_summary.tsv`. Generate the column tables from the code's column lists.
**Check.** A test per flag sentence; outputs.md complete, by a test.
**Code.** `interpret.py`, `output.py`, `sat/hor.py`.

## Speed, memory and reproducibility

### 5. A faster search

**Done in 0.11: the cost of low-complexity sequence.** A long model with a low-complexity tract spent 25-50 s per 2-Mb chunk without a
hit. Long models now search with those tracts neutral, 1.00-4.93x faster and never slower
([benchmark.md](benchmark.md#10-speed-memory-and-reproducibility); [design.md](design.md#engineering-notes), Search cost).
**Open: what that search does to results.** On 14 full reruns, none of the 143 chromosomes with truth changed outcome or confidence.
But only 8 of those genomes have truth, and the rest of the benchmark has not been searched again: its release reruns read stored
searches (#34). Without truth, overlapping families trade the lead in *Cuscuta* and *Panax*: 19 calls change main family and 2 change
place, one of them `high` -> `low`. A neutral model scores less, so a close family without a tract can win their overlaps. In cotton,
906 whole copies of a 512-bp family go to its 481-bp relative, which loses an rDNA tag. E-value parameters stay those of the whole model.
*Approach:* compare overlapping hits only on positions that can score in both models. *Check:* every benchmark arm rerun from the search
stage, with a split copy counted as lost; CENH3 in *Cuscuta* or *Panax* (#35).
**Open: the cost of real hits.** Chain families hit nearly everywhere (*Glycine dolichocarpa* gains 1.02x, *Cuscuta* 1.15x). Dense
short-unit arrays hold the rest (alpha satellite: 62% of CHM13's CPU), and in mouse 14 dispersed 1-2-kb models took 95% of a 2-Mb chunk
(v6). `search_cost.tsv` says which family costs what. Already tried, for the record:
- stopping a family's search early from a sample of chromosomes is v8's sample screen. On 216 stored genomes its second rule skipped 128
  families with no false skip, but real runs skipped none, for four reasons: fewer than 7 chromosomes cannot reach P ≤ 0.01; array-poor
  families cannot be measured (moss's fam009, 41% of moss's CPU); dense chaining families give too few breaks; and families interact,
  removing one shifting the others' arrays. Its code is saved (`v8/patches/v8-screen_screen.py`);
- tandem-copy models for short units are slow (4 copies: 43 s per 100 kb, against 0.1 s for 1);
- smaller pyhmmer blocks used 5x less memory and ran 3x faster on a hit-dense CEN180 job, but only about 40% of raw hits kept their
  position, since multimer models align at any phase of an array. design.md's "other block lengths: no gain" was measured on another
  workload;
- a BLASTN prefilter lost 12% of the hit bp; a minimap2 prefilter, proposed in 0.9, was never tested.

*Approach:* for dense short-unit arrays, a search that treats the seed arrays apart. For chain families, the screen's four failure
causes must be answered first. *Check:* monomer and array tables equal to the stored annotations, or every difference listed; monomer
counts against published ones (CEN180 on Col-CEN: 66,131, Naish 2021; cenrep found 67,241).
**Open: what a run costs.** The README's runtimes and memory predate 0.11, and GNU time's maximum counts only the largest process (the
seed pool's workers are not added). The runtime model behind the progress estimate (model length × genome size) fails on heavy
genomes: *Panax* 7.9 h predicted against more than 25 h actual, R² 0.72 over 104 genomes. *Check:* wall time, CPU and peak memory of a
fresh 0.13 run of 5-6 representative genomes, the pool included.
**Code.** `sat/annotate.py` (`search_cost.tsv`, `resolve_overlaps`), `sat/util.py` (progress).

### 6. Faster seed clustering

**Evidence.** BLAST takes 98% of seed clustering, driven by multimeric seeds: *Acanthus* took 76 min at 16 threads (the median over 106
genomes is 0.7 min; the next slowest, *Pontederia korsakowii*, 23 min). Only part of BLAST's output can ever be read by the clustering:
in *Acanthus*' heaviest round 7.0% of the HSPs (0.9% also pass the identity bar). Aligning only the pairs the clustering consults was 10×
faster with identical clusters on 8 genomes, but BLAST's hit lists differ for 0.2-0.5% of those pairs, so it is not exact by construction
(prototype `v9/bench/variant_blast.py`, `variant_cluster.py`). Single-thread query jobs on one database ran 1.6-1.9× faster per call,
with 3 of 873,000 hits different; `-mt_mode` gave identical hits in 12 of 13 reruns, with mixed times.
**Tested in 0.13, not built.** A fixed partition of the queries over single-thread BLAST processes, on *Acanthus*' stored rounds. The
heaviest round (200 x 200 seeds, 773,595 hits) ran slower (0.79x at 4 workers against BLAST's own 4 threads), and its hit lists were
not identical to one run's: with thousands of HSPs per pair, BLAST's output depends on how the queries are batched. Small rounds were
identical but gained little (0.65-1.39x). The script was lost (it lived in a session's scratch folder).
**Open: does the output depend on the thread count?** cenfuse runs `blastn -mt_mode 1` (queries split by thread) whenever queries
outnumber subjects. `sat/families.py` says the hits are identical and that "MAFFT --threadit 0 and BLAST give the same output at any
count", but the batching dependence above says otherwise. *Check:* `cluster_seeds` on *Acanthus*' stored seeds at 4 and 16 threads,
comparing the cluster assignments.
**Approach.** Fewer HSPs to report. A pair already known to match (or not) needs no more HSPs, so a two-pass search would cut the output
that dominates heavy rounds: first the best HSP per pair, then full tiling only for the pairs the clustering consults.
**Check.** Identical cluster assignments and representatives on every stored seed set.
**Code.** `sat/families.py` (`cluster_seeds`, `_blast`).

### 7. The FasTAN crash, reported upstream

**Evidence.** A 642.5-kb window of *Oryza minuta* (GWHESQK00000003:29,940,000-30,582,500) aborts FasTAN 0.8 every time (`malloc():
corrupted top size`, SIGABRT). The window has no N and no soft-masking, and about 90% of it is one satellite array (15-mer spacings of
243 and 342 bp); neither half crashes alone. *Sonchus oleraceus* crashed at another, uncharacterized site. Both sites lie inside
confident centromere calls. cenfuse recovers by splitting the piece: the 5-Mb seed chunk holding the window gives 757 rows after 4
splits, with nothing skipped. The report has not been filed.
**Approach.** Report it with the reproducer (`FAtoGDB`, then `FasTAN -m -p -T1` on that window), cenfuse's crash table as the second
example. Keep the splitting as a safety net. Filing is outward-facing, so the user decides. When a fixed release comes, the seed cache
must know the FasTAN version (#60).
**Check.** A fixed FasTAN release: no crash on the reproducer or on the 5-Mb chunk.

### 59. Checkpoints that miss changed inputs or code

**Evidence.** cenfuse has two checkpoint systems with different rules:
- **The satellite engine** keys its stages on hand-bumped stage versions and parameter values, not on its code or on its `genes/` and
  `tes/` inputs. A stale stage stops the run. Adding `tes/` or `genes/` to an earlier run needs `--force classify` by hand; otherwise the
  run resumes on stale satellite calls, the likeliest user error.
- **cenfuse's own stages** hash their source files and fingerprint their inputs. But they hash only the modules listed, not what those
  import (the 5S stage ignores the search code in `sat/annotate.py`), and whole files, comments included (a comment edited in `params.py`
  reruns Hi-C everywhere).
- **Pickled checkpoints** fail on a library upgrade with errors `_cached` does not catch. A moved or copied folder (`cp` without `-p`)
  changes every fingerprint. A corrupt checkpoint is recomputed without a log line.
- **The Hi-C matrix cache.** Deleting it drops the arm fold and alignment rules on resume, with one warning, because the hic checkpoint
  does not depend on it.
- **`--force STAGE`** recomputes later stages only inside the satellite engine, although README and outputs.md say every later one; and
  `--force all`, which the code supports, is rejected by the command line.
- **Concurrency.** There is no lock on the output folder, so two runs on one folder race. A seed-cache lock left by a job killed on
  another host makes later runs wait forever, with a single log line. `.done` stamps and the Hi-C engine's tables are written in place,
  not atomically.

**Approach.**
- One system: each stage keyed on its import closure's code (comments and docstrings stripped), its parameters, and its input
  fingerprints by content, not path and time, with `tes/` and `genes/` keyed for classify onward.
- Stale stages recomputed, or listed in one stop message.
- `--force all`.
- An output-folder lock (created exclusively, holding host and pid) and a stale-lock heartbeat.
- Atomic writes everywhere.

**Check.** Tests that a code or input change recomputes exactly the stages it should, and that a resume uses its checkpoints. Today's
resume test would pass if every checkpoint were ignored.
**Code.** `pipeline.py` (`code_key`, `_cached`), `sat/pipeline.py` (STAGE_VERSION, STAGE_PARAMS), `sat/seeds.py` (the lock).

### 60. Reproducibility across environments

**Evidence.** "Fresh runs are identical" holds within one Linux environment.
- **No record of the tools.** No output records the external tool or library versions (FasTAN, BLAST, MAFFT, cd-hit, Infernal,
  minimap2, TRF, numpy, pandas, pyhmmer), and the engine versions in the headers never move: 0.13's satellite outputs say "satellite
  engine 0.11.0".
- **Loose pins.** environment.yml pins some tools exactly but leaves numpy, pandas, scipy, Infernal, minimap2 and samtools open. The
  install uses `--no-deps`, there is no lock file, and the benchmark's versions are recorded nowhere.
- **macOS.** Without `setarch` about 1% of FasTAN's array cuts differ, with only a warning, and the seed-cache key records neither the
  `setarch` state nor the FasTAN version. A different family split moved pig's calls once (#41).
- **Thread counts.** Seed clustering is assumed to give the same hits at any thread count (#6).
- **Worker logs.** FasTAN crash lines are logged inside worker processes and are lost under `spawn` (macOS) or `forkserver` (Python
  3.14's Linux default), and one seeds test passes only through `fork`.
- **Old outputs.** Satellite outputs read with `--satellite-dir` from another version keep that version's rules, and nothing checks
  the version.

**Approach.** A run manifest (`run_info.json`: cenfuse version and git hash, command, host, threads, every tool's version string,
package versions, input fingerprints); a lock file per release; the FasTAN version and layout flag in the seed key; a version check on
`--satellite-dir`; worker log events returned to the parent process.
**Check.** A fresh run in a rebuilt environment equals the release; the macOS path gives a stated difference.
**Code.** `output.py`, `sat/seeds.py`, `sat/fastan.py`, `environment.yml`, `pipeline.py` (`stage_satellite`).

### 61. Memory

**Evidence.** The satellite search holds about one workspace per thread (up to ~2.4 GB for 2-kb units), and nothing caps it: `-t 32` on
a 64-GB node can be killed. Peak memory reaches 74 GB (rhesus at 12 threads), well above 12 workspaces, and the excess was never
explained. The CHG line holds about threads × 1-2 GB. The 3D feature-class test builds dense N × N matrices (1-2 GB per array at 15,000
bins) and is skipped above that size. Rejected in v5 (unrecorded in design.md):
- `malloc_trim` after each search: 229 GB;
- model-major job order: 125 GB;
- `MALLOC_ARENA_MAX=2` and freeing hit rows: no change.

**Approach.** A `--max-memory` resource option (a resource, not a threshold) that limits concurrent long-model searches, or the expected
peak logged before the search; a profile of the rhesus peak; the feature-class statistic from the sparse matrix.
**Check.** A run under a budget finishes with identical outputs.
**Code.** `sat/annotate.py`, `hicfuse.py` (`feature_contacts`).

### 62. Hi-C engine speed

**Evidence.** The Hi-C engine's own leave-one-partner-out test is cubic in the partners. The fast version, written for the anchored test
after a 25-partner genome took an hour, exists only in `hicfuse.py`. Null replicates, the cis simulations (200 maps per chromosome) and
the feature-class permutations run serially while the bootstrap is threaded. Each Hi-C library reruns the whole stage. Under load,
`scipy.sparse.load_npz` read the 32-MB members of the Hi-C caches at 20 s each. Hi-C takes minutes on the benchmark, so this matters for
birds, polyploids with 30-100 chromosomes and many libraries.
**Approach.** One leave-one-out and null implementation for both callers, mapped over the thread pool with per-replicate seeds.
**Check.** Identical tables, faster on chicken and on a 2n = 80 polyploid.
**Code.** `hic/confidence.py`, `hicfuse.py`.

## Code and docs

### 63. The integration cascade and other large functions

**Evidence.**
- **`fuse._decide`** is one 217-line ordered cascade (161 branches), to which each release added rules. Which rule wins is decided by
  its position in the cascade, and no output records which rule fired. Every rule proposed above (#14, #28, #32, #33, #45, #53) must be
  inserted into it.
- **Other large functions.** `interpret.genome_section`, `caveat` and `basis` need coordinated edits for every new flag. `Run.run` runs
  170 lines and builds the evidence twice. `hicfuse.py` holds four concerns in 877 lines.
- **Duplication.** Statistical helpers exist several times: binomial tails five ways, and the Hi-C confidence tests twice, nearly line
  for line. Small I/O helpers exist once per engine.
- **Coupling.** The chromosome-set rule lives in the methylation engine. cenfuse drives the Hi-C engine through its command line with a
  hand-built argparse Namespace. The LTR signature is computed twice per run, and `tes/` is parsed about six times.
- **Dead code**, such as `util.mcnemar_exact`.
- **The CLIs.** The four differ in exit codes, logging, the meaning of `-o` and the checks on `-t` (`-t 0` reaches a traceback in the
  satellite engine), and the engines' own CLIs are undocumented.

**Approach.** One function per rule, applied in a declared order and each returning its rule name (a `rule` column, with the order
documented in design.md); one `stats.py`; a Hi-C engine library call; the LTR line in its own package; one shared CLI helper.
**Check.** Byte-identical outputs on the benchmark before and after; `test_fuse_rules.py` as the per-rule spec.
**Code.** `fuse.py`, `interpret.py`, `pipeline.py`, `hicfuse.py`.

### 64. Constants without a reason in design.md

**Evidence.** design.md says every threshold is in the params files with its reason. In fact about 60 analysis constants live in 15
other modules (mostly with a reason in a comment), and the Constants table covers well under half of them. Missing are, among others:
- the Hi-C engine's strength tests, behind every "strong peak": N_BOOT 100, BOOT_IN_CORE 0.95, N_NULL 200, P_HIGH 0.01, MAPQ 30, END_WIN
  8, ICE_MIN_FRAC 0.2, and HUB_FOLD 2.0 and START_SMOOTH 5, both "heuristic";
- the CHG geometry: D_REF from maize, WIN_D and the ring, "heuristic" (#16);
- the satellite priors (#39);
- the LTR core resolution (TE_BINS 500, TE_MIN_BIN 1 kb, TE_SMOOTH 5), which on a 1.6-Mb chromosome gives 3-kb bins and 16-kb smoothing
  and splits small centromeres.

The code labels 19 constants "heuristic" and 8 "prior"; design.md names 9 heuristics. Bare literals also act as thresholds inside
functions: the 0.8 coverage rule appears four times in `sat/families.py`, and P 0.01 is hard-coded in `interpret.py`. LOCUS_FOLD = 2
copies the Hi-C hub fold (retrotransposon-derived arrays measured ≥ 2.5) and now decides both `dispersed` and `partial_chain`.
HALF_WIN's reason ("a 1-2 Mb CENH3 domain plus one bin") holds at 250-kb bins only, and the 250-kb standard bin floor was kept for
identity with 0.9.
**Approach.** Generate the Constants table from the code (each params module's reason comments and `sat/params.REASON`), with every
heuristic and prior marked. Add a test that fails on a module-level numeric constant without a reason, and name the bare literals.
**Check.** The table regenerated; no constant without a reason.
**Code.** `params.py`, `sat/params.py`, `hic/params.py`, `meth/params.py`.

### 65. Tests that pass without running

**Evidence.** The suite has 1,004 tests (collected), but those that need FasTAN, BLAST, MAFFT, cd-hit, minimap2 or TRF skip silently
without them. A machine without the bioconda tools therefore reports green while skipping the satellite engine's end-to-end tests. The
satellite tests run with test-scale values of 13 of the 51 parameters. These modules have no test of their own: `cli.py`, `linetests.py`,
`loci.py`, `methfuse.py`, `provenance.py`, `seqcalls.py`, `teline.py`, `util.py`. The checkpoint layer and the report's sentences are
covered only end to end, and no fixture holds IUPAC codes. The parallel family build is exact only while `frame_rotation` writes nothing
a later step reads, a premise guarded only by a docstring. The docs say the suite takes about 8 minutes; recorded runs took 137 s. There
is no CI.
**Approach.** A `CENFUSE_REQUIRE_TOOLS=1` switch that turns tool skips into failures (for CI and releases); unit tests for the checkpoint
layer, the summary row, each flag's sentences and the CLI parsing; one slow test at production values; CI from the lock file (#60).
**Check.** CI green with no skips.
**Code.** `tests/sat/sat_common.py`, `tests/fuse/test_fuse_e2e.py`.

### 66. Statements the 0.13 outputs no longer bear out

**Evidence.** To fix in the doc named; README.md and benchmark.md through their templates (#34).
- **design.md.**
  - The head-to-head test's cases (*B. rapa* A09, soybean chr1/chr6) no longer arise (entry 10).
  - "16-26 vs 14-32 elements" swaps the lower bounds: calls 14-26, cores 16-32.
  - The Known risks keep 0.12's Pa1A outcome; 0.13 moves it.
  - The gate text omits the majority of the whole reference set, which the code requires.
  - The pooled-LTR row gives only the first form's failure (#47).
  - The end-distance row carries the cis profiles' reason; its own was moss 23 → 25/26 but sheep 9 → 5, chicken 7 → 6, goat 12 → 11.
  - The not-goal "Hi-C choosing between satellite families" conflicts with `hic_chose_alt`, which can choose a rival acrocentric
    family's block. The intended meaning is "between the sequence line's two candidates on one chromosome".
  - It writes `gate_refs` where the output column is `hic_gate_refs`.
  - It says the 1,004 tests take "about 8 min".
- **benchmark.md.**
  - "459 of 2,355" disagrees with this file's 409 (entry 1).
  - Barley's fallback is still "3 confidently wrong" (0 since 0.12).
  - The annotation table is said to be from 0.11 (#48).
  - The sheep X dispute is called "unknown" without the karyotype check that supports the call.
  - The scoring section omits the 250-kb rule (#37).
  - Potato's dips lie on CENH3 on 8 of 11 chromosomes in the survey table and on 9 of 11 in the cis-trans text, with no stated
    difference between the two counts (entries 12 and 26).
  - Three figures date from 0.10.1, while the README says every number is from 0.13; the plants figure still shows tomato chr10 as
    "the failure".
- **README.**
  - "Generate the `tes/` file with exactly this command" includes `--proteins`, which no benchmark annotation used (#48).
  - Runtime and memory predate 0.11.
  - The workflow schematic still shows 0.12, the talk slides 0.10.1.
- **Code.**
  - Docstring numbers have drifted: `hicfuse.py` says mouse 20/20 (vs 19/19), lists soybean among the cis peaks, and gives tomato's P as
    6e-14 per sign; `fuse.py` says prior-free `high` calls are right on 45% (design.md: 14 of 63).
  - `sat/centromere.py` says "Decisive evidence (CRM) sets confidence; other clades only flag", but since 0.12 other lineages decide.
  - Six design.md section names cited in the code no longer exist.
  - `tests/fuse/test_align.py:165` cites "future.md 17", which was 0.10.1's entry (methylation at axis-pair calls) and is now the knob
    entry.
- **Spelling.** British spellings remain in user-visible messages (`meth/cli.py`: "analyse"; `meth/calibrate.py`) and in a name
  (`rdna_neighbourhood`).
- **The audit itself.** 0.12's docs were audited against its release outputs (20 stale statements found); 0.13's were not.

**Approach.** A release step that audits the docs against the release outputs, with docstrings citing benchmark sections instead of
numbers.
**Check.** The items above fixed; the audit written into the release PLAN.

## Leads for biology

Hypotheses that cenfuse's outputs could test. Each is tied to an entry, and none is a finding yet.

- **Young insertions mark where a centromere is now.** CRM copies are younger at Pa1A's Hi-C site than at its satellite (entry 31). The
  signature's youngest elements concentrate in the cores in nearly every CRM lineage and in Bryco, and the reverse holds in Athila,
  Retand and most Tekay (design.md). A general test: centromeres that moved within a species with CENH3 tables (soybean, tomato: #35)
  should carry young copies at the new site.
- **Which lineage CENH3 binds.** In tobacco, CENH3 binds Tekay and not the CRM clustered beside it; in tomato a single Tekay family marks
  the centromeres. Tekay may be the CENH3-bound lineage across Solanaceae (#46).
- **The Rabl configuration and condensin II.** The Rabl readouts pass in fungi, *Plasmodium* and barley, and in no territorial nucleus;
  condensin-II-deficient lineages keep the Rabl configuration (Hoencamp et al. 2021). Predictions: Diptera, which lack CAP-G2, pass;
  human cells without CAP-H2 gain the alignment (public Hi-C); declustering mutants (*S. pombe* csi1) lose it. G1 *C. neoformans* loses
  the cluster but keeps arm alignment on 8 of 14 chromosomes (#14).
- **What clusters in a mammal nucleus.** In CHM13 the off-target trans peaks sit on pericentromeric HSat2/3 and β-satellite rather than
  on the active array. Centromere clustering differs between mammals at equal depth: pig, horse, sheep and giraffe fail the test, goat
  and rhesus reach support. `hic_architecture` already labels about 50 genomes, with no compartment-matched null (#52).
- **The cis signal's sign.** It is a peak in tomato, rice and *C. neoformans*, a dip in potato and *Verticillium*, and absent in
  mammals. A fold of the flanks or an insulated domain would explain it, once a missing-length artifact is ruled out (#29).
- **Telomere repeats at bird centromeres.** Chicken's interstitial telomeric arrays lie at its CENP-A windows (P 4e-21), as described for
  its micro- and dot-chromosome cores (Huang et al. 2023). Other bird T2T assemblies would show whether this is general (#41).
- **5S rDNA beside HSA1q42 genes** in every primate, in horse ECA1 and in pig SSC14 (v6): a 5S locus that may be conserved across mammals.
- **Methylation and centromere position.** CHM13 clones whose methylation was changed moved CENP-A (PRJNA1270043). Whether the CpG
  kinetochore follows is a direct test of the CpG line's premise (#55).

## Tried and set aside (not in design.md)

Ideas measured in the development notes that design.md's "Choices tested and rejected" omits, or records without its reason or
numbers. Reopen one only with new evidence; most belong in that table.

| idea | result | where |
|---|---|---|
| `malloc_trim` after each search; model-major job order; `MALLOC_ARENA_MAX=2`; freeing hit rows | 229 GB; 125 GB; no change; no change (the fix was one pyhmmer pipeline per thread) | v4-v5 (#61) |
| smaller pyhmmer blocks | 5× less memory and 3× faster on a CEN180 job, but about 40% of raw hits keep their position | v5 (#5) |
| MAFFT `--thread 1` | reproducible, but 2.5× slower than `--threadit 0` | v5 |
| tandem-copy models for short units | 4 copies: 43 s per 100 kb (1 copy: 0.1 s) | cenrep (#5) |
| a sample screen to skip dispersed families | no family skipped on real genomes, for four reasons (entry 5) | v8 (#5) |
| BLAST `-mt_mode`; single-thread query jobs | identical in 12 of 13 reruns, mixed times; 1.6-1.9× per call, 3 of 873,000 hits differ | v9 (#6) |
| a minimap2 prefilter for long models | proposed, never tested (a BLASTN prefilter lost 12% of hit bp) | v9 |
| rebuilding a family on every length flag | 1,927 vs 1,755 s serial (33 of 40 *Panax* families flagged, 11 real) | v9 |
| recording the chromosome set after the engine runs | a silent reuse hole after a revert | v9 |
| an md5-only seed-cache key | a symlinked FASTA shared one entry between two sets | v9 |
| a pool of unclassified LTR elements | the guarded form placed *Cryptococcus* 10/14 and never lifted the cap, but qualified in 9 fish genomes; its 20-kb floor cost *O. punctata* a CRM call | v6 (#47) |
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
| cenhic 0.2's caps and re-located bootstrap | held-out H3b failed (`high` right 0 of 7); `telomere_dominated` capped horse's right acrocentric calls | cenhic |
| peaks restricted to B-compartment bins | gorilla 0 → 0.60, *P. annua* 0.93 → 0.36: the E1 sign at centromeres is opposite in plants and apes | cenhic (#14) |
| the best of many random starts; end windows excluded; a naive start | the telomere pole or a hub wins | cenhic |
| arm-scale segregation across the centromere | CHM13 6 of 23, plants at chance | v10 |
| the compartment eigenvector or the trans fraction as a locator | at or near chance; the trans fraction works on a minority, in mammals only | v10 |
| an effect-size test of nucleus type | separated holocentric and randomized libraries from yeast, not Type I nuclei from Type II | v10 (#52) |
| a confidence cap on `hic_sequence_conflict` | the flag fires where the sequence calls are right | v8 (#58) |
| the arms' depletion as a locator | tomato 0 of 12, potato 0/11, rice 0/12 per chromosome (a mean effect only) | v11 |
| the engine's peak levels required for core-anchored Hi-C | cleaned the controls, but left no real hits outside genomes that pass anyway | v11 |
| the Hi-C table swapped for the core-anchored one | potato 6 of 12, short of a majority | v11 |
| CUSUM on the GC deficit; the window by the largest trough gap | *M. sympodialis* 0 of 8; wrong windows in *C. neoformans* and *S. pombe* | v11 (#28) |
| lifting the honesty cap where one-per-chromosome passes | chicken 27 right, 5 wrong; pig 8 right, 9 wrong | v11 (#41) |
| interstitial telomeric arrays as a line | only chicken beats chance | v11 (#41) |
| a genome-level holocentric flag | also on 25 of 26 siamang and 25 of 32 EquCab3 calls | v4 (#44) |
| CpG: a cluster whose valley floor is half the call's | 6 changes in 7 primates: 1 right, 2 wrong | v4 (#38) |
| a CpG null from other satellite DNA; CpG bins as the mean of per-CpG fractions | infeasible (MFA8: 3.7 Mb outside the calls); MMU8 0 of 21 kinetochores at 1/16 depth | v4, v6 |
| a data-derived CHG domain size | rejected in v5, numbers not kept | v5 (#16) |
| 45S arrays by the Rfam SSU HMM; even spacing for 45S arrays | stalled on low-complexity sequence; only 100 of 272 arrays are regular | v3, v10.3 |
| 5S from the satellite families' Rfam tags | misses units over 2 kb (mammals, birds) | v6 |
| a MAPQ input check for Hi-C pairs | not built: the warning suffices (MAPQ ≥ 1 changed nothing) | v6 |
| sheep and goat CENP-A ChIP as truth | the antibody's epitope (phospho-Ser7) is absent in those species: at most 1.1-1.7× | v2 (#35) |
| loci split into dense runs; sparse clusters dropped | chimpanzee +5,757 loci; barley lost 108 of 156 subtelomeric loci | v8 |

## Closed entries

Kept for their numbers and for the evidence behind the rules. Each says what was done.

### 11. Hi-C anchors where the sequence line is wrong genome-wide

**Done in 0.12**, in two ways. Better calls make better anchors: with the LTR signature, barley's ordinary Hi-C test passes on 7 of 7
chromosomes and moss's on 26 of 26, with no new rule. And where the test still fails on most chromosomes, Hi-C is tested against the
LTR cores directly ([design.md](design.md#hi-c-and-the-ltr-line-test-each-other-hicfusepy)). Potato anchored on its cores finds 3 of 12
(P 0.02), which does not pass; its centromeres come from the cis signal instead (entry 12). Moving the windows to the best score near
each core found 6 of 12, but randomized Hi-C then passed too, so that was removed
([the benchmark, section 12](benchmark.md#12-how-the-benchmark-was-run)).

### 12. Cis signals of satellite-free centromeres

**Done in 0.12.** The contact crossing a bin, as observed over expected at each distance, has its extreme at CENH3 on all 12 tomato
chromosomes and 8 of 11 in potato. The LTR cores of the other chromosomes choose its sign and test it, with chance from simulated contact
maps. It moves `low` calls, and raises or lowers others by one level. See
[benchmark.md](benchmark.md#the-contact-crossing-the-centromere). Soybean, sorghum and the mammals have no such signal. Since 0.13 the
library's own trans contacts test it where there are no LTR cores (entry 26).

### 18. Organelle insertions called 45S arrays

**Done in 0.11.** Organelle rRNA holds the 18S motif too (rice chr4 9.23-9.37 Mb, reported as a NOR, is a chloroplast insertion). A
motif hit now counts only with nuclear context: the next base and the 19 after it. Hits followed by A fall in two groups: 110 with at
most 5 of 19 bases differing (nuclear) and 1,100 with at least 11 (organelle), 4 between. The base alone would have deleted real arrays
(barley 7H). Over the benchmark, arrays went from 272 to 261 (this entry had estimated 12 to go): rice chr4, *B. rapa* A02/A06, the
Sumatran orangutan. All 74 kept arrays sampled hold 5.8S genes, and removed ones hold at most 3. In the polyploids, 651 → 501.

### 19. Long-unit calls are mostly wrong

**Done in 0.11.** A long-unit array (units over 2 kb) with no support from Hi-C, `tes/` or methylation is reported as a feature, not as
the call (flag `long_unit_unsupported`; `seq_start`/`seq_end` keep it). Over all arms that removes 50 calls: 28 wrong, 8 right, 14 without
truth, all `low`. Hi-C's own call can then stand there: *M. sympodialis* LT671825.1 went from wrong to right. Costs: chicken 28/33 ->
27/33, *S. pombe* FASTA only 1/3 -> 0/3, one right call in sheep. Nothing yet tells a right long-unit array from a wrong one, so the 8
right calls go with the rest; another line can still keep them. Where centromeres are long-unit arrays of retrotransposon fragments
(*Coffea*, putatively), this leaves them uncalled from the FASTA.

### 21. Moved calls whose sentences disagree

**Done in 0.11.** `caveat` now reads the evidence `basis` reports. On moss's 4 moved calls it says the call rests on the 3D lines and
names the other evidence at that place, and it no longer says "no sequence or LTR evidence". No call changes.

### 23. A family's unit length can overstate its copies

**Done in 0.11.** Where they differ, `basis` gives the array's own spacing beside the family's unit ("165-bp unit, copies spaced 155 bp
here", column `seq_spacing`). In rice, 19 of the 28 large arrays of the 165-bp CentO family are not spaced 165 bp, 17 of them about 155
bp, because the 165-bp profile fits 155-bp copies with a 10-bp gap and scores higher. No call changes.

### 24. Telomeres outside the rules

**Done in 0.11.** Units of 13-29 bp are tried only where no 5-12-bp unit qualifies, and they must reach most ends: *C. tropicalis*' 23-bp
unit is found at 12 of 14 ends, with no spurious unit elsewhere. Once the unit is chosen, a stretch that reaches a sequence end counts
down to a density of 0.1 (61 of 6,336 ends lie under half density, none under 0.19). Ends went from 1,539 to 1,557 over the benchmark
and from 3,154 to 3,214 in the polyploids, none lost, interstitial arrays unchanged. Genomes where no unit is found at all are entry 56.

### 25. Non-CRM centromeric LTR retrotransposons

**Done in 0.12.** cenfuse learns each genome's signature at family level, leaving one chromosome out, and a lineage other than CRM
decides at 10-fold enrichment or more, unless a dominant satellite passes the one-per-chromosome test
([design.md](design.md#ltr-retrotransposons-tes-satsignaturepy)). Moss's Bryco, *Chlorella*'s Osser and tomato's Tekay family now
decide; *Arabidopsis*' Athila only flags, as CEN180 passes the test. What CRM's own standing rests on is entry 46.

### 26. The cis signal without `tes/`

**Evidence.** Tomato's cis peaks alone are right on 12 of 12 chromosomes, and potato's dips on 8 of 11, but in 0.12 the LTR cores had
to choose the sign and test it. Two ways to do without them were tried. One locus per chromosome: tomato passes (11 of 12 covered, P
0.003), while potato and *Verticillium* miss after correcting for two signs (P 0.044, 0.047), and chicken's dips pass (P 0.0008) with only
2 of 26 on CENP-A. Hi-C's own prior-free peaks: no agreement in tomato, potato or rice, and in *Ae. tauschii* 6 of 7 agree with 0 of 7 right.
**Done in 0.13**, with the library's own trans contacts as the test
([design.md](design.md#the-cis-and-trans-readouts-of-hi-c-test-each-other-hicfusecis_trans_line)). Each chromosome's cis extreme,
outside its terminal eighths, anchors the trans test of the LTR cores. Tomato: 10 of 12 chromosomes' trans peaks lie at their cis peaks
(P 1.2e-13, corrected for the two signs; 6e-14 for one), all on CENH3. *C. neoformans* in G1: 12 of 14. Read over whole chromosomes,
*Ae. tauschii* passed with 0 of 4 agreements on CENH3, all in terminal eighths. Randomized soybean Hi-C reached support (6/20), so the
line acts only at pass. Potato does not pass: its dips lie on CENH3 on 9 of 11 chromosomes here (8 of 11 in the survey of entry 12; the
two counts are not reconciled, #66), but its trans contacts at them do not meet. Chicken, the mammals and the holocentric genomes fail; no randomized library passes (3 reach support, where the line does
not act).
**What remains**: potato without `tes/` (#40 for its satellites), the calls in terminal eighths (#32), and the extremes the line trusts
on chromosomes where its own pair fails (#53).
