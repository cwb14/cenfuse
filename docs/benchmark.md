# Benchmarks

Every number here comes from criteria fixed before the run they score: plants 2026-10-16 and -19 (original text:
[`docs/benchmark.md` at 961a23d](https://github.com/cwb14/cenfuse/blob/961a23d/docs/benchmark.md)), mammals
2026-10-21, v4 2026-10-26 (CpG kinetochores, hard genomes, no-harm reruns: sections 6-8; the CpG criteria are
not blind, section 6), v5 2026-10-26 with dated amendments (memory, reproducibility, telomeres, new genomes:
section 9). Rules changed after a result was seen are listed in
[section 5](#5-rules-changed-after-seeing-a-result).

**Scoring.**
- *Right*: mammals, the call overlaps the CENP-A domain; plants, the call lies within 4 Hi-C bins of the CENH3
  domain (the method's resolution).
- *Confident* = `high` or `medium`; `low` = cenfuse says it is unsure.
- *P*: Poisson-binomial tail of the right count. A chromosome's chance is (call + truth + tolerance) / chromosome
  length.

## At a glance

![benchmark overview](img/benchmark_overview.png)

*[PDF with legends](img/benchmark_overview.pdf). **A** every chromosome with a mapped centromere: the call (box,
colored by outcome) against CENP-A/CENH3 (black bar); open boxes and arrows = barley calls moved by Hi-C. **B-F**
one chromosome per kind of outcome: right (B); uncertain, two similar satellite clusters (C); confidently wrong,
CENP-A on a small array (D); moved by the Hi-C arm fold (E, barley, with its contact map); no line has signal (F,
potato). **G** every chromosome as one square; black tick = chance. **H** each evidence line alone vs the fused
call. **I** the Hi-C test on real vs randomized Hi-C, and the arm-fold test.*

## 1. T2T mammals (14 genomes)

- **Truth**: CENP-A CDRs (Logsdon et al. 2024) for human CHM13 and five great apes. CENP-A ChIP for crab-eating
  macaque, pig and horse (pig and horse ChIP from other animals).
- **No functional truth** for siamang, rhesus, sheep, goat and giraffe. They are compared with the satellome locus,
  another sequence caller.
- **Inputs**: FASTA, + genes, + Hi-C (public libraries, 100 M pairs).
- **Controls**: the same Hi-C randomized two ways: ends moved, or partners re-paired with coverage kept.

**Overall (157 chromosomes with CENP-A truth):**

| inputs | right | confidently wrong | wrong, `low` | no call |
|---|---|---|---|---|
| FASTA only | **148** | 5 | 4 | 0 |
| FASTA + genes | 149 | 4 | 3 | 1 |
| FASTA + genes + Hi-C | **149** | 4 | 3 | 1 |
| *previous version (961a23d), FASTA only* | *112* | *19* | *26* | *0* |
| *previous version, all inputs* | *144* | *9* | *3* | *1* |

**Per genome (right / chromosomes; confidently wrong in brackets):**

| genome | truth | FASTA only | + genes + Hi-C | Hi-C test |
|---|---|---|---|---|
| human CHM13 | CENP-A CDR | 23/23 | 23/23 | pass |
| chimpanzee | CENP-A CDR | 18/18 | 18/18 | pass |
| bonobo | CENP-A CDR | 17/18 | 17/18 | pass |
| gorilla | CENP-A CDR | 9/9 | 9/9 | pass |
| Sumatran orangutan | CENP-A CDR | 8/8 | 8/8 | pass |
| Bornean orangutan | CENP-A CDR | 15/16 (1) | 15/16 (1) | pass |
| crab-eating macaque | CENP-A ChIP | 16/16 | 16/16 | pass |
| pig | CENP-A ChIP (other pigs) | 13/17 (4) | 13/17 (3) | fail |
| horse | CENP-A ChIP (other horse) | 29/32 | 30/32 | fail |
| rhesus | satellome locus | 22/22 | 22/22 | fail |
| goat | satellome locus | 30/30 | 30/30 | fail |
| giraffe | satellome locus | 15/15 | 15/15 | fail |
| sheep | satellome locus | 26/28 (1) | 26/28 (1) | fail |
| siamang | satellome locus | 18/26 | 20/26 (4) | pass |

- **Randomized Hi-C (28 runs):** the Hi-C test passed on no chromosome, and no call changed.
- **Arm-fold test:** failed in all 56 mammal Hi-C runs. Mammalian nuclei show no Rabl fold; the fold point hit
  CENP-A on 1/155 chromosomes.

**Pre-registered criteria:**

| criterion | result | numbers |
|---|---|---|
| B1 calls on the CENP-A domain (each ape > 50%, P ≤ 0.01) | PASS | 23/23, 18/18, 17/18, 9/9, 8/8, 15/16 |
| B2 final calls never worse than satellite calls | PASS | equal in all 9 genomes |
| B3 randomized Hi-C never passes the test | PASS | 0 of 28 runs |
| B4 the test passes only where Hi-C peaks lie within 4 bins of CENP-A | FAIL | in gorilla, orangutans and macaque the peaks sit beside the unmappable array, > 4 bins from CENP-A in bp |
| B5 `high` calls right more often than `medium`/`low` (P ≤ 0.01) | FAIL | 119/122 vs 30/34; only 3 `high` calls wrong, too few to test (P 0.25) |
| B6 Hi-C moves: right ≥ wrong | not tested | no call moved |
| B7 at least as good as the satellome caller on CENP-A | PASS | 90 vs 90 of 92 |

**Every miss (all inputs):**

| genome, chromosome | call | why |
|---|---|---|
| pig CM000817.6, CM000822.6, CM000828.6 | wrong, high/high/medium | the truth is weak ChIP from other pigs, with no centromeric-family array within 1 Mb; on the first two, Hi-C also favors the satellite array |
| pig CM001155.3 | wrong, low | same |
| B. orangutan CM054652.2 | wrong, high | CENP-A sits on 4 small alpha-satellite arrays (84 kb in all); cenfuse called the 3.4-Mb array cluster 14 Mb away. Hi-C points weakly to the right place |
| bonobo CM055489.2 | wrong, low | the true array is the reported competing cluster |
| horse CM128059.1 (chr11) | no call | satellite-free centromere; horse Hi-C failed its test, so it cannot place one |
| horse CM128070.1 | wrong, low | |

## 2. Plants and held-out genomes

Truth: CENH3 ChIP (rice: Shang et al. 2023; potato: Pham et al. 2020; barley), CENP-A ChIP (horse EquCab3),
satellome (gorilla). Inputs: FASTA + Hi-C. Each line scored on its own:

| genome | chromosomes | satellite line | prior-free Hi-C | anchored Hi-C | **final** | Hi-C test | fold test |
|---|---|---|---|---|---|---|---|
| rice (T2T-NIP) | 12 | 12 | 0 | 7 | **12** | pass | fail |
| barley (MorexV3) | 7 | 0 | 6 | 0 | **6** | fail | **pass** |
| potato (DM v6.1) | 11 | 2 | 1 | 0 | **2** | fail | fail |
| horse (EquCab3) | 5 | 3 | 3 | 3 | **3** | fail (on its `low` calls) | fail |
| gorilla (mGorGor1) | 25 | 25 | 0 | 19 | **25** | pass | fail |

- **Arm fold point vs CENH3** (exploration, not a criterion): rice 11/12 (500-kb bins), barley 4/7 within 8 Mb
  (its core holds CENH3 on 6/7), potato 0/11.
- **Rice**: the satellite line is right everywhere. Hi-C validates it, so all 12 calls are `high`.
- **Barley**: the satellite calls sit at the chromosome ends (CENH3 is on Cereba retrotransposons). Hi-C's
  prior-free peaks are right, but they can't be checked against wrong satellite calls.
  - The arm fold agrees with those peaks on 6/7 chromosomes (P 2e-5), so the two move 6 calls; all 6 are right.
  - The seventh fell from `high` to `low`.
  - Previous version: 0/7, 3 of them `high`.
- **Potato**: half its centromeres are satellite-free, and its Hi-C carries no centromere signal. Nothing can rescue
  them without `tes/`: 2 calls are confidently wrong, 7 are `low`.
- **Horse EquCab3**: right calls are `low` (no reference-grade call); the Hi-C test, run on those `low` calls, fails
  (4/24 agree), so Hi-C changes nothing; chr11 has no call.

**Poa annua / P. infirma (all four inputs; truth from earlier analyses, no ChIP):**

| criterion | result | numbers |
|---|---|---|
| P1 Pa1A on its functional centromere (167-170 Mb, 21 Mb from its satellite) | PASS | 163.0-169.4 Mb on the CRM core; Hi-C peak 168.0 Mb; `high` |
| P2 the other 20 calls unchanged | PASS | 20/20 |
| P3 Hi-C test passes | PASS | annua 14/14, infirma 7/7 |
| P4 methylation supports ≥ 10/13 annua calls | FAIL | 0/13 now. All 5 methylation plants are other plants, so methylation does not decide (12/13 before the provenance rule) |
| P5 per-plant candidates at Pa1A listed | FAIL | the J10 candidate lies inside the wide call |
| N1 permuted Hi-C: test fails, no change | PASS | 0/14, no change |
| demonstration: start from satellite-only calls | moved | Pa1A moved 20.7 Mb onto the functional centromere (Hi-C + CRM core), `high` |

**Other plant criteria:**

| criterion | result | numbers |
|---|---|---|
| H1 final ≥ satellite line (rice, potato, barley, horse) | PASS | 12 ≥ 12, 2 ≥ 2, 6 ≥ 0, 3 ≥ 3 |
| H2 changed calls: right ≥ wrong | PASS | 6 right, 0 wrong (barley) |
| M1 six imported T2T mammals: final ≥ satellite | PASS | 6/6 |
| M2 anchored Hi-C beats prior-free Hi-C (6 mammals pooled) | PASS | 99 vs 73 of 139, P 9e-7 (73 vs 73 in the previous version) |

**Several plants (P. annua reference + 5 other ONT plants; synthetic genome with planted plants):**

| criterion | result |
|---|---|
| S1-S4 synthetic: provenance, reference calls, planted shift, bad inputs | PASS |
| R1 all five plants recognized as other plants | PASS |
| R2 the five known missing segments found (within 5 Mb) | FAIL as written: all found, but the two Pa2A losses come out in 3 pieces |
| R3 the Hi-C library is `unknown` and decides | PASS |
| R4 reference calls unchanged | PASS |
| R5 no false shifts; P3 and R1 `call_absent` at Pa2A | PASS |

On 300 pure-noise genomes, the between-plant test gave 5.0% of null tests at P < 0.05 and no `shifted` verdict.

## 3. Satellite batch (185 plant and fish genomes, satellite line only)

This is a regression check of the satellite rules, not an accuracy test: only two complexes have CENH3 truth.

| change since the previous version | value |
|---|---|
| calls moved / confidence-only changes | 44 / 51 of 2,140 |
| Ae. geniculata (CENH3) | 10/14 overlapping, unchanged |
| S. verrucosum (CENH3) | 3/12 overlapping (7/11 within 1 Mb), unchanged |
| F. chiloensis D | 6 more chromosomes on the CENH3-validated Type I satellite |

## 4. What the controls show

| control | expected | observed |
|---|---|---|
| randomized Hi-C, ends moved (14 mammals, P. annua) | Hi-C test fails, no change | 0 passes, no change |
| randomized Hi-C, coverage kept (14 mammals) | same | 0 passes, no change |
| arm-fold test on mammals, potato, rice | fails where there is no Rabl fold, or where Hi-C's own peaks are wrong (rice) | fails in all |
| fold test, reverse conflicts (S. orangutan, giraffe: satellites right, prior-free Hi-C wrong) | must not side with Hi-C | 0/25, 0/16 agreement, no change |
| CpG values rotated along each chromosome (8 primates; values from half a chromosome away, on real CpGs) | CpG gate fails | 0/8 genomes |
| CpG positions shifted by half a chromosome (8 primates) | rejected or gate fails | 6/8 rejected as another assembly (rows off CpGs); 0/2 gate passes |

## 5. Rules changed after seeing a result

Each change has a reason independent of the case that prompted it, and was rechecked on cases not used to make it.

| date | change | prompted by | rechecked on |
|---|---|---|---|
| 2026-10-16 | Hi-C peaks in end bins never oppose a call | rice chr4 | potato, barley, horse, gorilla, mammals |
| 2026-10-16 | Hi-C opposes only beyond 4 bins (its resolution) | CHM13 | genomes scored later |
| 2026-10-19 | karyotype: skip only windows > half satellite; compare same-type coverage | P. annua R2 | R1 Pa5B, the Hi-C library |
| 2026-10-23 | satellite rules (tandem test, satellite kinds, acrocentric test, competing clusters); Hi-C in mappable bins; NOR detection; test majority | the previous mammal run | plant benchmark, 185-genome batch |
| 2026-10-24 | arm-fold line; centers only at centromeric index ≥ 12.5% | barley | rice, all mammals, randomized controls |
| 2026-10-26 | CpG line (new): within-call valley, shuffle test, positional gate | designed on the 7 primate tracks it is scored on (not blind) | MMU8 (R10), Arabidopsis, rotated and shifted controls |
| 2026-10-26 | CpG gate: positional (rotation) null instead of dips in other satellite clusters | the rotated control passed in gorilla | the other 7 rotated controls, real tracks |
| 2026-10-26 | telomeres: ≥ 100 bp counts when the array reaches the sequence end | *Neurospora* 0/14, *P. falciparum* 25/28 | 14 T2T mammals and mouse (no end lost or moved), 28 genomes |
| 2026-10-26 | Hi-C gate: `low` calls are references when no call is `high`/`medium` (clean peaks only) | siamang without genes (gate untestable) | 4 randomized Hi-C runs, chicken, EquCab3 |
| 2026-10-27 | acrocentric test: interior rival judged per satellite kind; a second unrelated acrocentric satellite competes with the call | T2T mouse (0/20; the v3 prediction failed) | the 72 earlier satellite runs (goat, the other all-acrocentric genome, among them) and 9 other v5 genomes: no change |

## 6. Kinetochores from CpG methylation (v4)

Criteria: `v4/CRITERIA.md` (2026-10-26). **Not blind**: the line was designed on these same primate tracks
(prototypes and rejected rules in `design.md`). Held out from the design: MMU8 and Arabidopsis.

- **Inputs**: public ONT CpG tracks of the assembled individuals. Apes: GenomeArk modkit 5mC-CpG bigWigs (R9, Guppy
  6.3.8; fraction only, as bedGraph). CHM13: nanopolish bigWig. MFA8: nanopolish with read counts (as bedMethyl). MMU8:
  R10 Dorado SUP + modkit bedMethyl. With FASTA + genes + Hi-C, as the mammal FULL arm.
- **Truth**: D = DiMeLo-seq CENP-A (CENP-A-directed 6mA, not CpG calls; 5 apes): domain = half-maximum run around the
  top 5-kb bin of each centromeric region. F = curated CDRs (same ONT data as the input: an implementation check). W =
  MFA8 CENP-A ChIP (testis, another animal; Mb-wide).

| criterion | result | numbers |
|---|---|---|
| C1 gate passes on real tracks (≥ 5/7) | PASS | 6/7: CHM13, bonobo, gorilla, both orangutans, MFA8. Fails: chimpanzee (positional dip in 1/25 references); also siamang (no truth) |
| C2 gate fails on control tracks | PASS | rotated 0/8; shifted: 6/8 rejected by the assembly check, 0/2 pass |
| C3 kinetochore on DiMeLo CENP-A | PASS | within 50 kb 35/39 (overlap 29/39); DiMeLo higher inside than over the call 39/39 (median 5.02 vs 0.81) |
| C4 kinetochore on the curated CDR | PASS | 63/63 |
| C5 MFA8 kinetochore inside the ChIP domain (≥ 90%) | FAIL | 12/15. The 3 others sit in the same array, 1-3 Mb from the domain of another animal's testis ChIP |
| C6 no call changed | PASS | 0 differences from v3 in calls, confidence, relation, support, flags (8 genomes) |
| C7 resolution | PASS | median kinetochore 60 kb vs call 5.1 Mb (85×); 132 kinetochores |

- **Positional dips** (share of tested calls; the gate counts them): CHM13 17/23, B. orangutan 19/25, S. orangutan
  16/25, bonobo 11/24, gorilla 7/25, MFA8 20/21; chimpanzee 1/25 and siamang 4/22 (gate fails). Rotated tracks: 1-2
  per genome (4-9%).
- **MMU8** (held out; R10 with counts; no functional truth): dips in 21/21 calls, positional in 20/21, gate passes,
  21 kinetochores (methylation 0.19-0.34 inside vs 0.52-0.61 over the array).
- **Wrong calls get kinetochores too**: B. orangutan CM054652.2 (the confident miss) gets one inside its inactive
  array (dip 0.15 vs 0.32). The line refines a call; it cannot check it. Tested and not adopted: CpG choosing between
  arrays (design.md).
- **Runtime** (final code): the CpG stage takes 1 min (CHM13 bedGraph) to 3 min (MMU8, 56 M per-strand
  bedMethyl rows), one thread.

## 7. Hard genomes (v4)

### Arabidopsis thaliana Col-CC (held out; criteria fixed before its data were processed)

Inputs: Col-CC (5 chromosomes); same-genotype ONT all-context modBAM (Dorado sup `5mC_5hmC`, ~56×) piled up with
modkit (CG combined, CHG and CHH per strand); `tes/` from LTRquest; no genes, no Hi-C. Truth: CENH3 ChIP of Col-0, by
the held-out plant rule (100-kb bins, half-maximum run, ≥ 4×): 1.3-1.6-Mb domains, 17-19×.

| criterion | result | numbers |
|---|---|---|
| A1 sequence calls on CENH3 (reported) | 5/5 | all `high`, on CEN180 (2.8-Mb calls around 1.4-Mb CENH3 domains); with `tes/` the Athila core supports all 5 |
| A2 CHG line does no harm | PASS | calls identical with and without methylation |
| A3 CHG line supports ≥ 3 right calls | FAIL | 0/5: no methylation signature (leave-one-out signatures on 1/5 chromosomes) |
| A4 CpG line in a plant | PASS | gate fails (positional dips 2/5); no kinetochore reported |
| A5 CpG line does no harm | PASS | CG-only calls identical to FASTA-only |

- **CHG line (first same-plant run with CENH3 truth):** inert. Its expected domain size comes from genome size (D =
  315 kb here), but CEN180 arrays span 2.2-3.5 Mb, so the calibration ring (1.5-4 D from the anchor) lies inside the
  centromeric arrays. A data-derived D is the likely fix; it was not adopted, because this genome prompted it and no
  second same-plant ChIP case exists to recheck it.
- **CpG:** the gate was right to fail. Exploratory: the strongest dips are shallow (0.62-0.74 vs 0.78-0.83 over
  CEN180) and lie off CENH3 on 4/5 chromosomes, so Arabidopsis CENH3 is not marked by a CpG dip in these data.

### Chicken GGswu (FASTA only; birds: 40 chromosomes, 2.5-200 Mb)

Truth: the authors' CENP-A windows on this assembly (Huang et al. 2023; DT40 ChIP of another bird; plotting windows of
0.1-2.6 Mb, 33 chromosomes). Chicken centromeres use chromosome-specific satellites, and chr5, chr27 and Z have none.

| result | numbers |
|---|---|
| calls on the CENP-A window | 28/33 (chance 6.3, P 6e-23); 2 more within 1 Mb (chr32, chr34, beside 0.3-0.8-Mb windows) |
| confidence | all 40 `low`: no satellite family dominates, so the satellite engine caps every call |
| misses (all `low`) | chr2 (a 1.4-kb-unit array 90 Mb from CENP-A), chr3 (3.4 Mb away), chr5 (satellite-free; called near its end) |

All 40 chromosomes, including the 2.5-7-Mb microchromosomes, are analyzed: the chromosome-set rule's only tenfold
drop is at the mitochondrion. The calls are right far more often than `low` suggests (28/33), so confidence is
conservative here. No confident call is wrong.

### Holocentric and CenH3-less genomes (FASTA only)

Any confident single centromere is wrong here.

| genome | why it is hard | calls |
|---|---|---|
| *C. elegans* CGC1 (T2T, 6 chr) | holocentric, no centromeric satellite | 6 `low`, all `no_dominant_satellite` |
| *Bombyx mori* P50T (29 chr) | holocentric; Lepidoptera lost CenH3 | 14 `low`, 15 none |
| *Rhynchospora breviuscula* (5 chr) | holocentric; Tyba arrays along every chromosome | 5 `low` |

0 confident calls in 40 chromosomes, with no code change. A genome-level "holocentric" flag was considered and not
added: `no_dominant_satellite` also marks most calls of monocentric siamang (FASTA only, 25/26) and horse EquCab3
(25/32).

### Telomeres (14 T2T mammals, v3 runs)

First scoring: cenfuse finds a telomere array at 687 of the 688 ends that hold one by the v2 truth rule (≥ 500 bp of
TTAGGG runs in the terminal 10 kb). The miss is sheep CP162267.1 (end). It also reports arrays at 4 ends the rule does
not count (it searches the terminal 100 kb).

## 8. No harm to the earlier benchmarks (v4 code)

Every earlier run was repeated with the v4 code (checkpoints resumed): the 14 T2T mammals with genes and Hi-C and the
16 plant-benchmark runs (rice, potato, barley, horse EquCab3, gorilla, the six imported mammals, P. annua with all
inputs and methylation, P. annua plants, the negative control, the Pa1A demonstration). Calls, loci and summaries are
identical in 30/30 runs; only the new `cpg_*` columns (empty) were added.

## 9. v5: memory, reproducibility, telomeres, TE formats, new lineages

Criteria: `v5/PLAN.md` (2026-10-26, before any v5 edit; dated amendments for tracks added later).

### Memory of the satellite search

Every pyhmmer `nhmmer()` call built a search pipeline that ended in a reference cycle. Its workspace (model length ×
the 256-kb search block: ~1.7 GB for an 1,853-bp model) waited for Python's full garbage collection, which a large
process runs rarely, so memory grew with every search done. One pipeline per thread, reused (as pyhmmer's own workers
do), gives the same hits and holds one workspace per thread.

| genome (threads) | before | after | outputs |
|---|---|---|---|
| chicken GGswu (24) | 277 GB | 58 GB | monomers, arrays, calls, loci identical |
| giraffe T2T (12) | 140 GB | 17 GB | monomers and arrays identical |
| *Arabidopsis* Col-CC (8, run side by side) | 7.8 GB, 23:47 | 4.8 GB, 23:59 | monomers identical |

Tried and not used (dead ends): glibc `malloc_trim` after each search (chicken 229 GB), job order (125 GB), smaller
pyhmmer blocks (hits change).

### Reproducibility of the satellite search

The families stage gave a different model for one *Arabidopsis* family on every run (3 runs, 3 models; one monomer
more or less downstream). MAFFT's multithreaded iterative refinement was the cause (4 runs, 4 alignments). With the
refinement single-threaded (`--threadit 0`), the alignment is the same at 1-32 threads at no time cost. Two
*Arabidopsis* runs from scratch at 16 and 8 threads now give byte-identical outputs (families, monomers, arrays, calls,
loci, summary).

### Telomeres (28 genomes: mammals, a bird, an insect, a nematode, plants, fungi, a protist)

A short array now counts when it reaches the sequence end: *Neurospora* NC12 0 → 13/14 ends (157-207-bp telomeres),
*P. falciparum* 25 → 27/28, horse EquCab3 19 → 27, yeast 0 → 7 (its irregular TG1-3 telomeres, as GGGTGT). The 14 T2T
mammals and mouse are unchanged (no end lost or moved); every added array lies within 98 bp of its end. Remaining
misses: ends without terminal repeats in the assembly (chicken 22), TTAGG broken by TRAS/SART retrotransposons
(*Bombyx* 15), fragmented NOR-arm telomeres (*Arabidopsis* 2).

**Interstitial telomeric arrays (new loci kind).** CHM13: 16 arrays ≥ 300 bp, most within 1 Mb of an end, one at the
2q13 fusion site (chr2 114.03 Mb). Chicken: 50 arrays, 22 of them on or within 100 kb of a CENP-A window (expected 1.4:
the padded windows cover 2.9% of the genome; binomial P 4e-21), holding 63% of all interstitial TTAGGG (up to 210 kb on
one chromosome): the telomere-repeat centromere cores Huang et al. 2023 described, genome-wide.

**NORs.** Every 45S array call that can be checked lies on a known NOR chromosome: *Arabidopsis* chr2/chr4, chicken
chr16, *C. elegans* chrI, pig chr8/chr10, human chr13/14/15/21/22, horse chr1/28/31; barley 1H, 2H, 5H, 6H, 7H (major
5H and 6H plus minor sites); mouse chr12, 15, 16, 18 and 19 (all five NOR chromosomes of the laboratory mouse).

### LTR retrotransposon input

DANTE_LTR GFF3 is read as it is (lineage from `Final_Classification`). *Arabidopsis* with the DANTE_LTR file: the same
5 calls, all `high`, as with LTRquest. DANTE classifies more lineages, so the flags-only marker clade can differ (Tekay
here, cored on 3/5 chromosomes; TE support 2/5 vs 5/5 with LTRquest's Athila). Ranking the marker by chromosomes first
was tested and rejected: over the 185-genome batch it agreed less with the independent satellite calls (te_discordant
111 → 129).

### Hi-C where every sequence call is `low`

When no satellite family dominates, every sequence call is capped `low`, so the Hi-C gate had no reference-grade call
and Hi-C could never act. The low calls are now the references, and only agreements at clean peaks (not at chromosome
ends, hubs or NORs) count.

| genome (inputs) | gate | effect |
|---|---|---|
| siamang (FASTA + Hi-C, no genes) | 0 → 26/26 chromosomes | 20 calls `low` → `medium`, one Hi-C-placed call; against the genes + Hi-C run, the `medium` calls match 18/20 and the remaining `low` ones 4/6 (before: 25 `low`, 21 matching) |
| siamang, randomized Hi-C (both modes) | 0/26, 0/26 | none |
| chicken (FASTA + 40 M Hi-C pairs) | fails on 30/30 | none; 25/30 anchored peaks sit in end windows: chicken Hi-C is dominated by telomere clustering (telomere arrays 4.0×, centromere calls 1.4×) |
| chicken, randomized Hi-C (both modes) | fails on 30/30 | none |
| horse EquCab3 (FASTA + Hi-C) | fails (4/24 agree) | none |

### New genomes

**Small genomes with point, regional or AT-rich centromeres** (FASTA only, seconds each): *S. cerevisiae* R64,
*S. pombe*, *P. falciparum* 3D7, *N. crassa* NC12. No confident call: yeast 0/16 calls (point centromeres), *S. pombe*
1 `low` call on cen3 (a long-unit tandem), *P. falciparum* 1 `low` call on a chromosome 4 var-gene cluster (the gene veto
removes it when genes are given), *Neurospora* 0/7 (RIP-mutated TE-relic centromeres). A genome with an IUPAC code every
5 kb (2,620 in all, inside an inserted 100-kb satellite array too) runs cleanly and finds the array exactly.

**Maize Mo17 T2T** (FASTA only; truth: same-line CENH3 ChIP, PRJNA751841, held-out plant rule). The 180-bp knob repeat
is the largest satellite (23.6 Mb, 20.6 Mb of it on chr8) and is classed `interstitial`; CentC is centromeric on all
10 chromosomes, all calls `high`, none holding any knob repeat. Calls on CENH3 (38-49× domains): 7/10 overlapping,
9/10 within 1 Mb (4 Hi-C bins, the plant rule; chance 0.2, P 1e-14). The misses (chr2 0.3 Mb, chr8 0.7 Mb, chr5 2.0 Mb
away; chr5 confidently wrong) have CENH3 on satellite-free DNA beside the CentC arrays, the case `tes/` (CRM cores) is
for.

**Mouse T2T (mhaESC; all 20 chromosomes acrocentric)** (FASTA only; truth: the minor-satellite array at the start of
each chromosome, 0.48-1.04 Mb, the CENP-A satellite). This tested the prediction that all-acrocentric karyotypes are
called without genes, and the prediction failed: 0/20 calls on the minor satellite, 14 confident calls 1.7-72 Mb
away (9 `high`), 2 chromosomes without a call. The major satellite (162 Mb, the largest satellite, pericentric) and
the minor satellite (16.5 Mb) were both judged subtelomeric. A 58-bp major-satellite subfamily at the inner edge of
the terminal block (beyond the 5% end zone on 10/20 chromosomes) counted as a strong interior rival and blocked the
acrocentric test. Two changes:
- the interior rival is judged per satellite kind (subfamilies pooled, as dominance already is), so the subfamily is
  part of its terminal block;
- both satellites then pass the acrocentric test. They are unrelated by sequence, and which one carries the
  kinetochore is not in the sequence: here the larger one is pericentric. In goat, the two acrocentric satellites are
  related layers and are pooled into one call. A second acrocentric satellite unrelated to the call now competes
  with it: the call is `low`, and that satellite's terminal block is reported (`seq_alt_*`).

Now all 20 calls sit at the centromeric end, on the major-satellite block, each within 45 kb of the minor satellite
(1/20 overlapping it), all `low`; the minor-satellite array is the reported alternative on 20/20. Of the criteria set
before the run, ≥ 15/20 calls overlapping the minor satellite still fails; no confident call off it now passes. Hi-C cannot
choose between the two arrays, because no mappable bins separate them. Telomeres are found at all 40 ends, plus 36
interstitial telomeric arrays.

### No harm (v5 code)

Every earlier run was repeated with the v5 code (checkpoints resumed): the 14 T2T mammals with genes and Hi-C and the
16 plant-benchmark runs. In 29/30 the calls, confidence, relation, support, against, flags and summaries are
identical; new columns (`hic_gate_refs`, `telomere_interstitial`) are added and the loci gain only interstitial
telomeric loci (+1 to +38 per genome; every earlier locus kept, or grown by one). Horse EquCab3: 8 short end telomeres
are found (19 → 27 ends), and its Hi-C gate, untested before (no reference-grade call), now tests the low calls and
fails; no call changed. The satellite search on the same families gives byte-identical monomers, arrays and calls
(chicken, *Arabidopsis*). The two acrocentric changes (mouse) act inside the satellite engine, so its classify and
call stages were replayed on all 72 earlier satellite runs (the 14 T2T mammals with and without genes, 6 plant runs,
38 polyploid complexes of the 185-genome batch: 3,134 chromosomes) and on the other 9 genomes of this round (104
chromosomes): no call, confidence or flag changed.
