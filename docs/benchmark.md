# Benchmarks

Every number here comes from criteria fixed before the run they score: plants 2026-10-16 and -19 (original text:
[`docs/benchmark.md` at 961a23d](https://github.com/cwb14/cenfuse/blob/961a23d/docs/benchmark.md)), mammals
2026-10-21. Rules changed after a result was seen are listed in [section 5](#5-rules-changed-after-seeing-a-result).

**Scoring.**
- *Right*: mammals, the call overlaps the CENP-A domain; plants, the call lies within 4 Hi-C bins of the CENH3
  domain (the method's resolution).
- *Confident* = `high` or `medium`; `low` = cenfuse says it is unsure.
- *P*: Poisson-binomial tail of the right count. A chromosome's chance is (call + truth + tolerance) / chromosome
  length.

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
| pig CM000817.6, CM000822.6, CM000828.6 | wrong, high/high/medium | the truth is weak ChIP from other pigs, far from any satellite; on the first two, Hi-C also favors the satellite array |
| pig CM001155.3 | wrong, low | same |
| B. orangutan CM054652.2 | wrong, high | CENP-A sits 14 Mb from the only satellite array (likely an evolutionary new centromere) |
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
| horse (EquCab3) | 5 | 3 | 3 | 3 | **3** | untested | fail |
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
- **Horse EquCab3**: right calls are `low` (no reference-grade call, so the Hi-C test is untested); chr11 has no call.

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

## 5. Rules changed after seeing a result

Each change has a reason independent of the case that prompted it, and was rechecked on cases not used to make it.

| date | change | prompted by | rechecked on |
|---|---|---|---|
| 2026-10-16 | Hi-C peaks in end bins never oppose a call | rice chr4 | potato, barley, horse, gorilla, mammals |
| 2026-10-16 | Hi-C opposes only beyond 4 bins (its resolution) | CHM13 | genomes scored later |
| 2026-10-19 | karyotype: skip only windows > half satellite; compare same-type coverage | P. annua R2 | R1 Pa5B, the Hi-C library |
| 2026-10-23 | satellite rules (tandem test, satellite kinds, acrocentric test, competing clusters); Hi-C in mappable bins; NOR detection; test majority | the previous mammal run | plant benchmark, 185-genome batch |
| 2026-10-24 | arm-fold line; centers only at centromeric index ≥ 12.5% | barley | rice, all mammals, randomized controls |
