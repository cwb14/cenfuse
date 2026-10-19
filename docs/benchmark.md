# cenfuse minimal benchmark: genomes and criteria

Written 2026-10-16, before cenfuse exists and before any cenfuse run. Nothing below may be changed after a run
except to add results. Budget: at most 6 h wall clock from the first benchmark launch, and at most 120 threads at
once. A run unfinished at 6 h is stopped and reported as postponed.

For publication, local file paths, file names and run-folder names were removed from this document (for example,
"the memo's bounds" now reads "the earlier analysis's bounds"); nothing else in the pre-registered text was changed. The run records (commands, scorers, per-chromosome scores) are kept outside the
repository.

Development status of the data. Every genome below was seen by at least one engine (the engines' own validation).
These are tests of the integration, not of the engines:
- P. annua and P. infirma: cenrep, cenhic and cenmeth development data.
- The 8 T2T mammals: cenhic development genomes.
- rice and potato: cenhic held-out genomes. Their cenhic calls are known: rice 0/12, potato 1/10.
- Anchored, leave-one-chromosome-out Hi-C with gate and fusion has never been run on any of them. That is what is
  tested.

## Genomes

| id | assembly | lines | truth | role |
|---|---|---|---|---|
| poa | P. annua + P. infirma reference assemblies | satellite, TE (LTRquest), genes, Hi-C; methylation for annua (5 ONT plants) | no ChIP. Earlier-analysis proxy: Pa1A functional at 167-170 Mb, the other 20 at their satellite | integration test, full stack |
| negctl | P. annua with Hi-C pairs permuted within chromosomes | as poa, minus infirma and methylation | Hi-C carries no centromere signal | negative control |
| rice | T2T-NIP GCA_034140825.1, 12 chromosomes | satellite, Hi-C (SRR25241089) | CENH3, Shang 2023 Table S1 | functional truth; Hi-C de novo fails (all ends) |
| potato | DM v6.1, chr01-chr12 | satellite, Hi-C (SRR11908553) | CENH3, Pham 2020 Table S5, 11 chromosomes | functional truth; half the centromeres satellite-free; Hi-C de novo fails (hotspots) |
| mammals | chm13, mPanTro3, mPonPyg2, mPonAbe1, mSymSyn1, mGirTip1 | satellite (cenrep 0.8.9 outputs, `--satellite-dir`), Hi-C (public libraries of each assembly) | CHM13: CENP-A active HOR (censat; independent). Others: satellome locus (sequence-based, not independent of the satellite line) | does anchoring fix the wrong-cluster failure? (de novo: 97/188 overall; mPonAbe1 4/25, mGirTip1 0/15) |
| stretch | barley MorexV3; gorilla mGorGor1 (cenrep engine run needed); horse EquCab3.0 | as rice; gorilla as mammals | ChIP truth (barley, horse); satellome (gorilla) | only if budget remains; else postponed |

Inputs are the files of the engines' earlier benchmarks, symlinked, never copied or modified.
Held-out FASTAs are reduced to chromosome-level sequences, as in the earlier benchmarks.

## Scoring

- **Tolerance:** 4 bins at the Hi-C bin size cenhic picks for that genome, as in cenhic's benchmark criteria
  (cenhic's partner half-window, the resolution the method claims).
- **Right:** the gap between the call interval and the truth interval is at most the tolerance (0 when they
  overlap).
  - For a Hi-C peak, the interval is the peak bin.
  - Each line is scored on its own: sequence call, de novo Hi-C peak, anchored Hi-C peak, integrated call.
  - The distance from the call's midpoint to the truth is reported too.
- **Chance per chromosome:** (call length + truth length + 2 x tolerance) / chromosome length, capped at 1. P is the
  Poisson-binomial tail P(K >= hits).
- **Change:** a chromosome whose integrated call does not agree with its sequence call (moved or functional_only).

## Criteria

poa (implementation; the answers come from earlier analyses of the same data):
- P1: the integrated Pa1A call overlaps 167.0-170.0 Mb.
- P2: on the other 20 chromosomes, the integrated call agrees with the sequence call on >= 19.
- P3: the Hi-C gate passes on >= 12/14 annua and >= 6/7 infirma chromosomes.
- P4: methylation (pooled, leave-one-out signature) is `centromere_like` at >= 10 of the 13 non-Pa1A annua
  integrated calls.
- P5: on Pa1A, the sample table holds a `centromere_like` or `sample_specific` C11 candidate inside 152.0-154.5 Mb,
  and a J10 candidate inside 163.2-165.7 Mb.

negctl:
- N1: the gate fails on every chromosome, no call changes, and every integrated call equals its sequence call.

rice, potato (held-out functional truth):
- H1, no harm: integrated calls are right on at least as many chromosomes as sequence calls.
- H2, changes: every change is listed with its outcome. Pass if right changes >= wrong changes. Zero changes passes
  and is reported as such.
- Reported, no pass/fail: each line's right count and P; the gate outcome; satellite-free potato chromosomes (Cen4,
  6, 10, 11, 12) by line.

mammals (import mode):
- M1, no harm: per genome, integrated right >= sequence right.
- M2, anchoring: anchored Hi-C peaks are right on more chromosomes than de novo peaks, pooled over the 6 genomes.
  McNemar exact P <= 0.01 on the discordant pairs. Truth is CHM13 active HOR, the other genomes satellome.
- M3, gate calibration, reported: per genome, gate pass share against the anchored-peak right rate.
- Caveat stated with every number: outside CHM13 the truth is satellite-based, like the sequence line.

stretch: the same scoring as rice and potato (barley, horse) or as mammals (gorilla); reported if finished.

## Verdict

Each criterion pass or fail, with its numbers. Every failure is investigated and reported, not explained away.
Results are appended below this line and never edit the text above.

## Results

Run 2026-10-16 22:38 to 2026-10-17 03:24. Calls come from code md5 3cbfadb519dc; two later reporting-only changes
(a6f32a623f61, then a8f8c0d889d2) left every call unchanged.

Rule changes made after a run, and the genome that prompted each (design.md, Changes after the first runs):
- `telomeric` peaks never oppose: rice. Rice is no longer a fair test of this rule.
- Hi-C opposes only beyond 4 bins: CHM13, a development genome.
- Code-review fixes: before scoring.

Every genome was rerun on the final code. All criteria are scored on those reruns.

### Summary by criterion

| criterion | result | numbers |
|---|---|---|
| P1 Pa1A integrated call on 167-170 Mb | PASS | 163.0-169.4 Mb (the CRM core; the sequence call already moved there through cenrep's CRM rule); concordant, high; supported by TE, Hi-C (peak 168.0 Mb) and methylation |
| P2 the other 20 Poa calls unchanged | PASS | 20/20 |
| P3 Poa gate | PASS | annua 14/14, infirma 7/7 |
| P4 methylation supports non-Pa1A annua calls | PASS | 12/13; Pa2A `dip` (as in cenmeth's own validation) |
| P5 Pa1A plant-specific candidates in the sample table | FAIL | C11 153.8 Mb (P 0.004) listed. J10 164.8 Mb (rank 2, P 0.006) exists in methylation/annua.candidates.tsv but is not in the sample table: it lies inside the reference's wide integrated call (163.0-169.4 Mb), and the table lists only candidates outside the call |
| N1 permuted Hi-C | PASS | gate 0/14, no change, 14/14 calls = sequence calls |
| H1 no harm (rice, potato; also horse) | PASS | rice 12/12 = 12/12; potato 2/11 = 2/11; horse 3/4 = 3/4 |
| H2 changes right >= wrong | PASS (no changes) | no call moved in any held-out genome |
| M1 no harm (6 mammals) | PASS 6/6 | integrated = sequence right count in every genome |
| M2 anchored beats prior-free Hi-C (mammals pooled) | FAIL | 73 vs 73 of 138 (15 vs 15 discordant; McNemar P = 1). Per genome: chm13 17 vs 16, chimpanzee 19 vs 15, S. orangutan 7 vs 4, B. orangutan 8 vs 16, siamang 22 vs 22, giraffe 0 vs 0 |
| M3 gate calibration (reported) | see below | gate passed in every genome where anchored Hi-C was mostly right, and failed in every genome where it was mostly wrong |

### Gate calibration across all scored genomes

| genome | gate passes | anchored Hi-C right (truth) | prior-free Hi-C right |
|---|---|---|---|
| chm13 | 24/24 | 17/23 | 16/23 |
| chimpanzee | 25/25 | 19/25 | 15/25 |
| siamang | 26/26 | 22/26 | 22/26 |
| rice | 12/12 | 8/12 | 0/12 |
| P. annua / P. infirma | 14/14, 7/7 | 14/14 and 7/7 on the sequence call (no ChIP) | 14/14, 7/7 |
| B. orangutan | 0/25 | 8/24 | 16/25 |
| S. orangutan | 0/25 | 7/25 | 4/25 |
| giraffe | 0/15 | 0/15 | 0/15 |
| potato | 0/11 | 0/11 | 1/10 |
| horse | 0/5 (all sequence calls low: no reference) | 3/5 | 3/5 |
| permuted Hi-C (negctl) | 0/14 | 0/14 on the sequence call | 0/14 |

- Where the gate passed, anchored Hi-C was right on 67-85% (truth) of chromosomes.
- Where it failed, anchored Hi-C was right on 0-33%. The exception is horse (3/5), where the gate could not be evaluated:
  no sequence call there was medium or high.
- Hi-C changed no call in any genome. Where it was validated, it added an independent supporting line (relation
  `concordant`): rice 8, chm13 14, chimpanzee 18, siamang 22, Poa 21 calls. Confidence rose only where the sequence
  call was below high.

### What failed, and why
- **Potato** (sequence 2/11).
  - Without TE input, cenrep calls short pericentromeric satellites. This is the Solanaceae limitation in cenrep's
    README; half of the potato centromeres are satellite-free.
  - Hi-C carries no centromere signal here: anchored 0/11, prior-free 1/10.
  - Neither line can rescue the other. The gate correctly kept Hi-C out.
- **Horse chr11** (satellite-free).
  - The sequence call is 13.6 Mb off; both Hi-C peaks are 32.7 Mb off, at the chromosome end. Missed by every line.
  - The acrocentric CM009178.1 has no sequence call; prior-free Hi-C was right there, but Hi-C was not validated.
- **M2.** Anchoring did not beat cenhic's prior-free location in mammals. In both orangutans and the giraffe, strongly
  self-associated heterochromatin captures the anchored peaks too (interstitial blocks 18x in S. orangutan). The
  advantage of anchoring is genome-specific (rice 8 vs 0). The pre-registered general claim fails.
- **CHM13.** 3 of 24 correct calls are `contested` (lowered to medium):
  - chr20, peak 2.05 Mb outside the HOR array;
  - chr14 and chr22, peaks on the acrocentric p-arm (NOR territory; 45S rDNA units are not detected).

### Added after the pre-registered runs (not criteria)

**The Hi-C move rule on real data: P. annua Pa1A.** No call moved in the pre-registered runs. With `tes/`, cenrep's CRM
rule had already placed Pa1A on its CRM core.

To exercise the integration's move rule:
1. Run the satellite engine on P. annua without TE input. Pa1A is then called at its old satellite array,
   146.12-146.76 Mb, satellite, high, as cenrep did before TE input.
2. Import those calls (`--satellite-dir`) into a run with the TE, Hi-C and methylation lines.

What the integration did on Pa1A:
- The gate passed: all 12 reference chromosomes agreed, P = 2.6e-20.
- The anchored peak was at 168.0 Mb, level high, core 167.0-170.0 Mb.
- The CRM core (163.0-169.4 Mb) lies within 4 bins of the peak and 16 Mb from the satellite: a second line.
- The call moved to 167.0-170.0 Mb: `offset` 20.2 Mb, high, flag `moved_with_te`.

The loci table keeps the old array as a `sequence_centromere` (145.8-146.8 Mb: satellite and a methylation dip).
The other 13 calls stayed (`concordant`).

Pooled methylation was a `dip` at both places, not centromere-like. The methylation plants are five other accessions,
whose Pa1A centromeres sit at different positions (C11 ~153.8 Mb, J10 ~164.8 Mb, P3 on the old satellite).

**3D self-association of feature classes** (fold over random places on the same chromosomes; * = P <= 0.01):

| genome | centromeres | chromosome ends | telomere arrays | interstitial blocks | subtelomeric blocks |
|---|---|---|---|---|---|
| P. annua | 5.05* | 5.47* | (2 ends) | 2.08* | 4.20* |
| P. infirma | 6.05* | 3.63* | . | 1.00 | 3.13* |
| rice | 1.28 | 5.98* | 7.41* | . | . |
| potato | 1.25 | 5.42* | . | 0.83 | 5.66* |
| horse | 2.90* | 7.89* | 13.76* | . | 9.34* |
| CHM13 | 3.59* | 2.67* | 3.22* | 5.88* | 16.71* |
| chimpanzee | 3.85* | 2.68* | 6.16* | 4.73* | 2.87* |
| B. orangutan | 4.15* | 4.99* | 8.58* | 5.27* | . |
| S. orangutan | 2.76* | 5.72* | 7.31* | 18.04* | . |
| siamang | 4.13* | 1.56* | . | 1.86* | 1.69* |
| giraffe | 1.33 | 7.45* | 11.28* | 1.52* | . |
| permuted Hi-C (P. annua) | 0.82 | 1.01 | . | 0.98 | 1.04* |

Caveats:
- Representative bins sit inside repeat arrays, where mis-mapped reads could add trans contacts.
- The null is not compartment-matched.
- On permuted pairs, 3 of 27 classes reached P <= 0.01 at folds 1.04-1.21. Only large folds mean clustering.

Taken at face value, chromosome ends cluster in every genome, and centromeres cluster in Poa, horse and all primates
but not in rice, potato or giraffe. Rice's and potato's cytology is non-Rabl. This is a lead to test with uniquely
mapped flanks and a compartment-matched null, not a finding.

**Reporting changes after the scored runs** (calls unchanged; code md5 a6f32a623f61):
- **Loci.** Dispersed and minor satellite families are no longer loci members, and only Hi-C members are padded. In
  satellite-dense genomes, single-linkage merging had chained P. annua Pa1A into 40-114-Mb loci.
- **Sample table.** It keeps candidates inside a wide call (`inside_call`). P5 stays FAIL as pre-registered; the table
  now lists J10's Pa1A candidate at 164.8 Mb (rank 2, plant-specific P = 0.006, inside the reference call).
- **Overview figure.** The TE marker density has its own lane.

Every finished run's outputs were rebuilt with these (checkpoints reused).

### Barley (stretch genome, finished within budget)

| genome | bin | sequence | prior-free Hi-C | anchored Hi-C | integrated | gate pass | changes | H1 | H2 |
|---|---|---|---|---|---|---|---|---|---|
| barley | 2000 kb | 0/7 (P 1) | 6/7 (P 5.7e-09) | 0/7 (P 1) | 0/7 | 0/7 | 0 | PASS (0 = 0) | PASS (no changes) |

- **The satellite line is wrong genome-wide.** Without TE input it picked non-centromeric satellites 185-356 Mb from
  every CENH3 centromere, three of them at high confidence. Barley centromeres are Cereba (CRM) retrotransposons
  plus a 6-bp repeat outside cenrep's search.
- **Hi-C has the answer but cannot use it.** cenhic's prior-free calls are right on 6/7 (barley nuclei are Rabl).
  The anchored test sits on the wrong anchors, and the gate correctly finds that Hi-C disagrees with the sequence
  line. Two lines that disagree everywhere cannot say which is wrong, so cenfuse keeps the sequence calls and is wrong
  on 7/7.
- **The fix is a third line.** `tes/` (CRM cores) would very likely move the sequence calls, through cenrep's CRM
  rule, and the gate would then validate Hi-C. Not run: no barley TE annotation on disk.
- **Added after this run (reporting only).** The genome-level diagnostic `hic_sequence_conflict` fires when the
  prior-free Hi-C calls are confident on most chromosomes but disagree with the sequence calls on most and the gate
  fails. It warns that one line is wrong genome-wide and asks for `tes/` or methylation. It fires on:
  - barley (satellite line wrong);
  - S. orangutan and giraffe (prior-free Hi-C wrong);
  - potato (both wrong);
  - horse.
  It does not fire on rice, CHM13, chimpanzee, siamang, Poa or the permuted control.

### Gorilla (stretch genome, finished within budget)

| genome | bin | sequence | prior-free Hi-C | anchored Hi-C | integrated | gate pass | changes | no harm |
|---|---|---|---|---|---|---|---|---|
| gorilla (mGorGor1, satellome truth) | 250 kb | 25/25 (P 9.6e-28) | 0/24 | 7/24 (P 3e-4) | 25/25 | 0/25 | 0 | PASS |

- **Sequence and Hi-C.** cenrep 0.9.0 (FasTAN seeds) called all 25 on the satellome locus. cenhic's prior-free Hi-C
  converges on interior clusters (0/24, as in the earlier benchmark).
- **Anchoring and the gate.** Anchoring to the sequence calls lifts Hi-C only to 7/24. The gate fails on all 25, so
  Hi-C changes nothing, and `hic_sequence_conflict` fires.
- **3D test.** Centromeres self-associate 3.5x (P = 0.001). The cluster exists but is not the dominant trans cluster.
- **Telomeres** (TTAGGG) at 50/50 ends.
- Wall time: satellite engine 280 min at 16 threads (annotate 4.4 h); Hi-C about 2 min.

### Budget
- First launch 2026-10-16 22:38:06; last run finished 2026-10-17 03:23:40, within the 6-h limit.
- At most 118 threads were allocated at once.
- Postponed: nothing. Every genome in the plan, the stretch genomes included, finished.

## Plants and provenance (pre-registered 2026-10-19, before any run of the new code)

Design: `docs/plants.md`.

Truth for *P. annua* comes from earlier analyses of the same data (not part of this repository):
- Genotyping: none of the five ONT plants is the reference plant (PI 595837 "MN234"). They carry 0-3 of the
  reference's 62 young Pa1A CRMs; the reference HiFi carries 62/62.
- The Hi-C library SRR19374716 comes from the reference plant (young-CRM presence, mixed heterozygous sites).
- Read depth along the genome (karyotype):
  - C11, J10 and R1 lack Pa1A 171-255 Mb;
  - P3 and R1 lack Pa2A ~162-212 Mb (both exchanges replaced by a B-subgenome copy);
  - I3 has neither.
- Per-plant positions from plant-private young-CRM clusters, inferred rather than measured (no CENH3 ChIP).
  - The other 13 chromosomes are stable in every plant.
  - Pa1A varies: C11 ~153.6 Mb, J10 ~164.3 Mb; R1 ~146-148, I3 ~145-148 and P3 ~147-150 (near the satellite array).
  - The reference's functional Pa1A centromere is 167.0-170.0 Mb (Hi-C).
- Known limits:
  - CHG dips alone are not diagnostic (C11 has comparable dips without CRM clusters elsewhere).
  - I3 and J10 see ~2% of satellite CHG sites.

**Synthetic genome** (`tests/fuse/fuse_synth.py`, planted plants):

| | criterion |
|---|---|
| S1 | Provenance: P1 is `reference`; P2 is `other` (genotype and its planted missing segment); P3 is `other`. P3's Hi-C library is `other` by name. The default library is `unknown` and decides. |
| S2 | The reference calls keep every 0.1 planted outcome (chr7 moved by methylation, chr8 by the TE core, chr9 contested, chr10 functional_only), with only P1's methylation deciding. |
| S3 | P3 chr3 is `shifted` to its planted site by methylation and by its Hi-C library. No other plant × chromosome is `shifted`. |
| S4 | A methylation file on another assembly stops the run; one methylation sample alone has genotype `unknown`. |

***P. annua*** (annua; the five plants; the reference Hi-C):

| | criterion |
|---|---|
| R1 | All five plants are `other` by genotype; none is `reference`. |
| R2 | Missing segments are found for the five known events (C11, J10, R1: Pa1A ~171-255; P3, R1: Pa2A ~162-212), each within 5 Mb of the earlier analysis's bounds. None is found for I3. Any other segment is reported and inspected. |
| R3 | The Hi-C library is `unknown` (not `other`) and decides. |
| R4 | Every reference call keeps its position from the scored run (methylation never moved one there). Confidence changes from losing the methylation support are reported. |
| R5a | P3 and R1 at Pa2A are `call_absent`. |
| R5b | No plant is `shifted` on any of the 13 stable chromosomes. `possible_shift` verdicts there are counted and reported. |
| R5c | At Pa1A, no plant is `shifted` or `possible_shift` to a site more than 2.5 Mb from its truth window. Hits and misses against the truth windows are reported; low power is expected. |

### Results (plants and provenance, 2026-10-19)

**Synthetic genome:** S1-S4 PASS (`tests/fuse/test_fuse_e2e.py`, `tests/fuse/test_plants.py`; 432 tests pass).
- **Calibration of the between-plant test** on 300 pure-noise genomes (4 plants x 10 chromosomes):
  - 5.0% of 24,000 null tests have P < 0.05;
  - 14/300 genomes have any q <= 0.05.

***P. annua*** (one scored results file per run, kept with the run record):

| | pre-registered run | after the karyotype fix | after the review fixes (final) |
|---|---|---|---|
| R1 all five plants `other` | PASS | PASS | PASS |
| R2 the five known missing segments | FAIL: 0/5 found | FAIL as written (see below) | FAIL as written |
| R3 Hi-C `unknown`, decides | PASS | PASS | PASS |
| R4 call positions unchanged | PASS (no confidence change) | PASS | PASS |
| R5a P3, R1 Pa2A `call_absent` | FAIL (`no_data`) | PASS | PASS |
| R5b no `shifted` on the 13 stable chromosomes | PASS | PASS | PASS (`possible_shift`: J10 Pa4B, Pa5A) |
| R5c no shift more than 2.5 Mb from truth | PASS (no shift called) | PASS (no shift called) | PASS (no shift called) |

R2 after the fix:
- **Pa1A events: found exactly.** C11 171-255, J10 171-255, R1 170-256 Mb.
- **Pa2A events (P3, R1): found in three pieces** (166-180, 195-205, 208-212 Mb), with coverage back at 180-194 and
  206 Mb. The earlier analysis also says "lost in pieces". The criterion asked for one segment within 5 Mb of
  162-212, so it fails as written.
- **I3: one unexpected segment**, Pa1A 199-203 Mb. I3's coverage drops from ~0.8 to 0 while P3 covers the region.
  Not in the earlier analysis; a lead.

Held-out checks (not criteria):
- R1, monosomic for Pa5B (~0.5x), has no segment there.
- The reference Hi-C library has no segment.
- P3's 2.5-Mb Pa1A loss is too short to detect.
- No segment is confined to an all-sample low-coverage patch. This check prompted the same-type rule, so it no
  longer tests it.

**Pa1A, plant by plant (final run) against the young-CRM positions:**

| plant | cenfuse | truth window (Mb) | relative to the reference call (163.0-169.4) |
|---|---|---|---|
| C11 | consistent with the call (0.68 at the call; 0.26 near 154.9 Mb) | 152.0-154.5 | shifted |
| I3 | unclear (dip at the call) | 145.2-147.7 | shifted |
| J10 | consistent with the call | 163.1-165.6 | at the call |
| P3 | consistent with the call; also centromere-like at 146.4 and 148.5 Mb | 147.3-149.8 | shifted |
| R1 | consistent with the call (0.6) | 146.0-148.5 | shifted |

- **Agreement.** Methylation agrees with the earlier insertion-based positions only for J10. C11, P3 and R1 look
  centromeric at the reference call although their insertion clusters lie elsewhere.
- **Why.** Some places that are not centromeres look centromeric in every plant (Pa4B 58.3 Mb: 2.4-3.4; Pa5A
  71.5 Mb: 1.1-2.6). The earlier analysis found the same ("CHG dips alone are not diagnostic").
- **What holds.** The between-plant test keeps these places from producing `shifted` verdicts. But on real ONT data
  at 2.5-8.6x, per-plant methylation does not resolve the Pa1A polymorphism.
- **Not validated on real data.** The per-plant Hi-C path has no other-plant library here.

**Time.** The first run took 7.7 min at 16 threads (methylation bins from scratch: 357 s); reruns from checkpoints
take under a minute.
