# cenfuse design

Written 2026-10-16, before any code. An earlier prototype with a different Hi-C rule, its benchmark and its rule
audit are used here as evidence.

## Intent

Satellite-based centromere callers exist (CentIER, quarTeT, TRASH, ...), so a satellite caller alone is not new. The
value is in combining independent evidence lines and placing their candidates on one map:
- satellite arrays (the satellite engine, `cenfuse.sat`);
- LTR-retrotransposon clades (the satellite engine with `tes/`);
- Hi-C centromere clustering (the Hi-C engine, `cenfuse.hic`);
- ONT methylation (the methylation engine, `cenfuse.meth`).

For every chromosome, cenfuse reports:
- each line's own call;
- where the lines agree and where they disagree;
- one integrated centromere call with its evidence;
- the other features the lines find (knobs and other interstitial blocks, rDNA/NORs, subtelomeric blocks, contact
  hubs, plant-specific methylation candidates), placed relative to it.

The sequence centromere (the satellite array) and the functional centromere (where the kinetochore sits) are kept
apart, because they can differ (Poa annua Pa1A: 21.8 Mb).

## Goals, non-goals, anti-goals

Goals:
- One command runs every evidence line present for a set of genomes and fuses them.
- Each line is optional except the assembly. The tool degrades to the satellite engine alone when only FASTA is given.
- Every call carries its supporting and opposing lines, a confidence and the relation between the sequence and
  functional calls.
- A loci table and GFF3 co-localize all candidate features, including the ones that are not centromeres.

Non-goals:
- Calling centromeres from methylation alone (the methylation engine's design: false dips are as deep as true ones).
- PacBio methylation (a design decision: CHG from ONT/modkit only).
- Changing how any engine works on its own. The engines are cenfuse's own, and fusion leaves their logic unchanged.

Anti-goals, which would make the tool a failure even if it runs:
- A line overriding another because of its own self-reported confidence. The Hi-C engine's prior-free `high` calls
  (version 0.1.1) were right on 32/71 in T2T mammals (45%), and accuracy fell as confidence rose (medium 51%, low 61%).
- Hi-C deciding satellite family choice. Gorilla would likely lose alpha satellite.
- Fixing Hi-C's failures from inside the Hi-C data, with more self-diagnostics on the Hi-C data alone (random
  restarts, caps for competing or telomere-dominated solutions). A trial of this left held-out accuracy unchanged.
- Tuning a threshold until the benchmark looks right.

## What the design rests on

Validation available before this design (the engines' own benchmarks):
- **The satellite engine.**
  - CHM13 against the CENP-A active HOR: 23/24.
  - 25 T2T mammal assemblies: 413/567 agree, and `high` is right on 318/397 (80%).
  - CENH3 tests: Ae. geniculata 10/14; S. verrucosum 7/12 within 1 Mb.
  - Satellite-free centromeres are called only with `tes/` (CRM cores).
- **The prior-free Hi-C calls (Hi-C engine 0.1.1).**
  - P. annua and P. infirma: 20/20 controls inside the satellite call, and the Pa1A functional core at 167.5-170 Mb.
  - 8 T2T mammals: 97/188 right (gorilla 0/24, giraffe 0/15, S. orangutan 4/25).
  - Held-out CENH3/CENP-A truth: barley 6/7, rice 0/12 (every peak in an end bin), potato 1/10 (1-bin hotspots),
    horse 3/5.
  - The Rabl test passed on every failing genome.
- **Root cause** (an analysis of where the prior-free Hi-C calls fail):
  - The prior-free iteration converges on whichever trans-contact cluster dominates and is self-consistent (the
    telomere pole, interior A-compartment clusters).
  - Its bootstrap, leave-one-partner-out and null tests keep the partner windows fixed. They measure how stable a
    peak is, not whether its cluster is the centromere cluster.
  - A trial that tried to fix this from inside the Hi-C data left held-out accuracy unchanged, and cost correct
    bonobo calls their `high` level.
- **Anchored Hi-C works where the iteration fails.** In a spike test, partner windows were placed at the other
  chromosomes' satellite loci; the peak scored if within ±1 bin:
  - CHM13 18/23 (P 8e-31).
  - Chimpanzee 17/25.
  - P. annua 13/14. The one miss is Pa1A, the real functional offset.
  - In gorilla, a state started at the true loci scores accuracy 0.72, but the wrong interior state has more total
    contact (106 vs 75). So the iteration leaves the centromere cluster; anchoring keeps it.
- **The prototype's rule audit:**
  - Letting prior-free Hi-C `high`/`medium` calls decide: 103 decisive calls, 45-48 wrong.
  - A leave-one-chromosome-out genome gate cut this to 16 wrong (in-sample). It was adopted together with two more
    rules: Hi-C moves a high sequence call only with a second line, and Hi-C stays out of family choice.
- **The methylation engine on P. annua.**
  - The signature is CRM CHG down (-0.16).
  - Among the top 3 dips, `centromere_like` marks 48/55 sites within 2 Mb of the Hi-C peak, against 2/140 elsewhere.
  - It needs anchors and never calls de novo.

## Principle: cross-validation between lines, not self-reported confidence

1. **Each line calls on its own first.**
   - Satellite and TE calls come from the satellite engine (`cenfuse.sat`), which already combines the two
     with its CRM rules.
   - Hi-C de novo calls come from the Hi-C engine (`cenfuse.hic`).
   - Methylation candidates come from the methylation engine (`cenfuse.meth`).
2. **Hi-C answers a question it can answer: does this place contact the other chromosomes' centromeres?**
   - Partner windows sit at the other chromosomes' sequence calls (anchored, leave one chromosome out).
   - For chromosome x, the score uses only the other chromosomes' calls, so x's own call never shapes its test.
   - Which trans cluster dominates the genome no longer matters.
3. **Whether Hi-C is informative in this genome is measured, not assumed (the gate).**
   - On the other chromosomes, does the anchored Hi-C peak land on the sequence call more often than chance?
   - A Rabl test or a bootstrap cannot answer this.
4. **Methylation is calibrated on the integrated calls of the other chromosomes** (leave one chromosome out). It
   supports calls and serves as a second line; it never moves a call alone.
5. **The integrated call is the sequence call** unless a validated functional line moves it under the explicit rules
   below. Every move and every disagreement is reported.

## Inputs

`GENOME_DIR/`: the GENOME_DIR layout the satellite engine reads, plus two optional folders.

| path | needed | used by |
|---|---|---|
| `<genome>.fa/.fasta/.fna[.gz]` + `.fai` (bgzip also `.gzi`) | yes | all lines |
| `genes/<genome>.gff3[.gz]` | no | satellites (subtelomeres, gene-rich veto), methylation (gene class), loci (gene density) |
| `tes/<genome>.gff3[.gz]` (`LTR_retrotransposon` with `clade=`) | no; needed for methylation | TE line (marker clade, cores), methylation classes |
| `hic/<genome>.pairs[.gz]` or `.cool`/`.mcool`, or `hic/<genome>.tsv` (one library per plant) | no | Hi-C line |
| `meth/<genome>.tsv` (the methylation engine's samples table) + optional `meth/<genome>.exclude.bed` | no | methylation line |

Command:
`cenfuse run GENOME_DIR [-o OUT] [-t T] [--satellite-dir DIR] [--seed-engine ...] [--seed-cache DIR]
[--classify-from T=S] [--force STAGE] [-v]`.
- `--satellite-dir` reads existing satellite-engine outputs (summary TSV + arrays GFF per genome) instead of running the
  satellite engine. It exists so a batch of satellite-engine runs need not be repeated.
- The other options pass through to the satellite engine.
- The default output dir is `GENOME_DIR/cenfuse`, so a batch `centromere/` dir is never overwritten.

## Method

Stages, checkpointed; a rerun resumes:

```
satellite -> genotype -> hic -> karyotype -> provisional -> meth -> fuse -> plants -> loci -> output -> plot
```

genotype, karyotype and plants are described in `docs/plants.md`. Inputs shown to come from another plant than the
assembled one never decide the reference calls, and every plant is scored on its own.

### satellite

The satellite engine's Pipeline (`cenfuse.sat`), or `--satellite-dir`. It produces:
- the sequence call per chromosome (interval, evidence, confidence, flags, main family);
- the TE marker clade and its core;
- every array with its family role (centromeric, cen_associated, subtelomeric, interstitial, dispersed, minor, rDNA).

TE cores of every concentrated clade, not only the marker, are recomputed from `tes/` with the satellite engine's
`te` module (they feed the loci).

### hic (per genome with `hic/`)

1. **Run the Hi-C engine de novo (prior-free), unchanged.** Its outputs are kept under `OUT/hic/`, and the matrix,
   balancing and block means are reused.
2. **Anchored score.** For every callable chromosome x:
   - Partner windows: each other chromosome's sequence-call midpoint bin ± 4 bins (the Hi-C engine's HALF_WIN).
   - Partners come only from chromosomes whose call is at least medium; all calls are used if fewer than 3 qualify.
   - Score: the Hi-C engine's centromere-cluster score (mean over partners of O/E with the partner window). x's own
     call is never a partner, so the score is leave-one-chromosome-out by construction.
   - Peak and core: the Hi-C engine's peak_core (half-maximum core).
3. **Support for every chromosome's anchored peak**, with the Hi-C engine's own tests:
   - 100 Poisson bootstraps (re-balanced, partner windows fixed);
   - leave one partner out;
   - 200 null replicates with the partner windows moved to random interior places;
   - the Hi-C engine's `level()` (high = all three) and flags (hub, rDNA, telomeric).
   - A peak is **strong** if its level is high and it carries no `hub`, `rDNA` or `telomeric` flag. Hubs and NORs
     contact everything, and end windows hold the chromosome-end cluster. `telomeric` was added after the first
     benchmark run (see Changes after the first runs).
4. **Agreement.** A Hi-C core agrees with an interval if the gap between them is at most 1 bin (placement error, as
   in the Hi-C engine's window rule).
5. **Gate, per chromosome x, leave one chromosome out.**
   - Reference set R: chromosomes other than x whose sequence call is high (high plus medium if fewer than 3 are
     high), with a finite anchored score.
   - k = the members of R whose anchored peak agrees with their own sequence call; n = |R|.
   - Chance of agreeing on chromosome y, p_y: the share of y's scored bins whose ±1-bin neighborhood touches y's call.
   - The gate passes if k > n/2 and the Poisson-binomial tail P(K >= k) <= 0.01. These are the benchmark's
     pre-registered C2 rule and the conventional level.
   - The same statistic for the de novo peaks is reported as a diagnostic, not used.

### provisional (per genome; Hi-C rules that need no methylation)

| sequence call | gate | anchored Hi-C on x | result |
|---|---|---|---|
| any | no Hi-C, or gate fails | any | sequence call; relation `sequence_only` (Hi-C reported) |
| present | passes | agrees | sequence call; `concordant`; Hi-C counts as a supporting line |
| none | passes | strong | Hi-C core; `functional_only` |
| low/medium | passes | strong, disagrees | Hi-C core; `offset` (sequence call kept as the sequence centromere) |
| high | passes | strong, disagrees | pending: moves only with a second line at the Hi-C core (fuse) |
| present | passes | not strong, disagrees | sequence call; flag `hic_disagrees_weak` |

### meth (per genome with `meth/` and `tes/`)

1. **Build the methylation engine's class bins** (cached). The depth filter and exclusions are unchanged.
2. **Anchors**: the provisional calls of confidence at least medium.
3. **Signatures.**
   - One calibration with all anchors gives the genome signature and the scan direction. The methylation
     engine's standard outputs are written under `OUT/methylation/`.
   - For every anchored chromosome x, a second calibration leaves x's anchor out and gives the leave-one-out
     signature sig_x.
   - The scan direction stays the genome's; it is one sign learned from all anchors, a negligible leak.
4. **Methylation support** for a candidate interval I on x:
   - Find the most extreme window of the adjusted CHG track (in the learned direction) inside I, widened by half a
     window, among the bins where the signature has data.
   - At that window, compute the methylation engine's signature score with sig_x, on class sums pooled over the
     samples that may decide (the reference plant's, or of unknown origin), with a signature calibrated on them alone
     (`docs/plants.md`).
   - `centromere_like` if score >= 0.5 with signature data (the methylation engine's rule).
   - The intervals evaluated are the sequence call, the anchored Hi-C core, the TE core and the provisional call.
5. **Each plant on its own.** Every sample is scored at every site of every chromosome and compared with the other
   plants there (`docs/plants.md`). The methylation engine's per-sample candidates, labeled with sig_x, add sites.

### fuse

- **Pending moves.** A high sequence call moves to the strong Hi-C core only if a second line supports the Hi-C
  core and not the sequence call:
  - the TE marker core agrees with the Hi-C core; or
  - methylation is `centromere_like` at the Hi-C core and not at the sequence call.
  Otherwise the sequence call stays, with relation `contested` and flag `hic_offset`.
- **Methylation.** Support only: it is the second line above, and it raises confidence (below). It never moves a
  call alone, and a `dip` is not evidence against (as in the methylation engine).
- **Supporting lines of the final call** (independent lines only):
  - `satellite`: the main centromeric family has arrays inside the call;
  - `te`: the marker core agrees;
  - `hic`: validated and agrees;
  - `meth`: `centromere_like`.
- **Confidence**, in the satellite engine's levels and adjusted one step at a time (its te_concordant convention):
  - A kept sequence call starts at the satellite engine's confidence. It goes up one level if validated Hi-C agrees,
    and up one level if methylation supports it; at most two raises, ending at high at most. It goes down one level if
    `contested`.
  - A moved or functional_only call is medium, and high with a second line.
- **Relation**: `concordant`, `offset` (with the distance), `contested`, `functional_only`, `sequence_only` or
  `none`.

### loci

All candidate intervals on a chromosome:
- satellite clusters, at least 20 kb (the satellite engine's S_MIN), of families with a feature role (centromeric,
  cen_associated, interstitial, subtelomeric, rDNA) or an rDNA tag. Dispersed and minor families are background, not
  members (their 3D association is still tested per family);
- the core of every concentrated TE clade;
- the anchored and de novo Hi-C cores, and Hi-C hub runs;
- methylation candidates per sample that are `centromere_like` or plant-specific;
- telomere arrays at the chromosome ends;
- the sequence call and the integrated call.

Overlapping intervals are merged into loci. Only Hi-C members are padded by 1 bin (their placement error) for the
merge, so dense satellite content cannot chain a locus along a whole chromosome. Each locus lists its members and the
value of every line at the locus:
- satellite bp by role, and telomere-array bp;
- LTR element count and the most common clades;
- the anchored, de novo and telomere-pole Hi-C scores over the chromosome median, and whether it is a hub;
- pooled methylation deviation and raw CHG level, and how many plants have a `centromere_like` or plant-specific
  candidate there;
- gene-density ratio and gap bp.

Its kind is the first that matches (the kind labels a locus; only the centromere call combines lines by rule):
1. `centromere` (holds the integrated call);
2. `sequence_centromere` (holds the sequence call when the integrated call moved);
3. `telomere` (a telomere array);
4. `rDNA` (a satellite family whose consensus holds an rRNA gene, or whose largest arrays sit among 45S rRNA genes);
5. `subtelomeric_block` (a family with most of its bp beyond the outermost genes, or at the ends of most chromosomes);
6. `interstitial_block` (a family concentrated at one locus per chromosome, away from the centromere: knobs and other
   heterochromatic blocks, which the satellite engine cannot tell apart);
7. `alt_centromere_candidate`: away from the call, a strong anchored Hi-C peak where the gate passes, a core of the
   marker TE clade, or `centromere_like` methylation in at least one plant;
8. `hic_hub`;
9. `hic_peak` (any other Hi-C core);
10. `te_core` (a core of a concentrated clade that is not the marker);
11. `satellite_block` (a cluster of a centromeric or cen_associated family away from the call);
12. `meth_candidate` (a plant-specific methylation candidate only);
13. `other`.

## Outputs (`OUT/`)

| file | content |
|---|---|
| `<genome>.centromeres.tsv` | one row per chromosome: integrated call, relation, confidence, supporting and opposing lines, flags; the sequence call; TE core; de novo and anchored Hi-C peak, core, level and gate; methylation scores at each interval |
| `<genome>.centromeres.gff3` / `.bed` | integrated calls (GFF3 attributes: relation, support, confidence) |
| `<genome>.loci.tsv`, `<genome>.features.gff3` | the loci and every member feature |
| `<genome>.provenance.tsv` | per input: reference plant, other plant or unknown, with the evidence (`docs/plants.md`) |
| `<genome>.plants.tsv`, `<genome>.plant_sites.tsv` | per plant and chromosome: its verdict relative to the reference call; every plant scored at every site |
| `<genome>.hic_anchored.bedgraph` | anchored Hi-C score per bin (over the chromosome median) |
| `cenfuse_summary.tsv` | one row per genome: what each line found and whether it was validated (gate pass count, de novo vs anchored agreement, Rabl, TE marker, methylation signature), calls by relation and confidence |
| `cenfuse_overview.pdf` (+ PNGs per genome) | per chromosome: satellite arrays by role, TE marker density, anchored and de novo Hi-C score, methylation track, sequence and integrated calls; with plants, a second page of plant verdicts |
| `satellite/`, `hic/`, `methylation/` | each engine's own outputs, unchanged in format |

## Constants (`cenfuse/params.py`)

| constant | value | reason |
|---|---|---|
| HIC_TOL_BINS | 1 | one bin of placement error (the Hi-C engine's window rule) |
| GATE_P | 0.01 | conventional; the Hi-C engine's P_HIGH; the benchmark's pre-registered C2 rule |
| GATE_MIN_REF | 3 | fewer than 3 reference chromosomes give no majority test |
| ANCHOR_LEVELS | high, medium | partner windows and methylation anchors from calls at least medium (wrong anchors weaken both) |
| STRONG_HIC | level high and no hub/rDNA/telomeric | the Hi-C engine's three tests all pass; hubs and NORs contact everything; end windows hold the chromosome-end cluster |
| METH_SCORE | 0.5 | the methylation engine's SIG_FRACTION |

The provenance and plant constants are listed in `docs/plants.md`.

## Validation

Tests (pytest):
- Unit tests of the gate statistic, interval agreement, the decision table, confidence arithmetic and loci merging.
- A synthetic end-to-end run at the engines' real scales. Sequence calls are given as satellite-engine-format tables;
  Hi-C pairs are simulated in a Rabl nucleus; bedMethyl is simulated with CRM CHG loss at the centromeres. Planted:
  - one chromosome whose functional centromere is 15 bins off its satellite, with CRM support (must move);
  - the same offset without a second line on a high call (must stay, `contested`);
  - a hub decoy;
  - a satellite-free chromosome (functional_only);
  - non-Rabl (permuted) Hi-C (no moves);
  - three plants: the reference plant, another plant missing a segment, and another plant whose chr3 centromere
    sits elsewhere, with its own Hi-C library (`docs/plants.md`).
- One small run of the full pipeline, satellite engine included, on the satellite engine's synthetic test genome.
- The engines (`cenfuse.sat`, `cenfuse.hic`, `cenfuse.meth`) keep their own test suites.

Benchmark: `docs/benchmark.md`. Criteria are fixed there before any real run, within 6 h wall and 120 threads.

## Known risks and open questions

- The gate trusts the sequence line as reference. When the sequence line is weak genome-wide (no dominant satellite,
  no `tes/`), Hi-C cannot be validated and only reports, even if it is right (barley without TE input). With `tes/`,
  the CRM cores give the gate a second reference. Not implemented now: the TE core as an alternative gate reference.
- Single wrong Hi-C peaks in a genome where Hi-C mostly works (the audit's acrocentric and Y residuals) pass the
  gate. The second-line rule protects high sequence calls; lower calls can still be moved wrongly.
- Methylation LOCO signatures depend on the other anchors; with fewer than 3 anchors with data, methylation is not
  used.
- Hi-C resolution is the bin (250 kb-2 Mb), so an offset smaller than 1 bin is not resolved.

## Changes after the first runs

**2026-10-16, after the rice run (held-out, CENH3 truth).** The rice run used the rule above as first written.

What it showed:
- The anchored Hi-C peak landed on the sequence call on 8/12 chromosomes and the gate passed on 12/12. The prior-free
  Hi-C peaks: 0/12, all at chromosome ends.
- On chr4 (CP132238.1) a strong anchored peak in the first bin made the correct CentO call `contested`. Its
  confidence dropped from high to medium; the call did not move.

The change: an anchored peak flagged `telomeric` (within the Hi-C engine's 8-bin end window) may support a call but
never opposes or moves one, as for hubs and NORs.

The reason is independent of this chromosome. The Hi-C engine's earlier benchmarks showed chromosome-end clusters to be
the main non-centromeric trans attractor: rice 12/12 and giraffe 15/15 prior-free peaks at ends, and 33 of the mammal
misses.

Effect on the test: rice chr4 no longer tests this rule fairly. The genomes not yet scored (potato, barley, horse,
gorilla, the six mammals) test it. Cost: a functional centromere near an end, on a chromosome whose sequence call is
elsewhere, cannot be moved there by Hi-C.

**2026-10-16, feature classes and telomeres (to label more features).** Added:
- telomere arrays found de novo (`telomere.py`);
- the 3D self-association test per feature class (`hicfuse.feature_contacts`);
- per-locus telomere-pole score, raw CHG level and LTR composition, and the per-kind signature table
  (`loci.feature_signatures`).

These report only; none changes a centromere call.

**2026-10-16, after the CHM13 run (a genome used to develop the Hi-C engine).** Five correct CHM13 calls were
`contested`:
- three by anchored peaks 0.5-2.0 Mb outside large HOR arrays;
- two (chr14, chr22) by p-arm peaks.

Change: Hi-C still supports a call within 1 bin, but opposes or moves one only beyond 4 bins (`HIC_FAR_BINS`, the Hi-C
engine's partner half-window, the resolution the method claims and the benchmark's tolerance). Between 1 and 4 bins it
does neither (flag `hic_near`).

The reason is the method's resolution, not the CHM13 answers. The genomes not yet scored test it.

The two p-arm peaks are NOR territory. The satellite engine does not detect 43-kb rDNA units, so those peaks carry no
`rDNA` flag and can still contest correct calls on acrocentric chromosomes (a limitation).

**2026-10-16, independent code review (before any scoring).** Fixes:
- **Second lines are symmetric.**
  - TE is a second line only for a decisive (CRM) marker whose core lies within 4 bins of the Hi-C peak and beyond
    4 bins of the sequence call. The satellite engine's width-scaled `agrees()` had let a wide Hi-C core claim a TE
    core that sat nearer the call.
  - Methylation is a second line only with `centromere_like` at the peak and a measured `dip` at the call. No data
    is not evidence.
  - Every place is judged in a window of one width, max(2 D, Hi-C bin), centered on it, so a wide core gets no more
    chances.
- **Support needs a call.** An anchored peak at Hi-C engine level `none` never supports a call.
- **Contested keeps methylation support.** A contested call with methylation support keeps its level (-1 +1), as
  written above.
- **A core that reaches the call is not an offset.** A Hi-C core within 1 bin of the call never moves or contests it.
- **Optional lines never abort the run.** An unusable Hi-C or methylation line is skipped with a warning. The TE
  reader for methylation uses the TE line's type rule (any type containing `LTR_retrotransposon`).
- **Bootstrap memory is bounded.** Hi-C bootstrap threads are capped at 8 (memory).
- **The matrix cache works.** The chromosome `.fai` given to the Hi-C engine is rewritten only when it changes, so
  its matrix cache hits.
- **Checkpoints see all their code.** Stage keys hash every module of the engine package they call, cenfuse's
  params and util, and (methylation) the assembly index.
- **3D smoothing stays inside chromosome pairs.** It no longer runs across chromosome boundaries, which had inflated
  the end classes.
- **Telomeres:** the end array holding the most unit copies is reported, not the nearest run.

**2026-10-17, after the barley run (held-out, CENH3 truth).**
- What it showed: the satellite line was wrong on 7/7 and the prior-free Hi-C calls were right on 6/7, but the gate
  could not validate Hi-C against wrong anchors. This is the known risk stated above.
- Added: a reporting-only genome-level diagnostic, `fuse.sequence_conflict`. It fires when the prior-free Hi-C calls
  are confident (Hi-C engine level high/medium) on a majority of chromosomes, agree with the sequence calls on fewer
  than half, and the gate fails on most. It is logged as a warning and recorded in `cenfuse_summary.tsv`.
- No call changes. Two lines cannot decide which is wrong; a third line (`tes/` CRM cores or methylation) can.

## Changes after the T2T mammal benchmark (2026-10-23)

The benchmark (14 T2T mammals, CENP-A truth for 157 chromosomes in 9 genomes; run record outside the repository)
exposed the weak points below. Each change is stated with its principle and the case that prompted it. Every
change was replayed on the plant benchmarks and on a 185-genome plant and fish satellite batch before it was kept;
the replay results are in docs/benchmark.md.

**What a satellite is (classify).** A family seeded by one tandem copy of an interspersed element (LINE, SINE or LTR
fragment) is found genome-wide by nhmmer and chained into "arrays" of partial hits. In five FASTA-only mammal runs
such a family (461-1,861-bp unit, 215-850 Mb, on every chromosome) became the centromeric family; only the gene
veto removed it. Three gene-free measures now judge every family:
- `tandem_frac`: share of its array bp inside de novo tandem repeats (the seed table, any period);
- `whole_frac`: share of its hit bp in whole-unit hits (>= FULL_MONOMER_COV of the model);
- `locus_fold`: its main cluster's density over its mean density on the chromosome (median over chromosomes).

A family with tandem_frac < 0.5, whole_frac < 0.5 (majority rules) and locus_fold < 2 (LOCUS_FOLD, a heuristic, as
the Hi-C engine's HUB_FOLD) is a dispersed repeat, never centromeric. Mammal interspersed families: tandem 5-14%,
whole <= 1%, fold 1.1-1.3. Every CENP-A-validated centromeric family: tandem >= 45%, whole >= 90%. The fold test
keeps retrotransposon-derived arrays that cluster in centromeres (Coffea, Glycine, Rosa: not tandem, fold >= 2.5),
which the first version of the rule removed on 160 batch chromosomes.

**Dominance per satellite kind (classify).** The families stage splits a satellite into subfamilies at 90%
identity (6 alpha subfamilies in bonobo). Dominance counted the largest cluster per family, so a 43-bp pericentric
satellite out-counted alpha and took 12/25 bonobo calls, 4 of them 1-5 Mb from CENP-A (3 high). The largest cluster
per chromosome is now judged per kind (families related by `same_sequence`, arrays pooled). Effect beyond
bonobo: F. chiloensis D, where no family was dominant, now calls the CENH3-validated 147-bp Type I satellite on 6
more chromosomes.

**Acrocentric karyotypes (classify).** A family at one end of most chromosomes was subtelomeric. In goat (all
autosomes acrocentric) that made satellite I subtelomeric and left 20/31 chromosomes without a call. Subtelomeric
repeats occupy chromosome ends independently; a centromere occurs once per chromosome. A family that would be
subtelomeric is the acrocentric centromere when it is strong, sits at an end of most chromosomes, and at both ends
of so few that the ends cannot be independent (binomial lower tail on the both-ends count, expected n q^2 from the
per-end occupancy q; P <= 0.01, ACRO_ALPHA), and no strong interior family competes. Power: with n chromosomes all
one-ended, P = (1 - q^2)^n, so about 16 chromosomes are needed; smaller karyotypes keep the old rule. Calls on such
a family are flagged `acrocentric` (and `near_end`, so at most medium).

**Classify's choice makes the call (centromeres).** The call took the largest cluster of any centromeric family,
even where classify had chosen another family for that chromosome (siamang: an 80-bp satellite chosen on 2
chromosomes took 11 calls; CHM13 chr21: an inactive alpha array 2 Mb from the active one). The call is now the
pooled-family cluster that holds the chosen families' main cluster (layered satellites still share one call); a
decisive CRM core outranks the choice in favour of the largest cluster (S. verrucosum chr07).

**Competing clusters (centromeres, fuse).** A second cluster at least half the main one costs two levels instead of
one (on CENP-A or CENH3 truth, 2 of 10 such calls were right: mammals 2/5, potato and barley 0/5), and its coordinates
are reported (`seq_alt_start/end`). Validated Hi-C that agrees with the second cluster and not the first chooses it
(`hic_chose_alt`): Hi-C chooses between two sequence candidates; it does not place the call.

**Hi-C in mappable-bin geometry (hicfuse, fuse).** With MAPQ >= 30, the CENP-A domain lay in a block of unmappable
bins on 139/157 chromosomes; partner windows (call midpoint +/- 4 bins) fell entirely inside such blocks (0/22 in
each macaque, 1/20 in pig), so Hi-C never ran there, and anchored peaks sat at the block edges, 1-5 Mb from CENP-A.
- Partner window: the 9 mappable bins nearest the call, at most 4 mappable bins outside it (unchanged when the call
  is mappable).
- Every distance involving a Hi-C interval counts mappable bins only (`visible_gap`), and the gate's chance uses the
  same geometry, so the gate stays calibrated: a peak at an array's edge is adjacent to a call inside it.

**Gate `untested`.** With fewer than GATE_MIN_REF reference chromosomes the gate is reported `untested`, not `fail`.

**Hi-C-only calls from clean medium peaks (fuse).** With no sequence call, a gate pass and a clean anchored peak of
Hi-C level medium, the call is `functional_only`, confidence low, flag `hic_medium_peak` (horse chr11, the
satellite-free centromere: the medium peak lay inside the CENP-A domain and was discarded).

**NOR peaks (hicfuse).** A peak that could oppose or place a call is searched for 18S rRNA genes (peak bin +/- 1
bin, Rfam SSU, cmsearch --hmmonly); a hit flags it rDNA, so it never opposes or places a call.
