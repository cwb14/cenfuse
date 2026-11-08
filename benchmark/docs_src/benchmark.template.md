# Benchmarks

cenfuse was scored on 55 assemblies (56 genomes and libraries), from T2T mammals to fungi, against centromeres mapped
functionally on the same assembly. Every number on this page comes from one release rerun of cenfuse 0.14.0 on every
benchmark arm and on the 40 arms held out from 0.12 (260 runs, 2026-11-08), scored with the rules fixed before each benchmark
was first run; 0.14 changed no rule, and its calls equal 0.13.0's on all 260 runs ([section 13](#13-what-changed-from-013)).
Each run read a replay of the satellite engine on its arm's stored search (`--satellite-dir`); the engine's rules are
unchanged since 0.11. Two genomes whose CENH3 truth was on disk but never scored, an allotetraploid among them, are held out
from 0.14 and run in full, search included ([section 3](#held-out-from-014-an-allotetraploid-and-a-wild-potato)), and so are 13 tomatoes ([section 3](#held-out-from-014-thirteen-tomatoes)); the genomes held
out from the 0.12 rules are in [section 4](#held-out-from-the-012-rules). The scorer, the truth manifest and this page's tables are
in the repository ([benchmark/](../benchmark/README.md)). The overview pages of every genome are in the [gallery](gallery.md); how
each rule works is in [design.md](design.md).

## How calls are scored

- **Truth.** CENP-A CDRs (Logsdon et al. 2024: human CHM13 and five great apes) or CENP-A ChIP (crab-eating macaque,
  pig, horse; pig and horse from other animals); CENH3 ChIP-seq of the assembled line (rice, potato, barley, maize,
  *Brachypodium*, *Ae. tauschii*, soybean, *B. rapa*, tomato, moss, *Chlorella*) or of a related line (sorghum: line V9
  on the BTx623 assembly; *Arabidopsis*: Col-0 on Col-CC), CUT&Tag (*Chlamydomonas*), the published centromeres of two
  *B. rapa* accessions whose ChIP is too weak, the CenSat active arrays of HG002, the authors' CENP-A windows of
  chicken, the CENP-A minor-satellite array of mouse, and CEN/CENP-A/CenH3 coordinates of the fungi and *Plasmodium*.
  ChIP domains were called with one frozen rule (100-kb bins, ChIP/input, half-maximum run, ≥ 4×). The scores do not hang on
  its two numbers: on the 18 genomes whose truth it made and whose reads are on disk (60 release runs), 3× or 5× instead of 4×,
  or 50-kb bins (which move 167 of 186 domain ends), change no outcome on the chromosome-runs both rules give truth (604-624);
  they add or drop chromosomes whose ChIP sits near 4×, scored right less often (3×: 26 more chromosome-runs, 21 right; 5×: 20
  fewer, 16 of them right; v13 `bench/truth_sensitivity.py`, since 0.14).
- **Right**, as the frozen scorer rules it ([benchmark/score.py](../benchmark/score.py), manifest
  [benchmark/truth.tsv](../benchmark/truth.tsv)). Mammals with CENP-A truth, *Arabidopsis* and HG002 (a CenSat active array):
  the call overlaps the domain by at least 1 bp; the mammals without functional truth are scored by the same overlap with
  another satellite caller's locus (agreement, not truth). Chicken and mouse: overlap or touch (gap 0). Rice, potato, barley and
  EquCab3: within 4 Hi-C bins of the run, and without Hi-C 1 Mb (barley 8 Mb). Maize, *Brachypodium* and the CENH3 truths
  added since 0.8 (*B. rapa*, soybean, sorghum, *Ae. tauschii*, tomato, grape, cucumber; 0.14's *Ae. geniculata* and
  *S. verrucosum*): 1 Mb whatever the Hi-C bin (*Ae. tauschii*'s 2-Mb bins: half a bin). Moss, *Chlamydomonas*,
  *Chlorella*: 100 kb. Fungi and *Plasmodium*: 1 Hi-C bin of the run (10-50 kb), and 250 kb in an arm without Hi-C, a
  lenient rule for genomes whose centromeres span 1-100 kb ([future.md #37](future.md#37-scoring-rules-and-controls-to-tighten)).
  Holocentric genomes have no single centromere, so any confident call there is wrong.
- **Confident** = `high` or `medium`; `low` = cenfuse says it is unsure.
- **Chance.** A chromosome's chance of a right call by luck is (call + truth + 2 × tolerance) / chromosome length; P
  is the Poisson-binomial tail of the right count.
- **How close a right call is** (since 0.14, [future.md #37](future.md#37-scoring-rules-and-controls-to-tighten)). Overlap
  and tolerance rule a call right, however wide it is. Every right call of the real arms (the randomized-Hi-C controls left out),
  against its domain ([benchmark/precision.py](../benchmark/precision.py)). In the animals a call is the satellite array (or
  a large part of it), about ten times the CENP-A domain inside it: right says where the centromere is, to about half a
  megabase, not which part of the array holds the kinetochore (the CpG line's `cpg_kt_*` columns do, where it passes).

{{TABLE_PRECISION}}

## Figures

New in 0.13 (from the 0.13 release):

- [Hi-C's two readouts test each other](img/benchmark_cis_trans.pdf) ![Cis extremes that trans contacts confirm](img/benchmark_cis_trans.png)
  *Two readouts of one Hi-C library find centromeres that no repeat marks, with no annotation. **A** tomato, all 12 chromosomes: the
  contact across each bin (purple) and the contact with the other chromosomes' cis extremes (blue) meet on CENH3 on 10 of 12, and a
  randomized library (gray) shows neither. **B** chr6: the contact map, the other chromosomes' trans contacts, CENH3 ChIP. **C** why
  chromosome ends are not read (schematic). **D** *Ae. tauschii* read with the ends kept: both readouts meet at the ends, off CENH3.
  **E** every Hi-C library against chance. **F** *C. neoformans* in G1: 12 of 14.*

New in 0.12 (from the 0.12 release):

- [The contact crossing the centromere](img/benchmark_cis.pdf) ![The two flanks of a centromere touch](img/benchmark_cis.png)
  *In tomato the two flanks of every centromere touch about twice as often as DNA at the same distance elsewhere. **A** all 12 tomato
  chromosomes: the cis profile, its randomized-Hi-C partner and simulated maps, CENH3, the LTR signature, the satellites and both
  releases' calls. **B** one chromosome's contact map. **C** maps piled up on CENH3 in tomato and potato, against random positions.
  **D** potato's dips. **E** every genome with Hi-C and truth against chance. **F** the profile around every centromere.*
- [One locus per chromosome](img/benchmark_one_per_chromosome.pdf) ![A count tells centromeric features from others](img/benchmark_one_per_chromosome.png)
  *A centromere occurs once per chromosome, so the genome's N largest candidates cover nearly every chromosome when they are
  centromeres. **A-C** chicken, rice and tomato: every cluster of the centromeric families, the chance draws beside them. **D, E** GC troughs
  in *P. falciparum*, real and shuffled. **F** confident errors from the FASTA alone, 0.11 and 0.12. **G** every genome, both tests.*
{{FIGURES_MORE}}

From the 0.10.1 release ([what changed since](#13-what-changed-from-013) applies to these, with [0.11 to 0.12](#from-011-to-012)):

- [Mammals and Hi-C](img/benchmark_overview.pdf) ![benchmark overview](img/benchmark_overview.png)
  *9 T2T mammals with CENP-A truth; rice, barley, potato and horse EquCab3 with ChIP truth. **A** every chromosome
  with a mapped centromere: the call (box, colored by outcome) against CENP-A/CENH3 (black bar); arrows = barley calls
  moved by Hi-C. **B-F** one chromosome per kind of outcome. **G** every chromosome as one square, FASTA only vs all
  inputs. **H** each line alone vs the integrated call. **I** the Hi-C tests on real vs randomized Hi-C.*
- [Plants and algae](img/benchmark_plants.pdf) ![plant benchmark](img/benchmark_plants.png)
  *14 plant and algal species with CENH3 truth, FASTA only vs + `tes/` (+ Hi-C), and 7 genomes without a single
  centromere. **A** every chromosome run with `tes/`: CENH3, the calls, every element of the centromeric LTR lineage.
  **B, C** *Ae. tauschii* Chr1 placed by its CRM core. **D** the failure: tomato chr10. **E** every conflict between a
  CRM core and a satellite call. **F** *Chlorella*'s Osser insertions. **G** every chromosome as a square.*
- [Kinetochores from CpG methylation](img/benchmark_cpg.pdf) ![CpG kinetochores](img/benchmark_cpg.png)
  *Primates with DiMeLo-seq CENP-A, curated CDRs or CENP-A ChIP, and HG002. **A, B** a typical bonobo centromere. **C**
  widths of call, kinetochore and CENP-A domain. **D** every ape kinetochore on its CENP-A domain. **E** HG002 from ONT
  and PacBio HiFi. **F** the CpG test on real vs rotated methylation. **G** reads thinned to 1/16.*
- [Hi-C in small genomes](img/benchmark_hic.pdf) ![Hi-C finds centromeres where nuclei cluster them](img/benchmark_hic.png)
  *Where nuclei cluster their centromeres, Hi-C contacts alone find every centromere. **A-D** budding yeast: contacts
  between chromosomes, real and shuffled. **E** every chromosome of four genomes with their centromeres and cenfuse's
  calls. **F** clustered centromeres vs territories. **G** chromosomes whose centromere is found.*
- [The Rabl axis in trans](img/benchmark_hic2.pdf) ![Chromosome arms align by distance from the centromere](img/benchmark_hic2.png)
  *Where a nucleus holds its centromeres at one pole, arms of different chromosomes align by distance from the
  centromere. **A** *C. neoformans* in anaphase vs G1 vs randomized. **B** the Rabl polarity (schematic). **C, D** one
  and all 14 chromosomes. **E** every axis test of the benchmark vs chance. **F** what cenfuse does with it.*

## 1. Every genome

Right / chromosomes with truth (confidently wrong in brackets), for each input set benchmarked; "all inputs" names the
richest set where it differs from the columns before it. Chance = the right calls expected by luck (FASTA only).

{{TABLE_EVERY_GENOME}}

## 2. Satellite centromeres

### T2T mammals

| genome | truth | FASTA only | + genes + Hi-C | chromosomes passing the Hi-C test |
|---|---|---|---|---|
{{TABLE_MAMMALS}}

- **Totals** (9 genomes with CENP-A, 157 chromosomes): FASTA only 143 right, 1 confidently wrong, 13 `low` or no call;
  with genes and Hi-C 148, 4, 5. `high` calls were right {{MAMMAL_HIGH}}.
- **Without functional truth** (siamang, rhesus, sheep, goat, giraffe; 121 chromosomes): the calls agree with an
  independent satellite caller (the satellome locus) on 111 from the FASTA, 113 with genes and Hi-C.
- **Pig, FASTA only (8/17):** the satellite families split so that no family dominates genome-wide, and the honesty cap
  makes every call `low`; 5 of them sit 0.6-2.0 Mb from CENP-A at the acrocentric ends. With genes and Hi-C: 13/17.

| every miss (all inputs) | call | why |
|---|---|---|
| pig CM000817.6, CM000822.6, CM000828.6 | wrong, confident | the truth is weak ChIP from other pigs, with no centromeric array within 1 Mb of its maximum (5.4-5.6×); on all three the ChIP's second peak (3.0-4.4×) lies inside the call, and the published karyotype agrees with the call on two (chr11, chr17). Second peaks ≥ 3× also occur on 3 of the 13 right pig calls, so this is evidence, not proof |
| pig CM001155.3 | wrong, `low` | the same |
| Bornean orangutan CM054652.2 | wrong, `high` | CENP-A sits on 4 small alpha-satellite arrays (73 kb in cenfuse's search), 14 Mb from the 3.4-Mb array cluster cenfuse called; Hi-C points weakly to the CENP-A site but cannot move a `high` call. The truth may be haplotype-specific: the CenSat annotation marks the called array as active, and the CDR region was redrawn in 2025 |
| bonobo CM055489.2 | wrong, `low` | the true array is the reported competing cluster |
| horse CM128059.1 (chr11) | no call | satellite-free centromere; horse Hi-C fails its test, so nothing can place it |
| horse CM128050.1, CM128070.1 | wrong, `low` | |

### Plants and other genomes with satellite centromeres

- **Rice** (T2T-NIP): 12/12 `high` from the FASTA; Hi-C passes its test on all 12 and changes nothing.
- ***B. rapa*** (11 T2T accessions): 110/110 from the FASTA. CX1 A05, the one miss until 0.11 (`low`, 2.8 Mb off), is now
  right and `high` (since 0.11).
- **Soybean** X26, **sorghum** BTx623, ***Brachypodium*** Bd21 (CENH3 truth on 3 of 5 chromosomes), ***Arabidopsis***
  Col-CC: 20/20, 10/10, 3/3, 5/5 from the FASTA.
- **Maize** Mo17 (T2T): the 180-bp knob repeat is the most abundant satellite family (23.6 Mb, 20.6 Mb on chr8) and
  gets the role `interstitial` (knobs and other blocks); CentC is centromeric on all 10 chromosomes and no call holds
  knob repeat. 9/10 right (7 overlap CENH3). The
  miss, chr5 (`high`), has CENH3 on satellite-free DNA 2 Mb from the CentC array.
- **HG002** (diploid human T2T, 46 chromosomes): 46/46 calls on the CenSat active array (45 `high`).
- **Chicken** GGswu (40 chromosomes; truth on 33): 27/33 on the CENP-A windows (chance 6.4, P 4e-22), 2 more within 1
  Mb; every call `low`, because no satellite family dominates (chicken centromeres use chromosome-specific satellites,
  and chr5, chr27 and Z have none). Confidence is conservative here: no call is confidently wrong.
- **Mouse** (T2T, all 20 chromosomes acrocentric): both the pericentric major satellite (162 Mb) and the CENP-A
  minor satellite (16.5 Mb) mark the centromeric end. All 20 calls sit there, on the major-satellite block, within
  45 kb of the minor satellite (1/20 overlapping it), all `low`, with the minor satellite reported as the alternative
  on 20/20. No mappable DNA separates the two arrays, so Hi-C cannot choose (`hic_alt_ambiguous`).
- **Horse EquCab3** (older assembly, CENP-A ChIP on 5 chromosomes): 3/5 right, all `low`; chr11 and chr31 have no call.

### One locus per chromosome

A chromosome has one centromere; anything else falls roughly in proportion to chromosome length. So the N largest clusters of the
families called centromeric (N = the number of chromosomes) should cover nearly every chromosome if they are centromeres, and far
fewer if they are not ([design.md](design.md#one-locus-per-chromosome-a-test-of-the-satellite-line-satkaryotypepy)). The test reads
the FASTA alone and asks nothing about any one locus:

{{TABLE_ONE_PER_CHROMOSOME}}

- **Every genome with satellite centromeres passes but one**, down to chicken's chromosome-specific satellites (36 of 40 chromosomes
  covered, 17 by chance) and mouse, whose two satellites both sit at each centromeric end. The exception is the siamang (18 of 26,
  P 0.2): on 8 chromosomes its alpha-satellite arrays fall into two clusters more than 500 kb apart.
- **Tomato, potato, barley and *Ae. tauschii* fail** (P 0.24-0.91): their centromeres are satellite-free or sit beside the satellite.
  These four genomes held 16 of the FASTA-only confident errors of 0.11; their satellite calls are now `low` unless an LTR core agrees
  with them, and from the FASTA alone they make none.
- **The holocentric genomes fail too** (P 0.13-0.73), as they should. The siamang's failure costs it confidence, not calls: its
  satellite calls without genes and Hi-C fall from 13 `high` to `low` (it has no CENP-A data; its calls are compared with another
  satellite caller).
- **Passing says little about a single call.** *Chlamydomonas* (LINE-cluster centromeres) and moss (Bryco clusters) pass on their
  satellites; their satellite calls stay `low`, because no family dominates.
- **Over all 56 genomes**, the FASTA-only confident errors fall from 19 to 3; the 3 left (maize, sheep, Bornean orangutan) are in genomes
  that pass.
- **Polyploids** (no truth, 106 allotetraploids): homoeologous chromosomes put their satellite calls (0.11's, stored) at the same relative position
  more often in genomes that pass (mean 0.62, 70 genomes) than in those that fail (0.43, 35; Mann-Whitney P 0.001; `bench/poly_one.py`).

## 3. Retrotransposon centromeres (`tes/`)

`tes/` = LTRquest annotations, each made with `ltrquest --genome G --max-rounds 1 --threads T --detect-args
ltr-timeout=300 --run-sdust --reboundary` (the README's recipe adds `--proteins` and `--no-plots`). cenfuse learns the signature in
each genome: the most specific centromere-concentrated family or lineage seeds it, and each chromosome's core comes from the families
enriched in the other chromosomes' cores ([design.md](design.md#ltr-retrotransposons-tes-satsignaturepy)). Right / chromosomes
(confidently wrong):

{{TABLE_SIGNATURE}}

- **Any lineage can mark centromeres, and the family says it better than the lineage.** Moss's Bryco (58.8-fold), *Chlorella*'s Osser
  (78.9-fold) and a single Tekay family in tomato (32.7-fold; 308 of the lineage's 11,657 elements, 1 in 38) now decide, as CRM does in the
  grasses. With `tes/` alone: moss 16 → 26/26, barley 4 → 7/7, potato 6 → 10/11, *Chlorella* 8 → 9/10, tomato 1 → 6/12.
- ***Arabidopsis*:** Athila concentrates beside the centromeres (13.3-fold): its cores flank the CENH3 domain, 0.46-1.38 Mb away, on all
  5 chromosomes. CEN180 passes the one-per-chromosome test, so Athila only flags and the satellite calls stand (5/5); had Athila decided,
  it would have placed all 5 beside the centromere.
- ***Ae. tauschii*:** satellite arrays lie beside CRM cores that hold CENH3 on all 7 chromosomes; 5 calls move onto the cores; Chr5 and
  Chr7 stay at a core's edge, 1.3 and 3.1 Mb from CENH3, and are confidently wrong (`high` and `medium` with `tes/`; both `medium`
  with Hi-C too).
- **The head-to-head test still protects right satellite calls:** in *B. rapa* A09 and soybean chr1/chr6 the densest window of a broad
  CRM plateau lay off CENH3, which sat on the satellite.
- **Tomato** with `tes/` alone: the Tekay signature has significant cores on all 12 chromosomes, 7 on CENH3. Two of the 5 misses stay
  `medium` (CM175074.1, CM175078.1): nothing in the LTR line tells them from the right ones ([future.md](future.md#27-wrong-ltr-cores-that-look-right)).
  With Hi-C, the cis signal places 10 of 12 and none is confidently wrong (section 4).
- **Controls:** on 157 annotations with element positions shuffled within chromosomes, no family or lineage qualifies as a seed. The
  signature does not hang on its internal choices: 152-155 of 172 right cores over 1-4 recruitment rounds, P from 1e-4 to 1e-2 and
  3-10 elements.

**The annotation matters** (measured on 0.11, the CRM line alone). The decisive CRM test needs CRM to concentrate on most chromosomes.
Right calls with each annotator, and with LTRquest's elements thinned at random (TE line alone: right / cores):

| genome | annotation | TE line alone | calls right (confidently wrong) |
|---|---|---|---|
| rice | LTRquest | 8/8 | 12/12 (0) |
| rice | EDTA intact + TEsorter; DANTE_LTR | 4/4; 3/8 | 12/12 (0) |
| rice | LTRquest at 50 / 25 / 10% | 3/3, 1/1, 1/1 | 12/12 (0) |
| potato | LTRquest | 6/7 | 6/11 (2) |
| potato | EDTA + TEsorter; DANTE_LTR; LTRquest at 50% | 2/2; 1/1; 3/3 | 2/11 (2) |
| maize | LTRquest, and at 50 / 25 / 10% | 10/10, 10/10, 10/10, 8/8 | 9/10 (1) |
| moss | LTRquest; at 50% | 25/25; 22/22 | 16/26 (0) |

Potato's 214 LTRquest CRM elements give exactly a majority of chromosomes; any less sensitive annotation loses the
gain from `tes/`. Every CRM core short of a majority still lay on CENH3 (24/24 in rice, potato, maize and
*Brachypodium*).

### Held out from 0.14: an allotetraploid and a wild potato

{{HELD14}}

### Held out from 0.14: thirteen tomatoes

{{HELD14_TOMATO}}

## 4. Hi-C

### The test, and Hi-C on its own

Before Hi-C may act on a chromosome, its peaks must agree with the other chromosomes' calls on most of them, beyond
chance ([design.md](design.md#hi-c-cenfusehic-and-hicfusepy)). Each Hi-C readout alone, on the genomes with truth:

| genomes (chromosomes) | satellite line | prior-free Hi-C | anchored Hi-C | integrated |
|---|---|---|---|---|
| 9 T2T mammals (157) | 148 | 32 | 44 | 148 |
| rice (12) | 12 | 0 | 7 | 12 |
| barley (7) | 0 | 6 | 0 | 6 |
| potato (11) | 2 | 1 | 0 | 2 |
| horse EquCab3 (5) | 3 | 3 | 3 | 3 |

Every line scored within 4 Hi-C bins of the truth (Hi-C's resolution); the mammals' satellite line is the one with
genes (FASTA alone: 143).

- In mammals the anchored peak sits beside the unmappable centromeric array, often more than 4 bins from CENP-A in bp;
  the satellite line decides, and Hi-C confirms (the test passes in all great apes, crab-eating macaque, CHM13 and
  siamang).
- Hi-C's prior-free peaks cannot be trusted where they can be tested (the ones the engine rates `high` lie within 4 bins
  of CENP-A on 14 of 63 T2T mammal chromosomes; potato's strongest mutual cluster is not its centromeres), so they
  never act unchecked.
- **Randomized Hi-C** (partners re-paired at random with coverage kept, or ends moved): in 57 runs (14 mammals
  randomized both ways, *P. annua*, 5 small genomes, 13 held-out genomes, and 10 genomes with `tes/`), no test passed, the support
  tier was never reached, and every call equals the call made without Hi-C (56 runs; *P. annua*'s has no run without Hi-C to compare).

### Hi-C and the LTR cores

Where the Hi-C test fails on most chromosomes because the sequence calls are wrong anchors, and `tes/` gives a signature, Hi-C is tested
against the LTR cores instead ([design.md](design.md#hi-c-and-the-ltr-line-test-each-other-hicfusepy)).

- **Better calls make better anchors.** With the LTR signature, the ordinary test passes in barley on 7 of 7 chromosomes (none
  without `tes/`) and in moss on 26 of 26: barley's calls are right on 7 of 7, all `high`, and moss's 26 of 26.
- **Hi-C anchored on the cores** passes in *Brachypodium* (3 of 5, P 0.0017) and *Arabidopsis* (P 0.0028), whose calls were right
  already: it flags 5 calls that were `high` already (`te_hic_pair`; 4 with truth, all right). It fails in potato (2 of 12, P 0.12), tomato (P 0.46) and sorghum.
- **Ends do not count.** In held-out *R. breviuscula* (holocentric), the Athila signature's cores lie at the ends of 4 of 5 chromosomes,
  and its nuclei cluster their chromosome ends (2.6-fold). Under the frozen 0.12 rules Hi-C anchored on those cores passed (3/5,
  P 5e-4) and moved 3 calls onto them at `medium`, in a genome with no single centromere. The test now counts a peak only where it is
  clean, the rule the Hi-C test already applies to `low` references: *R. breviuscula* fails (1/5, P 0.18) and its calls equal those
  without `tes/`; *Brachypodium* and *Arabidopsis* pass as before; of the 210 calls in the 20 runs where the test runs, no other
  changes ([section 12](#12-how-the-benchmark-was-run)).
- **The cis signal** passes in tomato, potato and *Ae. tauschii* (next section).

### The contact crossing the centromere

The contact between the 4 bins left and the 4 bins right of each bin, as observed over expected at those distances, has its extreme at
the centromere in many plants and fungi, and in no mammal. Which extreme differs: in tomato the two flanks of each centromere touch about
twice as often as DNA at the same distance elsewhere, in potato much less. Each genome with Hi-C and truth (the extreme within 4 bins of
truth; chance from contact maps simulated with the same coverage and distance decay; not corrected for the two signs), and the test
cenfuse runs against the LTR cores where `tes/` is given and the Hi-C test fails:

{{TABLE_CIS}}

- **Tomato:** the peak lies on CENH3 on all 12 chromosomes, where every other line fails on some. Against the LTR cores it passes (7 of
  12, P 1e-4), so it moves the `low` calls and lowers the confident calls it contradicts: 10 of 12 right, none confidently wrong (0.11: 1
  right, 7 confidently wrong).
- **Potato:** the dip, at a flank of the centromere, passes against the cores (7 of 12, P 1.4e-5): 9 of 11 right, 1 confidently wrong.
- **What the calls did** (all arms with truth): an extreme at the call raised it a level on 15 calls, all right; it moved 8 `low` calls
  (6 right after, 2 wrong); a confident call it contradicts fell a level on 9 (4 right `low`, 2 wrong `low`, 3 still confidently wrong).
- **Controls.** On randomized Hi-C (partners re-paired at random, every bin's coverage kept) neither test passes in any of 10 genomes. The
  closest is moss, whose re-paired map has its dip at the LTR cores on 11 of 26 chromosomes (P 0.011 against the level of 0.01); among 20
  tests (10 maps, two signs) one value this small is expected about one time in five. The first version of both tests did pass on such maps
  (the plain mean contact with the zone's share as chance; partner windows moved to the best score near each core), and both were changed
  before the release ([section 12](#12-how-the-benchmark-was-run)).

### Hi-C's two readouts test each other (0.13)

Without `tes/`, nothing chose the sign of the cis signal or tested it. Since 0.13 the library's own contacts between chromosomes do: each
chromosome's cis extreme outside its terminal eighths anchors the anchored test, and the cis signal acts where each chromosome's trans peak,
made from the other chromosomes' extremes, lies at its own extreme on most chromosomes ([design.md](design.md#the-cis-and-trans-readouts-of-hi-c-test-each-other-hicfusecis_trans_line)).
Every run where it ran (the anchored test failed on most chromosomes, no LTR cores):

{{TABLE_CIS_TRANS}}

{{CIS_TRANS_BULLETS}}

### Which end

A call in the terminal eighth of a chromosome sits where Hi-C's peaks cannot count, since chromosome ends cluster in 3D whatever the
centromere. But the anchored score is higher at the end holding the call than at the other end on most such chromosomes in goat, sheep and
horse:

{{TABLE_END}}

Of the 56 calls it raised (all `medium` → `high`), the 55 scored are right 54 times; the other is the sheep call that already disagreed with the other
satellite caller. In the T2T horse, 16 of them have CENP-A truth: 16 right. In mouse every call has a competing cluster at the same end, so
nothing is raised; chicken does not pass (12 of 19).

### The Rabl axis: barley, small genomes and the trans arm alignment

In a Rabl nucleus the arms run from the centromere pole to the telomere pole, so each chromosome's own map folds at
its centromere (the arm fold), and arms of different chromosomes align by distance from their centromeres (the trans
arm alignment). Either can test Hi-C without the sequence calls.

- **Barley:** the satellite calls sit toward the chromosome ends (centromeric index ≤ 0.14; CENH3, on Cereba
  retrotransposons, at 0.35-0.48). The arm fold agrees with Hi-C's prior-free peaks on 6/7 chromosomes (P 3e-5), so 6
  calls move, all onto CENH3 (0/7 → 6/7). On 5H the prior-free peak sits on CENH3 but outside the fold core: the call
  is not moved and stays `low` (from the FASTA it is `low` already, as barley's satellites fail the one-per-chromosome test).
- Chromosomes whose prior-free peak lies in the alignment core, in every 0.12 run with Hi-C ({{RABL_RUNS}} genome runs; 44 libraries
  thinned to 2-25% passed none in 0.11 and were not rerun):

{{TABLE_RABL}}

- One assembly in two nuclear states: *C. neoformans* in anaphase passes and gives 10/14 right calls; in G1, when its
  kinetochores are unclustered (Polisetty et al. 2025), nothing passes and nothing is called. Placed on the known
  centromeres, the G1 nucleus still aligns its arms on 8/14 chromosomes (P 0.0014): the alignment weakens without the
  cluster rather than vanishing.

### Point, regional and AT-rich centromeres: AT troughs and Hi-C

With no centromeric repeat, the satellite line gives no call. Two lines then remain: AT troughs from the FASTA, and Hi-C. Right =
within 1 Hi-C bin (250 kb without Hi-C); [high/medium/low]:

| genome | centromeres | Hi-C | bin | FASTA only | + Hi-C |
|---|---|---|---|---|---|
{{TABLE_SMALL}}

**AT troughs.** The genome's N deepest GC troughs must cover more chromosomes than chance, at one of the window sizes 1-500 kb:

{{TABLE_AT}}

- From the FASTA alone: *P. falciparum* 13/14, *N. crassa* 7/7 and *M. sympodialis* 8/8 right (all `low`; 0.11: 0 of 29).
- With Hi-C, the prior-free peaks test the troughs: where both agree, the call takes the trough's interval and rises a level (40 calls,
  all right).
- It fails in *Z. tritici*, *C. tropicalis*, *C. neoformans* and *V. dahliae*, and cannot say anything with 3 chromosomes (*S. pombe*);
  the holocentric *Bombyx* fails. GC blocks of 100 kb shuffled among chromosomes pass on 16 of 360 shuffles of 18 genomes (4%), the
  test's level being 5%.

**Hi-C.** Yeast and *Plasmodium*: the trans arm alignment passes and the arm fold leans the same way; *Neurospora*: the arm fold itself
passes. Where neither passes, Hi-C's clean prior-free peaks give `low` calls (*Z. tritici*, *C. neoformans*), raised to `medium` where an AT
  trough agrees (`at_pair`: 4 of *C. neoformans*' 10 calls in anaphase). *S. pombe* chr2's
prior-free peak is only `medium`; with the AT troughs it is now right (3/3). *C. tropicalis* (3C-seq) passes nothing. Randomized Hi-C
on each of these genomes: no test passes and no Hi-C-alone call is made.

### Held out: crop and mouse genomes newly given Hi-C

| genome | Hi-C (pairs kept) | FASTA only | + Hi-C | what the 3D data showed |
|---|---|---|---|---|
| *Brassica rapa* CC1 | 8.2 M | 10/10 [9/0/1] | 10/10 [9/0/1] | centromere_clustering; test passes on 10 |
| soybean | 10.7 M | 20/20 [18/0/2] | 20/20 [18/1/1] | centromere_clustering; test passes on 20 |
| sorghum | 3.9 M | 10/10 [10/0/0] | 10/10 [10/0/0] | not_detected; test passes on 0 |
| *Brachypodium distachyon* | 14.3 M | 3/3 [3/0/0] | 3/3 [3/0/0] | not_detected; test passes on 0 |
| *Aegilops tauschii* | 0.8 M | 1/7 [0/0/1] | 1/7 [0/0/1] | not_detected; test passes on 0 |
| mouse (T2T) | 15.7 M | 1/20 [0/0/1] | 1/20 [0/0/1] | weak_centromere_clustering; test passes on 0 |

No right call became wrong and no confident error was added in any of them. Three libraries are shallower than
downloaded (sorghum, *Ae. tauschii*, mouse: a mate file was damaged in transfer and only its intact part was mapped).

### Held out from the 0.12 rules

The 0.12 rules were frozen before any of these arms was scored: LTRquest annotations made for this release on genomes whose `tes/`
arms had never been run, and two genomes new to the benchmark (grape PN40024 and cucumber 9930, CENH3 ChIP-seq from their authors,
domains by the frozen rule).

{{TABLE_HELD12}}

{{HELD12_BULLETS}}

### Weak Hi-C: the support tier

Where Hi-C finds the centromeres beyond chance but not on most chromosomes (territorial nuclei: goat, rhesus, horse
EquCab3, moss without `tes/`), it may only support the calls it agrees with. At full depth the tier supported 23 calls in three genomes, all of them on CENH3 (moss without `tes/`, 10 calls raised to `medium`) or on the satellome locus (goat 8 raised to `high`, rhesus 5 already `high`). Randomized Hi-C never reached the tier.

### Depth, MAPQ, other pipelines' files, replicates

- **Depth.** Pairs thinned at random (nested; 2-50%). Chromosomes on which the Hi-C test passes:

| genome | pairs | 100% | 50% | 25% | 10% | 5% | 2% |
|---|---|---|---|---|---|---|---|
| CHM13 | 71 M | 24 | 24 | 24 | 24 | 24 | 0 |
| chimpanzee | 74 M | 25 | 25 | 25 | 25 | 25 | 0 |
| gorilla | 55 M | 25 | 25 | 25 | 0 | 0 | 0 |
| Sumatran; Bornean orangutan | 71 M; 67 M | 25; 25 | 0; 0 | 0; 0 | 0; 0 | 0; 0 | 0; 0 |
| crab-eating macaque | 70 M | 22 | 22 | 22 | 22 | 22 | 22 |
| pig | 54 M | 0 | 0 | 20 | 20 | 0 | 0 |
| horse (T2T); potato | 64 M; 66 M | 0 | 0 | 0 | 0 | 0 | 0 |
| rice | 53 M | 12 | 12 | 0 | 0 | 0 | 0 |
| barley (arm fold, or the trans alignment) | 9.2 M | pass | pass (alignment) | fail | fail | fail | fail |

  Tests flip at their threshold rather than fading: both orangutans fail at half depth (34-35 M pairs). Where a test
  fails, the calls are exactly those made without Hi-C (32 of 32 such runs), so thinner Hi-C never made a call worse
  than the FASTA alone; barley falls back to its FASTA-only calls (3 confidently wrong) below half depth. Give at
  least ~70 M pairs to a 3-Gb genome.
- **MAPQ floor.** The same reads kept at pairtools' default floor (MAPQ ≥ 1) instead of 30 changed no test and made no
  call confidently wrong in CHM13, orangutan, rice, moss or barley; near-identical repeat copies map at MAPQ 0, which
  pairtools drops at any floor.
- **Other pipelines.** Rice's 53 M pairs rewritten as HiC-Pro `allValidPairs` and Juicer `merged_nodups`, and its
  `.hic` converted with `hictk`: identical calls and tables.
- **Replicates.** Rice's pairs split into two files given as replicates of one plant: every table identical to the
  single file.
- **A territorial plant nucleus.** *Arabidopsis* Col-0 leaf Hi-C (15.4 M pairs): the five calls stay right and `high`;
  the prior-free peaks lie 1-4 Mb off CENH3 (chromocenter neighbors), as expected where centromeres do not dominate
  trans contacts.

## 5. Kinetochores from CpG methylation

Inputs: public ONT CpG tracks of the assembled individuals (apes: modkit 5mC bigWigs; CHM13 nanopolish; crab-eating
macaque with read counts; rhesus R10 modkit bedMethyl), with genes and Hi-C. Truth: DiMeLo-seq CENP-A (5 apes), the
curated CDRs (the same ONT data: an implementation check), and CENP-A ChIP of another crab-eating macaque. The line was
designed on these primate tracks (not blind); rhesus, *Arabidopsis* and HG002 were held out.

| test | result |
|---|---|
| the CpG test passes on real tracks | 6 of 7 primates with truth (CHM13, bonobo, gorilla, both orangutans, crab-eating macaque) and rhesus; not chimpanzee (positional dips at 1 of 25 calls) or siamang |
| it fails on misplaced methylation | rotated tracks 0/8 genomes; tracks shifted by half a chromosome rejected by the assembly check on 6/8 |
| kinetochore on DiMeLo-seq CENP-A (5 apes) | 35 of 39 within 50 kb (29 overlap); CENP-A signal higher inside the kinetochore than over the call on 39/39 |
| kinetochore on the curated CDR | 62/63 |
| crab-eating macaque kinetochore inside the ChIP domain | 12/15 (the other 3 lie 0.8-3.4 Mb from the domain of another animal's ChIP, inside the same call) |
| resolution | kinetochores a median 55 kb wide, calls 5.1 Mb (133 kinetochores, 6 primates) |
| calls changed | 0 |

- **HG002** (diploid human; no CENP-A data): CpG from ONT and from PacBio HiFi. The test passes on 46/46 with both; kinetochores on 40 (ONT, median 38 kb) and 43 (PacBio, median 25 kb), all inside the active array; where both place one (39 chromosomes), they lie within 50 kb on 31. PacBio tracks enter as bedGraphs, with no code change.
- **Low coverage.** Calls are pooled per 5-kb bin, so rhesus thinned to 1/16 of its reads (about one call per CpG
  strand) keeps 21/21 kinetochores within 50 kb of full depth; crab-eating macaque the same.
- **Wrong calls get kinetochores too**: Bornean orangutan CM054652.2 gets one inside its inactive array. The line
  refines a call; it cannot check it.
- ***Arabidopsis*:** the CpG test is right to fail (positional dips on 1/5; the strongest dips lie off CENH3 on 4/5).

## 6. CHG methylation (plants)

- ***Arabidopsis*** Col-CC, same-genotype ONT (~56×): calls identical with and without methylation; the CHG line finds
  no signature, because its expected domain size comes from genome size (315 kb) while CEN180 arrays span 2.2-3.5 Mb.
- ***Poa annua*** (all four inputs; five ONT plants, all other plants): methylation does not decide, since no sample is
  from the assembled plant.

## 7. Controls with no single centromere

| genome | why it is hard | calls (FASTA only; with 3D data) |
|---|---|---|
| *C. elegans* CGC1 | holocentric, no centromeric satellite | 6 `low`; with N2 Hi-C, the same |
| *C. nigoni* | holocentric | 6 `low`; with Hi-C, the same |
| *Bombyx mori* | holocentric; Lepidoptera lost CenH3 | 7 `low`, 22 none; with Arima Hi-C, the same |
| *Rhynchospora breviuscula* | holocentric; Tyba arrays along every chromosome | 5 `low`; with Omni-C, the same |

0 confident calls in 46 chromosomes, from the FASTA and with 3D data. Used without its test, Hi-C's own peaks would
have made 21 confident calls on these chromosomes (*C. elegans* 5, *C. nigoni* 5, *Bombyx* 10, *Rhynchospora* 1); every
Hi-C test failed, so nothing changed.

## 8. Several plants (*Poa annua*)

*P. annua* (14 chromosomes) and *P. infirma* (7), all four inputs; truth from earlier analyses (no ChIP).

- **Pa1A, a miss of 0.12.** Earlier analyses (no ChIP) place its functional centromere at 163.0-169.4 Mb, where Hi-C peaks
  (168.0 Mb), 21 Mb from its largest satellite; 0.11 called it there (`high`). 0.12 calls the satellite site instead (146.1-146.8 Mb,
  `medium`, `contested`, `hic_offset`): the CRM families learned from the other 13 chromosomes concentrate around the satellite
  (145.1-148.3 Mb, 30 of their 146 elements on Pa1A), so two lines agree there and Hi-C can only contest the call. The family
  signature marks where centromeric families accumulated; a centromere that has since moved keeps neither its satellite nor those
  families. Insertion age, an independent line, agrees with the earlier analyses: CRM copies at the Hi-C site are young (median age
  about 0, 89% under 100,000 years) and those at the satellite site old (325,000 years); within one family, 0 vs 157,000 years (P 2e-7),
  as recently active CRM colonize the present centromeres in maize ([Wolfgruber et al. 2009](https://pmc.ncbi.nlm.nih.gov/articles/PMC2776974/))
  ([future.md](future.md#31-a-centromere-that-left-its-satellite-and-its-ltr-families)).
- The Hi-C test passes on 14/14 and 7/7. With randomized *P. annua* Hi-C no test passes (0/14) and the calls rest on the sequence
  and LTR lines: Pa1A is then `high` at the satellite, so real Hi-C is what lowers it and marks the conflict.
- The five methylation plants are recognized as other plants (5/5), and the five known missing segments are found
  (the two Pa2A losses in 3 pieces each). On 300 pure-noise genomes the between-plant test gave 5.0% of null tests at
  P < 0.05 and no `shifted` verdict.

## 9. Other features

- **Telomeres.** Arrays at 687 of the 688 ends of 14 T2T mammals that hold one (≥ 500 bp of TTAGGG in the terminal 10
  kb). Short telomeres count when they reach the end: *Neurospora* 13/14 ends (157-207 bp), *P. falciparum* 27/28,
  yeast's irregular TG1-3 repeats found as GGGTGT. Misses: ends without telomere repeat in the assembly (chicken 22),
  TTAGG broken by retrotransposons (*Bombyx* 15). Units of 13-29 bp are tried where no 5-12-bp unit qualifies, and only if
  they reach most ends: *C. tropicalis*' 23-bp unit, at 12 of 14 ends, and no spurious unit in any other genome. A
  stretch that reaches the end counts down to a tenth of unit copies (the lowest seen at a real end is 0.19): 1,539 -> 1,557
  ends over the benchmark, none lost, interstitial arrays unchanged.
- **Telomere and 45S recall, every genome** (FASTA-only arms; since 0.14, [future.md #56](future.md#56-telomeres-and-rdna-arrays-the-benchmark-does-not-count)).
  Ends with a terminal telomere array, of two per chromosome; complete chromosomes (an array at both ends and no gap); 45S arrays and
  the chromosomes holding them. Where no unit is found, the assembly's ends lack telomere arrays: potato's DM assembly holds one
  (378 bp of TTTAGGG at the start of chr11) in 24 ends, and a unit is learned only from arrays at two ends or more; the
  *S. pombe* reference ends hold at most 8 copies of its GGTTAC-like repeat in 10 kb; *C. nigoni* one end 40 copies of TTAGGC
  (240 bp, below the 300-bp array floor). So the assemblies hold almost no telomere there; the one potato array goes unreported
  because a unit seen at a single end is not learned (#56):

{{TELOMERE_RECALL}}
- **Interstitial telomeric arrays.** CHM13: 16 arrays ≥ 300 bp, one at the 2q13 fusion site (chr2 114.03 Mb).
  Chicken: 50 arrays, 22 of them on or within 100 kb of a CENP-A window (expected 1.4; binomial P 4e-21), holding
  63% of all interstitial TTAGGG: the telomere-repeat centromere cores described by Huang et al. 2023, genome-wide.
- **45S rDNA (NORs).** Every checkable array lies on a known NOR chromosome: human chr13/14/15/21/22, *Arabidopsis*
  chr2/chr4, chicken chr16, *C. elegans* chrI, pig chr8/chr10, horse chr1/28/31, barley 1H/2H/5H/6H/7H, mouse
  chr12/15/16/18/19. Organelle rRNA holds the 18S motif too, so a hit counts only if the bases after it are nuclear (next
  base and the 19 after it; the two groups separate: 110 hits with at most 5 of 19 differing, 1,100 with at least 11, 4
  between). That removes 11 of 272 arrays (rice chr4, *B. rapa* A02/A06, Sumatran orangutan). The 5.8S genes confirm it:
  all 74 kept arrays sampled hold them, and removed arrays at most 3.
- **5S rDNA.** Every array of 100 copies or more lies on a chromosome known to carry 5S rDNA: CHM13 1q42 (126 copies of 2.24 kb), the HSA1
  homolog in every ape and macaque, mouse chr8, chicken chr9, *Arabidopsis* chr3/4/5, rice chr11, maize 2L, potato
  chr01, barley 2H/4H/7H; none in *S. pombe*, whose 36 5S genes are dispersed. Arrays of 3-24 copies also lie on *Arabidopsis* chr1 and
  barley 3H and 5H.
- **Morphology.** Acrocentric vs bi-armed agrees with the literature karyotype on 117 of 118 right calls (mouse 20/20,
  goat 30/30, horse 29/29, pig 13/13, sheep 25/26). By sequence, CHM13's acrocentric chr21 and chr22 come out
  submetacentric/subtelocentric: their T2T short arms, rDNA included, are 11-15 Mb.
- **Completeness.** Assembly gaps equal an independent scan on EquCab3, barley and potato (4,084 gaps).
- **GFF3.** Every features file passes `gt gff3validator -typecheck so`.

## 10. Speed, memory and reproducibility

- **The search since 0.11.** In long models (500 bp or more) low-complexity tracts are neutral ([design.md](design.md#engineering-notes)). Speed: the plain
  and the neutral model on the same three 2-Mb chunks, in one process, in thread CPU, since wall times on a shared server cannot be
  compared. Results: full reruns against a replay of the same code on the stored search, so only the search differs. Of 143
  chromosomes with truth in these genomes, none changed outcome and none changed confidence; the benchmark's other genomes were not
  searched again.

<!-- search table -->
| genome | models neutral | CPU per 3 chunks, plain -> neutral | gain | whole copies, stored | lost / gained | calls moved (at least 10 kb) | roles changed |
|---|---|---|---|---|---|---|---|
| *Lespedeza potaninii* (4x) | 13 of 21 | 3,401 -> 689 s | 4.93x | 26,716 | 9 / 407 | 12 (6) | 0 |
| *Chlamydomonas* | 5 of 20 | 273 -> 83 s | 3.28x | 20,971 | 4 / 203 | 6 (4) | 1 |
| *Panax japonicus* (4x) | 25 of 35 | 7,712 -> 3,868 s | 1.99x | 21,240 | 331 / 551 | 12 (3) | 5 |
| human CHM13 | 1 of 6 | 44 -> 23 s | 1.91x | 2,089,849 | 0 / 0 | 0 (0) | 0 |
| *Hibiscus tridactylites* (4x) | 13 of 23 | 493 -> 277 s | 1.78x | 248,392 | 15 / 164 | 0 (0) | 2 |
| moss | 3 of 10 | 329 -> 236 s | 1.40x | 16,365 | 7 / 6 | 0 (0) | 0 |
| *Arabidopsis* Col-CC | 2 of 10 | 320 -> 230 s | 1.39x | 72,186 | 3 / 148 | 0 (0) | 1 |
| potato | 4 of 24 | 273 -> 208 s | 1.31x | 48,334 | 9 / 20 | 0 (0) | 0 |
| *Cuscuta chinensis* (4x) | 15 of 30 | 2,199 -> 1,918 s | 1.15x | 115,817 | 2951 / 7322 | 20 (11) | 6 |
| chicken | 4 of 25 | 961 -> 906 s | 1.06x | 320,355 | 2 / 66 | 26 (1) | 1 |
| sheep | 2 of 16 | 74 -> 71 s | 1.05x | 206,701 | 122 / 158 | 0 (0) | 0 |
| cotton (4x) | 4 of 29 | 325 -> 314 s | 1.03x | 155,846 | 934 / 902 | 0 (0) | 2 |
| *Glycine dolichocarpa* (4x) | 2 of 13 | 791 -> 777 s | 1.02x | 5,638 | 0 / 1 | 0 (0) | 0 |
| *Chlorella sorokiniana* | 0 of 1 | 2 -> 2 s | 1.00x | 742 | 0 / 0 | 0 (0) | 0 |

Over 14 genomes: 3,349,152 whole copies, 4387 lost (0.131%; a copy split in two counts as lost), 9,948 gained; 1 confidence levels and 19 main families changed; 76 calls moved, 25 of them by at least 10 kb, of which 18 became narrower: hits of a tract alone no longer extend them.
<!-- end search table -->

  Whole copies "lost" are counted per family, so they include copies that only changed family. A neutral model scores less, so
  where two close families overlap, the one without a tract can take the other's copies: in cotton a 512-bp family's 906 whole copies
  now fall to its 481-bp relative, which also loses an rDNA tag (neither is centromeric; no call moved).

  Where calls stand on overlapping families with little between them, the new search can change which one leads; neither genome
  below has truth, so whether that is better is unknown. *Cuscuta chinensis*: 12 calls change main family while staying in place, 6
  families change role, and another family dominates. *Panax*: no family dominates and calls hold a median of 54 kb of satellite, so
  sites nearly tie; 7 calls change main family and 2 go to another site. On sequence 6 the lead passes between clusters of 33 and 35
  kb (31.8 -> 19.3 Mb, `high` -> `low`, with the old site now the reported competing cluster; neither overlaps the Tekay core); on
  sequence 8 it passes between clusters of 25 and 27 kb, both `low`.

- **Runtime** (measured before 0.11's search). The satellite search takes most of the time and grows with satellite content: T2T mammals (2.6-3.6 Gb)
  1-10 h at 12 threads (median 4.3 h), mouse (~6% satellite) 9 h at 16, maize 50 min at 16, *Ae. tauschii* (4.1 Gb)
  2 h at 16, *Arabidopsis* 15 min at 8. Hi-C and the integration take minutes.
- **Memory.** Peak of a whole run: rhesus 74 GB, horse 59 GB, giraffe 18 GB at 12 threads; chicken 42 GB at 16;
  *Arabidopsis* 5 GB at 8.
- **Reproducibility.** Two fresh moss runs gave 307 identical files apart from header dates; Col-CC reruns at 8 and 16
  threads are identical.

## 11. Every confident miss (all inputs)

{{TABLE_CONF_MISSES}}

{{CONF_MISSES_NOTE}}

## 12. How the benchmark was run

- **Pre-registration.** Each benchmark's criteria were written before its runs: plants 2026-10-16/19, mammals
  2026-10-21, then one plan per release (kinetochores, hard genomes, new lineages, the LTR line, a from-scratch rerun
  with new truths, Hi-C in small genomes, the Rabl axis and held-out genomes), each with dated amendments. The plans,
  run records and reports are kept beside the runs.
- **Criteria that failed as written** (reported, not reworded):
  - Hi-C peaks within 4 bins of CENP-A wherever the test passes: in four apes the peaks sit beside the unmappable
    array, more than 4 bins from CENP-A.
  - `high` calls right more often than `medium`/`low` at P ≤ 0.01: too few `high` errors to test.
  - Methylation supporting ≥ 10/13 *P. annua* calls: 0/13, since every sample is another plant.
  - The crab-eating macaque kinetochore inside the ChIP domain on ≥ 90%: 12/15; the other 3 lie inside the same call,
    0.8-3.4 Mb from the domain of another animal's ChIP.
  - The *Arabidopsis* CHG line supporting ≥ 3 calls: 0/5.
  - ≥ 15/20 mouse calls overlapping the minor satellite: 1/20.
  - ≥ 90% of the calls the Hi-C support tier raised being right: 88%.
  - The Rabl axis passing in the grasses held out: *Brachypodium* leaf Hi-C and *Ae. tauschii* (0.84 M pairs) passed
    nothing.
- **The 0.12 criteria, as written before its runs** (2026-11-06), with the results:

{{CRITERIA12}}

- **The 0.13 checks, as written before its runs** (PLAN.md, 2026-11-07), with the results:

{{CRITERIA13}}

- **The 0.14 checks, as written before its runs** (PLAN.md, 2026-11-07, with its amendment before the #45 measurement), with the
  results:

{{CRITERIA14}}

- **Rules changed after a result was seen.** Each change has a reason independent of the case that prompted it, and
  was rechecked on cases not used to make it:

| change | prompted by | rechecked on |
|---|---|---|
| Hi-C peaks in end bins never oppose a call; Hi-C opposes only beyond 4 bins | rice chr4; CHM13 | potato, barley, horse, gorilla, mammals |
| satellite rules (tandem test, satellite kinds, acrocentric test, competing clusters); Hi-C in mappable bins | the first mammal run | the plant benchmark, a 185-genome batch |
| the arm fold; fold centers only at centromeric index ≥ 12.5% | barley | rice, all mammals, randomized controls |
| the CpG line, and its positional (rotation) null | designed on the primate tracks; the rotated gorilla control passed the first null | the other rotated controls, rhesus, *Arabidopsis* |
| `low` calls as Hi-C references when no call is confident (clean peaks only) | siamang without genes | randomized Hi-C, chicken, EquCab3 |
| the acrocentric test judged per satellite kind; a second unrelated acrocentric satellite competes | T2T mouse | 72 earlier satellite runs, 9 other genomes |
| a dispersed-repeat family gets role `dispersed` | chromosome-scale loci in bonobo | 220-genome replay: calls identical |
| TE core floor of 3 independent insertion sites | *Chlorella* | 12 crop `tes/` arms, moss, *Chlamydomonas*, *Arabidopsis*: no call changed |
| a CRM core replaces a confident satellite call only if CRM marks it significantly more | *B. rapa* A09 | held-out *B. rapa*, soybean, sorghum, *Ae. tauschii*, tomato |
| the same test for a satellite call beside the core | *Ae. tauschii* Chr2/Chr4 | 21 `tes/` arms (only those 2 changed) |
| families only from seed clusters with a tandem array (speed) | T2T sheep's search time | siamang, barley, Bornean orangutan: calls identical |
| the support tier counts only significant peaks | randomized Sumatran orangutan reached the first form | 28 randomized controls |
| a peak is `telomeric` only when it cannot be told from the end | moss centromeres near the ends | the development replay (34 runs) and 13 held-out genomes |
| the trans arm alignment acts where the fold leans the same way | *C. elegans* nearly passed; barley at half depth could not be rescued | holocentric folds, held-out genomes |
| no Hi-C support where no mappable bin separates a call from its competing cluster | held-out mouse (7 confident calls on the major satellite) | development replay: costs 3 right ape calls (`medium` → `low`) |
| a long-unit array is no call without Hi-C, `tes/` or methylation behind it | 24 FASTA-only calls: 2 right, 16 wrong, 6 without truth | all 204 arms (section 13): 28 wrong `low` gone, 8 right `low` gone, 1 wrong replaced by a right Hi-C call |
| calls on chains of partial copies at most `medium`, `low` without an overlapping CRM core | 459 of 2,355 polyploid calls, in a replay; no benchmark genome has such a family dominant | benchmark replay: 17 calls flagged, 1 changed level, none moved. No truth tests the cap itself |
| the period test of `classify.related` works both ways | maize's knob family counting CentC's chromosomes, kept out only by a size floor | replay of 3,459 chromosomes, and all 204 arms: only *B. rapa* CX1 and 2 polyploids move. The knob stays out at floors of 50-100 kb |
| 45S motif hits must have nuclear context | rice chr4, a plastid insertion called a NOR | 5.8S genes at kept and removed arrays; a one-base rule was rejected, since it deleted a real barley 7H array |
| telomere units of 13-29 bp, loose stretches at ends | *C. tropicalis*; soybean chr3 | T2T assemblies: no end lost, no interstitial array gained |
| low-complexity tracts of long models neutral in the search | runaway search cost in *Lespedeza*, *Panax* | full reruns on 14 genomes (section 10) |
| 0.12: the one-per-chromosome test caps satellite calls; for a non-CRM marker it replaced a planned rule (more chromosomes marked than by the dominant satellite) | the FASTA-only confident errors of tomato, potato, barley, *Ae. tauschii* | every genome (section 2), held-out grape and *B. rapa*, 106 polyploids (homoeolog agreement) |
| 0.12: a core must be 10-fold enriched on its chromosome to overrule a confident call or reject a family | *Oryza minuta* (polyploid, no truth): homoeolog agreement 12 → 8 of 12 with weak CRM cores acting | polyploids; the `tes/` arms |
| 0.12: a thin core of a decisive signature takes a chromosome with no call or a failed satellite (`weak_core`, `low`) | potato chr04, chr09, chr12 | the other `tes/` arms: 3 such calls, all right |
| 0.12: a core too weak to replace a call does not reject a satellite family either | a code review before the release (consistency with the rule for replacing calls) | the release: tomato with `tes/` only loses one right `low` call (7 → 6 of 12) |
| 0.12: Hi-C partner windows stay at the LTR cores | randomized *Brachypodium* and soybean Hi-C passed with windows moved to the best score | 10 randomized controls; the real genomes |
| 0.12: the cis profile as observed/expected at each distance, chance from simulated contact maps | the first form passed on randomized sorghum, *Brachypodium*, moss and soybean Hi-C | 10 randomized controls; 38 genomes surveyed |
| 0.12: within 4 bins the cis signal neither supports nor opposes; at `pass` a confident call it contradicts loses a level | potato chr03, chr05, chr12 (a dip 1-2 bins off); tomato's two `medium` misses | Hi-C's existing tolerances (1 and 4 bins); potato, *Ae. tauschii* |
| 0.12: an AT trough places only what nothing else places, and pairs only with a Hi-C call it lies by | a code review before the release | the release: 40 paired calls, all right |
| 0.12: Hi-C anchored on the LTR cores counts a peak only where it is clean (not at a chromosome end, hub or NOR) | held-out *R. breviuscula*: 3 `medium` calls in a holocentric genome, its Athila cores at chromosome ends, which cluster | the 20 runs where the test runs, rerun (code/rel12b): only those 3 calls change, back to `low`; *Brachypodium* and *Arabidopsis* still pass; no untouched genome where the test runs was left to check |
{{RULES13}}

## 13. What changed from 0.13

{{CHANGES14}}

### From 0.12 to 0.13

{{CHANGES13}}

### From 0.11 to 0.12

{{CHANGES12}}
