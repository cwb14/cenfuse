# Plants and provenance (design, 2026-10-19)

## Intent

The data for one assembly come from many sources: the reference plant's own reads, other accessions, public
archives. Which data types exist differs by species. Centromeres can also sit elsewhere in other plants. In
*P. annua*, resequenced plants use Pa1A centromeres other than the reference's.

Without this, cenfuse would assume that every input describes the assembled plant, pool methylation over plants, and
make one call per chromosome. This adds three things:
1. **Provenance.** For every input, cenfuse measures whether it comes from the reference plant (its genotype),
   from another plant, or can't tell.
2. **Only reference-plant inputs and inputs of unknown origin decide the reference's centromere call.** An input
   shown to come from another plant never does.
3. **Every plant is scored on its own.** Does its centromere sit at the reference call or elsewhere (a shift)?

Non-goals:
- mapping calls between separate assemblies by synteny;
- genotyping Hi-C from read alignments (BAM);
- CENH3 ChIP as an input.

## Inputs

- `meth/<g>.tsv`: the cenmeth samples table. Its `sample` column names the plant.
- `hic/<g>.tsv` (new, optional): columns `sample` and `file` (pairs, cool or mcool; relative paths allowed). Rules:
  - One library per row, and `sample` names the plant.
  - A library and a methylation sample with the same name are the same plant.
  - Without a table, the single `hic/<g>.<ext>` is one library named `hic`.

## Provenance (per input)

1. **On this assembly (fail fast).** Methylation CHG rows must sit on reference CHG motifs: C, H, G on the + strand,
   or the reverse complement on the − strand. If fewer than half do, cenfuse stops with an error; below 99% it
   warns. A file mapped to another assembly falls far below half.
2. **Genotype (methylation samples).**
   - *Fixed differences.* modkit writes no row for a cytosine whose reads all show another base. At a symmetric CHG
     site (CAG/CTG: a CHG cytosine on each strand), a fixed genetic difference removes both strands' rows. ONT reads
     one strand per read, so read errors, failed calls and heterozygous sites remove one strand at a time.
   - *Tested sites.* A site is tested when reads span it on both strands: rows within 300 bp on both sides, at least
     3 reads per strand. A tested site with both rows missing is a fixed difference. The noise expected per site is
     q+ × q−, the single-strand loss rates.
   - *Regions.* Chromosome arms: per chromosome, up to two 2-Mb windows near 1/4 and 3/4 of its length. Windows avoid
     the sequence call (± 2 D), satellite clusters and gaps.
   - *Assembly errors are shared by every plant's reads.* So a difference that sample j has and sample k lacks, at
     sites tested in both, is a real genotype difference of j. j is **not the reference plant** if |F_j \ F_k|
     exceeds its noise (Poisson; P ≤ 0.01 after Bonferroni over the ordered pairs).
   - j is **the reference plant** (same genotype) if all three hold:
     - no such excess exists against any other sample;
     - at least one other sample is shown not to be the reference plant against j;
     - that test had ≥ 95% power to detect the smallest private-difference rate seen between plants shown to differ.
   - Otherwise **unknown**. This is always the case with one sample: assembly errors and another plant can't then be
     told apart.
3. **Karyotype (every input with coverage).**
   - *Coverage* is measured per window of max(1 Mb, the Hi-C bin). For methylation it is the summed valid calls
     (cenmeth class bins, before the depth filter); for Hi-C, contacts per bin.
   - *Missing segment.* At least 3 windows where the input has below 0.1 of its own median (a tenfold drop) while
     another input of the same data type has at least 0.5 of its median. Read length and mapping filters differ
     between ONT and Hi-C, so only same-type coverage shows that the DNA is mappable.
   - *Windows that cannot be judged* neither count toward a run nor break it:
     - windows overlapping a sequence call;
     - windows over 10% N, or over half satellite array (mapping, not DNA, sets their coverage);
     - chromosome ends;
     - windows no same-type input covers;
     - pile-ups above twice the input's median (collapsed repeats draw reads from elsewhere; cenmeth's DEPTH_HI).

     A window where the input has coverage ends a run.
   - *Meaning.* A missing segment shows that the input's plant lacks DNA the reference has, so it is not the
     reference plant. The absence of missing segments shows nothing.
4. **Verdict:**
   - `other` with genotype or karyotype evidence;
   - `reference` with genotype evidence;
   - `unknown` otherwise.

   An input of unknown provenance takes the verdict of a same-named input ("by name").

## Which inputs decide the reference call

- **Methylation.** The support for the integrated call pools only `reference` and `unknown` samples, and is judged
  with a signature calibrated on those samples alone, so other plants shape neither the evidence nor the rule. Without
  such samples, methylation neither supports nor helps move the reference call. The signature for scoring the plants
  is calibrated on all of them.
- **Hi-C.** The deciding library is a `reference` one; otherwise an `unknown` one. Among them, the library whose gate
  passes on the most chromosomes wins (its peaks land on the reference's sequence calls of the other chromosomes),
  then the one with the most reads. Libraries shown to come from other plants never decide.
  - An `unknown` library may decide, flagged `hic_provenance=unknown` in the calls table.
  - Reasons: Hi-C used to scaffold an assembly comes from the assembled plant, and the gate tests the library
    against the other chromosomes' sequence calls.
  - Known risk: an unknown library from another plant with one shifted chromosome can move that chromosome's call.
- **The other lines** (satellites, TEs, genes) describe the assembly itself.

## Per-plant scoring

Plants are the names in the `meth/` and `hic/` tables. Every plant is scored, including the reference plant when
present, as a positive control.

**Sites per chromosome**, fixed by the reference lines and the other plants rather than by one plant's data:
- the sequence call and the provisional call (the calls before methylation);
- the anchored and de novo Hi-C peaks;
- the TE marker core;
- centromeric-family satellite clusters away from the call;
- the peaks of the other libraries;
- every plant's `centromere_like` methylation candidates.

A site joins an earlier one (in that priority order) whose midpoint lies within W/2 of its own, and keeps its
origins. A site is *at the call* when its window overlaps the integrated call: calls can be wider than W (Pa1A:
6.4 Mb against W = 3.3 Mb). A plant's call site, s0, is the at-call site where it looks most centromeric.

**Methylation, per plant j and site s:**
- *Label and score:* cenmeth's support at s with sig_x, in a window of width W centered on s: `centromere_like`,
  `dip` (data, but not centromere-like), `no_signature_data` or `no_data`. The window's most extreme adjusted-CHG bin
  is chosen among the bins where the signature has data, so a low-depth plant is not labeled `no_signature_data`
  when its extreme falls on a bin without enough calls.
- *Between plants*, when at least 3 plants have data at s:
  - The statistic d_j(s) is j's signature score minus the median of the other plants' scores at the same site.
    The score is what defines `centromere_like`, and comparing plants on the same sequence cancels composition and
    mappability.
  - The null for chromosome x is d_j at the call sites of the other chromosomes, where plants mostly share the
    centromere (leave one chromosome out). It uses j's own values, or every plant's when j has fewer than 5. A new
    value is judged with Student's t prediction interval, t = (d − mean) / (sd √(1 + 1/n)) on n − 1 degrees of
    freedom, which carries the uncertainty of a spread estimated from about a dozen chromosomes.
  - **Loss** is tested at j's call site (the at-call site where j looks most centromeric): is j less centromeric
    there than the others?
  - **Gain** is tested at sites away from the call that do not come only from j's own candidates, since those were
    selected on j's data.
  - Benjamini-Hochberg over every test in the genome.

**Verdict per plant × chromosome:**

| verdict | condition |
|---|---|
| `call_absent` | the call lies in a missing segment of this plant |
| `no_data` | no signature data at any site at the call |
| `shared` (consistent with the call) | `centromere_like` at s0, and no significant loss there |
| `shifted` | `dip` at s0, `centromere_like` at an alternative, and q ≤ 0.05 for loss at s0 or gain at that alternative |
| `possible_shift` | `dip` at s0 and `centromere_like` at an alternative, without between-plant support (or with fewer than 3 plants) |
| `unclear` | any other combination |

**Hi-C, per library L that does not decide:** the anchored test with partner windows at the reference calls, and
L's own gate. Per chromosome:
- `shared`: the peak agrees with the call; or it is a strong peak beyond 4 bins that agrees with the deciding
  library's peak (the plant shares the reference plant's functional centromere, even where the reference call is
  contested);
- `shifted`: a strong peak beyond 4 bins, with the gate passing, where the deciding library peaks elsewhere;
- `unclear`: a near or weak peak;
- `uninformative`: the gate fails.

**Plant verdict** combines its lines:
- `call_absent` if the call lies inside a missing segment of the plant;
- `shifted` if a line says shifted and none says shared;
- `conflict` if one line says shifted and another says shared;
- `shared` if a line says shared;
- `possible_shift` if a line says so;
- `no_data` if no line has data;
- otherwise `unclear`.

`no_call` marks chromosomes without a reference call; they carry no tests.

**What `shared` means.** Methylation that looks centromeric at the call is consistent with the reference call, not
proof. On *P. annua*, some places that are not centromeres look centromeric in every plant:
- Pa4B 58.3 Mb scores 2.4-3.4 in all five plants, higher than the real Pa4B centromere;
- Pa5A 71.5 Mb scores 1.1-2.6 in all five plants.

The between-plant test keeps such places from producing `shifted` verdicts. By the same token, a centromere-like score
at a call does not show that a plant's centromere is there.

## Outputs

- `<g>.provenance.tsv`: one row per input. Columns:
  - kind, plant, verdict, reason, share of rows on reference motifs;
  - sites tested, fixed differences per kb, noise;
  - the strongest private-difference test (against which sample, count, expected, P);
  - missing segments, and whether the input decides.
- `<g>.plants.tsv`: one row per plant × chromosome. Columns:
  - provenance and verdict;
  - methylation at the call; the best alternative site with its origins, label and score; loss q and gain q;
  - the Hi-C verdict and peak;
  - distance from the call.
- `<g>.plant_sites.tsv`: one row per plant × site, with label, score, value, d, t, P, q and origins.
- `<g>.centromeres.tsv`, columns for this:
  - `meth_samples`: the samples whose methylation decided;
  - `hic_library`, `hic_provenance`;
  - `plants_shared`, `plants_shifted`, `plants_absent`, `plants_unresolved` (`.` when no plant is scored).
- `cenfuse_summary.tsv`: inputs by verdict, and the plants with a shift.
- **Figure:**
  - a plants × chromosomes panel with each cell's verdict;
  - per genome, the chromosome with the most shifts, drawn plant by plant along its length.

Coordinates in every TSV are 1-based inclusive. The provenance table's `read_ends` counts both ends of every Hi-C
pair (the matrix's row sums).

## Constants (`params.py`)

- `GENO_FLANK = 300` bp, `GENO_MIN_READS = 3` per strand. With ONT reads kb long, rows on both sides mean the reads
  span the site, and three reads all failing by chance has probability ~0.1³.
- `GENO_WINDOW = 2 Mb`, two per chromosome: enough sites at 2-3× coverage (*P. annua* probe: ~150,000).
- `GENO_ALPHA = 0.01` (conventional; Bonferroni over pairs). `GENO_POWER = 0.95` (conventional).
- `KARYO_ABSENT = 0.1`, a tenfold drop, the rule cenfuse uses for chromosome sets. `KARYO_PRESENT = 0.5`, half the
  median. `KARYO_MIN_WINDOWS = 3`. All three are heuristics, named as such.
- `PLANT_Q = 0.05` (Benjamini-Hochberg, conventional).
- Module constants, each with its reason in the code:
  - `provenance.py`: `MOTIF_MIN` 0.5, `MOTIF_WARN` 0.99, `MOTIF_SPAN` 200 kb, `END_MARGIN` 0.05, `GAP_SHARE` 0.1,
    `SAT_SHARE` 0.5;
  - `plants.py`: `MIN_NULL` 5.

## Validation

Pre-registered in `benchmark.md` before any real-data run: a synthetic genome with planted plants, and the five
*P. annua* plants.

## Changes during development (synthetic genomes only, before any real-data run)

- **The between-plant statistic.** At first it was the window's extreme adjusted-CHG value, with a null from every
  window of the chromosome. On the synthetic genome this made shared centromeres look lost (P1 chr10 q 1e-5,
  P2 chr8 q 0.049). The statistic varies most where the signal is strongest, so arm windows understate its spread,
  and one extreme bin is noisy. Changed to the signature-score difference, with a null from the other chromosomes'
  call sites.
- **The null's tail.** A normal tail on a median-absolute-deviation scale from nine values still gave false losses
  in pure-noise genomes. Changed to the t prediction interval. On 300 pure-noise genomes (4 plants × 10
  chromosomes):
  - 5.0% of 24,000 null tests have P < 0.05;
  - 14/300 genomes have any q ≤ 0.05, within the 5% Benjamini-Hochberg allows when every test is null;
  - 7 of 12,000 verdicts are not `shared`, and none is `shifted`.
- **At the call.** A site counts as at the call when its window overlaps the call, not only when its midpoint lies
  inside. A Hi-C peak one bin outside a call is the same centromere.
- **Hi-C shifts.** A far peak counts as a shift only if the deciding library's peak is elsewhere. On synthetic chr9
  (a contested call), the other plant's library peaked where the reference library did.

## Changes after the first *P. annua* run (karyotype only; the pre-registered results are kept)

The first run failed R2 (none of the five known missing segments found) and R5a
(P3 and R1 at Pa2A `no_data`, not `call_absent`). The coverage showed every event clearly. The detector failed for
reasons that do not depend on the answers:
- **Skipping.** The old rule skipped every window holding a satellite cluster of at least 20 kb. In this genome that is
  53% of all 1-Mb windows and 274 of Pa1A's 320, although satellites cover a median 2.6% of a window. Now only
  windows over half satellite are skipped (1.9% of windows).
- **Windows no one covers.** These broke runs, as if missing data were evidence of presence. Now they cannot be
  judged.

A rerun found the events but split Pa1A at 204 Mb. At 204 Mb, every ONT sample has 3-15x its
median (a collapsed repeat). It also flagged I3 and P3 across 203-211 Mb, where every ONT sample is near zero and only
the Hi-C library has reads. Two more rules, again independent of the answers:
- **Same-type coverage only.** Short Hi-C reads (no mapping-quality filter in that pairs file) and long ONT reads
  (MAPQ >= 20) map differently.
- **Pile-ups cannot be judged** (cenmeth's own depth rule).

The same-type rule was prompted by a check meant to be held out (no segment confined to a low-coverage patch shared
by all samples), so that check no longer tests it. Checks of the next rerun on cases not used for any change:
- R1, monosomic for Pa5B (~0.5x), has no segment there;
- the reference Hi-C library has none;
- P3's 2.5-Mb Pa1A loss is too short to detect.

## After the code review (2026-10-19; then a final rerun)

An independent review found these defects; all are fixed, with regression tests:
- **Label bin.** The site label took the most extreme bin and only then asked for signature data there. Low-depth
  plants were `no_signature_data` with usable bins in the window: P3, R1, I3 and J10 at Pa1A had 11-36 such bins.
  The label now uses the most extreme bin among bins with signature data.
- **Call absent.** `call_absent` accepted a call up to one window beyond a segment. The call must now lie inside it.
- **Deciding library.** Among unknown libraries, the most reads won whatever the gate said. In a synthetic case, a
  library from a plant with a shifted chromosome decided, its gate passing on 5/10 chromosomes against 10/10 for the
  reference library. Gate passes now come first.
- **Calibration.** Other plants still shaped the signature that judges the deciding samples. Now the deciders are
  calibrated alone.
- **Smaller fixes:**
  - chromosomes without a call carry no tests;
  - outputs are 1-based;
  - `contacts` is renamed `read_ends`;
  - the karyotype stage runs only when some data type has two inputs;
  - the assembly check reads 200 kb per window;
  - the between-plant statistic is vectorized;
  - the legend wording is corrected.

Two concerns remain untested:
- **The genotype noise model** (both-strand losses independent, at one rate across sites) has not been tested on two
  samples of one plant. A split of one sample's reads, or one plant sequenced with two chemistries, would test it.
  Chemistry-specific errors could make the reference plant's own reads look like another plant's.
- **Coverage at low depth.** Methylation coverage counts only cytosines with at least 3 calls (cenmeth's floor), so
  it is not proportional to depth. At ~2.5x, a threefold depth drop can fall below the tenfold "absent" level. The
  *P. annua* events sit at <= 0.03 of the median, far below it, and half depth (R1's monosomic Pa5B) stays above.
