# Several plants: provenance and per-plant scoring

Methylation samples and Hi-C libraries often come from plants other than the assembled one, and centromeres can sit
elsewhere in other plants (in *P. annua*, young CRM insertions suggest that other accessions use other Pa1A
centromeres; not yet confirmed). cenfuse therefore:

1. **tells where each input comes from**: the assembled (reference) plant, another plant, or unknown;
2. **lets only reference and unknown inputs decide the reference calls**;
3. **scores every plant on its own**: is its centromere at the reference call (`shared`) or elsewhere (`shifted`)?

Not goals: mapping calls between assemblies, genotyping Hi-C from BAMs, CENH3 ChIP as an input.

## Inputs

- `meth/asm.tsv`: one row per sample; `sample` names the plant.
- `hic/asm.tsv` (optional, for several libraries): columns `sample`, `file` (`.pairs[.gz]`, `.cool`, `.mcool`, HiC-Pro
  `.allValidPairs[.gz]` or Juicer `.merged_nodups.txt[.gz]`, relative to `hic/`). Rows with the same `sample` are replicate
  libraries of one plant, summed into one library (cool or mcool files not mixed with pair files in one sample). A lone
  `hic/asm.<ext>` beside the table is ignored, with a warning.
- The same `sample` name in both tables means the same plant; names are non-empty, without `/` or whitespace (they name
  output files). A lone `hic/asm.<ext>` is one library named `hic`.
- The CHG line needs a CHG file (modkit bedMethyl) for every sample: one CG-only sample, or one whose CHG file is in another
  layout, turns the CHG line off for every sample, and with it per-plant methylation scoring and the methylation karyotype
  test.

## Provenance (per input)

| step | rule |
|---|---|
| right assembly? | every methylation file's rows must sit on reference cytosines of its context (CHG, CG, CHH; a CG bedGraph or pb-CpG-tools file on CpGs): < 50% stops the run, < 99% warns; checked after the satellite search |
| genotype (methylation) | at symmetric CHG sites read on both strands (rows within 300 bp on both sides of each strand's C, averaging ≥ 3 reads), a fixed genetic difference removes both strands' rows; read errors remove one. A difference sample j has and sample k lacks, beyond its noise (Poisson, P ≤ 0.01, Bonferroni), makes j **not** the reference plant. Sites: two 2-Mb windows per chromosome on the arms. CHG files must be modkit bedMethyl (the genotype reads its counts of differing and failed calls); one in another layout is left out with a warning, and the sample's CG file is kept |
| reference plant | no such excess against any sample, another sample shown different, and ≥ 95% power. With one sample: always `unknown` (assembly errors and another plant look alike) |
| karyotype (Hi-C libraries, and methylation samples when the CHG line runs; two or more inputs of a type) | a **missing segment** = ≥ 3 windows (1 Mb, or the Hi-C bin if larger) below 0.1× the input's median while a same-type input has ≥ 0.5× (windows over calls, gaps, mostly-satellite windows, ends and pile-ups are not judged). A missing segment means another plant |
| verdict | `other` (genotype or karyotype evidence), `reference` (genotype), else `unknown`; an unknown input takes a same-named input's verdict |
| CG-only samples (CpG line) | not genotyped (genotypes are read at CHG sites) and not in the karyotype test: `unknown`, so they decide; their rows must still sit on reference CpGs |

## Which inputs decide the reference call

- **Methylation**: support pools only `reference` and `unknown` samples, with a signature calibrated on them alone.
  Without such samples, methylation does not decide.
- **Hi-C**: a `reference` library, else an `unknown` one. Among those, the one whose gate passes on the most
  chromosomes, then the most reads. Other-plant libraries never decide.
  - Risk: an unknown library from another plant with one shifted chromosome can move that call.

## Per-plant verdicts

**Sites per chromosome:**
- the sequence and provisional calls;
- the Hi-C peaks of every library;
- the TE core;
- centromeric satellite clusters away from the call;
- every plant's `centromere_like` methylation candidates.

Nearby sites merge.

**Methylation, plant j at site s:**
- the label `centromere_like` / `dip` / `no_signature_data` / `no_data`, and j's signature score;
- between plants (≥ 3 with data): d = j's score minus the other plants' median at s;
- the null for chromosome x comes from the other chromosomes' call sites (Student's t prediction interval);
- loss is tested at j's call site, gain at alternative sites;
- Benjamini-Hochberg over the genome (q ≤ 0.05).

| verdict | condition |
|---|---|
| `call_absent` | the call lies in a missing segment of this plant |
| `shared` | `centromere_like` at the call, no significant loss |
| `shifted` | `dip` at the call, `centromere_like` elsewhere, q ≤ 0.05 for the loss or the gain |
| `possible_shift` | the same without between-plant support (or < 3 plants) |
| `no_data` / `unclear` | no data / anything else |

**Hi-C, per non-deciding library**: `shared` when its peak agrees with the call (or with the deciding library's
peak); `shifted` when a strong far peak passes its own gate where the deciding library peaks elsewhere; `unclear`;
`uninformative` (gate fails).

A plant's verdict combines its lines: `shifted` when one line says so and none says `shared`; `conflict` when one
line says `shifted` and the other `shared`.

**`shared` is consistent with the call, not proof.** In *P. annua* some non-centromeric places look centromeric in
every plant. The between-plant test keeps them from producing false `shifted` verdicts.

## Outputs

| file | content |
|---|---|
| `asm.provenance.tsv` | per input: plant, verdict, reason, genotype test, missing segments, whether it decides (columns below) |
| `asm.plants.tsv` | per plant × chromosome: verdict, methylation at the call and the best alternative, Hi-C verdict (columns below) |
| `asm.plant_sites.tsv` | per plant × site: label, score, d, t, P, q (columns below) |
| `asm.centromeres.tsv` | `meth_samples`, `hic_library`, `hic_provenance`, `plants_shared/shifted/absent/unresolved` |
| `asm.plants.png` | plants × chromosomes verdicts, and the most-shifted chromosome plant by plant |

### `asm.provenance.tsv` (one row per input)

| column | meaning |
|---|---|
| `kind`, `input` | the input: a methylation sample (`meth`) or a Hi-C library (`hic`), by its `sample` name |
| `plant` | the plant it comes from: the same name |
| `verdict`, `reason` | `reference`, `other` or `unknown` (Provenance, above), and why, in words |
| `decides` | whether it decides the reference calls: a methylation sample that is not `other`, when the CHG line runs; the deciding Hi-C library (the CG files of samples that are not `other` decide the CpG line whatever this column says) |
| `on_motif` | methylation (CHG file): the share of its rows in the genotype windows that sit on reference CHG motifs |
| `sites_tested` | how many symmetric CHG sites (CAG, CTG) in those windows are read on both strands: each strand's C has rows within 300 bp on both sides, with ≥ 3 reads on average |
| `fixed_per_kb` | the fixed differences (both strands' rows missing) per kb of tested bases (3 per site) |
| `noise_per_site` | the chance that read errors remove both rows of a site: the share of tested sites missing only the + strand's row times the share missing only the − strand's (add-one smoothed) |
| `private_vs`, `private`, `private_expected`, `private_p`, `shared_sites` | its private-difference test with the lowest P: the other sample; this sample's fixed differences that the other lacks, at the sites both tested; how many its noise predicts there; the Poisson upper-tail P, uncorrected (a verdict needs P ≤ 0.01 over the number of ordered pairs of samples); those shared sites |
| `missing_segments` | its missing segments (`chrom:start-end`, `;`-separated) |
| `read_ends` | Hi-C libraries: its contacts counted at both ends (the per-bin totals of its contact matrix, summed). Among the libraries that may decide (`reference`, else `unknown`), the one whose gate passes on the most chromosomes decides; a tie goes to the most read ends |

### `asm.plants.tsv` (one row per plant and chromosome)

`alt_start`, `alt_end`, `hic_peak_start` and `hic_peak_end` are 1-based; `-1` = none.

| column | meaning |
|---|---|
| `plant`, `provenance` | the plant (a methylation sample, or a Hi-C library that does not decide) and its provenance verdict |
| `chrom` | the chromosome |
| `verdict` | the plant's verdict, its lines combined (above): `shared`, `shifted`, `possible_shift`, `conflict`, `call_absent`, `unclear`, `no_data`, or `no_call` (no reference call) |
| `meth_verdict` | the methylation verdict (above); `.` without a methylation sample |
| `call_site`, `call_label`, `call_score`, `loss_q` | the plant's site at the call (in `asm.plant_sites.tsv`): of the sites whose window overlaps the call, the one where it scores highest (without scores, the one nearest the call's middle); its label and score there, and the q of its loss test |
| `alt_site`, `alt_start`, `alt_end`, `alt_origins`, `alt_label`, `alt_score`, `gain_q` | the alternative: the highest-scoring site away from the call where the plant is `centromere_like` (among those whose gain has q ≤ 0.05, if any): its window, origins, label, score, and the q of its gain test |
| `shift_bp` | with `shifted`, `possible_shift` or `conflict`: bp from the call to the alternative, or to the plant's Hi-C peak when its methylation verdict is neither `shifted` nor `possible_shift` |
| `hic_library`, `hic_verdict`, `hic_peak_start`, `hic_peak_end` | the plant's own Hi-C library (one that does not decide), its verdict (`shared`, `shifted`, `unclear`, `uninformative`, `no_data`; above) and its anchored peak on this chromosome; `.` without one |

### `asm.plant_sites.tsv` (one row per plant and site)

Every methylation sample is scored at every site of its chromosome (Per-plant verdicts, above). `start` and `mid` are
1-based.

| column | meaning |
|---|---|
| `chrom`, `site` | the chromosome and the site (`chrom:S<n>`, numbered as the places merge into sites, the reference lines' places first) |
| `start`, `end`, `mid` | the window judged, centered on the site's midpoint (that of the first place merged into it): 2 D or one Hi-C bin wide, whichever is larger (D: the methylation engine's expected centromere size) |
| `origins` | the places merged into the site: `sequence_call`, `hic_peak` (the deciding library's anchored peak), `te_core`, `provisional` (the call before methylation), `satellite:<family>` (a cluster of a centromeric family), `denovo_hic` (the deciding library's prior-free peak), `hic:<library>` (another library's anchored peak), `meth:<sample>` (a sample's `centromere_like` candidate) |
| `sample` | the plant |
| `label`, `score`, `value`, `with_data` | at the plant's most extreme adjusted-CHG bin of the window: `centromere_like` (score ≥ 0.5), `dip`, `no_signature_data` or `no_data`; the signature score (1 = deviates like the other chromosomes' calls, 0 = like the background); the adjusted CHG value; the signature entries with data there |
| `inside` | the window overlaps the reference call |
| `plants_with_data` | how many plants have a score at this site |
| `d` | the plant's score minus the median of the other plants' (with ≥ 3 plants scored) |
| `call_site` | this is the plant's site at the call (`asm.plants.tsv`) |
| `t`, `df` | d against its null, the plant's d at the call sites of the other chromosomes (every plant's, when it has fewer than 5): t = (d − mean) / (sd × √(1 + 1/n)), on n − 1 degrees of freedom |
| `p_gain`, `p_loss` | the upper and the lower tail P of t |
| `test` | the test the row carries: `loss` (the plant's call site), `gain` (a site away from the call that does not come from the plant's own candidates alone), or empty (none: a test needs ≥ 3 plants scored and a reference call) |
| `q` | the Benjamini-Hochberg q, over every test of the genome; `.` without a test |

## Validation ([benchmark.md](benchmark.md#several-plants-poa-annua))

| test | result |
|---|---|
| synthetic genome with planted plants (provenance, planted shift, bad inputs) | pass |
| 300 pure-noise genomes | 5.0% of null tests at P < 0.05; no `shifted` verdict |
| *P. annua*: five other plants recognized as other | 5/5 |
| *P. annua*: known missing segments | all 5 found; the two Pa2A losses in 3 pieces each |
| *P. annua*: Pa1A per-plant positions | methylation matches the insertion-based position only for J10. At 2.5-8.6× ONT, per-plant methylation does not resolve this polymorphism |

## Constants (`cenfuse/params.py`)

| constant | value | reason |
|---|---|---|
| GENO_FLANK, GENO_MIN_READS | 300 bp, 3 reads | ONT reads span the site; three chance failures ~0.1³ |
| GENO_WINDOW | 2 Mb × 2 per chromosome | ~150,000 sites at 2-3× |
| GENO_ALPHA, GENO_POWER | 0.01, 0.95 | conventional |
| KARYO_WINDOW | 1 Mb (at least the Hi-C bin) | the coverage window for missing segments |
| KARYO_ABSENT, KARYO_PRESENT, KARYO_MIN_WINDOWS | 0.1, 0.5, 3 | heuristics, named as such |
| PLANT_Q | 0.05 | Benjamini-Hochberg, conventional |

## Untested

- **Genotype noise model**: never tested on two samples of one plant. Chemistry-specific errors could make the
  reference plant look like another.
- **Coverage at low depth**: below ~2.5×, a threefold drop can fall under the tenfold "absent" level.
- **Verdicts that depend on other chromosomes**: a plant's null comes from the other chromosomes' call sites and its q from every
  test of the genome, so a verdict can change when another chromosome's call moves (*P. annua* J10 Pa4B: `shifted` to
  `possible_shift` when Pa1A's call moved; [future.md #57](future.md#57-per-plant-verdicts-that-depend-on-other-chromosomes)).
