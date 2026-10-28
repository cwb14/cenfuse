# Several plants: provenance and per-plant scoring

Methylation samples and Hi-C libraries often come from plants other than the assembled one, and centromeres can sit
elsewhere in other plants (in *P. annua*, other accessions use other Pa1A centromeres). cenfuse therefore:

1. **tells where each input comes from**: the assembled (reference) plant, another plant, or unknown;
2. **lets only reference and unknown inputs decide the reference calls**;
3. **scores every plant on its own**: is its centromere at the reference call (`shared`) or elsewhere (`shifted`)?

Not goals: mapping calls between assemblies, genotyping Hi-C from BAMs, CENH3 ChIP as an input.

## Inputs

- `meth/asm.tsv`: one row per sample; `sample` names the plant.
- `hic/asm.tsv` (optional, for several libraries): columns `sample`, `file` (pairs, cool or mcool).
- The same `sample` name in both tables means the same plant. A lone `hic/asm.<ext>` is one library named `hic`.

## Provenance (per input)

| step | rule |
|---|---|
| right assembly? (fail fast) | methylation CHG rows must sit on reference CHG motifs: < 50% stops the run, < 99% warns |
| genotype (methylation) | at symmetric CHG sites read on both strands (≥ 3 reads each), a fixed genetic difference removes both strands' rows; read errors remove one. A difference sample j has and sample k lacks, beyond its noise (Poisson, P ≤ 0.01, Bonferroni), makes j **not** the reference plant. Sites: two 2-Mb windows per chromosome on the arms. CHG files must be modkit bedMethyl (the genotype reads its counts of differing and failed calls); one in another layout is left out with a warning, and the sample's CG file is kept |
| reference plant | no such excess against any sample, another sample shown different, and ≥ 95% power. With one sample: always `unknown` (assembly errors and another plant look alike) |
| karyotype (any input with coverage) | a **missing segment** = ≥ 3 windows below 0.1× the input's median while a same-type input has ≥ 0.5× (windows over calls, gaps, mostly-satellite windows, ends and pile-ups are not judged). A missing segment means another plant |
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

A plant's verdict combines its lines: `shifted` when one line says so and none says `shared`; `conflict` when the
two lines disagree.

**`shared` is consistent with the call, not proof.** In *P. annua* some non-centromeric places look centromeric in
every plant. The between-plant test keeps them from producing false `shifted` verdicts.

## Outputs

| file | content |
|---|---|
| `asm.provenance.tsv` | per input: plant, verdict, reason, genotype test, missing segments, whether it decides |
| `asm.plants.tsv` | per plant × chromosome: verdict, methylation at the call and the best alternative, Hi-C verdict |
| `asm.plant_sites.tsv` | per plant × site: label, score, d, t, P, q |
| `asm.centromeres.tsv` | `meth_samples`, `hic_library`, `hic_provenance`, `plants_shared/shifted/absent/unresolved` |
| `asm.plants.png` | plants × chromosomes verdicts, and the most-shifted chromosome plant by plant |

## Validation ([benchmark.md](benchmark.md))

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
| KARYO_ABSENT, KARYO_PRESENT, KARYO_MIN_WINDOWS | 0.1, 0.5, 3 | heuristics, named as such |
| PLANT_Q | 0.05 | Benjamini-Hochberg, conventional |

## Untested

- **Genotype noise model**: never tested on two samples of one plant. Chemistry-specific errors could make the
  reference plant look like another.
- **Coverage at low depth**: below ~2.5×, a threefold drop can fall under the tenfold "absent" level.
