# Slides: "cenfuse: an integrative and accurate approach to calling centromeres in eukaryotes"

Images for a talk, in order, one per slide (16:9, 300 dpi). The slide title is a suggestion; the takeaway is what the
image shows. Every number is from the cenfuse 0.10.1 benchmark ([benchmark.md](../../benchmark.md)).

Letter suffixes mark slides added between the numbered ones; `S` slides are supplements, for questions.

## What cenfuse does and how

| # | image | suggested slide title | takeaway |
|---|---|---|---|
| 1 | [01_centromere_diversity.png](01_centromere_diversity.png) | Centromeres are built from very different DNA | megabase satellite arrays, a few retrotransposons, a 118-bp point, AT-rich cores, whole holocentric chromosomes: no single sequence signal finds them all |
| 2 | [02_one_line_is_not_enough.png](02_one_line_is_not_enough.png) | No single line of evidence finds centromeres everywhere | which evidence works changes by lineage; tested and combined, cenfuse gets close to the best line picked in hindsight (79% vs 83% of chromosomes right) |
| 3 | [03_workflow.png](03_workflow.png) | How cenfuse calls a centromere | the sequence call stands unless evidence that passed a test in this genome moves it |
| 3a | [03a_tool_satellite_families.png](03a_tool_satellite_families.png) | Finding satellite families from scratch | a head-to-tail repeat matches itself one unit along, so no prior knowledge is needed |
| 3b | [03b_tool_satellite_roles.png](03b_tool_satellite_roles.png) | A satellite's place on the chromosomes says what it is | the call is the centromeric family's main block; it starts high, and warnings lower it |
| 3c | [03c_tool_retrotransposon_cores.png](03c_tool_retrotransposon_cores.png) | Which retrotransposon lineage may decide | CRM piles up once per chromosome, on CENH3; other lineages only flag |
| 3d | [03d_tool_hic_peaks.png](03d_tool_hic_peaks.png) | How cenfuse reads a centromere from Hi-C | centromeres gather, so each chromosome's contacts peak at its own: two ways of asking, three checks, one trusted peak |
| 4 | [04_lines_earn_a_vote.png](04_lines_earn_a_vote.png) | Hi-C gets a vote only if its peaks land on the other chromosomes' calls | on most other chromosomes the peak must land on that chromosome's own call; otherwise Hi-C is only reported (15 of 34 genomes pass) |
| 4a | [04a_tool_rabl_axis.png](04a_tool_rabl_axis.png) | The Rabl axis: two more readouts from the same Hi-C reads | where nuclei hold centromeres at one pole, side-by-side arms mark the centromere with no sequence call |
| 4b | [04b_tool_cpg_kinetochore.png](04b_tool_cpg_kinetochore.png) | How the CpG line finds the kinetochore | pool noisy CpGs over the arrays, rank the valleys, mark the top one's half-depth run; no call changes |
| 4c | [04c_tool_chg_methylation.png](04c_tool_chg_methylation.png) | CHG methylation: learned from the other chromosomes | it can support a call but never moves one alone; the signal is real in *Poa*, and it changed no call in this benchmark |
| 4d | [04d_integration_rules.png](04d_integration_rules.png) | No weights: a start level, one-level steps, moves only by rule | tested lines step a chromosome's level one rung at a time, and a confident call needs two lines to move; 678 of 800 chromosomes keep their level and place |
| 4e | [04e_integration_examples.png](04e_integration_examples.png) | Real chromosomes through the rules | seven real chromosomes, one through each rule, with the level stepping beside them |
| 5 | [05_example_output.png](05_example_output.png) | One page of output | calls, the evidence for and against each, Hi-C and methylation tracks, telomeres and rDNA (*Poa annua*, all four inputs) |

## How well it works

| # | image | suggested slide title | takeaway |
|---|---|---|---|
| 6 | [06_benchmark_scope.png](06_benchmark_scope.png) | A benchmark across eukaryotes | 55 assemblies across eukaryotes, 908 chromosomes; CENP-A/CENH3 truth on 633 chromosomes of 47 genomes |
| 7 | [07_accuracy.png](07_accuracy.png) | Right on most chromosomes, unsure rather than wrong | with all inputs 532 of 633 right; of the 101 misses, 84 are flagged `low` or left uncalled, 17 are confidently wrong |
| 8 | [08_calibration.png](08_calibration.png) | Confidence is honest | `high` calls right 329/343, `medium` 102/105, `low` 101/155; no confident call on 46 holocentric chromosomes |
| 9 | [09_mammals_t2t.png](09_mammals_t2t.png) | T2T mammals: 148 of 157 centromeres right | 143/157 CENP-A centromeres from satellites alone, 148 with genes and Hi-C; the 4 confident errors shown with their raw data |
| 10 | [10_retrotransposon_centromeres.png](10_retrotransposon_centromeres.png) | Retrotransposon centromeres: satellites alone miss them | potato 2 → 6/11, *Ae. tauschii* 1 → 5/7, *Chlorella* 0 → 8/10, moss with Hi-C 25/26; satellite centromeres stay right |
| 11 | [11_grasses_and_crops.png](11_grasses_and_crops.png) | Crops: which evidence places each centromere | satellites, CRM cores or Hi-C's arm fold, chromosome by chromosome; tomato fails with all |
| 12 | [12_hic_moves_calls.png](12_hic_moves_calls.png) | Tested Hi-C moves wrong satellite calls onto CENH3 | barley: the arm fold and Hi-C agree on CENH3, 0 → 6/7 right |
| 13 | [13_fungi_and_protists.png](13_fungi_and_protists.png) | Hi-C finds centromeres the sequence cannot place | point, regional and AT-rich centromeres: satellites alone 1/112, with Hi-C 82/112, none confidently wrong |
| 13a | [13a_callers_compared.png](13a_callers_compared.png) | Where satellite callers work, and where Hi-C is needed | TRASH2 and satellome find the centromere where a dominant tandem array marks it, and tie or beat cenfuse on some rows; on 87 chromosomes where cenfuse needs Hi-C they place 10 and 21 |
| 13b | [13b_callers_examples.png](13b_callers_examples.png) | Why: where the largest tandem arrays sit | four typical chromosomes: rice, barley, yeast, *C. neoformans* |
| 14 | [14_rabl_axis.png](14_rabl_axis.png) | Contacts between chromosomes mark centromeres in Rabl-like nuclei | where centromeres cluster at one pole, arms of different chromosomes align by distance from the centromere; absent in G1, territorial and holocentric nuclei |
| 15 | [15_cpg_kinetochore.png](15_cpg_kinetochore.png) | From a noisy CpG track to the kinetochore | the most prominent CpG valley places it within 50 kb of CENP-A on 35 of 39 ape centromeres |
| 16 | [16_controls.png](16_controls.png) | Where no single centromere exists, cenfuse makes no confident call | holocentric genomes: 0 confident calls of 46 (Hi-C used untested would make 21); every Hi-C test fails, so Hi-C changes nothing |
| 17 | [17_robustness.png](17_robustness.png) | Thinner data, same answers | Hi-C, CpG and LTR annotation thinned: a failed test falls back to the calls made without that data, never worse |
| 18 | [18_beyond_centromeres.png](18_beyond_centromeres.png) | One map of the chromosome | telomeres, interstitial telomeric arrays (the human 2q13 fusion at base-pair resolution), rDNA, knobs, completeness |
| 19 | [19_karyotype_morphology.png](19_karyotype_morphology.png) | Karyotypes from the sequence | acrocentric vs bi-armed matches literature karyotypes on 117 of 118 right calls |
| 20 | [20_where_it_fails.png](20_where_it_fails.png) | Where cenfuse fails: 17 of 448 confident calls | 17 confident errors in 4 kinds, each with the data it had and the data that would settle it |
| 21 | [21_runtime_memory.png](21_runtime_memory.png) | Run time is the satellite search | run time follows the satellite DNA in a genome, not its size; everything except the satellite search takes minutes |

## Supplements

| # | image | suggested slide title | takeaway |
|---|---|---|---|
| S1 | [S1_telomere_arrays.png](S1_telomere_arrays.png) | Telomeres: the unit is learned from the chromosome ends | no motif is given: the repeat far denser at the ends is the unit, and its copies map the arrays (a unit in 51 of 55 assemblies) |
| S2 | [S2_rdna_arrays.png](S2_rdna_arrays.png) | rDNA arrays: a run of copies, not a single hit | a conserved mark recurs once per unit; runs are kept, lone copies are ignored |
| S3 | [S3_rdna_satellites_and_knobs.png](S3_rdna_satellites_and_knobs.png) | Big satellite blocks that are not centromeres | a repeat among rRNA genes is tagged and never called; knobs stay listed beside the call |
| S4 | [S4_hor_and_long_units.png](S4_hor_and_long_units.png) | Two more parts of the satellite search | higher-order repeats are reported and change no call; long units are a rare last resort that starts unsure |
| S5 | [S5_several_plants.png](S5_several_plants.png) | Inputs from other plants are scored, and never decide | cenfuse tells where each input came from; only reference or unknown inputs decide |
| S6 | [S6_call_in_words.png](S6_call_in_words.png) | Every call explains itself in words | each call says what it stands on, what speaks against it and, where cenfuse can tell, which data would settle it |

Each image's build script, data and legend pages (a short legend that reads as speaker notes, and a long one that also
keeps the chance expectations and negative controls left off the slides) are kept with the benchmark runs, outside the
repository.
