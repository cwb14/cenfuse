# Slides: "cenfuse: an integrative and accurate approach to calling centromeres in eukaryotes"

Images for a talk, in order, one per slide (16:9, 300 dpi). The slide title is a suggestion; the takeaway is what the
image shows. Every number is from the cenfuse 0.10.1 benchmark ([benchmark.md](../../benchmark.md)).

| # | image | suggested slide title | takeaway |
|---|---|---|---|
| 1 | [01_centromere_diversity.png](01_centromere_diversity.png) | Centromeres are built from very different DNA | megabase satellite arrays, a few retrotransposons, a 118-bp point, AT-rich cores, whole holocentric chromosomes: no single sequence signal finds them all |
| 2 | [02_one_line_is_not_enough.png](02_one_line_is_not_enough.png) | Every kind of evidence fails somewhere | each line fails on some lineage, and which one works cannot be known in advance; testing each in the genome gets close to the best line picked in hindsight (79% vs 83% of chromosomes right) |
| 3 | [03_workflow.png](03_workflow.png) | How cenfuse works | each line calls on its own, is tested in the genome, and only then may confirm or move a call |
| 4 | [04_lines_earn_a_vote.png](04_lines_earn_a_vote.png) | A line is used only after it finds the other chromosomes' centromeres | real Hi-C passes the test in 15 of 34 genomes with calls to test against; 0 of 47 randomized libraries pass |
| 5 | [05_example_output.png](05_example_output.png) | One page of output | calls, the evidence for and against each, Hi-C and methylation tracks, telomeres and rDNA (*Poa annua*, all four inputs) |
| 6 | [06_benchmark_scope.png](06_benchmark_scope.png) | A benchmark across eukaryotes | 55 assemblies across eukaryotes, 908 chromosomes; CENP-A/CENH3 truth on 633 chromosomes of 47 genomes |
| 7 | [07_accuracy.png](07_accuracy.png) | Right on most chromosomes, unsure rather than confidently wrong | with all inputs 532 of 633 right; of the 101 misses, 84 are flagged `low` or left uncalled, 17 are confidently wrong |
| 8 | [08_calibration.png](08_calibration.png) | Confidence is honest | `high` calls right 329/343, `medium` 102/105, `low` 101/155; no confident call on 46 holocentric chromosomes |
| 9 | [09_mammals_t2t.png](09_mammals_t2t.png) | T2T mammals | 143/157 CENP-A centromeres from the FASTA alone, 148 with genes and Hi-C; the 4 confident errors shown with their raw data |
| 10 | [10_retrotransposon_centromeres.png](10_retrotransposon_centromeres.png) | Satellite-free centromeres from retrotransposons | potato 2 → 6/11, *Ae. tauschii* 1 → 5/7, *Chlorella* 0 → 8/10, moss with Hi-C 25/26; satellite centromeres stay right |
| 11 | [11_grasses_and_crops.png](11_grasses_and_crops.png) | Crops: which evidence places each centromere | satellites, CRM cores or Hi-C's arm fold, chromosome by chromosome; tomato fails with all |
| 12 | [12_hic_moves_calls.png](12_hic_moves_calls.png) | Tested Hi-C overrules a wrong satellite line | barley: the arm fold and Hi-C agree on CENH3, 0 → 6/7 right |
| 13 | [13_fungi_and_protists.png](13_fungi_and_protists.png) | Point, regional and AT-rich centromeres from Hi-C alone | yeasts, *Plasmodium*, *Neurospora* 1 → 39/40; five held-out fungi 0 → 43/72; no confident error |
| 14 | [14_rabl_axis.png](14_rabl_axis.png) | The Rabl axis in trans | where centromeres cluster at one pole, arms of different chromosomes align by distance from the centromere; absent in G1, territorial, holocentric and randomized nuclei |
| 15 | [15_cpg_kinetochore.png](15_cpg_kinetochore.png) | The kinetochore inside the array | a CpG methylation dip places it within 50 kb of CENP-A on 35 of 39 ape centromeres |
| 16 | [16_controls.png](16_controls.png) | Where any confident call is wrong, cenfuse makes none | holocentric genomes: 0 confident calls of 46 (Hi-C used untested would make 21); randomized Hi-C and misplaced methylation pass no test |
| 17 | [17_robustness.png](17_robustness.png) | Thinner data, same answers | Hi-C, CpG and LTR annotation thinned: a failed test falls back to the calls made without that data, never worse |
| 18 | [18_beyond_centromeres.png](18_beyond_centromeres.png) | One map of the chromosome | telomeres, interstitial telomeric arrays (the human 2q13 fusion at base-pair resolution), rDNA, knobs, completeness |
| 19 | [19_karyotype_morphology.png](19_karyotype_morphology.png) | Karyotypes from the sequence | acrocentric vs bi-armed matches literature karyotypes on 117 of 118 right calls |
| 20 | [20_where_it_fails.png](20_where_it_fails.png) | Where it fails, and what would settle it | 17 confident errors in 4 kinds, each with the data that would settle it |
| 21 | [21_runtime_memory.png](21_runtime_memory.png) | Runtime and memory | run time follows the satellite DNA in a genome, not its size; everything after the satellite search takes minutes |

Each image's build script, data and legend pages (a short legend that reads as speaker notes, and a long one) are
kept with the benchmark runs, outside the repository.
