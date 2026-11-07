# Outputs, options and running cenfuse

The reference for every file, column and flag cenfuse writes, and for running it. How the rules work:
[design.md](design.md). Coordinates are 1-based and inclusive unless said otherwise; `.` = none or not applicable.

## Files (`GENOME_DIR/cenfuse/`)

| file | content |
|---|---|
| `asm.centromeres.tsv` | one row per chromosome: the call and everything behind it (columns below) |
| `asm.centromeres.bed`, `.gff3` | the calls; the GFF3 also holds each CpG kinetochore (`regional_centromere_central_core`, Name=kinetochore) |
| `cenfuse_report.md` | per genome, in words: each line's result and test, what the centromeres are made of, the calls that are not `high` with what would settle them, morphology and completeness, the other features and their 3D association, notes on the inputs |
| `cenfuse_summary.tsv` | one row per genome (columns below); warnings are in the log and the report |
| `cenfuse_overview.pdf`, `asm.overview.png` | the overview figure with its legend pages; its first page at 300 dpi |
| `asm.loci.tsv`, `asm.features.gff3` | every feature cenfuse found (below), with each line's signal there |
| `asm.feature_signatures.tsv` | per kind of feature, each line's median signal |
| `asm.feature_contacts.tsv` | with Hi-C: how strongly each kind of feature contacts itself in 3D, against chance |
| `asm.hic_anchored.bedgraph` | the anchored Hi-C centromere score per bin, for a genome browser |
| `asm.line_tests.tsv` | do the lines agree with each other? One row per pair of lines: chromosomes where both have a locus, how many agree, how many chance would give, P, a verdict. It changes no call |
| `asm.te_signature.tsv` | with `tes/` from LTRquest: the LTR signature per chromosome: the lineage, the seed (lineage or family) and its fold enrichment, the verdict, the core, its elements, the families recruited for this chromosome, and the age of the signature's elements inside and outside the cores |
| `asm.provenance.tsv` | which plant each input comes from, and whether it decides ([plants.md](plants.md)) |
| `asm.plants.tsv`, `asm.plant_sites.tsv`, `asm.plants.png` | with several plants: each plant's verdict per chromosome and per site ([plants.md](plants.md)) |
| `satellite/` | the satellite engine: families (`satellite_families.tsv`, `.fa`, `.hmm`), arrays (`asm_satellite_arrays.gff`), its own calls (`asm_centromere_summary.tsv`), higher-order repeats, its overview, and `dropped_clusters.tsv`: the seed clusters of ≥ 20 kb left out by the cap of 40 families (largest first; span, period, seeds, chromosomes, the largest seed and its consensus) |
| `hic/`, `methylation/` | the Hi-C and methylation engines' tables and plots |
| `work/` | checkpoints (a rerun resumes from them) |

## `asm.centromeres.tsv`

**The call**

| column | meaning |
|---|---|
| `chrom`, `chrom_len` | the chromosome and its length |
| `start`, `end`, `length` | the final call |
| `call_source` | where the final call comes from: `sequence` (satellite or TE line), `hic`, `composition` (an AT trough) or `none` |
| `relation` | how the final call relates to the sequence call (below) |
| `confidence` | `high`, `medium`, `low` (unsure) or `none` |
| `support`, `against` | the lines that agree with the call / point elsewhere: `satellite` (≥ 5 kb of centromeric arrays inside the call), `te`, `hic`, `at`, `meth`, `fold`, `align`, `age` (younger centromeric LTR copies at a moved call) |
| `flags` | the integration's flags (below) |
| `offset_bp` | distance between the final and the sequence call (0 when they overlap) |
| `basis` | what the call stands on, in words |
| `caveat` | why it may be wrong, with where the other evidence points; after `->`, the input that would settle it. `.` = nothing known against it |

**The sequence line** (satellites and LTR retrotransposons)

| column | meaning |
|---|---|
| `seq_start`, `seq_end`, `seq_confidence`, `seq_flags` | the satellite engine's own call, confidence and flags |
| `seq_evidence` | what the sequence call stands on: `satellite`, `weak_satellite`, `long_unit_tandem` or `retrotransposon`. With the flag `long_unit_unsupported` the array is a feature, not the call |
| `main_family` | the satellite family of the call |
| `seq_spacing` | the median distance between neighboring copies of the main family in the call's arrays (bp); `basis` names it where it differs from the family's unit, since a profile can fit shorter copies with a gap |
| `seq_alt_start`, `seq_alt_end` | a competing cluster: the place tested Hi-C may move the call to |
| `seq_overruled_start`, `seq_overruled_end` | the satellite array a decisive LTR core replaced (reported only; Hi-C never moves a call there). `.` also when the satellite outputs were made before cenfuse 0.9 and read with `--satellite-dir` |
| `seq_rejected_family`, `seq_rejected_start`, `seq_rejected_end` | a family the LTR line rejected at classification (its arrays lie away from the decisive cores on most chromosomes), and its largest cluster on this chromosome, outside the call: the array the call would otherwise have stood on |
| `te_clade`, `te_seed`, `te_start`, `te_end` | the LTR marker's lineage, the lineage or family the signature was seeded by, and the core on this chromosome |
| `te_share`, `te_chrom_elements` | the share of this chromosome's signature elements inside the core, and how many there are: a lineage other than CRM needs most of them inside to act |
| `te_core_elements`, `te_call_elements` | the two counts of the head-to-head test: signature elements in the core, and in the call or the array it replaced |

**AT troughs** (only where most chromosomes have no sequence call)

| column | meaning |
|---|---|
| `at_gate` | the genome's test: `pass` (the deepest GC troughs fall one per chromosome), `pass_hic` (they agree with Hi-C's prior-free peaks on most chromosomes), `fail`, or `untested` (too few chromosomes) |
| `at_start`, `at_end`, `at_gc` | this chromosome's deepest trough at the window the test chose, and its GC fraction |

**Hi-C**

| column | meaning |
|---|---|
| `hic_gate`, `hic_gate_k`, `hic_gate_n`, `hic_gate_p` | the leave-one-out test for this chromosome: `pass` (Hi-C's peaks agree with the other chromosomes' calls on most of them, P ≤ 0.01), `support` (beyond chance on at least 3, but not most), `fail` or `untested` (fewer than 3 reference calls); k of n agree, with P |
| `hic_gate_refs` | the reference calls the test used: `reference` (the confident calls), `low` (no confident call anywhere, so the `low` calls at clean peaks), or `mixed` (fewer than 4 confident calls, joined by the `low` calls at clean peaks) |
| `hic_end_k`, `hic_end_n` | which end: of the other chromosomes with a call in a terminal eighth, how many score higher at the end holding the call |
| `cis_peak_start`, `cis_peak_end` | the extreme of the contact crossing each bin, where the cis signal passed its test (against the LTR cores, or without them by Hi-C's trans contacts) |
| `cis_test` | what tested the cis signal: `ltr` (the LTR cores) or `trans` (Hi-C's trans contacts anchored on the cis extremes, read outside the terminal eighths) |
| `age_peak`, `age_call`, `age_p`, `age_premise_p` | with `tes/` (LTRquest K2P_T) and Hi-C, only where age can break a tie (a `high` sequence call and a strong peak of validated Hi-C more than a bin from it): the median insertion age in years (and the number of copies) of the chromosome's signature families at Hi-C's anchored core and at the sequence call, each widened by 4 bins; the one-sided Mann-Whitney P that the copies at the peak are younger; and the one-sided Fisher P that, on the other chromosomes, the signature's younger copies lie in the cores more often than its older ones (the premise; both must be ≤ 0.01 for `moved_with_age`) |
| `hic_anch_peak_*`, `hic_anch_core_*`, `hic_anch_level`, `hic_anch_flags` | the anchored Hi-C peak (contacts with the other chromosomes' calls): position, core, strength (`high`, `medium`, `low`) and flags (`hub`, `rDNA`, `telomeric`) |
| `hic_anch_agrees_seq` | whether the anchored peak agrees with the sequence call |
| `hic_denovo_peak_*`, `hic_denovo_confidence`, `hic_denovo_flags` | Hi-C's own prior-free peak (no sequence call used) |
| `fold_gate`, `fold_point`, `fold_core_*` | the arm fold (the chromosome's own contact map folding at the centromere): the genome's test result, and this chromosome's fold point and core |
| `align_gate`, `align_point`, `align_core_*` | the trans arm alignment (arms of different chromosomes aligned by distance from their centromeres): the same |

**Methylation**

| column | meaning |
|---|---|
| `meth_seq`, `meth_hic`, `meth_call` (+ `_score`) | the CHG line's label at the sequence call, at the Hi-C peak and at the final call: `centromere_like`, `dip`, `no_signature_data` or `no_data` |
| `meth_samples` | the methylation samples that decided |
| `cpg_gate`, `cpg_gate_k`, `cpg_gate_n`, `cpg_gate_p` | the CpG test for this chromosome: do the other chromosomes' calls hold a CpG dip more often than methylation from elsewhere would put one there? |
| `cpg_kt_start`, `cpg_kt_end` | the kinetochore: the call's CpG dip (shown where the test passes and the dip is significant). It refines the call; it does not check it (an inactive array has dips too) |
| `cpg_dip_p` | the dip is longer than in the call's own 5-kb bins shuffled |
| `cpg_pos_p` | the dip is longer than methylation from elsewhere on the chromosome makes on the same bins (the test counts these) |
| `cpg_kt_meth`, `cpg_array_meth` | mean methylation in the kinetochore; median over the array |

**Several plants**: `hic_library`, `hic_provenance`, `plants_shared`, `plants_shifted`, `plants_absent`,
`plants_unresolved` ([plants.md](plants.md)).

**Morphology and completeness**

| column | meaning |
|---|---|
| `cen_point` | the centromere position the arms are measured from: the kinetochore's midpoint, else the call's |
| `short_arm` | the sequence end nearer the centromere (`.` when the call spans the middle) |
| `ci`, `arm_ratio` | centromeric index (short arm / length) and arm ratio (long / short) |
| `morphology` | Levan class: metacentric (arm ratio ≤ 1.7), submetacentric (≤ 3), subtelocentric (≤ 7), acrocentric; two classes when the call spans both. Measured from the sequence, not cytogenetic appearance |
| `telomeres` | telomere arrays at `both` ends, one (`start`, `end`) or `none` |
| `gaps`, `call_gaps` | assembly gaps (runs of ≥ 10 N) on the chromosome and inside the call |
| `t2t` | `yes` with telomere arrays at both ends and no gap |

## Relation and confidence

| `relation` | meaning |
|---|---|
| `concordant` | Hi-C or methylation agrees with the sequence call |
| `offset` | tested Hi-C (with a second line for a `high` call), or Hi-C with an LTR core, moved the call; the old one is in `seq_start`/`seq_end` |
| `contested` | a tested line points elsewhere without enough support to move the call: kept, one level lower |
| `functional_only` | no sequence call; Hi-C placed it |
| `sequence_only` | only the sequence line was usable (no Hi-C, Hi-C failed its test, or its peak was weak), including a call placed by an AT trough alone |
| `none` | no call |

The satellite engine sets the starting confidence. Each agreeing tested line raises a kept call one level (two at
most); each opposing line lowers it one. A call that Hi-C moved or placed is `medium`, `high` with a second line.

## Flags

| flag | meaning |
|---|---|
| `competing_cluster` | a second satellite cluster at least half the size, or a second satellite that also marks the chromosome end in an acrocentric genome: the call is `low`, the other is in `seq_alt_*` |
| `hic_chose_alt` | tested Hi-C chose the competing cluster |
| `hic_alt_ambiguous` | no mappable DNA separates the call from its competing cluster, so Hi-C cannot tell them apart: no Hi-C support |
| `acrocentric` | the family sits at one end of most chromosomes (an acrocentric karyotype) |
| `near_end`, `low_satellite`, `diffuse`, `no_dominant_satellite` | why the satellite engine lowered its confidence: the call is near a chromosome end; it holds little satellite; it is a broad satellite-sparse block; no family is the largest centromeric cluster on most chromosomes (every call `low`) |
| `not_one_per_chromosome` | the largest clusters of the centromeric families do not fall one per chromosome in this genome, so those arrays may not be centromeres: `low`, unless an LTR core agrees |
| `te_concordant`, `te_discordant` | the LTR marker's core agrees with the call / lies elsewhere; with a decisive core the call is raised to at least `medium` / lowered one level (otherwise they only flag) |
| `weak_core` | a decisive signature's core too thin to be significant took a chromosome without a call, or a satellite call of a family that failed the one-per-chromosome test: `low` |
| `satellite_away_from_core`, `satellite_beside_core`, `weak_satellite_beside_core` | a decisive core replaced the satellite call (kept in `seq_overruled_*`): a call elsewhere that was `low` or held significantly fewer of the signature's elements; a call beside the core that held significantly fewer; a weak call beside it |
| `partial_chain` | the main family is no tandem satellite: less than half of its bp in tandem repeats and less than half in whole-unit hits. The sequence call is `low`, or at most `medium` where a retrotransposon core overlaps it (for a lineage other than CRM, only if the family is not made of it); lines that pass their tests may raise the final call |
| `chain_to_core` | such a call narrowed to the decisive core inside it |
| `te_derived`, `gene_rich`, `gaps:N`, `contig_end` | the main family lies mostly inside LTR retrotransposons; the call is more gene-dense than its chromosome; N assembly gaps in the call; the call reaches a contig end (`--classify-from`) |
| `hic_offset` | strong tested Hi-C points elsewhere but cannot move a `high` call without a second line (`contested`) |
| `hic_disagrees_weak`, `hic_near` | a tested peak elsewhere too weak to oppose; a peak beside the call but within Hi-C's resolution (no effect) |
| `hic_support` | Hi-C passed only its support test: its own significant peak agrees, so the call rose one level (on this tier Hi-C never opposes or moves; with `hic_chose_alt`, a strong peak chose the competing cluster) |
| `hic_end` | the call lies in a terminal eighth and the anchored score is higher at its end than at the other, as on most such chromosomes of this genome: the call rose one level |
| `te_hic_pair` | anchored on the LTR cores, Hi-C's peak lies at this chromosome's core, the two agreeing beyond chance genome-wide: a call there rose one level, or a `low` call (with the test passed on most chromosomes, any call) moved there |
| `cis_agrees`, `cis_near`, `cis_only`, `cis_elsewhere` | the extreme of the contact crossing each bin, tested against the LTR cores (or, without them, by Hi-C's trans contacts on most chromosomes; column `cis_test`): at the call (one level up); within 4 bins of it (neither way); farther, a `low` call moved to it (still `low`), or a chromosome without a call placed there; farther, the call kept, one level lower if the test passed on most chromosomes (`contested`). The trans-tested signal is read outside the terminal eighths and leaves a call inside one as it is |
| `at_only`, `at_pair`, `at_elsewhere` | a call placed by an AT trough alone; a Hi-C call on a chromosome without a sequence call, with the trough agreeing: it takes the trough's interval, one level up; the trough lies elsewhere |
| `at_denovo_pair` | a chromosome nothing else places: its AT trough and Hi-C's prior-free peak (clean, placing no call of its own) agree, in a genome where the two agree beyond chance: the trough, `low` |
| `cis_blind` | the call lies in a terminal eighth, where the cis signal tested by trans contacts cannot see; its extreme lies elsewhere (the caveat says where), and the call is left as it is |
| `fold_pair`, `align_pair` | Hi-C's prior-free peak and the arm fold (or the trans arm alignment) agree here, away from the sequence call: the call moved, or was placed, at `medium` |
| `fold_offset`, `align_offset`, `hic_denovo_offset` | the arm fold / the alignment / Hi-C's prior-free peak points elsewhere: one level lower each |
| `moved_with_te`, `_meth`, `_fold`, `_align`, `_age` | a `high` call moved because a second line (TE core, methylation, arm fold, alignment; or copies of the centromeric LTR families younger at the peak than at the call, in a genome whose young copies mark its cores) agreed with Hi-C |
| `hic_medium_peak` | no sequence call; a clean medium-strength Hi-C peak placed a `low` call |
| `hic_denovo_only` | too few sequence calls to test Hi-C: Hi-C's own clean peak placed this `low` call (untested) |
| `hic_denovo_elsewhere` | in such a genome, Hi-C's clean peak lies away from the kept sequence call (the caveat says where) |
| `long_unit_unsupported` | the only sequence candidate is a long-unit tandem array that neither Hi-C, a retrotransposon core nor methylation supports: it is reported (`seq_start`, `seq_end`, the loci) and is not a call |

## `cenfuse_summary.tsv` (one row per genome)

| columns | meaning |
|---|---|
| `chromosomes`, `lines`, `called`, `high`, `medium`, `low` | chromosomes analyzed, the lines given, calls by confidence |
| `concordant` ... `none`, `seq_called`, `partial_chain_calls` | calls by relation; chromosomes with a sequence call; of those, calls on chains of partial hits |
| `sat_one_per_chromosome` | do the N largest clusters of the centromeric families fall one per chromosome? |
| `te_marker`, `te_decisive`, `te_seed`, `te_specificity`, `te_verdict` | the LTR marker lineage, whether it is decisive, the seed and its fold enrichment, and why |
| `at_gate`, `at_window`, `at_p`, `at_hic` | the AT-trough line: its test, the window it chose, P, and the agreement with Hi-C's prior-free peaks |
| `telomere_unit`, `telomere_ends`, `telomere_interstitial` | the telomere repeat found, ends with a telomere array, arrays inside chromosomes |
| `left_out_telomeric` | sequences left out of the chromosome set that hold telomere arrays at both ends (complete molecules: a small chromosome the length rule cut, a B chromosome); flagged in the report, never added |
| `rdna_arrays`, `rdna_chromosomes`, `rdna5s_arrays`, `rdna5s_chromosomes` | 45S (NOR) and 5S rDNA arrays |
| `karyotype`, `t2t_chromosomes`, `calls_with_gaps`, `short_arm_at_end` | Levan classes of the `high` and `medium` calls (`low` calls and chromosomes without a call counted apart); complete chromosomes; calls holding gaps; confident calls whose short arm is at the sequence end (the convention puts it first) |
| `hic_bin`, `hic_chromosomes` | the Hi-C bin size; chromosomes Hi-C analyzed |
| `hic_gate_pass`, `hic_gate_support`, `hic_gate_untested`, `hic_gate_refs` | chromosomes passing the Hi-C test, on the support tier, untested; which references it used |
| `hic_anch_on_seq`, `hic_denovo_on_seq`, `hic_denovo_on_anch`, `hic_seq_called` | agreement counts between anchored peaks, prior-free peaks and sequence calls |
| `hic_sequence_conflict` | Hi-C's prior-free peaks are confident but disagree with the sequence calls, and the test fails: one line is wrong genome-wide |
| `te_hic`, `cis_line` | where the Hi-C test failed: with `tes/`, Hi-C's peak at the LTR cores, and the cis signal at the LTR cores (k/n, P, verdict); without LTR cores, the cis signal tested by trans contacts (k/n, P, verdict, and the chromosomes it acts on; `ambiguous` when both signs pass) |
| `hic_rabl`, `fold_gate`, `fold_on_denovo`, `fold_p`, `align_gate`, `align_on_denovo`, `align_p`, `hic_axis` | the Rabl-axis tests: the arm fold and the trans arm alignment (k/n chromosomes whose prior-free peak lies in the core, P), and the axis that acted (`fold`, `align`) |
| `hic_3d_centromeres`, `hic_3d_ends`, `hic_architecture` | how strongly the calls and the chromosome ends contact each other in 3D; what the library showed: `rabl_axis`, `centromere_clustering` (test passed), `weak_centromere_clustering` (support tier), `untested`, `not_detected` |
| `meth_*`, `cpg_*` | the CHG line (samples, anchors, signature, deciding samples) and the CpG line (chromosomes passing, dips, kinetochores) |
| `inputs_reference`, `inputs_other`, `inputs_unknown`, `hic_library`, `plants*` | the inputs by provenance verdict, the deciding Hi-C library, the plants' verdicts ([plants.md](plants.md)) |

## Other features (`asm.loci.tsv`, `asm.features.gff3`)

Every candidate interval of every line (satellite clusters ≥ 20 kb, TE cores, Hi-C cores and hubs, methylation
candidates, 45S and 5S rDNA arrays, telomere arrays at the ends and inside chromosomes, the sequence and final calls)
becomes a member; overlapping members merge into a locus, which gets every line's value there and the first kind that
fits: `centromere` > `sequence_centromere` (where a call moved from) > `telomere` > `interstitial_telomeric` > `rDNA`
> `subtelomeric_block` > `interstitial_block` (knobs) > `long_unit_array` (an array that no line supports, so no call) >
`alt_centromere_candidate` > `hic_hub` > `hic_peak` > `te_core` > `satellite_block` > `meth_candidate` > `other`. GFF3
types are Sequence Ontology terms (`gt gff3validator
-typecheck so` passes); cenfuse's kind is the `kind=` attribute, and each member names its locus (`locus=`) and type
(`member=`).

## Which sequences are analyzed

An assembly also holds unplaced scaffolds, and every line needs whole chromosomes. cenfuse takes the sequences above
the first tenfold drop in length between neighbors (looked for once half the assembly is counted; all sequences when
there is none). A tapering tail of scaffolds has no such drop, so it is cut: from the first sequence after half the
assembly that is shorter than a tenth of the median length before it, if the sequences from there on span a tenfold
range. The log and the report's Notes say what set the chromosomes and how many sequences were left out.

The rule cannot tell a short chromosome from a scaffold by length. It drops a dot chromosome that stands alone after
a tenfold gap (*D. melanogaster* chromosome 4 without the Y), it cuts dot chromosomes followed by a tapering tail of
scaffolds, and it keeps large unplaced scaffolds (31 extra sequences in 12 of 121 polyploid assemblies with known
chromosome lists). For any of these, list the chromosomes in `GENOME_DIR/asm.chromosomes.txt`: one sequence name per
line, `#` starts a comment line (a `#` inside a name, as in `HG002#1#chr1`, is part of it); blank lines, Windows line
ends and trailing spaces are fine. List every chromosome: the CHG line sets its expected domain size from the number
listed. An unknown or repeated name, or an empty list, stops the run. `--classify-from` targets (contig assemblies)
keep the plain tenfold rule unless listed.

## Options

| option | meaning |
|---|---|
| `-o DIR` | output folder [`GENOME_DIR/cenfuse`] |
| `-t N` | threads [8] |
| `--satellite-dir DIR` | take the satellite outputs of an earlier run (`satellite/`) as they are, instead of running the satellite engine; cenfuse only warns if their chromosomes differ |
| `--seed-engine fastan\|trf` | the tandem-repeat finder that seeds the satellite families [fastan] |
| `--seed-cache DIR` | a folder of tandem-repeat results shared between runs |
| `--classify-from T=S[+S2]` | genome T (e.g. a contig assembly) takes the centromeric families of genome(s) S |
| `--force STAGE` | recompute STAGE and every later one: `satellite` (the whole satellite engine), `telomere`, `rdna`, `rdna5s`, `gaps`, `attrough`, `genotype`, `hic`, `fold`, `karyotype`, `meth`, `cpg`, or a satellite-engine stage (`seeds`, `families`, `annotate`, `classify`, `centromeres`, `hor`, `output`, `plot`) |
| `-v` | per-step progress and sanity checks |

`--force` is needed when the run stops and names a stage (after an upgrade that changed it), after adding or changing
`genes/` or `tes/` in an earlier run's folder (`--force classify`: the satellite engine tracks its code and
parameters, not its input files), and when the chromosome set changed (below).

## Running

**Resuming.** A rerun in the same output folder resumes from checkpoints. Each stage is keyed on its inputs and code,
so an upgrade recomputes only the stages it changed. The satellite engine's checkpoints also stand on the chromosome
set, recorded in `satellite/work/chromosome_set.json`: a rerun whose set changed (a genome added or removed, a new or
edited `chromosomes.txt`) stops, names the genomes, and says to rerun with `--force satellite`.

**The log.** Without `-v` the log reports the start, milestones, warnings and the end. The genome-wide satellite
search logs the share of its work done and a rough estimate of the time left at every 10%, and after an hour without
a line; when it ends, it names the three families that cost the most CPU (all families:
`satellite/work/annotate/<genome>.search_cost.tsv`) and, before it starts, how many models are searched with their
low-complexity positions neutral.

**Reproducibility.** On Linux, fresh runs give identical results apart from the dates and commands in file headers
(FasTAN runs with address-space randomization off, `setarch -R`). Two kinds of work files can differ: the search
cost table (CPU times), and the families stage's BLAST scratch under `satellite/work/families/`. Without `setarch`
(e.g. macOS) cenfuse warns, and about 1% of FasTAN's array cuts can differ between runs.

**Memory.** The satellite search holds about one search workspace (up to ~2.4 GB for 2-kb units) per thread.

**FasTAN crashes.** FasTAN can abort on a large satellite array (`malloc(): corrupted top size`). A piece of the seed
search on which it dies by a crash signal is split into two overlapping halves and retried, down to 20 kb; a piece
that still crashes is skipped with a warning, and the crashes are listed in `<seed cache file>.fastan_crashes.tsv`
(linked under `satellite/work/seeds/`). If FasTAN crashes on every piece of a chunk, the run stops and suggests
`--seed-engine trf`. A kill (out of memory, the scheduler, Ctrl-C) still stops the run.
