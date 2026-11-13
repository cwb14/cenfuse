# Outputs, options and running cenfuse

The reference for every file, column and flag cenfuse writes, and for running it. How the rules work:
[design.md](design.md). Coordinates are 1-based and inclusive unless said otherwise; `.` = none or not applicable.

## Files (`GENOME_DIR/cenfuse/`)

| file | content |
|---|---|
| `asm.centromeres.tsv` | one row per chromosome: the call and everything behind it (columns below) |
| `asm.centromeres.bed`, `.gff3` | the calls; the GFF3 also holds each CpG kinetochore (`regional_centromere_central_core`, Name=kinetochore; dip_depth, dip_z, pos_p, kt_meth, array_meth) |
| `cenfuse_report.md` | per genome, in words: each line's result and test, what the centromeres are made of, the calls that are not `high` with what would settle them, morphology and completeness, the other features and their 3D association, notes on the inputs |
| `cenfuse_summary.tsv` | one row per genome (columns below); warnings are in the log and the report |
| `run_info.json` | the run manifest: cenfuse's version and git commit, the command, host, threads and times, Python and package versions, every external tool's path, version and conda package, a content fingerprint of every input, how the run ended (`status`: `running`, `done`, `failed` with the error; `incomplete` outputs), and what it cost (`cost`: `wall_minutes`, `cpu_hours` of cenfuse and every process it started, `peak_memory_gb` of the whole process tree, sampled every 30 s, and `memory`, what that measured; also the run's last log line) |
| `.cenfuse.lock` | present while a run uses the folder (host, pid, start, command); removed at the end |
| `cenfuse_overview.pdf`, `asm.overview.png` | the overview figure with its legend pages; its first page at 300 dpi |
| `asm.loci.tsv`, `asm.features.gff3` | every feature cenfuse found (below), with each line's signal there |
| `asm.feature_signatures.tsv` | per kind of feature, each line's median signal (columns below) |
| `asm.feature_contacts.tsv` | with Hi-C: how strongly each kind of feature contacts itself in 3D, against chance (columns below) |
| `asm.hic_anchored.bedgraph` | the anchored Hi-C centromere score per bin, for a genome browser |
| `asm.line_tests.tsv` | do the lines agree with each other? One row per pair of lines: chromosomes where both have a locus, how many agree, how many chance would give, P, a verdict (columns below). It changes no call |
| `asm.te_signature.tsv` | with `tes/` from LTRquest: the LTR signature per chromosome: the lineage, the seed (lineage or family) and its fold enrichment, the verdict, the core, its elements, the families recruited for this chromosome, and the age of the signature's elements inside and outside the cores (columns below) |
| `asm.provenance.tsv` | which plant each input comes from, and whether it decides ([plants.md](plants.md); [columns](plants.md#asmprovenancetsv-one-row-per-input)) |
| `asm.plants.tsv`, `asm.plant_sites.tsv`, `asm.plants.png` | with several plants: each plant's verdict per chromosome and per site ([plants.md](plants.md); columns: [plants](plants.md#asmplantstsv-one-row-per-plant-and-chromosome), [sites](plants.md#asmplant_sitestsv-one-row-per-plant-and-site)) |
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
| `confidence` | `high`, `medium`, `low` (unsure) or `none` (no call) |
| `support`, `against` | the lines that agree with the call / point elsewhere: `satellite` (≥ 5 kb of centromeric arrays inside the call), `te`, `hic`, `at`, `meth`, `fold`, `align`, `age` (younger centromeric LTR copies at a moved call), `cpg` (against only: the call's arrays hold no CpG dip while another array does) |
| `flags` | the integration's flags ([Flags](#flags)) |
| `rule` | the rule of the integration that set the call's place or relation (below); `sequence` = the sequence line's call as it is |
| `offset_bp` | distance between the final and the sequence call when the call was moved (0 when they overlap); `.` when the call is the sequence call or its competing cluster |
| `basis` | what the call stands on, in words |
| `caveat` | why it may be wrong, with where the other evidence points; after `->`, the input that would settle it. `.` = nothing known against it |

**The sequence line** (satellites and LTR retrotransposons)

| column | meaning |
|---|---|
| `seq_start`, `seq_end`, `seq_confidence`, `seq_flags` | the sequence line's own call, confidence and flags ([Flags](#flags)) |
| `seq_evidence` | what the sequence call stands on: `satellite`, `weak_satellite`, `long_unit_tandem` or `retrotransposon`. With the flag `long_unit_unsupported` the array is a feature, not the call |
| `main_family` | the satellite family of the call |
| `seq_spacing` | the median distance between neighboring copies of the main family in the call's arrays (bp); `basis` names it where it differs from the family's unit, since a profile can fit shorter copies with a gap |
| `seq_alt_start`, `seq_alt_end` | a competing cluster: the place tested Hi-C may move the call to |
| `seq_overruled_start`, `seq_overruled_end` | the satellite array a decisive LTR core replaced (reported only; Hi-C never moves a call there) |
| `seq_rejected_family`, `seq_rejected_start`, `seq_rejected_end` | a family the LTR line rejected at classification (its arrays lie away from the decisive cores on most chromosomes), and its largest cluster on this chromosome, outside the call: the array the call would otherwise have stood on |
| `te_clade`, `te_seed`, `te_start`, `te_end` | the LTR marker's lineage, the lineage or family the signature was seeded by, and the core on this chromosome |
| `te_share`, `te_chrom_elements` | the share of this chromosome's signature elements inside the core, and how many there are: a lineage other than CRM needs most of them inside to act |
| `te_core_elements`, `te_call_elements` | the two counts of the head-to-head test: signature elements in the core, and in the call or the array it replaced |

**AT troughs** (only where most chromosomes have no sequence call or a `low` one with a competing cluster)

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
| `hic_anch_peak_*`, `hic_anch_core_*`, `hic_anch_level`, `hic_anch_flags` | the anchored Hi-C peak (contacts with the other chromosomes' calls): position, core, strength (`high`, `medium`, `low`) and flags (`hub`, `rDNA`, `telomeric`, `null_untestable`; as for the prior-free peak) |
| `hic_anch_agrees_seq` | whether the anchored peak agrees with the sequence call |
| `hic_denovo_peak_*`, `hic_denovo_confidence`, `hic_denovo_flags` | Hi-C's own prior-free peak (no sequence call used), its level and the Hi-C engine's flags: `hub` (2× coverage and 2× trans contacts: a contact hub, not a centromere), `rDNA` (at a 45S array or an rDNA-tagged satellite: NORs meet in the nucleolus), `telomeric` (in an end window: ends cluster too), `no_rabl` (the library shows no Rabl polarity, the engine's premise), `null_untestable` (too few chromosomes for the engine's null: judged on its other two tests, at most `medium`), `null_degenerate` (the null ran but could not score this chromosome: judged the same way, but never counted as tested), `not_converged` (the iterative estimate did not settle), `cycle_unstable` (it moves by more than a bin within its cycle), `no_contacts` (the chromosome has no contact with the others' calls). A prior-free peak with `hub`, `rDNA`, `telomeric`, `not_converged` or `cycle_unstable` is not clean, so it never places a call alone |
| `fold_gate`, `fold_point`, `fold_core_*` | the arm fold (the chromosome's own contact map folding at the centromere): the genome's test result, and this chromosome's fold point and core |
| `align_gate`, `align_point`, `align_core_*` | the trans arm alignment (arms of different chromosomes aligned by distance from their centromeres): the same |

**Methylation**

| column | meaning |
|---|---|
| `meth_seq`, `meth_hic`, `meth_call` (+ `_score`) | the CHG line's label at the sequence call, at the Hi-C peak and at the final call: `centromere_like`, `dip`, `no_signature_data` or `no_data` |
| `meth_samples` | the methylation samples that decided |
| `cpg_gate`, `cpg_gate_k`, `cpg_gate_n`, `cpg_gate_p` | the CpG test for this chromosome: do the other chromosomes' calls hold CpG dips deeper than methylation from elsewhere would put there, beyond chance? `k` of `n` references with `cpg_pos_p` ≤ 0.05; `p` combines all `n` P values (Fisher's method) |
| `cpg_kt_start`, `cpg_kt_end` | the kinetochore: the deepest CpG dip (CDR) in the call's centromeric arrays, shown where the test passes; the deepest spread dip where no dip falls below half its surroundings. It refines the call |
| `cpg_dip_depth` | the kinetochore candidate's methylation over its surroundings' (a level dip falls below 0.5; 1 = no dip; `.` = fewer than 50 kb of array with data) |
| `cpg_dip_z` | the same dip's methylation below its surroundings', in robust SD of the surroundings' own variation (a spread dip falls below -3 and at least 34% below its surroundings; 0 = no dip). Unlike `cpg_dip_depth`, it does not change when a track encodes methylation on another scale (a PacBio bigWig, nanopolish) |
| `cpg_pos_p` | the call's arrays hold dips deeper, or a valley wider, than methylation from elsewhere on the chromosome makes on the same array bins (level depth, spread depth and the most prominent valley's half-depth width, combined by minP over the same rotations; the test counts these). A call without a dip can have a small P, from its valley's width: it gets no kinetochore |
| `cpg_kt_meth`, `cpg_array_meth` | mean methylation in the kinetochore; median over the call's arrays |
| `cpg_array_data` | the share of the call's array bins with CpG data (a low share is a coverage hole: tracks piled up at MAPQ 60 lose most reads in arrays) |
| `cpg_alt_start`, `cpg_alt_end` | where the test passes and the call's arrays have data but no dip and no `cpg_pos_p` ≤ 0.05: the CDR-like dip on another satellite array of the chromosome (deeper than methylation from elsewhere makes there). The call is contested (`cpg_elsewhere`) and kept |

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

A call moved to its competing cluster (`rule` `hic_chose_alt` or `at_chose_alt`) keeps the relation `concordant` or
`sequence_only`: its `start` differs from `seq_start`.

The satellite engine sets the starting confidence. Each agreeing tested line raises a kept call one level (two at
most); each opposing line lowers it one. A call that Hi-C moved or placed is `medium`, `high` with a second line.

**`rule`.** `sequence` (kept as the sequence line made it); `hic_agrees` (validated Hi-C agrees: concordant); `hic_support` (the
support tier agrees); `hic_end` (the which-end test); `hic_chose_alt` (Hi-C chose the competing cluster); `hic_moved` (a strong
validated peak moved the call); `hic_offset` (it points elsewhere and cannot move a `high` call: contested); `hic_strong_peak`,
`hic_medium_peak`, `hic_denovo_only` (Hi-C placed a call where the sequence has none); `fold_pair`, `align_pair` (the Rabl-axis
line and the prior-free peak at one locus); `fold_offset`, `align_offset`, `hic_denovo_offset` (they oppose a kept call:
contested); `te_hic_pair` (Hi-C at the LTR core); `cis_agrees`, `cis_only`, `cis_elsewhere`, `cis_trans_pair`, `cis_end_pair` (the cis signal); `at_only`,
`at_denovo_pair`, `at_pair`, `at_chose_alt` (the AT trough); `cpg_elsewhere` (CpG finds the CDR on another array: contested); `meth_support` (methylation alone confirms the call); `none` (no call).

## Flags

Two columns hold flags. `seq_flags` holds the sequence line's flags: why the satellite engine set its confidence,
and what the LTR line did to its call. `flags` holds the integration's flags: what the other lines did. Several flags
are comma-separated.

**Sequence line (`seq_flags`)**

| flag | meaning |
|---|---|
| `competing_cluster` | a second satellite cluster at least half the size, or a second satellite that also marks the chromosome end in an acrocentric genome: the call is `low`, the other is in `seq_alt_*` |
| `acrocentric` | the family sits at one end of most chromosomes (an acrocentric karyotype) |
| `near_end`, `low_satellite`, `diffuse`, `no_dominant_satellite` | why the satellite engine lowered its confidence: the call is near a chromosome end; it holds little satellite; it is a broad satellite-sparse block; no family is the largest centromeric cluster on most chromosomes (every call `low`) |
| `not_one_per_chromosome` | the largest clusters of the centromeric families do not fall one per chromosome in this genome, so those arrays may not be centromeres: `low`, unless an LTR core agrees |
| `te_concordant`, `te_discordant` | the LTR marker's core agrees with the call / lies elsewhere; with a decisive core the call is raised to at least `medium` / lowered one level (otherwise they only flag) |
| `weak_core` | a decisive signature's core too thin to be significant took a chromosome without a call, or a satellite call of a family that failed the one-per-chromosome test: `low` |
| `minority_core` | a decisive signature's significant core holding under half of the signature's elements on the chromosome replaced a satellite call that was `low` on its own: `low` (the replaced array is in `seq_overruled_*`) |
| `te_age_core` | the chromosome's LTR core is the youngest of several comparable density peaks of the signature's elements (`te_signature.tsv`, `core_age_from`) |
| `satellite_away_from_core`, `satellite_beside_core`, `weak_satellite_beside_core` | a decisive core replaced the satellite call (kept in `seq_overruled_*`): a call elsewhere that was `low` or held significantly fewer of the signature's elements; a call beside the core that held significantly fewer; a weak call beside it |
| `partial_chain` | the main family is no tandem satellite: less than half of its bp in tandem repeats and less than half in whole-unit hits. The sequence call is `low`, or at most `medium` where a retrotransposon core overlaps it (for a lineage other than CRM, only if the family is not made of it); lines that pass their tests may raise the final call |
| `chain_to_core` | such a call narrowed to the decisive core inside it |
| `te_derived`, `gene_rich`, `gaps:N`, `contig_end` | the main family lies mostly inside LTR retrotransposons; the call is more gene-dense than its chromosome; N assembly gaps in the call; the call reaches a contig end (`--classify-from`) |

**Integration (`flags`)**

| flag | meaning |
|---|---|
| `hic_chose_alt` | tested Hi-C chose the competing cluster |
| `hic_alt_ambiguous` | no mappable DNA separates the call from its competing cluster, so Hi-C cannot tell them apart: no Hi-C support |
| `hic_offset` | strong tested Hi-C points elsewhere but cannot move a `high` call without a second line (`contested`) |
| `hic_disagrees_weak`, `hic_near` | a tested peak elsewhere too weak to oppose; a peak beside the call but within Hi-C's resolution (no effect) |
| `hic_support` | Hi-C passed only its support test: its own significant peak agrees, so the call rose one level (on this tier Hi-C never opposes or moves; with `hic_chose_alt`, a strong peak chose the competing cluster) |
| `hic_end` | the call lies in a terminal eighth and the anchored score is higher at its end than at the other, as on most such chromosomes of this genome: the call rose one level |
| `te_hic_pair` | anchored on the LTR cores, Hi-C's peak lies at this chromosome's core, the two agreeing beyond chance genome-wide: a call there rose one level, or a `low` call (with the test passed on most chromosomes, any call) moved there |
| `cis_agrees`, `cis_near`, `cis_only`, `cis_elsewhere` | the extreme of the contact crossing each bin, tested against the LTR cores (or, without them, by Hi-C's trans contacts on most chromosomes; column `cis_test`): at the call (one level up); within 4 bins of it (neither way); farther, a `low` call moved to it (still `low`), or a chromosome without a call placed there; farther, the call kept, one level lower if the test passed on most chromosomes (`contested`). The trans-tested signal is read outside the terminal eighths and leaves a call inside one as it is |
| `at_only`, `at_pair`, `at_elsewhere` | a call placed by an AT trough alone; a Hi-C call on a chromosome without a sequence call, with the trough agreeing (the troughs tested against the Hi-C readout that placed it: the prior-free peaks, or the cis extremes where those placed it): it takes the trough's interval, one level up; the trough lies elsewhere |
| `at_denovo_pair` | a chromosome nothing else places: its AT trough and Hi-C's prior-free peak (clean, placing no call of its own) agree, in a genome where the two agree beyond chance: the trough, `low` |
| `at_chose_alt` | the genome's AT troughs passed their test and this chromosome's deepest trough window lies on the competing cluster of a `low` sequence call and not on the call: the competing cluster, `low` unless another line agrees there; not where validated Hi-C can tell the two clusters apart. A heuristic on the evidence of one genome (design.md) |
| `cpg_elsewhere` | where the CpG test passes, the call's centromeric arrays hold no dip while another array holds a CDR-like one (`cpg_alt_*`): the call is kept, one level below what the agreeing lines give (`contested`; with another line against it, one level more). CpG never moves a call: in the truth genomes such a dip was the centromere once and elsewhere twice |
| `cis_blind` | the call lies in a terminal eighth, where the cis signal tested by trans contacts cannot see; its extreme lies elsewhere (the caveat says where), and the call is left as it is |
| `cis_trans_pair` | the call moved to the cis extreme, against the confident call it contradicts, where the cis signal passed two genome-wide tests (against the LTR cores and by Hi-C's trans contacts) and this chromosome's own trans contacts peak at the extreme; stays `low` (a rule shaped on one genome, tomato) |
| `cis_end_pair` | a `low` call in a terminal eighth moved to the interior cis extreme where this chromosome's own trans contacts peak (an acrocentric chromosome has no such pair and keeps its call); stays `low` (a rule shaped on one genome, tomato) |
| `fold_pair`, `align_pair` | Hi-C's prior-free peak and the arm fold (or the trans arm alignment) agree here, away from the sequence call: the call moved, or was placed, at `medium` |
| `fold_offset`, `align_offset`, `hic_denovo_offset` | the arm fold / the alignment / Hi-C's prior-free peak points elsewhere: one level lower each |
| `moved_with_te`, `_meth`, `_fold`, `_align`, `_age` | a `high` call moved because a second line (TE core, methylation, arm fold, alignment; or copies of the centromeric LTR families younger at the peak than at the call, in a genome whose young copies mark its cores) agreed with Hi-C; several are joined by `+` (`moved_with_te+age`) |
| `hic_medium_peak` | no sequence call; a clean medium-strength Hi-C peak placed a `low` call |
| `hic_denovo_only` | too few sequence calls to test Hi-C: Hi-C's own clean peak placed this `low` call (untested) |
| `hic_denovo_elsewhere` | in such a genome, Hi-C's clean peak lies away from the kept sequence call (the caveat says where) |
| `long_unit_unsupported` | the only sequence candidate is a long-unit tandem array that neither Hi-C, a retrotransposon core nor methylation supports: it is reported (`seq_start`, `seq_end`, the loci) and is not a call |

## `cenfuse_summary.tsv` (one row per genome)

| columns | meaning |
|---|---|
| `genome` | the genome, by the name of its FASTA (`asm`) |
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
| `meth_*`, `cpg_*` | the CHG line (samples, anchors, signature, deciding samples) and the CpG line (chromosomes passing, dips, kinetochores, calls contested by another array's dip, `cpg_array_data`: the arrays' share of bins with data over the genome's, below 0.8 a coverage hole) |
| `inputs_reference`, `inputs_other`, `inputs_unknown`, `hic_library`, `plants*` | the inputs by provenance verdict, the deciding Hi-C library, the plants' verdicts ([plants.md](plants.md)) |
| `incomplete` | outputs that failed for this genome (`report`, `figure`), or `.`: the tables are complete either way, the report says what is missing, and the log ends with a warning |

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

`asm.loci.tsv` has one row per locus:

| column | meaning |
|---|---|
| `locus` | the locus: `asm:chrom:L<n>`, numbered along the chromosome |
| `chrom`, `start`, `end`, `length` | its span, from the first member's start to the last member's end |
| `kind` | the first kind that fits (above) |
| `lines`, `n_lines` | the lines with a member in the locus, and how many: `satellite` (satellite clusters, rDNA arrays and the sequence call), `te`, `hic`, `meth`, `telomere`; the final call is not a line |
| `members` | every member, `;`-separated, as `type:label:start-end`: `integrated_call` (the final call; label: its relation), `sequence_call` (its evidence), `long_unit_array`, `satellite_cluster` (family[role]), `rdna_array` and `rdna5s_array` (`45S xN`, `5S xN`: N units), `telomere_array` and `telomere_interstitial` (the telomere unit), `te_core` (the lineage, starred for the marker), `hic_anchored_core` (`anchored:` and its level; `:strong` when `high` without a `hub`, `rDNA` or `telomeric` flag), `hic_denovo_core` (`denovo:` and its level), `hic_hub`, `meth_candidate` (sample:label; `:specific` when plant-specific) |
| `sat_bp_centromeric`, `sat_bp_cen_associated`, `sat_bp_interstitial`, `sat_bp_subtelomeric`, `sat_bp_dispersed`, `sat_bp_other` | bp of satellite arrays inside the locus, by the role of their family ([design.md](design.md#satellites-cenfusesat), step 4); `sat_bp_other`: the `rDNA`, `minor` and `unassessed` families |
| `rdna` | the rDNA there: the satellite engine's rDNA tags of its families (the Rfam rRNA model a family's consensus holds, such as `5S_rRNA`, or `45S_spacer` for arrays among 45S genes) and the 45S and 5S arrays (`45S_array:NxU`, `5S_array:NxU`: N units, U bp apart) |
| `telomere_bp` | bp of telomere-unit copies in its telomere arrays |
| `te_clades` | the LTR lineages whose core on this chromosome lies in the locus: the lineages with a significant core on most chromosomes, and the signature's lineage with the signature's core |
| `te_marker_elements` | elements of the chromosome's LTR signature (without one, of the marker lineage) in the bins the locus overlaps (bins of 1/500 of the chromosome, at least 1 kb) |
| `te_elements`, `te_top_clades` | LTR retrotransposons classified to a lineage whose midpoint lies in the locus, and the three most common lineages (`lineage:count`) |
| `hic_anch_max`, `hic_denovo_max`, `hic_telomere_max` | the highest Hi-C score in the locus' bins, over the chromosome's median: anchored (contacts with the other chromosomes' calls), prior-free, and telomere pole (contacts with the other chromosomes' end windows) |
| `hic_hub` | `True` when the locus holds a run of Hi-C hub bins (coverage and trans contacts both ≥ 2× the genome's median) |
| `meth_pooled_extreme` | the most extreme value, in the direction of the genome's centromeres (the lowest where they lose methylation), of the composition-adjusted CHG track in the locus' 50-kb bins, all samples pooled: the methylated share above or below each sequence class's own share in the flanks (1.5 to 4 expected centromere widths away) |
| `meth_raw_chg` | the share of methylated CHG calls in the locus' 50-kb bins, all samples pooled, not adjusted |
| `meth_samples_centromere_like`, `meth_samples_specific` | the samples with a methylation candidate in the locus that is `centromere_like` by the chromosome's leave-one-out signature, and that is plant-specific (with 3 or more samples: its value minus the other samples' median ranks among the most extreme 5% of the chromosome's windows); 0 without methylation |
| `gene_ratio` | with `genes/`: genes per bp in the locus over genes per bp on its chromosome (protein-coding genes, by midpoint) |
| `gap_bp` | N bases in the locus |
| `dist_to_call` | bp between the locus and its chromosome's final call, 0 when they overlap; `.` without a call |

## `asm.feature_signatures.tsv` (one row per kind of locus)

What each kind of locus in `asm.loci.tsv` looks like to every line: the median of its loci's values (`.` where none has
one), and, with Hi-C, the 3D self-association of the matching feature class (`asm.feature_contacts.tsv`). Kinds without
a locus have no row.

| column | meaning |
|---|---|
| `kind` | the kind, in the order of the kinds above |
| `loci`, `chromosomes` | how many loci are of this kind, and on how many chromosomes |
| `median_length` | their median length (bp) |
| `sat_bp_centromeric`, `sat_bp_interstitial`, `sat_bp_subtelomeric`, `telomere_bp` | the median of the loci's column of that name |
| `te_elements_per_mb`, `te_marker_per_mb` | the median per Mb of the loci's `te_elements` and `te_marker_elements` |
| `common_te_clade` | the lineage most often first in the loci's `te_top_clades` |
| `hic_anch_max`, `hic_telomere_max` | the median of the loci's column of that name |
| `hic_hub_share` | the share of the loci that hold Hi-C hub bins |
| `meth_pooled_extreme`, `meth_raw_chg`, `gene_ratio` | the median of the loci's column of that name |
| `gap_share` | the median share of N bases in the loci |
| `hic_3d_class`, `hic_3d_fold`, `hic_3d_p` | the feature class of `asm.feature_contacts.tsv` that matches the kind, with its `fold` and `p_enriched`: the sequence calls for `centromere` and `sequence_centromere`, the telomere arrays for `telomere`, the pooled blocks of all families for `rDNA`, `subtelomeric_block` and `interstitial_block`; `.` for the other kinds, where that class was not tested, and without Hi-C |

## `asm.feature_contacts.tsv` (with Hi-C: one row per feature class)

Does a class of features contact itself across chromosomes more than as many places drawn at random on the same
chromosomes? The sequence calls are the positive control: centromeres cluster in Rabl nuclei. Each instance stands as
one bin, the bin at its midpoint or the nearest usable bin inside it. The deciding library tests every class with at
least 3 instances on 2 or more chromosomes; in a genome of more than 15,000 bins nothing is tested (a warning says so).

| column | meaning |
|---|---|
| `feature_class` | the class: `centromeres (sequence calls)`; `chromosome ends` (the outermost usable bin at each end); `telomere arrays` (at the ends); `45S rDNA arrays (NORs)`; the satellite clusters of all families by role (`interstitial blocks (all families)`, `subtelomeric blocks (all families)`, `rDNA arrays (all families)`: role `rDNA` or an rDNA tag); each family's clusters of ≥ 20 kb (`satellite <family>`); each concentrated LTR lineage's cores (`LTR <lineage> cores`) |
| `kind`, `label` | the class's kind (`centromere`, `ends`, `telomere`, `nor`, `interstitial`, `subtelomeric`, `rDNA`, `satellite:<role>`, `te_core`) and label (`sequence calls`, `outermost bins`, the telomere unit, `18S motif`, `pooled`, the family or the lineage) |
| `instances`, `chromosomes` | how many instances have a usable bin, and on how many chromosomes the tested ones lie |
| `trans_pairs` | the pairs of tested instances on different chromosomes; a class of more than 400 instances is tested on a fixed random 400 |
| `observed_oe` | the mean trans observed/expected contact of those pairs: the balanced contacts of their two bins over the mean of their chromosome pair, smoothed over 3 × 3 bins |
| `null_median_oe` | the median of that mean over 1,000 draws that put each instance at a random usable bin of its own chromosome |
| `fold` | `observed_oe` / `null_median_oe` |
| `p_enriched` | one-sided: (1 + the draws at least as high as observed) / (1 + the draws), so at least 1/1,001 |

## `asm.line_tests.tsv` (one row per pair of lines)

Each line places the centromere on its own, so two lines that agree on most chromosomes beyond chance validate each
other without any truth, and two that do not show that at least one of them is wrong in this genome. The rules test
some pairs for themselves (the Hi-C test, the arm fold); this table reports every pair the same way and changes no
call. A pair is listed when both lines place a locus on at least one chromosome.

| column | meaning |
|---|---|
| `line_a`, `line_b` | the two lines, in alphabetical order: `satellite` (the sequence call on a `satellite` or `weak_satellite` array, or the array a decisive LTR core replaced, `seq_overruled_start`-`seq_overruled_end`), `ltr` (the LTR core, `te_start`-`te_end`), `hic_anchored` and `hic_prior_free` (the anchored and prior-free Hi-C peaks), `fold` and `alignment` (the cores of the arm fold and of the trans arm alignment), `at_trough` (the AT trough), `cis` (the cis extreme, `cis_peak_start`-`cis_peak_end`). `cis` with `ltr` is not listed: the LTR cores chose the cis signal's sign, and its own test is against them |
| `chromosomes` | the chromosomes where both lines place a locus |
| `agree` | of those, where the two loci overlap, or lie within one Hi-C bin when either line reads Hi-C (`hic_anchored`, `hic_prior_free`, `fold`, `alignment`, `cis`) |
| `expected` | the agreements chance would give: over those chromosomes, the sum of the share of placements of the two widths that would agree, (the two widths + twice the distance `agree` allows) / the chromosome length, at most 1 each |
| `p` | the Poisson-binomial P of at least `agree` agreements, with those shares |
| `verdict` | `agree on most, beyond chance` (P ≤ 0.01, on more than half), `agree beyond chance, on a minority` (P ≤ 0.01, on at least 3), `no agreement beyond chance`, or `too few chromosomes` (under 3) |

## `asm.te_signature.tsv` (with `tes/` from LTRquest: one row per chromosome with a core)

The LTR signature ([design.md](design.md#ltr-retrotransposons-tes-satsignaturepy)): the lineage or family that seeds
it and, on each chromosome, the families recruited on the other chromosomes and their densest region there (the core).
Without a signature (most elements lack `family=`, or no lineage or family qualifies as a seed) there is no table. The
seed, `decisive`, `why` and the ages are genome-wide: the same on every row.

| column | meaning |
|---|---|
| `chrom` | the chromosome |
| `marker_lineage` | the seed's lineage: the LTR marker |
| `seed_level`, `seed` | the seed: a `lineage` or a `family`, and its name |
| `seed_specificity` | the seed's fold enrichment in its cores: the share of its elements inside them over the share of those chromosomes they cover |
| `decisive` | `True` when the signature decides: a CRM seed with significant cores on most chromosomes; another lineage needs them too, at 10-fold or more, and no dominant satellite that passes the one-per-chromosome test (or that too few chromosomes leave untested) |
| `core_start`, `core_end` | the core: the densest region of the chromosome's signature elements (full width at half maximum) |
| `core_elements`, `chrom_elements` | the signature's elements in the core, and on the chromosome |
| `core_p`, `significant` | the binomial P of at least that many of the chromosome's signature elements in a core of that width, times the chromosome length over the core width (Bonferroni over its possible places), at most 1; significant at P < 0.001 |
| `source` | `signature` (the core of the families recruited on the other chromosomes) or `seed` (the seed's own core, where those families give no significant core) |
| `n_families`, `families` | the families recruited for this chromosome: how many, and the first 20; with `source` `seed`, the seed |
| `why` | the verdict on the signature, in words |
| `age_in_cores`, `age_outside` | with LTRquest's `K2P_T` (insertion time from the divergence of an element's two LTRs, years): the median age of the signature's elements inside the significant cores and outside them, on the chromosomes with one |
| `young_old` | the share of the younger half of those elements that lies inside the cores, over that of the older half: above 1, young copies mark the cores |
| `core_fold` | the core's own fold enrichment: the share of the chromosome's signature elements inside it over the share of the chromosome it covers |
| `core_age_from` | when the signature's elements form several comparable peaks on this chromosome and the youngest is significantly younger than the highest, the core is that younger peak and this column names the density maximum it replaced (1-based start-end); else `.` |
| `core_alt` | a second density peak of the signature's elements on this chromosome, at least half as high as the core's, with at least 3 independent insertions (1-based start-end; a second place to look; reported only), else `.` |

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
| `--satellite-dir DIR` | take the satellite outputs of an earlier run (`satellite/`) as they are, instead of running the satellite engine; cenfuse warns if their chromosomes differ, or if that run's `run_info.json` names another cenfuse version (its satellite rules are kept) |
| `--seed-engine fastan\|trf` | the tandem-repeat finder that seeds the satellite families [fastan] |
| `--seed-cache DIR` | a folder of tandem-repeat results shared between runs |
| `--classify-from T=S[+S2],...` | genome T (e.g. a contig assembly) takes the centromeric families of genome(s) S |
| `--force STAGE` | recompute STAGE: `all` (every stage), `satellite` (the whole satellite engine), `telomere`, `rdna`, `rdna5s`, `gaps`, `attrough`, `genotype`, `hic`, `fold`, `karyotype`, `meth`, `cpg`, or a satellite-engine stage (`seeds`, `families`, `annotate`, `classify`, `centromeres`, `hor`, `output`, `plot`: the engine from there on). A later stage reruns when what it reads changed, as on any rerun |
| `-v` | per-step progress and sanity checks |

`--force` is needed only when a run stops and names a stage (see Running, below).

Each engine also runs alone (`python -m cenfuse.sat run`, `python -m cenfuse.hic call`, `python -m cenfuse.meth run`;
`-h` for its options); the benchmark ran them through `cenfuse run`.

## Negative controls (`cenfuse control`)

`cenfuse control KIND IN OUT [--seed N] [--fai FAI] [-v]` writes a copy of one input with the signal one line reads destroyed
and the rest kept. Run cenfuse on a copy of GENOME_DIR with the control in place of the input: no test should pass and no call
should change against the run without it. The output is deterministic given `--seed` [1]; for a given seed the Hi-C controls
write the same bytes as the benchmark's (NEG, NEG2).

| KIND | randomized | kept |
|---|---|---|
| `hic-repair` (NEG2) | the partners of all read ends | each end's chromosome, position, strand and MAPQ, so every bin's coverage (about 40 bytes per pair in memory) |
| `hic-relocate` (NEG) | each end's position on its chromosome | each contact's chromosome pair |
| `gc-shuffle` | the order of 100-kb blocks across the chromosomes (laid end to end, shuffled, cut again at their lengths) | each chromosome's length, the set's base composition, other sequences; writes OUT's `.fai` |
| `ltr-shuffle --fai G.fai` | each LTR retrotransposon's position, with its LTRs and domains | lengths, clade, family, age |

The Hi-C controls check their own output: if contacts come out clumped (more than chance gives), the command stops and
leaves `OUT.part`.

## Running

**Resuming.** Rerun the same command in the same output folder and cenfuse resumes from its checkpoints. Each stage
is keyed on the code it runs, cenfuse's version and the content of the inputs it reads (not file names or dates, so a
moved or copied folder resumes), so an upgrade or a changed input recomputes only what it affects. Adding `tes/`,
`genes/`, Hi-C or methylation to a finished run needs no option. The run stops and asks for `--force` in two cases: a
stage of the satellite search (hours of work) is out of date after an upgrade, or the chromosome set changed (a
genome added or removed, a new or edited `chromosomes.txt`): rerun with `--force satellite`.

**One run per folder.** A run locks its output folder (`.cenfuse.lock`); a second run on the same folder stops. A lock
left by a run that died is taken over with a warning (on another host, after 15 minutes).

**The log.** Without `-v`: the start, milestones, warnings and the end. The satellite search reports its progress
with an estimate of the time left at every 10%.

**Reproducibility.** Every run writes `run_info.json`: the versions of cenfuse, Python, the packages and every
external tool, the fingerprints of the inputs, and the run's cost. On Linux, fresh runs give identical results apart
from dates and commands in file headers. On macOS (no `setarch`), about 1% of the satellite search's array cuts can
differ between runs; cenfuse warns.

**Memory.** The satellite search holds about one search workspace (up to ~2.4 GB for 2-kb units) per thread: lower
`-t` if memory is short.

**FasTAN crashes.** FasTAN can abort on a large satellite array. cenfuse splits the piece and retries, skips a piece
that still crashes (with a warning; listed in `<seed cache file>.fastan_crashes.tsv`), and stops if every piece of a
chunk crashes: rerun with `--seed-engine trf`.
