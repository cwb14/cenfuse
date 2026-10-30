# How cenfuse works

## Principles

1. **Each line calls on its own first**: satellites (+ TEs), Hi-C, methylation.
2. **Lines check each other.** No line is trusted on its own confidence: the Hi-C engine's prior-free `high` calls
   were right on only 45% of T2T mammal chromosomes.
3. **Whether Hi-C works in a genome is measured, not assumed.** Hi-C must find the *other* chromosomes' centromeres
   (the gate) before it may act on one. CpG dips must mark them before they place a kinetochore.
4. **The final call is the sequence call** unless a validated line moves it under the rules below. Every move and
   every disagreement is reported.
5. **Fixed constants, each with a reason** (a convention, a majority rule, a data-derived value, or a named
   heuristic). Nothing is tuned to a benchmark, and there are no options for thresholds.

Not goals: centromeres from methylation alone (false dips are as deep as true ones); Hi-C choosing the satellite
family.

**LTR annotations: LTRquest only** (decided 2026-10-28). The LTR line is developed and benchmarked with LTRquest
annotations only, and work on it assumes LTRquest input. The DANTE_LTR and EDTA + TEsorter readers stay as they are,
with their tests, but get no further work: no new rules, workarounds or benchmark arms for them. Their results in
[benchmark.md](benchmark.md) (LTR annotation) are the final record of how they compare.

## Stages

```
satellite -> telomere -> rdna -> rdna5s -> gaps -> genotype -> hic -> fold -> karyotype -> provisional -> meth ->
fuse -> cpg -> plants -> morphology -> interpret -> loci -> output -> report -> plot
```

Each stage is checkpointed, keyed on its inputs and code, so a rerun resumes. Every stage's key includes `params.py`: a
constant used by one new module lives in that module, so adding it does not recompute every stage. `genotype`,
`karyotype` and `plants` handle several plants ([plants.md](plants.md)).

## Lines

### Satellites and TEs (`cenfuse.sat`)

The pipeline: tandem-repeat seeds (FasTAN or TRF) → families (clustered, consensus) → genome-wide annotation
(nhmmer) → a role per family → one call per chromosome → HOR structure.

Roles: `centromeric`, `cen_associated`, `subtelomeric`, `interstitial` (knobs and other blocks), `dispersed`,
`minor`, `rDNA` (rRNA genes found with Rfam models).

| rule | test | why |
|---|---|---|
| family founding | a seed cluster becomes a family only if one of its seeds is an array by annotate's rule (>= 3 copies); then by size (>= 20 kb of seeds, at most 40) | a satellite is tandem. Clusters founded on two-copy duplications only (0.8-2-kb units, retrotransposon-like) were dispersed or absent in every genome of the 216-genome replay (12 families in 7 genomes) and took most of the genome-wide search where they occur (T2T sheep 93%, siamang 52%, orangutans 19-28%, barley 14%); without them calls and roles were unchanged (siamang, barley, B. orangutan); search CPU at equal machine load: sheep 11.7x, siamang 2.0x, B. orangutan 1.3x, barley 1.16x less |
| dispersed repeat | tandem_frac < 0.5 **and** whole-unit hits < 0.5 **and** locus fold < 2 → never centromeric, and role `dispersed` (never a satellite block) whatever its concentration | interspersed elements (LINE/SINE fragments) annotate as chains of partial hits. Satellites are tandem; retrotransposon-derived centromeric arrays are concentrated (fold ≥ 2.5). Before v7 the role ignored the test: one 0.5-1.9-kb family per genome (chimpanzee, bonobo, S. orangutan, rhesus, sheep; 215-719 Mb, 2-14% in tandem) was an `interstitial` block chained along whole chromosomes, making chromosome-scale loci |
| satellite kind | dominance judged per kind: subfamilies related by sequence are pooled | splitting alpha into 6 subfamilies let a pericentric satellite out-count it |
| acrocentric | a strong family at one end of most chromosomes, and at both ends of too few for independent ends (binomial, P ≤ 0.01), with no strong interior rival (judged per satellite kind) → centromeric, flag `acrocentric`, at most `medium` | subtelomeric repeats occupy ends independently; a centromere occurs once per chromosome. Needs ~16+ chromosomes. A subfamily at the inner edge of a terminal block is that block, not a rival (T2T mouse: a major-satellite subfamily beyond the 5% end zone had blocked the test) |
| the call | the cluster holding classify's chosen family (layered families pooled); a decisive CRM core may favor the largest cluster | |
| competing cluster | a second cluster ≥ half the main one, or the terminal block of a second acrocentric satellite unrelated to the call → two levels down (`low`), reported in `seq_alt_*` | such calls were right 2/10 times on CENP-A/CENH3 truth. The acrocentric test finds a centromere's end, not which of two satellites there holds the kinetochore (mouse: the 234-bp major satellite out-sizes the 120-bp minor satellite that carries CENP-A; related layers, as in goat, are pooled instead) |
| CRM cores (`tes/`) | the CRM clade's densest cores (full width at half maximum of its density; a core needs >= 3 independent insertion sites, overlapping elements counting once); a decisive core supports a call it agrees with, takes a chromosome without a call, replaces a `low` call, and replaces a `high` or `medium` call only if the clade marks it significantly more than the call (its elements in the core vs in the call, binomial at 1/2, P <= 0.01), and so does a confident call beside the core (agreeing without overlap: `satellite_beside_core`); otherwise the call stays, `te_discordant`, one level down | CRM retrotransposons mark grass centromeres. Site floor: *Chlorella* NS4-2's Osser (Copia) cores, 2-8 elements per chromosome, 53/54 inside CENH3: 0 -> 8/10 right (all `low`); no call changed in 12 crop `tes/` arms, B. rapa, soybean, sorghum. Head-to-head: in a broad CRM plateau the densest window is noise (B. rapa CC1/CC2 A09, soybean X26 chr1/chr6: CENH3 on the satellite, 16-26 vs 14-32 elements); a satellite on a smaller second cluster loses (P. annua Pa1A 19 vs 62; Ae. geniculata 7M/2U) |
| CRM cores short of a majority | where CRM is not the decisive marker, each significant CRM core replaces the marker's core on its chromosome (flags only; without a qualifying clade, the CRM cores alone are the marker) | CRM is the prior centromere-specific lineage. With less sensitive annotations CRM concentrates on only some chromosomes: rice DANTE_LTR 5/12, EDTA + TEsorter 4/12, LTRquest thinned to 50% 3/12; every such core lay on CENH3 (24/24 over 8 arms of rice, potato, maize and *Brachypodium* not used to make the rule), every rice Tekay core off it |
| honesty cap | centromeric families but none dominant, and no decisive (CRM) core: every call `low` (`no_dominant_satellite`); a flags-only marker (Athila, Tekay, Bryco) lifts the cap only for the calls it supports (`te_concordant`) | satellite evidence alone cannot tell a centromeric array from an interstitial one, and a line that contradicts a call, or says nothing about it, is no evidence for it. Before v6 any marker lifted the cap genome-wide: moss with `tes/` had 4 confident wrong calls (0 from FASTA alone), the 185-genome batch 66 confident calls the TE line did not support |
| `tes/` formats | an element's clade from `clade=` (type holding `LTR_retrotransposon`: LTRquest, EDTA + TEsorter) or from a DANTE_LTR `transposable_element`'s `Final_Classification` (its last REXdb lineage field; an element classified only to `Ty3/gypsy`, `chromovirus` etc. has none; rank D, domains without LTRs, is not an element). `.gff3` or `.gff`. A file with no classified element on the chromosomes is reported with the reason (no lineage, or other sequence names) | DANTE_LTR is a common plant classifier; Arabidopsis with its file: the same 5 calls as with LTRquest. Its rank-D clusters made 195 of potato's 223 CRM and put decisive cores up to 15 Mb off CENH3 (6 confidently wrong calls; 2 without them). EDTA alone gives superfamilies only, which no rule can use |
| genes (`genes/`) | gene-rich families are vetoed; subtelomeric = beyond the outermost genes. GFF3 or GTF; protein-coding `gene` lines (`gene_biotype`, `biotype`, GENCODE `gene_type`), or one locus per transcript `Parent` in a file without `gene` lines | potato DM v6.1 (Spud DB: mRNA lines only) 2 -> 4/11 on CENH3; barley: the 3 confident wrong satellite calls become `low` |

### rDNA and telomeres

- **NORs** (`rdna.py`): tandem runs of ≥ 3 copies of the universal 18S motif (GTACACACCGCCCGTC), ≤ 100 kb apart.
  This finds CHM13's 5 NORs in 32 s. Hi-C peaks beside a NOR never oppose a call (NORs cluster in the nucleolus).
- **5S rDNA** (`rdna.py`): tandem runs of >= 3 whole 5S genes, <= 20 kb apart. A gene is a hit of the Rfam 5S profile
  HMM (RF00001, the filter model inside the bundled covariance model) covering >= 90% of it, E <= 1e-5 per chromosome.
  Whole genes, because pseudogenes and 5S-derived fragments are mostly partial (CHM13: 126/126 hits of its 1q42 array
  whole, 108/434 scattered hits), and dispersed genes form no array (*S. pombe*: 36 genes, none). The satellite engine
  sees a 5S array only as one of its <= 40 families with units <= 2 kb, so it missed the 2.2-kb units of every T2T
  mammal. A loci feature only: 5S arrays do not gather in the nucleolus, so Hi-C does not use them. CHM13: 24 s.
- **Telomeres** (`telomere.py`): the unit is found de novo; arrays are reported at chromosome ends: >= 300 bp of unit
  copies anywhere in the terminal 100 kb, or >= 100 bp when the array reaches the sequence end (within 1 kb). Short
  telomeres are real (Neurospora NC12: 157-207 bp at 13/14 ends) and assemblers trim them; once the unit is the
  genome's (10-fold denser in tandem at the ends), 100 bp of it at an end is no chance match.
- **Interstitial telomeric sequences** (ITS, `telomere.py`): arrays of the same unit inside the chromosomes (beyond the
  terminal 100 kb; >= 300 bp of unit copies), as `interstitial_telomeric` loci. They mark old chromosome fusions (CHM13:
  16 arrays, most within 1 Mb of an end, one at the 2q13 fusion site) and the cores of some centromeres (chicken:
  78-210-kb TTAGGG blocks at the CENP-A domain of several chromosomes).

### Hi-C (`cenfuse.hic` engine + `hicfuse.py`)

1. **Prior-free call**: the engine's own iteration. It decides nothing by itself; only the fold test can validate it.
2. **Anchored score.** For chromosome x: its mean contact with partner windows at the *other* chromosomes' sequence
   calls (calls ≥ `medium`; all calls if fewer than 3). x's own call never shapes its score.
   - A partner window is the call midpoint ± 4 bins when most of those bins are mappable.
   - Otherwise it is the 9 mappable bins nearest the call: the array's flanks, since the array itself is unmappable.
3. **Peak tests** (the engine's): 100 bootstraps, leave one partner out, 200 null replicates.
   - A peak is *strong* when it passes all three and is not flagged `hub`, `rDNA` (incl. beside a NOR) or
     `telomeric` (end window).
   - A peak with any of those flags may support a call but never moves or opposes one.
4. **Distances in mappable bins.**
   - Agreement: within 1 bin, forgiving ≤ 4 unmappable bins (a peak at an array's edge supports a call inside it).
   - Opposition: beyond 4 *mappable* bins. Across a longer blind stretch Hi-C neither supports nor opposes.
5. **Gate** (per chromosome x).
   - References: the other chromosomes with `high` calls (`high` + `medium` if fewer than 3).
   - x passes when most references' own peaks agree with their calls, beyond chance (Poisson-binomial P ≤ 0.01),
     and also a majority of the whole reference set.
   - Fewer than 3 references: `untested`.
   - No call high or medium anywhere (no satellite family dominates, so every call is capped `low`: siamang without
     genes, chicken): the low calls are the references, and an agreement counts only at a clean peak (no hub, rDNA or
     telomeric flag); in that genome a peak supports a call only when clean. Chromosome ends, hubs and NORs cluster in
     3D whatever the centromere: in a Rabl nucleus, wrong calls at the ends would agree with Hi-C. Two independent
     lines that agree beyond chance on most chromosomes validate each other. `hic_gate_refs` says which references
     were used.

### Arm fold (`fold.py`)

In a Rabl nucleus both arms run from the centromere pole to the telomere pole. The chromosome's own (cis) contact
map then has an anti-diagonal ridge through the centromere.

- **Fold score of a center**: mean observed/expected over the pairs (c - d, c + d) spanning ≥ 9 bins, both
  mappable. A center inside an unmappable array is read from its flanks.
- **Scan only centers with centromeric index ≥ 12.5%.** Below that a chromosome is telocentric (Levan 1964) and the
  short arm cannot fold. Without the bound, centers near the ends won on noise.
- **Fold point** = top center; **core** = its half-maximum run. It is coarse: 5-10% of the chromosome.
- **Fold test**: the prior-free Hi-C peak lies inside the fold core on a majority of chromosomes, beyond chance
  (P ≤ 0.01). Trans clustering and cis folding are independent readouts, so passing validates both without the
  sequence line.
- Passes in barley and P. annua; fails in all mammals, potato and rice (rice's prior-free peaks are wrong).

### Methylation (`cenfuse.meth` engine + `methfuse.py`)

- **Anchors**: provisional calls ≥ `medium`.
- **Signature**: calibrated with x's anchor left out.
- **Candidate interval**: `centromere_like` when the signature score at its most extreme CHG window is ≥ 0.5.
- **Role**: support and second line only; it never moves a call alone. Only samples from the assembled plant (or of
  unknown origin) decide.
- Needs `tes/` and a CHG file for every sample.

### CpG kinetochores (`cdr.py`)

In primates CENP-A sits in a CpG-hypomethylated dip of the methylated active array (the CDR). The satellite and Hi-C
lines locate the array (0.8-5.6 Mb); the dip locates the kinetochore inside it. Inputs: the CG file of each deciding
methylation sample (modkit bedMethyl; ENCODE bedMethyl, whose coverage and percent are modkit's own columns 10-11,
read as counts; or a bedGraph of methylated fraction); `tes/` and CHG are not needed. The assembly check reads ENCODE
rows too. A CHG file must be modkit bedMethyl (the provenance test reads modkit's counts of failed and differing
calls); in any other layout it is not used (warning) and the sample's CG file is kept.

- **Profile**: CpG methylation per 5-kb bin (CDR-Finder's window) on bins that centromeric-role arrays cover for at
  least half, inside the final call: the bin's pooled calls (methylated over valid, every CpG with a call) from a
  bedMethyl, the mean fraction from a bedGraph. DNA between arrays (CpG islands, genes) is left out. Pooling is the
  binomial estimate at any coverage: with a 3-call floor per CpG, MMU8 thinned to 1/16 of its reads (~1 call per CpG
  strand) lost its gate; pooled, it keeps 21/21 kinetochores within 50 kb of full depth.
- **Dip**: the most prominent valley (higher methylation on both sides, so an array edge is never a dip), extended to
  its half-depth run (bins below the midpoint between its bottom and the array median).
- **Dip test**: the run is longer than in the call's own bins shuffled (200 shuffles; significant at P ≤ 0.05).
- **Positional test**: the run is longer than when the chromosome's CpG values are rotated by 200 random offsets onto
  the call's own CpGs and array bins (methylation from elsewhere; same size, bins and CpG density).
- **Gate** (per chromosome, leave one out): references = the other chromosomes' `high` calls. k of n have a positional
  dip; binomial tail at 0.05 (each reference's chance under the null) ≤ 0.01. No majority is required: the
  positional test is conservative (passing primates: 28-95% of tested calls; rotated tracks: 4-9%).
- **Role**: where the gate passes, every call with a significant dip test gets its kinetochore (`cpg_kt_start/end`,
  GFF3 `kinetochore`). It never moves a call or changes its confidence, and it cannot tell a wrong array.
- **Assembly check**: bedMethyl rows must sit on reference CpGs (as for CHG); bedGraph rows too, at the row start or
  one base off (public tracks differ by one base).

## Fusion rules (`fuse.py`)

| situation | result |
|---|---|
| no validated line | sequence call, `sequence_only` |
| validated Hi-C (or methylation) agrees | sequence call, `concordant`, +1 level per line (max +2) |
| validated strong Hi-C peak > 4 mappable bins away, call `low`/`medium` | moves to the Hi-C core, `offset`, `medium` |
| same, call `high` | moves only with a second line at the peak and not at the call: a CRM core within 4 bins, or methylation `centromere_like` at the peak and a measured dip at the call (`moved_with_te`/`_meth`, `high`). Otherwise `contested`, -1 |
| Hi-C peak 1-4 bins away | neither supports nor opposes (`hic_near`) |
| competing cluster, strong peak agrees with the second one, mappable bins between them | the second cluster is chosen (`hic_chose_alt`) |
| no sequence call, strong peak | Hi-C core, `functional_only`, `medium` |
| no sequence call, clean medium peak | Hi-C core, `functional_only`, `low` (`hic_medium_peak`) |
| gate not passed, fold test passed; fold core and prior-free peak agree, the sequence call at neither | moves to the prior-free core, `offset`, `medium` (`fold_pair`; `moved_with_fold` for a high call) |
| same, but the fold or the prior-free peak opposes the kept call | `contested`, -1 per opposing line (`fold_offset`, `hic_denovo_offset`) |

- **Confidence**: the satellite engine's level, then ±1 per line as above. A contest keeps the level when
  methylation supports the call. A moved or Hi-C-placed call is `medium`, `high` with a second line.
- **Support / against columns**: `satellite` (centromeric arrays inside the call, ≥ 5 kb), `te`, `hic`, `meth`,
  `fold`.
- **CpG kinetochores** are added after fusion (`cdr.py`): they change no call, confidence or flag.
- **`hic_sequence_conflict`** (per genome): Hi-C's prior-free peaks are confident but disagree with the sequence
  calls and the gate fails. One line is wrong genome-wide. The fold settles it in Rabl nuclei; otherwise add `tes/`
  or methylation.

## Reading the integration (`morphology.py`, `interpret.py`)

Nothing here changes a call; every statement is read from the tables cenfuse writes.

- **Morphology.** The integrated call fixes where the arms meet. Arms are measured from the CpG kinetochore's
  midpoint where the CpG line placed one, else from the call's midpoint: centromeric index (short arm / length), arm
  ratio (long / short), Levan class (Levan, Fredga & Sandberg 1964: metacentric <= 1.7 < submetacentric <= 3 <
  subtelocentric <= 7 < acrocentric; telocentric falls in the last, a sequence cannot show a missing short arm). A
  call spanning two classes gets both; `short_arm` is `.` when the call spans the chromosome's middle (CHM13 chr1's p
  arm is 50 kb the longer by sequence). The summary's karyotype counts `high` and `medium` calls only: a `low` call's
  class is a guess. Sequence classes are not cytogenetic descriptors: by sequence, CHM13's acrocentric chr21 and chr22
  are submetacentric/subtelocentric (their T2T p arms, rDNA included, are 11-15 Mb).
- **Completeness.** Telomere arrays at the ends (telomere.py), assembly gaps (runs of >= 10 N, the satellite
  engine's rule) on the chromosome and inside the call; `t2t` = both telomeres and no gap.
- **basis / caveat.** One sentence each per call: the sequence evidence and its size in the call, a move and its
  second line, every agreeing line, the kinetochore; then every flag that lowered the call, every opposing line with
  where it points, gaps in the call, and the input that would decide where cenfuse's own tests say one can (two arrays
  of one satellite: CENP-A/CENH3 data, since Hi-C needs mappable DNA between them and CpG methylation was tested and
  does not choose; satellite-free or capped calls: Hi-C or a CRM core when not given).
- **Report** (`cenfuse_report.md`): per genome, the lines and their tests, the calls by confidence and what they stand
  on, the calls that are not `high` with their caveats, morphology and completeness, the other features with their
  satellite share and 3D association. A failing report section is a warning; the results stay complete.

## Loci (`loci.py`)

Every candidate interval becomes a locus member:
- satellite clusters ≥ 20 kb of feature roles;
- TE cores;
- Hi-C cores and hubs;
- methylation candidates;
- telomere arrays at the ends and inside the chromosomes;
- the sequence and final calls.

Overlapping members merge; only Hi-C members are padded (1 bin). Each locus gets every line's value there and the
first kind that matches:

`centromere` > `sequence_centromere` (where the call moved from) > `telomere` > `interstitial_telomeric` > `rDNA` >
`subtelomeric_block` >
`interstitial_block` (knobs) > `alt_centromere_candidate` > `hic_hub` > `hic_peak` > `te_core` > `satellite_block` >
`meth_candidate` > `other`.

## Constants (`cenfuse/params.py`; engines keep their own)

| constant | value | reason |
|---|---|---|
| HIC_TOL_BINS | 1 | Hi-C placement error |
| HIC_FAR_BINS | 4 | Hi-C opposes only beyond its partner half-window (its resolution) |
| GATE_P | 0.01 | conventional evidence level |
| GATE_MIN_REF | 3 | fewer references give no majority test |
| ANCHOR_LEVELS | high, medium | wrong anchors blur Hi-C and methylation |
| NOT_STRONG_FLAGS | hub, rDNA, telomeric | hubs and NORs contact everything; chromosome ends are the main non-centromeric trans cluster |
| MAX_RAISE | 2 | low → high at most |
| METH_SCORE | 0.5 | the methylation engine's signature rule |
| SAT_SUPPORT_BP | 5 kb | the satellite engine's weak-satellite floor |
| TANDEM_FRAC, ACRO_ALPHA | 0.5, 0.01 | majority rule; conventional level |
| LOCUS_FOLD | 2 | heuristic (as the Hi-C engine's hub fold) |
| fold: S_MIN, CI_MIN | 9 bins, 12.5% | beyond the ± 4-bin centromeric domain; Levan's telocentric bound |
| CpG (`cdr.py`): BIN, MIN_SITES, MIN_BINS | 5 kb, 10 CpGs, 10 bins | CDR-Finder's window; heuristics (a dip needs shoulders) |
| CpG: N_PERM, N_ROT, DIP_P | 200, 200, 0.05 | P resolution 10× finer than the conventional level |
| 5S (`rdna.py`): FIVE_S_E, FIVE_S_UNIT_MAX, MIN_UNITS | 1e-5, 20 kb, 3 | the satellite engine's nhmmer E; twice the longest 5S-bearing unit (9.1-kb yeast rDNA); a tandem array |
| morphology (`morphology.py`): LEVAN, GAP_MIN_N | 1.7 / 3 / 7, 10 N | Levan et al. 1964 (convention); the satellite engine's gap rule (assemblers write gaps as >= 10 N) |

## Decisions and alternatives tried

| question | chosen | tried and rejected (why) |
|---|---|---|
| fix Hi-C from inside the Hi-C data | anchor to the sequence calls | random restarts and caps: held-out accuracy unchanged |
| Hi-C partner windows over unmappable arrays | engine window if mostly mappable, else the 9 nearest mappable bins | padding every window: lost agreements in CHM13, chimpanzee, S. orangutan. Replacing only empty windows: right peaks 97 → 86 of 139 |
| blind bins in Hi-C distances | ≤ 4 forgiven for agreement, none for opposition | full credit: let Hi-C support a wrong bonobo call |
| gate chance | uniform share + majority of the whole reference set | chance from the null replicates: removed 1 of 2 spurious passes on randomized sheep Hi-C (the majority removed both) |
| NOR detection | 18S motif scan | Rfam SSU HMM search: stalled for minutes on low-complexity sequence |
| 5S detection | Rfam 5S HMM scan, whole genes in tandem | satellite families tagged by Rfam (kept): miss units > 2 kb (mammals, birds) and arrays outside the 40 families |
| Hi-C from a Juicer `.hic` | stop with the `hictk convert` recipe | reading `.hic` natively: a new dependency for what one command does (rice: converted `.mcool` gives the same calls as the pairs) |
| LTR elements without a lineage (fungi) | not used | pooling them when no lineage concentrates (Cryptococcus JEC21 1 -> 10/14 on CENP-A with 20-kb TE bins): in fish and *Oryza minuta* B the pool lifted the no-dominant cap (`low` -> `high`). Again with the per-call cap, unclassified elements inside satellite arrays left out and a pool that never lifts the cap: JEC21 10/14 and the hybrid's D subgenome 0 -> 10/13 (all `low`, none wrong), but the pool still qualified in 9 fish genomes (18-46 flags each, no truth there) and the 20-kb bins lost a CRM call in *O. punctata* |
| dispersed-repeat test | tandem + whole-unit + locus fold | without the locus fold: 160 batch calls moved off TE-derived centromeric arrays (Coffea, Glycine, Rosa) |
| CRM core vs classify's choice | the core may favor only the largest cluster | any cluster the core touches: lost the CENH3 call of Ae. geniculata chr02U |
| decisive core vs a distant confident satellite call | the core replaces it only if the clade marks the core significantly more (binomial head-to-head) | the core always (v7): A09 wrong with `tes/` in B. rapa CC1/CC2 (9/10 vs 10/10 without), soybean X26 chr1/chr6 (17/20 vs 20/20). Bootstrap of the core position (keep the call unless <= 1% of resampled cores agree): fixed A09 but lost P. annua Pa1A, whose satellite sits on a real second CRM cluster (15% agreement vs A09's 6%). Element counts without a margin: A09 16 vs 18 |
| confident satellite call beside a decisive core | the core replaces it if the clade marks it significantly more (the same binomial head-to-head) | kept as agreeing (v7): Ae. tauschii Chr2/Chr4 arrays 5-7 Mb beside a 292-332-element CRM core holding CENH3 (3/7 with `tes/`; with the rule 5/7); replay: 17 calls in 15 genomes move to the core, with truth Ae. geniculata chr02M wrong -> right, S. verrucosum chr05 wrong -> wrong (closer); B. rapa A06 (call holds more CRB than the core) kept |
| TE core floor | 3 independent insertion sites | 5 elements (v7): *Chlorella* Osser cores on 5/12 chromosomes, no marker. Sites in the binomial test too: cores moved in 113/220 replay genomes |
| skip dispersed families on a sample of the genome (sample screen) | not adopted | every model searched on each chromosome's end chunks + every 5th chunk; a family skipped elsewhere when the sample shows at P <= 0.01 that classify would call it dispersed, and dropping it changes no other family's sampled arrays. Simulated on 216 stored annotations: 0 false skips. On real runs (sheep, mouse, S. orangutan, moss, *Brachypodium*): 0 skips. Dispersed families overlap each other on the same retrotransposon copies, and dropping them changes a few sampled arrays of small families (mouse: of the minor and major satellites); the check is needed (without it moss Chr08's call moved); < 7 chromosomes never reach the sign test; families whose hits rarely chain into arrays cannot be measured. Without the check mouse would gain <= 1.6x |
| FasTAN reruns (reproducibility) | FasTAN under `setarch -R` (address-space randomization off) where the OS allows | as is: FasTAN 0.8's cuts depend on memory layout (one moss 5-Mb chunk, 5 runs: 17/982 arrays cut differently), so consensuses, families and some calls varied between fresh runs (moss: 5 call boundaries, one by 0.9 Mb); off: 5/5 runs identical, a whole fresh moss run byte-identical twice |
| fold center selection | centromeric index ≥ 12.5% | lower confidence bound of the mean: biased to the chromosome middle (barley 2/7) |
| Hi-C peaks in end bins | support only, never oppose | (rice chr4: an end-bin peak contested a correct call) |
| numpy huge pages | off by default | on: Hi-C bootstrap 30× slower on fragmented memory |
| CRM on a minority of chromosomes | its cores replace the marker's per chromosome | CRM as the genome's marker when its cores beat chance (binomial, P <= 0.01): with 5 chromosomes one core is enough, and in TAIR10 one CRM core displaced Athila (te_concordant on 5/5), two right `high` calls fell to `low` |
| telomeric Hi-C peaks | never strong, so never move a call (ends cluster in 3D whatever the centromere: all 10 such peaks in rice, barley, CHM13 and orangutan wrong; at reduced depth plant ones 0/27 right, other peaks 27/120) | strong when the LTR marker's core lies within 4 bins: moss (11/26 CENH3 domains within 2 Mb of an end) 22 -> 26/26, all confident, and no other run changes; moss alone prompted it and is the only genome where it acts (lead) |
| Hi-C pairs without MAPQ columns | a warning (uniqueness not checked) | an input check reading pairtools' parse command, or refusing such pairs: not needed. The same reads kept at MAPQ >= 1 changed no gate or fold and made no call confidently wrong in 5 genomes (CHM13, orangutan, rice, barley, moss); multi-copy repeat reads map at MAPQ 0, which pairtools drops at any floor |
| Hi-C fragility | reported as depth guidance | a split-half check (gate and fold rerun on half the contacts, flag when the result flips): stable genomes agree with another random half (13/13), but where two halves disagree with each other (rice, barley's fold, both orangutans, *P. annua*'s fold) one half matched another in 1/5, and Hi-C time doubles |
| CpG bin estimate | pooled calls per bin (bedMethyl) | mean of per-CpG fractions with >= 3 calls (the CHG floor): kinetochores within 50 kb of full depth at 1/8 and 1/16 of the reads: MMU8 15.3 and 0 of 21 (gate fails), MFA8 21 and 14.3; pooled: 21 everywhere |
| CpG dip inside a call | most prominent valley on array bins, half-depth run | deepest bin: lands on array edges (CHM13). Longest run below half the median: on the CDR 52/90 (valley 88/90) |
| CpG dip significance | shuffle test per call; positional (rotation) test for the gate | equal-length windows of other satellite DNA as null: infeasible where most satellite DNA is centromeric (MFA8: 3.7 Mb outside calls vs 3.5-11 Mb calls). Dips in calls vs dips in satellite clusters away from calls: passed on rotated (misplaced) methylation in gorilla, because long calls hold coherent runs of any structured profile |
| pyhmmer search pipelines (memory) | one per thread, reused (`annotate.search_one`) | a new pipeline per search (pyhmmer's `nhmmer()`): each sits in a reference cycle and its workspace waits for Python's full collection (chicken 277 GB). glibc `malloc_trim` after each search: 229 GB. Model-major job order: 125 GB. Smaller pyhmmer blocks: less memory, but hits change |
| MAFFT threads (reproducibility) | `--thread N --threadit 0` | `--thread N` alone: the multithreaded refinement gave a different alignment on every run (4 runs, 4 HMMs). `--thread 1`: reproducible but 2.5× slower MAFFT |
| TE marker among several qualifying clades | most concentrated (median core share) | chromosomes with a core first: picked Retand over Tekay in *Arabidopsis* (DANTE_LTR), but agreed less with the satellite calls over the 185-genome batch (te_discordant 111 → 129) |
| centromeric region (call + satellite clusters chained through gaps <= 500 kb) and a `pericentromeric_block` kind | not adopted | held out (HG002 CenSat, 5 ape centromeric regions): apes median precision 0.07 (an interspersed 1.5-kb family classed `interstitial` chains whole chromosomes), HG002 median recall 0.60 (0.5-0.7-Mb non-satellite stretches stop the chain); a longer gap would be fitted to this truth |
| chromosome-scale loci from an interspersed family classed as a satellite | the family's role follows classify's own dispersed-repeat test (`dispersed`) | loci members that are mostly satellite (sparse clusters split into dense runs): removed the chromosome-scale loci of 5 genomes but made thousands of small loci of the same family (chimpanzee +5,757); dropping sparse clusters would erase diffuse real satellites (barley subtelomeric, 108/156 clusters) |
| GFF3 feature types | Sequence Ontology terms; cenfuse's kind as `kind=`, members linked by `locus=` | cenfuse's own kinds as types (`interstitial_block`, `hic_hub`): no SO validator or annotation pipeline accepts them. `Parent=` for members: GFF3's Parent means SO part_of, and a satellite block is no part of a centromere (gt gff3validator rejects it) |
| CpG choosing among satellite clusters | not done | least-methylated window anywhere: on the CDR 3/23 (CHM13). Most prominent valley across clusters: wrong cluster on 8-10/25 (orangutans). A cluster with half the call's valley floor: 6 changes in 7 primates, 1 right (B. orangutan), 2 wrong (CHM13) |

## Known risks

- The gate takes the sequence calls as references. If those are wrong genome-wide, only the fold (Rabl nuclei),
  `tes/` or methylation can rescue Hi-C.
- One wrong Hi-C peak in a genome where Hi-C mostly works passes the gate. Only `high` calls need a second line to move.
- Hi-C tests lose power gradually but flip at their threshold: the gate failed in both orangutans at half depth and
  barley's fold below its 9.2 M pairs; the pig gate passed only at 10-25% depth. A failed test reverts to the sequence
  line (no mammal call became confidently wrong at any depth). One random half does not flag them reliably (tested); the
  share of several halves passing might.
- The decisive CRM test counts chromosomes with a significant core: annotation sensitivity matters (potato: 214 CRM
  elements give exactly a majority; half of them, DANTE_LTR or EDTA's intact set do not).
- Centromeres within 8 Hi-C bins of a chromosome end: Hi-C cannot move a call there (`telomeric` flag), and without
  `tes/` its gate cannot pass where every sequence call is `low` (moss FASTA + Hi-C: 0/26). Moss keeps 4 wrong `low`
  calls where Hi-C and the Bryco core both sit on CENH3.
- Satellite-free centromeres beside large pericentromeric satellites are beyond the sequence lines: in tomato
  Heinz 1706, 11/12 CENH3 domains hold no satellite array and no LTR clade marks them; the dominant 53-bp satellite
  (largest cluster on 7/12 chromosomes, so no honesty cap) took 7 confident wrong calls (1/12 right).
- Methylation leave-one-out signatures need ≥ 3 anchors with data.
- Hi-C resolution is a mappable bin (250 kb-2 Mb). Offsets inside an unmappable array are invisible.
- CpG kinetochores need an array long enough for a dip with shoulders (≥ 10 bins with data). A small array that is
  hypomethylated as a whole (bonobo and B. orangutan minicentromeres) is not refined, and CpG never chooses it.
- The CpG line was designed on the same primates it is benchmarked on (docs/benchmark.md says which tests are held out).
- Satellite-search memory: about one pyhmmer workspace per thread (model length × the 256-kb search block, up to ~2.4
  GB for 2-kb units): chicken 58 GB at 24 threads (277 GB before v5), *Arabidopsis* 4.8 GB at 8.
- Satellite outputs made before v7 (`--satellite-dir`) keep the old roles: there an interspersed family classed
  `interstitial` (one 0.5-1.9-kb family each in chimpanzee, bonobo, S. orangutan, rhesus and sheep) chains into
  chromosome-wide clusters and chromosome-scale loci. Rerun the satellite engine from classify (`--force classify`) to
  apply the v7 role rule. Calls are not affected; the report gives each block class's satellite share.
- Morphology is only as right as the call: a wrong call gives a wrong class (`low` calls are kept out of the summary's
  karyotype). An incomplete chromosome (a missing telomere, gaps) shifts the index.

## Tests

`pytest -q -n 8`: 538 tests on small synthetic fixtures.
- Unit tests of every rule above.
- A synthetic end-to-end genome with planted cases: a functional offset with a CRM core (must move), the same
  without a second line (`contested`), a hub decoy, a satellite-free chromosome, permuted Hi-C (no moves), three
  plants.
- One full run of the satellite engine on its synthetic genome.
- CpG: planted dips (found, positional), shifted tracks and dips everywhere (gate fails), every input format (modkit
  and ENCODE bedMethyl, bedGraph), malformed rows, the bedGraph assembly check, and a CG-only table that changes no
  call.
- Inputs: `.gff`/`.gtf` annotations, two gene files for one genome, and a `.hic` or BAM in `hic/` (stops with the
  conversion recipe).
- Reading: Levan classes at their bounds, arms from the kinetochore, gaps across scan chunks; one sentence test per
  relation and flag; the synthetic run's report (every chromosome, the moved, contested and placed calls).
