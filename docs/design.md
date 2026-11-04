# How cenfuse works

cenfuse reads several kinds of evidence (*lines*) that each point to the centromere: satellite arrays, LTR
retrotransposons, Hi-C contacts and DNA methylation. Each line calls on its own first. A line may then change the
final call only after it has passed a test in the genome at hand, and every agreement and disagreement is reported.
This page describes each line, the tests, how the lines are combined, and the alternatives that were tried and
rejected. Benchmark numbers quoted here are from [benchmark.md](benchmark.md).

## Principles

1. **Each line calls on its own first.**
2. **Lines check each other.** No line is trusted on its own confidence: the Hi-C calls the Hi-C engine rates `high` on
   its own lie within 4 bins of CENP-A on only 14 of 63 T2T mammal chromosomes.
3. **Whether a line works in a genome is measured, not assumed.** Hi-C must find the *other* chromosomes'
   centromeres before it may act on one; CpG dips must mark them before they place a kinetochore.
4. **The final call is the sequence call** unless a tested line moves it under the rules below.
5. **Fixed constants, each with a reason** (a convention, a majority rule, a value measured in the data, or a named
   heuristic). Nothing is tuned to a benchmark, and there are no threshold options.

Not goals: centromeres from methylation alone (false dips are as deep as true ones); Hi-C choosing between satellite
families. The LTR line is developed and benchmarked with LTRquest annotations only; the readers for DANTE_LTR and
EDTA + TEsorter files stay as they are, with their tests, but get no further work.

## Stages

```
satellite -> telomere -> rdna -> rdna5s -> gaps -> genotype -> hic -> fold -> karyotype -> provisional -> meth ->
fuse -> cpg -> plants -> morphology -> interpret -> loci -> output -> report -> plot
```

The chromosome set is chosen once per genome ([outputs.md](outputs.md#which-sequences-are-analyzed)) and every line
reads the same set. Each stage is checkpointed on its inputs and code, so a rerun resumes. `genotype`, `karyotype`
and `plants` handle inputs from several plants ([plants.md](plants.md)).

## Satellites (`cenfuse.sat`)

**From tandem repeats to a call.**

1. *Seeds.* FasTAN (or TRF) finds tandem repeats in 5-Mb chunks that overlap by 10 kb, so no array is lost at a cut.
2. *Families.* Seeds are clustered by sequence (BLAST) into families, and each family gets a consensus and a profile
   HMM (MAFFT alignment).
3. *Genome-wide search.* nhmmer finds every copy of every family.
4. *Roles.* Each family gets a role from where and how it occurs across the chromosomes: `centromeric`,
   `cen_associated`, `subtelomeric`, `interstitial` (knobs and other blocks), `dispersed`, `minor` or `rDNA` (rRNA
   genes, found with Rfam models).
5. *The call.* On each chromosome, the array cluster of the chosen centromeric family; where no satellite call is
   possible, a long-unit tandem array or a retrotransposon core.
6. *Higher-order repeats* are reported for each call.

| rule | test | why |
|---|---|---|
| family founding | a seed cluster becomes a family only if one of its seeds is a tandem array (≥ 3 copies of its unit), then by size (≥ 20 kb of seeds, at most 40 families) | a satellite is tandem. Clusters founded on two-copy duplications were dispersed or absent in every genome of a 216-genome replay, and took most of the search time where they occur (T2T sheep: 93%) |
| dispersed repeat | less than half in tandem, less than half whole-unit hits, and no locus concentration (fold < 2) → never centromeric, role `dispersed` | interspersed elements (LINE/SINE fragments) annotate as chains of partial hits; satellites are tandem, and retrotransposon-derived centromeric arrays are concentrated |
| satellite kind | dominance is judged per kind: subfamilies related by sequence are pooled | splitting alpha satellite into 6 subfamilies let a pericentric satellite outcount it |
| acrocentric | a strong family at one end of most chromosomes, and at both ends of too few for independent ends (binomial, P ≤ 0.01), with no strong interior rival → centromeric, flag `acrocentric`, at most `medium` | subtelomeric repeats occupy ends independently; a centromere occurs once per chromosome. Needs about 16 or more chromosomes |
| competing cluster | a second cluster at least half the size of the main one, or a second, unrelated satellite that also marks the acrocentric end → `low`, reported in `seq_alt_*` | such calls were right only 2/10 times. Which of two satellites at one end holds the kinetochore is not in the sequence (mouse: the major satellite is larger than the CENP-A minor satellite) |
| honesty cap | centromeric families exist but none is the largest cluster on most chromosomes, and no CRM core decides → every call `low` (`no_dominant_satellite`) | satellite evidence alone cannot tell a centromeric array from an interstitial one |
| genes (`genes/`) | gene-rich families are vetoed; subtelomeric = beyond the outermost genes | potato: 2 → 4/11 right; barley: the 3 confident wrong satellite calls become `low` |

## LTR retrotransposons (`tes/`)

cenfuse finds, for each lineage of LTR retrotransposon, the chromosomes where its elements concentrate at one locus.

- **Cores.** A lineage's core on a chromosome is its densest region (the full width at half maximum of its density),
  and it needs at least 3 independent insertion sites (overlapping elements count once). *Chlorella* centromeres
  hold only 2-7 Copia insertions each; 53 of the genome's 54 sit in the CENH3 domains.
- **The marker.** The lineage most concentrated at one locus per chromosome. A CRM marker is **decisive** when its
  cores are significant on most chromosomes: CRM is the known centromeric lineage of grasses and many other plants.
  Any other marker (Athila, Tekay, Bryco...) only flags (`te_concordant`, `te_discordant`) and may lift the honesty
  cap only for the calls it supports.
- **What a decisive core does.** It raises a call it agrees with to at least `medium`, takes a chromosome without a
  call, and replaces a weak call beside it or a disagreeing call that is `low` on its own. A `high` or `medium` call,
  or a non-weak call beside the core, it replaces only when the lineage's elements sit significantly more in the core
  than in the call (binomial at 1/2, P ≤ 0.01); otherwise the call stays, one level lower (`te_discordant`). The
  replaced array is reported in `seq_overruled_*`. Reason: in a broad CRM plateau the densest window is noise (*B.
  rapa* A09 and soybean chr1/chr6 had CENH3 on the satellite, 16-26 vs 14-32 elements), while a satellite beside a
  core that holds hundreds of elements is pericentromeric (*Ae. tauschii*: 3/7 → 5/7).
- **CRM short of a majority.** Where CRM is not decisive, each significant CRM core replaces the marker's core on its
  chromosome (flags only): with less sensitive annotations CRM concentrates on only a few chromosomes, but every such
  core lay on CENH3 (24/24 in rice, potato, maize, *Brachypodium*).
- **Formats.** An element's lineage comes from `clade=` (LTRquest, EDTA + TEsorter) or from a DANTE_LTR element's
  `Final_Classification`. DANTE_LTR's rank-D clusters (protein domains without LTRs) are not elements: in potato they
  made 195 of 223 CRM and put cores up to 15 Mb off CENH3. A file without lineages (EDTA alone gives superfamilies
  only) is reported and the line is off.

## rDNA and telomeres

- **45S rDNA (NORs):** tandem runs of ≥ 3 copies of the universal 18S motif, ≤ 100 kb apart (CHM13's 5 NORs in 32 s).
  Hi-C peaks beside a NOR never oppose a call: NORs meet in the nucleolus.
- **5S rDNA:** tandem runs of ≥ 3 whole 5S genes (the Rfam 5S profile HMM, ≥ 90% of the gene, E ≤ 1e-5), ≤ 20 kb
  apart. Whole genes, because pseudogenes and fragments are mostly partial; dispersed 5S genes form no array.
- **Telomeres:** the repeat unit is found from scratch; an end holds a telomere array with ≥ 300 bp of unit copies in
  the terminal 100 kb, or ≥ 100 bp that reach the sequence end (assemblers trim short telomeres: *Neurospora* has
  157-207 bp).
- **Interstitial telomeric arrays:** ≥ 300 bp of the unit inside a chromosome. They mark old fusions (human 2q13) and
  the cores of some centromeres (chicken).

## Hi-C (`cenfuse.hic` and `hicfuse.py`)

Centromeres cluster in many nuclei, so the inter-chromosomal (trans) contacts of a chromosome's centromere are
strongest with the other chromosomes' centromeres. Chromosome ends, the nucleolus and other hubs cluster too, which
is why Hi-C is tested before it may act.

1. **Bins.** 250 kb or coarser, chosen by the library's depth (at least 25 expected centromere contacts per bin). A
   chromosome needs 20 bins; for genomes with small chromosomes the bins go down to the coarsest size that gives the
   shortest one 20 bins (yeast and *Plasmodium* 10 kb; *S. pombe* and *Neurospora* 50 kb), within memory caps (at
   most 4,000 bins on the largest chromosome). Replicate libraries of one plant are summed first.
2. **Prior-free peak.** The Hi-C engine's own iteration, using no sequence call. It decides nothing by itself.
3. **Anchored score.** For each chromosome, its mean contact with windows at the *other* chromosomes' sequence calls
   (its own call never shapes its score). Where the window falls in an unmappable satellite array, the 9 mappable
   bins nearest the call (the array's flanks) stand in for it.
4. **Peaks.** A peak is *strong* when it passes the engine's three tests (bootstrap, leave one partner out, null
   replicates) and is not flagged as a `hub`, `rDNA` (beside a NOR) or `telomeric` (it cannot be told from the
   chromosome end). A flagged peak may support a call but never moves or opposes one: chromosome ends cluster in 3D
   whatever the centromere.
5. **Distances in mappable bins.** A peak agrees with a call within 1 bin, forgiving up to 4 unmappable bins between
   them (a peak at an array's edge supports a call inside it); it opposes only beyond 4 mappable bins.
6. **The test (gate).** For chromosome x, the references are the other chromosomes' `high` calls (with `medium` if
   fewer than 4 are `high`). x passes when the references' own peaks agree with their calls on most of them, beyond
   chance (Poisson-binomial P ≤ 0.01). Fewer than 3 references: `untested`. Where no call is `high` or `medium`
   anywhere, the `low` calls are the references, and only agreements at clean peaks count.
7. **Support tier.** Where agreements are beyond chance (P ≤ 0.01, at least 3) but not on most chromosomes, Hi-C
   carries centromere information, yet its peak is often elsewhere. It may then only support: a call whose own
   significant peak agrees goes up one level (`hic_support`); it never opposes, moves or places a call. Only
   agreements at peaks of strength `medium` or `high` count, because randomized libraries put weak peaks on 3-6 calls
   per genome at the edges of unmappable arrays.
8. **Hi-C alone.** In a genome with fewer than 3 reference calls (point, regional, AT-rich centromeres) nothing can
   test Hi-C against the sequence. There the Rabl axis decides (below): where it passed and agrees with the
   prior-free peak, the call is `medium` (`fold_pair`, `align_pair`). Elsewhere a clean `high` prior-free peak gives
   a `low` call (`hic_denovo_only`): the best candidate, untested. A sequence call is kept; a clean prior-free peak
   away from it is reported (`hic_denovo_elsewhere`).

## The Rabl axis: the arm fold and the trans arm alignment (`fold.py`)

In a Rabl nucleus the centromeres sit at one pole and the arms extend toward the other, so a locus's position in
the nucleus follows its distance from its own centromere. Two readouts of that axis can test Hi-C without any
sequence call:

- **Arm fold (cis).** Both arms of a chromosome run side by side, so its own contact map has an anti-diagonal ridge
  through the centromere. The fold score of a center is the mean observed/expected contact over bin pairs (c - d,
  c + d) at least 9 bins apart. Only centers with a centromeric index ≥ 12.5% are scanned (below that a chromosome
  is telocentric and the short arm cannot fold). The fold core is the top center's half-maximum run.
- **Trans arm alignment.** Arms of different chromosomes also lie side by side, so loci at equal distances from
  their own centromeres meet (the trans "X" of Hoencamp et al. 2021). The score of a center on chromosome x is the
  contact between x's bins and the other chromosomes' bins at the same distance from their prior-free peaks, with
  the centromere cluster itself left out, so only the arms count. It pools every partner chromosome, so it works
  where one chromosome's map is too small for the fold (yeast: alignment 13/16, fold 6/16).
- **The test.** The prior-free Hi-C peak lies in the core on most chromosomes, beyond chance (P ≤ 0.005 for each
  readout). The alignment stands in for a failed fold only where the fold leans the same way (P ≤ 0.05): the
  alignment is read relative to the prior-free peaks, and a holocentric genome whose compartments are symmetric about
  each chromosome's middle can mimic it (*C. elegans*: 2/5), while the fold reads each chromosome's own map.
- **What it does.** Where the gate fails, an axis line that passed and agrees with the prior-free peak can move a call
  that both contradict (barley: 6 calls moved, all onto CENH3).

## Methylation: CHG (`cenfuse.meth` and `methfuse.py`)

The CHG line learns, leave-one-chromosome-out, what methylation looks like at the confident calls of the other
chromosomes, and labels each candidate interval `centromere_like` when its signature score is ≥ 0.5. It supports and
can be a second line for a Hi-C move; it never moves a call alone, and only samples of the assembled plant (or of
unknown origin) decide. It needs `tes/` and a CHG file for every sample.

## CpG kinetochores (`cdr.py`)

In primates CENP-A sits in a CpG-hypomethylated dip of the methylated active array. The satellite and Hi-C lines
locate the array (0.8-5.6 Mb); the dip locates the kinetochore inside it.

- **Profile.** CpG methylation per 5-kb bin over the centromeric arrays inside the final call; a bedMethyl's calls are
  pooled per bin, so low coverage works (rhesus thinned to 1/16 of its reads kept 21/21 kinetochores within 50 kb).
- **Dip.** The most prominent valley (higher methylation on both sides, so an array edge is never a dip), extended to
  its half-depth run.
- **Tests.** The dip is longer than in the call's own bins shuffled (P ≤ 0.05), and than methylation from elsewhere on
  the chromosome rotated onto the call's CpGs makes (the positional test).
- **The test (gate).** For each chromosome, the other chromosomes' confident calls hold a positional dip more often
  than chance (binomial at 0.05, P ≤ 0.01).
- **Role.** Where the gate passes, each call with a significant dip gets its kinetochore (`cpg_kt_*`). It changes no
  call or confidence, and it cannot tell a wrong array: an inactive array has dips too.
- **Inputs.** The CG file of each deciding sample (modkit or ENCODE bedMethyl, or a bedGraph of methylated fraction);
  rows must sit on reference CpGs.

## Integration (`fuse.py`)

| situation | result |
|---|---|
| no tested line | the sequence call, `sequence_only` |
| tested Hi-C (or methylation) agrees | the sequence call, `concordant`, +1 level per line (+2 at most) |
| tested strong Hi-C peak more than 4 mappable bins away, call `low`/`medium` | moves to the Hi-C core, `offset`, `medium` |
| same, call `high` | moves only with a second line at the peak (a CRM core within 4 bins, or methylation `centromere_like` at the peak and a dip at the call): `high`. Otherwise `contested` (`hic_offset`), one level lower |
| Hi-C peak 1-4 bins away | neither supports nor opposes (`hic_near`) |
| competing cluster: a strong peak agrees with the second cluster, mappable DNA between them | the second cluster is chosen (`hic_chose_alt`) |
| competing cluster with no mappable bin between it and the call | no Hi-C support (`hic_alt_ambiguous`): Hi-C cannot tell the two apart |
| no sequence call, strong peak | the Hi-C core, `functional_only`, `medium` (a clean medium peak: `low`, `hic_medium_peak`) |
| gate failed, a Rabl-axis test passed, its core and the prior-free peak agree away from the call | moves to the prior-free core, `offset`, `medium` (`fold_pair`/`align_pair`) |
| same, but the axis or the prior-free peak opposes the kept call | `contested`, one level lower per opposing line |
| support tier only, the call's own significant peak agrees | `concordant`, +1 (`hic_support`) |
| no references anywhere, no sequence call, clean `high` prior-free peak | the prior-free core, `functional_only`, `low` (`hic_denovo_only`) |

CpG kinetochores are added after the integration and change nothing. Per genome, `hic_sequence_conflict` says when
Hi-C's prior-free peaks are confident but disagree with the sequence calls while the gate fails: one line is wrong
genome-wide, which only the Rabl axis, `tes/` or methylation can settle.

## Reading the integration (`morphology.py`, `interpret.py`)

Nothing here changes a call; every statement is read from the tables cenfuse writes.

- **Morphology.** The arms are measured from the kinetochore's midpoint where the CpG line placed one, else from the
  call's: centromeric index, arm ratio and Levan class (Levan, Fredga & Sandberg 1964). The summary's karyotype counts
  only `high` and `medium` calls. The classes are sequence measurements, not cytogenetic descriptors: CHM13's
  acrocentric chr21 and chr22 come out submetacentric/subtelocentric, because their T2T short arms (rDNA included)
  are 11-15 Mb.
- **Completeness.** Telomere arrays at the ends and assembly gaps (runs of ≥ 10 N) on the chromosome and in the call.
- **`basis` and `caveat`.** One sentence each: what the call stands on (the satellite and its size in the call, a CRM
  core and the array it replaced, a move and its second line, each agreeing line, the kinetochore); then what speaks
  against it (the overruled array, each lowering flag, each opposing line and where it points, gaps) and the input
  that would settle it, where cenfuse's own tests say one can (two arrays of one satellite: CENP-A/CENH3 data;
  satellite-free or capped calls: Hi-C or `tes/` when not given).
- **The report.** Per genome, the lines and their tests, the calls by confidence, the uncertain calls with their
  caveats, morphology and completeness, and the other features with their 3D association.

## Constants

Every threshold is fixed, with its reason. cenfuse's are in `cenfuse/params.py`; each engine keeps its own
(`sat/params.py`, `hic/params.py`, `meth/params.py`).

| constant | value | reason |
|---|---|---|
| HIC_TOL_BINS, HIC_FAR_BINS | 1, 4 bins | Hi-C placement error; Hi-C opposes only beyond its partner half-window (its resolution) |
| GATE_P, GATE_MIN_REF | 0.01, 3 | conventional evidence level; fewer references give no majority test |
| ANCHOR_LEVELS | high, medium | wrong anchors blur Hi-C and methylation |
| NOT_STRONG_FLAGS | hub, rDNA, telomeric | hubs and NORs contact everything; chromosome ends are the main non-centromeric trans cluster |
| MAX_RAISE | 2 | `low` → `high` at most |
| SAT_SUPPORT_BP | 5 kb | the satellite engine's weak-satellite floor |
| METH_SCORE | 0.5 | the methylation engine's signature rule |
| fold: S_MIN, CI_MIN; axis P | 9 bins, 12.5%; 0.005 (0.05 to lean) | beyond the ± 4-bin centromeric window; Levan's telocentric bound; 0.01 shared by two readouts (Bonferroni) |
| Hi-C bins | 5 kb-2 Mb ladder, ≥ 20 bins per chromosome, ≤ 4,000 per chromosome, ≤ 50,000 per genome | each size a multiple of every finer one; a centromere needs flanks; memory |
| sat: TANDEM_FRAC, ACRO_ALPHA, LOCUS_FOLD | 0.5, 0.01, 2 | majority rule; conventional level; heuristic |
| sat: TE_MIN_SITES, TE_ALPHA, TE_CONTRA_ALPHA | 3 sites, 1e-3, 0.01 | three independent insertions make a locus; conventional levels |
| sat: min_array_monomers, LEN_TOL | 3 copies, 5% | a tandem array; unit-length tolerance |
| chromosome set: GAP_FOLD | 10 | a tenfold step in length ends the chromosomes (heuristic) |
| CpG: BIN, MIN_SITES, MIN_BINS; N_PERM, N_ROT, DIP_P | 5 kb, 10 CpGs, 10 bins; 200, 200, 0.05 | CDR-Finder's window; a dip needs shoulders; P resolution 10× finer than the level |
| 5S: E, unit max, copies | 1e-5, 20 kb, 3 | the satellite search's E; twice the longest 5S-bearing unit (yeast rDNA); a tandem array |
| morphology: Levan bounds, gap | 1.7 / 3 / 7, 10 N | Levan et al. 1964; assemblers write gaps as ≥ 10 N |

## Choices tested and rejected

| question | chosen | rejected, and why |
|---|---|---|
| fix Hi-C from inside the Hi-C data | anchor to the sequence calls | random restarts and caps: accuracy unchanged; in potato every start converges on a non-centromeric cluster 3.5× more mutually connected than the CENH3 set |
| Hi-C partner windows over unmappable arrays | the 9 nearest mappable bins | padding every window: lost agreements in three apes |
| Hi-C peaks at chromosome ends | flagged only where a peak cannot be told from the end; such peaks support, never move | an 8-bin end window: moss centromeres 3-6 bins from the end could not move wrong calls. Contact-based tests (the other end's height, mirrored windows) did not separate right from wrong end peaks |
| Hi-C that finds centromeres beyond chance but not on most chromosomes | support agreeing calls only | no role (the agreeing peaks were right 61/64 times); letting it oppose or move (a disagreeing peak is weak evidence) |
| genomes without sequence calls | Hi-C's clean prior-free peaks as `low` calls, `medium` with a Rabl-axis line | a Rabl polarity test: passed on randomized Hi-C too |
| cis profiles (insulation, compartments), an end-distance expected for trans contacts | not used | informative in single genomes only (tomato, potato), with opposite signs and no per-genome test |
| Hi-C depth fragility | reported as depth guidance | a split-half check: unstable halves did not predict each other, and Hi-C time doubles |
| `.hic` files | stop with the `hictk convert` recipe | reading `.hic` natively: a new dependency for what one command does |
| a CRM core vs a confident distant or adjacent satellite call | replace only if CRM marks the core significantly more | always replace: wrong on *B. rapa* A09 and soybean chr1/chr6. A bootstrap of the core position: lost *P. annua* Pa1A |
| TE core floor | 3 independent insertion sites | 5 elements: missed *Chlorella*'s 2-7-insertion centromeres |
| LTR elements without a lineage (fungi) | not used | pooling them found *Cryptococcus* centromeres (10/14, `low`) but lifted the honesty cap in fish and *Oryza* genomes |
| several qualifying TE lineages | the most concentrated | chromosomes with a core first: agreed less with independent satellite calls |
| CpG dip inside a call | most prominent valley, half-depth run | deepest bin: lands on array edges. Longest run below half the median: on the CDR 52/90 (valley 88/90) |
| CpG null | positional (rotation) test for the gate | dips in other satellite clusters: passed on misplaced methylation in gorilla |
| CpG choosing between arrays | not done | least-methylated window: on the CDR 3/23; most prominent valley across clusters: wrong cluster on 8-10/25 |
| the centromeric region (call + nearby satellite clusters) | not reported | held out: precision 0.07 in apes (an interspersed family chains whole chromosomes) |
| chromosome set | tenfold drop, tapering tail cut, optional list | the drop alone kept 5,614 extra scaffolds over 121 polyploid assemblies; other cuts dropped dot chromosomes or microchromosomes |
| a sample screen to skip dispersed families | not used | no family was skipped on real genomes |
| a faster, inexact seed clustering | not used | 10× faster, but BLAST's hit lists differ for 0.2-0.5% of pairs, so it is not exact by construction ([future.md](future.md)) |
| GFF3 feature types | Sequence Ontology terms, cenfuse's kind as `kind=` | cenfuse's own kinds as types: no validator or pipeline accepts them |

## Known risks

- The Hi-C test takes the sequence calls as references. Where those are wrong genome-wide, only the Rabl axis,
  `tes/` or methylation can rescue Hi-C (potato and barley Hi-C anchored at the true centromeres find 7/11 and 6/7 of
  them; anchored at the calls, none).
- Hi-C alone trusts the strongest trans cluster to be the centromeres. In a genome without sequence calls whose
  dominant cluster is something else, its `low` calls are wrong; a holocentric genome without sequence calls would get
  them too.
- One wrong Hi-C peak passes in a genome where Hi-C mostly works; only `high` calls need a second line to move.
- Hi-C tests lose power gradually but flip at their threshold (both orangutans at half depth); a failed test reverts
  to the sequence line.
- The decisive CRM test needs a sensitive annotation (potato: 214 CRM elements give exactly a majority).
- Satellite-free centromeres beside a large pericentromeric satellite are beyond the sequence lines (tomato).
- Hi-C resolution is a mappable bin (10 kb-2 Mb); offsets inside an unmappable array are invisible.
- The CpG line was designed on the primates it is benchmarked on; rhesus, *Arabidopsis* and HG002 were held out.
- Morphology is only as right as the call; an incomplete chromosome shifts the index.
- Satellite outputs of earlier versions read with `--satellite-dir` keep that version's rules; rerun with
  `--force satellite` to apply the current ones.

## Engineering notes

- **Reproducibility.** FasTAN's array cuts depend on memory layout, so it runs with address-space randomization off
  (`setarch -R`); MAFFT's refinement runs single-threaded (`--threadit 0`). Fresh runs are then identical.
- **Memory.** One pyhmmer search pipeline per thread, reused: a new pipeline per search held its workspace until
  Python's garbage collection (chicken: 277 → 58 GB).
- **Parallel family models.** A family's model is built in the frame of the most similar earlier family, so models
  are built concurrently but accepted in family order, each rebuilt if an earlier family it did not see could have
  changed it: the result is the serial one (*Panax*: families stage 3.4× faster at 16 threads).
- **FasTAN crashes** are recovered by splitting the piece ([outputs.md](outputs.md#running)).
- **Resume guard.** The satellite checkpoints stand on the chromosome set, which is recorded before the engine runs; a
  changed set stops the run instead of reusing seeds made on other sequences.

## Tests

`pytest -q -n 8`: 774 tests on small synthetic fixtures (about 3 min). They cover every rule above with unit tests; a
synthetic genome with planted cases (a Hi-C move with a CRM core, the same without a second line, a hub decoy, a
satellite-free chromosome, randomized Hi-C, three plants); a full satellite-engine run on its synthetic genome;
planted CpG dips and every methylation format; every Hi-C pairs layout; a yeast-like Rabl genome called at 10-kb
bins; the chromosome-set rule, the list file and the resume guard; FasTAN crash recovery; parallel family models
byte for byte equal to the serial build; the report's sentences.
