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

The satellite engine checkpoints its stages on stage version and parameters: a change stops the run and names the
`--force` stage, and changed input files (FASTA, `genes/`, `tes/`) are not detected (after adding `tes/` or `genes/`
to an earlier run: `--force classify`). Its checkpoints also stand on the chromosome set, which cenfuse checks itself
(below). telomere, rdna, rdna5s, gaps, genotype, hic, fold, karyotype, meth and cpg are
keyed on their inputs and code, so a rerun resumes; the later stages are cheap and always rerun. Their keys include
the code of their modules and `params.py`: a constant used by one new module lives in that module, so adding it does
not recompute every stage. `genotype`, `karyotype` and `plants` handle several plants ([plants.md](plants.md)).
In 0.9 the satellite engine keeps its version (0.9.0) and its stage versions, so nothing forces a rerun, although its
summary and region tables gained the `overruled_*` columns; the standalone methylation engine (`python -m
cenfuse.meth`, 0.1.0) now counts chromosomes with the 0.9 rule (`track.chromosome_count`).

### Chromosome set (`inputs.py`, `meth/track.py`)

An assembly also holds unplaced scaffolds, and every line needs whole chromosomes (the satellite engine's prevalence
tests count them; Hi-C and methylation need chromosome-scale sequences). The set is chosen once per genome, and all
three engines read it through a view of the genome folder (the FASTA is a symlink and its `.fai` lists the
chromosomes only).

- **Rule** (`track.chromosome_count`, GAP_FOLD = 10). Start from 0.8's set (`tenfold_count`): the sequences above the
  first >= 10x drop in sorted length between neighbors, looked for once half the assembly is accounted for; all
  sequences when there is none. A gradual scaffold tail has no tenfold step, so that set takes it whole. The tail test:
  take the first sequence, after half the assembly, shorter than a tenth of the median length of the sequences before
  it (`_median_cut`); if the sequences from there to the end of 0.8's set span a >= 10x range, they are a tail and the
  set ends before the first of them. A few small sequences of similar size before a clean gap (dot chromosomes,
  microchromosomes) span less and stay, as in 0.8.
- **List** (`<genome>.chromosomes.txt`, named after the FASTA, in GENOME_DIR): the chromosomes, one sequence name per
  line. A line starting with `#` is a comment (a `#` inside a name, as in `HG002#1#chr1`, is part of it); blank lines,
  CRLF line ends, trailing spaces and a UTF-8 byte-order mark are accepted. An unknown name, a duplicate or no names
  stops the run (a list is the user's decision, never guessed at); a list that matches no genome logs a warning. It
  replaces the rule for every line, also for a `--classify-from` target, which otherwise keeps 0.8's rule (a contig
  assembly has no chromosome/scaffold distinction for the tail test to find). It must name every chromosome: the
  methylation line's expected domain size is D = d_formula(G, number of chromosomes), so a partial list inflates D.
- **Source**: the length rule, the 0.8 rule of a `--classify-from` target, or the list's name is logged, with the number
  of sequences left out, and shown in the report's Notes. When the tail cut fired (`chromosome_count` below
  `tenfold_count`), the note says the set ends at a gradual tail of short sequences, with the number left out and the
  shortest one kept, and does not call them unplaced scaffolds (a dot chromosome can be among them): if chromosomes were
  left out, list every chromosome in `<genome>.chromosomes.txt`. The cut fires on none of the 234 benchmark and batch
  assemblies and on 16 of the 149 polyploid ones.
- **Evidence** (383 assemblies: 49 benchmark genomes, 185 batch genomes, 149 polyploid assemblies). The set is 0.8's on
  all 234 benchmark and batch assemblies (chicken's 40 chromosomes with 2.5-Mb microchromosomes included). On the 121
  polyploid assemblies that come with a list of their chromosomes, extra sequences fall from 5,614 to 31 (in 12
  assemblies) and none of the listed chromosomes is lost; an exact match rises from 102 to 109. The 31 left are
  1.5-117 Mb sequences outside the lists that a length rule cannot tell from chromosomes.
- **Resume guard** (`Run._check_chromosome_set`). The satellite engine's resume check covers its code and parameters,
  not the chromosomes, and 0.9 changes the set (gradual tails, a list). `satellite/work/chromosome_set.json` records
  {genome: sorted names}. cenfuse writes it once the engine is set up and before it runs, since every checkpoint,
  present or future, then stands on this set. A resumed folder whose set differs (a genome added, removed or changed)
  stops before anything is written, names the genomes with their chromosome counts, and says `rerun with --force
  satellite` (`--force seeds` also proceeds: it redoes the seeds, the first stage the set shapes). A folder made before
  0.9 has no record: each current genome is compared with what 0.8's rule gives on its `.fai`, so a genome added or
  removed since is detected only once a record exists; one without `seeds.done` has nothing to protect.
  `--satellite-dir` imports are taken as they are and only warn when the imported summary's chromosomes differ from
  the set.

## Lines

### Satellites and TEs (`cenfuse.sat`)

The pipeline: tandem-repeat seeds (FasTAN or TRF) → families (clustered, consensus) → genome-wide annotation
(nhmmer) → a role per family → one call per chromosome → HOR structure.

Roles: `centromeric`, `cen_associated`, `subtelomeric`, `interstitial` (knobs and other blocks), `dispersed`,
`minor`, `rDNA` (rRNA genes found with Rfam models).

How the stages run (none of this changes a call):

- **Seeds.** FasTAN (or TRF) runs on 5-Mb chunks that overlap by 10 kb (at least twice the longest seed unit, so no
  array is lost at a chunk edge). FasTAN 0.8 aborts on some large satellite arrays (a 642-kb array of *Oryza minuta*:
  SIGABRT, `malloc(): corrupted top size`, at the same place on every run). A piece on which it dies by SIGABRT,
  SIGSEGV, SIGBUS, SIGFPE or SIGILL is cut in two halves that overlap by 10 kb, and each is retried, recursively
  (`fastan.chunk_rows`); `seeds.dedupe` merges what both halves report, as for chunk overlaps. A piece of at most 20 kb
  (2 x the overlap) that still crashes is skipped with a warning. A chunk longer than that of which no piece was
  searched (FasTAN crashed on every piece down to 20 kb; pieces without bases never reach FasTAN) stops the run
  (`fastan.ChunkFailed`): FasTAN's known crash is local to one array (the *Oryza minuta* chunk lost 0 bp after 4
  splits), so a whole chunk failing means FasTAN itself is failing, and a run without the chunk's seeds could pass for
  a satellite-free genome. The error names the chunk, the signal and the crash table, written first with the genome's
  crashes so far, and suggests `--seed-engine trf`; no seed table is cached. Each crash (chromosome, start, end,
  signal, `split` or `skipped`; 0-based) goes to `<seed cache file>.fastan_crashes.tsv` beside the cached seeds, so a
  shared `--seed-cache` keeps it, linked as `satellite/work/seeds/<genome>.fastan_crashes.tsv`; the genome's summary
  line in the log counts the chunks crashed, the pieces split and the bp skipped. SIGKILL (the OOM killer), SIGTERM,
  SIGINT and every other failure still stop the run. A piece FasTAN does not crash on gives the rows cenfuse 0.8 gave,
  so FasTAN's cache key (`fastan.KEY`) and the seeds stage version are unchanged; the genome part of the key changed
  only for a FASTA indexed beside a symlink (seed cache key, below). An *Oryza minuta* 5-Mb piece that aborts FasTAN
  gives 757 seed rows after 4 splits.
- **Seed clustering** (`cluster_seeds`, `absorb_multimers`). BLAST is 97.7% of the wall time, driven by the HSP volume
  of multimeric seeds; Python is about 2%. *Acanthus* (5,456 seeds) took 76 min at 16 threads, against a median of 0.7 min
  over the 106 genomes of an allotetraploid plant set. No speedup that keeps the clusters exact was found, so 0.9 only
  logs each round and the seeds left at INFO ([future.md](future.md)).
- **Families.** A family's model is built in the monomer frame of the most similar earlier family (`frame_rotation`;
  references within `RELATED_LEN_TOL` = 15% in length), so models are not independent, and `build_models` builds them
  concurrently with a guarantee: the result is the serial one. A build starts with the families already accepted as
  references (final ones only) and records the length and the rotation of each `frame_rotation` call, the only steps
  that read the references. Families are accepted in family order, each with all earlier families final. A finished
  build is kept if no earlier family it started without is a frame candidate (within 15% in length) for one of its
  queries, or if replaying its recorded rotations with all earlier families gives the same rotations (the build is
  then the same computation step for step); otherwise it is rebuilt with all earlier families. Build k samples with
  `random.Random(k + 1)`, so a model never depends on another's sampling. Up to threads // `THREADS_PER_BUILD` (4)
  builds run at once, no more than there are families, and the threads are shared out among them (16 threads and 2
  families: 2 builds of 8; 16 and 40: 4 of 4); a failed build raises once the running builds end. families.tsv and
  .fa, families.hmm without its DATE and COM lines, and seeds.tsv.gz are byte-identical to the serial code on
  *Arabidopsis* Col-CC, rice T2T-NIP and *Panax japonicus*. Families stage, 0.8 (models one at a time) -> 0.9 at 16
  threads: *Panax* 2,290 s (0.8 at 16 threads) -> 674 s, 3.4x, the one like-for-like pair (model phase
  1,755 -> 402 s); Col-CC 230 s (0.8 at 8) -> 130 s; rice 138 s (0.8 at 12) -> 87 s. Which families a build started
  with depends on timing, so the BLAST scratch folders of the stage (`fam*/frame*/`, `fam*/replay*/`) can differ
  between runs and from 0.8's; the models, and so every results file, cannot.
- **Annotate.** nhmmer searches every (chunk x model) job on a thread pool. The search is the longest stage and hard to
  predict (*Panax*: 7.9 h predicted from family length x genome size, more than 25 h actual), so it logs at INFO each 10%
  of its work (chunk bp x model length, counted as searches finish; the time left assumes the rest costs per unit what
  the done part did, so it is an estimate) and, checked as each search finishes, after an hour without a line. Each
  finished search also adds the thread CPU time and wall time its nhmmer call took, and the hits kept, to its family's
  totals: `satellite/work/annotate/<genome>.search_cost.tsv` (family, model_len, searches, cpu_h, cpu_share, wall_h,
  hits_reported, monomers_kept; costliest first; wall_h sums the searches' own times, the threads' busy time). CPU
  times differ between runs, so this is a work table, never a results table: `satellite_families.tsv` is unchanged and
  this is the one run-dependent table (README, Reproducible). The log names the three costliest families when the
  search ends, and the output stage names them again with their final role (a dispersed family that cost much decides
  nothing, which only the roles, known after annotate, show). The cost sits in short-unit families with millions of
  hits in mammals (CHM13's 172-bp alpha satellite: 62% of the CPU) and in 1-2-kb families in plants and chicken.
- **Replaced calls.** Every call that a TE core took the place of (`centromere.te_call(core, call)`) is kept: the
  satellite summary gains `overruled_start`, `overruled_end` (1-based start) and `overruled_family` (its final name),
  and centromeres.tsv `seq_overruled_start`, `seq_overruled_end` right after `seq_alt_end`. That is, in a genome that
  is not borrowed (`--classify-from`): with a decisive (CRM) core, a weak call (weak-satellite or long-unit evidence, or
  less than `low_sat_bp` of satellite) that agrees with the core without overlapping it (`weak_satellite_beside_core`),
  a satellite call that agrees with it without overlapping it and holds significantly fewer of the clade's elements
  (`satellite_beside_core`), and a satellite or weak-satellite call that does not agree with it, unless the call is
  `high` or `medium` on its own and the clade does not mark the core significantly more (`satellite_away_from_core`);
  and, with a core of any clade, decisive or not, a long-unit call that does not agree with it (no flag). Agreeing is
  a gap of at most the wider of the two intervals; a long-unit call is made only on a chromosome without a satellite
  call, and with a core there the long-unit scan covers the core and `lu_flank` (500 kb) on each side. The report's
  `basis` sentence names the interval for the three flagged branches, and the `caveat` names the array for all of them
  ("the tandem array at X-Y Mb was overruled by the CRM core"). They are never written to `seq_alt_*`, which validated
  Hi-C can move a call to (`hic_chose_alt`), so they change no call. A kept call (`te_concordant`, `te_discordant`)
  and a chromosome without a call record nothing. Folders from 0.8 without a rerun of the centromeres stage read `.`,
  also on a chromosome whose flags say a core replaced a call. Replays of the centromeres stage on 3 stored runs (49
  chromosomes): every 0.8 column identical; the intervals equal an earlier audit's on all 16 chromosomes it covers.

| rule | test | why |
|---|---|---|
| family founding | a seed cluster becomes a family only if one of its seeds is a tandem array (>= 3 copies of its unit by FasTAN, 5% short allowed); then by size (>= 20 kb of seeds, at most 40) | a satellite is tandem. Clusters founded on two-copy duplications only (0.8-2-kb units, retrotransposon-like) were dispersed or absent in every genome of the 216-genome replay (12 families in 7 genomes) and took most of the genome-wide search where they occur (T2T sheep 93%, siamang 52%, orangutans 19-28%, barley 14%); without them calls and roles were unchanged (siamang, barley, B. orangutan); search CPU at equal machine load: sheep 11.7x, siamang 2.0x, B. orangutan 1.3x, barley 1.16x less |
| dispersed repeat | tandem_frac < 0.5 **and** whole-unit hits < 0.5 **and** locus fold < 2 → never centromeric, and role `dispersed` (never a satellite block) whatever its concentration | interspersed elements (LINE/SINE fragments) annotate as chains of partial hits. Satellites are tandem; retrotransposon-derived centromeric arrays are concentrated (fold ≥ 2.5). Before v7 the role ignored the test: one 0.5-1.9-kb family per genome (chimpanzee, bonobo, S. orangutan, rhesus, sheep; 215-719 Mb, 2-14% in tandem) was an `interstitial` block chained along whole chromosomes, making chromosome-scale loci |
| satellite kind | dominance judged per kind: subfamilies related by sequence are pooled | splitting alpha into 6 subfamilies let a pericentric satellite out-count it |
| acrocentric | a strong family at one end of most chromosomes, and at both ends of too few for independent ends (binomial, P ≤ 0.01), with no strong interior rival (judged per satellite kind) → centromeric, flag `acrocentric`, at most `medium` | subtelomeric repeats occupy ends independently; a centromere occurs once per chromosome. Needs ~16+ chromosomes. A subfamily at the inner edge of a terminal block is that block, not a rival (T2T mouse: a major-satellite subfamily beyond the 5% end zone had blocked the test) |
| the call | the cluster holding classify's chosen family (layered families pooled); a decisive CRM core may favor the largest cluster | |
| competing cluster | a second cluster ≥ half the main one, or the terminal block of a second acrocentric satellite unrelated to the call → two levels down (`low`), reported in `seq_alt_*` | such calls were right 2/10 times on CENP-A/CENH3 truth. The acrocentric test finds a centromere's end, not which of two satellites there holds the kinetochore (mouse: the 234-bp major satellite out-sizes the 120-bp minor satellite that carries CENP-A; related layers, as in goat, are pooled instead) |
| CRM cores (`tes/`) | the clade's densest cores (full width at half maximum of its density; any clade's core needs >= 3 independent insertion sites, overlapping elements counting once); any marker's core takes a chromosome without a call (`medium` if decisive, else `low`). A decisive (CRM) core supports a call it agrees with (`te_concordant`, at least `medium`) and replaces a weak call beside it (`weak_satellite_beside_core`) and a disagreeing call that is `low` on its own; a disagreeing `high`/`medium` call, or a non-weak call beside it, it replaces only if the clade marks the core significantly more than the call (its elements in the core vs in the call, binomial at 1/2, P <= 0.01; `satellite_away_from_core` / `satellite_beside_core`); otherwise a disagreeing call stays, `te_discordant`, one level down, and a call beside it stays `te_concordant` | CRM retrotransposons mark grass centromeres. Site floor: *Chlorella* NS4-2's Osser (Copia) cores, 2-8 elements per chromosome, 53/54 inside CENH3: 0 -> 8/10 right (all `low`); no call changed in 12 crop `tes/` arms, B. rapa, soybean, sorghum. Head-to-head: in a broad CRM plateau the densest window is noise (B. rapa CC1/CC2 A09, soybean X26 chr1/chr6: CENH3 on the satellite, 16-26 vs 14-32 elements); a satellite on a smaller second cluster loses (P. annua Pa1A 19 vs 62; Ae. geniculata 7M/2U) |
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
   - References: the other chromosomes with `high` calls (`high` + `medium` if fewer than 4 calls are `high`).
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
- **Gate** (per chromosome, leave one out): references = the other chromosomes' `high` calls (`high` + `medium` if
  fewer than 4 are `high`). k of n have a positional
  dip; binomial tail at 0.05 (each reference's chance under the null) ≤ 0.01. No majority is required: the
  positional test is conservative (passing primates: 28-95% of tested calls; rotated tracks: 4-9%).
- **Role**: where the gate passes, every call with a significant dip test gets its kinetochore (`cpg_kt_start/end`,
  GFF3 `regional_centromere_central_core`, Name=kinetochore). It never moves a call or changes its confidence, and it
  cannot tell a wrong array.
- **Assembly check**: bedMethyl rows must sit on reference CpGs (as for CHG); bedGraph rows too, at the row start or
  one base off (public tracks differ by one base).

## Fusion rules (`fuse.py`)

| situation | result |
|---|---|
| no validated line | sequence call, `sequence_only` |
| validated Hi-C (or methylation) agrees | sequence call, `concordant`, +1 level per line (max +2) |
| validated strong Hi-C peak > 4 mappable bins away, call `low`/`medium` | moves to the Hi-C core, `offset`, `medium` |
| same, call `high` | moves only with a second line at the peak and not at the call: a CRM core within 4 bins, or methylation `centromere_like` at the peak and a measured dip at the call (`moved_with_te`/`_meth`, `high`). Otherwise `contested` (`hic_offset`), -1 |
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
- **basis / caveat.** One sentence each per call: the sequence evidence and its size in the call, a CRM core that
  replaced the satellite call (with the interval of the call it replaced), a move and its second line, every agreeing
  line, the kinetochore; then the array a TE core overruled, with where it lies (`seq_overruled_*`), every flag that
  lowered the call, every opposing line with where it points, gaps in the call, and the input that would decide where
  cenfuse's own tests say one can (two arrays of one satellite: CENP-A/CENH3 data, since Hi-C needs mappable DNA
  between them and CpG methylation was tested and does not choose; satellite-free or capped calls: Hi-C or a CRM core
  when not given).
- **Report** (`cenfuse_report.md`): per genome, the lines and their tests, the calls by confidence and what they stand
  on, the calls that are not `high` with their caveats, morphology and completeness, the other features with their
  satellite share and 3D association. A failing report section is a warning; the results stay complete.

## Loci (`loci.py`)

Every candidate interval becomes a locus member:
- satellite clusters ≥ 20 kb of feature roles;
- TE cores;
- Hi-C cores and hubs;
- methylation candidates;
- 45S and 5S rDNA arrays;
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
| satellite engine (`sat/params.py`): TANDEM_FRAC, ACRO_ALPHA | 0.5, 0.01 | majority rule; conventional level |
| sat: LOCUS_FOLD | 2 | heuristic (as the Hi-C engine's hub fold) |
| sat: TE_MIN_SITES, TE_ALPHA, TE_CONTRA_ALPHA | 3 sites, 1e-3, 0.01 | nested copies are one insertion, three independent ones make a locus; conventional levels |
| sat: min_array_monomers, LEN_TOL | 3 copies, 5% | a tandem array (also founds a family); unit-length tolerance |
| chromosome set (`meth/params.py`): GAP_FOLD | 10 | a tenfold step in sorted length ends the chromosomes, and a tail is a run under a tenth of the median before it that spans a tenfold range: heuristic, the methylation engine's factor reused |
| sat: seed_overlap, FasTAN skip floor | 10 kb, 2 x overlap | the overlap is >= 2 x SEED_MAX_PERIOD, so no array of a seed-length unit is lost at a cut; at 20 kb the overlap is half the piece and each half three quarters of it, so the recursion stops there (heuristic) |
| sat: THREADS_PER_BUILD, RELATED_LEN_TOL | 4, 15% | the fewest threads per build while families remain: MAFFT is 3.6x faster at 4 threads and 5.6x at 16 (*Panax* fam006, identical alignment), so 4 builds of 4 beat 1 of 16, and spare threads go to the builds; length tolerance of related families (subfamilies with indels) |
| sat: progress step, heartbeat | 10%, 1 h | a line per tenth of the work; a line after an hour without one, checked only when a search finishes (so while every thread is busy with a long search the log can stay quiet for longer) |
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
| FasTAN reruns (reproducibility) | FasTAN under `setarch -R` (address-space randomization off) where the OS allows | as is: FasTAN 0.8's cuts depend on memory layout (one moss 5-Mb chunk, 5 runs: 17/982 arrays cut differently), so consensuses, families and some calls varied between fresh runs (moss: 5 call boundaries, one by 0.9 Mb); off: 5/5 runs identical; two whole fresh moss runs identical apart from header dates and commands |
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
| chromosome set | 0.8's set (the sequences above the first >= 10x drop between neighbors, after half the assembly), ended before a gradual scaffold tail (the sequences from the first one after half the assembly under a tenth of the median length of those before it, when they span >= 10x); a `<genome>.chromosomes.txt` replaces the rule | 0.8 alone: a gradual tail has no tenfold step, so it took the whole tail (*Cissus*: 646 extra sequences for 24 chromosomes; 5,614 extra over 121 polyploid assemblies), which dilutes the prevalence tests (rerun on the assemblies' own chromosome lists: *Cissus* 0 -> 14 `high` calls with all 24 positions unchanged, *Mimulus* 9 -> 25). The median cut alone: 27 extras, but on synthetic karyotypes it drops dot chromosomes that 0.8 keeps (*D. melanogaster* chromosome 4 beside the Y: 7 -> 6; a chicken with five 0.3-1.5-Mb dot chromosomes added: 45 -> 40). The mean instead of the median: 33 extras and a chicken microchromosome lost (40 -> 39). The largest neighbor ratio after half the assembly: 106 extras, 3 listed chromosomes lost, 11 of 234 benchmark and batch sets changed |
| `--classify-from` targets | keep 0.8's rule | the tail test: a target is a contig assembly, with no chromosome/scaffold distinction for it to find |
| chromosome-set record for the resume guard | written once the satellite engine is set up, before it runs | after the engine finishes: a run that fails after its seeds leaves seeds made on the new set beside the old record, and a return to the old set then passes the guard and resumes on them; and a first run that fails after its seeds leaves no record, so the next run is compared with 0.8's rule and can stop although nothing changed |
| seed cache key | the chromosome set joins the key whenever the real FASTA has no `.fai` of its own to compare with | the FASTA md5 alone, as before: when the user's FASTA is a symlink indexed beside the link, two chromosome sets of one file shared an entry, so `--force satellite` (the guard's own advice) reused the old set's seeds. Such folders compute their seeds once more; with the real file indexed (the standard layout) every key is unchanged (4 stored runs checked) |
| FasTAN crash (SIGABRT on a large array) | split the piece in two halves overlapping by 10 kb and retry, recursively, on a crash signal only (SIGABRT, SIGSEGV, SIGBUS, SIGFPE, SIGILL); skip a piece of <= 20 kb that still crashes, with a warning; record every crash; stop when no piece of a chunk longer than 20 kb was searched | 0.8: the first crash stops the run (two genomes of an allotetraploid plant set died, *Sonchus oleraceus* and *Oryza minuta*, both with the crash inside a high-confidence centromere call, so complete assemblies are hit hardest). Splitting on any failure: SIGKILL is the OOM killer and SIGTERM or SIGINT the scheduler or the user, not FasTAN's bug, and a retry on smaller pieces would hide the cause. Skipping without a ceiling: if FasTAN crashed on every piece after its self-test passed, the run would log a skip per piece and go on with almost no seeds, a satellite-free result that reads as biology; the known crash is local (*Oryza minuta*: 0 bp lost after 4 splits; the local patch that recovered both genomes stopped on a crash of any piece under 200 kb, and neither genome had one) |
| family models built in parallel | concurrently, at least 4 threads each (with fewer families than threads // 4, the threads are shared out among them); a build is kept if no earlier family it started without is a frame candidate for its queries (length within 15%), or if its recorded rotations replay unchanged with all earlier families; else rebuilt with all of them | rebuilding on every flag (the length test alone): 33 of 40 *Panax* families flagged, 11 for real and the rest length coincidences, and the rebuilds ran one after another in family order: model phase 1,927 s against 1,755 s serial. One thread per build (threads // workers): MAFFT is 3.6x faster at 4 threads and 5.6x at 16. Exactly 4 threads per build: with fewer families than builds the spare threads idled (16 threads, 2 families: 8 idle) |
| seed clustering speedup | none in 0.9; the rounds are logged | aligning only the pairs `seed_match` can consult (length-banded, N-padded databases, `-dbsize` pinned): 10x on *Acanthus* (266 + 69 s against 2,223 + 1,123 s) with identical clusters on 8 genomes, but BLAST's hit lists differ for 0.2-0.5% of those pairs (borderline hits near E 1e-5 depend on the query batch), so it is not exact by construction. Queries as single-thread jobs on the same database: 1.6-1.9x per call, hits differ (3 of 873k). `-mt_mode`: hits identical in 12 of 13 reruns, times mixed. [future.md](future.md) |
| CPU time per family | a work table (`search_cost.tsv`) and log lines | columns of `satellite_families.tsv`: a results file would differ between fresh runs |
| array replaced by a TE core | reported in `seq_overruled_*` (satellite summary: `overruled_*`) and named in the caveat | `seq_alt_*`: fusion lets validated Hi-C move a call to the alternative cluster (`hic_chose_alt`), so the replaced array would become a candidate and change calls |
| overview pages | drawn, saved and closed one at a time | all pages built first (0.8): a 670-sequence assembly held 28 figures open (matplotlib warns at 20) and 243 MB; one at a time: 1 figure and 144 MB, the same PDF and PNG |

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
- The chromosome-set rule cannot tell an unplaced scaffold from a short chromosome by length. The tail cut's own
  failures, where 0.8 kept the chromosomes: dot chromosomes followed by a gradual scaffold tail are cut with the tail (a
  bird-like karyotype of 37 chromosomes plus 300 scaffolds tapering from 0.4 Mb: 0.8 keeps all 337 sequences, 0.9
  keeps 31 and drops six 0.5-1.4-Mb dot chromosomes; with a clean gap before the scaffolds both keep 37), and small
  chromosomes that themselves span a tenfold range are cut as a tail. Shared with 0.8: a dot chromosome more than
  tenfold shorter than every other chromosome with nothing between stays out (*D. melanogaster* X, 2, 3 and 4 without
  the Y), and large unplaced scaffolds can stay in (31 sequences in 12 of 121 polyploid assemblies). A list in
  `<genome>.chromosomes.txt` settles any of them, and must name every chromosome.
- FasTAN crash recovery skips a piece of at most 20 kb that still crashes; only a chunk of which no piece was searched
  stops the run, so pieces skipped here and there in many chunks would pass with a warning each. The genome's summary
  line counts the bp skipped (an upper bound: skipped pieces overlap their neighbors).
- The time left in the progress lines is a linear extrapolation of the work done, and the work (bp x model length)
  ignores that families differ in cost per unit (*Arabidopsis* Col-CC: one 502-bp family took 75% of the search CPU),
  so it can be far off: at 8 threads the line at 90% said about 46 s left and the last tenth took 5 min. The hourly line
  is checked when a search finishes, so while all threads are busy with long searches the log can stay quiet for longer.
- Satellite outputs of earlier versions read with `--satellite-dir` keep their version's rules: before v7 an
  interspersed family classed `interstitial` (one 0.5-1.9-kb family each in chimpanzee, bonobo, S. orangutan, rhesus
  and sheep) chains into chromosome-wide clusters and chromosome-scale loci; before v8 the TE floor is 5 elements, a
  CRM core replaces every confident distant call, two-copy duplications can found families, and FasTAN's cuts vary
  between runs; before 0.9 the chromosome set was 0.8's and a replaced call was not recorded (`overruled_*` read as
  empty). A rerun in an earlier version's output folder stops at the first changed stage (`seeds` from v7): rerun the
  satellite engine (`--force satellite`) to apply the current rules. A 0.8 folder keeps its satellite checkpoints (0.9
  changed no satellite stage version) unless its chromosome set changed; cenfuse's own stages carry cenfuse's version
  in their keys and are computed once more after the upgrade.
- Morphology is only as right as the call: a wrong call gives a wrong class (`low` calls are kept out of the summary's
  karyotype). An incomplete chromosome (a missing telomere, gaps) shifts the index.

## Tests

`pytest -q -n 8`: 704 tests on small synthetic fixtures (about 2 min).
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
- Chromosome set: the length rule on synthetic length vectors (a gradual tail, a deep tail, dot chromosomes before a
  clean gap, the documented limit), equal to the rule written out in full on 600 random vectors, a million contigs
  counted in one pass; the list file (hand-written and Windows files, unknown and duplicate names, an empty list,
  `--classify-from` targets); the resume guard, with a stand-in for the satellite engine that writes only its seeds
  checkpoint (same set, changed set, no record, forced reruns, a run that fails after its seeds, `--satellite-dir`
  imports); cache keys of both FASTA layouts; the note when the tail cut fires, and the 0.8 wording otherwise.
- FasTAN crashes (a stub that dies by signal): the recursion, absolute offsets, the 20-kb floor, pieces at a
  sequence's edges and mostly N, kills and exit statuses that still stop the run, the crash table beside the cache; a
  chunk crashing on every piece stops the run with the chunk and the crash table named, while one crashing piece among
  searched ones and a short crashing sequence do not.
- Parallel family models: the acceptance and replay rules, output byte for byte equal at 1 (the serial build), 4 and
  16 threads, a failing build that raises instead of hanging, no reference that is not accepted; the thread split.
- Progress and cost: each tenth once, the hourly line, no line for no work, thread safety; one cost row per family and
  shares that sum to 1, no lost search under concurrent callbacks, `satellite_families.tsv` columns those of 0.8, and a
  malformed cost table that only warns while the output stage completes.
- Plots: a genome with more pages than matplotlib's warning limit ends with one figure open at a time, and none after
  an error.
- Replaced calls: each branch of the CRM replacement records its interval, `.` for outputs of 0.8, and the basis and
  the caveat name it.
