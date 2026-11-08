# How cenfuse works

cenfuse reads several kinds of evidence (*lines*) that each point to the centromere: satellite arrays, LTR
retrotransposons, AT-rich sequence, Hi-C contacts and DNA methylation. Each line calls on its own first. A line may then
change the final call only after it has passed a test in the genome at hand, and every agreement and disagreement is
reported. This page describes each line, the tests, how the lines are combined, and the alternatives that were tried and
rejected. Benchmark numbers quoted here are from [benchmark.md](benchmark.md).

## Principles

1. **Each line calls on its own first.**
2. **Lines check each other.** No line is trusted on its own confidence: the Hi-C calls the Hi-C engine rates `high` on
   its own lie within 4 bins of CENP-A on only 14 of 63 T2T mammal chromosomes.
3. **Whether a line works in a genome is measured, not assumed.** Hi-C must find the *other* chromosomes'
   centromeres before it may act on one; CpG dips must mark them before they place a kinetochore; a satellite family must
   fall one per chromosome before its calls may be confident.
4. **The final call is the sequence call** unless a tested line moves it under the rules below.
5. **Fixed constants, each with a reason** (a convention, a majority rule, a value measured in the data, or a named
   heuristic). Nothing is tuned to a benchmark, and there are no threshold options.

Not goals: centromeres from methylation alone (false dips are as deep as true ones); Hi-C choosing between satellite
families genome-wide (Hi-C may only choose between the sequence line's own candidates on one chromosome, its call and its
competing cluster: `hic_chose_alt`). The LTR line is developed and benchmarked with LTRquest annotations only; the readers for DANTE_LTR and
EDTA + TEsorter files stay as they are, with their tests, but get no further work.

## Stages

```
satellite -> telomere -> rdna -> rdna5s -> gaps -> attrough -> genotype -> hic -> fold -> karyotype -> provisional ->
meth -> fuse -> cpg -> plants -> morphology -> interpret -> loci -> output -> report -> plot
```

The chromosome set is chosen once per genome ([outputs.md](outputs.md#which-sequences-are-analyzed)) and every line
reads the same set. Each stage is checkpointed on its inputs and code, so a rerun resumes (Engineering notes, Checkpoints).
`genotype`, `karyotype` and `plants` handle inputs from several plants ([plants.md](plants.md)).

## Satellites (`cenfuse.sat`)

**From tandem repeats to a call.**

1. *Seeds.* FasTAN (or TRF) finds tandem repeats in 5-Mb chunks that overlap by 10 kb, so no array is lost at a cut.
2. *Families.* Seeds are clustered by sequence (BLAST) into families, and each family gets a consensus and a profile
   HMM (MAFFT alignment).
3. *Genome-wide search.* nhmmer finds every copy of every family. In a model of 500 bp or more, the positions of a
   low-complexity tract in the consensus are neutral: they emit the background, so they score nothing for any base. Such
   a tract otherwise matches every low-complexity stretch of the genome, so the model's filters pass nearly everything,
   and that was the runaway cost of the search (engineering notes). A copy usually aligns through a neutral tract; where
   it comes out in pieces they are joined, and a copy is whole when it covers 90% of the positions that can score.
4. *Roles.* Each family gets a role from where and how it occurs across the chromosomes: `centromeric`,
   `cen_associated`, `subtelomeric`, `interstitial` (knobs and other blocks), `dispersed`, `minor` or `rDNA` (rRNA
   genes, found with Rfam models).
5. *The call.* On each chromosome, the array cluster of the chosen centromeric family; where no satellite call is
   possible, a retrotransposon core, or a long-unit tandem array if another line agrees with it.
6. *Higher-order repeats* are reported for each call.

| rule | test | why |
|---|---|---|
| family founding | a seed cluster becomes a family only if one of its seeds is a tandem array (≥ 3 copies of its unit), then by size (≥ 20 kb of seeds, at most 40 families) | a satellite is tandem. Clusters founded on two-copy duplications were dispersed or absent in every genome of a 216-genome replay, and took most of the search time where they occur (T2T sheep: 93%) |
| dispersed repeat | less than half in tandem, less than half whole-unit hits, and no locus concentration (fold < 2) → never centromeric, role `dispersed` | interspersed elements (LINE/SINE fragments) annotate as chains of partial hits; satellites are tandem, and retrotransposon-derived centromeric arrays are concentrated |
| satellite kind | dominance is judged per kind: subfamilies related by sequence are pooled | splitting alpha satellite into 6 subfamilies let a pericentric satellite outcount it |
| acrocentric | a strong family at one end of most chromosomes, and at both ends of too few for independent ends (binomial, P ≤ 0.01), with no strong interior rival → centromeric, flag `acrocentric`, at most `medium` | subtelomeric repeats occupy ends independently; a centromere occurs once per chromosome. Needs about 16 or more chromosomes |
| competing cluster | a second cluster at least half the size of the main one, or a second, unrelated satellite that also marks the acrocentric end → `low`, reported in `seq_alt_*` | such calls were right only 2/10 times. Which of two satellites at one end holds the kinetochore is not in the sequence (mouse: the major satellite is larger than the CENP-A minor satellite) |
| chain of partial hits | less than half in tandem and less than half whole-unit hits, but concentrated at a locus (fold ≥ 2): the family is no tandem satellite, and a call standing on it (under 20 kb of real satellite inside) is flagged `partial_chain`: `low`, or at most `medium` where a CRM core overlaps it, or a core of another lineage that the family is not made of | one tandem copy of a retrotransposon can found a family whose hits are then fragments of that element (*Glycine dolichocarpa*: 44 whole hits of 25,732, 65% inside SIRE, calls 10 Mb wide). Only CRM may decide among retrotransposons, so a chain of Tekay or SIRE fragments must not be `high` on sequence evidence alone. It never moves a call; lines that pass their tests may raise it (moss with Hi-C and `tes/`: 3 calls, all right) |
| honesty cap | centromeric families exist but no real satellite (chains left out) is the largest cluster on most chromosomes, and no CRM core decides → every call `low` (`no_dominant_satellite`) | satellite evidence alone cannot tell a centromeric array from an interstitial one |
| genes (`genes/`) | gene-rich families are vetoed; subtelomeric = beyond the outermost genes | potato: 2 → 4/11 right; barley: the 3 confident wrong satellite calls become `low` |

## One locus per chromosome: a test of the satellite line (`sat/karyotype.py`)

A chromosome has one centromere; anything else falls roughly in proportion to length, so long chromosomes hold several and
short ones none. So take the N largest clusters of the families called centromeric, N = the number of chromosomes, and count
the chromosomes they cover. N clusters that fall in proportion to length cover far fewer than N (P from 10,000 draws; the test
passes at P ≤ 0.05). It asks nothing about any one locus, needs no known centromere, and works when the satellite differs from
chromosome to chromosome (chicken: 36 of 40 chromosomes, 17 expected).

Where the test fails, satellite evidence alone cannot say those arrays are centromeres: a satellite call is `low` unless an LTR
core agrees with it (flag `not_one_per_chromosome`). Tomato, potato, barley and *Ae. tauschii* fail, and hold most of the confident
errors the FASTA alone used to make; the holocentric genomes fail too; every genome with satellite centromeres passes but the siamang,
whose arrays split into two clusters on 8 chromosomes ([benchmark.md](benchmark.md#one-locus-per-chromosome)). Too few chromosomes give
no power (4 or fewer: `untested`), and then nothing is capped.

## LTR retrotransposons (`tes/`, `sat/signature.py`)

A lineage is a weak proxy for centromere targeting. Only part of the CRM clade is centromeric, and other lineages do it too:
Bryco in moss, Osser in *Chlorella*, and a single Tekay family in tomato, 308 of the lineage's 11,657 elements (1 in 38). So cenfuse
works at two levels, lineage and family, and learns each chromosome's signature from the other chromosomes. The input is an
LTRquest annotation, made with the [README's command](../README.md#inputs) (`--max-rounds 1 --run-sdust --reboundary --detect-args
"ltr-timeout=300"`); the tests below were benchmarked on those annotations only. A file without `family=` has no signature, and
the lineage rules below apply to its lineages alone.

- **Cores.** A set of elements has its core on a chromosome at its densest region (the full width at half maximum of its density),
  and it needs at least 3 independent insertion sites (overlapping elements count once). *Chlorella* centromeres
  hold only 2-7 Copia insertions each; 53 of the genome's 54 sit in the CENH3 domains. A core is significant at P < 0.001
  (binomial tail, Bonferroni over its possible places on the chromosome).
- **Candidates.** A lineage or a family with significant cores on most chromosomes (CRM from 3 cores on, since every CRM core
  short of a majority lay on CENH3). Only classified lineages qualify: a family without one can be a tandem array read as
  LTR pairs, which would not be independent of the satellite line. On 157 annotations with positions shuffled, nothing qualifies.
- **Specificity.** The share of a candidate's elements inside its cores divided by the share of those chromosomes the cores
  cover: a fold enrichment, comparable between a family of 50 elements and a lineage of 50,000.
- **The seed.** The most specific CRM candidate (on most chromosomes if any is), else the most specific candidate. CRM comes first
  because a non-centromeric family can be as specific (rice: a Tekay family at 18-fold, CRM at 16-fold).
- **Recruitment, leave one chromosome out.** For chromosome x, the families enriched in the current cores of the *other* chromosomes
  (binomial, Bonferroni over the families tested, P ≤ 0.001, at least 5 elements) make x's signature, and x's core is their
  densest locus. x never shapes its own signature, so a family recruited elsewhere that concentrates on x too is independent
  evidence there. The first round starts from the seed's cores, each next one from the significant cores of the round before, until
  they stop changing. Where recruitment leaves no significant core, the seed's own core stands (*Chlorella*'s 54 Osser elements
  fall in families too small to recruit).
- **What may decide.** A CRM seed with significant cores on most chromosomes decides. Another lineage decides at 10-fold or more,
  unless a dominant satellite passes the [one-per-chromosome test](#one-locus-per-chromosome-a-test-of-the-satellite-line-satkaryotypepy):
  tomato's dominant satellite fails it and its Tekay family decides, *Arabidopsis*' CEN180 passes and Athila only flags, moss has no
  dominant satellite and Bryco decides. The 10-fold level comes from 42 polyploids whose best marker is not CRM, against
  0.11's FASTA-only calls, with no truth: the share of cores within 1 Mb of the call has median 0.30 at or below 10-fold (15
  genomes) and 0.75 above (27; one-sided Mann-Whitney P = 2e-5; chance about 0.11). Counting only calls that engine did not place
  from retrotransposon clusters (843 of 934), the medians are the same (P = 4e-4). Six genomes above 10-fold agree on fewer than a
  quarter of those calls: three *Sphagnum* (Bryco, 54-80-fold; 0.04-0.07), *Lotus subbiflorus*, *Malus toringoides* and *Solanum
  stoloniferum*. Satellite calls are not truth (in tomato the satellite is wrong where the Tekay family is right), so which line
  errs there is unknown.
- **What a decisive core does.** It raises a call it agrees with to at least `medium`, takes a chromosome without a
  call, and replaces a weak call beside it or a disagreeing call that is `low` on its own (a satellite that failed the
  one-per-chromosome test counts as `low` here). A `high` or `medium` call, or a non-weak call beside the core, it replaces only
  when the core is at least 10-fold enriched on its chromosome and the signature's elements sit significantly more in the core
  than in the call (binomial at 1/2, P ≤ 0.01); otherwise the call stays, one level lower (`te_discordant`). The
  replaced array is reported in `seq_overruled_*`. Reasons: in a broad CRM plateau the densest window is noise (when the
  test was made, *B. rapa* A09 and soybean chr1/chr6 had CENH3 on the satellite: calls of 14-26 elements against cores of
  16-32; in 0.13 neither case arises any more: A09 is `high` and `te_concordant`, soybean chr1/chr6 carry no LTR flag, and
  the test acts on no benchmark call: [future.md #10](future.md#10-the-ltr-retrotransposon-line)), a satellite beside a
  core that holds hundreds of elements is pericentromeric (*Ae. tauschii*), and a core only 3-4-fold enriched says little
  (*Oryza minuta*: such cores took 9 CentO calls, and the homoeologous pairs agreeing fell from 12 of 12 to 8).
- **Weak cores.** A core of a lineage other than CRM must hold most of the signature's elements on its chromosome to act
  like this; otherwise it flags, and a call placed on it is `low` (a core is right 39 of 41 times with most of the
  elements, 3 of 11 without). A core too thin to be significant still has its place in a genome whose signature decides: on a chromosome
  without a call, or with a satellite call of a family that failed the one-per-chromosome test, it is the better candidate,
  `low` (`weak_core`).
- **The age of the elements.** The signature's youngest elements concentrate in the cores more than the oldest in nearly every CRM
  lineage and in Bryco, the reverse of lineages that accumulate near centromeres (Athila, Retand, most Tekay). It separates lineages
  but not families, so it is reported (`te_signature.tsv`) and places nothing. Since 0.13 the same contrast, measured on the other
  chromosomes, is the premise of the age tie-breaker (`ageline.py`): where it holds, younger copies at a strong validated Hi-C peak
  than at a `high` call are the second line that lets Hi-C move the call (Integration, below).
- **Formats.** An element's lineage comes from `clade=` (LTRquest, EDTA + TEsorter) or from a DANTE_LTR element's
  `Final_Classification`. DANTE_LTR's rank-D clusters (protein domains without LTRs) are not elements: in potato they
  made 195 of 223 CRM and put cores up to 15 Mb off CENH3. A file without lineages (EDTA alone gives superfamilies
  only) is reported and the line is off.

## AT troughs (`atline.py`)

Many regional centromeres sit on the most AT-rich stretch of their chromosome (*Plasmodium*, *Malassezia*, *Neurospora*, blast
fungi), whatever their sequence. That leaves a mark in the FASTA alone, but AT-rich stretches of other kinds exist in every
genome. So the one-per-chromosome test applies again: for each window of 1-500 kb (their width is not known in advance), the N deepest
independent GC troughs of the genome must cover more chromosomes than chance. The window with the smallest P is kept, corrected for
the windows tried, at P ≤ 0.05; 7 or fewer chromosomes can hardly reach that (*S. pombe*: `untested`). The call on a chromosome is
its deepest trough at that window, widened to the half-depth run around it.

The line runs only where most chromosomes have no sequence call, and never moves a satellite or LTR call. Alone it gives `low` calls
(`at_only`). With Hi-C, the prior-free peaks test the troughs directly: agreement beyond chance (Bonferroni over the windows, P ≤ 0.01)
on most chromosomes validates the line, on fewer lets the two support each other where they agree. A Hi-C call in agreement then
takes the trough's interval, which is finer than a Hi-C bin, and rises a level (`at_pair`); a chromosome Hi-C cannot place takes its
trough at `medium` if the line is validated.

## rDNA and telomeres

- **45S rDNA (NORs):** tandem runs of ≥ 3 copies of the universal 18S motif, ≤ 100 kb apart (CHM13's 5 NORs in 32 s).
  Organelle SSU rRNA holds the motif too, so a hit followed by A, whose next 19 bases also differ from nuclear 18S at
  most positions, is not counted: insertions of plastid DNA otherwise form arrays (rice chr4). Hi-C peaks beside a NOR
  never oppose a call: NORs meet in the nucleolus.
- **5S rDNA:** tandem runs of ≥ 3 whole 5S genes (the Rfam 5S profile HMM, ≥ 90% of the gene, E ≤ 1e-5), ≤ 20 kb
  apart. Whole genes, because pseudogenes and fragments are mostly partial; dispersed 5S genes form no array.
- **Telomeres:** the repeat unit is found from scratch (5-12 bp; where none qualifies, 13-29 bp, reaching the sequence
  end at most chromosome ends: *C. tropicalis* has 23 bp); an end holds a telomere array with ≥ 300 bp of unit copies in
  the terminal 100 kb at half density or more, or ≥ 100 bp that reach the sequence end (assemblers trim short telomeres:
  *Neurospora* has 157-207 bp). Once the unit is chosen, an end still without an array takes a looser stretch that reaches
  the end, ≥ 300 bp at a tenth density or more (soybean chr3: 48 kb, 42% unit copies).
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
   chance (Poisson-binomial P ≤ 0.01), and on most of the whole reference set too: with few references, leaving out one
   that disagrees can turn 2 of 4 into 2 of 3 (randomized sheep Hi-C passed twice that way; on real Hi-C the rule changed no
   gate in 14 mammals). Fewer than 3 references: `untested`. Where no call is `high` or `medium`
   anywhere, the `low` calls are the references, and only agreements at clean peaks count; where confident calls exist but
   are fewer than 4, the `low` calls join them on the same terms (0.13, `hic_gate_refs` `mixed`).
7. **Support tier.** Where agreements are beyond chance (P ≤ 0.01, at least 3) but not on most chromosomes, Hi-C
   carries centromere information, yet its peak is often elsewhere. It may then only support: a call whose own
   significant peak agrees goes up one level (`hic_support`); it never opposes, moves or places a call, but a strong
   peak on the sequence line's own competing cluster chooses that cluster (0.13; moss without `tes/`: 5 calls, all right). Only
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
  that both contradict (barley without `tes/`: 6 calls moved, all onto CENH3).

## Hi-C and the LTR line test each other (`hicfuse.py`)

Hi-C anchored on wrong calls finds nothing: potato's Hi-C finds 1 of 12 anchored on the sequence calls and 7 of 11 anchored on the true
centromeres. So where the gate fails on most chromosomes and `tes/` gives a signature, Hi-C gets two more readings, both tested against
the LTR cores, weak ones included. The sequence calls play no part.

- **The peak at the core.** Partner windows at the LTR cores. A chromosome agrees when its peak, made from the other chromosomes'
  windows only, lies in its own core's zone (the core and 4 bins each side), whose share of the chromosome is the chance. Beyond chance
  (P ≤ 0.01, at least 3) the status is `support`, on most chromosomes `pass`. A chromosome whose peak lies at its core then holds two
  independent lines at one locus (`te_hic_pair`): a call there rises a level; a `low` call elsewhere moves there (`medium`); with `pass`,
  any call does; a chromosome without a call takes the core. The windows stay where the LTR line put them: moved to the best score
  within their zones, they find each other in randomized Hi-C too (soybean: 1 of 20 chromosomes becomes 4, P 0.002). A peak counts
  only where it is clean, as for the gate's `low` references: chromosome ends, hubs and NORs attract trans contacts whatever the
  centromere. Without that, held-out holocentric *R. breviuscula*, whose Athila cores lie at the ends of 4 of 5 chromosomes, passed
  (3/5, P 5e-4) and got 3 `medium` calls; with it, 1/5 (P 0.18), and its calls equal those without `tes/`.
- **The contact crossing a bin (cis).** The two flanks of a centromere can touch each other more, or less, than DNA at the same distance
  elsewhere. So the profile is the contact between the 4 bins left and the 4 bins right of a bin, over what the chromosome's mean contact
  at those distances predicts. Its extreme lies at the centromere in many plants and fungi, not in mammals, and which extreme differs:
  a peak in tomato, rice and *C. neoformans*, a dip in potato and *Verticillium*
  ([benchmark.md](benchmark.md#the-contact-crossing-the-centromere)). So the LTR cores choose the sign: for chromosome x, the sign whose
  extreme lies at the cores of the other chromosomes most significantly, Bonferroni over the two signs. A stretch's share of the chromosome
  understates the chance here, because a profile is noisiest where coverage is low, and LTR cores lie there: randomized *Brachypodium* has
  extremes at the cores on 3 of 5 chromosomes, zones 7% of the chromosome. So the chance is the larger of that share and the rate in 200
  contact maps simulated with the same coverage and distance decay and no locus effect (55% there). Where the test passes (`support` or
  `pass`), an extreme at a call raises it a level (`cis_agrees`); within Hi-C's resolution of 4 bins it neither supports nor opposes
  (`cis_near`); farther away, a `low` call moves to it and stays `low` (`cis_only`), and a confident call stays, flagged
  (`cis_elsewhere`), one level lower at `pass`: one tested line against another, and neither more sure.
- **Which end.** A call in the terminal eighth of a chromosome cannot be supported by a peak: chromosome ends cluster in 3D whatever
  the centromere, so such a peak is flagged `telomeric`. But the two ends can be told apart: the anchored score is about twice as high
  at the end holding the call in goat (28 of 30), sheep (22/23), horse (18/19) and mouse (19/19). Where that holds on most of the
  other chromosomes with such a call, beyond the 1/2 of chance (binomial, P ≤ 0.01), a call at the higher end rises a level
  (`hic_end`), unless a competing cluster shares that end. It uses the sequence calls as anchors, and needs no `tes/`.

## The cis and trans readouts of Hi-C test each other (`hicfuse.cis_trans_line`)

Without `tes/` nothing chose the sign of the cis signal or tested it, although tomato's profile alone has its extreme on CENH3 on all
12 chromosomes. Two readouts of one library can do it: the contact crossing each bin (cis) and the contact between chromosomes
(trans). Where the gate fails on most chromosomes and there are no LTR cores:

- **Anchors.** For each sign, each chromosome's anchor is the extreme of its cis profile outside its terminal eighths (Levan's
  telocentric bound, as for the which-end test). Chromosome ends fold and cluster whatever the centromere, so there both readouts see
  the ends: read over whole chromosomes, *Ae. tauschii* passed (4/7, P 3e-5) and potato reached support (3/12) with every agreement
  in a terminal eighth and none on CENH3.
- **The test.** The anchored test of the LTR cores, on these anchors: each chromosome's trans peak, made from the other chromosomes'
  anchors only, lies at its own anchor (clean peaks; the zone's share is the chance), Bonferroni over the two signs. Both signs
  passing leaves the sign unknown (`ambiguous`): nothing acts.
- **What it may do.** Only at `pass`, and only on a chromosome whose other chromosomes pass too (leave one out, as the gate). The two
  readouts share the library's coverage, and agreement on a minority can come from coverage alone: randomized soybean Hi-C reached
  support (6/20, P 3e-5); no randomized library, mammal or holocentric genome reached a majority. For the same reason a sign at
  support never vetoes the other at pass; only both passing leaves the sign unknown. Where it acts, the cis rules apply
  (`cis_agrees`, `cis_near`, `cis_only` low, `cis_elsewhere`), except that a call inside a terminal eighth, which the line cannot
  see, is neither supported, opposed nor moved; where the extreme lies elsewhere it is flagged (`cis_blind`) and its caveat says so.
- **Result.** Tomato from the FASTA and Hi-C: 10/12 right (was 1/12); the other 2 calls lie in a terminal eighth. *C. neoformans*
  in G1, where the prior-free peaks sit at the chromosome ends: 12/14 (was 0/14), all `low`.

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
| same, call `high` | moves only with a second line at the peak (a CRM core within 4 bins; methylation `centromere_like` at the peak and a dip at the call; or, 0.13, copies of the chromosome's signature families significantly younger at the peak than at the call, in a genome whose young copies concentrate in the cores of the other chromosomes, `moved_with_age`): `high`. Otherwise `contested` (`hic_offset`), one level lower |
| Hi-C peak 1-4 bins away | neither supports nor opposes (`hic_near`) |
| competing cluster: a strong peak agrees with the second cluster, mappable DNA between them | the second cluster is chosen (`hic_chose_alt`); on the support tier too (0.13), the call rising one level (`hic_support`), unless Hi-C read on the LTR cores or the cis signal acts on that chromosome |
| competing cluster with no mappable bin between it and the call | no Hi-C support (`hic_alt_ambiguous`): Hi-C cannot tell the two apart |
| no sequence call, strong peak | the Hi-C core, `functional_only`, `medium` (a clean medium peak: `low`, `hic_medium_peak`) |
| a long-unit array that stays the call with no support from Hi-C, `tes/` or methylation | no call (`long_unit_unsupported`); the array stays in `seq_start`/`seq_end` and in the loci, and the chromosome is decided as one without a sequence call |
| gate failed, a Rabl-axis test passed, its core and the prior-free peak agree away from the call | moves to the prior-free core, `offset`, `medium` (`fold_pair`/`align_pair`) |
| same, but the axis or the prior-free peak opposes the kept call | `contested`, one level lower per opposing line |
| support tier only, the call's own significant peak agrees | `concordant`, +1 (`hic_support`) |
| no references anywhere, no sequence call, clean `high` prior-free peak | the prior-free core, `functional_only`, `low` (`hic_denovo_only`) |
| gate failed, Hi-C's peak at the LTR core, the two agreeing genome-wide | a call there +1 (`te_hic_pair`); a `low` call elsewhere moves to the core, `medium`; at `pass` any call does, and a `high` one stays `high` |
| gate failed, the cis signal tested against the LTR cores | at the call +1 (`cis_agrees`); within 4 bins nothing (`cis_near`); farther, a `low` call moves to it, still `low` (`cis_only`), and a confident call stays (`cis_elsewhere`), one level lower at `pass` |
| gate failed, no LTR cores, the cis signal tested by the trans contacts at `pass` (0.13) | the same rules, outside the terminal eighths only: a call inside one is left as it is |
| gate failed, no Rabl axis, a call in a terminal eighth at the end with the higher anchored score, as on most such chromosomes | +1 (`hic_end`) |
| no call from sequence, Hi-C or LTRs, AT troughs passed their test | the trough, `composition`: `low` (`at_only`), `medium` where Hi-C's prior-free peaks validated the troughs |
| a Hi-C call on a chromosome without a sequence call, a tested trough agreeing | takes the trough's interval, one level up (`at_pair`) |
| no call from any line, the trough and a clean prior-free Hi-C peak (placing no call of its own) agreeing, in a genome where the two agree beyond chance (0.13) | the trough, `low` (`at_denovo_pair`) |

The order of the table is the order of the rules in `fuse.py`, which since 0.14 is declared, not implied: `CASCADE` lists the rules
of which the first that applies sets the call (`_axis_pair`, `_hic_strong_peak`, `_hic_medium_peak`, `_hic_denovo_only`, `_no_call`,
`_hic_agrees_ambiguous`, `_hic_agrees`, `_hic_support`, `_hic_end`, `_hic_moved`, `_hic_offset`, `_kept`); `PHASES` then act on that
call in turn (`_te_hic_pair`, `_cis`, `_at_places`: the LTR-core and cis readings override a plain sequence call, the AT trough
places only what nothing else did); `_refine` takes an agreeing trough's finer interval and lets the Rabl-axis lines contest a kept
call; `_support`, `_against` and `_confidence` read the result. Every rule is a function of the chromosome's context (computed once)
and the call so far; a new rule is a new function in its place in the order. Each call names the rule that set its place or
relation in the `rule` column ([outputs.md](outputs.md#asmcentromerestsv)). The split was checked against the single function it
replaced on 200,000 random evidence combinations (every rule and flag reached; tests/fuse/test_fuse_equivalence.py) and on the
benchmark's 260 runs (byte-identical outputs).

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
| sat: TE_SPECIFIC_FOLD, TE_MAJOR_SHARE | 10-fold, 1/2 | agreement with independent FASTA-only calls in 42 polyploids whose marker is not CRM: median 0.30 at or below 10-fold, 0.75 above (6 exceptions above, named in [the LTR section](#ltr-retrotransposons-tes-satsignaturepy)); majority rule |
| sat: TE_SEED_MIN_CORES, TE_RECRUIT_MIN, TE_RECRUIT_ROUNDS | 3 cores, 5 elements, 5 rounds | fewer are no pattern; the results do not depend on them (right cores 152-155 of 172 across 3-10 elements and 1-4 rounds) |
| karyotype: ALPHA, N_NULL | 0.05, 10,000 draws | conventional level; P resolution 1e-4 |
| AT: BIN, SCALES, MIN_WINDOWS, APART, AT_ALPHA | 1 kb, 1-500 kb, 10 windows, 5 windows, 0.05 | Malassezia's 3-5 kb centromeres need 1 kb; known regional centromeres span the ladder; a trough needs flanks; one regional centromere gives several troughs closer than that; conventional level, Sidak over the sizes |
| Hi-C on LTR cores: zone | core ± HALF_WIN bins | the partner half-window, Hi-C's resolution |
| cis: half window, CIS_SIM | 4 bins, 200 maps | the partner half-window; the engine's null resolution (0.005) |
| which end: END_ZONE, END_SHARE, END_MIN_BINS | 1/8, 5%, 4 bins | Levan's telocentric bound; a named heuristic of about one half-window at the bins chosen; the half-window |
| cis tested by trans: the ends out, level | END_ZONE; `pass` only | the which-end test's bound; a minority agreement can come from the coverage the two readouts share |
| age: AGE_ALPHA, AGE_MIN, zones | 0.01 one-sided (the site test and the premise), 3 dated copies, ± 4 bins | the LTR line's level for overruling a call (TE_CONTRA_ALPHA); as the three insertions that make an LTR core; Hi-C's resolution at both sites |
| sat: min_array_monomers, LEN_TOL | 3 copies, 5% | a tandem array; unit-length tolerance |
| sat: TRACT_MIN_MODEL, TRACT_MIN_BP, TRACT_MAX_FRAC | 500 bp, 64 bp, 50% | what is left of a model must make a whole copy significant alone (250 positions clear the ~47 bits of E 1e-5 on a 500-Mb chromosome by far); one DUST window of low complexity in the consensus; a model mostly of low complexity stays as it is, with too little left to score |
| 45S: ORGANELLE_NEXT, NUCLEAR_AFTER | A; 19 bases, most differing | prokaryote-type SSU has A after the motif; 1,100 organelle hits differ at 11-17 of the next 19 bases, 110 nuclear hits with that base mutated at 1-5, 4 between |
| telomeres: K_LONG, MIN_DENSITY_END | 13-29 bp; 0.1 | up to the satellite engine's shortest unit; a heuristic floor that keeps out copies scattered every few kb (61 of 6,336 ends reach the end under half density, none below 0.19) |
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
| a cis profile of the centromere | observed / expected contact across each bin, its sign and its level tested against the LTR cores, chance from simulated maps | the plain mean contact with the zone's share as chance: randomized Hi-C passed (sorghum 5/10, P 3.5e-06), since beside a masked bin the box holds another mix of distances and noise is largest at low coverage. The same profile tested against one locus per chromosome, with no LTR cores: too little power (tomato passes, potato and *Verticillium* miss) and chicken's dips pass while 2/26 lie on CENP-A. Tested against Hi-C's prior-free peaks: no agreement in tomato, potato or rice, and *Ae. tauschii* agrees on 6 of 7 with 0 of 7 right |
| Hi-C partner windows on LTR cores | fixed at the cores | moved to the best score within reach: more real hits (potato 3 to 6 of 12), but randomized Hi-C passed too |
| an end-distance expected for trans contacts | not used | moss 23 → 25 of 26, but sheep 9 → 5, chicken 7 → 6, goat 12 → 11 (v10.1) |
| Hi-C depth fragility | reported as depth guidance | a split-half check: unstable halves did not predict each other, and Hi-C time doubles |
| `.hic` files | stop with the `hictk convert` recipe | reading `.hic` natively: a new dependency for what one command does |
| a CRM core vs a confident distant or adjacent satellite call | replace only if CRM marks the core significantly more | always replace: wrong on *B. rapa* A09 and soybean chr1/chr6. A bootstrap of the core position: lost *P. annua* Pa1A |
| TE core floor | 3 independent insertion sites | 5 elements: missed *Chlorella*'s 2-7-insertion centromeres |
| LTR elements without a lineage (fungi) | not used | pooling them, three forms (v6): the first lifted the honesty cap in fish and *Oryza* genomes; the guarded form placed *Cryptococcus* JEC21 10/14 and a hybrid's D subgenome 10/13, all `low`, none wrong, and never lifted the cap, but still qualified in 9 fish genomes without truth and in *A. suecica*, and its 20-kb bin floor cost *O. punctata* a CRM call ([future.md #47](future.md#47-centromeric-retroelements-the-ltr-line-cannot-read)) |
| several qualifying TE lineages | the most concentrated | chromosomes with a core first: agreed less with independent satellite calls |
| CpG dip inside a call | most prominent valley, half-depth run | deepest bin: lands on array edges. Longest run below half the median: on the CDR 52/90 (valley 88/90) |
| CpG null | positional (rotation) test for the gate | dips in other satellite clusters: passed on misplaced methylation in gorilla |
| CpG choosing between arrays | not done | least-methylated window: on the CDR 3/23; most prominent valley across clusters: wrong cluster on 8-10/25 |
| the centromeric region (call + nearby satellite clusters) | not reported | held out: precision 0.07 in apes (an interspersed family chains whole chromosomes) |
| chromosome set | tenfold drop, tapering tail cut, optional list | the drop alone kept 5,614 extra scaffolds over 121 polyploid assemblies; other cuts dropped dot chromosomes or microchromosomes |
| a sample screen to skip dispersed families | not used | no family was skipped on real genomes |
| long models with low-complexity tracts in the search | neutral positions in the model (background emissions at dustmasker's tracts of the consensus); pieces of a copy joined afterwards | masking the target and searching again around hits near a masked stretch: a copy split by the mask had to be searched twice, so hit-dense models cost more (0.45-1.42x overall, measured side by side). Turning off HMMER's composition correction for these models: it repairs some copies, but with 45% more significant hits in potato, none seen by the plain model. Tracts that keep 10-30% of their signal: they repair copies, but on sheep cost more than the plain model (26-31 vs 21 s). A BLASTN prefilter: lost 12% of the hit bp in chicken and sheep. Filtering harder in nhmmer (F1-F3): under 2x, and hits changed. Other block lengths: no gain |
| 45S arrays made of organelle DNA | the base after the motif and the next 19 | one strand only: large real NORs hold both (*Coffea*: 658 units, 336 and 322). The next base alone: 9% of the hits followed by A are nuclear copies with that base mutated (barley 7H, 10 5.8S genes in 120 kb) |
| which families count as one for dominance | same sequence, or period within 15% of the shorter | sequence alone (`same_sequence`): changes the main family on 90 of 3,459 chromosomes and removes dominance in 6 genomes, with no truth behind it |
| dominance within subgenomes (chromosome groups) | not used | real A/B groups lift 27 of 292 capped calls in 20 polyploids, random splits of the same chromosomes 9, and only 1 species beats its random splits; random splits of diploids lift up to 8 calls per genome |
| telomere stretches under half density in choosing the unit | only after the unit is chosen | decoy units in loose stretches could block the real unit |
| a faster, inexact seed clustering | not used | 10× faster, but BLAST's hit lists differ for 0.2-0.5% of pairs, so it is not exact by construction ([future.md](future.md)) |
| GFF3 feature types | Sequence Ontology terms, cenfuse's kind as `kind=` | cenfuse's own kinds as types: no validator or pipeline accepts them |
| the cis signal without `tes/` | tested by the library's own trans contacts, terminal eighths out, acting at `pass` | the whole chromosome read: *Ae. tauschii* passed with 0 of 4 agreements on CENH3, all in terminal eighths; acting at support: randomized soybean Hi-C reached it; a cis + trans pair as a `medium` call: 2 near misses in *C. neoformans* G1 would be confident errors |
| AT troughs of one centromere counted once | not changed | excluding each trough's half-depth run: no genome changes (*C. neoformans*' centromeres hold two sub-troughs separated by GC-normal sequence) |
| Hi-C on the support tier and a far strong peak | it chooses the competing cluster only | moving any `low` call to the peak: moss 7 right, 1 wrong, and nothing elsewhere |
| insertion age as a locator | a tie-breaker only (`moved_with_age`) | the core of the young copies: on CENH3 for 128 of 172 chromosomes, the density core 155 |

## Known risks

- The LTR line's `medium` calls outside CRM genomes. In 13 held-out tomatoes, whose centromeres a Tekay signature marks, 15 of the 69 confident calls with `tes/` are wrong
  (22%; benchmark.md, held-out tomatoes), and a call that merely agrees with the core is the least reliable kind (6 wrong of 17). Treat such calls as a location to check.
- The Hi-C test takes the sequence calls as references. Where those are wrong genome-wide, only the Rabl axis,
  `tes/` or methylation can rescue Hi-C (potato and barley Hi-C anchored at the true centromeres find 7/11 and 6/7 of
  them; anchored at the calls, none).
- Hi-C alone trusts the strongest trans cluster to be the centromeres. In a genome without sequence calls whose
  dominant cluster is something else, its `low` calls are wrong; a holocentric genome without sequence calls would get
  them too.
- One wrong Hi-C peak passes in a genome where Hi-C mostly works; only `high` calls need a second line to move.
- Hi-C tests lose power gradually but flip at their threshold (both orangutans at half depth); a failed test reverts
  to the sequence line.
- LTR families are learned from the genome's own cores, so a family that concentrates next to centromeres can be learned as if centromeric:
  tomato's cores miss CENH3 on 5 of 12 chromosomes, and without Hi-C two of those misses stay `medium`. *Arabidopsis*' Athila (13.3-fold)
  flanks CENH3 by 0.46-1.38 Mb on all 5 chromosomes; it only flags because CEN180 passes the one-per-chromosome test. In a genome whose
  satellite fails that test, such a lineage would decide and place calls beside the centromere. In 6 of 27 polyploids above 10-fold
  the cores and the FASTA-only calls mostly disagree (no truth).
- A call can agree with an LTR core it only touches: where the core is much wider than the centromere, a call at one edge is
  `te_concordant` while CENH3 sits at the other, and the agreement lifts the caps that would keep it `low`. 5 of the 7 confident
  errors of the `tes/` arms are such calls, and so is one of the two held-out genomes' errors (*Ae. geniculata* chr05U, in a 9.3-Mb
  core). Agreement at the core's maximum was measured as a remedy and lowers right calls too
  ([future.md #45](future.md#45-confident-calls-at-the-edge-of-a-wide-ltr-core)).
- A thin core of a decisive signature takes a chromosome whose satellite failed its test (`weak_core`). On the held-out
  *S. verrucosum* that replaced 2 right satellite calls (cores of 2 and 4 elements; both stay `low`), and a decisive core that
  replaced a weak array beside it was confidently wrong (chr06); the development genomes had shown neither
  ([future.md #36](future.md#36-no-held-out-test-of-the-013-rules)).
- The LTR signature marks where centromeric families accumulated. A centromere that has moved away from its satellite and those families
  is outvoted: the two sequence lines agree on the old site and Hi-C can only contest it, unless a second line at the Hi-C peak
  lets Hi-C move the call (*P. annua* Pa1A: earlier analyses put its centromere at the Hi-C peak, 21 Mb from the satellite; 0.12
  kept the satellite site at `medium`, `contested`; since 0.13 insertion age is that second line and moves it 20.7 Mb to the
  Hi-C site, `high`, without CENH3 truth: [future.md #31](future.md#31-a-centromere-that-left-its-satellite-and-its-ltr-families)).
- Without LTR cores the cis line is tested by the trans contacts of the same library: coverage they share can make them agree on a
  minority of chromosomes (randomized soybean: support), so it acts only at `pass`; a structure that both folds and clusters away
  from the centromere on most chromosomes would pass (none seen outside the terminal eighths). It cannot see the terminal eighths,
  so a wrong call there stays (tomato: 2 of 12).
- The age tie-breaker trusts that centromeric families still insert where the centromere is now. It tests that premise on the
  other chromosomes (their younger copies must lie in their cores more often than their older ones), so lineages that accumulate
  beside centromeres with age (Athila, Retand, most Tekay) cannot use it; within a genome that passes, recent dispersed insertions
  at a non-centromeric Hi-C peak could still tip it. It only acts where validated Hi-C already contests a `high` call (in the
  benchmark: *P. annua* Pa1A, without truth; its premise P 5e-10).
- The AT test passes on 4% of shuffled genomes (P ≤ 0.05 by design), and then gives `low` calls on chromosomes nothing else places. Point
  centromeres and AT-rich genomes without a clear trough (*Z. tritici*, *C. tropicalis*) do not pass.
- The one-per-chromosome test has no power with few chromosomes (`untested`) and passes on any family that happens to lie once per chromosome.
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
- **Search cost.** On the slowest genomes (*Panax*, *Lespedeza*) models of 1 kb or more took 87-96% of the CPU. A long
  model with a low-complexity tract spent 25-50 s per 2-Mb chunk with no hit: simple sequence passes its filters, and each
  window then costs model x window DP. A neutral position cannot trigger that. What neutral tracts cannot remove is the
  DP on real hits: long chain families hit nearly everywhere. A neutral position also earns nothing for a copy, which
  costs some alignments three ways, each seen on real families at 97-100% identity. A local alignment stops at the inner
  side of a tract at a model end (Col-CC's 502-bp rDNA family: a 98-bp tract at its start). nhmmer's correction for biased
  composition (null2) still charges a hit for the biased residues it holds, so a hit can stop at a tract with little
  beyond it (*Chlamydomonas*' 924-bp family: 59 bp, a 279-bp tract, 586 bp). Where a model holds an internal duplicate, a
  jump to it can outscore a walk over the tract (potato's 1,362-bp family: 9-514 and 608-1360). So a hit over 90% of the
  positions that can score counts as whole, and colinear pieces with a tract between them are joined if they lie as far
  apart on the target as on the model except for the tract (whose length varies between copies, so it may differ by its
  own length plus 20 bp) and miss no more scoring positions than a whole copy may. Such a hit reports the whole model in
  `hmm_from`, `hmm_to` and the GFF's `model_cov`, while its coordinates stay those of the alignment. E-value parameters stay
  those of the whole model. A known side effect: a neutral model scores less, so a close family without a tract can win their
  overlaps and take its copies (cotton's 512-bp family, to a 481-bp relative); no call has moved for it.
- **Parallel family models.** A family's model is built in the frame of the most similar earlier family, so models
  are built concurrently but accepted in family order, each rebuilt if an earlier family it did not see could have
  changed it: the result is the serial one (*Panax*: families stage 3.4× faster at 16 threads).
- **FasTAN crashes** are recovered by splitting the piece ([outputs.md](outputs.md#running)).
- **Resume guard.** The satellite checkpoints stand on the chromosome set, which is recorded before the engine runs; a
  changed set stops the run instead of reusing seeds made on other sequences.
- **Checkpoints (0.14).** One rule for what a checkpoint stands on: the code, the parameters and the content of the inputs.
  cenfuse's stages key their checkpoints on the syntax trees of the modules they run and of everything those import inside
  cenfuse, of their own `Run.stage_*` method, and of the package data they read (the Rfam model of the 5S scan), with comments
  and docstrings stripped (editing one recomputes nothing, a changed import is seen); on cenfuse's version; on the whole
  FASTA's md5 for every stage that reads the sequence (a same-length edit such as hard-masking keeps the `.fai`); and on
  content fingerprints of their other inputs (whole files up to 32 MiB, else the size and four 1-MiB samples; no names or
  time stamps, so a moved or copied folder resumes). The helpers that keep a run (fingerprints, atomic writes, locks, the run's
  cost, worker pools) live in `runtime.py`, and the manifest in `manifest.py`: both are kept out of every key (`pipeline.KEYLESS`),
  although every engine imports them, since they change no result. The satellite engine keys each stage on its results version,
  the parameters of it and the earlier stages, the FASTA's md5 (a stamp made before 0.14 lacks it: it is trusted only while the
  FASTA is older than it, then rewritten with the full key),
  and from classify on the content of `genes/` and `tes/`: adding `tes/` to a finished run, the commonest reason a run resumed
  on stale satellite calls, now recomputes classify onward with a warning. A stale search stage still stops the run (it may
  cost hours). Every output and checkpoint is written to a temporary file beside its target and moved over it, and a stage
  that recomputes first unlinks the files of its folder that have other hard links (its tools, BLAST or FasTAN, write in place):
  so a reader or a resume never sees a half-written file, and a folder copied with hard links (`cp -al`, `rsync --link-dest`) and
  recomputed leaves the original as it was (found by the 0.14 verification replay, whose held-out arms were such copies). An
  unreadable checkpoint is recomputed with a warning. `--force STAGE` recomputes that stage, `--force all` every stage; a
  stage that reads another's result has it in its key, so it follows. Tested: tests/fuse/test_checkpoints.py.
- **Locks (0.14).** A run locks its output folder (`.cenfuse.lock`, created exclusively, holding host, pid, start time and
  command, and touched every minute while the run holds it): a second run on the folder stops. A stale lock is taken over
  with a warning: on this host when its pid is gone or now belongs to another process (its kernel start time differs from the one
  the lock records: a reused pid); from
  another host, or when the lock cannot be read (a run killed between creating and writing it), after 15 minutes without a
  touch (a heuristic: 15 missed heartbeats). The takeover runs under a second, exclusively created guard file, so of two runs
  that find one stale lock only one replaces it and neither removes the other's live lock. SIGTERM (a scheduler's kill) ends
  a run as an exception, so the lock is released and the manifest says it failed. The shared seed cache uses the same rules
  per genome, so a lock left by a job killed on another node no longer makes later runs wait forever.
- **The run manifest (0.14).** `run_info.json` records cenfuse's version and git commit, the command, host, threads and times,
  Python and package versions, every external tool's path, version line and conda package (FasTAN and TRF print no
  version), the inputs' fingerprints, and how the run ended. cenfuse and its engines share one version, so the headers of
  the satellite, Hi-C and methylation outputs name it (the Hi-C contact-matrix cache keeps its own format version).
  `--satellite-dir` warns when the imported outputs' manifest names another version. A failed report section or figure is
  recorded (summary column `incomplete`, a placeholder in the report, a final warning) instead of one log line. The manifest
  also records what the run cost (`cost`: wall time, the CPU hours of cenfuse and of every process it started, and the peak
  memory of the whole process tree, sampled every 30 s as the sum of proportional set sizes: GNU time's maximum counts only the
  largest single process, not the seed pool's workers or the tools). `cost.memory` says what was measured: a kernel older than
  Linux 4.14 gives only resident sizes, which count shared pages in every worker (an overestimate).
- **Worker processes (0.14).** Under fork a worker inherits the log handlers and writes its own lines, as before. Under spawn
  (macOS) or forkserver (Python 3.14's default on Linux) it has none, and FasTAN's crash lines were lost: there the process
  pools (seeds, methylation bins, genotypes) send their workers' records over a queue to the main process, which writes them,
  and the queue's listener is stopped within 10 s whatever a worker did (one killed while it held the queue's lock would
  otherwise hang the run: found by the 0.14 code review). A write's temporary file left by a killed run is removed by the next
  write of that file.
- **Hi-C speed (0.14).** Profiled on the slowest arm (siamang, randomized Hi-C, 999 s of Hi-C), 95% of the stage was two Poisson
  bootstraps of 100 replicates: the engine's, on its own partner windows, and cenfuse's, on the windows of the sequence calls.
  They drew the same resamples (one seed, one matrix), and a resample, its balancing and its block means are about 97% of a
  replicate's cost; the windows differ only in the last, cheap step. So the engine's bootstrap scores cenfuse's windows too, on
  the same resamples (`hic.confidence.bootstrap`, `extra`; cenfuse passes its windows to the engine's call), and each resample
  is written straight into the matrix's symmetric layout (`_Resampler`: the same draws, the same arrays, without a sparse
  transpose and two sums per replicate). The null replicates (their random partner positions drawn first, in the former serial
  order) and the cis simulations (one random stream per chromosome) run on the thread pool too, though they cost little.
  Every number is the same, at any thread count.

## Tests

`pytest -q -n 8`: about 1,150 tests on small synthetic fixtures (5 min at 6 workers on a shared server).
`CENFUSE_REQUIRE_TOOLS=1` turns every skip for a missing external tool into a failure (tests/conftest.py), for CI and
releases. They cover every rule above with unit tests; a
synthetic genome with planted cases (a Hi-C move with a CRM core, the same without a second line, a hub decoy, a
satellite-free chromosome, randomized Hi-C, three plants); a full satellite-engine run on its synthetic genome, with and without
`tes/`; the LTR signature (seed, leave-one-out recruitment, shuffled positions); the one-per-chromosome test; GC troughs planted one
per chromosome against a decoy; Hi-C anchored on LTR cores, and the cis signal on maps with low coverage, masked bins and no locus
effect, which must fail; the cis signal tested by trans contacts (pass, support that does not act, both signs, extremes at the ends,
leave one out) and the integration rules of 0.13 (calls in a terminal eighth, the support tier's competing cluster, `mixed` references,
an AT trough with a prior-free peak, insertion age and its premise); planted CpG dips and every methylation format; every Hi-C pairs layout; a yeast-like Rabl genome called at
10-kb bins; the chromosome-set rule, the list file and the resume guard; FasTAN crash recovery; parallel family models byte for byte
equal to the serial build; the report's sentences.
