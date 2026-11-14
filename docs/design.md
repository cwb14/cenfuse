# How cenfuse works

cenfuse calls the centromere of each chromosome from several kinds of evidence, called *lines*: satellite arrays, LTR
retrotransposons, AT-rich sequence, repeat density, repeat families, Hi-C contacts and DNA methylation. Each line calls on its
own first. A line may move or raise a call only after it passes a test in the genome at hand; an untested line may only place a
`low` call where nothing else does (Hi-C alone; an LTR core that only flags). Every agreement and disagreement is reported. For
each line, this page says what it reads, how it calls, the test it must pass and what it may change; then how the lines are
combined. Numbers come from the benchmark ([benchmark.md](benchmark.md)); every benchmark genome is development data. Release
history: [CHANGELOG.md](../CHANGELOG.md); open problems: [future.md](future.md).

## Principles

1. **Each line calls on its own first.**
2. **Lines check each other.** No line is trusted on its own confidence: the Hi-C peaks that the Hi-C engine rates `high`
   lie within 4 bins of CENP-A on only 14 of 63 T2T mammal chromosomes. A chromosome's own agreement never counts toward
   the test that lets a line act on it: every test that validates one line by its agreement with another (Hi-C against the
   calls or the LTR cores, the cis signal, the arm fold and trans alignment, troughs and repeat-density peaks against Hi-C) acts on
   a chromosome only at the weaker of two verdicts, the genome's and the one on the other chromosomes alone (`loo.py`). A
   genome of three chromosomes therefore validates no such test (*S. pombe*), as the Hi-C gate never could there.
3. **Whether a line works is measured in each genome.** Most tests rest on one fact: a chromosome has one centromere. What
   marks centromeres falls once per chromosome; anything else falls roughly in proportion to chromosome length.
4. **The final call is the sequence call** unless a tested line moves it ([Integration](#integration-fusepy)).
5. **Fixed constants, each with a reason** (a convention, a majority rule, a value measured in the data, or a named
   heuristic). Nothing is tuned to a benchmark; there are no threshold options.

Not goals: centromeres from methylation alone (false dips are as deep as true ones); Hi-C choosing among satellite families
genome-wide (it only chooses between the sequence line's two candidates on one chromosome). The LTR line is developed on
LTRquest annotations; the DANTE_LTR and EDTA + TEsorter readers are kept, with their tests, but not developed further.

| term | meaning |
|---|---|
| sequence call | the call from the FASTA (satellites, with LTR cores when `tes/` is given); every other line starts from it |
| `high`, `medium`, `low` | confidence. *Confident* = `high` or `medium`; `low` = cenfuse is unsure |
| `pass`, `support`, `fail`, `untested` | a test's status. Agreement tests (Hi-C, and what Hi-C tests): agreement beyond chance on most chromosomes; on a minority only; neither; too few references. One-locus tests (satellites, AT troughs, repeats): P at or below the level; no `support`; above it; out of reach (even one locus per chromosome could not pass, or too many sequences) |
| most | more than half |
| bin | a Hi-C bin. Hi-C's resolution is 4 bins (the Hi-C engine's partner half-window) |
| terminal eighth | the outer 1/8 of a chromosome at either end (Levan's telocentric bound) |
| competing cluster | a second cluster of the centromeric families' arrays, at least half the called one, or the terminal block (≥ 20 kb, in the terminal 5%) of an acrocentric satellite unrelated by sequence |

## Stages

```
satellite -> telomere -> rdna -> rdna5s -> gaps -> attrough -> repeats -> repfam -> genotype -> hic -> karyotype -> fold ->
provisional -> meth -> cpg -> fuse -> plants -> morphology -> interpret -> output -> loci -> report -> plot
```

The satellite engine keeps its own checkpoints; `telomere` to `cpg` are checkpointed per genome (the `--force` names); the
other steps are cheap and always rerun. Some run only when needed: `attrough`, `repeats` and `repfam` where most chromosomes
lack a (confident) sequence call; `genotype` with CHG methylation; `karyotype` with two or more methylation samples or Hi-C
libraries; `fold` with a deciding Hi-C library. The LTR line acts inside the satellite stage (classify: the signature, its
verdict, rejected families; centromeres: cores replace or flag calls); `teline.py` reads `tes/` again for the integration,
loci and figures, and insertion ages are read before the provisional call. The CpG gate reads the provisional calls, its
contest acts in `fuse`, and its kinetochores are read on the final calls. `genotype`, `karyotype` (missing chromosome
segments, not the one-per-chromosome test) and `plants` handle inputs from several plants ([plants.md](plants.md)). The
chromosome set is chosen once per genome ([outputs.md](outputs.md#which-sequences-are-analyzed)), and every line reads the
same set.

## Satellites (`cenfuse.sat`)

**Reads** the FASTA, and the gene annotation (`genes/`) and the LTR annotation (`tes/`) if given.

**Calls.**
1. *Seeds.* FasTAN (or TRF) finds tandem arrays of 30-2,000-bp units, at least 2 kb long, in 5-Mb chunks that overlap by
   10 kb. Satellites of units under 30 bp are missed; longer units go to a long-unit detector.
2. *Families.* BLAST clusters the seeds into families; each gets a consensus and a profile HMM (MAFFT).
3. *Search.* nhmmer finds every copy. In models of 500 bp or more, low-complexity tracts score nothing, so they do not
   match every simple stretch of the genome; a copy is whole at 90% of the positions that can score.
4. *Roles,* from where and how each family occurs across the chromosomes: `centromeric`, `cen_associated`, `subtelomeric`,
   `interstitial` (knobs and other blocks), `dispersed`, `minor` or `rDNA` (rRNA genes, found with Rfam models).
5. *The call:* the cluster of centromeric-family arrays (all such families pooled, arrays < 500 kb apart) that holds the main
   cluster of the family chosen on that chromosome, else the largest (where only the largest agrees with a decisive LTR core,
   the largest). A cluster of 5-20 kb gives a `weak_satellite` call, `low` on its own (`medium` where a decisive LTR core
   agrees). With no cluster: an LTR core, or a long-unit tandem array that stays only if another line supports it.
6. *Confidence:* `high`, one level lower for each of `low_satellite` (under 50 kb of satellite), `near_end` (reaching into
   the terminal 5%) and `diffuse` (a call of 5 Mb or more, under 20% satellite), two lower for a competing cluster.
   Higher-order repeats are reported for each call.

| rule | test | why |
|---|---|---|
| founding a family | a seed cluster founds a family only if a seed has ≥ 3 copies of its unit; then by size (≥ 20 kb of seeds, at most 40 families) | a satellite is tandem; clusters founded on two-copy duplications were dispersed or absent in a 216-genome replay |
| dispersed repeat | under half of its bp in tandem, under half of its hit bp in whole units, main cluster < 2× its mean density → `dispersed`, never centromeric | interspersed elements annotate as chains of partial hits; retrotransposon-derived centromeric arrays are concentrated |
| dominance | a *strong* family (≥ 20 kb on at least half of the chromosomes, median main cluster ≥ 100 kb) whose kind is largest on at least half of them. A kind pools the families of one sequence; a chromosome counts when its largest kind holds the family or one related to it (same sequence, or period within 15%) | splitting alpha satellite into 6 subfamilies let a pericentric satellite outcount it |
| acrocentric | a strong family that would be subtelomeric, at one end of most chromosomes, at both ends of too few for independent ends (binomial, P ≤ 0.01; out of reach below about 17 chromosomes), no strong interior rival → eligible as centromeric, flag `acrocentric` | subtelomeric repeats occupy ends independently; a centromere occurs once per chromosome |
| competing cluster | see the terms above → two levels lower; reported in `seq_alt_*` | which array holds the kinetochore is not in the sequence (mouse: the major satellite outweighs the CENP-A minor satellite) |
| chain of partial hits | a concentrated family under half tandem and under half whole is no satellite: a call with under 20 kb of real satellite on it is `partial_chain`, `low`, at most `medium` where a core of a decisive signature (or of a lineage the family is not made of) overlaps it; a narrower strong core inside narrows it (`chain_to_core`) | one tandem copy of a retrotransposon can found a family of its fragments (*Glycine dolichocarpa*) |
| honesty cap | centromeric families exist, none dominates, no decisive CRM signature → satellite calls `low` (`no_dominant_satellite`), except those an LTR core agrees with | satellite evidence alone cannot tell a centromeric array from an interstitial one |
| genes | families whose main clusters hold ≥ half the genome's gene density are vetoed; subtelomeric = beyond the outermost genes; a call denser in genes than its chromosome is flagged `gene_rich` | potato from the FASTA: 2/11 right, with genes 4/11 |

## One locus per chromosome: a test of the satellite line (`sat/karyotype.py`)

Take the N largest clusters of the families called centromeric (N = the number of chromosomes) and count the chromosomes they
cover. Clusters that fall in proportion to length cover far fewer. P comes from 10,000 random draws; the test passes at
P ≤ 0.05. It needs no known centromere and works when the satellite differs between chromosomes (chicken: 36 of 40 covered,
17 by chance).

- **Fail:** a satellite call is `low` unless an LTR core agrees with it (`not_one_per_chromosome`). Tomato, potato, barley,
  *Ae. tauschii* and the holocentric genomes fail; every genome with satellite centromeres passes except the siamang, whose
  arrays split into two clusters on 8 chromosomes.
- **`untested`** where even covering all N could not reach P ≤ 0.05 (in practice 4 or fewer chromosomes), or above 5,000
  sequences; `none` without a centromeric family. Nothing is capped then.

## LTR retrotransposons (`tes/`, `sat/signature.py`)

**Reads** an LTRquest annotation made with the [recipe](inputs.md#ltr-retrotransposons-tes). Lineages come from `clade=` (LTRquest,
EDTA + TEsorter) or DANTE_LTR's `Final_Classification` (its rank-D domain clusters are not elements). A file without lineages
turns the line off. Only part of the CRM clade is centromeric, and other lineages do it too (Bryco in moss, Osser in
*Chlorella*, one Tekay family in tomato: 308 of the lineage's 11,657 elements), so cenfuse learns, at family level, each
chromosome's signature from the other chromosomes.

**Calls.**
- *Core:* the densest region of a set of elements on a chromosome (full width at half maximum), with ≥ 3 independent
  insertion sites (overlapping elements count once); significant at P < 0.001 (binomial, Bonferroni over its possible
  places). A call *agrees* with a core when the gap between them is at most the width of the wider one.
- *Seed:* among classified lineages and families with significant cores on most chromosomes (CRM: from 3 cores), the most
  specific CRM candidate with cores on most chromosomes, else the most specific CRM candidate, else the most specific
  candidate (where that one is under 10-fold, a pooled seed first). Specificity is the fold enrichment of elements in the
  cores. CRM comes first because a non-centromeric family, or a small CRM family on a few chromosomes, can be more specific.
  With element positions shuffled within chromosomes (270 shuffles of 54 genomes), nothing seeds.
- *Recruitment, one chromosome left out:* x's signature is the families enriched in the *other* chromosomes' cores (binomial,
  Bonferroni over families, P ≤ 0.001, ≥ 5 elements); x's core is their densest locus. Rounds repeat until the cores stop
  changing (at most 5). Where recruitment leaves no significant core, the seed's own core stands.
- *Pooled seed* (`signature.pooled_candidates`), for families with few copies per chromosome. Where no CRM candidate exists
  and none reaches 10-fold, a family seeds if its core has P ≤ 0.05 on more than half of the chromosomes, Fisher's combination
  of the core P values of every chromosome holding a core is P ≤ 0.001 after Bonferroni over the file's families, ≥ 3 cores
  are significant, and it is ≥ 10-fold; of several, the one with the most significant cores.
- *Lineages only* where under half of the elements name a family (`family=`) or nothing seeds: CRM decides with significant
  cores on most chromosomes, other lineages only flag; CRM's significant cores, where fewer, stand in for the marker's on their
  chromosomes, and overruling and rejecting need no fold.

**What may decide.** A *decisive* signature may change calls. Any other places only `low` calls (where there is no call, or
only a long-unit array away from its core) and otherwise flags; a call that agrees with any core escapes both caps. CRM
decides with significant cores on most chromosomes. Another lineage decides at ≥ 10-fold, with significant cores on most
chromosomes (or a pooled seed) and a median core fold ≥ 10 on their own chromosomes, unless the satellite line alone
(classified without the LTR line) has a dominant satellite whose centromeric families pass the
[one-per-chromosome test](#one-locus-per-chromosome-a-test-of-the-satellite-line-satkaryotypepy) or are `untested`. Tomato's
satellite fails and its Tekay family decides; *Arabidopsis*' CEN180 passes and Athila only flags.

**Satellites and TEs: what a decisive core does.**
- *Rejects a satellite kind* (the families of one sequence) that agrees with no strong ≥ 10-fold core on most of the
  chromosomes it is chosen on that hold one; the families are chosen again until none disagrees (`seq_rejected_family`, where
  a rejected cluster of ≥ 20 kb lies outside the call).
- *Agrees:* a call that agrees with the core is at least `medium` (`te_concordant`; a chain of partial hits beside the core
  stays `low`), and the two caps above do not apply to it. A chromosome without a call takes the core, `medium`.
- *Replaces* a weak call beside it (`weak_satellite`, long-unit, or under 50 kb of satellite: `weak_satellite_beside_core`),
  and a satellite call away from it that is `low` on its own. A `high` or `medium` call away from it, or a satellite array
  beside it, it replaces only when the core is ≥ 10-fold on its chromosome and the signature's elements sit significantly
  more in the core than in the call (binomial at 1/2, P ≤ 0.01; `satellite_away_from_core`, `satellite_beside_core`). A
  confident call away from it that survives stays, one level lower (`te_discordant`); an array beside it that survives stays
  `te_concordant`. The replaced array is in `seq_overruled_*`. In a broad CRM plateau the densest window is noise; a
  satellite beside a core of hundreds of elements is pericentromeric (*Ae. tauschii*).
- *Weak cores.* A non-CRM core acts so only when it holds at least half of the signature's elements on its chromosome. A
  significant core holding fewer replaces only a satellite call away from it that is `low` on its own, and stays `low`
  (`minority_core`; 13 tomatoes: right on 41 of 58); it places a `low` call where there is none; a call that agrees with it
  keeps its level (`te_concordant`, not capped). A core too thin to be significant takes a chromosome without a call, or
  whose satellite failed the one-per-chromosome test: `low` (`weak_core`).

**The age of the elements** (LTRquest's insertion times, `K2P_T=`). Young copies concentrate in the cores more than old ones
in nearly every CRM lineage and in Bryco; Athila, Retand and most Tekay show the reverse (reported in `te_signature.tsv`). Age
is used twice:
- *Choosing a core* (`signature.age_choice`): among comparable density peaks of the signature's elements (each ≥ half the
  highest, with 3 independent insertions, and itself a significant core), one significantly younger (one-sided Mann-Whitney
  P ≤ 0.01, Bonferroni over the peaks, ≥ 3 dated copies each) is the core (`te_age_core`). Copies insert into the active
  kinetochore, so the older cluster is where the centromere was (Wolfgruber et al. 2009). `core_age_from` gives the interval
  age replaced; `core_alt`, on every core, the highest other peak at least half as high, with 3 independent insertions.
- *A second line for Hi-C* (`ageline.py`): where young copies lie in the cores of the other chromosomes more often than old
  ones (the premise, one-sided Fisher), younger copies at a strong validated Hi-C peak than at a `high` call (one-sided
  Mann-Whitney), both at P ≤ 0.01, let Hi-C move it (`moved_with_age`, or `moved_with_te+age` with an LTR core too).

## AT troughs (`atline.py`)

**Reads** the FASTA (GC per 1-kb bin). Many regional centromeres sit on the most AT-rich stretch of their chromosome
(*Plasmodium*, *Malassezia*, *Neurospora*), whatever their sequence; AT-rich stretches of other kinds exist in every genome.

**The test.** For each window of 1-500 kb, the genome's independent troughs are ranked by depth, and the deepest class is read
off the depths alone (the largest step, the class holding N to 1.5 N troughs, N = the chromosomes long enough for the window).
One centromere per chromosome spreads that class over the chromosomes more evenly than random deals that place each trough in
proportion to the room a chromosome has left (Freeman-Tukey statistic, lower tail, mid-P of 10,000 deals). The window with the
smallest P is kept, corrected for the windows tried (Sidak); pass at P ≤ 0.05. `untested` where even one trough per
chromosome could not pass after the correction (about 6 or fewer chromosomes of similar length), where no window size has 10
windows on at least 3 chromosomes and on more than half of them, or above 200 sequences. GC blocks shuffled among
chromosomes pass 4% of the time; holocentric genomes fail.

**Calls.** The deepest trough at the chosen window, widened to its half-depth run (the window Hi-C's test chose where that
test passes, or supports while the FASTA test does not).

**What it may change.** It runs where most chromosomes have no sequence call (a long-unit array counts as none) or a `low`
one with a competing cluster. It never moves a confident satellite or LTR call. Where the troughs passed their test or Hi-C's
prior-free peaks validated them:
- A chromosome nothing else places takes its trough (`at_only`: `low`; `medium` where validated).
- On a `low` call no other rule changed, with a competing cluster, the trough's deepest window on that cluster and not on the
  call chooses the cluster (`at_chose_alt`), still `low` unless Hi-C or methylation agrees there; not where validated Hi-C
  can tell the two apart. **A named heuristic, on one genome**: two satellites differ in composition whatever their
  centromeres do. Only T2T mouse exercises it (from the FASTA: 1/20 → 18/20 right).

With Hi-C, the prior-free peaks test the troughs (Bonferroni over the windows, P ≤ 0.01): on most chromosomes they validate
the line, on a minority they support each other where they agree.
- A Hi-C call on a chromosome without a sequence call whose trough agrees (within 1 bin) takes the trough's finer interval and
  rises a level (`at_pair`). With mere support, a chromosome nothing else places takes its trough, `low`, where a clean
  prior-free peak that placed no call agrees (`at_denovo_pair`).
- Where the prior-free peaks neither validate nor support the troughs and Hi-C's cis signal acts (tested by the LTR cores or
  by trans contacts), the cis extremes test the troughs the same way (`atline.with_hic`): they are the readout that places
  Hi-C's calls there (*C. neoformans* in G1, whose prior-free peaks sit at the clustered telomeres: 5 of 14 agree, P 6e-5,
  `support`). A cis call whose trough agrees takes the trough's interval and is `medium` (`at_pair`), as a Hi-C-alone call
  with an agreeing trough is: two lines from two inputs at one locus. In G1 *C. neoformans* 5 calls rise, all on CENP-A.

## Repeat density (`repline.py`)

**Reads** the FASTA (every sequence counts) and the satellite engine's arrays. Where centromeres are built from transposons or
their remnants (*Cryptococcus*' Tcn clusters, *Chlorella*'s Osser insertions, *Chlamydomonas*' ZeppL clusters), each
chromosome's densest stretch of sequence repeated elsewhere in the genome marks its centromere, with no annotation; the LTR
line cannot read fungal Ty3 elements without a lineage ([future.md #47](future.md#47-centromeric-retroelements-the-ltr-line-cannot-read)).

**The readout.** Per window, the mean over its counted positions of log2 of the genome-wide count of the position's canonical
21-mer (of A, C, G and T): single-copy sequence adds 0, and no copy-number cutoff is chosen. K-mers that touch a satellite
array are left out of the readout (they still count copies): tandem arrays are the satellite line's, which has its own test;
with them in, six dense blocks gave holocentric *C. elegans* a near pass (P 0.052). Windows of 10-500 kb (one transposon copy
fills a smaller one); windows whose midpoint lies in a terminal eighth are left out (chromosome ends gather repeats whatever
the centromere: *Malassezia*'s densest class was the last 10 kb of every chromosome). Above 200 Mb of FASTA, 1 in 32 21-mers is
kept by a fixed hash (every copy of a kept 21-mer is kept, so its count is exact).

**The test** is the AT line's, on the negated readout, at P ≤ 0.01 after the correction for the window sizes: block shuffles
(2,100 of 8 genomes) passed it 6.5% of the time at 0.05. It runs where most chromosomes have no confident sequence call.
Development (rules fixed before computing): *C. neoformans* P 0.0004, all 14 peaks on CENP-A (median width
49 kb, centromeres 20-65 kb); *V. dahliae*, *Chlorella* and *Chlamydomonas* pass, every peak right; *C. elegans*,
*C. nigoni*, *B. mori*, *M. sympodialis* and tomato are silent; soybean and rice, whose centromeres are satellite arrays,
are silent once the arrays are left out.

**Calls.** The densest window at the chosen size, widened to its half-height run.

**What it may change.** Like a trough, it places only a chromosome nothing else places, `low` (`rep_only`). Where the line
passed its own test and its peaks agree beyond chance (`pass` or `support`) with Hi-C's prior-free peaks, or with the cis
extremes for a `cis_only` call, it is the second line of a Hi-C call on a chromosome without a sequence call that it agrees
with: the call takes its interval and rises a level (`rep_pair`). It never moves, opposes or raises a sequence call. Wherever
the AT trough may place or pair, it does so first.

## Repeat families (`famline.py`)

**Reads** the FASTA (every sequence counts copies) and the satellite engine's arrays. Where centromeres are clusters of
particular retrotransposon families (Tekay-type elements in *Solanum*, CRM in grasses) and no satellite family passes its
one-per-chromosome test, the FASTA still holds what the LTR line reads from an annotation: families concentrated at one
locus on most chromosomes. Repeat density cannot see them in large genomes, where the whole pericentromere is dense.

**The readout.** K-mers stand in for families. A fixed hash keeps 1 in 32 canonical 31-mers (every copy of a kept 31-mer is
kept, so its copy number is exact). Left out: 31-mers that touch a satellite array (a 31-mer with half its copies in arrays is
dropped), 31-mers with a period of at most 15 bp (tandem units), and 31-mers with fewer than 1.5 N copies (N chromosomes). A
31-mer is concentrated on a chromosome where at least 3 of its copies lie there and at least half of them in one window; a
family 31-mer is concentrated on at least half of the chromosomes. Each family 31-mer votes once per chromosome where it is
concentrated, at the center of its densest window. The profile is the votes per window (0.1, 0.3, 1 and 3 Mb).

**The test** is the AT line's, on the negated votes, at P ≤ 0.01 after the correction for the window sizes, with windows in a
terminal eighth left out. It runs where the repeat-density line runs: most chromosomes have no confident sequence call.
The test is not calibrated. Every family 31-mer votes at most once per chromosome, so one family that clusters once per
chromosome puts its votes one per chromosome wherever its clusters sit: a pass says that such a family exists, not where the
centromeres are. It passes in genomes whose repeat units were moved to random places, and holocentric genomes do not pass
it (Results). No calibrated null could be built from K-mers: grouping them into families either splits a diverged family, whose
parts then seem to coincide, or merges neighbors and leaves no power
([future.md #76](future.md#76-the-repeat-family-line-and-low-calls-it-contradicts)). Its peaks are leads.

**Calls.** The strongest window at the chosen size, widened to its half-height run.

**What it may change.** It places only a chromosome nothing else places, `low` (`fam_only`), after the AT trough and the
repeat-density peak; the caveat says the test is not calibrated. On a `low` sequence call that no other rule changed, outside
the terminal eighths, whose chromosome's family window lies elsewhere it says so (`fam_elsewhere`, in the caveat; `families` is
listed `against`) and the call stays: moving such calls onto the window would mostly have been right on the benchmark, but did
not bring the polyploid panel's homoeologous calls closer, so it waits for more truth where the line passes
([future.md #76](future.md#76-the-repeat-family-line-and-low-calls-it-contradicts)). It never pairs with Hi-C: a line whose test
is not calibrated cannot be one of the two tested lines that raise a call, so 0.16.1 removed 0.16.0's `fam_pair` (it had raised
no call in the release runs). It never moves, contests or raises a sequence call.

**Results** (every number of this section; rules fixed before computing). The 0.16.0 release runs the line where it runs in
cenfuse; the development scan ran its test on all 71 FASTA-only benchmark genomes.
- *Release* (FASTA-only runs, scored by the frozen scorer): the line runs in 39 genomes and passes in 18, its peaks right on
  183 of 211 chromosomes: 13 of the 16 *Solanum* genomes (the tomatoes LA1028, LA1777 and LA2134 fail, P 0.023-0.12, their
  peaks still right on 10 of 12 each), *Ae. geniculata* (14 of 14), *Chlamydomonas*, *Chlorella*, *N. crassa* and *V.
  dahliae*. Of the `low` calls it
  contradicts, its window is right on 89 of 110 (in 15 genomes: 13 *Solanum*, *Ae. geniculata*, *Chlamydomonas*) and the
  calls on 8; on the polyploid panel, moving them brought the homoeologous pairs no closer (64 pairs: 24 closer, 27 farther).
  Over all release arms, 189 of the 207 flagged calls with truth are wrong (16 genomes).
- *Development scan* (all 71): passes in 32 (also in the 10 *B. rapa* lines, soybean, maize, rice and sorghum, whose
  confident satellite calls keep the line from running), peaks right on 323 of 363 chromosomes. It fails in every mammal and
  chicken (their satellite arrays left out), in the four holocentric genomes (P 0.34-0.71) and in barley (P 0.20, though its
  peaks are right on 7 of 7); *Brachypodium*'s 5 chromosomes leave it `untested`; T2T mouse, whose centromeres lie in the
  terminal eighths it does not read, comes near (P 0.015) with every peak interstitial. Repeat density's peaks in the same
  scan: 0-2 of 12 right in the tomatoes, 0 of 7 in barley.
- *Calibration* (0.16.1). Rebuilt from each genome's own copies with every repeat unit moved as a whole to a random place on
  its chromosome, the test passes at its 1% level in 13-52% of draws in each of the 11 monocentric genomes tried (133 of 400:
  *Chlamydomonas*, *Chlorella*, *N. crassa*, *V. dahliae*, potato, *S. verrucosum* and five tomatoes). Synthetic genomes
  whose planted families cluster at independent random places pass in 1-10% (three settings of 100 genomes each). In tomato
  the 1-Mb peak on 7 of 12 chromosomes is one family (365 linked 31-mers).

## rDNA and telomeres

Reported features. A Hi-C peak beside a NOR never opposes or places a call: NORs meet in the nucleolus.

- **45S rDNA (NORs):** runs of ≥ 3 copies of the universal 18S motif, ≤ 100 kb apart. A hit followed by A whose next 19
  bases differ from nuclear 18S at most positions is organelle DNA, not counted.
- **5S rDNA:** runs of ≥ 3 whole 5S genes (Rfam 5S profile HMM, ≥ 90% of the gene, E ≤ 1e-5), ≤ 20 kb apart.
- **Telomeres:** the unit is found from scratch (5-12 bp, else 13-29 bp). An end holds an array with ≥ 300 bp of unit copies
  in the terminal 100 kb at half density, or ≥ 100 bp at half density ending within 1 kb of the sequence end (assemblers trim
  short telomeres); once the unit is known, ≥ 300 bp at a tenth density reaching the end. Interstitial arrays: ≥ 300 bp at
  half density outside the terminal 100 kb (old fusions; some centromere cores, as in chicken).

## Hi-C (`cenfuse.hic` and `hicfuse.py`)

Centromeres cluster in many nuclei, so a centromere's contacts with other chromosomes (*trans*) are strongest at the other
centromeres. Chromosome ends, the nucleolus and other hubs cluster too, which is why Hi-C is tested first.

1. **Bins.** 250 kb or coarser, by the library's depth (≥ 25 expected centromere contacts per bin). A chromosome needs 20
   bins; for small chromosomes bins go down to the coarsest size that gives the shortest one 20 (yeast 10 kb), while the
   largest keeps ≤ 4,000 bins and the genome ≤ 50,000. Replicate libraries of one plant are summed.
2. **Prior-free peak.** The Hi-C engine's own call, using no sequence call. It never acts unchecked: only where nothing can
   test Hi-C (step 8), or paired with a Rabl-axis line or an AT trough.
3. **Anchored score.** Each chromosome's mean contact with windows at the *other* chromosomes' `high` and `medium` calls (all
   calls where fewer than 3 are confident); its own call never shapes it. Over an unmappable array, the 9 mappable bins
   nearest the call stand in.
4. **Strong peak.** Passes the engine's three tests (bootstrap, leave one partner out, null replicates) and is not flagged
   `hub`, `rDNA` (beside a NOR) or `telomeric` (within 1 mappable bin of the end). A flagged peak may support, never move or
   oppose.
5. **Distances in mappable bins.** A peak agrees with a call within 1 bin, forgiving up to 4 unmappable bins (a peak at an
   array's edge supports a call inside it). It opposes only beyond 4 mappable bins, its core more than 1 bin from the call;
   in between it does neither (`hic_near`).
6. **The test (gate).** For chromosome x, the references are the other chromosomes' `high` calls (with `medium` if fewer than
   4 are `high`). x passes when the references' peaks agree with their calls on most of them, beyond chance
   (Poisson-binomial P ≤ 0.01), and on most of the whole reference set. Fewer than 3 references: `untested`. With no
   confident call the `low` calls are the references, and only clean peaks count; with fewer than 4 confident calls, the
   `low` calls join them on those terms (`hic_gate_refs` = `mixed`).
7. **Support tier.** Agreements at significant peaks (`medium` or `high`) beyond chance (P ≤ 0.01, ≥ 3) but on a minority:
   a call whose own significant peak agrees rises a level (`hic_support`), and a strong peak on the competing cluster chooses it. Hi-C never
   opposes, moves or places a call here. Randomized libraries put weak peaks on 3-6 calls per genome.
8. **Hi-C alone.** With fewer than 3 reference calls in the genome (point, regional, AT-rich centromeres), nothing tests Hi-C
   against the sequence. A passed Rabl-axis line that agrees with the prior-free peak gives a `medium` call (`fold_pair`,
   `align_pair`); else a clean `high` prior-free peak (or one held at `medium` only by `null_untestable`, too few partner
   chromosomes for the engine's null to run; not by `null_degenerate`, a null that ran but could not score it) gives a `low` call
   (`hic_denovo_only`; *C. tropicalis*: 6 of 7 right). A sequence call is kept: contested where a passed axis line, or under
   it the prior-free peak, opposes it (one level per opposing line), else flagged where a clean prior-free peak lies elsewhere
   (`hic_denovo_elsewhere`).

`hic_sequence_conflict` (per genome): Hi-C's prior-free peaks are confident but disagree with the sequence calls while the gate
fails. One line is wrong genome-wide; only a third line (the Rabl axis, `tes/`) can say which.

## The Rabl axis: the arm fold and the trans arm alignment (`fold.py`)

In a Rabl nucleus the centromeres sit at one pole and the arms extend toward the other. Two readouts test Hi-C without any
sequence call:
- **Arm fold (cis).** Both arms run side by side, so the chromosome's own map has an anti-diagonal ridge through the
  centromere: the mean observed/expected contact of bin pairs (c - d, c + d) at least 9 bins apart, scanned at centromeric
  index ≥ 12.5%. The fold core is the top center's half-maximum run.
- **Trans arm alignment.** Loci at equal distances from their own centromeres meet across chromosomes. It pools all partners,
  so it works where one map is too small for the fold (yeast: alignment 13/16, fold 6/16).
- **The test.** The prior-free peak overlaps or touches the core on most chromosomes, beyond chance (P ≤ 0.005 per readout). The alignment
  stands in for a failed fold only where the fold leans the same way (P ≤ 0.05): holocentric compartments can mimic the
  alignment (*C. elegans*), not the fold.
- **What it does.** Where the gate fails, a passed axis line that agrees with the prior-free peak moves a call both contradict
  (barley without `tes/`: 6 calls moved, all onto CENH3); a call they oppose is `contested`. The alignment pairs and
  opposes only through a clean prior-free peak, the fold through any. The support tier, the which-end test and
  `hic_denovo_elsewhere` are then off; the LTR-core and cis readings still act.

## Hi-C and the LTR line test each other (`hicfuse.py`)

Hi-C anchored on wrong calls finds nothing (potato: 0 of 11 centromeres when anchored on its sequence calls). So
where the gate passes on at most half of the chromosomes and `tes/` gives a signature, two more readings are tested against
the LTR cores (weak ones included), not the sequence calls:
- **The peak at the core.** Partner windows sit at the LTR cores. A chromosome agrees when its peak, made from the other
  chromosomes' windows, lies in its own core ± 4 bins (chance: that zone's share of the chromosome); only clean peaks count.
  Beyond chance (P ≤ 0.01, ≥ 3) the status is `support`, on most chromosomes `pass`. Where the peak lies at the core
  (`te_hic_pair`), a call on the core rises a level, a `low` call off it moves there (`medium`; at `pass` any call does, a
  `high` one staying `high`), and a chromosome without a call takes the core.
- **The contact crossing a bin (cis).** The two flanks of a centromere can touch more, or less, than DNA at the same distance
  elsewhere. The profile is the contact across each bin boundary (4 bins on each side), observed over expected. Its extreme
  marks the centromere in many plants and fungi, not in mammals, and its sign differs (a peak in tomato, a dip in potato), so
  the LTR cores choose the sign (Bonferroni over two). The chance is the larger of the zone's share and the rate in 200 maps
  simulated with the same coverage and distance decay. At `support` or `pass`: an extreme at the call raises it (`cis_agrees`);
  within 4 bins, nothing (`cis_near`); farther, at `pass` a `low` call moves to it and stays `low` (`cis_only`), and a confident
  call stays, flagged (`cis_elsewhere`), one level lower. At `support` a call it contradicts stays, flagged: a readout right
  on a minority of chromosomes raises the calls it agrees with and moves none (*S. verrucosum* with `tes/` and Omni-C, the
  first cis line at `support`: its 3 moves were all wrong). At either status a chromosome without a call takes the extreme,
  `low` (`cis_only`).
- **The cis extreme tested twice** (future.md #27). Where the LTR cores pass the cis signal, the library's trans contacts
  test it too (as below, without LTR cores). A confident call the extreme contradicts moves to it where the chromosome's own
  trans peak lies at that extreme (`cis_trans_pair`): the extreme has passed two genome-wide tests and is marked by the
  chromosome's own contacts, against one LTR core. **A named heuristic** (tomato with `tes/` and Hi-C: its 2 confident LTR
  misses move onto CENH3; TS-60, the one other genome whose cis signal passes both tests, 2 of 2): the moved call stays `low`.
- **Which end.** Where neither the gate, its support tier nor a Rabl-axis line speaks, Hi-C can still tell which end holds a
  call in a terminal eighth: the anchored score is higher there (goat 28 of 30). Where that holds on most other such
  chromosomes (at least 3), beyond the 1/2 of chance (binomial, P ≤ 0.01), the call rises a level (`hic_end`), unless no mappable bin
  separates it from its competing cluster. It needs no `tes/`.

## The cis and trans readouts of Hi-C test each other (`hicfuse.cis_trans_line`)

Without LTR cores, and where the gate passes on at most half of the chromosomes, the library tests its own cis signal:
- **Anchors.** For each sign, each chromosome's cis extreme outside its terminal eighths. Ends both fold and cluster, so over
  whole chromosomes the two readouts agree at the ends (*Ae. tauschii*: 4 of 7, none on CENH3).
- **The test.** The anchored test above, on these anchors (clean peaks; Bonferroni over the two signs). The lower-P sign is
  read; if the other sign passes too, the sign is unknown (`ambiguous`) and nothing acts.
- **What it may do.** Only at `pass`, on a chromosome whose other chromosomes pass too: the two readouts share the library's
  coverage, which alone can give agreement on a minority (randomized soybean reached `support`). The cis rules then apply,
  except that a call in a terminal eighth is neither supported, opposed nor moved (`cis_blind` where the extreme lies
  elsewhere), with one exception (future.md #32): a `low` call there moves to the interior extreme where the
  chromosome's own trans peak lies at that extreme (`cis_end_pair`); an acrocentric chromosome among metacentric ones has no
  such pair and keeps its call. **A named heuristic** (tomato; TS-60, the one other genome with the case, 2 of 2): the moved
  call stays `low`.
  Tomato from the FASTA and Hi-C: 12 of 12 right (from 1).

## Methylation: CHG (`cenfuse.meth` and `methfuse.py`)

The CHG line learns, one chromosome left out, what methylation looks like at the other chromosomes' confident calls, and
labels each candidate interval `centromere_like` at a signature score ≥ 0.5. It supports calls and can be a second line for a
Hi-C move; it never moves a call alone. Only samples of the assembled plant (or of unknown origin) decide. It needs `tes/` and
a modkit bedMethyl CHG file for every sample; otherwise it is skipped.

## CpG kinetochores (`cdr.py`)

In primates CENP-A sits in a CpG-hypomethylated dip of the methylated active array (the CDR); stickleback, sparrows, duck
and some reptiles and amphibians show it too, chicken, bats and many fish do not (VGP survey, benchmark.md). The other lines
find the array; the dip finds the kinetochore inside it. Input: the CG file of each deciding sample (modkit or ENCODE
bedMethyl, pb-CpG-tools, or a bedGraph of methylated fraction), on reference CpGs.
- **Profile.** CpG methylation per 5-kb bin, each pooling the calls of the three bins around it (so low coverage works);
  its background is the running median over ±500 kb.
- **Dips.** A run of at least 15 kb of array bins (one missing bin bridged), by two measures of the same depth: a *level*
  dip falls below half the background (depth: its methylation over the background's); a *spread* dip falls more than 3
  robust SD below the background, the SD being the running median of the bins' distance from the background over the
  same ±500 kb (depth: in those SD), and at least 34% below it (CDR-Finder's minimal valley height): more than one
  aberrant 5-kb bin can make through the pooling over three bins, so neither such a bin nor a wiggle of a saturated
  array is a dip. Tracks of the same DNA encode methylation on different scales (a + b × methylation:
  genome-wide median bin 0.36-0.92 across the benchmark tracks). The level depth moves with that scale: the same HG002
  CDRs fall to 0.14 of their background in the ONT track and only to 0.41 in the PacBio bigWig, a level that genic
  methylation rotated onto the arrays reaches too (the rotation test then fails on 46 of 46 chromosomes). The spread depth
  does not move. A call's kinetochore candidate is the deepest level dip in its centromeric arrays, or the deepest spread
  dip where there is none: the level dips' places agree best between two technologies on the same assembly.
- **Width.** The half-depth run of the arrays' most prominent valley (topographic prominence, so an
  array edge is none), the bins below the midpoint between its bottom and the arrays' median, read on each 5-kb bin's own
  methylation (bins with at least 10 CpG rows). Pooling over three bins leaves a CDR's half-depth run as it
  is but triples the width of noise valleys (simulated noise: a median 1 → 3 bins), and the rotations are made of those.
  Like the spread depth, it does not move with a track's scale. In the three macaques every call holds a dip, but methylation from elsewhere makes
  dips as deep (a median P of 0.96-1 by either depth); the crab-eating and rhesus macaques' CDRs stand out by width (a
  median 115 kb; P ≤ 0.05 on 18 and 15 of 21 calls). Width finds some valley in any array, so it never places a
  kinetochore.
- **Positional test.** The call's arrays hold dips deeper, or a valley wider, than the chromosome's methylation, rotated
  by 200 random offsets (any offset: the observed window is one of the rotations, so P is exact up to Monte Carlo error),
  makes on the same array bins. The three measures are combined by minP: each window's measures are ranked among the
  observed and the rotated windows, and the observed window's smallest rank is compared with the rotated windows'
  smallest ranks (exact, whatever the measures' correlation). Width alone is not enough: the partially methylated domains
  of cultured cells (the T2T chimpanzee track) and subtelomeric caps make wide shallow valleys anywhere, which failed the
  length test on every chimpanzee call; the depths carry those tracks. A call spanning much of its chromosome cannot reach
  a small P (a 3-Mb call on a 60-Mb chromosome: about 0.05 at best).
- **The test (gate).** References: the other chromosomes' `high` calls (with `medium` if fewer than 4 are `high`; fewer than
  3: `untested`). Their positional P values, combined by Fisher's method, must give P ≤ 0.01: many moderate dips count, as
  wide ape calls rarely reach P ≤ 0.05 one by one. It is computed on the provisional calls, before the final integration.
- **Role.** Where the gate passes, each call with a dip gets its kinetochore (`cpg_kt_*`). A sequence call whose arrays have
  data, no dip and no positional P (a valley wider than chance is CDR evidence in the call), while the deepest level dip on its
  chromosome's other satellite arrays, else the deepest spread dip, passes the positional test of its two depths there (the
  dip's own measures), is contested (`cpg_elsewhere`, the dip in `cpg_alt_*`): one level below what the agreeing lines give
  where it is the only contest, one level more where another line contests too. CpG never moves a call: in the truth genomes such a
  dip was the centromere once (stickleback) and elsewhere twice.
- **Coverage holes.** Where the calls' arrays (the median call) have data in under 0.8 of the genome's share of bins, the
  genome is reported as a coverage hole; each call's share is in `cpg_array_data`. Tracks piled up at MAPQ 60 give 0.36-0.68
  of the genome's share, at MAPQ ≥ 1 0.97 or more (most reads in arrays map at MAPQ 20-59).

## Integration (`fuse.py`)

The rules act in the code's order, which is the specification (`fuse._decide`); a later rule may override an earlier one:
1. `fuse._context`: a strong validated peak on the competing cluster, more than 1 mappable bin from the call, chooses that
   cluster (`hic_chose_alt`; on the support tier too, unless the LTR-core or cis reading acts there). The rules below read it
   as the sequence call.
2. `fuse.CASCADE`: the first rule that applies sets the call (first block of the table).
3. `fuse.PHASES`: the readings that act where the gate fails, each on the call the earlier ones left: Hi-C at the LTR cores,
   the cis signal, then the AT trough and the repeat lines (second block).
4. `fuse._refine` and `fuse._cpg`: an agreeing trough or repeat-density peak gives a Hi-C call its interval, the Rabl axis contests a
   kept call, and the CpG line contests last (third block). Support, what is against the call and its confidence follow.

A long-unit array that stays the call with no support from Hi-C, `tes/` or methylation is no call (`long_unit_unsupported`):
the chromosome is decided again as one without a sequence call, and the array stays in `seq_start`/`seq_end` and the loci. The
`rule` column names the last rule that set the call's place or relation ([outputs.md](outputs.md#relation-and-confidence)).

**Confidence.** A kept sequence call rises one level for agreeing tested Hi-C (any readout, counted once) and one for CHG
methylation (two at most); it falls one level when contested (one per opposing line under a Rabl-axis line; methylation support
at the call cancels one). A call that Hi-C moved or placed is
`medium`, `high` with a second line at its peak or an agreeing trough or repeat-density peak; one Hi-C readout alone gives `low`
(`hic_medium_peak`, `hic_denovo_only`, `cis_only`; `medium` with an agreeing trough or repeat-density peak), and `cis_trans_pair` and
`cis_end_pair` give `low` always. A trough alone is `low` (`medium` where Hi-C's prior-free peaks validated the troughs), a
repeat peak alone `low`.

| rule | situation | result |
|---|---|---|
| **`CASCADE`** | | |
| `fold_pair`, `align_pair` | gate failed; a passed Rabl-axis line and the prior-free peak agree, the sequence call at neither | the prior-free core, `offset` (`functional_only` without a sequence call), `medium`; a `high` call adds `moved_with_fold` or `moved_with_align` |
| `hic_strong_peak` | no sequence call; a strong validated peak | the Hi-C core, `functional_only` |
| `hic_medium_peak` | no sequence call; gate passed; a clean `medium` peak | the Hi-C core, `functional_only`, `low` |
| `hic_denovo_only` | no sequence call; fewer than 3 reference calls in the genome (no gate can be tested); no Rabl-axis pair; a clean prior-free peak (`high`, or `medium` held only by `null_untestable`) | the prior-free core, `functional_only`, `low` |
| `none` | no sequence call and none of the above | no call (the phases may still place one) |
| flag `hic_alt_ambiguous` | validated Hi-C agrees, but no mappable bin separates the call from its competing cluster | the call, without Hi-C support: Hi-C cannot tell the two apart |
| `hic_agrees` | validated Hi-C agrees with the call (within 1 bin) | `concordant`, one level up |
| `hic_support` | support tier; the call's own significant peak agrees | `concordant`, one level up |
| `hic_end` | gate failed, no Rabl axis; a call in a terminal eighth at the end with the higher anchored score, as on most such chromosomes | `concordant`, one level up |
| `hic_moved` | a strong validated peak more than 4 mappable bins away, its core not reaching the call; the call `low` or `medium`, or `high` with a second line at the peak and not at the call (`moved_with_*`): a strong core of a decisive LTR signature within 4 bins; methylation `centromere_like` at the peak and a `dip` at the call; or the signature's copies significantly younger at the peak, where the premise holds | the Hi-C core, `offset` |
| `hic_offset` | the same, a `high` call without a second line | kept, `contested`, one level lower |
| `sequence` | anything else | the sequence call (`sequence_only`); flags `hic_near` (a validated peak within Hi-C's resolution, neither way) and `hic_disagrees_weak` (one farther, too weak to oppose) |
| **`PHASES`** | | |
| `te_hic_pair` | gate failed; Hi-C's peak at the LTR core, tested genome-wide (`support` or `pass`) | a call on the core rises a level; a `low` call off it (at `pass` any call) moves to the core, `medium` (a `high` one stays `high`); a chromosome without a call takes the core, `medium` |
| `cis_agrees`, `cis_only`, `cis_elsewhere`; flag `cis_near` | gate failed; the cis signal tested against the LTR cores (`support` or `pass`) or by trans contacts (`pass`) | the extreme at the call: one level up; within 4 bins: nothing; farther: at `pass` a `low` call moves to it, still `low`, and a confident call stays, `contested`, one level lower; at `support` both stay, flagged; a chromosome without a call takes the extreme, `low` |
| `cis_trans_pair` | the LTR cores and the trans contacts both pass the cis signal; this chromosome's own trans peak lies at its extreme; a confident call elsewhere | moves to the extreme, `low` (a named heuristic) |
| `cis_end_pair`; flag `cis_blind` | the cis signal tested by trans contacts; a call in a terminal eighth, the extreme elsewhere | a `low` call moves to the interior extreme where its own trans peak lies, `low` (a named heuristic); else kept, flagged |
| `at_only`, `at_denovo_pair` | no call yet; AT troughs passed (or validated by Hi-C's prior-free peaks); or a trough and a clean prior-free peak that placed no call agree, in a genome where the two agree beyond chance | the trough, `low` (`at_only` with validated troughs: `medium`) |
| `rep_only` | no call yet; the repeat-density test passed | the densest repeat stretch, `low` |
| `fam_only` | no call yet; the repeat-family test passed | the strongest family window, `low` |
| `at_chose_alt` | AT troughs passed; a `low` call with a competing cluster; the deepest trough window on the cluster, not on the call; not where validated Hi-C can tell the two apart | that cluster, `low` unless another line agrees there (a named heuristic) |
| flag `fam_elsewhere` | the repeat-family test passed; a `low` sequence call left as it is, outside the terminal eighths; the family window elsewhere | the call stays; the caveat names the window |
| **`_refine`, `_cpg`** | | |
| `at_pair`; flag `at_elsewhere` | a Hi-C call on a chromosome without a sequence call; a tested trough agrees (tested against the prior-free peaks, or against the cis extremes for a `cis_only` call) | the trough's interval, one level up; a validated trough elsewhere is flagged |
| `rep_pair` | the same, no agreeing trough; a tested repeat-density peak agrees (the repeat-family window never pairs) | its interval, one level up |
| flag `hic_denovo_elsewhere` | Hi-C alone in the genome; a clean prior-free peak beyond 4 bins of a kept sequence call | the call stays; the caveat says where |
| `fold_offset`, `align_offset`, `hic_denovo_offset` | a passed Rabl-axis line, or under it the prior-free peak, opposes a kept call | `contested`, one level lower per opposing line |
| `cpg_elsewhere` | the CpG gate passed; the call's arrays have data and no dip, while another array holds a positional one | `contested` ([CpG kinetochores](#cpg-kinetochores-cdrpy)); never moved |
| `meth_support` | a kept call that methylation alone confirms (`centromere_like`) | `concordant`, one level up |

CpG kinetochores are added after the integration; the CpG gate is computed on the provisional calls. CpG never moves or raises a
call.

## Reading the integration (`morphology.py`, `interpret.py`)

Nothing here changes a call; every statement is read from the tables cenfuse writes.
- **`basis` and `caveat`.** One sentence each: what the call stands on (the satellite and its size in the call, an LTR core and
  the array it replaced, a move and its second line, each agreeing line, the kinetochore); then what speaks against it (the
  overruled array, each lowering flag, each opposing line and where it points, gaps) and the input that would settle it, where
  cenfuse's own tests say one can (two arrays of one satellite: CENP-A/CENH3 data; satellite-free or capped calls: Hi-C or
  `tes/` when not given).
- **Line tests** (`asm.line_tests.tsv`). Every pair of lines (but cis with the LTR cores, which chose its sign) scored the
  same way: how often they agree, against chance from each chromosome's geometry. Two lines that agree beyond chance on most
  chromosomes validate each other without truth.
- **Morphology.** Arms are measured from the kinetochore where the CpG line placed one, else from the call's midpoint:
  centromeric index, arm ratio and Levan class (Levan, Fredga & Sandberg 1964); the summary's karyotype counts only confident
  calls. These are sequence measurements: CHM13's acrocentric chr21 and chr22 come out submetacentric/subtelocentric, because
  their T2T short arms (rDNA included) are 11-15 Mb.
- **Completeness.** Telomere arrays at the ends, and assembly gaps (runs of ≥ 10 N) on the chromosome and in the call.
- **The report.** Per genome: the lines and their tests, the calls by confidence, the uncertain calls with their caveats,
  morphology and completeness, and the other features with their 3D association. Where no call is confident and no line
  passed a genome-level test, it says so first: each call is then one untested line's best candidate, as in holocentric
  genomes and satellite-free genomes from the FASTA alone.

## Constants

Every threshold is fixed, with its reason. cenfuse's are in `cenfuse/params.py`; each engine keeps its own (`sat/params.py`,
`hic/params.py`, `meth/params.py`). [constants.md](constants.md), generated from the code, lists all of them. The main ones:

| what | value | reason |
|---|---|---|
| test levels | 0.01: Hi-C gate and support tier, Hi-C on LTR cores, cis, which end, AT vs Hi-C, acrocentric, LTR overruling, age, repeat density and repeat families (after Sidak), CpG gate (Fisher). 0.05: one-per-chromosome, AT troughs (after Sidak), CpG dip, a pooled seed's cores. 0.001: LTR cores, recruitment and a pooled seed's Fisher combination | conventional evidence levels |
| most; shares | more than half; tandem, whole and competing-cluster shares 1/2; half of the signature's elements | majority rules |
| Hi-C agreement; opposition | within 1 bin; beyond 4 mappable bins | the Hi-C engine's placement error; its partner half-window |
| Hi-C references | ≥ 3; `high`, `medium` | fewer give no majority test; wrong anchors blur Hi-C |
| Hi-C bins | 5 kb-2 Mb; ≥ 20 per chromosome; below 250 kb, ≤ 4,000 per chromosome and ≤ 50,000 per genome | each size a multiple of the finer ones; a centromere needs flanks; memory |
| raise of a kept call | 2 levels at most | `low` → `high` |
| Rabl axis | pairs ≥ 9 bins apart; centromeric index ≥ 12.5%; P ≤ 0.005 per readout, 0.05 to lean | beyond the ± 4-bin window; Levan's telocentric bound; 0.01 shared by two readouts |
| chromosome ends | terminal eighth; end score over the terminal 5% (≥ 4 bins) | Levan's telocentric bound; a named heuristic of about one half-window |
| satellite sizes | 20 kb a centromere-sized cluster; 5 kb a weak one and the support floor; 50 kb `low_satellite`; 100 kb a strong family | priors |
| locus fold | 2 | a named heuristic: interspersed families 1.1-1.3, plant TE-derived arrays 2.5-3.4, satellites ≥ 15 |
| LTR core | ≥ 3 independent insertion sites | three independent insertions make a locus (*Chlorella*'s centromeres hold a few elements each) |
| LTR fold | 10 | measured without truth: in 42 polyploids whose best marker is not CRM, cores agree with independent FASTA-only calls at median 0.30 at or below 10-fold, 0.75 above |
| LTR recruitment | ≥ 5 elements; ≤ 5 rounds | fewer cannot reach P ≤ 0.001; results do not depend on them (152-155 of 172 right cores over 1-4 rounds, 3-10 elements) |
| age | 3 dated copies per site; peaks ≥ half the highest | as the three insertions of a core; the core's half-maximum |
| AT troughs | 1-kb bins; 1-500-kb windows; troughs ≥ 5 windows apart; class ≤ 1.5 N | *Malassezia*'s 3-5-kb centromeres; the known regional centromeres; one centromere gives several close troughs; a named heuristic |
| repeat density | 21-mers; windows 10-500 kb; a window half outside satellite arrays; 1 in 32 21-mers above 200 Mb; P ≤ 0.01 | odd K, chance recurrence ≤ 2e-4 in 0.5 Gb; one transposon copy fills a smaller window; a majority; memory; shuffles passed 6.5% at 0.05 |
| repeat families | 31-mers, 1 in 32 in every genome; windows 0.1, 0.3, 1 and 3 Mb; ≥ 3 copies to be concentrated on a chromosome, ≥ 1.5 N copies in all; periods ≤ 15 bp left out; P ≤ 0.01 | odd K, the longest a 64-bit integer holds (copies at 90% identity still share 4% of their 31-mers); the rate the line was developed at; a ×3 series from fungal to large-genome centromeres; the fewest that show a cluster beyond a pair, on half of N; a period of half the K-mer or less is a tandem unit; as the repeat-density line |
| random draws | 10,000 (one-per-chromosome, AT troughs, both repeat lines); 200 (cis maps, CpG rotations, the Hi-C engine's null) | P resolution 1e-4; 0.005, 10× finer than the CpG dip's 0.05 |
| telomeres | ≥ 300 bp at half density; ≥ 100 bp at the end; a tenth density at the end | ~40 copies of a 7-bp unit; trimmed telomeres; a named heuristic |
| chromosome set | the first tenfold drop in length after half the assembly; a gradual tail cut | a named heuristic |

## Choices tested and rejected

| question | chosen | rejected, and why |
|---|---|---|
| fixing Hi-C from inside the Hi-C data | anchor to the sequence calls | random restarts and caps: accuracy unchanged; in potato every start converges on a non-centromeric cluster |
| Hi-C windows over unmappable arrays | the 9 nearest mappable bins | padding every window: lost agreements in three apes |
| Hi-C peaks at chromosome ends | flag only peaks within a bin of the end; they support, never move | an 8-bin end window: moss centromeres 3-6 bins from the end could not move wrong calls |
| Hi-C beyond chance on a minority | support agreeing calls, choose the competing cluster | letting it oppose or move: a disagreeing peak is weak evidence; moving any `low` call to a strong peak: moss 7 right, 1 wrong |
| genomes without sequence calls | clean prior-free peaks as `low` calls, `medium` with a Rabl-axis line | a Rabl polarity test: passed on randomized Hi-C too |
| an end-distance expected for trans contacts | not used | moss 23 → 25 of 26, but sheep 9 → 5, chicken 7 → 6, goat 12 → 11 |
| a cis profile | observed/expected, sign and level tested against the LTR cores, chance from simulated maps | the plain mean contact with the zone's share as chance: randomized Hi-C passed. Against one locus per chromosome: too little power. Against the prior-free peaks: no agreement in tomato, potato or rice |
| Hi-C windows on LTR cores | fixed at the cores | moved to the best score nearby: randomized Hi-C passed too |
| the cis signal without `tes/` | tested by trans contacts, terminal eighths out, acting at `pass` | whole chromosomes: *Ae. tauschii* passed with no agreement on CENH3; acting at `support`: randomized soybean reached it; cis + trans as `medium`: 2 near misses in *C. neoformans* would be confident errors |
| Hi-C depth fragility | depth guidance ([benchmark.md](benchmark.md)) | a split-half check: unstable halves did not predict each other, and Hi-C time doubles |
| an LTR core vs a confident satellite call | replace only if the core is ≥ 10-fold and marked significantly more | always replace: wrong on *B. rapa* A09 and soybean chr1/chr6. A bootstrap of the core position: lost *P. annua* Pa1A |
| LTR core floor | 3 independent insertion sites | 5 elements: missed *Chlorella*'s centromeres |
| LTR elements without a lineage (fungi) | not used | pooled, they lifted the honesty cap in fish and *Oryza*, or qualified in 9 fish genomes without truth ([future.md #47](future.md#47-centromeric-retroelements-the-ltr-line-cannot-read)) |
| several qualifying LTR candidates | CRM first, then the most specific | chromosomes with a core first: agreed less with independent satellite calls |
| insertion age as a locator | only among comparable peaks, and as a second line for Hi-C | the young copies' core alone lies on CENH3 on 128 of 172 chromosomes, the density core on 155 |
| CpG dip inside a call | the deepest run below half the ±500-kb running median (level dip), else the deepest spread dip | the deepest bin: lands on array edges. The longest run below half the array median: on the CDR 52/90. The most prominent valley, its half-depth run (the earlier rule): on CENP-A 57/67 right calls in five truth genomes (deepest run 61/67, same width) |
| CpG rotations | any offset | offsets at least a tenth of the chromosome away (the earlier rule): P ≤ 0.05 for 7%, 15% and 32% of null calls spanning 5%, 15% and 30% of their chromosome |
| CpG gate pooling | Fisher's method over the references' P | the count of P ≤ 0.05 (binomial at 0.05, the earlier rule): with exact rotations it fails chimpanzee and gorilla (P 0.03), whose wide calls rarely reach 0.05 one by one |
| CpG measures | level depth, spread depth and the width, combined by minP; the level dip places the kinetochore, else the spread dip; the width never places | chosen on the six genomes with two tracks of the same assembly (HG002 ONT and PacBio; five apes, ONT R9 and HiFi) by rules fixed beforehand and blind to truth. Level alone fails HG002 PacBio and CHM13 (nanopolish); spread alone fails bonobo, whose small calls take their spread from the flanks; genome-wide quantile ranks and the dip's length score lower. Spread dips must fall at least 34% below their background: heavy-tailed noise and single aberrant bins pooled over three made spread dips. Level + spread swapped genomes against the earlier width test (chimpanzee, HG002 HiFi and siamang gained; gorilla, the three macaques, duck and two fish lost: a review after the truth had been seen, disclosed). With the width, read on each 5-kb bin alone: the gate passes at least as often on all 12 paired tracks and on both gorilla tracks (25 and 24 of 25), all 6 pairs pass on both tracks, simulated misplaced tracks pass on 0.07% of genomes (at most 0.2%); per-call verdicts agree less between the two tracks of a pair (0.66 vs 0.82). The pig-tailed macaque's dips off CENP-A (4 of 9) rest on cross-species ChIP; the crab-eating macaque's own ChIP agrees with its width-passed dips on 12 of 15. Placing by the level dip first: the pairs' kinetochores within 50 kb of each other 0.85 (spread first 0.81) |
| CpG null | positional (rotation) test of the dips' depths and the valley's width for the gate | dips in other satellite clusters: passed on misplaced methylation in gorilla. The dip's length (the earlier test): failed on every chimpanzee call (cell-line partially methylated domains) |
| CpG choosing between arrays | not done; a dipless call is contested where another array holds a positional dip | least-methylated window: on the CDR 3/23; most prominent valley across clusters: wrong cluster on 8-10/25 |
| a centromeric region (call + nearby clusters) | not reported | precision 0.07 in apes: an interspersed family chains whole chromosomes |
| the chromosome set | tenfold drop, tapering tail cut, optional list | the drop alone kept 5,614 extra scaffolds over 121 polyploid assemblies; other cuts dropped dot chromosomes or microchromosomes |
| families counted together for dominance | a kind pools the families of one sequence; a chromosome counts for a family when its largest kind holds a related one (same sequence, or period within 15% of the shorter) | sequence alone: changes the main family on 90 of 3,459 chromosomes, with no truth behind it |
| dominance within subgenomes | not used | real A/B groups lift 27 of 292 capped calls in 20 polyploids, random splits of the same chromosomes 9 |
| 45S arrays of organelle DNA | the base after the motif and the next 19 | one strand only: large real NORs hold both. The next base alone: 9% of the hits followed by A are nuclear |
| telomere stretches under half density | used only after the unit is chosen | decoy units in loose stretches could block the real unit |
| an AT window per chromosome | one window per genome | choosing it by the trough's prominence loses right calls |
| an LTR core under tied density maxima | the first maximum, as before | half of all core profiles have tied maxima in separate half-maximum runs; taking the run with more elements moved 24 of 634 final cores (2 wrong to right, 4 right to wrong) and changed 3 seeds |
| the repeat-density readout | mean log2 copy number of 21-mers, satellite arrays and terminal eighths left out, windows of 10-500 kb, P ≤ 0.01 | the share of 21-mers seen twice: low-copy duplications at chromosome ends outweighed *C. neoformans*' centromeres (P 0.18). Seen ten times: a cutoff with no reason (it passed, 14 of 14). Arrays in: holocentric *C. elegans* nearly passed (P 0.052) on six dense satellite blocks. Ends in: *Malassezia*'s chromosome ends formed the deepest class. P ≤ 0.05: block shuffles passed 6.5% of the time |
| a gene-poor competing cluster over a gene-rich call | not used | one chromosome each way in the benchmark (horse1 and siamang) |

## Known risks

- **LTR calls outside CRM genomes.** In 13 tomatoes, 14 of the 72 confident calls with `tes/` are wrong (19%); calls that
  merely agree with the core are the worst (5 wrong of 18). No statistic of the signature tells wrong cores from right ones.
  Treat such a `medium` call as a location to check.
- **A family learned beside the centromere.** *Arabidopsis*' Athila (13.3-fold) flanks CENH3 by 0.46-1.38 Mb and only flags
  because CEN180 passes its test. Where the satellite fails, such a lineage would decide and place calls beside the centromere.
- **Loose agreement.** A call agrees with a core up to a gap as wide as the wider of the two, so a call at the edge of a wide
  core is `te_concordant` while CENH3 can sit at the other edge (*Ae. geniculata* chr05U, in a 9.3-Mb core;
  [future.md #45](future.md#45-confident-calls-at-the-edge-of-a-wide-ltr-core)).
- **Thin cores** (`weak_core`) can replace right satellite calls that failed their test (*S. verrucosum*: 2, both `low`).
- **A centromere that left its satellite and LTR families** is found only where the centromeric LTR families form a second,
  significantly younger peak there: insertion age then makes that peak the core (`te_age_core`), and Hi-C confirms it
  (*P. annua* Pa1A, without CENH3 truth). Recent dispersed insertions at a non-centromeric site could tip that choice.
- **The Hi-C test takes the sequence calls as references.** Where they are wrong genome-wide, only the Rabl axis, the LTR cores
  or the cis-trans test let Hi-C act.
- **Hi-C alone** trusts the strongest trans cluster to be the centromeres; where it is something else (or the genome is
  holocentric), its `low` calls are wrong.
- **One wrong Hi-C peak** can pass in a genome where Hi-C mostly works; only `high` calls need a second line to move.
- **Hi-C depth.** Tests flip at their threshold (both orangutans fail at half depth); where every Hi-C test fails, the calls
  are those without Hi-C. Give at least ~70 M pairs to a 3-Gb genome. Hi-C cannot see inside an unmappable array.
- **Cis tested by trans** would pass a structure that both folds and clusters away from the centromere (none seen outside the
  terminal eighths), and cannot see the terminal eighths; a `low` call there moves only to an extreme that is the chromosome's
  own trans pair (`cis_end_pair`). Its two readouts share one library: the moves it makes from it, `cis_end_pair` and
  `cis_trans_pair`, were shaped on tomato (TS-60 repeats them) and stay `low`.
- **The AT test** passes on about 4% of shuffled genomes by design, and then places `low` calls and can move `low` calls to
  their competing cluster. It is a little liberal where chromosomes hold room for fewer than about 5 troughs.
- **The one-per-chromosome test** has no power with few chromosomes, and passes any family that happens to lie once per
  chromosome.
- **The repeat-density line** was developed on 23 genomes and passes in four, every peak right; block shuffles pass its
  test 6.5% of the time at 0.05, so it acts at 0.01. It cannot see a centromere in a terminal eighth (11 of moss's 26), reads
  only interspersed repeats (a satellite centromere is the satellite line's), and leaves low calls it contradicts
  (*Chlamydomonas*: 17 of 17 peaks right, 6 of 17 calls). Recent transposon bursts away from the centromeres could make it
  pass on the wrong clusters; alone it places only `low` calls, and where it agrees with a Hi-C call it raises that call a
  level.
- **The repeat-family line**'s test is not calibrated: with each genome's repeat units moved to random places it still passes
  13-52% of the time at its 1% level, because one family that clusters once per chromosome passes wherever it clusters
  (Repeat families, Results). Its peaks are leads: alone it places only `low` calls, it never pairs with Hi-C, and the `low`
  calls it contradicts stay, flagged.
- **The CpG line** was designed on the primate tracks it is benchmarked on; its depth measures were chosen on the paired
  tracks without truth, and the truth genomes test it; the width joined them after a review that had seen the truth, on
  grounds without it (gate passes and agreement on the paired tracks, calibration). Wide calls weaken its test: gorilla's
  3.5-8 Mb calls pass by the width alone, and the duck's calls, up to half of a microchromosome, cannot reach a small P
  under rotations within the chromosome. The width is less specific than the depths: in horse, whose CENP-A domains are as
  methylated as its arrays, 13 of 24 calls are positional by their width (none by the depths), and with Hi-C the test
  passes there.
- **Morphology** is only as right as the call; an incomplete chromosome shifts the index.
- **Old satellite outputs** read with `--satellite-dir` keep their version's rules; rerun without `--satellite-dir`.

## Engineering notes

Resuming, the output-folder lock and `run_info.json`, as users meet them: [outputs.md](outputs.md#running).
- **Reproducible.** On Linux, FasTAN runs with address-space randomization off (`setarch -R`) and MAFFT's refinement
  single-threaded, so fresh runs are identical at any thread count; without `setarch` (macOS) about 1% of FasTAN's array cuts
  can differ between runs, with a warning.
- **Checkpoints.** The satellite engine's stages are keyed on their results version, parameters and the FASTA's md5 (from
  classify on, `genes/` and `tes/` too). cenfuse's stages `telomere` to `cpg` are keyed on the code they run (the syntax trees
  of their import closure, without comments or docstrings; most also on cenfuse's version) and on the content of their inputs;
  the integration and outputs always rerun. A stale satellite-search stage (seeds, families, annotate) or a changed chromosome
  set stops the run and names the `--force` to use. Files are written whole or not at all.
- **Locks.** A lock on this host is stale when its process is gone or its pid was reused; one from another host, or one that
  cannot be read, after 15 minutes without a heartbeat.
- **Tests.** In the `cenfuse` environment, `python -m pytest -q -n 8` runs about 1,230 tests on small synthetic fixtures
  (`-n` needs pytest-xdist, in `environment.yml`). `CENFUSE_REQUIRE_TOOLS=1` turns skips for missing tools into failures.
