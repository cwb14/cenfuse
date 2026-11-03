# Future development

What cenfuse 0.9 leaves open. Each entry gives the evidence (numbers), a proposed approach and the check it needs: a
rerun, or truth data. Nothing here is implemented unless an entry says so.

The evidence comes from two stress tests: a run of 0.8 on a 106-genome allotetraploid plant set (2,320 calls on 2,355
sequences: 54% `high`, 23% `medium`, 23% `low`) and the CENH3/CENP-A benchmark ([benchmark.md](benchmark.md)). The
problems are general: the plant set only exposed them. 0.9 fixes what could be verified without rerunning either set or
collecting new CENH3 data; what needs a rerun or new truth is here. Denominators (n/N) count the species with the data a
check needs.

| # | open problem | needs |
|---|---|---|
| 1 | [weak (single-sequence) family models](#1-weak-single-sequence-family-models) | rerun from the families stage |
| 2 | [the 40-family cap](#2-the-40-family-cap) | rerun with a higher cap |
| 3 | [satellites specific to chromosome subsets vs the honesty cap](#3-satellites-specific-to-chromosome-subsets-vs-the-honesty-cap) | CENH3 truth, a null |
| 4 | [the CRM decisive test with subgenome-specific clades](#4-the-crm-decisive-test-with-subgenome-specific-clades) | CENH3 truth, a null |
| 5 | [a faster search](#5-a-faster-search) | output identity on stored annotations |
| 6 | [seed clustering](#6-seed-clustering) | cluster identity on every stored seed set |
| 7 | [the FasTAN upstream report](#7-the-fastan-upstream-report) | the user's decision, then a fixed release |
| 8 | [large unplaced scaffolds in the chromosome set](#8-large-unplaced-scaffolds-in-the-chromosome-set) | the 121 chromosome lists as truth |
| 9 | [CENH3 truths for allopolyploids](#9-cenh3-truths-for-allopolyploids) | new data |
| 10 | [chromosome-set limits](#10-chromosome-set-limits) | known karyotypes |
| 11 | [the LTR retrotransposon line](#11-the-ltr-retrotransposon-line) | held-out truth; for 11a, a rerun and a per-rule diagnostic switch |

## 1. Weak (single-sequence) family models

**Evidence.** A family whose sample has fewer than 20 monomers (`min_hmm_seqs`) gets a single-sequence model: its
consensus alone, with no profile of how the copies vary (`hmm_mode` = `single` in the families table). 214 such models
in 68/105 species; 11 are the main family of 131 calls. A consensus alone may find diverged copies less well; how far
the arrays of these families, and the calls that stand on them, are affected has not been measured.

**Approach.** Build the profile from more monomers: pool the unit of every seed in the family's cluster, search again
with the single model at a relaxed cutoff and realign what it finds, or borrow monomers of the same family from the
other genomes of a complex. Keep `single` for families that still have too few.

**Check.** Rerun the satellite engine from the families stage on the genomes with a weak main family and compare arrays
and calls; the quantity to watch is how many calls move and which way. Truth: a CENH3 call set for those genomes (9).

## 2. The 40-family cap

**Evidence.** At most `max_families` = 40 clusters per complex are modeled and searched, since annotate time grows with
models x genome size. In *Sonchus oleraceus*, 60 clusters of >= 20 kb were dropped, the largest 107 kb (the log reports
the count and the largest cluster dropped, once per complex and only when one was dropped). What the dropped clusters
are, and whether any is a centromeric satellite, is not recorded.

**Approach.** Replace the count by a cost budget: price a model by its length and the genome size, and now also by the
per-family CPU that `search_cost.tsv` records, and annotate the largest dropped clusters while the budget allows. Write
the dropped clusters (seeds, span, unit length) to a table, so that what was left out is visible.

**Check.** Rerun the genomes with dropped clusters >= 20 kb with a higher cap (start with *S. oleraceus*) and compare
roles, main families and calls with the 40-family run; the cost table prices the extra search.

## 3. Satellites specific to chromosome subsets vs the honesty cap

**Evidence.** A genome in which no family is the largest satellite cluster on a majority of the chromosomes gets every
call capped `low` (`no_dominant_satellite`), because satellite evidence alone cannot tell a centromeric array from an
interstitial one. In an allopolyploid each subgenome can have its own satellite, so no family is the largest on most
chromosomes although each chromosome's own family may be the right one. In 30/105 species the top families differ
between subgenomes A and B; 85 calls were capped `low`.

**Approach.** Judge dominance within chromosome subsets when the subsets are given (a table beside the FASTA, like
`<genome>.chromosomes.txt`: chromosome, group). Lift the cap for a group only when one family dominates inside it. With
no table nothing changes.

**Check.** A CENH3 call set for allopolyploids with subgenome-specific families (9) says whether the lifted calls are
right. The null: the same lift with a random split of a diploid's chromosomes into groups of the same sizes, to count
how often the cap lifts by chance. Without groups the benchmark outputs must be identical.

## 4. The CRM decisive test with subgenome-specific clades

**Evidence** (proxy only: CRM counts, no truth). A CRM marker is decisive only if its cores are significant on a
majority of chromosomes (annotation sensitivity matters: [design.md](design.md), Known risks). If the subgenomes carry
different CRM clades, a marker can fail that majority although it is centromeric in each. The CRM count differs
between subgenomes A and B in 3/87 species at P <= 0.01 (Fisher; at most 0.9 expected by chance, P = 0.06 for three
or more) and in 9 at P <= 0.05 (at most 4.4 expected, P = 0.03 for nine or more), so imbalance may be more common
than chance. But the marker is decisive in most of them anyway: 2 of the 3 (*Brassica napus*, A 8/8 vs B 2/8;
*Euphrasia confusa*, 10/11 vs 3/11) and 5 of the 9. Only the others can have lost a decisive marker this way: one of
the 3 (*Cerastium alpinum*, flags-only, 11/18 vs 0/18) and 4 of the 9. So there is no evidence yet that a decisive
marker is often lost this way.

**Approach.** Test CRM concentration within the chromosome subsets of entry 3 and pick the marker per subset.

**Check.** The three species against CENH3 truth (9), with the diploid null of entry 3. Not worth the code before truth
shows that a decisive marker is lost this way.

## 5. A faster search

**Evidence.** The genome-wide search is the longest stage and hard to predict (*Panax japonicus*: 7.9 h predicted,
more than 25 h actual). The cost sits in few families: the 172-bp alpha satellite took 62% of the CPU in CHM13, and in
plants and chicken 1-2-kb families do. Since 0.9, `satellite/work/annotate/<genome>.search_cost.tsv` names them for
every genome.

**Approach.** For long units, a BLASTN or minimap2 prefilter to find candidate windows, with nhmmer only on them. For
dense short-unit arrays (millions of hits per family in mammals), a search that treats the seed arrays apart from the
rest of the genome. Both are to be tested.

**Check.** The monomer and array tables must equal the stored annotations of the benchmark runs and the replay set, or
every difference must be listed hit by hit; the gain is read from the cost table (CPU per family before and after).

## 6. Seed clustering

**Evidence.** Seed clustering (`cluster_seeds` and `absorb_multimers`) took 76 min for *Acanthus tetraploideus* (5,456
seeds, 16 threads) and 2.3 min for *Hibiscus tridactylites* (7,631 seeds); the median over the 106 genomes is 0.7 min,
and two take >= 10 min. BLAST is 97.7% of the wall time, driven by the HSP volume of multimeric seeds (14.4 M and 4.6 M
HSPs in *Acanthus*); Python is about 2% (the stored hits replayed through the unmodified Python took 34 + 20 s). No
speedup that keeps the output exact was found, so 0.9 only logs each round at INFO.

**Approach.** Align only the pairs `seed_match` can consult: length-banded, N-padded databases with `-dbsize` pinned,
queries as single-thread jobs. Prototyped outside the repository: 10x on *Acanthus* (266 + 69 s against 2,223 + 1,123 s
profiled), identical clusters on 8 genomes. It is not exact by construction: BLAST's hit lists differ for 0.2-0.5% of
the consulted pairs (borderline hits near E 1e-5 depend on the batch of queries searched together). Other variants:
single-thread query jobs on the same database, 1.6-1.9x per call, hits differ (3 of 873k); `-mt_mode`, hits identical in
12 of 13 reruns, times mixed.

**Check.** Before it ships: a fixed, thread-independent partition of the queries, and identical cluster assignments and
representatives (after absorption) on every stored seed set: the benchmark runs and the 106 stored runs.

## 7. The FasTAN upstream report

**Evidence.** A 642.5-kb window of the public *Oryza minuta* assembly, GWHESQK00000003:29,940,000-30,582,500 (no N, no
soft-masking, about 90% one tandem satellite array with 15-mer spacings of 243 and 342 bp), aborts FasTAN
deterministically: `malloc(): corrupted top size`, SIGABRT; neither half of it crashes alone. A second genome
(*Sonchus oleraceus*) crashed at another site. Both crash sites lie inside a high-confidence centromere call, so
assemblies with complete centromeres are hit hardest. 0.9 recovers by splitting the piece, and stops when FasTAN
crashes on every piece of a chunk ([design.md](design.md)); the report to FasTAN's authors has not been filed.

**Approach.** Report it upstream with the reproducer (that window: `FAtoGDB`, then `FasTAN -m -p -T1`) and cenfuse's
crash table as the second example. Outward-facing, so the user decides.

**Check.** With a fixed FasTAN release, rerun the reproducer and the 5-Mb *O. minuta* piece (no crash, no crash table).
The splitting stays as a safety net for other crashes.

## 8. Large unplaced scaffolds in the chromosome set

**Evidence.** On the 121 polyploid assemblies that come with a list of their chromosomes, the 0.9 rule leaves 31
sequences outside the list in the set (5,614 in 0.8), in 12 assemblies (7 in one, 5 in another), 1.5-117 Mb, median
7 Mb. None of the listed chromosomes is lost. A length rule cannot tell such a sequence from a short chromosome; a few
may be real chromosomes that the list left out. Extra sequences dilute the per-chromosome prevalence tests.

**Approach.** A second signal for the sequences at the edge of the set: telomere arrays at the ends, a centromere-like
tandem array, gene density, Hi-C contacts. Use it to flag them in the report's Notes, not to drop them.

**Check.** The 121 chromosome lists as truth for the flag's precision and recall. The 234 benchmark and batch assemblies,
where the rule changes nothing, as a control: the flag should be rare there.

## 9. CENH3 truths for allopolyploids

**Evidence.** CENH3 truth for allopolyploids is scarce: *Ae. geniculata* (14 chromosomes) is the one allotetraploid with
CENH3 truth in the benchmark's satellite batch. The calls of the 106-genome set are otherwise checked by proxies, and
the proxies cannot arbitrate the question that matters most: where a CRM core replaced a satellite array (143 of 2,212
audited calls), homoeolog agreement favors the core (26 vs 14 pairs, P = 0.08) and gene density the replaced array (80
vs 34, P = 2e-5; an array holds no genes by construction, so that proxy is biased).

**Approach.** CENH3 ChIP-seq or CUT&Tag for a few allotetraploids chosen to settle open entries: calls where a core
replaced an array (the `seq_overruled_*` columns list them), subgenome-specific families (3, 4), weak main families (1).

**Check.** Score as in the benchmark (a call is right within 4 Hi-C bins of the CENH3 domain). No rerun is needed.

## 10. Chromosome-set limits

**Evidence.** The rule cannot tell a short chromosome from a scaffold by length ([README](../README.md#inputs)). A dot
chromosome more than tenfold shorter than every other chromosome, with nothing between, is dropped (as in 0.8):
*D. melanogaster* X, 2L, 2R, 3L, 3R and the 1.35-Mb chromosome 4, without the Y, keeps 5 of 6. The tail cut of 0.9
fails where 0.8 kept the chromosomes in two cases: dot chromosomes followed by a gradual scaffold tail are cut with the
tail (a bird-like karyotype of 37 chromosomes plus 300 scaffolds tapering from 0.4 Mb: 0.8 keeps all 337 sequences,
0.9 keeps 31 and drops six 0.5-1.4-Mb dot chromosomes; with a clean gap before the scaffolds both keep 37), and small
chromosomes that themselves span a tenfold range are cut as a tail. The workaround for each is
`<genome>.chromosomes.txt`; the run's notes say when the tail cut fired.

**Approach.** The signal of entry 8, applied from the other side: a sequence below the cut with telomere arrays at both
ends and a tandem array is a chromosome candidate, to be named in the report.

**Check.** Karyotypes with known dot chromosomes (*D. melanogaster*, the bird microchromosomes), and the 383 assemblies
of the chromosome-set evidence as a no-change control.

## 11. The LTR retrotransposon line

A check of the LTR line on stored results (17 genomes with `tes/` truth, 190 chromosomes; the 106-genome set) points to
four changes worth testing: 11a-11c below and the decisive test within subgenomes (entry 4). Caveat: the 17 genomes are
a development set, not a test set: all 20 gains are on genomes the rules were built or developed on, and the rest only
check that nothing breaks. So each of the four needs a `tes/` truth set made after the rules are fixed (new CENH3 data
with an LTRquest annotation, or genomes never cited as rule cases): the rule changes (11b, 11c, 4) before they can be
adopted, and 11a, which changes no call, before the gains its per-rule switch attributes can be credited to a rule.

**11a. Show which families the CRM line rejected.** 30 calls moved with `tes/`; 19 of them had a family in the
FASTA-only arm that the CRM line had rejected at classify (13 labeled family-level, and 6 of the 8 moves that carry a
replacement flag). The rejection leaves no flag, so the share of calls the TE line changed (at least 22%: 496 of the
2,212 calls audited in 105 of the 106 species) is a floor, and no gain can be credited to a flagged rule alone: all 5
flagged gains also had a rejected family. The family-level rule also left two `high` wrong calls. Potato chr04 has no TE
support (no CRM core on the chromosome): its competing cluster became the call and, with no rival left, went from `low`
to `high`. *Ae. tauschii* Chr5 moved to another family's array that overlaps the edge of a 7-Mb core, so it is
`te_concordant` (concordant by that overlap alone), which raised it from `low` to `high`; CENH3 lies 1-3 Mb away.
*Approach:* write the rejected families per call (a `te_rejected_family` flag or a caveat); output only. *Check:* rerun
the centromeres stage on the stored satellite folders and compare tables (no call changes); a diagnostic switch that
turns each TE rule off would let the gains be attributed.

**11b. The head-to-head bar for `low` calls.** A decisive core replaces a disagreeing call that is `low` on its own
(`satellite_away_from_core`) without the test a `high` or `medium` call has to pass (the clade marks the core
significantly more, P <= 0.01). That branch fired in 3 truth cases: potato chr10 (core right, array wrong: +1), soybean
chr13 (array right, core wrong: -1; 13 vs 7 elements, P = 0.13, so the bar would have kept the array) and *Ae. tauschii*
Chr6 (right in both arms: 0). Two of them (chr10, Chr6) also had a rejected family, so only chr13 is purely call-level.
The replaced call's confidence and the element counts are not stored, so the 83 audited `satellite_away_from_core` calls
of the 106-genome run cannot be split into tested and untested ones. *Approach:* one standard for replacing a call
whatever its confidence: the same bar for `low` calls. *Check:* the `seq_overruled_*` columns on a replay of the
106-genome run, plus new CRM-plateau genomes with truth; at least 10 cases, since 3 cannot decide.

**11c. A decisive `te_discordant` call with an inconclusive test.** When a decisive core disagrees with a `high` or
`medium` call and the head-to-head test is not significant (equal marking, or too few elements), the call stays, flagged
`te_discordant` and one level lower. In the truth set that lowered 6 calls from `high` to `medium`, and all 6 were right
in both arms: an under-confidence cost. They are 3 loci in 2 species (*B. rapa* A09, one locus in four accessions;
soybean chr1 and chr6), all of them conflict cases the head-to-head rule was built on. Of the 2,212 calls audited in
the 106-genome run, 16 were lowered this way (12 agree with the homoeolog's call). *Approach:* keep the flag and the
level when a decisive marker's test is inconclusive: a non-significant test is absence of evidence, and a clade that
marks the call as much as the core supports it. The counter-argument: two equally marked loci are an unresolved
conflict. A flags-only marker never lowers a call; its `te_discordant` calls (right 1/16) stay a flag, and lowering
them was considered and not proposed (see below). *Check:* store the element counts (`k_core`, `k_call`) per call
first, then new genomes with a broad CRM plateau and truth.

Considered, not proposed: lowering confidence under a flags-only `te_discordant` (15/16 wrong in moss, tomato and
Chlorella, but pericentromeric non-CRM clades are common, rice Tekay cores lay off CENH3, and the proxies favor the
satellite call); raising a flags-only core-only call above `low` (Chlorella 8/8 right, tomato 0/7).
