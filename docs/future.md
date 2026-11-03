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
| 11 | [the LTR retrotransposon line](#11-the-ltr-retrotransposon-line) | held-out truth |

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

**Evidence.** At most `max_families` = 40 clusters per complex are modelled and searched, since annotate time grows with
models x genome size. In *Sonchus oleraceus*, 60 clusters of >= 20 kb were dropped, the largest 107 kb (the log reports
the count and the largest cluster dropped in every genome). What the dropped clusters are, and whether any is a
centromeric satellite, is not recorded.

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

**Evidence.** A CRM marker is decisive only if its cores are significant on a majority of chromosomes (annotation
sensitivity matters: [design.md](design.md), Known risks). If the subgenomes carry different CRM clades, a marker can
fail that majority although it is centromeric in each. The CRM count differs between subgenomes A and B in 3/87
species at P <= 0.01 (Fisher): at most 0.9 expected by chance (P = 0.06 for three or more). So there is no evidence
yet that this happens often.

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

**Evidence.** A 642.5-kb window of *Oryza minuta* (no N, no soft-masking, about 90% one tandem satellite array with
15-mer spacings of 243 and 342 bp) aborts FasTAN deterministically: `malloc(): corrupted top size`, SIGABRT; neither
half of it crashes alone. A second genome (*Sonchus oleraceus*) crashed at another site. Both crash sites lie inside a
high-confidence centromere call, so assemblies with complete centromeres are hit hardest. 0.9 recovers by splitting the
piece ([design.md](design.md)); the report to FasTAN's authors has not been filed.

**Approach.** Report it upstream with the reproducer (`FAtoGDB`, then `FasTAN -m -p -T1`) and cenfuse's crash table as
the second example. Outward-facing, so the user decides.

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
*D. melanogaster* X, 2L, 2R, 3L, 3R and the 1.35-Mb chromosome 4, without the Y, keeps 5 of 6. Tiny chromosomes inside
an assembly with a gradual scaffold tail are lost (0.8 took the whole tail), and small chromosomes that themselves span
a tenfold range are cut as a tail. The workaround is `<genome>.chromosomes.txt`.

**Approach.** The signal of entry 8, applied from the other side: a sequence below the cut with telomere arrays at both
ends and a tandem array is a chromosome candidate, to be named in the report.

**Check.** Karyotypes with known dot chromosomes (*D. melanogaster*, the bird microchromosomes), and the 383 assemblies
of the chromosome-set evidence as a no-change control.

## 11. The LTR retrotransposon line

A check of the LTR line on stored results (17 genomes with `tes/` truth, 190 chromosomes; the 106-genome set) points to
three changes worth testing. Caveat: all the gains were on genomes that were used to build the rules, so these are
leads, not a test.

**11a. Show which families the CRM line rejected.** 30 calls moved with `tes/`; 19 of them had a family in the FASTA-only
arm that the CRM line had rejected at classify. The rejection leaves no flag, so the share of calls the TE line changed
(22% of the 106-genome audit) is a floor, and two high-confidence wrong calls (potato chr04, *Ae. tauschii* Chr5) had no
TE support and no flag to say so. *Approach:* write the rejected families per call (a column, or in `caveat`). *Check:*
`classify` rerun on the truth genomes and the 106; the count must reproduce 19/30. No call changes.

**11b. The head-to-head bar for `low` calls.** A decisive core replaces a disagreeing call that is `low` on its own
without the bar a `high` or `medium` call has to pass (the clade marks the core significantly more, P <= 0.01). Applying
the bar to `low` calls would change 3 truth cases (+1, -1 and 0 right calls). *Approach:* apply the same bar.
*Check:* 3 cases cannot decide; it needs truth beyond the development genomes.

**11c. A decisive `te_discordant` call with an inconclusive test.** When a decisive core disagrees with a `high` or
`medium` call and the head-to-head test is inconclusive, the call stays, flagged `te_discordant` and one level lower. In
the truth set that lowered 6 calls, and all 6 were right (*B. rapa* A09 in four accessions, soybean chr1 and chr6): an
under-confidence cost, since an inconclusive test is no evidence against a call. *Approach:* keep the flag and the
confidence when the marker is decisive and the test is inconclusive; flags-only markers keep lowering (their
`te_discordant` calls were right 1/16). *Check:* those 6 are the chromosomes the rule was built on; it needs held-out
truth.
