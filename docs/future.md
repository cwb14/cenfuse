# Future development

Open problems in cenfuse 0.10.1. Each gives the evidence, a proposed approach and the check it needs (a rerun, or
new truth data). Nothing here is implemented. The evidence comes from the benchmark ([benchmark.md](benchmark.md)) and
from a run on 106 allotetraploid plant genomes (2,320 calls: 54% `high`, 23% `medium`, 23% `low`).

| # | open problem | needs |
|---|---|---|
| 1 | [weak (single-sequence) family models](#1-weak-single-sequence-family-models) | rerun from the families stage |
| 2 | [the 40-family cap](#2-the-40-family-cap) | rerun with a higher cap |
| 3 | [satellites specific to subgenomes vs the honesty cap](#3-satellites-specific-to-subgenomes-vs-the-honesty-cap) | CENH3 truth, a null |
| 4 | [the decisive CRM test with subgenome-specific lineages](#4-the-decisive-crm-test-with-subgenome-specific-lineages) | CENH3 truth, a null |
| 5 | [a faster search](#5-a-faster-search) | identical outputs on stored annotations |
| 6 | [faster seed clustering](#6-faster-seed-clustering) | identical clusters on every stored seed set |
| 7 | [the FasTAN crash, reported upstream](#7-the-fastan-crash-reported-upstream) | the user's decision, then a fixed release |
| 8 | [large unplaced scaffolds and dot chromosomes](#8-large-unplaced-scaffolds-and-dot-chromosomes) | known chromosome lists and karyotypes |
| 9 | [CENH3 truths for allopolyploids](#9-cenh3-truths-for-allopolyploids) | new data |
| 10 | [the LTR retrotransposon line](#10-the-ltr-retrotransposon-line) | held-out `tes/` truth |
| 11 | [Hi-C anchors where the sequence line is wrong genome-wide](#11-hi-c-anchors-where-the-sequence-line-is-wrong-genome-wide) | CENH3 genomes with Hi-C and `tes/` |
| 12 | [cis signals of satellite-free centromeres](#12-cis-signals-of-satellite-free-centromeres) | more satellite-free CENH3 genomes with Hi-C |
| 13 | [bird microchromosomes in Hi-C](#13-bird-microchromosomes-in-hi-c) | a sparse fold scan; deep bird Hi-C |
| 14 | [the Rabl axis where the fold has no power](#14-the-rabl-axis-where-the-fold-has-no-power) | more small genomes with truth and Hi-C |
| 15 | [centromeres at the edge of a terminal unmappable array](#15-centromeres-at-the-edge-of-a-terminal-unmappable-array) | acrocentric genomes with truth |

## 1. Weak (single-sequence) family models

**Evidence.** A family with fewer than 20 sampled monomers gets a single-sequence model, with no profile of how its
copies vary: 214 such models in 68 of 105 allotetraploids, 11 of them the main family of 131 calls. How far their
arrays and calls are affected is not measured.
**Approach.** Build the profile from more monomers (every seed of the cluster, or the same family in the other
genomes of a complex); keep the single model where there are still too few.
**Check.** Rerun those genomes from the families stage and count how many calls move, and which way; CENH3 truth (9).

## 2. The 40-family cap

**Evidence.** At most 40 clusters per complex are modeled, since the search time grows with models × genome size. In
*Sonchus oleraceus* 60 clusters of ≥ 20 kb were dropped, the largest 107 kb; what they are is not recorded.
**Approach.** A cost budget instead of a count (model length × genome size, priced from `search_cost.tsv`), and a
table of the dropped clusters.
**Check.** Rerun the genomes with dropped clusters ≥ 20 kb and compare roles, main families and calls.

## 3. Satellites specific to subgenomes vs the honesty cap

**Evidence.** When no family is the largest centromeric cluster on most chromosomes, every call is capped `low`. In
an allopolyploid each subgenome can have its own satellite: in 30 of 105 species the top families differ between
subgenomes, and 85 calls were capped.
**Approach.** Judge dominance within chromosome groups when the user gives them (a table like
`<genome>.chromosomes.txt`); lift the cap for a group only where one family dominates inside it.
**Check.** CENH3 truth for allopolyploids (9); the null: random splits of a diploid's chromosomes into groups.

## 4. The decisive CRM test with subgenome-specific lineages

**Evidence** (CRM counts only, no truth). CRM counts differ between subgenomes in 3 of 87 species at P ≤ 0.01 (0.9
expected), but the marker is decisive in 2 of the 3 anyway: no evidence yet that a decisive marker is often lost.
**Approach.** Test CRM concentration within the chromosome groups of entry 3.
**Check.** Those species against CENH3 truth (9), with the diploid null of entry 3.

## 5. A faster search

**Evidence.** The genome-wide nhmmer search is the longest stage and hard to predict (*Panax japonicus*: 7.9 h
predicted, more than 25 h actual). A few families carry the cost: alpha satellite took 62% of CHM13's CPU; in plants
and chicken 1-2-kb families do (`satellite/work/annotate/<genome>.search_cost.tsv`).
**Approach.** A BLASTN or minimap2 prefilter for long units, with nhmmer only on candidate windows; for dense
short-unit arrays, a search that treats the seed arrays apart.
**Check.** Monomer and array tables equal to the stored annotations, or every difference listed.

## 6. Faster seed clustering

**Evidence.** BLAST is 98% of seed clustering, driven by multimeric seeds: *Acanthus* took 76 min at 16 threads
(median over 106 genomes: 0.7 min). Aligning only the pairs the clustering consults was 10× faster with identical
clusters on 8 genomes, but BLAST's hit lists differ for 0.2-0.5% of those pairs, so it is not exact by construction.
**Approach.** A fixed, thread-independent partition of the queries.
**Check.** Identical cluster assignments and representatives on every stored seed set.

## 7. The FasTAN crash, reported upstream

**Evidence.** A 642.5-kb window of *Oryza minuta* (GWHESQK00000003:29,940,000-30,582,500, about 90% one satellite
array) aborts FasTAN every time (`malloc(): corrupted top size`); *Sonchus oleraceus* crashed elsewhere. Both sites lie
inside confident centromere calls. cenfuse recovers by splitting the piece; the report has not been filed.
**Approach.** Report it with the reproducer (`FAtoGDB`, then `FasTAN -m -p -T1` on that window). Outward-facing, so
the user decides.
**Check.** A fixed FasTAN release: no crash on the reproducer and the 5-Mb piece.

## 8. Large unplaced scaffolds and dot chromosomes

**Evidence.** On 121 polyploid assemblies with known chromosome lists, the length rule keeps 31 sequences outside the
lists (1.5-117 Mb) in 12 assemblies; it drops a dot chromosome that stands alone after a tenfold gap (*D.
melanogaster* chromosome 4 without the Y) and dot chromosomes followed by a tapering scaffold tail.
**Approach.** A second signal for the sequences at the edge of the set (telomere arrays at both ends, a centromere-like
array, genes, Hi-C), to flag candidates in the report, not to drop or add them.
**Check.** The 121 lists as truth; *D. melanogaster* and bird karyotypes; the 383 assemblies of the rule's evidence
as a no-change control.

## 9. CENH3 truths for allopolyploids

**Evidence.** The benchmark's only allotetraploid with CENH3 truth is *Ae. geniculata* (14 chromosomes). Where a CRM
core replaced a satellite array (143 of 2,212 audited allotetraploid calls), homoeolog agreement favors the core (26
vs 14 pairs, P = 0.08) and gene density the array (a biased proxy).
**Approach.** CENH3 ChIP-seq or CUT&Tag for a few allotetraploids chosen to settle entries 1, 3, 4 and 10
(`seq_overruled_*` lists the replaced arrays).
**Check.** Score as in the benchmark; no rerun needed.

## 10. The LTR retrotransposon line

A check on 17 genomes with `tes/` and truth (190 chromosomes) points to three changes worth testing. These genomes are
a development set: each change needs a `tes/` truth set made after the rules are fixed.

- **Show which families the CRM line rejected.** 19 of 30 calls that moved with `tes/` had a family the CRM line
  rejected at classification, which leaves no flag. *Approach:* a flag or caveat naming it (output only). *Check:* a
  rerun of the centromeres stage changes no call.
- **One standard for replacing a `low` call.** A CRM core replaces a disagreeing `low` call without the head-to-head
  test a confident call must pass; in 3 truth cases that gave +1, -1 and 0. *Approach:* the same test for `low` calls.
  *Check:* at least 10 such cases with truth.
- **An inconclusive head-to-head test.** When a CRM core disagrees with a confident call and the test is not
  significant, the call drops one level; in the truth set all 6 such calls were right (*B. rapa* A09, soybean chr1 and
  chr6). *Approach:* keep the level when the test is inconclusive. *Check:* store the element counts per call, then
  new genomes with a broad CRM plateau and truth.

## 11. Hi-C anchors where the sequence line is wrong genome-wide

**Evidence.** Potato and barley Hi-C anchored at the true centromeres find 7/11 and 6/7 of them, but anchored at the
(wrong) sequence calls none, so the test fails; prior-free clustering does not help in potato.
**Approach.** Anchor and test Hi-C against another independent line's calls (decisive CRM cores, methylation
candidates).
**Check.** Potato and barley with `tes/` and Hi-C; no change where the satellite line is right.

## 12. Cis signals of satellite-free centromeres

**Evidence.** In tomato, the contact insulation within 1 Mb peaks on CENH3 on 9/12 chromosomes (chance 0.4), where
every sequence line fails; potato shows the opposite sign (6/11), mammals nothing, and no per-genome test is specific.
**Approach.** Test a cis profile like the anchored Hi-C, against an independent line, with its sign chosen
leave-one-out.
**Check.** More satellite-free CENH3 genomes with Hi-C (Solanaceae, legumes).

## 13. Bird microchromosomes in Hi-C

**Evidence.** Hi-C bins finer than 250 kb are capped by the largest chromosome (4,000 bins): chicken's 200-Mb chr1
needs 4,001 at 50 kb, so its 2.5-7-Mb microchromosomes stay out of the Hi-C line.
**Approach.** A sparse or banded fold scan, which lifts the per-chromosome cap.
**Check.** Chicken with deeper Hi-C; calls on the microchromosomes against the CENP-A windows.

## 14. The Rabl axis where the fold has no power

**Evidence.** The trans arm alignment passes in *C. neoformans* (anaphase, 12/14) and *M. sympodialis* (6/8), and
Hi-C's prior-free calls there are right on 10/10 and 6/6, but the arm fold, which must lean the same way before the
alignment may act (the holocentric guard), has no power in these small nuclei (1/14; 0/8), so the calls stay `low`.
In holocentric *C. elegans* the alignment reached 2/5 (P 0.095), and a position-matched null that rejects it also
rejects barley. In G1 *C. neoformans* the prior-free peaks sit at the chromosome ends, so the alignment's partners
cannot be placed (0/14), although anchored on CENP-A 8/14 chromosomes still align.
**Approach.** A compartment-adjusted alignment (trans contact explained by compartment identity removed) that can
stand alone; or a cis readout with power in small chromosomes; or partners placed on the sequence line's confident
calls, as a corroboration test.
**Check.** *Candida*, *Cryptococcus* relatives, *Toxoplasma*, *T. brucei* with truth and Hi-C; small-chromosome
holocentrics (nematodes, Lepidoptera) and randomized Hi-C must never pass.

## 15. Centromeres at the edge of a terminal unmappable array

**Evidence.** Chicken's acrocentric microchromosomes carry their centromeric array at the very end, unmappable for
1-1.25 Mb: their Hi-C peaks lie at its inner edge, which in mappable bins is the end, so they stay `telomeric`, and
chicken's calls stay `low`.
**Approach.** Where the sequence line finds an acrocentric karyotype, compare the two ends of each chromosome:
telomere clustering lifts both, a centromere one. It needs a per-genome test before it may count.
**Check.** Chicken with CENP-A truth (25/33 centromeres within 8 bins of an end), goat, sheep, horse; randomized Hi-C.
