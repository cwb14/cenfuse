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
| 16 | [the CHG line has changed no call](#16-the-chg-line-has-changed-no-call) | CHG from the assembled plant, with CENH3 truth |
| 17 | [knobs kept out only by a size floor](#17-knobs-kept-out-only-by-a-size-floor) | knob-rich lines with CENH3 truth |
| 18 | [organelle insertions called 45S arrays](#18-organelle-insertions-called-45s-arrays) | a rerun of the rDNA stage |
| 19 | [long-unit calls are mostly wrong](#19-long-unit-calls-are-mostly-wrong) | a rerun of the centromeres stage |
| 20 | [the rDNA spacer tag at its threshold](#20-the-rdna-spacer-tag-at-its-threshold) | *B. rapa* accessions run together |
| 21 | [moved calls whose sentences disagree](#21-moved-calls-whose-sentences-disagree) | a rerun of the report |
| 22 | [Pa1A in other *P. annua* plants](#22-pa1a-in-other-p-annua-plants) | CENH3 or deeper reads from those plants |
| 23 | [a family's unit length can overstate its copies](#23-a-familys-unit-length-can-overstate-its-copies) | output only |
| 24 | [telomeres outside the rules](#24-telomeres-outside-the-rules) | T2T assemblies as a no-change control |

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

## 16. The CHG line has changed no call

**Evidence.** Only *P. annua* and *Arabidopsis* have CHG data here, and across all 3,595 chromosomes of the release
runs the line supported no call and backed no move. *P. annua*'s five samples come from other plants, so they do not
decide; *Arabidopsis* learned a signature on 1 of 5 chromosomes. The signal is real in *P. annua*: the dip is the
chromosome's lowest point on 12 of 14.
**Approach.** None until the line has been seen to act.
**Check.** CHG from the assembled plant of a genome with CENH3 truth.

## 17. Knobs kept out only by a size floor

**Evidence.** Nothing recognizes a knob. In maize Mo17 the 180-bp knob family misses the centromeric role only because
its median largest block (74 kb) is under the 100-kb floor (`strong_cluster_bp`). At a 70-kb floor, FASTA only, it takes
chr4, chr6 and chr8 (none with `tes/`). `classify.related` is one-way by unit length (155 vs 180 bp passes, 180 vs 155
fails), so CentC's chromosomes count for the knob.
**Approach.** Relatedness by sequence only (`same_sequence`): in a scratch rerun at the 70-kb floor CentC takes all ten.
**Check.** Knob-rich lines with CENH3 truth; no call changes elsewhere.

## 18. Organelle insertions called 45S arrays

**Evidence.** A 45S array is 3 or more motif hits, each within 100 kb of the next; neither spacing nor strand is
tested, and organelle rRNA holds the motif too. Rice chr4 9.23-9.37 Mb, reported as a NOR, is a chloroplast insertion
(4 hits, mixed strands, 99.8% identical to the plastid genome). By the base after the motif (G in 18S, A in organelle
rRNA), 12 of 272 arrays look organelle-like. A tested Hi-C peak beside one is flagged `rDNA` and can no longer oppose or
place a call; the effect on calls is not measured.
**Approach.** Require one strand, or test that base. Even spacing is too strict: only 100 of 272 arrays are regular.
**Check.** Those 12 go; every other array stays.

## 19. Long-unit calls are mostly wrong

**Evidence.** A long-unit array (units over 2 kb) is the call on 24 of 908 chromosomes FASTA only: 2 right, 16 wrong (6
of them holocentric), 6 without truth, all `low`.
**Approach.** Report such arrays as features unless another line agrees.
**Check.** Right and wrong counts after a rerun of the centromeres stage.

## 20. The rDNA spacer tag at its threshold

**Evidence.** A family is tagged as an rDNA spacer when 45S genes flank at least half of its 10 largest arrays. *B.
rapa*'s 238-bp family reaches 5-6 of 10 in 5 accessions (tagged) and 3-4 in 6 (not). In CX1, at 4 of 10, it stays
eligible and leads A05 (wrong, `low`) and A06 (right).
**Approach.** Run the accessions together: a family tagged in one genome of a complex is tagged in all.
**Check.** The 11 accessions against CENH3.

## 21. Moved calls whose sentences disagree

**Evidence.** On moss's 4 moved calls `basis` says, for example, "21 kb of centromeric satellite inside", while `caveat`
says "no sequence or LTR evidence at this place": the caveat does not read the evidence the basis reports.
**Approach.** Word the caveat from the same evidence (`interpret.py`).
**Check.** A unit test with a moved call on a small array; no call changes.

## 22. Pa1A in other P. annua plants

**Evidence.** Young CRM insertions of four of the five other plants cluster 9-15 Mb from the assembled plant's Pa1A
call (163.0-169.4 Mb), which suggests other centromeres; J10's lie inside it. C11 and J10 show a plant-only methylation
dip at their clusters. No verdict confirms it: on Pa1A 4 `shared` and 1 `unclear`; over all chromosomes 0 `shifted` and
2 `possible_shift` (J10 Pa4B, Pa5A) of 70, at 4-8× reads.
**Approach.** None in the code until the data can settle it.
**Check.** CENH3 ChIP or deeper long reads from those plants.

## 23. A family's unit length can overstate its copies

**Evidence.** In rice, 19 of the 28 large arrays (≥ 20 kb) of the 165-bp CentO family are not spaced 165 bp, 17 of them
about 155 bp: the 165-bp profile fits 155-bp copies with a 10-bp gap and scores higher. The forms are pooled, so calls
are unaffected (12/12), but `basis` reports a 165-bp unit.
**Approach.** Report the array's own spacing beside the family's unit.
**Check.** Output only; no call changes.

## 24. Telomeres outside the rules

**Evidence.** Units of 13-29 bp fall between the telomere search (5-12 bp) and the satellite search (30 bp and up): *C.
tropicalis*'s 23-bp unit, at 13 of 14 ends, is missed. Of 208 ends without an array, 28 hold a real stretch under half
unit copies, which the density rule rejects (soybean chr3: 48 kb, 20.2 kb of copies).
**Approach.** Extend the unit range; relax the density rule for stretches that reach the end.
**Check.** T2T assemblies unchanged; those ends gain arrays, and no new interstitial array appears.
