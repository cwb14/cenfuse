# Future development

Open problems in cenfuse 0.11. Each gives the evidence, a proposed approach and the check it needs (a rerun, or new truth
data). Entries closed in 0.11 keep their numbers and say what was done; the rest is not implemented. The evidence comes from
the benchmark ([benchmark.md](benchmark.md)) and from a replay of 106 allotetraploid plant genomes (2,355 calls: 44% `high`, 30%
`medium`, 24% `low`).

| # | problem | needs |
|---|---|---|
| 1 | [calls on chains of partial copies](#1-calls-on-chains-of-partial-copies) | CENH3 truth where a chain dominates; a test of replacing the call by its core |
| 2 | [the 40-family cap](#2-the-40-family-cap) | a genome where it costs a call |
| 3 | [satellites specific to subgenomes vs the honesty cap](#3-satellites-specific-to-subgenomes-vs-the-honesty-cap) | tested in 0.11: no gain over random splits |
| 4 | [the decisive CRM test with subgenome-specific lineages](#4-the-decisive-crm-test-with-subgenome-specific-lineages) | CENH3 truth, a null |
| 5 | [a faster search](#5-a-faster-search) | low-complexity cost done in 0.11; the whole benchmark searched again; the cost of real hits |
| 6 | [faster seed clustering](#6-faster-seed-clustering) | identical clusters on every stored seed set |
| 7 | [the FasTAN crash, reported upstream](#7-the-fastan-crash-reported-upstream) | the user's decision, then a fixed release |
| 8 | [large unplaced scaffolds and dot chromosomes](#8-large-unplaced-scaffolds-and-dot-chromosomes) | known chromosome lists and karyotypes |
| 9 | [CENH3 truths for allopolyploids](#9-cenh3-truths-for-allopolyploids) | new data |
| 10 | [the LTR retrotransposon line](#10-the-ltr-retrotransposon-line) | held-out `tes/` truth (first bullet done in 0.11) |
| 11 | [Hi-C anchors where the sequence line is wrong genome-wide](#11-hi-c-anchors-where-the-sequence-line-is-wrong-genome-wide) | CENH3 genomes with Hi-C and `tes/` |
| 12 | [cis signals of satellite-free centromeres](#12-cis-signals-of-satellite-free-centromeres) | more satellite-free CENH3 genomes with Hi-C |
| 13 | [bird microchromosomes in Hi-C](#13-bird-microchromosomes-in-hi-c) | a sparse fold scan; deep bird Hi-C |
| 14 | [the Rabl axis where the fold has no power](#14-the-rabl-axis-where-the-fold-has-no-power) | more small genomes with truth and Hi-C |
| 15 | [centromeres at the edge of a terminal unmappable array](#15-centromeres-at-the-edge-of-a-terminal-unmappable-array) | acrocentric genomes with truth |
| 16 | [the CHG line has changed no call](#16-the-chg-line-has-changed-no-call) | CHG from the assembled plant, with CENH3 truth |
| 17 | [knobs and the size floor](#17-knobs-and-the-size-floor) | done in 0.11 for maize; other knob-rich lines with CENH3 truth |
| 18 | [organelle insertions called 45S arrays](#18-organelle-insertions-called-45s-arrays) | done in 0.11 |
| 19 | [long-unit calls are mostly wrong](#19-long-unit-calls-are-mostly-wrong) | done in 0.11 |
| 20 | [the rDNA spacer tag at its threshold](#20-the-rdna-spacer-tag-at-its-threshold) | *B. rapa* accessions run together |
| 21 | [moved calls whose sentences disagree](#21-moved-calls-whose-sentences-disagree) | done in 0.11 |
| 22 | [Pa1A in other *P. annua* plants](#22-pa1a-in-other-p-annua-plants) | CENH3 or deeper reads from those plants |
| 23 | [a family's unit length can overstate its copies](#23-a-familys-unit-length-can-overstate-its-copies) | done in 0.11 |
| 24 | [telomeres outside the rules](#24-telomeres-outside-the-rules) | done in 0.11 |

## 1. Calls on chains of partial copies

**Evidence.** A family can concentrate in one place without being a tandem satellite: under half of its bp in tandem and under
half in whole-unit hits, so it is a chain of pieces of a longer repeat, usually a retrotransposon. *Glycine dolichocarpa*'s 1,790-bp
family has 44 whole hits of 25,732, 65% inside SIRE elements, and calls 10 Mb wide. In the 106 allotetraploids, 41 such families in
29 species carry 409 of 2,355 calls, 209 of them `high` before 0.11. This was entry 1's "131 calls on single-sequence models" under
another name: those models are nearly all at least 1 kb, and a better profile would not make them tandem. Since 0.11 such a call
is at most `medium`, and `low` without an overlapping CRM core; no call moved. The new search changes these counts little: over its
14 reruns, 75 chain calls become 74 (entry 5).
That rule rests on principle (only CRM decides among retrotransposons), not on truth. No benchmark genome has a chain as its
dominant family; chain calls in the benchmark, all `low` or `medium`, are 15 right and 3 wrong. Proxies, which cannot settle it:
calls on CRM chains lie between TE cores and real satellites in gene density (0.135, 0.054, 0.009), and match their homoeolog's
place within 5% of the chromosome about as often (0.60, 0.70, 0.66).
**Approach.** Where a CRM core overlaps a chain call, make the core the call: such cores are tighter and more gene-poor than the
calls they overlap (gene density 0.021 vs 0.053), with about the same homoeolog agreement.
**Check.** CENH3 for a genome where a chain dominates (6 cotton species have one on 23-26 of 26 chromosomes); entry 9.

## 2. The 40-family cap

**Evidence.** At most 40 clusters per complex are modeled, since the search time grows with models × genome size; since 0.11 the
dropped ones are listed (`satellite/dropped_clusters.tsv`). Over 71 runs that dropped any: 2,376 clusters, none on half of the
chromosomes, 3 on three or more, 1,553 on one; 760 hold under 20 kb per chromosome. Only 4 reach 20 kb on a chromosome without a
call, all in *B. carinata*; 349 lie on a chromosome with a `low` call. So far nothing suggests the cap costs a centromere.
**Approach.** A cost budget instead of a count (model length × genome size, priced from `search_cost.tsv`), if a genome ever shows
the cap costing a call.
**Check.** Rerun *B. carinata* without the cap and compare roles, main families and calls.

## 3. Satellites specific to subgenomes vs the honesty cap

**Tested in 0.11, not implemented.** When no family is the largest centromeric cluster on most chromosomes, every call is capped
`low`, and in an allopolyploid each subgenome can have its own satellite. Judging dominance within the subgenomes of 20 capped
polyploids lifts 27 of their 292 capped calls, but random splits of the same chromosomes lift 9.3; one species of 20 beats its own
random splits (*O. schlechteri*, P 0.005). Random splits of diploids lift up to 7.9 calls (siamang), 5.7 (pig) and 4.7 (chicken), so
a lift this size is what splitting does anyway.
**Check, to reopen it.** CENH3 truth for an allopolyploid with subgenome-specific satellites (9).

## 4. The decisive CRM test with subgenome-specific lineages

**Evidence** (CRM counts only, no truth). CRM counts differ between subgenomes in 3 of 87 species at P ≤ 0.01 (0.9
expected), but the marker is decisive in 2 of the 3 anyway: no evidence yet that a decisive marker is often lost.
**Approach.** Test CRM concentration within subgenomes.
**Check.** Those species against CENH3 truth (9), with random splits of diploids as the null, as in entry 3.

## 5. A faster search

**Done in 0.11: the cost of low-complexity sequence.** A long model with a low-complexity tract spent 25-50 s per 2-Mb chunk with no
hit. Long models now search with those tracts neutral, 1.00-4.93x faster and never slower
([benchmark.md](benchmark.md#10-speed-memory-and-reproducibility); [design.md](design.md), search cost).
**Open: what that search does to results.** On 14 full reruns, none of the 143 chromosomes with truth changed outcome or confidence, but
only 8 of those genomes have truth, and the rest of the benchmark has not been searched again: its release rerun read stored searches.
Without truth, overlapping families trade the lead in *Cuscuta* and *Panax*: 19 calls change main family and 2 change place, one of them
`high` -> `low`. A neutral model scores less, so a close family without a tract can win their overlaps (cotton: 906 whole copies of a
512-bp family go to its 481-bp relative, which loses an rDNA tag). E-value parameters stay those of the whole model. *Approach:* compare
overlapping hits only on positions that can score in both models. *Check:* every benchmark arm rerun from the search stage, with a split
copy counted as lost; CENH3 in *Cuscuta* or *Panax* (9).
**Open: the cost of real hits.** Chain families hit nearly everywhere (*Glycine dolichocarpa* gains 1.02x, *Cuscuta* 1.15x), and dense
short-unit arrays hold the rest (alpha satellite: 62% of CHM13's CPU). `search_cost.tsv` says which family costs what. *Approach:* for
dense short-unit arrays, a search that treats the seed arrays apart; for chain families, stopping once a family's role is clear from a
sample of chromosomes. *Check:* monomer and array tables equal to the stored annotations, or every difference listed.

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
**Approach.** CENH3 ChIP-seq or CUT&Tag for a few allotetraploids chosen to settle entries 1, 3, 4, 5 and 10
(`seq_overruled_*` lists the replaced arrays).
**Check.** Score as in the benchmark; no rerun needed.

## 10. The LTR retrotransposon line

A check on 17 genomes with `tes/` and truth (190 chromosomes) points to the first three changes, and the polyploid runs of 0.11 to the
fourth. These genomes are a development set: each change needs a `tes/` truth set made after the rules are fixed.

- **Show which families the CRM line rejected. Done in 0.11:** columns `seq_rejected_family`, `seq_rejected_start` and
  `seq_rejected_end`, and a caveat naming it. 19 of 30 calls that moved with `tes/` had such a family.
- **One standard for replacing a `low` call.** A CRM core replaces a disagreeing `low` call without the head-to-head
  test a confident call must pass; in 3 truth cases that gave +1, -1 and 0. *Approach:* the same test for `low` calls.
  *Check:* at least 10 such cases with truth.
- **An inconclusive head-to-head test.** When a CRM core disagrees with a confident call and the test is not
  significant, the call drops one level; in the truth set all 6 such calls were right (*B. rapa* A09, soybean chr1 and
  chr6). *Approach:* keep the level when the test is inconclusive. *Check:* new genomes with a broad CRM plateau and truth
  (since 0.11 the element counts are stored per call: `te_core_elements`, `te_call_elements`).
- **`te_concordant` is a weak flag.** Two estimates agree if the gap between them is at most the wider of the two (`te.agrees`), so a wide
  core or call agrees with much of its chromosome: a call spanning 23% of a chromosome agrees with about 70% of random cores, and a 340-kb
  call on *Panax*'s sequence 6 agrees with a core 4.3 Mb away. The chain cap asks for overlap instead. *Approach:* score agreement against
  random placements of the same widths on that chromosome. *Check:* whether the calls that keep the flag are right more often, on these
  genomes and then a held-out set.

## 11. Hi-C anchors where the sequence line is wrong genome-wide

**Evidence.** Potato and barley Hi-C anchored at the true centromeres find 7/11 and 6/7 of them, but anchored at the
(wrong) sequence calls none, so the test fails; prior-free clustering does not help in potato. Anchored at potato's calls from CRM
cores the test fails too (0/3-0/4), so those cores alone do not rescue it there.
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
Hi-C's prior-free calls there are right on 10/10 and 7/7, but the arm fold, which must lean the same way before the
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

## 17. Knobs and the size floor

**Done in 0.11 for maize; one genome.** Nothing recognizes a knob as such. Until 0.11, maize Mo17's 180-bp knob family missed the centromeric
role only because its median largest block (74 kb) is under the 100-kb floor (`strong_cluster_bp`): at a 70-kb floor, FASTA only, it took
chr4, chr6 and chr8. `classify.related` was one-way by unit length (155 vs 180 bp passed, 180 vs 155 failed), so CentC's chromosomes counted
for the knob. The test now works both ways: at floors of 100, 70 and 50 kb the knob stays `interstitial` and CentC keeps all ten, `high`
(replay; the 0.10.1 code at 70 kb gives the three knob calls). The same change put *B. rapa* CX1 right ([20](#20-the-rdna-spacer-tag-at-its-threshold)).
Relatedness by sequence alone, proposed here before, stays out: it changes the main family on 90 of 3,459 chromosomes and removes dominance
in 6 genomes with no truth behind it, and it works against satellites specific to subgenomes.
**Check.** Other knob-rich lines with CENH3 truth.

## 18. Organelle insertions called 45S arrays

**Done in 0.11.** Organelle rRNA holds the 18S motif too (rice chr4 9.23-9.37 Mb, reported as a NOR, is a chloroplast insertion). A
motif hit now counts only with nuclear context: the next base and the 19 after it. A hit followed by A falls in two groups, 110 with
at most 5 of 19 bases differing (nuclear) and 1,100 with at least 11 (organelle), 4 between; the base alone would have deleted real
arrays (barley 7H). Over the benchmark 272 -> 261 arrays (this entry had estimated 12 to go): rice chr4, *B. rapa* A02/A06, the Sumatran
orangutan. All 74 kept arrays sampled hold 5.8S genes; removed ones hold at most 3. In the polyploids 651 -> 501.

## 19. Long-unit calls are mostly wrong

**Done in 0.11.** A long-unit array (units over 2 kb) with no support from Hi-C, `tes/` or methylation is reported as a feature, not as
the call (flag `long_unit_unsupported`; `seq_start`/`seq_end` keep it). Over all arms that removes 50 calls: 28 wrong, 8 right, 14 without
truth, all `low`. Hi-C's own call can then stand there: *M. sympodialis* LT671825.1 went from wrong to right. Costs: chicken 28/33 -> 27/33,
*S. pombe* FASTA only 1/3 -> 0/3, one right call in sheep. Nothing yet tells a right long-unit array from a wrong one, so the 8 right calls
go with the rest; another line can still keep them.

## 20. The rDNA spacer tag at its threshold

**Evidence.** A family is tagged as an rDNA spacer when 45S genes flank at least half of its 10 largest arrays. *B.
rapa*'s 238-bp family reaches 5-6 of 10 in 5 accessions (tagged) and 3-4 in 6 (not). In CX1, at 4 of 10, it stays
eligible; since 0.11 the two-way period test ([17](#17-knobs-and-the-size-floor)) takes it out of CX1's centromeric
class, which puts A05 right and narrows A06, but the tag still differs between accessions.
**Approach.** Run the accessions together: a family tagged in one genome of a complex is tagged in all.
**Check.** The 11 accessions against CENH3.

## 21. Moved calls whose sentences disagree

**Done in 0.11.** `caveat` now reads the evidence `basis` reports: on moss's 4 moved calls it says the call rests on the 3D lines and
names the other evidence at that place, and no longer says "no sequence or LTR evidence". No call changes.

## 22. Pa1A in other P. annua plants

**Evidence.** Young CRM insertions of four of the five other plants cluster 9-15 Mb from the assembled plant's Pa1A
call (163.0-169.4 Mb), which suggests other centromeres; J10's lie inside it. C11 and J10 show a plant-only methylation
dip at their clusters. No verdict confirms it: on Pa1A 4 `shared` and 1 `unclear`; over all chromosomes 0 `shifted` and
2 `possible_shift` (J10 Pa4B, Pa5A) of 70, at 4-8× reads.
**Approach.** None in the code until the data can settle it.
**Check.** CENH3 ChIP or deeper long reads from those plants.

## 23. A family's unit length can overstate its copies

**Done in 0.11.** `basis` gives the array's own spacing beside the family's unit where they differ ("165-bp unit, copies spaced 155 bp
here", column `seq_spacing`): in rice, 19 of the 28 large arrays of the 165-bp CentO family are not spaced 165 bp, 17 of them about 155
bp, because the 165-bp profile fits 155-bp copies with a 10-bp gap and scores higher. No call changes.

## 24. Telomeres outside the rules

**Done in 0.11.** Units of 13-29 bp are tried only where no 5-12-bp unit qualifies, and must reach most ends: *C. tropicalis*' 23-bp unit
is found at 12 of 14 ends, with no spurious unit elsewhere. Once the unit is chosen, a stretch that reaches a sequence end counts down
to a density of 0.1 (61 of 6,336 ends lie under half density, none under 0.19). Ends 1,539 -> 1,557 over the benchmark and 3,154 -> 3,214 in
the polyploids, none lost, interstitial arrays unchanged.
