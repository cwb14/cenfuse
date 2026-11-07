# Future development

Open problems in cenfuse 0.13. Each gives the evidence, a proposed approach and the check it needs (a rerun, or new truth
data). Entries closed in 0.11, 0.12 or 0.13 keep their numbers and say what was done; the rest is not implemented. The evidence comes
from the benchmark ([benchmark.md](benchmark.md)) and from a replay of 106 allotetraploid plant genomes (2,355 calls).

| # | problem | needs |
|---|---|---|
| 1 | [calls on chains of partial copies](#1-calls-on-chains-of-partial-copies) | narrowing to the core done in 0.12; CENH3 truth where a chain dominates |
| 2 | [the 40-family cap](#2-the-40-family-cap) | a genome where it costs a call |
| 3 | [satellites specific to subgenomes vs the honesty cap](#3-satellites-specific-to-subgenomes-vs-the-honesty-cap) | tested in 0.11: no gain over random splits |
| 4 | [the decisive CRM test with subgenome-specific lineages](#4-the-decisive-crm-test-with-subgenome-specific-lineages) | CENH3 truth, a null |
| 5 | [a faster search](#5-a-faster-search) | low-complexity cost done in 0.11; the whole benchmark searched again; the cost of real hits |
| 6 | [faster seed clustering](#6-faster-seed-clustering) | identical clusters on every stored seed set; a query partition tested in 0.13: not exact |
| 7 | [the FasTAN crash, reported upstream](#7-the-fastan-crash-reported-upstream) | the user's decision, then a fixed release |
| 8 | [large unplaced scaffolds and dot chromosomes](#8-large-unplaced-scaffolds-and-dot-chromosomes) | a flag for left-out complete molecules done in 0.13; known chromosome lists and karyotypes |
| 9 | [CENH3 truths for allopolyploids](#9-cenh3-truths-for-allopolyploids) | new data |
| 10 | [the LTR retrotransposon line](#10-the-ltr-retrotransposon-line) | measured in 0.12; `te_concordant` against random placements is reported, and the rules stay |
| 11 | [Hi-C anchors where the sequence line is wrong genome-wide](#11-hi-c-anchors-where-the-sequence-line-is-wrong-genome-wide) | done in 0.12 |
| 12 | [cis signals of satellite-free centromeres](#12-cis-signals-of-satellite-free-centromeres) | done in 0.12, with `tes/`; entry 26 for without |
| 13 | [bird microchromosomes in Hi-C](#13-bird-microchromosomes-in-hi-c) | a sparse fold scan; deep bird Hi-C |
| 14 | [the Rabl axis where the fold has no power](#14-the-rabl-axis-where-the-fold-has-no-power) | AT troughs help in 0.12; more small genomes with truth and Hi-C |
| 15 | [centromeres at the edge of a terminal unmappable array](#15-centromeres-at-the-edge-of-a-terminal-unmappable-array) | which-end test done in 0.12; chicken still `low` |
| 16 | [the CHG line has changed no call](#16-the-chg-line-has-changed-no-call) | CHG from the assembled plant, with CENH3 truth |
| 17 | [knobs and the size floor](#17-knobs-and-the-size-floor) | done in 0.11 for maize; other knob-rich lines with CENH3 truth |
| 18 | [organelle insertions called 45S arrays](#18-organelle-insertions-called-45s-arrays) | done in 0.11 |
| 19 | [long-unit calls are mostly wrong](#19-long-unit-calls-are-mostly-wrong) | done in 0.11 |
| 20 | [the rDNA spacer tag at its threshold](#20-the-rdna-spacer-tag-at-its-threshold) | *B. rapa* accessions run together |
| 21 | [moved calls whose sentences disagree](#21-moved-calls-whose-sentences-disagree) | done in 0.11 |
| 22 | [Pa1A in other *P. annua* plants](#22-pa1a-in-other-p-annua-plants) | CENH3 or deeper reads from those plants |
| 23 | [a family's unit length can overstate its copies](#23-a-familys-unit-length-can-overstate-its-copies) | done in 0.11 |
| 24 | [telomeres outside the rules](#24-telomeres-outside-the-rules) | done in 0.11 |
| 25 | [non-CRM centromeric LTR retrotransposons](#25-better-recognition-and-use-of-non-crm-centromeric-localized-ltr-rt) | done in 0.12 |
| 26 | [the cis signal without `tes/`](#26-the-cis-signal-without-tes) | done in 0.13: tested by trans contacts |
| 27 | [wrong LTR cores that look right](#27-wrong-ltr-cores-that-look-right) | more genomes with a non-CRM signature and truth |
| 28 | [where AT troughs do not work](#28-where-at-troughs-do-not-work) | a line for point centromeres and AT-rich genomes; merging troughs of one centromere tested in 0.13: no change |
| 29 | [what the contact crossing a centromere is](#29-what-the-contact-crossing-a-centromere-is) | Hi-C of mutants or of other stages |
| 30 | [a Hi-C test with too few confident references](#30-a-hi-c-test-with-too-few-confident-references) | done in 0.13 |
| 31 | [a centromere that left its satellite and its LTR families](#31-a-centromere-that-left-its-satellite-and-its-ltr-families) | age tie-breaker done in 0.13; CENH3 data for *P. annua* |
| 32 | [calls in a terminal eighth that an interior line contradicts](#32-calls-in-a-terminal-eighth-that-an-interior-line-contradicts) | a test that can see the ends |
| 33 | [cis and trans at one locus as a confident call](#33-cis-and-trans-at-one-locus-as-a-confident-call) | finer truth for small genomes |

## 1. Calls on chains of partial copies

**Evidence.** A family can concentrate in one place without being a tandem satellite: under half of its bp in tandem and under
half in whole-unit hits, so it is a chain of pieces of a longer repeat, usually a retrotransposon. *Glycine dolichocarpa*'s 1,790-bp
family has 44 whole hits of 25,732, 65% inside SIRE elements, and calls 10 Mb wide. In the 106 allotetraploids, 41 such families in
29 species carry 409 of 2,355 calls, 209 of them `high` before 0.11. This was entry 1's "131 calls on single-sequence models" under
another name: those models are nearly all at least 1 kb, and a better profile would not make them tandem. Since 0.11 such a call
is at most `medium`, and `low` without an overlapping CRM core; no call moved. The new search changes these counts little: over its
14 reruns, 75 chain calls become 74 (entry 5).
That rule rests on principle (only CRM decides among retrotransposons), not on truth. No benchmark genome has a chain as its
dominant family; in 0.12 the benchmark's 35 chain calls (chicken, potato, moss, tomato) are 29 right and 4 wrong of the 33 with truth,
all `low` or `medium` from the sequence line; 3 in moss rise to `high` with Hi-C and `tes/`, all right. Proxies, which cannot settle it:
calls on CRM chains lie between TE cores and real satellites in gene density (0.135, 0.054, 0.009), and match their homoeolog's
place within 5% of the chromosome about as often (0.60, 0.70, 0.66).
**Done in 0.12.** Where a decisive core lies inside a chain call and is narrower, it becomes the call (`chain_to_core`): such cores are
tighter and more gene-poor than the calls they overlap (gene density 0.021 vs 0.053), with about the same homoeolog agreement.
**Check still open.** CENH3 for a genome where a chain dominates (6 cotton species have one on 23-26 of 26 chromosomes); entry 9.

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
**Tested in 0.13, not built.** A fixed partition of the queries over single-thread BLAST processes, on *Acanthus*' stored rounds: the
heaviest round (200 x 200 seeds, 773,595 hits) ran slower (0.79x at 4 workers against BLAST's own 4 threads) and its hit lists were not
identical to one run's: with thousands of HSPs per pair, BLAST's output depends on how the queries are batched. Small rounds were identical
but gained little (0.65-1.39x).
**Approach.** Fewer HSPs to report: a pair already known to match (or not) needs no further HSPs, so a two-pass search (first the best HSP per
pair, then full tiling only for the pairs the clustering consults) would cut the output that dominates heavy rounds.
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
**Done in 0.13, in part.** The longest sequences left out of the set (up to 100, at least 200 kb) that hold arrays of the
genome's telomere unit at both ends are named in the report and in `left_out_telomeric`: complete molecules, such as a dot
chromosome after a tenfold gap or a B chromosome. They are flagged, never added.
**Approach, still open.** The other edge: large scaffolds kept in the set without a telomere, a centromere-like array or genes.
**Check.** The 121 lists as truth; *D. melanogaster* and bird karyotypes; the 383 assemblies of the rule's evidence
as a no-change control.

## 9. CENH3 truths for allopolyploids

**Evidence.** No allotetraploid in the benchmark has CENH3 truth. Where a CRM
core replaced a satellite array (143 of 2,212 audited allotetraploid calls), homoeolog agreement favors the core (26
vs 14 pairs, P = 0.08) and gene density the array (a biased proxy).
**Approach.** CENH3 ChIP-seq or CUT&Tag for a few allotetraploids chosen to settle entries 1, 3, 4, 5 and 10
(`seq_overruled_*` lists the replaced arrays).
**Check.** Score as in the benchmark; no rerun needed.

## 10. The LTR retrotransposon line

A check on 17 genomes with `tes/` and truth (190 chromosomes) pointed to these changes. Status in 0.12, on the same genomes:

- **Show which families the CRM line rejected. Done in 0.11:** columns `seq_rejected_family`, `seq_rejected_start` and
  `seq_rejected_end`, and a caveat naming it. 19 of 30 calls that moved with `tes/` had such a family.
- **One standard for replacing a `low` call. Measured, left as it is.** A core replaces a disagreeing `low` call without the head-to-head
  test a confident call must pass. Among such moves the test would have blocked, 2 are right and 1 wrong: no case for a stricter rule, and
  still too few to settle it.
- **An inconclusive head-to-head test. Left as it is.** When a core disagrees with a call and cannot overrule it, the call drops one level
  (`te_discordant`). The `tes/` arms now hold 3 such calls with truth, all in tomato, all `low` and all wrong: no case for keeping the level.
- **`te_concordant` is a weak flag. Measured, left as it is.** Two estimates agree if the gap between them is at most the wider of the two,
  so a wide core or call agrees with much of its chromosome. But the calls that keep the flag are right 125 of 130 times, and the widest
  are no worse (35 of 35 above a tenth of the chromosome), so a stricter rule has nothing to gain. What the agreement is worth against chance is
  now reported per genome (`asm.line_tests.tsv`).

## 11. Hi-C anchors where the sequence line is wrong genome-wide

**Done in 0.12**, in two ways. Better calls make better anchors: with the LTR signature, barley's ordinary Hi-C test passes on 7 of 7
chromosomes and moss's on 26 of 26, with no new rule. And where the test still fails on most chromosomes, Hi-C is tested against the LTR
cores directly ([design.md](design.md#hi-c-and-the-ltr-line-test-each-other-hicfusepy)). Potato anchored on its cores finds 3 of 12
(P 0.02), which does not pass; its centromeres come from the cis signal instead (entry 12). Moving the windows to the best score near each
core found 6 of 12, but randomized Hi-C then passed too, so that was removed ([the benchmark, section 12](benchmark.md#12-how-the-benchmark-was-run)).

## 12. Cis signals of satellite-free centromeres

**Done in 0.12.** The contact crossing a bin, as observed over expected at each distance, has its extreme at CENH3 on all 12 tomato
chromosomes and 8 of 11 in potato. The LTR cores of the other chromosomes choose its sign, and test it, with chance from simulated contact
maps. It moves `low` calls, and raises or lowers others by one level. See
[benchmark.md](benchmark.md#the-contact-crossing-the-centromere). Soybean, sorghum and the mammals have no such signal.

## 13. Bird microchromosomes in Hi-C

**Evidence.** Hi-C bins finer than 250 kb are capped by the largest chromosome (4,000 bins): chicken's 200-Mb chr1
needs 4,001 at 50 kb, so its 2.5-7-Mb microchromosomes stay out of the Hi-C line.
**Approach.** A sparse or banded fold scan, which lifts the per-chromosome cap.
**Check.** Chicken with deeper Hi-C; calls on the microchromosomes against the CENP-A windows.

## 14. The Rabl axis where the fold has no power

**Evidence.** The trans arm alignment passes in *C. neoformans* (anaphase, 12/14) and *M. sympodialis* (6/8), and
Hi-C's prior-free calls there are right on 10/10 and 7/7, but the arm fold, which must lean the same way before the
alignment may act (the holocentric guard), has no power in these small nuclei (1/14; 0/8), so the alignment alone raises nothing. In 0.12, AT troughs that agree with Hi-C raise
*M. sympodialis*' 8 calls and 4 of *C. neoformans*' to `medium` (`at_pair`); the rest stay `low`.
In holocentric *C. elegans* the alignment reached 2/5 (P 0.095), and a position-matched null that rejects it also
rejects barley. In G1 *C. neoformans* the prior-free peaks sit at the chromosome ends, so the alignment's partners
cannot be placed (0/14), although anchored on CENP-A 8/14 chromosomes still align.
**In 0.12.** AT troughs give a second line where they pass (*N. crassa*, *P. falciparum*, *M. sympodialis*). The cis profile lies on
CENP-A on 13 of 14 chromosomes in G1 *C. neoformans*, where Hi-C's other readouts fail, but nothing tests it there (entry 26).
**Approach.** A compartment-adjusted alignment (trans contact explained by compartment identity removed) that can
stand alone; or partners placed on the sequence line's confident calls, as a corroboration test.
**Check.** *Candida*, *Cryptococcus* relatives, *Toxoplasma*, *T. brucei* with truth and Hi-C; small-chromosome
holocentrics (nematodes, Lepidoptera) and randomized Hi-C must never pass.

## 15. Centromeres at the edge of a terminal unmappable array

**Evidence.** Chicken's acrocentric microchromosomes carry their centromeric array at the very end, unmappable for
1-1.25 Mb: their Hi-C peaks lie at its inner edge, which in mappable bins is the end, so they stay `telomeric`, and
chicken's calls stay `low`.
**Done in 0.12** (`hic_end`). Telomere clustering lifts both ends of a chromosome, a centromere one. The anchored score is about twice as high
at the end holding the call in goat (28 of 30), sheep (22/23), horse (18/19) and mouse (19/19), so a call at the higher end rises a level
where that holds on most of the other chromosomes (binomial, P ≤ 0.01). No randomized library passes.
**Still open.** Chicken does not pass (12 of 19), so its calls stay `low`.

## 16. The CHG line has changed no call

**Evidence.** Only *P. annua* and *Arabidopsis* have CHG data here, and across all 3,773 chromosomes of the release
runs the line supported no call and backed no move. *P. annua*'s five samples come from other plants, so they do not
decide; *Arabidopsis* finds no signature (the expected domain, from genome size, is 315 kb; its CEN180 arrays span 2.2-3.5 Mb). The signal is real in *P. annua*: the dip is the
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

**Evidence.** Young CRM insertions of four of the five other plants cluster 9-15 Mb from the functional centromere that earlier
analyses found on Pa1A (163.0-169.4 Mb), which suggests other centromeres; J10's lie inside it. C11 and J10 show a plant-only
methylation dip at their clusters. 0.12 calls Pa1A at its satellite instead (entry 31), and the plant verdicts there are 2 `shared`,
2 `unclear` and 1 `possible_shift` (R1, toward the Hi-C site 20 Mb away); over all chromosomes 1 `shifted` (J10 Pa4B) and 2
`possible_shift` (R1 Pa1A, J10 Pa5A) of 70, at 4-8× reads. No verdict settles it.
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

## 25. Better recognition and use of non-CRM centromeric-localized LTR-RT

**Done in 0.12.** cenfuse learns each genome's signature at family level, leave one chromosome out, and a lineage other than CRM decides at
10-fold enrichment or more, unless a dominant satellite passes the one-per-chromosome test
([design.md](design.md#ltr-retrotransposons-tes-satsignaturepy)). Moss's Bryco, *Chlorella*'s Osser and tomato's Tekay family now decide;
*Arabidopsis*' Athila only flags, as CEN180 passes the test.

## 26. The cis signal without `tes/`

**Evidence.** Tomato's cis peaks alone are right on 12 of 12 chromosomes, and potato's dips on 8 of 11, but today the LTR cores must choose the
sign and test it. Two ways to do without them were tried. One locus per chromosome: tomato passes (11 of 12 covered, P 0.003), potato and
*Verticillium* miss after correcting for two signs (P 0.044, 0.047), and chicken's dips pass (P 0.0008) while 2 of 26 lie on CENP-A. Hi-C's own
prior-free peaks: no agreement in tomato, potato or rice, and in *Ae. tauschii* 6 of 7 agree with 0 of 7 right.
**Done in 0.13**, with the library's own trans contacts as the test ([design.md](design.md#the-cis-and-trans-readouts-of-hi-c-test-each-other-hicfusecis_trans_line)):
each chromosome's cis extreme, outside its terminal eighths, anchors the trans test of the LTR cores. Tomato: 10 of 12 chromosomes'
trans peaks lie at their cis peaks (P 6e-14), all on CENH3; *C. neoformans* in G1: 12 of 14. Read over whole chromosomes, *Ae. tauschii*
passed with 0 of 4 agreements on CENH3, all in terminal eighths; randomized soybean Hi-C reached support (6/20), so the line acts only
at pass. Potato does not pass (its dips are right on 9 of 11, but its trans contacts at them are not). Chicken, the mammals, the
holocentric genomes and every randomized library fail. What remains: potato without `tes/`, and the ends (entry 32).

## 27. Wrong LTR cores that look right

**Evidence.** Tomato's signature has significant cores on all 12 chromosomes, and 5 of them miss CENH3. Without Hi-C two of those misses stay `medium`
(CM175074.1, CM175078.1), and with Hi-C they fall to `low` because the cis peak lies elsewhere; on both, the cis peak is the right one. Soybean chr13 is
similar: a CRM core 1.8 Mb from CENH3 at `medium`. *Arabidopsis*' Athila (13.3-fold, above the level at which a non-CRM lineage may decide) has
its cores beside CENH3 on all 5 chromosomes, 0.46-1.38 Mb from the domain; it only flags because CEN180 passes the one-per-chromosome test. In 6 of 27
polyploids above 10-fold (three *Sphagnum*, *Lotus subbiflorus*, *Malus toringoides*, *Solanum stoloniferum*) the cores and the FASTA-only calls mostly
disagree, with no truth to say which is wrong.
**Approach.** Nothing in the LTR line tells these apart so far: their fold and share look like the right cores'. Two leads: a core that sits beside
the dominant satellite array rather than over it (*Arabidopsis*), and the Hi-C test anchored on the cores wherever Hi-C exists.
**Check.** More genomes with a non-CRM signature and truth; *Arabidopsis* with CEN180 masked must not move its calls onto the Athila cores.

## 28. Where AT troughs do not work

**Evidence.** The test passes in *P. falciparum*, *M. sympodialis* and *N. crassa*. It fails in *Z. tritici*, *C. tropicalis*, *C. neoformans* and *V. dahliae*,
whose deepest troughs do not fall one per chromosome, and it can say nothing about point centromeres (yeast: 125 bp) or genomes with few chromosomes (*S. pombe*:
`untested`). It passes on about 4% of genomes with shuffled GC blocks, as its level of 5% allows, and then gives `low` calls.
**Tested in 0.13, not changed.** Counting the troughs of one AT-rich run once (each trough's half-depth run excluded before the next)
changes no genome: *C. neoformans*' 30-65-kb centromeres hold two sub-troughs with GC-normal sequence between them (CEN11: 0.865 and
0.925 Mb), so the run breaks, and *V. dahliae* has three deep troughs Mb apart on one chromosome. In *C. neoformans* the cis signal tested
by trans contacts now places them (entry 26), and an AT trough agreeing with a weak prior-free peak places 3 more (`at_denovo_pair`).
**Approach.** A second feature with the trough (a gene desert, without `genes/` a run of low k-mer complexity).
**Check.** *Z. tritici* and *C. tropicalis* pass; holocentric genomes and shuffled genomes do not.

## 29. What the contact crossing a centromere is

**Evidence.** Piled up on CENH3, tomato's two flanks contact each other 2.4-fold within 0.25 Mb and back to background by 1 Mb; potato shows a smaller central peak between deep
dips; both lose contact with the arms out to at least 6 Mb. *S. pombe* has inter-arm contacts close to its centromeres that depend on heterochromatin
([Mizuguchi et al. 2014](https://doi.org/10.1038/nature13833)). A fold of the two flanks, or a domain insulated from the arms, could give all of these, with the sign
depending on the domain's width against the 1 Mb the profile reads. That is a hypothesis: nothing here perturbs it.
The chance the test uses comes from contact maps simulated with each bin's coverage and the chromosome's mean contact at each distance.
On moss's re-paired control map the dip lies at the LTR cores on 11 of 26 chromosomes where the simulation expects about 5 (P 0.011,
just short of passing). One such value among 20 control tests is expected by chance, but if the simulation is short in some genomes,
balancing each simulated map again (instead of reusing the real map's balancing) would be the first thing to try.
**Approach.** Rerun the profile at several box sizes; look for genomes where the sign changes with the box.
**Check.** Hi-C of heterochromatin or condensin mutants, or of other cell-cycle stages, in a genome with this signal.

## 30. A Hi-C test with too few confident references

**Evidence.** The Hi-C test takes the `high` calls as references (with `medium` ones if too few), and falls back to the `low` calls at clean peaks
only when the genome has no `high` or `medium` call at all. Cucumber 9930 (held out) without `tes/` has 7 `low` calls, and Hi-C is tested on them;
with `tes/`, the Tekay signature raises 2 calls to `high`, so 2 references remain, too few, and the test is `untested` on all 7 chromosomes. No call
differs: on the low references Hi-C changed no level either.
**Done in 0.13.** When the confident references are fewer than the test needs (4), the `low` calls join them, agreements at clean peaks
only for the `low` ones (`hic_gate_refs` `mixed`). No benchmark call with truth changes; the randomized controls stay silent.
**Check still open.** A genome where this lets Hi-C act, with truth.

## 31. A centromere that left its satellite and its LTR families

**Evidence.** *P. annua* Pa1A: earlier analyses (no ChIP) place its functional centromere at 163.0-169.4 Mb, where Hi-C peaks (168.0 Mb),
21 Mb from its largest satellite, and 0.11 called it there. In 0.12 the CRM families learned from the other 13 chromosomes concentrate
around the satellite (145.1-148.3 Mb, 30 of their 146 elements on Pa1A), so the satellite and the LTR line agree on the old site and
the call stays there (`medium`, `contested`, `hic_offset`), although the Hi-C test passes on 14 of 14 chromosomes. With randomized Hi-C
the same call is `high`. A family that marks centromeres on the other chromosomes marks where centromeric families accumulated, which a
centromere that has moved leaves behind; only Hi-C sees the new site. Insertion age, an independent line, agrees: CRM copies at the Hi-C site are young (median age about 0, 89% under 100,000 years; 62 copies)
and those at the satellite site old (median 325,000 years, 39% young; 28 copies); within one family, 0 vs 157,000 years (P 2e-7); 61 of
Pa1A's CRM copies under 50,000 years old lie at 163-172 Mb and 2 at the satellite (`bench/pa1a_age.py`). CRM inserts into active
centromeric chromatin, so young copies mark where the centromere is now, as in maize, where recently active CRM subgroups colonize
the present centromeres and older ones their flanks ([Wolfgruber et al. 2009](https://pmc.ncbi.nlm.nih.gov/articles/PMC2776974/)).
No CENH3 data place Pa1A's centromere. In the scored benchmark, Hi-C contests a sequence call
(`hic_offset`) on one chromosome only, CHM13 CP068258.2 in two runs, and there the sequence call is right: one case each way.
**Done in 0.13** (`ageline.py`). Where validated Hi-C's strong peak lies far from a `high` call, copies of the chromosome's signature
families (learned from the other chromosomes) that are younger at the peak than at the call (one-sided Mann-Whitney, P ≤ 0.01, at least
3 dated copies at each site, both widened by 4 bins) are the second line that lets Hi-C move the call (`moved_with_age`), in a genome
where the premise holds: on the other chromosomes the signature's younger copies lie in the cores more often than its older ones
(Fisher, P ≤ 0.01; *P. annua* leaving Pa1A out: P 5e-10). As a locator on
its own, the core of the young copies does worse than the density core (on CENH3 for 128 vs 155 of 172 chromosomes with truth;
`bench/young_core.py`), so age only breaks a tie between sites that another line raised. In the benchmark it can act only on Pa1A.
**Check.** CENH3 ChIP or CUT&Tag of *P. annua*; genomes with known evolutionarily new centromeres and Hi-C; the randomized Hi-C controls stay silent.

## 32. Calls in a terminal eighth that an interior line contradicts

**Evidence.** The cis signal tested by trans contacts (entry 26) is read outside the terminal eighths, because there both readouts see the
chromosome ends. It therefore has no say on a call inside one: tomato without `tes/` keeps 2 wrong `low` satellite calls 8.6% and 11.8%
from an end (CM175073.1, CM175075.1), while its interior extremes lie on CENH3. Moving such a call would be right in tomato and wrong for an
acrocentric chromosome among metacentric ones, which the line cannot see.
**Approach.** A second line that can see the ends: the which-end test (anchored score at the two ends) where it passes, or an LTR core.
**Check.** Tomato moves its 2 calls; acrocentric mammals (goat, sheep, mouse) keep theirs.

## 33. Cis and trans at one locus as a confident call

**Evidence.** Where the cis signal passes its test by trans contacts, a chromosome whose own trans peak lies at its cis extreme holds two
readouts at one locus, as `te_hic_pair` holds an LTR core and Hi-C. Such calls stay `low` in 0.13: in *C. neoformans* G1, 10 of the 12 pairs
lie on CENP-A and the other 2 lie 2 bins (19-20 kb) from it, wrong by the 1-bin rule, so as `medium` they would be confident errors.
**Approach.** `medium` for a pair where the bin is coarse relative to the centromere, or a finer interval from the AT trough.
**Check.** Small genomes with truth at 10-kb bins; the randomized controls stay silent.
