# Future development

Open problems after cenfuse 0.19.0, for developers. Each entry gives the problem, its evidence, what would settle it
(**Settle**) and where to start (**Code**). Entry numbers are stable (the code, tests and docs cite `future.md #N`); a closed
entry keeps its number as one line under [Closed entries](#closed-entries). Release history: [CHANGELOG.md](../CHANGELOG.md).
The next new entry is 82.

**Status.** *open*; *partly done*; *built, unchecked* (in the code, untested by truth); *waiting for data*; *user's decision*.
Ideas measured and not adopted: [Tested, not adopted](#tested-not-adopted) and [Tried and set aside](#tried-and-set-aside-not-in-designmd).

**Words.** *Arm*: one benchmark run of a genome with one input set; *chromosome-run*: one chromosome in one arm. *Release runs*:
every arm rerun for a release, most replaying a stored satellite search (`--satellite-dir`). *Confident*: `high` or `medium`.
*Headroom*: chromosome-runs where some line's own locus is right but the final call is not (`benchmark/headroom.py`). *The panel*:
106 allotetraploids, replayed to check that homoeologous calls agree (#68); 4 of them have CENH3 truth since 0.17.

**Where things live.** The repo: `cenfuse/`, `tests/`, `docs/`, `benchmark/` (truth manifest, frozen scorer, the table chain
`benchmark/release.sh` that writes docs/benchmark.md). docs/constants.md is generated from the code (`python -m cenfuse.consts`):
edit the code. Truth files, stored runs and analysis scripts live outside the repo under `$CENFUSE_BENCH_DATA`, one nested
folder per release (v2-v7 = 0.1.0 development, v8 = 0.8, v9 = 0.9, v10 = 0.10.0, v10.1-v10.3 = 0.10.1, v10.4 = 0.11, v11 = 0.12,
v12 = 0.13, v13 = 0.14, v14 = 0.15, v15 = the 0.16 CpG line, v16 = 0.16, v17-v18 = 0.17, v19 = 0.18, v20 = 0.19); scripts are cited as `v14 bench/x.py`.

## Rules for any change

An entry is done only when these hold (design.md, Principles).

1. **Write the check first**: what should change, on which genomes, and what must not.
2. **Give every constant a reason**: a convention (P ≤ 0.01), a majority, a value measured in the data, or a named heuristic. No
   threshold options; never a value chosen because it matches known answers.
3. **A rule prompted by one case** needs a reason independent of that case, a recheck on cases not used to make it, and a row in
   benchmark.md's table of rules changed after a result was seen. A rule built on one genome waits for a second, or stays a named
   heuristic that makes no confident call (#70, #72).
4. **No data are held out** (user's decision, 2026-11-08): every benchmark genome, the 13 tomatoes, *S. verrucosum*, *Ae.
   geniculata* and the panel are development data. Rule 3 does the work that holding data back did.
5. **Replay on the panel too** (#68). It has no truth, but it caught a bad 0.15 rule and three 0.16 candidates the benchmark favored.
6. **Keep the controls silent** (`cenfuse control`; the methylation one has no generator, #34): randomized Hi-C (57 runs) passes no
   test and changes no call; shuffled LTR positions (270 shuffles of 54 genomes) give no signature; holocentric genomes (196
   chromosome-runs) get no confident call; misplaced methylation passes no test; GC-block shuffles pass the AT test at most at its
   5% level (3.9% of 3,336).
7. **Do no harm**: no right call becomes wrong and no confident error is added, or each change is listed. Right, then unsure
   (`low`), never confidently wrong.
8. **Report a missed criterion as written**; do not reword it.
9. **Test every rule** on a small synthetic fixture. Develop the LTR line on LTRquest annotations only.

## Where to start

0.19: a census of the 0.18.0 release runs found that the integration already uses nearly everything its tested lines offer:
of the 599 development chromosome-runs without a right call, a line is right in 439, but a line whose own test passed in few
(the repeat-family window 147 times, 70 of them on calls in a terminal eighth, where it cannot see, and 66 within both
readouts' resolution; repeat density 13, and 2 in *S. pombe* where only Hi-C validated it; Hi-C 10; the AT trough and the
Rabl axis about none). What remains needs a better
line or more truth, not another rule over the same readouts: a calibrated family test (#76; 0.19 lets the family line act
where it agrees with repeat density on most chromosomes, `fam_rep_pair`), a second line for LTR cores alone (#79) and
for centromeres away from the main satellite (#78), and the leads of #45 (an array inside a core holding few of its
elements: tested in 0.19, two criteria missed), #80 and #81.

0.17: on the 22 genomes added in 0.17, scored after 0.16.1's rules were fixed, `medium` calls were right 74% of the time against
88% on the development genomes; their confident misses come from centromeres away from the genome's main satellite (#78) and
from CRM cores (#79, #46). 0.16.0: 2,309 of 2,624 chromosome-runs with truth right in the real arms (88.0%), 35 confidently
wrong. The largest headroom was in retrotransposon genomes from the FASTA alone, where the repeat-family window was right on 89
of the 110 `low` calls it contradicted but only flagged them; since 0.18 it replaces an unsure satellite call it contradicts
beyond both readouts' resolution, where its window holds a cluster and nothing else speaks to the call (#76: on the
development genomes 94 wrong calls became right and 2 right ones wrong; on the genomes added in 0.17, 11 and 0), and 77 calls
remain flagged beside their window, within that resolution (window right 64 times, call 13). Next come the competing
clusters (40 chromosome-runs where the second array is right) and yeasts' AT troughs, "right" only at the 250-kb tolerance
(#37). Confident misses: [benchmark.md](benchmark.md#every-confident-miss).

1. **#79, #46, #78**: confidence on genomes cenfuse was not developed on: a CRM core alone, and centromeres away from the
   genome's main satellite (Hi-C or methylation needed).
2. **#76, #27, #72**: retrotransposon centromeres: truth beyond *Solanum* where the repeat-family line passes and the calls are
   `low`; Hi-C settles them only where the library shows the centromeres (TS-60).
3. **#35**: the public truth sets.
4. **#32, #27**: `cis_end_pair` and `cis_trans_pair` are right on 4 of 4 in two genomes each (tomato, TS-60); a third genome
   with the case would let them make confident calls.
5. **#45**: which locus holds CENH3 where an array sits at the edge of a wide LTR core (0.19 tested a form that lowers such
   calls: 8 confident errors and 4 right calls lowered on the development genomes, one right call on the added ones; it missed
   two of its four criteria, and the numbers are in #45).
6. **#37, #38, #50**: the mammal misses (6 scored against other animals' ChIP or another caller).
7. **#74, #76**: the repeat lines' genome-wide sort after reading runs on one thread (the reading is parallel since 0.18).
8. **#63, #60, #65**: the remaining large functions; CI's first result.
9. **#80** (the user's decision) and **#81**: confidence for a composition line alone; the family line tested against Hi-C.

## All entries

Grouped by area below; each entry opens with its status.
- **Open:** 2, 6, 16, 27, 28, 38, 41, 42, 43, 45, 48, 49, 50, 51, 54, 55, 57, 61, 70, 71, 72, 77, 78, 79, 81.
- **Partly done:** 1, 5, 8, 9, 10, 14, 15, 34, 35, 37, 39, 44, 46, 47, 52, 53, 56, 58, 60, 62, 63, 65, 66, 73, 76.
- **Built, unchecked by truth:** 30, 31, 32, 68, 74. **Waiting for data:** 4, 13, 22, 29. **User's decision:** 7, 80.
- **Tested, not adopted:** 3, 40, 67, 75 ([below](#tested-not-adopted)). **Closed:** 11, 12, 17-21, 23-26, 33, 36, 59, 64, 69
  ([below](#closed-entries)).

## Benchmark, truth and validation

### 9. CENH3 truths for allopolyploids

Partly done. One allotetraploid has CENH3 truth on its own assembly, *Ae. geniculata* (Yunusbaev 2026, Table S2, Mb
resolution): 0.15 placed 0 of 14 from the FASTA (all `low`) and 13 of 14 with `tes/`, 1 confidently wrong (chr05U, #45).
- **Genomes added in 0.17** ([benchmark.md](benchmark.md#genomes-added-in-017)): four panel allotetraploids now have CENH3 ChIP
  on their own assembly. Cotton ZM113: 24 of 26 right from the FASTA (all `low`), 25 of 26 with `tes/` (1 confidently wrong,
  #78). *Brachypodium hybridum* and *Erianthus rockii*: 15 of 15 each, from the FASTA and with `tes/`. *N. benthamiana*: 10 of
  18, with 6 confidently wrong from the FASTA and 8 with `tes/`, on chromosomes whose centromere sits away from the genome's
  main satellite (#78), not from polyploidy.

### 34. The benchmark cannot be rerun from the repository

Partly done. In the repo: the truth manifest (`benchmark/truth.tsv`), the frozen scorer, the table chain (`benchmark/release.sh`)
and the controls (`cenfuse control`). Outside: each release's doc chain after the repo chain (v16 `bench/docs16.sh`: the genomes
added since 0.14, the tomatoes), the panel check (v14 `bench/polyeval.py`), the truth-making scripts (ChIP mapping, domain calls,
fetches), the analysis scripts the docs cite, and the misplaced-methylation control's generator (never found). Release runs replay
searches mostly made before 0.11 (#5); fresh runs equal the replays where tried (0.15: *V. dahliae*, *M. sympodialis*).
**Settle:** a clean checkout regenerates benchmark.md from the stored runs; a fresh run of every arm reproduces the published
counts or lists each difference.

### 35. Truth on disk or public, not yet used

Partly done (*Ae. geniculata*, *S. verrucosum* and 13 tomatoes scored since 0.14; 22 genomes added in 0.17,
`benchmark/truth_017.tsv`). Confirm each antibody and target first (the
horse CUT&Tag target and the pig antibody are unconfirmed); score with the frozen ChIP rule (#37). Sheep and goat CENP-A ChIP
gives no truth ([Tried and set aside](#tried-and-set-aside-not-in-designmd)).

| data | what it would test |
|---|---|
| tobacco (*N. tabacum*), eggplant CENH3 (accessions to find; *N. benthamiana* was added in 0.17) | CRM that is not centromeric (#46) |
| Grevy's zebra mEquGre1 (GCA_057876975.1; CENP-A ChIP SRR17956804, input SRR17956803, other animals; Cappelletti et al. 2022, 10.1093/molbev/msac162) | 13 of 23 centromeres satellite-free, several beside inactivated satellite centromeres: #78 in a mammal |
| Malayan tapir mTapInd1 (GCA_031878655.1; CENP-A ChIP SRR35264170, input SRR35264171, the assembled animal; Biundo et al. 2026, 10.1038/s42003-026-10612-2) | 2-3 repositioned satellite-free centromeres (#50) |
| *N. benthamiana* NbT2T Hi-C (GSA only, 180-250 GB) | the 8 chromosomes whose centromere sits away from the main satellite (#78) |
| *Leishmania major* KKT1 ChIP (ERR2030267, input ERR2030268) | no truth at 4-fold (1.5-2.7-fold at best, 10 M reads): kinetoplastid kinetochores need their own rule |
| wheat Chinese Spring (T2T GCA_040256815.2; ChIP SRR9089557/8, DRR821424/5) | a hexaploid; CRM plateaus (#45) |
| *Avena sterilis* (GWHERDU00000000; CRR1339975/6); *A. maroccana*, *A. atlantica*, *A. eriantha*, 'Starter' (CRA007177) | oat's non-centromeric satellite on every chromosome (#42); allopolyploids |
| 6 accessions of Shi et al. 2026, Supp Table 8, published without deposited ChIP | Solanaceae LTR rules (#46, #72), once how their centromeres were found is confirmed |
| 4 more soybean lines (ChiHei23, QiHuang34, ZH13, Williams82); the pan-centromere of 27 accessions (PNAS 2023, Dataset S1; 0.1-0.5-Mb resolution; needs each accession's 2020 pan-genome assembly) | CRM plateaus (#45); centromeres that move within a species (#31) |
| cattle CENP-A CUT&RUN and input (SRR28462359/8) | an all-acrocentric ruminant; the which-end rule (#15) |
| hoolock gibbon CENP-A CUT&RUN | the siamang's failure mode |
| siamang same-cell-line CENP-A DiMeLo-seq (bigWig, 15.1 GB) | the siamang's 2 confident disagreements; `low` references (#30) |
| horse: same-animal assembly (GCA_052323815.1), Hi-C (SRR36069624), CUT&Tag/ChIP with IgG (SRR36096936-9, PRJNA1248178); CENP-A ChIP of about 15 horses (96 runs) | horse chr11 (#50); the horse Hi-C that fails; centromere sliding |
| pig CENP-A ChIP and input of the assembly's breed, Duroc liver (DRR1139646-9; GSA PRJCA032804; antibody not stated) | pig's 3 confident misses (#37) |
| the repaired ape assemblies (56 centromeres misassembled in v2.1, with repair BEDs), both haplotypes | about twice the ape truth; misassembly as a positive control (#43) |
| HGSVC3 CDRs (67 genomes, 2,156 centromeres; Gao 2026; figshare 10.6084/m9.figshare.29921990) | the CpG line (#55) |
| holocentric CENH3: *R. breviuscula* ChIP (SRR27965399/8), *Chionographis japonica* (GSM7120336), *Luzula sylvatica* (358 domains); the last two added in 0.17 as controls, their domains not yet used | the satellite line on satellite holocentromeres (#44) |
| same-plant CHG with CENH3: *S. verrucosum* CPC54, potato DM FAST5 (about 625 GB), rice T2T-NIP modBAMs | the CHG line (#16) |
| Hi-C: *Chlorella* NS4-2 (CRR559578), the primary Bd21 library (CRR611509, 289 M pairs), deeper libraries for the shallow controls (#37) | Hi-C on LTR cores in an alga; the Rabl axis in a grass (#14) |

### 37. Scoring rules and controls to tighten

Partly done (0.14: call and domain widths and midpoint distances reported, `benchmark/precision.py`; the truth rule at 3×, 5× and
50-kb bins moves no outcome). Open:
- **The chance model** draws calls uniformly along the chromosome, optimistic where centromeres sit near a known end (goat, sheep,
  horse, mouse). Condition it on end distance.
- **250 kb without Hi-C** for fungi and *Plasmodium*, whose centromeres span 1-100 kb: *S. cerevisiae*'s AT troughs are "right" on 12
  of 16 that way and overlap on 2. Score small genomes without Hi-C at a fine fixed tolerance. Barley is scored at 8 Mb without Hi-C.
- **Weak truth.** Pig and horse ChIP come from other animals; pig's 3 confident misses lie where that ChIP's second peak (3.0-4.4×)
  falls, and the karyotype agrees with the calls (chr11 bi-armed but its domain at CI 0.002; chr17 acrocentric but at CI 0.24).
  Siamang, sheep, goat, rhesus and giraffe are scored against another caller. Bornean orangutan CM054652.2 may be
  haplotype-specific (#38); *N. crassa* LG VI carries an inversion (Galazka 2016); *C. tropicalis* CENP-A comes from an older
  assembly; *V. dahliae*'s authors inferred its centromeres from this Hi-C; *C. neoformans* counts twice in the calibration. The
  ChIP truth rule (100-kb bins, ≥ 4×) has no stated reason for either number; in cucumber 1 of 7 chromosomes passes it.
- **Shallow controls** (library prefixes): *R. breviuscula* 3.3 M pairs, *B. mori* 5.9 M, *Ae. tauschii* 0.84 M (another
  accession), sorghum 3.9 M, mouse 15.7 M; "no holocentric genome passes" rests on them (#51). Deeper: *B. mori* SRR25988387 (369
  M), *R. breviuscula* SRR17429201 (863 M), soybean SRR37686999 (415 M), *B. rapa* CC1 SRR8633037 (368 M), Col-0 SRR2626163 (143 M).
- **Attribution** (user's decision): benchmark.md credits the ape CDRs to Logsdon et al. 2024; the development notes record an
  unpublished curated set (v1.3, "contact the lab before publication"), and the ape DiMeLo-seq and CUT&RUN tracks as unpublished.
- **Small genomes added in 0.17** are scored within 25 kb (fixed before any of their runs was scored); the development fungi
  stay at 250 kb so that releases compare. **Wide calls.** A call is right however wide it is: *Populus euphratica*'s FASTA
  calls span 5-12 Mb (19 of 19 right, a chance of 44% per chromosome), donkey chr1's 52.7 Mb (26%). Report width and chance
  beside each outcome, not only pooled.
- **Settle:** each scorer change reported beside the old numbers for one release; no README claim rests on overlap alone (an animal
  call is about ten times its CENP-A domain).

### 68. Homoeolog agreement on the polyploid panel

Built in 0.15 as a check, no truth. Each of the 106 allotetraploids is replayed from classify; a call lifted through the
collinear anchors agrees when it lies within 5% of the chromosome of its homoeolog's (v14 `bench/polyeval.py`). 0.15's rules: 663
of 1,103 pairs agree over 105 species (0.14: 657). It stopped a 0.15 rule that let a signature with unspecific cores decide in
*Arachis* (11.6-fold seed, 1.5-fold cores; agreement on the changed pairs 11 → 2 of 14), and three 0.16 candidates (#45, #74,
#76). A disagreement can be a right call on an evolutionarily new centromere. **Settle:** agreement calibrated against CENH3 on
one panel species (#9).

### 79. How far confidence holds on genomes cenfuse was not developed on

Open (0.17). 0.16.1's rules were fixed before the 22 genomes added in 0.17 were scored
([benchmark.md](benchmark.md#genomes-added-in-017)). On the 17 with truth, each with its run with the most inputs, 0.16.1's
`high` calls were right 104 of 110 times (95%), `medium` 23 of 31 (74%) and `low` 57 of 83 (69%), against 453 of 459, 137 of 155
and 149 of 197 on the development genomes. Most confident misses come from centromeres away from the genome's main satellite
(#78) and from `medium` calls on an LTR core alone.
- **An LTR core alone** (`support` te, always `medium`): right on 107 of 126 such calls in the 0.16.1 release (Tekay 81 of 99,
  CRM 21 of 22, Osser 5 of 5; 112 of the 126 in *Solanum*) and on 19 of 29 on the new genomes, all CRM (pepper 14 of 22,
  *Populus euphratica* 3 of 3, *O. australiensis* 2 of 2, *Sesbania* 0 of 2). Registered beforehand: below 80% across at least 3
  lineages held out from it (5 of 7; pepper was its basis), the property that best separated right from wrong calls in the
  development genomes was to be checked on the new ones. That property, a satellite array elsewhere on the chromosome
  (`satellite_away_from_core`: 32 of 41 right against 75 of 85), does not hold there (13 of 20 against 6 of 9; pepper reversed).
  No rule. CRM's standing (#46) is the clearest suspect: 21 of 22 in development, 19 of 29 here; but holding CRM to the other
  lineages' 10-fold would also send pepper's 14 right core calls back to its wrong satellite calls (1 of 12 right from the
  FASTA).
- **Approach.** A calibration that does not come from the genomes the rules were made with: leave one lineage out over the
  development benchmark and report it beside each level. The CRM core's own enrichment does not separate right from wrong core
  calls (pepper: wrong cores 4.0-18.1-fold on their chromosomes, right ones 3.5-12.4; *O. australiensis*'s two right cores 2.9
  and 4.4), so weighing CRM by its core fold (#46) would not fix pepper; Hi-C does (pepper with `tes/` and Hi-C: 11 and 12 of
  12). A second line, not a stricter LTR rule, is what the data point to.
- **Checked in 0.19, no rule.** *Core age:* the signature's copies inside each core against the same families' copies elsewhere
  on its chromosome (one-sided Mann-Whitney, P ≤ 0.01): where young copies concentrate in the cores genome-wide (the premise
  of `ageline`), CRM cores were right 167 of 174 times when young and 62 of 73 otherwise, but pepper's wrong cores are as
  young as its right ones (CRM LTR-core-alone calls where the premise holds, one per chromosome from its `tes/`-only run:
  young 14 of 18 right, not young 4 of 9). *The FASTA's family window* (the
  same genome's FASTA-only run): where it lies beyond both readouts' resolution from an LTR-core-alone call, the core was wrong
  4 of 4 times (*Capsicum* Zunla-1 Chr08, LA1777, both *Sesbania*); where it agrees, 90 of 107 right. A lead: the family line
  does not run where calls are confident, and its test may not oppose a call.
- **Settle:** `medium` on the added genomes as reliable as on the development genomes (88%), with no right confident call lost
  on either. **Code:** `sat/signature.py`, `sat/centromere.py`, `fuse.py` (`_confidence`).

## Satellites

### 1. Calls on chains of partial copies

Partly done. A family can concentrate at one place without being tandem (under half its bp in tandem, under half in whole-unit
hits): a chain of pieces of a longer repeat, usually a retrotransposon (*Glycine dolichocarpa*: 44 whole hits of 25,732, 65% inside
SIRE; 10-Mb calls). Such calls are `partial_chain` (design.md), a cap that rests on a principle, not on truth. On the panel 41 such
families in 29 species carry 409 of 2,355 calls. No benchmark genome has a dominant chain; its chain calls were 29 right, 4 wrong
of 33 with truth (0.12), `chain_to_core` 7 of 7; *Ae. geniculata* with `tes/` elects a 1,216-bp chain family inside CRM cores (13
of 14 right, none `high`). Never measured: whether short (< 1 kb) single-sequence families lose diverged copies (`min_hmm_seqs` 20).
**Settle:** CENH3 for a genome whose chain dominates without a CRM core (6 cotton species: 23-26 of 26 chromosomes). **Code:**
`sat/classify.py` (roles, LOCUS_FOLD), `sat/centromere.py` (`partial_chain`, `chain_to_core`).

### 2. The 40-family cap

Open. At most 40 clusters per complex are modeled; the dropped ones are listed (`satellite/dropped_clusters.tsv`). Over 71 runs that
dropped any: 2,376 clusters, 1,553 on one chromosome; only 4 reach 20 kb on a chromosome without a call (all *B. carinata*), 349 lie
on a chromosome with a `low` call. No case yet where the cap costs a centromere. Changing the family set is not neutral: overlapping
hits go to the higher-scoring family (`resolve_overlaps`). **Settle:** *B. carinata* without the cap (roles, main families,
arrays, calls); if a genome shows a cost, a budget priced from `search_cost.tsv` instead of a count. **Code:** `sat/families.py`
(`max_families`).

### 38. Which array holds the kinetochore

Open. Where a chromosome holds several candidate arrays, the sequence line takes the largest cluster; no tested line can choose
another.
- **Bornean orangutan CM054652.2** (`high`, 14.3 Mb off; the only confident ape miss): CENP-A on 4 small alpha-satellite arrays (73
  kb), the call on a 3.4-Mb cluster. Hi-C's peak lies 0.26 Mb from the CDR but is `low` (1 of 3 engine tests) and cannot oppose
  `high`; the called array dips in CpG too. The truth may be haplotype-specific (CenSat v2.0 marks the called array active; the CDR
  region was redrawn in 2025; hap2's DiMeLo-seq shows CENP-A at both sites).
- **Bonobo CM055489.2** (`low`, wrong): the true array is the reported competing cluster, but only a strong Hi-C peak may choose it
  (of the two weaker peaks that rule kept out, this one would have been right and CHM13 chr21's wrong).
- **Mouse**: the AT trough chooses the minor satellite (`at_chose_alt`, 18 of 20 from the FASTA, `low`; #70), a rule shaped on mouse
  alone. Chr2, 5 and 14 share one byte-identical end (9,259-807,536).
- **Approach.** Settle the orangutan truth (2025 CDR, CenSat v2.0, both haplotypes). A second line for a weak Hi-C peak against
  `high` (a CpG dip at the peak; CpG alone chose badly: design.md, Choices tested and rejected). For two satellites at one end, a
  sequence feature of the CENP-A satellite (CENP-B boxes in mouse; an idea).
- **Settle:** the orangutan, bonobo and mouse calls; no right ape call lost; CHM13 chr21 stays right. **Code:** `sat/centromere.py`
  (`satellite_call`, `rival_bp`), `fuse.py` (`hic_chose_alt`, `hic_alt_ambiguous`), `cdr.py`.

### 39. Satellite arrays below the engine's absolute sizes

Partly done. One prior, S_MIN = 20 kb, scales the satellite engine's size rules in every genome (the "prior" fields of
`sat/params.py`, 5 kb to 5 Mb), and units under 30 bp or over 30 kb are not searched (`min_period`, `LU_MAX_UNIT`), so the
satellite line places 0 of 14 *C. neoformans*, 0 of 8 *V. dahliae*, 0 of 21 *Z. tritici*. Since 0.16 the repeat-density line
places the first two from the FASTA (#74). In 0.10.1's comparison (readouts picked by the truth) TRASH2 placed 7 of 14 *C.
neoformans* with 0 wrong and satellome's largest tandem region 10 of 14 *C. neoformans* and 7 of 8 *V. dahliae*. **Approach.** A
family-free tandem readout tested like the AT troughs (each chromosome's densest tandem region, accepted where the N densest fall
one per chromosome beyond chance). **Settle:** *Z. tritici* placed; no change where satellites are large; shuffled and holocentric
genomes fail. **Code:** `sat/params.py`, `atline.py`, `sat/karyotype.py`.

### 41. Right calls held at low by the honesty cap

Open. Where no family dominates, every satellite call is `low`, even when mostly right: chicken 27 of 33 (chromosome-specific
satellites), moss 16 of 26 from the FASTA, *Chlamydomonas* 6 of 17. Lifting the cap where the one-per-chromosome test passes was
rejected (chicken +27 right, +5 wrong; pig 8 and 9). Dominance itself is unstable (a fresh 0.8 run split pig's satellites
differently: FASTA-only pig 13 → 8 of 17), and nothing measures how far a small change of the family set moves it. A bird lead: 22
of chicken's 50 interstitial telomeric arrays lie within 100 kb of a CENP-A window (1.4 expected, P 4e-21; Huang 2023); as a line
only chicken beats chance (6 right, 0 wrong), so it needs a second bird T2T genome. **Settle:** a stability test of dominance
(resampled seeds, families pooled by kind); chicken's right calls rise and no wrong one does. **Code:** `sat/centromere.py`
(`no_dominant_satellite`), `telomere.py`.

### 42. The one-per-chromosome test near its level

Open. The test passes on any family that happens to lie once per chromosome (design.md, Known risks), and genomes near its level flip:
- **Potato** fails from the FASTA (7 of 12, P 0.87) but passes with `genes/` (10 of 12, P 0.044), whose gene veto changes the family
  set, so two wrong calls become confident (chr04 3.4 Mb off, chr08 1.7 Mb off): an input that should only help adds 2 errors.
- ***S. chilense*** passes (10 of 12, P 0.045) on mostly non-centromeric arrays; only the honesty cap keeps its calls `low`.
- **Oat**: a 335-360-bp family on every chromosome is not centromeric yet would pass (#35).
- **Siamang** fails (18 of 26, P 0.2): its arrays form two clusters over 500 kb apart on 8 chromosomes.
- **Few chromosomes**: no power at 4 or fewer; the acrocentric test needs about 16; no small acrocentric karyotype tried.
- **Approach.** Report the test's P beside each call that depends on it; measure how the outcome moves with the family set (with and
  without the gene veto; families pooled by kind). **Settle:** no confident error in potato with `genes/`; oat's FASTA-only calls
  against CENH3. **Code:** `sat/karyotype.py`, `sat/genomes.py` (the gene veto).
- **Six chromosomes** (*Sesbania bispinosa*, 0.17): one 139-bp family lies on every chromosome and is centromeric on all 6, but
  two chromosomes hold a second large array, so the 6 largest clusters fall on 4 chromosomes (null 4, P 0.71) and every call is
  capped `low` (6 of 6 right in both lines). With 6 chromosomes even 5 of 6 gives P of about 0.2: the test cannot pass a small
  karyotype unless every chromosome is perfect. With `tes/`, weak CRM cores then replace 3 right `low` calls (SbTA02 chr1 and
  chr5, SbSA44 chr1; #46).

### 43. Assemblies that lack or misassemble the centromeric array

Open. Where the array is missing, a non-centromeric satellite elsewhere can be called with confidence (*B. napus* A06/A09, `high`);
gaps are reported (`call_gaps`) but change nothing. The 56 ape centromeres misassembled in NCBI v2.1, excluded from the truth and
with published repair intervals, are free positive controls. A collapsed array makes spanning pairs farther apart than their
coordinates, mimicking the cis dip (#29); mouse chr2, 5 and 14 share one templated end (#38). **Approach.** Report per call what an
incomplete centromere looks like (a gap in or beside the call, a short array, a call away from its karyotype position).
**Settle:** the ape repairs flagged; T2T assemblies not. **Code:** `morphology.py`, `interpret.py`.

### 44. Holocentric genomes: one low call per chromosome

Partly done. Holocentric genomes get no confident call (0 of 196 chromosome-runs) but still one `low` "centromere" per chromosome.
Since 0.16 the report says when no call is confident and no line passed a genome-level test; nothing yet says the genome may be
holocentric. A flag from `no_dominant_satellite` was rejected (it fires on 25 of 26 siamang and 25 of 32 EquCab3 calls too); the
one-per-chromosome test fails there (P 0.13-0.73) and in many monocentric genomes. Never scored: the satellite line on satellite
holocentromeres (*R. breviuscula*, *Chionographis japonica*, *Luzula sylvatica*: #35). An LTR signature can form at the ends (*R.
breviuscula*: Athila, 5.8-fold, cores within 1 Mb of an end on 4 of 5): a non-CRM signature above 10-fold or a CRM seed would place
calls at the ends from `tes/` alone. **Approach.** Several independent signs, each tested on the three genomes with CENH3: no
dominant satellite, a failed one-per-chromosome test, Hi-C's holocentric signature, CENH3-domain-like arrays along the chromosome.
- **Five more holocentric controls** (0.17, [benchmark.md](benchmark.md#genomes-added-in-017)): *Chionographis japonica*, whose
  CENH3 sits on Mb-scale satellite arrays, *Luzula sylvatica*, on satellite arrays of about 100-200 kb and satellite-free
  domains, *Leptidea sinapis*, *Ascaris suum* (germline) and *Cuscuta epithymum*: no confident call on 97 chromosome-runs. The
  satellite line on *Chionographis* and *Luzula* stays `low` (12 and 6 calls). The repeat-family test passes in *Leptidea* (the
  16 calls it places are `low`), so a pass does not rule out holocentry.

### 78. Centromeres away from the genome's main satellite

Open (0.17). In *N. benthamiana* (NbT2T, CENH3 ChIP of the assembled line; [benchmark.md](benchmark.md#genomes-added-in-017))
CENH3 sits on arrays of the 43-bp satellite, the genome's main centromeric family, on 8 of 18 chromosomes. On the other 10 it
sits on other satellites (a 181-bp array on chr4, 143 kb inside the CENH3 domain; a 1-kb-unit array on chr11, 93 kb; a 31-bp
array on chr15, 39 kb) or on at most 5 kb of any satellite (chr6, 7, 10, 12, 13, 14, 19). cenfuse calls arrays of the main
families, right on chr7 and chr14 (an array within 1 Mb) and wrong on the other 8: confidently on 6 chromosomes from the FASTA
and on 8 with `tes/`, 14 of the 29 confident misses on the genomes added in 0.17. On 5 of the 8 (chr4, 6, 10, 11, 19) the called
arrays are of the 43-bp satellite, and the CRM signature (55-fold) sits with them, so the two lines agree on the wrong place
(#46); on 3 (chr12, 13, 15) they are of 164- and 87-bp families. A rule that prefers the genome's dominant family cannot see a
centromere on a minor one or on none. Cotton CM104175.1 is a related case: satellite, CRM core and both Hi-C readouts place it
4.2 Mb from the line's CENH3 domain (86.8-fold), `high` with Hi-C.
- **Approach.** A line that reads function rather than sequence: Hi-C where its tests pass (NbT2T's Hi-C is in GSA only, 180-250
  GB), or CpG or CHG methylation of the same plant (#16). Checked in 0.19 from the FASTA: nothing in the satellite summary
  separates the right calls from the wrong ones (main-family bp, cluster size, monomer identity to the family model 0.835-0.925,
  HOR fraction, spacing); chr4 and chr10 carry 1.9-Mb arrays of the 43-bp family like the right chr8, chr14 and chr17, and
  the most model-like array of each chromosome lies within 1 Mb of CENH3 no more often than the largest (8 of 13 each).
- **Settle:** no confident call on those 8 *N. benthamiana* chromosomes, and no right satellite call of the development
  benchmark lost. **Code:** `sat/centromere.py`, `fuse.py`.

## LTR retrotransposons

### 4. The decisive CRM test with subgenome-specific lineages

Waiting for data. CRM core counts differ between subgenomes in 3 of 87 species at P ≤ 0.01 (0.9 expected), and the marker is decisive
in 2 of them anyway (*B. napus* A 8/8, B 2/8; *Euphrasia confusa* 10/11 vs 3/11); at P ≤ 0.05, 9 species (4.4 expected, P 0.03), 5
still decisive (*Cerastium alpinum* only flags: 11/18 vs 0/18). No evidence yet that a decisive marker is often lost. **Settle:**
those species against CENH3 (#35), random splits of diploids as the null; then test CRM concentration within subgenomes.

### 10. The LTR retrotransposon line

Partly done. Open questions on the branches that move, overrule or lower a call (design.md, What a decisive core does):
- **Which rule earned which gain.** No benchmark-only switch per rule attributes the gains.
- **The head-to-head test** (a confident call is replaced only if the core is 10-fold enriched and holds significantly more of the
  signature's elements than the call; otherwise `te_discordant`): its cited cases (*B. rapa* A09, soybean chr1/chr6) have not arisen
  since 0.13. In 0.15 it keeps two right calls with truth, grape chr11 (`high`, HT and HTH arms) and *Ae. geniculata* chr05M (`low`).
  Replacing a `low` call needs no such test (of the moves it would have blocked, 2 right, 1 wrong: soybean chr13, #27).
- **`weak_core`** cost *S. verrucosum* 2 right satellite calls (cores of 2 and 4 elements, chr10 and chr12) and is right on 7 of 15
  tomato chromosomes: no confident error, but right `low` calls lost. **`weak_satellite_beside_core`**: one confident error (*S.
  verrucosum* chr06, 1.4 Mb), tomatoes 7 of 7.
- **`te_concordant`** is a weak flag (#45).
- **Settle:** the per-rule switch replayed on every `tes/` arm and on the panel, each branch's right and wrong counts reported.
  **Code:** `sat/centromere.py`, `sat/signature.py`.

### 27. Wrong LTR cores that look right

Open. A signature learned from the other chromosomes can mark where centromeric families accumulated beside the centromere, and no
statistic of the signature tells those cores from the right ones.
- **The tomatoes' own Hi-C (0.16)**: the 13 tomatoes (CNCB CRA018974, 20-31 M mapped pairs each) and *S. verrucosum* (Omni-C, 56 M
  pairs). Hi-C's prior-free peaks are right on 21 of 161 chromosomes, and the cis and trans readouts meet on most chromosomes only
  in TS-60. With `tes/` the 14 *Solanum* genomes go from 123 to 130 of 168 right and 15 to 14 confident errors: TS-60's cis signal
  passes both tests and moves its 2 confident LTR errors onto CENH3 (`cis_trans_pair`), TS-9's lowers its 2, and Hi-C anchored on
  wrong cores adds 2 (LA0111, LA2134). An independent line exists only where the library shows the centromeres
  ([benchmark.md](benchmark.md#the-tomatoes-own-hi-c)).
- **Others.** Soybean chr13: a CRM core 1.8 Mb from CENH3, `medium`. Potato chr08 (#45). *Arabidopsis*: Athila (13.3-fold) has its
  cores 0.46-1.38 Mb beside CENH3 on all 5 chromosomes and only flags because CEN180 passes the one-per-chromosome test. In 6 of 27
  polyploids above 10-fold cores and FASTA-only calls mostly disagree (no truth).
- **Tried in 0.15, none separates them:** seed specificity, the core's share or fold, a second comparable peak (accuracy 0.85 with
  one, 0.83 without; reported as `core_alt`), a super-majority of the core's share; fragments as a locator (#67).
- **Leads.** A core beside the dominant satellite array rather than over it (*Arabidopsis*); Hi-C anchored on the cores.
- **Settle:** non-CRM signatures with truth beyond *Solanum* (#35, #72); *Arabidopsis* with CEN180 masked must not move onto the
  Athila cores; randomized tomato Hi-C silent. **Code:** `sat/signature.py`, `hicfuse.py`.

### 31. A centromere that left its satellite and its LTR families

Built, unchecked (no CENH3). *P. annua* Pa1A: Hi-C (peak 168.0 Mb; the test passes on 14 of 14) and earlier analyses put its
centromere at 163.0-169.4 Mb, 21 Mb from its largest satellite, around which the CRM families learned from the other chromosomes
concentrate (145.1-148.3 Mb). The signature's copies are younger at the Hi-C site (median 0 vs 620,000 years, P 2e-13; within one
CRM family 0 vs 157,000, P 2e-7), as young CRM colonize present maize centromeres and older ones their flanks (Wolfgruber et al.
2009). The age choice between density peaks makes the young peak Pa1A's core (`high`, Hi-C agreeing); `moved_with_age` (0.13) acts
on no benchmark run and has never been scored against truth. **Limits.** The premise is tested on the other chromosomes (Fisher P
5e-10 here); the young copies' core alone is a worse locator than the density core (128 vs 155 of 172 with truth); a younger
family can mimic a move (#71). **Settle:** CENH3 of *P. annua*; genomes with known evolutionarily new centromeres and Hi-C (horse
chr11, #50; the orangutan chr10 neocentromere polymorphism, de Gennaro 2026); within-species centromere tables (#35). **Code:**
`sat/signature.py` (`age_choice`), `ageline.py`, `fuse.py`.

### 45. Confident calls at the edge of a wide LTR core

Open. A call agrees with a core when the gap is at most the wider of the two (`te_concordant`), and an agreeing core lifts both the
honesty cap and the one-per-chromosome cap. Of 130 such calls with truth in the 0.13 FASTA + `tes/` arms, 125 are right and 5
confidently wrong, 5 of the 7 confident errors of those arms. In 3 of the 4 below the core holds CENH3 and the call sits at its
edge:

| chromosome | confidence | call (Mb) | LTR core (Mb) | CENH3 (Mb) |
|---|---|---|---|---|
| maize Mo17 CM039154.1 (chr5) | `high` | 107.19-107.40 | 103.54-108.83 | 103.5-105.2 |
| *Ae. tauschii* Chr5 | `high` | 199.36-201.42 | 194.02-201.03 | 195.4-198.1 |
| *Ae. tauschii* Chr7 | `medium` | 340.27-340.50 (45 kb of satellite) | 340.47-349.60 | 343.6-344.9 |
| potato chr08 | `high` | 15.54-16.92 | 16.46-17.41 (5 elements) | 18.6-18.7 |

The fifth is tomato CM175074.1 (its core misses CENH3: #27); *Ae. geniculata* chr05U has the same shape (2.1 Mb off, in a 9.3-Mb
core whose other edge reaches CENH3). Maize chr5 is the plants' only confident error from the FASTA alone. Where the core is at
least 3 times wider than the call, 11 right and 4 wrong; otherwise 114 and 1 (exploratory).
- **Tested, not adopted.** Agreement at the core's maximum (0.14): lowers 2 of the 5 errors and 3 right calls.
  `satellite_in_core` (0.16, an untested array that agrees with a decisive core gives way to the core): right on *Ae. tauschii*
  Chr5/Chr7 and *S. lycopersicoides* CM175231.1, but on the panel it changed 80 calls in 17 species and the pairs touched came
  apart (67 pairs: 26 closer, 40 farther), also when restricted to arrays overlapping the core.
- **Ties.** Most core profiles have tied maxima, so a core's position can depend on which equal maximum comes first (#75).
- **Approach.** Judge agreement by chance, not by gap (`linetests.py` computes each chromosome's chance of agreement).
- **Genomes added in 0.17.** *O. australiensis* chr2 is `high` from the FASTA and with `tes/`, 1.5 Mb off: a satellite array 1
  Mb outside the LTR core that holds CENH3 (the authors' neocentromere). Not acted on: where a call lies beside a decisive core
  (disjoint, within tolerance), the call is right more often than the core in the 0.16.1 release (13 against 10 of 15); the new
  genomes have too few such cases to say (2 against 3 of 3).
- **The 3-fold split holds out of sample (0.18 check).** On the genomes added in 0.17, which the 0.13 split never saw, satellite
  calls flagged `te_concordant` (every arm with `tes/`, 0.17 release runs) are right 17 of 24 times where the agreeing core is at
  least 3 times wider than the call, against 125 of 131 otherwise (development genomes in the same runs: 55 of 74 against 359
  of 370). A core merely wider than the call does not separate on the added genomes (98 of 107 against 44 of 48), and judging
  agreement by chance does not either (chance of agreement 0.02-1 for every such call; above 0.2 right 203 of 204), so no
  principled cut is in hand; the rule would sit in the satellite engine (`finalize`), which the stored FASTA-only replays
  read as is. A lead, not a rule.
- **A cut without a width: the call's share of the core's elements (0.19).** Over the 0.18.0 release runs, `te_concordant` calls
  holding more than half of their core's elements were right 322 of 329 times (development) and 118 of 125 (added genomes);
  those holding fewer, 92 of 115 and 24 of 30: the split replicates out of sample. Registered as `te_core_elsewhere` with
  criteria before its replay: the head-to-head test the engine already applies to an array beside a core (te.overrules:
  decisive, core ≥ 10-fold on its chromosome, the signature's elements significantly more in the rest of the core than in
  the array, binomial at 1/2, P ≤ 0.01) also judging an array inside or overlapping the core; where it contradicts, the array
  keeps its place and the confidence the satellite line gives it alone. Replayed (v20 c12c): on the development genomes it
  lowered 8 confident errors (*Ae. tauschii* Chr5 and Chr7, with and without Hi-C; LA2093 CM175086.1 twice; tomato
  CM175074.1; TS-60 CM175038.1) and 4 right confident calls (*Ae. geniculata* chr02M, chr03U; LA0426 CM175098.1 twice) to
  `low`, none to wrong, and a fifth right call went from `high` to `medium` through a changed Hi-C reference (LA0426
  CM175091.1 with Hi-C); on the added genomes only pepper G1-36576 Chr07, a right call (`medium` to `low` without Hi-C, `high`
  to `medium` with it): every other added genome's cores are under 10-fold; on the panel it lowered 16 confident calls in 10
  species (4 `high`, 12 `medium`), moving none. Two of its four criteria were missed as written, so it was not adopted: the
  share right among the calls it took out of the confident class had to fall below the confident class's in both sets
  (development 4 of 12 chromosome-runs; added genomes 1 of 1), and no call was to move, but a lowered call becomes movable by a
  tested Hi-C line (*Ae. tauschii* Chr5 with Hi-C, wrong to wrong) and changed Hi-C references moved 2 right pepper calls by a
  bin. Maize Mo17 chr5 is untouched either way (its satellite passes its own tests). A form that also judged long-unit arrays
  was replayed (v20 c12b) and narrowed after its results were read, for a reason of mechanism: without the core's support a
  long-unit array is no call, which left its chromosome without one beside its own decisive core (LA1982 CM175142.1 and
  LA2093 CM175082.1, with and without Hi-C: right `medium` calls became no call; with Hi-C, LA0103 CM175169.1 went from right
  `medium` to wrong `low` when the dropped call left Hi-C's references). The code is kept outside the repository (the v20
  development folder, patches/c2_engine.diff).
- **Settle:** the five calls fall to `low` or move onto CENH3; none of the 125 right calls drops below `medium`; CRM-plateau genomes
  not used to make the rule (wheat, the soybean lines: #35). **Code:** `sat/centromere.py` (`te_concordant`), `sat/te.py`
  (`agrees`), `linetests.py`.

### 46. CRM decides by its name

Partly done. CRM decides from 3 significant cores with no fold level; other lineages need 10-fold, specific cores (median core fold
≥ 10) and, to act on a confident call, a core with most of the signature's elements on its chromosome (`TE_MAJOR_SHARE`). CRM's
standing rests on every CRM core short of a majority lying on CENH3 (24 of 24), but those come from thinned annotations of genomes
whose CRM is centromeric (rice, potato, maize, *Brachypodium*). Non-centromeric CRM exists: tobacco (CENH3 binds Tekay; CRM clusters
beside), eggplant, soybean chr13 (#27); LTRquest labels 14 CHM13 ERVs CRM. The 3-core CRM seed was not retested against the
*Arabidopsis* case that sank an earlier form ([Tried and set aside](#tried-and-set-aside-not-in-designmd)). **Approach.** The same
test for CRM as for other lineages (a fold level, or the decisive test measured per genome), or CRM's standing only where its cores
pass the one-per-chromosome test. **Settle:** no right call lost in the CRM-decisive genomes; no tobacco or eggplant CENH3 domain
overruled by CRM (#35). **Code:** `sat/signature.py` (`seed_of`, `TE_SEED_MIN_CORES`).
- **Weak CRM cores** (0.17, *Sesbania bispinosa*): a CRM signature (11- and 21.5-fold in the two lines, decisive by its name)
  replaced right `low` satellite calls with cores 22 and 9.7 Mb away (SbTA02 chr1, chr5; SbSA44 chr1 likewise), `medium`; those
  cores are 2.8-, 3.2- and 2.9-fold on their chromosomes. Holding such cores to the 10-fold of other lineages was measured after
  this result and not adopted: cores under 10-fold that replaced a satellite were right 10 of 10 in the 0.16.1 release and 39
  against 8 for the replaced arrays in the new genomes. In *N. benthamiana* the CRM signature (55-fold) sits with the 43-bp
  arrays cenfuse calls on 5 of the 8 chromosomes whose centromere sits away from that satellite, and none of its copies lies
  inside CENH3 on any of the 8 (#78).

### 47. Centromeric retroelements the LTR line cannot read

Partly done. The LTR line reads classified LTR elements only. In *Cryptococcus* 89-94% of LTRquest's elements lie in the CENP-A
domains but none has a REXdb lineage (unknown, chromo-outgroup); since 0.16 the repeat-density line places *C. neoformans* from
the FASTA instead (#74). *Chlamydomonas*' centromeres are clusters of ZeppL, a LINE: LTRquest finds 669 elements, 613 without
lineage; 6 of 17 right, all `low` (12 of 17 since 0.18, where the repeat-family window takes the place of 6 wrong satellite
calls). Pooling unclassified elements qualified in fish genomes without truth (design.md); all LTR
elements pooled per chromosome (0.15) are right on 6 of 14 *C. neoformans* and 3 of 8 *V. dahliae* but never significant (too few
independent insertions). **Approach.** Lineage labels covering fungal Ty3 and chromoviruses; a non-LTR annotation read by the
family signature; each needs a second genome with truth. **Settle:** *Chlamydomonas* placed with `tes/`; fish and *Oryza*
unchanged; shuffled positions silent. **Code:** `sat/te.py` (reads `LTR_retrotransposon` only), `sat/signature.py`.

### 48. How much the LTR line depends on the annotation

Open. The decisive test needs significant cores on most chromosomes; potato sits at its threshold (8 of 12 in 0.13; a majority is
7). benchmark.md's annotation table holds an early lineage line's numbers (half of potato's CRM elements, EDTA or DANTE_LTR each
lose the gain from `tes/`; rice drops from 8 cores to 3, 1, 1 at 50, 25, 10%); the family signature (0.12 on) was never thinned.
LTRquest writes no version into its GFF3, so a run records only the `tes/` file's fingerprint. inputs.md marks `--proteins`
untested. **Settle:** the thinning rerun on current code with LTRquest only (v6 `bench/td_thin.py`), reporting where each signature
stops deciding; the recipe benchmarked with `--proteins`, or the flag dropped. **Code:** `sat/signature.py`, `teline.py`.

### 71. Age choice and family composition

Open (0.15 code review). `signature.age_choice` compares the ages of all recruited families' copies at two peaks; families burst at
different times, so a peak made of a younger family reads as younger with no move of the centromere. In 0.15 it moved 27 cores over
32 genomes with truth and `tes/` (12 wrong to right, 1 right to wrong); on the panel the 62 pairs it touches agree 32 times, as
before. **Settle:** the test within each peak's commonest family, where copies suffice; compare the moves. **Code:**
`sat/signature.py`.

### 72. Minority cores beyond *Solanum*

Open (0.15 code review). `minority_core` was built on the tomatoes (#46). Outside *Solanum* one scored chromosome: cucumber chr1 (HT,
HTH), right in 0.14, 2.1 Mb off since 0.15, still `low`. On the panel it acts on 14 calls in 3 species (*Arabidopsis suecica*,
*Malus toringoides*, *Sphagnum papillosum*); the 11 pairs they touch agree in 7 against 3 under 0.14's rules. It makes no
confident call. **Settle:** non-CRM genomes with truth outside *Solanum* (#35). **Code:** `sat/centromere.py`, `sat/signature.py`.

## Hi-C

### 13. Bird microchromosomes in Hi-C

Waiting for data. Chicken's 2.5-7-Mb microchromosomes stay out of the Hi-C line because its library (the first 40 M read pairs of
SRR15423486; 26.1 M unique) is too shallow for bins under 250 kb. The bin cap is not the limit: at `MAX_CHROM_BINS` 8,000 the depth
rule still chose 250 kb, the gate passed on 0 of 40 and no call changed. **Settle:** deeper chicken Hi-C, scored against the CENP-A
windows. **Code:** `hic/params.py`, `hic/bins.py`, `fold.py`.

### 14. The Rabl axis where the fold has no power

Partly done. The trans arm alignment passes in *C. neoformans* (anaphase, 12/14) and *M. sympodialis* (6/8), but the arm fold, which
must lean the same way first (the holocentric guard), has no power in these small nuclei (1/14; 0/8), so the alignment raises
nothing. In holocentric *C. elegans* the alignment reached 2/5 (P 0.095), and a position-matched null that rejects it also rejects
barley. G1 *C. neoformans* loses the cluster, yet anchored on CENP-A 8 of 14 chromosomes still align (P 0.0014). Other lines now
place these genomes (`at_pair`, the cis-trans line, `at_denovo_pair`, Hi-C alone: #52); with 6 chromosomes *C. tropicalis*' Rabl
test cannot reach P ≤ 0.01 (floor 0.0156). **Approach.** Partners placed on the confident sequence calls as a corroboration test
(a prototype: G1 8/14, anaphase 12/14, randomized 0/14), with a per-genome test like the gate's; a compartment-adjusted alignment
(the E1 sign at centromeres is opposite in plants and apes); the fold at the 5-kb rung. **Settle:** barley, *Neurospora* and the
yeasts pass; *C. elegans*, *Rhynchospora*, *Bombyx* and randomized Hi-C fail. More data: *Cryptococcus* relatives, *Toxoplasma*,
*T. brucei*, other *Candida*; Rabl tissue for grasses (*B. distachyon* root tips, not leaves; barley endosperm loses Rabl; Bd21
CRR611509, never used). The docs give no tissue or cell-cycle advice. **Code:** `fold.py` (`align_profile`), `fuse.py`
(`fold_pair`, `align_pair`).

### 15. Centromeres at the edge of a terminal unmappable array

Partly done. The which-end test (`hic_end`) raises a call at the end with the higher anchored score where that holds on most other
chromosomes (goat, sheep, horse, mouse; 54 of the 55 scored raises right). Open: chicken does not pass (12 of 19): its
microchromosomes' arrays are unmappable for 1-1.25 Mb at the very end, so their peaks are `telomeric` and the calls stay `low`. The
one raise not scored right, sheep CP162266.1 (X, `high`, 83.4 Mb from another caller's locus), agrees with the karyotype (CI 0.054,
acrocentric), so the other caller is probably wrong (#37). **Code:** `hicfuse.py` (`end_test`).

### 29. What the contact crossing a centromere is

Waiting for data. Piled up on CENH3, tomato's two flanks contact each other about twice as much as DNA at the same distance elsewhere
(within 0.25 Mb; 1.79 against 0.97 at random positions); potato shows a smaller central peak between deep dips. Both lose contact
with the arms out to 6 Mb, a mean effect only (per chromosome the lowest arm contact is not the centromere: tomato 0/12, potato
0/11, rice 0/12). *S. pombe* has heterochromatin-dependent inter-arm contacts near its centromeres (Mizuguchi et al. 2014). A fold
of the flanks or an insulated domain could explain all of it: a hypothesis.
- **An artifact first.** Pairs spanning a collapsed or missing array are farther apart than their coordinates, a short-range
  depletion that fades with distance; the dips (potato, *Verticillium*) were never compared with the gaps and collapsed arrays.
- **Chance** comes from simulated maps; on moss's re-paired control the dip lies at the LTR cores on 11 of 26 (about 5 expected, P
  0.011, just short; about one in five of 20 control tests gets this low). If the simulation runs short, rebalance each map.
- **Settle:** the profile at several box sizes; Hi-C of heterochromatin or condensin mutants, or other cell-cycle stages, in a
  genome with the signal. **Code:** `hicfuse.py` (`cis_profile`, `cis_chance`).

### 30. A Hi-C test with too few confident references

Built, unchecked. When the confident references are fewer than the gate needs (4), the `low` calls join them, agreements counted at
clean peaks only (`hic_gate_refs` `mixed`). It acts only on cucumber 9930 with `tes/` (Hi-C passes on 3 of 7, no call changes). The
rule for a genome without confident calls (`low`) acts only in the siamang (no CENP-A truth), where it raises the two calls that
confidently disagree with the satellome locus (CM054514.2 by 1.0 Mb, CM054529.2 by 15.0 Mb). **Settle:** a genome where it acts,
with truth; the siamang's DiMeLo-seq track (#35). **Code:** `hicfuse.py` (`gate`).

### 32. Calls in a terminal eighth that an interior line contradicts

Built in 0.16 as `cis_end_pair` (`low`, a named heuristic); right on 4 of 4 in two genomes. The cis-trans line is read outside
the terminal eighths, so it had no say on a call inside one: tomato without `tes/` kept 2 wrong `low` calls 8.6% and 11.8% from an
end (CM175073.1, CM175075.1; `cis_blind`) while their interior extremes lie on CENH3. The rule takes the chromosome's own trans
pair: the test's trans peak over the whole chromosome (`hicfuse.anchored_on_cores`, clean peaks only, so an end never counts) at
the interior extreme. Both tomato calls move onto CENH3 (extremes 27.00-27.25 in CENH3 26.1-27.5 Mb, 20.25-20.50 in 19.6-20.9
Mb), and TS-60 with its own Hi-C moves its 2 such calls onto CENH3 too; an acrocentric chromosome among metacentric ones has no
pair and keeps its call. **Settle:** a third genome with the case (rule 3), before the rule makes confident calls; goat, sheep
and mouse keep their calls. **Code:** `hicfuse.py` (`cis_trans_line`), `fuse.py`.

### 49. Hi-C cannot see inside megabase arrays

Open. At about 100 M pairs and MAPQ 30, the CENP-A domain lies in bins without a uniquely mapped contact on 139 of 157 T2T mammal
chromosomes (blind blocks of median 2.5 Mb in chimpanzee to 19 Mb in gorilla); the anchored peak overlaps CENP-A on 5 of 157, a
median 3.2 Mb away at the blind block's edge. A passing test certifies that Hi-C sees the array, not the kinetochore. MAPQ ≥ 1
changes nothing (near-identical copies map at MAPQ 0). **Settle:** remap one ape library keeping multi-mapped reads, or through
array-unique k-mers, and measure what is recoverable; Pore-C if public data exist for a T2T genome with CENP-A truth. **Code:**
`hic/io.py`, `hic/bins.py`, `hicfuse.py` (`visible_window`).

### 50. Satellite-free centromeres in mammals

Open. Horse chr11 (CM128059.1), the classic evolutionarily new centromere and the benchmark's only satellite-free mammal centromere
with truth, gets no call in any arm (nor EquCab3's chr11 and chr31). The rule for its case exists (`hic_medium_peak`, a `low` call),
but horse Hi-C fails its test at every depth; the library comes from another horse, and chr11's CENP-A slides about 500 kb between
horses. The cis signal works in no mammal. A gap-rank test (does the call's gap have the top flank score among its chromosome's
gaps?) reached P ≤ 0.01 in horse (9 of 11), pig (8/11) and rhesus (13/16), where the gate fails, and 0 of 14 on randomized Hi-C;
never adopted. **Settle:** the same-animal horse set and the multi-horse CENP-A union (#35), then the gap-rank test with a stated
role (a `low` call where there is no sequence call); other mammals and controls unchanged. **Code:** `fuse.py`, `hicfuse.py`
(`gate`).
- **Donkey** (0.17, T2T, CENP-A of two other donkeys, [benchmark.md](benchmark.md#genomes-added-in-017)): 14 of 31 right from
  the FASTA, all `low`, no confident error. The satellite calls on its satellite-free centromeres are wrong, but
  `no_dominant_satellite` and the failed one-per-chromosome test keep them `low`.

### 51. Hi-C tests at their threshold

Open. The tests flip rather than fade with depth: both orangutans pass at full depth and fail at half (the bin doubles to 500 kb); pig
fails at 100% and 50%, passes on 20 chromosomes at 25% and 10% (500-kb and 1-Mb bins), fails at 5%, its calls unchanged; rice passes
at 50% and fails at 25%. A verdict that improves when reads are discarded suggests the bin decides as much as the signal. The depth
series comes from older code; the share of about 8 random halves that pass was never tried (a split-half flag was rejected:
design.md). **Settle:** the depth series rerun on current code; the gate's chance model checked across bin sizes at fixed depth;
the share of passing halves reported as depth guidance, not a threshold. **Code:** `hic/bins.py`, `hicfuse.py` (`gate`).

### 52. Hi-C alone, untested

Partly done. In a genome without sequence calls, a clean `high` prior-free peak becomes a `low` call that nothing tests; the
strongest trans cluster may be something else, and a holocentric genome without sequence calls would get such calls (design.md,
Known risks). A peak held at `medium` only by `null_untestable` counts too (*C. tropicalis* 0 → 6 of 7; on randomized libraries the
flag never occurs without `not_converged`, which keeps a peak unclean).
- **Sharp single-bin peaks**: potato's 7 at 12-27 times the chromosome median are `high` and wrong; one carries `hub`.
- **An effect-size readout** (median peak over median, ≥ 2) separated *C. elegans* (1.21) and randomized libraries (1.05-1.51) from
  yeast (5.96) and *Plasmodium* (5.40); never tried as a guard.
- ***Z. tritici***: 19 of 21, all `low`; 8 of its 21 chromosomes are accessory and may dilute every "most chromosomes" test.
- **Dormant controls** that `cenfuse run` never runs: the read-class split (`--repeats`: trans contacts of reads in repeats against
  unique ones; passed in 6 of 6 mammals in the first benchmark) and `simulate` (cis-only pairs from the assembly).
- **Settle:** the effect-size guard, a single-bin flag and the read-class split by default (satellite arrays and LTR cores as the
  repeat class): holocentric and randomized libraries give no Hi-C-alone call, the yeasts and *Plasmodium* unchanged, *Z. tritici*
  scored on its core chromosomes. **Code:** `fuse.py` (`hic_denovo_only`), `hic/qc.py`, `hicfuse.py` (`repeats=None`),
  `hic/simulate.py`.

### 53. Where the cis-trans line passes, every extreme is trusted

Partly done; the diagnostic is open. Where it passes, the line acts on every chromosome, not only where that chromosome's own trans
peak meets its extreme. Some extremes are far off (*N. crassa* 911 and 82 kb on 2 of 7, *V. dahliae* 112 kb on 1 of 8); no call
changed only because those chromosomes had Hi-C or AT calls. Anchoring at the other chromosomes' extremes once agreed in rice on 8 of
12 (P 2e-10), none right; the line now skips terminal eighths (where rice's extremes lay) and runs only where the gate fails, so it
never ran on rice, soybean or the apes. Acting only at pairs (tested in 0.14) would lose 6 right calls and fix none. Since 0.16 a
cis line at `support` raises the calls it agrees with and moves none (*S. verrucosum* with its own Omni-C, the first at `support`:
its 3 moves under the old rule were all wrong). `cis_agrees` at `support` raised 2 wrong calls there: the cis sign and test come
from the same LTR cores, so at that tier the chromosomes that agree are the ones that made the test pass (a lead). **Settle:** the
line as a diagnostic on gate-passing genomes with truth: rice, soybean and the apes pass nothing; the off extremes would not act.
**Code:** `hicfuse.py` (`cis_trans_line`), `fuse.py`.

### 54. The support tier's evidence

Open. The tier lets Hi-C that finds centromeres beyond chance, but not on most chromosomes, raise the calls its significant peaks
agree with and choose a competing cluster. Its evidence when adopted (agreeing peaks right 61 of 64 times) counts agreement with
another caller as right in goat (12/12), rhesus (8/8) and giraffe (4/4); goat's 8 raises to `high` stand on that caller. The
pre-registered bar (≥ 90% of raises right) was missed: 14 of 16, 88% (pig's two `high` misses). Its risk, raising a wrong call on a
non-centromeric trans cluster (potato: 1 of 64), is not stated in design.md. **Settle:** its raises scored against functional truth
only (cattle and goat CENP-A first: #35), ≥ 90% right or its role narrowed; the risk stated in design.md. **Code:** `fuse.py`
(`hic_support`).

## AT troughs

### 28. Where AT troughs do not work

Open. The test passes in *P. falciparum*, *M. sympodialis*, *N. crassa*, *V. dahliae* (deepest trough right on 13/14, 8/8, 7/7, 8/8)
and mouse, on 3.9% of 3,336 GC-shuffled genomes and on no holocentric genome. Where it fails:
- ***C. neoformans*** (P 0.36): its deepest 10-kb troughs are right on 12 of 14 but tie with other AT-rich loci (GC 0.423-0.452 vs
  0.43-0.46), so nothing in the depths sets a class apart. With Hi-C the cis-trans line (G1) and `at_denovo_pair` (anaphase) place
  it, and since 0.16 the repeat-density line from the FASTA (#74). *S. pombe*: `untested` with 3 chromosomes, its troughs right.
- **Rightly failed**: *Z. tritici* (P 0.69; 0 of 21 troughs on the centromere), *C. tropicalis* (P 0.15; 0 of 7).
- ***Chlamydomonas*** no longer passes (P 0.09): its ends make two deep troughs per chromosome (21 of the 30 deepest in a terminal 5%).
- **Latent false passes**: *B. rapa* TC1 (P 0.002), tomato LA0426 (P 0.04); the line acts there only on a chromosome without a call.
- **Genomes added in 0.17**, scored within 25 kb ([benchmark.md](benchmark.md#genomes-added-in-017)): *C. deuterogattii* 14 of
  14 from troughs (the test passes). *P. carinii* passes and its troughs were right on 9 of 17 until `at_interior` (0.17)
  stopped troughs in a terminal 5% placing calls alone: 14 of 17. In *P. murina*, *K. pini*, *C. albicans* and *Toxoplasma* the
  test fails and their troughs lie on the centromere on 6 of 17, 0 of 11, 0 of 8 and 2 of 13 chromosomes, so cenfuse rightly
  places nothing; Hi-C places *C. albicans* 8 of 8 and *Toxoplasma* 12 of 13 (`low`).
- **Settle:** a second feature beside the trough (a gene desert; without `genes/`, low k-mer complexity) so that *Z. tritici*, *C.
  tropicalis*, *Chlamydomonas*, holocentric and shuffled genomes fail while the passing fungi still pass. **Code:** `atline.py`.

### 70. The AT trough choosing between two satellites rests on one genome

Open (0.15 code review); a named heuristic (design.md, AT troughs). `at_chose_alt` lets a tested trough choose a `low` call's
competing cluster, though two satellites differ in composition whatever their centromeres do. Only T2T mouse exercises it (17
chromosomes in 3 arms, all right, all `low`); it never fired on the panel. Of the 136 scored `low` calls with a competing cluster,
the 2 outside mouse where exactly one cluster is right and the trough lies on exactly one (*S. tuberosum*, *S. verrucosum*, whose AT
test fails, so the rule does not act) have the trough on the wrong one; outside mouse the first choice is right on 27 of 46 such
chromosomes. It stays because the call stays `low`, the flag names it and the first choice stays reported (`seq_start`, `seq_end`).
**Settle:** a genome with two satellites at one centromere per chromosome, CENH3 or CENP-A truth and a passing AT test. **Code:**
`fuse.py`, `atline.py`.

## Repeat density and repeat families

### 74. The repeat-density line

Built in 0.16 (`repline.py`; design.md). It passes in *C. neoformans* (14 of 14 on CENP-A), *V. dahliae*, *Chlorella* and
*Chlamydomonas*, every peak right, and places only chromosomes nothing else places (`rep_only`) or pairs with a tested Hi-C call
(`rep_pair`).
- **Calibration.** Block shuffles passed the test 6.5% of the time at 0.05 (2,100 of 8 genomes), 1.3% at 0.01 (3.3% in *B. mori*).
- **Large genomes.** No large benchmark genome passes; 24 of 106 polyploids do, 20 at 500-kb windows with pericentromere-scale
  peaks (0.7-3.5 Mb). Replacing the `low` calls it contradicts was right on 14 of 14 benchmark cases (13 in *Chlamydomonas*) but
  took the panel's homoeologous calls apart (36 pairs: agreeing 21 → 11): not adopted. In *Malus toringoides* it placed 4
  chromosomes that had no call; the 2 homoeologous pairs they form disagree.
- **Speed.** Since 0.18 both repeat lines read the FASTA on the run's threads (`repline.pieces`; outputs byte-identical at any
  thread count): on *E. rockii* (2.0 Gb), 0.17's single-thread code took 206 s (repeat density) and 233 s (repeat families),
  0.18 at 16 threads 34-42 s and 42-51 s, at no higher peak memory (one thread: about 0.4 GB more). What remains serial is
  the genome-wide sort after reading (`np.unique`, a stable `argsort`), about half of the 16-thread time. On a loaded shared
  server at 4 threads per run, stage times in the release replays followed the load more than the code.
- **Settle:** a large genome with truth where it passes. **Code:** `repline.py`, `fuse.py` (`rep_only`, `rep_pair`).

### 76. The repeat-family line and low calls it contradicts

Partly done (0.18: `fam_over_low`). Built in 0.16 (`famline.py`; design.md, where its Results paragraph holds the numbers). Its test
is not calibrated, so since 0.16.1 it places only `low` calls (`fam_pair` removed).
- **What it finds.** Release runs: the line runs in 39 FASTA-only genomes and passes in 18, peaks right on 183 of 211
  chromosomes (13 of the 16 *Solanum* genomes, *Ae. geniculata* 14 of 14, *Chlamydomonas*, *Chlorella*, *N. crassa*, *V.
  dahliae*). The development scan of all 71 FASTA-only genomes passes 32 (also *B. rapa*, soybean, maize, rice, sorghum,
  where confident satellite calls keep the line from running), peaks right on 323 of 363. Holocentric genomes, mammals and
  chicken fail.
- **The move, 0.18.** A `low` satellite call that no other rule changed, outside the terminal eighths, moves to the window and
  stays `low` (`fam_over_low`, `fuse._fam_alone_against`) when the satellite could not vouch for it (`not_one_per_chromosome`,
  `no_dominant_satellite`, or the weak-satellite fallback), the window holds a cluster (at least `famline.MIN_COPIES` family
  K-mers, `fuse._fam_speaks`, which also gates `fam_only` and `fam_elsewhere`), the window lies farther from the call than the
  wider of the two intervals (`sat/te.agrees`, the rule for an array and an LTR core), and no other line speaks to the call on
  the chromosome (a Hi-C reading acting there, an agreeing LTR core, tested AT trough or repeat-density peak, centromere-like
  methylation, a CpG line); any other window off a `low` call only flags it (`fam_elsewhere`). Registered with its criteria
  before the replay, after the counterfactual on the stored runs was known; the resolution rule answers C1's miss in *Populus*
  (below), so it is a rule changed after a result was seen. The guards came after the first replays: a code review found the
  first form also moved low LTR cores (`minority_core`, `weak_core`) and calls other tested lines spoke to (fixed from the
  rule's own reason); the satellite-test guard came from designing for genomes outside the benchmark (it changed no call with
  truth); the cluster guard from the holocentric *Leptidea*, whose test passes on windows of 0-2 family K-mers (a guard on
  the window's rank among the genome's strongest, tried first, was decided there by ties among empty windows and cost a right
  pepper move: withdrawn). Replayed (development, panel, added genomes): development 123 calls moved (122 with truth), 94 wrong
  -> right, 2 right -> wrong (LA0111 CM175153.1 without `tes/`, with and without Hi-C), no confident or control call changed;
  added genomes 18 moved, 11 wrong -> right (pepper 10, cotton 1), none lost, *Populus* untouched; *Leptidea*'s 16 `low`
  `fam_only` calls gone; panel 29 `low` calls moved and 10 `fam_only` calls gone (*Sesamum radiatum*, windows of 0-1 K-mers),
  homoeolog agreement 665 of 1,119 pairs to 666 of 1,110. Per lineage, window right against replaced call right: *Solanum* 76
  vs 4, *Aegilops* 8 vs 0, *Chlamydomonas* 12 vs 0, *Capsicum* 10 vs 0. The 140 moved calls with truth are right 107 times
  where about 11 would be right by chance (the scorer's chance); the calls they replaced, 4 times against about 8. Still
  open: 77 calls flagged beside their window, within the readouts' resolution (15 genomes), window right 64 times and call
  13; and 2 checks with no gain, run on the development and added genomes: running the line in arms where it does not run now
  (+3, -4 on the `low` calls it would reach, from 0.17's F-arm windows) and narrowing a wide call (at least twice the window's
  width) to a window inside it (both right in 64 of 65, but the window's midpoint closer to truth only in *Populus*, 16 of 16;
  farther in *N. crassa*, 24 of 24).
- **The pair test (0.19, `fam_rep_pair`).** The family windows that hold a cluster are tested against the repeat-density peaks
  (`famline.with_repeats`: overlap, chance over the span both lines read, P ≤ 0.01, most chromosomes, each left out in turn);
  where that passes, the line acts as where its own test passed. Registered with criteria after the counterfactual on the
  stored 0.18.0 runs; replayed: 21 `low` calls wrong to right (*Ae. tauschii* 12, LA2134 8, moss 1), none right to wrong, no
  confident call changed, panel agreement 666 to 667 of 1,110 pairs, no change on the genomes added in 0.17. A code review
  found the first form's chance spread over the whole chromosome (a quarter too small, since neither line reads a terminal
  eighth) and empty windows counted as agreements; fixed, which undid its one move in the holocentric *Luzula sylvatica*.
  Still: both lines read K-mer copy number, so the pair cannot tell a centromere from any repeat cluster the two share; it
  passes genome-wide in *Luzula* (acting nowhere) and on siamang's 4 clustered windows of 25. A test of the family line
  against an independent input (Hi-C: #81) is the stronger route.
- **The move not made (0.16-0.17).** Of the `low` calls it contradicts (outside the terminal eighths), its window is right on 89 of 110
  in the release (15 genomes: 13 *Solanum*, *Ae. geniculata*, *Chlamydomonas*) against 8 for the calls, but on the panel the
  homoeologous pairs touched came no closer (64 pairs: 24 closer, 27 farther; agreeing 20 → 17). Such calls are flagged
  `fam_elsewhere`. In *Solanum* FASTA-only arms this is the largest headroom left.
- **The test is not calibrated.** A family 31-mer votes at most once per chromosome, so one family that clusters once per
  chromosome passes wherever it clusters (tomato: the 1-Mb peak on 7 of 12 chromosomes is one family of 365 linked 31-mers).
  With each genome's own copies and every repeat unit moved as a whole to a random place on its chromosome, it passes at its
  1% level in 13-52% of draws in each of 11 monocentric genomes (133 of 400); synthetic genomes whose planted families
  cluster at independent random places pass in 1-10% (3 settings of 100). Holocentric genomes did not pass until 0.17's
  *Leptidea sinapis*, which passes (its calls stay `low`).
- **A family-level null, tried in 0.16.1, not adopted.** 31-mers were grouped into families by linkage (copies side by side
  at the same distance in at least 3 places) and moved as units, against criteria set beforehand (holocentric genomes fail;
  planted controls pass at most 2% at 0.01; seeded; no choice by truth). Moving each 31-mer alone passed 199 of 200 planted
  nulls; moving whole families had no power (planted positives 0 of 50, tomato P 0.93); counting each family once per
  chromosome (four linkage versions) passed none of 20 planted nulls and none of 10 planted positives, but 64% of semi-real
  *N. crassa* draws, where one ~7-kb unit split into several families. The grouping is the problem: a diverged family splits
  into parts that seem to coincide, and a looser linkage merges neighbors and leaves no power.
- **Tested on the genomes added in 0.17, not adopted.** The move (`fam_over_low`: the `low` call moves to the window and stays
  `low`) was registered with its criterion before any new truth was examined: in at least 3 lineages not used to make it, the
  moved calls right at least twice as often as the calls they replace. Held out, the line contradicts `low` calls in two
  lineages only: *Gossypium* (window right 1, call 0, of 1) and *Populus euphratica* (window 2, calls 3, of 3). Pepper, whose
  windows were part of the rule's basis, gives 10 against 0 of 17. On the development benchmark, replayed with the move, 153
  `low` calls became right and 8 wrong (all *Solanum*), no confident call changed, and the panel's homoeolog agreement stayed at
  665 of 1,119 pairs. Why *Populus* differs: its satellite calls are 5-12 Mb wide and right on 19 of 19 chromosomes (chance 44%
  each), while the 200-kb family window is right on 18 of 19; of the three windows the move would take, two lie 270-350 kb
  outside their wide right calls and are right, the third 3.6 Mb past its call and wrong. The window is the more precise readout
  where both are right; the move waits for lineages where calls are narrow and the line passes.
- **Settle:** more truth where the line passes and the calls are `low` (#35: oat, wheat, other Solanaceae with ChIP). For the
  test: a null that deals each family's clusters among chromosomes by length, which tests the line's own claim (the old test
  passes 8 of 100 such planted genomes); or family membership from `tes/` annotations, which removes the grouping problem.
  **Code:** `famline.py`, `fuse.py` (`fam_only`, `fam_over_low`, `_fam_speaks`, `_fam_alone_against`, `fam_elsewhere`).

### 77. Leave-one-out for the one-locus tests

Open. Since 0.16.1 a test that validates one line by its agreement with another acts on a chromosome only at the weaker of the
genome's verdict and the verdict on the other chromosomes (design.md, Principles). Not so yet: the one-per-chromosome tests of the
FASTA lines (`at_only`, `rep_only`, `fam_only`, `at_chose_alt`, `fam_elsewhere`) and the satellite engine's karyotype and
signature tests, where a chromosome's own peak counts toward the pass that then places or flags it; leaving each chromosome out
reruns the 10,000-deal null per window size per chromosome (up to 200 x 9). Choices made on the whole genome also stay: the AT and
repeat window size, the methylation scan direction, the Hi-C engine's partner windows. First measure how many passes depend on one
chromosome. **Code:** `loo.py`, `atline.scan_profiles`, `sat/karyotype.py`.

### 80. A composition line alone stays `low` however often it is right

User's decision (0.19). A trough or a repeat-density peak that places a chromosome nothing else places is `low` by design (design.md,
AT troughs: AT-rich stretches of other kinds exist in every genome; `medium` only where Hi-C validated the troughs), though its
genome's one-per-chromosome test passed at a calibrated level (AT: 3.9% of 3,336 GC-shuffled genomes at 0.05; repeat density
1.3% of block shuffles at 0.01). In the 0.18.0 release runs `at_only` calls left `low` were right 79 of 83 times (*M. sympodialis*,
*N. crassa*, *P. falciparum*, *V. dahliae*; added in 0.17: *C. deuterogattii* 14 of 14, *P. carinii* 14 of 17) and `rep_only`
calls 67 of 67 (*C. neoformans*, *Chlorella*): better than the `medium` class (92% on the
development genomes, 83% on the added ones). An LTR core alone is `medium` on a weaker record (#79). Raising them where the
test passes at P ≤ 0.01 (the gate's level, at which repeat density already acts) would give *M. sympodialis*, *P. falciparum*,
*C. deuterogattii*, both *C. neoformans* and *Chlorella* confident calls from the FASTA alone (1 more confident error, in
*P. falciparum*); *V. dahliae* (P 0.0104) sits at that level. It rests on 8 genomes, 2 of them not seen while the lines were
made, and changes what `low` means for a whole class of genomes, so it waits for the user. **Code:** `fuse.py`
(`_confidence`), `atline.py`, `repline.py`.

### 81. The repeat-family line tested against Hi-C

Open (0.19, a lead). The family windows are tested against Hi-C's prior-free peaks and cis extremes in every run with Hi-C
(`fam_hic`, reported only: the line never pairs; since 0.19 only windows that hold a cluster, with chance over the span the
windows can occupy). In the 0.19.0 release runs the test passes on most chromosomes in barley (prior-free 6 of 7, P 3.7e-9;
the family test alone fails, P 0.20), *C. neoformans* (prior-free 9 of 9; cis 9 of 9 in G1), *N. crassa* (3 of 4),
*V. dahliae* (8 of 8), *S. pombe* (3 of 3), and in tomato, TS-60 and LA0458 (cis 9 of 12, 8 of 12, 7 of 12). Two lines from two inputs validating each other
is principle 2, the standard by which a trough pairs with a Hi-C call (`at_pair`, which needs no passing trough test); 0.16.1
removed `fam_pair` because the family test alone is not calibrated, a reason that does not touch this test. A pairing at
`pass` (the Hi-C call takes the window's agreement as a second line, one level up) would raise 19 right `low` calls in the
0.18.0 and 0.19.0 release runs alike (21 counting `cis_end_pair`, which stays `low` by its own rule; tomato and TS-60 with
and without `tes/`, *S. pombe*) and make 2 confident errors (a *C. neoformans* G1 cis
call 1 bin off, in two arms). **Settle:** a pairing rule registered with its criteria; *C. neoformans*' two cis near-misses
(design.md: cis and trans as `medium` would make them confident errors) kept `low` by a reason of their own, or the rule
dropped. **Code:** `famline.py`, `repline.with_hic`, `fuse._refine`.

## Methylation

### 16. The CHG line has changed no call

Open. Only *P. annua*, *Arabidopsis* and rice have CHG data here, and the line has never supported or moved a call on scored data.
*P. annua*'s five samples are other plants, so they do not decide, though the signal is real (the dip is the chromosome's lowest
point on 12 of 14). *Arabidopsis*' expected domain (315 kb, from genome size) is far below its CEN180 arrays (2.2-3.5 Mb), so the
reference ring (1.5-4 times the domain) falls inside the arrays. Rice T2T-NIP with same-assembly ONT CHG (0.16): the signature
forms (leave-one-out on 7 of 12), `centromere_like` on 4 of 12 calls (all right), but every rice call is `high` already. Downloaded,
not run: barley Morex leaf WGBS (bigWig) and *S. verrucosum* CPC54 ONT CG/CHG (another genotype than the benchmark's).
- **D, the domain size**, is a maize prior for every genome (22.2 Mb over Mo17's 10 centromeres, scaled by genome size, slope 0.97,
  Plackova et al. 2021); it sets the window, ring and candidate spacing, depends on the chromosomes listed, and is computed twice
  (`pipeline.py`, `methfuse.py`). A D measured from the anchor calls was never tested on a second same-plant case.
- **A per-site floor** (MIN_CALLS_SITE 3), which the CpG line dropped for per-bin pooling.
- **Fraction-only input**: the reader takes modkit or ENCODE bedMethyl only, so the public *Solanum* bedGraphs cannot be read.
- **Settle:** same-plant CHG with CENH3 truth (#35; *S. verrucosum* CPC54 first), with D from the genome's confident calls (median
  width, else the prior, reported) and CHG pooled per bin. **Code:** `meth/params.py` (D_REF_BP, D_REF_G, D_SLOPE, MIN_CALLS_SITE),
  `meth/track.py`, `meth/io.py`, `methfuse.py`.

### 55. CpG kinetochores outside primates

Open. The CpG line was designed on the primates it is benchmarked on (rhesus, *Arabidopsis* and HG002 were not used to design it),
and its specificity rests on rotated tracks and one plant.
- **No negative**: mouse should fail, but its ONT reads carry no modification tags (a HiFi 5mC track would do). No grass, other
  vertebrate or holocentric genome with truth run. Unused: the HGSVC3 CDRs (#35), giraffe's and rhesus's ONT methylation.
- **A gate without a majority**: with five chromosomes the gate has little power (2 of 4 at chance 0.05: P 0.014); wide calls
  weaken it too (#73). CG-only samples decide without a provenance test.
- **Small arrays** hypomethylated as a whole (bonobo and Bornean orangutan minicentromeres) are never refined (10 bins with data
  needed); a wrong array cannot be told (#38).
- **Genomes without truth**: the duck passes on 2 chromosomes of 50 (0.15.0's own test, with every rotation offset, passes none):
  21 of its 47 calls span at least a sixth of their (micro)chromosome, where no exact test within the chromosome reaches a small P;
  the width alone would pass it (43 chromosomes). The house wren likewise (0.15.0 with every offset: 0 of 40), which the depths pass.
- **Chromosomes are not independent under a shared move**: the toad's track moved by the same share of every chromosome (a
  quarter, a half) still passes the gate, though moves drawn per chromosome pass on 0.05% of simulated genomes; Fisher's method
  assumes the chromosomes' P values independent. Lead: a gate calibrated against moves shared by all chromosomes.
- **A perturbation**: CHM13 clones whose methylation was changed by dox moved CENP-A (PRJNA1270043, fiber-seq and DiMeLo-seq).
- **Settle:** kinetochores within 50 kb of the HGSVC3 CDRs; a mouse HiFi 5mC track fails the gate; in the perturbed clones the
  called kinetochore follows CENP-A; rotated controls silent. **Code:** `cdr.py`.

### 73. The CpG test in wide calls

Partly done. The positional test compares a call's dips with those that methylation from anywhere on the chromosome makes on the
same array bins, so its power falls as calls widen. Since the width joined the depths (0.16), gorilla (calls 3.5-8 Mb) and the
macaques pass with both technologies, gorilla's kinetochores on DiMeLo-seq CENP-A on 9 of 9 truth chromosomes. Open for genomes of
small chromosomes (duck, #55). Tried: Stouffer's method for the gate (rescues nothing, lowers siamang); testing only the largest
array run (arrays fill 90-99% of calls). **Leads.** The active array alone, where the satellite line can name it (HOR structure); a
test of the dip's place within the array rather than its depth against the chromosome; for microchromosomes, rotations across
chromosomes (which changes the question: chromosomes differ in methylation as a whole). **Code:** `cdr.py`.

## Inputs, chromosome set and other features

### 8. Large unplaced scaffolds and dot chromosomes

Partly done. On 121 polyploid assemblies with known chromosome lists, the length rule keeps 31 sequences outside the lists (1.5-117
Mb) in 12 assemblies, and extra sequences dilute the prevalence tests (on the lists *Cissus* gains 14 `high` calls, *Mimulus* 16;
*Malus toringoides* moves 19 calls). It also drops a dot chromosome alone after a tenfold gap (*D. melanogaster* chr4 without the Y),
dots before a tapering scaffold tail (a simulated 37-chromosome bird karyotype keeps 31), and small chromosomes spanning a tenfold
range; it keeps organelles (yeast's mitochondrion, NC_001224.1, enters every *S. cerevisiae* denominator). MAX_CHROMOSOMES 50 (a
warning) and MAX_SEQUENCES 5,000 (no one-per-chromosome test) express one concept 100-fold apart. Complete molecules left out
(telomere arrays at both ends) are named (`left_out_telomeric`), never added; it fired in none of the 0.15 release runs and cannot
act without a telomere unit (#56). **Settle:** large scaffolds kept without a telomere, centromere-like array or genes flagged; an
organelle test from the sequence (the organelle rRNA context `rdna.py` reads, or size and copy number); the 121 lists, *D.
melanogaster*, bird karyotypes and the simulated karyotypes as regression tests, the rule's 383 assemblies unchanged. **Code:**
`inputs.py` (`chromosome_set`, `left_out`), `meth/track.py` (`chromosome_count`).

### 56. Telomeres and rDNA arrays the benchmark does not count

Partly done: benchmark.md ([Other features](benchmark.md#other-features)) lists every genome's telomere unit, ends with an array,
complete chromosomes and 45S arrays. Open: a unit is learned only from arrays at two ends or more (`telomere.MIN_ENDS`), so potato
(one 378-bp TTTAGGG array at chr11's start), *S. pombe* (at most 8 copies of a GGTTAC-like k-mer per 10 kb) and *C. nigoni* (one
240-bp TTAGGC array, under `MIN_ARRAY`) get none, and a single-end array is not reported. NOR recall was never measured (goat:
arrays on 3 of its 5 NOR chromosomes); the 45S scan needs 3 motif copies, so a collapsed 2-unit locus is missed (*Cryptococcus*
JEC21). Four T2T mammal ends hold an array the truth rule (terminal 10 kb) does not count, never inspected. **Settle:** single-end
arrays reported; a 2-copy rDNA rule from the organelle context the 45S scan reads; telomere and NOR recall in the benchmark tables.
**Code:** `telomere.py`, `rdna.py` (`MIN_UNITS`).

## Several plants

### 22. Pa1A in other P. annua plants

Waiting for data. Young CRM insertions of four of the five other plants cluster 9-15 Mb from Pa1A's centromere (163.0-169.4 Mb),
which suggests other centromeres; J10's lie inside it, and C11 and J10 show a plant-only methylation dip at their clusters. 0.15
calls Pa1A at 163.0-169.4 Mb (`high`, #31), with plant verdicts 4 `shared` (C11, J10, P3, R1) and 1 `unclear` (I3). Over all 70
plant-chromosomes: 59 `shared`, 5 `unclear`, 2 `possible_shift` (J10 Pa4B, Pa5A), 2 `call_absent`, 2 `no_data`, none `shifted`, at
4-8× reads. **Settle:** CENH3 ChIP or deeper long reads from those plants.

### 57. Per-plant verdicts that depend on other chromosomes

Open. A plant's verdicts change when another chromosome's call moves: when 0.13 moved Pa1A, J10 Pa4B went from `shifted` (loss q
0.015) to `possible_shift` (q 0.32) at the same site, and R1 Pa3B from `shared` to `unclear`. The cause (the multiple-testing pool
or the between-plant test) was not traced. Untested since 0.1 (plants.md, Untested): the noise model on two samples of one plant,
coverage below about 2.5×, the per-plant Hi-C path (pairs cannot be genotyped, so another plant's library decides as `unknown`).
**Settle:** a leave-one-call-out check of each verdict; one sample's reads split in two must come out as the reference plant.
**Code:** `plants.py`, `provenance.py`.

## Reporting and outputs

### 58. Report wording and undocumented outputs

Partly done (0.14: failed sections recorded in `incomplete`; every column documented, with a test; 0.16: `hic_sequence_conflict`
reworded as a note, `hic_rabl` documented apart from the Rabl-axis tests, a genome-level verdict, #44). Open:
- **One table of flag → basis phrase, caveat phrase and settling input**, read by `basis` and `caveat`, tested flag by flag (#63).
- **`tes/` advice where it cannot act**: caveats advise `tes/` in mammals and fungi too (`interpret._could`; the siamang's confident
  misses, every *C. tropicalis* no-call); cenfuse does not know the lineage.
- **The HOR summary**: the strongest class period equals the top tandem unit on 6 of 46 HG002 calls (`sat/hor.py`).
- **No per-run homolog check** (homoeologs, haplotypes), though cheap (one goat's assemblies chose different main families); the
  panel check (#68) runs outside cenfuse.
- **Replaced calls**: their own confidence is not stored; `basis` calls a replaced long-unit array a "weak satellite array".
- **A Hi-C call inside a decisive LTR core** (`cis_only`, `cis_trans_pair`) gets no `te` support (`fuse._support` counts the
  core only through a strong peak's second line), so its caveat says "no sequence or LTR evidence at this place"; 0.16.1 only
  stopped listing the core against it. No release call is affected; `support` also feeds the long-unit rule, so change both together.
- **Code:** `interpret.py`, `output.py`, `sat/hor.py`.

## Speed, memory and reproducibility

### 5. A faster search

Partly done: since 0.11 long models search with low-complexity tracts neutral, 1.00-4.93× faster on 14 genomes. Open:
- **What it does to results.** Over 14 full reruns no chromosome with truth (143) changed, but only 8 of those genomes have truth and
  the rest of the benchmark replays older searches (#34). Overlapping families trade the lead in *Cuscuta* and *Panax* (19 calls
  change main family, 2 change place): a neutral model scores less, so a close family without a tract wins the overlaps (cotton: 906
  copies move to a relative). **Settle:** every arm rerun from the search; overlaps compared only on positions both models can score.
- **The cost of real hits.** Chain families hit nearly everywhere; dense short-unit arrays hold the rest (alpha satellite: 62% of
  CHM13's CPU; mouse: 14 dispersed models took 95% of a 2-Mb chunk). Tried: [Tried and set aside](#tried-and-set-aside-not-in-designmd).
  **Settle:** a search that treats the seed arrays apart, with monomer and array tables equal to the stored ones (CEN180 on Col-CEN:
  66,131 in Naish 2021, 67,241 by the satellite engine).
- **The runtime model** of the progress estimate (model length × genome size) fails on heavy genomes (*Panax*: 7.9 h predicted, over
  25 h actual; R² 0.72 over 104 genomes). Runs record their cost (`run_info.json`).
- **Code:** `sat/annotate.py` (`search_cost.tsv`, `resolve_overlaps`), `sat/util.py` (progress).

### 6. Faster seed clustering

Open. BLAST takes 98% of seed clustering, driven by multimeric seeds (*Acanthus* 76 min at 16 threads; median over 106 genomes 0.7
min), and the clustering reads little of its output (*Acanthus*' heaviest round: 7.0% of the HSPs). The thread count changes no
family (*S. verrucosum* and *Acanthus* at 4 and 16 threads), though BLAST's hit lists depend on it below what the clustering reads
(`sat/families.py`'s comment says the hits are identical). Tried: aligning only the consulted pairs (10× faster, but hit lists
differ for 0.2-0.5% of them, so not exact); a fixed partition of queries over single-thread BLAST processes (0.79× on the heaviest
round, not identical). **Settle:** a two-pass search (the best HSP per pair, then full tiling only for consulted pairs) with
identical clusters and representatives on every stored seed set. **Code:** `sat/families.py` (`cluster_seeds`, `blast_pairs`).

### 7. The FasTAN crash, reported upstream

User's decision; not filed yet. A 642.5-kb window of *Oryza minuta* (GWHESQK00000003:29,940,000-30,582,500), about 90% one satellite
array (15-mer spacings of 243 and 342 bp), aborts FasTAN 0.8 every time (`malloc(): corrupted top size`); neither half crashes
alone. *Sonchus oleraceus* crashed at another site. cenfuse recovers by splitting the piece (757 rows after 4 splits, nothing
skipped). **To do:** file with the reproducer (`FAtoGDB`, then `FasTAN -m -p -T1` on the window) and cenfuse's crash table; keep
the splitting.

### 60. Reproducibility across environments

Partly done (design.md, Engineering notes): `run_info.json` records the versions of cenfuse, Python, packages and every tool, the
inputs' fingerprints and the cost; `environment.lock.txt` pins the exact environment; the seed cache keys on FasTAN's binary and
layout state; `--satellite-dir` warns on another version; workers log under every start method. Open: on macOS, without `setarch`,
about 1% of FasTAN's array cuts differ (a warning only): measure what that does to calls and state it. CI
(`.github/workflows/tests.yml`): record its first result. **Code:** `manifest.py`, `sat/seeds.py`, `sat/fastan.py`.

### 61. Memory

Open. The satellite search holds about one workspace per thread (about 2 GB for 2-kb units) with no cap, so `-t 32` on a 64-GB
node can be killed. Peaks reach 74 GB (rhesus, 12 threads), well above 12 workspaces, unexplained; the 13 tomatoes (0.8-1.2 Gb, 16
threads) peaked at 7-46 GB, so genome size and threads do not set it. The CHG line holds about threads × 1-2 GB; the 3D
feature-class test builds dense N × N matrices and is skipped above 15,000 bins. Tried: [Tried and set
aside](#tried-and-set-aside-not-in-designmd). **Settle:** profile the rhesus peak; then a `--max-memory` option limiting concurrent
long-model searches, with identical outputs; the feature-class statistic from the sparse matrix. **Code:** `sat/annotate.py`,
`hicfuse.py` (`feature_contacts`).

### 62. Hi-C engine speed

Partly done in 0.14 (the 139 Hi-C stages 4.6 h instead of 8.5 h, every number unchanged). Open: each Hi-C library reruns the
whole stage; a 2n = 80 polyploid was never timed on real data. **Settle:** identical tables, faster on chicken and a 2n = 80
polyploid. **Code:** `hic/confidence.py`, `hicfuse.py`.

## Code and docs

### 63. The integration cascade and other large functions

Partly done in 0.14: `fuse._decide` is one function per rule in a declared order (`CASCADE`, `PHASES`), each call names its rule
(`rule`), and an equivalence test checks it against 0.14.0's single function. Open:
- `interpret.genome_section`, `caveat` and `basis` need coordinated edits for every new flag (#58); `pipeline.Run._run` is long.
- `hicfuse.py` holds four concerns and drives the Hi-C engine through a hand-built argparse Namespace (a library call instead).
- Binomial tails written several ways and the Hi-C confidence tests twice (one `stats.py`); the chromosome-set rule lives in the
  methylation engine; the LTR signature is computed twice per run and `tes/` parsed about six times.
- The four CLIs differ in exit codes, logging and the meaning of `-o`; the engines' CLIs are undocumented.
- Bare numeric literals inside functions (the constants test covers module-level ones only: closed #64).
- **Settle:** byte-identical benchmark outputs before and after; `test_fuse_rules.py` as the per-rule spec. **Code:** `fuse.py`,
  `interpret.py`, `pipeline.py`, `hicfuse.py`.

### 65. Tests that pass without running

Partly done in 0.14: `CENFUSE_REQUIRE_TOOLS=1` turns a skip for a missing tool into a failure; unit tests cover the checkpoint
layer, locks, manifest, constants, `cli.py`, `linetests.py`, `seqcalls.py`, `teline.py` and parts of `loci.py`, `methfuse.py` and
`provenance.py`; 1,232 tests pass (0.16.0). Open: the satellite tests use test-scale values for 13 of the 51 parameters (one slow
test at production values); the parallel family build is exact only while `frame_rotation` writes nothing a later step reads, a
premise guarded by a docstring (make it a test); the report's sentences are tested only end to end (#58); no fixture holds IUPAC
codes; CI's first result (#60). **Code:** `tests/conftest.py`, `tests/sat/sat_common.py`, `tests/fuse/test_fuse_e2e.py`.

### 66. Statements the outputs no longer bear out

Partly done (each release's docs are audited against its runs). Open: the audit as a release step (docs checked against each
release's outputs; docstrings citing benchmark sections, not numbers); figures not redrawn for the release say which release drew
them; potato's dips on CENH3 count 8 of 11 in the cis survey and 9 of 11 in the cis-trans text; `sat/families.py` on BLAST hits at
any thread count (#6).

### 67. Inference from ethology of related species. 

Cenfuse may take related genomes and infer orthologous centromeres etc. 
If centromere call is weak or missing in one ortholog but strong in the other, maybe we could use the strong to place the weak or missing. 
We could also potentially report other weird biology, like centromere positional/composition shifts OR centromere replacement (activation of a new and deactivation of an old). 
Thus, with related species, we could potentially report other helpful comparative centromere features and/or use them to piggyback placement and confidence off each other. 
All of this would need tested. 

## Tested, not adopted
## Note
Reopen only with new evidence.

- **3. Satellites specific to subgenomes vs the honesty cap.** Dominance judged within subgenomes (design.md, Choices tested and
  rejected). Reopen with CENH3 for an allopolyploid with subgenome-specific satellites (*Fragaria × ananassa*, Jin 2025: #35).
- **40. A centromeric satellite smaller than a shared family's cluster.** Potato's centromeres sit on kb-unit families (St24, 981
  bp, 299 kb at Cen1) while a 157-bp family dominates. Taking the largest cluster of families on at most 3 chromosomes places 2 of
  11 potato, 0 of 12 tomato, 3 of 12 *S. verrucosum*, 16 of 33 chicken, 4 of 17 *Chlamydomonas* (0.15): these satellites are rarely
  the chromosome's largest array. Reopen with a readout that does not need the largest array (#39).
- **67. Fragments and solo LTRs as a core locator.** The signature families' LTRs mapped back with minimap2 (one chromosome left
  out) peak on CENH3 on 248 of 299 chromosomes with an element core, against 237 for the element core and 248 for the age choice
  (251 both; 0.15): +3 for a new external stage. Reopen with a larger gain.
- **75. Split ties in LTR core profiles.** Of 52,949 core profiles of the 0.15 `tes/` runs, 82% have tied maxima of the smoothed
  counts and 50% tied maxima in separate half-maximum runs, so the leftmost maximum decides half of all cores. A tie-break by the
  run holding more elements moved 24 of 634 final cores (2 wrong → right, 4 right → wrong) and 3 seeds (0.16). A smoother density
  estimator would change cores broadly; not tried. **Code:** `sat/te.py`, `sat/signature.py`.

## Leads for biology

Hypotheses cenfuse's outputs could test; none is a finding yet.

- **Young insertions mark where a centromere is now.** CRM copies are younger at Pa1A's Hi-C site than at its satellite (#31), and
  choosing the youngest comparable LTR peak moved 27 cores over 32 genomes with truth, 12 wrong to right and 1 right to wrong (#71).
  Test: centromeres that moved within a species with CENH3 tables (soybean, tomatoes: #35) should carry young copies at the new site.
- **Which lineage CENH3 binds.** Tobacco's CENH3 binds Tekay, not the CRM beside it; a Tekay signature marks all 13 tomatoes (18- to
  53-fold). Tekay may be the CENH3-bound lineage across Solanaceae (#46).
- **The Rabl configuration and condensin II.** The Rabl readouts pass in fungi, *Plasmodium* and barley, in no territorial nucleus;
  condensin-II-deficient lineages keep Rabl (Hoencamp et al. 2021). Predictions: Diptera (no CAP-G2) pass; human cells without
  CAP-H2 gain the alignment (public Hi-C); declustering mutants (*S. pombe* csi1) lose it (#14).
- **What clusters in a mammal nucleus.** CHM13's off-target trans peaks sit on pericentromeric HSat2/3 and β-satellite, not the
  active array; clustering differs between mammals at equal depth (pig, horse, sheep, giraffe fail; goat, rhesus reach support).
  `hic_architecture` labels about 50 genomes, without a compartment-matched null (#52).
- **The cis signal's sign**: a peak in tomato, rice and *C. neoformans*, a dip in potato and *Verticillium*, absent in mammals (#29).
- **Telomere repeats at bird centromeres**: chicken's interstitial telomeric arrays lie at its CENP-A windows (P 4e-21; Huang et al.
  2023); other bird T2T assemblies would show whether this is general (#41).
- **5S rDNA beside HSA1q42 genes** in every primate, horse ECA1 and pig SSC14: a 5S locus perhaps conserved across mammals.
- **Methylation and centromere position**: whether the CpG kinetochore follows CENP-A in CHM13 clones with changed methylation
  (PRJNA1270043) tests the CpG line's premise (#55).

## Tried and set aside (not in design.md)

Measured ideas that design.md's "Choices tested and rejected" omits. Reopen only with new evidence. *Where*: the development
folder (see Where things live; cenrep, cenhic: the engines' own projects).

| idea | result | where |
|---|---|---|
| `malloc_trim` after each search; model-major job order; `MALLOC_ARENA_MAX=2`; freeing hit rows | 229 GB; 125 GB; no change; no change (the fix was one pyhmmer pipeline per thread) | v4-v5 (#61) |
| smaller pyhmmer blocks | 5× less memory and 3× faster on a CEN180 job, but about 40% of raw hits keep their position (multimer models align at any phase) | v5 (#5) |
| tandem-copy models for short units | 4 copies: 43 s per 100 kb (1 copy: 0.1 s) | cenrep (#5) |
| a sample screen to skip dispersed families | 128 families skipped on 216 stored genomes with no false skip, none on real runs: under 7 chromosomes cannot reach P ≤ 0.01; array-poor families cannot be measured (moss's fam009, 41% of its CPU); dense chaining families give too few breaks; removing a family shifts the others' arrays | v8 (#5) |
| a BLASTN prefilter for long models; a minimap2 prefilter | lost 12% of the hit bp; never tested | v9 (#5) |
| BLAST `-mt_mode`; single-thread query jobs | identical in 12 of 13 reruns, mixed times; 1.6-1.9× per call, 3 of 873,000 hits differ | v9 (#6) |
| a CRM minority as the genome's marker | two right `high` *Arabidopsis* calls fell to `low`; 41 genomes changed | v6 (#46) |
| a share-in-core cutoff for specificity | a cut in a 0.18-0.24 gap fitted to five clades | v6 |
| sites in the core's binomial test too | cores moved in 113 of 220 replay genomes | v8 |
| element counts without a margin | failed its design case (*B. rapa* A09: call 16, core 18) | v8 |
| a core that may favor any cluster it touches | lost *Ae. geniculata* chr02U's CENH3 call | v3 (#9) |
| every candidate clade must disagree before a satellite is rejected | Triticeae chromosomes give 6-8 clades; nothing is ever rejected | cenrep |
| a two-fold concentration margin between clades | labels *Aegilops* genomes inconsistently | cenrep |
| agreement as "the midpoint of either interval inside the other" | 13 of 14 on *Ae. geniculata*, but overrides CENH3-bound satellite calls (CentO-C 4/8 in *O. minuta* C) | cenrep |
| a minimum fold for recruited families | right cores 155, 155, 148, 131 of 172 at none, 5, 10, 20-fold; barley 7 → 3 at 10-fold | v11 |
| an LTR core choosing the competing cluster | no benchmark case to build on (moss's 5 of 5 were solved another way) | v7, v11 |
| the dispersed test without the locus fold | moved 160 batch calls off retrotransposon-derived arrays (*Coffea*, *Glycine*, *Rosa*) | v3 |
| gate chance from the engine's null replicates | removed 1 of 2 spurious passes on randomized sheep; the whole-set majority removed both | v3 |
| full credit for blind bins; replacing only empty partner windows | supported a wrong bonobo call; right anchored peaks 97 → 86 of 139 | v3 |
| both flanks of an unmappable array rising | 75 of 82 centromeric gaps, but also 67 of 96 other gaps under 5 Mb | v3 |
| the fold center by the lower confidence bound | biased to the middle (barley 2/7) | v3 |
| numpy huge pages | the Hi-C bootstrap 30× slower on fragmented memory (now off at import) | v3 |
| peaks restricted to B-compartment bins | gorilla 0 → 0.60, *P. annua* 0.93 → 0.36: the E1 sign at centromeres is opposite in plants and apes | cenhic (#14) |
| arm-scale segregation across the centromere | CHM13 6 of 23, plants at chance | v10 |
| the compartment eigenvector or the trans fraction as a locator | at or near chance; the trans fraction works on a minority, in mammals only | v10 |
| an effect-size test of nucleus type | separated holocentric and randomized libraries from yeast, not Type I nuclei from Type II | v10 (#52) |
| a confidence cap on `hic_sequence_conflict` | the flag fires where the sequence calls are right | v8 (#58) |
| the arms' depletion as a locator | tomato 0 of 12, potato 0/11, rice 0/12 per chromosome (a mean effect only) | v11 |
| the engine's peak levels required for core-anchored Hi-C | cleaned the controls, but left no real hits outside genomes that pass anyway | v11 |
| the Hi-C table swapped for the core-anchored one | potato 6 of 12, short of a majority | v11 |
| CUSUM on the GC deficit; the window by the largest trough gap | *M. sympodialis* 0 of 8; wrong windows in *C. neoformans* and *S. pombe* | v11 (#28) |
| lifting the honesty cap where one-per-chromosome passes | chicken 27 right, 5 wrong; pig 8 right, 9 wrong | v11 (#41) |
| interstitial telomeric arrays as a line | only chicken beats chance (9 calls carry one, 3.9 expected, P 6e-4) | v11 (#41) |
| a genome-level holocentric flag | also on 25 of 26 siamang and 25 of 32 EquCab3 calls | v4 (#44) |
| CpG: a cluster whose valley floor is at most half the call's | 6 changes in 7 primates: 1 right (CM054652.2), 2 wrong | v4 (#38) |
| a CpG null from other satellite DNA; CpG bins as the mean of per-CpG fractions | infeasible (MFA8: 3.7 Mb outside the calls); MMU8 0 of 21 kinetochores at 1/16 depth | v4, v6 |
| a data-derived CHG domain size | rejected, numbers not kept | v5 (#16) |
| 45S arrays by the Rfam SSU HMM; even spacing for 45S arrays | stalled on low-complexity sequence; only 100 of 272 arrays are regular | v3, v10.3 |
| 5S from the satellite families' Rfam tags | misses units over 2 kb (mammals, birds) | v6 |
| sheep and goat CENP-A ChIP as truth | the antibody (CST #2187) recognizes phospho-Ser7 of CENP-A, absent in those species: at most 1.1-1.7× | v2 (#35) |
| loci split into dense runs; sparse clusters dropped | chimpanzee +5,757 loci; barley lost 108 of 156 subtelomeric loci | v8 |
| the family window and the repeat-density peak together over any `low` call they agree against | interior +44/-10, in a terminal eighth +29/-10 (0.18.0 runs): siamang loses 6 right calls, mouse and *Z. tritici* their right end calls; adopted only as a test (`fam_rep_pair`: on most chromosomes, windows that hold a cluster, interior calls), under which siamang's 6 are spared because their windows hold under 3 family K-mers | v20 (#76) |
| `hic_disagrees_weak` (a validated peak elsewhere, too weak to move) as a contest | such calls right 37 of 39 `high`, 7 of 10 `medium` | v20 |
| AT and repeat-density calls only where the chromosome's window ranks among the genome's N deepest (`first`) | calls off the deepest class right 30 of 32 | v20 (#77) |
| a right call narrowed to its overlap with the lines that support it | the domain stays in the narrowed call: LTR core 87%, Hi-C 75% (plants) and 44% (animals, beside the arrays) | v20 (#37) |
| the array closest to its family model as the active one | *N. benthamiana*: 8 of 13 within 1 Mb of CENH3, as for the largest | v20 (#78) |

## Closed entries

One line each; the full text is in git history.

- **11. Hi-C anchors where the sequence line is wrong genome-wide.** 0.12: the LTR signature makes better anchors, and where the
  test still fails Hi-C is tested against the LTR cores.
- **12. Cis signals of satellite-free centromeres.** 0.12: the cis profile, its sign and level tested against the LTR cores; without
  cores, #26.
- **17. Knobs and the size floor.** 0.11: `classify.related`'s period test works both ways, so maize's knob family no longer counts
  CentC's chromosomes (size floors: #39).
- **18. Organelle insertions called 45S arrays.** 0.11: a 45S motif hit counts only with nuclear context.
- **19. Long-unit calls are mostly wrong.** 0.11: an unsupported long-unit array is a feature, not the call (`long_unit_unsupported`).
- **20. The rDNA spacer tag at its threshold.** 0.14: the 11 *B. rapa* accessions run as one complex tag the family in all 11, 110
  of 110 calls right.
- **21. Moved calls whose sentences disagree.** 0.11: `caveat` reads the evidence `basis` reports.
- **23. A family's unit length can overstate its copies.** 0.11: `basis` gives the array's own spacing (`seq_spacing`).
- **24. Telomeres outside the rules.** 0.11: 13-29-bp units and loose stretches at sequence ends; genomes with no unit: #56.
- **25. Non-CRM centromeric LTR retrotransposons.** 0.12: the family-level signature, in which another lineage decides at 10-fold;
  CRM's own standing: #46.
- **26. The cis signal without `tes/`.** 0.13: the library's own trans contacts test the cis extremes (`cis_trans_line`); what
  remains: #32, #53.
- **33. Cis and trans at one locus as a confident call.** 0.16: the AT troughs are tested against the cis extremes too, and a cis
  call whose tested trough agrees is `medium` (`at_pair`; G1 *C. neoformans*: 5 calls, all on CENP-A); a tested repeat-density
  peak does the same (`rep_pair`; the repeat-family peak's `fam_pair` was removed in 0.16.1, #76).
- <a name="36-no-held-out-test-of-the-013-rules"></a>**36. No held-out test of the 0.13 rules.** 0.15, the user's decision
  (2026-11-08): no data are held out; rules made on few genomes are checked under Rules for any change, 3-5.
- **59. Checkpoints that miss changed inputs or code.** 0.14: one checkpoint rule (code, parameters, input content), `--force all`,
  an output-folder lock, atomic writes (design.md, Engineering notes); a changed search stage stops for `--force` by design.
- **64. Constants without a reason in design.md.** 0.14: every module-level constant has a reason, enforced by
  tests/fuse/test_constants.py, and docs/constants.md is generated from the code; bare literals inside functions: #63.
- **69. Why a null test could not be run.** 0.16: a null that ran but could not score a chromosome is flagged `null_degenerate` (P
  written NaN); only `null_untestable` (too few partner chromosomes) counts as a test the peak could not take.
