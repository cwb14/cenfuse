# Changelog

What changed for users in each release, newest first, with the benchmark headline of the release. Current results:
[docs/benchmark.md](docs/benchmark.md). A *chromosome-run* is one chromosome in one benchmark run; *confident* = `high`
or `medium`.

## Docs update (2026-11-23, no new version)

- **LTRquest is now [LTR_peel](https://github.com/cwb14/LTR_peel).** The `tes/` recipe is LTR_peel at its defaults
  (`ltr_peel --genome asm.fa --threads 32`); the GFF3 has the same fields, so files made with LTRquest work as before.
  `python -m cenfuse.meth run` takes `--ltr_peel` (`--ltrquest` still works). The benchmark annotations were made under the old
  name with one round of detection; cenfuse was not re-benchmarked on annotations made with today's defaults.
- **The performance figure and the README's numbers now cover every benchmark genome**: the 22 genomes added in 0.17
  (0.19.0's results, `benchmark/results_017.tsv`) are pooled with the others under their centromere type
  (`benchmark/canonical_017.tsv`); before, both drew only `benchmark/results.tsv`. Over all 80 genomes with a mapped
  centromere `high` calls were right 557 of 569 times, `medium` 160 of 186, `low` 217 of 280. The Hi-C control is now
  panel C, beside B. docs/benchmark.md still tabulates the two sets apart.
- No call changed: results and tables are those of 0.19.0. Three gallery pages (siamang, both *C. neoformans*
  runs) are now the ones 0.19.0 draws, with the repeat-family line's new verdict; the other 44 were already identical.
  The gallery's note for *Ae. tauschii* gives 0.19.0's FASTA-only result (7 of 7, `low`).

## 0.19.0 (2026-11-18)

- **The repeat-family line is tested against the repeat-density line.** The family windows and the densest-repeat peaks of
  the FASTA are compared chromosome by chromosome, as the AT troughs are compared with Hi-C: where they overlap on most
  chromosomes, more often than each chromosome's own geometry gives by chance (P ≤ 0.01; a chromosome's own agreement never
  counts toward the test that lets the line act on it), the family line acts on a chromosome whose window and peak overlap
  as it acts where its own one-per-chromosome test passed (flag `fam_rep_pair`, with `repeats` added to the call's support;
  summary columns `fam_rep` and `fam_rep_acts`). The family line's own test can fail where its windows are right: in
  *Ae. tauschii* from the FASTA it gives P 0.015 (level 0.01) with the windows on CENH3 on 7 of 7 chromosomes, while the pair
  test passes there at P 1.9e-12. Only windows that hold a cluster take part, and chance is taken over the span both lines
  read (neither reads a terminal eighth). Benchmark: 21 `low` calls went from wrong to right
  (*Ae. tauschii* 12, *S. juglandifolium* LA2134 8, moss 1), none from right to wrong, and no confident call changed; no call
  changed on the genomes added in 0.17 (the holocentric *Luzula sylvatica* passes the pair test genome-wide but on no
  chromosome with itself left out); the polyploid panel's homoeologous calls agree once more (666 to 667 of 1,110 pairs). Both
  lines read repeat copy number, so a pass says where a centromere made of repeats would sit, not that one sits there: the
  line still makes only `low` calls and never pairs with Hi-C; where a window lies off a `low` call but within both
  readouts' resolution of it, or another line speaks to the call, the call only gains a flag (`fam_elsewhere`; LA2134, moss,
  siamang, and on the panel *O. schlechteri*) (future.md #76).
- **Fixed: the repeat lines' tests against Hi-C took chance over the whole chromosome**, though their windows skip the
  terminal eighths, so for a Hi-C peak in the interior chance was a quarter too small and the test too easy to pass (in
  simulation, 1.3-1.6% false passes at the 0.01 level; 0.2-0.4% now). The repeat-density test is what lets a Hi-C call with
  an agreeing peak rise to `medium` (`rep_pair`). Chance is now taken over the span the windows can occupy, as in the family
  line's test against repeat density; the family line's test against Hi-C (reported only) reads only windows that hold a
  cluster; `line_tests.tsv` does both for pairs with a repeat line. No call changed: `rep_pair` acts only in
  *C. neoformans*, where both tests still pass (P ≤ 1e-13). *C. tropicalis*' repeat density no longer agrees with Hi-C beyond
  chance (4 of 6, P 0.022; it moved no call before), two weak family agreements with Hi-C no longer pass (siamang, LA1982),
  and where no family window holds a cluster (*Bombyx*, *Z. tritici*, *P. falciparum*, mouse, horse and others) the
  family test against Hi-C is now untested.
- **The report's verdict for a genome without a confident call names the tests that passed but cannot make a call
  confident**: the repeat families' own test, their agreement with the densest repeats, and either repeat line's agreement
  with Hi-C, each with how many chromosomes agree, and it counts the repeat clusters among the possibilities. Before, it said
  no test had passed where one of these had (the tomatoes, *S. verrucosum*, *S. pombe*).
- **Tested and not adopted:** an LTR core no longer vouching for an array inside it whose elements lie significantly more in
  the rest of the core (a binomial test; the edge-of-core errors, future.md #45). On the development genomes it lowered 8
  confident errors and 4 right confident calls to `low` and made none wrong (a fifth right call went from `high` to `medium`
  through a changed Hi-C reference); on the genomes added in 0.17 it lowered only one call, a right one (pepper G1-36576
  Chr07); on the polyploid panel it lowered 16 confident calls. It missed two of its four registered criteria: the calls it
  lowered were not less reliable than confident calls on the added genomes, and three calls moved as side effects (*Ae.
  tauschii* Chr5 with Hi-C, lowered and then moved by the cis Hi-C line, wrong to wrong; two pepper calls with Hi-C, right
  to right, when Hi-C's references changed). The numbers and the code's place are in future.md #45.
- Docs: future.md's "Where to start" opens with where the integration stands (it uses nearly everything its tested lines
  offer; what remains needs a better line or more truth) and adds two entries, #80 (a passing AT trough or repeat-density
  peak alone stays `low` though its calls were right 146 of 150 times; a decision left open) and #81 (the repeat-family line
  tested against Hi-C).

## 0.18.0 (2026-11-18)

- **A repeat-family window now takes the place of an unsure satellite call it contradicts.** Where the repeat-family line
  passed its test, a `low` satellite call that no other rule changed, outside the terminal eighths, moves to its
  chromosome's family window and stays `low` (`fam_over_low`; `seq_start`/`seq_end` keep the call it replaced) when four
  things hold: the satellite could not vouch for the call in this genome (its families failed the one-per-chromosome test or
  none dominates, or the call is the weak-satellite fallback); the window holds a cluster, at least 3 family K-mers; the
  window lies farther from the call than the wider of the two intervals, the rule cenfuse already uses to judge whether a
  satellite array and an LTR core agree; and no other line speaks to the call on that chromosome (no Hi-C reading acting
  there, no agreeing LTR core, tested AT trough or repeat-density peak, no centromere-like methylation, no CpG line). Each
  condition is a principle, not a fitted value: a tested line outranks the family line, whose test is not calibrated.
  Benchmark: on the development genomes 94 calls went from wrong to right and 2 from right to wrong (one tomato
  chromosome, without `tes/`), and no confident call or control changed; on the genomes added in 0.17, 11 went from wrong
  to right (pepper 10, cotton 1) and none the other way; the polyploid panel's homoeologous calls agree a little more
  often (665 of 1,119 pairs to 666 of 1,110). The 140 moved calls with truth are right 107 times, where about 11 would be
  right by chance; the calls they replaced were right 4 times, against about 8 by chance. From the FASTA alone,
  retrotransposon and satellite-free centromeres go from 60 to 114 of 272 chromosomes right (13 tomatoes 20 to 51 of 156).
  0.17 had tested the move on any gap and not adopted it: in *Populus euphratica* the windows lay 0.3-3.6 Mb outside wide
  satellite calls that were right. Requiring the window to point beyond both readouts' resolution answers that miss, so
  the rule is listed among the rules changed after a result was seen (future.md #76).
- **The repeat-family line speaks only where its window holds a cluster** (at least 3 family K-mers, the line's own
  smallest cluster beyond a pair), for the calls it places too. In the holocentric control *Leptidea sinapis*, whose
  family test passes on windows of 0-2 family K-mers, it no longer places its 16 `low` calls; no benchmark call with truth
  changes.
- **Faster repeat lines.** The repeat-density and repeat-family lines read the FASTA on the run's threads (`-t`) instead of
  one, with byte-identical outputs at any thread count: on a 2.0-Gb genome, 0.17's code took 206 s and 233 s, 0.18 at 16
  threads 34-42 s and 42-51 s, at no higher peak memory (at one thread, about 0.4 GB more) (future.md #74).
- Report: the genome section counts the calls a repeat-family cluster placed apart from those it took over.

## 0.17.0 (2026-11-18)

- **22 genomes added to the benchmark** to test cenfuse where it was weakest (`benchmark/truth_017.tsv`, every scored
  chromosome in `benchmark/results_017.tsv`, [docs/benchmark.md](docs/benchmark.md#genomes-added-in-017)): two peppers,
  *Populus euphratica*, *Oryza australiensis*, two *Sesbania* lines, four allotetraploids of the polyploid panel
  (cotton, *N. benthamiana*, *B. hybridum*, *E. rockii*), the T2T donkey, six fungi and parasites, and five holocentric
  controls, with CENH3 or CENP-A ChIP-seq or CUT&Tag of the assembled line (the donkey's from other donkeys,
  *Toxoplasma*'s from another strain) or published centromeres. 0.16.1's rules were fixed before any of them was scored,
  so its results on them measure how cenfuse does on genomes it was not developed on: `high` right 104 of 110 times,
  `medium` 23 of 31, `low` 57 of 83 (development genomes: 453 of 459, 137 of 155, 149 of 197). Most confident misses
  come from centromeres that sit away from the genome's main satellite (*N. benthamiana*) and from `medium` calls on a
  CRM core (pepper, *Sesbania*). The holocentric controls got no confident call.
- AT troughs at chromosome ends no longer place a call alone: subtelomeres are AT-rich in many genomes, and a trough
  cannot tell an acrocentric centromere from a subtelomere. Where a chromosome's deepest trough is centred in its
  terminal 5%, its deepest trough whose run lies clear of the ends places the call (`at_interior`; `at_start`, `at_end`,
  `at_gc` give that trough), or none (`at_terminal`, and the caveat says why), unless Hi-C tested against the troughs
  agrees with the terminal trough on that chromosome. *Pneumocystis carinii*: 9 to 14 of 17 right; no call of the
  development benchmark changes.
- Tested on the added genomes and not adopted: moving a `low` call to the repeat-family window that contradicts it
  (registered beforehand; it held in one of the two lineages it touched), and holding weak CRM cores to the 10-fold of
  other lineages (measured after *Sesbania*'s result: such cores were right more often than the arrays they replaced)
  (future.md #76, #46, #79).
- The overview figure is rewritten for biologists: one centromere per chromosome with its confidence, every signal on
  its own track with the place where it alone puts the centromere, each marked used, points elsewhere or not used, and
  one plain sentence per close-up on how the call was made. Each signal's genome test is shown in words (trusted, not
  trusted here, a lead only). New `asm.closeups.png`.

## 0.16.1 (2026-11-14)

- No evidence confirms itself: every test that validates one line by its agreement with another (Hi-C on the LTR cores,
  the cis signal, the arm fold and trans alignment, AT troughs and repeat-density peaks tested against Hi-C) now acts on a
  chromosome only at the weaker of the genome's verdict and the verdict on the other chromosomes alone, as the Hi-C gate
  always did (`loo.py`). Per-chromosome `at_gate`, `fold_gate` and `align_gate` show the status the line had on that
  chromosome. A genome of three chromosomes can no longer validate such a test (*S. pombe*).
- The repeat-family test is not calibrated: with each genome's repeat units moved to random places on each chromosome
  it still passes 13-52% of the time at its 1% level (11 genomes, 133 of 400 draws), because one family that clusters
  once per chromosome passes wherever its clusters lie; no calibrated null could be built from K-mers. Its peaks are now
  leads: it no longer pairs with a Hi-C call (`fam_pair` removed; it had raised no call in the 0.16.0 benchmark), still
  places `low` calls nothing else places and flags `low` calls it contradicts, and the report says the test is not
  calibrated and no longer counts it as a passed test.
- Inputs fail clearly instead of silently: 4DN `.pairs` with `chr1`/`chr2` headers are read; genes as GTF (or GFF3
  without `gene` lines) feed the methylation line; empty, duplicated, truncated or mis-compressed Hi-C and methylation
  files stop with a message naming the file (exit 2); a methylation exclusion naming no sample stops the line; sequence
  names such as `NA` are kept; 5hmC rows of unfiltered modkit files are no longer counted as cytosines (methylation
  caches re-bin once); a moved output folder resumes; the satellite search deletes its scratch (GBs per genome). The
  conda environment excludes Anaconda's `defaults` channel.
- Reports and figures: repeat-line calls get their own band and `rep`/`fam` columns on the overview page and their own
  caveat sentence; a Hi-C call raised by a second line names it; GFF3 sequence names are percent-encoded; `cen_point` is
  1-based; the AT and repeat lines do not run on assemblies of more than 200 sequences (their test cannot).
- Benchmark: the 51 runs these changes can reach (and controls of every input type) were rerun; 26 chromosome-runs
  changed, all right calls that lost confidence (20 one level) or their call (*S. pombe* chr II, two runs); no wrong call
  changed, and the tomato confident errors the leave-one-out rule was aimed at remain (their tests pass without them).
  Over all 63 genomes with truth, each with its fullest inputs: `high` 453 of 459 right, `medium` 137 of 155, `low` 149
  of 197. Docs: a short README, `docs/inputs.md`, this changelog, the per-chromosome results in `benchmark/results.tsv`
  with the script that draws the new performance figure, and figure scripts in `docs/img/src/`.

## 0.16.0 (2026-11-13)

- Two new lines read repeats from the FASTA alone, where most chromosomes have no confident sequence call. Repeat
  density (`repline.py`, stage `repeats`): each chromosome's densest stretch of 21-mers repeated elsewhere. Repeat
  families (`famline.py`, stage `repfam`): sampled 31-mers that cluster at one locus on most chromosomes. Each is tested
  one per chromosome (P ≤ 0.01); a passed line places a chromosome nothing else places (`rep_only`, `fam_only`, `low`)
  and is the second line of a tested Hi-C call, raising it a level (`rep_pair`, `fam_pair`). A `low` call the family window
  contradicts stays and its caveat names the window (`fam_elsewhere`). New columns `rep_*`, `fam_*`.
- CpG kinetochores: a dip is measured by two depths (level, and spread in robust SD: new column `cpg_dip_z`) and the
  valley's width, combined by minP over rotations at every offset; the gate combines the references' P values by
  Fisher's method. A call whose arrays hold no dip while another array holds one is contested, never moved
  (`cpg_elsewhere`, `cpg_alt_*`). PacBio pb-CpG-tools BEDs and `.csi`-indexed files are read; coverage holes are
  reported (`cpg_array_data`).
- Hi-C: AT troughs are tested against the cis extremes too, and a cis call with an agreeing trough is `medium`
  (`at_pair`). A cis extreme that passes both the LTR-core and the trans test moves a confident call it contradicts
  (`cis_trans_pair`, `low`); a `low` call in a terminal eighth moves to the interior extreme where its own trans peak
  lies (`cis_end_pair`, `low`). The cis line at `support` raises the calls it agrees with and moves none. A Hi-C null
  that ran but could not score a chromosome is flagged `null_degenerate`.
- The report says when no call is confident and no line passed a genome-level test.
- Benchmark: on the 260 arms scored for both releases (real inputs; the 13 tomatoes, *Ae. geniculata* and *S.
  verrucosum* apart), 2,309 of 2,624 chromosome-runs with truth right (0.15.0: 2,236), 35 confidently wrong (36). The
  repeat-family line runs in 39 FASTA-only genomes and passes in 18, its peaks right on 183 of 211 chromosomes. CpG (the
  0.16 line's own runs): the test passes in 7 of 7 primates with truth; kinetochores within 50 kb of DiMeLo-seq CENP-A on
  53 of 63.

## 0.15.0 (2026-11-09)

- LTR line: insertion age chooses between comparable density peaks (`te_age_core`; columns `core_age_from`,
  `core_alt`). A significant core holding under half of the signature's elements replaces a call that is `low` on its
  own and stays `low` (`minority_core`). A family concentrated across the chromosomes can seed the signature (pooled
  seed). A signature other than CRM decides only if its cores are specific (median core fold ≥ 10).
- AT troughs: a new test (the deepest class of troughs, read off the depths, spread more evenly over the chromosomes
  than a room-aware random deal gives; Freeman-Tukey, mid-P). A tested trough on a `low` call's competing cluster
  chooses that cluster (`at_chose_alt`, a named heuristic from mouse).
- Hi-C alone counts a prior-free peak held at `medium` only by `null_untestable` (too few chromosomes for its null).
- No benchmark data are held out: every genome is development data.
- Benchmark: on the 258 release arms, 89 more right calls than 0.14.0, confident errors unchanged; mouse from the FASTA
  1 → 18 of 20 (`low`); *C. tropicalis* with Hi-C 0 → 6 of 7.

## 0.14.0 (2026-11-08)

- No rule changed: calls equal 0.13.0's on all 260 benchmark runs.
- Safe resumes: checkpoints keyed on code, parameters and input content (adding `tes/` or `genes/` recomputes the
  satellite engine from classify on); `--force all`; an output-folder lock; atomic writes; SIGTERM ends a run cleanly.
- `run_info.json`: versions, tools, input fingerprints, CPU time and peak memory. Column `rule`: the integration rule
  that set the call.
- Faster Hi-C engine with identical numbers (the benchmark's 139 Hi-C stages: 4.6 h instead of 8.5 h).
- `cenfuse control`: randomized Hi-C, GC blocks and LTR positions.
- The benchmark's truth manifest, frozen scorer and table chain in `benchmark/`; docs/constants.md generated from the
  code; `environment.lock.txt`.

## 0.13.0 (2026-11-07)

- Without `tes/`, Hi-C's cis signal is tested by the library's own trans contacts (`cis_test` = `trans`), read
  outside the terminal eighths (`cis_blind`).
- Insertion age is a second line that lets validated Hi-C move a `high` call (`moved_with_age`, columns `age_*`).
- The support tier chooses a competing cluster with a strong clean peak; an AT trough and a clean prior-free peak
  together place a `low` call (`at_denovo_pair`).
- With fewer than 4 confident references, `low` calls at clean peaks join the Hi-C test (`hic_gate_refs` = `mixed`).
- Left-out sequences with telomere arrays at both ends are named (`left_out_telomeric`).
- Benchmark: against 0.12.0 on 258 arms, 44 more right calls and no new confident error; tomato from the FASTA and
  Hi-C 1 → 10 of 12; *C. neoformans* G1 0 → 12 of 14; moss from the FASTA and Hi-C 16 → 21 of 26.

## 0.12.0 (2026-11-07)

- LTR line at family level (`sat/signature.py`): each chromosome's signature is learned from the other chromosomes'
  cores; a lineage other than CRM decides at ≥ 10-fold (Tekay in tomato, Bryco in moss); `weak_core`.
- The satellite line's one-per-chromosome test (`sat/karyotype.py`): where the centromeric families' largest clusters
  do not fall one per chromosome, satellite calls are `low` (`not_one_per_chromosome`).
- AT troughs, a new line (`atline.py`): `at_only`, `at_pair`.
- Where the Hi-C test fails: Hi-C anchored on the LTR cores (`te_hic_pair`), the cis signal tested against them
  (`cis_*`), and which end holds a call in a terminal eighth (`hic_end`).
- `asm.line_tests.tsv`: every pair of lines tested for agreement.
- Benchmark: 400 of the 3,062 calls with truth changed and none became confidently wrong; 49 confident errors became
  `low` and 12 became right; 60 chromosomes without a call were placed right (AT troughs).

## 0.11.0 (2026-11-06)

- Faster satellite search: in models of 500 bp or more, low-complexity tracts score nothing (1.00-4.93× faster on 14
  genomes); a copy is whole at 90% of the positions that can score.
- A call on a chain of partial copies of a longer repeat is at most `medium`, `low` without an overlapping CRM core
  (`partial_chain`); a long-unit array with no support from Hi-C, `tes/` or methylation is a feature, not a call
  (`long_unit_unsupported`).
- 45S hits need nuclear context (no organelle insertions); telomere units of 13-29 bp.
- New outputs: `seq_rejected_*`, `te_core_elements`, `te_call_elements`, `seq_spacing`, `satellite/dropped_clusters.tsv`.
- Benchmark: 58 of 3,595 calls changed (mostly unsupported long-unit calls removed: 28 wrong, 8 right); confident
  errors unchanged (104).

## 0.10.1 (2026-11-03)

- Hi-C engine 0.3.0: a peak is `telomeric` only within one mappable bin of the chromosome end.
- The trans arm alignment, a second Rabl-axis readout (`align_*`, `align_pair`).
- No Hi-C support where no mappable bin separates a call from its competing cluster (`hic_alt_ambiguous`).
- Replicate libraries (rows of `hic/asm.tsv` with one `sample`) are summed. Summary: `hic_axis`, `hic_architecture`.
- Benchmark: every number rerun on 204 arms; yeasts, *Plasmodium* and *Neurospora* with Hi-C 39 of 40; five more fungi
  with Hi-C 43 of 72.

## 0.10.0 (2026-11-03)

- Hi-C bins from the chromosome lengths (5 kb to 2 Mb), so yeasts, protists and other fungi are callable.
- HiC-Pro `.allValidPairs` and Juicer `merged_nodups` read as they are.
- Support tier: Hi-C beyond chance on a minority of chromosomes raises the calls its own significant peaks agree with
  (`hic_support`); it never opposes, moves or places a call.
- Hi-C alone: in a genome with fewer than 3 reference calls a clean `high` prior-free peak places a `low` call
  (`hic_denovo_only`, `hic_denovo_elsewhere`).
- Benchmark: yeast 16 of 16, *Plasmodium* 13 of 14, *Neurospora* 6 of 7 from Hi-C, where the FASTA placed none.

## 0.9.0 (2026-11-03)

- The chromosome set ignores a gradual tail of scaffolds; `<genome>.chromosomes.txt` lists the chromosomes. A run
  stops when the set changed since its satellite checkpoints.
- FasTAN crashes: the piece is split and retried, crashes are recorded beside the seed cache, and a chunk that cannot
  be searched stops the run.
- Family models built in parallel, identical to the serial build; a per-family search cost table; progress with a
  time estimate.
- `seq_overruled_start`, `seq_overruled_end`: the satellite array a CRM core replaced.
- Benchmark: chromosome sets unchanged on all 234 benchmark and batch assemblies; on 121 polyploid assemblies with known
  chromosome lists, extra sequences 5,614 → 31.

## 0.8.0 (2026-10-30)

- Satellites: interspersed repeats are never satellites (role `dispersed`); dominance judged per satellite kind;
  acrocentric karyotypes; competing clusters reported (`seq_alt_*`) and `low`; a family needs a seed array of ≥ 3
  copies. Reproducible (FasTAN with `setarch -R`, single-threaded MAFFT refinement) and lighter (one pyhmmer pipeline
  per thread).
- Hi-C: distances in mappable bins; the gate also needs a majority of the whole reference set; the arm fold, a
  Rabl-axis line that can move calls where the gate fails (`fold_pair`).
- LTR cores need ≥ 3 independent insertion sites; a CRM core replaces a confident satellite call only where its
  elements sit there significantly more (otherwise the call stays, `te_discordant`; new flag `satellite_beside_core`);
  DANTE_LTR read.
- CpG methylation: the kinetochore inside each call. 45S rDNA arrays (NORs); short telomeres at sequence ends;
  interstitial telomeric arrays.
- The call in words (`basis`, `caveat`, `cenfuse_report.md`); morphology and completeness; GFF3 with Sequence Ontology
  types.
- Benchmark: every arm rerun from scratch (156 runs); `high` calls right on 112 of 113 with CENP-A truth; no confident
  call on 40 holocentric chromosomes.

## 0.1.0 (2026-10-19)

- First release: satellite arrays and CRM retrotransposon cores, Hi-C anchored on the other chromosomes' calls and
  tested before it acts (the gate), CHG methylation calibrated on the other chromosomes' calls, one call per chromosome
  with the lines for and against it, other features (telomeres, knobs, subtelomeric blocks, rDNA, contact hubs), and
  provenance and per-plant scoring for inputs from several plants.
- Benchmark (criteria fixed before the first run): the integrated calls right at least as often as the sequence calls
  in every genome; the gate passed where anchored Hi-C was right on 67-85% of chromosomes and failed where it was right
  on 0-33%; rice: anchored Hi-C 8 of 12 on CENH3, prior-free Hi-C 0 of 12.
