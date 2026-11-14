# Changelog

What changed for users in each release, newest first, with the benchmark headline of the release. Current results:
[docs/benchmark.md](docs/benchmark.md). A *chromosome-run* is one chromosome in one benchmark run; *confident* = `high`
or `medium`.

## Unreleased

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
