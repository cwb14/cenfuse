# Constants

Every numeric constant cenfuse decides with, as the code states it: module, name, value, reason, and a mark for the named heuristics, the absolute priors and figure layout. Generated from the code by `python -m cenfuse.consts docs/constants.md` (cenfuse/consts.py); do not edit by hand. The curated summary is [design.md, Constants](design.md#constants).


## `atline.py`

| constant | value | reason | mark |
|---|---|---|---|
| `BIN` | `1_000` | GC per 1 kb: the finest window (Malassezia's 3-5 kb centromeres need it), and 1 MB of memory per 500 Mb of sequence at 2 x 2 bytes per bin |  |
| `SCALES` | `(1_000, 2_000, 5_000, 10_000, 20_000, 50_000, 100_000, 200_000, 500_000)` | the known regional centromeres span them; each step about doubles |  |
| `MIN_WINDOWS` | `10` | a chromosome takes part at a size when it holds at least 10 windows (a trough needs flanks) |  |
| `MIN_CALLED` | `0.8` | a window counts with >= 80% of its bases called (not N) |  |
| `APART` | `5` | independent troughs lie >= 5 windows apart: one regional centromere gives several troughs closer than that (C. neoformans at 10 kb: 13 of the top 14 lie on centromeres, on 9 chromosomes) |  |
| `CLASS_MAX` | `1.5` | the deepest class of troughs holds at most 1.5 N of them: a chromosome's one centromere and a minority of others (chromosome ends, AT-rich at both ends, make two per chromosome: Chlamydomonas; heuristic) | heuristic |
| `PACKING` | `0.7476` | Renyi's parking constant: independent troughs >= r windows apart found deepest first fill a stretch at 0.7476 of the closest possible packing (checked: 1.4-1.7 x n / (2r + 1) on real chromosomes) |  |
| `AT_DRAWS` | `10_000` | null deals per window: P resolution 1e-4, as the karyotype test's |  |
| `MAX_CHROMS` | `200` | more sequences than this are contigs, not chromosomes: no test (the deal of many short sequences passed noise genomes of 300 contigs: no room to tell one trough per chromosome from one per contig) |  |
| `AT_ALPHA` | `0.05` | conventional level, after the correction for the sizes tried. The line alone gives `low` calls |  |
| `MIN_AGREE` | `3` | agreements needed to speak of a test at all (cenfuse's GATE_MIN_REF) |  |
| `CHUNK` | `10_000_000` | bp read at a time: streaming, never a whole chromosome string above this |  |

## `cdr.py`

| constant | value | reason | mark |
|---|---|---|---|
| `BIN` | `5_000` | CDR-Finder's window (Mastrorosa et al. 2024): >= 4 bins across the smallest CDRs (~20 kb) |  |
| `SMOOTH` | `3` | bins pooled into each bin's methylation (15 kb): a lone 5-kb bin of a sparse track is noise |  |
| `MIN_SITES` | `5` | CpG rows in a SMOOTH window for its methylation to count (heuristic; alpha satellite holds ~30-60 per 5 kb, so only bins with almost no data are dropped) | heuristic |
| `SURR` | `500_000` | the background: running median over +-SURR, wider than CDRs (20-200 kb), narrower than arrays |  |
| `FRAC` | `0.5` | a dip: methylation below half its background (CDRs fall to 10-35% of their arrays; survey) |  |
| `MIN_RUN` | `3` | bins in a dip (15 kb): the narrowest CDRs are ~20 kb |  |
| `ARRAY_COVER` | `0.5` | a bin is array sequence when arrays cover most of it (majority rule) |  |
| `MIN_BINS` | `10` | array bins with data for a call to be tested (50 kb: a dip needs shoulders; heuristic) | heuristic |
| `N_ROT` | `200` | rotations per call: P resolution 0.005, ten times finer than DIP_P (every offset when fewer) |  |
| `DIP_P` | `0.05` | a significant dip (conventional level); each reference's chance under the gate's null |  |
| `BRIDGE` | `2` | bins of non-array sequence an array run spans (an alternative dip's array; heuristic) | heuristic |
| `HOLE` | `0.8` | arrays with data in < HOLE x the genome's share of bins: a coverage hole (MAPQ-60 tracks: 0.36-0.68, MAPQ >= 1: >= 0.97, four truth genomes) |  |
| `SEED` | `7` | fixed seed: reruns give identical P values |  |
| `PBCPG_COLS` | `9` | pb-CpG-tools: score, type (Total/hap1/hap2), coverage (6), estimated modified count (7) |  |

## `controls.py`

| constant | value | reason | mark |
|---|---|---|---|
| `CHUNK` | `2_000_000` | pairs read and written at a time. hic-repair numbers the read ends by these blocks (each block's first ends, then its second ends), as permute_pairs.py did: a seed gives its controls only with this value |  |
| `DISP_BINS` | `1024` | the finest genome-wide bins of the dispersion check: a power of 2, so that halving keeps the bins nested; 1024 bins make 524,800 bin pairs |  |
| `DISP_MIN_E` | `5` | a bin pair is tested when its count's null variance is >= 5 (Cochran's rule for Pearson's statistic: >= 5 expected for Poisson counts) |  |
| `DISP_MIN_CELLS` | `10` | fewer testable bin pairs: untested, with a warning |  |
| `DISP_Z` | `10` | clumped above 1 + 10 null sd (the sd of the index for independent counts; the spread measured is 0.85-1.1x it). Random controls sit at 1: 40 seeds of C. tropicalis Hi-C (11 M pairs, 1024 bins) gave 0.9999 +- 0.0017 re-paired and 1.0002 +- 0.0047 relocated, at most +2.2 sd; 200 seeds of a 3,000-pair fixture at most +3.0 sd; the 14 stored NEG2 controls -1.4 to +0.5 sd, the stored NEG of mouse MMU8 +1.7 sd. The clumped `shuf` partners (v11/data/neg) scored 176-1,712, over 88,000 sd |  |
| `GC_BLOCK` | `100_000` | at_shuffle.py's blocks (100 bins of 1 kb): the benchmark's GC control |  |
| `FASTA_WIDTH` | `60` | bases per FASTA line written (samtools faidx's) |  |
| `COPY` | `10_000_000` | bp read at a time when a sequence is copied |  |

## `fold.py`

| constant | value | reason | mark |
|---|---|---|---|
| `CI_MIN` | `0.125` | Levan et al. 1964: centromeric index < 12.5% = telocentric ('t'), no second arm to fold |  |
| `N_MIN` | `3` | heuristic: fewest anti-diagonal pairs whose mean scores a center | heuristic |

## `hic/params.py`

| constant | value | reason | mark |
|---|---|---|---|
| `MAPQ` | `30` | conventional unique-mapping threshold (applied when the pairs carry MAPQ columns) |  |
| `BIN_LADDER` | `(5_000, 10_000, 50_000, 250_000, 500_000, 1_000_000, 2_000_000)` | bin sizes, finest first (bp); each a multiple of every finer one, so coarsening the base matrix is exact (bins.coarsen); which ones a genome may use: bins.candidate_sizes |  |
| `STD_FLOOR` | `250_000` | finest bin while every chromosome has MIN_BINS bins at it (cenfuse 0.9's ladder): bins finer than this only make small chromosomes (yeasts, protists, bird microchromosomes) callable, and are never used where they are not needed |  |
| `MAX_CHROM_BINS` | `4_000` | a finer candidate must keep the largest chromosome <= 4,000 bins (its dense cis block, fold.py: 128 MB) ... |  |
| `MAX_GENOME_BINS` | `50_000` | ... and the genome <= 50,000 bins (sparse matrix, bootstrap copies; 0.9 ran 3-4-Gb genomes at 250 kb, 12,000-16,000 bins) |  |
| `MIN_EXPECTED_CEN` | `25` | expected centromere contacts per bin: a 2x rise is then >= 5 SD over Poisson noise |  |
| `HALF_WIN` | `4` | partner window = estimate +/- 4 bins: a 1-2 Mb CENH3 domain plus one bin of error |  |
| `END_WIN` | `8` | telomere-pole window at each chromosome end (bins), about one partner window |  |
| `END_TOL` | `1` | a call is 'telomeric' when its peak cannot be told from a chromosome end: within one mappable bin of the terminal bin (the placement error; cenfuse's HIC_TOL_BINS) ... |  |
| `MIN_BINS` | `20` | callable sequence: two end windows plus an interior |  |
| `MIN_CALLABLE` | `3` | fewer than 3 chromosomes leave a call fewer than 2 partners (no mean over partners, no leave-one-out) |  |
| `FANIN` | `16` | per-chunk matrices held before a partial tree-sum (bounds memory while reducing pairs) |  |
| `ICE_MIN_FRAC` | `0.2` | ICE masks bins below 0.2x median coverage (standard practice) |  |
| `ICE_TOL` | `1e-3` | ICE stops when the largest bias change is below this |  |
| `ICE_MAX_ITER` | `200` | ICE iteration cap |  |
| `START_SMOOTH` | `5` | running mean (bins) of the telomere-pole score for the start (heuristic) | heuristic |
| `MAX_ROUNDS` | `50` | cap on location rounds |  |
| `CYCLE_SPREAD` | `1` | a cycle may alternate between neighboring bins of a flat-topped peak; a larger jump of a chromosome's estimate within the cycle = its location is unstable (flag) |  |
| `N_BOOT` | `100` | Poisson bootstrap replicates |  |
| `N_NULL` | `200` | null replicates (P resolution 0.005) |  |
| `NULL_EXCL` | `20` | null partner windows >= 20 bins from that partner's centromere and its ends |  |
| `NULL_REACH` | `20` | the null may reach the observed height within 20 bins of the call |  |
| `P_HIGH` | `0.01` | conventional evidence level |  |
| `P_NONE` | `0.05` | conventional evidence level; above it no call is made |  |
| `BOOT_IN_CORE` | `0.95` | share of bootstrap peaks inside the core for 'stable' |  |
| `RABL_P` | `0.01` | Rabl sign test across chromosomes |  |
| `HUB_FOLD` | `2.0` | hub flag: coverage and trans contacts >= 2x the genome median (heuristic) | heuristic |
| `SEED` | `7` | fixed seed: reruns give identical confidence numbers |  |
| `RABL_MIN_BINS` | `5` | fewest finite bins for a chromosome's Spearman rho (fewer cannot rank-correlate) |  |
| `PAIRS_CHUNK` | `2_000_000` | pair lines parsed at once (bounds memory while streaming the pairs) |  |
| `CACHE_FP_BYTES` | `1 << 20` | cache key reads the first and last MiB of each input (a fast fingerprint, not a checksum) |  |
| `MIN_CIS_DIST` | `1_000` | cis pairs < 1 kb are dropped (usual Hi-C filter: self-ligation, undigested); simulated mates are therefore >= 1 kb apart (log-uniform from 1 kb to the chromosome length) |  |
| `SIM_MAXLEN` | `150` | simulated reads are capped at 150 bp (typical Hi-C read; bounds the read matrix) |  |
| `SIM_DEFAULT_LEN` | `150` | read length assumed when no FASTQ is given (stated in the output) |  |
| `SIM_CHUNK` | `200_000` | pairs simulated at once: bounds memory (a pairs x 150-byte read matrix) at any -n |  |
| `SIM_ERR` | `0.005` | substitution rate of simulated reads (0.5%, a typical Illumina error rate) |  |
| `SIM_TRUTH_TOL` | `1_000` | a simulated mate counts as correctly placed within 1 kb of its true position |  |
| `SIM_FASTQ_READS` | `200_000` | FASTQ reads sampled for the read-length distribution (ample for an empirical one) |  |

## `hicfuse.py`

| constant | value | reason | mark |
|---|---|---|---|
| `END_ZONE` | `0.125` | end_test: Levan's telocentric bound (fold.CI_MIN) |  |
| `END_SHARE` | `0.05` | ... the terminal 5% of a chromosome's scored bins stands for an end (a heuristic, about one | heuristic |
| `END_MIN_BINS` | `4` | partner half-window at the bins chosen); at least the half-window itself |  |
| `CIS_SIM` | `200` | cis_chance: contact maps simulated per chromosome (resolution 0.005, as the engine's null) |  |
| `FC_PERM` | `1000` | permutations: P resolution 0.001 |  |
| `FC_MIN` | `3` | instances needed, on >= 2 chromosomes (fewer give no trans pair to average) |  |
| `FC_MAX_INST` | `400` | larger classes are tested on a fixed random sample of this many (cost; the mean over ~80,000 trans pairs is already precise) |  |
| `FC_MAX_BINS` | `15_000` | dense genome-wide O/E above this many bins is too large; the test is skipped |  |

## `inputs.py`

| constant | value | reason | mark |
|---|---|---|---|
| `LEFT_OUT_MAX` | `100` | the longest sequences outside the chromosome set read for left_out (two telomere.TERM_SEARCH reads each) |  |

## `meth/bins.py`

| constant | value | reason | mark |
|---|---|---|---|
| `CACHE_VERSION` | `1` | the class-bin cache's format: bump when its content changes (fingerprint) |  |

## `meth/io.py`

| constant | value | reason | mark |
|---|---|---|---|
| `MODKIT_COLS` | `18` | modkit pileup bedMethyl: N_valid (column 10), N_mod (12), N_other_mod (14) |  |
| `ENCODE_COLS` | `11` | ENCODE bedMethyl (WGBS pipelines, Megalodon): coverage (10), percent methylated (11); modkit's first 11 columns are the same layout |  |

## `meth/params.py`

| constant | value | reason | mark |
|---|---|---|---|
| `MIN_CALLS_SITE` | `3` | min N_valid per cytosine: per-site usability floor for ~4x samples |  |
| `BIN` | `50_000` | bin (bp): >= 20 bins per CENH3 domain; validated on P. annua |  |
| `DEPTH_LO` | `0.5` | drop a bin below half the sample's median depth: one-copy loss, missing DNA (heuristic) | heuristic |
| `DEPTH_HI` | `2.0` | drop a bin above double the median depth: duplication, collapsed repeat (heuristic) | heuristic |
| `WIN_D` | `0.4` | sliding window = 0.4 D: narrower than the domain (heuristic) | heuristic |
| `RING_LO_D` | `1.5` | class reference ring starts 1.5 D away: beyond the domain (heuristic) | heuristic |
| `RING_HI_D` | `4.0` | ...and ends 4 D away: wide enough for stable class levels (heuristic) | heuristic |
| `MIN_LEN_D` | `10` | sequences >= 10 D are analyzed: a 4-D ring each side of a 2-D center |  |
| `GAP_FOLD` | `10` | chromosomes vs unplaced scaffolds (heuristic): the first tenfold drop in sorted sequence length after half the assembly ends the set (cenfuse 0.8, track.tenfold_count); in 0.9 (track.chromosome_count) a sequence under 1/GAP_FOLD of the median length before it is small, and small ones spanning a GAP_FOLD range form a gradual scaffold tail, cut off | heuristic |
| `D_REF_BP` | `22.2e6` | summed CENH3 domain length of maize Mo17 (2.22 Mb x 10 centromeres; Chen et al. 2023) |  |
| `D_REF_G` | `2.18e9` | maize Mo17 genome size (bp) |  |
| `D_SLOPE` | `0.97` | CENH3 area scales with genome size^0.97 (Plackova et al. 2021 Sci Rep) |  |
| `MIN_ANCHORS` | `3` | calibration needs >= 3 anchored chromosomes with data |  |
| `N_PERM` | `10_000` | random placements: P resolution 1e-4, so Benjamini-Hochberg over up to ~500 class x context tests can still pass a lone true class |  |
| `FDR` | `0.05` | Benjamini-Hochberg level (conventional) |  |
| `MIN_CALLS` | `500` | calls needed for a class (or the track) in a window: deviation SE ~0.02 |  |
| `MIN_EFFECT` | `0.04` | a signature entry's mean deviation at the anchors must be >= 2 x the ~0.02 SE of one candidate window at MIN_CALLS calls, so one window can tell it from zero |  |
| `SEP_D` | `4` | candidates, random placements, 'near an anchor': 4 D apart, beyond the ring's inner edge |  |
| `MAX_CAND` | `5` | candidates per sample x chromosome: enough to hold the true centromere without flooding (heuristic) | heuristic |
| `SIG_FRACTION` | `0.5` | nearest of two references, used twice: a signature class is concentrated at the anchors (genome share <= 0.5 x anchor share), and a centromere_like candidate is closer to the anchors' signature (score 1) than to the background (score 0) |  |
| `SPECIFIC_P` | `0.05` | plant-specific flag (conventional level; per chromosome, not corrected for multiple tests) |  |
| `LOW_DEPTH` | `8.0` | warn below 8x median depth: few repeat cytosines reach MIN_CALLS_SITE |  |
| `CHUNK` | `10_000_000` | bedMethyl region read at once (bp): bounds memory. Each worker holds one chunk of rows (~0.5 GB CHG, ~1 GB CHH at plant densities), so a run needs about threads x 1-2 GB |  |
| `SEED` | `7` | fixed seed: reruns give identical P values |  |
| `CACHE_FP_BYTES` | `1 << 20` | cache key reads the first and last MiB of each input (a fast fingerprint, not a checksum) |  |

## `meth/plot.py`

| constant | value | reason | mark |
|---|---|---|---|
| `PER_PAGE` | `8` | layout: chromosome panels per page | layout |

## `morphology.py`

| constant | value | reason | mark |
|---|---|---|---|
| `GAP_MIN_N` | `10` | an assembly gap is a run of >= 10 N: the satellite engine's min_gap_n (assemblers write gaps as runs of >= 10 N; shorter runs are ambiguity codes) |  |
| `SCAN_CHUNK` | `10_000_000` | bp read at once in the gap scan (memory) |  |

## `params.py`

| constant | value | reason | mark |
|---|---|---|---|
| `HIC_TOL_BINS` | `1` | a Hi-C peak agrees with an interval within one bin: the Hi-C engine's placement error |  |
| `HIC_FAR_BINS` | `4` | ... and points to another locus only beyond the Hi-C engine's partner half-window (HALF_WIN, the resolution the method claims): in between, Hi-C neither supports nor opposes the call |  |
| `GATE_P` | `0.01` | Hi-C gate: Poisson-binomial tail of agreements with the other chromosomes' sequence calls (conventional level; Hi-C engine P_HIGH; the benchmark's pre-registered C2 rule) |  |
| `FC_LEVEL` | `0.01` | feature classes in 3D: a class self-associates at P <= this (conventional level, as GATE_P) |  |
| `GATE_MIN_REF` | `3` | fewer reference chromosomes give no majority test (a call has < 2 partners, Hi-C engine) |  |
| `AXIS_LEAN_P` | `0.05` | the trans arm alignment stands in for a failed arm fold only where the fold's own test leans the same way at the conventional evidence level (Hi-C engine P_NONE): a Rabl axis shows in both readouts, a symmetric compartment pattern (holocentrics) only in the trans one |  |
| `HIC_MAX_THREADS` | `8` | Hi-C bootstrap threads: each holds ~3 copies of the contact matrix (memory, not speed) |  |
| `METH_SCORE` | `0.5` | centromere_like = signature score >= 0.5 (methylation engine SIG_FRACTION: nearer the anchors than the background) |  |
| `MAX_RAISE` | `2` | supporting lines raise a kept sequence call at most two levels (low -> high) |  |
| `SAT_SUPPORT_BP` | `5_000` | the satellite line supports a call holding >= 5 kb of centromeric-role arrays: the satellite engine's weak_satellite floor (fallback_min_bp = S_MIN / 4) |  |
| `LOCI_MIN_SAT_BP` | `20_000` | satellite clusters entered as loci: the satellite engine's S_MIN (centromere-sized cluster prior) | prior |
| `LOCI_CLUSTER_GAP` | `500_000` | arrays of one family closer than this are one cluster: the satellite engine's cluster_gap |  |
| `GENO_FLANK` | `300` | bp each side of a tested CHG site: ONT reads are kb long, so rows on both sides mean reads span the site |  |
| `GENO_MIN_READS` | `3` | mean reads per neighboring row on each strand: all of them failing at a true C ~0.1^3 |  |
| `GENO_WINDOW` | `2_000_000` | bp per genotype window, two per chromosome on the arms: ~150,000 tested sites at 2-3x coverage (P. annua probe) |  |
| `GENO_ALPHA` | `0.01` | private-difference test: conventional level, Bonferroni over ordered sample pairs |  |
| `GENO_POWER` | `0.95` | a 'reference' verdict needs this power against the smallest private rate seen between plants |  |
| `KARYO_WINDOW` | `1_000_000` | coverage window for missing segments (at least the Hi-C bin) |  |
| `KARYO_ABSENT` | `0.1` | missing: below a tenth of the input's median coverage, a tenfold drop (heuristic) | heuristic |
| `KARYO_PRESENT` | `0.5` | ... where another input of the same data type has at least half its median (heuristic) | heuristic |
| `KARYO_MIN_WINDOWS` | `3` | consecutive windows for a missing segment (heuristic) | heuristic |
| `PLANT_Q` | `0.05` | Benjamini-Hochberg level of the between-plant tests (conventional) |  |

## `plants.py`

| constant | value | reason | mark |
|---|---|---|---|
| `MIN_NULL` | `5` | call-site values needed for the null spread (4 degrees of freedom); heuristic | heuristic |

## `plot.py`

| constant | value | reason | mark |
|---|---|---|---|
| `MM` | `1 / 25.4` | layout: inches per millimeter | layout |
| `B_COL_MM` | `4.6` | panel B: spacing of the line columns | layout |
| `ROWS_PER_PAGE` | `24` | layout: chromosome rows per page | layout |
| `OVERLAY_BINS` | `40` | panel D spans +/- 40 Hi-C bins, at most 10 Mb (10 Mb at 250-kb bins, as before 0.10) | layout |

## `provenance.py`

| constant | value | reason | mark |
|---|---|---|---|
| `MOTIF_MIN` | `0.5` | rows on reference motifs below this share: the file belongs to another assembly (error) |  |
| `MOTIF_WARN` | `0.99` | ... below this: a warning (a few rows off-motif are possible at assembly edits) |  |
| `END_MARGIN` | `0.05` | genotype windows stay out of the outer 5% of each chromosome (subtelomeric repeats) |  |
| `MAX_TRIES` | `12` | windows tried per target position before giving up (gap-rich chromosomes) |  |
| `GAP_SHARE` | `0.1` | a karyotype window with more N than this is not judged |  |
| `SAT_SHARE` | `0.5` | ... nor one mostly satellite array (mapping, not DNA, sets its coverage); heuristic | heuristic |
| `KEY` | `1 << 34` | site key = chromosome index x KEY + position (chromosomes < 17 Gb) |  |
| `MOTIF_SPAN` | `200_000` | bp of each window read for the assembly check (CHH files hold ~10x the rows of CHG) |  |

## `rdna.py`

| constant | value | reason | mark |
|---|---|---|---|
| `UNIT_MAX` | `100_000` | prior: 45S units are ~8-45 kb (plants, fungi ~8-11; mammals 43-45); hits farther apart than twice the longest unit are separate loci | prior |
| `MIN_UNITS` | `3` | a tandem array has >= 3 copies (the satellite engine's min_array_monomers) |  |
| `CHUNK` | `10_000_000` | bp fetched at once (bounds memory on Gb chromosomes) |  |
| `FIVE_S_E` | `1e-5` | per-gene E-value scaled to the chromosome (the satellite engine's nhmmer_evalue) |  |
| `FIVE_S_UNIT_MAX` | `20_000` | prior: 5S units are 0.2-2.5 kb (plants 0.3-1 kb; mammals and birds 2.2 kb), and a 5S gene sits in each 9.1-kb rDNA unit of budding yeast; genes farther apart than twice that are separate loci | prior |

## `runtime.py`

| constant | value | reason | mark |
|---|---|---|---|
| `WHOLE_FILE` | `32 << 20` | files up to 32 MiB are hashed whole (an .fai, a GFF of genes, most tables: well under a second) |  |
| `SAMPLE` | `1 << 20` | larger ones (pairs, LTR annotations, bedMethyl): size + four 1-MiB samples at fixed offsets |  |
| `HEARTBEAT` | `60` | s between the owner's touches of a lock it holds |  |
| `STALE_AFTER` | `15 * 60` | a lock whose owner cannot be checked (another host, a garbled lock) is stale after 15 missed heartbeats (heuristic) | heuristic |
| `GUARD_STALE` | `60` | s: a takeover guard older than this was left by a process that died while taking a lock over |  |
| `SAMPLE_EVERY` | `30` | s between samples of the process tree's memory (a peak shorter than this can be missed) |  |
| `LISTENER_STOP` | `10` | s to wait for the log forwarding to stop: a worker killed while it held the queue's lock (SIGTERM, the OOM killer) leaves it blocked for good, and a run must not hang on that |  |

## `sat/families.py`

| constant | value | reason | mark |
|---|---|---|---|
| `UNIT_COVER` | `0.8` | heuristic: a BLAST alignment matches two units when it spans >= 80% of the shorter (frayed ends allowed); the same share of the longer unit's dimer must be covered | heuristic |
| `THREADS_PER_BUILD` | `4` | the fewest per build while families remain: MAFFT (Panax fam006, 65 monomers) is 3.6x faster at 4 threads and 5.6x at 16, same alignment at 1-16, so 4 builds of 4 beat 1 of 16 (model_workers) |  |

## `sat/karyotype.py`

| constant | value | reason | mark |
|---|---|---|---|
| `ALPHA` | `0.05` | conventional level |  |
| `N_NULL` | `10_000` | null draws: P resolution 1e-4 |  |
| `SEED` | `7` | fixed: reruns give identical P values |  |
| `MAX_SEQUENCES` | `5_000` | more sequences than this are contigs, not chromosomes: no test (and the null grows with them) |  |

## `sat/params.py`

| constant | value | reason | mark |
|---|---|---|---|
| `SEED_MAX_PERIOD` | `2000` | longest seed unit = TRF's maximum period (FasTAN's longer units are dropped); longer units are left to the long-unit detector |  |
| `FAMILY_IDENT` | `0.90` | same satellite family: seed consensuses are near-homogeneous, 90% splits variants >10% apart |  |
| `SAME_VARIANT_IDENT` | `0.97` | same monomer variant: Naish et al. 2021 HOR rule, <= 5 differences per 178 bp |  |
| `RELATED_IDENT` | `0.70` | related families (superfamily, multimer folding, shared monomer frame) |  |
| `FULL_MONOMER_COV` | `0.90` | a hit covering >= 90% of the model is a whole monomer |  |
| `S_MIN` | `20_000` | prior: smallest satellite cluster treated as centromere-sized (bp); other sizes are rungs | prior |
| `LEN_TOL` | `0.05` | length tolerance for "k copies of a unit" and "same-length" comparisons |  |
| `RELATED_LEN_TOL` | `0.15` | length tolerance for related families (subfamilies with indels) |  |
| `FULL_PERIOD_TOL` | `0.01` | a seed consensus within 1% (>= 1 bp) of its period is one whole repeat unit |  |
| `MAX_CHROMOSOMES` | `50` | more sequences than this: probably not a chromosome-level assembly (warning, plot) |  |
| `LU_MAX_UNIT` | `30_000` | longest long-unit period searched (no principled upper bound; 2 copies fit a window) |  |
| `CONS_UNITS` | `30` | copies per seed consensus: the majority base of a column is wrong with p < 1e-3 when copies are 20% divergent (binomial tail); more copies cost time, not accuracy |  |
| `BLAST_BATCH` | `200` | seed representatives per BLAST round (bounds output size) |  |
| `ANNOT_CHUNK` | `2_000_000` | nhmmer chunk (bp): keeps all threads busy; E-values rescaled to chromosome length |  |
| `MULTIMER_RATIO` | `1.9` | a model >= 1.9 x longer than another is a multimer of it (k >= 2, minus indels) |  |
| `MULTIMER_MIN_COV` | `0.5` | drop a multimer hit when shorter-unit hits cover the majority of it |  |
| `OVERLAP_MAX_COV` | `0.2` | accept a hit only if <= 20% of it is already covered by better hits |  |
| `OVERHANG_MIN_BP` | `20` | terminal redundancy shorter than this is a chance match, not a run-through |  |
| `MSA_MAX_GAP` | `0.5` | consensus keeps alignment columns occupied in the majority of sequences |  |
| `NHMMER_REPORT_E` | `10.0` | hits kept up to HMMER's default reporting E 10 per search (pyhmmer returns weaker ones too); significance is then decided per hit or per run of copies (annotate.run_significant) |  |
| `RDNA_TAG_EVALUE` | `1e-5` | cmsearch E for an rRNA gene inside a family consensus dimer |  |
| `RDNA_FLANK_EVALUE` | `1e-10` | cmsearch --hmmonly E for 45S genes in array flanks (stricter: long flank queries) |  |
| `TRACT_MIN_MODEL` | `500` | a model's other half (TRACT_MAX_FRAC) must make a whole copy significant alone: 250 scoring positions at about 1 bit each clear E 1e-5 on a 500-Mb chromosome (~47 bits) by far |  |
| `TRACT_MIN_BP` | `64` | low-complexity bp in the consensus (one DUST window) to match such stretches of the genome |  |
| `TRACT_MAX_FRAC` | `0.5` | a model mostly of low complexity stays as it is: little of it would be left to score |  |
| `TE_BINS` | `500` | density profile resolution: 1/500 of the chromosome |  |
| `TE_MIN_BIN` | `1_000` | ... but bins of at least 1 kb |  |
| `TE_SMOOTH` | `5` | smoothing window, in bins |  |
| `TE_MIN_SITES` | `3` | fewer independent insertion sites of a clade on a chromosome: no peak estimate. Nested elements (overlapping or abutting) are one site, so one nested cluster cannot pass alone; three independent insertions at one locus can (Chlorella NS4-2: 2-8 Osser elements per chromosome, 53/54 inside CENH3; a floor of 5 elements left 5/12 chromosomes) |  |
| `TE_DERIVED_FRAC` | `0.5` | a call's main family with most array bp inside LTR retrotransposons: flag te_derived |  |
| `TE_ALPHA` | `1e-3` | core significance: binomial tail x (L / core length) < 0.001 (conventional level) |  |
| `TE_CONTRA_ALPHA` | `0.01` | a core contradicts a call if the clade marks it more at P <= 0.01 (conventional) |  |
| `TE_SPECIFIC_FOLD` | `10.0` | a lineage other than CRM may decide only at >= 10-fold enrichment in its cores, and no core, CRM's included, may overrule a confident call or a satellite family below it on its chromosome. Measured on 42 polyploids whose best marker is not CRM, against 0.11's FASTA-only calls (an independent line, no truth): cores within 1 Mb of the call, median 0.30 at <= 10-fold, 0.75 above |  |
| `TE_MAJOR_SHARE` | `0.5` | a core of a lineage other than CRM is strong where it holds most of the signature's elements on its chromosome (majority rule): the rest of them say the locus is not unique there |  |
| `TE_SEED_MIN_CORES` | `3` | a CRM seed needs this many cores: fewer leave a chromosome under 2 others to learn from |  |
| `TE_RECRUIT_MIN` | `5` | elements of a family on the other chromosomes: fewer cannot reach P <= 1e-3 at a core share of a few percent (0.05^5 = 3e-7, before the Bonferroni correction) |  |
| `TE_RECRUIT_ROUNDS` | `5` | recruitment repeats until its significant cores stop changing, this many rounds at most (a bound on the time: the development genomes settle within 3) |  |
| `TE_POOLED_BREADTH` | `0.05` | a pooled family's core P on more than half of the chromosomes, each by itself: conventional level (its concentration is broad, not the work of a few chromosomes, as Fisher's P alone could be) |  |
| `TE_POOLED_MIN_CORES` | `3` | a family seeded by its pooled concentration needs this many significant cores to recruit from, as a CRM seed (TE_SEED_MIN_CORES) |  |
| `TE_PEAK_SHARE` | `0.5` | a second density peak of the signature's elements on a chromosome is comparable to the highest when at least half as high (the half-maximum convention of te.core) |  |
| `TE_AGE_ALPHA` | `0.01` | copies younger at a comparable peak than at the highest: one-sided Mann-Whitney at this level (the LTR line's level for overruling a call, TE_CONTRA_ALPHA; cenfuse.ageline) |  |
| `TE_AGE_MIN` | `3` | dated copies needed at each peak (as the three insertions that make an LTR core, TE_MIN_SITES) |  |
| `TE_PEAK_MAX` | `8` | density peaks of one chromosome examined (more than any genome's real clusters) |  |
| `TANDEM_FRAC` | `0.5` | a satellite is tandem: a family with neither the majority of its array bp inside de novo tandem repeats (seeds) nor the majority of its hit bp in whole units is a dispersed (interspersed) repeat annotated as chains of partial hits (majority rule) |  |
| `LOCUS_FOLD` | `2.0` | heuristic (as the Hi-C engine's HUB_FOLD): a non-tandem family whose main cluster is less than twice as dense as the family's average along the chromosome is spread along it, not concentrated at a locus. T2T mammal interspersed families: 1.1-1.3; plant TE-derived centromeric arrays (Coffea, Glycine CRM-derived): 2.5-3.4; satellites >= 15 | heuristic |
| `ACRO_ALPHA` | `0.01` | acrocentric test: so few chromosomes carry a terminal block at both ends that ends cannot be independent (binomial lower tail; conventional level) |  |

## `sat/plot.py`

| constant | value | reason | mark |
|---|---|---|---|
| `BIN` | `100_000` | layout: satellite bp are drawn per 100-kb bin | layout |
| `H` | `0.8` | layout: bar height (track units) | layout |
| `MAX_TRACKS` | `20` | layout: family tracks drawn per page at most | layout |

## `sat/seeds.py`

| constant | value | reason | mark |
|---|---|---|---|
| `POLL` | `15` | seconds between checks while another process holds a genome's seed lock |  |

## `sat/signature.py`

| constant | value | reason | mark |
|---|---|---|---|
| `FAMILY_MIN_SHARE` | `0.5` | a file holds families when most of its LTR retrotransposons name one (majority rule) |  |

## `sat/trf.py`

| constant | value | reason | mark |
|---|---|---|---|
| `DUP_TOL` | `0` | the same array seen from two chunks keeps its TRF period (within 1 bp) |  |

## `sat/util.py`

| constant | value | reason | mark |
|---|---|---|---|
| `_FLOAT_NOISE` | `1e-12` | 3 x 0.1 != 0.3, and sums of float work are off by a few ulps: a tenth is not missed |  |

## `telomere.py`

| constant | value | reason | mark |
|---|---|---|---|
| `TERM_SCAN` | `20_000` | bp at each end used to learn the unit (telomere arrays are kb-long at assembled ends) |  |
| `TERM_SEARCH` | `100_000` | bp at each end searched for the array (assemblies can leave a short non-telomeric tail) |  |
| `MAX_BREAK` | `30` | copies of the unit closer than this are one run (sequencing errors, variant copies) |  |
| `CHAIN_GAP` | `1_000` | runs closer than this form one array (variant-repeat stretches split real arrays) |  |
| `MIN_ARRAY` | `300` | bp of unit copies in an array: ~40 copies of a 7-bp unit, beyond chance tandem matches |  |
| `MIN_TERMINAL` | `100` | ... or this many bp when the array reaches the sequence end (within TERMINAL): a short telomere (Neurospora NC12: 157-207 bp at 13/14 ends) or one the assembler trimmed. The unit is already the genome's telomere unit (END_ENRICH), and 100 bp of one exact 5-12-bp unit is no chance match |  |
| `TERMINAL` | `1_000` | an array this close to the sequence end reaches it |  |
| `MIN_COPIES` | `4` | runs shorter than this (chance matches every ~100 bp in random sequence) are not chained |  |
| `SCAN_CHUNK` | `10_000_000` | bp read at once in the interstitial scan (memory) |  |
| `MIN_DENSITY` | `0.5` | an array is mostly unit copies: covered bp >= half of its span |  |
| `MIN_DENSITY_END` | `0.1` | ... or, once the unit is known, a tenth for a stretch that reaches the sequence end: its position is the evidence (stretches this rich occur in 1.1% of 100-kb windows inside chromosomes, 61 of 6,336 ends reach the end under half density, none below 0.19). The floor is a heuristic that keeps out copies scattered every few kb, which CHAIN_GAP would chain | heuristic |
| `MIN_ENDS` | `2` | a telomere unit forms arrays at two or more chromosome ends ... |  |
| `END_ENRICH` | `10` | ... and is >= 10-fold denser in tandem at the ends than inside the chromosomes |  |
| `INTERIOR` | `5` | interior windows (TERM_SCAN bp each) sampled per chromosome for that comparison |  |
| `SEED` | `7` | fixed: the interior sample is the same on every rerun |  |

## `sat/params.py (Params)`

| constant | value | reason | mark |
|---|---|---|---|
| `seed_engine` | `'fastan'` | tandem-repeat finder for the seeds: fastan (FasTAN + abPOA unit consensus) or trf (the --seed-engine option) |  |
| `seed_chunk` | `5000000` | seed-search chunk (bp); speed only (parallel chunks), bounded by TRF -l 6 (6 Mb) |  |
| `seed_overlap` | `10000` | overlap between seed chunks (bp); >= 2 x SEED_MAX_PERIOD so no array is lost at a chunk edge |  |
| `min_period` | `30` | shortest seed period (bp): satellites, not microsatellites (misses 5-12-bp satellites, e.g. Drosophila) |  |
| `max_period` | `2000` | longest seed period (bp) = SEED_MAX_PERIOD |  |
| `min_seed_span` | `2000` | shortest tandem array that seeds a family (bp) = one longest unit; bounds BLAST work only |  |
| `fam_ident` | `0.9` | identity for one family (seeds of equal period); see FAMILY_IDENT |  |
| `multi_ident` | `0.8` | identity when one seed is a k-mer of another: a HOR-period consensus averages divergent monomers |  |
| `superfam_ident` | `0.7` | related-family identity: superfamily table, shared monomer frame (so the HMMs), co-centromeric and dominance tests |  |
| `absorb_ident` | `0.7` | identity for folding a k-mer family into its monomer family = RELATED_IDENT |  |
| `dedup_ident` | `0.97` | refined families this identical (same length) are one variant; see SAME_VARIANT_IDENT |  |
| `chain_ident` | `0.9` | subfamilies this identical are chained into shared arrays = FAMILY_IDENT |  |
| `min_fam_span` | `20000` | prior: union seed span in >= 1 genome for a family to be annotated (bp) = S_MIN | prior |
| `max_families` | `40` | families annotated per complex; cost ceiling (annotate time grows with models x genome size) |  |
| `hmm_sample` | `300` | monomers sampled per family for MSA/HMM (MAFFT cost) |  |
| `hmm_sample_bp` | `500000` | seed-array bp scanned per genome when sampling monomers (cost) |  |
| `min_hmm_seqs` | `20` | below ~20 sampled monomers, profile emissions are noise: single-sequence HMM |  |
| `nhmmer_evalue` | `1e-05` | nhmmer E-value for monomer hits (conventional) |  |
| `array_gap_factor` | `2.0` | max gap inside an array, in monomer lengths (scales with the unit) |  |
| `min_array_monomers` | `3` | a tandem array has >= 3 copies (annotate's arrays; families: a seed cluster founds a family only with a seed of >= 3 copies, LEN_TOL short allowed) |  |
| `evalue_ref_len` | `0` | score sequences shorter than this as if this long (bp; 0 = off), so contigs get a chromosome's nhmmer cutoff |  |
| `chrom_min_bp` | `20000` | prior: min family bp on a chromosome / min main cluster (bp) = S_MIN | prior |
| `cluster_gap` | `500000` | prior: arrays closer than this form one cluster; also the co-centromeric distance (bp) | prior |
| `min_concentration` | `0.5` | the main cluster holds the majority of a chromosome's family bp (median over chromosomes) |  |
| `max_subtel_frac` | `0.5` | majority share: of a family's bp beyond the outermost genes (terminal) or near ends, and of genomes in the subtelomeric vote |  |
| `subtel_frac_len` | `0.05` | end zone (near_end flag; near-end grade; subtelomere when a sequence has < 2 gene loci), fraction of length |  |
| `max_gene_ratio` | `0.5` | centromeres are gene-poor: main clusters at >= half the genome's gene density are not centromeric |  |
| `strong_prevalence` | `0.5` | 'strong' and dominant families are centromeric on the majority of chromosomes |  |
| `strong_cluster_bp` | `100000` | prior: min median main-cluster bp for 'strong' evidence = 5 x S_MIN | prior |
| `rdna_flank` | `15000` | window each side of a tested array searched for 45S rRNA genes (bp); fits plant 45S units |  |
| `rdna_arrays` | `10` | largest arrays per eligible family tested for flanking 45S genes (cost) |  |
| `rdna_array_frac` | `0.5` | the majority of tested arrays flanked by 45S genes -> rDNA spacer repeat, not centromeric |  |
| `classify_from` | `''` | TARGET=SRC[+SRC2][,...]: TARGET (e.g. a contig assembly) takes the centromeric families of its SRC genomes (the --classify-from option) |  |
| `low_sat_bp` | `50000` | prior: centromere satellite bp below this -> low_satellite flag = S_MIN x 2.5 | prior |
| `competing_ratio` | `0.5` | a second cluster >= half the main one -> competing_cluster flag |  |
| `fallback_min_bp` | `5000` | prior: min cluster bp for a weak_satellite fallback call = S_MIN / 4 | prior |
| `cen_assoc_frac` | `0.5` | a non-primary family with the majority of its bp inside centromeres is cen_associated |  |
| `min_gap_n` | `10` | assemblers write gaps as runs of >= 10 N; shorter runs are ambiguity codes |  |
| `diffuse_min_len` | `5000000` | prior: min centromere length for the diffuse flag (bp) | prior |
| `diffuse_max_frac` | `0.2` | diffuse flag when the satellite fraction of a long region is below this |  |
| `lu_flank` | `500000` | long-unit scan flank around each centromere (bp) |  |
| `lu_window` | `1000000` | long-unit self-alignment window (bp); cost, >= 4 x lu_max_unit |  |
| `lu_window_overlap` | `60000` | overlap between long-unit windows = 2 x LU_MAX_UNIT, so no unit pair is split |  |
| `lu_min_unit` | `2000` | shortest long-unit period = SEED_MAX_PERIOD (hand-off from the seeds) |  |
| `lu_max_unit` | `30000` | longest long-unit period = LU_MAX_UNIT |  |
| `lu_min_ident` | `0.8` | min self-alignment identity for long-unit arrays: the divergence minimap2 asm20 aligns through |  |
| `contig_end_bp` | `10000` | borrowed genomes: a call this close to a contig end is flagged contig_end (bp) |  |
| `hor_ident` | `0.97` | cd-hit-est identity for monomer variant classes = SAME_VARIANT_IDENT |  |
| `hor_min_cov` | `0.9` | only whole monomers enter the HOR search = FULL_MONOMER_COV |  |
| `hor_max_len` | `200` | longest HOR unit searched, in monomers (runtime bound) |  |
| `tandem_min_len` | `500` | shortest other tandem array reported inside centromeres (bp; reporting filter) |  |
