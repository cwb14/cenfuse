"""cenfuse's own constants: one place, one reason each (design: docs/design.md). The engines keep theirs
(cenfuse.sat.params, cenfuse.hic.params, cenfuse.meth.params); none is a command-line option."""

HIC_TOL_BINS = 1          # a Hi-C peak agrees with an interval within one bin: cenhic's placement error
HIC_FAR_BINS = 4          # ... and points to another locus only beyond cenhic's partner half-window (HALF_WIN, the
                          # resolution the method claims): in between, Hi-C neither supports nor opposes the call
GATE_P = 0.01             # Hi-C gate: Poisson-binomial tail of agreements with the other chromosomes' sequence calls
                          # (conventional level; cenhic P_HIGH; the benchmark's pre-registered C2 rule)
GATE_MIN_REF = 3          # fewer reference chromosomes give no majority test (a call has < 2 partners, cenhic)
ANCHOR_LEVELS = ("high", "medium")   # calls good enough to anchor Hi-C partner windows and methylation calibration
                          # (wrong anchors blur both); low calls are used only when fewer than GATE_MIN_REF qualify
STRONG_HIC_LEVEL = "high"            # an anchored Hi-C peak that may move a call passes all three cenhic tests
                                     # (bootstrap in core >= 0.95, leave-one-partner-out stable, null P <= 0.01)
NOT_STRONG_FLAGS = ("hub", "rDNA", "telomeric")   # peaks on contact hubs or NORs (they contact everything), or in
                          # an end window (8 bins): chromosome-end clusters are the main non-centromeric trans attractor
                          # (cenhic benchmarks: rice 12/12 and giraffe 15/15 prior-free peaks at ends; 33 mammal misses).
                          # They may support a call, never oppose or move one. Added 2026-10-16 after the rice run
                          # (docs/design.md, Changes after the first runs)
HIC_MAX_THREADS = 8       # Hi-C bootstrap threads: each holds ~3 copies of the contact matrix (memory, not speed)
METH_SCORE = 0.5          # centromere_like = signature score >= 0.5 (cenmeth SIG_FRACTION: nearer the anchors than
                          # the background)
CONF_LEVELS = ("low", "medium", "high")   # cenrep's confidence levels, in order
MAX_RAISE = 2             # supporting lines raise a kept sequence call at most two levels (low -> high)
SAT_SUPPORT_BP = 5_000    # the satellite line supports a call holding >= 5 kb of centromeric-role arrays: cenrep's
                          # weak_satellite floor (fallback_min_bp = S_MIN / 4)
LOCI_MIN_SAT_BP = 20_000  # satellite clusters entered as loci: cenrep's S_MIN (centromere-sized cluster prior)
LOCI_CLUSTER_GAP = 500_000  # arrays of one family closer than this are one cluster: cenrep's cluster_gap

# Provenance and plants (docs/plants.md)
GENO_FLANK = 300          # bp each side of a tested CHG site: ONT reads are kb long, so rows on both sides mean reads
                          # span the site
GENO_MIN_READS = 3        # mean reads per neighboring row on each strand: all of them failing at a true C ~0.1^3
GENO_WINDOW = 2_000_000   # bp per genotype window, two per chromosome on the arms: ~150,000 tested sites at 2-3x
                          # coverage (P. annua probe)
GENO_ALPHA = 0.01         # private-difference test: conventional level, Bonferroni over ordered sample pairs
GENO_POWER = 0.95         # a 'reference' verdict needs this power against the smallest private rate seen between plants
KARYO_WINDOW = 1_000_000  # coverage window for missing segments (at least the Hi-C bin)
KARYO_ABSENT = 0.1        # missing: below a tenth of the input's median coverage, a tenfold drop (heuristic)
KARYO_PRESENT = 0.5       # ... where another input of the same data type has at least half its median (heuristic)
KARYO_MIN_WINDOWS = 3     # consecutive windows for a missing segment (heuristic)
PLANT_Q = 0.05            # Benjamini-Hochberg level of the between-plant tests (conventional)
