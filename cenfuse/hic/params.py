"""Fixed internal constants, one place, one reason each (rationale in docs/design.md)."""
from __future__ import annotations

MAPQ = 30                       # conventional unique-mapping threshold (applied when the pairs carry MAPQ columns)
BIN_LADDER = (5_000, 10_000, 50_000, 250_000, 500_000, 1_000_000, 2_000_000)  # bin sizes, finest first (bp); each a
                                # multiple of every finer one, so coarsening the base matrix is exact (bins.coarsen);
                                # which ones a genome may use: bins.candidate_sizes
STD_FLOOR = 250_000             # finest bin while every chromosome has MIN_BINS bins at it (cenfuse 0.9's ladder): bins
                                # finer than this only make small chromosomes (yeasts, protists, bird
                                # microchromosomes) callable, and are never used where they are not needed
MAX_CHROM_BINS = 4_000          # a finer candidate must keep the largest chromosome <= 4,000 bins (its dense cis block,
                                # fold.py: 128 MB) ...
MAX_GENOME_BINS = 50_000        # ... and the genome <= 50,000 bins (sparse matrix, bootstrap copies; 0.9 ran 3-4-Gb
                                # genomes at 250 kb, 12,000-16,000 bins)
MIN_EXPECTED_CEN = 25           # expected centromere contacts per bin: a 2x rise is then >= 5 SD over Poisson noise
HALF_WIN = 4                    # partner window = estimate +/- 4 bins: a 1-2 Mb CENH3 domain plus one bin of error
END_WIN = 8                     # telomere-pole window at each chromosome end (bins), about one partner window
END_TOL = 1                     # a call is 'telomeric' when its peak cannot be told from a chromosome end: within one
                                # mappable bin of the terminal bin (the placement error; cenfuse's HIC_TOL_BINS) ...
END_FREE = HALF_WIN             # ... counting at most HALF_WIN unmappable bins between them as no distance (cenfuse's
                                # HIC_FAR_BINS); a longer blind stretch separates the peak from the end. 0.10 flagged
                                # every peak within END_WIN bins in bp: true centromeres lie there on acrocentric
                                # chromosomes (moss 11/26, goat 29/30, horse 17/32), while the prior-free end clusters
                                # it was made for lie in the terminal 0-1 bins (rice 12/12, giraffe 14/16)
MIN_BINS = 20                   # callable sequence: two end windows plus an interior
MIN_CALLABLE = 3                # fewer than 3 chromosomes leave a call fewer than 2 partners (no mean over
                                # partners, no leave-one-out)
FANIN = 16                      # per-chunk matrices held before a partial tree-sum (bounds memory while reducing pairs)
ICE_MIN_FRAC = 0.2              # ICE masks bins below 0.2x median coverage (standard practice)
ICE_TOL = 1e-3                  # ICE stops when the largest bias change is below this
ICE_MAX_ITER = 200              # ICE iteration cap
START_SMOOTH = 5                # running mean (bins) of the telomere-pole score for the start (heuristic)
MAX_ROUNDS = 50                 # cap on location rounds
CYCLE_SPREAD = 1                # a cycle may alternate between neighbouring bins of a flat-topped peak; a larger jump
                                # of a chromosome's estimate within the cycle = its location is unstable (flag)
N_BOOT = 100                    # Poisson bootstrap replicates
N_NULL = 200                    # null replicates (P resolution 0.005)
NULL_EXCL = 20                  # null partner windows >= 20 bins from that partner's centromere and its ends
NULL_REACH = 20                 # the null may reach the observed height within 20 bins of the call
P_HIGH = 0.01                   # conventional evidence level
P_NONE = 0.05                   # conventional evidence level; above it no call is made
BOOT_IN_CORE = 0.95             # share of bootstrap peaks inside the core for 'stable'
RABL_P = 0.01                   # Rabl sign test across chromosomes
HUB_FOLD = 2.0                  # hub flag: coverage and trans contacts >= 2x the genome median (heuristic)
SEED = 7                        # fixed seed: reruns give identical confidence numbers
RABL_MIN_BINS = 5               # fewest finite bins for a chromosome's Spearman rho (fewer cannot rank-correlate)
PAIRS_CHUNK = 2_000_000         # pair lines parsed at once (bounds memory while streaming the pairs)
CACHE_FP_BYTES = 1 << 20        # cache key reads the first and last MiB of each input (a fast fingerprint, not a checksum)

# simulate (mis-mapping null)
MIN_CIS_DIST = 1_000            # cis pairs < 1 kb are dropped (usual Hi-C filter: self-ligation, undigested); simulated
                                # mates are therefore >= 1 kb apart (log-uniform from 1 kb to the chromosome length)
SIM_MAXLEN = 150                # simulated reads are capped at 150 bp (typical Hi-C read; bounds the read matrix)
SIM_DEFAULT_LEN = 150           # read length assumed when no FASTQ is given (stated in the output)
SIM_MIN_CHROM = max(MIN_CIS_DIST, 2 * SIM_MAXLEN)   # shorter chromosomes cannot hold a >= 1 kb pair plus read margins
SIM_CHUNK = 200_000             # pairs simulated at once: bounds memory (a pairs x 150-byte read matrix) at any -n
SIM_ERR = 0.005                 # substitution rate of simulated reads (0.5%, a typical Illumina error rate)
SIM_TRUTH_TOL = 1_000           # a simulated mate counts as correctly placed within 1 kb of its true position
SIM_FASTQ_READS = 200_000       # FASTQ reads sampled for the read-length distribution (ample for an empirical one)
