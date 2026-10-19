"""Fixed internal constants, one place, one reason each (rationale in docs/design.md)."""
from __future__ import annotations

MAPQ = 30                       # conventional unique-mapping threshold (applied when the pairs carry MAPQ columns)
BIN_LADDER = (250_000, 500_000, 1_000_000, 2_000_000)  # candidate bin sizes, finest first (bp)
MIN_EXPECTED_CEN = 25           # expected centromere contacts per bin: a 2x rise is then >= 5 SD over Poisson noise
HALF_WIN = 4                    # partner window = estimate +/- 4 bins: a 1-2 Mb CENH3 domain plus one bin of error
END_WIN = 8                     # telomere-pole window at each chromosome end (bins), about one partner window
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
