"""Fixed internal constants, one place, one reason each (rationale in docs/design.md)."""
from __future__ import annotations

MIN_CALLS_SITE = 3          # min N_valid per cytosine: per-site usability floor for ~4x samples
BIN = 50_000                # bin (bp): >= 20 bins per CENH3 domain; validated on P. annua
DEPTH_LO = 0.5              # drop a bin below half the sample's median depth: one-copy loss, missing DNA (heuristic)
DEPTH_HI = 2.0              # drop a bin above double the median depth: duplication, collapsed repeat (heuristic)
WIN_D = 0.4                 # sliding window = 0.4 D: narrower than the domain (heuristic)
RING_LO_D = 1.5             # class reference ring starts 1.5 D away: beyond the domain (heuristic)
RING_HI_D = 4.0             # ...and ends 4 D away: wide enough for stable class levels (heuristic)
MIN_LEN_D = 10              # sequences >= 10 D are analysed: a 4-D ring each side of a 2-D centre
GAP_FOLD = 10               # chromosomes vs unplaced scaffolds: the first tenfold drop in sorted sequence length after
                            # half the assembly (heuristic)
D_REF_BP = 22.2e6           # summed CENH3 domain length of maize Mo17 (2.22 Mb x 10 centromeres; Chen et al. 2023)
D_REF_G = 2.18e9            # maize Mo17 genome size (bp)
D_SLOPE = 0.97              # CENH3 area scales with genome size^0.97 (Plackova et al. 2021 Sci Rep)
MIN_ANCHORS = 3             # calibration needs >= 3 anchored chromosomes with data
N_PERM = 10_000             # random placements: P resolution 1e-4, so Benjamini-Hochberg over up to ~500 class x
                            # context tests can still pass a lone true class
FDR = 0.05                  # Benjamini-Hochberg level (conventional)
MIN_CALLS = 500             # calls needed for a class (or the track) in a window: deviation SE ~0.02
MIN_EFFECT = 0.04           # a signature entry's mean deviation at the anchors must be >= 2 x the ~0.02 SE of one
                            # candidate window at MIN_CALLS calls, so one window can tell it from zero
SEP_D = 4                   # candidates, random placements, 'near an anchor': 4 D apart, beyond the ring's inner edge
MAX_CAND = 5                # candidates per sample x chromosome: enough to hold the true centromere without flooding
                            # (heuristic)
SIG_FRACTION = 0.5          # nearest of two references, used twice: a signature class is concentrated at the anchors
                            # (genome share <= 0.5 x anchor share), and a centromere_like candidate is closer to the
                            # anchors' signature (score 1) than to the background (score 0)
SPECIFIC_P = 0.05           # plant-specific flag (conventional level; per chromosome, not corrected for multiple tests)
LOW_DEPTH = 8.0             # warn below 8x median depth: few repeat cytosines reach MIN_CALLS_SITE
CHUNK = 10_000_000          # bedMethyl region read at once (bp): bounds memory. Each worker holds one chunk of rows
                            # (~0.5 GB CHG, ~1 GB CHH at plant densities), so a run needs about threads x 1-2 GB
SEED = 7                    # fixed seed: reruns give identical P values
CACHE_FP_BYTES = 1 << 20    # cache key reads the first and last MiB of each input (a fast fingerprint, not a checksum)
