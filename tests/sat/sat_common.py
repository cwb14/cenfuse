"""Helpers shared by the sat tests."""
import shutil
import pytest
from cenfuse.sat.params import Params
from sat_synth import build_complex
TEST_PARAMS = dict(
    seed_chunk=150_000, seed_overlap=5_000, min_fam_span=5_000, hmm_sample=100, hmm_sample_bp=100_000,
    chrom_min_bp=5_000, cluster_gap=50_000, strong_cluster_bp=20_000,
    low_sat_bp=10_000, fallback_min_bp=2_000, lu_flank=20_000, lu_window=100_000, lu_window_overlap=20_000,
)
def have(*tools) -> bool:
    return all(shutil.which(t) for t in tools)
