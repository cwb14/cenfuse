"""Helpers shared by the sat tests."""
import shutil
import numpy as np
TEST_PARAMS = dict(
    seed_chunk=150_000, seed_overlap=5_000, min_fam_span=5_000, hmm_sample=100, hmm_sample_bp=100_000,
    chrom_min_bp=5_000, cluster_gap=50_000, strong_cluster_bp=20_000,
    low_sat_bp=10_000, fallback_min_bp=2_000, lu_flank=20_000, lu_window=100_000, lu_window_overlap=20_000,
)
def have(*tools) -> bool:
    return all(shutil.which(t) for t in tools)


ONE_ARRAY = "M 0 100 300\nL 3 50\nX 0\n"   # ONEview text: one array of 50-bp units at [100, 300) of the piece


def dna(rng, n: int) -> str:
    """n random bases; fast enough for megabases."""
    return np.frombuffer(b"ACGT", dtype=np.uint8)[rng.integers(0, 4, n)].tobytes().decode()


def crashing_fastan(monkeypatch, seqs: dict, dies):
    """Replace fastan.run_fastan by a FasTAN stand-in (no FasTAN needed) and return the list of pieces it finished.
    A piece is located in seqs ({name: sequence}, random enough that a piece occurs once) by its content.
    dies(name, start, end) says whether the real FasTAN would abort on that piece: it then dies by SIGABRT as FasTAN
    0.8 does (util.run reports 'exit -6'); otherwise it finds ONE_ARRAY and the piece is recorded as (name, start,
    end). Like FAtoGDB it refuses a sequence without bases."""
    from cenfuse.sat import fastan
    from cenfuse.sat.util import SatelliteError
    done = []

    def run(piece, tmpdir):
        assert any(b in piece for b in "ACGT"), "FAtoGDB rejects a sequence without bases"
        at = [(name, s.find(piece)) for name, s in seqs.items() if piece in s]
        assert len(at) == 1, "the piece must occur once in the test sequences"
        (name, start), end = at[0], at[0][1] + len(piece)
        if dies(name, start, end):
            raise SatelliteError("command failed (exit -6): FasTAN -m -p -T1 s.1gdb\nmalloc(): corrupted top size")
        done.append((name, start, end))
        return ONE_ARRAY

    monkeypatch.setattr(fastan, "run_fastan", run)
    return done
