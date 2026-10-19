from __future__ import annotations

import numpy as np
import pandas as pd
import scipy.sparse as sp

from cenfuse.hic.bins import Bins
from cenfuse.hic.io import read_intervals
from cenfuse.hic.qc import call_flags, hub_bins, inside, merge_intervals, rdna_intervals, read_class_split, satellite_centres


def test_merge_and_inside():
    df = pd.DataFrame(dict(chrom=["a", "a", "a"], start=[10, 15, 40], end=[20, 30, 50]))
    iv = merge_intervals(df)
    assert list(iv["a"][0]) == [10, 40] and list(iv["a"][1]) == [30, 50]
    assert list(inside(iv, "a", np.array([9, 10, 29, 30, 45]))) == [False, True, True, False, True]
    assert not inside(iv, "b", np.array([10])).any()


def test_hub_bins():
    b = Bins({"a": 1000, "b": 1000}, 100)
    A = np.ones((b.N, b.N))
    A[3, 10:] = A[10:, 3] = 6.0          # bin 3: trans x6
    hub = hub_bins(sp.csr_matrix(A), b, np.ones(b.N, bool))
    assert hub[3] and hub.sum() == 1


def test_satellite_centres_and_rdna(tmp_path):
    p = tmp_path / "s.gff"
    p.write_text("##gff-version 3\nchr1\tc\tcentromere\t1001\t2000\t.\t.\t.\tID=x\n"
                 "chr1\tc\tsatellite_DNA\t1001\t1100\t.\t+\t.\tID=y;family=F;rdna=45S_spacer\n"
                 "chr1\tc\tsatellite_DNA\t5001\t5100\t.\t+\t.\tID=z;family=G\n")
    df = read_intervals(p)
    b = Bins({"chr1": 10_000}, 1000)
    assert satellite_centres(df, b) == {0: 1500.0}
    iv = rdna_intervals(df)
    assert list(iv["chr1"][0]) == [1000]


def test_flags():
    b = Bins({"a": 100 * 100}, 100)
    pc = dict(peak=3, core_lo=2, core_hi=4)
    hub = np.zeros(b.N, bool)
    hub[3] = True
    iv = {"a": (np.array([250]), np.array([260]))}
    f = call_flags(0, pc, b, hub, iv, rabl_passed=False, p_null=float("nan"), converged=False, cycle_unstable=True)
    assert set(f) == {"hub", "rDNA", "telomeric", "no_rabl", "null_untestable", "not_converged", "cycle_unstable"}
    pc = dict(peak=50, core_lo=49, core_hi=51)
    assert call_flags(0, pc, b, np.zeros(b.N, bool), {}, rabl_passed=True, p_null=0.01, converged=True,
                      cycle_unstable=False) == []


def test_read_class_split_counts(tmp_path):
    b = Bins({"a": 2000, "b": 2000}, 100)
    pairs = tmp_path / "p.txt"
    pairs.write_text("a\t1050\t+\tb\t550\t-\n"     # core end (non-repeat) -> b window: centromere contact
                     "a\t1010\t+\tb\t1950\t-\n"    # core end inside a repeat -> b outside window
                     "a\t150\t+\tb\t560\t-\n")     # background end -> b window
    rep = {"a": (np.array([1000]), np.array([1020]))}
    t = read_class_split(pairs, {"a": 2000, "b": 2000}, b, np.ones(2, bool),
                         {0: dict(peak=10, core_lo=10, core_hi=10)}, {0: 10, 1: 5}, rep)
    r = t.set_index(["zone", "cls"])
    assert r.loc[("core", "non_repeat"), "cen"] == 1 and r.loc[("core", "repeat"), "trans"] == 1
    assert r.loc[("core", "repeat"), "cen"] == 0 and r.loc[("background", "non_repeat"), "cen"] == 1
