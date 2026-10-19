from __future__ import annotations

import h5py
import numpy as np
import pytest

from cenfuse.hic.io import InputError, cool_resolutions, read_cool

EXPECTED = np.zeros((5, 5))
EXPECTED[2, 2] = 5                                  # a0-a0
EXPECTED[2, 0] = EXPECTED[0, 2] = 2                 # a0-b0
EXPECTED[4, 1] = EXPECTED[1, 4] = 7                 # a2-b1


def _fill(g, res):
    g["chroms/name"] = np.array([b"a", b"b", b"z"])
    g["chroms/length"] = np.array([300, 200, 100])
    g["bins/chrom"] = np.array([0, 0, 0, 1, 1, 2])
    g["bins/start"] = np.array([0, 100, 200, 0, 100, 0])
    g["bins/end"] = np.array([100, 200, 300, 100, 200, 100])
    g["pixels/bin1_id"] = np.array([0, 0, 2, 4])
    g["pixels/bin2_id"] = np.array([0, 3, 4, 5])
    g["pixels/count"] = np.array([5, 2, 7, 9])


def make_mcool(path, res=100):
    with h5py.File(path, "w") as f:
        _fill(f.create_group(f"resolutions/{res}"), res)


def make_cool(path, res=100):
    with h5py.File(path, "w") as f:
        f.attrs["bin-size"] = res
        _fill(f, res)


def test_read_mcool_maps_to_fai_order(tmp_path):
    p = tmp_path / "x.mcool"
    make_mcool(p)
    assert cool_resolutions(p) == [100]
    M = read_cool(p, 100, {"b": 200, "a": 300}).toarray()    # .fai order: b first; 'z' not in the .fai -> dropped
    # our bins: b0 b1 a0 a1 a2
    assert M[2, 2] == 5                     # a0-a0
    assert M[2, 0] == 2 and M[0, 2] == 2    # a0-b0 (cool bins 0,3)
    assert M[4, 1] == 7 and M[1, 4] == 7    # a2-b1 (cool bins 2,4)
    assert M.sum() == 5 + 2 * 2 + 2 * 7     # the b1-z pixel is dropped


def test_absent_resolution_in_mcool_lists_available(tmp_path):
    p = tmp_path / "x.mcool"
    make_mcool(p, res=100)
    with pytest.raises(InputError, match=r"resolution 200.*available: \[100\]"):
        read_cool(p, 200, {"b": 200, "a": 300})


def test_single_resolution_cool_reads_correct_matrix(tmp_path):
    p = tmp_path / "x.cool"
    make_cool(p, res=100)
    assert cool_resolutions(p) == [100]
    M = read_cool(p, 100, {"b": 200, "a": 300}).toarray()
    assert np.array_equal(M, EXPECTED)


def test_single_resolution_cool_mismatch_raises(tmp_path):
    p = tmp_path / "x.cool"
    make_cool(p, res=100)
    with pytest.raises(InputError, match="bin-size 100, requested 50"):
        read_cool(p, 50, {"b": 200, "a": 300})


def test_non_hdf5_file_raises(tmp_path):
    p = tmp_path / "x.cool"
    p.write_text("chr1\t300\n")
    with pytest.raises(InputError, match="cannot open"):
        cool_resolutions(p)
    with pytest.raises(InputError, match="cannot open"):
        read_cool(p, 100, {"a": 300})


def test_hdf5_without_cool_groups_raises(tmp_path):
    p = tmp_path / "other.h5"
    with h5py.File(p, "w") as f:
        f.attrs["bin-size"] = 100
        f["something"] = np.arange(3)
    with pytest.raises(InputError, match="not a .cool/.mcool"):
        read_cool(p, 100, {"b": 200, "a": 300})
