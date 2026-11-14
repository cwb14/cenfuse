import os
from pathlib import Path

import pandas as pd
import pysam
import pytest

from cenfuse.meth import bins, classes, io, params


def _bm(path, rows):
    with open(path, "w") as f:
        for c, p, nv, nm in rows:
            f.write(f"{c}\t{p}\t{p + 1}\tm\t{nv}\t+\t{p}\t{p + 1}\t255,0,0\t{nv}\t0\t{nm}\t{nv - nm}\t0\t0\t0\t0\t0\n")
    pysam.tabix_index(str(path), preset="bed", force=True)
    return f"{path}.gz"


def test_bin_chromosome_counts_by_bin_and_class(tmp_path):
    B = params.BIN
    f = _bm(tmp_path / "a.bed", [("c1", 10, 4, 1), ("c1", 20, 5, 5), ("c1", B + 5, 6, 3), ("c1", B + 7, 2, 2)])
    cm = classes.build([pd.DataFrame([("c1", 0, 15, "gene", 0)], columns=["chrom", "start", "end", "cls", "depth"])],
                       ["c1"])
    cb, bad = bins.bin_chromosome(f, "c1", 3 * B, cm.segs["c1"], cm.other, len(cm.names))
    n, v, m = cb.dense(len(cm.names), 3)
    assert bad == 0
    assert n[0].tolist() == [1, 0, 0] and v[0].tolist() == [4, 0, 0] and m[0].tolist() == [1, 0, 0]   # gene
    assert n[1].tolist() == [1, 1, 0] and v[1].tolist() == [5, 6, 0] and m[1].tolist() == [5, 3, 0]   # other


def _samples(tmp_path, rows_by_sample):
    out = []
    for s, rows in rows_by_sample.items():
        out.append(io.Sample(s, {"CHG": _bm(tmp_path / f"{s}.bed", rows)}))
    return out


def test_collect_uses_cache_and_reruns_on_new_annotation(tmp_path, monkeypatch):
    smp = _samples(tmp_path, {"S1": [("c1", 10, 4, 1)]})
    cm = classes.build([pd.DataFrame([("c1", 0, 15, "gene", 0)], columns=["chrom", "start", "end", "cls", "depth"])],
                       ["c1"])
    a = bins.collect(smp, {"c1": 1000}, cm, "k1", tmp_path / "cache", 1)

    def boom(*args, **kw):
        raise AssertionError("re-binned although cached")
    monkeypatch.setattr(bins, "bin_chromosome", boom)
    b = bins.collect(smp, {"c1": 1000}, cm, "k1", tmp_path / "cache", 1)
    assert b[("S1", "CHG")]["c1"].v.tolist() == a[("S1", "CHG")]["c1"].v.tolist() == [4]
    with pytest.raises(AssertionError, match="re-binned"):
        bins.collect(smp, {"c1": 1000}, cm, "k2", tmp_path / "cache", 1)


def test_depth_filter_drops_outlier_bins_and_sample_exclusions(tmp_path):
    B = params.BIN
    rows = [("c1", i * B + 7, 10, 5) for i in range(20)]
    rows[3] = ("c1", 3 * B + 7, 40, 5)                     # 4x median depth
    smp = _samples(tmp_path, {"S1": rows, "S2": [("c1", i * B + 7, 10, 5) for i in range(20)]})
    cm = classes.build([pd.DataFrame([("c1", 0, 1, "gene", 0)], columns=["chrom", "start", "end", "cls", "depth"])],
                       ["c1"])
    raw = bins.collect(smp, {"c1": 20 * B}, cm, "k", None, 1)
    excl = pd.DataFrame(dict(chrom=["c1"], start=[10 * B], end=[12 * B], sample=["S2"]))
    filt, st = bins.depth_filter(raw, smp, {"c1": 20 * B}, excl)
    s1, s2 = st
    assert (s1.dropped_depth, s1.dropped_excluded, s2.dropped_depth, s2.dropped_excluded) == (1, 0, 0, 2)
    assert 3 not in filt[("S1", "CHG")]["c1"].bin.tolist()
    assert {10, 11} & set(filt[("S2", "CHG")]["c1"].bin.tolist()) == set()
    assert {10, 11} <= set(filt[("S1", "CHG")]["c1"].bin.tolist())


def _one_class(chroms):
    genes = pd.DataFrame([(chroms[0], 0, 1, "gene", 0)], columns=["chrom", "start", "end", "cls", "depth"])
    return classes.build([genes], chroms)


def test_fingerprint_sees_a_change_past_the_first_mib(tmp_path):
    """A 1.5-MiB file is read whole: a same-size edit at 1.2 MiB with the old mtime changes the key."""
    f = tmp_path / "a.gff3"
    f.write_bytes(b"x" * 1_500_000)
    st = f.stat()
    k = bins.fingerprint([f])
    with open(f, "r+b") as h:
        h.seek(1_200_000)
        h.write(b"y")
    os.utime(f, ns=(st.st_atime_ns, st.st_mtime_ns))
    assert bins.fingerprint([f]) != k


def test_cache_file_of_a_sequence_name_with_a_slash(tmp_path):
    smp = _samples(tmp_path, {"S1": [("chr1/hap1", 10, 4, 1)]})
    raw = bins.collect(smp, {"chr1/hap1": 1000}, _one_class(["chr1/hap1"]), "k", tmp_path / "cache", 1)
    assert raw[("S1", "CHG")]["chr1/hap1"].v.tolist() == [4]
    assert len(list((tmp_path / "cache").iterdir())) == 1


def test_cache_is_rebinned_when_a_chromosome_length_changes(tmp_path, monkeypatch):
    smp = _samples(tmp_path, {"S1": [("c1", 10, 4, 1)]})
    cm = _one_class(["c1"])
    bins.collect(smp, {"c1": 1000}, cm, "k", tmp_path / "cache", 1)
    lengths, real = [], bins.bin_chromosome
    monkeypatch.setattr(bins, "bin_chromosome", lambda *a: lengths.append(a[2]) or real(*a))
    raw = bins.collect(smp, {"c1": 3 * params.BIN}, cm, "k", tmp_path / "cache", 1)
    assert lengths == [3 * params.BIN] and raw[("S1", "CHG")]["c1"].v.tolist() == [4]


def test_failed_cache_write_leaves_no_temporary_file(tmp_path, monkeypatch):
    smp = _samples(tmp_path, {"S1": [("c1", 10, 4, 1)]})

    def disk_full(f, **kw):
        Path(f).write_bytes(b"part")
        raise OSError(28, "No space left on device")
    monkeypatch.setattr(bins.np, "savez", disk_full)
    with pytest.raises(OSError, match="No space"):
        bins.collect(smp, {"c1": 1000}, _one_class(["c1"]), "k", tmp_path / "cache", 1)
    assert list((tmp_path / "cache").iterdir()) == []


def test_exclusion_with_a_negative_start(tmp_path):
    """BED [-5, 60000) excludes bins 0 and 1, not the slice from the chromosome's last bin."""
    B = params.BIN
    smp = _samples(tmp_path, {"S1": [("c1", i * B + 7, 10, 5) for i in range(12)]})
    raw = bins.collect(smp, {"c1": 12 * B}, _one_class(["c1"]), "k", None, 1)
    excl = pd.DataFrame(dict(chrom=["c1", "c1"], start=[-5, -3 * B], end=[B + 10_000, -B], sample=["", ""]))
    filt, st = bins.depth_filter(raw, smp, {"c1": 12 * B}, excl)
    assert st[0].dropped_excluded == 2
    assert sorted(filt[("S1", "CHG")]["c1"].bin.tolist()) == list(range(2, 12))


def test_bin_chromosome_restores_htslib_verbosity(tmp_path):
    f = _bm(tmp_path / "a.bed", [("c1", 10, 4, 1)])
    old = pysam.set_verbosity(3)
    try:
        bins.bin_chromosome(f, "c1", 1000, None, 0, 1)
        assert pysam.set_verbosity(3) == 3
    finally:
        pysam.set_verbosity(old)
