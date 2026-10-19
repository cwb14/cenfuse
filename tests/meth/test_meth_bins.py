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
