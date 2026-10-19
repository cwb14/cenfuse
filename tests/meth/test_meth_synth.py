import numpy as np
import pysam
import meth_synth as synth

from cenfuse.meth import classes, io


def test_files_and_indexes(syn):
    for k in ("fai", "samples", "ltrquest", "satellites", "genes", "anchors", "anchors_wrong", "exclude"):
        assert syn[k].exists(), k
    s = io.read_samples(syn["samples"])
    assert [x.name for x in s] == ["S1", "S2", "S3"] and list(s[0].files) == ["CHG", "CG"]
    assert all(f.with_name(f.name + ".tbi").exists() for x in s for f in x.files.values())


def test_classes_match_truth(syn):
    cm = classes.build([io.read_satellites(syn["satellites"]), io.read_ltrquest(syn["ltrquest"]),
                        io.read_genes(syn["genes"])], synth.CHROMS)
    pos = np.arange(15_000_000, 25_000_000, 97)
    truth = synth.site_classes(syn["ann"], "chr1", pos)
    got = np.array(cm.names)[cm.classify("chr1", pos)]
    want = np.array([{"satA": "sat:satA", "CRM": "LTR:CRM", "Tekay": "LTR:Tekay"}.get(x, x) for x in truth])
    ok = (got == want) | ((got == "LTR:Gypsy_unclassified") & np.isin(want, ["LTR:CRM", "LTR:Tekay"]))
    assert ok.all()


def test_planted_crm_loss_in_raw_sites(syn):
    with pysam.TabixFile(str(syn["dir"] / "S1.CHG.bedmethyl.gz")) as tb:
        pos, v, m, _ = io.read_bedmethyl(tb, "chr1", 15_000_000, 25_000_000, 3)
    k = synth.site_classes(syn["ann"], "chr1", pos)
    inside = np.abs(pos - syn["cen"]["chr1"]) <= synth.DOMAIN / 2
    crm_in = m[(k == "CRM") & inside].sum() / v[(k == "CRM") & inside].sum()
    crm_out = m[(k == "CRM") & ~inside].sum() / v[(k == "CRM") & ~inside].sum()
    sat = m[k == "satA"].sum() / v[k == "satA"].sum()
    assert abs(crm_in - 0.20) < 0.03 and abs(crm_out - 0.45) < 0.02 and abs(sat - 0.68) < 0.02


def test_depth_filter_on_synthetic(prep):
    st = {d.sample: d for d in prep["stats"]}
    assert st["S1"].dropped_depth == 2 and st["S2"].dropped_excluded == 40 and st["S3"].dropped_depth == 0
