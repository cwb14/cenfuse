import numpy as np
import pandas as pd

from cenfuse.meth import classes


def _f(rows):
    return pd.DataFrame(rows, columns=["chrom", "start", "end", "cls", "depth"])


def test_priority_satellite_then_transposon_then_gene():
    cm = classes.build([_f([("c1", 0, 100, "gene", 0)]), _f([("c1", 20, 80, "LTR:Tekay", 0)]),
                        _f([("c1", 40, 60, "sat:f1", 0)])], ["c1"])
    assert cm.names == ["sat:f1", "LTR:Tekay", "gene", "other"]
    got = [cm.names[k] for k in cm.classify("c1", np.array([5, 25, 50, 75, 99, 100, 500]))]
    assert got == ["gene", "LTR:Tekay", "sat:f1", "LTR:Tekay", "gene", "other", "other"]


def test_innermost_element_wins_by_nesting_depth():
    te = _f([("c1", 0, 1000, "LTR:Tekay", 0), ("c1", 300, 600, "LTR:CRM", 1), ("c1", 400, 450, "LTR:SIRE", 2)])
    cm = classes.build([te], ["c1"])
    got = [cm.names[k] for k in cm.classify("c1", np.array([100, 350, 420, 500, 700]))]
    assert got == ["LTR:Tekay", "LTR:CRM", "LTR:SIRE", "LTR:CRM", "LTR:Tekay"]


def test_same_depth_shorter_element_wins():
    te = _f([("c1", 0, 1000, "LTR:Tekay", 0), ("c1", 100, 200, "LTR:CRM", 0)])
    cm = classes.build([te], ["c1"])
    assert [cm.names[k] for k in cm.classify("c1", np.array([150, 250]))] == ["LTR:CRM", "LTR:Tekay"]


def test_paint_merges_adjacent_same_label_and_keeps_gaps():
    s, e, lab = classes.paint(np.array([0, 10, 30]), np.array([10, 20, 40]), np.array([0, 1, 2]),
                              np.array([0, 0, 0]))
    assert s.tolist() == [0, 30] and e.tolist() == [20, 40] and lab.tolist() == [0, 0]


def test_unannotated_chromosome_is_other():
    cm = classes.build([_f([("c1", 0, 10, "gene", 0)])], ["c1", "c2"])
    assert cm.classify("c2", np.array([5])).tolist() == [cm.other]
    assert cm.classify("c1", np.zeros(0, int)).size == 0


def test_all_empty_features_gives_only_other():
    cm = classes.build([_f([])], ["c1"])
    assert cm.names == ["other"] and cm.segs == {}
    assert cm.classify("c1", np.array([5])).tolist() == [cm.other]


def test_features_outside_chroms_are_dropped():
    cm = classes.build([_f([("c2", 0, 10, "gene", 0)])], ["c1"])
    assert "c2" not in cm.segs
    assert cm.names == ["other"]
