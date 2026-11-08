"""benchmark/precision.py: widths and midpoints of right calls against their domain."""
import pandas as pd

from benchmark import precision


def test_widths_midpoints_and_the_classes():
    rows = [dict(arm="S", kind="CENP-A", start="1000001", end="5000000", truth="2000000-2400000", confidence="high",
                 outcome="right_high"),
            dict(arm="NEG", kind="CENP-A", start="1", end="10", truth="2000000-2400000", confidence="high", outcome="right_high"),
            dict(arm="S", kind="agreement", start="1", end="10", truth="1-20", confidence="high", outcome="right_high"),
            dict(arm="S", kind="CENH3", start="100001", end="200000", truth="300000-400000", confidence="low", outcome="right_low"),
            dict(arm="S", kind="CENH3", start="1", end="10", truth="1-20", confidence="low", outcome="wrong_low")]
    d = precision.calls(pd.DataFrame(rows))
    assert len(d) == 2                                   # controls, agreement and wrong calls left out
    a = d[d.kind == "CENP-A"].iloc[0]
    assert a.call_bp == 4_000_000 and a.truth_bp == 400_000 and a.mid_dist == 800_000 and not a.mid_inside
    assert abs(a.on_domain - 0.1) < 1e-9
    b = d[d.kind == "CENH3"].iloc[0]
    assert b.on_domain == 0 and not b.confident          # right within the tolerance, no overlap
    t = precision.table(d)
    assert "animals" in t and "plants" in t and t.count("\n") == 2 + 3   # header, rule, animals x2, plants x1
