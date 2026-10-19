"""Helpers shared by the meth tests (moved out of conftest.py by tools/vendor.py)."""
import sys
from pathlib import Path
import pytest
sys.path.insert(0, str(Path(__file__).parent))
import meth_synth as synth  # noqa: E402
def run_scan(prep, anchor_file, direction):
    """Anchors (or none) -> (anchors, scan, calibration table or None)."""
    from cenfuse.meth import calibrate, io
    from cenfuse.meth.scan import scan
    anchors = calibrate.anchor_table(io.read_anchors(anchor_file)[0], prep["chroms"]) if anchor_file else {}
    abin = {c: a["bin"] for c, a in anchors.items()}
    s = scan(prep["bins"], prep["samples"], prep["chroms"], len(prep["cmap"].names), prep["geo"], abin, direction)
    cal = calibrate.calibrate(s.pooled, s.totals, prep["cmap"].names) if anchors else None
    return anchors, s, cal
