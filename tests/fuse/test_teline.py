"""teline.py (future.md #65): the LTR line read for the loci table, and what it says about an unusable tes/ file."""
import logging

import numpy as np

from cenfuse import teline


def _gff(path, rows):
    path.write_text("##gff-version 3\n" + "".join(f"{c}\tLTRquest\tLTR_retrotransposon\t{s}\t{e}\t.\t+\t.\tID=te{i};clade={cl}\n"
                                                  for i, (c, s, e, cl) in enumerate(rows)))
    return path


def test_no_file_means_no_line():
    d = teline.te_line(None, {"c1": 1000})
    assert d["marker"] is None and d["clades"] == [] and d["signature"] is None


def test_a_file_on_other_sequence_names_is_said(tmp_path, caplog):
    gff = _gff(tmp_path / "g.gff3", [("other", 100, 5000, "CRM")] * 3)
    lg = logging.getLogger("cenfuse")
    keep = lg.propagate
    lg.propagate = True
    try:
        with caplog.at_level(logging.WARNING, logger="cenfuse"):
            d = teline.te_line(gff, {"c1": 10_000_000})
    finally:
        lg.propagate = keep
    assert d["marker"] is None and "other sequence names" in caplog.text


def test_composition_counts_midpoints_in_the_interval():
    mids = {"CRM": {"c1": np.array([10, 20, 30])}, "Tekay": {"c1": np.array([25])}, "Athila": {"c2": np.array([5])}}
    assert teline.composition(mids, "c1", 15, 31) == (3, "CRM:2,Tekay:1")
    assert teline.composition(mids, "c1", 100, 200) == (0, ".")
