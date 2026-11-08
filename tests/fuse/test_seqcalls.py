"""seqcalls.py (future.md #65): the satellite engine's outputs read back, 1-based inclusive to 0-based half-open."""
import pytest

from cenfuse import seqcalls
from cenfuse.util import CenfuseError

HEAD = ("chrom\tchrom_len\tcen_start\tcen_end\tcen_len\tevidence\tconfidence\tflags\tmain_family\tmain_family_bp\tcen_sat_bp\t"
        "main_spacing\tte_clade\tte_core_start\tte_core_end\n")


def _write(d, summary_rows, arrays):
    (d / "g_centromere_summary.tsv").write_text(HEAD + "".join(summary_rows))
    (d / "g_satellite_arrays.gff").write_text("##gff-version 3\n" + "".join(arrays))


def test_calls_and_arrays_become_half_open(tmp_path):
    _write(tmp_path, ["c1\t1000\t101\t200\t100\tsatellite\thigh\t.\tfamA\t90\t90\t180\tCRM\t151\t300\n",
                      "c2\t1000\t.\t.\t.\tnone\tnone\t.\t.\t.\t.\t.\t.\t.\t.\n"],
           ["c1\tcenfuse\tsatellite_DNA\t101\t200\t.\t+\t.\tID=a;family=famA;role=centromeric;monomer_len=180;"
            "n_monomers=5;in_centromere=true\n",
            "c9\tcenfuse\tsatellite_DNA\t1\t50\t.\t+\t.\tfamily=famA\n"])                   # not a chromosome: left out
    s = seqcalls.read(tmp_path, "g", {"c1": 1000, "c2": 1000, "c3": 1000})
    c1 = s.call("c1")
    assert (c1.start, c1.end, c1.te_start, c1.te_end, c1.sat_bp, c1.spacing) == (100, 200, 150, 300, 90, 180)
    assert s.call("c2") is None and s.call("c3") is None                     # no call; no row at all
    assert list(s.arrays.chrom) == ["c1"] and s.arrays.iloc[0].start == 100 and bool(s.arrays.iloc[0].in_centromere)
    assert seqcalls.bp_by_role(s.arrays, "c1", 150, 400) == {"centromeric": 50}
    assert s.test is None and s.families.empty


def test_a_file_that_is_no_summary_or_is_missing_stops(tmp_path):
    with pytest.raises(CenfuseError, match="missing satellite-engine output"):
        seqcalls.read(tmp_path, "g", {"c1": 10})
    (tmp_path / "g_centromere_summary.tsv").write_text("a\tb\n1\t2\n")
    (tmp_path / "g_satellite_arrays.gff").write_text("")
    with pytest.raises(CenfuseError, match="not a satellite-engine centromere summary"):
        seqcalls.read(tmp_path, "g", {"c1": 10})


def test_gff_attributes_are_unescaped():
    assert seqcalls.gff_attrs("ID=x;family=a%3Bb; role = centromeric") == {"ID": "x", "family": "a;b", "role": "centromeric"}
