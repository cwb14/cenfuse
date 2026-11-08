"""Every constant has a reason, and docs/constants.md lists them as the code states them (future.md #64, the rule 'give
every constant a reason')."""
from pathlib import Path

from cenfuse import consts

DOC = Path(__file__).resolve().parents[2] / "docs" / "constants.md"


def test_every_numeric_constant_has_a_reason():
    bare = [f"{m}:{line} {name}" for m, line, name, _, reason in consts.scan() if not reason]
    assert not bare, "constants without a reason comment: " + ", ".join(bare)


def test_every_satellite_parameter_has_a_reason():
    assert all(reason for *_, reason in consts.params_rows())


def test_the_constants_doc_is_current():
    want = consts.markdown(consts.scan() + consts.params_rows())
    assert DOC.read_text() == want, "docs/constants.md is stale: python -m cenfuse.consts docs/constants.md"


def test_reason_reading(tmp_path):
    pkg = tmp_path / "pk"
    pkg.mkdir()
    (pkg / "m.py").write_text('# reason above\nA = 1\nB = (1, 2)  # same line\n'
                              'C = 3       # first part\n            # second part\nD = "x"  # not numeric\nlow = 4\n')
    rows = {name: reason for _, _, name, _, reason in consts.scan(pkg)}
    assert rows == {"A": "reason above", "B": "same line", "C": "first part second part"}
    assert consts.mark("a named heuristic", "x.py") == "heuristic" and consts.mark("prior: 20 kb", "x.py") == "prior"
