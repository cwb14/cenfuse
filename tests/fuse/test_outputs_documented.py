"""Every column of every table cenfuse writes at the top of its output folder is documented (future.md #58).

A column counts as documented when its name is in backticks in docs/outputs.md or docs/plants.md, alone or as an item
of a comma list (`chrom`, `chrom_len`), or matches a backticked wildcard (`meth_*`). The tables come from one run on
the synthetic genome of fuse_synth.py (Hi-C from two plants, tes/, genes/, CHG methylation of three plants), whose LTR
retrotransposons are given LTRquest's family= so that the LTR signature table is written too."""
import fnmatch
import re
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent))
import fuse_synth  # noqa: E402

from cenfuse.pipeline import Run  # noqa: E402

DOCS = Path(__file__).parents[2] / "docs"
DOC_FILES = ("outputs.md", "plants.md")
GENOME = "Syn"
# every table the run must write: one it stopped writing would escape the column check
TABLES = ("centromeres", "loci", "feature_signatures", "feature_contacts", "line_tests", "te_signature", "provenance",
          "plants", "plant_sites")
# tables whose own section names every column: the rule above also accepts a name documented for another table
# (the loci's `rdna` passes on the --force stage `rdna`)
OWN_SECTION = {"loci": "outputs.md", "feature_signatures": "outputs.md", "feature_contacts": "outputs.md",
               "line_tests": "outputs.md", "te_signature": "outputs.md", "provenance": "plants.md", "plants": "plants.md",
               "plant_sites": "plants.md"}
# a code span, as in CommonMark: it may hold a line break (outputs.md wraps one), not a blank line
CODE_SPAN = re.compile(r"`((?:[^`\n]|\n(?![ \t]*\n))+)`")


def documented(text: str) -> tuple:
    """(names, wildcards) of the backticked items of text: each code span split at commas; an item holding * is a
    wildcard."""
    names, wild = set(), set()
    for span in CODE_SPAN.findall(text):
        for item in span.replace("\n", " ").split(","):
            item = item.strip()
            if item:
                (wild if "*" in item else names).add(item)
    return names, wild


def undocumented(columns: list, names: set, wild: set) -> list:
    return [c for c in columns if c not in names and not any(fnmatch.fnmatchcase(c, w) for w in wild)]


def header(path: Path) -> list:
    with open(path) as f:
        for line in f:
            if not line.startswith("#"):
                return line.rstrip("\r\n").split("\t")
    return []


def section(text: str, table: str) -> str:
    """The section whose heading names `asm.<table>.tsv`, to the next heading of its level or above; '' if none."""
    lines, fence, start, level = text.splitlines(), False, None, 0
    for i, x in enumerate(lines):
        if x.lstrip().startswith("```"):
            fence = not fence
        m = None if fence else re.match(r"(#+) ", x)
        if m is None:
            continue
        if start is not None and len(m.group(1)) <= level:
            return "\n".join(lines[start:i])
        if start is None and f"`asm.{table}.tsv`" in x:
            start, level = i, len(m.group(1))
    return "\n".join(lines[start:]) if start is not None else ""


def with_families(gff: Path) -> None:
    """Name each element's family (LTRquest's family=; fuse_synth gives the lineage only), one family per lineage, so
    the LTR signature is built and its table written."""
    rows = []
    for x in gff.read_text().splitlines():
        if not x.startswith("#"):
            x += f";family={re.search(r'clade=([^;]+)', x).group(1)}_fam1"
        rows.append(x)
    gff.write_text("\n".join(rows) + "\n")


@pytest.fixture(scope="module")
def out(tmp_path_factory):
    root = tmp_path_factory.mktemp("documented")
    truth = fuse_synth.make(root)
    with_families(truth["dir"] / "tes" / f"{GENOME}.gff3")
    o = root / "out"
    Run(truth["dir"], o, 4, satellite_dir=truth["sat"]).run()
    return o


def test_the_run_writes_every_table(out):
    found = {p.name for p in out.glob("*.tsv")}
    missing = [f for f in [f"{GENOME}.{t}.tsv" for t in TABLES] + ["cenfuse_summary.tsv"] if f not in found]
    assert not missing, f"the synthetic run wrote no {', '.join(missing)}: their columns would go unchecked"


def test_every_column_is_documented(out):
    names, wild = set(), set()
    for f in DOC_FILES:
        n, w = documented((DOCS / f).read_text())
        names, wild = names | n, wild | w
    assert not [w for w in wild if not w.strip("*")], "a backticked `*` in the docs would document every column"
    missing = {p.name: undocumented(header(p), names, wild) for p in sorted(out.glob("*.tsv"))}
    missing = {k: v for k, v in missing.items() if v}
    assert not missing, ("columns documented in neither docs/outputs.md nor docs/plants.md (name each in backticks):\n"
                         + "\n".join(f"  {k}: {', '.join(v)}" for k, v in missing.items()))


@pytest.mark.parametrize("table", sorted(OWN_SECTION))
def test_every_column_is_in_its_table_section(out, table):
    doc = OWN_SECTION[table]
    text = section((DOCS / doc).read_text(), table)
    assert text, f"docs/{doc} has no section whose heading names `asm.{table}.tsv`"
    missing = undocumented(header(out / f"{GENOME}.{table}.tsv"), *documented(text))
    assert not missing, f"`asm.{table}.tsv` columns missing from its section of docs/{doc}: {', '.join(missing)}"


def test_the_rule_reads_backticks_as_the_docs_write_them():
    names, wild = documented("| `a`, `b_c` | x |\n`d, e` and `f_*`; a span `g\nh` across a line; `i` `j`\n\n`k`")
    assert names == {"a", "b_c", "d", "e", "g h", "i", "j", "k"} and wild == {"f_*"}
    assert undocumented(["a", "e", "f_1", "f", "h", "z"], names, wild) == ["f", "h", "z"]
    text = "# T\n## `asm.x.tsv` (rows)\n`a`\n### sub\n`b`\n## next\n`c`\n"
    assert documented(section(text, "x"))[0] == {"asm.x.tsv", "a", "b"} and section(text, "y") == ""
