#!/usr/bin/env python3
"""Fill the {{NAME}} placeholders of the doc templates: README.md and docs/benchmark.md under OUT_ROOT.

  python -m benchmark.fill_docs TABLES OUT_ROOT [--src DIR] [--check]
A placeholder takes TABLES/NAME.md (the tables computed from a release's scores), else SRC/NAME.md (texts written by
hand from the same tables; SRC = benchmark/docs_src, which also holds README.template.md and benchmark.template.md).
A filled text may hold placeholders of its own (three rounds). Unknown names stay in place and are listed. OUT_ROOT =
the repository root rewrites the repo's README.md and docs/benchmark.md: edit the templates, not those files.
--check: write nothing. Port of v12 bench/fill_docs.py."""
import argparse
import re
from pathlib import Path

from benchmark import paths

DOCS = (("benchmark.template.md", "docs/benchmark.md"), ("README.template.md", "README.md"))


def fill(text, dirs, left, doc):
    def sub(m):
        name = m.group(1)
        for d in dirs:
            f = d / f"{name}.md"
            if f.exists():
                return f.read_text().rstrip("\n")
        if f"{doc}:{name}" not in left:
            left.append(f"{doc}:{name}")
        return m.group(0)

    for _ in range(3):
        text = re.sub(r"\{\{([A-Z0-9_]+)\}\}", sub, text)
    return text


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("tables", type=Path, help="folder of computed tables (md_tables, md_013, readme_table, rabl_table)")
    ap.add_argument("out", type=Path, help="root to write README.md and docs/benchmark.md under")
    ap.add_argument("--src", type=Path, default=paths.DOCS_SRC, help="templates and texts [benchmark/docs_src]")
    ap.add_argument("--check", action="store_true", help="list unfilled placeholders, write nothing")
    a = ap.parse_args(argv)
    left = []
    for tpl, doc in DOCS:
        new = fill((a.src / tpl).read_text(), [a.tables, a.src], left, Path(doc).name)
        if not a.check:
            (a.out / doc).parent.mkdir(parents=True, exist_ok=True)
            (a.out / doc).write_text(new)
    print("unfilled:", left or "none")


if __name__ == "__main__":
    main()
