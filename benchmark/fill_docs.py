#!/usr/bin/env python3
"""Fill the {{NAME}} placeholders of docs_src/benchmark.template.md and write OUT_ROOT/docs/benchmark.md.

  python -m benchmark.fill_docs TABLES OUT_ROOT [--src DIR] [--check]
A placeholder takes TABLES/NAME.md (a table computed from a release's scores), else SRC/NAME.md (a text written by hand;
SRC = benchmark/docs_src). A filled text may hold placeholders of its own (three rounds). Unknown names stay in place
and are listed. OUT_ROOT = the repository root rewrites the repo's docs/benchmark.md: edit the template, not that file.
--check: write nothing."""
import argparse
import re
from pathlib import Path

from benchmark import paths

TEMPLATE, DOC = "benchmark.template.md", "docs/benchmark.md"


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
    ap.add_argument("tables", type=Path, help="folder of computed tables (release.sh: OUT/mdTAG)")
    ap.add_argument("out", type=Path, help="root to write docs/benchmark.md under")
    ap.add_argument("--src", type=Path, default=paths.DOCS_SRC, help="template and hand-written texts [benchmark/docs_src]")
    ap.add_argument("--check", action="store_true", help="list unfilled placeholders, write nothing")
    a = ap.parse_args(argv)
    left = []
    new = fill((a.src / TEMPLATE).read_text(), [a.tables, a.src], left, Path(DOC).name)
    if not a.check:
        (a.out / DOC).parent.mkdir(parents=True, exist_ok=True)
        (a.out / DOC).write_text(new)
    print("unfilled:", left or "none")


if __name__ == "__main__":
    main()
