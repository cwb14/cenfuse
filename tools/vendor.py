#!/usr/bin/env python3
"""Vendor the three engines into cenfuse, logic unchanged: cenrep -> cenfuse.sat, cenhic -> cenfuse.hic,
cenmeth -> cenfuse.meth, with their test suites under tests/<sub>/. Only imports and test-helper module names are
rewritten (so three suites can share one pytest run). Writes VENDORED.md with the source commits.

usage: python tools/vendor.py SRC_DIR   (SRC_DIR holds git checkouts cenrep/, cenhic/, cenmeth/)"""
import ast
import re
import shutil
import subprocess
import sys
from pathlib import Path

ENGINES = {"sat": "cenrep", "hic": "cenhic", "meth": "cenmeth"}
HELPERS = {"synth.py": "{sub}_synth.py", "gffcheck.py": "{sub}_gffcheck.py"}
SKIP = shutil.ignore_patterns("__pycache__", "*.pyc", "*.egg-info", ".pytest_cache")


def rewrite_pkg(text: str, eng: str) -> str:
    text = re.sub(rf"^(\s*)from {eng}\.(\w+) import", r"\1from .\2 import", text, flags=re.M)
    return re.sub(rf"^(\s*)from {eng} import", r"\1from . import", text, flags=re.M)


def rewrite_test(text: str, eng: str, sub: str) -> str:
    pkg = f"cenfuse.{sub}"
    text = re.sub(rf"^(\s*)from {eng}(\.\w+)? import", lambda m: f"{m[1]}from {pkg}{m[2] or ''} import", text, flags=re.M)
    text = re.sub(rf"^(\s*)import {eng}\.(\w+)", rf"\1import {pkg}.\2", text, flags=re.M)
    text = re.sub(rf"(?<![\w.]){eng}\.(\w+)\.", rf"{pkg}.\1.", text)          # cenhic.cli.x -> cenfuse.hic.cli.x
    text = re.sub(rf"([\"']){eng}\.", rf"\1{pkg}.", text)                      # "cenrep.params.X" strings
    text = text.replace(f'"-m", "{eng}"', f'"-m", "{pkg}"')
    text = re.sub(r"^(\s*)from synth import", rf"\1from {sub}_synth import", text, flags=re.M)
    text = re.sub(r"^(\s*)import synth\b", rf"\1import {sub}_synth as synth", text, flags=re.M)
    text = re.sub(r"^(\s*)from gffcheck import", rf"\1from {sub}_gffcheck import", text, flags=re.M)
    return re.sub(r"^(\s*)from conftest import", rf"\1from {sub}_common import", text, flags=re.M)


def split_conftest(tdir: Path, sub: str) -> None:
    """Names other test modules import from conftest move to <sub>_common.py; conftest re-imports them."""
    conf = tdir / "conftest.py"
    wanted = set()
    for f in tdir.glob("test_*.py"):
        for m in re.finditer(rf"from {sub}_common import ([\w, ]+)", f.read_text()):
            wanted |= {x.strip() for x in m[1].split(",")}
    if not wanted or not conf.exists():
        return
    src = conf.read_text()
    lines = src.splitlines()
    moved, kept = [], []
    for node in ast.parse(src).body:                 # top-level statements, decorators included
        first = min([node.lineno] + [d.lineno for d in getattr(node, "decorator_list", [])]) - 1
        block = "\n".join(lines[first:node.end_lineno])
        names = {node.name} if isinstance(node, (ast.FunctionDef, ast.ClassDef)) else \
            {t.id for t in getattr(node, "targets", []) if isinstance(t, ast.Name)}
        (moved if names & wanted else kept).append(block)
    head = [b for b in kept if b.startswith(("import ", "from ", "sys.path"))]
    (tdir / f"{sub}_common.py").write_text(f'"""Helpers shared by the {sub} tests (moved out of conftest.py by '
                                           f'tools/vendor.py)."""\n' + "\n".join(head) + "\n\n\n" + "\n\n\n".join(moved) + "\n")
    names = ", ".join(sorted(wanted))
    conf.write_text("\n\n".join(kept) + f"\n\nfrom {sub}_common import {names}  # noqa: E402,F401\n")


def main(src: Path) -> None:
    root = Path(__file__).resolve().parents[1]
    rows = []
    for sub, eng in ENGINES.items():
        repo = src / eng
        commit = subprocess.run(["git", "-C", str(repo), "log", "-1", "--format=%h %s"], capture_output=True,
                                text=True, check=True).stdout.strip()
        dirty = subprocess.run(["git", "-C", str(repo), "status", "--porcelain"], capture_output=True, text=True,
                               check=True).stdout.strip()
        if dirty:
            sys.exit(f"{repo} has uncommitted changes; vendor a clean checkout")
        dst = root / "cenfuse" / sub
        shutil.rmtree(dst, ignore_errors=True)
        shutil.copytree(repo / eng, dst, ignore=SKIP)
        for f in dst.glob("*.py"):
            f.write_text(rewrite_pkg(f.read_text(), eng))
        tdst = root / "tests" / sub
        shutil.rmtree(tdst, ignore_errors=True)
        shutil.copytree(repo / "tests", tdst, ignore=SKIP)
        for old, new in HELPERS.items():
            if (tdst / old).exists():
                (tdst / old).rename(tdst / new.format(sub=sub))
        for f in list(tdst.glob("test_*.py")):
            if not f.name.startswith(f"test_{sub}_"):
                f = f.rename(tdst / f"test_{sub}_{f.name[5:]}")
        for f in tdst.glob("*.py"):
            f.write_text(rewrite_test(f.read_text(), eng, sub))
        split_conftest(tdst, sub)
        rows.append(f"| `cenfuse.{sub}` | {eng} | {commit} |")
    (root / "VENDORED.md").write_text(
        "# Vendored engines\n\nCopied by `tools/vendor.py` from clean git checkouts. Logic is unchanged; only "
        "imports (absolute to relative) and test-helper module names were rewritten. Local changes made after "
        "vendoring are listed below the table.\n\n| package | source | commit |\n|---|---|---|\n"
        + "\n".join(rows) + "\n\n## Local changes after vendoring\n\n- none\n")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    main(Path(sys.argv[1]).resolve())
