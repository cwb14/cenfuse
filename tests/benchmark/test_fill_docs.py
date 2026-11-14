"""benchmark/fill_docs.py: computed tables win over hand-written texts, nested placeholders fill, unknown names stay."""
from benchmark import fill_docs, paths


def test_fill_order_nesting_and_unknown(tmp_path, capsys):
    tables, src, out = tmp_path / "tables", tmp_path / "src", tmp_path / "out"
    tables.mkdir()
    src.mkdir()
    (tables / "T.md").write_text("computed {{NOTE}}\n")
    (src / "T.md").write_text("hand-written\n")              # shadowed by the computed table
    (src / "NOTE.md").write_text("(see {{FIG}})\n")
    (src / "FIG.md").write_text("figure 1\n")
    (src / "benchmark.template.md").write_text("# B\n{{T}}\n{{MISSING}}\n")
    fill_docs.main([str(tables), str(out), "--src", str(src)])
    assert (out / "docs/benchmark.md").read_text() == "# B\ncomputed (see figure 1)\n{{MISSING}}\n"
    assert [p.name for p in out.rglob("*") if p.is_file()] == ["benchmark.md"]     # the README is written by hand
    assert "benchmark.md:MISSING" in capsys.readouterr().out
    fill_docs.main([str(tables), str(tmp_path / "none"), "--src", str(src), "--check"])
    assert not (tmp_path / "none").exists()


def test_repository_template_present():
    assert "{{" in (paths.DOCS_SRC / "benchmark.template.md").read_text()
    assert [p.name for p in paths.DOCS_SRC.glob("*.template.md")] == ["benchmark.template.md"]
