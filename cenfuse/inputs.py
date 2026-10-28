"""GENOME_DIR discovery: one FASTA per genome plus optional genes/ (GFF3 or GTF), tes/ (GFF3), hic/ and meth/ inputs,
and the chromosome set every evidence line uses.

The chromosome set: an assembly also holds unplaced scaffolds, and every line needs whole chromosomes (the satellite
engine's prevalence tests count chromosomes; Hi-C and methylation need chromosome-scale sequences). cenfuse takes the
methylation engine's rule (the sequences above the first tenfold drop in sorted length once half the assembly is accounted for; all
sequences when there is no such drop) and hands the same set to all three engines through a view of the genome
directory: the FASTA is a symlink and its .fai lists the chromosomes only, so nothing is copied.

Hi-C libraries: hic/<g>.tsv lists them (columns sample, file; sample names the plant, as in the methylation table),
or the lone hic/<g>.<ext> is one library named 'hic' (docs/plants.md)."""
from __future__ import annotations

import os
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

from .meth.track import chromosome_count
from .sat.genomes import GENE_EXT, _fasta_name, load_fai
from .sat.te import TE_EXT
from .sat.util import SatelliteError
from .util import CenfuseError, log

HIC_EXT = (".pairs.gz", ".pairs", ".mcool", ".cool")
# Hi-C files cenfuse does not read, and how to make one it does (README recipe): better a clear stop than a silent run
# without the Hi-C the user gave
HIC_CONVERT = {".hic": "convert it: hictk convert {f} {stem}.mcool (or hic2cool convert {f} {stem}.mcool)",
               ".bam": "make pairs: pairtools parse ... | pairtools sort | pairtools dedup (README recipe)",
               ".cram": "make pairs: pairtools parse ... | pairtools sort | pairtools dedup (README recipe)",
               ".sam": "make pairs: pairtools parse ... | pairtools sort | pairtools dedup (README recipe)",
               ".merged_nodups.txt": "convert it to 4DN pairs (columns chrom1 pos1 chrom2 pos2), or use the .hic "
                                     "Juicer made from it (hictk convert)",
               ".merged_nodups.txt.gz": "convert it to 4DN pairs (columns chrom1 pos1 chrom2 pos2), or use the .hic "
                                        "Juicer made from it (hictk convert)"}


@dataclass
class HicLibrary:
    name: str                       # the plant: a methylation sample of the same name is the same plant
    path: Path


@dataclass
class GenomeInputs:
    name: str
    fasta: Path                     # the user's FASTA
    fai: dict                       # every sequence in its .fai (name -> length), file order
    chroms: dict                    # the chromosome set (name -> length), .fai order
    genes: Path | None = None
    tes: Path | None = None
    hic: list = field(default_factory=list)     # HicLibrary per library
    meth: Path | None = None        # methylation samples table
    exclude: Path | None = None     # methylation exclusions BED
    notes: list = field(default_factory=list)

    @property
    def lines(self) -> list:
        return ["satellite"] + (["te"] if self.tes else []) + (["hic"] if self.hic else []) + \
            (["meth"] if self.meth else [])


def _one(folder: Path, name: str, exts) -> Path | None:
    hits = [folder / f"{name}{e}" for e in exts if (folder / f"{name}{e}").exists()]
    if len(hits) > 1:
        raise CenfuseError(f"{folder}: several inputs for genome {name}: {', '.join(h.name for h in hits)}")
    return hits[0] if hits else None


def hic_libraries(folder: Path, name: str) -> list:
    """The Hi-C libraries of genome `name`: the rows of hic/<name>.tsv, or the lone hic/<name>.<ext> as 'hic'."""
    table = folder / f"{name}.tsv"
    if not table.exists():
        lone = _one(folder, name, HIC_EXT)
        if lone is None:
            for ext, how in HIC_CONVERT.items():
                f = folder / f"{name}{ext}"
                if f.exists():
                    raise CenfuseError(f"{f}: cenfuse reads Hi-C as .pairs[.gz], .cool or .mcool; "
                                       f"{how.format(f=f, stem=name)}")
        return [HicLibrary("hic", lone)] if lone else []
    try:
        d = pd.read_csv(table, sep="\t", dtype=str).fillna("")
    except Exception as e:
        raise CenfuseError(f"{table}: cannot read ({e})") from None
    if not {"sample", "file"} <= set(d.columns):
        raise CenfuseError(f"{table}: needs columns 'sample' and 'file' (found {list(d.columns)})")
    out = []
    for r in d.to_dict("records"):
        n, f = r["sample"].strip(), r["file"].strip()
        if not n or "/" in n or any(ch.isspace() for ch in n):
            raise CenfuseError(f"{table}: bad sample name {n!r} (non-empty, without '/' or whitespace)")
        p = Path(f) if Path(f).is_absolute() else table.parent / f
        if not p.name.endswith(HIC_EXT):
            raise CenfuseError(f"{table}: {n}: {p.name} is not .pairs[.gz], .cool or .mcool")
        if not p.exists():
            raise CenfuseError(f"{table}: {n}: no file {p}")
        out.append(HicLibrary(n, p))
    if not out:
        raise CenfuseError(f"{table}: no rows")
    names = [x.name for x in out]
    if len(set(names)) < len(names):
        raise CenfuseError(f"{table}: duplicate sample names")
    listed = {x.path.resolve() for x in out}
    for e in HIC_EXT:
        lone = folder / f"{name}{e}"
        if lone.exists() and lone.resolve() not in listed:
            log.warning(f"{lone}: not listed in {table.name}; not used")
    return out


def chromosome_set(fai: dict) -> dict:
    """The n longest sequences, n = the methylation engine's chromosome count (tenfold length gap), in .fai order."""
    n = chromosome_count(fai.values())
    keep = set(sorted(fai, key=lambda c: -fai[c])[:n])
    return {c: L for c, L in fai.items() if c in keep}


def discover(genome_dir) -> list[GenomeInputs]:
    gdir = Path(genome_dir)
    if not gdir.is_dir():
        raise CenfuseError(f"{gdir}: not a directory")
    fas = sorted((f for f in gdir.iterdir() if f.is_file() and _fasta_name(f)), key=lambda f: f.name)
    if not fas:
        raise CenfuseError(f"{gdir}: no genome FASTA (*.fa/*.fasta/*.fna, plain or bgzip .gz, with .fai)")
    out, seen = [], set()
    for f in fas:
        name = _fasta_name(f)
        if name in seen:
            raise CenfuseError(f"{gdir}: genome {name} is given twice")
        seen.add(name)
        if f.name.endswith(".gz") and not Path(f"{f}.gzi").exists():
            raise CenfuseError(f"{f}: a compressed genome must be bgzip with .fai + .gzi (run: samtools faidx {f})")
        try:
            fai = load_fai(f)
        except SatelliteError as e:
            raise CenfuseError(str(e)) from None
        g = GenomeInputs(name, f, fai, chromosome_set(fai))
        g.genes = _one(gdir / "genes", name, GENE_EXT)
        g.tes = _one(gdir / "tes", name, TE_EXT)
        g.hic = hic_libraries(gdir / "hic", name)
        g.meth = _one(gdir / "meth", name, (".tsv",))
        g.exclude = _one(gdir / "meth", name, (".exclude.bed",))
        dropped = len(fai) - len(g.chroms)
        if dropped:
            g.notes.append(f"{dropped} shorter sequences (unplaced scaffolds; {sum(fai.values()) - sum(g.chroms.values()):,} bp)"
                           f" are not analyzed")
        if g.meth and not g.tes:
            g.notes.append("meth/ without tes/: the CHG line needs LTR retrotransposon classes and is skipped; CG "
                           "files still locate kinetochores (CpG line)")
        out.append(g)
    return out


def make_view(genomes: list[GenomeInputs], src_dir: Path, view: Path) -> Path:
    """view/<genome>.<ext> -> the user's FASTA (symlink; .gzi symlinked too) with a .fai of the chromosome set;
    view/genes and view/tes -> the user's folders. Rebuilt on every run (cheap)."""
    view.mkdir(parents=True, exist_ok=True)
    for g in genomes:
        link = view / g.fasta.name
        for p in (link, Path(f"{link}.fai"), Path(f"{link}.gzi")):
            if p.is_symlink() or p.exists():
                p.unlink()
        os.symlink(g.fasta.resolve(), link)
        if Path(f"{g.fasta}.gzi").exists():
            os.symlink(Path(f"{g.fasta}.gzi").resolve(), f"{link}.gzi")
        with open(f"{g.fasta}.fai") as src, open(f"{link}.fai", "w") as dst:
            for line in src:
                if line.split("\t", 1)[0] in g.chroms:
                    dst.write(line)
    for sub in ("genes", "tes"):
        link = view / sub
        if link.is_symlink() or link.exists():
            link.unlink()
        if (src_dir / sub).is_dir():
            os.symlink((src_dir / sub).resolve(), link)
    log.debug(f"genome view {view}: {', '.join(f'{g.name} ({len(g.chroms)} chromosomes)' for g in genomes)}")
    return view
