"""GENOME_DIR discovery: one FASTA per genome plus optional genes/ (GFF3 or GTF), tes/ (GFF3), hic/ and meth/ inputs,
and the chromosome set every evidence line uses.

The chromosome set: an assembly also holds unplaced scaffolds, and every line needs whole chromosomes (the satellite
engine's prevalence tests count chromosomes; Hi-C and methylation need chromosome-scale sequences). cenfuse takes the
sequences GENOME_DIR/<genome>.chromosomes.txt lists (one name per line), else the methylation engine's rule
(track.chromosome_count: the sequences above the first tenfold drop in sorted length once half the assembly is
accounted for, all sequences when there is no such drop, unless the small ones among them are a gradual scaffold tail;
the --classify-from targets keep cenfuse 0.8's rule, which has no tail test). It hands the same set to all three
engines through a view of the genome directory: the FASTA is a symlink and its .fai lists the chromosomes only, so
nothing is copied.

Hi-C libraries: hic/<g>.tsv lists them (columns sample, file; sample names the plant, as in the methylation table),
or the lone hic/<g>.<ext> is one library named 'hic' (docs/plants.md)."""
from __future__ import annotations

import os
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

from .meth.track import chromosome_count, tenfold_count
from .sat.genomes import GENE_EXT, _fasta_name, load_fai
from .sat.te import TE_EXT
from .sat.util import SatelliteError
from .util import CenfuseError, log

HIC_MAIN = (".pairs.gz", ".pairs", ".mcool", ".cool")          # 4DN pairs and cooler: preferred when present
HIC_OTHER = (".allValidPairs", ".allValidPairs.gz", ".merged_nodups.txt", ".merged_nodups.txt.gz")   # HiC-Pro, Juicer
HIC_EXT = HIC_MAIN + HIC_OTHER        # hic.io.pairs_format reads the layout
CHROM_LIST = ".chromosomes.txt"      # GENOME_DIR/<genome>.chromosomes.txt: the user's chromosomes, one name per line
RULE_SOURCE = "length rule"          # GenomeInputs.chrom_source: what set the chromosomes
LEGACY_SOURCE = "length rule of cenfuse 0.8 (--classify-from target)"
READS = ".pairs[.gz], .cool, .mcool, .allValidPairs[.gz] (HiC-Pro) or .merged_nodups.txt[.gz] (Juicer)"
# Hi-C files cenfuse does not read, and how to make one it does (README recipe): better a clear stop than a silent run
# without the Hi-C the user gave
HIC_CONVERT = {".hic": "convert it: hictk convert {f} {stem}.mcool (or hic2cool convert {f} {stem}.mcool)",
               ".bam": "make pairs: pairtools parse ... | pairtools sort | pairtools dedup (README recipe)",
               ".cram": "make pairs: pairtools parse ... | pairtools sort | pairtools dedup (README recipe)",
               ".sam": "make pairs: pairtools parse ... | pairtools sort | pairtools dedup (README recipe)",
               ".validPairs": "use HiC-Pro's deduplicated .allValidPairs (or 4DN pairs)"}


@dataclass
class HicLibrary:
    name: str                       # the plant: a methylation sample of the same name is the same plant
    path: Path
    more: tuple = ()                # replicate libraries of the same plant (rows of hic/<g>.tsv with the same sample
                                    # name): their contacts are summed before the Hi-C engine (0.10.1)

    @property
    def files(self) -> list:
        return [self.path, *self.more]


@dataclass
class GenomeInputs:
    name: str
    fasta: Path                     # the user's FASTA
    fai: dict                       # every sequence in its .fai (name -> length), file order
    chroms: dict                    # the chromosome set (name -> length), .fai order
    chrom_source: str = RULE_SOURCE     # what set it: RULE_SOURCE, LEGACY_SOURCE or the list file's name
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
        lone = _one(folder, name, HIC_MAIN)
        others = [folder / f"{name}{e}" for e in HIC_OTHER if (folder / f"{name}{e}").exists()]
        if lone is not None and others:     # 0.9 asked to convert Juicer files to pairs: the original may stay beside
            log.warning(f"{folder}: {', '.join(o.name for o in others)} not used: {lone.name} is")
        elif lone is None:
            lone = _one(folder, name, HIC_OTHER)
        if lone is None:
            for ext, how in HIC_CONVERT.items():
                f = folder / f"{name}{ext}"
                if f.exists():
                    raise CenfuseError(f"{f}: cenfuse reads Hi-C as {READS}; {how.format(f=f, stem=name)}")
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
            raise CenfuseError(f"{table}: {n}: {p.name} is not one of {READS}")
        if not p.exists():
            raise CenfuseError(f"{table}: {n}: no file {p}")
        out.append(HicLibrary(n, p))
    if not out:
        raise CenfuseError(f"{table}: no rows")
    merged = {}                     # rows of one sample = replicate libraries of one plant, merged (0.10.1)
    for x in out:
        merged.setdefault(x.name, []).append(x.path)
    for n, ps in merged.items():
        if len({p.resolve() for p in ps}) < len(ps):
            raise CenfuseError(f"{table}: {n}: a file is listed twice (its contacts would count twice)")
        if len({p.name.endswith((".cool", ".mcool")) for p in ps}) > 1:
            raise CenfuseError(f"{table}: {n}: replicate libraries mix .cool/.mcool and pairs files (give one kind)")
        if len(ps) > 1:
            log.info(f"{table.name}: {n}: {len(ps)} replicate libraries merged")
    out = [HicLibrary(n, ps[0], tuple(ps[1:])) for n, ps in merged.items()]
    listed = {p.resolve() for x in out for p in x.files}
    for e in HIC_EXT:
        lone = folder / f"{name}{e}"
        if lone.exists() and lone.resolve() not in listed:
            log.warning(f"{lone}: not listed in {table.name}; not used")
    return out


def chromosome_set(fai: dict, legacy: bool = False) -> dict:
    """The n longest sequences, in .fai order; n = the methylation engine's chromosome count (track.chromosome_count),
    or cenfuse 0.8's (track.tenfold_count) when legacy: a --classify-from target is a contig assembly, with no
    chromosome/scaffold distinction for the tail test to find."""
    n = (tenfold_count if legacy else chromosome_count)(fai.values())
    keep = set(sorted(fai, key=lambda c: -fai[c])[:n])
    return {c: L for c, L in fai.items() if c in keep}


def _quoted_names(names: list, k: int = 5) -> str:
    """The first k names, quoted, and how many more."""
    return ", ".join(repr(n) for n in names[:k]) + (f" and {len(names) - k} more" if len(names) > k else "")


def read_chromosome_list(path: Path, fai: dict) -> dict:
    """The chromosome set the user lists in <genome>.chromosomes.txt: one sequence name per line; blank lines, lines
    starting with '#', CRLF, a byte-order mark and blanks around a name are accepted (written by hand or on Windows).
    Returns name -> length in .fai order. A name absent from the .fai, a name listed twice and a list without names
    stop the run: a list is the user's decision, never guessed at."""
    try:
        text = Path(path).read_text(encoding="utf-8-sig")
    except UnicodeDecodeError as e:                  # e.g. PowerShell's default UTF-16
        raise CenfuseError(f"{path}: not UTF-8 text ({e}); save the list as plain UTF-8 text") from None
    except OSError as e:
        raise CenfuseError(f"{path}: cannot read ({e})") from None
    names = [x for x in (line.strip() for line in text.splitlines()) if x and not x.startswith("#")]
    if not names:
        raise CenfuseError(f"{path}: no sequence names (one per line; blank lines and # comments are ignored)")
    twice = [n for n, k in Counter(names).items() if k > 1]
    if twice:
        raise CenfuseError(f"{path}: listed more than once: {_quoted_names(twice)}")
    unknown = [n for n in names if n not in fai]
    if unknown:
        raise CenfuseError(f"{path}: not sequences of the genome's .fai: {_quoted_names(unknown)}")
    keep = set(names)
    return {c: L for c, L in fai.items() if c in keep}


def _chromosomes(gdir: Path, name: str, fai: dict, legacy: bool) -> tuple[dict, str, list]:
    """A genome's chromosome set, what set it (GenomeInputs.chrom_source), and the notes that say so for the log and
    the report: the user's <genome>.chromosomes.txt, else the length rule (cenfuse 0.8's when legacy)."""
    listed = gdir / f"{name}{CHROM_LIST}"
    from_list = listed.is_file()
    if from_list:
        chroms, source = read_chromosome_list(listed, fai), listed.name
    else:
        chroms, source = chromosome_set(fai, legacy), LEGACY_SOURCE if legacy else RULE_SOURCE
    dropped, left_bp = len(fai) - len(chroms), sum(fai.values()) - sum(chroms.values())
    if from_list:
        note = f"chromosome set from {source}: {len(chroms)} of {len(fai)} sequences"
        if dropped:
            note += f"; the other {dropped} sequences ({left_bp:,} bp) are not analyzed"
        return chroms, source, [note]
    notes = [f"chromosome set by the {source}"]
    if dropped and not legacy and len(chroms) < tenfold_count(fai.values()):
        # the 0.9 tail cut fired: length cannot tell a small chromosome cut with the tail from a scaffold, so say so
        notes.append(f"the set ends at a gradual tail of short sequences: {dropped} shorter sequences ({left_bp:,} bp) "
                     f"are not analyzed, and the shortest one kept is {min(chroms.values()):,} bp; if chromosomes were "
                     f"left out, list every chromosome in {name}{CHROM_LIST}")
    elif dropped:
        notes.append(f"{dropped} shorter sequences (unplaced scaffolds; {left_bp:,} bp) are not analyzed")
    return chroms, source, notes


def discover(genome_dir, legacy: frozenset = frozenset()) -> list[GenomeInputs]:
    """The genomes of GENOME_DIR with their inputs and chromosome sets. legacy: the genomes that keep cenfuse 0.8's
    chromosome rule (the --classify-from targets); a <genome>.chromosomes.txt overrides either rule."""
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
        chroms, source, notes = _chromosomes(gdir, name, fai, name in legacy)
        g = GenomeInputs(name, f, fai, chroms, source, notes=notes)
        g.genes = _one(gdir / "genes", name, GENE_EXT)
        g.tes = _one(gdir / "tes", name, TE_EXT)
        g.hic = hic_libraries(gdir / "hic", name)
        g.meth = _one(gdir / "meth", name, (".tsv",))
        g.exclude = _one(gdir / "meth", name, (".exclude.bed",))
        if g.meth and not g.tes:
            g.notes.append("meth/ without tes/: the CHG line needs LTR retrotransposon classes and is skipped; CG "
                           "files still locate kinetochores (CpG line)")
        out.append(g)
    for f in sorted(gdir.iterdir()):                 # a list under a misspelled genome name would be ignored silently
        if f.is_file() and f.name.endswith(CHROM_LIST) and f.name[:-len(CHROM_LIST)] not in seen:
            log.warning(f"{f}: no genome {f.name[:-len(CHROM_LIST)]} in {gdir} (genomes: {', '.join(sorted(seen))}); "
                        f"not used")
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
