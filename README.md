# cenfuse

cenfuse finds the centromere of every chromosome in a genome assembly. It needs only the FASTA; an LTR
retrotransposon annotation, Hi-C and DNA methylation make it better. Each kind of evidence must first pass a test in
your genome (does it mark one place on each chromosome?) before it may confirm or move a call; untested evidence offers
at most a `low` best candidate. Every call comes with a confidence (`high`, `medium` or `low`) and a sentence saying
what supports it and what would settle it.

![How cenfuse works](docs/img/workflow.png)

*How cenfuse works ([PDF with legend](docs/img/workflow.pdf)). **A** six lines of evidence, three of them from the
FASTA alone, each tested in your genome before it may confirm or move a call; lines that pass are combined into one
call per chromosome, graded high, medium or low. **B** the test on real satellites: rice passes, tomato fails. **C** how
tested lines change a call.*

## What you give it

One folder per genome, every file named after the FASTA (`asm`). Only the FASTA is required.

```
mygenome/
├── asm.fa + asm.fa.fai    a chromosome-level assembly
├── tes/asm.gff3           LTR retrotransposons (LTRquest)        strongly recommended for plants
├── hic/asm.pairs.gz       Hi-C contacts (also .cool/.mcool, HiC-Pro, Juicer)
├── meth/asm.tsv           methylation samples (ONT or PacBio HiFi)
└── genes/asm.gff3         gene annotation
```

| your genome | its centromeres are usually | give cenfuse |
|---|---|---|
| mammals, birds; plants with a centromeric satellite (rice, *Brassica*, soybean, maize, *Arabidopsis*) | satellite arrays | the FASTA |
| plants without one (barley, *Aegilops*, potato, tomato), moss, green algae | clusters of LTR retrotransposons, or no repeat | the FASTA + `tes/`; Hi-C if you have it |
| yeasts, other fungi, *Plasmodium* | point, regional or AT-rich | the FASTA + Hi-C |
| primates, to place the kinetochore inside the array | alpha-satellite arrays | + CpG methylation |
| holocentric species (nematodes, moths, *Rhynchospora*) | no single centromere | anything: expect no confident call |

Hi-C or methylation that fails its test in your genome is reported but changes no call. How to make each input
(LTRquest, Hi-C pairs, methylation), which sequences count as chromosomes, and inputs from several plants:
[docs/inputs.md](docs/inputs.md).

## What you get

In `mygenome/cenfuse/`:

- **`cenfuse_report.md`**: start here. In words: what each kind of evidence found and whether it passed its test,
  what the centromeres are made of, the uncertain calls and the data that would settle them.
- **`asm.centromeres.tsv`** (also `.bed`, `.gff3`): one row per chromosome with the call, its `confidence`, the evidence
  for and against it, and in words what it stands on (`basis`) and why it may be wrong (`caveat`).
- **`cenfuse_overview.pdf`** (and `asm.overview.png`, `asm.closeups.png`): the figure. One centromere per chromosome,
  how confident the call is, and which evidence it rests on (examples below). Per genome: a whole-genome page,
  close-ups, then legends. Telomeres, rDNA arrays and other features are in `asm.loci.tsv`.

`high` and `medium` calls are confident; `low` is the best candidate, with another place not excluded. Read the
`caveat` of any call you rely on. Every file and column: [docs/outputs.md](docs/outputs.md).

<img src="docs/img/example_overview.png" width="640" alt="cenfuse genome page for Poa annua">

*Page 1 for* Poa annua *with all four inputs ([PDF with legends](docs/img/example_overview.pdf)). **A** every
chromosome to scale, with its satellite arrays, retrotransposon elements and Hi-C contact, and the centromere as one
triangle (solid black = high confidence, gray = medium, open = low). **B** for each chromosome, which evidence was used
(●), points elsewhere (×) or was not used (–), then the confidence. **C** whether each kind of evidence marks the
centromeres in this genome at all, against chance. **D** each signal's average around all centromeres, against random
places, with how many of them that signal placed (it rises there by construction).*

<img src="docs/img/example_closeups.png" width="640" alt="cenfuse close-ups for Poa annua">

*The close-up page* (**E**)*: each centromere is one band running through its signal tracks, shaded by confidence
(open with dashed edges when low). A bar in a track marks where that signal alone would put the centromere; a colored
bar well off the band is a disagreement, a gray bar a signal that was not used. The sentence under each name says how
the call was made.*

## How good it is

Numbers are from cenfuse 0.17.0 on every genome benchmarked before 0.17, scored against centromeres mapped on the
same assembly by CENP-A/CENH3 ChIP, CUT&Tag or DiMeLo-seq; the genomes added in 0.17 have their own bullet.

- **Satellite centromeres** (mammals, birds, most crops): from the FASTA alone, right on 425 of 450 chromosomes in 32
  genomes, 2 of them confidently wrong.
- **Retrotransposon or satellite-free centromeres** (barley, *Aegilops*, potato, 14 tomato genomes, moss, green
  algae): the FASTA alone finds 60 of 272; with `tes/` 207 of 272 (22 genomes), with `tes/` and Hi-C 189 of 231 (19
  genomes). Most confident errors are here (21 with `tes/`): treat a `medium` call that rests on an LTR core alone
  with caution.
- **Fungi and *Plasmodium*** (9 genomes): 50 of 98 chromosomes from the FASTA alone (AT-rich and transposon-rich
  centromeres), 94 of 98 with Hi-C.
- **Confidence means what it says.** Over all 63 genomes, each with its fullest set of inputs, `high` calls were right
  453 of 459 times (99%), `medium` 137 of 155 (88%), `low` 149 of 197 (76%).
- **Controls stay silent.** Randomized Hi-C changed no call (55 runs); holocentric genomes got no confident call (0 of
  46 chromosomes, with any input).
- **On genomes it was not developed on** (22 added in 0.17 and scored with 0.16.1, whose rules had never seen them):
  `high` calls were right 104 of 110 times (95%), `medium` 23 of 31 (74%), `low` 57 of 83 (69%). Most confident
  misses were centromeres that sit away from the genome's main satellite (*Nicotiana benthamiana*) and `medium` calls
  resting on one retrotransposon core; the 5 holocentric genomes got no confident call
  ([details](docs/benchmark.md#genomes-added-in-017)).

Right means overlapping the mapped domain in mammals, within 1 Mb of it in plants (8 Mb in barley), within 100 kb in
moss and algae, and within one Hi-C bin in fungi (250 kb without Hi-C, a lenient rule); every rule:
[benchmark/README.md](benchmark/README.md#scoring-rules). These genomes are not an independent test: cenfuse's rules
were developed while looking at them; the genomes added in 0.17 are an independent test of 0.16.1.

![cenfuse on every benchmark genome, by the data it was given](docs/img/performance.png)

*How often cenfuse finds the centromere, by the data it was given (every benchmark genome;
[PDF with legend](docs/img/performance.pdf)). **A** one genome with each set of inputs: calls, colored by outcome,
beside the mapped centromeres (black). **B** every chromosome with a mapped centromere, by inputs (dots beneath) and
centromere type; gray bars: chance; circles: single genomes. **C** every genome and every input set it was run with.
**D** how often `high`, `medium` and `low` calls were right. **E** Hi-C changes calls only when it carries real signal.*

Every genome, control and failure: [docs/benchmark.md](docs/benchmark.md); every genome's overview page:
[docs/gallery.md](docs/gallery.md); every result, one row per chromosome:
[benchmark/results.tsv](benchmark/results.tsv) (the genomes added in 0.17:
[benchmark/results_017.tsv](benchmark/results_017.tsv)).

## Install

```bash
git clone https://github.com/cwb14/cenfuse.git && cd cenfuse
mamba env create -f environment.yml && mamba activate cenfuse
pip install --no-deps .
python -m pytest -q -n 8      # optional: ~1,300 tests on tiny synthetic genomes, a few minutes
```

On linux-64, `mamba create -n cenfuse --file environment.lock.txt` gives the benchmark's exact environment. The tools
that make the optional inputs (LTRquest, bwa-mem2, pairtools, dorado, modkit) are installed separately.

## Run

```bash
mkdir mygenome && cp asm.fa mygenome/ && samtools faidx mygenome/asm.fa
cenfuse run mygenome -t 32
```

To add `tes/`, Hi-C or methylation later, put the files in place and run the same command: it reuses the finished
work. From the FASTA, a 1-Gb plant genome takes 13-50 min at 16 threads and up to 46 GB of memory, a 3-Gb T2T mammal
1-10 h at 12 threads and up to 74 GB; adding an input to a finished run takes minutes.

| option | meaning |
|---|---|
| `-t N` | threads [8] |
| `-o DIR` | output folder [`mygenome/cenfuse`] |
| `--force STAGE` | recompute a stage (`all`: every stage) |
| `-v` | per-step progress |

There are no thresholds to tune: each is a fixed constant with a stated reason. All options, resuming and negative
controls (`cenfuse control`): [docs/outputs.md](docs/outputs.md#options).

## More

[docs/design.md](docs/design.md) how each kind of evidence is tested and combined ·
[CHANGELOG.md](CHANGELOG.md) what changed in each version · [docs/future.md](docs/future.md) open problems ·
[benchmark/](benchmark/README.md) truth, scorer and results · License: MIT
