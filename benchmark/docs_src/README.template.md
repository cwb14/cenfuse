# cenfuse

cenfuse finds the centromere of every chromosome in a genome assembly. It needs only the FASTA; LTR
retrotransposons, Hi-C, DNA methylation and genes make it better. Each kind of evidence calls the centromere on its
own, and cenfuse tests, in your genome, whether that evidence can be trusted before it may confirm or move a call.
Every call comes with a confidence and, in words, what speaks for and against it.

![cenfuse overview of the Poa annua assembly](docs/img/example_overview.png)

*cenfuse's overview of the* Poa annua *assembly, from all four inputs ([PDF with legends](docs/img/example_overview.pdf)).
**A** one row per chromosome: Hi-C score (blue), CHG methylation (black), satellite arrays colored by role,
retrotransposon density (orange); black box = the sequence call, gray band = the final call; ▶ ▼ telomere arrays at
the ends and inside; ◆ ◇ 45S and 5S rDNA arrays. **B** the evidence that supports (●) or opposes (×) each call.
**C** where each kind of evidence puts the centromere. **D** the Hi-C score around each call. On Pa1A validated Hi-C peaks 21 Mb
from the satellite that the LTR signature also marks (TE ×), and the copies of the centromeric LTR families there are young (median
0 years) while those at the satellite are old (620,000 years; column `age`): the centromere moved, and the call (arrow) moves with
it, where earlier analyses put it ([docs/future.md #31](docs/future.md#31-a-centromere-that-left-its-satellite-and-its-ltr-families)).*

## How it works

![How cenfuse works](docs/img/slides/03_workflow.png)

Each kind of evidence calls the centromere on its own, and each must pass a test in the genome at hand before it may act.
Most tests rest on one fact: a chromosome has one centromere. So whatever marks centromeres falls once per chromosome,
while anything else falls in proportion to chromosome length.

1. **Satellites** (always). cenfuse finds the genome's tandem-repeat families from scratch and gives each a role
   (centromeric, subtelomeric, knob-like, dispersed, rDNA) from where and how often it occurs across the
   chromosomes. Its call is the main array of the centromeric family on each chromosome. The test: do the N largest
   arrays of those families cover more of the N chromosomes than chance? Where they do not (tomato, potato, barley,
   *Ae. tauschii*), a satellite call cannot be confident on its own.
2. **AT troughs** (always, where most chromosomes have no other call). Many regional centromeres sit on the most AT-rich
   stretch of their chromosome, whatever their sequence. The same test, on the genome's deepest GC troughs:
   *Plasmodium*, *Malassezia* and *Neurospora* pass from the FASTA alone.
3. **LTR retrotransposons** (`tes/`). cenfuse learns which families mark centromeres in this genome, at family level and
   for each chromosome from the other chromosomes, so any lineage can do it: CRM in grasses, Bryco in moss, Osser in
   *Chlorella*, a single Tekay family in tomato. A decisive signature can replace a satellite call it contradicts.
4. **Hi-C** (`hic/`). Centromeres meet in the nucleus. Before Hi-C may act on a chromosome, it has to find the
   centromeres of the *other* chromosomes in your library. Where those calls are wrong genome-wide, Hi-C is tested
   against the LTR cores instead, and so is a cis signal: in tomato the two flanks of each centromere touch about twice
   as often as expected. Without LTR cores, Hi-C's two readouts test each other: the cis signal's extremes anchor the
   contacts between chromosomes, and where those find the other chromosomes' extremes on most chromosomes, the cis
   signal may act (tomato from the FASTA and Hi-C: 10/12 right, from 1/12). Where nuclei hold their centromeres at one
   pole (a Rabl configuration), the way chromosome arms fold and align tests Hi-C, and Hi-C tells which end of an
   acrocentric chromosome holds the centromere. In genomes whose centromeres carry no repeat, Hi-C places the calls
   itself.
5. **Methylation** (`meth/`). In primates, CpG methylation dips where the kinetochore sits: where cenfuse shows that
   the dips mark centromeres in your genome, it places the kinetochore inside each call. In plants, CHG methylation
   can be a second line of evidence, though it [changed no call](docs/future.md#16-the-chg-line-has-changed-no-call)
   in the benchmark.
6. **Integration.** The sequence call stands unless tested evidence moves it. Agreement raises its confidence,
   disagreement lowers it, and two tested lines at one locus outrank one: a `high` call moves only when a second line
   sits at Hi-C's peak (an LTR core, methylation, or younger copies of the centromeric LTR families than at the call: a
   centromere that moved). No line speaks where it cannot see: contact readouts never judge the chromosome ends, which
   cluster in the nucleus whatever the centromere. `low` means "unsure": cenfuse prefers right, then uncertain, over
   wrong. How well the lines agree with each other, against chance, is reported for every genome.

Details: [docs/design.md](docs/design.md).

## How well it works

Calls were scored against centromeres mapped by CENP-A/CENH3 ChIP (or CUT&Tag, DiMeLo-seq) on the same assembly.
A call is right when it overlaps the CENP-A domain (mammals), lies within 4 Hi-C bins of the CENH3 domain (plants:
about 1 Mb; moss and algae: 100 kb), or within 1 Hi-C bin (fungi, *Plasmodium*); *confident* = `high` or `medium`.
Every number is from cenfuse 0.14.0, whose calls equal 0.13.0's on all 260 benchmark runs. Per genome, the controls and every failure are in
[docs/benchmark.md](docs/benchmark.md); every genome's overview page is in the [gallery](docs/gallery.md).

{{README_TABLE}}

{{README_BULLETS}}

## Install

```bash
git clone https://github.com/cwb14/cenfuse.git && cd cenfuse
mamba env create -f environment.yml && mamba activate cenfuse      # or, on linux-64, the benchmark's exact environment:
                                                                    # mamba create -n cenfuse --file environment.lock.txt
pip install --no-deps .
pytest -q -n 8          # optional: about 1,150 tests, a few minutes; CENFUSE_REQUIRE_TOOLS=1 fails instead of skipping without the tools
```

## Quick start

```bash
mkdir mygenome && cp asm.fa mygenome/ && samtools faidx mygenome/asm.fa
cenfuse run mygenome -t 32
```

Then read `mygenome/cenfuse/cenfuse_report.md` (what was found, in words), `asm.centromeres.tsv` (one row per
chromosome) and `cenfuse_overview.pdf`.

**Runtime and memory.** The satellite search takes most of the time and grows with genome size and satellite content. Fresh
runs from the FASTA (2026-11-07; the search code is unchanged in 0.14): the allotetraploid *Ae. geniculata* (7.9 Gb) 64 min at 32 threads, the wild potato *S. verrucosum*
(0.7 Gb) 16 min at 16, 13 tomatoes (0.8-1.2 Gb) 13-50 min at 16; adding `tes/` to the finished run, 1-12 min. From before 0.11 (whose search is up to 4.9x faster
where long repeat models hold low-complexity tracts): T2T mammals (2.6-3.6 Gb) 1-10 h at 12 threads, mouse 9 h at 16. Hi-C and
the integration take minutes. Peak memory (the largest process): 17.5 GB for *Ae. geniculata* at 32 threads, 7-46 GB for the tomatoes at 16, 8-74 GB for T2T
mammals at 12. A rerun resumes where it stopped, fresh runs give identical results (Linux), and every run records its versions,
tools, inputs and cost (CPU hours; the peak memory of all its processes) in `run_info.json`.

## Inputs

One folder per genome. Every file is named after the FASTA (`asm`):

```
mygenome/
├── asm.fa + asm.fa.fai    required: a chromosome-level assembly
├── asm.chromosomes.txt    optional: the chromosomes to analyze, one name per line
├── tes/asm.gff3           LTR retrotransposons    strongly recommended for plants (LTRquest recipe below)
├── hic/asm.pairs.gz       Hi-C contacts           (or .cool/.mcool, HiC-Pro .allValidPairs, Juicer merged_nodups)
├── meth/asm.tsv           methylation samples     (ONT or PacBio HiFi)
└── genes/asm.gff3         gene annotation         (or .gff/.gtf)
```

| input | what it adds |
|---|---|
| FASTA | satellite families and their roles, telomeres, rDNA arrays, one sequence call per chromosome |
| `tes/` | the LTR retrotransposon families that mark centromeres in this genome, learned from it (from the FASTA alone barley gets 0/7 right, potato 2/11, tomato 1/12; with `tes/` 7/7, 10/11, 6/12); with Hi-C, they also test Hi-C where the satellite calls cannot |
| `hic/` | a test of every call; tested Hi-C can move a call, and where the FASTA gives no call (point, regional, AT-rich centromeres) it places one; a cis signal at satellite-free centromeres, tested by the LTR cores or, without `tes/`, by Hi-C's own contacts between chromosomes (tomato: 10/12 either way) |
| `meth/` | CpG: the kinetochore inside each call (primates); CHG (with `tes/`): a possible second line of evidence in plants |
| `genes/` | gene density at every locus; rejects gene-rich repeats as centromeres |

**Which sequences are chromosomes.** cenfuse keeps the sequences above the first tenfold drop in length and leaves
out a tapering tail of unplaced scaffolds. If that picks the wrong set (a tiny dot chromosome, large unplaced
scaffolds), list every chromosome in `asm.chromosomes.txt` ([details](docs/outputs.md#which-sequences-are-analyzed)).
Several FASTAs in one folder are analyzed together and share satellite family names.

<details>
<summary><b>Recipes: making each input</b> (click)</summary>

**LTR retrotransposons** with [LTRquest](https://github.com/cwb14/LTRquest), the annotator cenfuse is developed and
benchmarked with. Generate the `tes/` file with this command (the benchmark's annotations used it without `--proteins`):

```bash
ltrquest --genome asm.fa --proteins related.pep.fa --max-rounds 1 --threads 32 \
  --detect-args "ltr-timeout=300" --run-sdust --no-plots --reboundary
mkdir -p mygenome/tes && cp asm_LTRs_all_depth_LTR_cleaned.gff3 mygenome/tes/asm.gff3
```

| flag | why |
|---|---|
| `--proteins` | optional: a protein FASTA of this or any related species, to mask genes during detection (untested with cenfuse: no benchmark annotation used it) |
| `--max-rounds 1` | one detection pass, as in every benchmark run; LTRquest's default keeps searching for nested elements for up to 10 rounds, which cenfuse's tests were not benchmarked on |
| `--detect-args "ltr-timeout=300"` | caps the slowest detectors at 5 min per chunk: they stall on centromeric satellites |
| `--run-sdust` | drops candidates made mostly of low-complexity sequence |
| `--reboundary` | repairs element ends that lack a target-site duplication, and merges duplicate calls of one element |
| `--no-plots` | skips LTRquest's plots: cenfuse reads only the GFF3 |

All 17 `tes/` annotations of the benchmark, and those made for the genomes held out from 0.12, come from this command, except
that they were made without `--proteins`. cenfuse reads LTRquest's `family=` labels to learn the centromeric families; a file without
them uses lineages only. The LTR line needs a sensitive annotation: with less sensitive ones (DANTE_LTR, EDTA's intact set, half of
LTRquest's elements) CRM concentrated on too few rice chromosomes to decide in 0.11, and potato lost its gain from `tes/`.
Other GFF3 files are read too: each element's type must contain `LTR_retrotransposon` and carry `clade=` (e.g. EDTA +
TEsorter), or be a DANTE_LTR GFF3 as it is. A file without lineages, or on other sequence names, is reported and the
LTR line is left out.

**Hi-C** to uniquely mapped, deduplicated pairs (about 100 M read pairs; map to this exact assembly):

```bash
bwa-mem2 index asm.fa && cut -f1,2 asm.fa.fai > asm.chrom.sizes && mkdir -p mygenome/hic
bwa-mem2 mem -5SP -T0 -t 32 asm.fa hic_R1.fq.gz hic_R2.fq.gz \
  | pairtools parse --min-mapq 30 --walks-policy 5unique --max-inter-align-gap 30 \
      --add-columns mapq --drop-sam --drop-seq -c asm.chrom.sizes \
  | pairtools sort --nproc 8 --memory 16G --tmpdir . | pairtools dedup \
  | pairtools select '(pair_type=="UU") or (pair_type=="UR") or (pair_type=="RU")' -o mygenome/hic/asm.pairs.gz
```

Omni-C, Micro-C and Arima libraries take the same recipe. Give at least ~70 M pairs for a 3-Gb genome: below that,
Hi-C's tests can fail and the calls fall back to the sequence evidence. HiC-Pro `asm.allValidPairs[.gz]` and Juicer
`asm.merged_nodups.txt[.gz]` files are read as they are. Convert a `.hic` file with
[hictk](https://github.com/paulsengroup/hictk): `hictk convert asm.hic mygenome/hic/asm.mcool`.

**Methylation** from ONT reads, one sample per plant (CHG needs `tes/` too; CpG alone needs nothing else):

```bash
dorado basecaller sup pod5/ --modified-bases 5mC_5hmC --min-qscore 10 > plant1.ubam     # not the CpG-only model
samtools fastq -T MM,ML,MN plant1.ubam | minimap2 -t 32 -ax lr:hq -y asm.fa - \
  | samtools view -u -q 20 -F 0x904 - | samtools sort -@ 8 -o plant1.bam && samtools index plant1.bam
M=mygenome/meth; mkdir -p $M
for ctx in CHG CHH; do                                  # no --combine-strands for CHG: strands are read separately
  modkit pileup plant1.bam - --ref asm.fa --modified-bases 5mC 5hmC --motif $ctx 0 -t 32 \
    | awk -F'\t' '$4=="m"' | bgzip > $M/plant1.$ctx.bed.gz && tabix -p bed $M/plant1.$ctx.bed.gz
done
modkit pileup plant1.bam - --ref asm.fa --modified-bases 5mC 5hmC --motif CG 0 --combine-strands -t 32 \
  | awk -F'\t' '$4=="m"' | bgzip > $M/plant1.CG.bed.gz && tabix -p bed $M/plant1.CG.bed.gz
printf 'sample\tCHG\tCG\tCHH\nplant1\tplant1.CHG.bed.gz\tplant1.CG.bed.gz\tplant1.CHH.bed.gz\n' > $M/asm.tsv
```

For CpG only (animals; Dorado `5mCG_5hmCG` calls), one CG file per sample:

```bash
modkit pileup sample1.bam - --ref asm.fa --cpg --combine-strands -t 32 \
  | awk -F'\t' '$4=="m"' | bgzip > $M/sample1.CG.bed.gz && tabix -p bed $M/sample1.CG.bed.gz
printf 'sample\tCG\nsample1\tsample1.CG.bed.gz\n' > $M/asm.tsv
```

A CG file can also be an ENCODE-style bedMethyl (e.g. from bisulfite sequencing) or a bedGraph of methylated fraction
(a public bigWig through `bigWigToBedGraph`, or PacBio HiFi scores from pb-CpG-tools), sorted, bgzipped and tabixed.
Prefer files with read counts: cenfuse pools the calls in 5-kb bins, so low coverage works. CHG files must be modkit
bedMethyl. An optional `meth/asm.exclude.bed` lists regions the CHG line ignores (column 4: a sample, or empty for all).

**Genes:** any GFF3 or GTF of this assembly, with sequence names as in the FASTA (cenfuse warns when they do not
match).

</details>

## Outputs (`mygenome/cenfuse/`)

| file | content |
|---|---|
| **`asm.centromeres.tsv`** (+ `.bed`, `.gff3`) | **the answer**: per chromosome, the call, its `confidence` and `relation`, the evidence for and against it, the call in words (`basis`, `caveat`), the CpG kinetochore, the chromosome's morphology (centromeric index, Levan class) and completeness (telomeres, gaps) |
| **`cenfuse_report.md`** | per genome, in words: what each kind of evidence found and whether it passed its test, what the centromeres are made of, the uncertain calls and what would settle them, the karyotype, the other features |
| `cenfuse_overview.pdf`, `asm.overview.png` | the figure above, with its legends; the first page as a PNG |
| `cenfuse_summary.tsv` | per genome: counts of calls, test results, the nucleus's 3D organization |
| `asm.loci.tsv`, `asm.features.gff3` | every other feature: telomeres (also inside chromosomes), 45S and 5S rDNA arrays, knobs, subtelomeric blocks, retrotransposon cores, Hi-C hubs |
| `asm.line_tests.tsv` | how well each pair of lines agrees in this genome, against chance (it changes no call) |
| `asm.te_signature.tsv` | with `tes/`: the LTR signature per chromosome: seed family, fold enrichment, the families learned from the other chromosomes, the core |
| `asm.hic_anchored.bedgraph` | the Hi-C centromere score along the genome, for a genome browser |
| `asm.provenance.tsv`, `asm.plants.tsv` | which plant each input comes from; with several plants, each plant's centromeres ([plants.md](docs/plants.md)) |
| `satellite/`, `hic/`, `methylation/` | each line's detailed outputs |

Every column, flag and option: [docs/outputs.md](docs/outputs.md).

## Reading a call

**`confidence`**: `high`, `medium`, or `low` = unsure (the best candidate, but another place is not excluded).

**`relation`**: how the final call relates to the sequence call.

| `relation` | meaning |
|---|---|
| `concordant` | Hi-C or methylation agrees with the sequence call |
| `offset` | tested Hi-C moved the call (the sequence call is kept in `seq_start`/`seq_end`) |
| `contested` | tested evidence points elsewhere but could not move the call: kept, one level lower |
| `functional_only` | no sequence call; Hi-C placed it |
| `sequence_only` | only the sequence call was usable |
| `none` | no call |

`basis` says what the call stands on and `caveat` why it may be wrong, ending with the input that would settle it.
`flags` lists the rules that acted ([all flags](docs/outputs.md#flags)).

**Several plants.** Methylation and Hi-C often come from plants other than the assembled one, whose centromeres can
sit elsewhere. Give one row per plant (`meth/asm.tsv`; `hic/asm.tsv` with columns `sample`, `file`); the same
`sample` name means the same plant, and several Hi-C rows of one plant are summed. cenfuse tells which plant each
input comes from, lets only the assembled plant's inputs (or those of unknown origin) decide, and scores every plant
on its own ([docs/plants.md](docs/plants.md)).

## Options

| option | meaning |
|---|---|
| `-o DIR` | output folder [`GENOME_DIR/cenfuse`] |
| `-t N` | threads [8] |
| `--satellite-dir DIR` | reuse the `satellite/` folder of an earlier run |
| `--seed-engine fastan\|trf` | tandem-repeat finder for the satellite search [fastan] |
| `--seed-cache DIR` | share tandem-repeat results between runs |
| `--classify-from T=S` | genome T (e.g. a contig assembly) takes its satellite families from genome S |
| `--force STAGE` | recompute a stage (`all`: every stage); a later stage reruns when what it reads changed. Adding `tes/` or `genes/` to an earlier run needs none: the satellite engine sees the new files and recomputes from classify |
| `-v` | per-step progress |

`cenfuse control KIND IN OUT` writes a negative control of one input: Hi-C pairs re-paired or relocated at random, GC blocks or
LTR positions shuffled ([outputs.md](docs/outputs.md#negative-controls-cenfuse-control)). Run cenfuse on the control: no test should pass.

There are no tuning knobs: every threshold is a fixed constant with its reason ([docs/design.md](docs/design.md#constants)).

## Documentation

| file | read it for |
|---|---|
| [`docs/design.md`](docs/design.md) | how each kind of evidence is tested and combined, the constants, the alternatives that were rejected |
| [`docs/benchmark.md`](docs/benchmark.md) | every result: per genome, per input, the controls, every failure |
| [`docs/gallery.md`](docs/gallery.md) | cenfuse's overview page for every benchmark genome |
| [`docs/outputs.md`](docs/outputs.md) | every output column and flag, the options, running and resuming |
| [`docs/plants.md`](docs/plants.md) | inputs from several plants: provenance and per-plant scoring |
| [`docs/future.md`](docs/future.md) | open problems, each with its evidence and the check it needs |
| [`docs/constants.md`](docs/constants.md) | every constant in the code with its reason, generated from the code |
| [`benchmark/`](benchmark/README.md) | the truth manifest, the scorer and the scripts that make this page's and benchmark.md's tables |
| [`.github/workflows/tests.yml`](.github/workflows/tests.yml) | the test suite in the benchmark's exact environment, every tool required |
| [`docs/img/slides/`](docs/img/slides/README.md) | figures for a talk on cenfuse |

## License

MIT (see `LICENSE`).
