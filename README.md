# cenfuse

Finds the centromere of every chromosome in a genome assembly. It combines satellite repeats, LTR
retrotransposons, Hi-C and DNA methylation, lets them check each other, and says when it is unsure. Only the FASTA
is required.

![cenfuse overview of the Poa annua assembly](docs/img/example_overview.png)

*Poa annua, all four inputs ([PDF with legends](docs/img/example_overview.pdf)). **A** one row per chromosome:
Hi-C score (blue), CHG methylation (black), satellite arrays by role, TE density; black box = sequence call, gray
band = final call. **B** lines that support (●) or oppose (×) each call. **C** each line's call relative to the
sequence call. **D** Hi-C score around each call.*

## At a glance

| | |
|---|---|
| **Needs** | a chromosome-level FASTA. Optional: LTR retrotransposons, Hi-C, ONT methylation, genes |
| **Gives** | one centromere per chromosome, with a confidence and the evidence for and against it; telomeres, rDNA, knobs and other repeat blocks on the same map |
| **Priority** | right > uncertain > wrong. `low` means "uncertain", and cenfuse uses it rather than guess |
| **Runtime** | the satellite search dominates (3-Gb genome: ~5 h at 16 threads); Hi-C and fusion take minutes; reruns resume |

## How well it works

Right = the call overlaps the CENP-A/CENH3 domain (mammals) or lies within 4 Hi-C bins of it (plants).
Figure, details and every failure: [`docs/benchmark.md`](docs/benchmark.md). Every genome's overview plot:
[`docs/gallery.md`](docs/gallery.md).

| genomes | truth | inputs | chromosomes | right | confidently wrong | wrong but `low`, or no call |
|---|---|---|---|---|---|---|
| 9 T2T mammals (human, great apes, macaque, pig, horse) | CENP-A | FASTA only | 157 | 148 | 5 | 4 |
| same | CENP-A | + genes + Hi-C | 157 | 149 | 4 | 4 |
| rice | CENH3 | FASTA + Hi-C | 12 | 12 | 0 | 0 |
| barley | CENH3 | FASTA + Hi-C | 7 | 6 | 0 | 1 |
| potato | CENH3 | FASTA + Hi-C | 11 | 2 | 2 | 7 |
| horse (EquCab3) | CENP-A | FASTA + Hi-C | 5 | 3 (all `low`) | 0 | 2 |

Controls: on randomized Hi-C (29 runs) the Hi-C test never passed and no call changed.

| strong | weak |
|---|---|
| Satellite centromeres from the FASTA alone (mammals 148/157, rice 12/12) | Satellite-free centromeres without `tes/` (potato 2/11; horse chr11 has no call) |
| `high` is reliable: 119/122 `high` calls on CENP-A truth were right | Confidence is conservative: many `medium`/`low` calls are right too (horse 17/18) |
| Hi-C, satellites and TEs are tested against each other, not trusted on their own word | When a small array is the active one and a large array of the same satellite sits elsewhere, the large one is called, confidently (one orangutan chromosome) |
| Hi-C must pass a test in each genome before it can act; it never passed on randomized Hi-C | Hi-C cannot see inside repeat arrays: calls span the array (median 0.8-5.6 Mb), not the 0.1-0.5 Mb kinetochore |
| In grasses, Hi-C can overrule a wrong satellite line (barley 0 → 6/7) | Hi-C failed its test in 6/14 mammals (pig, horse, rhesus, sheep, goat, giraffe) and then only reports |
| Rejects interspersed repeats and handles acrocentric karyotypes without annotation | Acrocentric detection needs ≥ ~16 chromosomes |

## Install

```bash
git clone https://github.com/cwb14/cenfuse.git && cd cenfuse
mamba env create -f environment.yml && mamba activate cenfuse
pip install --no-deps .
pytest -q -n 8          # optional: 460 tests, ~2 min
```

## Quick start

```bash
mkdir mygenome && cp asm.fa mygenome/ && samtools faidx mygenome/asm.fa
cenfuse run mygenome -t 32
```

Read `mygenome/cenfuse/asm.centromeres.tsv` (one row per chromosome), and look at `cenfuse_overview.pdf`.

## Inputs

One folder per genome. Every file uses the FASTA's name (`asm`):

```
mygenome/
├── asm.fa + asm.fa.fai    required: chromosome-level assembly (unplaced scaffolds are skipped)
├── tes/asm.gff3           LTR retrotransposons    strongly recommended for plants; needed for methylation
├── hic/asm.pairs.gz       Hi-C contacts           (or .cool/.mcool)
├── meth/asm.tsv           ONT methylation samples table
└── genes/asm.gff3         gene annotation
```

| input | adds |
|---|---|
| FASTA | satellite families found de novo, each array's role, telomeres, rDNA, one sequence call per chromosome |
| `tes/` | CRM retrotransposon cores, which mark grass centromeres (without them, satellite calls were wrong on barley 7/7 and potato 9/11) |
| `hic/` | a functional test of every call; validated Hi-C can move a call |
| `meth/` | a second functional line; per-plant verdicts when samples come from several plants |
| `genes/` | gene density at every locus; helps reject gene-rich and subtelomeric repeats |

Several FASTAs in one folder are analyzed together and share satellite family names.

<details>
<summary><b>Recipes: making each input</b> (click)</summary>

**LTR retrotransposons** with [LTRquest](https://github.com/cwb14/LTRquest):

```bash
ltrquest --genome asm.fa --proteins related.pep.fa --threads 32      # uncompressed FASTA; --proteins optional
mkdir -p mygenome/tes && cp asm_LTRs_all_depth_LTR_cleaned.gff3 mygenome/tes/asm.gff3
```
Set `--mutation-rate` for your organism (default 3e-8, plants). Any GFF3 works if each element's type contains
`LTR_retrotransposon` and it carries `clade=` (e.g. EDTA + TEsorter).

**Hi-C** to uniquely mapped, deduplicated pairs (~100 M read pairs is enough; map to this exact assembly):

```bash
bwa-mem2 index asm.fa && cut -f1,2 asm.fa.fai > asm.chrom.sizes && mkdir -p mygenome/hic
bwa-mem2 mem -5SP -T0 -t 32 asm.fa hic_R1.fq.gz hic_R2.fq.gz \
  | pairtools parse --min-mapq 30 --walks-policy 5unique --max-inter-align-gap 30 \
      --add-columns mapq --drop-sam --drop-seq -c asm.chrom.sizes \
  | pairtools sort --nproc 8 --memory 16G --tmpdir . | pairtools dedup \
  | pairtools select '(pair_type=="UU") or (pair_type=="UR") or (pair_type=="RU")' -o mygenome/hic/asm.pairs.gz
```

**Methylation** from ONT reads, one sample per plant (CHG required; needs `tes/` too):

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
Optional `meth/asm.exclude.bed`: regions to ignore (column 4 = a sample name, or empty for all).

**Genes:** any GFF3 of this assembly with `gene` lines.

</details>

## Outputs (`mygenome/cenfuse/`)

| file | content |
|---|---|
| **`asm.centromeres.tsv`** (+ `.bed`, `.gff3`) | **the answer**: per chromosome, the call, `relation`, `confidence`, supporting and opposing lines, and each line's own call |
| `cenfuse_overview.pdf` | the figure above |
| `cenfuse_summary.tsv` | per genome: what each line found, whether Hi-C passed its tests, warnings |
| `asm.loci.tsv`, `asm.features.gff3` | every other feature (telomeres, rDNA, knobs, subtelomeric blocks, Hi-C hubs, LTR cores, former centromeres) with each line's value there |
| `asm.hic_anchored.bedgraph` | the Hi-C centromere score per bin, for a browser |
| `asm.provenance.tsv`, `asm.plants.tsv` | with several plants: which plant each input comes from; each plant's centromere |
| `satellite/`, `hic/`, `methylation/` | each line's detailed outputs |

## Reading a call

| `relation` | meaning |
|---|---|
| `concordant` | Hi-C or methylation agrees with the satellite/TE call (best case) |
| `offset` | validated Hi-C (with a second line for a high call) moved the call; the old one is kept in `seq_start`/`seq_end` |
| `contested` | a validated line points elsewhere without enough support to move: call kept, one level lower |
| `functional_only` | no satellite/TE call; Hi-C placed it |
| `sequence_only` | only the satellite/TE call was usable (no Hi-C, Hi-C failed its test, or its peak was weak) |
| `none` | no call |

**`confidence`**: `high`, `medium`, or `low` = uncertain (the best candidate, but another place is not excluded).
Agreeing independent lines raise a kept call (at most two levels); an opposing line lowers it one. A call that Hi-C
moved or placed is `medium` (`high` with a second line).

| flag | meaning |
|---|---|
| `competing_cluster` | a second satellite cluster at least half the size; call is `low`, the other is in `seq_alt_start`/`seq_alt_end` |
| `hic_chose_alt` | validated Hi-C picked the competing cluster |
| `acrocentric` | the family sits at one end of most chromosomes (acrocentric karyotype) |
| `fold_pair` | Hi-C's centromere cluster and the chromosome's arm fold agree here, away from the satellite call: the call moved |
| `fold_offset`, `hic_denovo_offset` | the arm fold / Hi-C's prior-free peak points elsewhere (one level lower each) |
| `hic_medium_peak` | no sequence call; a clean medium Hi-C peak placed a `low` call |
| `moved_with_te`, `_meth`, `_fold` | a high call moved: a second line (TE core, methylation, arm fold) agreed with Hi-C |
| `near_end`, `low_satellite`, `diffuse`, `no_dominant_satellite` | why the satellite engine lowered its confidence |

Per-genome tests (`*.centromeres.tsv` and `cenfuse_summary.tsv`): `hic_gate` (`pass`/`fail`/`untested`: does Hi-C
find the other chromosomes' centromeres?) and `fold_gate` (does the arm fold agree with Hi-C's own peaks?).

## Options

| option | meaning |
|---|---|
| `-o DIR` | output folder [`GENOME_DIR/cenfuse`] |
| `-t N` | threads [8] |
| `--satellite-dir DIR` | reuse the `satellite/` folder of an earlier run |
| `--seed-engine fastan\|trf` | tandem-repeat finder for satellite seeds [fastan] |
| `--seed-cache DIR` | share tandem-repeat results between runs |
| `--classify-from T=S` | genome T (e.g. a contig assembly) takes its satellite families from genome S |
| `--force STAGE` | recompute one stage |
| `-v` | per-step progress |

There are no tuning knobs: every threshold is a fixed constant with its reason in a `params.py`.

## Several plants

Methylation and Hi-C often come from plants other than the assembled one, whose centromeres can sit elsewhere. Give
one row per plant (`meth/asm.tsv`; `hic/asm.tsv` with columns `sample`, `file`); the same `sample` name means the
same plant. cenfuse genotypes each input, lets only the assembled plant's (or unknown) inputs decide the calls, and
scores every plant on its own (`shared` or `shifted`). Details: [`docs/plants.md`](docs/plants.md).

## Limitations

- **Satellite-free centromeres** need `tes/` (CRM cores) or validated Hi-C; the satellite line alone misses them.
- **A genome-wide wrong satellite line** can be overruled by Hi-C only through the arm fold, which needs Rabl nuclei
  (grasses). Elsewhere `cenfuse_summary.tsv` flags `hic_sequence_conflict`; add `tes/` or methylation.
- **Resolution**: Hi-C locates the unmappable array, not the kinetochore inside it.
- **One wrong Hi-C peak** in a genome where Hi-C mostly works passes the test; only high calls need a second line to move.
- **Methylation** needs ONT CHG calls and `tes/`; plants below ~4x often lack data at centromeres.
- **Hi-C provenance**: pairs or cool files cannot be genotyped, so their plant of origin is usually `unknown`.

## Documentation

| file | read it for |
|---|---|
| [`docs/benchmark.md`](docs/benchmark.md) | all results as a figure and tables: per genome, per line, controls, failures |
| [`docs/gallery.md`](docs/gallery.md) | cenfuse's overview plot for every benchmark genome |
| [`docs/design.md`](docs/design.md) | how each line works, the fusion rules, constants, and alternatives that were rejected |
| [`docs/plants.md`](docs/plants.md) | provenance and per-plant scoring |

## License

MIT (see `LICENSE`).
