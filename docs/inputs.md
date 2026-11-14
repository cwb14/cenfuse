# Inputs

One folder per genome. Every file is named after the FASTA (`asm`):

```
mygenome/
├── asm.fa + asm.fa.fai    required: a chromosome-level assembly
├── asm.chromosomes.txt    optional: the chromosomes to analyze, one name per line
├── tes/asm.gff3           LTR retrotransposons (LTRquest)     strongly recommended for plants
├── hic/asm.pairs.gz       Hi-C contacts (or .cool/.mcool, HiC-Pro .allValidPairs, Juicer merged_nodups)
├── meth/asm.tsv           methylation samples (ONT or PacBio HiFi)
└── genes/asm.gff3         gene annotation (or .gff/.gtf)
```

| input | what it adds |
|---|---|
| FASTA | satellite families and their roles, telomeres, rDNA arrays, AT troughs, repeat density and repeat families, one sequence call per chromosome |
| `tes/` | the retrotransposon families that mark centromeres in this genome, learned from it; essential where centromeres have no satellite |
| `hic/` | a test of every call; tested Hi-C can move a call, or place one where the sequence gives none |
| `meth/` | CpG: the kinetochore inside each call (primates); CHG: a second line of evidence in plants (needs `tes/`) |
| `genes/` | gene density at every locus; rules out gene-rich repeats as centromeres |

Also read: `hic/asm.tsv` (several Hi-C libraries: [plants.md](plants.md#inputs)) and `meth/asm.exclude.bed` (regions the
methylation line drops; a 4th column names one sample, blank for all). The FASTA may be bgzipped (with `.fai` and `.gzi`), and
so may the `tes/` and `genes/` files; methylation files need a tabix index (`.tbi`, or `.csi` for sequences over 512 Mb).

## Which sequences are chromosomes

cenfuse keeps the sequences above the first tenfold drop in length and leaves out the tail of unplaced scaffolds. If that picks
the wrong set (a tiny dot chromosome, large unplaced scaffolds), list the chromosomes in `asm.chromosomes.txt`
([details](outputs.md#which-sequences-are-analyzed)). Several FASTAs in one folder are analyzed together and share satellite
family names.

## Inputs from other plants

Hi-C and methylation often come from a plant other than the assembled one, whose centromeres can sit elsewhere. cenfuse tells
which plant each input comes from and lets only the assembled plant's inputs decide ([plants.md](plants.md)).

## Recipes: making each input

### LTR retrotransposons (`tes/`)

With [LTRquest](https://github.com/cwb14/LTRquest), the annotator cenfuse is developed and benchmarked with:

```bash
ltrquest --genome asm.fa --max-rounds 1 --threads 32 \
  --detect-args "ltr-timeout=300" --run-sdust --no-plots --reboundary
mkdir -p mygenome/tes && cp asm_LTRs_all_depth_LTR_cleaned.gff3 mygenome/tes/asm.gff3
```

`--max-rounds 1` is how every benchmark annotation was made; `ltr-timeout=300` stops the slowest detectors from stalling on
centromeric satellites. `--proteins related.pep.fa` (genes of this or a related species, masked during detection) is allowed but
untested with cenfuse. The LTR line needs a sensitive annotation: with less sensitive ones (DANTE_LTR, EDTA's intact set) it lost
its gains in rice and potato. cenfuse reads LTRquest's `family=` labels; other GFF3 files work if each element's type contains
`LTR_retrotransposon` and carries `clade=` (e.g. EDTA + TEsorter), or if they are DANTE_LTR output as is.

### Hi-C (`hic/`)

As uniquely mapped, deduplicated pairs (about 100 M read pairs; map to this exact assembly):

```bash
bwa-mem2 index asm.fa && cut -f1,2 asm.fa.fai > asm.chrom.sizes && mkdir -p mygenome/hic
bwa-mem2 mem -5SP -T0 -t 32 asm.fa hic_R1.fq.gz hic_R2.fq.gz \
  | pairtools parse --min-mapq 30 --walks-policy 5unique --max-inter-align-gap 30 \
      --add-columns mapq --drop-sam --drop-seq -c asm.chrom.sizes \
  | pairtools sort --nproc 8 --memory 16G --tmpdir . | pairtools dedup \
  | pairtools select '(pair_type=="UU") or (pair_type=="UR") or (pair_type=="RU")' -o mygenome/hic/asm.pairs.gz
```

Omni-C, Micro-C and Arima libraries take the same recipe. Give at least ~70 M pairs for a 3-Gb genome: below that, Hi-C's tests
can fail and the calls fall back to the sequence evidence. HiC-Pro `asm.allValidPairs[.gz]` and Juicer
`asm.merged_nodups.txt[.gz]` files are read as they are. Convert a `.hic` file with [hictk](https://github.com/paulsengroup/hictk):
`hictk convert asm.hic mygenome/hic/asm.mcool`. Several libraries (replicates, or other plants) go in `hic/asm.tsv`
([plants.md](plants.md#inputs)).

### Methylation (`meth/`)

From ONT reads, one sample per plant (CHG needs `tes/` and a CHG file for every sample; CpG alone needs nothing else):

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
M=mygenome/meth; mkdir -p $M
modkit pileup sample1.bam - --ref asm.fa --cpg --combine-strands -t 32 \
  | awk -F'\t' '$4=="m"' | bgzip > $M/sample1.CG.bed.gz && tabix -p bed $M/sample1.CG.bed.gz
printf 'sample\tCG\nsample1\tsample1.CG.bed.gz\n' > $M/asm.tsv
```

A CG file can also be an ENCODE-style bedMethyl (e.g. from bisulfite sequencing), the BED of PacBio HiFi pb-CpG-tools (made with
`--min-mapq 1`: at 60, centromere arrays lose most of their data, and cenfuse warns of the hole) or a bedGraph of methylated
fraction (a public bigWig through `bigWigToBedGraph`), sorted, bgzipped and tabixed. Tracks encode methylation on different
scales (a PacBio bigWig, a nanopolish-derived bedGraph, ONT): the CpG test reads dips relative to each track's own variation and
by their width as well, so the ONT and PacBio tracks of HG002 give the same verdict. cenfuse pools calls in 5-kb bins, so low
coverage works. CHG files must be modkit bedMethyl. For sequences over 512 Mb, index with `samtools index -c` and
`tabix -C -p bed`.

### Genes (`genes/`)

Any GFF3 or GTF of this assembly, with sequence names as in the FASTA.
