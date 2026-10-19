# cenfuse

Centromeres from four independent lines of evidence: satellite arrays, LTR retrotransposons, Hi-C and DNA
methylation. Each line makes its own call; cenfuse tests the lines against each other, places one centromere per
chromosome with the lines that support or oppose it, and maps the other features the lines find:
- telomeres;
- knobs and other interstitial satellite blocks;
- subtelomeric blocks;
- rDNA (NORs);
- contact hubs;
- plant-specific methylation candidates.

All on one coordinate system. When the data come from several plants, it also tells which inputs come from the
assembled plant, and scores every other plant chromosome by chromosome.

## Why

Satellite-based centromere finders exist (cenrep, CentIER, quarTeT, TRASH), but a satellite array is not a centromere.
The functional centromere can sit beside the array, far from it, or on no satellite at all (neocentromeres,
satellite-free centromeres in potato, horse chr11; the *Poa annua* Pa1A centromere 21.8 Mb from its satellite). Each
line fails in its own way:

| line | fails when |
|---|---|
| satellite arrays | the array is pericentromeric, inactive or a knob; the centromere has no satellite |
| LTR retrotransposons (CRM clade) | the CRM cluster flanks rather than marks the centromere (some Solanaceae) |
| Hi-C centromere clustering | another trans cluster dominates the nucleus. cenhic's prior-free calls: rice 0/12 (all at chromosome ends), gorilla 0/24, giraffe 0/15; its own confidence did not predict accuracy |
| methylation | it cannot find centromeres de novo: false dips are as deep as true ones |

The lines fail independently, so they can check each other. Two ideas carry the integration:

1. **Hi-C tests the other lines instead of competing with them.**
   - For each chromosome, cenfuse asks how strongly every place contacts the other chromosomes' centromeres, as the
     other chromosomes' sequence calls place them. A chromosome's own call is never used for its own test.
   - Which trans cluster dominates the nucleus then no longer matters.
   - In rice, where cenhic's prior-free calls are all wrong, this anchored test lands on the CENH3-verified centromere
     on 8/12 chromosomes.
2. **How far a line can be trusted is measured in each genome, not assumed.**
   - Hi-C is allowed to move a call only when, on the other chromosomes, its anchored peaks agree with the sequence
     calls more often than chance (majority, P <= 0.01).
   - Methylation is calibrated on the integrated calls of the other chromosomes.

## What it does

1. **Satellite arrays and LTR retrotransposons.** cenrep 0.9.0 (vendored): satellite families de novo, every array
   with its family's role, and one sequence call per chromosome. With `tes/`, the CRM retrotransposon cores also
   place calls.
2. **Hi-C.**
   - cenhic's prior-free call, unchanged and reported.
   - The anchored test, with cenhic's bootstrap, leave-one-partner-out and null tests.
   - The gate.
   - The 3D self-association of every feature class: do centromeres, chromosome ends, telomere arrays, each
     satellite family's blocks, each LTR clade's cores contact each other across chromosomes more than random places
     do?
3. **Methylation (ONT, CHG required).** cenmeth's composition-adjusted track. The centromere signature (which
   sequence classes gain or lose which methylation context at centromeres) is learned from the other chromosomes'
   integrated calls. Each place is scored against it, using only the samples that may decide the reference calls
   (5).
4. **Integration.**
   - The integrated call is the sequence call unless validated Hi-C places a strong peak elsewhere.
   - A high-confidence sequence call moves only when a second line (the CRM core or methylation) supports the Hi-C
     peak and not the sequence call. Without one, the call stays and is marked `contested`.
   - Methylation never moves a call on its own.
5. **Which plant each input comes from, and each plant on its own** (`docs/plants.md`). Data for one assembly
   often come from several sources, and centromeres can sit elsewhere in other plants.
   - **Provenance.** Every methylation sample is genotyped against the assembly. A site whose reads all show
     another base on both strands is a fixed difference; read errors strike one strand. A difference that one
     sample has and another lacks cannot be an assembly error, which every plant's reads share. Every input is
     also checked for long segments it lacks while another input covers them. Verdicts: `reference` (the
     assembled genotype), `other`, or `unknown` (always with one sample, or for Hi-C given as pairs).
   - **Only reference-plant inputs and inputs of unknown origin decide the reference calls.** An input shown to
     come from another plant never does.
   - **Every plant is scored on its own**, chromosome by chromosome. Is its centromere at the reference call, or
     elsewhere? It is compared with the other plants at the same places, against their spread at the other
     chromosomes' calls. A plant's own Hi-C library is tested the same way.
6. **Telomeres.** The genome's telomere unit is found de novo: the 5-12-bp unit enriched in tandem arrays at the
   chromosome ends, such as TTTAGGG or TTAGGG. Each end then gets its array.
7. **Loci.** Every candidate feature of every line is merged into loci. Each locus lists its members, every line's
   value there, and a kind:
   - centromere; former (sequence) centromere; alternative centromere candidate;
   - telomere; rDNA; subtelomeric block; interstitial block (knob-like);
   - Hi-C hub; Hi-C peak; LTR core; satellite block.

   A per-kind table gives each kind's joint signature across the four lines. Only the centromere call combines the
   lines by rule. Every other kind is labeled by the line that finds it (mostly cenrep's satellite roles), and the
   other lines are measured there and reported, not used to decide the label.

## Install

```bash
mamba env create -f environment.yml && mamba activate cenfuse
pip install --no-deps .        # the `cenfuse` command
pytest -q -n 8                 # 432 tests (cenfuse's own plus the three engines' suites), ~1-2 min
```

## Run

```bash
cenfuse run GENOME_DIR -t 32                 # outputs in GENOME_DIR/cenfuse/; a rerun resumes
cenfuse run GENOME_DIR --satellite-dir old/centromere   # reuse existing cenrep outputs (no satellite engine run)
```

`GENOME_DIR/` holds one FASTA per genome. Genomes analyzed together share one satellite-family namespace, as in
cenrep.

| file | needed | content |
|---|---|---|
| `<g>.fa` / `.fasta` / `.fna` (`.gz` bgzip) + `.fai` | yes | chromosome-level assembly; unplaced scaffolds are dropped by a length-gap rule and reported |
| `genes/<g>.gff3` | no | `gene` lines: subtelomere and gene-density tests, gene class for methylation |
| `tes/<g>.gff3` | no (needed for methylation) | LTR retrotransposons with `clade=` (LTRquest, or EDTA + TEsorter) |
| `hic/<g>.pairs.gz` (4DN or 6-column) or `.cool`/`.mcool` | no | uniquely mapped, deduplicated Hi-C pairs on this assembly |
| `hic/<g>.tsv` | no | several Hi-C libraries: `sample`, `file`, one row per library. A library and a methylation sample with the same name are the same plant |
| `meth/<g>.tsv` (+ `meth/<g>.exclude.bed`) | no | samples table: `sample`, `CHG` [, `CG`, `CHH`] = modkit bedMethyl (bgzip + tabix); one row per plant |

Options: `-o` output dir, `-t` threads, `--satellite-dir`, `--seed-engine fastan|trf`, `--seed-cache`,
`--classify-from T=S` (a contig assembly takes its satellite families from S), `--force STAGE`, `-v`. Every
threshold is an internal constant with its reason: `cenfuse/params.py` and the engines' `params.py`.

## Outputs (`OUT/`)

| file | content |
|---|---|
| `<g>.centromeres.tsv` / `.gff3` / `.bed` | one row per chromosome: integrated call, relation (`concordant`, `offset`, `contested`, `functional_only`, `sequence_only`, `none`), confidence, supporting and opposing lines, flags; the sequence call, TE core, anchored and prior-free Hi-C peaks with the gate, methylation at each place |
| `<g>.loci.tsv`, `<g>.features.gff3` | every locus with its kind, members and the value of every line (satellite bp by role, telomere bp, LTR composition, anchored / prior-free / telomere-pole Hi-C score, hub, methylation deviation and raw CHG, gene density, gaps) |
| `<g>.feature_signatures.tsv` | per locus kind: the median signal of every line and the kind's 3D self-association |
| `<g>.feature_contacts.tsv` | per feature class: trans contact among its members vs 1,000 random placements (fold, P) |
| `<g>.provenance.tsv` | one row per input (methylation sample, Hi-C library): `reference`, `other` or `unknown`, the reason, fixed differences per kb, the strongest private-difference test, missing segments, whether it decides the reference calls |
| `<g>.plants.tsv` | one row per plant x chromosome: `shared`, `shifted`, `possible_shift`, `call_absent`, `conflict`, `unclear` or `no_data`; methylation at the call and at the best site away from it, with q values; the plant's Hi-C verdict |
| `<g>.plant_sites.tsv` | every plant scored at every site of every chromosome: label, signature score, between-plant difference, t, P, q |
| `<g>.hic_anchored.bedgraph` | anchored Hi-C score per bin (x chromosome median) |
| `cenfuse_summary.tsv` | one row per genome: what each line found and whether it was validated |
| `cenfuse_overview.pdf`, `<g>.overview.png`, `<g>.plants.png` | per genome: every chromosome with all lines, the evidence matrix, offsets, and the Hi-C score around the calls; with plants, a second page: every plant x chromosome verdict, and the chromosome where the plants differ most drawn plant by plant; short and long legends |
| `satellite/`, `hic/`, `methylation/` | the engines' own outputs (cenrep, cenhic, cenmeth formats) |

## How to read a call

- `concordant`: a validated functional line (Hi-C or methylation) agrees with the satellite or TE call.
- `offset`: validated Hi-C, with a second line where needed, places the centromere away from the sequence call. The
  sequence call stays in the loci table as `sequence_centromere`.
- `contested`: validated Hi-C places a strong peak elsewhere, but nothing else supports it. The call stays one level
  lower.
- `functional_only`: no satellite or TE call; validated Hi-C (with support) places it.
- `sequence_only`: no functional line could be used (no data, gate failed, or Hi-C weak or at a hub, NOR or end).
- Confidence uses cenrep's levels. Each agreeing independent line raises a kept call one level (at most two); a
  contest lowers it one.

## Validation

- **Tests.** 432 pytest tests:
  - cenfuse's own: the decision table, the gate, loci and telomere units, provenance and per-plant scoring, and a
    synthetic genome at the engines' real scales with every planted case (concordant, moved by TE, moved by
    methylation, contested, satellite-free; a reference plant, another plant missing a segment, another plant with a
    shifted centromere and its own Hi-C library) plus a permuted-Hi-C control;
  - the three engines' suites.
- **Benchmark.** Criteria were fixed before any run, within 6 h and 120 threads; full numbers in
  `docs/benchmark.md`. Truth: CENH3/CENP-A for rice, potato, horse and CHM13; satellome for the apes and giraffe; earlier
  analyses for Poa.
  - **No harm.** In every genome the integrated calls are right at least as often as the sequence calls: rice 12/12,
    CHM13 23/24, chimpanzee 25/25, siamang 22/26, gorilla 25/25, potato 2/11, horse 3/4, barley 0/7.
  - **The gate is calibrated.**
    - It passed where anchored Hi-C was right on 67-85% of chromosomes: CHM13, chimpanzee, siamang, rice, Poa.
    - It failed where anchored Hi-C was right on 0-33%: both orangutans, gorilla, giraffe, potato, barley.
    - It also failed on permuted Hi-C, and on horse, whose sequence calls were all low.
    - Hi-C therefore changed no call in a genome where it was misleading.
  - **Anchoring helps in some genomes, not in general.**
    - In rice, cenhic's prior-free calls are 0/12 (all at chromosome ends), while the anchored test finds the CENH3
      centromere on 8/12.
    - Pooled over six mammals the two tie, 73 vs 73 of 138 (pre-registered criterion M2: fail).
  - **A real offset recovered.** *P. annua* Pa1A:
    - From satellite-only calls (146.1 Mb, high), validated Hi-C (peak 168.0 Mb) with the CRM core as second line
      moves the call 20.2 Mb to 167.0-170.0 Mb, the functional centromere of the earlier analysis.
    - With `tes/` given, cenrep's CRM rule already places it there, and Hi-C and methylation confirm it.
  - **Failures.**
    - Potato: without TE input the satellite line is wrong on 9/11, and Hi-C has no centromere signal.
    - Barley: without TE input the satellite line is wrong on 7/7, and cenfuse keeps those calls, although
      cenhic's prior-free Hi-C is right on 6/7. It flags the genome-wide conflict.
    - Horse chr11: satellite-free, missed by every line.
    - CHM13 chr14 and chr22: acrocentric p-arm (NOR) Hi-C peaks contest correct calls.
    - The pre-registered sample-table criterion missed one plant-specific candidate that lies inside a wide call.
- **Provenance and plants** (*P. annua*, five resequenced ONT plants; truth from the earlier analyses;
  `docs/benchmark.md`).
  - **Provenance is right.** All five plants are `other`, none the reference plant: the earlier work showed each
    carries few of the reference's young insertions. The reference plant's own Hi-C stays `unknown`, as a pairs
    file allows, and decides.
  - **The karyotype test finds the known exchanges once fixed.**
    - Pa1A 171-255 Mb absent in C11, J10 and R1, to within 1 Mb.
    - Pa2A in P3 and R1, in the three pieces the earlier work described.
    - So P3 and R1 lack the Pa2A call region (`call_absent`).
    - The first, pre-registered run found none of these events: its skip rule threw out half the genome. That run's
      failures stay on record.
  - **No false shift** on the 13 stable chromosomes. A low-depth plant got two `possible_shift` verdicts.
  - **Per-plant methylation does not resolve Pa1A.** Four plants have insertion evidence that their centromere lies
    off the reference call. Three of them look centromeric at the call anyway, and the fourth is unclear. Some places
    that are not centromeres look centromeric in every plant.

## Limitations

- **Per-plant methylation verdicts are weak evidence.**
  - A plant shift is called only with between-plant support, and none is called falsely.
  - A plant that looks centromeric at the call is only consistent with it: some non-centromeric places look
    centromeric in every plant.
  - Plants at 2-4x often lack data at centromeres.
  - Stronger per-plant lines would be each plant's own Hi-C, CENH3 ChIP, or plant-private young-CRM insertions
    (from reads).
- **Provenance.**
  - Hi-C given as pairs or cool cannot be genotyped. Only missing segments can show it comes from another plant, so
    it is usually `unknown` and decides, which is right when it is the assembled plant's scaffolding library.
  - With one methylation sample, its genotype is `unknown`: assembly errors and another plant cannot be told apart.
  - The genotype noise model is untested on two samples of one plant (e.g. two chemistries).
- Hi-C is validated against the sequence line. When the sequence line is weak or wrong everywhere, the gate cannot pass
  and Hi-C only reports, even when it is right. Barley without TE input: satellite calls 0/7, prior-free Hi-C 6/7,
  integrated 0/7. `cenfuse_summary.tsv` then flags `hic_sequence_conflict`; give `tes/` (CRM cores) or methylation as
  a third line. The flag says that one line is wrong, not which one. It fired on 6 of the 13 real benchmark genomes:
  the sequence line was wrong in barley and potato (potato: both lines), Hi-C was wrong in gorilla, Sumatran orangutan
  and giraffe, and in horse the truth covers too few chromosomes to tell.
- A wrong Hi-C peak on one chromosome, in a genome where Hi-C mostly works, passes the gate. The second-line rule
  protects high sequence calls only.
- Resolution is the Hi-C bin (250 kb-2 Mb); offsets smaller than a bin are not resolved.
- Methylation needs CHG (ONT, modkit) and `tes/`. Wrong anchors weaken the signature, and plants below ~4x often
  lack signature data at centromeres.
- The 3D self-association null moves each feature along its own chromosome. It does not match compartments, so
  heterochromatic classes can look self-associated through the B compartment alone.
- Telomeres are reported only where the assembly holds them. Most non-T2T assemblies lack them at most ends.

## Layout

| path | content |
|---|---|
| `cenfuse/cli.py`, `pipeline.py` | command, stages, checkpoints |
| `cenfuse/inputs.py`, `seqcalls.py`, `teline.py` | genome discovery and chromosome set; the sequence and TE lines read back |
| `cenfuse/hicfuse.py` | anchored test, gate, feature-class 3D association |
| `cenfuse/methfuse.py` | leave-one-chromosome-out calibration, support, candidates |
| `cenfuse/fuse.py` | integration rules and confidence |
| `cenfuse/provenance.py`, `plants.py` | which plant each input comes from; each plant scored on its own |
| `cenfuse/loci.py`, `telomere.py` | loci, kinds and signatures; telomere arrays |
| `cenfuse/output.py`, `plot.py` | files and the overview figure |
| `cenfuse/sat`, `hic`, `meth` | vendored cenrep 0.9.0, cenhic 0.1.1, cenmeth 0.1.0 (`VENDORED.md`) |
| `docs/design.md`, `docs/plants.md`, `docs/benchmark.md` | design and its reasons (plants and provenance in their own note); benchmark criteria and results |

## License

MIT (see `LICENSE`).
