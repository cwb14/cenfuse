# Benchmarks

cenfuse was scored on 55 assemblies (56 genomes and libraries), from T2T mammals to fungi, against centromeres mapped
functionally on the same assembly. Every number on this page comes from one release rerun of cenfuse 0.10.1 on every
benchmark arm (204 runs, 2026-11-04), scored with the rules fixed before each benchmark was first run. Each run read
the stored satellite search of its arm (`--satellite-dir`; the satellite rules are unchanged since), so the column
naming an array that a CRM core replaced (`seq_overruled_*`) is empty in these runs. The overview pages of every
genome are in the [gallery](gallery.md); how each rule works is in [design.md](design.md).

## How calls are scored

- **Truth.** CENP-A CDRs (Logsdon et al. 2024: human CHM13 and five great apes) or CENP-A ChIP (crab-eating macaque,
  pig, horse; pig and horse from other animals); CENH3 ChIP-seq of the assembled line (rice, potato, barley, maize,
  *Brachypodium*, *Ae. tauschii*, soybean, *B. rapa*, tomato, moss, *Chlorella*) or of a related line (sorghum: line V9
  on the BTx623 assembly; *Arabidopsis*: Col-0 on Col-CC), CUT&Tag (*Chlamydomonas*), the published centromeres of two
  *B. rapa* accessions whose ChIP is too weak, the CenSat active arrays of HG002, the authors' CENP-A windows of
  chicken, the CENP-A minor-satellite array of mouse, and CEN/CENP-A/CenH3 coordinates of the fungi and *Plasmodium*.
  ChIP domains were called with one frozen rule (100-kb bins, ChIP/input, half-maximum run, ≥ 4×).
- **Right.** Mammals: the call overlaps the CENP-A domain. Plants: within 4 Hi-C bins of the CENH3 domain (1 Mb
  without Hi-C; barley 8 Mb, its Hi-C bin being 2 Mb). Moss, *Chlamydomonas*, *Chlorella*: within 100 kb. Chicken and
  mouse: overlap. Fungi and *Plasmodium*: within 1 Hi-C bin (10-50 kb). Holocentric genomes have no single centromere,
  so any confident call there is wrong.
- **Confident** = `high` or `medium`; `low` = cenfuse says it is unsure.
- **Chance.** A chromosome's chance of a right call by luck is (call + truth + 2 × tolerance) / chromosome length; P
  is the Poisson-binomial tail of the right count.

## Figures

- [Mammals and Hi-C](img/benchmark_overview.pdf) ![benchmark overview](img/benchmark_overview.png)
  *9 T2T mammals with CENP-A truth; rice, barley, potato and horse EquCab3 with ChIP truth. **A** every chromosome
  with a mapped centromere: the call (box, colored by outcome) against CENP-A/CENH3 (black bar); arrows = barley calls
  moved by Hi-C. **B-F** one chromosome per kind of outcome. **G** every chromosome as one square, FASTA only vs all
  inputs. **H** each line alone vs the integrated call. **I** the Hi-C tests on real vs randomized Hi-C.*
- [Plants and algae](img/benchmark_plants.pdf) ![plant benchmark](img/benchmark_plants.png)
  *14 plant and algal species with CENH3 truth, FASTA only vs + `tes/` (+ Hi-C), and 7 genomes without a single
  centromere. **A** every chromosome run with `tes/`: CENH3, the calls, every element of the centromeric LTR lineage.
  **B, C** *Ae. tauschii* Chr1 placed by its CRM core. **D** the failure: tomato chr10. **E** every conflict between a
  CRM core and a satellite call. **F** *Chlorella*'s Osser insertions. **G** every chromosome as a square.*
- [Kinetochores from CpG methylation](img/benchmark_cpg.pdf) ![CpG kinetochores](img/benchmark_cpg.png)
  *Primates with DiMeLo-seq CENP-A, curated CDRs or CENP-A ChIP, and HG002. **A, B** a typical bonobo centromere. **C**
  widths of call, kinetochore and CENP-A domain. **D** every ape kinetochore on its CENP-A domain. **E** HG002 from ONT
  and PacBio HiFi. **F** the CpG test on real vs rotated methylation. **G** reads thinned to 1/16.*
- [Hi-C in small genomes](img/benchmark_hic.pdf) ![Hi-C finds centromeres where nuclei cluster them](img/benchmark_hic.png)
  *Where nuclei cluster their centromeres, Hi-C contacts alone find every centromere. **A-D** budding yeast: contacts
  between chromosomes, real and shuffled. **E** every chromosome of four genomes with their centromeres and cenfuse's
  calls. **F** clustered centromeres vs territories. **G** chromosomes whose centromere is found.*
- [The Rabl axis in trans](img/benchmark_hic2.pdf) ![Chromosome arms align by distance from the centromere](img/benchmark_hic2.png)
  *Where a nucleus holds its centromeres at one pole, arms of different chromosomes align by distance from the
  centromere. **A** *C. neoformans* in anaphase vs G1 vs randomized. **B** the Rabl polarity (schematic). **C, D** one
  and all 14 chromosomes. **E** every axis test of the benchmark vs chance. **F** what cenfuse does with it.*

## 1. Every genome

Right / chromosomes with truth (confidently wrong in brackets), for each input set benchmarked; "all inputs" names the
richest set where it differs from the columns before it. Chance = the right calls expected by luck (FASTA only).

| genome | truth | FASTA only | + genes | + `tes/` | + Hi-C | all inputs | chance (FASTA only) |
|---|---|---|---|---|---|---|---|
| **mammals** | | | | | | | |
| human (CHM13) | CENP-A | 23/23 | 23/23 |  |  | 23/23 (genes, Hi-C) | 0.8 |
| chimpanzee | CENP-A | 18/18 | 18/18 |  |  | 18/18 (genes, Hi-C) | 0.5 |
| bonobo | CENP-A | 17/18 | 17/18 |  |  | 17/18 (genes, Hi-C) | 0.3 |
| gorilla | CENP-A | 9/9 | 9/9 |  |  | 9/9 (genes, Hi-C) | 0.3 |
| Sumatran orangutan | CENP-A | 8/8 | 8/8 |  |  | 8/8 (genes, Hi-C) | 0.3 |
| Bornean orangutan | CENP-A | 15/16 (1) | 15/16 (1) |  |  | 15/16 (1) (genes, Hi-C) | 0.7 |
| crab-eating macaque | CENP-A | 16/16 | 16/16 |  |  | 16/16 (genes, Hi-C) | 1.8 |
| pig | CENP-A | 8/17 | 13/17 (3) |  |  | 13/17 (3) (genes, Hi-C) | 1.8 |
| horse (T2T) | CENP-A | 29/32 | 29/32 |  |  | 29/32 (genes, Hi-C) | 4.0 |
| HG002 (diploid human) | CenSat active array | 46/46 |  |  |  | 46/46 (ONT CpG methylation) | 2.6 |
| horse (EquCab3) | CENP-A | 3/5 |  |  | 3/5 |  | 0.1 |
| mouse (T2T) | minor-satellite array | 1/20 |  |  | 1/20 |  | 1.4 |
| **mammals (no CENP-A data)** | | | | | | | |
| siamang | none (satellome locus) | 18/26 | 19/26 (3) |  |  | 20/26 (3) (genes, Hi-C) | 1.0 |
| rhesus macaque | none (satellome locus) | 22/22 | 22/22 |  |  | 22/22 (genes, Hi-C) | 4.2 |
| sheep | none (satellome locus) | 26/28 (1) | 27/28 (1) |  |  | 27/28 (1) (genes, Hi-C) | 3.0 |
| goat | none (satellome locus) | 30/30 | 30/30 |  |  | 30/30 (genes, Hi-C) | 4.2 |
| giraffe | none (satellome locus) | 15/15 | 15/15 |  |  | 15/15 (genes, Hi-C) | 0.3 |
| **birds** | | | | | | | |
| chicken | CENP-A windows | 28/33 |  |  | 28/33 |  | 6.6 |
| **plants** | | | | | | | |
| rice | CENH3 | 12/12 | 12/12 | 12/12 | 12/12 | 12/12 (`tes/`, genes, Hi-C) | 1.2 |
| barley | CENH3 | 0/7 (3) | 0/7 | 4/7 (1) | 6/7 | 6/7 (`tes/`, Hi-C) | 0.2 |
| potato | CENH3 | 2/11 (2) | 4/11 (2) | 6/11 (2) | 2/11 (2) | 6/11 (2) (`tes/`, genes, Hi-C) | 0.6 |
| maize (Mo17) | CENH3 | 9/10 (1) | 9/10 (1) | 9/10 (1) |  | 9/10 (1) (`tes/`, genes) | 0.2 |
| sorghum | CENH3 | 10/10 |  | 10/10 | 10/10 |  | 1.4 |
| *Brachypodium distachyon* | CENH3 | 3/3 |  | 3/3 | 3/3 |  | 0.2 |
| *Aegilops tauschii* | CENH3 | 1/7 (4) |  | 5/7 (2) | 1/7 (4) |  | 0.1 |
| soybean | CENH3 | 20/20 |  | 19/20 (1) | 20/20 |  | 2.2 |
| *Brassica rapa* CC1 | CENH3 | 10/10 |  | 10/10 | 10/10 |  | 2.4 |
| *Brassica rapa* CC2 | CENH3 | 10/10 |  | 10/10 |  |  | 2.4 |
| *Brassica rapa* CC3 | CENH3 | 10/10 |  | 10/10 |  |  | 2.5 |
| *Brassica rapa* CCA03 | CENH3 | 10/10 |  | 10/10 |  |  | 2.5 |
| *Brassica rapa* PC1 | CENH3 | 10/10 |  |  |  |  | 2.4 |
| *Brassica rapa* TC1 | CENH3 | 10/10 |  |  |  |  | 2.6 |
| *Brassica rapa* TU1 | CENH3 (published domains) | 10/10 |  |  |  |  | 2.3 |
| *Brassica rapa* OI1 | CENH3 | 10/10 |  |  |  |  | 2.5 |
| *Brassica rapa* CX1 | CENH3 (published domains) | 9/10 |  |  |  |  | 3.2 |
| *Brassica rapa* MIZ1 | CENH3 | 10/10 |  |  |  |  | 2.5 |
| *Brassica rapa* WTC1 | CENH3 | 10/10 |  |  |  |  | 2.3 |
| *Arabidopsis thaliana* (Col-CC) | CENH3 | 5/5 |  | 5/5 | 5/5 | 5/5 (`tes/`, ONT methylation) | 0.8 |
| tomato | CENH3 | 1/12 (7) |  | 1/12 (7) | 1/12 (7) | 1/12 (7) (`tes/`, Hi-C) | 0.9 |
| moss (*Physcomitrium patens*) | CENH3 | 16/26 |  | 16/26 | 16/26 | 25/26 (`tes/`, Hi-C) | 1.4 |
| **algae** | | | | | | | |
| *Chlamydomonas reinhardtii* | CENH3 | 6/17 |  | 6/17 |  |  | 2.6 |
| *Chlorella sorokiniana* | CENH3 | 0/10 |  | 8/10 |  |  | 0.7 |
| **fungi** | | | | | | | |
| *Saccharomyces cerevisiae* | CEN / CENP-A / CenH3 | 0/16 |  |  | 16/16 |  | 0.0 |
| *Schizosaccharomyces pombe* | CEN / CENP-A / CenH3 | 1/3 |  |  | 2/3 |  | 0.2 |
| *Neurospora crassa* | CEN / CENP-A / CenH3 | 0/7 |  |  | 7/7 |  | 0.0 |
| **fungi (held out)** | | | | | | | |
| *Verticillium dahliae* | CEN / CENP-A / CenH3 | 0/8 |  |  | 8/8 |  | 0.0 |
| *Zymoseptoria tritici* | CEN / CENP-A / CenH3 | 0/21 |  |  | 19/21 |  | 0.0 |
| *Cryptococcus neoformans* (anaphase Hi-C) | CEN / CENP-A / CenH3 | 0/14 |  |  | 10/14 |  | 0.0 |
| *Cryptococcus neoformans* (G1 Hi-C) | CEN / CENP-A / CenH3 | 0/14 |  |  | 0/14 |  | 0.0 |
| *Malassezia sympodialis* | CEN / CENP-A / CenH3 | 0/8 |  |  | 6/8 |  | 0.8 |
| *Candida tropicalis* | CEN / CENP-A / CenH3 | 0/7 |  |  | 0/7 |  | 0.2 |
| **protists** | | | | | | | |
| *Plasmodium falciparum* | CEN / CENP-A / CenH3 | 0/14 |  |  | 14/14 |  | 0.4 |
| **nematodes (holocentric control)** | | | | | | | |
| *Caenorhabditis elegans* | none (holocentric) | 0 confident / 6 |  |  | 0 confident / 6 |  |  |
| *Caenorhabditis nigoni* | none (holocentric) | 0 confident / 6 |  |  | 0 confident / 6 |  |  |
| **insects (holocentric control)** | | | | | | | |
| *Bombyx mori* | none (holocentric) | 0 confident / 29 |  |  | 0 confident / 29 |  |  |
| **plants (holocentric control)** | | | | | | | |
| *Rhynchospora breviuscula* | none (holocentric) | 0 confident / 5 |  |  | 0 confident / 5 |  |  |

## 2. Satellite centromeres

### T2T mammals

| genome | truth | FASTA only | + genes + Hi-C | chromosomes passing the Hi-C test |
|---|---|---|---|---|
| human (CHM13) | CENP-A | 23/23 | 23/23 | 24 |
| chimpanzee | CENP-A | 18/18 | 18/18 | 25 |
| bonobo | CENP-A | 17/18 | 17/18 | 25 |
| gorilla | CENP-A | 9/9 | 9/9 | 25 |
| Sumatran orangutan | CENP-A | 8/8 | 8/8 | 25 |
| Bornean orangutan | CENP-A | 15/16 (1) | 15/16 (1) | 25 |
| crab-eating macaque | CENP-A | 16/16 | 16/16 | 22 |
| pig | CENP-A | 8/17 | 13/17 (3) | 0 |
| horse (T2T) | CENP-A | 29/32 | 29/32 | 0 |
| siamang | none (satellome locus) | 18/26 | 20/26 (3) | 26 |
| rhesus macaque | none (satellome locus) | 22/22 | 22/22 | 0 (support tier on 22) |
| sheep | none (satellome locus) | 26/28 (1) | 27/28 (1) | 0 |
| goat | none (satellome locus) | 30/30 | 30/30 | 0 (support tier on 31) |
| giraffe | none (satellome locus) | 15/15 | 15/15 | 0 |

- **Totals** (9 genomes with CENP-A, 157 chromosomes): FASTA only 143 right, 1 confidently wrong, 13 `low` or no call;
  with genes and Hi-C 148, 4, 5. `high` calls were right 117 of 120 times.
- **Without functional truth** (siamang, rhesus, sheep, goat, giraffe; 121 chromosomes): the calls agree with an
  independent satellite caller (the satellome locus) on 111 from the FASTA, 114 with genes and Hi-C.
- **Pig, FASTA only (8/17):** the satellite families split so that no family dominates genome-wide, and the honesty cap
  makes every call `low`; 5 of them sit 0.6-2.0 Mb from CENP-A at the acrocentric ends. With genes and Hi-C: 13/17.

| every miss (all inputs) | call | why |
|---|---|---|
| pig CM000817.6, CM000822.6, CM000828.6 | wrong, confident | the truth is weak ChIP from other pigs, with no centromeric array within 1 Mb of its maximum (5.4-5.6×); on all three the ChIP's second peak (3.0-4.4×) lies inside the call, and the published karyotype agrees with the call on two (chr11, chr17). Second peaks ≥ 3× also occur on 3 of the 13 right pig calls, so this is evidence, not proof |
| pig CM001155.3 | wrong, `low` | the same |
| Bornean orangutan CM054652.2 | wrong, `high` | CENP-A sits on 4 small alpha-satellite arrays (73 kb in cenfuse's search), 14 Mb from the 3.4-Mb array cluster cenfuse called; Hi-C points weakly to the CENP-A site but cannot move a `high` call. The truth may be haplotype-specific: the CenSat annotation marks the called array as active, and the CDR region was redrawn in 2025 |
| bonobo CM055489.2 | wrong, `low` | the true array is the reported competing cluster |
| horse CM128059.1 (chr11) | no call | satellite-free centromere; horse Hi-C fails its test, so nothing can place it |
| horse CM128050.1, CM128070.1 | wrong, `low` | |

### Plants and other genomes with satellite centromeres

- **Rice** (T2T-NIP): 12/12 `high` from the FASTA; Hi-C passes its test on all 12 and changes nothing.
- ***B. rapa*** (11 T2T accessions): 109/110 from the FASTA; the miss (CX1 A05) is `low`, 2.8 Mb off.
- **Soybean** X26, **sorghum** BTx623, ***Brachypodium*** Bd21 (CENH3 truth on 3 of 5 chromosomes), ***Arabidopsis***
  Col-CC: 20/20, 10/10, 3/3, 5/5 from the FASTA.
- **Maize** Mo17 (T2T): the 180-bp knob repeat is the most abundant satellite family (23.6 Mb, 20.6 Mb on chr8) and
  gets the role `interstitial` (knobs and other blocks); CentC is centromeric on all 10 chromosomes and no call holds
  knob repeat. 9/10 right (7 overlap CENH3). The
  miss, chr5 (`high`), has CENH3 on satellite-free DNA 2 Mb from the CentC array.
- **HG002** (diploid human T2T, 46 chromosomes): 46/46 calls on the CenSat active array (45 `high`).
- **Chicken** GGswu (40 chromosomes; truth on 33): 28/33 on the CENP-A windows (chance 6.6, P 9e-23), 2 more within 1
  Mb; every call `low`, because no satellite family dominates (chicken centromeres use chromosome-specific satellites,
  and chr5, chr27 and Z have none). Confidence is conservative here: no call is confidently wrong.
- **Mouse** (T2T, all 20 chromosomes acrocentric): both the pericentric major satellite (162 Mb) and the CENP-A
  minor satellite (16.5 Mb) mark the centromeric end. All 20 calls sit there, on the major-satellite block, within
  45 kb of the minor satellite (1/20 overlapping it), all `low`, with the minor satellite reported as the alternative
  on 20/20. No mappable DNA separates the two arrays, so Hi-C cannot choose (`hic_alt_ambiguous`).
- **Horse EquCab3** (older assembly, CENP-A ChIP on 5 chromosomes): 3/5 right, all `low`; chr11 has no call.

## 3. Retrotransposon centromeres (`tes/`)

`tes/` = LTRquest annotations. Right / chromosomes (confidently wrong):

| genome | centromeres | FASTA only | + `tes/` | marker lineage |
|---|---|---|---|---|
| rice | CentO satellite + CRM | 12/12 | 12/12 | CRM (decisive) |
| potato | satellite on half, satellite-free on half | 2/11 (2) | 6/11 (2) | CRM (decisive) |
| maize (Mo17) | CentC satellite + CRM | 9/10 (1) | 9/10 (1) | CRM (decisive) |
| barley | CRM (Cereba) retrotransposons | 0/7 (3) | 4/7 (1) | CRM (decisive) |
| sorghum | satellite + CRM | 10/10 | 10/10 | CRM (decisive) |
| *Brachypodium distachyon* | satellite + CRM | 3/3 | 3/3 | CRM (decisive) |
| *Aegilops tauschii* | CRM beside satellites | 1/7 (4) | 5/7 (2) | CRM (decisive) |
| soybean | satellite | 20/20 | 19/20 (1) | CRM (decisive) |
| *Brassica rapa* CC1 | satellite | 10/10 | 10/10 | CRM (decisive) |
| *Brassica rapa* CC2 | satellite | 10/10 | 10/10 | CRM (decisive) |
| *Brassica rapa* CC3 | satellite | 10/10 | 10/10 | CRM (decisive) |
| *Brassica rapa* CCA03 | satellite | 10/10 | 10/10 | CRM (decisive) |
| *Arabidopsis thaliana* (Col-CC) | CEN180 satellite | 5/5 | 5/5 | Athila (flags only) |
| moss (*Physcomitrium patens*) | Copia (Bryco) retrotransposon clusters | 16/26 | 16/26 | Bryco (flags only) |
| *Chlamydomonas reinhardtii* | LINE (ZeppL) clusters | 6/17 | 6/17 | none |
| *Chlorella sorokiniana* | 2-7 Copia (Osser) insertions | 0/10 | 8/10 | Osser (flags only) |
| tomato | satellite-free, beside a dominant satellite | 1/12 (7) | 1/12 (7) | Ale (flags only) |
| moss, + `tes/` + Hi-C | | | 25/26 | Bryco (flags only) |

- **CRM decides** in rice, potato, maize, barley, sorghum, *Brachypodium*, *Ae. tauschii*, soybean and *B. rapa*: its
  cores sit on CENH3 (rice 8/8, potato 6/7, maize 10/10 within 1 Mb).
- **Potato:** `tes/` places the satellite-free centromeres (2 → 6/11 right); genes move chr01 and chr03 by vetoing
  gene-rich satellite families.
- ***Ae. tauschii*:** satellite arrays lie beside CRM cores that hold CENH3 on all 7 chromosomes. The head-to-head
  test (the core holds significantly more CRM than the call) moves 5 calls onto the cores; Chr5 and Chr7 stay at a
  core's edge, 1.3 and 3.1 Mb from CENH3. From the FASTA alone 4 calls are confidently wrong.
- **The head-to-head test protects right satellite calls:** in *B. rapa* A09 and soybean chr1/chr6 the densest window
  of a broad CRM plateau lay off CENH3, which sat on the satellite (16-26 vs 14-32 elements); without the test, `tes/`
  made those calls wrong.
- ***Chlorella*:** centromeres of 2-7 Osser (Copia) insertions (53 of the genome's 54 lie inside the CENH3 domains),
  found from 3 independent insertion sites: 0 → 8/10 (all `low`).
- **Moss:** Bryco (Copia) is the only lineage concentrated at one locus per chromosome; its cores sit on CENH3 on 25/26,
  but a non-CRM lineage only flags. `tes/` alone: 16/26 (11 confident); with Hi-C as well: 25/26, all confident.
- ***Chlamydomonas*:** centromeres are LINE (ZeppL) clusters; no LTR lineage concentrates, so `tes/` changes nothing.
- **Tomato** (a reported failure): 9 of 12 CENH3 domains hold no satellite array (11 hold under 10 kb), no LTR lineage marks them, and the
  dominant 53-bp satellite (the largest cluster on 7 chromosomes, so no honesty cap) and a 181-bp satellite take 7
  confident wrong calls (5 and 2).
  43 M Hi-C pairs give prior-free peaks 1-37 Mb off on all 12.

**The annotation matters.** The decisive CRM test needs CRM to concentrate on most chromosomes. Right calls with each
annotator, and with LTRquest's elements thinned at random (TE line alone: right / cores):

| genome | annotation | TE line alone | calls right (confidently wrong) |
|---|---|---|---|
| rice | LTRquest | 8/8 | 12/12 (0) |
| rice | EDTA intact + TEsorter; DANTE_LTR | 4/4; 3/8 | 12/12 (0) |
| rice | LTRquest at 50 / 25 / 10% | 3/3, 1/1, 1/1 | 12/12 (0) |
| potato | LTRquest | 6/7 | 6/11 (2) |
| potato | EDTA + TEsorter; DANTE_LTR; LTRquest at 50% | 2/2; 1/1; 3/3 | 2/11 (2) |
| maize | LTRquest, and at 50 / 25 / 10% | 10/10, 10/10, 10/10, 8/8 | 9/10 (1) |
| moss | LTRquest; at 50% | 25/25; 22/22 | 16/26 (0) |

Potato's 214 LTRquest CRM elements give exactly a majority of chromosomes; any less sensitive annotation loses the
gain from `tes/`. Every CRM core short of a majority still lay on CENH3 (24/24 in rice, potato, maize and
*Brachypodium*).

## 4. Hi-C

### The test, and Hi-C on its own

Before Hi-C may act on a chromosome, its peaks must agree with the other chromosomes' calls on most of them, beyond
chance ([design.md](design.md#hi-c-cenfusehic-and-hicfusepy)). Each Hi-C readout alone, on the genomes with truth:

| genomes (chromosomes) | satellite line | prior-free Hi-C | anchored Hi-C | integrated |
|---|---|---|---|---|
| 9 T2T mammals (157) | 148 | 32 | 44 | 148 |
| rice (12) | 12 | 0 | 7 | 12 |
| barley (7) | 0 | 6 | 0 | 6 |
| potato (11) | 2 | 1 | 0 | 2 |
| horse EquCab3 (5) | 3 | 3 | 3 | 3 |

Every line scored within 4 Hi-C bins of the truth (Hi-C's resolution); the mammals' satellite line is the one with
genes (FASTA alone: 143).

- In mammals the anchored peak sits beside the unmappable centromeric array, often more than 4 bins from CENP-A in bp;
  the satellite line decides, and Hi-C confirms (the test passes in all great apes, crab-eating macaque, CHM13 and
  siamang).
- Hi-C's prior-free peaks cannot be trusted where they can be tested (the ones the engine rates `high` lie within 4 bins
  of CENP-A on 14 of 63 T2T mammal chromosomes; potato's strongest mutual cluster is not its centromeres), so they
  never act unchecked.
- **Randomized Hi-C** (partners re-paired at random with coverage kept, or ends moved): in 47 runs (14 mammals
  randomized both ways, *P. annua*, 5 small genomes and 13 held-out genomes), no test passed, the support tier was never
  reached, and every call equals the call made without Hi-C.

### The Rabl axis: barley, small genomes and the trans arm alignment

In a Rabl nucleus the arms run from the centromere pole to the telomere pole, so each chromosome's own map folds at
its centromere (the arm fold), and arms of different chromosomes align by distance from their centromeres (the trans
arm alignment). Either can test Hi-C without the sequence calls.

- **Barley:** the satellite calls sit toward the chromosome ends (centromeric index ≤ 0.14; CENH3, on Cereba
  retrotransposons, at 0.35-0.48). The arm fold agrees with Hi-C's prior-free peaks on 6/7 chromosomes (P 3e-5), so 6
  calls move, all onto CENH3 (0/7 → 6/7). On 5H the prior-free peak sits on CENH3 but outside the fold core: the call
  is not moved and drops from `high` to `low`.
- Chromosomes whose prior-free peak lies in the alignment core, in every run of the benchmark (159 runs):

| nuclei | runs | alignment passes | examples (k/n, P) |
|---|---|---|---|
| Type I fungi and *Plasmodium* | 4 | 3 | yeast 13/16 (3e-5), *Plasmodium* 12/14 (4e-8), *Neurospora* 6/7 (3e-4); *S. pombe* 1/3 |
| barley, full and half depth | 3 | 3 | 7/7 (3e-5); half depth 5/7 (0.002) |
| Type I fungi, held out | 5 | 3 | *V. dahliae* 8/8 (2e-5), *C. neoformans* anaphase 12/14 (1.5e-7), *M. sympodialis* 6/8 (4e-4); not *Z. tritici* 6/20 (0.055), *C. tropicalis* 3/6 (0.049) |
| *C. neoformans* in G1 (centromeres unclustered) | 1 | 0 | 0/14 |
| territorial mammals | 23 | 0 | best CHM13 7/23 (P 0.0013: beyond chance, but a minority) |
| other plants, chicken, mouse, horse EquCab3 | 19 | 0 | |
| holocentric | 4 | 0 | *C. elegans* 2/5 (0.095) |
| randomized Hi-C | 46 | 0 | |
| thinned libraries (2-25%) | 44 | 0 | |

- One assembly in two nuclear states: *C. neoformans* in anaphase passes and gives 10/14 right calls; in G1, when its
  kinetochores are unclustered (Polisetty et al. 2025), nothing passes and nothing is called. Placed on the known
  centromeres, the G1 nucleus still aligns its arms on 8/14 chromosomes (P 0.0014): the alignment weakens without the
  cluster rather than vanishing.

### Point, regional and AT-rich centromeres: Hi-C alone

With no centromeric repeat, the FASTA gives (almost) no call. Bins go down to 10-50 kb for small chromosomes. Right =
within 1 Hi-C bin:

| genome | centromeres | Hi-C | bin | FASTA only | + Hi-C |
|---|---|---|---|---|---|
| *Saccharomyces cerevisiae* | point (125 bp) | 14.9 M | 10 kb | 0/16 | 16/16 [0/13/3] |
| *Plasmodium falciparum* | AT-rich core | 15.9 M | 10 kb | 0/14 | 14/14 [0/12/2] |
| *Schizosaccharomyces pombe* | regional | 22.1 M | 50 kb | 1/3 | 2/3 [0/0/2] |
| *Neurospora crassa* | regional, RIP-mutated TE relics | 29.5 M | 50 kb | 0/7 | 7/7 [0/6/1] |
| *Verticillium dahliae* | regional | 17.7 M | 50 kb | 0/8 | 8/8 [0/7/1] |
| *Zymoseptoria tritici* | regional | 23.8 M | 10 kb | 0/21 | 19/21 [0/0/19] |
| *Cryptococcus neoformans* (anaphase Hi-C) | regional, retrotransposon-rich | 20.3 M | 10 kb | 0/14 | 10/14 [0/0/10] |
| *Cryptococcus neoformans* (G1 Hi-C) | regional, retrotransposon-rich | 20.0 M | 10 kb | 0/14 | 0/14 |
| *Malassezia sympodialis* | short regional | 15.0 M | 10 kb | 0/8 | 6/8 [0/0/6] |
| *Candida tropicalis* | regional (3C-seq) | 11.4 M | 50 kb | 0/7 | 0/7 |

- Yeast and *Plasmodium*: the trans arm alignment passes and the arm fold leans the same way, so the calls are
  `medium`; *Neurospora*: the arm fold itself passes. Where neither passes, Hi-C's clean prior-free peaks give `low`
  calls (*Z. tritici*, *C. neoformans*, *M. sympodialis*).
- *S. pombe* chr2: the prior-free peak is only `medium`, so no call. *C. tropicalis* (3C-seq) passes nothing.
- Randomized Hi-C on each of these genomes: no test passes and no Hi-C-alone call is made.

### Held out: crop and mouse genomes newly given Hi-C

| genome | Hi-C (pairs kept) | FASTA only | + Hi-C | what the 3D data showed |
|---|---|---|---|---|
| *Brassica rapa* CC1 | 8.2 M | 10/10 [9/0/1] | 10/10 [9/0/1] | centromere_clustering; test passes on 10 |
| soybean | 10.7 M | 20/20 [18/0/2] | 20/20 [18/1/1] | centromere_clustering; test passes on 20 |
| sorghum | 3.9 M | 10/10 [10/0/0] | 10/10 [10/0/0] | not_detected; test passes on 0 |
| *Brachypodium distachyon* | 14.3 M | 3/3 [3/0/0] | 3/3 [3/0/0] | not_detected; test passes on 0 |
| *Aegilops tauschii* | 0.8 M | 1/7 (4) [0/0/1] | 1/7 (4) [0/0/1] | not_detected; test passes on 0 |
| mouse (T2T) | 15.7 M | 1/20 [0/0/1] | 1/20 [0/0/1] | weak_centromere_clustering; test passes on 0 |

No right call became wrong and no confident error was added in any of them. Three libraries are shallower than
downloaded (sorghum, *Ae. tauschii*, mouse: a mate file was damaged in transfer and only its intact part was mapped).

### Weak Hi-C: the support tier

Where Hi-C finds the centromeres beyond chance but not on most chromosomes (territorial nuclei: goat, rhesus, horse
EquCab3, moss without `tes/`), it may only support the calls it agrees with. At full depth the tier supported 23 calls in three genomes, all of them on CENH3 (moss without `tes/`, 10 calls raised to `medium`) or on the satellome locus (goat 8 raised to `high`, rhesus 5 already `high`). Randomized Hi-C never reached the tier.

### Depth, MAPQ, other pipelines' files, replicates

- **Depth.** Pairs thinned at random (nested; 2-50%). Chromosomes on which the Hi-C test passes:

| genome | pairs | 100% | 50% | 25% | 10% | 5% | 2% |
|---|---|---|---|---|---|---|---|
| CHM13 | 71 M | 24 | 24 | 24 | 24 | 24 | 0 |
| chimpanzee | 74 M | 25 | 25 | 25 | 25 | 25 | 0 |
| gorilla | 55 M | 25 | 25 | 25 | 0 | 0 | 0 |
| Sumatran; Bornean orangutan | 71 M; 67 M | 25; 25 | 0; 0 | 0; 0 | 0; 0 | 0; 0 | 0; 0 |
| crab-eating macaque | 70 M | 22 | 22 | 22 | 22 | 22 | 22 |
| pig | 54 M | 0 | 0 | 20 | 20 | 0 | 0 |
| horse (T2T); potato | 64 M; 66 M | 0 | 0 | 0 | 0 | 0 | 0 |
| rice | 53 M | 12 | 12 | 0 | 0 | 0 | 0 |
| barley (arm fold, or the trans alignment) | 9.2 M | pass | pass (alignment) | fail | fail | fail | fail |

  Tests flip at their threshold rather than fading: both orangutans fail at half depth (34-35 M pairs). Where a test
  fails, the calls are exactly those made without Hi-C (32 of 32 such runs), so thinner Hi-C never made a call worse
  than the FASTA alone; barley falls back to its FASTA-only calls (3 confidently wrong) below half depth. Give at
  least ~70 M pairs to a 3-Gb genome.
- **MAPQ floor.** The same reads kept at pairtools' default floor (MAPQ ≥ 1) instead of 30 changed no test and made no
  call confidently wrong in CHM13, orangutan, rice, moss or barley; near-identical repeat copies map at MAPQ 0, which
  pairtools drops at any floor.
- **Other pipelines.** Rice's 53 M pairs rewritten as HiC-Pro `allValidPairs` and Juicer `merged_nodups`, and its
  `.hic` converted with `hictk`: identical calls and tables.
- **Replicates.** Rice's pairs split into two files given as replicates of one plant: every table identical to the
  single file.
- **A territorial plant nucleus.** *Arabidopsis* Col-0 leaf Hi-C (15.4 M pairs): the five calls stay right and `high`;
  the prior-free peaks lie 1-4 Mb off CENH3 (chromocenter neighbors), as expected where centromeres do not dominate
  trans contacts.

## 5. Kinetochores from CpG methylation

Inputs: public ONT CpG tracks of the assembled individuals (apes: modkit 5mC bigWigs; CHM13 nanopolish; crab-eating
macaque with read counts; rhesus R10 modkit bedMethyl), with genes and Hi-C. Truth: DiMeLo-seq CENP-A (5 apes), the
curated CDRs (the same ONT data: an implementation check), and CENP-A ChIP of another crab-eating macaque. The line was
designed on these primate tracks (not blind); rhesus, *Arabidopsis* and HG002 were held out.

| test | result |
|---|---|
| the CpG test passes on real tracks | 6 of 7 primates with truth (CHM13, bonobo, gorilla, both orangutans, crab-eating macaque) and rhesus; not chimpanzee (positional dips at 1 of 25 calls) or siamang |
| it fails on misplaced methylation | rotated tracks 0/8 genomes; tracks shifted by half a chromosome rejected by the assembly check on 6/8 |
| kinetochore on DiMeLo-seq CENP-A (5 apes) | 35 of 39 within 50 kb (29 overlap); CENP-A signal higher inside the kinetochore than over the call on 39/39 |
| kinetochore on the curated CDR | 62/63 |
| crab-eating macaque kinetochore inside the ChIP domain | 12/15 (the other 3 lie 0.8-3.4 Mb from the domain of another animal's ChIP, inside the same call) |
| resolution | kinetochores a median 55 kb wide, calls 5.1 Mb (133 kinetochores, 6 primates) |
| calls changed | 0 |

- **HG002** (diploid human; no CENP-A data): CpG from ONT and from PacBio HiFi. The test passes on 46/46 with both; kinetochores on 40 (ONT, median 38 kb) and 43 (PacBio, median 25 kb), all inside the active array; where both place one (39 chromosomes), they lie within 50 kb on 31. PacBio tracks enter as bedGraphs, with no code change.
- **Low coverage.** Calls are pooled per 5-kb bin, so rhesus thinned to 1/16 of its reads (about one call per CpG
  strand) keeps 21/21 kinetochores within 50 kb of full depth; crab-eating macaque the same.
- **Wrong calls get kinetochores too**: Bornean orangutan CM054652.2 gets one inside its inactive array. The line
  refines a call; it cannot check it.
- ***Arabidopsis*:** the CpG test is right to fail (positional dips on 2/5; the strongest dips lie off CENH3 on 4/5).

## 6. CHG methylation (plants)

- ***Arabidopsis*** Col-CC, same-genotype ONT (~56×): calls identical with and without methylation; the CHG line finds
  no signature, because its expected domain size comes from genome size (315 kb) while CEN180 arrays span 2.2-3.5 Mb.
- ***Poa annua*** (all four inputs; five ONT plants, all other plants): methylation does not decide, since no sample is
  from the assembled plant.

## 7. Controls with no single centromere

| genome | why it is hard | calls (FASTA only; with 3D data) |
|---|---|---|
| *C. elegans* CGC1 | holocentric, no centromeric satellite | 6 `low`; with N2 Hi-C, the same |
| *C. nigoni* | holocentric | 6 `low`; with Hi-C, the same |
| *Bombyx mori* | holocentric; Lepidoptera lost CenH3 | 13 `low`, 16 none; with Arima Hi-C, the same |
| *Rhynchospora breviuscula* | holocentric; Tyba arrays along every chromosome | 5 `low`; with Omni-C, the same |

0 confident calls in 46 chromosomes, from the FASTA and with 3D data. Used without its test, Hi-C's own peaks would
have made 21 confident calls on these chromosomes (*C. elegans* 5, *C. nigoni* 5, *Bombyx* 10, *Rhynchospora* 1); every
Hi-C test failed, so nothing changed.

## 8. Several plants (*Poa annua*)

*P. annua* (14 chromosomes) and *P. infirma* (7), all four inputs; truth from earlier analyses (no ChIP).

- Pa1A is called on its functional centromere (163.0-169.4 Mb, on the CRM core; Hi-C peak 168.0 Mb, `high`), 21 Mb
  from its largest satellite; starting from satellite-only calls, Hi-C and the CRM core move it there.
- The Hi-C test passes on 14/14 and 7/7; randomized *P. annua* Hi-C fails 0/14 and changes nothing.
- The five methylation plants are recognized as other plants (5/5), and the five known missing segments are found
  (the two Pa2A losses in 3 pieces each). On 300 pure-noise genomes the between-plant test gave 5.0% of null tests at
  P < 0.05 and no `shifted` verdict.

## 9. Other features

- **Telomeres.** Arrays at 687 of the 688 ends of 14 T2T mammals that hold one (≥ 500 bp of TTAGGG in the terminal 10
  kb). Short telomeres count when they reach the end: *Neurospora* 13/14 ends (157-207 bp), *P. falciparum* 27/28,
  yeast's irregular TG1-3 repeats found as GGGTGT. Misses: ends without telomere repeat in the assembly (chicken 22),
  TTAGG broken by retrotransposons (*Bombyx* 15).
- **Interstitial telomeric arrays.** CHM13: 16 arrays ≥ 300 bp, one at the 2q13 fusion site (chr2 114.03 Mb).
  Chicken: 50 arrays, 22 of them on or within 100 kb of a CENP-A window (expected 1.4; binomial P 4e-21), holding
  63% of all interstitial TTAGGG: the telomere-repeat centromere cores described by Huang et al. 2023, genome-wide.
- **45S rDNA (NORs).** Every checkable array lies on a known NOR chromosome: human chr13/14/15/21/22, *Arabidopsis*
  chr2/chr4, chicken chr16, *C. elegans* chrI, pig chr8/chr10, horse chr1/28/31, barley 1H/2H/5H/6H/7H, mouse
  chr12/15/16/18/19.
- **5S rDNA.** Every array lies on a chromosome known to carry 5S rDNA: CHM13 1q42 (126 copies of 2.24 kb), the HSA1
  homolog in every ape and macaque, mouse chr8, chicken chr9, *Arabidopsis* chr3/4/5, rice chr11, maize 2L, potato
  chr01, barley 2H/4H/7H; none in *S. pombe*, whose 36 5S genes are dispersed.
- **Morphology.** Acrocentric vs bi-armed agrees with the literature karyotype on 117 of 118 right calls (mouse 20/20,
  goat 30/30, horse 29/29, pig 13/13, sheep 25/26). By sequence, CHM13's acrocentric chr21 and chr22 come out
  submetacentric/subtelocentric: their T2T short arms, rDNA included, are 11-15 Mb.
- **Completeness.** Assembly gaps equal an independent scan on EquCab3, barley and potato (4,084 gaps).
- **GFF3.** Every features file passes `gt gff3validator -typecheck so`.

## 10. Speed, memory and reproducibility

- **Runtime.** The satellite search takes most of the time and grows with satellite content: T2T mammals (2.6-3.6 Gb)
  1-10 h at 12 threads (median 4.3 h), mouse (~6% satellite) 9 h at 16, maize 50 min at 16, *Ae. tauschii* (4.1 Gb)
  2 h at 16, *Arabidopsis* 15 min at 8. Hi-C and the integration take minutes.
- **Memory.** Peak of a whole run: rhesus 74 GB, horse 59 GB, giraffe 18 GB at 12 threads; chicken 42 GB at 16;
  *Arabidopsis* 5 GB at 8.
- **Reproducibility.** Two fresh moss runs gave 307 identical files apart from header dates; Col-CC reruns at 8 and 16
  threads are identical.

## 11. Every confident miss (all inputs)

| genome (inputs) | chromosome | confidence | distance from truth | the call's caveat |
|---|---|---|---|---|
| Bornean orangutan (genes, Hi-C) | CM054652.2 | high | 14.3 Mb | Hi-C's peak lies elsewhere (25.75-26.00 Mb) but is not strong enough to oppose |
| pig (genes, Hi-C) | CM000817.6 | high | 73.3 Mb | none given |
| pig (genes, Hi-C) | CM000822.6 | high | 37.4 Mb | none given |
| pig (genes, Hi-C) | CM000828.6 | medium | 58.2 Mb | near a chromosome end; Hi-C failed its test here |
| potato (`tes/`, genes, Hi-C) | chr04 | high | 3.4 Mb | none given |
| potato (`tes/`, genes, Hi-C) | chr08 | high | 1.7 Mb | none given |
| maize (Mo17) (`tes/`, genes) | CM039154.1 | high | 2.0 Mb | none given |
| *Aegilops tauschii* (`tes/`) | Chr5 | high | 1.3 Mb | none given |
| *Aegilops tauschii* (`tes/`) | Chr7 | medium | 3.1 Mb | another satellite cluster (339.52-339.62 Mb) competes: the sequence cannot tell which is active; little centromeric satellite -> CENP-A/CENH3 data would decide, or Hi-C when mappable DNA separates them |
| soybean (`tes/`) | chr13 | medium | 1.9 Mb | none given |
| tomato (`tes/`, Hi-C) | CM175067.1 | high | 2.3 Mb | the Ale core lies elsewhere (5.83-7.89 Mb) |
| tomato (`tes/`, Hi-C) | CM175070.1 | high | 41.8 Mb | the Ale core lies elsewhere (51.96-53.46 Mb) |
| tomato (`tes/`, Hi-C) | CM175072.1 | high | 3.5 Mb | none given |
| tomato (`tes/`, Hi-C) | CM175074.1 | high | 11.3 Mb | the Ale core lies elsewhere (56.47-58.26 Mb) |
| tomato (`tes/`, Hi-C) | CM175076.1 | high | 4.4 Mb | the Ale core lies elsewhere (3.89-5.37 Mb) |
| tomato (`tes/`, Hi-C) | CM175077.1 | high | 1.9 Mb | none given |
| tomato (`tes/`, Hi-C) | CM175078.1 | high | 17.4 Mb | none given |

In the mammals without CENP-A data, 4 confident calls disagree with the satellome locus (siamang 3, sheep 1); which of the two callers is right there is unknown.

## 12. How the benchmark was run

- **Pre-registration.** Each benchmark's criteria were written before its runs: plants 2026-10-16/19, mammals
  2026-10-21, then one plan per release (kinetochores, hard genomes, new lineages, the LTR line, a from-scratch rerun
  with new truths, Hi-C in small genomes, the Rabl axis and held-out genomes), each with dated amendments. The plans,
  run records and reports are kept beside the runs.
- **Criteria that failed as written** (reported, not reworded):
  - Hi-C peaks within 4 bins of CENP-A wherever the test passes: in four apes the peaks sit beside the unmappable
    array, more than 4 bins from CENP-A.
  - `high` calls right more often than `medium`/`low` at P ≤ 0.01: too few `high` errors to test.
  - Methylation supporting ≥ 10/13 *P. annua* calls: 0/13, since every sample is another plant.
  - The crab-eating macaque kinetochore inside the ChIP domain on ≥ 90%: 12/15; the other 3 lie inside the same call,
    0.8-3.4 Mb from the domain of another animal's ChIP.
  - The *Arabidopsis* CHG line supporting ≥ 3 calls: 0/5.
  - ≥ 15/20 mouse calls overlapping the minor satellite: 1/20.
  - ≥ 90% of the calls the Hi-C support tier raised being right: 88%.
  - The Rabl axis passing in the grasses held out: *Brachypodium* leaf Hi-C and *Ae. tauschii* (0.84 M pairs) passed
    nothing.
- **Rules changed after a result was seen.** Each change has a reason independent of the case that prompted it, and
  was rechecked on cases not used to make it:

| change | prompted by | rechecked on |
|---|---|---|
| Hi-C peaks in end bins never oppose a call; Hi-C opposes only beyond 4 bins | rice chr4; CHM13 | potato, barley, horse, gorilla, mammals |
| satellite rules (tandem test, satellite kinds, acrocentric test, competing clusters); Hi-C in mappable bins | the first mammal run | the plant benchmark, a 185-genome batch |
| the arm fold; fold centers only at centromeric index ≥ 12.5% | barley | rice, all mammals, randomized controls |
| the CpG line, and its positional (rotation) null | designed on the primate tracks; the rotated gorilla control passed the first null | the other rotated controls, rhesus, *Arabidopsis* |
| `low` calls as Hi-C references when no call is confident (clean peaks only) | siamang without genes | randomized Hi-C, chicken, EquCab3 |
| the acrocentric test judged per satellite kind; a second unrelated acrocentric satellite competes | T2T mouse | 72 earlier satellite runs, 9 other genomes |
| a dispersed-repeat family gets role `dispersed` | chromosome-scale loci in bonobo | 220-genome replay: calls identical |
| TE core floor of 3 independent insertion sites | *Chlorella* | 12 crop `tes/` arms, moss, *Chlamydomonas*, *Arabidopsis*: no call changed |
| a CRM core replaces a confident satellite call only if CRM marks it significantly more | *B. rapa* A09 | held-out *B. rapa*, soybean, sorghum, *Ae. tauschii*, tomato |
| the same test for a satellite call beside the core | *Ae. tauschii* Chr2/Chr4 | 21 `tes/` arms (only those 2 changed) |
| families only from seed clusters with a tandem array (speed) | T2T sheep's search time | siamang, barley, Bornean orangutan: calls identical |
| the support tier counts only significant peaks | randomized Sumatran orangutan reached the first form | 28 randomized controls |
| a peak is `telomeric` only when it cannot be told from the end | moss centromeres near the ends | the development replay (34 runs) and 13 held-out genomes |
| the trans arm alignment acts where the fold leans the same way | *C. elegans* nearly passed; barley at half depth could not be rescued | holocentric folds, held-out genomes |
| no Hi-C support where no mappable bin separates a call from its competing cluster | held-out mouse (7 confident calls on the major satellite) | development replay: costs 3 right ape calls (`medium` → `low`) |
