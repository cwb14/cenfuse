# Gallery: cenfuse's overview figure for every benchmark genome

Each image is the first page of `cenfuse_overview.pdf`, unedited, from cenfuse 0.8.0 on the benchmark inputs
(2026-10-30). Each run reads the satellite search of the v8 from-scratch rerun (`--satellite-dir`) and recomputes
everything else; its calls are identical to the v8 run's. *Poa infirma*'s comes from a fresh 0.8.0 run of both *Poa*
genomes with all four inputs, whose *P. annua* page is in the [README](../README.md). Numbers:
[benchmark.md](benchmark.md).

**How to read one** (panel letters as in the images):

| panel | shows |
|---|---|
| **A** | every chromosome to scale: satellite arrays colored by role, the Hi-C centromere score (blue line), with `tes/` the
marker retrotransposon clade's density (orange), the satellite call (black box) and the final call (band; blue band = placed by Hi-C, arrow = moved); telomere arrays at the ends (▶ ◀) and inside the chromosome (▼); 45S (◆) and 5S (◇) rDNA arrays |
| **B** | which evidence lines support (●) or oppose (×) each final call, then its relation and confidence (a `fold` column appears only where the arm fold acted: barley) |
| **C** | where each line puts the centromere, relative to the satellite call |
| **D** | Hi-C score around every final call: a peak at 0 means the calls sit where centromeres meet in the nucleus |

## Mammals with CENP-A truth (right / chromosomes, all inputs)

| genome | right | Hi-C test | note |
|---|---|---|---|
| [human (CHM13)](#human-chm13) | 23/23 | pass | |
| [chimpanzee](#chimpanzee) | 18/18 | pass | |
| [bonobo](#bonobo) | 17/18 | pass | one miss, flagged uncertain (the true array is the reported competing cluster) |
| [gorilla](#gorilla) | 9/9 | pass | |
| [Sumatran orangutan](#sumatran-orangutan) | 8/8 | pass | |
| [Bornean orangutan](#bornean-orangutan) | 15/16 | pass | one confident miss: CENP-A on a small alpha array, 14 Mb from the large one that was called |
| [crab-eating macaque](#crab-eating-macaque) | 16/16 | pass | |
| [pig](#pig) | 13/17 | fail | truth is weak ChIP from other pigs |
| [horse (T2T)](#horse-t2t) | 29/32 | fail | chr11 (satellite-free) has no call; CM128050.1 is a `low` call 2.7 Mb from CENP-A |

## Mammals without functional truth (agreement with the satellome locus)

| genome | agree | Hi-C test |
|---|---|---|
| [siamang](#siamang) | 20/26 | pass |
| [rhesus macaque](#rhesus-macaque) | 22/22 | fail |
| [sheep](#sheep) | 27/28 | fail |
| [goat](#goat) (all acrocentric) | 30/30 | fail |
| [giraffe](#giraffe) | 15/15 | fail |

## Plants and horse EquCab3 (CENH3 / CENP-A ChIP; right = within 4 Hi-C bins)

| genome | right | Hi-C test | note |
|---|---|---|---|
| [rice](#rice) | 12/12 | pass | |
| [barley](#barley) | 6/7 | fail; arm fold passes | satellite calls at the chromosome ends; Hi-C and the arm fold moved 6 |
| [potato](#potato) | 2/11 | fail | half the centromeres are satellite-free; no TE input |
| [horse (EquCab3)](#horse-equcab3) | 3/5 | fail (on its `low` calls) | right calls are uncertain |
| [*Poa infirma*](#poa-infirma) | 7/7 | pass | truth from earlier analyses |

## Other genomes (no Hi-C)

| genome | truth | inputs | right | note |
|---|---|---|---|---|
| [*Arabidopsis thaliana* (Col-CC)](#arabidopsis-thaliana-col-cc) | CENH3 | FASTA, `tes/`, ONT methylation | 5/5 | all `high` |
| [chicken (GGswu)](#chicken-ggswu) | CENP-A windows (33 of 40 chromosomes) | FASTA | 28/33 | all `low`: no satellite family dominates |
| [maize (Mo17 T2T)](#maize-mo17-t2t) | CENH3 | FASTA | 9/10 within 1 Mb (7 overlap) | the knob repeat (23.6 Mb) is never called; chr5 confidently wrong: CENH3 on satellite-free DNA 2 Mb from CentC |
| [mouse (mhaESC T2T)](#mouse-mhaesc-t2t) | minor-satellite (CENP-A) array | FASTA | 1/20 (20/20 at the centromeric end) | all acrocentric: calls are `low` on the larger pericentric major satellite, with the minor satellite as the reported alternative |
## Plants and algae with CENH3 truth and LTR retrotransposons (`tes/`)

| genome | inputs | right | note |
|---|---|---|---|
| [moss (*Physcomitrium* V6)](#moss-physcomitrium-v6) | FASTA, `tes/`, Hi-C | 22/26 | Copia (Bryco) cluster centromeres without satellite; the 4 misses are `low` |
| [*Chlamydomonas* (UL-1690)](#chlamydomonas-ul-1690) | FASTA, `tes/` | 6/17 | LINE-cluster centromeres; every call `low` |
| [*Brachypodium distachyon* Bd21](#brachypodium-distachyon-bd21) | FASTA, `tes/` | 3/3 | CENH3 truth on 3 of 5 chromosomes; all `high` |
| [*Brassica rapa* CC1](#brassica-rapa-cc1) | FASTA, `tes/` | 10/10 | one of 11 accessions (FASTA alone: 109/110 in all 11) |
| [soybean X26](#soybean-x26) | FASTA, `tes/` | 19/20 | chr13: a CRM core replaced a `low` satellite call and is `medium`, 1.9 Mb from CENH3 (FASTA alone 20/20) |
| [sorghum BTx623](#sorghum-btx623) | FASTA, `tes/` | 10/10 | CENH3 ChIP of line V9; all `high` |
| [*Aegilops tauschii* AL8/78](#aegilops-tauschii-al878) | FASTA, `tes/` | 5/7 | calls on CRM cores; Chr5 and Chr7 sit at a core's edge, 1.3 and 3.1 Mb from CENH3 (FASTA alone 1/7) |
| [*Chlorella sorokiniana* NS4-2](#chlorella-sorokiniana-ns4-2) | FASTA, `tes/` | 8/10 | centromeres of 2-7 Osser insertions; every call `low` (FASTA alone 0/10) |
| [tomato Heinz 1706](#tomato-heinz-1706) | FASTA, `tes/`, Hi-C | 1/12 | satellite-free centromeres beside a dominant pericentromeric satellite: 7 confidently wrong, a reported failure |

## A diploid human T2T genome

| genome | truth | inputs | calls on truth | note |
|---|---|---|---|---|
| [HG002](#hg002) | CenSat active alpha-satellite arrays | FASTA, ONT CpG methylation | 46/46 | both haplotypes; 45 `high`; CpG kinetochores on 40, all inside the active array |

## Controls: no single regional centromere (any confident call is wrong)

| genome | centromeres | chromosomes | confident calls | note |
|---|---|---|---|---|
| [*Caenorhabditis elegans*](#caenorhabditis-elegans) | holocentric | 6 | 0 | 6 `low` |
| [*Bombyx mori*](#bombyx-mori) | holocentric, no CenH3 | 29 | 0 | 13 `low`, 16 none |
| [*Rhynchospora breviuscula*](#rhynchospora-breviuscula) | holocentric (Tyba repeat) | 5 | 0 | 5 `low` |
| [*Saccharomyces cerevisiae*](#saccharomyces-cerevisiae) | point | 17 | 0 | none |
| [*Schizosaccharomyces pombe*](#schizosaccharomyces-pombe) | regional | 3 | 0 | 1 `low`, 2 none |
| [*Neurospora crassa*](#neurospora-crassa) | regional, AT-rich | 7 | 0 | none |
| [*Plasmodium falciparum*](#plasmodium-falciparum) | short, AT-rich | 14 | 0 | 1 `low`, 13 none |

---

### human (CHM13)
![human CHM13](img/gallery/chm13.png)

### chimpanzee
![chimpanzee](img/gallery/mPanTro3.png)

### bonobo
![bonobo](img/gallery/mPanPan1.png)

### gorilla
![gorilla](img/gallery/mGorGor1.png)

### Sumatran orangutan
![Sumatran orangutan](img/gallery/mPonAbe1.png)

### Bornean orangutan
![Bornean orangutan](img/gallery/mPonPyg2.png)

### crab-eating macaque
![crab-eating macaque](img/gallery/MFA8.png)

### pig
![pig](img/gallery/Sscrofa.png)

### horse (T2T)
![horse T2T](img/gallery/horse1.png)

### siamang
![siamang](img/gallery/mSymSyn1.png)

### rhesus macaque
![rhesus macaque](img/gallery/MMU8.png)

### sheep
![sheep](img/gallery/sheep1.png)

### goat
![goat](img/gallery/goat1.png)

### giraffe
![giraffe](img/gallery/mGirTip1.png)

### rice
![rice](img/gallery/rice.png)

### barley
![barley](img/gallery/barley.png)

### potato
![potato](img/gallery/potato.png)

### horse (EquCab3)
![horse EquCab3](img/gallery/horse_EquCab3.png)

### Poa infirma
![Poa infirma](img/gallery/infirma.png)

### Arabidopsis thaliana (Col-CC)
![Arabidopsis thaliana Col-CC](img/gallery/colcc.png)

### chicken (GGswu)
![chicken GGswu](img/gallery/ggswu.png)

### maize (Mo17 T2T)
![maize Mo17](img/gallery/mo17.png)

### mouse (mhaESC T2T)
![mouse mhaESC](img/gallery/mmT2T.png)

### moss (Physcomitrium V6)
![moss Physcomitrium V6](img/gallery/ppat.png)

### Chlamydomonas (UL-1690)
![Chlamydomonas UL-1690](img/gallery/chlamy.png)

### Brachypodium distachyon Bd21
![Brachypodium distachyon Bd21](img/gallery/brachy.png)

### Brassica rapa CC1
![Brassica rapa CC1](img/gallery/CC1.png)

### soybean X26
![soybean X26](img/gallery/gmX26.png)

### sorghum BTx623
![sorghum BTx623](img/gallery/sbBTx623.png)

### Aegilops tauschii AL8/78
![Aegilops tauschii AL8/78](img/gallery/aetAL878.png)

### Chlorella sorokiniana NS4-2
![Chlorella sorokiniana NS4-2](img/gallery/csNS42.png)

### tomato Heinz 1706
![tomato Heinz 1706](img/gallery/tomato.png)

### HG002
![HG002](img/gallery/hg002.png)

### Caenorhabditis elegans
![Caenorhabditis elegans](img/gallery/celegans.png)

### Bombyx mori
![Bombyx mori](img/gallery/bmori.png)

### Rhynchospora breviuscula
![Rhynchospora breviuscula](img/gallery/rbrev.png)

### Saccharomyces cerevisiae
![Saccharomyces cerevisiae](img/gallery/scer.png)

### Schizosaccharomyces pombe
![Schizosaccharomyces pombe](img/gallery/spom.png)

### Neurospora crassa
![Neurospora crassa](img/gallery/ncra.png)

### Plasmodium falciparum
![Plasmodium falciparum](img/gallery/pfal.png)
