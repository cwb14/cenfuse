# Gallery: cenfuse's overview figure for every benchmark genome

Each image is the first page of `cenfuse_overview.pdf` from a default run on the benchmark inputs, unedited
(the *Poa annua* one is in the [README](../README.md)). Numbers: [benchmark.md](benchmark.md).

**How to read one** (panel letters as in the images):

| panel | shows |
|---|---|
| **A** | every chromosome to scale: satellite arrays colored by role, the Hi-C centromere score (blue line), the satellite call (black box) and the final call (band; blue band = placed by Hi-C, arrow = moved); telomere arrays at the ends (▶ ◀) and inside the chromosome (▼) |
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
| [horse (T2T)](#horse-t2t) | 30/32 | fail | chr11 (satellite-free) has no call |

## Mammals without functional truth (agreement with the satellome locus)

| genome | agree | Hi-C test |
|---|---|---|
| [siamang](#siamang) | 20/26 | pass |
| [rhesus macaque](#rhesus-macaque) | 22/22 | fail |
| [sheep](#sheep) | 26/28 | fail |
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
