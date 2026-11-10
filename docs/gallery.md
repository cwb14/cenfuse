# Gallery: cenfuse's overview page for every benchmark genome

Each image is the first page of `cenfuse_overview.pdf`, unedited, from the release rerun of cenfuse 0.15.0
with the genome's richest benchmarked inputs. Counts: right / chromosomes with truth (for mammals without
CENP-A data: agreement with an independent satellite caller). Every number: [benchmark.md](benchmark.md).

**How to read one** (panel letters as in the images):

| panel | shows |
|---|---|
| **A** | every chromosome to scale: satellite arrays colored by role, the Hi-C centromere score (blue line), with `tes/` the density of the marker retrotransposon lineage (orange), the sequence call (black box) and the final call (band; blue band = placed by Hi-C, arrow = moved); telomere arrays at the ends (▶ ◀) and inside (▼); 45S (◆) and 5S (◇) rDNA arrays |
| **B** | which lines support (●) or oppose (×) each final call, then its relation and confidence (`fold` and `align` columns appear only where those Hi-C tests acted) |
| **C** | where each line puts the centromere, relative to the sequence call |
| **D** | the Hi-C score around every final call (prior-free score where no sequence call anchors it): a peak at 0 means the calls sit where centromeres meet in the nucleus |

## Mammals with CENP-A truth

| genome | inputs | right | note |
|---|---|---|---|
| [human (CHM13)](#human-chm13) | FASTA, genes, Hi-C | 23/23 |  |
| [chimpanzee](#chimpanzee) | FASTA, genes, Hi-C | 18/18 |  |
| [bonobo](#bonobo) | FASTA, genes, Hi-C | 17/18 | the miss is `low`: the true array is the reported competing cluster |
| [gorilla](#gorilla) | FASTA, genes, Hi-C | 9/9 |  |
| [Sumatran orangutan](#sumatran-orangutan) | FASTA, genes, Hi-C | 8/8 |  |
| [Bornean orangutan](#bornean-orangutan) | FASTA, genes, Hi-C | 15/16 (1 confidently wrong) | one confident miss: CENP-A on 84 kb of small arrays, 14 Mb from the large array called |
| [crab-eating macaque](#crab-eating-macaque) | FASTA, genes, Hi-C | 16/16 |  |
| [pig](#pig) | FASTA, genes, Hi-C | 13/17 (3 confidently wrong) | truth is weak ChIP from other pigs |
| [horse (T2T)](#horse-t2t) | FASTA, genes, Hi-C | 29/32 | chr11 is satellite-free and has no call |

## Mammals without CENP-A data (agreement with an independent satellite caller)

| genome | inputs | agree | note |
|---|---|---|---|
| [siamang](#siamang) | FASTA, genes, Hi-C | 20/26 agree (2 confident disagreements) |  |
| [rhesus macaque](#rhesus-macaque) | FASTA, genes, Hi-C | 22/22 agree |  |
| [sheep](#sheep) | FASTA, genes, Hi-C | 26/28 agree (1 confident disagreement) |  |
| [goat](#goat) | FASTA, genes, Hi-C | 30/30 agree | all acrocentric |
| [giraffe](#giraffe) | FASTA, genes, Hi-C | 15/15 agree |  |

## Other mammals and a bird

| genome | inputs | right | note |
|---|---|---|---|
| [HG002 (diploid human)](#hg002-diploid-human) | FASTA, ONT CpG methylation | 46/46 | both haplotypes; CpG kinetochores on 40, all inside the active array |
| [horse (EquCab3)](#horse-equcab3) | FASTA, Hi-C | 3/5 | older assembly; right calls are `low`; chr11 has no call |
| [mouse (T2T)](#mouse-t2t) | FASTA, Hi-C | 18/20 | all acrocentric: on 17 chromosomes the AT trough moves the call from the major satellite to the CENP-A minor satellite beside it; every call `low` |
| [chicken](#chicken) | FASTA, Hi-C | 27/33 | CENP-A windows on 33 of 40 chromosomes; every call `low` (no satellite family dominates) |

## Plants with satellite centromeres

| genome | inputs | right | note |
|---|---|---|---|
| [rice](#rice) | FASTA, `tes/`, genes, Hi-C | 12/12 |  |
| [*Brassica rapa* CC1](#brassica-rapa-cc1) | FASTA, `tes/`, Hi-C | 10/10 | one of 11 accessions (FASTA alone: 110/110 in all 11) |
| [soybean](#soybean) | FASTA, `tes/`, Hi-C | 19/20 (1 confidently wrong) | chr13: a CRM core replaced a `low` satellite call, 1.9 Mb from CENH3 (FASTA alone 20/20) |
| [sorghum](#sorghum) | FASTA, `tes/`, Hi-C | 10/10 | CENH3 ChIP of line V9 |
| [*Brachypodium distachyon*](#brachypodium-distachyon) | FASTA, `tes/`, Hi-C | 3/3 | CENH3 truth on 3 of 5 chromosomes |
| [*Arabidopsis thaliana* (Col-CC)](#arabidopsis-thaliana-col-cc) | FASTA, `tes/`, ONT methylation | 5/5 |  |
| [maize (Mo17)](#maize-mo17) | FASTA, `tes/`, genes | 9/10 (1 confidently wrong) | the knob repeat (23.6 Mb) is never called; chr5: CENH3 on satellite-free DNA 2 Mb from CentC |
| [*Poa infirma*](#poa-infirma) | all four inputs | 7/7 | truth from earlier analyses; its sister *P. annua* is the README example |

## Plants and algae with retrotransposon or satellite-free centromeres

| genome | inputs | right | note |
|---|---|---|---|
| [barley](#barley) | FASTA, `tes/`, Hi-C | 7/7 | satellites fail the one-per-chromosome test; the CRM signature places all 7, and Hi-C's test then passes on 7 |
| [potato](#potato) | FASTA, `tes/`, genes, Hi-C | 9/11 (1 confidently wrong) | half the centromeres are satellite-free; the CRM signature places them; Hi-C's cis signal, tested against its cores, moves 3 `low` calls |
| [*Aegilops tauschii*](#aegilops-tauschii) | FASTA, `tes/`, Hi-C | 5/7 (2 confidently wrong) | calls on CRM cores; Chr5 and Chr7 at a core's edge (FASTA alone 1/7) |
| [moss (*Physcomitrium patens*)](#moss-physcomitrium-patens) | FASTA, `tes/`, Hi-C | 26/26 | Copia (Bryco) cluster centromeres without satellite; the Bryco signature decides |
| [*Chlorella sorokiniana*](#chlorella-sorokiniana) | FASTA, `tes/` | 9/10 | centromeres of 2-7 Osser (Copia) insertions; the Osser signature decides |
| [*Chlamydomonas reinhardtii*](#chlamydomonas-reinhardtii) | FASTA, `tes/` | 6/17 | LINE-cluster centromeres; every call `low` |
| [tomato](#tomato) | FASTA, `tes/`, Hi-C | 10/12 | satellite-free centromeres beside a dominant satellite: a Tekay family and Hi-C's cis signal place 10 of 12, none confidently wrong |

## Point, regional and AT-rich centromeres (Hi-C)

| genome | inputs | right | note |
|---|---|---|---|
| [*Saccharomyces cerevisiae*](#saccharomyces-cerevisiae) | FASTA, Hi-C | 16/16 | Hi-C alone at 10-kb bins; the trans arm alignment makes 13 calls `medium` |
| [*Schizosaccharomyces pombe*](#schizosaccharomyces-pombe) | FASTA, Hi-C | 3/3 | chr2's Hi-C peak is too weak for a call; its AT trough, which Hi-C's peaks validate, places it |
| [*Neurospora crassa*](#neurospora-crassa) | FASTA, Hi-C | 7/7 | the arm fold validates Hi-C; AT troughs alone place all 7 |
| [*Plasmodium falciparum*](#plasmodium-falciparum) | FASTA, Hi-C | 14/14 | AT troughs from the FASTA alone (13/14); Hi-C at 10-kb bins agrees with them |
| [*Verticillium dahliae*](#verticillium-dahliae) | FASTA, Hi-C | 8/8 |  |
| [*Zymoseptoria tritici*](#zymoseptoria-tritici) | FASTA, Hi-C | 19/21 | 8 accessory chromosomes |
| [*Cryptococcus neoformans* (anaphase Hi-C)](#cryptococcus-neoformans-anaphase-hi-c) | FASTA, Hi-C | 13/14 | Hi-C in anaphase (centromeres clustered); an AT trough with a weak prior-free peak places 3 more |
| [*Cryptococcus neoformans* (G1 Hi-C)](#cryptococcus-neoformans-g1-hi-c) | FASTA, Hi-C | 12/14 | Hi-C in G1 (centromeres unclustered): the cis signal, tested by trans contacts, places all 14, `low` |
| [*Malassezia sympodialis*](#malassezia-sympodialis) | FASTA, Hi-C | 8/8 | truth is a kinetochore protein (Mtw1) |
| [*Candida tropicalis*](#candida-tropicalis) | FASTA, Hi-C | 6/7 | 3C-seq: no test passes; Hi-C's own peaks place 6 `low` calls |

## Holocentric controls (any confident call is wrong)

| genome | inputs | right | note |
|---|---|---|---|
| [*Caenorhabditis elegans*](#caenorhabditis-elegans) | FASTA, Hi-C | 0 confident / 6 | every call `low` |
| [*Caenorhabditis nigoni*](#caenorhabditis-nigoni) | FASTA, Hi-C | 0 confident / 6 | every call `low` |
| [*Bombyx mori*](#bombyx-mori) | FASTA, Hi-C | 0 confident / 29 | no CenH3; 13 `low`, 16 without a call |
| [*Rhynchospora breviuscula*](#rhynchospora-breviuscula) | FASTA, Hi-C | 0 confident / 5 | Tyba arrays along every chromosome; every call `low` |

---

### human (CHM13)
![human (CHM13)](img/gallery/chm13.png)

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
![horse (T2T)](img/gallery/horse1.png)

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

### HG002 (diploid human)
![HG002 (diploid human)](img/gallery/hg002.png)

### horse (EquCab3)
![horse (EquCab3)](img/gallery/horse_EquCab3.png)

### mouse (T2T)
![mouse (T2T)](img/gallery/mmT2T.png)

### chicken
![chicken](img/gallery/ggswu.png)

### rice
![rice](img/gallery/rice.png)

### *Brassica rapa* CC1
![Brassica rapa CC1](img/gallery/CC1.png)

### soybean
![soybean](img/gallery/gmX26.png)

### sorghum
![sorghum](img/gallery/sbBTx623.png)

### *Brachypodium distachyon*
![Brachypodium distachyon](img/gallery/brachy.png)

### *Arabidopsis thaliana* (Col-CC)
![Arabidopsis thaliana (Col-CC)](img/gallery/colcc.png)

### maize (Mo17)
![maize (Mo17)](img/gallery/mo17.png)

### *Poa infirma*
![Poa infirma](img/gallery/infirma.png)

### barley
![barley](img/gallery/barley.png)

### potato
![potato](img/gallery/potato.png)

### *Aegilops tauschii*
![Aegilops tauschii](img/gallery/aetAL878.png)

### moss (*Physcomitrium patens*)
![moss (Physcomitrium patens)](img/gallery/ppat.png)

### *Chlorella sorokiniana*
![Chlorella sorokiniana](img/gallery/csNS42.png)

### *Chlamydomonas reinhardtii*
![Chlamydomonas reinhardtii](img/gallery/chlamy.png)

### tomato
![tomato](img/gallery/tomato.png)

### *Saccharomyces cerevisiae*
![Saccharomyces cerevisiae](img/gallery/scer.png)

### *Schizosaccharomyces pombe*
![Schizosaccharomyces pombe](img/gallery/spom.png)

### *Neurospora crassa*
![Neurospora crassa](img/gallery/ncra.png)

### *Plasmodium falciparum*
![Plasmodium falciparum](img/gallery/pfal.png)

### *Verticillium dahliae*
![Verticillium dahliae](img/gallery/vdah.png)

### *Zymoseptoria tritici*
![Zymoseptoria tritici](img/gallery/ztri.png)

### *Cryptococcus neoformans* (anaphase Hi-C)
![Cryptococcus neoformans (anaphase Hi-C)](img/gallery/cneoA.png)

### *Cryptococcus neoformans* (G1 Hi-C)
![Cryptococcus neoformans (G1 Hi-C)](img/gallery/cneoG1.png)

### *Malassezia sympodialis*
![Malassezia sympodialis](img/gallery/msym.png)

### *Candida tropicalis*
![Candida tropicalis](img/gallery/ctro.png)

### *Caenorhabditis elegans*
![Caenorhabditis elegans](img/gallery/celegans.png)

### *Caenorhabditis nigoni*
![Caenorhabditis nigoni](img/gallery/cnig.png)

### *Bombyx mori*
![Bombyx mori](img/gallery/bmori.png)

### *Rhynchospora breviuscula*
![Rhynchospora breviuscula](img/gallery/rbrev.png)
