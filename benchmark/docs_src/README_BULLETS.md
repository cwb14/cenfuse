- **Confidence means something.** Over the 47 genomes with truth, each with all its inputs, `high` calls were right 406 of 410 times,
  `medium` 64/69, `low` 118/145 (0.14: 406/410, 62/67, 97/141).
- **Lines that fail their test do not act.** Of the 16 confident errors the FASTA alone made in tomato, potato, barley and *Ae. tauschii*,
  none is left: their satellites do not fall one per chromosome, so those calls stay `low`. Against 0.14, on the same 258 runs: {{DEV15_RIGHT}} more
  right calls and not one more confident error ({{DEV15_TOTAL}}).
- **Hi-C's two readouts test each other.** Without `tes/`, the contacts between chromosomes test the cis signal: from the FASTA and Hi-C,
  tomato goes from 1 to 10 of 12 right, *C. neoformans* in G1 from 0 to 12 of 14; no randomized library, mammal or holocentric genome passes.
- **Controls.** On randomized Hi-C (57 runs) and misplaced methylation (8 primates), no test passed, and no call changed against
  the same genome's run without them (56 of the 57 Hi-C runs have one). Shuffling LTR
  positions within chromosomes gives no signature (157 annotations, and 270 further shuffles of the pooled-seed rule of 0.15). The AT test passes on 3.9% of 3,336 genomes with shuffled GC blocks, under its 5% level. Five more mammals without CENP-A data agree with an independent
  satellite caller on 113/121 chromosomes.
- **LTR retrotransposons at satellite-free centromeres (0.15).** A significant LTR core that holds under half of the signature's elements now replaces a call that is `low` on its own, insertion age picks between comparable peaks, and a family concentrated across the chromosomes can seed the signature where no other can. In 13 tomatoes with `tes/`, 115 of 156 chromosomes are right (70 in 0.14), and 14 of the 72 confident calls are wrong (19%): a `medium` call from a Tekay core still deserves caution. The wild potato *S. verrucosum*: 5/12 from the FASTA, 8/12 with `tes/` (7 confident, 1 wrong); the allotetraploid *Ae. geniculata*: 13/14 with `tes/`, 1 confidently wrong, and 0/14 from the FASTA, where none is confident ([benchmark.md](docs/benchmark.md#an-allotetraploid-and-a-wild-potato), [tomatoes](docs/benchmark.md#thirteen-tomatoes)).
- **AT troughs test the karyotype, and choose between a call's own candidates (0.15).** *V. dahliae* is 8/8 right from the FASTA, and in mouse the trough picks the CENP-A minor satellite on 18 of 20 chromosomes; *C. neoformans* still does not pass the test.
- **Polyploids.** On 106 allotetraploids without truth, the agreement between homoeologous chromosomes is as high after the 0.15 LTR rules as before ({{CRIT15_PANEL}}).
- **Earlier genomes.** The genomes first scored after the 0.12 rules were written: {{HELD_README}} 0.13 reran them with the others: unchanged, but for *C. neoformans* with `tes/` and Hi-C, which gains as it does without `tes/`.
- **Kinetochores from CpG methylation** (5 great apes with DiMeLo-seq CENP-A): 39 kinetochores, a median of 45 kb
  wide inside calls of 3.9 Mb; 35 within 50 kb of CENP-A. With human, 62 of 63 lie on curated CENP-A domains.

**Where it is strong:** satellite centromeres from the FASTA alone; satellite-free centromeres from LTR retrotransposon families learned in
the genome itself (CRM in grasses, potato and *Ae. tauschii*; Bryco in moss; Osser in *Chlorella*; Tekay families in the tomato clade); AT-rich
centromeres from the FASTA alone (*Plasmodium*, *Malassezia*, *Neurospora*, *V. dahliae*); point and regional centromeres from Hi-C (yeasts, fungi);
satellite-free centromeres that Hi-C sees in cis, tested by the LTR cores or by Hi-C's own contacts between chromosomes (tomato 10/12 either
way, none confidently wrong); a centromere that moved away from its satellite, where the younger of two peaks of the centromeric LTR families and validated Hi-C agree (*P. annua* Pa1A, no ChIP yet); no confident call on holocentric chromosomes; knobs and other large satellites were not
mistaken for centromeres.

**Where it is weak:** an LTR signature can have cores beside the centromere that look like the right ones (tomato without Hi-C: 2 of the
5 such cores stay `medium`), and a call at one edge of a wide LTR core counts as agreeing with it (the *Ae. geniculata*'s one
confident error; [future.md #45](docs/future.md#45-confident-calls-at-the-edge-of-a-wide-ltr-core)); the cis Hi-C signal works in some plants and fungi but in no mammal, and without `tes/` it cannot judge calls
in the terminal eighth of a chromosome (tomato keeps 2 wrong `low` calls there); a small active array beside a large inactive one of the
same satellite is called confidently on the large one (one orangutan chromosome); right calls are often left `low` (chicken 27/33,
*Z. tritici* 19/21, *C. neoformans* 12/14); Hi-C helps only where nuclei cluster their centromeres, fold their arms or mark them in cis,
and with enough reads. In allopolyploids, a family called centromeric can be a chain of partial retrotransposon copies:
such calls are flagged and the sequence line holds them to `medium` at most (lines that pass their tests may raise them), a
[choice no genome with truth has tested](docs/future.md#1-calls-on-chains-of-partial-copies).
