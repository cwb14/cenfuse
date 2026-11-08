- **Confidence means something.** Over the 47 genomes with truth, each with all its inputs, `high` calls were right 406 of 410 times,
  `medium` 62/67, `low` 97/141 (0.12: 406/410, 62/67, 82/124).
- **Lines that fail their test do not act.** Of the 16 confident errors the FASTA alone made in tomato, potato, barley and *Ae. tauschii*,
  none is left: their satellites do not fall one per chromosome, so those calls stay `low`. Against 0.12, on the same 258 runs: 44 more
  right calls and not one more confident error.
- **Hi-C's two readouts test each other.** Without `tes/`, the contacts between chromosomes test the cis signal: from the FASTA and Hi-C,
  tomato goes from 1 to 10 of 12 right, *C. neoformans* in G1 from 0 to 12 of 14; no randomized library, mammal or holocentric genome passes.
- **Controls.** On randomized Hi-C (57 runs) and misplaced methylation (8 primates), no test passed, and no call changed against
  the same genome's run without them (56 of the 57 Hi-C runs have one). Shuffling LTR
  positions within chromosomes (157 annotations) gives no signature. Five more mammals without CENP-A data agree with an independent
  satellite caller on 113/121 chromosomes.
- **Held out from 0.14: the first allotetraploid with CENH3 truth.** *Ae. geniculata* (UUMM, 14 chromosomes) from the FASTA: 0 of 14
  right, but none confident, because its satellites fail the one-per-chromosome test; with `tes/`, 13 of 14 right (all `medium`) and 1
  confidently wrong. The wild potato *S. verrucosum*: 5/12 (all `low`) from the FASTA, 8/12 with `tes/` (7 of them `high` or `medium`,
  1 confidently wrong) ([benchmark.md](docs/benchmark.md#held-out-from-014-an-allotetraploid-and-a-wild-potato)).
- **Held out from 0.14: 13 tomatoes of one genus.** From the FASTA 20 of 156 chromosomes are right and none is confident: their centromeres hold little
  satellite, and cenfuse says it is unsure. With `tes/` 70 of 156 are right, but 15 of the 69 confident calls are wrong (22%), so in this genus a `medium` call
  from a Tekay core deserves caution ([benchmark.md](docs/benchmark.md#held-out-from-014-thirteen-tomatoes)).
- **Held out from the 0.12 rules** (never scored before they were frozen): {{HELD_README}} 0.13 reran them with the others (it was
  developed with them in view): unchanged, but for *C. neoformans* with `tes/` and Hi-C, which gains as it does without `tes/`.
- **Kinetochores from CpG methylation** (5 great apes with DiMeLo-seq CENP-A): 39 kinetochores, a median of 45 kb
  wide inside calls of 3.9 Mb; 35 within 50 kb of CENP-A. With human, 62 of 63 lie on curated CENP-A domains.

**Where it is strong:** satellite centromeres from the FASTA alone; satellite-free centromeres from LTR retrotransposon families learned in
the genome itself (CRM in grasses, potato and *Ae. tauschii*; Bryco in moss; Osser in *Chlorella*; a Tekay family in tomato); AT-rich
centromeres from the FASTA alone (*Plasmodium*, *Malassezia*, *Neurospora*); point and regional centromeres from Hi-C (yeasts, fungi);
satellite-free centromeres that Hi-C sees in cis, tested by the LTR cores or by Hi-C's own contacts between chromosomes (tomato 10/12 either
way, none confidently wrong); a centromere that moved away from its satellite, where validated Hi-C and young copies of the centromeric
LTR families agree (*P. annua* Pa1A, no ChIP yet); no confident call on holocentric chromosomes; knobs and other large satellites were not
mistaken for centromeres.

**Where it is weak:** an LTR signature can have cores beside the centromere that look like the right ones (tomato without Hi-C: 2 of the
5 such cores stay `medium`), and a call at one edge of a wide LTR core counts as agreeing with it (the held-out *Ae. geniculata*'s one
confident error; [future.md #45](docs/future.md#45-confident-calls-at-the-edge-of-a-wide-ltr-core)); the cis Hi-C signal works in some plants and fungi but in no mammal, and without `tes/` it cannot judge calls
in the terminal eighth of a chromosome (tomato keeps 2 wrong `low` calls there); a small active array beside a large inactive one of the
same satellite is called confidently on the large one (one orangutan chromosome); right calls are often left `low` (chicken 27/33,
*Z. tritici* 19/21, *C. neoformans* 12/14); Hi-C helps only where nuclei cluster their centromeres, fold their arms or mark them in cis,
and with enough reads. In allopolyploids, a family called centromeric can be a chain of partial retrotransposon copies:
such calls are flagged and the sequence line holds them to `medium` at most (lines that pass their tests may raise them), a
[choice no genome with truth has tested](docs/future.md#1-calls-on-chains-of-partial-copies).
