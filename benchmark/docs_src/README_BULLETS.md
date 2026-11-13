- **Confidence is calibrated.** Over the 47 genomes with truth, each with all its inputs, `high` calls were right 406
  of 410 times, `medium` 64 of 69 and `low` 118 of 145.
- **Evidence that fails its test changes nothing.** On randomized Hi-C (57 runs) and misplaced methylation (8
  primates) no test passed and no call changed. Shuffled LTR positions give no centromeric signature; shuffled GC
  blocks pass the AT test 3.9% of the time, under its 5% level. On holocentric genomes no call is confident.
- **Satellite-free centromeres are harder.** In 13 tomato genomes with `tes/`, 115 of 156 chromosomes are right, but
  14 of the 72 confident calls are wrong: treat a `medium` call from a Tekay core with caution. The allotetraploid
  *Ae. geniculata*: 13 of 14 right with `tes/` (1 confidently wrong); the wild potato *S. verrucosum*: 8 of 12.
- **Kinetochores from CpG methylation** (5 great apes with DiMeLo-seq CENP-A): a median 45 kb wide inside calls of
  3.0 Mb; 53 of 63 lie within 50 kb of CENP-A. All five pass with ONT and with PacBio.

**Where it is weak**
- An LTR core beside the centromere can look like the right one, and a call at one edge of a wide core counts as
  agreeing with it.
- Hi-C helps only where nuclei cluster their centromeres, fold their arms or mark them in cis, and only with enough
  reads (about 70 M pairs for a 3-Gb genome). Its cis signal works in some plants and fungi but in no mammal.
- A small active array beside a large inactive one of the same satellite is called on the large one (one orangutan
  chromosome).
- Right calls are often left `low` (chicken 27 of 33, *Z. tritici* 19 of 21).
- In allopolyploids, a family that looks like a satellite can be a chain of partial retrotransposon copies. Such calls
  are flagged (`partial_chain`) and held to `medium` at most, a choice no genome with truth has tested yet.
