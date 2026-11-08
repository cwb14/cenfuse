  1. Randomized Hi-C controls (NEG, NEG2, TNEG): no test passes and no call changes: **met** (57 runs; the cis test by trans contacts
     reaches support on 3, where it does not act).
  2. Holocentric genomes: no confident call: **met** (0 of 196 chromosome-runs, with and without `tes/` and Hi-C).
  3. The expected effects of each rule (PLAN.md): tomato from the FASTA and Hi-C 10 of 12 (all right), *C. neoformans* G1 12 of 14
     (10 of the 12 pairs on CENP-A), moss from the FASTA and Hi-C 5 calls wrong to right, *C. neoformans* in anaphase 3 more:
     **met**, each as predicted.

  Checked as well, not written beforehand: no call that was right became wrong, and there is no new confident error (0 of 3,676
  chromosome-runs with truth); 4 chromosome-runs without a call got a wrong `low` one (*C. neoformans* G1, 19-20 kb from CENP-A).
