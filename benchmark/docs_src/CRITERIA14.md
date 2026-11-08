  1. No call may change on the benchmark (0.14 changed no rule): {{CRIT14_SAME}}.
  2. The held-out genomes (*Ae. geniculata*, *S. verrucosum*, each from the FASTA and with `tes/`): calls equal to 0.13's: **met**
     ({{HELD14_SAME}}). Their scores, first reported here, are in [section 3](#held-out-from-014-an-allotetraploid-and-a-wild-potato).
  3. Randomized Hi-C controls: no test passes and no call changes; holocentric genomes: no confident call: {{CRIT14_CONTROLS}}.
  4. A code or input change recomputes exactly the stages it should, and a resume uses its checkpoints: **met** (unit and end-to-end
     tests, tests/fuse/test_checkpoints.py and test_fuse_e2e.py).
  5. The scorer in the repository equals the frozen one on every 0.13 release arm: **met** (byte for byte, and the whole table chain).
  6. The one rule considered (#45: agreement at the LTR core's maximum, written down before its measurement), adopted only if none
     of the 125 right `te_concordant` calls drops below `medium`: **not met** (3 would), so not built.
  7. Thirteen held-out tomatoes (named before any run; no rule changed on them), expected as the benchmark's tomato: from the FASTA most calls `low` and wrong: **met** (20 of
     156 right, none confident); with `tes/`, about half right with confident errors: **as expected, and reported as a limit** (70 of 156 right, 15 of 69 confident calls wrong;
     the benchmark's tomato: 6 of 12, 2 of 8) ([section 3](#held-out-from-014-thirteen-tomatoes)).
  8. The changes made after the release pass (atomic writes, the merged Hi-C bootstrap and resampler, the fast leave-one-out) change no number: **met** (every table of the 260
     runs and the four held-out arms, and again of the 138 runs with Hi-C).
  9. A folder copied with hard links and recomputed leaves the original as it was: **met** (four held-out arms recomputed from classify in such copies: every file of the originals
     unchanged, md5).
