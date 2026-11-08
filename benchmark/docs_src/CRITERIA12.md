  1. Over the arms with `tes/`, more right calls and no more confidently wrong than 0.11 (and the 104 confident errors over all arms):
     **met**. Arms with `tes/`: 248 → 292 right of 324, 29 → 12 confidently wrong. The 202 arms run by both releases: 104 → 43
     (the 16 arms new in 0.12, 10 of them randomized-Hi-C controls that repeat their partners' calls, add 9).
  2. Potato ≥ 8/11: **met** (10/11 with `tes/`, 9/11 with Hi-C too). Barley with `tes/` only ≥ 6/7: **met** (7/7). Moss with `tes/`
     only ≥ 24/26: **met** (26/26). Tomato ≥ 6/12 with fewer than 7 confident errors: **met** (6/12 and 2 with `tes/`; 10/12 and 0
     with Hi-C too).
  3. Held-out *B. rapa*: right on ≥ 90% of chromosomes, confident errors ≤ 2%: {{HELD_BRAPA}}.
  4. Holocentric genomes: 0 confident calls: {{HELD_HOLO}}. CHM13 with `tes/`: calls unchanged: **met** (0 of 24 differ; LTRquest
     gives its elements no lineage, so the LTR line is off). Randomized or shuffled controls, no test passes and no call changes: **met** on
     57 randomized Hi-C runs (no test passes; no call differs from the run without Hi-C in the 56 that have one) and on 157 annotations with LTR positions
     shuffled (no seed); **failed as written** on GC shuffled among chromosomes, where the AT test passes on 16 of 360 shuffles (4%): a test at
     the 5% level cannot meet "no pass"; those shuffles are not runs and change no call.
  5. `high` right ≥ 95% and `medium` ≥ 90% (47 genomes, each with all its inputs): **met**: `high` 406/410 (99.0%), `medium` 62/67
     (92.5%), `low` 82/124 (0.11: 331/345, 102/105, 100/149).
