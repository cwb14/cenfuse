Between the two releases, 400 of the 3,062 chromosome calls that have truth changed: 189 moved, 211 changed confidence only. No call
became confidently wrong that was not already. Each release scored with the same rules ([how calls are scored](#how-calls-are-scored)):

| from (0.11) | to (0.12) | calls | where |
|---|---|---|---|
| confidently wrong | `low` and wrong | 49 | tomato 17, *Ae. tauschii* 12, siamang 11, potato 6, barley 3: the one-per-chromosome test (FASTA only), and the cis signal |
| confidently wrong | right | 12 | tomato 9, potato 2, barley 1, with `tes/` (and Hi-C) |
| no call | right | 60 | AT troughs: *P. falciparum* 13, *N. crassa* 7, *M. sympodialis* 8 (FASTA only and their randomized-Hi-C runs), *S. pombe* chr2; barley with `tes/` 2 |
| no call | wrong `low` | 2 | *P. falciparum* chr10, an AT trough 0.5 Mb off (FASTA only and its randomized run) |
| right `medium` | right `high` | 105 | the which-end test (goat 20, sheep 17, horse 16), moss with `tes/` 13, AT troughs paired with Hi-C (*P. falciparum* 10, *V. dahliae* 7, *N. crassa* 6) |
| wrong `low` | right | 31 | the LTR signature (potato 12, moss 11, tomato 6, barley 1, *Chlorella* 1), with the cis signal in tomato |
| right `high` | right `low` | 48 | siamang 43 (its satellites fail the one-per-chromosome test), potato and tomato without `tes/` |
| right `high` | wrong `low` | 1 | tomato CM175069.1 with `tes/` alone: a long-unit array on CENH3 gives way to a Tekay core 3.3 Mb off; with Hi-C the cis signal puts it back |
| right | right, at another level or moved | 80 | |
| confidently wrong | confidently wrong, at another level or moved | 10 | siamang 4, potato 2, tomato 2, soybean chr13, sheep |
| wrong `low` | wrong `low`, moved | 2 | tomato with `tes/` alone (CM175068.1, CM175073.1) |

On randomized Hi-C (57 runs) no test passes and no call differs from the same genome's run without Hi-C. Each line alone and the
morphology check are unchanged (117 of 118 karyotype classes agree with the literature). Over 106 allotetraploids (no truth), 636 of
2,355 calls moved, 244 by more than 5% of the chromosome: those land in more gene-poor DNA (median gene density 0.23 → 0.19) and homoeologous
pairs agree slightly more often (642 → 653 of 1,128).
