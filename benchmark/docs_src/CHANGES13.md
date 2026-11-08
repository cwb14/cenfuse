Between the two releases, 49 of the 3,676 chromosome calls that have truth changed (48 moved or placed, 1 raised), all from the
cis signal tested by trans contacts, the AT trough paired with a prior-free Hi-C peak, and the support tier choosing the competing
cluster. No right call became wrong and no confident call changed; 4 chromosomes without a call got a wrong `low` one (2 *C.
neoformans* near misses, in two arms). Each release scored with the same rules ([how calls are scored](#how-calls-are-scored)):

| from (0.12) | to (0.13) | calls | where |
|---|---|---|---|
| no call | right `low` | 30 | *C. neoformans* G1 12 and 12 (FASTA + Hi-C, and with `tes/`): the cis signal tested by trans contacts; *C. neoformans* in anaphase 3 and 3: the AT trough and a prior-free peak (`at_denovo_pair`) |
| no call | wrong `low` | 4 | *C. neoformans* G1 CP003824.1 and CP003829.1 (both arms): the cis extreme 2 bins (19-20 kb) from CENP-A, wrong by the 1-bin rule |
| wrong `low` | right `low` | 9 | tomato from the FASTA and Hi-C: the cis signal moved 9 satellite calls onto CENH3 |
| wrong `low` | right `medium` | 5 | moss from the FASTA and Hi-C: the support tier's strong peaks chose the competing satellite cluster |
| right `low` | right `medium` | 1 | tomato CM175069.1: the cis extreme at the satellite call |

On randomized Hi-C (57 runs) no test passes and no call differs from the same genome's run without Hi-C; the cis test by trans contacts
reaches support on 3 of them and does not act. The holocentric genomes keep 0 confident calls (196 chromosome-runs). Without truth: *P.
annua* Pa1A moves 20.7 Mb with validated Hi-C, the second line being younger copies of its centromeric LTR families at the new site
(median 0 vs 620,000 years, P 2e-13), where earlier analyses put its functional centromere ([future.md #31](future.md#31-a-centromere-that-left-its-satellite-and-its-ltr-families));
cucumber's Hi-C, untested in 0.12 for want of confident references, is now tested on its `low` calls too and changes no call. The
satellite engine is unchanged, so are the 106 allotetraploid replays (no Hi-C).
