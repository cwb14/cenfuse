# How cenfuse works

## Principles

1. **Each line calls on its own first**: satellites (+ TEs), Hi-C, methylation.
2. **Lines check each other.** No line is trusted on its own confidence: the Hi-C engine's prior-free `high` calls
   were right on only 45% of T2T mammal chromosomes.
3. **Whether Hi-C works in a genome is measured, not assumed.** Hi-C must find the *other* chromosomes' centromeres
   (the gate) before it may act on one. CpG dips must mark them before they place a kinetochore.
4. **The final call is the sequence call** unless a validated line moves it under the rules below. Every move and
   every disagreement is reported.
5. **Fixed constants, each with a reason** (a convention, a majority rule, a data-derived value, or a named
   heuristic). Nothing is tuned to a benchmark, and there are no options for thresholds.

Not goals: centromeres from methylation alone (false dips are as deep as true ones); PacBio methylation; Hi-C
choosing the satellite family.

## Stages

```
satellite -> rdna -> genotype -> hic -> fold -> karyotype -> provisional -> meth -> fuse -> cpg -> plants -> loci -> output -> plot
```

Each stage is checkpointed, keyed on its inputs and code, so a rerun resumes. `genotype`, `karyotype` and `plants`
handle several plants ([plants.md](plants.md)).

## Lines

### Satellites and TEs (`cenfuse.sat`)

The pipeline: tandem-repeat seeds (FasTAN or TRF) → families (clustered, consensus) → genome-wide annotation
(nhmmer) → a role per family → one call per chromosome → HOR structure.

Roles: `centromeric`, `cen_associated`, `subtelomeric`, `interstitial` (knobs and other blocks), `dispersed`,
`minor`, `rDNA` (rRNA genes found with Rfam models).

| rule | test | why |
|---|---|---|
| dispersed repeat | tandem_frac < 0.5 **and** whole-unit hits < 0.5 **and** locus fold < 2 → never centromeric | interspersed elements (LINE/SINE fragments) annotate as chains of partial hits. Satellites are tandem; retrotransposon-derived centromeric arrays are concentrated (fold ≥ 2.5) |
| satellite kind | dominance judged per kind: subfamilies related by sequence are pooled | splitting alpha into 6 subfamilies let a pericentric satellite out-count it |
| acrocentric | a strong family at one end of most chromosomes, and at both ends of too few for independent ends (binomial, P ≤ 0.01), with no strong interior rival → centromeric, flag `acrocentric`, at most `medium` | subtelomeric repeats occupy ends independently; a centromere occurs once per chromosome. Needs ~16+ chromosomes |
| the call | the cluster holding classify's chosen family (layered families pooled); a decisive CRM core may favor the largest cluster | |
| competing cluster | a second cluster ≥ half the main one → two levels down (`low`), reported in `seq_alt_*` | such calls were right 2/10 times on CENP-A/CENH3 truth |
| CRM cores (`tes/`) | the CRM clade's densest cores; a decisive core places or supports the call | CRM retrotransposons mark grass centromeres |
| genes (`genes/`) | gene-rich families are vetoed; subtelomeric = beyond the outermost genes | |

### rDNA and telomeres

- **NORs** (`rdna.py`): tandem runs of ≥ 3 copies of the universal 18S motif (GTACACACCGCCCGTC), ≤ 100 kb apart.
  This finds CHM13's 5 NORs in 32 s. Hi-C peaks beside a NOR never oppose a call (NORs cluster in the nucleolus).
- **Telomeres** (`telomere.py`): the unit is found de novo; arrays are reported at chromosome ends.

### Hi-C (`cenfuse.hic` engine + `hicfuse.py`)

1. **Prior-free call**: the engine's own iteration. It decides nothing by itself; only the fold test can validate it.
2. **Anchored score.** For chromosome x: its mean contact with partner windows at the *other* chromosomes' sequence
   calls (calls ≥ `medium`; all calls if fewer than 3). x's own call never shapes its score.
   - A partner window is the call midpoint ± 4 bins when most of those bins are mappable.
   - Otherwise it is the 9 mappable bins nearest the call: the array's flanks, since the array itself is unmappable.
3. **Peak tests** (the engine's): 100 bootstraps, leave one partner out, 200 null replicates.
   - A peak is *strong* when it passes all three and is not flagged `hub`, `rDNA` (incl. beside a NOR) or
     `telomeric` (end window).
   - A peak with any of those flags may support a call but never moves or opposes one.
4. **Distances in mappable bins.**
   - Agreement: within 1 bin, forgiving ≤ 4 unmappable bins (a peak at an array's edge supports a call inside it).
   - Opposition: beyond 4 *mappable* bins. Across a longer blind stretch Hi-C neither supports nor opposes.
5. **Gate** (per chromosome x).
   - References: the other chromosomes with `high` calls (`high` + `medium` if fewer than 3).
   - x passes when most references' own peaks agree with their calls, beyond chance (Poisson-binomial P ≤ 0.01),
     and also a majority of the whole reference set.
   - Fewer than 3 references: `untested`.

### Arm fold (`fold.py`)

In a Rabl nucleus both arms run from the centromere pole to the telomere pole. The chromosome's own (cis) contact
map then has an anti-diagonal ridge through the centromere.

- **Fold score of a center**: mean observed/expected over the pairs (c - d, c + d) spanning ≥ 9 bins, both
  mappable. A center inside an unmappable array is read from its flanks.
- **Scan only centers with centromeric index ≥ 12.5%.** Below that a chromosome is telocentric (Levan 1964) and the
  short arm cannot fold. Without the bound, centers near the ends won on noise.
- **Fold point** = top center; **core** = its half-maximum run. It is coarse: 5-10% of the chromosome.
- **Fold test**: the prior-free Hi-C peak lies inside the fold core on a majority of chromosomes, beyond chance
  (P ≤ 0.01). Trans clustering and cis folding are independent readouts, so passing validates both without the
  sequence line.
- Passes in barley and P. annua; fails in all mammals, potato and rice (rice's prior-free peaks are wrong).

### Methylation (`cenfuse.meth` engine + `methfuse.py`)

- **Anchors**: provisional calls ≥ `medium`.
- **Signature**: calibrated with x's anchor left out.
- **Candidate interval**: `centromere_like` when the signature score at its most extreme CHG window is ≥ 0.5.
- **Role**: support and second line only; it never moves a call alone. Only samples from the assembled plant (or of
  unknown origin) decide.
- Needs `tes/` and a CHG file for every sample.

### CpG kinetochores (`cdr.py`)

In primates CENP-A sits in a CpG-hypomethylated dip of the methylated active array (the CDR). The satellite and Hi-C
lines locate the array (0.8-5.6 Mb); the dip locates the kinetochore inside it. Inputs: the CG file of each deciding
methylation sample (modkit bedMethyl, or a bedGraph of methylated fraction); `tes/` and CHG are not needed.

- **Profile**: mean CpG methylation per 5-kb bin (CDR-Finder's window) on bins that centromeric-role arrays cover
  for at least half, inside the final call. DNA between arrays (CpG islands, genes) is left out.
- **Dip**: the most prominent valley (higher methylation on both sides, so an array edge is never a dip), extended to
  its half-depth run (bins below the midpoint between its bottom and the array median).
- **Dip test**: the run is longer than in the call's own bins shuffled (200 shuffles; significant at P ≤ 0.05).
- **Positional test**: the run is longer than when the chromosome's CpG values are rotated by 200 random offsets onto
  the call's own CpGs and array bins (methylation from elsewhere; same size, bins and CpG density).
- **Gate** (per chromosome, leave one out): references = the other chromosomes' `high` calls. k of n have a positional
  dip; binomial tail at 0.05 (each reference's chance under the null) ≤ 0.01. No majority is required: the
  positional test is conservative (passing primates: 28-95% of tested calls; rotated tracks: 4-9%).
- **Role**: where the gate passes, every call with a significant dip test gets its kinetochore (`cpg_kt_start/end`,
  GFF3 `kinetochore`). It never moves a call or changes its confidence, and it cannot tell a wrong array.
- **Assembly check**: bedMethyl rows must sit on reference CpGs (as for CHG); bedGraph rows too, at the row start or
  one base off (public tracks differ by one base).

## Fusion rules (`fuse.py`)

| situation | result |
|---|---|
| no validated line | sequence call, `sequence_only` |
| validated Hi-C (or methylation) agrees | sequence call, `concordant`, +1 level per line (max +2) |
| validated strong Hi-C peak > 4 mappable bins away, call `low`/`medium` | moves to the Hi-C core, `offset`, `medium` |
| same, call `high` | moves only with a second line at the peak and not at the call: a CRM core within 4 bins, or methylation `centromere_like` at the peak and a measured dip at the call (`moved_with_te`/`_meth`, `high`). Otherwise `contested`, -1 |
| Hi-C peak 1-4 bins away | neither supports nor opposes (`hic_near`) |
| competing cluster, strong peak agrees with the second one, mappable bins between them | the second cluster is chosen (`hic_chose_alt`) |
| no sequence call, strong peak | Hi-C core, `functional_only`, `medium` |
| no sequence call, clean medium peak | Hi-C core, `functional_only`, `low` (`hic_medium_peak`) |
| gate not passed, fold test passed; fold core and prior-free peak agree, the sequence call at neither | moves to the prior-free core, `offset`, `medium` (`fold_pair`; `moved_with_fold` for a high call) |
| same, but the fold or the prior-free peak opposes the kept call | `contested`, -1 per opposing line (`fold_offset`, `hic_denovo_offset`) |

- **Confidence**: the satellite engine's level, then ±1 per line as above. A contest keeps the level when
  methylation supports the call. A moved or Hi-C-placed call is `medium`, `high` with a second line.
- **Support / against columns**: `satellite` (centromeric arrays inside the call, ≥ 5 kb), `te`, `hic`, `meth`,
  `fold`.
- **CpG kinetochores** are added after fusion (`cdr.py`): they change no call, confidence or flag.
- **`hic_sequence_conflict`** (per genome): Hi-C's prior-free peaks are confident but disagree with the sequence
  calls and the gate fails. One line is wrong genome-wide. The fold settles it in Rabl nuclei; otherwise add `tes/`
  or methylation.

## Loci (`loci.py`)

Every candidate interval becomes a locus member:
- satellite clusters ≥ 20 kb of feature roles;
- TE cores;
- Hi-C cores and hubs;
- methylation candidates;
- telomere arrays;
- the sequence and final calls.

Overlapping members merge; only Hi-C members are padded (1 bin). Each locus gets every line's value there and the
first kind that matches:

`centromere` > `sequence_centromere` (where the call moved from) > `telomere` > `rDNA` > `subtelomeric_block` >
`interstitial_block` (knobs) > `alt_centromere_candidate` > `hic_hub` > `hic_peak` > `te_core` > `satellite_block` >
`meth_candidate` > `other`.

## Constants (`cenfuse/params.py`; engines keep their own)

| constant | value | reason |
|---|---|---|
| HIC_TOL_BINS | 1 | Hi-C placement error |
| HIC_FAR_BINS | 4 | Hi-C opposes only beyond its partner half-window (its resolution) |
| GATE_P | 0.01 | conventional evidence level |
| GATE_MIN_REF | 3 | fewer references give no majority test |
| ANCHOR_LEVELS | high, medium | wrong anchors blur Hi-C and methylation |
| NOT_STRONG_FLAGS | hub, rDNA, telomeric | hubs and NORs contact everything; chromosome ends are the main non-centromeric trans cluster |
| MAX_RAISE | 2 | low → high at most |
| METH_SCORE | 0.5 | the methylation engine's signature rule |
| SAT_SUPPORT_BP | 5 kb | the satellite engine's weak-satellite floor |
| TANDEM_FRAC, ACRO_ALPHA | 0.5, 0.01 | majority rule; conventional level |
| LOCUS_FOLD | 2 | heuristic (as the Hi-C engine's hub fold) |
| fold: S_MIN, CI_MIN | 9 bins, 12.5% | beyond the ± 4-bin centromeric domain; Levan's telocentric bound |
| CpG (`cdr.py`): BIN, MIN_SITES, MIN_BINS | 5 kb, 10 CpGs, 10 bins | CDR-Finder's window; heuristics (a dip needs shoulders) |
| CpG: N_PERM, N_ROT, DIP_P | 200, 200, 0.05 | P resolution 10× finer than the conventional level |

## Decisions and alternatives tried

| question | chosen | tried and rejected (why) |
|---|---|---|
| fix Hi-C from inside the Hi-C data | anchor to the sequence calls | random restarts and caps: held-out accuracy unchanged |
| Hi-C partner windows over unmappable arrays | engine window if mostly mappable, else the 9 nearest mappable bins | padding every window: lost agreements in CHM13, chimpanzee, S. orangutan. Replacing only empty windows: right peaks 97 → 86 of 139 |
| blind bins in Hi-C distances | ≤ 4 forgiven for agreement, none for opposition | full credit: let Hi-C support a wrong bonobo call |
| gate chance | uniform share + majority of the whole reference set | chance from the null replicates: removed 1 of 2 spurious passes on randomized sheep Hi-C (the majority removed both) |
| NOR detection | 18S motif scan | Rfam SSU HMM search: stalled for minutes on low-complexity sequence |
| dispersed-repeat test | tandem + whole-unit + locus fold | without the locus fold: 160 batch calls moved off TE-derived centromeric arrays (Coffea, Glycine, Rosa) |
| CRM core vs classify's choice | the core may favor only the largest cluster | any cluster the core touches: lost the CENH3 call of Ae. geniculata chr02U |
| fold center selection | centromeric index ≥ 12.5% | lower confidence bound of the mean: biased to the chromosome middle (barley 2/7) |
| Hi-C peaks in end bins | support only, never oppose | (rice chr4: an end-bin peak contested a correct call) |
| numpy huge pages | off by default | on: Hi-C bootstrap 30× slower on fragmented memory |
| CpG dip inside a call | most prominent valley on array bins, half-depth run | deepest bin: lands on array edges (CHM13). Longest run below half the median: on the CDR 52/90 (valley 88/90) |
| CpG dip significance | shuffle test per call; positional (rotation) test for the gate | equal-length windows of other satellite DNA as null: infeasible where most satellite DNA is centromeric (MFA8: 3.7 Mb outside calls vs 3.5-11 Mb calls). Dips in calls vs dips in satellite clusters away from calls: passed on rotated (misplaced) methylation in gorilla, because long calls hold coherent runs of any structured profile |
| CpG choosing among satellite clusters | not done | least-methylated window anywhere: on the CDR 3/23 (CHM13). Most prominent valley across clusters: wrong cluster on 8-10/25 (orangutans). A cluster with half the call's valley floor: 6 changes in 7 primates, 1 right (B. orangutan), 2 wrong (CHM13) |

## Known risks

- The gate takes the sequence calls as references. If those are wrong genome-wide, only the fold (Rabl nuclei),
  `tes/` or methylation can rescue Hi-C.
- One wrong Hi-C peak in a genome where Hi-C mostly works passes the gate. Only `high` calls need a second line to move.
- Methylation leave-one-out signatures need ≥ 3 anchors with data.
- Hi-C resolution is a mappable bin (250 kb-2 Mb). Offsets inside an unmappable array are invisible.
- CpG kinetochores need an array long enough for a dip with shoulders (≥ 10 bins with data). A small array that is
  hypomethylated as a whole (bonobo and B. orangutan minicentromeres) is not refined, and CpG never chooses it.
- The CpG line was designed on the same primates it is benchmarked on (docs/benchmark.md says which tests are held out).
- Satellite-search memory: 40-160 GB for T2T mammals, 260 GB for chicken. Freeing each chromosome's raw hits once it
  is resolved did not lower the peak (dead_ends of the v4 benchmark): it lies inside one chromosome's hits.

## Tests

`pytest -q -n 8`: 475 tests on small synthetic fixtures.
- Unit tests of every rule above.
- A synthetic end-to-end genome with planted cases: a functional offset with a CRM core (must move), the same
  without a second line (`contested`), a hub decoy, a satellite-free chromosome, permuted Hi-C (no moves), three
  plants.
- One full run of the satellite engine on its synthetic genome.
- CpG: planted dips (found, positional), shifted tracks and dips everywhere (gate fails), both input formats,
  malformed rows, the bedGraph assembly check, and a CG-only table that changes no call.
