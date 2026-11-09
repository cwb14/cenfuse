Two genomes whose CENH3 truth was on disk and first scored in 0.14: *Aegilops geniculata* (UUMM allotetraploid, 14 chromosomes, 7.9 Gb, GCA_053567775.1; CENH3 ChIP of
Yunusbaev et al. 2026 on the same assembly, their Table S2 at Mb resolution; cenrep, cenfuse's predecessor, used it to shape two LTR
rules) and *Solanum verrucosum* (a diploid wild potato, 12
chromosomes, 0.7 Gb, GCA_031230405.1; CENH3 ChIP SRR18548893 and SRR18548892 by the frozen rule, with WGS of another accession,
SRR5349572, standing in for the input that was never deposited). Both were run in full, search included, from the FASTA; the `tes/`
arm (LTRquest annotations of the cenrep batch) reuses the FASTA arm's search. Scored within 1 Mb.

| genome | inputs | right / chromosomes | right `high` / `medium` / `low` | confidently wrong | chance expected | P |
|---|---|---|---|---|---|---|
{{HELD14_ROWS}}

- ***Ae. geniculata* from the FASTA:** its centromeric satellites fail the one-per-chromosome test (9 of 14 chromosomes hold the 14
  largest clusters, null 9, P 0.67), so every call stays `low`, and all 14 are wrong: they lie on other arrays, 6-395 Mb from CENH3.
  The rule that a failed test leaves no confident call held on a genome it was not made on.
- ***Ae. geniculata* with `tes/`:** a CRM signature (38.1-fold in its cores, significant cores on 14 of 14 chromosomes) decides and
  rejects the satellite families that lie away from its cores, the FASTA arm's 86-bp family among them (away on 11 of 12). The
  centromeric family becomes a 1,216-bp family made of chains of partial copies of a longer repeat (`partial_chain`), whose arrays
  lie in the cores (`te_concordant`): 13 of 14 right, 12 `medium` and 1 `low` (a chain call is never `high`; 0.15: chr05M went from `medium` to `low` when the age rule moved its core). The one confident error, chr05U
  (167.1-175.2 Mb, 2.1 Mb off), lies in the right part of a 9.3-Mb CRM core (164.5-173.9 Mb) whose left edge reaches CENH3
  (160-165 Mb): the pattern of [future.md #45](future.md#45-confident-calls-at-the-edge-of-a-wide-ltr-core).
- ***S. verrucosum* from the FASTA:** the test fails too (9/12, P 0.24): 5 of 12 right, all `low`. **With `tes/`** (a CRM signature,
  42-fold, decisive): 8 of 12 right, 7 of them confident (3 `high`, 4 `medium`), and one confident error, chr06 (`medium`, 1.4 Mb
  off), where the CRM core replaced a weak satellite array beside it (`weak_satellite_beside_core`). Two right satellite calls were
  lost to thin CRM cores of 2 and 4 elements (chr10, chr12, `weak_core`, `low`): the 0.12 rule that a thin decisive core takes a
  chromosome whose satellite failed its test costs right calls here
  ([future.md #36](future.md#36-no-held-out-test-of-the-013-rules)).
- **Earlier.** The satellite batch of cenfuse 0.10.1 (development data there, its docs since rewritten) placed 10 of 14 *Ae. geniculata*
  centromeres and 7 of 11 *S. verrucosum* ones within 1 Mb, with its own processing of the same ChIP; 0.15 with `tes/` places 13 of 14
  and 8 of 12.
- **Truth caveats.** *Ae. geniculata*'s domains have Mb resolution. *S. verrucosum*'s input comes from another accession, and its
  chr04 holds two CENH3 peaks of nearly equal height (28.4 and 27.3-fold in 50-kb bins); the frozen rule takes the one at 23.2-23.5 Mb.
- **Cost of a fresh run** (search included; GNU time's maximum counts the largest process only): *Ae. geniculata* 64 min at 32
  threads, 17.5 GB (the search 60 min: seeds 3, families 13, annotate 38); *S. verrucosum* 16 min at 16 threads, 3.2 GB. Adding
  `tes/` to the finished run (classify onward): 12 and 1.4 min.
