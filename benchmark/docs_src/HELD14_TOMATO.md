Thirteen tomatoes held out from 0.14, named with their truth before any of their runs: the telomere-to-telomere assemblies of Shi et al. 2026
(BioProject PRJNA1201608; 12 chromosomes each, 0.8-1.2 Gb), wild species from every clade of the genus and three cultivated lines. Truth: the CENH3
centromeres of their Supplementary Table 8 for the 11 accessions with deposited CENH3 ChIP, and for the two with ChIP but no row there (LA2093, TS-9) the
frozen ChIP rule. Left out when the set was named: Heinz 1706 (the benchmark's own tomato), six accessions whose centromeres were published without
deposited ChIP, and two without an assembly. Each was run from the FASTA, search included, and with `tes/` (LTRquest, the benchmark tomato's command).
Scored within 1 Mb. They are one genus, and three are the species of the benchmark's tomato, so they are 13 tests of one lineage's centromeres, not of 13
lineages.

| accession | truth | FASTA: right | `high` / `medium` / `low` | confidently wrong | P | with `tes/`: right | `high` / `medium` / `low` | confidently wrong | P |
|---|---|---|---|---|---|---|---|---|---|
{{TOM14_ROWS}}

Pooled, {{TOM14_POOLED}}.

- **From the FASTA, no confident call.** The CENH3 domains hold a median of 4.5 kb of annotated satellite (0.25% of the domain; 100 kb or more on 24 of 156), and cenfuse says
  it is unsure: 20 of 156 right against 11.3 expected by chance, every call `low`. In *S. chilense* the one-per-chromosome test passes (10 of 12 chromosomes, P 0.045) on arrays that are mostly not centromeric, and
  the calls stay `low` only because no satellite dominates: the two safeguards are not redundant.
- **With `tes/`, more than three times as many right, and confident errors.** A Tekay signature decides in 12 of the 13 accessions (18- to 53-fold; *S. chmielewskii*'s
  best is SIRE at 9.5-fold, under the 10-fold level, and flags only). 54 calls are right and confident, but 15 of the 69 confident calls are wrong (22%; 12 of 60 without
  LA2951, below), as in the benchmark's own tomato (2 of 8). They lie a median 4.0 Mb from CENH3, two of them just past the 1-Mb tolerance. This is the
  held-out estimate of how far a `medium` call from a Tekay core can be trusted in this genus.
- **What each branch did** (future.md [#36](future.md#36-no-held-out-test-of-the-013-rules)): the branches that move a call to a strong core are right on 33 of 39, and the
  thin cores of `weak_core` on 3 of 8, all `low`. Calls that merely agree with the core are the worst confident ones: 6 wrong of 17.

{{TOM14_BRANCHES}}

- **What the line leaves unused.** {{TOM14_DISAGREE}}. That rule (a lineage other than CRM replaces a call only where its core holds most of the signature's elements on the
  chromosome) has a reason, that the signature's other elements say the locus is not unique; here it costs right calls. Seen on these genomes, it can only be a lead
  ([future.md #46](future.md#46-crm-decides-by-its-name)): no rule was changed on them.

{{TOM14_TE}}

- **Truth caveats.** {{TOM14_CROSS}}. *S. lycopersicoides*' ChIP is weak by the frozen rule (at most 5.2-fold, and the rule's peak inside Table 8's interval on 5 of 12 chromosomes,
  against 12 of 12 in *S. habrochaites*); Table 8 stayed its truth, as written in advance. Without it the `tes/` arms are 62 of 144 right and 12 confidently wrong. The LTRquest annotations
  come from the environment as of 2026-11-05, later than the benchmark's (2026-10-27).
- **Cost of a fresh run** (search included, 16 threads, GNU time's maximum, the largest process only): 13-50 min and 7-46 GB (*S. lycopersicoides*, 1.2 Gb); LTRquest, 12 threads:
  1.8-4.1 h. Adding `tes/` to the finished run: 1-3 min.
