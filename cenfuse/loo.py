"""Leave one chromosome out: the status at which a tested line may act on chromosome c (0.16.1).

A line tested against another (Hi-C against the LTR cores, the cis extremes against the LTR cores or the trans contacts, the arm
fold and the trans arm alignment against the prior-free Hi-C peaks, the AT troughs and the repeat peaks against Hi-C) passes on
the agreements of all chromosomes, c's among them. Where the line then raises, moves or places c's call, c's own agreement would
count twice: once in the test, and again as the evidence the rule acts on. So the line acts on c at the stricter of two statuses:
the genome's, and the one the same test (same statistic, level and minimums) reaches on the other chromosomes alone. A
chromosome that disagrees keeps the genome's status, so its disagreement still counts against the line: the rule of the Hi-C
gate, whose majority is taken over the whole reference set (hicfuse.gate). Genome-level verdicts (summary, report) stay
genome-wide."""
from __future__ import annotations

RANK = {"pass": 2, "support": 1}          # statuses at which a line acts; any other (fail, untested, ambiguous) lets it do nothing


def stricter(genome: str, others: str) -> str:
    """The status of lower rank; the genome's on a tie."""
    return others if RANK.get(others, 0) < RANK.get(genome, 0) else genome


def by_chromosome(genome: str, names, test) -> dict:
    """{c: stricter(genome, test(the names other than c))} for every c of names; test(names) -> a status."""
    names = list(names)
    return {c: stricter(genome, test([x for x in names if x != c])) for c in names}


def status_on(res: dict | None, chrom: str):
    """The status a test's result gives chromosome chrom: its own entry (status_by_chrom), or the genome's status for a chromosome
    the test did not read (it holds no agreement to leave out); None without a result."""
    if not res:
        return None
    return (res.get("status_by_chrom") or {}).get(chrom, res.get("status"))
