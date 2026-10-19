# Vendored engines

Copied by `tools/vendor.py` from clean git checkouts. Logic is unchanged; only imports (absolute to relative) and test-helper module names were rewritten. Local changes made after vendoring are listed below the table.

| package | source | commit |
|---|---|---|
| `cenfuse.sat` | cenrep | eaa58a1 Add FastTAN |
| `cenfuse.hic` | cenhic | e47cbf0 cenhic 0.1.1 |
| `cenfuse.meth` | cenmeth | dba6c53 cenmeth 0.1.0 |

## Local changes after vendoring

- `cenfuse/sat/seeds.py`: `genome_key()`. The seed-cache key adds the chromosome set when a genome uses a subset of
  its FASTA's sequences (cenfuse's chromosome-only view), so a subset and the whole file never share a cache entry.
  Results unchanged; chromosome-only FASTAs keep cenrep's plain md5 key.
- `cenfuse/hic/cli.py`: `run_call()` also returns its in-memory state (matrix, balancing, bins, location, peaks,
  hubs), so cenfuse scores anchored partner windows on the same matrix. Outputs unchanged.
