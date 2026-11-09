0.14 changed no rule. On all 260 release runs, every table cenfuse writes (calls, summary, loci, feature tables, line tests,
LTR signature, provenance, plants; and the calls' BED and GFF3) equals 0.13.0's cell by cell, apart from the new columns
(v13 `bench/cmp_runs.py`: {{CMP14}}). What changed is how a run can be trusted and repeated:

- **Checkpoints** stand on the code a stage runs (its modules and their imports, without comments or docstrings), its parameters
  and the content of its inputs. Adding `tes/` or `genes/` to a finished run now recomputes the satellite engine from classify
  by itself; a moved or copied folder resumes; an unreadable checkpoint is recomputed with a warning; `--force all` exists.
- **One run per output folder** (a lock with a heartbeat; stale locks are taken over), and the shared seed cache no longer waits
  forever on a lock left by a job killed on another host.
- **A run manifest** (`run_info.json`): versions of cenfuse, Python, the packages and every external tool, the inputs'
  fingerprints, how the run ended, and what it cost (CPU hours and the peak memory of every process the run started). The
  engines report cenfuse's version; the seed cache knows FasTAN's binary.
- **Atomic writes:** every output and checkpoint is written beside its target and moved over it, so a folder copied with hard
  links and recomputed leaves the original untouched (the satellite engine rewrote its tables in place before).
- **Hi-C faster, same numbers:** the engine's bootstrap scores cenfuse's anchored windows on its own resamples instead of a
  second bootstrap, and builds each resample straight into the matrix's layout ({{HIC_SPEED14}}).
- **Worker processes** send their log lines to the main process under every start method (spawn on macOS, forkserver in
  Python 3.14); the suite passes in an environment rebuilt from `environment.lock.txt`, and a CI workflow runs it.
- **New columns:** `rule` (which integration rule set each call) in the calls table and GFF3; `incomplete` in the summary (a
  report section or figure that failed, instead of a silent warning). Every column of every table is documented
  ([outputs.md](outputs.md), [plants.md](plants.md)).
- **The benchmark in the repository** ([benchmark/](../benchmark/README.md)): the truth manifest, the frozen scorer as a tested
  module, and the chain that writes this page's tables and fills it from templates, reproducing 0.13's tables byte for byte.
- **New truth, held out:** *Ae. geniculata* and *S. verrucosum* ([section 3](#an-allotetraploid-and-a-wild-potato)), and 13 tomatoes
  ([section 3](#thirteen-tomatoes)).
- **The benchmark's scoring, measured:** how close a right call is ([How calls are scored](#how-calls-are-scored)), and the truth rule's two numbers changed: no outcome moves.
