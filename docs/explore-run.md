# The exploration run (`scripts/explore/`)

Status: untrusted tooling.  Every number the exploration run produces is numerical (double precision) and unverified, and no proof may cite it.  What the run finds are *candidates*: critical configurations of N points on a sphere that the existing search code (`random_ascent` and `surgery` of `kisscrit/search.py`, unchanged) reaches, merged into one census and prepared for a later refinement stage.  The run does not prove that a census is complete; the surgery step in particular is not exhausted by any run of bounded length.

## What it does

`run_explore.sh` runs `explore_worker.py` on several cores in batches until a wall-clock budget ends:

1. **Worker** (`explore_worker.py`): random ascents, then contact-graph surgery, writing a census file of the format of `Census.dump`.  A worker is reproducible: its generator is `numpy.random.default_rng(SeedSequence([seed, index]))`, so the same seed, index and counts give the same census when the budget is not the limit.  Numerics are single-threaded (`OMP_NUM_THREADS` and the BLAS variables are set to 1 before numpy is imported).  It dumps its census at least every three minutes and on SIGTERM; each dump is written to a temporary name and renamed.
2. **Merge** (`explore_merge.py`, run after every batch): merges the worker files by the type key of `kisscrit.search.type_key` (the key of `Census.add` and `Census.load`; no new key).  It keeps the earliest copy of each type, sums `hits`, adds the max-min balanced stress `sigma` where a record has none, and reports types before and after, the yield per worker and the Euler line.  Merging always starts from the worker files, so running it again gives the same file; a file that cannot be read (for example a truncated copy) and a merge output are skipped and listed.
3. **Export** (`explore_export.py`, run at the end): writes `refine-input.json`, `refine-queue.tsv` and `summary.txt`.

From the second batch on, each worker is given a slice of the merged census: the types that no surgery call has used yet (marked `surg` once a call that used them returned), taken in the order of decreasing r, position i of every W, at most `--surgery-types` of them.  Workers do surgery on their slice first (at most 70 % of the batch) and spend the rest of the batch on ascents.  A batch cut short by its time limit marks nothing, so those types come back in a later batch.

## Running it

```
bash scripts/explore/run_explore.sh start            # N = 12, d = 3, 16 workers (at most nproc), 10 hours, nice 19
bash scripts/explore/run_explore.sh status
bash scripts/explore/run_explore.sh stop             # ends cleanly, then merges and exports
bash scripts/explore/run_explore.sh start --cores 8  # a second start with the same run home resumes the run
bash scripts/explore/run_explore.sh merge | export   # by hand, on the files found so far
```

`start` pins the code: it takes a snapshot of the repository's HEAD (`git archive`) into the run home, so a later `git pull` does not change a run in progress; uncommitted changes are not in the snapshot (the script says so).  `--tree DIR` runs the code of DIR as it is.  The run needs the Python packages of `scripts/full/run.sh`'s conda environment `kisscrit-run` (numpy, scipy, networkx); `KISSCRIT_PYTHON` names another interpreter.

| flag | default | meaning |
|---|---|---|
| `--n`, `--d` | 12, 3 | N points on S^(d-1); d = 3 is the sphere S^2 |
| `--cores` | 16, capped at `nproc` | parallel workers; an explicit value above `nproc` is refused |
| `--nice` | 19 | nice level of every worker |
| `--hours` | 10 | wall-clock budget of the whole run (decimals allowed) |
| `--runs` | 60 | random ascents per worker and batch |
| `--rounds` | 2 | surgery rounds per surgery call |
| `--near` | 8 | nearest non-contacts that surgery tries |
| `--batch-minutes` | 30 | length of a batch |
| `--surgery-types` | 60 | surgery slice per worker and batch |
| `--seed` | 20261027 | base seed; worker j is seeded from (seed, j) |
| `--run-home` | `$KISSCRIT_RUN_HOME` or `$HOME/kisscrit-run` | as in `scripts/full/run.sh`; the run lives in its `explore/` subdirectory |
| `--conda-env` | `kisscrit-run` | environment under `$KISSCRIT_MINIFORGE` (default `$HOME/miniforge3`) |

The defaults of `--runs`, `--rounds`, `--batch-minutes` and `--surgery-types` are guesses to be tuned after the first batches; `status` and `merge-report.txt` show the yield per worker.

**Priority rule.** The exploration never preempts the completeness runs: its workers run at nice 19 and the number of workers may not exceed `nproc`.  While another job runs on the server, lower the number of workers: `stop`, then `start --cores K`.  A resume refuses any other changed parameter.

**Layout** (`<run-home>/explore/`): `tree/` (the code snapshot), `batches/bNNNN/wNNNNNN.json` (worker censuses), `logs/`, `state/`, `loop.log`, and `result/` with the only files meant to leave the machine: `manifest.json` (commit, parameters, versions, start and deadline; no host, user or path), `merged-census.json`, `merge-report.txt`, `refine-input.json`, `refine-queue.tsv`, `summary.txt`.

## Reading the result

* `merged-census.json`: `Census.load` reads it.  Each record has the fields of `data/n12_above1.json` (`r`, `c`, `k`, `coindex`, `kernel`, `sdim`, `tmin`, `sym`, `sym_proper`, `pairs`, `sigma`, `X`) and the census fields (`nsupport`, `rattlers`, `nzero`, `ev`, `labeled`, `hits`, `src`), plus `prov` (kisscrit commit, base seed, worker index, batch, run id, UTC time at which the worker first saw the type) and `surg`.  `hits` counts finds; surgery also counts a rediscovery of a type as a hit, so a type reached by many surgery moves has many hits.
* `refine-input.json`: the layout of `data/n12_above1.json` (`N`, `d`, `note`, `types`), every type of the census, in the order of decreasing r.  A type with a rattler has `sym` and `sym_proper` null.  The refine scripts load it unchanged; `index` of the queue is the position in its `types`.
* `refine-queue.tsv`: one row per type, sorted by priority (new types first, then fewer hits, lower co-index, larger r).  Columns: `index`, `priority`, `flag`, `new`, `hits`, `coindex`, `kernel`, `sdim`, `rattlers`, `r`, `n_contacts`, `kkt_resid`, `min_gap`, `stress_min`, `sym`, `first_found_utc`, `run_id`, `commit`, `worker`.
* Flags, recomputed from X and the contact pairs in double precision, first rule that applies: **suspect** (KKT residual above 1e-9, or minimum non-contact gap below 1e-6); **rattler** (a vertex outside the stressed graph); **regular** (kernel 0 and a max-min balanced stress above 1e-10, no rattler: the types the refine stage solves); **degenerate** (everything else; listed, never solved here).  `kkt_resid` is the largest of | |x_i| - 1 |, | <x_i, x_j> - c | over contacts and the norm of the tangential balance of the stress; `min_gap` is c minus the largest non-contact cosine; `stress_min` is the max-min value of the stress LP.  `new` is `yes` when the type's key is in none of the known censuses given to the export (the run passes `data/n12_above1.json` for N = 12 and `results/census_n{N}.json` where it exists).
* The Euler line compares the signed count of `euler_check` with chi(M_{0,N}).  For d = 3 it uses `results/census_n{n}.json` for n < N where those files exist; it is complete only if the census is, so a difference is expected while the run is not exhaustive, and equality would not prove completeness.  On the repository's own census files the count equals chi for N = 4, 6, 7 and differs for N = 5, whose census holds a degenerate type that `euler_check` skips.
* Every output file says that its values are numerical (double precision), that nothing in it is verified, and what the thresholds are.

## Limits

* The type key is a heuristic invariant (degree sequence, Weisfeiler-Lehman hash, rounded radius); two different types can in principle share a key and be merged.
* Many types of a family with a continuous kernel (`kernel` > 0) are listed as degenerate even when their stress is unique; the refine stage may handle some of them with other means.
* The ascent and the surgery are the existing single-process Python code; speeding them up is a separate step.

## Tests

`tests/test_explore.py` (run with `python -m pytest -q tests/test_explore.py -s` to see the printed results): worker reproducibility, budget and SIGTERM; surgery slices; merge idempotence, file-order independence and a truncated file; the layout of `refine-input.json` against `data/n12_above1.json`; the flags; the refine contract (each `regular` type of short N = 6 and N = 7 runs converges to 60 digits under the Gauss-Newton of `scripts/exact/hiprec.py`); parameter validation; a 60-second two-worker run for N = 5; stop and resume.
