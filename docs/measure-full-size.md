# Sampled size of the `full` run

Tool: `scripts/measure/sample_full_size.py` (version 1.0). Status: measurement tooling, no part of the proof. Its output is a projection from a random sample of triangulations. It is never the count of record: the full enumeration of all 7,595 triangulations replaces it.

## What it does

The `full` route enumerates candidates from all 7,595 triangulations with 12 vertices (`kisscrit/completeness/candidates.py`). The four `cond` triangulations (indices 0, 1, 5, 7) are already enumerated in `results/completeness/gen`: 45,260 distinct codes. The script draws a seeded simple random sample of `n` of the other 7,591 triangulations, runs the unchanged enumerator on each with a wall-clock cap, and scales up. Per triangulation it records the tilings produced, the distinct canonical codes before and after the linear relaxation, wall and CPU seconds, and whether the cap was reached. A capped triangulation keeps what was enumerated before the cap and is never extrapolated.

## Running it

    python scripts/measure/sample_full_size.py --out sample.json --seed 1 --n 150 --cap-seconds 1200 --include-known

Defaults: `--n 150`, `--cap-seconds 1200`, `--workers` = number of cores. `--include-known` first re-enumerates indices 1, 5 and 7 (cap `--calib-cap-seconds`, default 3600) and compares them with the gen files code by code. A completed index must reproduce the file; a capped one must find only codes the file lists. Anything else stops the run with exit status 2, before any sampling. Exit status 3 means an enumeration raised an exception. The output contains no host name, user name or path.

It writes `sample.json` (header with commit, seed, sampled indices and library versions; per-triangulation records; estimates; overlap; projection) and `sample.summary.txt`, which ends in one sentence of the form "The full run generates about N candidates, about H CPU-hours, about D days on the available servers (projection from a sample of n triangulations)."

## Reading the numbers

- **A** = 45,260 + 7,591 x (sample mean of the codes per triangulation). It adds the codes of each triangulation separately, so a code that arises from several triangulations is counted several times: A estimates a sum that is an upper bound for the number of distinct candidates (the sampling error comes on top), provided no sampled triangulation was capped. In the `cond` gen files the 45,260 distinct codes occur in 1, 2 or 3 of the four triangulations (9,140, 26,664 and 9,456 codes), so the overlap is not small.
- **B** = the same with the codes already in the `cond` set removed from the sample. It is a tighter sum, but still an upper bound on the distinct candidates outside `cond`, because codes shared between two non-`cond` triangulations are still counted twice.
- The overlap measured inside the sample (codes found in more than one sampled triangulation) is a lower bound for the overlap of the whole population: a code's other triangulations are each in the sample only with probability n/7,591. It is reported, not extrapolated.
- If a sampled triangulation is capped, A is also given with the capped ones dropped and the rest reweighted. The as-recorded value estimates a sum that is at most the full sum; the reweighted value is biased low when the capped triangulations are the large ones. With any capped triangulation the total is bracketed from neither side.
- The standard error uses the finite-population correction; the 95% interval is a percentile bootstrap of the sample. With heavy-tailed counts and a few hundred triangulations both tend to be too short on the upper side; the summary reports the median, mean, maximum and the share carried by the five largest.
- **H** = search CPU-hours at the two pilot rates (3.29 CPU-hours per 45,260 candidates on a 64-core server, `results/pilot/cond-2026-10-22/summary-cond.json`; 9.28 for the older `cond` run, `results/completeness/summary-cond.json`) plus checker CPU-hours at 0.24 to 0.83 of the search time (`docs/decisions.md`, 2026-10-20). **D** = H / 64 / 24, assuming perfect parallel efficiency on 64 cores. Both are crude: the cost per candidate varies widely (on `cond`, 54 budget candidates took 41% of the boxes).
- The script adds no tolerance. The counts are integers produced by the enumerator, whose floating-point steps (the cell masks with slack 1e-9, the HiGHS discard) are the assumptions G2 and G1 in the table of `STATUS.md`.

## Tests

`tests/test_sample_full_size.py` checks the constants against the repository files, the estimator on a toy population with a known total, the tiling counter against a direct count of `tilings`, the cap with fake enumerators, calibration mutations, and that no output carries a host, user or path.
