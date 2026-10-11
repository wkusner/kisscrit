# Completeness above r = 1: runbook

This is the operating manual for the computation that is meant to show the census of critical configurations of 12 points is complete for r ≥ r_min. It covers what the code does, how to run it on one machine, and how to spread it across several machines using only GitHub.

Version 0.2, 2026-10-11. Revise this when the certifier changes, or when the first full generation run gives the real candidate count.

## 1. Status: read this first

- **What a finished run gives.** It is a computer-assisted *survey*, not yet a proof. Every candidate stressed graph gets one of three outcomes:
  - `excluded` — no critical configuration with this stressed graph and r ≥ r_min;
  - `terminal` — small boxes are left that may contain a critical point, and these still have to be matched against the census;
  - `budget` — not resolved within the box budget.

  The completeness theorem needs zero `budget` results and every `terminal` box identified.
- **Pilot results.**
  - On a sample of T0 candidates and one generic triangulation: about 98% `excluded`; median under 1 s; a typical mean of about 3 CPU-s per candidate; 1–2% unresolved at a budget of 100k–300k boxes.
  - The unresolved ones are mostly graphs with few edges (E = 17–19). Their realization space has dimension 22 − E = 3–5. They need a stronger exclusion test, a first-order (Taylor-model) Gordan certificate, which is not written yet.
- **The low end r → 1 is excluded by design.** The default is r_min = 1.0001, i.e. θ ≥ 60.0033°. The smallest known critical value is 1.0002111. Below r_min, many graphs become nearly realizable, because all the kissing configurations live at r = 1, and box counts blow up. The range (1, r_min) needs a separate perturbation argument from r = 1.
- **Two routes.**
  - `cond` uses only triangulations 0, 1, 5, 7, which assumes Conjecture 7.3(a) (hull triangulations).
  - `full` uses all 7595 triangulations and is unconditional.

  Run `cond` first. Its generation phase gives the candidate count, which sets the real cost.

## 2. What the code does

There are two phases, and both are sharded, resumable and append-only.

**Phase A, `gen`.** For each triangulation T of the sphere with 12 vertices (`data/completeness/triangulations12.json`):
1. Enumerate the partitions of its 20 triangles into faces of 1, 2 or 3 triangles (3-, 4- and 5-gons), pruned by the vertex angle conditions.
2. Build the plane graph and compute its canonical code (`kisscrit/completeness/canon.py`; the code is invariant under relabelling and reflection).
3. Keep each new code that passes a linear relaxation built from facts valid for every critical configuration.

The output is `results/completeness/gen/gen-<route>-<k>-of-<K>.jsonl.gz`.

**Phase B, `certify`.** This phase deduplicates the candidates from all gen files and takes those whose code hashes to shard k. For each candidate, `kisscrit/completeness/runner.py` runs an interval branch and bound over face-angle parameters (θ, one corner per quadrilateral, two per pentagon), with three regions per pentagon. A box is discarded when it is:
- `invalid` — not a convex equilateral polygon in this chart;
- `vsum` — some vertex angle sum cannot equal 360°;
- `linear` — infeasible for the linearized vertex equations; the LP multipliers are re-checked with interval arithmetic;
- `gordan` — a Gordan or Motzkin certificate shows that no strictly positive balanced stress exists on the whole box.

Otherwise the box is contracted and bisected. The output is `results/completeness/cert/cert-<route>-<k>-of-<K>.jsonl`, one line per candidate with status, statistics, host and seconds. The line also carries the terminal boxes, when there are any.

**Rigour.**
- All exclusions use outward-rounded interval arithmetic: one ulp per arithmetic operation, and ≥ 4 ulp plus 1e-15 for each libm call.
- Floating-point LP solutions are used only to choose multipliers.
- `kisscrit/completeness/model.py` is a slower reference implementation of the same mathematics. `selftest.py` checks the compiled kernels against it.
- Lemma P (no two adjacent near-straight pentagon corners) is proved by `scripts/completeness/lemma_pentagon.py`.

## 3. One-time setup on each machine

Requirements: Python ≥ 3.9 and git. Use Linux, macOS, or Windows with WSL. On native Windows, everything except `claim.py` works.

```
git clone https://github.com/<you>/<repo>.git
cd <repo>
python3 -m venv .venv && source .venv/bin/activate      # Windows: .venv\Scripts\activate
pip install numpy scipy networkx mpmath numba
python3 scripts/completeness/selftest.py                 # 1-3 min; must end with "selftest OK"
```

For `claim.py` to work, git must be able to push to the repository from this machine without a password prompt. Use the GitHub CLI login, a credential helper, or an SSH key.

Find the core count with `nproc` (Linux), `sysctl -n hw.ncpu` (macOS), or `echo %NUMBER_OF_PROCESSORS%` (Windows).

## 4. Running on one machine

```
# Phase A: the cond route is only 4 triangulations; 4 shards = one triangulation each
for k in 0 1 2 3; do python3 scripts/completeness/run_shard.py gen --route cond --shard $k --of 4 & done; wait

# Phase B: one shard covering everything, all cores
python3 scripts/completeness/run_shard.py certify --route cond --shard 0 --of 1 --workers $(nproc)

python3 scripts/completeness/merge.py --route cond
```

To stop a run, interrupt it. To resume, rerun the same command: finished work is skipped.

## 5. Distributing across machines with GitHub

The repository is the only coordination point. Each shard has its own output file, so two machines never write the same file, and merges never conflict. A shard is claimed by committing a one-line file under `results/completeness/claims/`.

**Choose K once, before starting.** K is the number of shards, and it must stay fixed for the whole run. Pick roughly 5–10 × the number of machines, e.g. K = 40, so that fast machines simply take more shards. Record K in the commit message of the first claim.

### Phase A (gen)

For `cond` there are only 4 shards (`--of 4`). On one machine:

```
for k in 0 1 2 3; do python3 scripts/completeness/run_shard.py gen --route cond --shard $k --of 4 & done; wait
git add results/completeness/gen && git commit -m "gen cond" && git push
```

For `full`, distribute gen exactly like certify below, using `--phase gen` and `--of 40`.

### Phase B (certify): the loop to run on each machine

```
NAME=lawrence-laptop          # any unique name per machine
K=40
while K_SHARD=$(python3 scripts/completeness/claim.py --phase certify --route cond --of $K --name $NAME); do
  python3 scripts/completeness/run_shard.py certify --route cond --shard $K_SHARD --of $K --workers $(nproc)
  git pull --rebase --quiet
  git add results/completeness/cert && git commit -m "certify cond shard $K_SHARD/$K ($NAME)" && git push
done
```

`claim.py` pulls first, takes the lowest unclaimed shard, and commits and pushes the claim. If another machine pushed first, it backs out and tries again. When no shards are left it exits with status 1, which ends the loop.

To check progress from anywhere:

```
git pull && python3 scripts/completeness/merge.py --route cond
```

**If a machine dies mid-shard:**
1. Commit whatever partial `cert-…jsonl` it left, if you can reach it.
2. Delete its claim file and push.
3. The next machine to claim that shard resumes from the partial file.

**To push long-running partial results** without waiting for the shard to finish, run this from a second terminal now and then:

```
git add results/completeness/cert && git commit -m wip && git pull --rebase && git push
```

**Re-running `budget` cases with a larger budget.** Rerun the same shard with, e.g., `--budget 3000000`. Only `budget` records with a smaller budget are redone. `merge.py` keeps the best record for each candidate.

## 6. Cost estimates

Measured here, on 2 cores:
- **Generation:** about 7 ms per tiling. A generic triangulation yields 10³–10⁵ tilings. T0 yields 1.2 × 10⁶, so it takes about 2.3 h on one core.
- **Certification:** about 3 CPU-s per candidate on average, with a long tail.

Total certify cost ≈ (number of unique candidates) × 3 CPU-s. That gives about 1 CPU-hour per 1,200 candidates.

The unique candidate count for `cond` is not known until Phase A finishes. Expect 10⁵–10⁶, so the certify phase is roughly 80–800 CPU-hours. On an 8-core laptop that is 10–100 h of wall-clock. A 64-core server takes 1.3–13 h.

The `full` route is at least an order of magnitude larger and is not practical until the certifier improves.

**Claude cloud sessions are a poor fit for this batch work.**
- The container in this session has 2 cores.
- The container is reclaimed when the session goes idle.
- The credits pay for model usage, not compute.

Use them for development (improving the certifier, identifying `terminal` boxes, writing the proof), not for grinding shards.

**The departmental LLM is not needed.** The run is plain CPU work with no model calls. The machine it runs on is useful only as CPU.

## 7. After the run: what to send back

Push everything, then run `merge.py`. The file `results/completeness/summary-<route>.json` lists every `terminal` and `budget` candidate. Those are the inputs to the next step:
- identify each terminal box with a census point (the existing X-space Krawczyk certification);
- attack the `budget` graphs with the stronger certificate.
