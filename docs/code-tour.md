# Code tour: from the paper to the code

Purpose: a reading guide for a referee. For each statement of the paper that the computation supports, the file and function that establish it, in the order a reader would follow. Only the **trusted base** is listed; the search is untrusted and proposes what the checker verifies.

Version 1, 2026-10-16 (checklist F9). Update this file in the same commit as any change to trusted code.

## 0. What is trusted

| Component | Files | Status |
|---|---|---|
| Interval arithmetic of the checker | `kisscrit/core/interval.py` (Python reference), `csrc/kc_arith.c` (C kernel) | trusted |
| Face model (corners, vertex sums, Jacobians) | `kisscrit/core/model.py` (Python reference), `csrc/kc_model.c` (C kernel) | trusted |
| Certificate replay | `kisscrit/core/certificate.py` (Python reference), `csrc/kc_replay.c` (C kernel) | trusted |
| Exclusion tests and structural checks | `kisscrit/core/checker.py` | trusted |
| Leaf tests of the C kernel (non-LP exclusions, per-chart driver) | `csrc/kc_leaf.c`, `csrc/kc_leaf.h` | trusted |
| Front end of the C kernel (ctypes, Python fallback) | `kisscrit/ckernel.py` | untrusted for now: not used by `check_certificates.py`. It aggregates the verdict (`verified`), runs the structural checks and the starting-box checks, matches kernel units to reference units and decides the fallback, so it must enter the trusted base, or be re-derived by trusted code, before a later item wires it into the checker. |
| Front end: manifest and summary | `scripts/completeness/check_certificates.py` | trusted |
| Computer-assisted lemmas P, R, U | `scripts/completeness/lemmas.py` | trusted |
| Candidate generation and its pruning | `kisscrit/completeness/candidates.py`, `canon.py` | trusted, with the known gap G1 (floating-point LP discard; checklist F6) |
| Second-order data at FCC and HCP | `scripts/r1/exact_lambda.py` | trusted (exact arithmetic) |
| The search (branch and bound, certificates) | `kisscrit/completeness/fast.py`, `runner.py`, `model.py`, `structure.py`, `scripts/completeness/run_shard.py` | **untrusted** |

Every Python file says which in its module docstring (`Status:` line, enforced by `tests/test_code_standards.py`).

## 0.1 The arithmetic contract and the structural checks (restated from the checker specification v0.2 and the pipeline specification v0.3, internal documents outside this repository, available on request)

Arithmetic: `+ − × ÷ √` on doubles are IEEE correctly rounded and are stepped one ulp outward; `sin, cos, atan, asin, acos` and π are evaluated as Arb balls at 64 bits and rounded outward to doubles; no verdict depends on the system libm, on `-ffast-math`, on FMA contraction, on a floating-point LP result, or on the search's rounding. Minima and maxima of endpoints follow Python's `min` and `max`: the first argument is kept unless the next one compares strictly smaller (larger), so ties between 0.0 and −0.0 and comparisons with a NaN give the same bits in C as in the reference. Read for `kc_hull`, `kc_mul` and the model's `intersect` (2026-10-20, F4 session); `kc_hull` was brought in line then, and `tests/test_carith.py::test_hull_ties_and_nan` compares it by bit pattern on ties and NaNs.

Order of the leaf tests of the C kernel (pipeline spec v0.3 §4.2, an internal document available on request, restated): per leaf or slab, (1) the witness test without tightening, (2) invalid and vertex sums on all vertices without tightening, (3) the linear inclusion with the parent's context and the hint, (4) the tests of (2) with monotonicity tightening, (5) Gordan or Motzkin with a witness y, (6) bisection to the reference's depth and piece limits, (7) if still open, a fallback box. In certificate format v2 the only witness is the multiplier hint, so (1) and (3) are the one hint test, run where `check_certificate` runs it; (5) has no witness until format v3 (checklist F5) and is left to the fallback, as is the linear test with an LP-chosen μ.

Structural checks V1–V5 and V1b (L31–L33) are listed in `docs/certificate-format.md` §4; V6 (the checker imports search code only as an untrusted oracle) is `tests/test_checker_imports.py`; V7 (the checker reports its own bisections and fallbacks) is `check_certificates.py --summary`; V8 (check results are pushed incrementally) is `scripts/completeness/mac_run.sh` phase C.

## 1. Arithmetic (paper Appendix D)

Two implementations, bit-for-bit identical on every tested input (`tests/test_carith.py`): the Python reference and the C kernel's.

- `csrc/kc_arith.c`, `csrc/kc_arith.h` (checklist F2): the same contract in C11 on FLINT/Arb. `kc_init` asserts round-to-nearest and the absence of FMA contraction (the product 0.1·10 − 1 is 0 only without contraction) and computes the constants; `kc_down`/`kc_up` round a ball to doubles through `arb_get_lbound_arf`/`arb_get_ubound_arf` at 64 bits and `arf_get_d` with floor/ceil; `op_down`/`op_up` step IEEE results one ulp outward; `kc_cos` enumerates the multiples of π that may lie inside the interval with the same margins as the reference and decides containment only by Arb comparisons (for |x| > 10⁹, which the model never produces, it returns [−1, 1]); `kc_sin` is cos shifted by π/2 on an outward-rounded interval; `kc_atan`, `kc_asin`, `kc_acos` are monotone ranges from the endpoints. Built by `csrc/Makefile` with `-O2 -fno-fast-math -ffp-contract=off`. `csrc/test_arith.c` checks containment against Arb at 256 bits on 10⁶ random intervals per function plus the special endpoints (measured 2026-10-16: 0 failures).
- `kisscrit/core/carith.py`: the ctypes loader (untrusted; it only moves numbers).


- `kisscrit/core/interval.py`: an interval is a pair of doubles. `+ − × ÷` and `sqrt` use IEEE double arithmetic (correctly rounded) and step one ulp outward (`_op_down`, `_op_up`). `sin`, `cos`, `atan`, `asin`, `acos` and π go through Arb at 64 bits and are rounded outward (`down`, `up`, `_mono`). `cos` widens to ±1 where a multiple of π may lie inside the interval.
- Nothing in the trusted base calls the system libm.

## 2. The face model (paper Appendix A, Lemmas A.2, A.6–A.8; `docs/lemmas.md` L3, L6–L8, L10–L11)

- `kisscrit/core/model.py`: forward-mode AD over intervals (`var`, `add`, `mul`, `div`, `chain`, …), the clamped functions (`clip`, `asin_angle`) whose interval Jacobians enclose the Clarke generalized Jacobian, and the faces:
  - `alpha_range`: α(θ) = arccos(cos θ/(1 + cos θ)) (equilateral corner, A.2);
  - `beta`, `Graph._beta_pt`: base angle of the isosceles triangle with legs θ (rhombus corners, A.7);
  - `pentagon`: the corner formulas of Lemma A.8 in the chart (θ, f₁, f₂), with the pentagon regions a, b, c (`SHIFT`, `LSTRAIGHT`);
  - `_tighten`: monotonicity tightening of pentagon corners (optional; `tighten=False` skips it);
  - `Graph.corners`: all corners on a box with the validity flag `ok` (false only when the box certainly contains no valid configuration: corners outside [α, π], invalid faces);
  - `Graph.vertex_sums`: Σ corners − 2π per vertex with gradients; `Graph.directions`: edge directions at each vertex.
  - The domain: `THETA_FLOOR`, `THETA_CEIL` (1.0 and 1.5): every enclosure above is valid only for θ in this range, where sin θ, cos θ > 0 and α, β are monotone in θ (L17; the places that use it are listed there). The checker enforces it in V3.
- `csrc/kc_model.c`, `csrc/kc_model.h` (checklist F3a): the same model in C11 on `kc_arith`, a line-for-line reimplementation of `model.py` (same function names, same interval operations in the same order, so the same doubles). `kc_graph_new` rebuilds the rotation systems and edges from the oriented faces (`Graph.__init__`, with the same errors); `kc_corners`, `kc_vertex_sums`, `kc_directions` are `Graph.corners`, `Graph.vertex_sums`, `Graph.directions`; `alpha_arb` and `beta_pt` repeat the Arb operations of `_alpha_arb` and `_beta_pt` at 64 bits; `KC_LSTRAIGHT` is the reference's double. Fixed sizes: at most 64 vertices, 128 faces, degree 32, faces of size 3–5 (outside these, `KC_ERR_LIMIT`). `tests/test_cmodel.py` compares the IEEE bit patterns of every corner value and gradient, the validity flag, the vertex sums with their Jacobian rows and the edge directions, on charts of the `cond` certificates: the a-priori boxes, random sub-boxes, point boxes with the clamp bound from an enclosing box, and the first replay boxes. The replay in C is section 3.

## 3. Certificates (`docs/certificate-format.md`)

- `kisscrit/core/certificate.py:replay`: rebuilds every box of the proof tree from the event stream, bit for bit (midpoint splits, grid contractions). It uses no other floating-point computation.
- `csrc/kc_replay.c`, `csrc/kc_replay.h` (checklist F3b): the same replay in C11, as an iterator: `kc_replay_init` takes one chart's decoded streams (events, q, mu), its starting box and caller-owned workspace (`kc_replay_stack_need`, `kc_replay_work_doubles`); each `kc_replay_next` returns the next node in preorder (leaf, split, or a contraction round with its new box, removed slabs and hints), then `KC_REPLAY_DONE` or an error code at the point where the reference raises (each ValueError message, and IndexError for a q read past the end or a coordinate k ≥ n). No globals, no allocation, no elementary functions: `grid_value` and `midpoint` are the format's two IEEE expressions, written as separate statements (built with `-ffp-contract=off`). `tests/test_creplay.py` compares node kinds, the IEEE bit patterns of every box, slab coordinates and sides, hints and the way the replay ends, against `certificate.replay`: on 1,000 sampled `cond` certificates, on hand-made charts for every raise site, on truncated, extended and randomly edited real streams, and with compiled mutants of `kc_replay.c` (an `fma()` grid value, another grid formula, two other midpoints, the upper slab first), each of which the comparison detects.

## 4. The checker (paper Appendix D; `docs/lemmas.md` L15–L16, L31–L33; checker spec v0.2 V1–V7, and V1b)

`kisscrit/core/checker.py`, in the order of a check:

1. `structural_check`:
   - V3 `theta_domain_ok` and `theta_claim_ok`: θ_lo, θ_hi finite with THETA_FLOOR ≤ θ_lo < θ_hi ≤ THETA_CEIL, the model's domain (L17; outside it the a-priori box and `alpha_range` can be inverted and a box excluded on that), then θ_lo ≤ 2 arcsin(r_min/(1 + r_min)), θ_hi ≥ arctan 2, all decided in Arb;
   - V4 `canonical_code` (an independent reimplementation of the canonical code) and `faces_match`;
   - V2 `required_modes`: the pentagon regions that L10′ and L11′ require (the region definitions and the cover proof are L10′);
   - V1 `apriori_box`: the a-priori box of each chart from L3, L4, L7 and L10′, rounded outward; the certificate's box must contain it.
   - V1b (`docs/lemmas.md` L31–L33): every start-box coordinate is finite, a rhombus coordinate lies in (0, π] and a pentagon coordinate in [0, 2π] (hygiene). For a rhombus box outside (0, π] the model's range `ex` of the corner O is false (L32), which let a forged certificate verify at cc375bc (reported by DES-F2, reproduced by the DOM-1 session). L31–L33 are proved by DES-F2 by reading; AUD-2 reviewed them on 2026-10-23, and they are listed again for the next audit because the status lines in this file and in `STATUS.md` disagreed (this line was brought in line on 2026-10-23, DOCS-4).
2. `check_certificate`: replays each chart; every excluded leaf and every contraction slab goes to `Checker.exclude` / `exclude_slab`; terminal leaves are collected; an open leaf fails the certificate.
3. `Checker.exclude`: cheap pass without tightening (`test_invalid_vsum`, `test_linear`), then with tightening, then `test_gordan`, `test_motzkin`, then bisection (bounded by `MAX_DEPTH`, `MAX_PIECES`). A box no test excludes is reported, never accepted.
   - `test_invalid_vsum`: the model's validity flag (L3, L6–L8, L10′) and the vertex-sum intervals (L15);
   - `test_linear`, `test_linear_combo`, `test_linear_mu`, `test_linear_bound`: the mean-value inclusion 0 ∉ μ·F(m) + Σ_j (μᵀJ(B))_j (B_j − m_j); floating point only chooses μ (`combination` is the interval proof);
   - `test_gordan`, `test_motzkin`: Gordan's and Motzkin's alternatives (L16) with `_edge_terms` the interval edge terms in the per-vertex frames of `model.Graph.directions` (L16, "edge directions and frames"); floating point only chooses y.
   - One import from untrusted code: `test_linear_bound` imports `kisscrit.completeness.fast` and uses its simplex only to propose a multiplier vector μ for a contraction slab; the proof is the interval `combination` with that μ, and any real μ is admissible there, so a wrong μ can only make the test fail (`tests/test_checker_imports.py` keeps this the only such import).
   - Terminal leaves are not exclusions: `check_certificate` records them and checks that they are no wider than the certificate's `tol`; completeness needs them identified (checklist C4).
4. `scripts/completeness/check_certificates.py`: `manifest` (V5: gen = cert = certificates = checks, no duplicates, every code verified or terminal-and-identified) and `summary` (V7: the checker's own bisections and fallbacks).

5. `csrc/kc_leaf.c`, `csrc/kc_leaf.h` (checklist F4, spec v0.3 §4.2): the exclusion tests of `Checker` that need no linear program, in C, bit for bit:
   - `kc_leaf_invalid_vsum` is `test_invalid_vsum` on `G.corners(lo, hi, tighten=0|1)`;
   - `kc_leaf_context` and `kc_leaf_hint` are `linear_context` and `test_linear_mu` with `combination` (the hint μ from the certificate, rows in increasing order, the same `Y` and summation order);
   - `kc_leaf_exclude` is `exclude` without `test_linear`, `test_gordan`, `test_motzkin`: untightened and tightened invalid/vsum, then bisection with the same `MAX_DEPTH`, `MAX_PIECES`, split coordinate (first maximum width) and order of the halves;
   - `kc_check_chart` is the per-chart loop of `check_certificate`: it replays the chart with `kc_replay` and gives every leaf and slab a verdict (closed by a named test, unclosed, terminal within or beyond `tol`, open).
   The left-out tests choose μ or y by a floating-point LP (HiGHS, or the search's simplex) and format v2 has no Gordan/Motzkin witness, so a box they would be needed for is returned unclosed; `kisscrit/ckernel.py` (untrusted) hands it to the unchanged reference calls. On every box where the reference returns without raising, the kernel closes it only if the reference does (an argument about the code, in `kc_leaf.h`): the reference runs the same tests in the same order with others interleaved, and visits a subset of the kernel's bisection nodes in the same order, so its piece counter is never larger. The reference also evaluates the model at the centre of every box it tests, and that evaluation can raise (`model.intersect`, a signal of an enclosure bug); the kernel does not, so it could close a box on which the reference raises. `diff_kernel.py` counts a chart on which the reference raises as rejecting all its units, so such a case appears as "C accepts, reference rejects". That the kernel's tests equal the reference's bit for bit is computed, not proved; a referee reads soundness from `kc_model.c` and `kc_arith.c` directly. `scripts/completeness/diff_kernel.py` compares kernel and reference unit by unit (verdicts, and each test's result as bit patterns), and `tests/test_ckernel.py` runs a small sample of it plus compiled mutants of `kc_leaf.c`.

Mutation tests: `tests/test_certificates.py` (each named for its gap); `tests/test_checker_imports.py` (V6).

## 5. Lemmas by computation (paper Appendix A, Lemmas A.8 and A.9)

`scripts/completeness/lemmas.py`, all in the checker's arithmetic, by exhaustive bisection (`run`):
- `lemma_P`: no two adjacent corners of a face pentagon are ≥ 155° (L10);
- `lemma_R`: no rattlers — reduction to circumcenters (`no_deep_hole`, `triple_excluded`), the vertex chart (`pentagon_by_vertices`), test E near the regular pentagon (L9);
- `U_problem` / `lemma_U`: Fejes Tóth's area bound at θ_max with its equality case (L9, uniqueness of the icosahedron).

`tests/test_lemmas.py` checks that these computations are not vacuous.

## 6. Candidates (paper Appendix A, proof of Corollary 6.3; L12–L14)

- `kisscrit/completeness/candidates.py`: `tilings` merges the triangles of each triangulation into faces of 1–3 triangles, pruned by `_vmask` (L3, L4, L7, L13); `lp_feasible` is the linear relaxation (L14) — its floating-point discard is gap G1 until exact Farkas certificates are logged (F6).
- `kisscrit/completeness/canon.py`: canonical codes for deduplication (L12); the checker recomputes them independently (V4).

## 7. The range 1 < r < r_min (paper Appendix B)

- `scripts/r1/exact_lambda.py`: Proposition B.5 — rank 20 contact graphs, one nontrivial flex, per-edge second-order values along the flex: 1/6 at FCC (λ = 1/6), {1/12, 8/33, 3/22} at HCP (λ ≥ 1/12), exactly in Q(√2) and Q(√2, √3). Test: `tests/test_r1.py`.
- `scripts/r1/exact_lambda_stress_space.py` (untrusted, control item C5a): an independent exact check of the Proposition B.5 sentence on the whole self-stress space. Everything is decided in K = Q(√2, √3) (sympy `QQ.algebraic_field`, DomainMatrix; no `simplify`). It recomputes the contact graph, R (24 × 24), a basis of the self-stresses (left null space of R) and the unit flex; gets the per-edge values q_e two ways (the stress Laplacian of Lemma B.2(b), and the second t-derivative of the closed form of φ_e(exp(t f))); and exhibits a witness y with R y = q − λ(1, …, 1), verified by exact multiplication, so q − λ·1 is orthogonal to every self-stress. It also enumerates the extreme rays of the nonnegative stress cone, with an exact sign test. Part D holds the mpmath cross-checks (floating point, not proof). Test: `tests/test_r1_stress_space.py` (exact claims at FCC and HCP, the uniform weight is a stress at FCC and not at HCP, a strictly positive HCP stress, five mutants that must fail).
- The rest of Appendix B (constants, δ per stratum, the band run) is checklist F7 and not yet in code.

## 8. What the search is for

`kisscrit/completeness/fast.py` and `runner.py` find the proof trees quickly and write certificates (format v2). Their rounding (libm with ulp padding) and their LPs affect only how often the checker has to bisect. Soundness never depends on them.

## 9. Exact data of the critical configurations found so far (EXACT-1; not part of the proof)

- `scripts/exact/` (all untrusted): `kfield.py` (K = Q[t]/(P), Sturm isolation, exact signs, linear algebra, inertia), `hiprec.py` (guesses only: Newton and LLL), `verify.py` (Gram checks, stress space, motions and co-index, symmetry, canonical graph, s and r), `build_e2.py` and `build_e3.py` (the eleven configurations), `elim_system.py` (symmetric-subspace KKT system for msolve), `weights_e4.py` (stress cones from the stored data), `table_md.py`. Data: `data/exact/critical-configs-n12.json` and `.md`. Description and conventions: `docs/exact-configs.md`. Test: `tests/test_exact_configs.py` (re-verifies every exact claim from the stored strings; five mutants must fail).

## 10. The run scripts (checklist F8; untrusted tooling, not part of the proof)

`scripts/full/run.sh` (version 1, F8 milestone M1) drives the completeness run on a Linux server without root. It runs the untrusted search (`scripts/completeness/run_shard.py`, unchanged) and then the configured checker, by default the reference checker `scripts/completeness/check_certificates.py` (unchanged). It decides nothing itself. Every verdict it reports is labelled "checker verdict at commit X; not certified", where X is the last commit at the pinned commit that touched `kisscrit/core/` or the checker script. The C kernel front end `kisscrit/ckernel.py` is not used (audit finding B-8). A different checker can be named with `start --checker SCRIPT`. It must be a tracked file of the pinned tree with the same command line, and the run record names it.

- **Phases:** gen, then certify (shards run concurrently), then check (one certificate file at a time), then manifest (V5) and summary (V7).
- **Resume:** each phase skips finished work: gen markers, the `done-*` files of `run_shard.py`, and check markers.
- **One run at a time:** `flock` on `~/kisscrit-run/run.lock`, inherited by the shard processes, so a new copy waits until stragglers of a killed one exit.
- **Pinned code:** the code is a tag or full hash, checked out on the branch `run/<route>-<run-id>`. A tree whose code differs from that commit (`git status` outside `results/`) is refused at start and at every resume. Commits on the run branch touch only `results/`; the certificate files are never committed.
- **Scratch:** under `~/kisscrit-run` (`TMPDIR`, numba cache, bytecode cache). Hashes use Python `hashlib`.
- **Survival:** a `systemd --user` unit when lingering is on, else a detached tmux session plus a tagged `@reboot` crontab line (`start --boot`; a run stopped by hand stays stopped).

`scripts/full/run_state.py` is the helper for the run's JSON files: configuration, run record, progress. It also writes `results/full/progress.json` (rates only, no time estimate), computes the worker cap from the self-test's peak RSS (`workers()`), measures that peak (`rss`) and records the package versions.

`scripts/full/run_repair.py` runs before every resume. It handles what a kill can leave behind:
- a partial last line in a JSON-lines result, certificate or check file;
- a truncated gzip member at the end of a gen file;
- a certificate without its result line, or a result line without its certificate. `run_shard.py` writes the certificate first, so a kill between the two writes leaves this pair, and the resumed run would otherwise append a duplicate certificate.

It removes those lines into a quarantine directory and logs each action. It never adds or changes a result.

Tests: `tests/test_run_scripts.py` (dry runs of every subcommand, refusals, a route-test run killed and resumed, tmux start and stop, the systemd unit through a fake `systemctl`, the self-test with a pytest subset). Manual steps: `docs/setup-guide.md` Part B.

## 11. Moved to `attic/` (not part of the tour; nothing here is trusted or used)

`attic/` keeps files that no code imports or runs, moved with `git mv` so their history stays. The tour does not describe them.

- `attic/kisscrit/completeness/bnb.py`, `contract.py`, `ad.py` (CLEAN-6, 2026-10-22): superseded search prototypes (checklist B1 "Delete"). Before the move, an AST scan of the 96 tracked Python files found no importer of `bnb`, `contract` or `ad` except `contract.py` importing `ad.py`; the full suite passes after the move. `bnb2.py` (read by `scripts/completeness/selftest.py`) and `faces.py` (imported by `structure.py`) stay in `kisscrit/completeness/`.
- `attic/scripts/completeness/lpfilter.py`, `rattler_check.py` (CLEAN-7, 2026-10-22): unused scripts. The AST scan of the tracked Python files and a grep of all tracked text outside `results/` and `data/` find no importer and no reference by name. Every other file in `scripts/completeness/` stays; `candgen2.py` waits for checklist F6.
