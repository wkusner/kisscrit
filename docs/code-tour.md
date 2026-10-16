# Code tour: from the paper to the code

Purpose: a reading guide for a referee. For each statement of the paper that the computation supports, the file and function that establish it, in the order a reader would follow. Only the **trusted base** is listed; the search is untrusted and proposes what the checker verifies.

Version 1, 2026-10-16 (checklist F9). Update this file in the same commit as any change to trusted code.

## 0. What is trusted

| Component | Files | Status |
|---|---|---|
| Interval arithmetic of the checker | `kisscrit/core/interval.py` (Python reference), `csrc/kc_arith.c` (C kernel) | trusted |
| Face model (corners, vertex sums, Jacobians) | `kisscrit/core/model.py` | trusted |
| Certificate replay | `kisscrit/core/certificate.py` | trusted |
| Exclusion tests and structural checks | `kisscrit/core/checker.py` | trusted |
| Front end: manifest and summary | `scripts/completeness/check_certificates.py` | trusted |
| Computer-assisted lemmas P, R, U | `scripts/completeness/lemmas.py` | trusted |
| Candidate generation and its pruning | `kisscrit/completeness/candidates.py`, `canon.py` | trusted, with the known gap G1 (floating-point LP discard; checklist F6) |
| Second-order data at FCC and HCP | `scripts/r1/exact_lambda.py` | trusted (exact arithmetic) |
| The search (branch and bound, certificates) | `kisscrit/completeness/fast.py`, `runner.py`, `model.py`, `structure.py`, `scripts/completeness/run_shard.py` | **untrusted** |

Every Python file says which in its module docstring (`Status:` line, enforced by `tests/test_code_standards.py`).

## 0.1 The arithmetic contract and the structural checks (restated from the specifications, which are Project documents outside this repository)

Arithmetic: `+ − × ÷ √` on doubles are IEEE correctly rounded and are stepped one ulp outward; `sin, cos, atan, asin, acos` and π are evaluated as Arb balls at 64 bits and rounded outward to doubles; no verdict depends on the system libm, on `-ffast-math`, on FMA contraction, on a floating-point LP result, or on the search's rounding.

Structural checks V1–V5 are listed in `docs/certificate-format.md` §4; V6 (the checker imports search code only as an untrusted oracle) is `tests/test_checker_imports.py`; V7 (the checker reports its own bisections and fallbacks) is `check_certificates.py --summary`; V8 (check results are pushed incrementally) is `scripts/completeness/mac_run.sh` phase C.

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

## 3. Certificates (`docs/certificate-format.md`)

- `kisscrit/core/certificate.py:replay`: rebuilds every box of the proof tree from the event stream, bit for bit (midpoint splits, grid contractions). It uses no other floating-point computation.

## 4. The checker (paper Appendix D; `docs/lemmas.md` L15–L16; checker spec v0.2 V1–V7)

`kisscrit/core/checker.py`, in the order of a check:

1. `structural_check`:
   - V3 `theta_claim_ok`: θ_lo ≤ 2 arcsin(r_min/(1 + r_min)), θ_hi ≥ arctan 2, decided in Arb;
   - V4 `canonical_code` (an independent reimplementation of the canonical code) and `faces_match`;
   - V2 `required_modes`: the pentagon regions that L10′ and L11′ require (the region definitions and the cover proof are L10′);
   - V1 `apriori_box`: the a-priori box of each chart from L3, L4, L7 and L10′, rounded outward; the certificate's box must contain it.
2. `check_certificate`: replays each chart; every excluded leaf and every contraction slab goes to `Checker.exclude` / `exclude_slab`; terminal leaves are collected; an open leaf fails the certificate.
3. `Checker.exclude`: cheap pass without tightening (`test_invalid_vsum`, `test_linear`), then with tightening, then `test_gordan`, `test_motzkin`, then bisection (bounded by `MAX_DEPTH`, `MAX_PIECES`). A box no test excludes is reported, never accepted.
   - `test_invalid_vsum`: the model's validity flag (L3, L6–L8, L10′) and the vertex-sum intervals (L15);
   - `test_linear`, `test_linear_combo`, `test_linear_mu`, `test_linear_bound`: the mean-value inclusion 0 ∉ μ·F(m) + Σ_j (μᵀJ(B))_j (B_j − m_j); floating point only chooses μ (`combination` is the interval proof);
   - `test_gordan`, `test_motzkin`: Gordan's and Motzkin's alternatives (L16) with `_edge_terms` the interval edge terms in the per-vertex frames of `model.Graph.directions` (L16, "edge directions and frames"); floating point only chooses y.
   - One import from untrusted code: `test_linear_bound` imports `kisscrit.completeness.fast` and uses its simplex only to propose a multiplier vector μ for a contraction slab; the proof is the interval `combination` with that μ, and any real μ is admissible there, so a wrong μ can only make the test fail (`tests/test_checker_imports.py` keeps this the only such import).
   - Terminal leaves are not exclusions: `check_certificate` records them and checks that they are no wider than the certificate's `tol`; completeness needs them identified (checklist C4).
4. `scripts/completeness/check_certificates.py`: `manifest` (V5: gen = cert = certificates = checks, no duplicates, every code verified or terminal-and-identified) and `summary` (V7: the checker's own bisections and fallbacks).

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
- The rest of Appendix B (constants, δ per stratum, the band run) is checklist F7 and not yet in code.

## 8. What the search is for

`kisscrit/completeness/fast.py` and `runner.py` find the proof trees quickly and write certificates (format v2). Their rounding (libm with ulp padding) and their LPs affect only how often the checker has to bisect. Soundness never depends on them.
