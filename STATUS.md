# Status of the results in this repository

What is proved, what is only numerical, and what the completeness computation assumes. Each claim names the file or script that supports it.

Version 3, 2026-10-16 (v3: the `full` route replaces hypothesis H; `cond` is a pilot; the candidate count resolved. v2: L6, L8 and L9 closed). This version covers the completeness computation (checklist item A3). The full inventory of the other results (census, co-indices, joins, hull survey, containers, Edwards pilot) is checklist item B2 and will be added here.

## Completeness of the critical census above r = 1 (in progress)

**Claim being worked towards.** Every critical configuration of 12 points on S² with r ≥ r_min = 1.0001 is one of the census configurations.

**What a finished run will prove.** If every candidate is `excluded` or `terminal`, every certificate verifies with `scripts/completeness/check_certificates.py`, and every terminal box is identified with a census point, then the claim holds under the assumptions below.

**Not yet true.** The `cond` run is done but conditional on H (below); terminal boxes are not yet identified (C4), and budget cases have no strong test yet (C3). The `full` run has not been done (checklist F).

**State of the `cond` run (measured, 2026-10-15).** 45,260 candidate codes (the gen files hold 90,836 code lines, 90,840 with four end markers; computed, `zcat` and a count over `results/completeness/gen/gen-cond-*-of-4.jsonl.gz` at e4ae530. A code is written once for each of the triangulations T0, T1, T5, T7 that produces it: 9,140 codes appear once, 26,664 twice and 9,456 three times (computed, same files). An earlier version of this sentence said every code is generated twice across the four shards; that was wrong. 45,260 is the count of distinct codes, and it is the number to use). Results: 44,740 `excluded`, 466 `terminal`, 54 `budget`; 39.5 M boxes and 9.3 search CPU-hours, of which the 54 budget candidates used 16.2 M boxes and 3.9 h. Checker results from the iMac's phase C (old Python checker, labelled `checker-1`) are informational; the claim will rest on `checker-2` (the C kernel with Python fallback and V1–V5, checklist F1–F4).

**Decision (2026-10-15).** The `full` route over all 7,595 triangulations will be run, which removes H. `cond` is kept as a pilot and a regression set. See `docs/decisions.md`.

### Assumptions

The conclusion of the run is conditional on these. Each is a checklist item.

| | assumption | used by | to remove it |
|---|---|---|---|
| H | **Hull conjecture (a).** For θ > 60°, the hull triangulation of a critical configuration is one of T0, T1, T5, T7. | the `cond` route only (`run_shard.py --route cond`) | the `full` route over all 7,595 triangulations (checklist F); C1 is dropped |
| R0 | **Range.** Only r ≥ r_min = 1.0001, i.e. θ ≥ 60.0033°, is covered. | every route | the paper's Appendix B (compactness and second-order isolation at the balanced 60° configurations, exact at FCC and HCP in `scripts/r1/exact_lambda.py`) plus a band run with δ-ball exclusion (C5, F7) |

### Proved

Every rule the search and the checker use, other than H and R0, is proved in `docs/lemmas.md`, by hand or by computer:

- `scripts/completeness/lemmas.py P` — Lemma P (pentagon regions cover all pentagons);
- `scripts/completeness/lemmas.py R` — Lemma R (no rattlers for θ < θ_max; about 27 CPU-minutes);
- `scripts/completeness/lemmas.py U` — Lemma U (Fejes Tóth's area bound at θ_max with its equality case: the icosahedron is the only configuration at θ_max);
- `tests/test_lemmas.py` — checks that these computations are not vacuous.

All three use the checker's arithmetic (`kisscrit/core/interval.py`).
