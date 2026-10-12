# Status of the results in this repository

What is proved, what is only numerical, and what the completeness computation assumes. Each claim names the file or script that supports it.

Version 2, 2026-10-12 (v2: L6, L8 and L9 closed). This version covers the completeness computation (checklist item A3). The full inventory of the other results (census, co-indices, joins, hull survey, containers, Edwards pilot) is checklist item B2 and will be added here.

## Completeness of the critical census above r = 1 (in progress)

**Claim being worked towards.** Every critical configuration of 12 points on S² with r ≥ r_min = 1.0001 is one of the census configurations.

**What a finished run will prove.** If every candidate is `excluded` or `terminal`, every certificate verifies with `scripts/completeness/check_certificates.py`, and every terminal box is identified with a census point, then the claim holds under the assumptions below.

**Not yet true.** The run has not been done (checklist A4, A5), terminal boxes are not yet identified (C4), and budget cases have no strong test yet (C3).

### Assumptions

The conclusion of the run is conditional on these. Each is a checklist item.

| | assumption | used by | to remove it |
|---|---|---|---|
| H | **Hull conjecture (a).** For θ > 60°, the hull triangulation of a critical configuration is one of T0, T1, T5, T7. | the `cond` route (`run_shard.py --route cond`) | prove it (C1), or run the `full` route over all 7595 triangulations |
| R0 | **Range.** Only r ≥ r_min = 1.0001, i.e. θ ≥ 60.0033°, is covered. | every route | an argument for 1 < r < r_min (C5) |

### Proved

Every rule the search and the checker use, other than H and R0, is proved in `docs/lemmas.md`, by hand or by computer:

- `scripts/completeness/lemmas.py P` — Lemma P (pentagon regions cover all pentagons);
- `scripts/completeness/lemmas.py R` — Lemma R (no rattlers for θ < θ_max; about 27 CPU-minutes);
- `scripts/completeness/lemmas.py U` — Lemma U (Fejes Tóth's area bound at θ_max with its equality case: the icosahedron is the only configuration at θ_max);
- `tests/test_lemmas.py` — checks that these computations are not vacuous.

All three use the checker's arithmetic (`kisscrit/core/interval.py`).
