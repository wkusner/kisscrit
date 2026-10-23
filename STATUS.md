# Status of the results in this repository

What is proved, what is only numerical, and what the completeness computation assumes. Each claim names the file or script that supports it.

Version 3.1, 2026-10-22 (v3.1: sections added at the end for the server pilot of 2026-10-22, the DES-F2 findings and the EXACT-1 data; the gen count corrected, CLEAN-8. v3, 2026-10-16: the `full` route replaces hypothesis H; `cond` is a pilot; the candidate count resolved. v2: L6, L8 and L9 closed). This version covers the completeness computation (checklist item A3). The full inventory of the other results (census, co-indices, joins, hull survey, containers, Edwards pilot) is checklist item B2 and will be added here.

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

## Pilot of the `cond` route on the server, 2026-10-22 (computed; not a completeness claim)

Provenance (from the pilot's own report, Project document `claude/code-pilot-report-2026-10-22.txt`, written 2026-10-22 12:47 CDT; the CLEAN-A session did not rerun it): route `cond`; tree commit 322f8bc (before V3-FIX, BALL-KILL and DOM-1); server W227-2, 64 cores; Python 3.12.15, numpy 2.5.4, numba 0.68.0, python-flint 0.8.0; 40 certify shards, budget 300,000; wall seconds: gen 490, certify 360, check 8,042.

Results (computed, same report): generated 45,260; certificates written 45,260 (a count of candidates for which the search wrote a certificate, not a checker verdict); excluded 44,737; terminal 467; budget 56; search 3.29 CPU-hours. The reference checker (`check_certificates.py` at that tree) verified 45,204 and did not verify 56, all 56 structural failures; 0 failed boxes; 740 open leaves; 6,079 terminal boxes. The manifest says `complete: false`, because terminal-box identification (C4) is not built. The report's manifest also gives `terminal_unidentified` 468 against `terminal` 467 in the counts line; the report does not explain the difference of one (unchecked).

Against the iMac run above (44,740 excluded, 466 terminal, 54 budget; 9.3 search CPU-hours) the pilot differs by three candidates (excluded 3 fewer, terminal 1 more, budget 2 more). The cause is not established. One possible cause, unchecked: commit f0a32cb changed the search's start boxes after the iMac run. The CPU-hours are not comparable (other machine, other tree, other Python); not analysed.

## DES-F2 design audit (reported by DES-F2, not reproduced by the control chat)

Correction first: the sentence under "Proved" that every rule the search and the checker use, other than H and R0, is proved in `docs/lemmas.md` was wrong at cc375bc: the checker did not test the domain of a rhombus start-box coordinate. Fix: check V1b (lemmas L31-L33, `docs/lemmas.md`), commit f3e8d0d. L32 was re-derived by the DOM-1 session; L31 and L33 are proved by DES-F2 by reading and are not yet reviewed by the control chat.

- The checker at cc375bc accepted a forged certificate with a rhombus start box that leaves (0, π] (DES-F2 document `claude/admin-twelve-fable-design-domain-2026-10-21.md`; reproduced by the DOM-1 session at 5da4bba, `docs/session-log.md`, computed there; not reproduced by the control chat). V1b at f3e8d0d rejects it (DOM-1 session, computed there).
- The premise of checklist C3 is reported wrong: for 53 of the 54 budget candidates of the iMac run, the 300,000 boxes are spent on coarse `vsum` and `invalid` leaves (width 0.24-0.48), not on leaves a stronger Gordan test would close; only the E = 19 census candidate (r = 1.0002111636) is a near-critical cluster (DES-F2 document `claude/admin-twelve-fable-design-C3-2026-10-21.md`, scratch computations in floating point; not reproduced by the control chat). C3 is open.

## Exact data of the census configurations (EXACT-1; not part of the proof)

`data/exact/critical-configs-n12.json` and `.md`, described in `docs/exact-configs.md` (EXACT-1 session, 2026-10-21/22; reported there, and `tests/test_exact_configs.py` re-verifies the stored exact claims, passing in the CLEAN-A baseline, computed). Eleven configurations are exactly verified in K = Q(c): DOD, M, r = 1.0010073939, 1.0008452029, B, r1, 1.0003062544, A (eight of the twelve census values above r = 1), and FCC, HCP and the hexagonal antiprism band. Four census values have no exact polynomial: 1.0003062547, 1.0002374095 and 1.0002115891 (no integer relation of degree at most 150 at 8000 digits: high precision, not proof), and 1.0002111636 (a candidate of degree 126, high precision only, not exactly verified). Not a completeness claim. Open question for the stored symmetry of `data/n12_above1.json` entry 6: `docs/open-questions.md`.
