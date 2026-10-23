# Status of the results in this repository

What is proved, what is only numerical, and what the completeness computation assumes. Each claim names the file or script that supports it.

Version 3.2, 2026-10-23 (v3.2, DOCS-1: the assumptions table gains G1, G2, G3, C4, tol and V4(b); "Proved" names each file and section; the hypotheses H1–H4; the pilot's count differences explained; the check records of the iMac run. v3.1, 2026-10-22: sections added at the end for the server pilot of 2026-10-22, the DES-F2 findings and the EXACT-1 data; the gen count corrected, CLEAN-8. v3, 2026-10-16: the `full` route replaces hypothesis H; `cond` is a pilot; the candidate count resolved. v2: L6, L8 and L9 closed). This version covers the completeness computation (checklist item A3). The full inventory of the other results (census, co-indices, joins, hull survey, containers, Edwards pilot) is checklist item B2 and will be added here.

## Completeness of the critical census above r = 1 (in progress)

**Claim being worked towards.** Every critical configuration of 12 points on S² with r ≥ r_min = 1.0001 is one of the census configurations.

**What a finished run would prove, under hypotheses H1–H4.** The claim follows from four hypotheses (the conditional statement is lemma L27 of the DES-F2 C4 design document, 2026-10-21, available on request; proved there as a conditional statement and re-derived by the AUD-2 review, 2026-10-23):

- **H1.** The candidate list is complete.
- **H2.** Every candidate has a certificate that `scripts/completeness/check_certificates.py` accepts with no undecided leaf, and every unexcluded leaf is either excluded by Krawczyk emptiness or identified with a census configuration.
- **H3.** The census ball certificates and the isolation constants at DOD and M verify.
- **H4.** The manifest ties every check and every identification to its certificate.

These are hypotheses, not results. Today none of H1–H4 is met: the `full` candidate list does not exist (checklist F6), terminal boxes are not identified (C4), 54 (iMac run) or 56 (server pilot) candidates are over budget (C3), and the manifest ties codes, not file contents (DOM-3). The assumptions table below lists what the present computation rests on.

**Not yet true.** The `cond` run is done but conditional on H (below); terminal boxes are not yet identified (C4), and budget cases have no strong test yet (C3). The `full` run has not been done (checklist F).

**State of the `cond` run (measured, 2026-10-15).** 45,260 candidate codes (the gen files hold 90,836 code lines, 90,840 with four end markers; computed, `zcat` and a count over `results/completeness/gen/gen-cond-*-of-4.jsonl.gz` at e4ae530. A code is written once for each of the triangulations T0, T1, T5, T7 that produces it: 9,140 codes appear once, 26,664 twice and 9,456 three times (computed, same files). An earlier version of this sentence said every code is generated twice across the four shards; that was wrong. 45,260 is the count of distinct codes, and it is the number to use). Results: 44,740 `excluded`, 466 `terminal`, 54 `budget`; 39.5 M boxes and 9.3 search CPU-hours, of which the 54 budget candidates used 16.2 M boxes and 3.9 h. Check records (computed at 111ca4b in the DOCS-1 session, a count over `results/completeness/checks/checker-2-reference-*.jsonl` and `results/completeness/cert/`): the reference checker's `checker-2-reference` records cover all 45,260 `cond` certificates: 0 verified, every one a structural failure (V3 in all 45,260, V1 in 39,608; the one-ulp start-box defect of the search that f0a32cb fixed after this run, checklist F1), 0 failed boxes, 4,673 terminal boxes on 467 codes (466 `terminal` codes plus the one `budget` code with terminal boxes, the E = 19 census candidate with 3,761). These records are a box-level regression set, not a verdict on the run. The 35 older `checks-certificates-cond-*` files are `checker-1` records (the old Python checker; informational). The key `certified` (45,260) in `results/completeness/summary-cond.json` counts certificates written by the search, not a checker verdict; the key is not renamed. The claim will rest on `checker-2` (the C kernel with Python fallback and V1–V5, checklist F1–F4).

**Decision (2026-10-15).** The `full` route over all 7,595 triangulations will be run, which removes H. `cond` is kept as a pilot and a regression set. See `docs/decisions.md`.

### Assumptions

The conclusion of the run is conditional on these. Each is a checklist item.

| | assumption | used by | to remove it |
|---|---|---|---|
| H | **Hull conjecture (a).** For θ > 60°, the hull triangulation of a critical configuration is one of T0, T1, T5, T7. | the `cond` route only (`run_shard.py --route cond`) | the `full` route over all 7,595 triangulations (checklist F); C1 is dropped |
| R0 | **Range.** Only r ≥ r_min = 1.0001, i.e. θ ≥ 60.0033°, is covered. | every route | the paper's Appendix B (compactness and second-order isolation at the balanced 60° configurations, exact at FCC and HCP in `scripts/r1/exact_lambda.py`) plus a band run with δ-ball exclusion (C5, F7) |
| G1 | **LP discard in floating point.** `kisscrit/completeness/candidates.py:lp_feasible` discards a tiling when the HiGHS LP (`scipy.optimize.linprog`) reports it infeasible; that decision is floating point. | candidate generation, every route | exact rational Farkas certificates for every discard (checklist F6, G1) |
| G2 | **Vertex masks in floating point.** `candidates.py:_vmask` decides the admissible vertex types in floating point. | candidate generation, every route | an exact vertex-mask table (F6, G2) |
| G3, L12 | **Triangulation list.** `data/completeness/triangulations12.json` holds 7,595 triangulations (counted); its content is not re-derived in this repository (no generator here), and `candidates.py` imports the untrusted `kisscrit/completeness/structure.py`. | the `full` route (all triangulations); the `cond` route for T0, T1, T5, T7 | the G3 set comparison and the L12 cross-check (F6) |
| C4 | **Terminal identification not built.** The checker accepts a terminal leaf as "left for identification"; nothing yet identifies terminal leaves with census configurations. | every route | C4a and C4b (per-leaf identification, lemmas L18–L26 of the DES-F2 design, not yet in `docs/lemmas.md`) |
| tol | **Terminal width trusted.** The checker accepts a terminal leaf whose width is at most the certificate's own `tol` (`kisscrit/core/checker.py`, the terminal-leaf branch of `check_certificate`); `tol` is read from the certificate. | the checker | per-leaf identification that does not use `tol` (DES-F2 lemma L26, with C4) |
| V4(b) | **Gen files needed.** Check V4(b) (gen faces match the certificate's faces) runs only when the gen files are present; without them the code identity is not checked (`checker.py:check_certificate`). The manifest ties codes, not file contents. | the checker, the manifest | gen files with every check, and manifest hashes (DOM-3) |

### Proved

The exclusion tests and the structural checks V1–V4 and V1b are proved in `docs/lemmas.md`, by hand or in the checker's arithmetic:

- structure of the stressed graph: `docs/lemmas.md` §A, L1–L9 (L9, no rattlers, uses Lemma R and Lemma U below);
- pentagon charts: §B, L10 (Lemma P), L10′, L11, L11′;
- the rules that candidate generation applies: §C, L12–L14 (the rules are proved; their implementation is not: G1, G2, G3 in the table);
- the exclusion tests and the domain of the model's enclosures: §D, L15–L17;
- the domain check V1b: §D, L31–L33.

The computer-assisted lemmas, in the checker's arithmetic (`kisscrit/core/interval.py`):

- `scripts/completeness/lemmas.py P` — Lemma P (pentagon regions cover all pentagons; `docs/lemmas.md` §B, L10);
- `scripts/completeness/lemmas.py R` — Lemma R (no rattlers for θ < θ_max; §A, L9; about 27 CPU-minutes);
- `scripts/completeness/lemmas.py U` — Lemma U (Fejes Tóth's area bound at θ_max with its equality case: the icosahedron is the only configuration at θ_max; §A, L9, step 3);
- `tests/test_lemmas.py` — checks that these computations are not vacuous.

Not proved: the candidate list (G1, G2, G3), the identification of terminal boxes (C4), and the resolution of budget candidates (C3). See the table.

## Pilot of the `cond` route on the server, 2026-10-22 (computed; not a completeness claim)

Provenance (from the pilot's own report, `code-pilot-report-2026-10-22.txt` in the maintainer's project folder `04_Compute`, available on request; written 2026-10-22 12:47 CDT; the CLEAN-A session did not rerun it): route `cond`; tree commit 322f8bc (before V3-FIX, BALL-KILL and DOM-1); server W227-2, 64 cores; Python 3.12.15, numpy 2.5.4, numba 0.68.0, python-flint 0.8.0; 40 certify shards, budget 300,000; wall seconds: gen 490, certify 360, check 8,042.

Results (computed, same report): generated 45,260; certificates written 45,260 (a count of candidates for which the search wrote a certificate, not a checker verdict); excluded 44,737; terminal 467; budget 56; search 3.29 CPU-hours. The reference checker (`check_certificates.py` at that tree) verified 45,204 and did not verify 56, all 56 structural failures; 0 failed boxes; 740 open leaves; 6,079 terminal boxes. The manifest says `complete: false`, because terminal-box identification (C4) is not built. The report's manifest also gives `terminal_unidentified` 468 against `terminal` 467 in the counts line. This is bookkeeping, not a discrepancy: `terminal_unidentified` counts every code whose check record holds terminal boxes, whatever its search status (`scripts/completeness/check_certificates.py`, `manifest`); the one extra code is the E = 19 census candidate (status `budget`, 5,165 terminal boxes; 6,079 terminal boxes in all) (computed by the AUD-2 audit, 2026-10-23, from the pilot's files; not reproduced here). The iMac run has the same structure: 466 `terminal` codes plus that `budget` code, 467 codes with terminal boxes (computed here, see the check records above). The 717 open leaves of the iMac run (54 budget codes) and the 740 of the pilot (56 budget codes) are measurements of different runs on different candidate sets, not a discrepancy (AUD-2).

Against the iMac run above (44,740 excluded, 466 terminal, 54 budget; 9.3 search CPU-hours) the pilot differs by three candidates (excluded 3 fewer, terminal 1 more, budget 2 more). The cause is commit f0a32cb, which changed the search's start boxes after the iMac run (they now contain π and θ(r_min) exactly). Established for one candidate: the E = 23 candidate (6 triangles, 7 rhombi) that became `terminal` with one box was excluded on the iMac by a Gordan leaf after three contraction rounds at the root; once arctan 2 lies inside the start box it carries a nonnegative stress at θ_max (the icosahedron), so the Gordan test cannot close it (computed by AUD-2 by replaying its iMac certificate; the mechanism inferred). Inferred for the other two: the two E = 19/20 candidates that became `budget` were excluded on the iMac at 273,094 and 271,419 boxes, 91 % of the 300,000 cap, and the wider start boxes need more (which candidates moved is AUD-2's comparison with the pilot's files; the two box counts are in `results/completeness/cert/`, checked here). The CPU-hours are not comparable (other machine, other tree, other Python); not analysed.

## DES-F2 design audit (reported by DES-F2; reviewed by the AUD-2 audit, 2026-10-23)

Correction first: the sentence under "Proved" that every rule the search and the checker use, other than H and R0, is proved in `docs/lemmas.md` was wrong at cc375bc: the checker did not test the domain of a rhombus start-box coordinate. Fix: check V1b (lemmas L31-L33, `docs/lemmas.md`), commit f3e8d0d. L31–L33 are proved (DES-F2; reviewed by AUD-2, 2026-10-23; L32 also re-derived by the DOM-1 session). The "Proved" section above was rewritten in v3.2 and no longer makes that claim.

- The checker at cc375bc accepted a forged certificate with a rhombus start box that leaves (0, π] (DES-F2 document `admin-twelve-fable-design-domain-2026-10-21.md`, maintainer's project folder `06_Admin`, available on request; reproduced by the DOM-1 session at 5da4bba, `docs/session-log.md`, computed there, and by AUD-2 at cc375bc and 5da4bba, scratch). V1b at f3e8d0d rejects it in the reference and the kernel front end (DOM-1 session, computed there; reproduced by AUD-2 at 111ca4b, scratch).
- The premise of checklist C3 is reported wrong: for 53 of the 54 budget candidates of the iMac run, the 300,000 boxes are spent on coarse `vsum` and `invalid` leaves (width 0.24-0.48), not on leaves a stronger Gordan test would close; only the E = 19 census candidate (r = 1.0002111636) is a near-critical cluster (DES-F2 document `admin-twelve-fable-design-C3-2026-10-21.md`, maintainer's project folder `06_Admin`, available on request, scratch computations in floating point; reproduced by AUD-2 by replaying all 54 budget certificates, scratch). C3 is open; decision D-031 replaces its text by R1–R4 (`docs/checklist.md`).

## Exact data of the census configurations (EXACT-1; not part of the proof)

`data/exact/critical-configs-n12.json` and `.md`, described in `docs/exact-configs.md` (EXACT-1 session, 2026-10-21/22; reported there, and `tests/test_exact_configs.py` re-verifies the stored exact claims, passing in the CLEAN-A baseline, computed). Eleven configurations are exactly verified in K = Q(c): DOD, M, r = 1.0010073939, 1.0008452029, B, r1, 1.0003062544, A (eight of the twelve census values above r = 1), and FCC, HCP and the hexagonal antiprism band. Four census values have no exact polynomial: 1.0003062547, 1.0002374095 and 1.0002115891 (no integer relation of degree at most 150 at 8000 digits: high precision, not proof), and 1.0002111636 (a candidate of degree 126, high precision only, not exactly verified). Not a completeness claim. Open question for the stored symmetry of `data/n12_above1.json` entry 6: `docs/open-questions.md`.
