# Decisions

One line each, dated. Newest last. Only decisions Wöden stated or accepted; proposals go in `docs/open-questions.md`.

- 2026-10-15 — Run the `full` route (all 7,595 triangulations). Hypothesis H is removed by it; `cond` becomes a pilot and a regression set.
- 2026-10-15 — Drop the hull lemma (checklist C1, C2) from the proof path; keep it in the paper as a conjecture with the survey, marked unused.
- 2026-10-15 — Rebuild the checker as a C kernel on FLINT/Arb, with every error bound from a certified source (spec v0.3); the Python reference checker, patched with V1–V5, stays as the oracle.
- 2026-10-15 — Report progress to GitHub during the run (`results/full/progress.json`, hourly, rates only, no ETA). Archive certificates on Zenodo as a draft the script never publishes; only Wöden publishes.
- 2026-10-15 — The range 1 < r < r_min is handled by compactness plus second-order isolation at the balanced 60° configurations (paper, Appendix B), with a band run outside the certified balls.
- 2026-10-15 — Every manual step lives in `docs/setup-guide.md`, updated in the same commit as the code that changes it.
- 2026-10-16 — The `cond` candidate count is 45,260 (distinct codes); 90,836 was the number of gen lines including duplicates and is not to be used.
