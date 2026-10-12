# Project checklist: completeness of the N = 12 census and the Experimental Mathematics paper

Purpose: the working checklist for this project, in the order we intend to do the items. Each item says what "done" means. After an item is closed, Claude gives a short progress report and review, and asks before starting the next one.

Version 1, 2026-10-11. Update this file as items close; don't keep a separate list.

Status keys: [ ] open · [~] in progress · [x] done · [-] dropped

## A. Certification core — must be done before the long run starts

The long run on the iMac should produce data that can be checked independently. These items come first so that weeks of compute are not spent on output nobody can audit.

- [x] **A1. Certificate format.** (2026-10-11: format v2 in `docs/certificate-format.md` — contraction rounds and multiplier hints; `kisscrit/core/certificate.py` replays it)
  - **Change:** the branch and bound writes, for each candidate, a compact proof tree: starting boxes per chart, bisections, contractions (removed slabs), and leaves with their exclusion reason (invalid, vertex sum, linear, Gordan or Motzkin, terminal).
  - **Spec:** `docs/certificate-format.md`, written so the checker needs nothing else.
  - **Done when:** the run writes certificates, and their size per candidate has been measured.
- [x] **A2. Independent checker.** (2026-10-11: `kisscrit/core/{interval,model,checker}.py`, `scripts/completeness/check_certificates.py`, `tests/test_certificates.py`; verifies all pilot certificates; rejects tampered ones)
  - **What it is:** a small, readable program built on a verified interval library (Arb, via python-flint).
  - **What it does:**
    - re-evaluates the face-angle model;
    - checks that the leaves of each tree partition each starting box exactly;
    - re-derives every exclusion, recomputing multipliers and Gordan vectors itself.
  - **Rigour:** it does not use the search code's rounding scheme, and floating-point LPs appear in it only as a source of multipliers. A leaf it cannot verify is subdivided and retried, never accepted.
  - **Done when:** the checker verifies the certificates from the pilot (A4) end to end.
- [x] **A3. Lemma inventory.** (2026-10-11: `docs/lemmas.md` L1–L16, `STATUS.md`, `scripts/completeness/lemmas.py`, `tests/test_lemmas.py`. Lemma P and Lemma R proved in checker arithmetic. Three citations or write-ups remain and are listed in `STATUS.md`: L6 convexity, the L8 derivation, and L9 Tammes uniqueness at θ_max. They go to D2/D3.) Every rule used to prune candidates or boxes gets a numbered entry, with a proof, a reference, or a certificate script, and its status. The rules:
  - G ⊆ hull edges;
  - faces are convex 3-, 4- or 5-gons;
  - degrees are 2..5;
  - all corners lie in [α, π];
  - quadrilateral corners are ≤ 2α;
  - degree-2 vertices have straight corners;
  - a near-straight corner forces degree ≤ 3;
  - Lemma P;
  - no rattlers;
  - the facts used in the linear relaxation.

  **Done when:** every rule the run relies on is either proved or is explicitly listed as an assumption in `STATUS.md`.
- [ ] **A4. Pilot with certificates.**
  - **Run:** phase A (candidate generation) of the conditional route, plus one certify shard, here.
  - **Measure:** certificate size, time, and the fraction of candidates that are excluded, left terminal, or out of budget.
  - **Done when:** the checker (A2) verifies the pilot, and the run settings — budget, number of shards, where certificates are stored (GitHub, or Zenodo at the end) — are fixed from the measurements.
- [ ] **A5. Freeze and start the run.**
  - **Steps:** tag the run version, rebuild `kisscrit-imac-run.sh` from the tag, and start the iMac.
  - **Done when:** the iMac is pushing certified shards.

## B. Repository and first dissemination

- [ ] **B1. Restructure.**
  - **Layout:**
    - `core/`: interval model, certificate checker, lemma scripts — the audited part;
    - `search/`: candidate generation, fast branch and bound, runners;
    - `explore/`: census, containers, Edwards work — no rigour claims.
  - **Delete:** superseded prototypes (`bnb.py`, `contract.py`, the old `ad.py`, `faces.py` once nothing uses them).
  - **Done when:** there is no dead code and the tests pass.
- [ ] **B2. `STATUS.md`.** Say what is certified (the twelve configurations, their co-indices, the index-1 joins), what is numerical (the census, the hull survey, the containers, the Edwards pilot), and what is in progress (completeness). Each claim links to the script that supports it.
- [ ] **B3. Reproducibility.**
  - **Contents:** a pinned environment (conda lock file), `make verify` (re-checks every stored certificate and certified configuration), and tests tying each statement in the paper to a script.
- [ ] **B4. v0.1 research snapshot with a DOI.**
  - **Steps:** make the repository public, add `CITATION.cff` and `.zenodo.json`, enable the Zenodo–GitHub integration, and make a GitHub release.
  - **Label:** "research snapshot; completeness not established".
  - **Versioning:** the paper's version will be v1.0 under the same concept DOI.

## C. Mathematics needed for completeness

- [ ] **C1. Hull lemma (a).**
  - **To prove:** a certified bound θ*(T) ≤ 60° for each of the 7591 triangulations other than T0, T1, T5 and T7.
  - **How:** first a degree and area argument to dispose of most of them, then interval branch and bound for the rest.
  - **Known difficulty:** about 8 triangulations reach exactly 60°.
- [ ] **C2. Hull lemma (b).** In T1, T5 and T7 with minimum angle above 60°, the flipped diagonal is never a contact.
- [ ] **C3. Budget cases.** Stronger exclusion for low-edge graphs (E = 17–19, where V_G has dimension 3–5): first-order (Taylor-model) Gordan certificates, so that exclusion along V_G works at coarse widths.
- [ ] **C4. Terminal boxes.** Identify every terminal box with a census configuration, using X-space Krawczyk uniqueness. Here a box is never declared to be "a new critical point" without certification.
- [ ] **C5. The range 1 < r < r_min.**
  - **Why it is open:** the computation starts at r_min = 1.0001.
  - **Known:** the lowest critical points lie about 25° away from FCC and HCP, so this is not a local problem around FCC/HCP.
  - **Needed:** an argument covering this range.
- [ ] **C6. (Optional for the paper.) Higher attaching maps.** The cells at the co-index-2 and co-index-3 saddles, and b1, b2 and b3 just above r = 1.

## D. The paper (Experimental Mathematics)

- [ ] **D1. Mathematical fixes in Sections 3–5.**
  - SO(3) acts freely on labelled configurations (N ≥ 3), so the deformation lemma is applied on the manifold BConf(N) instead of claiming equivariance.
  - "The largest r with BConf connected" becomes a supremum, with a proof that it is a critical radius.
  - Precise theorem numbers from Bochnak–Coste–Roy for the components, semialgebraic arcs, piecewise C¹ arcs, and algebraic points.
- [ ] **D2. References.** Check the CGKM journal version, the Holmes-Cerfon–Theran–Gortler volume and pages, the Musin–Tarasov volume and pages, and that `tfs` is the style the journal requires.
- [ ] **D3. Certification appendix.** Written from the checker (A2): the model, the exclusion tests, the certificate format, and how to re-verify.
- [ ] **D4. Completeness section.** Rewrite it from the results of A and C: the theorem if everything closes, otherwise a precise statement of what is proved and what remains.
- [ ] **D5. Final pass.** Read-through, the paragraph on authorship and AI use, the cover letter, and submission. No submission before D4 is settled.

## E. Future projects (not blocking; see `open-questions.md`)

- Edwards entropy and basin weights: free clusters, containers, channels.
- Sticky growth in flat geometry.
- Exactness of the remaining critical values.
