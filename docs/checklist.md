# Project checklist: completeness of the N = 12 census and the Experimental Mathematics paper

Purpose: the working checklist for this project, in the order we intend to do the items. Each item says what "done" means. After an item is closed, Claude gives a short progress report and review, and asks before starting the next one.

Version 2, 2026-10-16 (v2: section F for the `full` route and the compiled pipeline; C1/C2 dropped; C5 now Appendix B). Update this file as items close; don't keep a separate list.

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
- [x] **A3. Lemma inventory.** (2026-10-11: `docs/lemmas.md` L1–L16, `STATUS.md`, `scripts/completeness/lemmas.py`, `tests/test_lemmas.py`. Lemma P and Lemma R proved in checker arithmetic. 2026-10-12: L6 proved, L8 derived, uniqueness at θ_max proved via Lemma U; the LaTeX write-up goes to D3.) Every rule used to prune candidates or boxes gets a numbered entry, with a proof, a reference, or a certificate script, and its status. The rules:
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
  - **Run:** the pilot of the conditional route on the run server (2026-10-22, tree 322f8bc: gen, 40 certify shards, check; see `STATUS.md`). (2026-10-23, S-04: the text said "plus one certify shard, here"; the server replaced the iMac.)
  - **Measure:** certificate size, time, and the fraction of candidates that are excluded, left terminal, or out of budget.
  - **Done when:** the checker (A2) verifies the pilot, and the run settings — budget, number of shards, where certificates are stored (GitHub, or Zenodo at the end) — are fixed from the measurements.
- [ ] **A5. Freeze and start the run.**
  - **Steps:** tag the run version and start `scripts/full/run.sh` from the tag on the run server (checklist F8). The iMac runner `scripts/completeness/mac_run.sh` and its instructions are a record of the old run. (2026-10-23, S-04: the text said "rebuild `kisscrit-imac-run.sh`", a file that does not exist.)
  - **Done when:** the run is under way on the server and pushes its progress and check results (F8 M2).

## B. Repository and first dissemination

- [ ] **B1. Restructure.**
  - **Layout:**
    - `core/`: interval model, certificate checker, lemma scripts — the audited part;
    - `search/`: candidate generation, fast branch and bound, runners;
    - `explore/`: census, containers, Edwards work — no rigour claims.
  - **Retire:** superseded prototypes (`bnb.py`, `contract.py`, the old `ad.py`, `faces.py` once nothing uses them) move to `attic/`, with their history kept; nothing is deleted (2026-10-23, S-04: the text said "Delete").  *Status 2026-10-22 (CLEAN-A, CLEAN-B): `bnb.py`, `contract.py`, `ad.py`, `lpfilter.py` and `rattler_check.py` moved to `attic/` (not deleted; nothing imports them); `faces.py` stays: `structure.py` imports it.*
  - **Done when:** there is no dead code and the tests pass.
- [ ] **B2. `STATUS.md`.** Say what is verified in interval arithmetic or exactly (the twelve configurations, their co-indices, the index-1 joins), what is numerical (the census, the hull survey, the containers, the Edwards pilot), and what is in progress (completeness). Each claim links to the script that supports it.
- [ ] **B3. Reproducibility.**
  - **Contents:** a pinned environment (conda lock file), `make verify` (re-checks every stored certificate and certified configuration), and tests tying each statement in the paper to a script.
- [ ] **B4. v0.1 research snapshot with a DOI.**
  - **Steps:** make the repository public, add `CITATION.cff` and `.zenodo.json`, enable the Zenodo–GitHub integration, and make a GitHub release.
  - **Label:** "research snapshot; completeness not established".
  - **Versioning:** the paper's version will be v1.0 under the same concept DOI.

## C. Mathematics needed for completeness

- [-] **C1. Hull lemma (a).** (2026-10-16: dropped from the proof path; the `full` route removes H. Kept in the paper as Conjecture 6.3 with the survey, marked unused.)
  - **To prove:** a rigorous bound θ*(T) ≤ 60° for each of the 7591 triangulations other than T0, T1, T5 and T7.
  - **How:** first a degree and area argument to dispose of most of them, then interval branch and bound for the rest.
  - **Known difficulty:** about 8 triangulations reach exactly 60°.
- [-] **C2. Hull lemma (b).** (2026-10-16: dropped with C1.) In T1, T5 and T7 with minimum angle above 60°, the flipped diagonal is never a contact.
- [ ] **C3. Budget cases.** (Rewritten 2026-10-23 per D-031. The earlier text asked for first-order (Taylor-model) Gordan certificates; its premise was found wrong by the DES-F2 design document `admin-twelve-fable-design-C3-2026-10-21.md` (maintainer's project folder `06_Admin`, available on request): for 53 of the 54 `cond` budget candidates the budget is spent on coarse `vsum` and `invalid` leaves; reproduced by AUD-2, scratch.) In order:
  - **R1.** Measure the budget first: rerun the 53 non-census budget candidates at 3,000,000 boxes (30,000,000 for any still open), about 38 CPU-hours; needs a `--codes` option in `run_shard.py` (session C3-1).
  - **R2.** Identification-based termination of near-critical clusters (C4-5).
  - **R3.** A bisection rule aimed at the vertex-sum width; untrusted, no proof obligation (session C3-2).
  - **R4.** Fallback, only if R1–R3 leave candidates open: the realization chart with the Krawczyk contractor of DES-F2 lemma L28 (session C3-3).
  - Wording for C3 and H5: session C3-4. If candidates stay unresolved, the paper's statement falls back to H5, an explicit list of unresolved codes (D-032).
  - **Done when:** every candidate has zero open leaves, or the unresolved codes are listed for H5.
- [ ] **C4. Terminal boxes.** Identify every terminal box with a census configuration, per leaf (C4a, C4b below). A box is never declared to be "a new critical point" without a certificate accepted by the checker; such a point is recorded and the run continues (D-027). (2026-10-22, D-021: split into C4a and C4b; 2026-10-22, D-025: the per-leaf trichotomy replaces the X-space Krawczyk uniqueness balls; neither is done.)
- [ ] **C4a. Regular terminal boxes.** (Rewritten 2026-10-23 per D-025; the D-021 text named X-space Krawczyk uniqueness balls.) Decide every unexcluded leaf (terminal or open) at a point where the square KKT system is nonsingular by the per-leaf trichotomy of DES-F2 lemma L23: the leaf is excluded by Krawczyk emptiness, or contains exactly one critical configuration, which is compared with the census point, or is undecided. L23 is proposed in the DES-F2 design document `admin-twelve-fable-design-C4-2026-10-21.md` (maintainer's project folder `06_Admin`, available on request) and was reviewed by AUD-2 (proved from L18–L22, modulo the orientation convention C1 of the developing map, L20); it is not yet in `docs/lemmas.md`. Sessions C4-1, C4-2, C4-4.
- [ ] **C4b. Singular terminal boxes.** An isolation argument for the terminal boxes at DOD and M, where the square KKT system is singular: first-order isolation with an explicit radius δ_C (DES-F2 lemmas L24, L25; reviewed by AUD-2: proved; the constants computed in floating point only, not yet verified in interval arithmetic); M stored as 80-digit intervals (D-026); the √N loss of the Procrustes step accepted (D-028). Session C4-3. 458 singular boxes in the `cond` run, of which 456 are DOD and 2 are M (reported by the DES-F2 documents, the DES-F2 design document `admin-twelve-fable-design-C4-2026-10-21.md` (maintainer's project folder `06_Admin`, available on request) and the DES-F2 session prompt, from BALL-KILL's floating developing map; not reproduced by the control chat, by the CLEAN-A session or by AUD-2, which did not open the ball file), plus the lattice rule D5 for the DOD sub-support codes (DES-F2 idea D5, proposed, not implemented).
- [~] **C5. The range 1 < r < r_min.** (2026-10-16: the argument is Appendix B of the paper — compactness, second-order isolation with explicit radius, exact lambda at FCC (1/6) and HCP (>= 1/12) in `scripts/r1/exact_lambda.py`. Open: the Bott form at the KKLS M_5 4-simplex, the list of strata of Z_1, the band run with the delta-ball exclusion; see F7.) C5a (2026-10-20, computed with `scripts/r1/exact_lambda_stress_space.py` at f339a3b, tests in `tests/test_r1_stress_space.py`): at HCP f^T Q_w f = (3/22) sum w for every self-stress w. The witness is y with R y = q - (3/22)(1,...,1), exact in Q(√2, √3), and two routes give the same 24 per-edge values. At FCC the value is 1/6. The uniform weight is a self-stress at FCC but not at HCP, where a strictly positive self-stress exists.
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

## F. The `full` route with the compiled pipeline (every error bound from Arb) (2026-10-16)

Decisions of 2026-10-15 (`docs/decisions.md`): run `full` (all 7,595 triangulations, no hypothesis H); rebuild the checker in C on FLINT/Arb with every error bound from a rigorous source, Arb (the decision of 2026-10-15 says "a certified source"; pipeline specification v0.3, `code-completeness-compiled-pipeline-spec-v0.3-2026-10-15.md`, an internal document, available on request); report progress to GitHub during the run; archive certificates on Zenodo, not in git. `cond` becomes a pilot and a regression set.

- [x] **F0. Checklist, status, decisions.** (2026-10-16: this section; `STATUS.md` v3; `docs/decisions.md`; the 90,836/45,260 count resolved — 90,836 is the number of gen lines including duplicates, 45,260 the number of distinct codes.)
- [x] **F1. Complete reference checker (Python).** (2026-10-16: V1–V5 in `kisscrit/core/checker.py` and `scripts/completeness/check_certificates.py` (manifest, summary counters); V6 `tests/test_checker_imports.py`; V8 in `mac_run.sh`; lazy tightening and the combination-form multiplier LP; mutation tests 1–8 in `tests/test_certificates.py`. Measured here: the 3,654-box certificate 5.0 s → 3.0 s with identical verdict (the spec's 42.7 → 20.6 s was a slower machine); 25 random `cond` certificates give identical verdicts and counts under the old and new checker. One-ulp search fixes for V1/V3 (`f0a32cb`); the `cond` certificates fail V1/V3 by those ulps and are a pilot. Readability review by a fresh agent done; its 21 findings fixed in the same milestone: pentagon regions written out as L10′/L11′ with the cover proof, the straight-corner box [π⁻, π⁺], the Lipschitz exception at the collapsing diagonal, terminal-box width check, manifest refuses completeness without gen, references to untrusted code labelled.) V1–V5 of spec v0.2, lazy tightening, and the mutation tests.
- [x] **F2. Elementary functions in C on Arb (rigorous enclosures).** (2026-10-16: `csrc/kc_arith.{c,h}`, `csrc/Makefile` (`-O2 -fno-fast-math -ffp-contract=off`), `csrc/test_arith.c`; `kc_init` asserts round-to-nearest and no FMA contraction; containment against Arb at 256 bits on 10⁶ random intervals per function plus special endpoints: 0 failures, measured 0.9–1.8 µs per call unoptimised; `tests/test_carith.py`: bit-for-bit identical to `kisscrit/core/interval.py` on random and special inputs (the one documented divergence: cos/sin of |x| > 10⁹ return [−1, 1] in C). FLINT 3.0.1 headers here came from Ubuntu's `libflint-dev`; the server path (conda-forge `libflint` or a source build) is settled by the preflight.) FLINT/Arb on interval endpoints; IEEE basic operations stepped outward; rounding-environment asserts; the FMA test; containment against Arb at 256 bits on ≥ 10⁶ random intervals per function; build with `-O2 -fno-fast-math -ffp-contract=off`.
- [x] **F3a. C model.** (Split from F3 by Wöden, 2026-10-17 20:16 CT. Done 2026-10-17: `tests/test_cmodel.py` bit-pattern diff, seed 20261017, 19,348 boxes on 691 `cond` charts, 13,348,384 doubles, 0 mismatches; full suite 60 passed.) `csrc/kc_model.{c,h}`: corners, vertex sums, Jacobians, tightening, mirroring `kisscrit/core/model.py`; bit-pattern diff against the Python model on random boxes (`tests/test_cmodel.py`).
- [x] **F3b. C replay.** (Split from F3 by Wöden, 2026-10-17 20:16 CT; a separate session. Done 2026-10-20: `tests/test_creplay.py` bit-pattern comparison, seed 20261020, 1,000 `cond` certificates, 1,560 charts, 550,688 nodes, 0 mismatches; every raise site of `certificate.replay` covered; five compiled mutants detected; full suite 116 passed.) `csrc/kc_replay.{c,h}`: bit-exact replay diff against `certificate.replay` (spec v0.3 section 4.2).
- [x] **F4. C kernel and fallback.** (2026-10-20, F4 session; ticked 2026-10-20 22:40 CT after Wöden accepted the performance gate (below). Done: `csrc/kc_leaf.{c,h}` (hint test, invalid/vsum untightened and tightened, bisection; LP tests to the fallback), `kisscrit/ckernel.py` (fallback to the unchanged reference), `tests/test_ckernel.py` (7 compiled mutants, all detected), `scripts/completeness/diff_kernel.py`, `bench_kernel.py`; `kc_hull` brought in line with `interval.hull`; F9 readability review done and fixed. Differential, computed at b06b2f5: 500 `cond` certificates (seed 20261021), 228,236 units, C closes 228,215, C unclosed 10 (all Gordan leaves, the reference accepts all), C accepts where Python rejects 0, 603,714 per-test comparisons 0 mismatches; a fresh sample of 242 certificates (shard 174 of 200, seed 20261022, regenerated at 77711f5), 127,165 units, C unclosed 1 (Gordan, accepted), C accepts where Python rejects 0. Performance, computed at b3f2f8f on a 2-core Xeon 2.1 GHz VM, gcc 13.3.0 `-O2 -fno-fast-math -ffp-contract=off`, 7 repeats: kernel 317 µs per unit on the fresh sample (min 292, max 334), 285 µs on the 500-sample (min 270, max 318), i.e. 0.33–0.36 of the 0.87 ms per box quoted below; fallback 0.06 s and 0.23 s in all. Spec v0.3 §4.4 target, quoted from Project knowledge: "checker CPU time ≤ search CPU time on the same candidates". On the fresh sample, same machine: checker 40.4 s against search 169.6 s (ratio 0.24), or 48.5 s (ratio about 0.83) without the two candidates whose search time is numba's compile time. Gate: accepted by Wöden, 2026-10-20 22:38 CT, as the spec v0.3 §4.4 target above (checker CPU time ≤ search CPU time on the same candidates), measured ratio 0.24 on the fresh sample and about 0.83 without the two numba-compile candidates, on one 2-core VM; a stricter threshold, if needed, belongs to the F8 run plan.) Leaf tests in the v0.3 §4.2 order, `ctypes` interface, Python fallback for unclosed boxes; differential run on the `cond` certificates (stop if C accepts where Python rejects); the performance gate (kernel ns per box against the search's 0.87 ms per box).
- [ ] **F5. Format v3 witnesses in the search.** Decoder still reads v2; 1,000 `cond` candidates re-run with v3, checked with the C kernel; certificate bytes per candidate measured and recorded here.
- [ ] **F6. Gen for `full`.** G2 exact vertex-mask table; compiled enumerator written fresh (so `candidates.py` stays an independent second enumerator); G1 exact rational Farkas certificates for every LP discard; G3 set comparison on T0, T1, T5, T7 and 2% of the rest, counts on all 7,595; the L12 cross-check. **Done when:** the measured `full` candidate count is reported.
- [ ] **F7. Algorithmic items** (the G2–G10 handoff; rewritten 2026-10-23 per D-025 and D-031, the earlier text named "C3 Taylor-model Gordan for budget cases"):
  - C4 terminal identification, sessions C4-1 to C4-6: C4-1 the trusted square system and Krawczyk (L18, L19); C4-2 developing map, gauge and stress enclosure (L20–L22); C4-3 isolation at DOD and M (L24, L25); C4-4 the per-leaf trichotomy and checker integration (L23, L26, L27); C4-5 search-side termination by identification (untrusted); C4-6 paper and `STATUS.md` wording. Order: C4-1, C4-2, C4-4; C4-3 independent; C4-5 after C4-4; C4-6 last.
  - C3 budget cases, sessions C3-1 to C3-4 (R1 budget rerun, R3 bisection rule, R4 killing test and contractor, wording; see C3).
  - The r1 package (constants, δ per stratum, band-run hook, M_5 stratum, margin diagnostic).
  - Attribution of the budget-box blow-up: done (the DES-F2 design document `admin-twelve-fable-design-C3-2026-10-21.md` (maintainer's project folder `06_Admin`, available on request) §0.1, computed; reproduced by AUD-2, scratch).
  - Priority: whatever cuts `budget` boxes (on `cond`, 54 of 45,260 candidates took 41% of the boxes).
- [ ] **F8. Run scripts (Linux first).** `scripts/full/run.sh`, Linux the primary target (a dedicated Ubuntu 26.04 server, no root, no scheduler), macOS the fallback; it replaces `mac_run.sh` for the new run. Setup without root: Miniforge, a conda env with `libflint`, `c-compiler`, git and gh; `gh auth login` (device code over SSH); a preflight that reuses `scripts/full/kisscrit-preflight.sh` (cores, cgroup limits, memory, disk, noexec, tmux and cron, lingering, HTTPS to GitHub, Miniforge, conda-forge, Zenodo and the sandbox), stopping with a clear message when a host is unreachable. Zenodo tokens from `~/.config/kisscrit/zenodo.token` and `zenodo-sandbox.token` (mode 600; looser permissions refused), on macOS from the Keychain; never from the repo, an environment file or the logs. Self-test: unit tests, C kernel smoke test, sandbox Zenodo round trip when a sandbox token exists, peak RSS per worker measured and the worker count capped to fit RAM (default `nproc`). Phases gen → certify (shards) → check (C kernel + fallback) → manifest, all resumable; scratch files under `~/kisscrit-run`, never `/tmp` (tmpfs on 26.04); hashes with Python `hashlib`. Survival: `systemd --user` service when lingering is on, else a detached `tmux` session plus a `@reboot` crontab line, launchd on macOS; `status` reports which. Progress to GitHub: a commit per shard and `results/full/progress.json` at least hourly (rates, no ETA); check results pushed incrementally. Zenodo: one draft per run (id in `results/full/zenodo-state.json`), each certify shard gzipped and PUT to the bucket, `results/full/zenodo-manifest.json` with bytes and SHA-256 pushed; never publish; keep the record at ≤ 90 files; **[ask]** above a projected 45 GB. `status` and `stop` subcommands. The setup guide `docs/setup-guide.md` changes in the same commit as any manual step. M1 done (2026-10-22): scripts/full/run.sh v1 for Linux with setup, preflight, selftest, start, status, stop; M2 GitHub progress, M3 Zenodo, M4 macOS not done.
- [ ] **F9. Code standards (referee-readable).** From F1 on: file headers (purpose, trusted/untrusted, lemmas and appendix sections), function comments with the mathematical claim, a source for every inequality and rounding direction, plain code with paper-notation names, `docs/code-tour.md` kept in the same commit as any trusted-code change, `tests/test_code_standards.py` enforcing the headers. A fresh agent reviews readability at the end of F1, F4 and F6; its findings are fixed before the item is ticked.

## G. Done since the F section (catch-up, 2026-10-22, D-021)

- [x] **EXACT-1. Exact data of the critical configurations found so far.** (2026-10-21/22, EXACT-1 session, origin 5da4bba and before; reported in `docs/session-log.md`; `tests/test_exact_configs.py` passes in the CLEAN-A baseline, computed.) Exact arithmetic in K = Q(c) for 11 configurations: the census values DOD, M, 1.0010073939, 1.0008452029, B, r1, 1.0003062544 and A (eight of the twelve), and FCC, HCP and the hexagonal antiprism band. Data `data/exact/critical-configs-n12.{json,md}`, description `docs/exact-configs.md`. Not done: exact data for the other four census values: 1.0003062547, 1.0002374095 and 1.0002115891 (no integer relation of degree at most 150 at 8000 digits, high precision only) and 1.0002111636 (a candidate minimal polynomial of c of degree 126, noted in the data as "high precision only", not exactly verified). Not a completeness claim. (The CLEAN-A prompt said "the 12 census configurations"; the session log and `docs/exact-configs.md` say 11 exact and four census values without an exact polynomial, and this line follows them.)
- [x] **DOM-1. Domain check V1b.** (2026-10-22, DOM-1 session, f3e8d0d on origin.) V1b in `structural_check` (`kisscrit/core/checker.py`): every start-box coordinate finite, rhombus coordinates in (0, π], pentagon coordinates in [0, 2π]; lemmas L31–L33 in `docs/lemmas.md`. Reported by the DOM-1 session (computed there, not reproduced here): the DES-F2 forged rhombus start box that the checker accepted at cc375bc is rejected by the reference and by the kernel front end; the 45,260 committed `cond` lines contain 0 violations. L31–L33 reviewed by AUD-2 (2026-10-23): proved; L32 also re-derived by the DOM-1 session. Not done: DOM-2 (stop after a structural failure) and DOM-3 (manifest hashes).
