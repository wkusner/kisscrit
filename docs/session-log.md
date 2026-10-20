# Session log

One dated entry per checkpoint of a controlled session (commit hash, what changed, tests and results, deviations, seeds). Times are US Central.

## 2026-10-17, F3 code session (Opus 5.5, cloud clone), prompt v1.3

### 20:02 CT, start
- Start state: HEAD 1fe331f ("gen cond", iMac) on top of 0ff26f3. `git merge-base --is-ancestor 0ff26f3 HEAD` succeeded; no change since 0ff26f3 to the read-first files.
- Baseline: `python3 -m pytest -q` 55 passed; `make -C csrc test` 0 containment failures.
- Versions: gcc 13.3.0, FLINT 3.0.1 (Ubuntu `libflint-dev`), Python 3.13.16, python-flint 0.9.0, numpy 2.5.3.

### 20:16 CT, scope change by Wöden (verbatim)
> SCOPE CHANGE: split checklist item F3. F3a = the C model (csrc/kc_model.{c,h}, with the bit-pattern diff test against model.py). F3b = the bit-exact replay diff against certificate.replay (spec v0.3 section 4.2), which will be a separate session. This session does F3a only. In docs/checklist.md replace the F3 line with two lines, F3a (tick it when done) and F3b (unticked). "DONE" in the prompt now means F3a is done and the F3a line is ticked. Do not touch HCP lambda or Appendix B; the control chat handles those. Record this message verbatim in docs/session-log.md. Everything else in the prompt stands.

### Milestone: model compiles, diff test passes (commit: see the commit that adds this file)
- Changed: `csrc/kc_model.c`, `csrc/kc_model.h` (new), `csrc/Makefile` (the model is built into `libkisscheck.so`), `tests/test_cmodel.py` (new), `docs/checklist.md` (F3 split into F3a and F3b, both unticked), `docs/code-tour.md` (section 2 and the trust table), this file.
- Diff test (computed, `tests/test_cmodel.py`, seed 20261017, tolerance 0, comparison by IEEE bit pattern): default size 1,166 boxes on 53 charts, 806,432 doubles, 0 mismatches; `KISSCRIT_CMODEL_BOXES=20000`: 19,348 boxes on 691 charts (ok=1 2,953, ok=0 16,395; tightening narrowed a corner on 5,653), 13,348,384 doubles, 0 mismatches, max abs diff 0.
- Test power (computed): with the exact-zero rule of `sum0` disabled, the test reports 9,008 mismatches on 200 boxes; restored, 0.
- Full suite: 59 passed, 1 failed: `tests/test_code_standards.py::test_c_files_have_status_headers`, because its list `TRUSTED_C` names only `kc_arith.{c,h}` and the new trusted files are not on it. That file is not on this session's allowlist: stop rule. F3a is not ticked.
- Deviations: none. Seeds: none.

### 20:28 CT, scope change by Wöden (verbatim)
> SCOPE CHANGE: you may edit tests/test_code_standards.py, only to add csrc/kc_model.c and csrc/kc_model.h to TRUSTED_C. Change nothing else in that file and do not alter the check logic. Make it a separate commit with the "Control: F3" trailer, run the full suite, record the result in docs/session-log.md, then deliver one patch set for 1fe331f..HEAD (both commits) and send a new report block. F3a is done only when the suite is fully green and the F3a line is ticked.

### 20:32 CT, F3a done
- Commits: f6385e9 (the model, test and docs), f371d08 (`TRUSTED_C` gains `csrc/kc_model.c`, `csrc/kc_model.h`; nothing else in that file changed), and the commit that adds this entry (F3a ticked).
- Full suite after a clean rebuild (computed, HEAD f371d08): `python3 -m pytest -q` 60 passed (baseline 55 + 5 new); `make -C csrc test` 0 containment failures.
- F3a ticked in `docs/checklist.md`; F3b unticked (a separate session).
- Deviation: three commits rather than two. The tick and this entry can only be committed after the green run.
- Seeds: none.

## 2026-10-20, F3b code session (Opus 5.5, cloud clone), prompt v1.1

### 09:37 CT, start
- The repo was not attached when the session began; it was attached (read access) and cloned fresh from origin.
- Start state: START_HEAD c318385 ("gen cond", iMac; the three commits after 2ef4b4c touch only `results/`). `git merge-base --is-ancestor 2ef4b4c HEAD` succeeded; no change since 2ef4b4c to the read-first files.
- Baseline: `python3 -m pytest -q` 60 passed; `make -C csrc test` 0 containment failures.
- Versions: gcc 13.3.0, FLINT 3.0.1 (Ubuntu `libflint-dev`), Python 3.13.16, numpy 2.5.3, numba 0.68.0, python-flint 0.9.0, sympy 1.14.0; flags `-O2 -std=c11 -Wall -Wextra -fno-fast-math -ffp-contract=off -fPIC`.
- Measured on all 45,260 `cond` certificates (70,235 run charts): box dimension n from 1 to 14; largest chart 1,387,325 events (cond 39-of-40), q stream up to 5,198,740 entries, mu up to 2,501,465; proof-tree stack at most 85 entries; every run chart has a mu stream.
- Finding: `certificate.replay`'s docstring says it raises ValueError, but a q read past the end of q, or a coordinate k >= n, raises IndexError. `scripts/completeness/check_certificates.py` treats any exception as a failed certificate, so a malformed chart is still rejected. Not changed (off the allowlist); the C replay reports both, as `KC_RERR_INDEX` for the IndexError sites.

### 09:56 CT, milestone: replay compiles and matches
- Commits: a975fa4 (`csrc/kc_replay.{c,h}`; `csrc/Makefile` builds it into `libkisscheck.so`; `TRUSTED_C` gains the two files, nothing else in that file changed), 4c3fc6c (`tests/test_creplay.py`).
- `objdump -d build/libkisscheck.so` shows no FMA instruction.

### 10:09 CT, F3b done
- Replay test (computed, `tests/test_creplay.py`, tolerance 0, IEEE bit patterns): 1,000 certificates drawn with seed 20261020, 1,560 charts, 550,688 nodes, 0 mismatches, all complete.
- Extra runs, not in the suite (computed with the same functions): 5,000 certificates, seed 7: 7,746 charts, 6,360,896 nodes, 5,051,007 slabs, 0 mismatches; the heaviest chart (1,387,325 events, 503,007 nodes): 0 mismatches.
- Malformed streams: 47 hand-made charts reach each of the 5 `raise ValueError` lines and 8 IndexError reads, plus edge cases (no mu stream, empty, short and repeated hints, a coordinate contracted twice in a round, point boxes, -0.0, infinities, NaN, subnormals, n = 0); truncated and extended events, q and mu on 12 real charts (72 cases); 720 random edits of real charts (every end occurs). C and Python agree on every node and on where and how the replay ends.
- Test power (computed, mutant copies of `kc_replay.c` compiled with the same flags; hand-made cases plus 60 certificates): `fma()` grid value 985 mismatching nodes; grid as (1 - t) lo + t hi 1,047; midpoint one ulp low 1,159; upper slab before lower 124; midpoint lo + 0.5 (hi - lo) 2, only on the wide-box case (it equals 0.5 (lo + hi) whenever hi - lo is exact).
- Full suite after a clean rebuild: `python3 -m pytest -q` 116 passed (baseline 60 + 56 new); `make -C csrc test` 0 containment failures.
- F3b ticked in `docs/checklist.md`; `docs/code-tour.md` sections 0, 2 and 3 updated.
- Correction: the milestone report block sent in chat was stamped 10:24 CT; the clock shows the milestone at 09:56 CT and this entry at 10:09 CT.
- Deviations: none. Seeds: none.

## 2026-10-20, C5a code session (Opus 5.5, cloud clone), prompt v1.0

### 12:48 CT, start
- The repo was not attached when the session began; it was attached (read access) and cloned from origin (depth 50). Paused by Wöden 12:39 CT before the clone; resumed 12:41 CT on his message that the control chat approved.
- START_HEAD 77711f5. `git merge-base --is-ancestor 77711f5 HEAD` succeeded; no change since 77711f5 to scripts/r1, tests/test_r1.py, docs/checklist.md, docs/code-tour.md.
- Baseline: `python3 -m pytest -q` 116 passed (after installing `libflint-dev`; without it 42 passed, 74 skipped); `make -C csrc test` all containment tests passed.
- Commit 3f917fb: not found (`git cat-file -t`, `git fetch origin 3f917fb`, `git log --all --oneline | grep`).
- Versions: Python 3.13.16, sympy 1.14.0, mpmath 1.3.0, numpy 2.5.3, python-flint 0.9.0, gcc 13.3.0, FLINT 3.0.1.
- Noted: the prompt header says v1.0 of 11:35 CT; the message that carried it was sent at 11:17 CT.

### 12:58 CT, milestone: parts A-D computed, routes agree (commit 8faad96)
- New `scripts/r1/exact_lambda_stress_space.py` (untrusted: `tests/test_code_standards.py` rejects "Status: trusted" on files outside its list, which this session may not edit). Code tour section 7 entry in the same commit.
- Computed (script at 8faad96, exact in QQ<sqrt2+sqrt3>): FCC and HCP each have 24 contacts, rank R = 20, a 4-dimensional self-stress space and a 1-dimensional flex space in V. Route 1 (stress Laplacian) and route 2 (closed-form second derivative) agree on all 24 edges. FCC: all values are 1/6. HCP: 3/22 on the 6 equator edges, 1/12 on the 12 equator-triangle edges, 8/33 on the 6 triangle edges. At HCP the witness y with R y = q - (3/22)1 exists and is verified exactly, so f^T Q_w f = (3/22) sum w for every self-stress w. Each basis stress has ratio 3/22; HCP has 5 extreme rays, each with ratio 3/22.
- Finding (computed, three checks): at HCP the uniform weight is NOT a self-stress. Exact `is_stress` gives False. By hand, c1 + c5 + c6 + c9 = (2, sqrt3/3, 0) is not parallel to c0 = (1, 0, 0). A float cross product gives the same. A strictly positive self-stress does exist: the sum of the extreme rays, verified exactly. This contradicts `scripts/r1/exact_lambda.py` lines 16, 132 and 135 (which call the uniform weight at HCP a stress) and the Proposition B.5 sentence "K(C) is 4-dimensional and contains the uniform stress" for HCP. The uniform weight still gives 3/22, because q - (3/22)1 sums to 0. Off-limits file: stop rule; not fixed.
- Tests: `tests/test_code_standards.py` 6 passed. Full suite not rerun at this checkpoint. `tests/test_r1_stress_space.py` not yet written.
- Deviations: the new script is marked untrusted (see above). Seeds: none.
