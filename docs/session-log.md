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

### 13:08 CT, control ruling (relayed by Wöden, verbatim)
> Control ruling: continue. The mislabel in exact_lambda.py is a finding, not a blocker; do not change that file. Record it in the session log with the vertex-0 computation. In tests/test_r1_stress_space.py also assert: the uniform weight is a self-stress at FCC and not at HCP; a strictly positive self-stress exists at HCP (print it exactly). Then go on with the plan: at least 3 mutants, the full suite, the C5 sentence, the final log entry, the patches.

### 13:17 CT, C5a done (outcome CONFIRMED, computed)
- Finding, recorded per the ruling: `scripts/r1/exact_lambda.py` line 16 (module docstring), line 132 (docstring of `hcp_stresses`) and line 135 (the label 'uniform') call the uniform weight at HCP a stress. It is not one. Vertex-0 computation, exact in Q(sqrt2, sqrt3): c_0 = (1, 0, 0); its neighbours are c_1 = (1/2, sqrt3/2, 0), c_5 = (1/2, -sqrt3/2, 0), c_6 = (1/2, sqrt3/6, sqrt6/3) and c_9 = (1/2, sqrt3/6, -sqrt6/3). Their sum is (2, sqrt3/3, 0), and its cross product with c_0 is (0, 0, -sqrt3/3), which is nonzero. So with equal weights the tangential force at c_0 does not vanish. `tests/test_r1.py` still passes, because the uniform weight's value is the average of the per-edge values, (6 x 3/22 + 12 x 1/12 + 6 x 8/33)/24 = 3/22. `exact_lambda.py` also decides contacts with `sp.simplify(...) == 0` and per-edge values with `nsimplify(simplify(...))`, which are heuristics; the new script decides both by field equality in K. Neither file was changed.
- Commits: 8faad96 (script and code tour), 91b9e88 (log, start and milestone), 59e39ef (`tests/test_r1_stress_space.py`, code tour test line), and the commit that adds this entry (the C5 sentence in `docs/checklist.md`, appended).
- Computed (script at 8faad96): at HCP f^T Q_w f = (3/22) sum w for every self-stress w. The witness y with R y = q - (3/22)1 is verified by exact multiplication. The basis stresses have sums -18, 12, 18 and -6, each with ratio 3/22. There are 5 extreme rays of the nonnegative cone, each with ratio 3/22. A strictly positive stress exists: 1/20 on (0,1), (2,3), (4,5); 11/150 on (0,5), (1,2), (3,4); 7/200 on the other 18 edges. FCC: all 24 per-edge values are 1/6, the uniform weight is a stress, and there are 4 extreme rays (the great hexagons). Part D (mpmath 1.3.0 at 50 digits; floating point, not proof): max |finite difference - exact| 3.5e-21 with h = 1e-15; eigenvalue minus (3/22) sum w at most 4.3e-50 for three random combinations.
- Tests: `python3 -m pytest -q tests/test_r1_stress_space.py` 14 passed. Full suite 130 passed (baseline 116 + 14). `make -C csrc test` all containment tests passed.
- Mutants (5 run, 5 detected): (1) HCP point 6, x + 1/7, fails the set-up; (2) per-edge value on edge (0,6) + 1/1000 fails both the route comparison and the stress claim; (3) basis stress 0 replaced by e_0 fails both the set-up and the stress claim; (4) lambda = 1/12 fails the stress claim; (5) route 2 with the linear chart c + tF disagrees with route 1.
- Deviations: the new script and its tests are marked untrusted (see the milestone entry). Seeds: none.

### 2026-10-20, control note on the C5a hashes
- The hashes cited in the C5a entries above (8faad96, 91b9e88, 59e39ef, 87335e7) are those of the session's clone; after `git am` and the push they are f339a3b (script), b824ff9 (session log, start and milestone), 05ecc23 (tests) and 5e749ec (checklist and log). Use the origin hashes. The `docs/checklist.md` C5 sentence now cites f339a3b.
- Restored in `docs/code-tour.md`: the sentence "The rest of Appendix B (constants, δ per stratum, the band run) is checklist F7 and not yet in code.", which the C5a patch had replaced.


## 2026-10-20, F4 code session (Opus 5.5, cloud clone), prompt v1.0

### 11:01 CT, start
- The repo was not attached when the session began; it was attached (read access) and cloned fresh from origin.
- Start state: START_HEAD 77711f5 (= origin/main). `git merge-base --is-ancestor 77711f5 HEAD` succeeded; no change since 77711f5 to the read-first files.
- Baseline: `python3 -m pytest -q` 116 passed; `make -C csrc test` 0 containment failures.
- Versions: gcc 13.3.0, FLINT 3.0.1 (Ubuntu `libflint-dev`, installed this session), Python 3.13.16, python-flint 0.9.0, numpy 2.5.3, numba 0.68.0, sympy 1.14.0; flags `-O2 -std=c11 -Wall -Wextra -fno-fast-math -ffp-contract=off -fPIC`.
- Pipeline spec read from Project knowledge (`code-completeness-compiled-pipeline-spec-v0.3`, header 0.3.1), not pasted.
- Plan: sample (a) seed 20261021; fresh sample (b) seed 20261022 picks shard 174 of `--of 200` (242 candidates; the iMac's recorded search time for them is 59.6 s).

### 11:24 CT, stopped before M1 coding (stop rule: a bug in a file outside the allowlist)
- Finding: `csrc/kc_arith.c` lines 161-165, `kc_hull`, is not bit-for-bit `kisscrit/core/interval.py:hull` (lines 145-147), contrary to the contract in `csrc/kc_arith.h`. Python's `min`/`max` keep the first argument on a tie or a NaN comparison; `kc_hull` keeps the second. Failing inputs (computed, `carith.binary('hull', ...)` against `interval.hull`, HEAD 77711f5): hull([-0.0, 1], [0.0, 1]) gives lo -0.0 in Python and 0.0 in C; hull([-1, 0.0], [-1, -0.0]) gives hi 0.0 and -0.0; hull([nan, 1], [0.5, 1]) gives lo nan and 0.5; hull([0.5, 1], [nan, 2]) gives lo 0.5 and nan.
- Impact: no C code calls `kc_hull` at 77711f5 (grep); `tests/test_carith.py::test_clamp_and_hull` uses a strictly wider second argument, so no tie occurs. F4's linear-hint test is the first caller (the `Y` of `checker.test_linear_mu`, line 263). There its arguments come from `kc_sub`, whose endpoints are never exactly ±0 or NaN (`op_down`/`op_up` map 0 to ∓5e-324 and NaN to ∓inf), so the divergence cannot reach an F4 verdict (inferred from the code, not tested). Mathematically both results enclose the same set for the ±0 cases.
- Not changed. Waiting for Wöden.

### 12:33 CT, scope change by Wöden (verbatim)
> SCOPE CHANGE: F4 may also change csrc/kc_arith.c (kc_hull only), tests/test_carith.py (add the four failing inputs and tie and NaN cases) and docs/code-tour.md section 0.1. Make kc_hull follow Python’s min and max exactly: lo = b.lo < a.lo ? b.lo : a.lo; hi = b.hi > a.hi ? b.hi : a.hi. Do it as its own commit, with a test that fails on the old code. Re-run the F2 tests (make -C csrc test and tests/test_carith.py), then continue with M1.

### 12:40 CT, kc_hull fixed (commit: the one that adds this entry)
- `csrc/kc_arith.c` `kc_hull`: now `lo = b.lo < a.lo ? b.lo : a.lo; hi = b.hi > a.hi ? b.hi : a.hi`, with a comment stating the rule. Nothing else in the file changed.
- `tests/test_carith.py`: new `test_hull_ties_and_nan`, comparing bit patterns (the existing `_same` treats 0.0 and -0.0 as equal) on the four failing inputs and all 2,401 interval pairs with endpoints in {-1, -0.0, 0.0, 1, nan, -inf, inf}. On the old `kc_hull` it failed with 1,180 mismatches (computed, HEAD 323dbd8 plus the test); on the new one it passes.
- F2 tests rerun (computed): `make -C csrc test` 0 containment failures; `python3 -m pytest -q tests/test_carith.py` 14 passed; `objdump -d build/libkisscheck.so` shows no FMA instruction.
- `docs/code-tour.md` section 0.1: the min/max tie rule, read for `kc_hull`, `kc_mul` and the model's `intersect`.

### 13:05 CT, milestone M1: C leaf tests compile and agree on a trial
- Correction: the previous entry is stamped 12:40 CT; the clock read 12:35 CT when it was committed (76dd5d4).
- Commit: the one that adds this entry. `csrc/kc_leaf.{c,h}` (new: `kc_leaf_invalid_vsum`, `kc_leaf_context`, `kc_leaf_hint`, `kc_leaf_exclude`, `kc_check_chart`); `csrc/Makefile` builds it into `libkisscheck.so` (nothing else changed); `TRUSTED_C` gains the two files (nothing else changed); `docs/code-tour.md` sections 0 and 4.
- Build: `make -C csrc` without warnings; `objdump -d build/libkisscheck.so` shows no FMA instruction.
- Order of tests (spec v0.3 §4.2, read from Project knowledge): the hint test is step 1 (the only witness in format v2), where `check_certificate` runs it; then untightened invalid/vsum, tightened (charts with a pentagon), bisection with the reference's limits; LP tests left to the fallback.
- Trial differential (computed with the working-tree `scripts/completeness/diff_kernel.py` and `kisscrit/ckernel.py`, committed later; a check, not a reported result): 30 certificates, seed 1, 12,580 units; C closed 12,579, the reference accepted all of them; 1 unclosed (a Gordan leaf), which the reference accepts; C accepts where the reference rejects: 0; per-test comparisons 36,574 (invalid/vsum untightened 12,580, tightened 12,540, hint interval s as bit patterns 11,454), 0 mismatches.

### 13:10 CT, milestone M2: front end and fallback, kernel tests in the suite
- Correction: the M1 entry above is stamped 13:05 CT; the clock read 12:49 CT when it was committed (1a0ba67). From here on each stamp is read from the clock when the entry is written.
- Commit: the one that adds this entry. `kisscrit/ckernel.py` (new, untrusted: ctypes binding of kc_model and kc_leaf, `check_chart`, and `check_certificate` with the Python fallback through the unchanged `Checker` calls, or the whole chart to `checker.check_certificate` when the kernel stops on an error); `tests/test_ckernel.py` (new); `scripts/completeness/diff_kernel.py` (new; the differential driver, committed here because the test imports its `check_one`).
- `tests/test_ckernel.py` (computed): `test_differential_on_cond_sample` 120 certificates, seed 20261023, at most 3,000 events each; no disagreement. My first default (40 certificates) gave 876 units, below the test's own coverage floor of 1,000; I enlarged the sample rather than lower the floor. `test_power`: 7 compiled mutants of `kc_leaf.c`, all detected: hint lower endpoint one ulp inward; Y without outward rounding; hint rows in decreasing order; `excludes_zero` testing one side; `excludes_zero` with a 1e-9 slack; the tightening flag inverted; the validity flag ignored. The tightening mutant was not detected on the first 15 certificates (tightening acts only on pentagon boxes wider than 0.01), so the test runs every mutant over the whole sample, stopping at the first detection; it is caught at certificate 26. Hand-made cases: hint edge cases (repeated row, order, zero, ±inf, NaN, all rows) by bit pattern; malformed streams and a hint row out of range raise in the front end exactly when the reference raises; terminal width and open leaves; the depth and piece limits.
- Full suite after a clean rebuild (computed, working tree at this commit): `python3 -m pytest -q` 123 passed (baseline 116 + 1 hull test + 6 kernel tests); `make -C csrc test` 0 containment failures; no FMA instruction in `objdump -d build/libkisscheck.so`.
- Scratch worktree for the fresh sample: `git worktree add --detach /home/claude/f4-scratch 77711f5` (outside the committed tree).

### 13:43 CT, readability review (F9) and its fixes
- A fresh agent reviewed `csrc/kc_leaf.{c,h}`, `kc_hull`, `kisscrit/ckernel.py` and the code-tour entries at c355ca5 (reading only). It found two blockers, both wrong claims rather than wrong code: (1) the "kernel closes a subset of what the reference closes" argument fails where the reference raises (it evaluates the model at every box's centre, which can raise `model.intersect`'s ArithmeticError; the kernel does not); (2) `ckernel.py` was called "untrusted, moves arrays" although it aggregates `verified`, runs the structural and starting-box checks, aligns units and decides the fallback; and its unit check was an `assert`. It also gave should-fixes: the sum-first mean-value form needs a scalar argument (Lebourg on mu.F plus Clarke's chain rule), not the row-wise inclusion; m in B; the bisection cover with finite endpoints; the counts are the kernel's attribution; "bit-for-bit" is computed, not proved; spec step numbers point to a Project document; L16 tightening citation; and nits.
- Fixed in 7f0dfb8 (comments and docs only; `objdump -d` of a rebuild identical to the running build) and in the commit that adds this entry (`ckernel.py`: docstring states its role, that it must become trusted or be re-derived before it is wired into the checker, and that `counts` are exempt; the `assert` is now a `raise`). `tests/test_ckernel.py` 6 passed, `tests/test_code_standards.py` 6 passed (computed).
- Not fixed, off the allowlist, open for Wöden: the same sum-first gap is in the wording of `docs/lemmas.md` L16 and the `checker.combination` docstring; a finiteness check on lo0, hi0 (proposed for the front end and the reference together). Declined nit: the two capacity exits in `kc_check_chart` (`goto done` and `break`) behave the same; changing one would change the object code that the runs measured.
- Container restart at about 13:45 CT: the first differential runs (b, then a) were lost (b had written 241 of 242 certificates); both are rerun from the next commit.
