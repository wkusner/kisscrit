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
