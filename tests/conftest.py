"""pytest configuration for the suite (CLEAN-9; tooling only, no part of the proof).

KISSCRIT_REQUIRE_ALL=1 turns every skipped test into a failure.  Several tests skip silently when something is
missing (the C kernel not built, libflint absent, the `cond` certificates or gen files absent, no C compiler);
on the run server a skip would hide the loss of a check, so the server sets this variable.  Two skips stay allowed: the opt-in test that regenerates the DOM-1 forgery fixture (`DOM1_REGEN=1`; its reason
names DOM1_REGEN) and the shellcheck lint of run.sh (reason 'shellcheck not installed').  Without the variable nothing changes.
"""
import os

import pytest

# Skips that stay allowed under KISSCRIT_REQUIRE_ALL=1: the opt-in forgery-fixture regeneration, and the shellcheck lint
# of run.sh (a lint of run tooling, not a check of the mathematics; `bash -n` runs before it and is not skipped; the
# server's conda environment has no shellcheck).
ALLOWED_SKIP_MARKERS = ('DOM1_REGEN', 'shellcheck not installed')


@pytest.hookimpl(hookwrapper=True)
def pytest_runtest_makereport(item, call):
    outcome = yield
    report = outcome.get_result()
    if os.environ.get('KISSCRIT_REQUIRE_ALL') != '1' or not report.skipped or hasattr(report, 'wasxfail'):
        return
    lr = report.longrepr
    reason = lr[2] if isinstance(lr, tuple) and len(lr) == 3 else str(lr)
    if any(m in reason for m in ALLOWED_SKIP_MARKERS):
        return
    report.outcome = 'failed'
    report.longrepr = 'KISSCRIT_REQUIRE_ALL=1: a skipped test counts as a failure. Skip reason: %s' % reason
