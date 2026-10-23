"""pytest configuration for the suite (CLEAN-9; tooling only, no part of the proof).

KISSCRIT_REQUIRE_ALL=1 turns every skipped test into a failure.  Several tests skip silently when something is
missing (the C kernel not built, libflint absent, the `cond` certificates or gen files absent, no C compiler);
on the run server a skip would hide the loss of a check, so the server sets this variable.  The one skip that
stays allowed is the opt-in test that regenerates the DOM-1 forgery fixture (`DOM1_REGEN=1`): its reason names
DOM1_REGEN.  Without the variable nothing changes.
"""
import os

import pytest

ALLOWED_SKIP_MARKER = 'DOM1_REGEN'


@pytest.hookimpl(hookwrapper=True)
def pytest_runtest_makereport(item, call):
    outcome = yield
    report = outcome.get_result()
    if os.environ.get('KISSCRIT_REQUIRE_ALL') != '1' or not report.skipped or hasattr(report, 'wasxfail'):
        return
    lr = report.longrepr
    reason = lr[2] if isinstance(lr, tuple) and len(lr) == 3 else str(lr)
    if ALLOWED_SKIP_MARKER in reason:
        return
    report.outcome = 'failed'
    report.longrepr = 'KISSCRIT_REQUIRE_ALL=1: a skipped test counts as a failure. Skip reason: %s' % reason
