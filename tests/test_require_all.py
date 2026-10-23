"""KISSCRIT_REQUIRE_ALL=1 (tests/conftest.py): a skipped test fails, except the two allowed skips (CLEAN-B2).

Status: untrusted (a test of the test tooling, no part of the proof).  Each case runs a throw-away pytest in a
temporary directory with a copy of tests/conftest.py, so the real suite is not affected.
"""
import os, pathlib, shutil, subprocess, sys, textwrap

CONFTEST = pathlib.Path(__file__).resolve().parent / 'conftest.py'


def _run(tmp, reason, require):
    shutil.copy(CONFTEST, tmp / 'conftest.py')
    (tmp / 'test_x.py').write_text(textwrap.dedent('''
        import pytest
        def test_x():
            pytest.skip(%r)
    ''' % reason))
    env = {k: v for k, v in os.environ.items() if k != 'KISSCRIT_REQUIRE_ALL'}
    if require: env['KISSCRIT_REQUIRE_ALL'] = '1'
    return subprocess.run([sys.executable, '-m', 'pytest', '-q', '-p', 'no:cacheprovider', str(tmp / 'test_x.py')],
                          cwd=tmp, env=env, capture_output=True, text=True).returncode


def test_skip_fails_under_require_all(tmp_path):
    assert _run(tmp_path, 'C kernel not built', True) != 0


def test_skip_is_a_skip_without_the_variable(tmp_path):
    assert _run(tmp_path, 'C kernel not built', False) == 0


def test_shellcheck_and_regen_skips_stay_allowed(tmp_path):
    assert _run(tmp_path, 'shellcheck not installed', True) == 0
    assert _run(tmp_path, 'set DOM1_REGEN=1 (about 35 s)', True) == 0
