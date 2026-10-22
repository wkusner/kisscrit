"""Smoke test of an EXPLORATORY script: scripts/balls/ball_kill_test.py (BALL-KILL, 2026-10-21).

Status: untrusted (a test of an exploratory measurement script, not of the proof or the trusted base).

Checks, at the E = 19 saddle r = 1.0002111636 (coindex 3): (1) the floating-point developing map, fed the census
configuration's own face parameters in the chart of its terminal boxes (cond shard 32, line 1187, chart 2),
closes and lands back on the census configuration; (2) the Krawczyk test of kisscrit/certify.py's square KKT
system, run in Arb by the script, accepts a small radius and rejects a large one.  Under 60 s.
"""
import json, os, subprocess, sys

import numpy as np
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'scripts', 'balls'))
import ball_kill_test as B  # noqa: E402

SHARD, LINE, CHART = 'certificates-cond-32-of-40.jsonl', 1187, 2


@pytest.fixture(autouse=True)
def _restore_global_precision():
    """The script sets the process-wide precision of mpmath (`mp.mp.dps`, in load_census and KrawczykAt) and of
    python-flint (`flint.ctx.dps`, in KrawczykAt: scripts/balls/ball_kill_test.py, line 392) and never restores it,
    and `accepts()` relies on it staying set, so the script itself is left as it is.  This fixture saves both
    precisions (in bits) before each test of this module and restores them afterwards, so that tests which run
    later in the same process (tests/test_carith.py, test_cmodel.py, test_ckernel.py, which compare the Python
    reference with the C kernel at the default precision) are not affected (DOM-1, part 5)."""
    import flint, mpmath
    mp_prec, fl_prec = mpmath.mp.prec, flint.ctx.prec
    try:
        yield
    finally:
        mpmath.mp.prec = mp_prec; flint.ctx.prec = fl_prec


def _census():
    P = B.load_census(only=B.SMOKE_R)
    assert len(P) == 1 and P[0]['k'] == 19 and P[0]['coindex'] == 3
    return P[0]


def test_ball_kill_smoke():
    P = _census()
    with open(os.path.join(B.CERTS, SHARD)) as fh:
        for ln, line in enumerate(fh, 1):
            if ln == LINE:
                cert = json.loads(line)['certificate']; break
    ch = cert['charts'][CHART]
    S, M = B.structure_and_model(cert['faces'], ch['modes'])
    assert len(S.edges) == 19
    # label map and orientation: every graph isomorphism, both orientations; keep the best fit
    from networkx.algorithms import isomorphism
    GM = isomorphism.GraphMatcher(B._graph(P['A']), B._graph(S.edges))
    best = None
    for m in GM.subgraph_monomorphisms_iter():
        perm = np.array([m[i] for i in range(12)])
        inv = np.empty(12, int); inv[perm] = np.arange(12)
        for mirror in (False, True):
            Xc = P['X'][inv] * (np.array([1.0, 1.0, -1.0]) if mirror else 1.0)
            p = B.census_face_parameters(S, M, ch['modes'], Xc)
            C, F, _, _ = B.corners_at(M, p[:, None], p, p)
            Y, clo = B.develop(S, C, p[:1])
            m1, m2 = B.m1_m2(B.kabsch(Y[:, perm], P['X'], mirror), p[:1], P)
            if best is None or m1[0] < best[0]:
                best = (m1[0], m2[0], clo[0], np.abs(F).max())
    m1, m2, closure, vsum = best
    assert vsum < 1e-12, vsum            # the census parameters satisfy the vertex sums
    assert closure < 1e-9, closure       # the stop rule of the BALL-KILL prompt
    assert m1 < 1e-9 and m2 < 1e-9, (m1, m2)
    K = B.KrawczykAt(*P['mid'], P['A'])
    assert K.accepts(1e-6)               # a positive radius is accepted
    assert not K.accepts(1e-2)           # and a large one is not (the measured bracket is ~2.7e-4)


def test_ball_kill_does_not_leak_precision_into_later_tests():
    """Run the smoke test above and then one test of tests/test_carith.py in ONE pytest process: before the
    fixture, the carith test failed (5 parametrizations of test_unary_bit_identical, 7 failures over carith,
    cmodel and ckernel in the F8 session) because the global precision stayed at 80 digits.  The check is made in a
    subprocess so that it does not depend on the order of this module's tests in the running process."""
    from kisscrit.core import carith
    if not carith.available(): pytest.skip('C kernel not built (make -C csrc)')
    here = os.path.join(ROOT, 'tests')
    cmd = [sys.executable, '-m', 'pytest', '-q', '-p', 'no:cacheprovider',
           os.path.join(here, 'test_ball_kill.py') + '::test_ball_kill_smoke',
           os.path.join(here, 'test_carith.py') + '::test_unary_bit_identical']
    r = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True, timeout=600)
    assert r.returncode == 0 and ' passed' in r.stdout and 'failed' not in r.stdout, r.stdout[-1500:]
