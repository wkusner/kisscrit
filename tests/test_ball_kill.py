"""Smoke test of an EXPLORATORY script: scripts/balls/ball_kill_test.py (BALL-KILL, 2026-10-21).

Status: untrusted (a test of an exploratory measurement script, not of the proof or the trusted base).

Checks, at the E = 19 saddle r = 1.0002111636 (coindex 3): (1) the floating-point developing map, fed the census
configuration's own face parameters in the chart of its terminal boxes (cond shard 32, line 1187, chart 2),
closes and lands back on the census configuration; (2) the Krawczyk test of kisscrit/certify.py's square KKT
system, run in Arb by the script, accepts a small radius and rejects a large one.  Under 60 s.
"""
import json, os, sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'scripts', 'balls'))
import ball_kill_test as B  # noqa: E402

SHARD, LINE, CHART = 'certificates-cond-32-of-40.jsonl', 1187, 2


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
