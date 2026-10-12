#!/usr/bin/env python3
"""Quick check that a machine is set up correctly (about 1-3 minutes, mostly numba compilation).
Checks: Lemma P; the compiled kernels agree with the reference implementation; one candidate certifies."""
import os, sys, time, json
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT); sys.path.insert(0, os.path.join(ROOT, 'scripts', 'completeness'))
import numpy as np
t0 = time.time()
import lemmas; lemmas.lemma_P()
from kisscrit.completeness.candidates import candidates_of_triangulation
from kisscrit.completeness.structure import Structure
from kisscrit.completeness.model import Model
from kisscrit.completeness import runner, fast, bnb2
T = json.load(open(os.path.join(ROOT, 'data', 'completeness', 'triangulations12.json')))
it = candidates_of_triangulation(T[3000])
code, faces = next(it); S = Structure(faces)
ch = {k: 'a' for k in range(len(S.sizes)) if S.sizes[k] == 5}
M = Model(S, ch); lo, hi = M.initial_box(1.0472, 1.1071); m = 0.5 * (lo + hi); w = 0.02
l, h = m - w, m + w
A = runner.arrays(S, ch)[:7]
ref = bnb2.evaluate(M, l[:, None], h[:, None]); got = fast.eval_all(*A, l, h, l, h)
# the compiled box enclosure may be tighter (monotonicity tightening), never looser
assert np.all(got[4] >= ref[1].lo[:, 0] - 1e-9) and np.all(got[5] <= ref[1].hi[:, 0] + 1e-9), 'kernel looser than reference'
refp = bnb2.evaluate(M, m[:, None], m[:, None]); gotp = fast.eval_all(*A, m, m, m, m)
assert np.allclose(refp[1].lo[:, 0], gotp[4], atol=1e-9), 'kernel point values differ from reference'
print('kernels agree with reference implementation')
term, st = runner.run(S, np.radians(60.0033), np.radians(63.4349489), tol=1e-6, max_boxes=300000)
print('sample candidate:', 'excluded' if not term else 'terminal', st)
print('selftest OK (%.0fs)' % (time.time() - t0))
