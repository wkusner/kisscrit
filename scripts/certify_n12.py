"""Interval-certify the critical configurations of 12 points above r=1 (data/n12_above1.json).

Status: untrusted (exploration or tooling; no part of the proof).

"""
import sys, os, json, time
import numpy as np, mpmath as mp
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from kisscrit.exact import refine_kkt
from kisscrit.kkt import solve_kkt
from kisscrit.certify import certify
D = json.load(open(os.path.join(os.path.dirname(__file__), '..', 'data', 'n12_above1.json')))
which = [float(a) for a in sys.argv[1:]]
out = {}
for t in D['types']:
    if which and not any(abs(t['r'] - w) < 1e-6 for w in which): continue
    A = [tuple(p) for p in t['pairs']]; X = np.array(t['X'])
    Xs, c, s, nf = solve_kkt(X, A)
    t0 = time.time()
    x, cc, sig, res = refine_kkt(Xs, A, s, dps=90)
    mp.mp.dps = 90
    C = [[sum(x[i][a] * x[j][a] for a in range(3)) for j in range(12)] for i in range(12)]
    mu = [sum(sig[e] * C[i][j if i == a else a] for e, (a, j) in enumerate(A) if i in (a, j)) for i in range(12)]
    try:
        R = certify(x, cc, sig, mu, A, dps=80)
        R['time'] = round(time.time() - t0, 1); R['k'] = len(A)
        print('r=%.10f k=%d coindex(certified)=%d (numerical %d)  %.0fs' % (t['r'], len(A), R['coindex'], t['coindex'], R['time']), flush=True)
        out['%.10f' % t['r']] = R
    except (AssertionError, ZeroDivisionError) as e:
        print('r=%.10f FAILED: %s' % (t['r'], e), flush=True)
        out['%.10f' % t['r']] = dict(failed=str(e))
prev = {}
pth = os.path.join(os.path.dirname(__file__), '..', 'results', 'certified_n12.json')
if os.path.exists(pth): prev = json.load(open(pth))
prev.update(out); out = prev
json.dump(out, open(os.path.join(os.path.dirname(__file__), '..', 'results', 'certified_n12.json'), 'w'), indent=1)
