"""Local maxima of the packing radius by random ascent; best value, number of distinct maxima, and an LLL
minimal polynomial of the best radius.  usage: python scripts/pack_maxima.py CONTAINER n_lo n_hi ascents"""
import sys, os, json, time, warnings
import numpy as np
import mpmath as mp, mpmath as mp
warnings.filterwarnings('ignore')
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from kisscrit.packing import *
from kisscrit.exact import minpoly
CONT = {'sq': Cube(2), 'cube': Cube(3), 'tor': Torus(np.eye(2)), 'hex': Torus([[1, 0], [0.5, 3 ** 0.5 / 2]], exact=lambda: [[mp.mpf(1), mp.mpf(0)], [mp.mpf(1) / 2, mp.sqrt(3) / 2]]),
        'disk': Ball(2), 'ball': Ball(3), 'tor3': Torus(np.eye(3))}
w, lo, hi, asc = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
cont = CONT[w]; out = {}
path = os.path.join(os.path.dirname(__file__), '..', 'results', 'maxima_%s.json' % w)
if os.path.exists(path): out = json.load(open(path))
for n in range(lo, hi + 1):
    rng = np.random.default_rng(1000 + n); C = PCensus(cont, n); t0 = time.time()
    for _ in range(asc):
        X, r = ascend(cont, cont.random(n, rng))
        for tol in (1e-9, 1e-7, 1e-5):
            A, t = active(cont, X, tol)
            Xs, t2, s, nf = kkt_solve(cont, X, A)
            if kkt_check(cont, Xs, t2, s, A, nf):
                C.add(cont.wrap(Xs), 'ascent'); break
    maxima = [t for t in C.sorted() if t['coindex'] == 0]
    if not maxima:
        print(w, n, 'no maxima polished'); continue
    best = maxima[0]
    A = [tuple(c) if c[0] != 'pair' else ('pair', c[1], c[2], tuple(c[3]) if c[3] is not None else None) for c in best['active']]
    A = [('wall', c[1], c[2], c[3]) if c[0] == 'wall' else (('sphere', c[1]) if c[0] == 'sphere' else c) for c in A]
    poly = None
    try:
        Xh, th, sh, res = refine_mp(cont, np.array(best['X']), A, np.array(best['sigma']), dps=120)
        if res < mp.mpf(10) ** -100:
            mp.mp.dps = 120; poly = minpoly(th, maxdeg=24, digits=110)
    except Exception as e:
        print('refine failed', e)
    out[str(n)] = dict(best_rho=best['rho'], k=best['k'], rattlers=best['rattlers'], nzero=best['nzero'],
                       n_maxima=len(maxima), hits=[m['hits'] for m in maxima], second=maxima[1]['rho'] if len(maxima) > 1 else None,
                       poly=poly, labeled=best.get('labeled'), X=best['X'])
    print('%s n=%2d best rho=%.12f  k=%d ratt=%s nzero=%d  distinct maxima=%d  hits(best)=%d  poly(deg %s)=%s  %.0fs' % (
        w, n, best['rho'], best['k'], best['rattlers'], best['nzero'], len(maxima), best['hits'],
        (len(poly) - 1) if poly else '-', poly, time.time() - t0), flush=True)
    json.dump(out, open(path, 'w'), indent=1)
