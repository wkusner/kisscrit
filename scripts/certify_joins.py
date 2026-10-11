"""Certify the joins of the index-1 saddles of N=12 above r=1 (results/joins_n12.json)."""
import sys, os, json, time
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
import kisscrit as kc
from kisscrit.joins import push_on_variety, ascent_path, certify_side, max_cos
from kisscrit.geometry import to_ambient
R = os.path.join(os.path.dirname(__file__), '..', 'results')
D = json.load(open(os.path.join(os.path.dirname(__file__), '..', 'data', 'n12_above1.json')))['types']
cert = json.load(open(os.path.join(R, 'certified_n12.json')))
DOD = np.array(D[0]['X']); M = np.array(D[1]['X'])
Ms, *_ = kc.solve_kkt(M, [tuple(p) for p in D[1]['pairs']])
refs = [('DOD', kc.normalize(DOD)), ('M', kc.normalize(Ms))]
out = {}
t0 = float(sys.argv[1]) if len(sys.argv) > 1 else 0.02
for t in D:
    if t['coindex'] != 1: continue
    key = '%.10f' % t['r']
    lo = float(cert[key]['c'].strip('[]').split(',')[0]) - 1e-15
    A = [tuple(p) for p in t['pairs']]
    Xs, c, s, nf = kc.solve_kkt(np.array(t['X']), A)
    a = kc.analyze(Xs)
    w = to_ambient(a['B'], a['asc'][:, 0])
    res = dict(r=t['r'], c_lower=lo)
    for sign, side in ((1, 'plus'), (-1, 'minus')):
        tt = t0
        while True:
            Y = push_on_variety(Xs, A, sign * w, tt); cy = max_cos(Y)
            if cy + 1e-13 < lo or tt < 1e-4: break
            tt /= 2
        assert cy + 1e-13 < lo, 'start point not above the saddle level'
        assert sign * ((Y - Xs).ravel() @ w.ravel()) > 0
        t1 = time.time(); path = ascent_path(Y)
        r = certify_side(path, lo, refs)
        r.update(start_margin=float(lo - cy), t0=tt, time=round(time.time() - t1, 1))
        res[side] = r
        print(key, side, '->', r['end'], 'segments', r['segments'], 'checks', r['checks'],
              'start margin %.2e' % r['start_margin'], flush=True)
    out[key] = res
json.dump(out, open(os.path.join(R, 'joins_n12.json'), 'w'), indent=1)
