"""Euler-characteristic check for the packing censuses in results/pack_*.json, with one-rattler strata."""
import sys, os, json, glob, re, warnings
import numpy as np
import mpmath as mp
warnings.filterwarnings('ignore')
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from kisscrit.packing import Cube, Torus, Ball, PCensus, euler_with_rattlers, euler_target
CONT = {'sq': Cube(2), 'cube': Cube(3), 'tor': Torus(np.eye(2)), 'hex': Torus([[1, 0], [0.5, 3 ** 0.5 / 2]], exact=lambda: [[mp.mpf(1), mp.mpf(0)], [mp.mpf(1) / 2, mp.sqrt(3) / 2]]),
        'disk': Ball(2), 'ball': Ball(3), 'tor3': Torus(np.eye(3))}
R = os.path.join(os.path.dirname(__file__), '..', 'results')
data = {}
for p in glob.glob(os.path.join(R, 'pack_*_n*.json')):
    w, n = re.match(r'.*pack_(\w+?)_n(\d+)\.json', p).groups(); n = int(n)
    C = PCensus(CONT[w], n)
    for t in json.load(open(p))['types']:
        C.types[C.key(t)] = t
    data.setdefault(w, {})[n] = C
out = {}
for w in sorted(data):
    for n in sorted(data[w]):
        tot, sk = euler_with_rattlers(data[w], n)
        tgt = euler_target(CONT[w], n)
        out['%s_%d' % (w, n)] = dict(signed=tot, target=tgt, skipped=sk, types=len(data[w][n].types))
        print('%-5s n=%d  types=%2d  signed=%6s  target=%6s  %s  skipped=%s' % (
            w, n, len(data[w][n].types), tot, tgt, 'OK' if tot == tgt else '--', sk), flush=True)
json.dump(out, open(os.path.join(R, 'pack_euler.json'), 'w'), indent=1)
