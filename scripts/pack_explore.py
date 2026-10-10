"""Explore critical configurations of equal-ball packings in a container.

usage: python scripts/pack_explore.py CONTAINER n ascents kkt_starts
CONTAINER: sq (disks in unit square), cube (balls in unit cube), tor (square flat torus),
           hex (hexagonal flat torus), tor3 (cubic 3-torus), disk (disks in a disk), ball (balls in a ball)
Writes results/pack_CONTAINER_n{n}.json.
"""
import sys, os, time
import numpy as np
import mpmath as mp
import warnings
warnings.filterwarnings('ignore')
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from kisscrit.packing import Cube, Torus, Ball, explore, euler_target, surgery

which, n, asc, ks = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
cont = {'sq': Cube(2), 'cube': Cube(3), 'tor': Torus(np.eye(2)), 'hex': Torus([[1, 0], [0.5, 3 ** 0.5 / 2]], exact=lambda: [[mp.mpf(1), mp.mpf(0)], [mp.mpf(1) / 2, mp.sqrt(3) / 2]]),
        'disk': Ball(2), 'ball': Ball(3), 'tor3': Torus(np.eye(3))}[which]
t0 = time.time(); C = explore(cont, n, np.random.default_rng(n), asc, ks); surgery(C)
print(which, 'n', n, 'types', len(C.types), 'signed', C.signed_count(), 'target', euler_target(cont, n),
      '%.0fs' % (time.time() - t0), flush=True)
for t in C.sorted():
    print('  rho=%.10f k=%d coidx=%d nzero=%d ratt=%s labeled=%s family=%s hits=%d' % (
        t['rho'], t['k'], t['coindex'], t['nzero'], t['rattlers'], t['labeled'], t.get('family_count'), t['hits']))
C.dump(os.path.join(os.path.dirname(__file__), '..', 'results', 'pack_%s_n%d.json' % (which, n)))
