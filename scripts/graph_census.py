"""Graph-mode census for n points on S^2: solve the square KKT system for every graph on n vertices
(minimum degree >= 2, n <= 7) from random and graph-seeded starts.

usage: python scripts/graph_census.py n starts out.json
"""
import sys, os, time
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
import kisscrit as kc
from kisscrit.search import graphs_on, solve_graph

n, starts, out = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
rng = np.random.default_rng(0); C = kc.Census(n, 3); t0 = time.time()
Gs = graphs_on(n); print(len(Gs), 'graphs', flush=True)
for gi, G in enumerate(Gs):
    for X in solve_graph(G, 3, rng, starts):
        C.add(X, 'graph')
    if gi % 25 == 0:
        print(gi, 'graphs', len(C.types), 'types %.0fs' % (time.time() - t0), flush=True)
        C.dump(out)
C.dump(out)
print('types', len(C.types), 'time %.0fs' % (time.time() - t0))
for t in C.sorted():
    print('  r=%.9f k=%d coidx=%d nzero=%d kern=%d rattlers=%d sym=%s labeled=%s hits=%d' % (
        t['r'], t['k'], t['coindex'], t['nzero'], t['kernel'], t['rattlers'], t.get('sym'), t.get('labeled'), t['hits']))
