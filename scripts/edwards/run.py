"""Basin statistics for free planar clusters modulo similarity.  Usage: python run.py N runs {morse|slsqp} seed"""
import sys, json, time, numpy as np, collections, networkx as nx
import os; HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, '..', '..', 'results', 'edwards')
sys.path.insert(0, HERE); from shape import *
N = int(sys.argv[1]); M = int(sys.argv[2]); proto = sys.argv[3]; seed = int(sys.argv[4])
rng = np.random.default_rng(seed); types = {}; counts = collections.Counter(); t0 = time.time()
for s in range(M):
    x = rng.normal(size=(N, 2))
    y = ascend(x) if proto == 'morse' else x
    y = polish(y)
    r = rho(y); C = contacts(y)
    G = nx.Graph(C); G.add_nodes_from(range(N))
    key = None
    for k, v in types.items():
        if abs(v['rho'] - r) < 1e-6 and v['E'] == len(C) and nx.is_isomorphic(G, v['G']): key = k; break
    if key is None:
        key = len(types)
        types[key] = dict(rho=float(r), aut=aut_order(y), E=len(C), G=G, x=normalize(y).tolist(), ismax=is_max(y))
    counts[key] += 1
out = []
for k, v in sorted(types.items(), key=lambda kv: -kv[1]['rho']):
    v['count'] = counts[k]; v.pop('G'); out.append(v)
json.dump(dict(N=N, M=M, proto=proto, types=out, seconds=time.time() - t0), open(os.path.join(OUT, 'res_%s_%d.json' % (proto, N)), 'w'))
print(N, proto, '%.0fs' % (time.time() - t0))
for v in out: print('  rho %.6f  |Aut| %2d  E %2d  max %s  count %d' % (v['rho'], v['aut'], v['E'], v['ismax'], v['count']))
