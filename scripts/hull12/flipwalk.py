"""Explore triangulations reachable above 60 deg by single flips, seeding theta* with the parent realization.

Status: untrusted (exploration or tooling; no part of the proof).

"""
import json, numpy as np, networkx as nx, collections, itertools, sys, pathlib
# Repository root = two directories above this file (scripts/hull12/), found from __file__ so that the
# script does not depend on where the repository is cloned.
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[2]))
import warnings; warnings.filterwarnings('ignore')
from thetastar import theta_star
from scipy.spatial import ConvexHull
T=json.load(open('tri12.json'))
def G_of(F):
    G=nx.Graph()
    for t in F: a,b,c=tuple(t); G.add_edges_from([(a,b),(b,c),(a,c)])
    return G
key=lambda G:(tuple(sorted(d for _,d in G.degree())),nx.weisfeiler_lehman_graph_hash(G,iterations=3))
buck=collections.defaultdict(list)
for k,F in enumerate(T): G=G_of(F); buck[key(G)].append((k,G))
def ident(F):
    G=G_of(F); return next(k for k,H in buck[key(G)] if nx.is_isomorphic(G,H))
def hull_faces(X): return [tuple(sorted(map(int,s))) for s in ConvexHull(X).simplices]
def flips(F):
    F=[frozenset(t) for t in F]; G=G_of(F); ef=collections.defaultdict(list)
    for t in F:
        for e in itertools.combinations(sorted(t),2): ef[e].append(t)
    for (a,b),(t1,t2) in ef.items():
        c=next(iter(t1-{a,b})); d=next(iter(t2-{a,b}))
        if G.has_edge(c,d) or G.degree(a)<=3 or G.degree(b)<=3: continue
        F2=set(F)-{t1,t2}|{frozenset((a,c,d)),frozenset((b,c,d))}
        yield [tuple(sorted(t)) for t in F2]
rng=np.random.default_rng(5)
seeds=json.load(open('seeded.json'))
best={}; queue=[]
ico=json.load(open('ts_a.json'))['0']
for k,X in [(0,ico['X'])]+[(int(k),v[1]) for k,v in seeds.items()]:
    X=np.array(X); F=hull_faces(X); kk=ident(F)
    b=theta_star(F,rng,starts=1,X0=X)
    if b: best[kk]=max(best.get(kk,0),b[0]); queue.append((F,b[1]))
print('seeds',best,flush=True)
seen=set()
while queue:
    F,X=queue.pop(0)
    for F2 in flips(F):
        b=theta_star(F2,rng,starts=1,X0=X)
        if b is None or b[0]<=60.0+1e-6: continue
        kk=ident(F2)
        sig=(kk,round(b[0],4))
        if sig in seen: continue
        seen.add(sig)
        new = kk not in best or b[0]>best[kk]+1e-5
        best[kk]=max(best.get(kk,0),b[0])
        print('T%d theta*>=%.5f'%(kk,b[0]),'NEW' if new else '',flush=True)
        if new: queue.append((F2,b[1]))
print('final',{k:round(v,5) for k,v in sorted(best.items())})
json.dump(best,open('flipwalk.json','w'))
