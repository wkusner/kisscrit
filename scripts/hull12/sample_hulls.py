"""Ascend random 12-point configurations to just above r=1 (min angle > 60 deg) and record hull types.

Status: untrusted (exploration or tooling; no part of the proof).

"""
import sys, json, collections, numpy as np, networkx as nx, pathlib
# Repository root = two directories above this file (scripts/hull12/), found from __file__ so that
# `import kisscrit` below works from any clone, wherever the repository is.
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[2]))
import warnings; warnings.filterwarnings('ignore')
import kisscrit as kc
from scipy.spatial import ConvexHull
T=json.load(open('tri12.json'))
def G_of(F):
    G=nx.Graph()
    for a,b,c in F: G.add_edges_from([(a,b),(b,c),(a,c)])
    return G
ref={k:G_of(T[k]) for k in range(60)}
seed=int(sys.argv[1]); n=int(sys.argv[2]); thr=float(sys.argv[3])
rng=np.random.default_rng(seed); cnt=collections.Counter(); ex={}
ct=np.cos(np.radians(thr))
for it in range(n):
    X,info=kc.ascend(kc.random_config(12,3,rng),maxit=6000,target_c=ct)
    if info['c']>ct: cnt['not reached']+=1; continue
    h=ConvexHull(X); H=nx.Graph()
    for s in h.simplices: H.add_edges_from([(s[0],s[1]),(s[1],s[2]),(s[0],s[2])])
    k=next((k for k,g in ref.items() if nx.is_isomorphic(g,H)),'other')
    cnt[k]+=1; ex.setdefault(str(k),X.tolist())
    print(it,k,round(info['r'],6),flush=True)
print(dict(cnt))
json.dump(dict(counts={str(k):v for k,v in cnt.items()},examples=ex),open('sample_hulls_%d.json'%seed,'w'))
