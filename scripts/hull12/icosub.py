"""icosub.py: see the module.

Status: untrusted (exploration or tooling; no part of the proof).

"""
import networkx as nx, itertools, numpy as np, collections, sys
from networkx.algorithms import isomorphism as iso
G=nx.icosahedral_graph(); E=sorted(tuple(sorted(e)) for e in G.edges()); idx={e:i for i,e in enumerate(E)}
auts=list(iso.GraphMatcher(G,G).isomorphisms_iter()); print('aut',len(auts))
perms=[[idx[tuple(sorted((a[u],a[v])))] for (u,v) in E] for a in auts]
# enumerate deletion sets D (bitmask) with max degree <= 2 in the deleted subgraph, |D| <= 11
res=collections.Counter(); reps=collections.defaultdict(set)
deg=[0]*12
def canon(mask):
    best=None
    for p in perms:
        m=0
        for i in range(30):
            if mask>>i&1: m|=1<<p[i]
        if best is None or m<best: best=m
    return best
cnt=collections.Counter()
def dfs(i,mask,size):
    if size>=7: cnt[size]+=1; reps[size].add(canon(mask)) if size<=11 else None
    if size==11 or i==30: return
    for j in range(i,30):
        u,v=E[j]
        if deg[u]<2 and deg[v]<2:
            deg[u]+=1; deg[v]+=1
            dfs(j+1,mask|1<<j,size+1)
            deg[u]-=1; deg[v]-=1
dfs(0,0,0)
for s in sorted(cnt): print('deleted',s,'edges kept',30-s,'labeled',cnt[s],'up to symmetry',len(reps[s]))
