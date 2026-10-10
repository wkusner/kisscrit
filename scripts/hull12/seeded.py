import json,numpy as np,networkx as nx,sys
from scipy.spatial import ConvexHull
from thetastar import theta_star
T=json.load(open('tri12.json'))
def G_of(F):
    G=nx.Graph()
    for a,b,c in F: G.add_edges_from([(a,b),(b,c),(a,c)])
    return G
d=json.load(open('../crit12/explore2.json')); rng=np.random.default_rng(0)
out={}
for e in sorted(d,key=lambda e:-e['r'])[1:]:
    X=np.array(e['X']); h=ConvexHull(X); H=nx.Graph()
    for s in h.simplices: H.add_edges_from([(s[0],s[1]),(s[1],s[2]),(s[0],s[2])])
    for k in [1,5,7]:
        m=nx.vf2pp_isomorphism(H,G_of(T[k]))
        if m: break
    X0=np.zeros((12,3))
    for i,j in m.items(): X0[j]=X[i]
    b=theta_star(T[k],rng,starts=1,X0=X0)
    print(round(e['r'],8),'T',k,'seeded theta*', None if b is None else round(b[0],5),flush=True)
    if b and (k not in out or b[0]>out[k][0]): out[k]=(b[0],b[1].tolist())
json.dump(out,open('seeded.json','w'))
