"""For a triangulation T of the sphere (n vertices), maximize the minimum edge angle over convex
realizations inscribed in the unit sphere whose boundary complex is T (faces may become coplanar in the limit).
Variables: X (n x 3), c = max edge cosine.  minimize c  s.t.  |x_i|=1,  <x_i,x_j> <= c on edges,
det(x_j-x_i, x_k-x_i, x_l-x_i) <= 0 for every oriented face (i,j,k) and every other vertex l."""
import numpy as np, networkx as nx
from scipy.optimize import minimize
def oriented_faces(F, X0):
    """Consistent combinatorial orientation from a planar embedding; global sign chosen by the start X0."""
    G=nx.Graph()
    for a,b,c in F: G.add_edges_from([(a,b),(b,c),(a,c)])
    ok,emb=nx.check_planarity(G); faces=set(); out=[]
    for u,v in emb.edges():
        f=tuple(emb.traverse_face(u,v))
        if len(f)==3:
            m=f.index(min(f)); f=f[m:]+f[:m]
            if f not in faces: faces.add(f); out.append(f)
    vol=sum(np.cross(X0[j]-X0[i],X0[k]-X0[i])@(X0[i]) for i,j,k in out)
    if vol<0: out=[(i,k,j) for i,j,k in out]
    return out
def tutte_sphere(F, rng):
    G=nx.Graph()
    for a,b,c in F: G.add_edges_from([(a,b),(b,c),(a,c)])
    n=G.number_of_nodes(); outer=F[rng.integers(len(F))]
    pos={}
    for m,v in enumerate(outer): pos[v]=np.array([np.cos(2*np.pi*m/3),np.sin(2*np.pi*m/3)])*3.0
    inner=[v for v in G if v not in outer]; idx={v:i for i,v in enumerate(inner)}
    A=np.zeros((len(inner),len(inner))); B=np.zeros((len(inner),2))
    for v in inner:
        A[idx[v],idx[v]]=G.degree(v)
        for w in G[v]:
            if w in idx: A[idx[v],idx[w]]-=1
            else: B[idx[v]]+=pos[w]
    P=np.linalg.solve(A,B)
    for v in inner: pos[v]=P[idx[v]]
    X=np.zeros((n,3))
    for v,p in pos.items():
        s=p@p; X[v]=[2*p[0]/(1+s),2*p[1]/(1+s),(s-1)/(1+s)]
    # relax: spring toward equal edge lengths + repulsion
    E=list(G.edges())
    for it in range(400):
        D=np.zeros_like(X)
        for i,j in E:
            d=X[j]-X[i]; D[i]+=0.1*d; D[j]-=0.1*d
        for i in range(n):
            diff=X[i]-X; dist=np.linalg.norm(diff,axis=1)+1e-9; dist[i]=np.inf
            D[i]+=0.02*(diff/dist[:,None]**3).sum(0)
        X=X+0.5*D; X/=np.linalg.norm(X,axis=1,keepdims=True)
    return X, E
def theta_star(F, rng, starts=3, X0=None):
    """Largest min edge angle (degrees) over convex realizations inscribed in S^2 with boundary complex F;
    returns (theta, X) or None. Local optimization: a LOWER bound for theta*(F)."""
    F=[tuple(t) for t in F]; best=None; n=1+max(max(t) for t in F)
    for s in range(starts):
        if X0 is not None and s==0: X=X0.copy(); G=nx.Graph(); [G.add_edges_from([(a,b),(b,c),(a,c)]) for a,b,c in F]; E=list(G.edges())
        else: X,E=tutte_sphere(F,rng)
        OF=oriented_faces(F,X)
        others=[(f,l) for f in OF for l in range(n) if l not in f]
        I=np.array([f[0] for f,l in others]); J=np.array([f[1] for f,l in others]); K=np.array([f[2] for f,l in others]); L=np.array([l for f,l in others])
        Ei=np.array([e[0] for e in E]); Ej=np.array([e[1] for e in E])
        def obj(z): return z[-1]
        def obj_j(z): g=np.zeros_like(z); g[-1]=1; return g
        def sph(z): x=z[:-1].reshape(n,3); return (x*x).sum(1)-1
        def sph_j(z):
            x=z[:-1].reshape(n,3); J_=np.zeros((n,3*n+1))
            for i in range(n): J_[i,3*i:3*i+3]=2*x[i]
            return J_
        def edg(z): x=z[:-1].reshape(n,3); return z[-1]-(x[Ei]*x[Ej]).sum(1)
        def edg_j(z):
            x=z[:-1].reshape(n,3); J_=np.zeros((len(Ei),3*n+1))
            for r,(i,j) in enumerate(zip(Ei,Ej)): J_[r,3*i:3*i+3]=-x[j]; J_[r,3*j:3*j+3]=-x[i]
            J_[:,-1]=1; return J_
        def cvx(z):
            x=z[:-1].reshape(n,3); a,b,c,d=x[I],x[J],x[K],x[L]
            return -np.einsum('ij,ij->i',np.cross(b-a,c-a),d-a)
        def cvx_j(z):
            x=z[:-1].reshape(n,3); a,b,c,d=x[I],x[J],x[K],x[L]
            gb=np.cross(c-a,d-a); gc=np.cross(d-a,b-a); gd=np.cross(b-a,c-a); ga=-(gb+gc+gd)
            J_=np.zeros((len(I),3*n+1))
            for r in range(len(I)):
                J_[r,3*I[r]:3*I[r]+3]-=ga[r]; J_[r,3*J[r]:3*J[r]+3]-=gb[r]; J_[r,3*K[r]:3*K[r]+3]-=gc[r]; J_[r,3*L[r]:3*L[r]+3]-=gd[r]
            return J_
        # phase 1: push toward strict convexity
        def p1obj(w): return -w[-1]
        def p1j(w): g=np.zeros_like(w); g[-1]=-1; return g
        def p1c(w): return cvx(np.concatenate([w[:-1],[0]]))-w[-1]
        def p1cj(w):
            J_=np.zeros((len(I),3*n+1)); J_[:,:3*n]=cvx_j(np.concatenate([w[:-1],[0]]))[:,:3*n]; J_[:,-1]=-1; return J_
        w0=np.concatenate([X.ravel(),[min(cvx(np.concatenate([X.ravel(),[0]])).min(),0)]])
        r1=minimize(p1obj,w0,jac=p1j,constraints=[dict(type='eq',fun=lambda w:sph(np.concatenate([w[:-1],[0]])),jac=lambda w:np.hstack([sph_j(np.concatenate([w[:-1],[0]]))[:,:3*n],np.zeros((n,1))])),dict(type='ineq',fun=p1c,jac=p1cj),dict(type='ineq',fun=lambda w:0.05-w[-1])],method='SLSQP',options=dict(maxiter=300))
        if r1.x[-1]<=1e-9: continue
        X=r1.x[:-1].reshape(n,3); X=X/np.linalg.norm(X,axis=1,keepdims=True)
        z0=np.concatenate([X.ravel(),[max((X[Ei]*X[Ej]).sum(1))]])
        cons=[dict(type='eq',fun=sph,jac=sph_j),dict(type='ineq',fun=edg,jac=edg_j),dict(type='ineq',fun=cvx,jac=cvx_j)]
        res=minimize(obj,z0,jac=obj_j,constraints=cons,method='SLSQP',options=dict(maxiter=500,ftol=1e-12))
        z=res.x; x=z[:-1].reshape(n,3)
        feas=np.abs(sph(z)).max()<1e-7 and edg(z).min()>-1e-7 and cvx(z).min()>-1e-7
        if feas:
            th=np.degrees(np.arccos(min(1,z[-1])))
            if best is None or th>best[0]: best=(th,x.copy())
    return best


# ------------------------------------------------------------------ enumeration of triangulations
import itertools as _it, collections as _co

def tri_faces(G):
    ok, emb = nx.check_planarity(G); faces = set()
    for u, v in emb.edges():
        f = emb.traverse_face(u, v)
        if len(f) == 3: faces.add(frozenset(f))
    return faces

def graph_of(F):
    G = nx.Graph()
    for t in F:
        a, b, c = tuple(t); G.add_edges_from([(a, b), (b, c), (a, c)])
    return G

def flips(F):
    G = graph_of(F); ef = _co.defaultdict(list)
    for t in F:
        for e in _it.combinations(sorted(t), 2): ef[e].append(t)
    for (a, b), (t1, t2) in ef.items():
        c = next(iter(t1 - {a, b})); d = next(iter(t2 - {a, b}))
        if G.has_edge(c, d) or G.degree(a) <= 3 or G.degree(b) <= 3: continue
        F2 = set(F); F2 -= {t1, t2}; F2 |= {frozenset((a, c, d)), frozenset((b, c, d))}
        yield frozenset(F2)

def start_triangulation(n, seed=0):
    from scipy.spatial import ConvexHull
    X = np.random.default_rng(seed).normal(size=(n, 3)); X /= np.linalg.norm(X, axis=1, keepdims=True)
    return frozenset(frozenset(map(int, s)) for s in ConvexHull(X).simplices)

def enumerate_triangulations(n, start=None, verbose=False):
    """All combinatorial types of triangulations of S^2 with n vertices, by BFS over edge flips
    (the flip graph is connected: Wagner 1936). Order = BFS order from `start` (default: random hull).
    Known counts (OEIS A000109): n=4..12 -> 1,1,2,5,14,50,233,1249,7595."""
    if start is None:
        start = start_triangulation(n)
    key = lambda G: (tuple(sorted(d for _, d in G.degree())), nx.weisfeiler_lehman_graph_hash(G, iterations=3))
    buckets = _co.defaultdict(list); G0 = graph_of(start); buckets[key(G0)].append(G0)
    reps = [start]; queue = _co.deque([start])
    while queue:
        F = queue.popleft()
        for F2 in flips(F):
            G2 = graph_of(F2); k = key(G2)
            if any(nx.is_isomorphic(G2, H) for H in buckets[k]): continue
            buckets[k].append(G2); reps.append(F2); queue.append(F2)
        if verbose and len(reps) % 500 == 0: print(len(reps), flush=True)
    return [[sorted(t) for t in F] for F in reps]
