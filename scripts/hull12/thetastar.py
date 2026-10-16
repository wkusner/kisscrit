"""For a triangulation T of the sphere (12 vertices), maximize the minimum edge angle over convex
realizations inscribed in the unit sphere whose boundary complex is T (faces may become coplanar in the limit).
Variables: X (12x3), c = max edge cosine.  minimize c  s.t.  |x_i|=1,  <x_i,x_j> <= c on edges,
det(x_j-x_i, x_k-x_i, x_l-x_i) <= 0 for every oriented face (i,j,k) and every other vertex l.

Status: untrusted (exploration or tooling; no part of the proof).
"""
import numpy as np, networkx as nx, json, sys, time
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
    F=[tuple(t) for t in F]; best=None
    for s in range(starts):
        if X0 is not None and s==0: X=X0.copy(); G=nx.Graph(); [G.add_edges_from([(a,b),(b,c),(a,c)]) for a,b,c in F]; E=list(G.edges())
        else: X,E=tutte_sphere(F,rng)
        OF=oriented_faces(F,X)
        others=[(f,l) for f in OF for l in range(12) if l not in f]
        I=np.array([f[0] for f,l in others]); J=np.array([f[1] for f,l in others]); K=np.array([f[2] for f,l in others]); L=np.array([l for f,l in others])
        Ei=np.array([e[0] for e in E]); Ej=np.array([e[1] for e in E])
        def obj(z): return z[-1]
        def obj_j(z): g=np.zeros_like(z); g[-1]=1; return g
        def sph(z): x=z[:-1].reshape(12,3); return (x*x).sum(1)-1
        def sph_j(z):
            x=z[:-1].reshape(12,3); J_=np.zeros((12,37))
            for i in range(12): J_[i,3*i:3*i+3]=2*x[i]
            return J_
        def edg(z): x=z[:-1].reshape(12,3); return z[-1]-(x[Ei]*x[Ej]).sum(1)
        def edg_j(z):
            x=z[:-1].reshape(12,3); J_=np.zeros((len(Ei),37))
            for r,(i,j) in enumerate(zip(Ei,Ej)): J_[r,3*i:3*i+3]=-x[j]; J_[r,3*j:3*j+3]=-x[i]
            J_[:,-1]=1; return J_
        def cvx(z):
            x=z[:-1].reshape(12,3); a,b,c,d=x[I],x[J],x[K],x[L]
            return -np.einsum('ij,ij->i',np.cross(b-a,c-a),d-a)
        def cvx_j(z):
            x=z[:-1].reshape(12,3); a,b,c,d=x[I],x[J],x[K],x[L]
            gb=np.cross(c-a,d-a); gc=np.cross(d-a,b-a); gd=np.cross(b-a,c-a); ga=-(gb+gc+gd)
            J_=np.zeros((len(I),37))
            for r in range(len(I)):
                J_[r,3*I[r]:3*I[r]+3]-=ga[r]; J_[r,3*J[r]:3*J[r]+3]-=gb[r]; J_[r,3*K[r]:3*K[r]+3]-=gc[r]; J_[r,3*L[r]:3*L[r]+3]-=gd[r]
            return J_
        # phase 1: push toward strict convexity
        def p1obj(w): return -w[-1]
        def p1j(w): g=np.zeros_like(w); g[-1]=-1; return g
        def p1c(w): return cvx(np.concatenate([w[:-1],[0]]))-w[-1]
        def p1cj(w):
            J_=np.zeros((len(I),37)); J_[:,:36]=cvx_j(np.concatenate([w[:-1],[0]]))[:,:36]; J_[:,-1]=-1; return J_
        w0=np.concatenate([X.ravel(),[min(cvx(np.concatenate([X.ravel(),[0]])).min(),0)]])
        r1=minimize(p1obj,w0,jac=p1j,constraints=[dict(type='eq',fun=lambda w:sph(np.concatenate([w[:-1],[0]])),jac=lambda w:np.hstack([sph_j(np.concatenate([w[:-1],[0]]))[:,:36],np.zeros((12,1))])),dict(type='ineq',fun=p1c,jac=p1cj),dict(type='ineq',fun=lambda w:0.05-w[-1])],method='SLSQP',options=dict(maxiter=300))
        if r1.x[-1]<=1e-9: continue
        X=r1.x[:-1].reshape(12,3); X=X/np.linalg.norm(X,axis=1,keepdims=True)
        z0=np.concatenate([X.ravel(),[max((X[Ei]*X[Ej]).sum(1))]])
        cons=[dict(type='eq',fun=sph,jac=sph_j),dict(type='ineq',fun=edg,jac=edg_j),dict(type='ineq',fun=cvx,jac=cvx_j)]
        res=minimize(obj,z0,jac=obj_j,constraints=cons,method='SLSQP',options=dict(maxiter=500,ftol=1e-12))
        z=res.x; x=z[:-1].reshape(12,3)
        feas=np.abs(sph(z)).max()<1e-7 and edg(z).min()>-1e-7 and cvx(z).min()>-1e-7
        if feas:
            th=np.degrees(np.arccos(min(1,z[-1])))
            if best is None or th>best[0]: best=(th,x.copy())
    return best
if __name__=='__main__':
    T=json.load(open('tri12.json')); lo,hi=int(sys.argv[1]),int(sys.argv[2]); out=sys.argv[3]
    rng=np.random.default_rng(lo); res={}
    t0=time.time()
    for k in range(lo,min(hi,len(T))):
        b=theta_star(T[k],rng,starts=1)
        if b is None or b[0]>56.0:
            b2=theta_star(T[k],rng,starts=4)
            if b2 is not None and (b is None or b2[0]>b[0]): b=b2
        res[k]=None if b is None else dict(theta=b[0],X=b[1].tolist())
        if (k-lo)%100==0:
            json.dump(res,open(out,'w')); print(k,round(time.time()-t0),'s', None if b is None else round(b[0],3),flush=True)
    json.dump(res,open(out,'w'))
