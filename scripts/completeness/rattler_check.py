import numpy as np
exec(open('pent.py').read().split('def pentagons')[0])
def pent_list(th,n):
    out=[]
    P0=np.array([0,0,1.]); P1=rot(np.array([0,1.,0]),th)@P0
    for f1 in np.linspace(0.3,np.pi,n):
        P2=rot(P1,f1)@P0
        for f2 in np.linspace(0.3,np.pi,n):
            P3=rot(P2,f2)@P1
            c=np.cos(th); n_=np.cross(P3,P0)
            if np.linalg.norm(n_)<1e-9: continue
            x0=np.linalg.lstsq(np.array([P3,P0]),np.array([c,c]),rcond=None)[0]; nn=n_/np.linalg.norm(n_); r2=1-x0@x0
            if r2<0: continue
            for s in (1,-1):
                P4=x0+s*np.sqrt(r2)*nn; V=[P0,P1,P2,P3,P4]
                sg=[np.sign(np.dot(np.cross(V[i]-V[i-1],V[(i+1)%5]-V[i]),V[i])) for i in range(5)]
                if not (all(x>0 for x in sg) or all(x<0 for x in sg)): continue
                D=[np.arccos(np.clip(V[i]@V[j],-1,1)) for i in range(5) for j in range(i+2,5) if not (i==0 and j==4)]
                if min(D)<th-1e-12: continue
                out.append(np.array(V))
    return out
rng=np.random.default_rng(0)
W=rng.dirichlet(np.ones(5)*0.7,size=4000)
for thd in [60.0001,61,62,63,63.4349]:
    th=np.radians(thd); worst=0
    for V in pent_list(th,70):
        P=W@V; P/=np.linalg.norm(P,axis=1,keepdims=True)
        md=np.degrees(np.arccos(np.clip(P@V.T,-1,1))).min(1).max()
        worst=max(worst,md)
    print('theta %.2f  max over pentagons of max-min distance to vertices: %.3f (need < theta)'%(thd,worst))
