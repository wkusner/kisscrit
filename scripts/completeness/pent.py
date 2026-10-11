import numpy as np
def rot(axis, ang):
    axis=axis/np.linalg.norm(axis); K=np.array([[0,-axis[2],axis[1]],[axis[2],0,-axis[0]],[-axis[1],axis[0],0]])
    return np.eye(3)+np.sin(ang)*K+(1-np.cos(ang))*K@K
def corner(P,Q,R):  # angle at Q between QP and QR
    u=P-(P@Q)*Q; v=R-(R@Q)*Q
    return np.arccos(np.clip(u@v/np.linalg.norm(u)/np.linalg.norm(v),-1,1))
def pentagons(th, n=400):
    out=[]
    P0=np.array([0,0,1.]); P1=rot(np.array([0,1.,0]),th)@P0
    for f1 in np.linspace(0.3,np.pi,n):
        # P2: rotate P0 about P1 by f1 (interior angle at P1)
        P2=rot(P1,f1)@P0
        for f2 in np.linspace(0.3,np.pi,n):
            P3=rot(P2,f2)@P1
            # P4 at distance th from P3 and P0
            c=np.cos(th); A=np.array([P3,P0]); 
            # solve x: P3.x=c, P0.x=c, |x|=1
            n_=np.cross(P3,P0); 
            if np.linalg.norm(n_)<1e-9: continue
            M=np.array([P3,P0,n_]); 
            x0=np.linalg.lstsq(np.array([P3,P0]),np.array([c,c]),rcond=None)[0]
            nn=n_/np.linalg.norm(n_); r2=1-x0@x0
            if r2<0: continue
            for s in (1,-1):
                P4=x0+s*np.sqrt(r2)*nn
                V=[P0,P1,P2,P3,P4]
                # orientation: all corners measured, require convex: turning same sign
                sg=[np.sign(np.dot(np.cross(V[i]-V[i-1],V[(i+1)%5]-V[i]),V[i])) for i in range(5)]
                if not (all(x>0 for x in sg) or all(x<0 for x in sg)): continue
                D=[np.arccos(np.clip(V[i]@V[j],-1,1)) for i in range(5) for j in range(i+2,5) if not (i==0 and j==4)]
                if min(D)<th-1e-12: continue
                ang=[corner(V[i-1],V[i],V[(i+1)%5]) for i in range(5)]
                out.append(ang)
    return np.degrees(np.array(out))
for thd in [60.0001,61,62,63,63.4349]:
    th=np.radians(thd); a=np.degrees(np.arccos(np.cos(th)/(1+np.cos(th))))
    A=pentagons(th,260)
    S=A.sum(1)
    print('th %.2f alpha %.3f  n %d  sum in [%.2f, %.2f]  sum-9alpha in [%.2f, %.2f]  corner in [%.2f, %.2f]'%(thd,a,len(A),S.min(),S.max(),S.min()-9*a,S.max()-9*a,A.min(),A.max()))
