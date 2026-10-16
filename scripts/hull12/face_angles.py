"""face_angles.py: see the module.

Status: untrusted (exploration or tooling; no part of the proof).

"""
import numpy as np, sys
th=float(sys.argv[1]); rho_max=float(sys.argv[2]); step=float(sys.argv[3])
g=np.radians(np.arange(th,120.0+1e-9,step))
a,b,c=np.meshgrid(g,g,g,indexing='ij'); a,b,c=a.ravel(),b.ravel(),c.ravel()
# vertices: P1 north, P2 at polar angle c, P3 polar angle b with azimuth from cos a
cosA=(np.cos(a)-np.cos(b)*np.cos(c))/(np.sin(b)*np.sin(c))
ok=np.abs(cosA)<1
a,b,c,cosA=a[ok],b[ok],c[ok],cosA[ok]
P1=np.stack([0*a,0*a,1+0*a],1); P2=np.stack([np.sin(c),0*a,np.cos(c)],1)
sA=np.sqrt(1-cosA**2); P3=np.stack([np.sin(b)*cosA,np.sin(b)*sA,np.cos(b)],1)
n=np.cross(P2-P1,P3-P1); d=np.abs(np.einsum('ij,ij->i',n,P1))/np.linalg.norm(n,axis=1)
R=np.degrees(np.arccos(np.clip(d,-1,1)))
keep=R<rho_max
A=np.degrees(np.arccos(cosA[keep]))
i=np.argmin(A); j=np.argmax(A)
print('theta',th,'rho<',rho_max,'triangles',keep.sum(),'min angle %.3f at sides'%A[i],np.degrees([a[keep][i],b[keep][i],c[keep][i]]).round(2),'R %.2f'%R[keep][i])
print('   max angle %.3f at sides'%A[j],np.degrees([a[keep][j],b[keep][j],c[keep][j]]).round(2),'R %.2f'%R[keep][j])
