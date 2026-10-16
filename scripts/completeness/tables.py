"""tables.py: see the module.

Status: untrusted (exploration or tooling; no part of the proof).

"""
import numpy as np, json
from scipy.optimize import minimize
exec(open('pent.py').read().split('for thd in')[0])
def theta_of_alpha(ad):
    ca=np.cos(np.radians(ad)); c=ca/(1-ca); return np.arccos(c)
def g(th,A):
    cb=np.cos(th)**2+np.sin(th)**2*np.cos(A); b=np.arccos(cb)
    return 2*np.arccos((np.cos(th)-np.cos(th)*cb)/(np.sin(th)*np.sin(b)))
A60=np.degrees(np.arccos(1/3)); NC=16
edges=np.linspace(A60,72.0,NC+1)
rows=[]
for i in range(NC+1):
    a=edges[i]; th=theta_of_alpha(a); ar=np.radians(a)
    As=np.linspace(ar,2*ar,4001); s4=np.degrees((As+g(th,As)).max())-3*a
    P=pentagons(th,110); S=P.sum(1)
    rows.append(dict(alpha=a, theta=np.degrees(th), s4max=s4, s5min=S.min()-9*a, s5max=S.max()-9*a))
    print('%.4f %.4f s4max %.3f s5 [%.2f, %.2f]'%(a,np.degrees(th),s4,S.min()-9*a,S.max()-9*a),flush=True)
json.dump(rows,open('slack_tables.json','w'),indent=1)
