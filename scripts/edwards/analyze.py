import json, numpy as np, sys
from scipy.optimize import minimize_scalar
from scipy.stats import chi2
for proto in ['morse','slsqp']:
  for N in range(6,11):
    try: d=json.load(open(__import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)),'..','..','results','edwards','res_%s_%d.json'%(proto,N))))
    except: continue
    T=[t for t in d['types'] if t['ismax']]
    n=np.array([t['count'] for t in T]); a=np.array([t['aut'] for t in T]); r=np.array([t['rho'] for t in T])
    tot=n.sum(); obs=n/tot
    flat=(1/a)/(1/a).sum()
    def nll(b):
        w=np.exp(b*(r-r.max())*100)/a; p=w/w.sum(); return -(n*np.log(p)).sum()
    b=minimize_scalar(nll,bounds=(-50,50),method='bounded').x
    w=np.exp(b*(r-r.max())*100)/a; fit=w/w.sum()
    g=lambda p: 2*(n*np.log(np.maximum(n,1e-300)/(tot*p))).sum()
    k=len(n)
    print('%s N=%d samples=%d types=%d'%(proto,N,tot,k))
    print('   observed ', np.round(obs,3)); print('   flat 1/|Aut|', np.round(flat,3), ' G=%.0f p=%.1e'%(g(flat), chi2.sf(g(flat),k-1)))
    if k>2: print('   compactivity fit beta=%.2f per 0.01 rho'%b, np.round(fit,3), ' G=%.1f p=%.2g'%(g(fit), chi2.sf(g(fit),k-2)))
