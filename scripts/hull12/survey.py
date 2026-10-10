"""Resumable theta*(T) survey over all 7595 triangulations; worker w of W handles indices k with k%W==w.
Results appended line-by-line to ts_w{w}.jsonl."""
import json, sys, os, numpy as np, time
from thetastar import theta_star
w,W=int(sys.argv[1]),int(sys.argv[2]); out='ts_w%d.jsonl'%w
T=json.load(open('tri12.json'))
done=set()
for f in ['ts_a.json','ts_b.json']:
    for k in json.load(open(f)): done.add(int(k))
if os.path.exists(out):
    for line in open(out): done.add(json.loads(line)['k'])
rng=np.random.default_rng(1000+w); t0=time.time()
with open(out,'a') as fh:
    for k in range(len(T)):
        if k%W!=w or k in done: continue
        b=theta_star(T[k],rng,starts=1)
        if b is None or b[0]>58.0:
            b2=theta_star(T[k],rng,starts=4)
            if b2 is not None and (b is None or b2[0]>b[0]): b=b2
        fh.write(json.dumps(dict(k=k,theta=None if b is None else b[0],X=None if b is None else b[1].tolist()))+'\n'); fh.flush()
