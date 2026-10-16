"""Ascent + surgery census for N points on S^{d-1} (any N).
usage: python scripts/ascent_census.py N d runs rounds out.json

Status: untrusted (exploration or tooling; no part of the proof).
"""
import sys, os, time
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
import kisscrit as kc
N, d, runs, rounds, out = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), sys.argv[5]
C = kc.Census.load(out) if os.path.exists(out) else kc.Census(N, d); t0 = time.time()
kc.random_ascent(C, runs, np.random.default_rng(N * 100 + d)); C.dump(out)
print('ascent: %d types (%.0fs)' % (len(C.types), time.time() - t0), flush=True)
kc.surgery(C, max_rounds=rounds, near=8, verbose=True); C.dump(out)
print('surgery: %d types (%.0fs)' % (len(C.types), time.time() - t0))
for t in C.sorted():
    print('  r=%.9f k=%d coidx=%d nzero=%d ratt=%d labeled=%s hits=%d' % (t['r'], t['k'], t['coindex'], t['nzero'], t['rattlers'], t.get('labeled'), t['hits']))
