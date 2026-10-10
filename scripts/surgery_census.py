"""Extend a census by contact-graph surgery (and optional random ascent).
usage: python scripts/surgery_census.py census.json rounds [ascent_runs]"""
import sys, os, time
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
import kisscrit as kc
path, rounds = sys.argv[1], int(sys.argv[2]); runs = int(sys.argv[3]) if len(sys.argv) > 3 else 0
C = kc.Census.load(path); n0 = len(C.types); t0 = time.time()
if runs:
    kc.random_ascent(C, runs, np.random.default_rng(7))
kc.surgery(C, max_rounds=rounds, near=8, verbose=True)
C.dump(path)
print('types %d -> %d  (%.0fs)' % (n0, len(C.types), time.time() - t0))
