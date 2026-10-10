"""Signed counts vs chi(M_{0,N}) from results/census_n*.json."""
import sys, os, glob, warnings
warnings.filterwarnings('ignore')
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from kisscrit import Census
from kisscrit.search import euler_check, chi_M0N
D = {}
for p in glob.glob(os.path.join(os.path.dirname(__file__), '..', 'results', 'census_n*.json')):
    C = Census.load(p); D[C.N] = C
for N in sorted(D):
    if N < 4: continue
    tot, sk = euler_check({n: D[n] for n in D if n <= N}, N)
    print('N=%d  signed count %d  chi(M_0N) %d  skipped %s' % (N, tot, chi_M0N(N), sk))
    for t in D[N].sorted():
        print('     r=%.9f k=%d coidx=%d nzero=%d ratt=%d labeled=%s' % (t['r'], t['k'], t['coindex'], t['nzero'], t['rattlers'], t.get('labeled')))
