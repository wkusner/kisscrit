"""From the joins (results/joins_n12.json; double precision, assumed margin): component counts of BConf(12; theta(r)) for r > 1,
assuming the census of critical configurations above r = 1 is complete.

Status: untrusted (exploration or tooling; no part of the proof).
"""
import sys, os, json, math
import numpy as np
from sympy.combinatorics import Permutation, PermutationGroup
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
import kisscrit as kc
from kisscrit.joins import align
R = os.path.join(os.path.dirname(__file__), '..', 'results')
D = json.load(open(os.path.join(os.path.dirname(__file__), '..', 'data', 'n12_above1.json')))['types']
J = json.load(open(os.path.join(R, 'joins_n12.json')))
DOD = kc.normalize(np.array(D[0]['X'])); Ms, *_ = kc.solve_kkt(np.array(D[1]['X']), [tuple(p) for p in D[1]['pairs']]); Ms = kc.normalize(Ms)
rot = [Permutation(p.tolist()) for _, p in align(DOD, DOD)]          # DOD[p(i)] = R DOD[i]
msym = [p for _, p in align(Ms, Ms)]
print('|Rot(DOD)| =', len(rot), ' |Rot(M)| =', len(msym))
inv = lambda p: np.argsort(p)
out = {}
moves_M = {}; moves_D = {}
for key, v in sorted(J.items(), key=lambda kv: -float(kv[0])):
    ends = {v['plus']['end']: np.array(v['plus']['perm']), v['minus']['end']: np.array(v['minus']['perm'])}
    if 'M' in ends:   # DOD end in M-reference labels: q[m[i]] = d[i]
        m, d = ends['M'], ends['DOD']; q = np.empty(12, int); q[m] = d
        S = []
        for tau in msym:            # M symmetry relabels the reference: q o tau
            S.append(q[tau])
        moves_M[key] = S
        print(key, 'M-DOD saddle; DOD ends around a labeled M (mod rotation) distinct:',
              len({tuple(sorted([0])) for _ in S}), len(S))
    else:
        p1, p2 = np.array(v['plus']['perm']), np.array(v['minus']['perm'])
        g = Permutation(p2[inv(p1)].tolist())
        par = 'even' if g.is_even else 'odd'
        same = any((g * r).is_Identity for r in rot) or any(Permutation(p2[inv(p1)].tolist()) == r for r in rot)
        moves_D[key] = g
        print(key, 'DOD-DOD saddle; relabeling parity', par, '; same labeled DOD:', same)
        out[key] = dict(type='DOD-DOD', parity=par, same=bool(same))
keys = sorted(moves_M, key=float, reverse=True)
def group_from(Slists, extra):
    mv = []
    S = [s for L in Slists for s in L]
    for i in range(len(S)):
        for j in range(len(S)):
            if i != j: mv.append(Permutation(S[j][inv(S[i])].tolist()))
    return PermutationGroup(rot + mv + extra)
N = math.factorial(12)
H1 = group_from([moves_M[keys[0]]], [])
H2 = group_from([moves_M[k] for k in keys], [])
H3 = group_from([moves_M[k] for k in keys], list(moves_D.values()))
for name, H in [('after %s' % keys[0], H1), ('after %s' % keys[1], H2), ('after all index-1 saddles', H3)]:
    print(name, '|H| =', H.order(), ' DOD classes of labelings joined into', N // H.order(), 'component(s)')
out['groups'] = {'H_top': int(H1.order()), 'H_both_M_saddles': int(H2.order()), 'H_all': int(H3.order())}
json.dump(out, open(os.path.join(R, 'joins_groups_n12.json'), 'w'), indent=1)
