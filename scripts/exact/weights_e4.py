"""EXACT-1 milestone E4: weight (self-stress) spaces of the configurations in data/exact/critical-configs-n12.json.

Status: untrusted (exploration or tooling; no part of the completeness proof).

Input is the stored exact data only (P, the isolating interval, the contacts and the Gram entries as strings);
nothing is re-derived numerically.  The Gram checks of scripts/exact/verify.py are re-run first.  Then, in K:

  balance system: for every point i and every frame point b (the frame B of the Gram checks),
        sum_{j ~ i} w_ij (G_bj - c G_bi) = 0,
    which is sum_j w_ij (u_j - c u_i) = 0 tested against a basis of R^3 (contacts have G_ij = c); it is
    linear in w with coefficients in K;
  K(C) = its solution space: dimension and an exact basis (reduced echelon form);
  the cone K(C) intersected with {w >= 0}: its extreme rays by the double description method (exact sign tests in K;
    adjacency by the combinatorial test: two rays are adjacent iff no third ray vanishes on every constraint
    on which both vanish), each ray normalized to sum 1;
  a strictly positive stress, if one exists: the normalized sum of the extreme rays (it is positive on an
    edge iff some extreme ray is), otherwise the edges on which every nonnegative stress vanishes;
  whether the uniform weight is a stress;
  the orbits of the contacts under the symmetry group (permutations preserving G exactly).

Writes a "weights" entry into each configuration of the JSON file.  Usage:
  python3 scripts/exact/weights_e4.py [labels...]   (default: every configuration in the file)
"""
import json
import os
import subprocess
import sys
import time

sys.set_int_max_str_digits(0)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from kfield import Field, KE, rref, rank, nullspace
import verify as V

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..')
OUT = os.path.join(ROOT, 'data', 'exact', 'critical-configs-n12.json')


def load(entry):
    """Field, contacts and Gram dict of a stored configuration (exact, from strings)."""
    F = Field([int(a) for a in entry['c']['minpoly']], lo=entry['c']['interval'][0], hi=entry['c']['interval'][1])
    A = [tuple(e) for e in entry['contacts']]
    G = {}
    for key, co in entry['gram']['entries'].items():
        i, j = map(int, key.split(','))
        G[(i, j)] = F(co)
    return F, A, G


def extreme_rays(F, basis):
    """Extreme rays (as alpha-coordinates in the basis) of {alpha : sum_i alpha_i basis_i >= 0 on every edge}."""
    m, k = len(basis), len(basis[0])
    a = [[basis[i][e] for i in range(m)] for e in range(k)]  # w_e = a_e . alpha
    dot = lambda u, v: sum((x * y for x, y in zip(u, v)), F.const(0))
    if m == 1:
        s = [x.sign() for x in basis[0]]
        if all(x >= 0 for x in s):
            return [[F.const(1)]]
        if all(x <= 0 for x in s):
            return [[F.const(-1)]]
        return []
    live = [e for e in range(k) if not all(x.is_zero() for x in a[e])]
    # initial simplicial cone from m independent constraints
    S = []
    for e in live:
        if rank([a[f] for f in S + [e]]) == len(S) + 1:
            S.append(e)
        if len(S) == m:
            break
    assert len(S) == m
    # rays: columns of the inverse of the matrix with rows a_e (e in S)
    M = [a[e] + [F.const(1) if r == i else F.const(0) for r in range(m)] for i, e in enumerate(S)]
    R, piv = rref(M)
    inv = [row[m:] for row in R]  # inverse of A_S
    rays = []
    for j in range(m):
        r = [inv[i][j] for i in range(m)]
        rays.append((r, frozenset(S[i] for i in range(m) if i != j)))
    done = list(S)
    for e in live:
        if e in S:
            continue
        s = [dot(a[e], r).sign() for r, _ in rays]
        val = [dot(a[e], r) for r, _ in rays]
        P = [i for i in range(len(rays)) if s[i] > 0]
        Z = [i for i in range(len(rays)) if s[i] == 0]
        N = [i for i in range(len(rays)) if s[i] < 0]
        new = [rays[i] for i in P] + [(rays[i][0], rays[i][1] | {e}) for i in Z]
        for p in P:
            for n in N:
                common = rays[p][1] & rays[n][1]
                if len(common) < m - 2:
                    continue
                if any(t != p and t != n and common <= rays[t][1] for t in range(len(rays))):
                    continue
                vp, vn = val[p], val[n]
                r = [vp * y - vn * x for x, y in zip(rays[p][0], rays[n][0])]
                new.append((r, common | {e}))
        rays = new
        done.append(e)
    return [r for r, _ in rays]


def edge_orbits(A, syms):
    """Orbits of the contacts under the permutations in syms."""
    idx = {tuple(sorted(e)): n for n, e in enumerate(A)}
    seen, orbits = set(), []
    for n, (i, j) in enumerate(A):
        if n in seen:
            continue
        orb = sorted({idx[tuple(sorted((p[i], p[j])))] for p, _ in syms})
        seen.update(orb)
        orbits.append(orb)
    return orbits


def weights(label, entry):
    """The E4 record for one configuration."""
    t0 = time.time()
    F, A, G = load(entry)
    checks, frame = V.gram_checks(F, A, G)
    B = frame[0]
    rows = V.balance_matrix(F, A, G, frame)
    basis = V.stress_space(F, A, G, frame)
    m = len(basis)
    rec = {'balance_system': ('for every point i and every frame point b in %s: sum over contacts ij of '
                              'w_ij (G_bj - c G_bi) = 0 (linear in w, coefficients in K = Q(c)); equivalent to '
                              'sum_j w_ij (u_j - c u_i) = 0' % list(B)),
           'balance_rank': len(A) - m, 'dim': m,
           'basis': [[F.coeffs(x) for x in b] for b in basis]}
    one = [F.const(1)] * len(A)
    rec['uniform_is_stress'] = V.is_stress(F, A, G, frame, one)
    raysA = extreme_rays(F, basis)
    rays = []
    for r in raysA:
        w = [sum((r[i] * basis[i][e] for i in range(m)), F.const(0)) for e in range(len(A))]
        s = sum(w[1:], w[0])
        w = [x / s for x in w]
        assert all(x.sign() >= 0 for x in w) and V.is_stress(F, A, G, frame, w)
        rays.append(w)
    rec['extreme_rays'] = [{'weights': [F.coeffs(x) for x in w],
                            'support': [A[e] for e in range(len(A)) if not w[e].is_zero()]} for w in rays]
    rec['n_extreme_rays'] = len(rays)
    if rays:
        tot = [sum((w[e] for w in rays), F.const(0)) for e in range(len(A))]
        s = sum(tot[1:], tot[0])
        tot = [x / s for x in tot]
        sg = [x.sign() for x in tot]
        rec['strictly_positive_stress_exists'] = all(x > 0 for x in sg)
        if rec['strictly_positive_stress_exists']:
            rec['strictly_positive_stress'] = {'weights': [F.coeffs(x) for x in tot],
                                               'how': 'normalized sum of the extreme rays (exact)'}
        else:
            rec['edges_zero_on_every_nonnegative_stress'] = [A[e] for e in range(len(A)) if sg[e] == 0]
    else:
        rec['strictly_positive_stress_exists'] = False
        rec['nonnegative_cone'] = '{0}'
    syms = V.symmetries(F, A, G, frame)
    orbs = edge_orbits(A, syms)
    rec['edge_orbits'] = [[A[e] for e in o] for o in orbs]
    rec['edge_orbit_sizes'] = sorted((len(o) for o in orbs), reverse=True)
    if m == 1 and rays:
        rec['the_stress'] = 'unique up to scale; the single extreme ray is the balanced stress (normalized sum 1)'
    rec['status'] = 'exactly verified (exact linear algebra and signs in K; double description with exact signs)'
    rec['time_s'] = round(time.time() - t0, 1)
    return rec


def main():
    """Compute and store the E4 records."""
    data = json.load(open(OUT))
    labels = sys.argv[1:] or list(data['configs'])
    h = subprocess.run(['git', '-C', ROOT, 'rev-parse', '--short', 'HEAD'], capture_output=True, text=True).stdout.strip()
    dirty = bool(subprocess.run(['git', '-C', ROOT, 'status', '--porcelain', 'scripts/exact'],
                                capture_output=True, text=True).stdout.strip())
    for lab in labels:
        rec = weights(lab, data['configs'][lab])
        rec['script'] = 'scripts/exact/weights_e4.py at %s%s' % (h, ' (DIRTY)' if dirty else '')
        data['configs'][lab]['weights'] = rec
        print('%-13s dim=%d rays=%d positive=%s uniform=%s orbits=%s %ss' % (
            lab, rec['dim'], rec['n_extreme_rays'], rec['strictly_positive_stress_exists'], rec['uniform_is_stress'],
            rec['edge_orbit_sizes'], rec['time_s']), flush=True)
    json.dump(data, open(OUT, 'w'), indent=1)


if __name__ == '__main__':
    main()
