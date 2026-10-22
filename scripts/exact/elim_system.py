"""EXACT-1 E3: the KKT system of a critical configuration, reduced to its symmetric subspace, for msolve.

Status: untrusted (exploration or tooling; no part of the completeness proof).

By the principle of symmetric criticality, a configuration with a symmetry s (here a reflection or a
half-turn) is critical iff it is critical for the restriction of the problem to s-symmetric configurations,
and its (unique) stress is then s-invariant.  This script writes that restricted KKT system with rational
coefficients, for elimination of every variable except c (msolve -e).  Numerical data from
data/n12_above1.json are used only for combinatorial choices: which points are fixed by s, the pairing of
the others, the orbits of contacts, and the gauge point.  The polynomial system itself is exact.

Frame: a reflection is y -> -y (mirror plane y = 0); a half-turn is (x, y, z) -> (-x, -y, z).
Gauge: for a reflection, a point on the mirror plane is put at (0, 0, 1) (rotation about the y-axis);
for a half-turn, one representative point has y = 0 (rotation about the z-axis).
Unknowns: coordinates of orbit representatives, one weight per contact orbit (one of them set to 1),
and c (last).  Equations: unit length of each representative, <x_i, x_j> = c for one contact per orbit,
and the balance sum_j w_ij (x_j - c x_i) = 0 at each representative (contacts have <x_i, x_j> = c).

Usage: python3 scripts/exact/elim_system.py LABEL OUT.ms   (LABEL as in build_e2/build_e3)
"""
import json
import os
import sys
import itertools
import numpy as np
import networkx as nx

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import hiprec as H
try:
    import build_e3  # noqa: F401  (registers the E3 labels in hiprec.ABOVE1)
except Exception:
    pass


def symmetry(X, A):
    """The nontrivial isometry of the configuration: (perm, Q) with Q x_i = x_perm(i), from the contact-graph
    automorphisms (numerical Procrustes fit, residual < 1e-9)."""
    G = nx.Graph()
    G.add_nodes_from(range(len(X)))
    G.add_edges_from(A)
    for iso in nx.algorithms.isomorphism.GraphMatcher(G, G).isomorphisms_iter():
        p = [iso[i] for i in range(len(X))]
        if p == list(range(len(X))):
            continue
        Y = X[p]
        U, _, Vt = np.linalg.svd(Y.T @ X)
        Q = U @ Vt
        if np.abs(X @ Q.T - Y).max() < 1e-9:
            return p, Q
    return None, None


def build(label):
    """Variables, equations (msolve syntax) and the numerical c, for the configuration `label`."""
    A, X0, sig, t = H.load_above1(label)
    X = np.array(X0, float)
    p, Q = symmetry(X, A)
    assert p is not None, 'no nontrivial symmetry: the symmetric reduction does not apply'
    det = np.linalg.det(Q)
    # rotate so the mirror normal (reflection) or the axis (half-turn) is a coordinate axis
    w, V = np.linalg.eigh((Q + Q.T) / 2)
    if det < 0:  # reflection: eigenvalue -1 eigenvector -> y axis
        n = V[:, np.argmin(w)]
        e2 = n
        e3 = np.cross([1, 0, 0] if abs(n[0]) < 0.9 else [0, 1, 0], n)
        e3 /= np.linalg.norm(e3)
        e1 = np.cross(e2, e3)
        kind = 'reflection'
    else:  # half-turn: eigenvalue +1 eigenvector -> z axis
        n = V[:, np.argmax(w)]
        e3 = n
        e1 = np.cross([1, 0, 0] if abs(n[0]) < 0.9 else [0, 1, 0], n)
        e1 /= np.linalg.norm(e1)
        e2 = np.cross(e3, e1)
        kind = 'half-turn'
    R = np.array([e1, e2, e3])
    Xr = X @ R.T
    N = len(X)
    fixed = [i for i in range(N) if p[i] == i]
    reps = sorted(set(min(i, p[i]) for i in range(N)))
    # coordinates of every point as expressions in the representative's variables
    coords = {}
    var = []
    gauge = None
    if kind == 'reflection':
        gauge = fixed[0] if fixed else None
        assert gauge is not None, 'reflection without a point on the mirror plane (gauge not implemented)'
    else:
        gauge = [i for i in reps if p[i] != i][0]
    for i in reps:
        if kind == 'reflection':
            if p[i] == i:
                if i == gauge:
                    coords[i] = ('0', '0', '1')
                else:
                    var += ['x%d' % i, 'z%d' % i]
                    coords[i] = ('x%d' % i, '0', 'z%d' % i)
            else:
                var += ['x%d' % i, 'y%d' % i, 'z%d' % i]
                coords[i] = ('x%d' % i, 'y%d' % i, 'z%d' % i)
                coords[p[i]] = ('x%d' % i, '(-y%d)' % i, 'z%d' % i)
        else:
            if p[i] == i:  # a point on the axis: (0, 0, +-1)
                sgn = '1' if Xr[i, 2] > 0 else '(-1)'
                coords[i] = ('0', '0', sgn)
            elif i == gauge:
                var += ['x%d' % i, 'z%d' % i]
                coords[i] = ('x%d' % i, '0', 'z%d' % i)
                coords[p[i]] = ('(-x%d)' % i, '0', 'z%d' % i)
            else:
                var += ['x%d' % i, 'y%d' % i, 'z%d' % i]
                coords[i] = ('x%d' % i, 'y%d' % i, 'z%d' % i)
                coords[p[i]] = ('(-x%d)' % i, '(-y%d)' % i, 'z%d' % i)
    # contact orbits
    E = [tuple(sorted(e)) for e in A]
    orb = {}
    for e in E:
        f = tuple(sorted((p[e[0]], p[e[1]])))
        key = min(e, f)
        orb.setdefault(key, set()).update({e, f})
    orbits = sorted(orb)
    # weights: one per orbit; the orbit carrying the largest stored weight gets w = 1
    sigma = dict(zip(E, sig))
    big = max(orbits, key=lambda o: sigma[o])
    wname = {}
    for o in orbits:
        wname[o] = '1' if o == big else 'w%d_%d' % o
        if o != big:
            var.append('w%d_%d' % o)
    weight = {e: wname[o] for o in orbits for e in orb[o]}
    var.append('c')
    eqs = []

    def dot(a, b):
        return '+'.join('%s*%s' % (a[k], b[k]) for k in range(3))
    for i in reps:
        if i == gauge and kind == 'reflection':
            continue
        if kind == 'half-turn' and p[i] == i:
            continue
        eqs.append('%s-1' % dot(coords[i], coords[i]))
    for o in orbits:
        i, j = o
        eqs.append('%s-c' % dot(coords[i], coords[j]))
    for i in reps:
        nb = [(j, weight[tuple(sorted((i, j)))]) for j in range(N) if tuple(sorted((i, j))) in weight]
        for k in range(3):
            terms = ['%s*(%s-c*%s)' % (wt, coords[j][k], coords[i][k]) for j, wt in nb]
            expr = '+'.join(terms)
            # drop components that vanish identically by symmetry (all coordinates '0')
            if all(coords[j][k] == '0' for j, _ in nb) and coords[i][k] == '0':
                continue
            eqs.append(expr)
    return {'label': label, 'kind': kind, 'perm': p, 'fixed': fixed, 'reps': reps, 'gauge': gauge,
            'vars': var, 'eqs': eqs, 'c_numerical': float(t['c']), 'orbits': len(orbits)}


def write_msolve(S, path):
    """Write the system in msolve's input format (variables, characteristic 0, expanded polynomials with
    integer coefficients; '^' for powers)."""
    import sympy as sp
    syms = sp.symbols(S['vars'])
    loc = dict(zip(S['vars'], syms))
    polys = []
    for e in S['eqs']:
        q = sp.Poly(sp.expand(sp.sympify(e, locals=loc)), *syms)
        if q.is_zero:
            continue
        q = q.clear_denoms()[1]
        polys.append(str(q.as_expr()).replace('**', '^').replace(' ', ''))
    with open(path, 'w') as f:
        f.write(','.join(S['vars']) + '\n0\n')
        f.write(',\n'.join(polys) + '\n')
    return len(polys)


if __name__ == '__main__':
    S = build(sys.argv[1])
    S['n_written'] = write_msolve(S, sys.argv[2])
    print(json.dumps({k: S[k] for k in ('label', 'kind', 'fixed', 'reps', 'gauge', 'c_numerical', 'orbits')}))
    print('variables %d, equations %d' % (len(S['vars']), S['n_written']))
