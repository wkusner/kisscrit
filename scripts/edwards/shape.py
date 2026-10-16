"""N points in the plane modulo similarity (Kendall shape space CP^{N-2}).
rho(x) = min_{i<j} |x_i - x_j| / sqrt(sum |x_i - xbar|^2).  Local maxima = jammed free clusters.

Status: untrusted (exploration or tooling; no part of the proof).
"""
import numpy as np, itertools, networkx as nx
from scipy.optimize import minimize, linprog
from networkx.algorithms.isomorphism import GraphMatcher

def normalize(x):
    x = x - x.mean(0); return x / np.sqrt((x**2).sum())

def pdist(x):
    i, j = np.triu_indices(len(x), 1); return np.linalg.norm(x[i] - x[j], axis=1), i, j

def rho(x):
    x = normalize(x); return pdist(x)[0].min()

def polish(x):
    """SLSQP: max t s.t. d_ij^2 >= t, centered, sum |x|^2 = 1."""
    n = len(x); x = normalize(x); d, I, J = pdist(x)
    z0 = np.concatenate([x.ravel(), [d.min()**2]])
    cons = [dict(type='ineq', fun=lambda z: ((z[:-1].reshape(n, 2)[I] - z[:-1].reshape(n, 2)[J])**2).sum(1) - z[-1]),
            dict(type='eq', fun=lambda z: np.concatenate([z[:-1].reshape(n, 2).sum(0), [(z[:-1]**2).sum() - 1]]))]
    r = minimize(lambda z: -z[-1], z0, constraints=cons, method='SLSQP', options=dict(maxiter=500, ftol=1e-14))
    return normalize(r.x[:-1].reshape(n, 2))

def ascend(x, steps=3000, h=0.02):
    """Morse protocol: steepest ascent of rho on the shape sphere (LP direction on the eps-active set)."""
    n = len(x); x = normalize(x)
    for k in range(steps):
        d, I, J = pdist(x); m = d.min(); eps = max(2 * h, 1e-7) * m
        act = np.nonzero(d <= m + eps)[0]
        G = np.zeros((len(act), 2 * n))
        for r, e in enumerate(act):
            u = (x[I[e]] - x[J[e]]) / d[e]
            G[r, 2 * I[e]:2 * I[e] + 2] = u; G[r, 2 * J[e]:2 * J[e] + 2] = -u
            G[r] -= d[e] * x.ravel()            # project out the radial (scale) part: grad of d/sqrt(I) at I=1
        # max s : G v >= s, v centered, |v|_inf <= 1
        A_ub = np.hstack([-G, np.ones((len(act), 1))]); b_ub = np.zeros(len(act))
        A_eq = np.zeros((2, 2 * n + 1)); A_eq[0, 0:2 * n:2] = 1; A_eq[1, 1:2 * n:2] = 1
        res = linprog(np.r_[np.zeros(2 * n), -1], A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=[0, 0],
                      bounds=[(-1, 1)] * (2 * n) + [(None, 1)], method='highs')
        s = res.x[-1]
        if s < 1e-9:
            h *= 0.5
            if h < 1e-6: break
            continue
        v = res.x[:-1].reshape(n, 2)
        x = normalize(x + h * m * v)
    return x

def invariant(x, dig=5):
    d = pdist(normalize(x))[0]; return tuple(np.round(np.sort(d), dig))

def aut_order(x, tol=1e-5):
    x = normalize(x); n = len(x); D = np.linalg.norm(x[:, None] - x[None], axis=2)
    G = nx.complete_graph(n)
    for a, b in G.edges: G[a][b]['d'] = D[a, b]
    gm = GraphMatcher(G, G, edge_match=lambda e1, e2: abs(e1['d'] - e2['d']) < tol)
    return sum(1 for _ in gm.isomorphisms_iter())

def contacts(x, tol=1e-6):
    d, I, J = pdist(normalize(x)); m = d.min(); return [(int(I[e]), int(J[e])) for e in np.nonzero(d < m * (1 + tol))[0]]


def is_max(y, trials=20, seed=0):
    """Numerical local-max test: re-ascend from small perturbations; a saddle lets rho climb."""
    rng = np.random.default_rng(seed); r0 = rho(y)
    for t in range(trials):
        z = polish(ascend(normalize(y) + 1e-4 * rng.normal(size=y.shape), steps=400, h=1e-3))
        if rho(z) > r0 + 1e-7: return False
    return True
