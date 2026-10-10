"""Monotone ascent of r (descent of the max cosine) by LP directions on the delta-active set.

At each step: take the pairs within delta of the max cosine, find a tangent direction v with
|v|_inf <= 1 decreasing all of them at rate >= t (maximize t by LP), and step with backtracking so
that the max cosine strictly decreases.  If t ~ 0 the delta-active set carries a nonnegative balanced
stress (LP duality) and delta is shrunk; when delta is tiny the configuration is (numerically) Clarke
critical.  This flow is monotone, unlike sticky growth.
"""
import numpy as np
from scipy.optimize import linprog
from .geometry import normalize, tangent_bases, exp_move, r_of_c
from .stress import grad_rows


def ascend(X, maxit=5000, delta0=1e-3, target_c=None, record=False):
    X = normalize(X); N, d = X.shape; delta = delta0
    iu = np.triu_indices(N, 1); pairs_all = list(zip(*iu)); path = []
    it = 0
    for it in range(maxit):
        C = X @ X.T; cm = C[iu].max()
        if target_c is not None and cm <= target_c:
            break
        act = [p for p in pairs_all if C[p] >= cm - delta]
        B = tangent_bases(X); G = grad_rows(X, B, act); k, n = G.shape
        cost = np.zeros(n + 1); cost[-1] = -1
        res = linprog(cost, A_ub=np.hstack([G, np.ones((k, 1))]), b_ub=np.zeros(k),
                      bounds=[(-1, 1)] * n + [(None, 1)], method='highs')
        t = res.x[-1]; v = res.x[:-1]
        if t < 1e-9:
            if delta < 1e-12:
                break
            delta *= 0.1; continue
        h = min(0.05, delta / t); ok = False
        while h > 1e-14:
            Y = exp_move(X, B, v, h); cY = (Y @ Y.T)[iu].max()
            if cY < cm - 0.3 * h * t:
                ok = True; break
            h *= 0.5
        if not ok:
            delta *= 0.1
            if delta < 1e-13:
                break
            continue
        X = Y
        if record:
            path.append(cm)
        delta = min(delta0, max(delta * 2, 1e-6))
    cm = (X @ X.T)[iu].max()
    info = dict(c=float(cm), r=float(r_of_c(cm)), it=it)
    if record:
        info['path'] = path
    return X, info
