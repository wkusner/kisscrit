"""Certified joins: for an index-1 critical configuration s on S^2, certify which local maxima the two sides
of its ascending direction lie with, in the superlevel set {rho > rho(s)}.

Status: untrusted (exploration or tooling; no part of the proof).

For each side we build a polygonal path
    y = y_+-  ->  ascent points  ->  a copy of the maximum (aligned)
and verify, rigorously up to an explicit floating-point error bound, that the max cosine along every
segment stays below c_s (so rho > rho(s) on the whole path).  Segments are the normalized straight segments
p_i(tau) = u_i(tau)/|u_i(tau)|, u_i(tau) = P_i + tau (Q_i - P_i).  On [a,b] each cosine <p_i,p_j> changes by
at most L (b - a) with L = (d_i + d_j) / (1 - max_k d_k), d_k = |Q_k - P_k|, since |p_k'| <= |u_k'| / |u_k|.
Point values are computed in double precision; the rounding error of a cosine of normalized 3-vectors is
below 1e-14, and we subtract 1e-13 from every margin.

The start points y_+- lie on the variety of configurations with all contacts of s equal, at distance t0 on
either side of s along the ascending direction w, with c(y) < c_s certified the same way.  That y_+ and y_-
lie in the two local components of {rho > rho(s)} near s is the min-type Morse lemma at a nondegenerate
critical point (Baryshnikov-Bubenik-Kahle), using the certified Hessian; it is not checked quantitatively here.
"""
import numpy as np
from scipy.optimize import linprog
from .geometry import normalize, tangent_bases, to_ambient, exp_move
from .stress import grad_rows, analyze
from .kkt import solve_kkt

FP_ERR = 1e-13


def max_cos(X):
    C = X @ X.T; np.fill_diagonal(C, -2.0); return C.max()


def segment_ok(P, Q, c_bound, depth=0, max_depth=40):
    """True if max cosine along the normalized segment P->Q is < c_bound (with error bound)."""
    stack = [(0.0, 1.0, 0)]
    D = np.linalg.norm(Q - P, axis=1); dmax = D.max()
    if dmax >= 0.5:
        return False
    L = (D[:, None] + D[None, :]) / (1 - dmax)
    np.fill_diagonal(L, 0)
    nchecks = 0
    while stack:
        a, b, dep = stack.pop()
        m = 0.5 * (a + b)
        Xm = normalize(P + m * (Q - P))
        C = Xm @ Xm.T; np.fill_diagonal(C, -2.0)
        bound = (C + L * (b - a) / 2).max() + FP_ERR
        nchecks += 1
        if bound < c_bound:
            continue
        if dep >= max_depth or C.max() + FP_ERR >= c_bound:
            return False
        stack += [(a, m, dep + 1), (m, b, dep + 1)]
    return nchecks


def push_on_variety(X, A, w, t0, iters=50):
    """Solve: |x_i| = 1, contacts in A equal, gauge sum x0_i x x_i = 0, <Y - X, w> = t0 (Newton)."""
    N = len(X); X0 = X.copy(); Y = normalize(X + t0 * w / np.dot(w.ravel(), w.ravel()))
    wf = w.ravel()
    for _ in range(iters):
        F = []; J = []
        for i in range(N):
            F.append(Y[i] @ Y[i] - 1); r = np.zeros(3 * N); r[3 * i:3 * i + 3] = 2 * Y[i]; J.append(r)
        (i0, j0) = A[0]
        for (i, j) in A[1:]:
            F.append(Y[i] @ Y[j] - Y[i0] @ Y[j0]); r = np.zeros(3 * N)
            r[3 * i:3 * i + 3] += Y[j]; r[3 * j:3 * j + 3] += Y[i]
            r[3 * i0:3 * i0 + 3] -= Y[j0]; r[3 * j0:3 * j0 + 3] -= Y[i0]; J.append(r)
        for a, b in ((0, 1), (0, 2), (1, 2)):
            F.append(np.sum(X0[:, a] * Y[:, b] - X0[:, b] * Y[:, a])); r = np.zeros(3 * N)
            for i in range(N):
                r[3 * i + b] += X0[i, a]; r[3 * i + a] -= X0[i, b]
            J.append(r)
        F.append((Y - X).ravel() @ wf - t0); J.append(wf.copy())
        F = np.array(F); J = np.array(J)
        dz, *_ = np.linalg.lstsq(J, -F, rcond=None)
        Y = Y + dz.reshape(N, 3)
        if np.linalg.norm(F) < 1e-15: break
    return Y


def ascent_path(X, maxit=20000, delta0=1e-3, hmax=0.02):
    """Monotone LP ascent recording every accepted configuration (steps <= hmax)."""
    X = normalize(X); N = len(X); iu = np.triu_indices(N, 1); pairs = list(zip(*iu))
    delta = delta0; path = [X.copy()]
    for it in range(maxit):
        C = X @ X.T; cm = C[iu].max()
        act = [p for p in pairs if C[p] >= cm - delta]
        B = tangent_bases(X); G = grad_rows(X, B, act); k, n = G.shape
        cost = np.zeros(n + 1); cost[-1] = -1
        res = linprog(cost, A_ub=np.hstack([G, np.ones((k, 1))]), b_ub=np.zeros(k),
                      bounds=[(-1, 1)] * n + [(None, 1)], method='highs')
        if res.x is None: break
        t = res.x[-1]; v = res.x[:-1]
        if t < 1e-9:
            if delta < 1e-12: break
            delta *= 0.1; continue
        h = min(hmax, delta / t); ok = False
        while h > 1e-14:
            Y = exp_move(X, B, v, h)
            if (Y @ Y.T)[iu].max() < cm - 0.3 * h * t:
                ok = True; break
            h *= 0.5
        if not ok:
            delta *= 0.1
            if delta < 1e-13: break
            continue
        X = Y; path.append(X.copy()); delta = min(delta0, max(2 * delta, 1e-6))
    return path


def align(Y, Ref, tol=1e-4):
    """Proper rotations R and permutations p with Y[i] ~ R Ref[p[i]]; returns list of (R, p)."""
    N = len(Y); out = []
    CY = Y @ Y.T; CR = Ref @ Ref.T
    i1 = int(np.argmin(np.abs(CY[0, 1:]))) + 1
    def frame(u, v):
        e2 = v - (v @ u) * u; e2 /= np.linalg.norm(e2); return np.column_stack([u, e2, np.cross(u, e2)])
    for a in range(N):
        for b in range(N):
            if a == b or abs(CR[a, b] - CY[0, i1]) > tol: continue
            R = frame(Y[0], Y[i1]) @ frame(Ref[a], Ref[b]).T
            if np.linalg.det(R) < 0: continue
            Z = Ref @ R.T; D = Y @ Z.T; p = D.argmax(1)
            if np.all(D[np.arange(N), p] > 1 - tol) and len(set(p.tolist())) == N:
                out.append((R, p))
    return out


def certify_side(Path, c_s, Ref_list):
    """Certify a path (list of configurations) and close it at an aligned copy of one of the references.
    Ref_list: list of (name, Ref configuration). Returns dict or raises AssertionError."""
    end = Path[-1]
    hit = None
    for name, Ref in Ref_list:
        al = align(end, Ref, tol=1e-3)
        if al:
            R, p = al[0]; hit = (name, R, p, Ref); break
    assert hit is not None, 'endpoint not identified'
    name, R, p, Ref = hit
    T = (Ref @ R.T)[p]                       # T[i] = R Ref[p[i]]
    pts = list(Path) + [T]
    tot = 0; worst = -1
    for P, Q in zip(pts[:-1], pts[1:]):
        r = segment_ok(P, Q, c_s)
        assert r, 'segment not certified'
        tot += r
        worst = max(worst, max_cos(P))
    return dict(end=name, perm=p.tolist(), segments=len(pts) - 1, checks=tot, max_cos_on_vertices=float(worst))
