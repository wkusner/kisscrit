"""Stresses, Clarke criticality, rigidity and the co-index of a critical configuration.

For the max-cosine function c(X) = max c_ij (equivalently the min-type function r), Clarke
criticality at X means: 0 lies in the convex hull of the tangential gradients of the active c_ij,
i.e. there is sigma >= 0, sum sigma = 1, with  sum_j sigma_ij P_i x_j = 0  at every vertex i
(a balanced nonnegative stress, P_i the projection to x_i^perp).
"""
import numpy as np
from scipy.optimize import linprog
from .geometry import tangent_bases, rotation_motions, contacts


def grad_rows(X, B, pairs):
    """Row e = tangential gradient of c_ij, (i,j) = pairs[e], in tangent coordinates."""
    N, dm1, _ = B.shape
    G = np.zeros((len(pairs), N * dm1))
    for e, (i, j) in enumerate(pairs):
        G[e, dm1 * i:dm1 * (i + 1)] = B[i] @ X[j]
        G[e, dm1 * j:dm1 * (j + 1)] = B[j] @ X[i]
    return G


def null_space(A, rtol=1e-9):
    if A.shape[0] == 0:
        return np.eye(A.shape[1])
    U, S, Vt = np.linalg.svd(A)
    rank = int((S > rtol * max(S[0], 1e-300)).sum()) if len(S) else 0
    return Vt[rank:].T


def rank(A, rtol=1e-9):
    if A.size == 0:
        return 0
    S = np.linalg.svd(A, compute_uv=False)
    return int((S > rtol * max(S[0], 1e-300)).sum())


def balanced_stress(X, pairs, B=None):
    """Max-min nonnegative balanced stress.

    Returns dict(sigma, tmin, support, resid, sdim):
      sigma  : stress with sum 1 maximizing min_e sigma_e (None if no nonnegative stress exists)
      tmin   : that minimum (>0 means a strictly positive stress exists)
      support: boolean mask of edges carrying weight in some nonnegative stress
      sdim   : dimension of the space of all (signed) self-stresses
    """
    if B is None:
        B = tangent_bases(X)
    G = grad_rows(X, B, pairs)
    k, n = G.shape
    sdim = k - rank(G)
    A_eq = np.zeros((n + 1, k + 1)); A_eq[:n, :k] = G.T; A_eq[n, :k] = 1
    b_eq = np.zeros(n + 1); b_eq[n] = 1
    cost = np.zeros(k + 1); cost[-1] = -1
    A_ub = np.zeros((k, k + 1)); A_ub[:, :k] = -np.eye(k); A_ub[:, -1] = 1
    res = linprog(cost, A_ub=A_ub, b_ub=np.zeros(k), A_eq=A_eq, b_eq=b_eq,
                  bounds=[(0, None)] * k + [(None, None)], method='highs')
    if res.status != 0:
        return dict(sigma=None, tmin=-np.inf, support=np.zeros(k, bool), resid=np.inf, sdim=sdim, G=G)
    sigma = res.x[:k]; tmin = res.x[-1]
    support = np.ones(k, bool) if tmin > 1e-10 else np.zeros(k, bool)
    if tmin <= 1e-10:
        for e in range(k):
            c2 = np.zeros(k); c2[e] = -1
            r2 = linprog(c2, A_eq=A_eq[:, :k], b_eq=b_eq, bounds=[(0, None)] * k, method='highs')
            support[e] = r2.status == 0 and -r2.fun > 1e-9
        # a stress with maximal support: average of per-edge maximizers is in the cone; use LP again
        if support.any():
            idx = np.where(support)[0]
            sub = balanced_stress_sub(G, idx)
            sigma = np.zeros(k); sigma[idx] = sub
    resid = float(np.linalg.norm(G.T @ sigma))
    return dict(sigma=sigma, tmin=float(tmin), support=support, resid=resid, sdim=sdim, G=G)


def balanced_stress_sub(G, idx):
    Gs = G[idx]; k, n = Gs.shape
    A_eq = np.zeros((n + 1, k + 1)); A_eq[:n, :k] = Gs.T; A_eq[n, :k] = 1
    b_eq = np.zeros(n + 1); b_eq[n] = 1
    cost = np.zeros(k + 1); cost[-1] = -1
    A_ub = np.zeros((k, k + 1)); A_ub[:, :k] = -np.eye(k); A_ub[:, -1] = 1
    res = linprog(cost, A_ub=A_ub, b_ub=np.zeros(k), A_eq=A_eq, b_eq=b_eq,
                  bounds=[(0, None)] * k + [(None, None)], method='highs')
    return res.x[:k]


def colspace(A, tol=1e-8):
    """Orthonormal basis of the column space of A."""
    if A.shape[1] == 0:
        return A
    U, S, _ = np.linalg.svd(A, full_matrices=False)
    return U[:, S > tol * max(1.0, S[0] if len(S) else 1.0)]


def quotient(K, R):
    """Orthonormal basis of the projection of span(K) onto the complement of span(R)."""
    Rb = colspace(R)
    return colspace(K - Rb @ (Rb.T @ K))


def lagrangian_Q(X, B, pairs, sigma):
    """Second variation of sum_e sigma_e c_e along great-circle motions, in tangent coordinates.

    r increases to second order along contact-preserving directions where Q < 0."""
    N, dm1, _ = B.shape
    Q = np.zeros((N * dm1, N * dm1))
    I = np.eye(dm1)
    for s, (i, j) in zip(sigma, pairs):
        if s == 0:
            continue
        cij = X[i] @ X[j]
        Bij = B[i] @ B[j].T
        si, sj = slice(dm1 * i, dm1 * (i + 1)), slice(dm1 * j, dm1 * (j + 1))
        Q[si, sj] += s * Bij
        Q[sj, si] += s * Bij.T
        Q[si, si] -= s * cij * I
        Q[sj, sj] -= s * cij * I
    return Q


def analyze(X, tol_contact=1e-9, tol_ev=1e-7):
    """Classify X as a critical configuration of the max-cosine function.

    Returns dict with: pairs (contacts), c, critical (positive-weight balance exists),
    sdim, tmin, support, rattlers, kernel (dimension of first-order contact-preserving motions on the
    stressed graph modulo rotations), ev (eigenvalues of Q on that kernel), coindex (# directions in
    which r increases to second order), nzero (degenerate directions), asc (ascending directions).
    """
    X = np.asarray(X, float); N, d = X.shape
    pairs, cm = contacts(X, tol_contact)
    B = tangent_bases(X)
    st = balanced_stress(X, pairs, B)
    out = dict(pairs=pairs, c=cm, k=len(pairs), sdim=st['sdim'], tmin=st['tmin'],
               critical=st['sigma'] is not None and st['resid'] < 1e-7)
    if not out['critical']:
        return out
    sup = st['support']; sp = [p for p, s in zip(pairs, sup) if s]; sig = st['sigma'][sup]
    deg = np.zeros(N, int)
    for i, j in sp:
        deg[i] += 1; deg[j] += 1
    out.update(support=sup.tolist(), nsupport=int(sup.sum()), rattlers=int((deg == 0).sum()), sigma=st['sigma'].tolist())
    Gs = grad_rows(X, B, sp)
    K = null_space(Gs)                       # motions keeping stressed contacts to first order
    Rm = rotation_motions(X, B)
    if Rm.size and K.shape[1]:
        K = quotient(K, Rm.T)
    Q = lagrangian_Q(X, B, sp, sig)
    if K.shape[1]:
        ev, V = np.linalg.eigh(K.T @ Q @ K)
        asc = K @ V[:, ev < -tol_ev]
    else:
        ev = np.array([]); asc = np.zeros((Q.shape[0], 0))
    out.update(kernel=int(K.shape[1]), ev=ev.tolist(), coindex=int((ev < -tol_ev).sum()),
               nzero=int((np.abs(ev) <= tol_ev).sum()), asc=asc, B=B)
    return out
