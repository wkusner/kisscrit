"""The square KKT system of a contact graph, and a Levenberg-Marquardt solver for it.

For a contact graph A on N points in S^{d-1}, critical configurations with contact graph A are the
solutions of

    |x_i|^2 = 1                                   (N equations)
    <x_i, x_j> = c                  (ij in A)     (|A| equations)
    sum_{j: ij in A} s_ij x_j - mu_i x_i = 0      (N d equations)
    sum s_ij = 1                                  (1 equation)

in the unknowns X (N d), c (1), s (|A|), mu (N), with s > 0 and every non-contact cosine < c.
The balance equations satisfy d(d-1)/2 torque identities (sum_i x_i ^ (balance_i) = 0), and the
solution set is invariant under O(d); after a gauge fixing the system is square:
(N d + 1 + |A| + N) unknowns minus d(d-1)/2 gauge = equations minus d(d-1)/2 identities.
We impose the linear gauge  sum_i x0_i ^ x_i = 0  (x0 = the starting configuration) and solve the
overdetermined-but-consistent system by Levenberg-Marquardt.
"""
import numpy as np
from .geometry import normalize, r_of_c, tangent_bases
from .stress import grad_rows


def _unpack(z, N, d, k):
    X = z[:N * d].reshape(N, d); c = z[N * d]
    s = z[N * d + 1:N * d + 1 + k]; mu = z[N * d + 1 + k:]
    return X, c, s, mu


def residual_and_jacobian(z, N, d, A, X0):
    k = len(A)
    X, c, s, mu = _unpack(z, N, d, k)
    ng = d * (d - 1) // 2
    m = N + k + N * d + 1 + ng
    n = len(z)
    F = np.zeros(m); J = np.zeros((m, n))
    ix = lambda i: slice(d * i, d * (i + 1))
    ic = N * d; isg = N * d + 1; imu = N * d + 1 + k
    r = 0
    for i in range(N):                              # sphere
        F[r] = X[i] @ X[i] - 1; J[r, ix(i)] = 2 * X[i]; r += 1
    for e, (i, j) in enumerate(A):                  # equal contacts
        F[r] = X[i] @ X[j] - c
        J[r, ix(i)] = X[j]; J[r, ix(j)] = X[i]; J[r, ic] = -1; r += 1
    b0 = r                                          # balance
    Bal = -mu[:, None] * X
    for i in range(N):
        J[b0 + d * i:b0 + d * (i + 1), ix(i)] -= mu[i] * np.eye(d)
        J[b0 + d * i:b0 + d * (i + 1), imu + i] = -X[i]
    for e, (i, j) in enumerate(A):
        Bal[i] += s[e] * X[j]; Bal[j] += s[e] * X[i]
        J[b0 + d * i:b0 + d * (i + 1), ix(j)] += s[e] * np.eye(d)
        J[b0 + d * j:b0 + d * (j + 1), ix(i)] += s[e] * np.eye(d)
        J[b0 + d * i:b0 + d * (i + 1), isg + e] = X[j]
        J[b0 + d * j:b0 + d * (j + 1), isg + e] = X[i]
    F[b0:b0 + N * d] = Bal.ravel(); r = b0 + N * d
    F[r] = s.sum() - 1; J[r, isg:isg + k] = 1; r += 1
    for a in range(d):                              # gauge: sum_i x0_i ^ x_i = 0
        for b in range(a + 1, d):
            F[r] = np.sum(X0[:, a] * X[:, b] - X0[:, b] * X[:, a])
            for i in range(N):
                J[r, d * i + b] += X0[i, a]; J[r, d * i + a] -= X0[i, b]
            r += 1
    return F, J


def initial_stress(X, A):
    B = tangent_bases(X); G = grad_rows(X, B, A); k = len(A)
    M = np.vstack([G.T, np.ones(k)]); b = np.zeros(M.shape[0]); b[-1] = 1
    s, *_ = np.linalg.lstsq(M, b, rcond=None)
    return s


def solve_kkt(X, A, s0=None, iters=100, tol=1e-13):
    """Solve the KKT system of contact graph A starting from X. Returns (X, c, s, residual_norm)."""
    X = normalize(X); N, d = X.shape; k = len(A)
    if s0 is None:
        s0 = initial_stress(X, A)
    C = X @ X.T
    mu = np.zeros(N)
    for se, (i, j) in zip(s0, A):
        mu[i] += se * C[i, j]; mu[j] += se * C[i, j]
    z = np.concatenate([X.ravel(), [np.mean([C[p] for p in A])], s0, mu])
    X0 = X.copy(); lam = 1e-8
    F, J = residual_and_jacobian(z, N, d, A, X0); nf = np.linalg.norm(F)
    for _ in range(iters):
        if nf < tol:
            break
        JTJ = J.T @ J; g = J.T @ F; D = np.diag(np.diag(JTJ)) + 1e-12 * np.eye(len(z))
        improved = False
        for _ in range(14):
            try:
                dz = np.linalg.solve(JTJ + lam * D, -g)
            except np.linalg.LinAlgError:
                lam *= 10; continue
            F2, J2 = residual_and_jacobian(z + dz, N, d, A, X0); n2 = np.linalg.norm(F2)
            if n2 < nf:
                z, F, J, nf = z + dz, F2, J2, n2; lam = max(lam / 10, 1e-15); improved = True; break
            lam *= 10
        if not improved:
            break
    Xs, c, s, mu = _unpack(z, N, d, k)
    return Xs, float(c), s, float(nf)


def check_solution(X, c, s, A, nf, res_tol=1e-10, gap_tol=1e-8):
    """Accept a KKT solution as a critical configuration with contact graph exactly A."""
    N = len(X)
    if nf > res_tol or not (-0.999 < c < 0.999) or s.min() < -1e-10:
        return None
    C = X @ X.T; Aset = set(map(tuple, A))
    others = [C[i, j] for i in range(N) for j in range(i + 1, N) if (i, j) not in Aset]
    gap = c - max(others) if others else 1.0
    if gap < gap_tol:
        return None
    return dict(c=float(c), r=float(r_of_c(c)), gap=float(gap), smin=float(s.min()))
