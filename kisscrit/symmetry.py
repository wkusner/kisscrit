"""Symmetry group of a point configuration on S^{d-1}: all (R in O(d), permutation p) with R x_i = x_p(i)."""
import numpy as np


def _frame(X, tol=1e-6):
    """Greedy choice of d indices whose points are linearly independent (or as many as possible)."""
    idx = []
    for i in range(len(X)):
        M = X[idx + [i]]
        if np.linalg.matrix_rank(M, tol) == len(idx) + 1:
            idx.append(i)
        if len(idx) == X.shape[1]:
            break
    return idx


def symmetries(X, tol=1e-6):
    """Return list of (det R, perm) for every orthogonal symmetry of the (spanning) configuration."""
    X = np.asarray(X, float); N, d = X.shape
    F = _frame(X, tol)
    if len(F) < d:
        raise ValueError('configuration does not span R^d; symmetry group is infinite')
    C = X @ X.T; out = []

    def extend(assign):
        m = len(assign)
        if m == d:
            Xi = X[F]; Yi = X[assign]
            R = np.linalg.solve(Xi, Yi).T          # R x_F = y
            if np.abs(R @ R.T - np.eye(d)).max() > 1e-5:
                return
            Y = X @ R.T
            D = np.linalg.norm(Y[:, None, :] - X[None, :, :], axis=2)
            p = np.argmin(D, axis=1)
            if np.all(D[np.arange(N), p] < tol) and len(set(p.tolist())) == N:
                out.append((int(round(np.linalg.det(R))), p.tolist()))
            return
        for j in range(N):
            if j in assign:
                continue
            if all(abs(C[F[m], F[q]] - C[j, assign[q]]) < tol for q in range(m)):
                extend(assign + [j])

    extend([])
    return out


def group_orders(X, tol=1e-6):
    """(|G|, |G ∩ SO(d)|).  A configuration spanning a hyperplane is handled in that hyperplane:
    each symmetry of the hyperplane extends to O(d) in two ways, exactly one of them proper."""
    X = np.asarray(X, float); N, d = X.shape
    rk = np.linalg.matrix_rank(X, 1e-8)
    if rk == d - 1:
        U, S_, Vt = np.linalg.svd(X)
        Y = X @ Vt[:d - 1].T
        g = len(symmetries(Y / np.linalg.norm(Y, axis=1, keepdims=True), tol))
        return 2 * g, g
    S = symmetries(X, tol)
    return len(S), sum(1 for det, _ in S if det == 1)
