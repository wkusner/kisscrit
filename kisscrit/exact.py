"""High-precision refinement of critical configurations and minimal polynomials of c = cos(theta).

refine_kkt  : Newton (Gauss-Newton) on the KKT system at arbitrary precision (mpmath).
minpoly     : LLL search (python-flint) for an integer polynomial vanishing at c, with an acceptance
              test comparing the height to the precision used, which rejects the spurious high-degree
              relations LLL always returns at the edge of the precision budget.
"""
import numpy as np
import mpmath as mp

try:
    from flint import fmpz_mat, fmpz_poly
except ImportError:  # pragma: no cover
    fmpz_mat = fmpz_poly = None


def refine_kkt(X, A, s0, dps=200, verbose=False):
    """Refine (X, c, s, mu) of the KKT system of contact graph A to `dps` digits.

    Works in any dimension d. Gauge: sum_i x0_i ^ x_i = 0. Returns (X, c, s, residual) as mpmath objects."""
    X0 = np.array(X, float); N, d = X0.shape; k = len(A)
    C = X0 @ X0.T
    mu0 = np.zeros(N)
    for se, (i, j) in zip(s0, A):
        mu0[i] += se * C[i, j]; mu0[j] += se * C[i, j]
    mp.mp.dps = 30
    z = mp.matrix(list(X0.ravel()) + [float(np.mean([C[p] for p in A]))] + list(map(float, s0)) + list(mu0))
    n = len(z); ng = d * (d - 1) // 2

    def FJ(z):
        x = [[z[d * i + t] for t in range(d)] for i in range(N)]; c = z[N * d]
        s = [z[N * d + 1 + e] for e in range(k)]; mu = [z[N * d + 1 + k + i] for i in range(N)]
        m = N + k + N * d + 1 + ng
        f = mp.matrix(m, 1); J = mp.matrix(m, n); r = 0
        for i in range(N):
            f[r] = mp.fsum(t * t for t in x[i]) - 1
            for t in range(d): J[r, d * i + t] = 2 * x[i][t]
            r += 1
        for (i, j) in A:
            f[r] = mp.fsum(x[i][t] * x[j][t] for t in range(d)) - c
            for t in range(d): J[r, d * i + t] = x[j][t]; J[r, d * j + t] = x[i][t]
            J[r, N * d] = -1; r += 1
        Bal = [[-mu[i] * x[i][t] for t in range(d)] for i in range(N)]
        for i in range(N):
            for t in range(d):
                J[r + d * i + t, d * i + t] -= mu[i]; J[r + d * i + t, N * d + 1 + k + i] = -x[i][t]
        for e, (i, j) in enumerate(A):
            for t in range(d):
                Bal[i][t] += s[e] * x[j][t]; Bal[j][t] += s[e] * x[i][t]
                J[r + d * i + t, d * j + t] += s[e]; J[r + d * j + t, d * i + t] += s[e]
                J[r + d * i + t, N * d + 1 + e] = x[j][t]; J[r + d * j + t, N * d + 1 + e] = x[i][t]
        for i in range(N):
            for t in range(d): f[r + d * i + t] = Bal[i][t]
        r += N * d
        f[r] = mp.fsum(s) - 1
        for e in range(k): J[r, N * d + 1 + e] = 1
        r += 1
        for a in range(d):
            for b in range(a + 1, d):
                f[r] = mp.fsum(mp.mpf(X0[i, a]) * x[i][b] - mp.mpf(X0[i, b]) * x[i][a] for i in range(N))
                for i in range(N):
                    J[r, d * i + b] += X0[i, a]; J[r, d * i + a] -= X0[i, b]
                r += 1
        return f, J

    p = 30
    while True:
        mp.mp.dps = p
        for _ in range(12):
            f, J = FJ(z)
            dz = mp.lu_solve(J.T * J, -(J.T * f))
            z = z + dz
            if mp.norm(dz) < mp.mpf(10) ** (-p + 5):
                break
        if verbose:
            print('dps', p, 'resid', mp.nstr(mp.norm(FJ(z)[0]), 5), flush=True)
        if p >= dps:
            break
        p = min(2 * p, dps)
    x = [[z[d * i + t] for t in range(d)] for i in range(N)]
    s = [z[N * d + 1 + e] for e in range(k)]
    return x, z[N * d], s, mp.norm(FJ(z)[0])


def minpoly(c, maxdeg=60, digits=None, mindeg=1):
    """Smallest-degree integer polynomial vanishing at c (to `digits` digits), or None.

    A relation of degree d with height H is accepted only if log10(H) < digits/(d+1) - 3
    (a random relation found by LLL has height about 10^(digits/(d+1)))."""
    if fmpz_mat is None:
        raise ImportError('python-flint is required for minpoly')
    if digits is None:
        digits = mp.mp.dps - 10
    for d in range(mindeg, maxdeg + 1):
        co = minpoly_at(c, d, digits)
        if co is None:
            continue
        h = max(abs(t) for t in co)
        if mp.log10(h) < digits / (d + 1) - 3:
            return co
    return None


def minpoly_at(c, d, digits):
    scale = mp.mpf(10) ** digits
    pw = [mp.mpf(1)]
    for _ in range(d):
        pw.append(pw[-1] * c)
    M = []
    for k in range(d + 1):
        row = [0] * (d + 1); row[k] = 1
        row.append(int(mp.nint(scale * pw[k]))); M.append(row)
    L = fmpz_mat(M).lll()
    for rr in range(min(3, d + 1)):
        v = [int(L[rr, k]) for k in range(d + 1)]
        P = fmpz_poly(v)
        if P.degree() < 1:
            continue
        for f, _ in P.factor()[1]:
            co = [int(f[k]) for k in range(f.degree() + 1)]
            if abs(mp.polyval([mp.mpf(t) for t in reversed(co)], c)) < mp.mpf(10) ** (-(digits // 2)):
                return co
    return None


def poly_in_r(co_c):
    """Given integer coefficients (low->high) of a polynomial in c, return one in r (low->high) via
    c = 1 - 2 r^2/(1+r)^2, cleared of denominators, as sympy Poly factors."""
    import sympy as sp
    r, c = sp.symbols('r c')
    P = sum(int(a) * c ** i for i, a in enumerate(co_c))
    Q = sp.together(P.subs(c, 1 - 2 * r ** 2 / (1 + r) ** 2))
    num = sp.numer(Q)
    return sp.factor_list(sp.Poly(num, r))
