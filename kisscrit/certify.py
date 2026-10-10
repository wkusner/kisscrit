"""Interval certification of a nondegenerate critical configuration on S^2 (Krawczyk operator).

Given a contact graph A and a numerical solution of the KKT system, we certify:

  (1) existence and uniqueness, in an explicit box, of a solution z* = (X, c, sigma, mu) of the
      square system  S(z) = 0  obtained from the KKT system by adding the rotation gauge
      sum_i x0_i x x_i = 0 and dropping three balance rows (two at vertex i0, one at vertex i1);
  (2) that z* solves the full KKT system: by the torque identity sum_i x_i x B_i = 0, the dropped
      balance components satisfy a 3x3 linear system whose matrix is certified nonsingular;
  (3) sigma_e > 0 for all e, and c - <x_i,x_j> > 0 for all non-contacts (so the contact graph is exactly
      A and X is Clarke critical for the max-cosine function);
  (4) the co-index: the bordered Hessian
            [ H   J^T ]      H = Hessian in (X, c) of  L = sum sigma_e <x_i,x_j> - sum mu_i/2 (|x_i|^2-1)
            [ J   0   ]      J = Jacobian of (sphere, contact, gauge) constraints
      is enclosed over the box; Weyl's inequality applied to an enclosure of its spectrum shows it is
      nonsingular with exactly (rank J) + coindex negative eigenvalues.  Since sum sigma = 1, c agrees with L
      on the contact-preserving variety, so the coindex counts the directions in which r increases.

All enclosures use mpmath interval arithmetic (mp.iv); the midpoint linear algebra is non-rigorous but only
used to build the preconditioner and the box (as Krawczyk's theorem allows).
"""
import itertools
import mpmath as mp
from mpmath import iv


def _F_rows(z, N, A, X0, keep_rows):
    """KKT residual as a list (works for mpf and iv.mpf entries).  Row order:
    sphere(N), contacts(k), balance(3N), sum(1), gauge(3)."""
    k = len(A)
    x = [[z[3 * i + t] for t in range(3)] for i in range(N)]; c = z[3 * N]
    s = z[3 * N + 1:3 * N + 1 + k]; mu = z[3 * N + 1 + k:]
    F = []
    for i in range(N):
        F.append(x[i][0] * x[i][0] + x[i][1] * x[i][1] + x[i][2] * x[i][2] - 1)
    for (i, j) in A:
        F.append(x[i][0] * x[j][0] + x[i][1] * x[j][1] + x[i][2] * x[j][2] - c)
    B = [[-mu[i] * x[i][t] for t in range(3)] for i in range(N)]
    for e, (i, j) in enumerate(A):
        for t in range(3):
            B[i][t] = B[i][t] + s[e] * x[j][t]; B[j][t] = B[j][t] + s[e] * x[i][t]
    for i in range(N):
        F.extend(B[i])
    tot = s[0]
    for e in range(1, k):
        tot = tot + s[e]
    F.append(tot - 1)
    for a, b in ((0, 1), (0, 2), (1, 2)):
        g = 0
        for i in range(N):
            g = g + X0[i][a] * x[i][b] - X0[i][b] * x[i][a]
        F.append(g)
    return [F[r] for r in keep_rows]


def _J_rows(z, N, A, X0, keep_rows, zero):
    """Jacobian of the KKT residual (same row order), entries of the type of z."""
    k = len(A); n = 3 * N + 1 + k + N
    x = [[z[3 * i + t] for t in range(3)] for i in range(N)]
    s = z[3 * N + 1:3 * N + 1 + k]; mu = z[3 * N + 1 + k:]
    rows = []
    for i in range(N):
        r = [zero] * n
        for t in range(3): r[3 * i + t] = 2 * x[i][t]
        rows.append(r)
    for (i, j) in A:
        r = [zero] * n
        for t in range(3): r[3 * i + t] = x[j][t]; r[3 * j + t] = x[i][t]
        r[3 * N] = zero - 1
        rows.append(r)
    bal = [[[zero] * n for _ in range(3)] for _ in range(N)]
    for i in range(N):
        for t in range(3):
            bal[i][t][3 * i + t] = bal[i][t][3 * i + t] - mu[i]
            bal[i][t][3 * N + 1 + k + i] = zero - x[i][t]
    for e, (i, j) in enumerate(A):
        for t in range(3):
            bal[i][t][3 * j + t] = bal[i][t][3 * j + t] + s[e]
            bal[j][t][3 * i + t] = bal[j][t][3 * i + t] + s[e]
            bal[i][t][3 * N + 1 + e] = x[j][t]
            bal[j][t][3 * N + 1 + e] = x[i][t]
    for i in range(N):
        rows.extend(bal[i])
    r = [zero] * n
    for e in range(k): r[3 * N + 1 + e] = zero + 1
    rows.append(r)
    for a, b in ((0, 1), (0, 2), (1, 2)):
        r = [zero] * n
        for i in range(N):
            r[3 * i + b] = r[3 * i + b] + X0[i][a]; r[3 * i + a] = r[3 * i + a] - X0[i][b]
        rows.append(r)
    return [rows[q] for q in keep_rows]


def _choose_dropped(xm, N, A):
    """Pick (i0, p, q, i1, s): drop balance components p,q at i0 and s at i1, maximizing |det M|,
    M = [x_i0 x e_p, x_i0 x e_q, x_i1 x e_s]."""
    E = [mp.matrix([1, 0, 0]), mp.matrix([0, 1, 0]), mp.matrix([0, 0, 1])]
    def cross(u, v): return mp.matrix([u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]])
    best = None
    for i0 in range(min(N, 4)):
        for i1 in range(N):
            if i1 == i0: continue
            for p, q in ((0, 1), (0, 2), (1, 2)):
                for s_ in range(3):
                    xi0 = mp.matrix(xm[i0]); xi1 = mp.matrix(xm[i1])
                    cols = [cross(xi0, E[p]), cross(xi0, E[q]), cross(xi1, E[s_])]
                    M = mp.matrix(3, 3)
                    for a in range(3):
                        for b in range(3): M[a, b] = cols[b][a]
                    d = abs(mp.det(M))
                    if best is None or d > best[0]:
                        best = (d, i0, p, q, i1, s_)
    return best[1:]


def certify(Xmid, cmid, smid, mumid, A, dps=60, rad=None, verbose=False):
    """Certify a critical configuration. Inputs: high-precision midpoint (mp numbers or strings).
    Returns a dict of certified facts (interval enclosures as strings) or raises AssertionError."""
    mp.mp.dps = dps; iv.dps = dps
    N = len(Xmid); k = len(A); n = 3 * N + 1 + k + N
    zm = [mp.mpf(v) for row in Xmid for v in row] + [mp.mpf(cmid)] + [mp.mpf(v) for v in smid] + [mp.mpf(v) for v in mumid]
    xm = [[zm[3 * i + t] for t in range(3)] for i in range(N)]
    # gauge reference: the midpoint rounded to 30 digits, as exact rationals (so iv and mp agree)
    X0 = [[mp.mpf(mp.nstr(xm[i][t], 30)) for t in range(3)] for i in range(N)]
    i0, p, q, i1, s_ = _choose_dropped(xm, N, A)
    nrow_full = N + k + 3 * N + 1 + 3
    drop = {N + k + 3 * i0 + p, N + k + 3 * i0 + q, N + k + 3 * i1 + s_}
    keep = [r for r in range(nrow_full) if r not in drop]
    assert len(keep) == n
    # Newton polish at the midpoint (square system)
    for _ in range(8):
        Fm = mp.matrix(_F_rows(zm, N, A, X0, keep))
        Jm = mp.matrix(_J_rows(zm, N, A, X0, keep, mp.mpf(0)))
        dz = mp.lu_solve(Jm, -Fm)
        zm = [zm[a] + dz[a] for a in range(n)]
        if mp.norm(dz) < mp.mpf(10) ** (-dps + 8): break
    Jm = mp.matrix(_J_rows(zm, N, A, X0, keep, mp.mpf(0)))
    C = mp.inverse(Jm)
    if rad is None:
        rad = mp.mpf(10) ** (-(dps // 2))
    # Krawczyk
    Y = [iv.mpf([zm[a] - rad, zm[a] + rad]) for a in range(n)]
    zi = [iv.mpf(v) for v in zm]
    X0i = [[iv.mpf(v) for v in row] for row in X0]
    Fp = _F_rows(zi, N, A, X0i, keep)
    JY = _J_rows(Y, N, A, X0i, keep, iv.mpf(0))
    Ci = [[iv.mpf(C[a, b]) for b in range(n)] for a in range(n)]
    CF = [sum((Ci[a][b] * Fp[b] for b in range(n)), iv.mpf(0)) for a in range(n)]
    # M = I - C*J(Y)
    Mi = [[(iv.mpf(1) if a == b else iv.mpf(0)) - sum((Ci[a][l] * JY[l][b] for l in range(n)), iv.mpf(0))
           for b in range(n)] for a in range(n)]
    D = [Y[a] - zi[a] for a in range(n)]
    K = [zi[a] - CF[a] + sum((Mi[a][b] * D[b] for b in range(n)), iv.mpf(0)) for a in range(n)]
    inside = all(K[a].a > Y[a].a and K[a].b < Y[a].b for a in range(n))
    assert inside, 'Krawczyk inclusion failed'
    Z = K  # the unique zero lies in K
    # (2) dropped rows: torque system nonsingular on the box
    x = [[Z[3 * i + t] for t in range(3)] for i in range(N)]
    def cross_e(v, t):
        e = [iv.mpf(0)] * 3; e[t] = iv.mpf(1)
        return [v[1] * e[2] - v[2] * e[1], v[2] * e[0] - v[0] * e[2], v[0] * e[1] - v[1] * e[0]]
    cols = [cross_e(x[i0], p), cross_e(x[i0], q), cross_e(x[i1], s_)]
    det = (cols[0][0] * (cols[1][1] * cols[2][2] - cols[1][2] * cols[2][1])
           - cols[1][0] * (cols[0][1] * cols[2][2] - cols[0][2] * cols[2][1])
           + cols[2][0] * (cols[0][1] * cols[1][2] - cols[0][2] * cols[1][1]))
    assert not (det.a <= 0 <= det.b), 'torque determinant not certified nonzero'
    # (3) positivity of stress and gaps
    c = Z[3 * N]; sig = Z[3 * N + 1:3 * N + 1 + k]; mu = Z[3 * N + 1 + k:]
    assert all(se.a > 0 for se in sig), 'stress not certified positive'
    Aset = set(A); mingap = None
    for i in range(N):
        for j in range(i + 1, N):
            if (i, j) in Aset: continue
            g = c - (x[i][0] * x[j][0] + x[i][1] * x[j][1] + x[i][2] * x[j][2])
            assert g.a > 0, 'gap not certified positive for %s' % ((i, j),)
            mingap = g if mingap is None or g.a < mingap.a else mingap
    # (4) co-index from the bordered Hessian
    nv = 3 * N + 1
    H = [[iv.mpf(0)] * nv for _ in range(nv)]
    for e, (i, j) in enumerate(A):
        for t in range(3):
            H[3 * i + t][3 * j + t] = H[3 * i + t][3 * j + t] + sig[e]
            H[3 * j + t][3 * i + t] = H[3 * j + t][3 * i + t] + sig[e]
    for i in range(N):
        for t in range(3):
            H[3 * i + t][3 * i + t] = H[3 * i + t][3 * i + t] - mu[i]
    Jc = []
    for i in range(N):
        r = [iv.mpf(0)] * nv
        for t in range(3): r[3 * i + t] = 2 * x[i][t]
        Jc.append(r)
    for (i, j) in A:
        r = [iv.mpf(0)] * nv
        for t in range(3): r[3 * i + t] = x[j][t]; r[3 * j + t] = x[i][t]
        r[3 * N] = iv.mpf(-1); Jc.append(r)
    for a, b in ((0, 1), (0, 2), (1, 2)):
        r = [iv.mpf(0)] * nv
        for i in range(N):
            r[3 * i + b] = r[3 * i + b] + X0i[i][a]; r[3 * i + a] = r[3 * i + a] - X0i[i][b]
        Jc.append(r)
    m = len(Jc); nb = nv + m
    Bd = [[iv.mpf(0)] * nb for _ in range(nb)]
    for a in range(nv):
        for b in range(nv): Bd[a][b] = H[a][b]
    for r in range(m):
        for a in range(nv):
            Bd[nv + r][a] = Jc[r][a]; Bd[a][nv + r] = Jc[r][a]
    Bmid = mp.matrix(nb, nb); radF = mp.mpf(0)
    for a in range(nb):
        for b in range(nb):
            lo, hi = mp.mpf(Bd[a][b].a), mp.mpf(Bd[a][b].b)
            Bmid[a, b] = (lo + hi) / 2
            radF += (hi - lo) ** 2          # (2 * radius)^2: a safe overestimate of the radius^2
    ev = mp.eigsy(Bmid, eigvals_only=True)
    ev = sorted([ev[a] for a in range(nb)])
    # error of the midpoint eigen-solver at this precision is far below 10^-(dps-10); add it
    tol = mp.sqrt(radF) + mp.mpf(10) ** (-(dps - 10))
    gap0 = min(abs(e) for e in ev)
    assert gap0 > tol, 'bordered Hessian not certified nonsingular'
    nneg = sum(1 for e in ev if e < 0)
    coindex = nneg - m
    r_iv = (1 - c) / 2
    s_iv = iv.sqrt(r_iv); r_enc = s_iv / (1 - s_iv)
    out = dict(c=str(c), r=str(r_enc), coindex=coindex, constraint_rows=m,
               min_sigma=str(min(sig, key=lambda v: v.a)), min_gap=str(mingap),
               min_abs_bordered_eig=mp.nstr(gap0, 6), weyl_radius=mp.nstr(tol, 6),
               box_radius=mp.nstr(rad, 3), dropped=dict(i0=i0, comps=(p, q), i1=i1, comp=s_))
    if verbose:
        print(out)
    return out
