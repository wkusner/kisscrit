"""High-precision configurations and integer-relation guesses for EXACT-1.

Status: untrusted (exploration or tooling; no part of the completeness proof).

Nothing computed here is a proof.  This file produces *guesses*: a high-precision approximation of each
critical configuration (Gauss-Newton on its KKT system, written independently of kisscrit/exact.py), the
minimal polynomial of c = cos(theta) found by LLL, and each Gram entry written in the basis 1, c, ..., c^(d-1)
found by LLL.  Every guess is then checked by exact arithmetic in K = Q(c) (scripts/exact/verify.py); a wrong
guess makes that check fail, so the guessing step plays no role in correctness.

Starting points: data/n12_above1.json (double precision) for DOD, M, B, r1, A; closed forms for the
hexagonal antiprism band, FCC (cuboctahedron) and HCP (triangular orthobicupola).
"""
import json
import os
import mpmath as mp
from flint import fmpz_mat, fmpq, fmpq_poly

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..')

# index into data/n12_above1.json "types", and whether the KKT stress is unique (used by Newton)
ABOVE1 = {'DOD': (0, False), 'M': (1, False), 'B': (4, True), 'r1': (5, True), 'A': (9, True)}


def load_above1(label):
    """Contacts, double-precision coordinates and stored stress of one configuration of data/n12_above1.json."""
    D = json.load(open(os.path.join(ROOT, 'data', 'n12_above1.json')))['types']
    i, _ = ABOVE1[label]
    t = D[i]
    return [tuple(p) for p in t['pairs']], t['X'], t['sigma'], t


def closed_form(label, dps):
    """Coordinates (mpmath, at dps digits) of FCC, HCP and the hexagonal antiprism band, from closed forms."""
    with mp.workdps(dps + 20):
        pts = []
        if label in ('FCC', 'HCP'):
            # equator hexagon, a triangle above at azimuths 30, 150, 270 degrees and elevation arcsin(sqrt(2/3)),
            # and a triangle below: mirror image (HCP, orthobicupola) or turned by 60 degrees (FCC, cuboctahedron)
            for k in range(6):
                a = mp.pi * k / 3
                pts.append([mp.cos(a), mp.sin(a), mp.mpf(0)])
            ca, sa = 1 / mp.sqrt(3), mp.sqrt(mp.mpf(2) / 3)
            for k in range(3):
                a = mp.pi / 6 + 2 * mp.pi * k / 3
                pts.append([ca * mp.cos(a), ca * mp.sin(a), sa])
            off = 0 if label == 'HCP' else mp.pi / 3
            for k in range(3):
                a = mp.pi / 6 + 2 * mp.pi * k / 3 + off
                pts.append([ca * mp.cos(a), ca * mp.sin(a), -sa])
        elif label == 'antiprism':
            # two hexagons at heights +-h, turned by 30 degrees; all 24 edges equal iff h^2 = (2 sqrt3 - 3)/3
            h = mp.sqrt((2 * mp.sqrt(3) - 3) / 3)
            rho = mp.sqrt(1 - h * h)
            for k in range(6):
                a = mp.pi * k / 3
                pts.append([rho * mp.cos(a), rho * mp.sin(a), h])
            for k in range(6):
                a = mp.pi / 6 + mp.pi * k / 3
                pts.append([rho * mp.cos(a), rho * mp.sin(a), -h])
        else:
            raise KeyError(label)
        return [[+x for x in p] for p in pts]


def dot(a, b):
    """Euclidean inner product of two 3-vectors."""
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def contacts_of(X, tol):
    """Pairs whose inner product is within tol of the largest inner product (numerical contact detection)."""
    N = len(X)
    G = [[dot(X[i], X[j]) for j in range(N)] for i in range(N)]
    m = max(G[i][j] for i in range(N) for j in range(i + 1, N))
    return [(i, j) for i in range(N) for j in range(i + 1, N) if abs(G[i][j] - m) < tol], m


def newton(X0, A, dps, stress=True, w0=None, iters_per_level=8, verbose=False):
    """Gauss-Newton on the KKT system of contact graph A, refined to `dps` digits.

    Unknowns: X (12 x 3), c, and (if stress) the weights w (one per contact).  Equations: |x_i|^2 = 1,
    <x_i, x_j> = c on contacts, and (if stress) the balance sum_j w_ij (x_j - <x_i, x_j> x_i) = 0 at every i,
    sum w = 1; plus the rotation gauge sum_i x0_i x x_i = 0.  The system is overdetermined but consistent;
    at a nondegenerate solution (unique stress, or no stress needed because the graph is rigid modulo a common
    change of c) the Jacobian has full column rank and the iteration converges quadratically.
    Returns (X, c, w, residual)."""
    N, k = len(X0), len(A)
    mp.mp.dps = 30
    X = [[mp.mpf(v) for v in p] for p in X0]
    G0 = [dot(X[i], X[j]) for (i, j) in A]
    c = mp.fsum(G0) / k
    if stress:
        w = [mp.mpf(v) for v in w0] if w0 is not None else [mp.mpf(1) / k] * k
        s = mp.fsum(w)
        w = [v / s for v in w]
    else:
        w = []
    X0r = [[mp.mpf(v) for v in p] for p in X0]
    nv = 3 * N + 1 + (k if stress else 0)
    level = 30
    res = None
    while True:
        mp.mp.dps = level + 20
        X0r = [[mp.mpf(v) for v in p] for p in X0r]
        for _ in range(iters_per_level):
            rows, Fv = [], []

            def add(f, entries):
                Fv.append(f)
                rows.append(entries)
            for i in range(N):
                add(dot(X[i], X[i]) - 1, {3 * i + t: 2 * X[i][t] for t in range(3)})
            for (i, j) in A:
                e = {3 * i + t: X[j][t] for t in range(3)}
                for t in range(3):
                    e[3 * j + t] = e.get(3 * j + t, 0) + X[i][t]
                e[3 * N] = -1
                add(dot(X[i], X[j]) - c, e)
            if stress:
                Bv = [[mp.mpf(0)] * 3 for _ in range(N)]
                JB = [[{} for _ in range(3)] for _ in range(N)]
                for e_, (i, j) in enumerate(A):
                    g = dot(X[i], X[j])
                    for (a, b) in ((i, j), (j, i)):
                        for t in range(3):
                            Bv[a][t] += w[e_] * (X[b][t] - g * X[a][t])
                            r = JB[a][t]
                            r[3 * N + 1 + e_] = r.get(3 * N + 1 + e_, 0) + (X[b][t] - g * X[a][t])
                            for u in range(3):
                                # d/dx_b of w (x_b - (x_a.x_b) x_a)_t = w (delta_tu - x_a,t x_a,u)
                                r[3 * b + u] = r.get(3 * b + u, 0) + w[e_] * ((1 if t == u else 0) - X[a][t] * X[a][u])
                                # d/dx_a of -w (x_a.x_b) x_a,t = -w (x_b,u x_a,t + (x_a.x_b) delta_tu)
                                r[3 * a + u] = r.get(3 * a + u, 0) - w[e_] * (X[b][u] * X[a][t] + (g if t == u else 0))
                for a in range(N):
                    for t in range(3):
                        add(Bv[a][t], JB[a][t])
                add(mp.fsum(w) - 1, {3 * N + 1 + e_: 1 for e_ in range(k)})
            # gauge: sum_i x0_i x x_i = 0
            for t in range(3):
                t1, t2 = (t + 1) % 3, (t + 2) % 3
                f = mp.fsum(X0r[i][t1] * X[i][t2] - X0r[i][t2] * X[i][t1] for i in range(N))
                e = {}
                for i in range(N):
                    e[3 * i + t2] = X0r[i][t1]
                    e[3 * i + t1] = -X0r[i][t2]
                add(f, e)
            m = len(Fv)
            J = mp.matrix(m, nv)
            for r, e in enumerate(rows):
                for col, v in e.items():
                    J[r, col] = v
            Fm = mp.matrix(Fv)
            res = max(abs(v) for v in Fv)
            JT = J.T
            dz = mp.lu_solve(JT * J, -(JT * Fm))
            for i in range(N):
                for t in range(3):
                    X[i][t] += dz[3 * i + t]
            c += dz[3 * N]
            for e_ in range(len(w)):
                w[e_] += dz[3 * N + 1 + e_]
            step = max(abs(v) for v in dz)
            if verbose:
                print('  dps %d residual %s step %s' % (level, mp.nstr(res, 3), mp.nstr(step, 3)), flush=True)
            if step < mp.mpf(10) ** (-(level + 5)):
                break
        if level >= dps:
            break
        level = min(2 * level, dps)
    return X, c, w, res


def relation(vals, digits):
    """An integer relation sum a_k vals[k] ~ 0 found by LLL (FLINT) at `digits` digits, or None.

    The relation is returned only if it is small compared to the precision: its residual is below
    10^(-digits/2) relative to the scale and its height uses less than half of the digits per entry.
    (These are heuristics for accepting a *guess*; correctness is decided later, exactly.)"""
    n = len(vals)
    scale = mp.mpf(10) ** digits
    rows = []
    for k in range(n):
        r = [0] * n + [int(mp.nint(vals[k] * scale))]
        r[k] = 1
        rows.append(r)
    L = fmpz_mat(rows).lll()
    best = None
    for i in range(n):
        a = [int(L[i, j]) for j in range(n)]
        if all(x == 0 for x in a):
            continue
        h = max(abs(x) for x in a)
        resid = abs(mp.fsum(a[j] * vals[j] for j in range(n)))
        if resid < mp.mpf(10) ** (-(digits // 2)) and len(str(h)) * n < 0.8 * digits:
            if best is None or h < best[0]:
                best = (h, a)
    return best[1] if best else None


def minpoly_guess(c, dps, dmax=40):
    """The integer polynomial of least degree <= dmax vanishing at c (LLL guess; constant term first)."""
    with mp.workdps(dps):
        for d in range(1, dmax + 1):
            a = relation([c ** k for k in range(d + 1)], dps - 10)
            if a is not None and a[-1] != 0:
                if a[-1] < 0:
                    a = [-x for x in a]
                return a
    return None


def in_basis_guess(v, c, d, dps):
    """Rational coefficients q_0..q_{d-1} with v ~ sum q_k c^k (LLL guess), as fmpq_poly, or None."""
    with mp.workdps(dps):
        a = relation([c ** k for k in range(d)] + [v], dps - 10)
    if a is None or a[-1] == 0:
        return None
    b = a[-1]
    return fmpq_poly([fmpq(-x, b) for x in a[:-1]])
