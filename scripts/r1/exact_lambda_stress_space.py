"""Second-order values on the whole self-stress space at FCC and HCP, exact, by two routes (control item C5a).

Status: untrusted (an independent exact check written for C5a; not on the trusted list of
tests/test_code_standards.py, which only the control chat may change).
Implements: paper Appendix B, Lemma B.2(b) (quadratic part of phi_ij in the exponential chart) and
Proposition B.5 (the sentence "For every self-stress w (of either sign) f^T Q_w f = lambda * sum_ij w_ij with
lambda = 1/6 at FCC and lambda = 3/22 at HCP").

Question.  At a configuration C = (c_1..c_12) on S^2 with contact graph G (pairs at exactly 60 degrees), let f be
the unit nontrivial infinitesimal flex (the jitterbug) in V = (sum of tangent planes) minus the rotations, with
the Euclidean norm.  For an edge e = ij let q_e = -(d^2/dt^2) phi_e(exp(t f)) at t = 0, where
phi_e(v) = <exp_{c_i}(v_i), exp_{c_j}(v_j)>.  Then f^T Q_w f = sum_e w_e q_e for every weight w.  Does
q - lambda (1,...,1) lie in the column space of the rigidity matrix R (equivalently, is it orthogonal to every
self-stress)?

Exact arithmetic.  Every number lives in K = QQ.algebraic_field(sqrt 2, sqrt 3) of sympy, i.e. Q(theta) with
theta = sqrt 2 + sqrt 3, minimal polynomial x^4 - 10 x^2 + 1.  Elements are polynomials in theta of degree < 4
reduced modulo that polynomial (sympy ANP), so equality and zero tests are exact comparisons of four rational
coefficients.  Ranks, null spaces and solves use sympy DomainMatrix over K (fraction-free and field
Gaussian elimination, exact).  Signs (needed only for the nonnegative cone in part C) are decided exactly by
`sign`, which writes x = a + b sqrt2 + c sqrt3 + d sqrt6 with rational a..d and compares squares.  No
sp.simplify is used as a decision procedure.  Floating point appears only in `float_crosschecks` (part D,
mpmath at 50 digits), which is labeled as a cross-check, not proof.

Two routes to q.
  Route 1 (the formula of Lemma B.2(b)): q_e = f^T ((cD_e - W_e) (x) I_3) f, built as a 36 x 36 matrix.
  Route 2 (independent of that formula): symbolic second t-derivative at 0 of the closed form
  phi_e(exp(tF)) = cos(t a_i) cos(t a_j) <c_i,c_j> + cos(t a_i) sin(t a_j)/a_j <c_i,F_j>
                  + sin(t a_i)/a_i cos(t a_j) <F_i,c_j> + sin(t a_i) sin(t a_j)/(a_i a_j) <F_i,F_j>,
  a_i = |F_i|, done by sympy differentiation (and, separately, by series) on symbols; the result is checked to
  depend on a_i only through a_i^2, then evaluated exactly in K and divided by |F|^2.

Run: python3 scripts/r1/exact_lambda_stress_space.py   (prints parts A-D for FCC and HCP).
"""
import itertools
import random
import sys
from fractions import Fraction

import sympy as sp
from sympy.polys.matrices import DomainMatrix

S2, S3, S6 = sp.sqrt(2), sp.sqrt(3), sp.sqrt(6)
K = sp.QQ.algebraic_field(S2, S3)
ZERO, ONE = K.zero, K.one


# ---------------------------------------------------------------- exact field helpers

def el(x):
    """The element of K equal to the sympy number x (exact: sympy's to_number_field via minimal polynomials)."""
    return K.from_sympy(sp.sympify(x))


def coords(x):
    """Rationals (a, b, c, d) with x = a + b sqrt2 + c sqrt3 + d sqrt6.  Exact: with theta = sqrt2 + sqrt3,
    theta^2 = 5 + 2 sqrt6 and theta^3 = 11 sqrt2 + 9 sqrt3; the identity is re-verified by `_check_coords`."""
    c = [Fraction(int(r.numerator), int(r.denominator)) for r in x.to_list()]
    c = [Fraction(0)] * (4 - len(c)) + c          # highest degree first: c3 theta^3 + c2 theta^2 + c1 theta + c0
    c3, c2, c1, c0 = c
    return (c0 + 5 * c2, c1 + 11 * c3, c1 + 9 * c3, 2 * c2)


def _from_coords(a, b, c, d):
    """The element a + b sqrt2 + c sqrt3 + d sqrt6 of K (exact)."""
    return el(sp.Rational(a) + sp.Rational(b) * S2 + sp.Rational(c) * S3 + sp.Rational(d) * S6)


def _check_coords(x):
    """Exact check that `coords` is right for x (field equality)."""
    return _from_coords(*coords(x)) == x


def _sign_q2(p, q):
    """Exact sign of p + q sqrt2 for rationals p, q (compare p^2 with 2 q^2)."""
    if q == 0:
        return (p > 0) - (p < 0)
    if p == 0:
        return (q > 0) - (q < 0)
    if (p > 0) == (q > 0):
        return 1 if p > 0 else -1
    s = (p * p > 2 * q * q) - (p * p < 2 * q * q)      # sign of |p| - |q| sqrt2
    return s if p > 0 else -s


def sign(x):
    """Exact sign of x in K.  Write x = P + Q sqrt3 with P = a + b sqrt2, Q = c + d sqrt2; signs of P and Q are
    exact by `_sign_q2`; if they differ, compare P^2 with 3 Q^2, an element of Q(sqrt2), again by `_sign_q2`."""
    a, b, c, d = coords(x)
    sP, sQ = _sign_q2(a, b), _sign_q2(c, d)
    if sQ == 0:
        return sP
    if sP == 0:
        return sQ
    if sP == sQ:
        return sP
    # P^2 - 3 Q^2 = (a^2 + 2 b^2 - 3 c^2 - 6 d^2) + (2ab - 6cd) sqrt2
    s = _sign_q2(a * a + 2 * b * b - 3 * c * c - 6 * d * d, 2 * a * b - 6 * c * d)   # sign of |P| - |Q| sqrt3
    return s * sP


def show(x):
    """Readable sympy form of x in K (display only; never used for a decision)."""
    a, b, c, d = coords(x)
    return str(sp.Rational(a.numerator, a.denominator) + sp.Rational(b.numerator, b.denominator) * S2
               + sp.Rational(c.numerator, c.denominator) * S3 + sp.Rational(d.numerator, d.denominator) * S6)


def is_rational(x):
    """Exact: x lies in Q (its sqrt2, sqrt3, sqrt6 coordinates are 0)."""
    return coords(x)[1:] == (0, 0, 0)


def frac(x):
    """x as a Fraction; requires is_rational(x)."""
    assert is_rational(x)
    return coords(x)[0]


# ---------------------------------------------------------------- vectors and matrices over K

def dot(u, v):
    """Euclidean inner product of two vectors over K (exact)."""
    s = ZERO
    for a, b in zip(u, v):
        s += a * b
    return s


def cross(u, v):
    """Cross product in R^3 over K (exact)."""
    return [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]]


def dm(rows, nrows=None, ncols=None):
    """A DomainMatrix over K from a list of rows of K elements."""
    nrows = len(rows) if nrows is None else nrows
    ncols = (len(rows[0]) if rows else 0) if ncols is None else ncols
    return DomainMatrix([list(r) for r in rows], (nrows, ncols), K)


def rows_of(M):
    """The rows of a DomainMatrix over K as lists of K elements."""
    return [list(r) for r in M.to_list()]


def col(v):
    """Column DomainMatrix from a list of K elements."""
    return dm([[x] for x in v], len(v), 1)


def flat(M):
    """Entries of a one-column or one-row DomainMatrix as a list."""
    R = rows_of(M)
    return [x for r in R for x in r]


# ---------------------------------------------------------------- the two configurations

def hcp_points():
    """The triangular orthobicupola (HCP kissing configuration), exactly in K: the great hexagon
    (cos k pi/3, sin k pi/3, 0) and two aligned triangles at height +-sqrt6/3 and azimuths 30, 150, 270 degrees
    (radius 1/sqrt3, so cos of the latitude is 1/sqrt3: latitude arctan sqrt2).  Points 0-5 hexagon, 6-8 top,
    9-11 bottom, in the order of scripts/r1/exact_lambda.py:hcp()."""
    h = sp.Rational(1, 2)
    hexa = [(1, 0, 0), (h, S3 / 2, 0), (-h, S3 / 2, 0), (-1, 0, 0), (-h, -S3 / 2, 0), (h, -S3 / 2, 0)]
    xy = [(h, S3 / 6), (-h, S3 / 6), (0, -S3 / 3)]
    tri = lambda z: [(x, y, z) for (x, y) in xy]
    return [[el(t) for t in p] for p in hexa + tri(S6 / 3) + tri(-S6 / 3)]


def fcc_points():
    """The cuboctahedron (FCC kissing configuration): the 12 vectors (+-1, +-1, 0)/sqrt2 and permutations, in the
    order of scripts/r1/exact_lambda.py:fcc() (itertools.product over {-1,0,1}^3)."""
    return [[el(sp.Integer(x) * S2 / 2) for x in p] for p in itertools.product([-1, 0, 1], repeat=3)
            if sum(abs(x) for x in p) == 2]


# ---------------------------------------------------------------- part A: set-up

def check_unit(U):
    """Exact: every point has |c_i|^2 = 1 (field equality)."""
    return all(dot(u, u) == ONE for u in U)


def contact_edges(U):
    """Pairs ij (i < j) with <c_i, c_j> = 1/2 exactly (field equality).  Also returns whether every other pair
    has <c_i, c_j> < 1/2 (exact sign), i.e. the configuration is a packing at 60 degrees."""
    half = el(sp.Rational(1, 2))
    E, packing = [], True
    for i, j in itertools.combinations(range(len(U)), 2):
        g = dot(U[i], U[j])
        if g == half:
            E.append((i, j))
        elif sign(g - half) > 0:
            packing = False
    return E, packing


def tangent_basis(c):
    """Two vectors a = c x e_k (first k with c x e_k != 0) and b = c x a spanning the tangent plane at c, with
    entries in K.  Claim (checked in `setup`): <a,c> = <b,c> = 0 and a, b independent.  Not orthonormal; the
    Euclidean metric is used through ambient coordinates, never through these bases."""
    for k in range(3):
        e = [ZERO] * 3; e[k] = ONE
        a = cross(c, e)
        if any(x != ZERO for x in a):
            return [a, cross(c, a)]
    raise ValueError('zero vector')


def rigidity(U, E, T):
    """R: one row per edge ij, columns 2i, 2i+1 equal to <T_i[k], c_j> and 2j, 2j+1 equal to <T_j[k], c_i>.  For a
    tangent motion x (F_i = sum_k x_{2i+k} T_i[k]) row ij of R x is d/dt <c_i + tF_i, c_j + tF_j> at 0."""
    R = [[ZERO] * (2 * len(U)) for _ in E]
    for r, (i, j) in enumerate(E):
        for k in range(2):
            R[r][2 * i + k] = dot(T[i][k], U[j])
            R[r][2 * j + k] = dot(T[j][k], U[i])
    return dm(R, len(E), 2 * len(U))


def ambient(T, x):
    """The ambient vectors F_i = x_{2i} T_i[0] + x_{2i+1} T_i[1] (a list of 12 vectors in K^3)."""
    return [[x[2 * i] * T[i][0][m] + x[2 * i + 1] * T[i][1][m] for m in range(3)] for i in range(len(T))]


def rotation_fields(U):
    """The three infinitesimal rotations rho_a with rho_a,i = e_a x c_i (ambient)."""
    out = []
    for a in range(3):
        e = [ZERO] * 3; e[a] = ONE
        out.append([cross(e, u) for u in U])
    return out


def setup(U):
    """Part A.  Exact set-up: contact edges, tangent bases, R (24 x 24), rank, self-stresses (left null space of R,
    an rref basis), null space of R, the flex x in tangent coordinates (null space of R stacked with the three
    rows 'ambient inner product with rho_a = 0'), its ambient form F and |F|^2.  The unit flex is f = F/|F|."""
    N = len(U)
    E, packing = contact_edges(U)
    T = [tangent_basis(u) for u in U]
    tang_ok = all(dot(T[i][k], U[i]) == ZERO for i in range(N) for k in range(2))
    R = rigidity(U, E, T)
    rank = R.rank()
    stresses = [r for r in rows_of(R.transpose().nullspace().rref()[0])]
    nullR = rows_of(R.nullspace())
    rho = rotation_fields(U)
    # each rotation is an infinitesimal motion of every edge: <rho_i, c_j> + <c_i, rho_j> = 0, and tangent
    rot_ok = all(dot(r[i], U[j]) + dot(U[i], r[j]) == ZERO for r in rho for (i, j) in E) and \
        all(dot(r[i], U[i]) == ZERO for r in rho for i in range(N))
    # rows rho_a . (Tb x): coefficient of x_{2i+k} is <rho_a,i, T_i[k]>
    rot_rows = [[dot(r[i], T[i][k]) for i in range(N) for k in range(2)] for r in rho]
    A = dm(rows_of(R) + rot_rows, len(E) + 3, 2 * N)
    flexes = rows_of(A.nullspace())
    out = dict(U=U, E=E, packing=packing, T=T, tangent_ok=tang_ok, R=R, rank=rank, stresses=stresses,
               null_dim=len(nullR), rot_ok=rot_ok, rho=rho, flex_dim=len(flexes))
    if len(flexes) == 1:
        x = flexes[0]
        F = ambient(T, x)
        out.update(x=x, F=F, normF2=sum((dot(v, v) for v in F), ZERO),
                   F_tangent=all(dot(F[i], U[i]) == ZERO for i in range(N)),
                   F_perp_rot=all(sum((dot(r[i], F[i]) for i in range(N)), ZERO) == ZERO for r in rho),
                   F_is_flex=all(dot(F[i], U[j]) + dot(U[i], F[j]) == ZERO for (i, j) in E))
    return out


# ---------------------------------------------------------------- part B: per-edge values by two routes

def stress_laplacian_form(U, E, w, F):
    """Route 1.  F^T ((cD - W) (x) I_3) F in ambient coordinates, with W the weighted adjacency matrix of w on E and
    cD = diag(sum_j w_ij <c_i, c_j>), as a 36 x 36 DomainMatrix (the formula of Lemma B.2(b), sign reversed)."""
    N = len(U)
    W = [[ZERO] * N for _ in range(N)]
    for wk, (i, j) in zip(w, E):
        W[i][j] = W[j][i] = wk
    cD = [sum((W[i][j] * dot(U[i], U[j]) for j in range(N)), ZERO) for i in range(N)]
    L = [[ZERO] * (3 * N) for _ in range(3 * N)]
    for i in range(N):
        for j in range(N):
            v = (cD[i] if i == j else ZERO) - W[i][j]
            if v != ZERO:
                for m in range(3):
                    L[3 * i + m][3 * j + m] = v
    Fv = col([F[i][m] for i in range(N) for m in range(3)])
    return flat(Fv.transpose() * dm(L, 3 * N, 3 * N) * Fv)[0]


def route1(S):
    """Per-edge values q_e = f^T Q_{e} f, f = F/|F|, by route 1 (exact in K)."""
    U, E, F, n2 = S['U'], S['E'], S['F'], S['normF2']
    out = []
    for k in range(len(E)):
        w = [ZERO] * len(E); w[k] = ONE
        out.append(stress_laplacian_form(U, E, w, F) / n2)
    return out


_t = sp.Symbol('t')
_ai, _aj, _Ai, _Aj = sp.symbols('a_i a_j A_i A_j', positive=True)
_gcc, _gcF, _gFc, _gFF = sp.symbols('g_cc g_cF g_Fc g_FF')


def second_derivative_closed_form():
    """Route 2, symbolic part.  With u_i(t) = cos(t a_i) c_i + sin(t a_i)/a_i F_i (the exponential chart along tF),
    phi(t) = <u_i(t), u_j(t)> is expanded in the symbols g_cc = <c_i,c_j>, g_cF = <c_i,F_j>, g_Fc = <F_i,c_j>,
    g_FF = <F_i,F_j>.  Returns the polynomial P(A_i, A_j, g..) with P = -(d^2/dt^2) phi at t = 0 and
    A = a^2.  Claims, each decided exactly by sympy on polynomials with rational coefficients: the second
    derivative computed by diff equals 2 x (the t^2 coefficient of the series), and after a = sqrt(A) it is a
    polynomial in A_i, A_j, g.. over Q."""
    ci, cj = sp.cos(_t * _ai), sp.cos(_t * _aj)
    si, sj = sp.sin(_t * _ai) / _ai, sp.sin(_t * _aj) / _aj
    phi = ci * cj * _gcc + ci * sj * _gcF + si * cj * _gFc + si * sj * _gFF
    d2 = sp.expand(sp.diff(phi, _t, 2).subs(_t, 0))
    ser = sp.expand(sp.series(phi, _t, 0, 3).removeO().coeff(_t, 2) * 2)
    assert sp.expand(d2 - ser) == 0                 # polynomial identity in a, g (expanded form is canonical)
    P = sp.expand(-d2.subs({_ai: sp.sqrt(_Ai), _aj: sp.sqrt(_Aj)}))
    poly = sp.Poly(P, _Ai, _Aj, _gcc, _gcF, _gFc, _gFF, domain=sp.QQ)   # raises if not a polynomial over Q
    return poly


def eval_poly(poly, vals):
    """Evaluate a polynomial over Q at a tuple of K elements, exactly (sum of rational multiples of monomials)."""
    s = ZERO
    for mon, c in poly.terms():
        term = K.convert(sp.QQ(int(c.numerator), int(c.denominator)))
        for v, e in zip(vals, mon):
            for _ in range(e):
                term = term * v
        s += term
    return s


def route2(S, poly=None):
    """Per-edge values q_e = -(d^2/dt^2) phi_e(exp(t f)) at 0, f = F/|F|, by route 2: the closed-form polynomial of
    `second_derivative_closed_form` evaluated exactly in K at the data of F, divided by |F|^2 (because
    phi_e(exp(t f)) = phi_e(exp((t/|F|) F)))."""
    poly = second_derivative_closed_form() if poly is None else poly
    U, E, F, n2 = S['U'], S['E'], S['F'], S['normF2']
    out = []
    for (i, j) in E:
        vals = (dot(F[i], F[i]), dot(F[j], F[j]), dot(U[i], U[j]), dot(U[i], F[j]), dot(F[i], U[j]), dot(F[i], F[j]))
        out.append(eval_poly(poly, vals) / n2)
    return out


# ---------------------------------------------------------------- part C: the whole stress space

def is_stress(R, w):
    """Exact: w^T R = 0 (w is a self-stress)."""
    return all(x == ZERO for x in flat(dm([w], 1, len(w)) * R))


def solve_witness(R, b):
    """A vector y with R y = b, by exact Gaussian elimination over K (rref of [R | b]); None if inconsistent.  The
    returned y is re-verified by `verify_witness`."""
    m, n = R.shape
    Ab = dm([r + [bb] for r, bb in zip(rows_of(R), b)], m, n + 1)
    rref, piv = Ab.rref()
    rr = rows_of(rref)
    if n in piv:
        return None
    y = [ZERO] * n
    for r, p in enumerate(piv):
        y[p] = rr[r][n]
    return y


def verify_witness(R, y, b):
    """Exact: R y = b, entry by entry (field equality)."""
    return flat(R * col(y)) == list(b)


def stress_space_report(S, q, lam):
    """Part C.  For each basis stress s: s is a stress, sum_e s_e q_e, sum_e s_e, and whether
    sum s q = lam sum s.  Whether the uniform stress is a stress.  The witness y with R y = q - lam 1 (exists iff
    q - lam 1 is orthogonal to every self-stress; covers the whole stress space at once)."""
    R, B = S['R'], S['stresses']
    rows = []
    for s in B:
        sq = sum((a * b for a, b in zip(s, q)), ZERO)
        ss = sum(s, ZERO)
        rows.append(dict(s=s, is_stress=is_stress(R, s), sq=sq, ss=ss, holds=(sq == lam * ss)))
    b = [qe - lam for qe in q]
    y = solve_witness(R, b)
    return dict(basis=rows, uniform_is_stress=is_stress(R, [ONE] * len(q)), b=b, y=y,
                witness_ok=(y is not None and verify_witness(R, y, b)), all_hold=all(r['holds'] for r in rows))


def extreme_rays(S):
    """Extreme rays of the cone {w in the stress space : w >= 0}.  The stress space is spanned by the basis B
    (rows); w = alpha^T B.  A ray is extreme iff the constraints w_e = 0 that it meets have rank dim - 1.  We solve
    every (dim - 1)-subset of the 24 constraints exactly, keep the solutions with all w_e >= 0 (exact sign) up to
    sign, and deduplicate by exact proportionality.  Returns the list of rays, each normalized to sum w = 1."""
    B = S['stresses']
    d, m = len(B), len(B[0])
    cols = [[B[r][e] for r in range(d)] for e in range(m)]          # constraint e: alpha . cols[e] >= 0
    rays = []
    for sub in itertools.combinations(range(m), d - 1):
        M = dm([cols[e] for e in sub], d - 1, d)
        if M.rank() != d - 1:
            continue
        a = rows_of(M.nullspace())[0]
        w = [sum((a[r] * B[r][e] for r in range(d)), ZERO) for e in range(m)]
        sg = [sign(x) for x in w]
        if all(s >= 0 for s in sg):
            pass
        elif all(s <= 0 for s in sg):
            w = [-x for x in w]
        else:
            continue
        tot = sum(w, ZERO)
        w = [x / tot for x in w]
        if all(any(a != b for a, b in zip(w, r)) for r in rays):
            rays.append(w)
    return rays


def strictly_positive_stress(r):
    """A self-stress with every entry > 0, or None: the sum of the extreme rays (normalized to sum 1) is checked to
    be a stress (exact) with every entry of exact sign +1.  If it fails, None (the cone may still have one only if
    some edge is zero on every ray, in which case none exists)."""
    rays = r.get('rays') or []
    if not rays:
        return None
    m = len(rays[0])
    w = [sum((ray[e] for ray in rays), ZERO) for e in range(m)]
    tot = sum(w, ZERO)
    w = [x / tot for x in w]
    if is_stress(r['S']['R'], w) and all(sign(x) > 0 for x in w):
        return w
    return None


def ratio(w, q):
    """f^T Q_w f / sum w = (sum w_e q_e) / (sum w_e), exactly."""
    return sum((a * b for a, b in zip(w, q)), ZERO) / sum(w, ZERO)


# ---------------------------------------------------------------- part D: floating-point cross-checks (not proof)

def float_crosschecks(S, q, dps=50, seed=20261020):
    """Part D, floating point (mpmath, `dps` digits), labeled cross-check, not proof.
    (1) per-edge q_e by a central second difference of phi_e(exp(t f)) with h = 1e-15, computed from the closed
    form of the exponential map; (2) for three random integer combinations w of the basis stresses, the eigenvalue
    of the stress Laplacian on the (one-dimensional) flex space, f^T ((cD - W) (x) I_3) f with |f| = 1, compared
    with lam sum w.  Returns the maximum absolute deviations."""
    import mpmath as mp
    mp.mp.dps = dps
    tofl = lambda x: mp.mpf(1) * sum(mp.mpf(Fraction(c).numerator) / mp.mpf(Fraction(c).denominator) * r
                                     for c, r in zip(coords(x), (1, mp.sqrt(2), mp.sqrt(3), mp.sqrt(6))))
    U = [[tofl(x) for x in u] for u in S['U']]
    F = [[tofl(x) for x in v] for v in S['F']]
    nF = mp.sqrt(sum(sum(x * x for x in v) for v in F))
    f = [[x / nF for x in v] for v in F]
    norm_f = mp.sqrt(sum(sum(x * x for x in v) for v in f))

    def expmap(c, v, t):
        a = mp.sqrt(sum(x * x for x in v))
        if a == 0:
            return c
        return [mp.cos(t * a) * ci + mp.sin(t * a) / a * vi for ci, vi in zip(c, v)]

    h = mp.mpf('1e-15')
    fd = []
    for (i, j) in S['E']:
        ph = lambda t: sum(x * y for x, y in zip(expmap(U[i], f[i], t), expmap(U[j], f[j], t)))
        fd.append(-(ph(h) - 2 * ph(0) + ph(-h)) / h ** 2)
    dev_fd = max(abs(a - tofl(b)) for a, b in zip(fd, q))
    rnd = random.Random(seed)
    ev = []
    for _ in range(3):
        coef = [rnd.randint(-9, 9) for _ in S['stresses']]
        w = [sum(c * tofl(s[e]) for c, s in zip(coef, S['stresses'])) for e in range(len(S['E']))]
        Qf = mp.mpf(0)
        for wk, (i, j) in zip(w, S['E']):
            gij = sum(x * y for x, y in zip(U[i], U[j]))
            Qf += wk * (gij * (sum(x * x for x in f[i]) + sum(x * x for x in f[j])) - 2 * sum(x * y for x, y in zip(f[i], f[j])))
        ev.append((coef, Qf, sum(w)))
    return dict(norm_f=norm_f, dev_fd=dev_fd, fd=fd, eig=ev, dps=dps, h=h)


# ---------------------------------------------------------------- driver

def analyze(name, U, lam):
    """Run parts A-C for one configuration and return the data (used by the tests and by `main`)."""
    S = setup(U)
    out = dict(name=name, S=S, unit=check_unit(U), lam=lam)
    if S.get('flex_dim') == 1:
        q1, q2 = route1(S), route2(S)
        out.update(q1=q1, q2=q2, routes_agree=(q1 == q2))
        out['C'] = stress_space_report(S, q1, lam)
        out['rays'] = extreme_rays(S)
    return out


def main():
    """Print parts A-D for FCC (lambda = 1/6) and HCP (lambda = 3/22)."""
    import mpmath, platform
    print('versions: Python %s, sympy %s, mpmath %s' % (platform.python_version(), sp.__version__, mpmath.__version__))
    print('field K = %s, minimal polynomial of theta = sqrt2 + sqrt3: %s' % (K, K.mod.to_list()))
    for name, U, lam in [('FCC', fcc_points(), el(sp.Rational(1, 6))), ('HCP', hcp_points(), el(sp.Rational(3, 22)))]:
        r = analyze(name, U, lam)
        S = r['S']
        print('\n==== %s  (lambda tested = %s)' % (name, show(lam)))
        print('A  unit points: %s; contacts: %d; packing (all other pairs < 1/2): %s; tangent bases ok: %s'
              % (r['unit'], len(S['E']), S['packing'], S['tangent_ok']))
        print('A  edges:', S['E'])
        print('A  rank R = %d (R is %d x %d); dim left null space (self-stresses) = %d; dim null space = %d; '
              'rotations are motions: %s; flex dim in V = %d'
              % (S['rank'], S['R'].shape[0], S['R'].shape[1], len(S['stresses']), S['null_dim'], S['rot_ok'], S['flex_dim']))
        if S['flex_dim'] != 1:
            print('STOP: flex space is not one-dimensional'); continue
        print('A  flex F (ambient, unnormalized):', [[show(x) for x in v] for v in S['F']])
        print('A  |F|^2 = %s; f = F/|F| so |f| = 1; F tangent: %s; F perp rotations: %s; F is a flex: %s'
              % (show(S['normF2']), S['F_tangent'], S['F_perp_rot'], S['F_is_flex']))
        print('B  per-edge values (route 1 = stress Laplacian, route 2 = closed-form second derivative):')
        for e, a, b in zip(S['E'], r['q1'], r['q2']):
            print('     edge %-8s  %-10s %-10s %s' % (str(e), show(a), show(b), 'agree' if a == b else 'DISAGREE'))
        print('B  routes agree on all edges: %s; set of values: %s' % (r['routes_agree'], sorted(set(show(x) for x in r['q1']))))
        C = r['C']
        print('C  uniform stress is a self-stress: %s' % C['uniform_is_stress'])
        for k, row in enumerate(C['basis']):
            rr = show(row['sq'] / row['ss']) if row['ss'] != ZERO else 'n/a (sum 0)'
            print('C  basis stress %d: %s' % (k, [show(x) for x in row['s']]))
            print('     is stress: %s; sum s = %s; sum s q = %s; ratio = %s; sum s q == lam sum s: %s'
                  % (row['is_stress'], show(row['ss']), show(row['sq']), rr, row['holds']))
        print('C  witness y with R y = q - lam 1: %s' % (None if C['y'] is None else [show(x) for x in C['y']]))
        print('C  witness verified by exact multiplication: %s' % C['witness_ok'])
        rats = [ratio(w, r['q1']) for w in r['rays']]
        print('C  extreme rays of the nonnegative stress cone: %d; ratios f^T Q_w f / sum w: %s'
              % (len(r['rays']), sorted(set(show(x) for x in rats))))
        for w in r['rays']:
            print('     ray (sum 1):', [show(x) for x in w])
        pos = strictly_positive_stress(r)
        print('C  strictly positive self-stress exists: %s; witness (sum of the extreme rays, then sum 1): %s'
              % (pos is not None, None if pos is None else [show(x) for x in pos]))
        D = float_crosschecks(S, r['q1'])
        print('D  [floating point, mpmath %d digits, cross-check not proof] |f| = %s' % (D['dps'], mpmath.nstr(D['norm_f'], 30)))
        print('D  finite differences (h = %s): max |fd - q_exact| = %s' % (mpmath.nstr(D['h'], 3), mpmath.nstr(D['dev_fd'], 5)))
        lamf = mpmath.mpf(frac(lam).numerator) / frac(lam).denominator
        for coef, Qf, sw in D['eig']:
            print('D  w = %s . basis: eigenvalue on flex space %s; lam sum w %s; diff %s'
                  % (coef, mpmath.nstr(Qf, 25), mpmath.nstr(lamf * sw, 25), mpmath.nstr(Qf - lamf * sw, 5)))
        verdict = r['routes_agree'] and C['witness_ok'] and C['all_hold']
        print('VERDICT %s: f^T Q_w f = %s sum w on the whole stress space: %s' % (name, show(lam), verdict))
    return 0


if __name__ == '__main__':
    sys.exit(main())
