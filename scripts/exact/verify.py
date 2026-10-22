"""Exact verification of a critical configuration of 12 points on S^2 given by its Gram matrix over K = Q(c).

Status: untrusted (exploration or tooling; no part of the completeness proof).

Input: the field K (scripts/exact/kfield.py: P irreducible, root isolated exactly), the contact list, and a
Gram matrix whose entries are elements of K.  Every check below is decided in K (exact zero tests, exact
sign tests on the isolating interval); floating point is used only for the face structure, which is labelled
numerical.

Frame.  Pick three points B = (a, b, e) whose 3x3 Gram block M is positive definite.  If the Schur complement
of M in G vanishes, G is the Gram matrix of 12 unit vectors in R^3, unique up to O(3), and in the oblique
basis (u_a, u_b, u_e) with metric M the point u_i has coordinates y_i = M^{-1} g_i, where g_i = (G_ai, G_bi,
G_ei); inner products are x^T M y and M y_i = g_i.  Everything that follows (balance, rigidity, flexes, the
second-order form, symmetries) is linear algebra over K in this frame; no square roots are needed.

Conventions (EXACT-1 prompt and control reply of 2026-10-21):
  - stress (weight) space K(C): weights w on the contacts with sum_j w_ij (u_j - G_ij u_i) = 0 at every i;
  - rigidity rank: rank of the spherical rigidity matrix (contact rows on the tangent space (TS^2)^12);
    rank = k - dim K(C) is checked;
  - flex dimension: first-order motions keeping every contact length, modulo rotations;
  - kernel ("ker"): first-order motions keeping all contacts equal to a common, varying c, modulo rotations
    (the space of the certify.py bordered Hessian and of the 2026-10-02 note);
  - q(v) = 2 sum_e w_e <v_i, v_j> - sum_i lambda_i |v_i|^2, lambda_i = c sum_j w_ij, sum w = 1: the second
    derivative of c along the variety of equal contacts (c = sum w_e <u_i, u_j> there, and the ambient
    gradient of the Lagrangian vanishes by balance);
  - co-index = number of negative directions of q on ker (directions in which c decreases, i.e. r increases);
    descending index = number of positive directions (r decreases); nondegenerate iff q has exactly the
    3 rotation directions as null space on ker.
"""
import itertools
import mpmath as mp
import networkx as nx
from flint import fmpq_poly, fmpz_poly, fmpq, arb, ctx
from kfield import Field, KE, qstr, rref, rank, nullspace, det3, inv3, inertia, sturm_count


# ----- Gram matrix ------------------------------------------------------------------------------------------

def gram_checks(F, A, G, N=12):
    """Exact checks of the Gram matrix.  G: dict (i, j) -> KE for all i <= j.  Returns (report, frame)."""
    c = F.gen()
    Aset = set(A)
    g = lambda i, j: G[(min(i, j), max(i, j))]
    out = {}
    out['diagonal_is_one'] = all(g(i, i) == 1 for i in range(N))
    out['contacts_equal_c'] = all(g(i, j) == c for (i, j) in A)
    margins = []
    ok = True
    for i in range(N):
        for j in range(i + 1, N):
            if (i, j) in Aset:
                continue
            s = (c - g(i, j)).sign()
            ok &= s > 0
    out['noncontacts_below_c'] = ok
    # frame: first triple (in lexicographic order) with a positive definite 3x3 block
    B = None
    for tri in itertools.combinations(range(N), 3):
        a, b, e = tri
        if (1 - g(a, b) * g(a, b)).sign() <= 0:
            continue
        M = [[g(x, y) for y in tri] for x in tri]
        if det3(M).sign() > 0:
            B = tri
            break
    assert B is not None, 'no positive definite 3x3 block'
    M = [[g(x, y) for y in B] for x in B]
    Mi = inv3(M)
    gv = [[g(x, i) for x in B] for i in range(N)]
    Y = [[sum((Mi[r][s] * gv[i][s] for s in range(3)), F.const(0)) for r in range(3)] for i in range(N)]
    schur = all(g(i, j) == sum((gv[i][r] * Y[j][r] for r in range(3)), F.const(0))
                for i in range(N) for j in range(i, N))
    out['frame'] = list(B)
    out['frame_block_positive_definite'] = True  # leading minors 1, 1 - G_ab^2, det M all > 0 (checked above)
    out['schur_complement_zero'] = schur
    out['rank3_psd'] = schur
    assert out['diagonal_is_one'] and out['contacts_equal_c'] and out['noncontacts_below_c'] and schur, out
    return out, (B, M, Y, gv)


def ip(M, x, y):
    """<x, y> = x^T M y in the oblique frame."""
    return sum((x[r] * M[r][s] * y[s] for r in range(3) for s in range(3)), x[0].F.const(0))


# ----- stresses, rigidity, kernel, second-order form ------------------------------------------------------

def balance_matrix(F, A, G, frame, N=12):
    """Rows: the 3N balance equations sum_j w_ij (y_j - G_ij y_i) = 0 (oblique coordinates); columns: contacts."""
    B, M, Y, gv = frame
    g = lambda i, j: G[(min(i, j), max(i, j))]
    rows = [[F.const(0) for _ in A] for _ in range(3 * N)]
    for e, (i, j) in enumerate(A):
        for (a, b) in ((i, j), (j, i)):
            for t in range(3):
                rows[3 * a + t][e] = rows[3 * a + t][e] + Y[b][t] - g(a, b) * Y[a][t]
    return rows


def stress_space(F, A, G, frame):
    """Exact basis of the self-stress space K(C) (as lists of KE, one per contact)."""
    return nullspace(balance_matrix(F, A, G, frame))


def is_stress(F, A, G, frame, w):
    """Exact test that w balances every vertex."""
    rows = balance_matrix(F, A, G, frame)
    return all(sum((r[e] * w[e] for e in range(len(A))), F.const(0)).is_zero() for r in rows)


def motion_rows(F, A, frame, with_dc, N=12):
    """Linear constraints on tangent motions v (3N oblique coordinates, plus dc if with_dc):
    <u_i, v_i> = 0 and <v_i, u_j> + <u_i, v_j> = dc (or 0) on every contact."""
    B, M, Y, gv = frame
    n = 3 * N + (1 if with_dc else 0)
    rows = []
    for i in range(N):
        r = [F.const(0)] * n
        for t in range(3):
            r[3 * i + t] = gv[i][t]  # (M y_i)_t
        rows.append(r)
    for (i, j) in A:
        r = [F.const(0)] * n
        for t in range(3):
            r[3 * i + t] = r[3 * i + t] + gv[j][t]
            r[3 * j + t] = r[3 * j + t] + gv[i][t]
        if with_dc:
            r[3 * N] = F.const(-1)
        rows.append(r)
    return rows


def rotation_fields(F, frame, N=12):
    """The three infinitesimal rotations v_i = M^{-1} S y_i (S skew), which are M-orthogonal motions."""
    B, M, Y, gv = frame
    Mi = inv3(M)
    out = []
    for (p, q) in ((0, 1), (0, 2), (1, 2)):
        S = [[F.const(0)] * 3 for _ in range(3)]
        S[p][q], S[q][p] = F.const(1), F.const(-1)
        MS = [[sum((Mi[r][s] * S[s][u] for s in range(3)), F.const(0)) for u in range(3)] for r in range(3)]
        v = []
        for i in range(N):
            v += [sum((MS[r][u] * Y[i][u] for u in range(3)), F.const(0)) for r in range(3)]
        out.append(v)
    return out


def motions(F, A, frame, w_normalized=None, N=12):
    """Rigidity rank, flex dimension, kernel dimension, and (if a stress is given) the inertia of q on ker."""
    rows = motion_rows(F, A, frame, False)
    flex = nullspace(rows, 3 * N)
    rows_dc = motion_rows(F, A, frame, True)
    ker = nullspace(rows_dc, 3 * N + 1)
    rot = rotation_fields(F, frame)
    # the rotations satisfy every constraint and are independent
    for v in rot:
        assert all(sum((r[k] * v[k] for k in range(3 * N)), F.const(0)).is_zero() for r in rows)
    assert rank(rot) == 3
    out = {'rigidity_rank': 2 * N - len(flex), 'flex_dim_mod_rotations': len(flex) - 3,
           'kernel_dim_mod_rotations': len(ker) - 3}
    if w_normalized is not None:
        B, M, Y, gv = frame
        lam = [F.const(0)] * N
        c = F.gen()
        for e, (i, j) in enumerate(A):
            lam[i] = lam[i] + c * w_normalized[e]
            lam[j] = lam[j] + c * w_normalized[e]
        V = [[v[3 * i:3 * i + 3] for i in range(N)] for v in ker]

        def q(a, b):
            s = F.const(0)
            for e, (i, j) in enumerate(A):
                s = s + w_normalized[e] * (ip(M, a[i], b[j]) + ip(M, a[j], b[i]))
            for i in range(N):
                s = s - lam[i] * ip(M, a[i], b[i])
            return s
        Q = [[q(V[x], V[y]) for y in range(len(V))] for x in range(len(V))]
        pos, neg, zero = inertia(Q)
        out['q_inertia_on_ker'] = {'positive': pos, 'negative': neg, 'zero': zero}
        out['nondegenerate'] = zero == 3
        out['coindex'] = neg
        out['descending_index'] = pos
    return out


# ----- symmetry ------------------------------------------------------------------------------------------

def symmetries(F, A, G, frame, N=12):
    """Permutations of the points induced by isometries: contact-graph automorphisms that preserve G exactly.
    Returns a list of (perm, det) with det = +1 (rotation) or -1, decided exactly."""
    B, M, Y, gv = frame
    g = lambda i, j: G[(min(i, j), max(i, j))]
    H = nx.Graph()
    H.add_nodes_from(range(N))
    H.add_edges_from(A)
    out = []
    for iso in nx.algorithms.isomorphism.GraphMatcher(H, H).isomorphisms_iter():
        p = [iso[i] for i in range(N)]
        if all(g(p[i], p[j]) == g(i, j) for i in range(N) for j in range(i, N)):
            T = [[Y[p[B[x]]][r] for x in range(3)] for r in range(3)]  # columns y_{p(a)}, y_{p(b)}, y_{p(e)}
            d = det3(T)
            assert d == 1 or d == -1, 'induced map is not orthogonal'
            out.append((p, 1 if d == 1 else -1))
    return out


def _order(p):
    """Order of a permutation."""
    from math import gcd
    seen, o = set(), 1
    for i in range(len(p)):
        if i in seen:
            continue
        L, j = 0, i
        while j not in seen:
            seen.add(j)
            j = p[j]
            L += 1
        o = o * L // gcd(o, L)
    return o


def group_name(syms, G, N=12):
    """Schoenflies name of the symmetry group from its order, its rotation subgroup, the element orders and
    whether it contains -I (identified by G_{i, p(i)} = -1 for all i).  Only the groups met here are named."""
    order = len(syms)
    proper = [p for p, d in syms if d == 1]
    g = lambda i, j: G[(min(i, j), max(i, j))]
    has_inv = any(all(g(i, p[i]) == -1 for i in range(N)) for p, d in syms)
    rot_orders = sorted({_order(p) for p in proper})
    improper_orders = sorted({_order(p) for p, d in syms if d == -1})
    key = (order, len(proper), has_inv, max(rot_orders))
    names = {(120, 60, True, 5): 'I_h', (48, 24, True, 4): 'O_h', (12, 6, False, 3): 'D_3h',
             (24, 12, False, 6): 'D_6d', (2, 2, False, 2): 'C_2', (2, 1, True, 1): 'C_i',
             (2, 1, False, 1): 'C_s', (1, 1, False, 1): 'C_1'}
    name = names.get(key, 'unnamed')
    if name == 'D_6d':
        assert 12 in improper_orders, 'expected an S_12 element in D_6d'
    if name == 'D_3h':
        assert 6 in improper_orders, 'expected an S_3 (order 6) element in D_3h'
    return {'order': order, 'proper_order': len(proper), 'contains_inversion': has_inv,
            'rotation_orders': rot_orders, 'name': name}


# ----- canonical form of the contact graph --------------------------------------------------------------

def canonical_form(A, N=12):
    """Canonical edge list of the contact graph: the lexicographically least sorted edge list over the
    labelings reached by individualization-refinement (colour refinement is label-independent, so the
    result depends only on the isomorphism class)."""
    adj = [set() for _ in range(N)]
    for i, j in A:
        adj[i].add(j)
        adj[j].add(i)

    def refine(col):
        while True:
            sig = [(col[v], tuple(sorted(col[u] for u in adj[v]))) for v in range(N)]
            keys = sorted(set(sig))
            new = [keys.index(s) for s in sig]
            if len(set(new)) == len(set(col)):
                return new
            col = new

    best = [None]

    def search(col):
        col = refine(col)
        if len(set(col)) == N:
            lab = col  # vertex v gets label col[v]
            code = sorted(tuple(sorted((lab[i], lab[j]))) for i, j in A)
            if best[0] is None or code < best[0]:
                best[0] = code
            return
        cells = {}
        for v in range(N):
            cells.setdefault(col[v], []).append(v)
        target = min((c for c in cells if len(cells[c]) > 1), key=lambda c: (len(cells[c]), c))
        for v in cells[target]:
            new = [2 * x + 1 for x in col]
            new[v] = 2 * col[v]
            search(new)
    search([0] * N)
    return [list(e) for e in best[0]]


# ----- faces (numerical) --------------------------------------------------------------------------------

def faces(X, A, N=12):
    """Face sizes of the geodesic embedding of the contact graph on S^2, from the rotation system at each
    vertex (numerical, from coordinates; returns sizes and the smallest angular gap between neighbours)."""
    adj = {i: [] for i in range(N)}
    for i, j in A:
        adj[i].append(j)
        adj[j].append(i)
    order, gap = {}, mp.inf
    for i in range(N):
        u = X[i]
        a = [1, 0, 0] if abs(u[0]) < 0.9 else [0, 1, 0]
        e1 = [a[1] * u[2] - a[2] * u[1], a[2] * u[0] - a[0] * u[2], a[0] * u[1] - a[1] * u[0]]
        n1 = mp.sqrt(sum(x * x for x in e1))
        e1 = [x / n1 for x in e1]
        e2 = [u[1] * e1[2] - u[2] * e1[1], u[2] * e1[0] - u[0] * e1[2], u[0] * e1[1] - u[1] * e1[0]]
        ang = sorted((mp.atan2(sum(X[j][t] * e2[t] for t in range(3)), sum(X[j][t] * e1[t] for t in range(3))), j)
                     for j in adj[i])
        order[i] = [j for _, j in ang]
        for k in range(len(ang)):
            d = (ang[(k + 1) % len(ang)][0] - ang[k][0]) % (2 * mp.pi)
            gap = min(gap, d)
    used, sizes = set(), []
    for i in range(N):
        for j in adj[i]:
            if (i, j) in used:
                continue
            f, a, b = 0, i, j
            while (a, b) not in used:
                used.add((a, b))
                f += 1
                L = order[b]
                c = L[(L.index(a) - 1) % len(L)]
                a, b = b, c
            sizes.append(f)
    return sorted(sizes), gap


# ----- s and r --------------------------------------------------------------------------------------------

def _int_primitive(p):
    """Primitive integer polynomial with positive leading coefficient, from a rational polynomial."""
    num = fmpz_poly([int((p[k] * p.denom()).p) for k in range(p.degree() + 1)])
    c = num.content()
    num = fmpz_poly([int(x) // int(c) for x in num.coeffs()])
    if num.leading_coefficient() < 0:
        num = -num
    return num


def _select_factor(fac, x):
    """The unique factor in `fac` that vanishes at the Arb ball x (others must exclude 0)."""
    hits = []
    for f, m in fac:
        v = arb(0)
        for k in range(f.degree(), -1, -1):
            v = v * x + arb(int(f[k]))
        if not (v > 0 or v < 0):
            hits.append(f)
    return hits


def isolate(fz, x, approx_digits=40):
    """Rational interval around the Arb ball x whose Sturm count for fz is exactly 1."""
    fq = fmpq_poly([int(a) for a in fz.coeffs()])
    # rational endpoints from the ball: [mid - 2 rad - 10^-D, mid + 2 rad + 10^-D]
    m = _arb_to_fmpq(x.mid())
    r = _arb_to_fmpq(x.rad()) * 2 + fmpq(1, 10 ** approx_digits)
    lo, hi = m - r, m + r
    assert fq(lo) != 0 and fq(hi) != 0
    assert sturm_count(fq, lo, hi) == 1
    return lo, hi


def _arb_to_fmpq(a):
    """The exact rational value of an exact Arb (a midpoint or radius)."""
    m, e = a.man_exp()
    m, e = int(m), int(e)
    return fmpq(m * 2 ** e) if e >= 0 else fmpq(m, 2 ** (-e))


def s_and_r(F):
    """Minimal polynomials over Q (integer, primitive, constant first) and isolating intervals for
    s^2 = (1 - c)/2, s = sqrt(s^2) > 0 and r = s/(1 - s)."""
    d = F.d
    # s^2 = x  <=>  c = 1 - 2x
    P2 = F.P(fmpq_poly([1, -2]))
    Q2 = _int_primitive(P2)
    # s: Q2(x^2)
    Qs = fmpz_poly(sum(([int(a), 0] for a in Q2.coeffs()), [])[:-1])
    # r: s = r/(1+r), so x = r^2/(1+r)^2; homogenize with (1+r)^(2d)
    one_r = fmpz_poly([1, 1])
    Qr = fmpz_poly([0])
    co = Q2.coeffs()
    for k, a in enumerate(co):
        Qr += int(a) * fmpz_poly([0, 0, 1]) ** k * one_r ** (2 * (len(co) - 1) - 2 * k)
    old = ctx.prec
    try:
        ctx.prec = 512
        while F.hi - F.lo > fmpq(1, 2 ** 300):
            F.refine()
        cb = F.ball()
        s2 = (1 - cb) / 2
        s = s2.sqrt()
        r = s / (1 - s)
        out = {}
        for name, poly, x in (('s2', Q2, s2), ('s', Qs, s), ('r', Qr, r)):
            content, fac = poly.factor()
            hits = _select_factor(fac, x)
            assert len(hits) == 1, (name, hits)
            f = hits[0]
            if f.leading_coefficient() < 0:
                f = -f
            lo, hi = isolate(f, x)
            out[name] = {'minpoly': [str(int(a)) for a in f.coeffs()], 'degree': f.degree(),
                         'interval': [qstr(lo), qstr(hi)], 'arb': x.str(30, radius=True)}
        return out
    finally:
        ctx.prec = old
