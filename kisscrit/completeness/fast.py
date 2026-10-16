"""Compiled (numba) kernels for the completeness branch and bound.

Status: untrusted (the search (numba branch and bound). Soundness does not depend on it: it proposes proof trees that the checker re-verifies).

Same mathematics as model.py / bnb2.py (which remain the readable reference implementation and are used by the
tests to cross-check this module), reorganised for speed: one box at a time, all in nopython mode.

Rigour.  Interval results are rounded outward by one ulp per arithmetic operation (np.nextafter); libm results
(cos, sin, atan, asin, acos, sqrt) are widened by 4 ulp and by an absolute 1e-15.  Linear programs are solved in
floating point only to choose multipliers; every exclusion or contraction is then re-derived from interval
enclosures with those multipliers (preconditioned interval Gauss-Seidel / Gordan check), so LP round-off cannot
make the result unsound.
"""
import math
import numpy as np
import numba as nb
from numba import njit

INF = np.inf
PI = np.pi
PI_LO = PI; PI_HI = float(np.nextafter(PI, INF))
HPI_LO = PI / 2; HPI_HI = float(np.nextafter(PI / 2, INF))
TWOPI_LO = 2 * PI; TWOPI_HI = float(np.nextafter(2 * PI, INF))
DELTA = float(np.radians(25.0))
LSTRAIGHT = PI - DELTA          # region threshold for near-straight pentagon corners
EPS_ABS = 1e-15
MONO = True        # monotonicity tightening of pentagon corners on wide boxes
MONO_W = 0.01


# ============================================================================ interval scalars
# Outward rounding without nextafter: for finite x, x - (|x| * 2.3e-16 + 1e-300) computed in round-to-nearest is
# at most the float below x (2.3e-16 > 2^-52 >= ulp(x)/|x|, so the exact difference lies at least one ulp below
# x and rounding to nearest cannot return x).  Infinite values are left unchanged.
@njit(cache=True, inline='always')
def dn(x):
    if x == INF: return 1.7976931348623157e308
    return x - (abs(x) * 2.3e-16 + 1e-300)


@njit(cache=True, inline='always')
def up(x):
    if x == -INF: return -1.7976931348623157e308
    return x + (abs(x) * 2.3e-16 + 1e-300)


@njit(cache=True, inline='always')
def dnl(x):  # libm widening: >= 4 ulp and >= 1e-15
    if x == INF: return 1.7976931348623157e308
    return x - (abs(x) * 9.0e-16 + EPS_ABS)


@njit(cache=True, inline='always')
def upl(x):
    if x == -INF: return -1.7976931348623157e308
    return x + (abs(x) * 9.0e-16 + EPS_ABS)


@njit(cache=True, inline='always')
def iadd(al, ah, bl, bh): return dn(al + bl), up(ah + bh)


@njit(cache=True, inline='always')
def isub(al, ah, bl, bh): return dn(al - bh), up(ah - bl)


@njit(cache=True)
def imul(al, ah, bl, bh):
    p1 = al * bl; p2 = al * bh; p3 = ah * bl; p4 = ah * bh
    if math.isnan(p1) or math.isnan(p2) or math.isnan(p3) or math.isnan(p4):
        # 0 * inf: be conservative unless one factor is exactly [0, 0]
        if (al == 0.0 and ah == 0.0) or (bl == 0.0 and bh == 0.0): return 0.0, 0.0
        return -INF, INF
    lo = min(min(p1, p2), min(p3, p4)); hi = max(max(p1, p2), max(p3, p4))
    return dn(lo), up(hi)


@njit(cache=True)
def idiv(al, ah, bl, bh):
    if bl <= 0.0 <= bh: return -INF, INF
    rl = dn(1.0 / bh); rh = up(1.0 / bl)
    return imul(al, ah, rl, rh)


@njit(cache=True)
def icos(l, h):
    if not (h - l < 2 * PI) or math.isnan(l) or math.isnan(h): return -1.0, 1.0
    ca = math.cos(l); cb = math.cos(h)
    lo = min(ca, cb); hi = max(ca, cb)
    k1 = math.ceil(l / (2 * PI))
    if k1 * 2 * PI <= h: hi = 1.0
    k2 = math.ceil((l - PI) / (2 * PI))
    if k2 * 2 * PI + PI <= h: lo = -1.0
    return max(dnl(lo), -1.0), min(upl(hi), 1.0)


@njit(cache=True)
def isin(l, h):
    # sin x = cos(x - pi/2); the float pi/2 differs from pi/2 by < 1e-16, absorbed by the 1e-15 widening
    return icos(dn(l - HPI_LO), up(h - HPI_LO))


@njit(cache=True)
def iatan(l, h):
    return max(dnl(math.atan(l)) if l > -INF else -HPI_HI, -HPI_HI), min(upl(math.atan(h)) if h < INF else HPI_HI, HPI_HI)


@njit(cache=True)
def iasin(l, h):
    l = min(max(l, -1.0), 1.0); h = min(max(h, -1.0), 1.0)
    return max(dnl(math.asin(l)), -HPI_HI), min(upl(math.asin(h)), HPI_HI)


@njit(cache=True)
def iacos(l, h):
    l = min(max(l, -1.0), 1.0); h = min(max(h, -1.0), 1.0)
    return max(dnl(math.acos(h)), 0.0), min(upl(math.acos(l)), PI_HI)


# ============================================================================ AD over intervals
# Local forward AD: a value is an 8-tuple (vl, vh, g0l, g0h, g1l, g1h, g2l, g2h); the three gradient components
# are with respect to theta and the (at most two) parameters of the face being evaluated.

@njit(cache=True, inline='always')
def V3(l, h, k):
    if k == 0: return (l, h, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0)
    if k == 1: return (l, h, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0)
    return (l, h, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0)


@njit(cache=True, inline='always')
def ad_add(a, b):
    v = iadd(a[0], a[1], b[0], b[1]); g0 = iadd(a[2], a[3], b[2], b[3])
    g1 = iadd(a[4], a[5], b[4], b[5]); g2 = iadd(a[6], a[7], b[6], b[7])
    return (v[0], v[1], g0[0], g0[1], g1[0], g1[1], g2[0], g2[1])


@njit(cache=True, inline='always')
def ad_sub(a, b):
    v = isub(a[0], a[1], b[0], b[1]); g0 = isub(a[2], a[3], b[2], b[3])
    g1 = isub(a[4], a[5], b[4], b[5]); g2 = isub(a[6], a[7], b[6], b[7])
    return (v[0], v[1], g0[0], g0[1], g1[0], g1[1], g2[0], g2[1])


@njit(cache=True, inline='always')
def ad_chain(a, vl, vh, dl, dh):
    g0 = imul(a[2], a[3], dl, dh); g1 = imul(a[4], a[5], dl, dh); g2 = imul(a[6], a[7], dl, dh)
    return (vl, vh, g0[0], g0[1], g1[0], g1[1], g2[0], g2[1])


@njit(cache=True, inline='always')
def ad_scale(a, cl, ch):
    v = imul(a[0], a[1], cl, ch)
    return ad_chain(a, v[0], v[1], cl, ch)


@njit(cache=True, inline='always')
def ad_shift(a, cl, ch):
    v = iadd(a[0], a[1], cl, ch)
    return (v[0], v[1], a[2], a[3], a[4], a[5], a[6], a[7])


@njit(cache=True, inline='always')
def ad_rsub(cl, ch, a):
    v = isub(cl, ch, a[0], a[1])
    return (v[0], v[1], -a[3], -a[2], -a[5], -a[4], -a[7], -a[6])


@njit(cache=True, inline='always')
def _mulg(ga, gb, va, vb):
    x = imul(ga[0], ga[1], vb[0], vb[1]); y = imul(gb[0], gb[1], va[0], va[1])
    return iadd(x[0], x[1], y[0], y[1])


@njit(cache=True)
def ad_mul(a, b):
    v = imul(a[0], a[1], b[0], b[1])
    va = (a[0], a[1]); vb = (b[0], b[1])
    g0 = _mulg((a[2], a[3]), (b[2], b[3]), va, vb)
    g1 = _mulg((a[4], a[5]), (b[4], b[5]), va, vb)
    g2 = _mulg((a[6], a[7]), (b[6], b[7]), va, vb)
    return (v[0], v[1], g0[0], g0[1], g1[0], g1[1], g2[0], g2[1])


@njit(cache=True, inline='always')
def _divg(ga, gb, q, vb):
    x = imul(gb[0], gb[1], q[0], q[1]); y = isub(ga[0], ga[1], x[0], x[1])
    return idiv(y[0], y[1], vb[0], vb[1])


@njit(cache=True)
def ad_div(a, b):
    q = idiv(a[0], a[1], b[0], b[1]); vb = (b[0], b[1])
    g0 = _divg((a[2], a[3]), (b[2], b[3]), q, vb)
    g1 = _divg((a[4], a[5]), (b[4], b[5]), q, vb)
    g2 = _divg((a[6], a[7]), (b[6], b[7]), q, vb)
    return (q[0], q[1], g0[0], g0[1], g1[0], g1[1], g2[0], g2[1])


@njit(cache=True)
def ad_cos(a):
    v = icos(a[0], a[1]); s = isin(a[0], a[1])
    return ad_chain(a, v[0], v[1], -s[1], -s[0])


@njit(cache=True)
def ad_sin(a):
    v = isin(a[0], a[1]); c = icos(a[0], a[1])
    return ad_chain(a, v[0], v[1], c[0], c[1])


@njit(cache=True)
def ad_atan(a):
    v = iatan(a[0], a[1])
    s = imul(a[0], a[1], a[0], a[1])
    d = iadd(1.0, 1.0, max(s[0], 0.0), s[1])
    r = idiv(1.0, 1.0, d[0], d[1])
    return ad_chain(a, v[0], v[1], r[0], r[1])


@njit(cache=True)
def ad_clip(a, lo, hi):
    """clamp to [lo, hi]; returns (value, nonempty flag).  Gradient times [0, 1] where the clamp may act."""
    ok = (a[1] >= lo) and (a[0] <= hi)
    active = (a[0] < lo) or (a[1] > hi)
    vl = min(max(a[0], lo), hi); vh = min(max(a[1], lo), hi)
    if active:
        return ad_chain(a, vl, vh, 0.0, 1.0), ok
    return (vl, vh, a[2], a[3], a[4], a[5], a[6], a[7]), ok


@njit(cache=True)
def ad_asin_angle(a, lo, hi):
    """y = clamp(asin(clamp(x, -1, 1)), lo, hi), dy/dx = 1/cos(y) times clamp factors."""
    ok = (a[0] <= 1.0) and (a[1] >= -1.0)
    y = iasin(a[0], a[1])
    ok = ok and (y[1] >= lo) and (y[0] <= hi)
    active = (a[0] < -1.0) or (a[1] > 1.0) or (y[0] < lo) or (y[1] > hi)
    yl = min(max(y[0], lo), hi); yh = min(max(y[1], lo), hi)
    c = icos(yl, yh)
    d = idiv(1.0, 1.0, c[0], c[1])
    if active: d = imul(d[0], d[1], 0.0, 1.0)
    return ad_chain(a, yl, yh, d[0], d[1]), ok


@njit(cache=True)
def ad_acos(a):
    ok = (a[0] <= 1.0) and (a[1] >= -1.0)
    y = iacos(a[0], a[1])
    s = isin(y[0], y[1])
    d = idiv(-1.0, -1.0, s[0], s[1])
    return ad_chain(a, y[0], y[1], d[0], d[1]), ok


# ============================================================================ face corners
@njit(cache=True)
def alpha_lo(t):
    c = math.cos(t); return dnl(math.acos(c / (1 + c)))


@njit(cache=True)
def alpha_hi(t):
    c = math.cos(t); return upl(math.acos(c / (1 + c)))


@njit(cache=True)
def ad_beta(ct, A):
    half = ad_scale(A, 0.5, 0.5)
    return ad_atan(ad_div(ad_cos(half), ad_mul(ct, ad_sin(half))))


@njit(cache=True)
def pentagon(ct, st, hmin, hmax, f1, f2):
    """Formula-order corners (c0, f1, f2, c3, f4) and validity."""
    half1 = ad_scale(f1, 0.5, 0.5)
    hd2, ok = ad_asin_angle(ad_mul(st, ad_sin(half1)), 0.0, HPI_LO)
    d2 = ad_scale(hd2, 2.0, 2.0)
    b1 = ad_beta(ct, f1)
    g, o = ad_clip(ad_sub(f2, b1), 0.0, PI_LO); ok = ok and o
    cg = ad_cos(g); sg = ad_sin(g)
    cd2 = ad_cos(d2); sd2 = ad_sin(d2)
    cd3 = ad_add(ad_mul(cd2, ct), ad_mul(ad_mul(sd2, st), cg))
    ok = ok and (cd3[0] <= ct[1])
    cd3, o = ad_clip(cd3, -1.0, 1.0)
    d3, o = ad_acos(cd3); ok = ok and o
    x0 = ad_sub(ad_mul(sd2, ct), ad_mul(ad_mul(cd2, st), cg))
    a0 = ad_rsub(HPI_LO, HPI_HI, ad_atan(ad_div(x0, ad_mul(sg, st))))
    x3 = ad_sub(ad_mul(st, cd2), ad_mul(ad_mul(ct, sd2), cg))
    a3 = ad_rsub(HPI_LO, HPI_HI, ad_atan(ad_div(x3, ad_mul(sg, sd2))))
    h, o = ad_asin_angle(ad_div(ad_sin(ad_scale(d3, 0.5, 0.5)), st), hmin, hmax); ok = ok and o
    f4 = ad_scale(h, 2.0, 2.0)
    b3 = ad_beta(ct, f4)
    c0 = ad_add(ad_add(b1, a0), b3)
    c3 = ad_add(a3, b3)
    return c0, f1, f2, c3, f4, ok


@njit(cache=True)
def _pent_at(thl, thh, f1l, f1h, f2l, f2h, hmin, hmax):
    th = V3(thl, thh, 0)
    ct = ad_cos(th); st = ad_sin(th)
    return pentagon(ct, st, hmin, hmax, V3(f1l, f1h, 1), V3(f2l, f2h, 2))


@njit(cache=True)
def _mono_pick(g_lo, g_hi, lo, hi, want_min):
    """Endpoint (or the whole interval) of a variable at which a function with derivative in [g_lo, g_hi]
    attains its minimum (want_min) or maximum over [lo, hi]."""
    if g_lo >= 0.0: return (lo, lo) if want_min else (hi, hi)
    if g_hi <= 0.0: return (hi, hi) if want_min else (lo, lo)
    return lo, hi


@njit(cache=True)
def _tighten(c, j, thl, thh, f1l, f1h, f2l, f2h, hmin, hmax):
    """Monotonicity tightening of the value enclosure of computed pentagon corner j (0, 3 or 4)."""
    lo, hi = c[0], c[1]
    for want_min in (True, False):
        a = _mono_pick(c[2], c[3], thl, thh, want_min)
        b = _mono_pick(c[4], c[5], f1l, f1h, want_min)
        d = _mono_pick(c[6], c[7], f2l, f2h, want_min)
        if a[0] == thl and a[1] == thh and b[0] == f1l and b[1] == f1h and d[0] == f2l and d[1] == f2h:
            continue
        r = _pent_at(a[0], a[1], b[0], b[1], d[0], d[1], hmin, hmax)
        cc = r[0] if j == 0 else (r[3] if j == 3 else r[4])
        if want_min: lo = max(lo, cc[0])
        else: hi = min(hi, cc[1])
    if lo > hi:      # cannot happen in exact arithmetic; keep the original enclosure
        return c
    return (lo, hi, c[2], c[3], c[4], c[5], c[6], c[7])


@njit(cache=True, inline='always')
def _store(C, f, q, c):
    for t in range(8): C[f, q, t] = c[t]


@njit(cache=True)
def corners_box(fsize, fp0, fp1, fmode, lo, hi, blo, bhi, C):
    """Corner enclosures with local gradients (theta, face params) for the box [lo, hi], clamp bounds from
    [blo, bhi].  C has shape (nF, 5, 8).  Returns validity."""
    nF = fsize.shape[0]
    th = V3(lo[0], hi[0], 0)
    ct = ad_cos(th); st = ad_sin(th)
    alad, _ = ad_acos(ad_div(ct, ad_shift(ct, 1.0, 1.0)))
    al = (alpha_lo(lo[0]), alpha_hi(hi[0]), alad[2], alad[3], 0.0, 0.0, 0.0, 0.0)
    ok = True
    hmin = 0.5 * alpha_lo(blo[0])
    for f in range(nF):
        k = fsize[f]
        if k == 3:
            for q in range(3): _store(C, f, q, al)
        elif k == 4:
            i = fp0[f]
            A = V3(lo[i], hi[i], 1)
            O = ad_scale(ad_beta(ct, A), 2.0, 2.0)
            exl = dnl(2 * math.atan2(math.cos(hi[i] / 2), math.cos(lo[0]) * math.sin(hi[i] / 2)))
            exh = upl(2 * math.atan2(math.cos(lo[i] / 2), math.cos(hi[0]) * math.sin(lo[i] / 2)))
            O = (max(O[0], exl), min(O[1], exh), O[2], O[3], O[4], O[5], O[6], O[7])
            _store(C, f, 0, A); _store(C, f, 1, O); _store(C, f, 2, A); _store(C, f, 3, O)
        else:
            i1 = fp0[f]; i2 = fp1[f]
            f1 = V3(lo[i1], hi[i1], 1); f2 = V3(lo[i2], hi[i2], 2)
            c0, a1, a2, c3, f4, o = pentagon(ct, st, hmin, 0.5 * LSTRAIGHT, f1, f2)
            ok = ok and o
            if MONO and max(hi[0] - lo[0], max(hi[i1] - lo[i1], hi[i2] - lo[i2])) > MONO_W:
                args = (lo[0], hi[0], lo[i1], hi[i1], lo[i2], hi[i2], hmin, 0.5 * LSTRAIGHT)
                c0 = _tighten(c0, 0, *args); c3 = _tighten(c3, 3, *args); f4 = _tighten(f4, 4, *args)
            mode = fmode[f]
            s = 0 if mode == 0 else (2 if mode == 1 else 3)
            if mode == 2: ok = ok and (c3[1] >= LSTRAIGHT)
            _store(C, f, s % 5, c0); _store(C, f, (s + 1) % 5, a1); _store(C, f, (s + 2) % 5, a2)
            _store(C, f, (s + 3) % 5, c3); _store(C, f, (s + 4) % 5, f4)
    # every corner lies in [alpha(theta), pi] (corner >= alpha: the two neighbours are >= theta apart)
    for f in range(nF):
        for q in range(fsize[f]):
            if C[f, q, 0] > PI_HI or C[f, q, 1] < al[0]: ok = False
    return ok


@njit(cache=True)
def vertex_sums(fp0, fp1, vstart, vface, vpos, C, n, FL, FH, JL, JH):
    """Sums minus 2 pi with global gradients.  FL/FH (V,), JL/JH (V, n)."""
    V = vstart.shape[0] - 1
    vok = True
    for v in range(V):
        sl, sh = isub(0.0, 0.0, TWOPI_LO, TWOPI_HI)
        for j in range(n): JL[v, j] = 0.0; JH[v, j] = 0.0
        for t in range(vstart[v], vstart[v + 1]):
            f = vface[t]; q = vpos[t]
            sl, sh = iadd(sl, sh, C[f, q, 0], C[f, q, 1])
            JL[v, 0], JH[v, 0] = iadd(JL[v, 0], JH[v, 0], C[f, q, 2], C[f, q, 3])
            if fp0[f] >= 0:
                k = fp0[f]; JL[v, k], JH[v, k] = iadd(JL[v, k], JH[v, k], C[f, q, 4], C[f, q, 5])
            if fp1[f] >= 0:
                k = fp1[f]; JL[v, k], JH[v, k] = iadd(JL[v, k], JH[v, k], C[f, q, 6], C[f, q, 7])
        FL[v] = sl; FH[v] = sh
        if not (sl <= 0.0 <= sh): vok = False
    return vok


# ============================================================================ small dense LP (simplex)
# Tableau for the linearised vertex equations of one box.  Unknowns z = y - ylo in [0, w] (y = p - m) and t >= 0:
#     Jc y - t <= R - fc ,   -Jc y - t <= R + fc        (2V rows)        z_j <= w_j     (n rows)
# Columns: z (n), t (1), slacks (2V + n), rhs.  The initial pivot on t makes the slack basis feasible.
# The reduced cost of the slack of row i is the multiplier lambda_i >= 0 of that row.

@njit(cache=True)
def _pivot(T, basis, r, e):
    piv = T[r, e]
    T[r] /= piv
    for i in range(T.shape[0]):
        if i != r:
            f = T[i, e]
            if f != 0.0: T[i] -= f * T[r]
    basis[r] = e


@njit(cache=True)
def _simplex(T, basis, maxit):
    """Primal simplex on a feasible tableau.  Dantzig pricing; after 50 iterations Bland's rule (smallest
    entering index, smallest leaving basis index on ratio ties), which cannot cycle.  Status 0 optimal."""
    M = T.shape[0] - 1; ncol = T.shape[1] - 1
    for it in range(maxit):
        bland = it > 50
        e = -1; best = -1e-10
        for j in range(ncol):
            if T[M, j] < best:
                e = j
                if bland: break
                best = T[M, j]
        if e < 0: return 0
        r = -1; bratio = INF
        for i in range(M):
            if T[i, e] > 1e-12:
                ratio = T[i, ncol] / T[i, e]
                if ratio < bratio - 1e-13 or (bland and ratio <= bratio + 1e-13 and r >= 0 and basis[i] < basis[r]):
                    if ratio < bratio: bratio = ratio
                    r = i
        if r < 0: return 3
        _pivot(T, basis, r, e)
    return 4


@njit(cache=True)
def lp_setup(Jc, fc, R, ylo, w):
    V, n = Jc.shape
    M = 2 * V + n; ncol = n + 1 + M
    T = np.zeros((M + 1, ncol + 1)); basis = np.empty(M, np.int64)
    for v in range(V):
        # row v:  Jc (ylo + z) - t <= R - fc   ->  Jc z - t <= R - fc - Jc ylo
        a = 0.0
        for j in range(n): a += Jc[v, j] * ylo[j]
        for j in range(n): T[v, j] = Jc[v, j]; T[V + v, j] = -Jc[v, j]
        T[v, n] = -1.0; T[V + v, n] = -1.0
        T[v, ncol] = R[v] - fc[v] - a
        T[V + v, ncol] = R[v] + fc[v] + a
    for j in range(n):
        T[2 * V + j, j] = 1.0; T[2 * V + j, ncol] = w[j]
    for i in range(M):
        T[i, n + 1 + i] = 1.0; basis[i] = n + 1 + i
    r = 0
    for i in range(2 * V):
        if T[i, ncol] < T[r, ncol]: r = i
    if T[r, ncol] < 0: _pivot(T, basis, r, n)
    return T, basis


@njit(cache=True)
def lp_set_objective(T, basis, cost):
    M = T.shape[0] - 1; ncol = T.shape[1] - 1
    T[M, :] = 0.0
    for j in range(cost.shape[0]): T[M, j] = cost[j]
    for i in range(M):
        bj = basis[i]
        if bj < cost.shape[0] and cost[bj] != 0.0:
            T[M] -= cost[bj] * T[i]


@njit(cache=True)
def lp_multipliers(T, V, n):
    """mu = lambda(+) - lambda(-) for the vertex equations (length V)."""
    M = T.shape[0] - 1
    mu = np.empty(V)
    for v in range(V):
        mu[v] = T[M, n + 1 + v] - T[M, n + 1 + V + v]
    return mu


# ============================================================================ rigorous steps
@njit(cache=True)
def combo(mu, FL, FH, JL, JH):
    """Interval s0 = mu.F(m) and row r = mu^T J(box)."""
    V, n = JL.shape
    s0l, s0h = 0.0, 0.0
    rl = np.zeros(n); rh = np.zeros(n)
    for v in range(V):
        a, b = imul(mu[v], mu[v], FL[v], FH[v]); s0l, s0h = iadd(s0l, s0h, a, b)
        for j in range(n):
            a, b = imul(mu[v], mu[v], JL[v, j], JH[v, j]); rl[j], rh[j] = iadd(rl[j], rh[j], a, b)
    return s0l, s0h, rl, rh


@njit(cache=True)
def excluded_by(mu, FL, FH, JL, JH, ylo, yhi):
    s0l, s0h, rl, rh = combo(mu, FL, FH, JL, JH)
    sl, sh = s0l, s0h
    for j in range(ylo.shape[0]):
        a, b = imul(rl[j], rh[j], ylo[j], yhi[j]); sl, sh = iadd(sl, sh, a, b)
    return sl > 0.0 or sh < 0.0


@njit(cache=True)
def gs_update(mu, k, FL, FH, JL, JH, ylo, yhi):
    """New bounds for y_k from 0 in mu.F(m) + sum_j (mu^T J)_j y_j; returns (lo, hi, empty)."""
    s0l, s0h, rl, rh = combo(mu, FL, FH, JL, JH)
    tl, th = -s0h, -s0l
    for j in range(ylo.shape[0]):
        if j == k: continue
        a, b = imul(rl[j], rh[j], ylo[j], yhi[j]); tl, th = isub(tl, th, a, b)
    if rl[k] <= 0.0 <= rh[k]: return ylo[k], yhi[k], False
    ql, qh = idiv(tl, th, rl[k], rh[k])
    nl = max(ylo[k], ql); nh = min(yhi[k], qh)
    if nl > nh: return ylo[k], yhi[k], True
    return nl, nh, False


# ============================================================================ per-box pipeline
@njit(cache=True)
def eval_all(fsize, fp0, fp1, fmode, vstart, vface, vpos, lo, hi, blo, bhi):
    nF = fsize.shape[0]; n = lo.shape[0]; V = vstart.shape[0] - 1
    C = np.zeros((nF, 5, 8))
    FL = np.zeros(V); FH = np.zeros(V); JL = np.zeros((V, n)); JH = np.zeros((V, n))
    ok = corners_box(fsize, fp0, fp1, fmode, lo, hi, blo, bhi, C)
    vok = vertex_sums(fp0, fp1, vstart, vface, vpos, C, n, FL, FH, JL, JH)
    return ok, vok, C[:, :, 0].copy(), C[:, :, 1].copy(), FL, FH, JL, JH


# Contracted bounds are rounded outward to a grid of 2^16 steps across the box the contraction started from, so
# that a certificate records each changed bound in two bytes.  The grid value is computed exactly this way by
# the checker:  q == 0 -> lo0;  q == QMAX -> hi0;  otherwise lo0 + (hi0 - lo0) * (q / QMAX).
QMAX = 65536


@njit(cache=True)
def grid_value(lo0, hi0, q):
    if q <= 0: return lo0
    if q >= QMAX: return hi0
    return lo0 + (hi0 - lo0) * (q / 65536.0)


@njit(cache=True)
def grid_down(lo0, hi0, x):
    """Largest grid index q with grid_value(q) <= x (x in [lo0, hi0])."""
    if hi0 <= lo0: return 0
    q = int(math.floor((x - lo0) / (hi0 - lo0) * QMAX))
    q = min(max(q, 0), QMAX)
    while q > 0 and grid_value(lo0, hi0, q) > x: q -= 1
    while q < QMAX and grid_value(lo0, hi0, q + 1) <= x: q += 1
    return q


@njit(cache=True)
def grid_up(lo0, hi0, x):
    """Smallest grid index q with grid_value(q) >= x."""
    if hi0 <= lo0: return QMAX
    q = int(math.ceil((x - lo0) / (hi0 - lo0) * QMAX))
    q = min(max(q, 0), QMAX)
    while q < QMAX and grid_value(lo0, hi0, q) < x: q += 1
    while q > 0 and grid_value(lo0, hi0, q - 1) >= x: q -= 1
    return q


@njit(cache=True)
def contract(fsize, fp0, fp1, fmode, vstart, vface, vpos, lo, hi, rounds, use_bounds):
    """Returns (status, lo, hi, CL, CH, RQL, RQH, MUL, MUH, nr, MU3): status 0 alive, 1 invalid, 2 vsum,
    3 linear exclusion of the returned box with multipliers MU3.
    Each contraction round r < nr maps the box at the start of the round to the grid box with indices
    RQL[r], RQH[r] relative to it; MUL[r, j] / MUH[r, j] are the multipliers that proved the new lower / upper
    bound of coordinate j (hints for the checker).  The returned box is the result of all rounds."""
    n = lo.shape[0]; V = vstart.shape[0] - 1
    RQL = np.zeros((rounds, n), np.int64); RQH = np.full((rounds, n), QMAX, np.int64); nr = 0
    MUL = np.zeros((rounds, n, V)); MUH = np.zeros((rounds, n, V)); MU3 = np.zeros(V)
    qlo = np.zeros(n, np.int64); qhi = np.full(n, QMAX, np.int64)
    lo = lo.copy(); hi = hi.copy()
    ok, vok, CL, CH, FL, FH, JL, JH = eval_all(fsize, fp0, fp1, fmode, vstart, vface, vpos, lo, hi, lo, hi)
    if not ok: return 1, lo, hi, CL, CH, RQL, RQH, MUL, MUH, nr, MU3
    if not vok: return 2, lo, hi, CL, CH, RQL, RQH, MUL, MUH, nr, MU3
    for rnd in range(rounds):
        m = 0.5 * (lo + hi)
        okm, _, _, _, FmL, FmH, _, _ = eval_all(fsize, fp0, fp1, fmode, vstart, vface, vpos, m, m, lo, hi)
        Jc = np.zeros((V, n)); R = np.zeros(V); fc = np.zeros(V)
        mask = np.ones(V, np.bool_)
        for v in range(V):
            good = np.isfinite(FmL[v]) and np.isfinite(FmH[v])
            for j in range(n):
                if not (np.isfinite(JL[v, j]) and np.isfinite(JH[v, j])): good = False
            mask[v] = good
            if good:
                fc[v] = 0.5 * (FmL[v] + FmH[v]); r = 0.5 * (FmH[v] - FmL[v])
                for j in range(n):
                    Jc[v, j] = 0.5 * (JL[v, j] + JH[v, j])
                    r += 0.5 * (JH[v, j] - JL[v, j]) * 0.5 * (hi[j] - lo[j])
                R[v] = r * (1 + 1e-9) + 1e-14
            else:
                R[v] = 1e30
        if not mask.any(): break
        ylo = lo - m; yhi = hi - m
        # zero multipliers on bad rows: mask F and J there to exact zeros for the rigorous combination
        FL2 = FmL.copy(); FH2 = FmH.copy(); JL2 = JL.copy(); JH2 = JH.copy()
        for v in range(V):
            if not mask[v]:
                FL2[v] = 0.0; FH2[v] = 0.0
                for j in range(n): JL2[v, j] = 0.0; JH2[v, j] = 0.0
        T, basis = lp_setup(Jc, fc, R, ylo, yhi - ylo)
        cost = np.zeros(n + 1); cost[n] = 1.0
        lp_set_objective(T, basis, cost)
        st = _simplex(T, basis, 400)
        if st != 0: break
        tstar = -T[T.shape[0] - 1, T.shape[1] - 1]
        if tstar > 1e-12:
            mu = lp_multipliers(T, V, n)
            for v in range(V):
                if not mask[v]: mu[v] = 0.0
            if excluded_by(mu, FL2, FH2, JL2, JH2, ylo, yhi):
                MU3[:] = mu
                return 3, lo, hi, CL, CH, RQL, RQH, MUL, MUH, nr, MU3
            break
        if not use_bounds: break
        w0 = hi - lo
        for k in range(n):
            for sgn in (1.0, -1.0):
                cost = np.zeros(n + 1); cost[k] = sgn; cost[n] = 1e3
                lp_set_objective(T, basis, cost)
                st = _simplex(T, basis, 400)
                if st != 0: continue
                mu = lp_multipliers(T, V, n)
                for v in range(V):
                    if not mask[v]: mu[v] = 0.0
                nl, nh, emp = gs_update(mu, k, FL2, FH2, JL2, JH2, ylo, yhi)
                if emp:
                    # the box is empty; record it as a leaf only if this mu excludes the whole round box
                    if excluded_by(mu, FL2, FH2, JL2, JH2, lo - m, hi - m):
                        MU3[:] = mu
                        return 3, lo, hi, CL, CH, RQL, RQH, MUL, MUH, nr, MU3
                    continue
                if nl > ylo[k]: MUL[nr, k, :] = mu
                if nh < yhi[k]: MUH[nr, k, :] = mu
                ylo[k] = nl; yhi[k] = nh
        newlo = np.empty(n); newhi = np.empty(n); changed = False
        for j in range(n):
            a = max(lo[j], dn(m[j] + ylo[j])); b = min(hi[j], up(m[j] + yhi[j]))
            ql = 0; qh = QMAX
            if hi[j] > lo[j]:
                # leave a margin of 1e-4 of the width (at least 1e-12) between the proved bound and the
                # recorded one, so that an independent check with slightly different rounding still succeeds
                mg = 1e-4 * (hi[j] - lo[j]) + 1e-12
                ql = grid_down(lo[j], hi[j], max(lo[j], a - mg)); qh = grid_up(lo[j], hi[j], min(hi[j], b + mg))
                # keep both indices encodable in 16 bits: ql <= QMAX - 1, qh >= ql + 1 (outward, sound)
                ql = min(ql, QMAX - 1)
                if qh < ql + 1: qh = ql + 1
            RQL[nr, j] = ql; RQH[nr, j] = qh
            if ql > 0 or qh < QMAX: changed = True
            newlo[j] = grid_value(lo[j], hi[j], ql); newhi[j] = grid_value(lo[j], hi[j], qh)
        if not changed: break
        nr += 1
        lo = newlo; hi = newhi
        ok, vok, CL, CH, FL, FH, JL, JH = eval_all(fsize, fp0, fp1, fmode, vstart, vface, vpos, lo, hi, lo, hi)
        if not ok: return 1, lo, hi, CL, CH, RQL, RQH, MUL, MUH, nr, MU3
        if not vok: return 2, lo, hi, CL, CH, RQL, RQH, MUL, MUH, nr, MU3
        shrink = 0.0
        for j in range(n):
            if w0[j] > 0: shrink = max(shrink, (w0[j] - (hi[j] - lo[j])) / w0[j])
        if shrink < 0.2: break
    return 0, lo, hi, CL, CH, RQL, RQH, MUL, MUH, nr, MU3


# ============================================================================ Gordan
@njit(cache=True)
def gordan_vector(fsize, vstart, vface, vpos, eu, eus, ew, ews, Cpt):
    """LP at point corners Cpt: max t s.t. y_u.t_u + y_w.t_w >= t for all edges, |y| <= 1.  Returns y or empty."""
    V = vstart.shape[0] - 1; E = eu.shape[0]
    ang = np.zeros(vface.shape[0])
    for v in range(V):
        a = 0.0
        for s in range(vstart[v], vstart[v + 1]):
            ang[s] = a; a += Cpt[vface[s], vpos[s]]
    # variables yp (2V), ym (2V), t ; rows: E edge rows, 4V bounds y<=1, t <= 1
    nv = 4 * V + 1; M = E + 4 * V + 1; ncol = nv + M
    T = np.zeros((M + 1, ncol + 1)); basis = np.empty(M, np.int64)
    for e in range(E):
        for side in range(2):
            v = eu[e] if side == 0 else ew[e]; s = eus[e] if side == 0 else ews[e]
            c = math.cos(ang[s]); sn = math.sin(ang[s])
            T[e, 2 * v] -= c; T[e, 2 * v + 1] -= sn
            T[e, 2 * V + 2 * v] += c; T[e, 2 * V + 2 * v + 1] += sn
        T[e, 4 * V] = 1.0
    for j in range(4 * V + 1):
        T[E + j, j] = 1.0; T[E + j, ncol] = 1.0
    for i in range(M):
        T[i, nv + i] = 1.0; basis[i] = nv + i
    T[M, 4 * V] = -1.0
    st = _simplex(T, basis, 1000)
    y = np.zeros(2 * V)
    if st != 0: return False, y
    for i in range(M):
        b = basis[i]
        if b < 2 * V: y[b] += T[i, ncol]
        elif b < 4 * V: y[b - 2 * V] -= T[i, ncol]
    t = T[M, ncol]
    return t > 1e-12, y


@njit(cache=True)
def gordan_check(vstart, vface, vpos, eu, eus, ew, ews, CL, CH, y):
    V = vstart.shape[0] - 1; E = eu.shape[0]
    cl = np.zeros(vface.shape[0]); ch = np.zeros(vface.shape[0]); sl = np.zeros(vface.shape[0]); sh = np.zeros(vface.shape[0])
    for v in range(V):
        al, ah = 0.0, 0.0
        for s in range(vstart[v], vstart[v + 1]):
            if s == vstart[v]:
                cl[s] = 1.0; ch[s] = 1.0; sl[s] = 0.0; sh[s] = 0.0
            else:
                cl[s], ch[s] = icos(al, ah); sl[s], sh[s] = isin(al, ah)
            al, ah = iadd(al, ah, CL[vface[s], vpos[s]], CH[vface[s], vpos[s]])
    for e in range(E):
        tl, th = 0.0, 0.0
        for side in range(2):
            v = eu[e] if side == 0 else ew[e]; s = eus[e] if side == 0 else ews[e]
            a, b = imul(cl[s], ch[s], y[2 * v], y[2 * v]); tl, th = iadd(tl, th, a, b)
            a, b = imul(sl[s], sh[s], y[2 * v + 1], y[2 * v + 1]); tl, th = iadd(tl, th, a, b)
        if not (tl > 0.0): return False
    return True


@njit(cache=True)
def motzkin_vector(vstart, vface, vpos, eu, eus, ew, ews, Cpt):
    """LP at point corners: max sum_e s_e, s_e <= a_e.y, 0 <= s_e <= 1, |y| <= 1.  Small components of y are then
    set to exactly 0 so that edges away from the support of y have exactly zero terms (Motzkin certificate
    against strictly positive stresses)."""
    V = vstart.shape[0] - 1; E = eu.shape[0]
    ang = np.zeros(vface.shape[0])
    for v in range(V):
        a = 0.0
        for s in range(vstart[v], vstart[v + 1]):
            ang[s] = a; a += Cpt[vface[s], vpos[s]]
    nv = 4 * V + E; M = E + nv; ncol = nv + M
    T = np.zeros((M + 1, ncol + 1)); basis = np.empty(M, np.int64)
    for e in range(E):
        for side in range(2):
            v = eu[e] if side == 0 else ew[e]; s = eus[e] if side == 0 else ews[e]
            c = math.cos(ang[s]); sn = math.sin(ang[s])
            T[e, 2 * v] -= c; T[e, 2 * v + 1] -= sn
            T[e, 2 * V + 2 * v] += c; T[e, 2 * V + 2 * v + 1] += sn
        T[e, 4 * V + e] = 1.0
    for j in range(nv):
        T[E + j, j] = 1.0; T[E + j, ncol] = 1.0
    for i in range(M):
        T[i, nv + i] = 1.0; basis[i] = nv + i
    for e in range(E): T[M, 4 * V + e] = -1.0
    st = _simplex(T, basis, 2000)
    y = np.zeros(2 * V)
    if st != 0: return False, y
    for i in range(M):
        b = basis[i]
        if b < 2 * V: y[b] += T[i, ncol]
        elif b < 4 * V: y[b - 2 * V] -= T[i, ncol]
    for v in range(V):
        if abs(y[2 * v]) + abs(y[2 * v + 1]) < 1e-7: y[2 * v] = 0.0; y[2 * v + 1] = 0.0
    return T[M, ncol] > 1e-9, y


@njit(cache=True)
def motzkin_check(vstart, vface, vpos, eu, eus, ew, ews, CL, CH, y):
    """True if every edge term y_u.t_u + y_w.t_w is exactly 0 (both y_u, y_w = 0) or > 0 on the whole box,
    and at least one is > 0.  Then no strictly positive balanced stress exists on the box."""
    V = vstart.shape[0] - 1; E = eu.shape[0]
    cl = np.zeros(vface.shape[0]); ch = np.zeros(vface.shape[0]); sl = np.zeros(vface.shape[0]); sh = np.zeros(vface.shape[0])
    for v in range(V):
        al, ah = 0.0, 0.0
        for s in range(vstart[v], vstart[v + 1]):
            if s == vstart[v]:
                cl[s] = 1.0; ch[s] = 1.0; sl[s] = 0.0; sh[s] = 0.0
            else:
                cl[s], ch[s] = icos(al, ah); sl[s], sh[s] = isin(al, ah)
            al, ah = iadd(al, ah, CL[vface[s], vpos[s]], CH[vface[s], vpos[s]])
    npos = 0
    for e in range(E):
        u = eu[e]; w = ew[e]
        zu = y[2 * u] == 0.0 and y[2 * u + 1] == 0.0
        zw = y[2 * w] == 0.0 and y[2 * w + 1] == 0.0
        if zu and zw: continue
        tl, th = 0.0, 0.0
        for side in range(2):
            v = u if side == 0 else w; s = eus[e] if side == 0 else ews[e]
            if y[2 * v] == 0.0 and y[2 * v + 1] == 0.0: continue
            a, b = imul(cl[s], ch[s], y[2 * v], y[2 * v]); tl, th = iadd(tl, th, a, b)
            a, b = imul(sl[s], sh[s], y[2 * v + 1], y[2 * v + 1]); tl, th = iadd(tl, th, a, b)
        if not (tl > 0.0): return False
        npos += 1
    return npos > 0


# ============================================================================ driver
# Certificate events (see docs/certificate-format.md), one byte each, in depth-first preorder:
#   1 invalid  2 vertex-sum  3 linear  4 gordan (inherited vector)  5 gordan  6 motzkin  7 terminal  8 open
#   16 + k     bisect coordinate k at the midpoint; the UPPER child's subtree follows, then the LOWER child's
#   15         start of a contraction round; the 64 + k events that follow are relative to the box at this point
#   64 + k     contract coordinate k to grid indices (qlo, QMAX - qhi), two uint16 in the q stream
EV_ROUND = 15
EV_SPLIT = 16
EV_CONTRACT = 64


@njit(cache=True)
def _push_ev(ev, nev, code):
    if nev == ev.shape[0]:
        e2 = np.empty(2 * nev, np.uint8); e2[:nev] = ev; ev = e2
    ev[nev] = code
    return ev, nev + 1


@njit(cache=True)
def _push_q(qs, nq, a, b):
    if nq + 2 > qs.shape[0]:
        q2 = np.empty(2 * qs.shape[0] + 2, np.uint16); q2[:nq] = qs[:nq]; qs = q2
    qs[nq] = a; qs[nq + 1] = b
    return qs, nq + 2


@njit(cache=True)
def _push_mu(qs, nq, ms, nm, mu):
    """Sparse multiplier hint: count and row indices in the q stream, values (float32) in the m stream.
    The vector is scaled to max |mu_v| = 1 and entries below 1e-9 are dropped; any real vector is admissible
    in the checker's inclusion, so this changes nothing about soundness."""
    mx = 0.0
    for v in range(mu.shape[0]): mx = max(mx, abs(mu[v]))
    cnt = 0
    if mx > 0:
        for v in range(mu.shape[0]):
            if abs(mu[v]) > 1e-9 * mx: cnt += 1
    if nq + cnt + 1 > qs.shape[0]:
        q2 = np.empty(2 * qs.shape[0] + cnt + 1, np.uint16); q2[:nq] = qs[:nq]; qs = q2
    if nm + cnt > ms.shape[0]:
        m2 = np.empty(2 * ms.shape[0] + cnt, np.float32); m2[:nm] = ms[:nm]; ms = m2
    qs[nq] = cnt; nq += 1
    if mx > 0:
        for v in range(mu.shape[0]):
            if abs(mu[v]) > 1e-9 * mx:
                qs[nq] = v; nq += 1
                ms[nm] = mu[v] / mx; nm += 1
    return qs, nq, ms, nm


@njit(cache=True)
def run_chart(fsize, fp0, fp1, fmode, vstart, vface, vpos, eu, eus, ew, ews, lo0, hi0, tol, max_boxes,
              w_lp, w_gordan, rounds, use_bounds):
    """Depth-first branch and bound over one chart.
    Returns (stats, terminal_lo, terminal_hi, aborted, events, qstream, mustream)."""
    n = lo0.shape[0]; V = vstart.shape[0] - 1
    cap = 1024
    SL = np.empty((cap, n)); SH = np.empty((cap, n)); SY = np.zeros((cap, 2 * V)); SHY = np.zeros(cap, np.bool_)
    SL[0] = lo0; SH[0] = hi0; top = 1
    TL = np.empty((64, n)); TH = np.empty((64, n)); nt = 0
    ev = np.empty(4096, np.uint8); nev = 0
    qs = np.empty(4096, np.uint16); nq = 0
    ms = np.empty(1024, np.float32); nm = 0
    stats = np.zeros(8, np.int64)   # boxes, invalid, vsum, linear, gordan_inherit, gordan, terminal, lpcalls
    aborted = False
    while top > 0:
        top -= 1
        if aborted:                       # budget exhausted: the rest of the stack is recorded as open
            ev, nev = _push_ev(ev, nev, 8); continue
        lo = SL[top].copy(); hi = SH[top].copy(); y = SY[top].copy(); hy = SHY[top]
        stats[0] += 1
        if stats[0] > max_boxes:
            aborted = True
            ev, nev = _push_ev(ev, nev, 8); continue
        wmax = 0.0
        for j in range(n): wmax = max(wmax, hi[j] - lo[j])
        if wmax <= w_lp:
            # near terminal size, only the exclusion test is used (no bound tightening): slabs cut from tiny
            # boxes next to a critical point would have no margin left for an independent check
            ub = use_bounds and wmax > 1000.0 * tol
            status, lo2, hi2, CL, CH, RQL, RQH, MUL, MUH, nr, MU3 = contract(fsize, fp0, fp1, fmode, vstart,
                                                                             vface, vpos, lo, hi, rounds, ub)
            for r in range(nr):
                ev, nev = _push_ev(ev, nev, EV_ROUND)
                for j in range(n):
                    if RQL[r, j] > 0 or RQH[r, j] < QMAX:
                        ev, nev = _push_ev(ev, nev, EV_CONTRACT + j)
                        qs, nq = _push_q(qs, nq, RQL[r, j], QMAX - RQH[r, j])
                        if RQL[r, j] > 0:
                            qs, nq, ms, nm = _push_mu(qs, nq, ms, nm, MUL[r, j])
                        if RQH[r, j] < QMAX:
                            qs, nq, ms, nm = _push_mu(qs, nq, ms, nm, MUH[r, j])
            lo = lo2; hi = hi2
        else:
            ok, vok, CL, CH, FL, FH, JL, JH = eval_all(fsize, fp0, fp1, fmode, vstart, vface, vpos, lo, hi, lo, hi)
            status = 0 if (ok and vok) else (1 if not ok else 2)
        if status > 0:
            stats[status] += 1; ev, nev = _push_ev(ev, nev, status)
            if status == 3:
                qs, nq, ms, nm = _push_mu(qs, nq, ms, nm, MU3)
            continue
        if hy and gordan_check(vstart, vface, vpos, eu, eus, ew, ews, CL, CH, y):
            stats[4] += 1; ev, nev = _push_ev(ev, nev, 4); continue
        wmax = 0.0
        for j in range(n): wmax = max(wmax, hi[j] - lo[j])
        if wmax <= w_gordan:
            m = 0.5 * (lo + hi)
            okm, _, PL, PH, _, _, _, _ = eval_all(fsize, fp0, fp1, fmode, vstart, vface, vpos, m, m, lo, hi)
            Cpt = 0.5 * (PL + PH)
            stats[7] += 1
            found, y2 = gordan_vector(fsize, vstart, vface, vpos, eu, eus, ew, ews, Cpt)
            if found:
                y = y2; hy = True
                if gordan_check(vstart, vface, vpos, eu, eus, ew, ews, CL, CH, y):
                    stats[5] += 1; ev, nev = _push_ev(ev, nev, 5); continue
            found, y3 = motzkin_vector(vstart, vface, vpos, eu, eus, ew, ews, Cpt)
            if found and motzkin_check(vstart, vface, vpos, eu, eus, ew, ews, CL, CH, y3):
                stats[5] += 1; ev, nev = _push_ev(ev, nev, 6); continue
        k = 0
        for j in range(n):
            if hi[j] - lo[j] > hi[k] - lo[k]: k = j
        if hi[k] - lo[k] < tol:
            if nt == TL.shape[0]:
                TL2 = np.empty((2 * nt, n)); TH2 = np.empty((2 * nt, n)); TL2[:nt] = TL; TH2[:nt] = TH; TL = TL2; TH = TH2
            TL[nt] = lo; TH[nt] = hi; nt += 1; stats[6] += 1
            ev, nev = _push_ev(ev, nev, 7)
            continue
        if top + 2 > cap:
            cap *= 2
            SL2 = np.empty((cap, n)); SH2 = np.empty((cap, n)); SY2 = np.zeros((cap, 2 * V)); SHY2 = np.zeros(cap, np.bool_)
            SL2[:top] = SL[:top]; SH2[:top] = SH[:top]; SY2[:top] = SY[:top]; SHY2[:top] = SHY[:top]
            SL = SL2; SH = SH2; SY = SY2; SHY = SHY2
        ev, nev = _push_ev(ev, nev, EV_SPLIT + k)
        mid = 0.5 * (lo[k] + hi[k])
        # pushed lower first, so the upper child is popped (and logged) first
        SL[top] = lo; SH[top] = hi; SH[top, k] = mid; SY[top] = y; SHY[top] = hy; top += 1
        SL[top] = lo; SH[top] = hi; SL[top, k] = mid; SY[top] = y; SHY[top] = hy; top += 1
    return stats, TL[:nt], TH[:nt], aborted, ev[:nev].copy(), qs[:nq].copy(), ms[:nm].copy()
