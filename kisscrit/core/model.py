"""Face-angle model of a candidate stressed graph, for the checker: inf-sup intervals with rounding through Arb
(kisscrit/core/interval.py).

This is the checker's own implementation of the model.  It shares no numerical code with the search
(kisscrit/completeness/fast.py); the formulas are the same, and they are restated here so that this file can be
read on its own.

Conventions.  A face is a convex spherical polygon with all sides theta.  At a vertex v of a face listed
counterclockwise (seen from outside) as (..., a, v, b, ...), the corner is the angle at v from the edge v->b
counterclockwise to the edge v->a.  The parameter vector p is:
  p[0] = theta;
  one parameter per quadrilateral: the corner at boundary positions 0 and 2 (the others are 2 beta(theta, A));
  two per pentagon: the corners at formula positions 1 and 2 of its chart (see `pentagon`).
Pentagon charts ('a', 'b', 'c') rotate the formula positions: formula position j is boundary position
(s + j) mod 5 with s = 0, 2, 3, and the apex corner (formula position 4) is clamped to at most LSTRAIGHT.

Clamped functions.  Where the formulas need a quantity in a range that every valid configuration satisfies
(corners at least alpha(theta), the apex corner at most LSTRAIGHT, the angle g in [0, pi]), the value is clamped
to that range and derivatives are multiplied by [0, 1] where the clamp may act.  The clamped function agrees
with the true corner function at every valid point and is Lipschitz, so its interval Jacobian encloses the
Clarke generalized Jacobian and the mean value inclusion F(p) in F(m) + J(B)(p - m) holds for p, m in a box B.

Every function returns balls; `ok` is False only when the box certainly contains no valid point.
"""
import math
from flint import arb
from . import interval as I
from .interval import Iv, PI, HALF_PI, TWO_PI, down, up

ZERO_ONE = Iv(0.0, 1.0)
LSTRAIGHT = math.pi - math.radians(25.0)      # the same double as in the search: the region threshold


def hull(lo, hi): return Iv(lo, hi)


def intersect(x, y):
    lo = max(x.lo, y.lo); hi = min(x.hi, y.hi)
    return x if lo > hi else Iv(lo, hi)


def lower(x): return x.lo
def upper(x): return x.hi
def finite(x): return x.finite()


# ------------------------------------------------------------------ forward AD with 3 local components
# a value is (v, g0, g1, g2): value and gradient with respect to theta and the face's two parameters

def var(lo, hi, k):
    g = [Z, Z, Z]; g[k] = Iv(1.0)
    return (Iv(lo, hi), g[0], g[1], g[2])


def m0(x, y): return x * y          # Iv keeps exact zeros exact
def d0(x, y): return x / y


Z = Iv(0.0)


def const(c): return (c if isinstance(c, Iv) else Iv(c), Z, Z, Z)
def add(a, b): return tuple(x + y for x, y in zip(a, b))
def sub(a, b): return tuple(x - y for x, y in zip(a, b))
def neg(a): return tuple(-x for x in a)
def scale(a, c): return tuple(x * Iv(c) for x in a)
def shift(a, c): return (a[0] + Iv(c), a[1], a[2], a[3])


def mul(a, b):
    return (a[0] * b[0],) + tuple(_sum0(m0(a[i], b[0]), m0(b[i], a[0])) for i in (1, 2, 3))


def _sum0(x, y):
    if x.is_zero(): return y
    if y.is_zero(): return x
    return x + y


def div(a, b):
    q = a[0] / b[0]
    return (q,) + tuple(d0(_sum0(a[i], -m0(b[i], q)), b[0]) for i in (1, 2, 3))


def chain(a, v, d): return (v, m0(a[1], d), m0(a[2], d), m0(a[3], d))
def cos_(a): return chain(a, I.cos(a[0]), -I.sin(a[0]))
def sin_(a): return chain(a, I.sin(a[0]), I.cos(a[0]))


def atan_(a):
    x = a[0]
    if not finite(x):
        return chain(a, I.atan(x), Iv(0.0, 1.0))
    return chain(a, I.atan(x), Iv(1.0) / (Iv(1.0) + I.square(x)))


def clip(a, lo, hi):
    """Clamp to [lo, hi].  Returns (value, possibly_nonempty): the value is the image of the clamp; the
    gradient is multiplied by [0, 1] unless the clamp certainly does not act on the box."""
    x = a[0]
    ok = not (x.hi < lo) and not (x.lo > hi)
    inside = x.lo > lo and x.hi < hi
    v = I.clamp(x, lo, hi)
    if inside: return (v, a[1], a[2], a[3]), ok
    return (v, a[1] * ZERO_ONE, a[2] * ZERO_ONE, a[3] * ZERO_ONE), ok


def asin_angle(a, lo, hi):
    """y = clamp(asin(clamp(x, -1, 1)), lo, hi) with [lo, hi] inside [-pi/2, pi/2); dy/dx = 1/cos(y) times
    clamp factors.  Returns (y, possibly_valid)."""
    x = a[0]
    ok = not (x.hi < -1) and not (x.lo > 1)
    y0 = I.asin(x)
    ok = ok and not (y0.hi < lo) and not (y0.lo > hi)
    inside = x.lo > -1 and x.hi < 1 and y0.lo > lo and y0.hi < hi
    y = I.clamp(y0, lo, hi)
    d = Iv(1.0) / I.cos(y)
    if not inside: d = d * ZERO_ONE
    return (y, a[1] * d, a[2] * d, a[3] * d), ok


def acos_(a):
    x = a[0]
    ok = not (x.hi < -1) and not (x.lo > 1)
    y = I.acos(x)
    return chain(a, y, Iv(-1.0) / I.sin(y)), ok


# ------------------------------------------------------------------ faces
def _alpha_arb(t):
    c = arb(t).cos()
    return (c / (1 + c)).acos()


def alpha_range(tlo, thi):
    """alpha(theta) = arccos(cos t / (1 + cos t)) is increasing on (0, pi/2): exact range from the ends."""
    return Iv(down(_alpha_arb(tlo)), up(_alpha_arb(thi)))


def beta(ct, A):
    """Base angle of the isosceles spherical triangle with legs theta and apex angle A."""
    h = scale(A, 0.5)
    return atan_(div(cos_(h), mul(ct, sin_(h))))


def pentagon(ct, st, hmin, hmax, f1, f2):
    """Corners (formula positions 0..4) of the convex equilateral pentagon P0..P4 with corners f1, f2 at P1, P2.
    P0 P2 is the base of the isosceles triangle P0 P1 P2; P3 follows from the angle g at P2 between P2 P0 and
    P2 P3; P4 is the apex of the isosceles triangle on the diagonal P0 P3, on the convex side."""
    hd2, ok = asin_angle(mul(st, sin_(scale(f1, 0.5))), 0.0, HALF_PI.lo)
    d2 = scale(hd2, 2)
    b1 = beta(ct, f1)
    g, o = clip(sub(f2, b1), 0.0, PI.lo); ok = ok and o
    cg, sg = cos_(g), sin_(g)
    cd2, sd2 = cos_(d2), sin_(d2)
    cd3 = add(mul(cd2, ct), mul(mul(sd2, st), cg))
    ok = ok and not (lower(cd3[0]) > upper(ct[0]))          # diagonal P0 P3 >= theta somewhere
    cd3, _ = clip(cd3, -1.0, 1.0)
    d3, o = acos_(cd3); ok = ok and o
    a0 = sub(const(HALF_PI), atan_(div(sub(mul(sd2, ct), mul(mul(cd2, st), cg)), mul(sg, st))))
    a3 = sub(const(HALF_PI), atan_(div(sub(mul(st, cd2), mul(mul(ct, sd2), cg)), mul(sg, sd2))))
    h, o = asin_angle(div(sin_(scale(d3, 0.5)), st), hmin, hmax); ok = ok and o
    f4 = scale(h, 2)
    b3 = beta(ct, f4)
    return (add(add(b1, a0), b3), f1, f2, add(a3, b3), f4), ok


def _pick(g, lo, hi, want_min):
    """Endpoint (or the whole interval) of a variable where a function with derivative enclosure g attains
    its minimum (want_min) or maximum."""
    if g.lo >= 0: return (lo, lo) if want_min else (hi, hi)
    if g.hi <= 0: return (hi, hi) if want_min else (lo, lo)
    return (lo, hi)


def _tighten(c, j, tlo, thi, alo, ahi, blo, bhi, hmin, hmax):
    """Monotonicity tightening of the computed pentagon corner j: where a partial derivative enclosure has a
    sign on the box, the corner function is monotone in that variable, so its extreme values are attained on
    the corresponding face of the box; evaluate there (with the same clamp bounds)."""
    lo, hi = c[0].lo, c[0].hi
    for want_min in (True, False):
        a = _pick(c[1], tlo, thi, want_min); b = _pick(c[2], alo, ahi, want_min); d = _pick(c[3], blo, bhi, want_min)
        if a == (tlo, thi) and b == (alo, ahi) and d == (blo, bhi): continue
        th = var(a[0], a[1], 0); ct, st = cos_(th), sin_(th)
        cs, _ = pentagon(ct, st, hmin, hmax, var(b[0], b[1], 1), var(d[0], d[1], 2))
        v = cs[j][0]
        if want_min: lo = max(lo, v.lo)
        else: hi = min(hi, v.hi)
    if lo > hi: return c
    return (Iv(lo, hi), c[1], c[2], c[3])


SHIFT = {'a': 0, 'b': 2, 'c': 3}


class Graph:
    """A candidate: oriented faces, pentagon modes, parameter indices, rotation systems, edges."""

    def __init__(self, faces, modes):
        self.faces = [list(f) for f in faces]
        self.sizes = [len(f) for f in self.faces]
        self.modes = {int(k): v for k, v in modes.items()}
        self.pidx = []; n = 1
        for k in self.sizes:
            if k == 3: self.pidx.append(())
            elif k == 4: self.pidx.append((n,)); n += 1
            else: self.pidx.append((n, n + 1)); n += 2
        self.npar = n
        nv = 1 + max(v for f in self.faces for v in f)
        succ = [dict() for _ in range(nv)]
        for fi, f in enumerate(self.faces):
            k = len(f)
            for pos in range(k):
                v, a, b = f[pos], f[pos - 1], f[(pos + 1) % k]
                if b in succ[v]: raise ValueError('inconsistent orientation')
                succ[v][b] = (a, fi, pos)
        self.rot = []
        for v in range(nv):
            start = min(succ[v]); seq = []; w = start
            for _ in range(len(succ[v])):
                a, fi, pos = succ[v][w]; seq.append((w, fi, pos)); w = a
            if w != start: raise ValueError('rotation does not close')
            self.rot.append(seq)
        self.n = nv
        E = set()
        for f in self.faces:
            for i in range(len(f)): E.add(tuple(sorted((f[i], f[(i + 1) % len(f)]))))
        self.edges = sorted(E)
        if self.n - len(self.edges) + len(self.faces) != 2: raise ValueError('not a sphere')

    # -------------------------------------------------------------- evaluation
    def corners(self, lo, hi, blo=None, bhi=None):
        """Corner values with local gradients on the box [lo, hi]; clamp bounds from the box [blo, bhi].
        Returns (C, ok): C[f][q] = (value, g_theta, g_p1, g_p2) in boundary order."""
        blo = lo if blo is None else blo
        th = var(lo[0], hi[0], 0)
        ct, st = cos_(th), sin_(th)
        al_v = alpha_range(lo[0], hi[0])
        al = (al_v,) + acos_(div(ct, shift(ct, 1)))[0][1:]
        hmin = 0.5 * down(_alpha_arb(blo[0]))
        hmax = LSTRAIGHT / 2
        ok = True; C = []
        for fi, k in enumerate(self.sizes):
            if k == 3:
                C.append([al, al, al]); continue
            if k == 4:
                i = self.pidx[fi][0]
                A = var(lo[i], hi[i], 1)
                O = scale(beta(ct, A), 2)
                # 2 beta is increasing in theta and decreasing in A: exact range from two corners of the box
                ex = Iv(down(2 * self._beta_pt(lo[0], hi[i])), up(2 * self._beta_pt(hi[0], lo[i])))
                Ov = intersect(ex, O[0])
                O = (Ov, O[1], O[2], O[3])
                C.append([A, O, A, O]); continue
            i1, i2 = self.pidx[fi]
            f1 = var(lo[i1], hi[i1], 1); f2 = var(lo[i2], hi[i2], 2)
            cs, o = pentagon(ct, st, hmin, hmax, f1, f2); ok = ok and o
            if max(hi[0] - lo[0], hi[i1] - lo[i1], hi[i2] - lo[i2]) > 0.01:
                cs = tuple(_tighten(cs[j], j, lo[0], hi[0], lo[i1], hi[i1], lo[i2], hi[i2], hmin, hmax)
                           if j in (0, 3, 4) else cs[j] for j in range(5))
            mode = self.modes[fi]; s = SHIFT[mode]
            if mode == 'c': ok = ok and not (upper(cs[3][0]) < LSTRAIGHT)
            row = [None] * 5
            for j in range(5): row[(s + j) % 5] = cs[j]
            C.append(row)
        # every corner lies in [alpha(theta), pi]
        for f in C:
            for c in f:
                if c[0].lo > PI.hi or c[0].hi < al_v.lo: ok = False
        return C, ok

    @staticmethod
    def _beta_pt(t, A):
        t = arb(t); A = arb(A)
        return ((A / 2).cos() / (t.cos() * (A / 2).sin())).atan()

    def vertex_sums(self, C):
        """For each vertex: (sum - 2 pi, gradient over the global parameters)."""
        out = []
        for v in range(self.n):
            s = -TWO_PI; g = [Z] * self.npar
            for (w, fi, pos) in self.rot[v]:
                c = C[fi][pos]; s = s + c[0]
                g[0] = g[0] + c[1]
                if len(self.pidx[fi]) >= 1: g[self.pidx[fi][0]] = g[self.pidx[fi][0]] + c[2]
                if len(self.pidx[fi]) == 2: g[self.pidx[fi][1]] = g[self.pidx[fi][1]] + c[3]
            out.append((s, g))
        return out

    def directions(self, C):
        """Direction angle of each edge v->w in the frame of v (first neighbour at angle 0): dict (v, w) -> ball."""
        d = {}
        for v in range(self.n):
            acc = Z
            for (w, fi, pos) in self.rot[v]:
                d[(v, w)] = acc; acc = acc + C[fi][pos][0]
        return d
