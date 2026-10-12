"""Batched interval evaluation of corner angles, vertex angle sums and their Jacobians.

A batch is B boxes in face-parameter space, given as arrays lo, hi of shape (npar, B).  Everything is computed
by forward-mode automatic differentiation over intervals (class D: value I of shape (B,), gradient I of
shape (npar, B)).  Division by an interval containing 0 gives the whole line (interval.STRICT_DIV = False), so
one bad box never stops a batch; such boxes simply get infinite enclosures and are split later.

Pentagon charts.  The face-angle chart (f1, f2) of an equilateral pentagon (corners at boundary positions
1, 2) is singular where the corner at position 4 (the apex of the last isosceles triangle) is straight.  Each
pentagon is treated in one of three regions, each with its own chart (formula position j is boundary
position (s + j) mod 5; c_q = corner at boundary position q; D = DELTA):
  'a':  c4 <= pi - D;                 s = 0, apex c4 clamped to <= pi - D;
  'b':  c4 >= pi - D, c1 <= pi - D;   s = 2, parameters c3, c4 (c4 in [pi - D, pi]), apex c1 <= pi - D;
  'c':  c4 >= pi - D, c1 >= pi - D;   s = 3, parameters c4, c0 (c4 in [pi - D, pi]), apex c2 <= pi - D,
                                       and boxes with c1 < pi - D everywhere are rejected.
The regions cover every pentagon, given Lemma P (no two adjacent corners are >= pi - D), which is checked
separately (scripts/completeness/lemmas.py P).  Every chart clamps its apex to <= pi - D, so all
formulas are smooth on their regions.

Value clipping.  Where a quantity is known to lie in a range on every valid point of the box (for instance the
half-apex angle h in [alpha/2, min(f1, f2)/2]), its value enclosure is intersected with that range and
derivatives that depend on it are evaluated on the clipped value.  This is sound: a box point violating the
range is not a valid configuration in this chart, so nothing is lost by ignoring it.
"""
import numpy as np
from . import interval as iv
from .interval import I

iv.STRICT_DIV = False
np.seterr(all='ignore')
PI = np.pi
DELTA = np.radians(25.0)


class D:
    __slots__ = ('v', 'g')

    def __init__(self, v, g): self.v, self.g = v, g

    def __add__(self, o):
        if isinstance(o, D): return D(self.v + o.v, self.g + o.g)
        return D(self.v + o, self.g)
    __radd__ = __add__

    def __sub__(self, o):
        if isinstance(o, D): return D(self.v - o.v, self.g - o.g)
        return D(self.v - o, self.g)

    def __rsub__(self, o): return D(o - self.v, -self.g)
    def __neg__(self): return D(-self.v, -self.g)

    def __mul__(self, o):
        if isinstance(o, D): return D(self.v * o.v, self.g * o.v + o.g * self.v)
        return D(self.v * o, self.g * o)
    __rmul__ = __mul__

    def __truediv__(self, o):
        if isinstance(o, D):
            q = self.v / o.v
            return D(q, (self.g - o.g * q) / o.v)
        return D(self.v / o, self.g / o)


def _chain(x, val, dval): return D(val, x.g * dval)
def cos(x): return _chain(x, iv.cos(x.v), -iv.sin(x.v))
def sin(x): return _chain(x, iv.sin(x.v), iv.cos(x.v))
def atan(x): return _chain(x, iv.arctan(x.v), 1.0 / (1.0 + x.v * x.v))


def _clamp_factor(active):
    """Generalized derivative of a clamp: 1 where inactive on the whole box, [0, 1] where it may be active."""
    return I(np.where(active, 0.0, 1.0), np.ones_like(active, float))


def clip(x, lo, hi):
    """Clamp to [lo, hi]: value = image of the clamp, gradient times the clamp's generalized derivative.
    Returns (x', mask where the unclamped range meets [lo, hi])."""
    ok = (x.v.hi >= lo) & (x.v.lo <= hi)
    active = (x.v.lo < lo) | (x.v.hi > hi)
    v = I(np.clip(x.v.lo, lo, hi), np.clip(x.v.hi, lo, hi))
    return D(v, x.g * _clamp_factor(active)), ok


def asin_angle(x, lo, hi):
    """y = clamp(arcsin(clamp(x, -1, 1)), lo, hi) with [lo, hi] inside [-pi/2, pi/2]; dy/dx = 1/cos(y) times
    clamp factors.  Returns (y, mask where some point of the box has |x| <= 1 and arcsin(x) in [lo, hi])."""
    ok = (x.v.lo <= 1) & (x.v.hi >= -1)
    y0 = iv.arcsin(x.v)
    ok &= (y0.hi >= lo) & (y0.lo <= hi)
    active = (x.v.lo < -1) | (x.v.hi > 1) | (y0.lo < lo) | (y0.hi > hi)
    y = I(np.clip(y0.lo, lo, hi), np.clip(y0.hi, lo, hi))
    return D(y, x.g / iv.cos(y) * _clamp_factor(active)), ok


def acos_val(x):
    ok = (x.v.lo <= 1) & (x.v.hi >= -1)
    y = iv.arccos(x.v)
    return D(y, -x.g / iv.sin(y)), ok


# ---------------------------------------------------------------------- faces
def alpha_value(th_lo, th_hi):
    """Exact monotone enclosure of alpha(theta) = arccos(cos t / (1 + cos t)), increasing in t."""
    f = lambda t: np.arccos(np.cos(t) / (1 + np.cos(t)))
    return I(iv._dn(f(th_lo), 4), iv._up(f(th_hi), 4))


def beta(ct, A):
    """Base angle of the isosceles triangle with legs theta and apex A:  atan(cos(A/2) / (cos t sin(A/2)))."""
    return atan(cos(A * 0.5) / (ct * sin(A * 0.5)))


class Model:
    """Corner functions for a Structure S with pentagon chart offsets `charts` (dict face -> s)."""

    def __init__(self, S, charts=None):
        self.S = S; self.n = S.npar
        self.charts = dict(charts or {})
        for fi, k in enumerate(S.sizes):
            if k == 5: self.charts.setdefault(fi, 'a')

    def variables(self, lo, hi):
        n, B = lo.shape
        P = []
        for k in range(n):
            g = np.zeros((n, B)); g[k] = 1.0
            P.append(D(I(lo[k], hi[k]), I(g, g.copy())))
        return P

    def corners(self, lo, hi, box=None):
        """Per face: tuple of D corner values in boundary order, and a validity mask (B,).
        The clamp bounds are taken from `box` = (blo, bhi) (default: the evaluated box itself), so that a
        point evaluation inside a box uses the same clamped function as the box evaluation."""
        S = self.S; P = self.variables(lo, hi); B = lo.shape[1]
        blo, bhi = (lo, hi) if box is None else box
        th = P[0]; ct, st = cos(th), sin(th)
        al_v = alpha_value(lo[0], hi[0])
        c = ct; al = acos_val(c / (1 + c))[0]; al = D(al_v, al.g)       # exact value, AD derivative
        ok = np.ones(B, bool)
        out = []
        for fi, k in enumerate(S.sizes):
            if k == 3:
                out.append((al, al, al)); continue
            if k == 4:
                A = P[S.pidx[fi][0]]
                O = beta(ct, A) * 2.0
                # exact range: 2 beta increasing in theta, decreasing in A
                ex_lo = 2 * np.arctan2(np.cos(hi[S.pidx[fi][0]] / 2), np.cos(lo[0]) * np.sin(hi[S.pidx[fi][0]] / 2))
                ex_hi = 2 * np.arctan2(np.cos(lo[S.pidx[fi][0]] / 2), np.cos(hi[0]) * np.sin(lo[S.pidx[fi][0]] / 2))
                O = D(I(np.maximum(O.v.lo, iv._dn(ex_lo, 4)), np.minimum(O.v.hi, iv._up(ex_hi, 4))), O.g)
                out.append((A, O, A, O)); continue
            i1, i2 = S.pidx[fi]; f1, f2 = P[i1], P[i2]
            hmin = 0.5 * alpha_value(blo[0], blo[0]).lo
            mode = self.charts[fi]
            hmax = 0.5 * (PI - DELTA)
            cs, o = self.pentagon(th, ct, st, hmin, hmax, f1, f2)
            ok &= o
            s = {'a': 0, 'b': 2, 'c': 3}[mode]
            if mode == 'c': ok &= cs[3].v.hi >= PI - DELTA          # c1 = formula position 3
            out.append(tuple(cs[(q - s) % 5] for q in range(5)))
        return out, ok

    @staticmethod
    def pentagon(th, ct, st, hmin, hmax, f1, f2):
        """Corners (formula order 0..4) of the equilateral pentagon P0..P4 with corners f1, f2 at P1, P2,
        with the apex corner at P4 clamped to [2 hmin, 2 hmax]."""
        hd2, ok = asin_angle(st * sin(f1 * 0.5), 0.0, PI / 2)
        d2 = hd2 * 2.0
        b1 = beta(ct, f1)
        g, o = clip(f2 - b1, 0.0, PI); ok &= o
        cg, sg = cos(g), sin(g)
        cd2, sd2 = cos(d2), sin(d2)
        cd3 = cd2 * ct + sd2 * st * cg
        ok &= cd3.v.lo <= ct.v.hi                      # diagonal P0 P3 >= theta
        cd3, _ = clip(cd3, -1.0, 1.0)
        d3, o = acos_val(cd3); ok &= o
        a0 = PI / 2 - atan((sd2 * ct - cd2 * st * cg) / (sg * st))
        a3 = PI / 2 - atan((st * cd2 - ct * sd2 * cg) / (sg * sd2))
        h, o = asin_angle(sin(d3 * 0.5) / st, hmin, hmax); ok &= o
        f4 = h * 2.0
        b3 = beta(ct, f4)
        c0 = b1 + a0 + b3; c3 = a3 + b3
        return (c0, f1, f2, c3, f4), ok

    # ------------------------------------------------------------------ vertex sums
    def vertex_sums(self, C):
        S = self.S; out = []
        for v in range(S.n):
            s = None
            for (_, fi, pos) in S.rot[v]:
                s = C[fi][pos] if s is None else s + C[fi][pos]
            out.append(s)
        return out

    def initial_box(self, th_lo, th_hi):
        S = self.S
        a_lo = float(alpha_value(np.array(th_lo), np.array(th_lo)).lo)
        a_hi = float(alpha_value(np.array(th_hi), np.array(th_hi)).hi)
        lo = np.zeros(S.npar); hi = np.zeros(S.npar)
        lo[0], hi[0] = th_lo, th_hi
        for fi, k in enumerate(S.sizes):
            for i in S.pidx[fi]:
                if k == 4: lo[i], hi[i] = a_lo, 2 * a_hi
                else: lo[i], hi[i] = a_lo, PI
            if k == 5 and self.charts[fi] == 'b':
                lo[S.pidx[fi][1]] = PI - DELTA          # formula position 2 = boundary position 4
            if k == 5 and self.charts[fi] == 'c':
                lo[S.pidx[fi][0]] = PI - DELTA          # formula position 1 = boundary position 4
        return lo, hi
