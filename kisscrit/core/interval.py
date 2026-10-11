"""Inf-sup interval arithmetic with float endpoints and rigorous rounding through Arb (python-flint).

An interval is a pair of doubles [lo, hi] (possibly infinite).  The basic operations (+, -, *, /) use IEEE 754
double arithmetic, which is correctly rounded, and step one ulp outward.  Elementary functions (sin, cos, atan,
asin, acos) are evaluated as Arb balls at 64-bit precision and rounded outward to doubles (`down`, `up`).  So
every result contains the exact real result.  Inf-sup intervals are used (rather than Arb's midpoint-radius balls) because
products and quotients of wide intervals are much tighter in inf-sup form, and the checker works on wide boxes.

Conventions shared with the search (kisscrit/completeness/fast.py):
  * division by an interval containing 0 gives (-inf, inf), unless the numerator is exactly [0, 0];
  * a product with an exactly-zero factor [0, 0] is exactly [0, 0] (structural zeros of gradients);
  * sin, cos are evaluated at the endpoints, and the interval is widened to +-1 where an extremum
    2 k pi or (2k+1) pi may lie inside (decided with Arb's pi, conservatively).
"""
import math
from flint import arb, ctx

ctx.prec = 64
INF = math.inf
_PI = arb.pi()


def down(x):
    """Largest double <= every point of the ball x (x finite)."""
    l = x.lower()
    f = float(l)
    while arb(f) > l: f = math.nextafter(f, -INF)
    return f


def up(x):
    u = x.upper()
    f = float(u)
    while arb(f) < u: f = math.nextafter(f, INF)
    return f


class Iv:
    __slots__ = ('lo', 'hi')

    def __init__(self, lo, hi=None):
        self.lo = float(lo); self.hi = float(lo if hi is None else hi)

    def __repr__(self): return '[%r, %r]' % (self.lo, self.hi)
    def is_zero(self): return self.lo == 0.0 and self.hi == 0.0
    def finite(self): return math.isfinite(self.lo) and math.isfinite(self.hi)
    def contains0(self): return self.lo <= 0.0 <= self.hi
    def width(self): return self.hi - self.lo

    # arithmetic --------------------------------------------------------
    def __add__(self, o):
        o = _iv(o)
        return Iv(_op_down(lambda a, b: a + b, self.lo, o.lo), _op_up(lambda a, b: a + b, self.hi, o.hi))
    __radd__ = __add__

    def __sub__(self, o):
        o = _iv(o)
        return Iv(_op_down(lambda a, b: a - b, self.lo, o.hi), _op_up(lambda a, b: a - b, self.hi, o.lo))

    def __rsub__(self, o): return _iv(o) - self
    def __neg__(self): return Iv(-self.hi, -self.lo)

    def __mul__(self, o):
        o = _iv(o)
        if self.is_zero() or o.is_zero(): return Iv(0.0)
        los = []; his = []
        for a in (self.lo, self.hi):
            for b in (o.lo, o.hi):
                if (a == 0.0 and math.isinf(b)) or (b == 0.0 and math.isinf(a)):
                    return Iv(-INF, INF)
                los.append(_op_down(lambda x, y: x * y, a, b)); his.append(_op_up(lambda x, y: x * y, a, b))
        return Iv(min(los), max(his))
    __rmul__ = __mul__

    def __truediv__(self, o):
        o = _iv(o)
        if self.is_zero(): return Iv(0.0)
        if o.lo <= 0.0 <= o.hi: return Iv(-INF, INF)
        inv = Iv(_op_down(lambda a, b: a / b, 1.0, o.hi), _op_up(lambda a, b: a / b, 1.0, o.lo))
        return self * inv

    def __rtruediv__(self, o): return _iv(o) / self


def _iv(x): return x if isinstance(x, Iv) else Iv(x)


# The four basic operations on doubles are correctly rounded (IEEE 754), so the exact result lies within one
# ulp of the computed one; stepping one ulp outward gives a rigorous bound.  (Arb is used for the elementary
# functions, where no such guarantee exists.)
_nx = math.nextafter


_MAX = 1.7976931348623157e308


def _op_down(f, a, b):
    r = f(a, b)
    if math.isnan(r): return -INF
    if r == INF: return INF if (math.isinf(a) or math.isinf(b)) else _MAX     # finite overflow
    if r == -INF: return -INF
    return _nx(r, -INF)


def _op_up(f, a, b):
    r = f(a, b)
    if math.isnan(r): return INF
    if r == -INF: return -INF if (math.isinf(a) or math.isinf(b)) else -_MAX
    if r == INF: return INF
    return _nx(r, INF)


def hull(*xs):
    return Iv(min(x.lo for x in xs), max(x.hi for x in xs))


def clamp(x, lo, hi):
    """Image of [x.lo, x.hi] under t -> min(max(t, lo), hi)."""
    return Iv(min(max(x.lo, lo), hi), min(max(x.hi, lo), hi))


def _mono(fun, a, increasing=True, lo_lim=-INF, hi_lim=INF):
    """Rigorous range of a monotone function on [a.lo, a.hi] from its values at the endpoints."""
    def ev(t, d):
        if math.isinf(t):
            return hi_lim if (t > 0) == increasing else lo_lim
        v = fun(arb(t))
        return down(v) if d else up(v)
    if increasing:
        return Iv(max(ev(a.lo, True), lo_lim), min(ev(a.hi, False), hi_lim))
    return Iv(max(ev(a.hi, True), lo_lim), min(ev(a.lo, False), hi_lim))


def atan(a): return _mono(lambda t: t.atan(), a, True, down(-_PI / 2), up(_PI / 2))


def asin(a):
    """asin on [a.lo, a.hi] clamped to [-1, 1] (callers check validity separately)."""
    return _mono(lambda t: t.asin(), clamp(a, -1.0, 1.0), True)


def acos(a):
    return _mono(lambda t: t.acos(), clamp(a, -1.0, 1.0), False)


def cos(a):
    if not a.finite() or not (a.hi - a.lo < 6.0):
        return Iv(-1.0, 1.0)
    ca, cb = arb(a.lo).cos(), arb(a.hi).cos()
    lo = min(down(ca), down(cb)); hi = max(up(ca), up(cb))
    # an even multiple of pi inside -> max 1; an odd multiple -> min -1 (conservative where undecidable)
    k = math.floor(a.lo / math.pi) - 1
    while k * math.pi <= a.hi + 1.0:
        p = k * _PI
        possibly_inside = not (p < arb(a.lo)) and not (p > arb(a.hi))
        if possibly_inside:
            if k % 2 == 0: hi = 1.0
            else: lo = -1.0
        k += 1
    return Iv(max(lo, -1.0), min(hi, 1.0))


def sin(a):
    if not a.finite(): return Iv(-1.0, 1.0)
    # sin t = cos(t - pi/2), computed on the shifted interval with outward rounding
    s = Iv(down(arb(a.lo) - _PI / 2), up(arb(a.hi) - _PI / 2))
    return cos(s)


def square(a):
    """Range of t^2."""
    if a.lo >= 0: return a * a
    if a.hi <= 0: return (-a) * (-a)
    m = max(-a.lo, a.hi)
    return Iv(0.0, (Iv(m) * Iv(m)).hi)


PI = Iv(down(_PI), up(_PI))
HALF_PI = Iv(down(_PI / 2), up(_PI / 2))
TWO_PI = Iv(down(2 * _PI), up(2 * _PI))
