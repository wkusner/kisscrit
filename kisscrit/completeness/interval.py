"""Minimal interval arithmetic on numpy arrays (lo, hi), with outward rounding.

Status: untrusted (the search's interval arithmetic with libm padding; no verdict depends on it).

Every arithmetic result is widened by one ulp in each direction (np.nextafter); every libm function
(cos, sin, arccos, sqrt) is widened by 4 ulps, which covers the documented accuracy of glibc's
double-precision implementations (< 1 ulp for these functions).  This is the usual standard for
interval libraries built on a non-correctly-rounded libm; the final certification steps use mpmath.
"""
import numpy as np

INF = np.inf
STRICT_DIV = True     # False: division by an interval containing 0 gives (-inf, inf) (batch mode)


def _dn(x, k=1):
    for _ in range(k): x = np.nextafter(x, -INF)
    return x


def _up(x, k=1):
    for _ in range(k): x = np.nextafter(x, INF)
    return x


class I:
    __slots__ = ('lo', 'hi')
    __array_priority__ = 1000

    def __init__(self, lo, hi=None):
        lo = np.asarray(lo, float)
        self.lo = lo
        self.hi = lo.copy() if hi is None else np.asarray(hi, float)

    @staticmethod
    def pt(x):
        return I(x, x)

    def mid(self): return 0.5 * (self.lo + self.hi)
    def width(self): return self.hi - self.lo
    def __repr__(self): return 'I(%s, %s)' % (self.lo, self.hi)

    def __add__(self, o):
        o = _as(o); return I(_dn(self.lo + o.lo), _up(self.hi + o.hi))
    __radd__ = __add__

    def __sub__(self, o):
        o = _as(o); return I(_dn(self.lo - o.hi), _up(self.hi - o.lo))

    def __rsub__(self, o):
        return _as(o) - self

    def __neg__(self): return I(-self.hi, -self.lo)

    def __mul__(self, o):
        o = _as(o)
        with np.errstate(invalid='ignore', over='ignore'):
            c = np.stack(np.broadcast_arrays(self.lo * o.lo, self.lo * o.hi, self.hi * o.lo, self.hi * o.hi))
        bad = np.isnan(c).any(0)
        lo = np.where(bad, -INF, c.min(0)); hi = np.where(bad, INF, c.max(0))
        return I(_dn(lo), _up(hi))
    __rmul__ = __mul__

    def __truediv__(self, o):
        o = _as(o)
        z = (o.lo <= 0) & (o.hi >= 0)
        if not STRICT_DIV or not np.any(z):
            with np.errstate(divide='ignore', invalid='ignore'):
                r = self * I(_dn(1 / o.hi), _up(1 / o.lo))
            if np.any(z):
                if STRICT_DIV: raise ZeroDivisionError('interval division by an interval containing 0')
                r = I(np.where(z, -INF, r.lo), np.where(z, INF, r.hi))   # whole line where divisor meets 0
            return r
        raise ZeroDivisionError('interval division by an interval containing 0')

    def __rtruediv__(self, o):
        return _as(o) / self

    def __getitem__(self, k): return I(self.lo[k], self.hi[k])


def _as(x):
    return x if isinstance(x, I) else I(x, x)


def sqrt(x):
    lo = np.maximum(x.lo, 0.0)
    return I(_dn(np.sqrt(lo), 4), _up(np.sqrt(np.maximum(x.hi, 0.0)), 4))


def arccos(x):
    """arccos of x clipped to [-1, 1] (caller must check that x meets [-1, 1])."""
    lo = np.clip(x.lo, -1, 1); hi = np.clip(x.hi, -1, 1)
    return I(_dn(np.arccos(hi), 4), _up(np.arccos(lo), 4))


def cos(x):
    a, b = np.nan_to_num(x.lo, nan=-INF), np.nan_to_num(x.hi, nan=INF)
    wide = ~((b - a) < 2 * np.pi)                      # includes infinite endpoints
    a = np.where(wide, 0.0, a); b = np.where(wide, 0.0, b)
    ca, cb = np.cos(a), np.cos(b)
    lo = np.minimum(ca, cb); hi = np.maximum(ca, cb)
    # interior extrema: cos = 1 at 2k pi, -1 at (2k+1) pi
    k1 = np.ceil(a / (2 * np.pi)); has_max = k1 * 2 * np.pi <= b
    k2 = np.ceil((a - np.pi) / (2 * np.pi)); has_min = k2 * 2 * np.pi + np.pi <= b
    hi = np.where(has_max | wide, 1.0, hi); lo = np.where(has_min | wide, -1.0, lo)
    return I(np.maximum(_dn(lo, 4), -1.0), np.minimum(_up(hi, 4), 1.0))


def sin(x):
    return cos(x - np.pi / 2)


def contains(x, v):
    return (x.lo <= v) & (v <= x.hi)


def hull(xs):
    return I(np.min([x.lo for x in xs], 0), np.max([x.hi for x in xs], 0))


def arctan(x):
    return I(_dn(np.arctan(x.lo), 4), _up(np.arctan(x.hi), 4))


def arcsin(x):
    lo = np.clip(x.lo, -1, 1); hi = np.clip(x.hi, -1, 1)
    return I(_dn(np.arcsin(lo), 4), _up(np.arcsin(hi), 4))
