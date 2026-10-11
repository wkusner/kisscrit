"""Forward-mode automatic differentiation over intervals (or floats).

AD(v, g): value v (interval I or float) and gradient g (I with arrays, or numpy array) with respect to the
face parameters.  Supports +, -, *, /, cos, sin, arctan, arcsin, arccos and constants.
"""
import numpy as np
from . import interval as iv
from .interval import I


def _isI(x): return isinstance(x, I)


class AD:
    __slots__ = ('v', 'g')

    def __init__(self, v, g): self.v, self.g = v, g

    @staticmethod
    def var(v, k, n):
        g = np.zeros(n); g[k] = 1.0
        return AD(v, I(g, g) if _isI(v) else g)

    @staticmethod
    def const(c, n, like):
        return AD(I(c, c) if _isI(like) else c, I(np.zeros(n), np.zeros(n)) if _isI(like) else np.zeros(n))

    def _c(self, o):
        if isinstance(o, AD): return o
        n = len(self.g.lo) if _isI(self.g) else len(self.g)
        return AD.const(o, n, self.v)

    def __add__(self, o): o = self._c(o); return AD(self.v + o.v, self.g + o.g)
    __radd__ = __add__
    def __sub__(self, o): o = self._c(o); return AD(self.v - o.v, self.g - o.g)
    def __rsub__(self, o): return self._c(o) - self
    def __neg__(self): return AD(-self.v, -self.g)

    def __mul__(self, o):
        o = self._c(o); return AD(self.v * o.v, self.g * o.v + o.g * self.v)
    __rmul__ = __mul__

    def __truediv__(self, o):
        o = self._c(o); q = self.v / o.v
        return AD(q, (self.g - o.g * q) / o.v)

    def __rtruediv__(self, o): return self._c(o) / self


def _f(x, val, dval):
    return AD(val, x.g * dval)


def cos(x): return _f(x, iv.cos(x.v) if _isI(x.v) else np.cos(x.v), -(iv.sin(x.v) if _isI(x.v) else np.sin(x.v)))
def sin(x): return _f(x, iv.sin(x.v) if _isI(x.v) else np.sin(x.v), iv.cos(x.v) if _isI(x.v) else np.cos(x.v))


def arctan(x):
    if _isI(x.v): return _f(x, iv.arctan(x.v), 1.0 / (1.0 + x.v * x.v))
    return _f(x, np.arctan(x.v), 1.0 / (1.0 + x.v * x.v))


def arcsin(x):
    if _isI(x.v):
        return _f(x, iv.arcsin(x.v), 1.0 / iv.sqrt(1.0 - x.v * x.v))
    return _f(x, np.arcsin(np.clip(x.v, -1, 1)), 1.0 / np.sqrt(max(1 - x.v * x.v, 1e-300)))


def arccos(x):
    if _isI(x.v):
        return _f(x, iv.arccos(x.v), -1.0 / iv.sqrt(1.0 - x.v * x.v))
    return _f(x, np.arccos(np.clip(x.v, -1, 1)), -1.0 / np.sqrt(max(1 - x.v * x.v, 1e-300)))
