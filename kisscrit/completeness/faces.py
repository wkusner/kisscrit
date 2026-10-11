"""Corner angles of the faces of a candidate stressed graph, as functions of the face parameters.

All faces are convex spherical polygons with all sides equal to theta.  Formulas use half-angle (Napier)
and cotangent four-part forms, which are smooth on the whole parameter range (including straight corners):
  base angle of the isosceles triangle with legs t and apex A:   beta(t, A) = atan2(cos(A/2), cos t sin(A/2))
  base of that triangle:                                          d = 2 asin(sin t sin(A/2))
  triangle:      corners alpha(t) = arccos(cos t / (1 + cos t))
  quadrilateral: parameter A = corner at boundary positions 0, 2; corners at 1, 3 are 2 beta(t, A)
  pentagon:      parameters f1, f2 = corners at boundary positions 1, 2; diagonals from position 0.
Works on floats and on intervals (kisscrit.completeness.interval.I); returns (corners, ok) where ok is False
where the parameters cannot belong to a convex equilateral polygon.
"""
import numpy as np
from . import interval as iv


def _is_iv(x): return isinstance(x, iv.I)
def _cos(x): return iv.cos(x) if _is_iv(x) else np.cos(x)
def _sin(x): return iv.sin(x) if _is_iv(x) else np.sin(x)


def _atan_pos(Y, X):
    """atan2(Y, X) for Y >= 0, X > 0 (value in [0, pi/2])."""
    if _is_iv(X):
        return iv.arctan(Y / X)
    return np.arctan2(Y, X)


def _acot_posY(X, Y):
    """atan2(Y, X) for Y > 0 (value in (0, pi)) = pi/2 - atan(X / Y)."""
    if _is_iv(Y):
        return np.pi / 2 - iv.arctan(X / Y)
    return np.arctan2(Y, X)


def _asin(x):
    if _is_iv(x):
        return iv.arcsin(x), (x.lo <= 1) & (x.hi >= -1)
    return np.arcsin(np.clip(x, -1, 1)), (x <= 1 + 1e-12) & (x >= -1 - 1e-12)


def alpha(th):
    if _is_iv(th):
        lo = alpha(float(th.lo)); hi = alpha(float(th.hi))           # increasing in theta
        return iv.I(iv._dn(np.asarray(lo), 4), iv._up(np.asarray(hi), 4))
    c = np.cos(th); return np.arccos(c / (1 + c))


def beta(th, A):
    return _atan_pos(_cos(A * 0.5), _cos(th) * _sin(A * 0.5))


def quad_other(th, A):
    if _is_iv(th) or _is_iv(A):
        th = th if _is_iv(th) else iv.I(th, th); A = A if _is_iv(A) else iv.I(A, A)
        # 2 beta is increasing in theta and decreasing in A on (0, pi): exact range from two corners
        lo = 2 * beta(float(th.lo), float(A.hi)); hi = 2 * beta(float(th.hi), float(A.lo))
        return iv.I(iv._dn(np.asarray(lo), 4), iv._up(np.asarray(hi), 4)), True
    return 2 * beta(th, A), True


def pentagon(th, f1, f2):
    ct, st = _cos(th), _sin(th)
    d2, ok = _asin(st * _sin(f1 * 0.5)); d2 = 2 * d2
    b1 = beta(th, f1)
    g = f2 - b1
    if _is_iv(g):
        ok = ok & (g.hi > 0) & (g.lo < np.pi)
        if not ok: return None, False
        # the part of the box with gamma outside (0, pi) contains no convex pentagon: clip (sound)
        g = iv.I(np.maximum(g.lo, 1e-9), np.minimum(g.hi, np.pi - 1e-9))
    else:
        ok = ok & (g > 0) & (g < np.pi)
    cg, sg = _cos(g), _sin(g)
    cd2, sd2 = _cos(d2), _sin(d2)
    cd3 = cd2 * ct + sd2 * st * cg
    if _is_iv(cd3):
        d3 = iv.arccos(cd3); ok = ok & (cd3.hi >= -1) & (cd3.lo <= 1)
        ok = ok & (cd3.lo <= ct.hi)          # diagonal P0 P3 >= theta
    else:
        d3 = np.arccos(np.clip(cd3, -1, 1))
    a0 = _acot_posY(sd2 * ct - cd2 * st * cg, sg * st)       # angle at P0 in triangle P0 P2 P3
    a3 = _acot_posY(st * cd2 - ct * sd2 * cg, sg * sd2)      # angle at P3 in triangle P0 P2 P3
    h, o = _asin(_sin(d3 * 0.5) / st); ok = ok & o           # half of the corner at P4
    if _is_iv(h):
        # the corner at P4 is at least alpha(theta) > 0 in a convex pentagon with diagonals >= theta
        amin = 0.5 * float(alpha(float(th.lo)) if _is_iv(th) else alpha(th))
        if h.hi < amin: return None, False
        h = iv.I(np.maximum(h.lo, amin), np.minimum(h.hi, np.pi / 2))
    f4 = 2 * h
    b3 = beta(th, f4)
    return (b1 + a0 + b3, f1, f2, a3 + b3, f4), ok
