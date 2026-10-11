"""Linear-relaxation contractor for the vertex angle-sum equations S_v(p) = 2 pi.

For a box B with midpoint m and an interval enclosure J(B) of the Jacobian (automatic differentiation over
intervals), every p in B satisfies S(p) in S(m) + J(B)(p - m)  (mean value theorem).  Writing J(B) = Jc + [-Jr, Jr],
any zero p in B satisfies the linear inequalities
    | S_v(m) - 2 pi + sum_j Jc_vj (p_j - m_j) | <= sum_j Jr_vj w_j / 2 + rad(S_v(m)),
and we tighten each coordinate of B by linear programming over these inequalities.  This removes no zero of S.
"""
import numpy as np
from scipy.optimize import linprog
from . import ad as AD_
from .ad import AD
from . import interval as iv
from .interval import I


def _corners_ad(S, P):
    """Corner AD values per face (P = list of AD parameters)."""
    th = P[0]; c = AD_.cos(th)
    al = AD_.arccos(c / (1 + c))
    def beta(t, A):
        return AD_.arctan(AD_.cos(A * 0.5) / (AD_.cos(t) * AD_.sin(A * 0.5)))
    out = []
    for fi, k in enumerate(S.sizes):
        if k == 3:
            out.append((al, al, al))
        elif k == 4:
            A = P[S.pidx[fi][0]]; O = beta(th, A) * 2
            out.append((A, O, A, O))
        else:
            i1, i2 = S.pidx[fi]; f1, f2 = P[i1], P[i2]
            ct, st = AD_.cos(th), AD_.sin(th)
            d2 = AD_.arcsin(st * AD_.sin(f1 * 0.5)) * 2
            b1 = beta(th, f1)
            g = f2 - b1
            cg, sg = AD_.cos(g), AD_.sin(g)
            cd2, sd2 = AD_.cos(d2), AD_.sin(d2)
            cd3 = cd2 * ct + sd2 * st * cg
            d3 = AD_.arccos(cd3)
            a0 = np.pi / 2 - AD_.arctan((sd2 * ct - cd2 * st * cg) / (sg * st))
            a3 = np.pi / 2 - AD_.arctan((st * cd2 - ct * sd2 * cg) / (sg * sd2))
            h = AD_.arcsin(AD_.sin(d3 * 0.5) / st)
            f4 = h * 2
            b3 = beta(th, f4)
            out.append((b1 + a0 + b3, f1, f2, a3 + b3, f4))
    return out


def sums_and_jacobian(S, lo, hi):
    n = S.npar
    P = [AD.var(I(lo[k], hi[k]), k, n) for k in range(n)]
    C = _corners_ad(S, P)
    rows_lo = []; rows_hi = []; vals = []
    for v in range(S.n):
        s = None
        for (_, fi, pos) in S.rot[v]:
            s = C[fi][pos] if s is None else s + C[fi][pos]
        vals.append(s.v); rows_lo.append(s.g.lo); rows_hi.append(s.g.hi)
    return vals, np.array(rows_lo), np.array(rows_hi)


def contract(S, lo, hi, rounds=2):
    """Return contracted (lo, hi), or None if the box contains no zero of the vertex equations."""
    for _ in range(rounds):
        m = 0.5 * (lo + hi); w = hi - lo
        try:
            Fm, _, _ = sums_and_jacobian(S, m, m)
            _, Jlo, Jhi = sums_and_jacobian(S, lo, hi)
        except (ZeroDivisionError, FloatingPointError):
            return lo, hi
        if not (np.all(np.isfinite(Jlo)) and np.all(np.isfinite(Jhi))):
            return lo, hi
        Fc = np.array([0.5 * (f.lo + f.hi) for f in Fm]) - 2 * np.pi
        Fr = np.array([0.5 * (f.hi - f.lo) for f in Fm]) + 1e-13
        Jc = 0.5 * (Jlo + Jhi); Jr = 0.5 * (Jhi - Jlo)
        R = Jr @ (w / 2) + Fr
        # Jc (p - m) <= R - Fc  and  -Jc (p - m) <= R + Fc
        A = np.vstack([Jc, -Jc]); b = np.concatenate([R - Fc + Jc @ m, R + Fc - Jc @ m])
        bounds = list(zip(lo, hi))
        new_lo = lo.copy(); new_hi = hi.copy()
        for k in range(S.npar):
            c = np.zeros(S.npar); c[k] = 1
            r1 = linprog(c, A_ub=A, b_ub=b, bounds=bounds, method='highs')
            if r1.status == 2: return None
            if r1.status != 0: return lo, hi
            r2 = linprog(-c, A_ub=A, b_ub=b, bounds=bounds, method='highs')
            if r2.status != 0: return lo, hi
            new_lo[k] = max(lo[k], r1.x[k] - 1e-12 - 1e-9 * w[k])
            new_hi[k] = min(hi[k], r2.x[k] + 1e-12 + 1e-9 * w[k])
        if np.any(new_lo > new_hi): return None
        shrink = np.max((w - (new_hi - new_lo)) / np.maximum(w, 1e-300))
        lo, hi = new_lo, new_hi
        if shrink < 0.05: break
    return lo, hi
