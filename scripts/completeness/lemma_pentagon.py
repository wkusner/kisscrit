"""Lemma P.  Let P be a convex spherical pentagon with all sides theta in [60 deg, 63.4350 deg], all corners
>= alpha(theta) and all diagonals >= theta.  Then no two adjacent corners of P are >= 180 - 25 deg.

Proof by interval branch and bound over the chart (theta, f1, f2) (corners at P1, P2) with the apex corner at
P4 unclamped (every such pentagon is represented; only values are used, no derivatives).  A box is discarded
when the parameters are invalid, when some corner is < alpha on the whole box, or when every adjacent pair of
corners has a member < 155 deg on the whole box.
"""
import sys, time
import numpy as np
sys.path.insert(0, '.')
from kisscrit.completeness.model import Model, D, I, cos, sin, alpha_value, PI, DELTA


def corners(lo, hi):
    B = lo.shape[1]; z = np.zeros((3, B))
    var = lambda k: D(I(lo[k], hi[k]), I(z, z.copy()))
    th, f1, f2 = var(0), var(1), var(2)
    ct, st = cos(th), sin(th)
    hmin = 0.5 * alpha_value(lo[0], lo[0]).lo
    cs, ok = Model.pentagon(th, ct, st, hmin, 0.5 * PI, f1, f2)
    return [c.v for c in cs], ok


def main(tol=1e-7, maxboxes=10 ** 8):
    L = PI - DELTA
    th0, th1 = np.radians(60.0), np.radians(63.435)
    a0 = float(alpha_value(np.array(th0), np.array(th0)).lo)
    Q = [np.array([th0, a0, a0]), np.array([th1, PI, PI])]
    Qlo = [Q[0]]; Qhi = [Q[1]]; n = 0; t0 = time.time()
    while Qlo:
        lo = np.stack(Qlo[-4096:], 1); hi = np.stack(Qhi[-4096:], 1); del Qlo[-4096:]; del Qhi[-4096:]
        n += lo.shape[1]
        cs, ok = corners(lo, hi)
        al = alpha_value(lo[0], hi[0])
        for c in cs: ok &= c.hi >= al.lo
        bad = np.zeros(lo.shape[1], bool)
        for i in range(5):
            bad |= (cs[i].hi >= L) & (cs[(i + 1) % 5].hi >= L)
        keep = ok & bad
        for b in np.nonzero(keep)[0]:
            w = (hi[:, b] - lo[:, b]); k = int(np.argmax(w))
            if w[k] < tol: raise SystemExit('Lemma P not proved: box %s %s' % (lo[:, b], hi[:, b]))
            m = 0.5 * (lo[k, b] + hi[k, b])
            h1 = hi[:, b].copy(); h1[k] = m; l2 = lo[:, b].copy(); l2[k] = m
            Qlo += [lo[:, b].copy(), l2]; Qhi += [h1, hi[:, b].copy()]
        if n > maxboxes: raise SystemExit('too many boxes')
    print('Lemma P proved: %d boxes, %.1fs' % (n, time.time() - t0))


if __name__ == '__main__':
    main()
