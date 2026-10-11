"""Branch and bound over the face-parameter box of a candidate stressed graph.

A box is excluded if
  (V) some vertex angle sum cannot equal 2 pi, or some face parameters cannot form a convex polygon; or
  (G) a Gordan certificate exists: vectors y_v in the tangent plane at each vertex (in the local frame given
      by the first edge of the rotation) with  y_u . t_u(e) + y_w . t_w(e) > 0  for every edge e = uw and every
      point of the box.  Then no nonnegative nonzero balanced stress exists on the graph anywhere in the box
      (Gordan's alternative), so the box contains no critical configuration with this stressed graph.
Boxes are split along their widest coordinate.  Boxes narrower than `tol` that are not excluded are returned.
"""
import numpy as np
from scipy.optimize import linprog
from . import interval as iv
from .interval import I


def eval_box(S, lo, hi):
    p = [I(lo[k], hi[k]) for k in range(S.npar)]
    try:
        C, ok = S.corners(p)
    except ZeroDivisionError:
        return None
    if not np.all(ok): return 'invalid'
    for f in C:
        for c in f:
            if c.lo > np.pi + 1e-12: return 'invalid'
    sums = S.vertex_sums(C)
    for s in sums:
        if not (s.lo <= 2 * np.pi <= s.hi): return 'vsum'
    return C


def directions_iv(S, C):
    out = []
    for v in range(S.n):
        d = {}; acc = None
        for (w, fi, pos) in S.rot[v]:
            if acc is None:
                d[w] = (I(1.0, 1.0), I(0.0, 0.0)); acc = C[fi][pos]
            else:
                d[w] = (iv.cos(acc), iv.sin(acc)); acc = acc + C[fi][pos]
        out.append(d)
    return out


def gordan_ok(S, D, y):
    for (u, w) in S.edges:
        cu, su = D[u][w]; cw, sw = D[w][u]
        val = cu * y[2 * u] + su * y[2 * u + 1] + cw * y[2 * w] + sw * y[2 * w + 1]
        if val.lo <= 0: return False
    return True


def gordan_lp(S, p):
    C, ok = S.corners(list(p))
    if not np.all(ok): return None
    phis = S.directions(C)
    n = S.n; E = len(S.edges)
    A = np.zeros((E, 2 * n + 1))
    for k, (u, w) in enumerate(S.edges):
        A[k, 2 * u] = -np.cos(phis[u][w]); A[k, 2 * u + 1] = -np.sin(phis[u][w])
        A[k, 2 * w] = -np.cos(phis[w][u]); A[k, 2 * w + 1] = -np.sin(phis[w][u])
        A[k, -1] = 1
    cost = np.zeros(2 * n + 1); cost[-1] = -1
    res = linprog(cost, A_ub=A, b_ub=np.zeros(E), bounds=[(-1, 1)] * (2 * n) + [(None, 1)], method='highs')
    if res.status != 0 or res.x[-1] < 1e-10: return None
    return res.x[:-1]


def run(S, th_lo, th_hi, tol=1e-6, max_boxes=2_000_000, verbose=False):
    lo0, hi0 = S.initial_box(th_lo, th_hi)
    stack = [(lo0, hi0, None)]
    stats = dict(boxes=0, vsum=0, invalid=0, gordan=0, lp=0)
    terminal = []
    while stack:
        lo, hi, y = stack.pop()
        stats['boxes'] += 1
        if stats['boxes'] > max_boxes:
            stats['aborted'] = True; terminal += [(l, h) for l, h, _ in stack] + [(lo, hi)]; break
        r = eval_box(S, lo, hi)
        if isinstance(r, str):
            stats[r] += 1; continue
        y2 = None
        if r is not None:                 # r is None: enclosure failed (division by 0); just split
            C = r; D = directions_iv(S, C)
            if y is not None and gordan_ok(S, D, y):
                stats['gordan'] += 1; continue
            mid = 0.5 * (lo + hi)
            y2 = gordan_lp(S, mid); stats['lp'] += 1
            if y2 is not None and gordan_ok(S, D, y2):
                stats['gordan'] += 1; continue
        else:
            stats['noeval'] = stats.get('noeval', 0) + 1
        w = hi - lo; k = int(np.argmax(w))
        if w[k] < tol:
            terminal.append((lo, hi)); continue
        m = 0.5 * (lo[k] + hi[k])
        hi1 = hi.copy(); hi1[k] = m; lo2 = lo.copy(); lo2[k] = m
        yy = y2 if y2 is not None else y
        stack.append((lo, hi1, yy)); stack.append((lo2, hi, yy))
        if verbose and stats['boxes'] % 20000 == 0:
            print(stats, len(stack), len(terminal), flush=True)
    return terminal, stats
