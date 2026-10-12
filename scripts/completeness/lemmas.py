#!/usr/bin/env python3
"""Computer-assisted lemmas used by the completeness computation (statements in docs/lemmas.md).

All computations use the checker's interval arithmetic (kisscrit/core/interval.py: IEEE-rounded basic
operations stepped outward, Arb for elementary functions) and the checker's face model (kisscrit/core/model.py).

  python scripts/completeness/lemmas.py P      Lemma P: no two adjacent corners of a face pentagon are >= 155 deg
  python scripts/completeness/lemmas.py R [i K] Lemma R: a face polygon contains no point at distance >= theta
                                               from all its vertices (so there are no rattlers)
Each prints the number of boxes and ends with "proved" or raises SystemExit with the box it could not exclude.
Theta ranges over [60 deg, 63.4350 deg] (the closed range, a little more than needed).
"""
import math, sys, os, time
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
from kisscrit.core import interval as I
from kisscrit.core.interval import Iv, PI, HALF_PI
from kisscrit.core import model as M

TH0, TH1 = math.radians(60.0), math.radians(63.435)
LSTRAIGHT = M.LSTRAIGHT                       # pi - 25 deg, the region threshold of the search


def alpha_lo(t): return M.alpha_range(t, t).lo


def pentagon_values(box, unclamped=True):
    """Corner value intervals (formula order) of the pentagon chart (theta, f1, f2) on a box, apex unclamped."""
    (t0, t1), (a0, a1), (b0, b1) = box
    th = M.var(t0, t1, 0); ct, st = M.cos_(th), M.sin_(th)
    hmin = 0.5 * alpha_lo(t0)
    hmax = HALF_PI.lo if unclamped else 0.5 * LSTRAIGHT
    cs, ok = M.pentagon(ct, st, hmin, hmax, M.var(a0, a1, 1), M.var(b0, b1, 2))
    return [c[0] for c in cs], ok


def bisect(box):
    k = max(range(len(box)), key=lambda j: box[j][1] - box[j][0])
    lo, hi = box[k]; m = 0.5 * (lo + hi)
    b1 = list(box); b1[k] = (lo, m); b2 = list(box); b2[k] = (m, hi)
    return b1, b2


def run(start, excluded, name, min_width=1e-9, max_boxes=10 ** 8):
    stack = [start]; n = 0; t = time.time()
    vol = lambda b: math.prod(h - l for l, h in b)
    v0, vdone = vol(start), 0.0
    while stack:
        box = stack.pop(); n += 1
        if excluded(box):
            vdone += vol(box); continue
        if max(b - a for a, b in box) < min_width or n > max_boxes:
            raise SystemExit('%s NOT proved: cannot exclude box %s' % (name, box))
        stack.extend(bisect(box))
        if n % 20000 == 0: print('  %s: %d boxes, %d open, %.4f%% of the volume done (%.0fs)' % (name, n, len(stack), 100 * vdone / v0, time.time() - t), flush=True)
    print('%s proved: %d boxes, %.0fs' % (name, n, time.time() - t))


# ---------------------------------------------------------------- Lemma P
def lemma_P():
    """Every convex equilateral pentagon with side theta in [60, 63.435] deg, all corners >= alpha(theta) and
    the diagonal P0 P3 >= theta (any pentagon, in the unclamped chart (theta, f1, f2)) has no two adjacent
    corners >= LSTRAIGHT."""
    a = alpha_lo(TH0)
    def excluded(box):
        cs, ok = pentagon_values(box)
        if not ok: return True
        al = M.alpha_range(box[0][0], box[0][1]).lo
        if any(c.hi < al for c in cs): return True          # a corner below alpha(theta)
        return all(cs[i].hi < LSTRAIGHT or cs[(i + 1) % 5].hi < LSTRAIGHT for i in range(5))
    run([(TH0, TH1), (a, math.pi), (a, math.pi)], excluded, 'Lemma P')


# ---------------------------------------------------------------- Lemma R
def rot(v, axis, c, s):
    """Rotate the interval vector v about the unit interval vector axis by the angle with cos c, sin s."""
    ax, ay, az = axis; vx, vy, vz = v
    cx = (ay * vz - az * vy, az * vx - ax * vz, ax * vy - ay * vx)
    d = ax * vx + ay * vy + az * vz
    one_c = Iv(1.0) - c
    return tuple(v[i] * c + cx[i] * s + axis[i] * d * one_c for i in range(3))


def dot(u, v): return u[0] * v[0] + u[1] * v[1] + u[2] * v[2]


def det(u, v, w):
    return (u[0] * (v[1] * w[2] - v[2] * w[1]) - u[1] * (v[0] * w[2] - v[2] * w[0])
            + u[2] * (v[0] * w[1] - v[1] * w[0]))


def polygon_vertices(th, corners):
    """Vertices P0..P_{k-1} of the counterclockwise equilateral spherical polygon with side th and the given
    corners (intervals): P0 = north pole, P1 at bearing 0, each next vertex from the corner at the previous
    one.  Returns the list of interval vectors."""
    ct, st = I.cos(th), I.sin(th)
    P = [(Iv(0.0), Iv(0.0), Iv(1.0)), (st, Iv(0.0), ct)]
    k = len(corners)
    for i in range(1, k - 1):
        # the corner at P_i spans counterclockwise from P_i -> P_{i+1} to P_i -> P_{i-1}: rotate P_{i-1}
        # about P_i by -corner
        c = corners[i]
        P.append(rot(P[i - 1], P[i], I.cos(c), -I.sin(c)))
    return P


def cross(u, v): return (u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0])


def sub(u, v): return (u[0] - v[0], u[1] - v[1], u[2] - v[2])


def triple_excluded(d, close, dets, gaps):
    """The decision of no_deep_hole for one triple, from enclosures of d = det(Pa, Pb, Pc) (= n.Pa),
    close = (d^2 - cos^2(th) |n|^2 > 0 proved), dets[m] = det(P_m, P_m+1, n), gaps = n.P_m - d (other m)."""
    for s in (1, -1):
        sd = d if s > 0 else -d
        if sd.hi <= 0: continue                                          # A
        if sd.lo > 0 and close: continue                                 # B
        if any((e.hi if s > 0 else -e.lo) < 0 for e in dets): continue   # C
        if any((g.lo if s > 0 else -g.hi) > 0 for g in gaps): continue   # D
        return False
    return True


TRIPLES = {k: [(a, b, c) for a in range(k) for b in range(a + 1, k) for c in range(b + 1, k)] for k in (3, 4, 5)}


def no_deep_hole(th, P):
    """True if it is proved that no point of the convex polygon P (interval vertices) is at distance >= th from
    every vertex.  By the argument in docs/lemmas.md (L9), such a point would be the circumcenter c of three
    vertices, with circumradius R in [th, 90 deg), inside the polygon, and no vertex closer than R.  So for every
    triple and both orientations s of its normal n = (Pj - Pi) x (Pk - Pi) (c = s n / |n|) we show one of:
      A  s n.Pi <= 0                 (R >= 90 deg)
      B  s n.Pi > 0 and (n.Pi)^2 > cos^2(th) |n|^2     (R < th)
      C  s det(P_m, P_m+1, n) < 0 for some edge         (c outside the polygon)
      D  s (n.P_m - n.Pi) > 0 for some vertex m         (a vertex closer than R)."""
    k = len(P); ct2 = I.square(I.cos(th))
    for a, b, c in TRIPLES[k]:
        n = cross(sub(P[b], P[a]), sub(P[c], P[a]))
        d = dot(n, P[a])
        close = (I.square(d) - ct2 * dot(n, n)).lo > 0
        dets = [det(P[m], P[(m + 1) % k], n) for m in range(k)]
        gaps = [dot(n, P[m]) - d for m in range(k) if m not in (a, b, c)]
        if not triple_excluded(d, close, dets, gaps): return False
    return True


def pentagon_by_vertices(th, c0, f1, cosL):
    """Vertices of the counterclockwise equilateral pentagon with side th, corner c0 at P0 and f1 at P1, with
    the corner at P3 at most LSTRAIGHT.  P0 = north pole, P1 at bearing 0, P4 at bearing c0, P2 from f1 (one
    rotation each, which keeps the interval enclosures tight), and P3 the apex over the diagonal P2 P4:
    P3 = l (P2 + P4) + m (P4 x P2) with l = cos th / (1 + c), m^2 = (1 - 2 cos^2 th / (1 + c)) / (1 - c^2),
    c = P2.P4; the sign of m puts P3 on the far side of P2 P4 from P0 (the convex side).
    Every convex equilateral pentagon has such a labelling: by Lemma P some corner is < LSTRAIGHT; call it P3.
    Returns (P, possibly_valid); possibly_valid is False when P3 certainly does not exist or its corner is
    certainly > LSTRAIGHT (by the law of cosines, corner(P3) <= L iff c >= cos^2 th + sin^2 th cos L)."""
    ct, st = I.cos(th), I.sin(th)
    one = Iv(1.0)
    P0 = (Iv(0.0), Iv(0.0), one); P1 = (st, Iv(0.0), ct)
    P4 = (st * I.cos(c0), st * I.sin(c0), ct)
    P2 = rot(P0, P1, I.cos(f1), -I.sin(f1))
    c = dot(P2, P4)
    if c.hi < (I.square(ct) + I.square(st) * cosL).lo: return None, False
    num = one - Iv(2.0) * I.square(ct) / (one + c)
    if num.hi < 0: return None, False
    l = ct / (one + c); m = I.sqrt(num / (one - I.square(c)))
    X = cross(P4, P2)
    P3 = tuple(l * (P2[i] + P4[i]) + m * X[i] for i in range(3))
    return [P0, P1, P2, P3, P4], True


# ---- pentagons: the same tests in mean-value form.  Every quantity q is a function of the box parameters
# (theta, c0, f1); it is evaluated with forward differentiation on the box (value q(B), gradient g(B)) and at the
# midpoint m, and enclosed by q(B) intersected with q(m) + g(B) . (B - m) (mean value theorem).  This removes
# most of the dependency blow-up of plain interval evaluation, which otherwise needs boxes ~100x smaller.
def _adot(u, v): return M.add(M.add(M.mul(u[0], v[0]), M.mul(u[1], v[1])), M.mul(u[2], v[2]))
def _asub(u, v): return tuple(M.sub(x, y) for x, y in zip(u, v))
def _across(u, v): return (M.sub(M.mul(u[1], v[2]), M.mul(u[2], v[1])), M.sub(M.mul(u[2], v[0]), M.mul(u[0], v[2])),
                           M.sub(M.mul(u[0], v[1]), M.mul(u[1], v[0])))
def _adet(u, v, w): return _adot(u, _across(v, w))
def _asqrt(a):
    v = I.sqrt(a[0])
    return M.chain(a, v, Iv(0.5) / v)


def _arot(v, axis, c, s):
    cx = _across(axis, v); d = _adot(axis, v); one_c = M.sub(M.const(1.0), c)
    return tuple(M.add(M.add(M.mul(v[i], c), M.mul(cx[i], s)), M.mul(M.mul(axis[i], d), one_c)) for i in range(3))


def pentagon_quantities(t, c0, f1, cosL):
    """All quantities the pentagon test needs, as AD values (see pentagon_by_vertices for the construction).
    Returns None when the pentagon certainly does not exist in this chart."""
    ct, st = M.cos_(t), M.sin_(t)
    one = M.const(1.0); zero = M.const(0.0)
    P0 = (zero, zero, one); P1 = (st, zero, ct)
    P4 = (M.mul(st, M.cos_(c0)), M.mul(st, M.sin_(c0)), ct)
    P2 = _arot(P0, P1, M.cos_(f1), M.neg(M.sin_(f1)))
    c = _adot(P2, P4); ct2 = M.mul(ct, ct)
    if c[0].hi < (ct2[0] + M.mul(st, st)[0] * cosL).lo: return None          # corner at P3 > LSTRAIGHT
    num = M.sub(one, M.div(M.scale(ct2, 2.0), M.add(one, c)))
    if num[0].hi < 0: return None                                           # P3 does not exist
    l = M.div(ct, M.add(one, c)); m = _asqrt(M.div(num, M.sub(one, M.mul(c, c))))
    X = _across(P4, P2)
    P3 = tuple(M.add(M.mul(l, M.add(P2[i], P4[i])), M.mul(m, X[i])) for i in range(3))
    P = [P0, P1, P2, P3, P4]
    q = {'turn': [_adet(P[i - 1], P[i], P[(i + 1) % 5]) for i in range(5)],
         'diag': [M.sub(_adot(P[i - 1], P[(i + 1) % 5]), ct) for i in range(5)]}
    for a, b, cc in TRIPLES[5]:
        n = _across(_asub(P[b], P[a]), _asub(P[cc], P[a]))
        d = _adet(P[a], P[b], P[cc])
        q[(a, b, cc)] = (d, M.sub(M.mul(d, d), M.mul(ct2, _adot(n, n))),
                         [_adet(P[i], P[(i + 1) % 5], n) for i in range(5)],
                         [M.sub(_adot(n, P[i]), d) for i in range(5) if i not in (a, b, cc)])
    return q


def _mv(qb, qm, rad):
    """Mean-value enclosure of one quantity: qb on the box (with gradient), qm at the midpoint."""
    v = qm[0] + qb[1] * rad[0] + qb[2] * rad[1] + qb[3] * rad[2]
    lo, hi = max(v.lo, qb[0].lo), min(v.hi, qb[0].hi)
    return Iv(lo, hi) if lo <= hi else qb[0]


def pentagon_excluded(box, cosL):
    """Lemma R test for the pentagon chart (theta, c0, f1) on a box; True if no pentagon of the box is valid or
    none has a deep hole."""
    vb = [M.var(box[k][0], box[k][1], k) for k in range(3)]
    mid = [0.5 * (box[k][0] + box[k][1]) for k in range(3)]
    vm = [M.var(mid[k], mid[k], k) for k in range(3)]
    rad = [Iv(box[k][0]) - Iv(mid[k]) for k in range(3)]
    rad = [Iv(r.lo, (Iv(box[k][1]) - Iv(mid[k])).hi) for k, r in enumerate(rad)]
    qb = pentagon_quantities(vb[0], vb[1], vb[2], cosL)
    if qb is None: return True
    qm = pentagon_quantities(vm[0], vm[1], vm[2], cosL)
    if qm is None:                                         # the midpoint is invalid: fall back to plain values
        enc = lambda x: x[0]
    else:
        enc = None
    def E(key_b, key_m): return key_b[0] if enc else _mv(key_b, key_m, rad)
    for i in range(5):
        if E(qb['turn'][i], qm and qm['turn'][i]).hi < 0: return True           # a corner > pi (L3)
        if E(qb['diag'][i], qm and qm['diag'][i]).lo > 0: return True           # a diagonal < theta (L8)
    for t in TRIPLES[5]:
        db, Qb, Db, Gb = qb[t]
        dm, Qm, Dm, Gm = qm[t] if qm else (None, None, [None] * 5, [None] * 2)
        if not triple_excluded(E(db, dm), E(Qb, Qm).lo > 0, [E(x, y) for x, y in zip(Db, Dm)],
                               [E(x, y) for x, y in zip(Gb, Gm)]):
            return False
    return True


def lemma_R(part=0, parts=1):
    """Lemma R: a face polygon (triangle, rhombus or pentagon with side theta, corners in [alpha(theta), pi],
    diagonals >= theta) contains no point at distance >= theta from all its vertices.  Checked with no_deep_hole
    over the face parameters only; the reduction to circumcenters is proved in docs/lemmas.md (L9)."""
    a = alpha_lo(TH0)

    def excl_tri(box):
        th = Iv(*box[0]); al = M.alpha_range(*box[0])
        return no_deep_hole(th, polygon_vertices(th, [al, al, al]))
    run([(TH0, TH1)], excl_tri, 'Lemma R (triangles)')

    # rhombi: corners (A, O, A, O) with O = 2 beta(theta, A), A in [alpha, 2 alpha]
    def excl_quad(box):
        C, ok = M.Graph([[0, 1, 2, 3], [3, 2, 1, 0]], {}).corners([box[0][0], box[1][0], box[1][0]],
                                                                  [box[0][1], box[1][1], box[1][1]])
        if not ok: return True
        cs = [c[0] for c in C[0]]
        al = M.alpha_range(*box[0]).lo
        if any(c.hi < al or c.lo > math.pi for c in cs): return True    # corners in [alpha, pi] (L3)
        th = Iv(*box[0])
        return no_deep_hole(th, polygon_vertices(th, cs))
    run([(TH0, TH1), (a, 2 * M.alpha_range(TH1, TH1).hi)], excl_quad, 'Lemma R (rhombi)')

    # pentagons: vertex chart (theta, c0, f1), corners c0 at P0 and f1 at P1 (see pentagon_by_vertices)
    cosL = I.cos(Iv(LSTRAIGHT))

    def excl_pent(box): return pentagon_excluded(box, cosL)
    t0 = TH0 + (TH1 - TH0) * part / parts; t1 = TH1 if part == parts - 1 else TH0 + (TH1 - TH0) * (part + 1) / parts
    run([(t0, t1), (a, math.pi), (a, math.pi)], excl_pent,
        'Lemma R (pentagons%s)' % ('' if parts == 1 else ', theta part %d of %d' % (part + 1, parts)))


if __name__ == '__main__':
    which = sys.argv[1:] or ['P', 'R']
    if 'P' in which: lemma_P()
    if 'R' in which:
        # optional: R <part> <parts> splits the pentagon theta range, to run the parts in parallel
        rest = which[which.index('R') + 1:]
        if len(rest) >= 2 and rest[0].isdigit(): lemma_R(int(rest[0]), int(rest[1]))
        else: lemma_R()
