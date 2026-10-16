"""Batched branch and bound over face-parameter boxes of one candidate stressed graph.

Status: untrusted (superseded prototype (checklist B1)).

For one Structure S, every pentagon chart assignment (see model.py) is a separate starting box.  A box is
discarded when one of these sound tests succeeds on the whole box:
  invalid  the parameters cannot describe convex equilateral faces in this chart;
  vsum     some vertex angle sum cannot equal 2 pi;
  gordan   a Gordan certificate y exists (no nonzero nonnegative balanced stress on the graph);
and is otherwise contracted with a preconditioned interval Gauss-Seidel step on the vertex equations
(contraction keeps every zero of the equations in the box), then bisected.  Boxes whose width (relative to the
starting box) drops below `tol` are returned as terminal boxes, with their chart.

Soundness notes are in docs/completeness.md.  Every exclusion uses interval enclosures only; floating point
enters only in choosing preconditioners, Gordan vectors and split coordinates.
"""
import itertools, time
import numpy as np
from scipy.optimize import linprog
from . import interval as iv
from .interval import I
from .model import Model, PI


TWO_PI = I(2 * np.pi, 2 * np.pi)
TWO_PI = I(iv._dn(np.array(2 * np.pi)), iv._up(np.array(2 * np.pi)))


def _stack(vals):
    return I(np.stack([x.lo for x in vals]), np.stack([x.hi for x in vals]))


def evaluate(M, lo, hi, box=None):
    """Corners, sums - 2pi (I of shape (V, B)), Jacobian (I of shape (V, npar, B)), validity mask."""
    C, ok = M.corners(lo, hi, box)
    vs = M.vertex_sums(C)
    F = _stack([s.v for s in vs]) - TWO_PI
    J = I(np.stack([s.g.lo for s in vs]), np.stack([s.g.hi for s in vs]))
    for f in C:
        for c in f:
            ok &= c.v.lo <= PI + 1e-12
    vok = np.all((F.lo <= 0) & (F.hi >= 0), axis=0)
    return C, F, J, ok, vok


def gauss_seidel(M, lo, hi, J, nrank=9):
    """One preconditioned interval Gauss-Seidel sweep.  Returns new (lo, hi) and an emptiness mask."""
    n, B = lo.shape
    m = 0.5 * (lo + hi)
    _, Fm, Jm, okm, _ = evaluate(M, m, m.copy(), box=(lo, hi))
    Jc = 0.5 * (J.lo + J.hi)                                   # (V, n, B)
    fin = np.all(np.isfinite(J.lo) & np.isfinite(J.hi), axis=(0, 1)) & np.all(np.isfinite(Fm.lo), axis=0)
    newlo, newhi = lo.copy(), hi.copy(); empty = np.zeros(B, bool)
    for b in np.nonzero(fin)[0]:
        A = 0.5 * (Jm.lo[:, :, b] + Jm.hi[:, :, b])            # point Jacobian at the midpoint
        U, s, Vt = np.linalg.svd(A)
        r = min(nrank, int(np.sum(s > 1e-9 * s[0])), n)
        if r == 0: continue
        G = U[:, :r].T @ A                                        # (r, n)
        # pivot columns by QR with column pivoting (greedy)
        piv = []; R = G.copy()
        for i in range(r):
            norms = np.linalg.norm(R, axis=0); norms[piv] = -1
            k = int(np.argmax(norms)); piv.append(k)
            q = R[:, k] / np.linalg.norm(R[:, k]); R = R - np.outer(q, q @ R)
        Gk = G[:, piv]
        if np.linalg.cond(Gk) > 1e8: continue
        Pc = np.linalg.solve(Gk, U[:, :r].T)                     # (r, V)
        # interval M = Pc J(box), f = Pc F(m)
        Jb = I(J.lo[:, :, b], J.hi[:, :, b])
        Mlo = np.zeros((r, n)); Mhi = np.zeros((r, n))
        Pp, Pn = np.maximum(Pc, 0), np.minimum(Pc, 0)
        Mlo = iv._dn(Pp @ Jb.lo + Pn @ Jb.hi, 2 + Pc.shape[1])
        Mhi = iv._up(Pp @ Jb.hi + Pn @ Jb.lo, 2 + Pc.shape[1])
        flo = iv._dn(Pp @ Fm.lo[:, b] + Pn @ Fm.hi[:, b], 2 + Pc.shape[1])
        fhi = iv._up(Pp @ Fm.hi[:, b] + Pn @ Fm.lo[:, b], 2 + Pc.shape[1])
        l, h = newlo[:, b].copy(), newhi[:, b].copy()
        for i, k in enumerate(piv):
            mk = m[k, b]
            # sum_{j != k} M_ij (p_j - m_j)
            ylo = l - m[:, b]; yhi = h - m[:, b]
            Mi = I(Mlo[i], Mhi[i]); y = I(ylo, yhi)
            t = Mi * y
            tlo = t.lo.copy(); thi = t.hi.copy(); tlo[k] = 0; thi[k] = 0
            rhs = I(-fhi[i] - iv._up(np.sum(thi), n), -flo[i] - iv._dn(np.sum(tlo), n))
            rhs = I(iv._dn(rhs.lo, n), iv._up(rhs.hi, n))
            d = I(Mlo[i, k], Mhi[i, k])
            if d.lo <= 0 <= d.hi: continue
            q = rhs / d
            a = max(l[k], float(iv._dn(np.array(mk + q.lo), 2))); c = min(h[k], float(iv._up(np.array(mk + q.hi), 2)))
            if a > c: empty[b] = True; break
            l[k], h[k] = a, c
        newlo[:, b], newhi[:, b] = l, h
    return newlo, newhi, empty


# ---------------------------------------------------------------------- Gordan
def directions(M, C):
    """Per vertex dict w -> (cos, sin) intervals of the edge direction in the local frame."""
    S = M.S; out = []
    for v in range(S.n):
        d = {}; acc = None
        for (w, fi, pos) in S.rot[v]:
            if acc is None:
                d[w] = None; acc = C[fi][pos].v
            else:
                d[w] = (iv.cos(acc), iv.sin(acc)); acc = acc + C[fi][pos].v
        out.append(d)
    return out


def gordan_check(M, Dd, Y):
    """Y: (2V, B) vectors.  True where y_u.t_u(e) + y_w.t_w(e) > 0 for all edges, on the whole box."""
    S = M.S; B = Y.shape[1]; good = np.ones(B, bool)
    for (u, w) in S.edges:
        tot = None
        for (a, b) in ((u, w), (w, u)):
            if Dd[a][b] is None:
                term = I(Y[2 * a], Y[2 * a])
            else:
                c, s = Dd[a][b]; term = c * Y[2 * a] + s * Y[2 * a + 1]
            tot = term if tot is None else tot + term
        good &= tot.lo > 0
    return good


def gordan_lp(M, p, C=None, col=0):
    """Gordan vector at the point p (or from precomputed point corners C, batch column col)."""
    S = M.S
    if C is None:
        P = p[:, None]; C, ok = M.corners(P, P.copy())
        if not ok[0]: return None
    n = S.n; E = len(S.edges)
    phi = []
    for v in range(n):
        d = {}; acc = 0.0
        for (w, fi, pos) in S.rot[v]:
            d[w] = acc; acc += C[fi][pos].v.lo[col]
        phi.append(d)
    A = np.zeros((E, 2 * n + 1))
    for k, (u, w) in enumerate(S.edges):
        A[k, 2 * u] = -np.cos(phi[u][w]); A[k, 2 * u + 1] = -np.sin(phi[u][w])
        A[k, 2 * w] = -np.cos(phi[w][u]); A[k, 2 * w + 1] = -np.sin(phi[w][u])
        A[k, -1] = 1
    cost = np.zeros(2 * n + 1); cost[-1] = -1
    res = linprog(cost, A_ub=A, b_ub=np.zeros(E), bounds=[(-1, 1)] * (2 * n) + [(None, 1)], method='highs')
    if res.status != 0 or res.x[-1] < 1e-12: return None
    return res.x[:-1]


# ---------------------------------------------------------------------- driver
def charts_of(S):
    pents = [fi for fi, k in enumerate(S.sizes) if k == 5]
    for combo in itertools.product('abc', repeat=len(pents)):
        yield dict(zip(pents, combo))


def run(S, th_lo, th_hi, tol=1e-6, batch=512, max_boxes=10 ** 7, w_lp=0.06, w_gordan=0.15, lp_rounds=3,
        verbose=False, charts=None, time_limit=None, lp_mode='full'):
    """Branch and bound for one structure.  Widths are absolute (radians)."""
    stats = dict(boxes=0, invalid=0, vsum=0, gordan=0, gordan_inherit=0, lp=0, contract=0, contract_empty=0,
                 terminal=0)
    terminal = []
    t0 = time.time()
    for ch in (charts if charts is not None else charts_of(S)):
        M = Model(S, ch)
        lo0, hi0 = M.initial_box(th_lo, th_hi)
        Q = [(lo0, hi0, None)]
        while Q:
            Bk = Q[-batch:]; del Q[-batch:]
            lo = np.stack([q[0] for q in Bk], 1); hi = np.stack([q[1] for q in Bk], 1)
            ys = [q[2] for q in Bk]
            stats['boxes'] += lo.shape[1]
            C, F, J, ok, vok = evaluate(M, lo, hi)
            stats['invalid'] += int(np.sum(~ok)); stats['vsum'] += int(np.sum(ok & ~vok))
            alive = ok & vok
            # LP contraction for narrow boxes
            for b in np.nonzero(alive)[0]:
                if np.max(hi[:, b] - lo[:, b]) > w_lp: continue
                l, h = lo[:, b:b + 1], hi[:, b:b + 1]; Jb = J[(slice(None), slice(None), slice(b, b + 1))] \
                    if False else I(J.lo[:, :, b:b + 1], J.hi[:, :, b:b + 1])
                for rnd in range(lp_rounds):
                    w0 = h - l
                    r = lp_contract(M, l, h, Jb, 0, coords=[] if lp_mode == 'feas' else None); stats['contract'] += 1
                    if r is None: alive[b] = False; stats['contract_empty'] += 1; break
                    l, h = r[0][:, None], r[1][:, None]
                    _, Fb, Jb, okb, vokb = evaluate(M, l, h)
                    if not (okb[0] and vokb[0]): alive[b] = False; stats['contract_empty'] += 1; break
                    if np.max((w0 - (h - l)) / np.maximum(w0, 1e-300)) < 0.2: break
                lo[:, b], hi[:, b] = l[:, 0], h[:, 0]
            if not alive.any(): continue
            C, F, J, ok, vok = evaluate(M, lo, hi)
            alive &= ok & vok
            Dd = directions(M, C)
            Y = np.stack([y if y is not None else np.zeros(2 * S.n) for y in ys], 1)
            has = np.array([y is not None for y in ys])
            g = gordan_check(M, Dd, Y) & has & alive
            stats['gordan_inherit'] += int(g.sum()); alive &= ~g
            cand = [b for b in np.nonzero(alive)[0] if np.max(hi[:, b] - lo[:, b]) <= w_gordan]
            if cand:
                mids = 0.5 * (lo[:, cand] + hi[:, cand])
                Cm, okm = M.corners(mids, mids.copy(), box=(lo[:, cand], hi[:, cand]))
                for j, b in enumerate(cand):
                    yb = gordan_lp(M, None, Cm, j); stats['lp'] += 1
                    if yb is not None:
                        ys[b] = yb; Y[:, b] = yb; has[b] = True
            g = gordan_check(M, Dd, Y) & has & alive
            stats['gordan'] += int(g.sum()); alive &= ~g
            for b in np.nonzero(alive)[0]:
                w = hi[:, b] - lo[:, b]; k = int(np.argmax(w))
                if w[k] < tol:
                    terminal.append((ch, lo[:, b].copy(), hi[:, b].copy())); stats['terminal'] += 1; continue
                mid = 0.5 * (lo[k, b] + hi[k, b])
                l1, h1 = lo[:, b].copy(), hi[:, b].copy(); h1[k] = mid
                l2, h2 = lo[:, b].copy(), hi[:, b].copy(); l2[k] = mid
                Q.append((l1, h1, ys[b])); Q.append((l2, h2, ys[b]))
            if stats['boxes'] > max_boxes or (time_limit and time.time() - t0 > time_limit):
                stats['aborted'] = True; stats['queue'] = len(Q); stats['time'] = time.time() - t0
                return terminal, stats
        if verbose:
            print(ch, stats, '%.1fs' % (time.time() - t0), flush=True)
    stats['time'] = time.time() - t0
    return terminal, stats


# ---------------------------------------------------------------------- LP contraction
def lp_contract(M, lo, hi, J, b, coords=None):
    """Linear-relaxation contraction of box b (mean value form, all finite rows of the vertex equations).
    Returns (lo_b, hi_b) or None if no zero of the vertex equations lies in the box."""
    n = lo.shape[0]
    l, h = lo[:, b], hi[:, b]; m = 0.5 * (l + h)
    _, Fm, _, _, _ = evaluate(M, m[:, None], m[:, None].copy(), box=(l[:, None], h[:, None]))
    Jl, Jh = J.lo[:, :, b], J.hi[:, :, b]
    rows = np.all(np.isfinite(Jl) & np.isfinite(Jh), axis=1) & np.isfinite(Fm.lo[:, 0]) & np.isfinite(Fm.hi[:, 0])
    if not rows.any(): return l, h
    Jl, Jh = Jl[rows], Jh[rows]; flo, fhi = Fm.lo[rows, 0], Fm.hi[rows, 0]
    Jc = 0.5 * (Jl + Jh); Jr = 0.5 * (Jh - Jl)
    # for p in box: F(p) in [flo, fhi] + J(box)(p - m); with y = p - m, |y| <= w/2:
    #   J(box) y  subset of  Jc y + [-Jr |y|, Jr |y|]  subset of  Jc y + [-R, R],  R = Jr w/2
    R = Jr @ ((h - l) / 2)
    slack = 1e-12 * (1 + np.abs(flo) + np.abs(fhi) + R + np.abs(Jc).sum(1))
    # zero:  Jc y <= -flo + R ,  -Jc y <= fhi + R
    A = np.vstack([Jc, -Jc]); bb = np.concatenate([-flo + R + slack, fhi + R + slack])
    bounds = list(zip(l - m, h - m))
    nl, nh = l.copy(), h.copy()
    if coords is not None and len(coords) == 0:
        r0 = linprog(np.zeros(n), A_ub=A, b_ub=bb, bounds=bounds, method='highs')
        return None if r0.status == 2 else (l, h)
    for k in (range(n) if coords is None else coords):
        c = np.zeros(n); c[k] = 1
        r1 = linprog(c, A_ub=A, b_ub=bb, bounds=bounds, method='highs')
        if r1.status == 2: return None
        if r1.status != 0: continue
        r2 = linprog(-c, A_ub=A, b_ub=bb, bounds=bounds, method='highs')
        if r2.status != 0: continue
        pad = 1e-9 * (h[k] - l[k]) + 1e-13
        nl[k] = max(l[k], m[k] + r1.fun - pad); nh[k] = min(h[k], m[k] - r2.fun + pad)
        bounds[k] = (nl[k] - m[k], nh[k] - m[k])
    return nl, nh
