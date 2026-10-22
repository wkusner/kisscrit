#!/usr/bin/env python3
"""BALL-KILL: does one Krawczyk uniqueness ball per census point cover the `cond` terminal boxes around it?

Status: untrusted (exploration or tooling; no part of the proof).  NOT trusted code: an exploratory
measurement (Project prompt BALL-KILL v1.0, 2026-10-21).  The developing map and every extent are floating
point ("computed (exploratory)"); only the Krawczyk inclusion is interval arithmetic with outward rounding
("computed (interval)").  Nothing here is a proof or a certificate, and nothing here may enter the trusted base.

Steps (one command runs all of them):
  1. read every terminal leaf of results/completeness/certificates/*.jsonl (replay by kisscrit.core.certificate);
  2. develop terminal-box points into configurations of 12 points on S^2 (face model kisscrit.completeness.model,
     combinatorics kisscrit.completeness.structure), align them to the census configurations
     (data/n12_above1.json, refined as in scripts/certify_n12.py) and measure the extents m1 and m2;
  3. bisect the largest radius rho at which the Krawczyk operator of kisscrit/certify.py (its _F_rows, _J_rows,
     _choose_dropped, gauge and dropped rows, imported unchanged) maps the cube of radius rho into its interior;
  4. write the verdict table to results/balls/ball-kill-2026-10-21.json.

  python3 scripts/balls/ball_kill_test.py              # all steps
  python3 scripts/balls/ball_kill_test.py --smoke      # the E = 19 saddle (r = 1.0002112) only
"""
import argparse, glob, json, math, os, sys, time
from multiprocessing import Pool

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
CERTS = os.path.join(ROOT, 'results', 'completeness', 'certificates')
OUT = os.path.join(ROOT, 'results', 'balls', 'ball-kill-2026-10-21.json')
EV_TERMINAL = 7          # docs/certificate-format.md section 3: the only event byte meaning `terminal`


# ---------------------------------------------------------------------------------------------- step 1
def _terminal_leaves_of_file(path):
    """All terminal leaves of one certificate file: list of dicts (code, E, sizes, faces, chart, modes, lo, hi,
    shard, line, theta, tol).  Charts whose event stream has no byte 7 have no terminal leaf and are skipped."""
    from kisscrit.core.certificate import replay, unpack
    out = []; shard = os.path.basename(path)
    with open(path) as fh:
        for ln, line in enumerate(fh, 1):
            if '"events"' not in line:
                continue
            rec = json.loads(line); cert = rec['certificate']
            faces = cert['faces']; E = sum(len(f) for f in faces) // 2
            for ci, ch in enumerate(cert['charts']):
                if ch.get('not_run'):
                    continue
                if not np.any(unpack(ch['events'], np.uint8) == EV_TERMINAL):
                    continue
                for kind, lo, hi, _ in replay(ch):
                    if kind == 'terminal':
                        out.append(dict(code=rec['code'], E=E, sizes=[len(f) for f in faces], faces=faces,
                                        chart=ci, modes=ch['modes'], lo=lo.tolist(), hi=hi.tolist(),
                                        shard=shard, line=ln, theta=cert['theta'], tol=cert['tol']))
    return out


def step1(workers=2):
    files = sorted(glob.glob(os.path.join(CERTS, 'certificates-cond-*.jsonl')))
    with Pool(workers) as pool:
        parts = pool.map(_terminal_leaves_of_file, files)
    leaves = [x for p in parts for x in p]
    leaves.sort(key=lambda d: (d['shard'], d['line'], d['chart'], d['lo']))
    return leaves, len(files)


# ---------------------------------------------------------------------------------------------- census
SMOKE_R = 1.0002111636          # the E = 19 saddle of coindex 3 named by audit D for the killing test


def load_census(refine=True, only=None):
    """The 12 census types of data/n12_above1.json.  For nondegenerate types (sdim 1) the KKT midpoint is
    re-derived exactly as scripts/certify_n12.py does: kkt.solve_kkt, then exact.refine_kkt at 90 digits, mu from
    the balance equations.  DOD and M (sdim 9 and 2) keep the solve_kkt floats (their square system is singular)."""
    import mpmath as mp
    from kisscrit.kkt import solve_kkt
    from kisscrit.exact import refine_kkt
    D = json.load(open(os.path.join(ROOT, 'data', 'n12_above1.json')))
    out = []
    for t in D['types']:
        if only is not None and abs(t['r'] - only) > 1e-9:
            continue
        A = [tuple(p) for p in t['pairs']]; X = np.array(t['X'])
        Xs, c, s, nf = solve_kkt(X, A)
        P = dict(r=t['r'], k=t['k'], coindex=t['coindex'], sdim=t['sdim'], A=A, c=float(c),
                 X=np.array(Xs, float), solve_residual=float(nf), mid=None)
        if refine and t['sdim'] == 1:
            x, cc, sig, res = refine_kkt(Xs, A, s, dps=90)
            mp.mp.dps = 90
            C = [[sum(x[i][a] * x[j][a] for a in range(3)) for j in range(12)] for i in range(12)]
            mu = [sum(sig[e] * C[i][j if i == a else a] for e, (a, j) in enumerate(A) if i in (a, j))
                  for i in range(12)]
            P['mid'] = (x, cc, sig, mu)
            P['X'] = np.array([[float(v) for v in row] for row in x]); P['c'] = float(cc)
            P['refine_residual'] = float(res)
        P['theta'] = math.acos(P['c'])
        out.append(P)
    return out


# ---------------------------------------------------------------------------------------------- step 2
def structure_and_model(faces, modes):
    from kisscrit.completeness.structure import Structure
    from kisscrit.completeness.model import Model
    S = Structure(faces)
    return S, Model(S, {int(k): v for k, v in modes.items()})


def corners_at(M, P, box_lo, box_hi):
    """Face corners (midpoints of the point enclosures) and vertex-sum Jacobian at the points P (npar, B),
    with the clamps of the leaf box (box_lo, box_hi), each of shape (npar,) or (npar, B), as in Model.corners.
    Returns (C, F, J, ok): C[face][pos] arrays (B,), F = vertex sums - 2 pi (B, n), J (B, n, npar), ok (B,)."""
    B = P.shape[1]
    blo = box_lo if np.ndim(box_lo) == 2 else np.repeat(np.asarray(box_lo)[:, None], B, 1)
    bhi = box_hi if np.ndim(box_hi) == 2 else np.repeat(np.asarray(box_hi)[:, None], B, 1)
    Cd, ok = M.corners(P, P.copy(), box=(blo, bhi))
    C = [[0.5 * (c.v.lo + c.v.hi) for c in face] for face in Cd]
    V = M.vertex_sums(Cd)
    F = np.stack([0.5 * (v.v.lo + v.v.hi) for v in V], 1) - 2 * np.pi
    J = np.stack([0.5 * (v.g.lo + v.g.hi) for v in V], 1).transpose(2, 1, 0)        # (B, n, npar)
    return C, F, J, ok


def develop(S, C, theta, root=0):
    """Develop corners C (per face, arrays (B,)) and contact angle theta (B,) into points on S^2, by breadth-first
    search over the edges from `root` (root at e3, its first rotation neighbour in direction e1).  The edge v->w
    leaves v at angle phi[v][w] counterclockwise (about the outward normal x_v) from v's first neighbour.
    Returns (X (B, n, 3), closure): closure = max over directed edges v->w of |x_w predicted from v - x_w|."""
    n = S.n; B = theta.shape[0]
    phi = []
    for v in range(n):
        d = {}; acc = np.zeros(B)
        for (w, fi, pos) in S.rot[v]:
            d[w] = acc.copy(); acc = acc + C[fi][pos]
        phi.append(d)
    ct, st = np.cos(theta)[:, None], np.sin(theta)[:, None]
    X = np.full((n, B, 3), np.nan); T = np.full((n, B, 3), np.nan)
    X[root] = [0.0, 0.0, 1.0]; T[root] = [1.0, 0.0, 0.0]
    done = {root}; queue = [root]

    def step(v, w):
        p = phi[v][w][:, None]
        t = np.cos(p) * T[v] + np.sin(p) * np.cross(X[v], T[v])
        return ct * X[v] + st * t

    while queue:
        v = queue.pop(0)
        for (w, _, _) in S.rot[v]:
            if w in done:
                continue
            X[w] = step(v, w)
            u = (X[v] - ct * X[w]) / st
            p = phi[w][v][:, None]
            T[w] = np.cos(p) * u - np.sin(p) * np.cross(X[w], u)
            done.add(w); queue.append(w)
    clo = np.zeros(B)
    for v in range(n):
        for (w, _, _) in S.rot[v]:
            clo = np.maximum(clo, np.linalg.norm(step(v, w) - X[w], axis=1))
    return X.transpose(1, 0, 2), clo


def project_to_variety(M, P, box_lo, box_hi, fixed, iters=12, rcond=1e-7):
    """Damped minimal-norm Gauss-Newton on the vertex sums: move the points P (npar, B) towards
    {vertex sums = 2 pi}.  `fixed` (bool (npar, B)) marks coordinates that are not moved (straight corners).
    Near a critical point the vertex-sum Jacobian is nearly rank-deficient (the self-stress is a dependency), so
    the pseudo-inverse uses a relative cut-off `rcond` and a step is halved until |F| decreases (else not taken).
    Only the columns still moving are re-evaluated.  Returns (P', |F|inf before, |F|inf after, |P' - P|inf)."""
    B = P.shape[1]
    blo = box_lo if np.ndim(box_lo) == 2 else np.repeat(np.asarray(box_lo)[:, None], B, 1)
    bhi = box_hi if np.ndim(box_hi) == 2 else np.repeat(np.asarray(box_hi)[:, None], B, 1)
    fixed = fixed if np.ndim(fixed) == 2 else np.repeat(np.asarray(fixed)[:, None], B, 1)
    P0 = P.copy(); P = P.copy()
    _, F, J, _ = corners_at(M, P, blo, bhi)
    nF = np.where(np.isfinite(F).all(1), np.abs(F).max(1), np.inf); F0 = nF.copy()
    for _ in range(iters):
        act = np.where((nF > 1e-14) & np.isfinite(J).all((1, 2)) & np.isfinite(nF))[0]
        if not len(act): break
        Ja = J[act] * ~fixed.T[act][:, None, :]
        dp = -np.einsum('bij,bj->bi', np.linalg.pinv(Ja, rcond=rcond), F[act]).T     # (npar, len(act))
        lam = 1.0; todo = act; moved = 0
        for _ in range(8):
            Pt = P[:, todo] + lam * dp[:, np.searchsorted(act, todo)]
            _, Ft, Jt, _ = corners_at(M, Pt, blo[:, todo], bhi[:, todo])
            nt = np.where(np.isfinite(Ft).all(1), np.abs(Ft).max(1), np.inf)
            good = nt < nF[todo]; g = todo[good]
            P[:, g] = Pt[:, good]; F[g] = Ft[good]; J[g] = Jt[good]; nF[g] = nt[good]; moved += len(g)
            todo = todo[~good]
            if not len(todo): break
            lam *= 0.5
        if not moved: break
    return P, F0, nF, np.abs(P - P0).max(0)


def kabsch(Y, X0, mirror=False):
    """Best proper rotation of Y (B, n, 3) (mirrored in z first if `mirror`) onto X0 (n, 3).  The optimum
    satisfies sum_i x0_i x (R y_i) = 0, which is the gauge of kisscrit/certify.py.  Returns R Y."""
    if mirror:
        Y = Y * np.array([1.0, 1.0, -1.0])
    H = np.einsum('bni,nj->bij', Y, X0)
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(np.einsum('bij,bjk->bik', U, Vt)))
    Dg = np.zeros((len(Y), 3, 3)); Dg[:, 0, 0] = 1; Dg[:, 1, 1] = 1; Dg[:, 2, 2] = d
    R = np.einsum('bji,bjk,bkl->bil', Vt, Dg, U.transpose(0, 2, 1))      # V D U^T
    return np.einsum('bij,bnj->bni', R, Y)


def m1_m2(Yal, theta, P):
    """m1: max over points of the angle to the census image; m2: sup norm of (X, c) - (X0, c0)."""
    Yn = Yal / np.linalg.norm(Yal, axis=2, keepdims=True)
    X0 = np.broadcast_to(P['X'], Yn.shape)              # atan2(|x y|, x.y): accurate for small angles
    m1 = np.arctan2(np.linalg.norm(np.cross(Yn, X0), axis=2), np.einsum('bnj,bnj->bn', Yn, X0)).max(1)
    m2 = np.maximum(np.abs(Yal - P['X']).max((1, 2)), np.abs(np.cos(theta) - P['c']))
    return m1, m2


def _graph(edges):
    import networkx as nx
    G = nx.Graph(); G.add_nodes_from(range(12)); G.add_edges_from(edges); return G


def identify_code(S, Y, theta, census, cap=200000):
    """Best census point, label map and orientation for one developed configuration Y (n, 3) of a code.
    Tries every subgraph monomorphism of the code graph into each census contact graph with k >= E."""
    from networkx.algorithms import isomorphism
    Gc = _graph(S.edges); res = []
    for P in census:
        if P['k'] < len(S.edges):
            continue
        GM = isomorphism.GraphMatcher(_graph(P['A']), Gc)
        best = None; cnt = 0
        for m in GM.subgraph_monomorphisms_iter():
            perm = np.array([m[i] for i in range(12)])
            for mirror in (False, True):
                Yal = kabsch(Y[perm][None], P['X'], mirror)
                m1, _ = m1_m2(Yal, np.array([theta]), P)
                if best is None or m1[0] < best[0]:
                    best = (float(m1[0]), perm, mirror)
            cnt += 1
            if cnt >= cap: break
        if best is not None:
            res.append((best[0], P, best[1], best[2], cnt))
    res.sort(key=lambda t: t[0])
    return res


def census_face_parameters(S, M, modes, Xc):
    """Face parameters of a census configuration Xc (labels of the code): theta, the quadrilateral corner at
    boundary position 0, and the pentagon chart corners (region a: positions 1, 2; b: 3, 4; c: 4, 0)."""
    def corner(f, pos):
        v, a, b = f[pos], f[pos - 1], f[(pos + 1) % len(f)]
        x = Xc[v]; ta = Xc[a] - (Xc[a] @ x) * x; tb = Xc[b] - (Xc[b] @ x) * x
        return math.atan2(np.cross(tb, ta) @ x, tb @ ta) % (2 * np.pi)
    th = float(np.mean([math.acos(np.clip(Xc[i] @ Xc[j], -1, 1)) for (i, j) in S.edges]))
    p = np.zeros(S.npar); p[0] = th
    pos_of = {'a': (1, 2), 'b': (3, 4), 'c': (4, 0)}
    for fi, f in enumerate(S.faces):
        if len(f) == 4: p[S.pidx[fi][0]] = corner(f, 0)
        if len(f) == 5:
            q1, q2 = pos_of[M.charts[fi]]
            p[S.pidx[fi][0]] = corner(f, q1); p[S.pidx[fi][1]] = corner(f, q2)
    return p


def box_samples(lo, hi, rng, max_corners=16):
    """Centre plus all corners (if 2^n' <= max_corners, n' = number of coordinates of positive width) or
    max_corners seeded random corners.  Returns (P (npar, B), description)."""
    lo = np.asarray(lo); hi = np.asarray(hi); n = len(lo)
    free = [i for i in range(n) if hi[i] > lo[i]]
    pts = [0.5 * (lo + hi)]
    if 2 ** len(free) <= max_corners:
        for bits in range(2 ** len(free)):
            p = lo.copy()
            for t, i in enumerate(free):
                if bits >> t & 1: p[i] = hi[i]
            pts.append(p)
        how = 'all %d corners' % 2 ** len(free)
    else:
        for _ in range(max_corners):
            p = lo.copy(); s = rng.integers(0, 2, len(free))
            for t, i in enumerate(free):
                if s[t]: p[i] = hi[i]
            pts.append(p)
        how = '%d random corners' % max_corners
    return np.array(pts).T, how


def _fixed_mask(lo, hi):
    """Coordinates of zero width (straight corners, [pi, pi] or [pi, pi+]) are held fixed by the push."""
    return (np.asarray(hi) - np.asarray(lo)) < 1e-14


def step2(leaves, census, max_corners=16, seed=20261021, chunk=8192, near=10.0, log=print):
    """Per code: identify the census point from the first box's centre; per box: extents m1 and m2 over the
    centre and corners (box_samples), for (raw) the spanning-tree development of the sample itself and (push) the
    development of the sample pushed onto {vertex sums = 2 pi}.  Pushed samples are used only if they converged
    (|F| <= 1e-12) and moved by at most `near` box widths.  Box extent = max of raw and the usable pushed values."""
    from collections import defaultdict
    rng = np.random.default_rng(seed)
    by = defaultdict(list)
    for d in leaves: by[d['code']].append(d)
    codes = []; boxes = []
    for code, L in sorted(by.items(), key=lambda kv: -len(kv[1])):
        t0 = time.time(); d0 = L[0]
        S, M = structure_and_model(d0['faces'], d0['modes'])
        lo = np.array(d0['lo']); hi = np.array(d0['hi']); c0 = 0.5 * (lo + hi)
        Pc, _, _, _ = project_to_variety(M, c0[:, None], lo, hi, _fixed_mask(lo, hi)[:, None])
        C, F, _, _ = corners_at(M, Pc, lo, hi)
        Y, clo = develop(S, C, Pc[0])
        ids = identify_code(S, Y[0], Pc[0, 0], census)
        if not ids:
            log('code %s: no census graph contains it' % code[:30]); continue
        m1b, P, perm, mirror, nmono = ids[0]
        # self-check: the census configuration's own face parameters (in this chart) must develop back onto it
        inv = np.empty(12, int); inv[perm] = np.arange(12)
        Xc = P['X'][inv] * (np.array([1.0, 1.0, -1.0]) if mirror else 1.0)
        pc = census_face_parameters(S, M, d0['modes'], Xc)
        Cc, Fc, _, _ = corners_at(M, pc[:, None], pc, pc)
        Yc, clc = develop(S, Cc, pc[:1])
        sc1, sc2 = m1_m2(kabsch(Yc[:, perm], P['X'], mirror), pc[:1], P)
        codes.append(dict(census_selfcheck=dict(closure=float(clc[0]), vsum=float(np.abs(Fc).max()),
                                                m1=float(sc1[0]), m2=float(sc2[0])),code=code, E=d0['E'], sizes=d0['sizes'], n_boxes=len(L),
                          theta_lo=min(x['lo'][0] for x in L), theta_hi=max(x['hi'][0] for x in L),
                          shard=d0['shard'], line=d0['line'], census_r=P['r'], census_k=P['k'],
                          mirror=bool(mirror), extra_contacts=P['k'] - d0['E'], m1_centre_first_box=m1b,
                          monomorphisms=nmono,
                          runner_up=(dict(r=ids[1][1]['r'], m1=ids[1][0]) if len(ids) > 1 else None),
                          centre_closure_after_push=float(clo[0])))
        # all samples of all boxes of this code, with per-column clamp boxes
        cols = []; owner = []; hows = []
        for ib, d in enumerate(L):
            Ps, how = box_samples(d['lo'], d['hi'], rng, max_corners)
            cols.append(Ps); owner += [ib] * Ps.shape[1]; hows.append(how)
        Pall = np.concatenate(cols, 1); owner = np.array(owner)
        BLO = np.array([L[i]['lo'] for i in owner]).T; BHI = np.array([L[i]['hi'] for i in owner]).T
        W = (BHI - BLO).max(0); FIX = _fixed_mask(BLO, BHI)
        res = {k: np.empty(Pall.shape[1]) for k in ('r1', 'r2', 'p1', 'p2', 'clr', 'clp', 'F0', 'F1', 'push')}
        okr = np.empty(Pall.shape[1], bool)
        for s0 in range(0, Pall.shape[1], chunk):
            sl = slice(s0, s0 + chunk); Ps = Pall[:, sl]
            Craw, _, _, ok = corners_at(M, Ps, BLO[:, sl], BHI[:, sl]); okr[sl] = ok
            Yr, res['clr'][sl] = develop(S, Craw, Ps[0])
            Pp, res['F0'][sl], res['F1'][sl], res['push'][sl] = project_to_variety(M, Ps, BLO[:, sl], BHI[:, sl], FIX[:, sl])
            Cp, _, _, _ = corners_at(M, Pp, BLO[:, sl], BHI[:, sl])
            Yp, res['clp'][sl] = develop(S, Cp, Pp[0])
            res['r1'][sl], res['r2'][sl] = m1_m2(kabsch(Yr[:, perm], P['X'], mirror), Ps[0], P)
            res['p1'][sl], res['p2'][sl] = m1_m2(kabsch(Yp[:, perm], P['X'], mirror), Pp[0], P)
        usable = (res['F1'] <= 1e-12) & (res['push'] <= near * W) & np.isfinite(res['p1'])
        for ib, d in enumerate(L):
            j = owner == ib; u = j & usable
            r1, r2 = res['r1'][j].max(), res['r2'][j].max()
            p1 = res['p1'][u].max() if u.any() else -np.inf; p2 = res['p2'][u].max() if u.any() else -np.inf
            boxes.append(dict(code=code, census_r=P['r'], E=d['E'], shard=d['shard'], line=d['line'],
                              chart=d['chart'], lo=d['lo'], hi=d['hi'], samples=hows[ib],
                              width=float(np.max(np.array(d['hi']) - np.array(d['lo']))),
                              m1_raw=float(r1), m2_raw=float(r2),
                              m1_push=float(p1) if u.any() else None, m2_push=float(p2) if u.any() else None,
                              m1=float(max(r1, p1)), m2=float(max(r2, p2)),
                              m1_centre_raw=float(res['r1'][np.where(j)[0][0]]),
                              closure_raw_max=float(res['clr'][j].max()), vsum_raw_max=float(res['F0'][j].max()),
                              closure_push_max=float(res['clp'][u].max()) if u.any() else None,
                              vsum_push_max=float(res['F1'][u].max()) if u.any() else None,
                              push_max=float(res['push'][u].max()) if u.any() else None,
                              pushed_usable=int(u.sum()), pushed_rejected=int((j & ~usable).sum()),
                              invalid_samples=int((~okr[j]).sum())))
        log('  code E=%d boxes=%d -> r=%.10f%s (m1 %.2e; runner-up %s) %.1f s' % (
            d0['E'], len(L), P['r'], ' mirror' if mirror else '', m1b,
            '%.10f %.2e' % (ids[1][1]['r'], ids[1][0]) if len(ids) > 1 else '-', time.time() - t0))
    return codes, boxes


# ---------------------------------------------------------------------------------------------- step 3
def _to_arb(x):
    """Exact conversion of an mpmath mpf to an Arb ball of radius 0 (no rounding through a double)."""
    import mpmath as mp, flint
    s, m, e, _ = mp.mpf(x)._mpf_
    if m == 0: return flint.arb(0)
    a = flint.arb(flint.fmpz(m)) * flint.arb(2) ** e
    return -a if s else a


class KrawczykAt:
    """The Krawczyk test of kisscrit/certify.py (lines `# Krawczyk`), at a variable radius rho, in Arb.
    Same square system (_F_rows, _J_rows), same gauge reference X0 (midpoint rounded to 30 digits), same dropped
    rows (_choose_dropped), same Newton polish and midpoint inverse C.  Only the enclosures change from mp.iv to
    Arb (python-flint), which is outward rounded as well; matrix products use arb_mat."""

    def __init__(self, x, cc, sig, mu, A, dps=80):
        import mpmath as mp, flint
        from kisscrit.certify import _F_rows, _J_rows, _choose_dropped
        self._F, self._J = _F_rows, _J_rows
        mp.mp.dps = dps; flint.ctx.dps = dps
        N = len(x); k = len(A); n = 3 * N + 1 + k + N
        zm = [mp.mpf(v) for row in x for v in row] + [mp.mpf(cc)] + [mp.mpf(v) for v in sig] + [mp.mpf(v) for v in mu]
        xm = [[zm[3 * i + t] for t in range(3)] for i in range(N)]
        X0 = [[mp.mpf(mp.nstr(xm[i][t], 30)) for t in range(3)] for i in range(N)]
        i0, p, q, i1, s_ = _choose_dropped(xm, N, A)
        drop = {N + k + 3 * i0 + p, N + k + 3 * i0 + q, N + k + 3 * i1 + s_}
        keep = [r for r in range(N + k + 3 * N + 1 + 3) if r not in drop]
        assert len(keep) == n
        for _ in range(8):
            Fm = mp.matrix(_F_rows(zm, N, A, X0, keep)); Jm = mp.matrix(_J_rows(zm, N, A, X0, keep, mp.mpf(0)))
            dz = mp.lu_solve(Jm, -Fm); zm = [zm[a] + dz[a] for a in range(n)]
            if mp.norm(dz) < mp.mpf(10) ** (-dps + 8): break
        Jm = mp.matrix(_J_rows(zm, N, A, X0, keep, mp.mpf(0)))
        C = mp.inverse(Jm)
        self.N, self.k, self.n, self.A, self.keep = N, k, n, A, keep
        self.zm = zm; self.X0 = X0; self.dropped = dict(i0=i0, comps=(p, q), i1=i1, comp=s_)
        self.zi = [_to_arb(v) for v in zm]
        self.X0i = [[_to_arb(v) for v in row] for row in X0]
        self.Ci = flint.arb_mat([[_to_arb(C[a, b]) for b in range(n)] for a in range(n)])
        Fp = _F_rows(self.zi, N, A, self.X0i, keep)
        self.CF = self.Ci * flint.arb_mat([[f] for f in Fp])
        self.residual = max(float(abs(f).mid()) + float(f.rad()) for f in Fp)

    def accepts(self, rho):
        """True iff K(Y) lies in the interior of Y, Y = the cube of radius rho about the polished midpoint."""
        import flint
        n = self.n
        Y = [z + flint.arb(0, rho) for z in self.zi]        # radius rounded up: contains [z - rho, z + rho]
        JY = flint.arb_mat(self._J(Y, self.N, self.A, self.X0i, self.keep, flint.arb(0)))
        Mi = flint.arb_mat(n, n)
        for a in range(n): Mi[a, a] = 1
        Mi = Mi - self.Ci * JY
        Dv = flint.arb_mat([[y - z] for y, z in zip(Y, self.zi)])
        K = flint.arb_mat([[z] for z in self.zi]) - self.CF + Mi * Dv
        return all(Y[a].contains_interior(K[a, 0]) for a in range(n))

    def bisect(self, lo=1e-40, hi=1.0, rel=1.01, maxit=200):
        """Largest accepted rho (log bisection).  Returns (accepted, rejected, evaluations)."""
        assert self.accepts(lo), 'Krawczyk fails at the starting radius'
        if self.accepts(hi): return hi, None, 2
        it = 2
        while hi / lo > rel and it < maxit:
            m = math.sqrt(lo * hi); it += 1
            if self.accepts(m): lo = m
            else: hi = m
        return lo, hi, it


def step3(census, log=print):
    out = []
    for P in census:
        row = dict(r=P['r'], k=P['k'], coindex=P['coindex'], sdim=P['sdim'])
        if P['sdim'] != 1:
            row.update(status='not applicable', reason='stress space of dimension %d: the square KKT system '
                       '(one stress normalisation) is singular at this point (certified_n12.json: "matrix is '
                       'numerically singular")' % P['sdim'])
            out.append(row); continue
        t0 = time.time()
        Kz = KrawczykAt(*P['mid'], P['A'])
        acc, rej, it = Kz.bisect()
        row.update(status='computed (interval)', rho_accepted=acc, rho_rejected=rej, evaluations=it,
                   midpoint_residual_enclosure=Kz.residual, dropped=Kz.dropped, n=Kz.n,
                   radius_m1=acc, m1_conversion='m1-ball of radius rho (Kabsch gauge) lies in the X-part of the '
                   'z-cube: |x_i - x0_i|_inf <= 2 sin(phi/2) <= phi', seconds=round(time.time() - t0, 1))
        log('r=%.10f k=%d coindex=%d: rho accepted %.4g, rejected %s (%d evals, %.0f s)'
            % (P['r'], P['k'], P['coindex'], acc, '%.4g' % rej if rej else '-', it, time.time() - t0))
        out.append(row)
    return out


# ---------------------------------------------------------------------------------------------- step 4
def verdict(ratio):
    if ratio is None: return 'not applicable (D1(c))'
    if ratio >= 2: return 'one ball covers the cluster'
    if ratio >= 1: return 'marginal'
    return 'one ball does not cover'


def step4(codes, boxes, radii, census):
    rad = {round(x['r'], 9): x for x in radii}
    from collections import defaultdict
    cl = defaultdict(list)
    for b in boxes: cl[round(b['census_r'], 9)].append(b)
    table = []
    for P in sorted(census, key=lambda P: P['r']):
        key = round(P['r'], 9); B = cl.get(key, [])
        R = rad.get(key, {})
        m1 = np.array([b['m1'] for b in B]) if B else None
        m2 = np.array([b['m2'] for b in B]) if B else None
        row = dict(r=P['r'], k=P['k'], coindex=P['coindex'], sdim=P['sdim'], n_boxes=len(B),
                   codes=sorted({b['code'] for b in B}),
                   extent_m1_max=float(m1.max()) if B else None, extent_m1_median=float(np.median(m1)) if B else None,
                   extent_m2_max=float(m2.max()) if B else None, extent_m2_median=float(np.median(m2)) if B else None,
                   radius_m1=R.get('radius_m1'), rho_accepted=R.get('rho_accepted'), rho_rejected=R.get('rho_rejected'))
        if not B:
            row.update(ratio_m1=None, ratio_m2=None, verdict='no terminal boxes')
        elif P['sdim'] != 1:
            row.update(ratio_m1=None, ratio_m2=None, verdict=verdict(None))
        else:
            row['ratio_m1'] = R['radius_m1'] / row['extent_m1_max']
            row['ratio_m2'] = R['rho_accepted'] / row['extent_m2_max']
            row['verdict'] = verdict(row['ratio_m1'])
            row['verdict_m2'] = verdict(row['ratio_m2'])
        if P['sdim'] == 9:
            ex = defaultdict(list)
            for b in B: ex[30 - b['E']].append(b['m1'])
            row['by_extra_contacts'] = {str(e): dict(boxes=len(v), m1_max=float(max(v))) for e, v in sorted(ex.items())}
        table.append(row)
    covered = sum(r['n_boxes'] for r in table if r['verdict'] == 'one ball covers the cluster')
    return table, covered


def versions():
    import mpmath, scipy, flint, networkx
    return dict(python=sys.version.split()[0], numpy=np.__version__, scipy=scipy.__version__,
                mpmath=mpmath.__version__, python_flint=flint.__version__, networkx=networkx.__version__)


def git_head():
    import subprocess
    try:
        h = subprocess.run(['git', '-C', ROOT, 'rev-parse', 'HEAD'], capture_output=True, text=True).stdout.strip()
        dirty = subprocess.run(['git', '-C', ROOT, 'status', '--porcelain', 'scripts/balls'],
                               capture_output=True, text=True).stdout.strip()
        return h, bool(dirty)
    except OSError:
        return None, None


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--smoke', action='store_true', help='the E = 19 saddle r = 1.0002112 only')
    ap.add_argument('--workers', type=int, default=2)
    ap.add_argument('--max-corners', type=int, default=16)
    ap.add_argument('--out', default=None, help='output JSON (default: results/balls/ball-kill-2026-10-21.json; '
                    'with --smoke, nothing is written unless --out is given)')
    a = ap.parse_args()
    t0 = time.time(); head, dirty = git_head()
    leaves, nfiles = step1(a.workers)
    print('step 1: %d terminal leaves, %d codes, %d files (%.1f s)' % (len(leaves), len({d['code'] for d in leaves}),
                                                                      nfiles, time.time() - t0))
    census = load_census(only=SMOKE_R if a.smoke else None)
    if a.smoke:
        Pk = census[0]
        leaves = [d for d in leaves if d['E'] == Pk['k'] and d['lo'][0] <= Pk['theta'] + 1e-5
                  and d['hi'][0] >= Pk['theta'] - 1e-5]
    codes, boxes = step2(leaves, census, a.max_corners)
    print('step 2: %d boxes developed (%.1f s)' % (len(boxes), time.time() - t0))
    radii = step3(census)
    table, covered = step4(codes, boxes, radii, census)
    for r in table:
        print('r=%.10f k=%d ci=%d boxes=%4d ext_m1=%s rad_m1=%s ratio=%s %s' % (
            r['r'], r['k'], r['coindex'], r['n_boxes'], '%.3g' % r['extent_m1_max'] if r['extent_m1_max'] else '-',
            '%.3g' % r['radius_m1'] if r['radius_m1'] else '-', '%.3g' % r['ratio_m1'] if r.get('ratio_m1') else '-',
            r['verdict']))
    out = dict(
        task='BALL-KILL', prompt_version='1.0 (2026-10-21)', script='scripts/balls/ball_kill_test.py',
        commit=head, script_dirty=dirty, command=' '.join(['python3'] + sys.argv), versions=versions(),
        seconds=round(time.time() - t0, 1),
        status_of_claims=dict(
            step1='computed (replay of the committed cond certificates by kisscrit.core.certificate.replay)',
            step2='computed (exploratory): floating-point developing map, IEEE double; not a proof',
            step3='computed (interval): Krawczyk inclusion in Arb (python-flint), outward rounded; the midpoint '
                  'polish and preconditioner are non-rigorous as Krawczyk allows; not a certificate',
            step4='computed (exploratory): ratios of an interval radius to floating-point extents',
            nothing_certified=True),
        definitions=dict(
            m1='max over the 12 points of the angle (radians) between a point and its census image, after the best '
               'proper rotation (Kabsch; mirror first if the code is the mirror image)',
            m2='sup norm of (X, c) - (X0, c0) in the z-coordinates of kisscrit/certify.py after its gauge '
               'sum x0_i x x_i = 0 (the Kabsch optimum satisfies it); sigma and mu are not defined at a '
               'non-critical box point and are not included',
            extent='per box: max over the samples of the larger of the raw development and the development of '
                   'the sample pushed onto the vertex-sum variety; per cluster: max and median over boxes',
            samples='box centre plus all corners when 2^(free coordinates) <= %d, else %d seeded random corners '
                    '(seed 20261021)' % (a.max_corners, a.max_corners),
            radius='largest rho for which the Krawczyk operator maps the cube of radius rho (in z = (X, c, sigma, '
                   'mu)) into its interior; bisection bracket [accepted, rejected], ratio 1.01'),
        terminal_totals=dict(leaves=len(leaves), codes=len({d['code'] for d in leaves})),
        codes=codes, radii=radii, verdict_table=table, boxes_covered_by_single_balls=covered,
        boxes=[b | dict(theta=[b['lo'][0], b['hi'][0]]) for b in boxes])
    path = a.out or (None if a.smoke else OUT)
    if path:
        os.makedirs(os.path.dirname(path), exist_ok=True)
        json.dump(out, open(path, 'w'), indent=1)
        print('wrote', os.path.relpath(path, ROOT))
    return out


if __name__ == '__main__':
    main()
