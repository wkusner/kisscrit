"""Python front end for the compiled branch and bound (fast.py): one Structure, all pentagon regions."""
import itertools, time
import numpy as np
from . import fast
from .model import Model


def arrays(S, chart):
    nF = len(S.sizes)
    fsize = np.array(S.sizes, np.int64)
    fp0 = np.full(nF, -1, np.int64); fp1 = np.full(nF, -1, np.int64); fmode = np.zeros(nF, np.int64)
    for fi, k in enumerate(S.sizes):
        if k >= 4: fp0[fi] = S.pidx[fi][0]
        if k == 5: fp1[fi] = S.pidx[fi][1]; fmode[fi] = 'abc'.index(chart[fi])
    vstart = [0]; vface = []; vpos = []; slot = {}
    for v in range(S.n):
        for (w, fi, pos) in S.rot[v]:
            slot[(v, w)] = len(vface); vface.append(fi); vpos.append(pos)
        vstart.append(len(vface))
    eu = np.array([u for u, w in S.edges], np.int64); ew = np.array([w for u, w in S.edges], np.int64)
    eus = np.array([slot[(u, w)] for u, w in S.edges], np.int64); ews = np.array([slot[(w, u)] for u, w in S.edges], np.int64)
    return (fsize, fp0, fp1, fmode, np.array(vstart, np.int64), np.array(vface, np.int64), np.array(vpos, np.int64),
            eu, eus, ew, ews)


def oriented(S):
    """Rotate each pentagon's boundary so that position 4 is a vertex of largest degree (then position 1)."""
    from .structure import Structure
    deg = {v: len(S.rot[v]) for v in range(S.n)}
    faces = []
    for f in S.faces:
        if len(f) == 5:
            # a degree-2 vertex (two straight corners) goes to position 4, a second one to position 1;
            # otherwise position 4 gets a vertex of largest degree
            def score(s):
                d4, d1 = deg[f[(s + 4) % 5]], deg[f[(s + 1) % 5]]
                return (d4 == 2, d1 == 2, d4, d1)
            best = max(range(5), key=score)
            f = [f[(best + j) % 5] for j in range(5)]
        faces.append(list(f))
    return Structure(faces)


def charts_of(S):
    """Pentagon regions needed.  A corner >= 180 - 25 deg can only sit at a vertex of degree <= 3, because the
    other corners there are each >= alpha > 70.52 deg and 155 + 3 * 70.52 > 360.  So region 'b' (corner at
    position 4 near straight) is needed only if that vertex has degree <= 3, and region 'c' only if the vertices
    at positions 4 and 1 both do."""
    pents = [fi for fi, k in enumerate(S.sizes) if k == 5]
    opts = []
    for fi in pents:
        f = S.faces[fi]; d4 = len(S.rot[f[4]]); d1 = len(S.rot[f[1]])
        if d4 == 2:            # corner at position 4 is exactly straight
            opts.append('c' if d1 == 2 else ('b' if d1 >= 4 else 'bc'))
        else:
            opts.append('a' if d4 >= 4 else ('ab' if d1 >= 4 else 'abc'))
    for combo in itertools.product(*opts):
        yield dict(zip(pents, combo))


def run(S, th_lo, th_hi, tol=1e-6, max_boxes=10 ** 8, w_lp=10.0, w_gordan=0.2, rounds=3, use_bounds=True,
        charts=None, time_limit=None):
    if charts is None: S = oriented(S)
    t0 = time.time(); tot = np.zeros(8, np.int64); terminal = []; aborted = False
    for ch in (charts if charts is not None else charts_of(S)):
        A = arrays(S, ch)
        lo0, hi0 = Model(S, ch).initial_box(th_lo, th_hi)
        for fi, mode in ch.items():          # straight corners at degree-2 vertices are fixed at pi
            if len(S.rot[S.faces[fi][4]]) == 2:
                k = S.pidx[fi][1] if mode == 'b' else S.pidx[fi][0]
                lo0[k] = hi0[k] = np.pi
        st, TL, TH, ab = fast.run_chart(*A, lo0, hi0, tol, max_boxes, w_lp, w_gordan, rounds, use_bounds)
        tot += st; terminal += [(ch, l, h) for l, h in zip(TL, TH)]
        if ab: aborted = True; break
        if time_limit and time.time() - t0 > time_limit: aborted = True; break
    names = ['boxes', 'invalid', 'vsum', 'linear', 'gordan_inherit', 'gordan', 'terminal', 'gordan_lp']
    stats = dict(zip(names, map(int, tot))); stats['time'] = time.time() - t0; stats['aborted'] = aborted
    stats['faces'] = [list(map(int, f)) for f in S.faces]
    return terminal, stats
