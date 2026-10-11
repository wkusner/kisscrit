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


def charts_of(S):
    pents = [fi for fi, k in enumerate(S.sizes) if k == 5]
    for combo in itertools.product('abc', repeat=len(pents)):
        yield dict(zip(pents, combo))


def run(S, th_lo, th_hi, tol=1e-6, max_boxes=10 ** 8, w_lp=10.0, w_gordan=0.2, rounds=3, use_bounds=True,
        charts=None, time_limit=None):
    t0 = time.time(); tot = np.zeros(8, np.int64); terminal = []; aborted = False
    for ch in (charts if charts is not None else charts_of(S)):
        A = arrays(S, ch)
        lo0, hi0 = Model(S, ch).initial_box(th_lo, th_hi)
        st, TL, TH, ab = fast.run_chart(*A, lo0, hi0, tol, max_boxes, w_lp, w_gordan, rounds, use_bounds)
        tot += st; terminal += [(ch, l, h) for l, h in zip(TL, TH)]
        if ab: aborted = True; break
        if time_limit and time.time() - t0 > time_limit: aborted = True; break
    names = ['boxes', 'invalid', 'vsum', 'linear', 'gordan_inherit', 'gordan', 'terminal', 'gordan_lp']
    stats = dict(zip(names, map(int, tot))); stats['time'] = time.time() - t0; stats['aborted'] = aborted
    return terminal, stats
