"""LP relaxation test for a candidate face structure.
Variables: alpha; per quad: D (corner at the endpoints of the deleted diagonal), O (other corners);
per pentagon: one variable per corner.  Constraints:
  vertex sums = 360; triangle corner = alpha; alpha <= D,O <= 2 alpha; 3 alpha <= D+O <= 3 alpha + S4;
  alpha <= pentagon corner <= 180; 9 alpha + S5L <= pentagon sum <= 9 alpha + S5U;
  alpha in [A60, 72].  Returns min alpha (None if infeasible)."""
import numpy as np
from scipy.optimize import linprog
A60 = float(np.degrees(np.arccos(1 / 3)))
S4 = 9.0; S5L = 10.0; S5U = 75.0

def lp_min_alpha(F, faces):
    """F: triangles of T; faces: list of tuples of triangle indices (1,2,3 triangles)."""
    var = {'alpha': 0}; nv = 1
    vrow = [dict() for _ in range(12)]       # vertex -> {var: coeff}
    cons_ub = []; cons_eq = []
    def addv(name):
        nonlocal nv; var[name] = nv; nv += 1; return nv - 1
    for j, tris in enumerate(faces):
        if len(tris) == 1:
            for v in F[tris[0]]: vrow[v][0] = vrow[v].get(0, 0) + 1
        elif len(tris) == 2:
            a, b = tris; sh = set(F[a]) & set(F[b]); others = (set(F[a]) | set(F[b])) - sh
            D = addv(('D', j)); O = addv(('O', j))
            for v in sh: vrow[v][D] = vrow[v].get(D, 0) + 1
            for v in others: vrow[v][O] = vrow[v].get(O, 0) + 1
            cons_ub += [({D: -1, 0: 1}, 0), ({D: 1, 0: -2}, 0), ({O: -1, 0: 1}, 0), ({O: 1, 0: -2}, 0),
                        ({D: -1, O: -1, 0: 3}, 0), ({D: 1, O: 1, 0: -3}, S4)]
        else:
            vs = set(v for t in tris for v in F[t])
            ids = []
            for v in vs:
                c = addv(('P', j, v)); ids.append(c); vrow[v][c] = vrow[v].get(c, 0) + 1
                cons_ub += [({c: -1, 0: 1}, 0), ({c: 1}, 180.0)]
            row = {c: -1 for c in ids}; row[0] = 9; cons_ub.append((row, -S5L))
            row = {c: 1 for c in ids}; row[0] = -9; cons_ub.append((row, S5U))
    A_ub = np.zeros((len(cons_ub), nv)); b_ub = np.zeros(len(cons_ub))
    for r, (row, b) in enumerate(cons_ub):
        for k, v in row.items(): A_ub[r, k] = v
        b_ub[r] = b
    A_eq = np.zeros((12, nv)); b_eq = np.full(12, 360.0)
    for v in range(12):
        for k, c in vrow[v].items(): A_eq[v, k] = c
    bounds = [(A60, 72.0)] + [(0, 180)] * (nv - 1)
    cost = np.zeros(nv); cost[0] = 1
    res = linprog(cost, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=b_eq, bounds=bounds, method='highs')
    return None if res.status != 0 else res.fun
