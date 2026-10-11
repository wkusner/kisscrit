"""Candidate stressed contact graphs for 12 points with minimum angle in (60, 63.4349] degrees.

Faces of a candidate are unions of 1, 2 or 3 triangles of a triangulation T (triangle, quadrilateral, pentagon;
a pentagon is a path of three triangles in the dual).  We enumerate the partitions of the 20 triangles of each T
into such faces by depth-first search, pruning as soon as a vertex has all its triangles assigned:
the corners at v are the faces containing v (a face counts once per vertex), and they must be able to sum to
360 degrees with triangle corners = alpha(theta) in [70.53, 72], quadrilateral corners in [alpha, 2 alpha],
pentagon corners in [alpha, 180], for a common theta (alpha is cut into cells; a mask of feasible cells is kept).
Degree-2 vertices need two straight pentagon corners.
"""
import json, sys, itertools, collections, time
import numpy as np
import networkx as nx

A60 = float(np.degrees(np.arccos(1 / 3))); A63 = 72.0
NC = 64
CELLS = [(A60 + (A63 - A60) * i / NC, A60 + (A63 - A60) * (i + 1) / NC) for i in range(NC)]
FULL = (1 << NC) - 1


def vmask(t, q, p):
    d = t + q + p
    if d < 2 or d > 5: return 0
    if d == 2: return FULL if (t == 0 and q == 0 and p == 2) else 0
    m = 0
    for i, (lo, hi) in enumerate(CELLS):
        if d * lo <= 360 + 1e-9 and t * hi + q * 2 * hi + p * 180 >= 360 - 1e-9:
            m |= 1 << i
    return m
VMASK = {(t, q, p): vmask(t, q, p) for t in range(6) for q in range(6) for p in range(6)}


def process_T(F, out, buckets):
    F = [tuple(sorted(t)) for t in F]; nF = len(F)
    adj = collections.defaultdict(list)          # triangle adjacency via shared edges
    edge_faces = collections.defaultdict(list)
    for fi, t in enumerate(F):
        for e in itertools.combinations(t, 2): edge_faces[e].append(fi)
    for e, (a, b) in edge_faces.items():
        adj[a].append(b); adj[b].append(a)
    vtris = collections.defaultdict(list)
    for fi, t in enumerate(F):
        for v in t: vtris[v].append(fi)
    # order triangles: greedy so vertices complete early
    order = []; seen = set()
    vorder = sorted(range(12), key=lambda v: -len(vtris[v]))
    # BFS by vertex completion
    cur = [vorder[0]]; done_v = set()
    while len(order) < nF:
        best = None
        for v in range(12):
            if v in done_v: continue
            rem = [f for f in vtris[v] if f not in seen]
            if best is None or len(rem) < best[0] or (len(rem) == best[0] and len(rem) > 0):
                if rem or best is None: best = (len(rem), v)
        v = best[1]; done_v.add(v)
        # add remaining triangles around v in cyclic-ish order
        for f in vtris[v]:
            if f not in seen: seen.add(f); order.append(f)
    pos = {f: i for i, f in enumerate(order)}
    # vertex completion step: after assigning triangles order[:k], vertex v complete when all its tris assigned
    vlast = {v: max(pos[f] for f in vtris[v]) for v in range(12)}
    face_of = [-1] * nF; faces = []          # faces: list of (tri tuple)
    count = 0

    def corners_at(v):
        fs = set(face_of[f] for f in vtris[v])
        t = q = p = 0
        for fid in fs:
            k = len(faces[fid]) + 2
            if k == 3: t += 1
            elif k == 4: q += 1
            else: p += 1
        return t, q, p

    def record():
        nonlocal count
        cl = [corners_at(v) for v in range(12)]
        H = nx.Graph()
        for v in range(12): H.add_node(('v', v), kind='v')
        flist = []
        for j, tris in enumerate(faces):
            vs = sorted(set(v for f in tris for v in F[f]))
            k = len(tris) + 2
            H.add_node(('f', j), kind='f%d' % k)
            for v in vs: H.add_edge(('v', v), ('f', j))
            flist.append(vs)
        key = (tuple(sorted(sum(c) for c in cl)), tuple(sorted(len(t) + 2 for t in faces)),
               nx.weisfeiler_lehman_graph_hash(H, node_attr='kind', iterations=4))
        nm = lambda a, b: a['kind'] == b['kind']
        for H2 in buckets[key]:
            if nx.is_isomorphic(H, H2, node_match=nm):
                return
        buckets[key].append(H)
        # kept edges = triangulation edges not interior to a merged face
        interior = set()
        for tris in faces:
            for a, b in itertools.combinations(tris, 2):
                sh = set(F[a]) & set(F[b])
                if len(sh) == 2: interior.add(tuple(sorted(sh)))
        edges = [list(e) for e in edge_faces if e not in interior]
        out.append(dict(faces=flist, sizes=[len(t) + 2 for t in faces], edges=edges, corners=cl))
        count += 1

    def check_vertices(k, mask):
        # vertices completed exactly at step k
        for v in range(12):
            if vlast[v] == k:
                c = corners_at(v)
                mask &= VMASK.get(c, 0)
                if not mask: return 0
        return mask

    def rec(k, mask):
        # advance to next unassigned triangle; check vertices completed along the way
        while k < nF and face_of[order[k]] != -1:
            mask = check_vertices(k, mask)
            if not mask: return
            k += 1
        if k == nF:
            record(); return
        f = order[k]
        options = [(f,)]
        for g in adj[f]:
            if face_of[g] == -1:
                options.append((f, g))
                for h in adj[g]:
                    if h != f and face_of[h] == -1 and h not in adj[f]:
                        options.append((f, g, h))
                for h in adj[f]:
                    if h != g and face_of[h] == -1 and h not in adj[g] and g < h:
                        options.append((g, f, h))
        for opt in options:
            # a 3-triangle face must be a pentagon: the three triangles share no common vertex
            if len(opt) == 3:
                common = set(F[opt[0]]) & set(F[opt[1]]) & set(F[opt[2]])
                if any(len(vtris[v]) == 3 for v in common):   # common vertex would be interior
                    continue
            fid = len(faces); faces.append(opt)
            for g in opt: face_of[g] = fid
            m2 = check_vertices(k, mask)
            if m2: rec(k + 1, m2)
            for g in opt: face_of[g] = -1
            faces.pop()

    rec(0, FULL)
    return count


if __name__ == '__main__':
    T = json.load(open('tri12.json'))
    lo, hi, outp = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
    out = []; buckets = collections.defaultdict(list); t0 = time.time()
    for k in range(lo, min(hi, len(T))):
        process_T(T[k], out, buckets)
        if k % 100 == 0:
            print(k, len(out), '%.0fs' % (time.time() - t0), flush=True)
    json.dump(out, open(outp, 'w'))
    print('total', len(out), '%.0fs' % (time.time() - t0))
