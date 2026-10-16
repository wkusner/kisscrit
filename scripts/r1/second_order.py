"""Second-order isolation test at r=1 balanced configurations (FCC, HCP).
For C with contact graph G (pairs at 60 deg):
  K = cone of nonneg stresses (extreme rays via nullspace + LP vertex enumeration),
  Flex = first-order flexes of G on S^2 modulo rotations,
  CC   = critical cone: flexes of supp(w) with other G-edges nondecreasing (here = Flex(G) if a full-support stress exists),
  Q_w  = (cD - W) x I_3 restricted to tangent space,
report eigenvalues of Q_w on Flex(G)/SO(3) for extreme rays w and the uniform-combination stress.

Status: untrusted (floating-point exploration; the exact version is exact_lambda.py).

"""
import numpy as np, itertools
from scipy.optimize import linprog
from scipy.linalg import null_space

def cubocta():
    P = [p for p in itertools.product([-1, 0, 1], repeat=3) if sum(abs(x) for x in p) == 2]
    return np.array(P, float) / np.sqrt(2)

def anticubocta():
    # triangular orthobicupola: hexagon at equator, triangles above/below aligned (eclipsed)
    t = np.arange(6) * np.pi / 3
    hexa = np.c_[np.cos(t), np.sin(t), np.zeros(6)]
    phi = np.arctan(np.sqrt(2))              # latitude 54.74 deg, so distance to hexagon = 60 deg
    tt = np.radians([30, 150, 270])
    top = np.c_[np.cos(phi) * np.cos(tt), np.cos(phi) * np.sin(tt), np.sin(phi) * np.ones(3)]
    bot = np.c_[np.cos(phi) * np.cos(tt), np.cos(phi) * np.sin(tt), -np.sin(phi) * np.ones(3)]
    return np.vstack([hexa, top, bot])

def analyze(name, U):
    N = len(U); G = np.abs(U @ U.T - 0.5) < 1e-9
    E = [(i, j) for i in range(N) for j in range(i + 1, N) if G[i, j]]
    print('==', name, 'N', N, 'contacts', len(E), 'min pair angle %.3f' % np.degrees(np.arccos(max(U[i] @ U[j] for i in range(N) for j in range(i + 1, N)))))
    # tangent bases
    T = []
    for u in U:
        a = np.cross(u, [1, 0, 0]) if abs(u[0]) < 0.9 else np.cross(u, [0, 1, 0]); a /= np.linalg.norm(a); T.append(np.c_[a, np.cross(u, a)])
    # rigidity matrix on S^2 (tangent coords): d/dt <u_i,u_j>
    R = np.zeros((len(E), 2 * N))
    for k, (i, j) in enumerate(E):
        R[k, 2 * i:2 * i + 2] = T[i].T @ U[j]; R[k, 2 * j:2 * j + 2] = T[j].T @ U[i]
    # stresses: left kernel of R, then nonneg part
    S = null_space(R.T); print('  stress space dim', S.shape[1], ' rank R', np.linalg.matrix_rank(R), ' flex dim (incl. 3 rotations)', 2 * N - np.linalg.matrix_rank(R))
    # extreme rays of {w >= 0, w = S a}: enumerate via LP vertices on sum w = 1
    rays = []
    for sub in itertools.combinations(range(len(E)), S.shape[1] - 1):
        A = S[list(sub), :]
        ns = null_space(A)
        if ns.shape[1] != 1: continue
        w = S @ ns[:, 0]
        if (w >= -1e-9).all() or (w <= 1e-9).all():
            w = np.abs(w); w /= w.sum()
            if not any(np.allclose(w, r, atol=1e-7) for r in rays): rays.append(w)
    print('  extreme rays of stress cone:', len(rays), ' supports:', [int((r > 1e-9).sum()) for r in rays])
    # rotations in tangent coords
    rot = np.zeros((2 * N, 3))
    for a in range(3):
        e = np.zeros(3); e[a] = 1
        for i in range(N): rot[2 * i:2 * i + 2, a] = T[i].T @ np.cross(e, U[i])
    F = null_space(R)                       # flexes incl rotations
    # remove rotations
    Fn = F - rot @ np.linalg.lstsq(rot, F, rcond=None)[0]
    Fn = null_space(np.vstack([R, rot.T]))  # flexes orthogonal to rotations
    print('  nontrivial flex dim', Fn.shape[1])
    def Q_of(w):
        W = np.zeros((N, N))
        for k, (i, j) in enumerate(E): W[i, j] = W[j, i] = w[k]
        c = np.array([sum(W[i, j] * (U[i] @ U[j]) for j in range(N)) for i in range(N)])
        Q = np.kron(np.diag(c) - W, np.eye(3))
        # restrict to tangent coords
        Tb = np.zeros((3 * N, 2 * N))
        for i in range(N): Tb[3 * i:3 * i + 3, 2 * i:2 * i + 2] = T[i]
        return Tb.T @ Q @ Tb
    full = np.mean(rays, axis=0)
    for lab, w in [('ray%d' % k, r) for k, r in enumerate(rays)] + [('mean', full)]:
        Qt = Q_of(w); ev_all = np.linalg.eigvalsh(Qt)
        Qf = Fn.T @ Qt @ Fn; ev = np.linalg.eigvalsh(Qf)
        print('  %-5s support %2d  Q on tangent: neg %d zero %d pos %d | Q on Flex/SO(3): %s' % (
            lab, (w > 1e-9).sum(), (ev_all < -1e-9).sum(), (abs(ev_all) <= 1e-9).sum(), (ev_all > 1e-9).sum(), np.round(ev, 4)))
    return E, rays, Fn

analyze('FCC (cuboctahedron)', cubocta())
analyze('HCP (anticuboctahedron)', anticubocta())
