"""EXACT-1 milestone E2: exact data for DOD, M, B, r1, A, the hexagonal antiprism band, FCC and HCP.

Status: untrusted (exploration or tooling; no part of the completeness proof).

For each configuration:
  1. a high-precision approximation (scripts/exact/hiprec.py: Gauss-Newton from data/n12_above1.json, or a
     closed form), the minimal polynomial P of c = cos(theta) and every Gram entry in Q(c) guessed by LLL;
  2. exact checks in K = Q[t]/(P) (scripts/exact/verify.py): P irreducible, the root isolated by a Sturm
     count, unit diagonal, contacts = c, non-contacts < c, a positive definite 3x3 block with zero Schur
     complement (rank 3, PSD), the self-stress space, the balance of the stated stress and the signs of its
     weights, the rigidity rank, flex and kernel dimensions, the co-index and descending index (exact
     inertia of the second-order form on the kernel), the symmetry group (permutations preserving G exactly,
     with the determinant of the induced map), and the minimal polynomials of s^2, s and r;
  3. comparisons with the 2026-10-01 and 2026-10-08 notes and the stored double-precision data (labelled).

Writes data/exact/critical-configs-n12.json (exact values as strings; floats are labelled and decide
nothing).  "Exactly verified" in that file means: an exact computation in K accepted every listed identity
and inequality.  It is not the completeness checker.

Usage: python3 scripts/exact/build_e2.py [labels...]   (default: all eight)
"""
import json
import os
import subprocess
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mpmath as mp
import flint
import sympy as sp
from flint import fmpq_poly, fmpq
from kfield import Field, KE, qstr, nullspace
import hiprec as H
import verify as V

ROOT = H.ROOT
OUT = os.path.join(ROOT, 'data', 'exact', 'critical-configs-n12.json')
LABELS = ['DOD', 'M', 'B', 'r1', 'A', 'antiprism', 'FCC', 'HCP']
DPS = {'DOD': 120, 'M': 300, 'B': 700, 'r1': 700, 'A': 1600, 'antiprism': 120, 'FCC': 80, 'HCP': 80}

# Minimal polynomials as written in the notes (opened this session), converted to constant term first.
NOTES = {
    'DOD': {'c': [-1, 0, 5], 'r': [1, 4, 1, -6, 1], 'src': '02_Census/research-n12-critical-algebraic-2026-10-01.md'},
    'antiprism': {'c': [-1, 0, 3], 'r': [1, 4, 0, -8, 1], 'src': '02_Census/research-n12-critical-algebraic-2026-10-01.md'},
    'M': {'c': [-1, -8, -9, 40, 45, -32, 29], 'r': [1, 12, 57, 130, 121, -40, -158, -100, -45, -4, 38, -12, 1],
          'src': '02_Census/research-n12-above1-critical-data-2026-10-08.md (c); paper v6 Theorem thm:M (r)'},
    'B': {'c': [-1, -15, -72, -32, 478, 498, -792, -768, 683, 581, -152, -408, -16, 144],
          'src': '02_Census/research-n12-above1-critical-data-2026-10-08.md; 10-01 lists the same, highest degree first'},
    'r1': {'c': [1, -2, -4, 22, 73, -40, -316, 204, 1715, 2318, -2932, -10818, -16717, -9412, -140, 9344, 7888, 2240, 192],
           'src': '02_Census/research-n12-above1-critical-data-2026-10-08.md; 10-01 lists the same, highest degree first'},
    'A': {'c': [1, 8, 20, -56, -658, -2512, -3956, 6984, 58105, 160472, 199624, -166432, -1318992, -2985600, -3543296,
                -570368, 6978816, 16748544, 23597056, 23404544, 16855040, 8683520, 3047424, 655360, 65536],
          'src': '02_Census/research-n12-above1-critical-data-2026-10-08.md; 10-01 lists -1 times it, highest degree first'},
}
# Gram entries of r1 in the 2026-10-02 tex (spot check; same labels as data/n12_above1.json)
R1_TEX = {(1, 3): [-1, 0, 2], (1, 4): [0, -1, 2], (1, 5): [0, -1, 2]}


def git_state():
    """Commit hash of HEAD and whether scripts/exact has uncommitted changes."""
    h = subprocess.run(['git', '-C', ROOT, 'rev-parse', '--short', 'HEAD'], capture_output=True, text=True).stdout.strip()
    dirty = subprocess.run(['git', '-C', ROOT, 'status', '--porcelain', 'scripts/exact'],
                           capture_output=True, text=True).stdout.strip()
    return h, bool(dirty)


def same_up_to_sign(a, b):
    """True if integer coefficient lists a and b agree up to a global sign."""
    return [int(x) for x in a] == [int(x) for x in b] or [int(x) for x in a] == [-int(x) for x in b]


def approximate(label):
    """(contacts, mp coordinates, c, stored float stress or None, r_float) at DPS[label] digits."""
    dps = DPS[label]
    if label in H.ABOVE1:
        A, X0, sig, t = H.load_above1(label)
        stress = H.ABOVE1[label][1]
        X, c, w, res = H.newton(X0, A, dps, stress=stress, w0=sig if stress else None)
        assert res < mp.mpf(10) ** (-(dps - 30)), ('newton did not converge', label, res)
        return A, X, c, sig, t['r'], 'data/n12_above1.json types[%d], Newton to %d digits (residual < 1e-%d)' % (
            H.ABOVE1[label][0], dps, dps - 30)
    mp.mp.dps = dps + 20
    X = H.closed_form(label, dps)
    A, c = H.contacts_of(X, mp.mpf(10) ** (-(dps // 2)))
    s = mp.sqrt((1 - c) / 2)
    return A, X, c, None, float(s / (1 - s)), 'closed form at %d digits' % dps


def build(label):
    """Exact data and checks for one configuration (dict, JSON-ready)."""
    t0 = time.time()
    A, X, c, sig, r_float, how = approximate(label)
    dps = DPS[label]
    mp.mp.dps = dps + 20
    N = len(X)
    P = H.minpoly_guess(c, dps, 40)
    assert P is not None, 'no minimal polynomial guess'
    F = Field(P, approx=c)
    d = F.d
    # Gram entries: group equal values, guess each class in the basis 1, c, ..., c^(d-1)
    G, classes = {}, {}
    for i in range(N):
        for j in range(i, N):
            v = H.dot(X[i], X[j])
            key = mp.nstr(v, dps // 2)
            classes.setdefault(key, (v, []))[1].append((i, j))
    for key, (v, pairs) in classes.items():
        guess = H.in_basis_guess(v, c, d, dps)
        assert guess is not None, ('no Gram guess', label, pairs[0])
        e = F([guess[k] for k in range(d)])
        for p in pairs:
            G[p] = e
    rep = {'label': label, 'r_numerical': '%.13f (double precision, does not decide anything)' % r_float,
           'approximation': how, 'contacts': [list(e) for e in A], 'n_contacts': len(A)}
    # (1) identity
    rep['canonical_edges'] = V.canonical_form(A)
    sizes, gap = V.faces(X, A)
    from collections import Counter
    rep['faces'] = {'sizes': {str(k): v for k, v in sorted(Counter(sizes).items())},
                    'status': 'numerical (rotation system from %d-digit coordinates; smallest angular gap %s rad)'
                              % (dps, mp.nstr(gap, 4))}
    # (2) c and K
    rep['c'] = {'minpoly': [str(a) for a in P], 'degree': d, 'irreducible_over_Q': True,
                'interval': [qstr(F.lo), qstr(F.hi)], 'arb': F.eval_arb(F.gen(), 256).str(40, radius=True),
                'K': 'Q[t]/(P), t = c; elements are coefficient lists in 1, t, ..., t^(d-1), constant first'}
    # (4) Gram matrix
    checks, frame = V.gram_checks(F, A, G)
    rep['gram'] = {'entries': {'%d,%d' % p: F.coeffs(e) for p, e in sorted(G.items())}, 'checks': checks,
                   'distinct_values': len({tuple(F.coeffs(e)) for e in G.values()})}
    # symmetry (needed for the symmetric stress below)
    syms = V.symmetries(F, A, G, frame)
    rep['symmetry'] = V.group_name(syms, G)
    rep['symmetry']['status'] = 'exact (graph automorphisms preserving every Gram entry in K; det of the induced map in K)'
    # (5) stress space
    basis = V.stress_space(F, A, G, frame)
    st = {'dim': len(basis), 'basis': [[F.coeffs(x) for x in b] for b in basis]}
    one = [F.const(1)] * len(A)
    st['uniform_is_stress'] = V.is_stress(F, A, G, frame, one)
    # symmetric stresses: the group average (Reynolds operator) of the basis
    eidx = {e: n for n, e in enumerate(A)}
    sym_vecs = []
    for b in basis:
        acc = [F.const(0)] * len(A)
        for perm, _ in syms:
            for n, (i, j) in enumerate(A):
                acc[n] = acc[n] + b[eidx[tuple(sorted((perm[i], perm[j])))]]
        sym_vecs.append([x / len(syms) for x in acc])
    R, piv = __import__('kfield').rref(sym_vecs) if sym_vecs else ([], [])
    sym_basis = R[:len(piv)]
    st['symmetric_stress_dim'] = len(sym_basis)
    sym_pos = None
    if len(sym_basis) == 1:
        b = sym_basis[0]
        tot = sum(b[1:], b[0])
        if not tot.is_zero():
            w = [x / tot for x in b]
            sg = [x.sign() for x in w]
            st['symmetric_stress_normalized'] = [F.coeffs(x) for x in w]
            st['symmetric_stress_signs'] = sg
            if all(x > 0 for x in sg):
                sym_pos = w
    w_use = None
    if len(basis) == 1:
        b = basis[0]
        s = sum(b[1:], b[0])
        w = [x / s for x in b]
        signs = [x.sign() for x in w]
        st['normalized'] = [F.coeffs(x) for x in w]
        st['all_positive'] = all(x > 0 for x in signs)
        st['normalization'] = 'sum of weights = 1'
        w_use = w
    else:
        st['all_positive'] = 'see E4 (stress space of dimension > 1)'
    if sig is not None and len(basis) == 1:
        mp.mp.dps = 30
        st['stored_sigma_max_abs_difference'] = '%.2e (double precision vs exact; numerical)' % max(
            abs(float(F.eval_arb(x, 128).mid()) - s0) for x, s0 in zip(w_use, sig))
    # a nonnegative witness when the stress is not unique
    if len(basis) > 1:
        if sym_pos is not None:
            st['nonnegative_witness'] = 'the symmetric stress (group average), all weights > 0 (exact)'
            st['strictly_positive_stress_exists'] = True
            w_use = sym_pos
        elif st['uniform_is_stress']:
            st['nonnegative_witness'] = 'uniform weight (exact)'
            w_use = [F.const(fmpq(1, len(A)))] * len(A)
        elif label == 'M':
            # the stored stress has weight 0 on edge (5,9): the stress with w_(5,9) = 0, normalized
            e59 = A.index((5, 9))
            rows = V.balance_matrix(F, A, G, frame) + [[F.const(1) if e == e59 else F.const(0) for e in range(len(A))]]
            ray = nullspace(rows, len(A))
            assert len(ray) == 1
            r0 = ray[0]
            s = sum(r0[1:], r0[0])
            w = [x / s for x in r0]
            sg = [x.sign() for x in w]
            st['ray_with_zero_on_(5,9)'] = {'normalized': [F.coeffs(x) for x in w],
                                           'signs': sg, 'nonnegative': all(x >= 0 for x in sg)}
            mp.mp.dps = 30
            st['ray_vs_stored_sigma_max_abs_difference'] = '%.2e (numerical)' % max(
                abs(float(F.eval_arb(x, 128).mid()) - s0) for x, s0 in zip(w, sig))
            if all(x >= 0 for x in sg):
                st['nonnegative_witness'] = 'the ray with w_(5,9) = 0 (exact; one zero weight)'
                w_use = w
        elif label == 'HCP':
            eq = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (0, 5)]
            w = [F.const(1) if e in eq else F.const(0) for e in A]
            assert V.is_stress(F, A, G, frame, w)
            st['nonnegative_witness'] = 'weight 1 on the six equator edges, 0 elsewhere (exact; great-hexagon stress)'
            w_use = [x / 6 for x in w]
    rep['stress'] = st
    # rigidity, flexes, kernel, co-index
    unique = len(basis) == 1
    mo = V.motions(F, A, frame, w_use)
    assert mo['rigidity_rank'] == len(A) - len(basis), 'rank R != k - dim stress'
    rep['motions'] = mo
    if 'coindex' in mo and not unique and mo['kernel_dim_mod_rotations'] > 0:
        # q depends on the choice of stress when the stress is not unique and the kernel is not trivial
        for key in ('q_inertia_on_ker', 'nondegenerate', 'descending_index'):
            mo.pop(key, None)
        mo['coindex'] = 'not defined by a single stress (stress not unique and kernel > 0); see E4'
    elif 'coindex' in mo:
        mo['convention'] = ('co-index: negative directions of q on the kernel (r increases); descending index: '
                            'positive directions (r decreases); q has the 3 rotations as its only null directions '
                            'iff nondegenerate')
    else:
        mo['coindex'] = 'no nonnegative stress witness found in E2; see E4'
    # (3) s and r
    sr = V.s_and_r(F)
    rep['s2'], rep['s'], rep['r'] = sr['s2'], sr['s'], sr['r']
    # closed forms for degree <= 2
    if d <= 2:
        tt = sp.Symbol('t')
        roots = sp.solve(sp.Poly([int(a) for a in reversed(P)], tt).as_expr(), tt)
        cval = [z for z in roots if abs(float(z) - float(F.lo)) < 1e-6][0]
        sv = sp.sqrt((1 - cval) / 2)
        rv = sp.radsimp(sv / (1 - sv))
        x = sp.Symbol('x')
        mr = sp.Poly(sp.minimal_polynomial(rv, x), x)
        ours = [int(a) for a in sr['r']['minpoly']]
        ok = same_up_to_sign([int(a) for a in reversed(mr.all_coeffs())], ours)
        rep['closed_forms'] = {'c': str(cval), 's': str(sv), 'r': str(rv),
                               'r_minpoly_matches_sympy_minimal_polynomial': ok}
        assert ok
    # comparisons with the notes
    cmp = {}
    if label in NOTES:
        n = NOTES[label]
        cmp['c_minpoly_vs_note'] = 'equal up to sign' if same_up_to_sign(n['c'], P) else 'DIFFERENT'
        if 'r' in n:
            cmp['r_minpoly_vs_note'] = 'equal up to sign' if same_up_to_sign(n['r'], sr['r']['minpoly']) else 'DIFFERENT'
        cmp['source'] = n['src']
    if label == 'r1':
        cmp['gram_spot_check_vs_10-02_tex'] = all(G[p] == F(v) for p, v in R1_TEX.items())
    rep['comparisons'] = cmp
    rep['time_s'] = round(time.time() - t0, 1)
    return rep


def main():
    """Build the requested configurations and write the JSON file (merging with an existing one)."""
    labels = sys.argv[1:] or LABELS
    commit, dirty = git_state()
    data = {}
    if os.path.exists(OUT):
        data = json.load(open(OUT))
    data.setdefault('configs', {})
    for lab in labels:
        rep = build(lab)
        rep['status'] = {
            'identity': 'contact graph exact (contacts = c, non-contacts < c in K); symmetry group exact; '
                        'canonical form computed; faces numerical',
            'c': 'exactly verified (irreducible by factoring; root isolated by an exact Sturm count)',
            'r_and_s': 'exactly verified (factored over Z; factor selected by Arb; Sturm count 1)',
            'gram': 'exactly verified (identities and inequalities in K; rank 3 and PSD by exact Schur complement)',
            'stress': 'exactly verified (dimension, basis and signs in K)',
            'index': ('exactly verified (inertia in K)' if isinstance(rep['motions'].get('coindex'), int)
                      else 'not applicable (degenerate)'),
            'script': 'scripts/exact/build_e2.py at %s%s' % (commit, ' (DIRTY)' if dirty else ''),
        }
        data['configs'][lab] = rep
        print('%-9s d=%-3s k=%-2d stress_dim=%s  coindex=%s  sym=%s  %ss' % (
            lab, rep['c']['degree'], rep['n_contacts'], rep['stress']['dim'], rep['motions'].get('coindex'),
            rep['symmetry']['name'], rep['time_s']), flush=True)
    data['about'] = ('EXACT-1: exact descriptions of the N = 12 critical configurations found so far (not a '
                     'completeness claim). Exact values are strings: integers, rationals p/q, and elements of '
                     'K = Q[t]/(P) as coefficient lists in 1, t, ..., t^(d-1), constant term first. Floats are '
                     'labelled and decide nothing. "Exactly verified" = an exact computation in K accepted every '
                     'listed identity and inequality; this is not the completeness checker.')
    data['generated_by'] = {'script': 'scripts/exact/build_e2.py', 'commit': commit, 'dirty': dirty,
                            'versions': {'python-flint': flint.__version__, 'sympy': sp.__version__,
                                         'mpmath': mp.__version__, 'networkx': __import__('networkx').__version__,
                                         'python': sys.version.split()[0]}}
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    json.dump(data, open(OUT, 'w'), indent=1, sort_keys=False)


if __name__ == '__main__':
    main()
