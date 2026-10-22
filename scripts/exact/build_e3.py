"""EXACT-1 milestone E3: the remaining critical values above r = 1 (data/n12_above1.json).

Status: untrusted (exploration or tooling; no part of the completeness proof).

Same pipeline as scripts/exact/build_e2.py (high-precision Newton and LLL guesses, then exact checks in
K = Q[t]/(P)), applied to the seven configurations that had no exactly verified minimal polynomial:
three with a high-precision relation in the 2026-10-08 note (degrees 28, 46, 59) and four with no relation
of degree <= 90 at 1600 digits (docs/open-questions.md section 3).

  python3 scripts/exact/build_e3.py LABEL [DPS]     build one configuration and merge it into the JSON file
  python3 scripts/exact/build_e3.py --probe LABEL DPS DMAX
        refine to DPS digits and run one LLL on 1, c, ..., c^DMAX (timing and a relation if one is found);
        prints, writes nothing.

Labels are the ten-digit radii used as keys of results/certified_n12.json.
"""
import json
import os
import sys
import time

sys.set_int_max_str_digits(0)  # high-precision values have more than 4300 digits

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mpmath as mp
import flint
import sympy as sp
import hiprec as H
import build_e2 as B2

# label -> index in data/n12_above1.json types
E3 = {'1.0010073939': 2, '1.0008452029': 3, '1.0003062547': 6, '1.0003062544': 7,
      '1.0002374095': 8, '1.0002115891': 10, '1.0002111636': 11}
# high-precision relations of the 2026-10-08 note (constant term first), for comparison only
NOTE_0810 = {
    '1.0003062544': [1, -2, -33, 52, 527, -514, -5151, 1744, 32368, 7872, -126784, -89536, 286976, 257536, -471552,
                     -233472, 1261568, 648192, -2473984, -2580480, 1507328, 3342336, 294912, -3932160, -6553600,
                     -6029312, -2621440, 0, 262144],
    '1.0010073939': [1, -9, -62, 124, 2754, 3222, -9936, 97868, -843739, -5610029, -2895174, 13966308, -33019248,
                     478757592, 2222442672, -129070496, -4476594622, 3201426862, -66584234764, -244958627512,
                     203394325980, 1164327808532, 147066116624, 3120035394840, 12514334373618, -9438611164226,
                     -39375740290508, 26316542962872, 27394480710744, -199476657206536, -175331810845808,
                     81038097563776, -388363839647299, -803583807548069, -321051022957814, -668406066707268,
                     -1159466815150414, -713832071096186, -748937309628672, -619650104933444, -283005582104343,
                     -164373430871185, -76890853580462, -14896160966940, -202190308920, 207954079488, 14668221312],
    '1.0008452029': [4, -71, 191, 1942, -5600, -54034, -75028, 2620744, 4200456, -67506350, -184492458, 1610295412,
                     7539666240, -49311474880, -220706637804, 1367973801824, 5447180442652, -27504004096601,
                     -129002907354119, 382040091865026, 2500233224004832, -3436832945979706, -34972333786030016,
                     15643126902919096, 344147259774908720, 37677772759474876, -2399133223143880620,
                     -1149324984015654856, 11981984331711930496, 8993705379722419792, -43355680907638020792,
                     -41883171037626117312, 115918847654357606140, 132538093268484554583, -238451714637420297823,
                     -305021747841565128582, 401103682510760074976, 538861542188064624634, -581232990867322769396,
                     -758152760700760039688, 733196215454463592456, 865044365616295679506, -789162123915131619402,
                     -811877187451468296300, 707069368008736857664, 630144349365798998064, -515809728982865991388,
                     -405174842960679478048, 294937873883444076900, 219571116472747359401, -128446915369015911609,
                     -97415787155885905938, 37362044217859712544, 35345904886787687058, -6097116176439216568,
                     -8797824588068002680, -178670960810440128, 1104928774081296000, 266409929136496640,
                     20449956625776640],
}
SRC = '02_Census/research-n12-above1-critical-data-2026-10-08.md (high-precision relation, not verified there)'

def minpoly_guess_steps(c, dps, dmax=40, steps=(20, 40, 60, 80, 100, 130, 160, 200)):
    """LLL guess of the minimal polynomial of c: one LLL per trial degree bound in `steps` (up to dmax),
    then the irreducible factor (over Z) that is smallest at c.  A guess only; checked exactly later."""
    from flint import fmpz_poly
    for D in [s for s in steps if s <= max(dmax, steps[0])]:
        with mp.workdps(dps):
            a = H.relation([c ** k for k in range(D + 1)], dps - 10)
        if a is None:
            continue
        _, fac = fmpz_poly(a).factor()
        with mp.workdps(dps):
            best = min(fac, key=lambda f: abs(mp.polyval([int(x) for x in reversed(f[0].coeffs())], c)))
        f = [int(x) for x in best[0].coeffs()]
        return [-x for x in f] if f[-1] < 0 else f
    return None


_LAST_P = {}


def _minpoly_and_remember(c, dps, dmax=40):
    """minpoly_guess_steps with degree bound 200; remembers P for the ratio guess below."""
    P = minpoly_guess_steps(c, dps, 200)
    _LAST_P['P'] = P
    return P


def omega_basis(P):
    """The integral elements 1, w_1, ..., w_{d-1} of Q(c), w_k = a_d c^k + a_{d-1} c^(k-1) + ... + a_{d-k+1} c
    for P = sum a_i x^i (constant first), as integer coefficient lists in powers of c."""
    d = len(P) - 1
    return [[1]] + [[0] + [P[d - k + 1 + j] for j in range(k)] for k in range(1, d)]


def ratio_guess(v, c, d, dps):
    """Guess v = -a(c)/b(c) with a, b of small height: one LLL on w_k(c) and v w_k(c) (dimension 2d).
    Gram entries of configurations whose P is far from monic are quotients of small elements of Q(c) whose
    power-basis form has very large coefficients; the power-basis guess cannot see them at practical
    precision.  Returns the power-basis form a b^(-1) mod P (exact division in Q[t]/(P)), or None."""
    P = _LAST_P['P']
    om = omega_basis(P)
    with mp.workdps(dps):
        ov = [mp.polyval(list(reversed(o)), c) for o in om]
        a = H.relation(ov + [x * v for x in ov], dps - 10)
    if a is None:
        return None
    A = fmpq_poly([0])
    B = fmpq_poly([0])
    for k, o in enumerate(om):
        A += a[k] * fmpq_poly(o)
        B += a[d + k] * fmpq_poly(o)
    if B.is_zero():
        return None
    Pq = fmpq_poly(P)
    g, s, _ = (B % Pq).xgcd(Pq)
    assert g.degree() == 0
    return (-A * s / g[0]) % Pq


def in_basis_or_ratio(v, c, d, dps):
    """Power-basis guess (hiprec.in_basis_guess) first; if it finds nothing, the ratio guess."""
    g = _in_basis_power(v, c, d, dps)
    return g if g is not None else ratio_guess(v, c, d, dps)


from flint import fmpq_poly  # noqa: E402
_in_basis_power = H.in_basis_guess
H.minpoly_guess = _minpoly_and_remember
H.in_basis_guess = in_basis_or_ratio
for lab, idx in E3.items():
    H.ABOVE1[lab] = (idx, True)
    if lab in NOTE_0810:
        B2.NOTES[lab] = {'c': NOTE_0810[lab], 'src': SRC}


def probe(label, dps, dmax):
    """Refine to dps digits and run a single LLL on the powers of c up to dmax (timing; prints only)."""
    A, X0, sig, t = H.load_above1(label)
    t0 = time.time()
    X, c, w, res = H.newton(X0, A, dps, stress=True, w0=sig)
    t1 = time.time()
    with mp.workdps(dps):
        a = H.relation([c ** k for k in range(dmax + 1)], dps - 10)
    t2 = time.time()
    print('%s dps=%d newton %.1fs residual %s  LLL(dim %d) %.1fs  relation: %s' % (
        label, dps, t1 - t0, mp.nstr(res, 3), dmax + 1, t2 - t1,
        'none' if a is None else 'degree %d, height %d digits' % (max(i for i, x in enumerate(a) if x), len(str(max(abs(x) for x in a))))), flush=True)
    return a


def main():
    """Build one configuration (or probe) and merge it into data/exact/critical-configs-n12.json."""
    if sys.argv[1] == '--probe':
        probe(sys.argv[2], int(sys.argv[3]), int(sys.argv[4]))
        return
    label = sys.argv[1]
    if len(sys.argv) > 2:
        B2.DPS[label] = int(sys.argv[2])
    commit, dirty = B2.git_state()
    rep = B2.build(label)
    rep['status'] = {
        'identity': 'contact graph exact (contacts = c, non-contacts < c in K); symmetry group exact; '
                    'canonical form computed; faces numerical',
        'c': 'exactly verified (irreducible by factoring; root isolated by an exact Sturm count)',
        'r_and_s': 'exactly verified (factored over Z; factor selected by Arb; Sturm count 1)',
        'gram': 'exactly verified (identities and inequalities in K; rank 3 and PSD by exact Schur complement)',
        'stress': 'exactly verified (dimension, basis and signs in K)',
        'index': ('exactly verified (inertia in K)' if isinstance(rep['motions'].get('coindex'), int)
                  else 'not applicable (degenerate)'),
        'script': 'scripts/exact/build_e3.py at %s%s (pipeline of build_e2.py)' % (commit, ' (DIRTY)' if dirty else ''),
    }
    data = json.load(open(B2.OUT))
    data['configs'][label] = rep
    data.setdefault('generated_by_e3', {})[label] = {'script': 'scripts/exact/build_e3.py', 'commit': commit,
                                                     'dirty': dirty, 'dps': B2.DPS[label]}
    json.dump(data, open(B2.OUT, 'w'), indent=1)
    print('%-13s d=%-3s k=%-2d stress_dim=%s coindex=%s sym=%s %ss cmp=%s' % (
        label, rep['c']['degree'], rep['n_contacts'], rep['stress']['dim'], rep['motions'].get('coindex'),
        rep['symmetry']['name'], rep['time_s'], rep['comparisons'].get('c_minpoly_vs_note')), flush=True)


if __name__ == '__main__':
    main()
