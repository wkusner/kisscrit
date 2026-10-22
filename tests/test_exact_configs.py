"""Exact re-verification of data/exact/critical-configs-n12.json (EXACT-1, milestone E5).

Status: untrusted (a test of stored exact data; no part of the completeness proof).

Every claim is re-checked from the stored strings with exact arithmetic in K = Q[t]/(P) (python-flint; signs
on a rational isolating interval checked by an exact Sturm count, decided by Arb balls that exclude zero):
the field (P irreducible, one root in the interval), the Gram matrix (unit diagonal, contacts = c, non-contacts
< c, a positive definite frame block with zero Schur complement), the canonical contact graph, the stress
space (each stored basis vector balances, dimension = nullity, rank = k - dim), the stored unique stress
(balanced, all weights > 0), the cone (each stored extreme ray is a nonnegative stress of sum 1 whose zero set
has rank dim - 1, and the set is closed under the symmetry group), the stored strictly positive stress, the
uniform-weight flag, the rigidity rank, flex and kernel dimensions, the co-index and descending index, the
symmetry group, and the minimal polynomials and intervals of s^2, s, r.  Mutants of the stored data must fail.
No floating point decides anything.  Runtime is dominated by the degree-46 and degree-59 configurations.
"""
import copy
import json
import os
import sys

import pytest

sys.set_int_max_str_digits(0)
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
sys.path.insert(0, os.path.join(ROOT, 'scripts', 'exact'))
from flint import fmpz_poly, fmpq_poly, arb, ctx  # noqa: E402
from kfield import Field, rank, sturm_count, _q  # noqa: E402
import verify as V  # noqa: E402
import weights_e4 as W  # noqa: E402

DATA = json.load(open(os.path.join(ROOT, 'data', 'exact', 'critical-configs-n12.json')))
LABELS = list(DATA['configs'])
_CACHE = {}


def setup(label):
    """Field, contacts, Gram matrix, frame and stress basis of a stored configuration (cached)."""
    if label not in _CACHE:
        F, A, G = W.load(DATA['configs'][label])
        checks, frame = V.gram_checks(F, A, G)
        basis = V.stress_space(F, A, G, frame)
        _CACHE[label] = (F, A, G, frame, checks, basis)
    return _CACHE[label]


def vec(F, coeff_lists):
    """A stored vector of K elements."""
    return [F(c) for c in coeff_lists]


def test_inventory_and_wording():
    """The eleven configurations of E2 and E3 are present; nothing calls itself certified or complete."""
    assert set(LABELS) == {'DOD', 'M', 'B', 'r1', 'A', 'antiprism', 'FCC', 'HCP', '1.0003062544', '1.0010073939',
                           '1.0008452029'}
    for f in ('data/exact/critical-configs-n12.json', 'data/exact/critical-configs-n12.md', 'docs/exact-configs.md'):
        text = open(os.path.join(ROOT, f)).read().lower()
        assert 'certified' not in text, f
    assert 'found so far' in DATA['about'] and 'not a completeness claim' in DATA['about']


@pytest.mark.parametrize('label', LABELS)
def test_field_and_gram(label):
    """P irreducible with exactly one root in the stored interval; every Gram identity and inequality."""
    d = DATA['configs'][label]
    F, A, G, frame, checks, basis = setup(label)
    assert F.d == d['c']['degree']
    content, fac = fmpz_poly([int(a) for a in d['c']['minpoly']]).factor()
    assert len(fac) == 1 and fac[0][1] == 1
    assert sturm_count(fmpq_poly([int(a) for a in d['c']['minpoly']]), _q(d['c']['interval'][0]),
                       _q(d['c']['interval'][1])) == 1
    for key in ('diagonal_is_one', 'contacts_equal_c', 'noncontacts_below_c', 'schur_complement_zero'):
        assert checks[key] and d['gram']['checks'][key]
    assert len(G) == 78  # all pairs i <= j
    assert V.canonical_form(A) == d['canonical_edges']


@pytest.mark.parametrize('label', LABELS)
def test_stress_space_and_cone(label):
    """Stress basis, unique stress, extreme rays, positive stress and uniform flag, from the stored data."""
    d = DATA['configs'][label]
    w = d['weights']
    F, A, G, frame, checks, basis = setup(label)
    k = len(A)
    assert len(basis) == w['dim'] == d['stress']['dim']
    stored = [vec(F, b) for b in w['basis']]
    assert all(V.is_stress(F, A, G, frame, b) for b in stored)
    assert rank(stored) == w['dim'] and w['balance_rank'] == k - w['dim']
    assert V.is_stress(F, A, G, frame, [F.const(1)] * k) == w['uniform_is_stress']
    if w['dim'] == 1:
        u = vec(F, d['stress']['normalized'])
        assert V.is_stress(F, A, G, frame, u) and all(x.sign() > 0 for x in u)
        assert sum(u[1:], u[0]) == 1
    m = w['dim']
    a = [[basis[i][e] for i in range(m)] for e in range(k)]
    rays = [vec(F, r['weights']) for r in w['extreme_rays']]
    assert len(rays) == w['n_extreme_rays'] >= 1
    keys = set()
    for r in rays:
        assert V.is_stress(F, A, G, frame, r) and all(x.sign() >= 0 for x in r) and sum(r[1:], r[0]) == 1
        assert rank([a[e] for e in range(k) if r[e].is_zero()]) == m - 1  # extreme
        keys.add(tuple(tuple(F.coeffs(x)) for x in r))
    assert len(keys) == len(rays)
    idx = {tuple(sorted(e)): n for n, e in enumerate(A)}
    for p, _ in V.symmetries(F, A, G, frame):
        for r in rays:
            img = [None] * k
            for n, (i, j) in enumerate(A):
                img[idx[tuple(sorted((p[i], p[j])))]] = r[n]
            assert tuple(tuple(F.coeffs(x)) for x in img) in keys
    assert w['strictly_positive_stress_exists']
    pos = vec(F, w['strictly_positive_stress']['weights'])
    assert V.is_stress(F, A, G, frame, pos) and all(x.sign() > 0 for x in pos)


@pytest.mark.parametrize('label', LABELS)
def test_motions_index_and_symmetry(label):
    """Rigidity rank, flex and kernel dimensions, co-index and descending index, symmetry group."""
    d = DATA['configs'][label]
    F, A, G, frame, checks, basis = setup(label)
    pos = vec(F, d['weights']['strictly_positive_stress']['weights'])
    mo = V.motions(F, A, frame, pos)
    st = d['motions']
    for key in ('rigidity_rank', 'flex_dim_mod_rotations', 'kernel_dim_mod_rotations'):
        assert mo[key] == st[key], key
    if isinstance(st.get('coindex'), int):
        assert mo['coindex'] == st['coindex'] and mo['descending_index'] == st['descending_index']
        assert mo['nondegenerate']
    else:  # FCC, HCP: stress not unique and kernel > 0; the form depends on the stress
        assert d['weights']['dim'] > 1 and st['kernel_dim_mod_rotations'] > 0
    g = V.group_name(V.symmetries(F, A, G, frame), G)
    for key in ('order', 'proper_order', 'name', 'contains_inversion'):
        assert g[key] == d['symmetry'][key], key
    assert sorted(len(o) for o in W.edge_orbits(A, V.symmetries(F, A, G, frame))) == \
        sorted(d['weights']['edge_orbit_sizes'])


@pytest.mark.parametrize('label', LABELS)
def test_s_and_r(label):
    """Minimal polynomials of s^2, s, r: irreducible, one root in the stored interval, and that root is the
    value computed from c (Arb enclosure from the c interval inside the stored interval)."""
    d = DATA['configs'][label]
    F = setup(label)[0]
    old = ctx.prec
    try:
        ctx.prec = 512
        while F.hi - F.lo > _q('1/%d' % 2 ** 300):
            F.refine()
        cb = F.ball()
        s2 = (1 - cb) / 2
        vals = {'s2': s2, 's': s2.sqrt(), 'r': s2.sqrt() / (1 - s2.sqrt())}
        for name, x in vals.items():
            e = d[name]
            f = fmpz_poly([int(a) for a in e['minpoly']])
            _, fac = f.factor()
            assert len(fac) == 1 and fac[0][1] == 1 and f.degree() == e['degree']
            lo, hi = _q(e['interval'][0]), _q(e['interval'][1])
            assert sturm_count(fmpq_poly([int(a) for a in e['minpoly']]), lo, hi) == 1
            assert x > arb(lo) and x < arb(hi), name
    finally:
        ctx.prec = old


def test_mutants_fail():
    """Altered data must be rejected: a Gram entry, P, an interval, a contact, a stress weight."""
    d0 = DATA['configs']['r1']

    def fails(d, stage):
        try:
            F, A, G = W.load(d)
            checks, frame = V.gram_checks(F, A, G)
            if stage == 'stress':
                u = vec(F, d['stress']['normalized'])
                assert V.is_stress(F, A, G, frame, u) and all(x.sign() > 0 for x in u)
        except AssertionError:
            return True
        return False
    m1 = copy.deepcopy(d0)  # one Gram entry changed by 1/1000 in its constant term
    co = m1['gram']['entries']['0,2']
    co[0] = str(_q(co[0]) + _q('1/1000'))
    assert fails(m1, 'gram')
    m2 = copy.deepcopy(d0)  # P changed
    m2['c']['minpoly'][0] = str(int(m2['c']['minpoly'][0]) + 1)
    assert fails(m2, 'gram')
    m3 = copy.deepcopy(d0)  # interval moved away from the root
    m3['c']['interval'] = ['1/3', '2/5']
    assert fails(m3, 'gram')
    m4 = copy.deepcopy(d0)  # a contact dropped: that pair has G = c, not < c
    m4['contacts'] = m4['contacts'][1:]
    assert fails(m4, 'gram')
    m5 = copy.deepcopy(d0)  # a stress weight with the wrong sign
    m5['stress']['normalized'][3] = ['-' + x if not x.startswith('-') else x[1:] for x in m5['stress']['normalized'][3]]
    assert fails(m5, 'stress')
    assert not fails(copy.deepcopy(d0), 'stress')  # the unaltered data pass


def test_high_precision_candidates():
    """Stored high-precision candidates are labelled as such; the polynomial itself is irreducible (an exact
    statement about the polynomial only; that c is its root is high precision, not proved)."""
    for label, cd in DATA.get('high_precision_candidates', {}).items():
        assert cd['status'].startswith('high precision only')
        f = fmpz_poly([int(a) for a in cd['minpoly']])
        _, fac = f.factor()
        assert f.degree() == cd['degree'] and len(fac) == 1 and fac[0][1] == 1
        assert label not in DATA['configs']
