"""Exact second-order values on the whole self-stress space at FCC and HCP (scripts/r1/exact_lambda_stress_space.py).

Status: untrusted (tests of an untrusted exact check; control item C5a).
Implements: paper Appendix B, Lemma B.2(b) and the Proposition B.5 sentence "for every self-stress w (of either
sign) f^T Q_w f = lambda sum w_ij, lambda = 1/6 at FCC and 3/22 at HCP".

Every assertion is an exact equality in K = Q(sqrt2, sqrt3) (sympy algebraic field elements, compared
coefficient by coefficient) or an exact sign; nothing here is floating point.  The mutant tests run the same
assertion helpers on deliberately wrong data and require them to fail (test power).
"""
import copy
import importlib.util
import pathlib

import pytest
import sympy as sp

_p = pathlib.Path(__file__).resolve().parents[1] / 'scripts' / 'r1' / 'exact_lambda_stress_space.py'
_spec = importlib.util.spec_from_file_location('exact_lambda_stress_space', _p)
X = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(X)

R = sp.Rational
LAM = {'FCC': X.el(R(1, 6)), 'HCP': X.el(R(3, 22))}
PER_EDGE = {'FCC': {R(1, 6): 24}, 'HCP': {R(3, 22): 6, R(1, 12): 12, R(8, 33): 6}}


@pytest.fixture(scope='module')
def fcc():
    """The full analysis at FCC with lambda = 1/6."""
    return X.analyze('FCC', X.fcc_points(), LAM['FCC'])


@pytest.fixture(scope='module')
def hcp():
    """The full analysis at HCP with lambda = 3/22."""
    return X.analyze('HCP', X.hcp_points(), LAM['HCP'])


# ---------------------------------------------------------------- assertion helpers (also used on mutants)

def assert_setup(U, S):
    """Part A: unit points, 24 contacts at exactly 60 degrees and none closer, rank 20, a 4-dimensional stress
    space whose basis vectors are stresses, null space of dimension 4 containing the rotations, and one unit flex
    F/|F| that is tangent, a flex of every edge and orthogonal to the rotations."""
    assert X.check_unit(U)
    assert len(S['E']) == 24 and S['packing'] and S['tangent_ok']
    assert S['R'].shape == (24, 24) and S['rank'] == 20
    assert len(S['stresses']) == 4 and all(X.is_stress(S['R'], s) for s in S['stresses'])
    assert S['null_dim'] == 4 and S['rot_ok']
    assert S['flex_dim'] == 1 and S['F_tangent'] and S['F_is_flex'] and S['F_perp_rot']
    assert S['normF2'] != X.ZERO


def assert_routes(q1, q2, expected):
    """Part B: the two routes agree edge by edge, every value is rational, and the multiset is `expected`."""
    assert len(q1) == len(q2) == 24
    assert all(a == b for a, b in zip(q1, q2))
    assert all(X.is_rational(a) for a in q1)
    counts = {}
    for a in q1:
        v = R(X.frac(a).numerator, X.frac(a).denominator)
        counts[v] = counts.get(v, 0) + 1
    assert counts == expected


def assert_stress_claim(S, q, lam):
    """Part C: sum_e s_e q_e = lam sum_e s_e for every basis stress, and a witness y with R y = q - lam 1 exists and
    passes exact multiplication (so the identity holds on the whole stress space)."""
    C = X.stress_space_report(S, q, lam)
    assert all(row['is_stress'] for row in C['basis'])
    assert all(row['holds'] for row in C['basis'])
    assert C['y'] is not None and C['witness_ok']
    assert X.verify_witness(S['R'], C['y'], [a - lam for a in q])
    return C


# ---------------------------------------------------------------- the exact claims

@pytest.mark.parametrize('name', ['FCC', 'HCP'])
def test_setup_routes_and_stress_space(name, fcc, hcp):
    """Parts A-C at both configurations."""
    r = fcc if name == 'FCC' else hcp
    assert_setup(r['S']['U'], r['S'])
    assert_routes(r['q1'], r['q2'], PER_EDGE[name])
    assert_stress_claim(r['S'], r['q1'], LAM[name])


@pytest.mark.parametrize('name', ['FCC', 'HCP'])
def test_extreme_rays_all_have_ratio_lambda(name, fcc, hcp):
    """Every extreme ray of the nonnegative stress cone is a stress and has ratio exactly lambda."""
    r = fcc if name == 'FCC' else hcp
    assert r['rays']
    for w in r['rays']:
        assert X.is_stress(r['S']['R'], w) and all(X.sign(x) >= 0 for x in w)
        assert X.ratio(w, r['q1']) == LAM[name]


def test_uniform_weight_is_a_stress_at_fcc_but_not_at_hcp(fcc, hcp):
    """Control ruling 2026-10-20: the uniform weight (1,...,1) is a self-stress at FCC and is not one at HCP."""
    assert X.is_stress(fcc['S']['R'], [X.ONE] * 24)
    assert fcc['C']['uniform_is_stress']
    assert not X.is_stress(hcp['S']['R'], [X.ONE] * 24)
    assert not hcp['C']['uniform_is_stress']


def test_hcp_vertex_0_is_not_in_equilibrium_under_the_uniform_weight():
    """By hand: the neighbours of c_0 = (1,0,0) at HCP are c_1, c_5 (equator), c_6 (top), c_9 (bottom); their sum is
    (2, sqrt3/3, 0), which is not parallel to c_0, so the uniform weight is not a stress at vertex 0."""
    U = X.hcp_points()
    E, _ = X.contact_edges(U)
    nb = sorted(j for (i, j) in E if i == 0) + sorted(i for (i, j) in E if j == 0)
    assert sorted(nb) == [1, 5, 6, 9]
    s = [sum((U[j][m] for j in nb), X.ZERO) for m in range(3)]
    assert s == [X.el(2), X.el(sp.sqrt(3) / 3), X.ZERO]
    assert X.cross(s, U[0]) != [X.ZERO] * 3


def test_hcp_strictly_positive_stress(hcp, capsys):
    """A self-stress with every entry > 0 exists at HCP: the sum of the extreme rays, normalized to sum 1.  It is
    printed exactly and pinned."""
    w = X.strictly_positive_stress(hcp)
    assert w is not None
    assert X.is_stress(hcp['S']['R'], w) and all(X.sign(x) == 1 for x in w) and sum(w, X.ZERO) == X.ONE
    with capsys.disabled():
        print('\nHCP strictly positive self-stress (edge: weight):',
              ', '.join('%s: %s' % (e, X.show(x)) for e, x in zip(hcp['S']['E'], w)))
    a, b, c = R(1, 20), R(11, 150), R(7, 200)
    expect = {(0, 1): a, (2, 3): a, (4, 5): a, (0, 5): b, (1, 2): b, (3, 4): b}
    for e, x in zip(hcp['S']['E'], w):
        assert x == X.el(expect.get(e, c)), e
    assert X.ratio(w, hcp['q1']) == LAM['HCP']


def test_points_agree_with_exact_lambda():
    """The points built here equal those of scripts/r1/exact_lambda.py, compared in K (sympy expressions are
    converted by to_number_field, which is exact)."""
    old = importlib.util.spec_from_file_location('exact_lambda', _p.parent / 'exact_lambda.py')
    m = importlib.util.module_from_spec(old); old.loader.exec_module(m)
    for mine, theirs in [(X.hcp_points(), m.hcp()), (X.fcc_points(), m.fcc())]:
        assert [[X.el(x) for x in p] for p in theirs] == mine


def test_exact_sign_and_coords():
    """`coords` round-trips and `sign` is right on hand cases, including near-cancellations."""
    cases = [(sp.sqrt(2) - R(141421, 100000), 1), (sp.sqrt(3) - sp.sqrt(2) - R(3178, 10000), 1),
             (sp.sqrt(6) - sp.sqrt(2) - sp.sqrt(3) + R(57, 100), -1), (5 - 2 * sp.sqrt(6), 1),
             (sp.sqrt(2) + sp.sqrt(3) - sp.sqrt(6) - R(7, 10000), 1), (sp.Integer(0), 0), (-sp.sqrt(3), -1)]
    for expr, sg in cases:
        x = X.el(expr)
        assert X._check_coords(x)
        assert X.sign(x) == sg, expr


# ---------------------------------------------------------------- test power: deliberate mutants

def test_mutant_perturbed_point_is_detected():
    """Mutant 1: the x coordinate of HCP point 6 moved by 1/7 (exact rational).  The set-up assertions fail."""
    U = X.hcp_points()
    U[6][0] = U[6][0] + X.el(R(1, 7))
    S = X.setup(U)
    with pytest.raises(AssertionError):
        assert_setup(U, S)


def test_mutant_wrong_per_edge_value_is_detected(hcp):
    """Mutant 2: one per-edge value off by 1/1000 (edge (0, 6), value 1/12).  Both the route comparison and the
    whole-stress-space claim fail."""
    q = list(hcp['q1'])
    k = hcp['S']['E'].index((0, 6))
    q[k] = q[k] + X.el(R(1, 1000))
    with pytest.raises(AssertionError):
        assert_routes(q, hcp['q2'], PER_EDGE['HCP'])
    with pytest.raises(AssertionError):
        assert_stress_claim(hcp['S'], q, LAM['HCP'])


def test_mutant_non_stress_in_basis_is_detected(hcp):
    """Mutant 3: the first basis stress replaced by the edge indicator e_0, which is not a stress."""
    S = dict(hcp['S'])
    S['stresses'] = copy.copy(S['stresses'])
    S['stresses'][0] = [X.ONE] + [X.ZERO] * 23
    with pytest.raises(AssertionError):
        assert_setup(S['U'], S)
    with pytest.raises(AssertionError):
        assert_stress_claim(S, hcp['q1'], LAM['HCP'])


def test_mutant_wrong_lambda_is_detected(hcp):
    """Mutant 4: lambda = 1/12 (the old lower bound) instead of 3/22."""
    with pytest.raises(AssertionError):
        assert_stress_claim(hcp['S'], hcp['q1'], X.el(R(1, 12)))


def test_mutant_wrong_chart_is_detected(hcp):
    """Mutant 5: route 2 with the linear chart u_i = c_i + t F_i instead of the exponential chart, whose second
    derivative is -2 <F_i, F_j>.  The two routes disagree."""
    t = X._gFF
    lin = sp.Poly(-2 * t, X._Ai, X._Aj, X._gcc, X._gcF, X._gFc, X._gFF, domain=sp.QQ)
    with pytest.raises(AssertionError):
        assert_routes(hcp['q1'], X.route2(hcp['S'], poly=lin), PER_EDGE['HCP'])
