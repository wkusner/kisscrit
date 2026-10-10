import json, os, math
import numpy as np
import networkx as nx
import kisscrit as kc
from kisscrit.search import graphs_on, solve_graph, chi_M0N
from kisscrit.hull import enumerate_triangulations, tri_faces, theta_star

HERE = os.path.dirname(__file__)


def ico():
    phi = (1 + 5 ** 0.5) / 2
    V = []
    for s1 in (-1, 1):
        for s2 in (-1, 1):
            V += [(0, s1, s2 * phi), (s1, s2 * phi, 0), (s2 * phi, 0, s1)]
    return kc.normalize(np.array(V, float))


def test_radius_conversions():
    assert abs(kc.r_of_c(0.5) - 1.0) < 1e-15                   # 60 degrees: r = 1
    assert abs(kc.c_of_r(kc.r_of_c(0.3)) - 0.3) < 1e-14


def test_icosahedron_is_strict_max():
    a = kc.analyze(ico())
    assert a['critical'] and a['k'] == 30 and a['coindex'] == 0 and a['kernel'] == 0
    assert abs(kc.r_of_c(a['c']) - 1.1085085) < 1e-6
    assert kc.group_orders(ico()) == (120, 60)


def test_n4_census():
    rng = np.random.default_rng(0); C = kc.Census(4, 3)
    for G in graphs_on(4):
        for X in solve_graph(G, 3, rng, 10):
            C.add(X)
    T = C.sorted()
    assert len(T) == 2
    assert abs(T[0]['r'] - (2 + 6 ** 0.5)) < 1e-9 and T[0]['coindex'] == 0 and T[0]['labeled'] == 2
    assert abs(T[1]['r'] - (1 + 2 ** 0.5)) < 1e-9 and T[1]['coindex'] == 1 and T[1]['labeled'] == 3
    assert 2 - 3 == chi_M0N(4)


def test_n12_data_coindices():
    D = json.load(open(os.path.join(HERE, '..', 'data', 'n12_above1.json')))
    for t in D['types']:
        a = kc.analyze(np.array(t['X']))
        assert a['critical'] and a['coindex'] == t['coindex']


def test_kkt_roundtrip_trigonal_prism():
    # N=6 trigonal prism saddle: r = (3 + sqrt 21)/4
    rng = np.random.default_rng(3); C = kc.Census(6, 3)
    G = nx.circular_ladder_graph(3)   # prism graph, 9 edges
    for X in solve_graph(G, 3, rng, 30):
        C.add(X)
    rs = [t['r'] for t in C.sorted()]
    assert any(abs(r - (3 + 21 ** 0.5) / 4) < 1e-9 for r in rs)


def test_triangulation_counts():
    assert [len(enumerate_triangulations(n)) for n in range(4, 10)] == [1, 1, 2, 5, 14, 50]


def test_theta_star_icosahedron():
    F = [sorted(t) for t in tri_faces(nx.icosahedral_graph())]
    th, X = theta_star(F, np.random.default_rng(0), starts=2)
    assert abs(th - 63.43494882) < 1e-5


def test_minpoly_tetrahedron():
    import mpmath as mp
    from kisscrit.exact import minpoly
    mp.mp.dps = 60
    assert minpoly(mp.mpf(-1) / 3, maxdeg=4) in ([1, 3], [-1, -3])
