"""Exact second-order data at FCC and HCP (scripts/r1/exact_lambda.py).  Status: exact."""
import importlib.util, pathlib
import sympy as sp

_p = pathlib.Path(__file__).resolve().parents[1] / 'scripts' / 'r1' / 'exact_lambda.py'
_spec = importlib.util.spec_from_file_location('exact_lambda', _p)
exact_lambda = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(exact_lambda)


def test_fcc_lambda_is_one_sixth():
    r = exact_lambda.analyze(exact_lambda.fcc(), exact_lambda.fcc_stresses)
    assert (r['contacts'], r['rank'], r['flex_dim']) == (24, 20, 1)
    assert r['per_edge'] == {sp.Rational(1, 6)}
    for ev in r['eigen'].values():
        assert ev == {sp.Rational(1, 6): 1}


def test_hcp_lambda_at_least_one_twelfth():
    r = exact_lambda.analyze(exact_lambda.hcp(), exact_lambda.hcp_stresses)
    assert (r['contacts'], r['rank'], r['flex_dim']) == (24, 20, 1)
    assert r['per_edge'] == {sp.Rational(1, 12), sp.Rational(8, 33), sp.Rational(3, 22)}
    assert min(r['per_edge']) == sp.Rational(1, 12)
    for ev in r['eigen'].values():
        assert ev == {sp.Rational(3, 22): 1}
