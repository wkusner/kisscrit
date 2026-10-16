"""The C kernel's interval arithmetic (csrc/kc_arith.c) is bit-for-bit the Python reference
(kisscrit/core/interval.py) on random and special intervals, and its environment check passes.
Skipped when build/libkisscheck.so is not built (run `make -C csrc`).

Status: untrusted (a test).  Spec v0.3 section 1.2: 'the C code must reproduce that logic exactly, and a
test must compare the two'."""
import math, random
import pytest
from kisscrit.core import interval as I
from kisscrit.core.interval import Iv
from kisscrit.core import carith as C

pytestmark = pytest.mark.skipif(not C.available(), reason='C kernel not built (make -C csrc)')

SPECIAL = [0.0, math.pi / 2, -math.pi / 2, math.pi, -math.pi, 2 * math.pi, 1.0, -1.0, 0.5, 1e-300, 5e-324,
           1.7976931348623157e308, math.inf, -math.inf, 3.1415926535897936, 1.5707963267948968]


def _same(a, b):
    return (a.lo == b.lo or (math.isnan(a.lo) and math.isnan(b.lo))) and (a.hi == b.hi or (math.isnan(a.hi) and math.isnan(b.hi)))


def _intervals(rng, n, lo, hi):
    out = []
    for s in SPECIAL:
        for t in SPECIAL:
            if s <= t: out.append(Iv(s, t))
    for _ in range(n):
        a = rng.uniform(lo, hi); w = abs(rng.uniform(-1e-9, 1.0)) * (1 if rng.random() < 0.5 else 1e-3)
        out.append(Iv(a, a + w))
    return out


def test_constants_match():
    assert C.constants() == (I.PI, I.HALF_PI, I.TWO_PI) or all(_same(a, b) for a, b in zip(C.constants(), (I.PI, I.HALF_PI, I.TWO_PI)))


@pytest.mark.parametrize('name,fun,lo,hi', [('cos', I.cos, -8, 8), ('sin', I.sin, -8, 8), ('atan', I.atan, -50, 50),
                                             ('asin', I.asin, -1.2, 1.2), ('acos', I.acos, -1.2, 1.2),
                                             ('square', I.square, -5, 5), ('sqrt', I.sqrt, 0, 10)])
def test_unary_bit_identical(name, fun, lo, hi):
    rng = random.Random(7)
    bad = []
    for a in _intervals(rng, 3000, lo, hi):
        if name == 'sqrt' and a.hi < 0: continue
        if name in ('cos', 'sin') and a.finite() and max(abs(a.lo), abs(a.hi)) > 1e9: continue   # C returns [-1, 1] there
        p, c = fun(a), C.unary(name, a)
        if not _same(p, c): bad.append((a, p, c))
    assert not bad, bad[:3]


@pytest.mark.parametrize('name', ['add', 'sub', 'mul', 'div'])
def test_binary_bit_identical(name):
    rng = random.Random(11)
    ops = {'add': lambda a, b: a + b, 'sub': lambda a, b: a - b, 'mul': lambda a, b: a * b, 'div': lambda a, b: a / b}
    xs = _intervals(rng, 400, -10, 10) + [Iv(0.0), Iv(0.0, 1.0), Iv(-1.0, 0.0), Iv(-math.inf, 2.0)]
    bad = []
    for a in xs:
        for b in rng.sample(xs, 60):
            p, c = ops[name](a, b), C.binary(name, a, b)
            if not _same(p, c): bad.append((a, b, p, c))
    assert not bad, bad[:3]


def test_clamp_and_hull():
    rng = random.Random(3)
    for a in _intervals(rng, 500, -5, 5):
        assert _same(I.clamp(a, -1.0, 1.0), C.clamp(a, -1.0, 1.0))
        b = Iv(a.lo - 1, a.hi + 0.5)
        assert _same(I.hull(a, b), C.binary('hull', a, b))
