"""ctypes interface to the C kernel's interval arithmetic (csrc/kc_arith.c), for the differential tests and,
later, for the kernel's front end.

Status: untrusted (a thin loader; the arithmetic it calls is the trusted C code, and the verdicts of the C
kernel are checked against the Python reference in tests/test_carith.py).
Implements: nothing mathematical; it maps kc_iv structs to kisscrit.core.interval.Iv.

The library is build/libkisscheck.so, built by `make -C csrc`.  `available()` says whether it is present and
its kc_init() accepted the floating-point environment.
"""
import ctypes, os, pathlib
from .interval import Iv

ROOT = pathlib.Path(__file__).resolve().parents[2]
LIB = ROOT / 'build' / 'libkisscheck.so'
_lib = None


class KcIv(ctypes.Structure):
    """The C struct kc_iv {double lo, hi}."""
    _fields_ = [('lo', ctypes.c_double), ('hi', ctypes.c_double)]


def load():
    """Load the library once and run kc_init; raises OSError if it is missing or the environment check fails."""
    global _lib
    if _lib is not None: return _lib
    lib = ctypes.CDLL(str(LIB))
    lib.kc_init.restype = ctypes.c_int
    rc = lib.kc_init()
    if rc != 0: raise OSError('kc_init failed with code %d (1 rounding mode, 2 FMA contraction, 3 double format)' % rc)
    for name in ('add', 'sub', 'mul', 'div', 'hull'):
        f = getattr(lib, 'kc_' + name); f.restype = KcIv; f.argtypes = [KcIv, KcIv]
    for name in ('neg', 'square', 'sqrt', 'cos', 'sin', 'atan', 'asin', 'acos'):
        f = getattr(lib, 'kc_' + name); f.restype = KcIv; f.argtypes = [KcIv]
    lib.kc_clamp.restype = KcIv; lib.kc_clamp.argtypes = [KcIv, ctypes.c_double, ctypes.c_double]
    _lib = lib
    return lib


def available():
    """True iff the library can be loaded and initialised."""
    try:
        load(); return True
    except OSError:
        return False


def _to_c(x): return KcIv(x.lo, x.hi)
def _to_py(x): return Iv(x.lo, x.hi)


def unary(name, a):
    """Apply kc_<name> to the interval a and return an Iv."""
    return _to_py(getattr(load(), 'kc_' + name)(_to_c(a)))


def binary(name, a, b):
    """Apply kc_<name> to the intervals a, b and return an Iv."""
    return _to_py(getattr(load(), 'kc_' + name)(_to_c(a), _to_c(b)))


def clamp(a, lo, hi):
    return _to_py(load().kc_clamp(_to_c(a), lo, hi))


def constants():
    """The C kernel's PI, HALF_PI, TWO_PI."""
    lib = load()
    return tuple(_to_py(KcIv.in_dll(lib, n)) for n in ('KC_PI', 'KC_HALF_PI', 'KC_TWO_PI'))
