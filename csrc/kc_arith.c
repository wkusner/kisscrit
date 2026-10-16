/* kc_arith.c: certified interval arithmetic for the C checker kernel (see kc_arith.h for the contract).
 *
 * Status: trusted.  Implements docs/code-tour.md section 0.1; a bit-for-bit reimplementation of
 * kisscrit/core/interval.py, whose docstrings carry the justification of every rule restated here.
 * Build with -O2 -fno-fast-math -ffp-contract=off (csrc/Makefile); kc_init checks the environment.
 */
#include <math.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846   /* the float seed only; never decides a verdict */
#endif
#include <fenv.h>
#include <float.h>
#include <flint/arb.h>
#include <flint/arf.h>
#include "kc_arith.h"

kc_iv KC_PI, KC_HALF_PI, KC_TWO_PI;
static arb_t ARB_PI;                      /* Arb's pi at KC_PREC bits, as interval.py's _PI */
static int initialized = 0;

/* --------------------------------------------------------------------- rounding of balls to doubles */
double kc_down(const arb_t x)
{   /* the lower bound of the ball at KC_PREC bits, then the largest double not above it: this is what
     * interval.py's down() computes by stepping float(x.lower()) downward until arb(f) <= lower. */
    arf_t b; arf_init(b);
    arb_get_lbound_arf(b, x, KC_PREC);
    double d = arf_get_d(b, ARF_RND_FLOOR);
    arf_clear(b);
    return d;
}

double kc_up(const arb_t x)
{
    arf_t b; arf_init(b);
    arb_get_ubound_arf(b, x, KC_PREC);
    double d = arf_get_d(b, ARF_RND_CEIL);
    arf_clear(b);
    return d;
}

/* --------------------------------------------------------------------- environment */
static int fma_contraction_present(void)
{   /* 0.1 * 10 rounds to exactly 1.0 in binary64, so a * b + c with c = -1 is 0 without contraction; a fused
     * multiply-add would keep the exact product 1.0000000000000000555 and return 5.55e-17.  volatile stops
     * constant folding. */
    volatile double a = 0.1, b = 10.0, c = -1.0;
    double r = a * b + c;
    return r != 0.0;
}

int kc_init(void)
{
    if (fegetround() != FE_TONEAREST) return 1;      /* the ulp-stepping argument assumes round-to-nearest */
    if (fma_contraction_present()) return 2;          /* -ffp-contract=off was not honoured */
    if (sizeof(double) != 8 || DBL_MANT_DIG != 53) return 3;
    arb_init(ARB_PI); arb_const_pi(ARB_PI, KC_PREC);
    arb_t t; arb_init(t);
    KC_PI.lo = kc_down(ARB_PI); KC_PI.hi = kc_up(ARB_PI);
    arb_mul_2exp_si(t, ARB_PI, -1); KC_HALF_PI.lo = kc_down(t); KC_HALF_PI.hi = kc_up(t);
    arb_mul_2exp_si(t, ARB_PI, 1);  KC_TWO_PI.lo = kc_down(t);  KC_TWO_PI.hi = kc_up(t);
    arb_clear(t);
    initialized = 1;
    return 0;
}

/* --------------------------------------------------------------------- basic operations */
/* The four IEEE operations are correctly rounded, so the exact result is within one ulp of the computed
 * double; stepping one ulp outward (nextafter) is a rigorous bound.  NaN (undefined, e.g. inf - inf) widens
 * to the whole line; a finite overflow is bounded by the largest double.  These are _op_down/_op_up. */
static double op_down(double r, double a, double b)
{
    if (isnan(r)) return -INFINITY;
    if (r == INFINITY) return (isinf(a) || isinf(b)) ? INFINITY : DBL_MAX;
    if (r == -INFINITY) return -INFINITY;
    return nextafter(r, -INFINITY);
}

static double op_up(double r, double a, double b)
{
    if (isnan(r)) return INFINITY;
    if (r == -INFINITY) return (isinf(a) || isinf(b)) ? -INFINITY : -DBL_MAX;
    if (r == INFINITY) return INFINITY;
    return nextafter(r, INFINITY);
}

kc_iv kc_point(double x) { kc_iv r = {x, x}; return r; }
int kc_is_zero(kc_iv a) { return a.lo == 0.0 && a.hi == 0.0; }
int kc_finite(kc_iv a) { return isfinite(a.lo) && isfinite(a.hi); }

kc_iv kc_add(kc_iv a, kc_iv b)
{   /* [a.lo + b.lo rounded down, a.hi + b.hi rounded up] */
    kc_iv r = { op_down(a.lo + b.lo, a.lo, b.lo), op_up(a.hi + b.hi, a.hi, b.hi) };
    return r;
}

kc_iv kc_sub(kc_iv a, kc_iv b)
{   /* [a.lo - b.hi rounded down, a.hi - b.lo rounded up] */
    kc_iv r = { op_down(a.lo - b.hi, a.lo, b.hi), op_up(a.hi - b.lo, a.hi, b.lo) };
    return r;
}

kc_iv kc_neg(kc_iv a) { kc_iv r = { -a.hi, -a.lo }; return r; }     /* exact */

kc_iv kc_mul(kc_iv a, kc_iv b)
{   /* the outward-rounded extremes of the four endpoint products; an exactly-zero factor gives exactly
     * [0, 0] (structural zeros, see the header); 0 * inf is undefined and widens to the whole line. */
    kc_iv whole = { -INFINITY, INFINITY };
    if (kc_is_zero(a) || kc_is_zero(b)) return kc_point(0.0);
    double ends_a[2] = { a.lo, a.hi }, ends_b[2] = { b.lo, b.hi };
    double lo = INFINITY, hi = -INFINITY;
    for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++) {
            double x = ends_a[i], y = ends_b[j];
            if ((x == 0.0 && isinf(y)) || (y == 0.0 && isinf(x))) return whole;
            double p = x * y;
            double d = op_down(p, x, y), u = op_up(p, x, y);
            if (d < lo) lo = d;
            if (u > hi) hi = u;
        }
    kc_iv r = { lo, hi };
    return r;
}

kc_iv kc_div(kc_iv a, kc_iv b)
{   /* (-inf, inf) if b may be 0 (unless a is exactly [0, 0]); else a * [1/b.hi down, 1/b.lo up] */
    kc_iv whole = { -INFINITY, INFINITY };
    if (kc_is_zero(a)) return kc_point(0.0);
    if (b.lo <= 0.0 && 0.0 <= b.hi) return whole;
    kc_iv inv = { op_down(1.0 / b.hi, 1.0, b.hi), op_up(1.0 / b.lo, 1.0, b.lo) };
    return kc_mul(a, inv);
}

kc_iv kc_square(kc_iv a)
{   /* range of t^2: monotone on each side of 0; across 0 it is [0, max(|lo|, |hi|)^2 rounded up] */
    if (a.lo >= 0) return kc_mul(a, a);
    if (a.hi <= 0) return kc_mul(kc_neg(a), kc_neg(a));
    double m = (-a.lo > a.hi) ? -a.lo : a.hi;
    kc_iv r = { 0.0, kc_mul(kc_point(m), kc_point(m)).hi };
    return r;
}

kc_iv kc_sqrt(kc_iv a)
{   /* range of sqrt on [max(lo, 0), hi]; IEEE sqrt is correctly rounded, so one ulp outward is rigorous.
     * The caller checks a.hi >= 0.  (Not used by the kernel at present; kept to match the reference.) */
    double lo = a.lo > 0.0 ? a.lo : 0.0;
    double l = (lo == 0.0) ? 0.0 : nextafter(sqrt(lo), -INFINITY);
    double h = isinf(a.hi) ? INFINITY : nextafter(sqrt(a.hi > 0.0 ? a.hi : 0.0), INFINITY);
    kc_iv r = { l > 0.0 ? l : 0.0, h };
    return r;
}

kc_iv kc_clamp(kc_iv a, double lo, double hi)
{   /* image of [a.lo, a.hi] under t -> min(max(t, lo), hi) */
    double l = a.lo < lo ? lo : a.lo; if (l > hi) l = hi;
    double h = a.hi < lo ? lo : a.hi; if (h > hi) h = hi;
    kc_iv r = { l, h };
    return r;
}

kc_iv kc_hull(kc_iv a, kc_iv b)
{
    kc_iv r = { a.lo < b.lo ? a.lo : b.lo, a.hi > b.hi ? a.hi : b.hi };
    return r;
}

/* --------------------------------------------------------------------- elementary functions */
/* mono: rigorous range of a monotone function from its values at the endpoints (interval.py _mono).
 * An infinite endpoint maps to the limit value lo_lim / hi_lim. */
typedef void (*arb_fun)(arb_t, const arb_t, slong);

static double ev(arb_fun f, double t, int round_down, int increasing, double lo_lim, double hi_lim)
{
    if (isinf(t)) return ((t > 0) == increasing) ? hi_lim : lo_lim;
    arb_t x, v; arb_init(x); arb_init(v);
    arb_set_d(x, t); f(v, x, KC_PREC);
    double d = round_down ? kc_down(v) : kc_up(v);
    arb_clear(x); arb_clear(v);
    return d;
}

static kc_iv mono(arb_fun f, kc_iv a, int increasing, double lo_lim, double hi_lim)
{
    kc_iv r;
    if (increasing) { r.lo = ev(f, a.lo, 1, 1, lo_lim, hi_lim); r.hi = ev(f, a.hi, 0, 1, lo_lim, hi_lim); }
    else            { r.lo = ev(f, a.hi, 1, 0, lo_lim, hi_lim); r.hi = ev(f, a.lo, 0, 0, lo_lim, hi_lim); }
    if (r.lo < lo_lim) r.lo = lo_lim;
    if (r.hi > hi_lim) r.hi = hi_lim;
    return r;
}

kc_iv kc_atan(kc_iv a)
{   /* arctan is increasing; its range is within [-pi/2, pi/2], so the limits are the outward-rounded ends */
    return mono(arb_atan, a, 1, -KC_HALF_PI.hi, KC_HALF_PI.hi);
}

kc_iv kc_asin(kc_iv a)
{   /* asin on the argument clamped to [-1, 1] (the caller raises its own validity flag), increasing */
    return mono(arb_asin, kc_clamp(a, -1.0, 1.0), 1, -INFINITY, INFINITY);
}

kc_iv kc_acos(kc_iv a)
{   /* acos on the argument clamped to [-1, 1], decreasing */
    return mono(arb_acos, kc_clamp(a, -1.0, 1.0), 0, -INFINITY, INFINITY);
}

kc_iv kc_cos(kc_iv a)
{   /* the values at the endpoints (Arb, rounded outward), widened to 1 where an even multiple of pi may lie
     * inside and to -1 where an odd multiple may.  Candidate multiples k are enumerated with the float
     * M_PI from one below floor(a.lo / pi) to one unit past a.hi, margins that exceed any rounding error of
     * the float seed; whether k pi lies inside is decided only by the Arb comparisons, conservatively
     * (undecidable counts as inside).  A non-finite or very wide interval (>= 6 > 2 pi) gives [-1, 1]. */
    kc_iv whole = { -1.0, 1.0 };
    if (!kc_finite(a) || !(a.hi - a.lo < 6.0)) return whole;
    /* arguments beyond 1e9 in magnitude never occur in the face model (all angles lie in [-2 pi, 4 pi]);
     * returning [-1, 1] there is sound and keeps the multiple-of-pi enumeration in a safe integer range
     * (this is the one place where the C enclosure is wider than the Python reference's) */
    if (fabs(a.lo) > 1e9 || fabs(a.hi) > 1e9) return whole;
    arb_t x, ca, cb, p, lo_b, hi_b; arb_init(x); arb_init(ca); arb_init(cb); arb_init(p); arb_init(lo_b); arb_init(hi_b);
    arb_set_d(x, a.lo); arb_cos(ca, x, KC_PREC);
    arb_set_d(x, a.hi); arb_cos(cb, x, KC_PREC);
    double da = kc_down(ca), db = kc_down(cb), ua = kc_up(ca), ub = kc_up(cb);
    double lo = da < db ? da : db, hi = ua > ub ? ua : ub;
    arb_set_d(lo_b, a.lo); arb_set_d(hi_b, a.hi);
    long k = (long)floor(a.lo / M_PI) - 1;
    while (k * M_PI <= a.hi + 1.0) {
        arb_mul_si(p, ARB_PI, k, KC_PREC);
        int possibly_inside = !arb_lt(p, lo_b) && !arb_gt(p, hi_b);
        if (possibly_inside) { if (k % 2 == 0) hi = 1.0; else lo = -1.0; }
        k++;
    }
    arb_clear(x); arb_clear(ca); arb_clear(cb); arb_clear(p); arb_clear(lo_b); arb_clear(hi_b);
    kc_iv r = { lo < -1.0 ? -1.0 : lo, hi > 1.0 ? 1.0 : hi };
    return r;
}

kc_iv kc_sin(kc_iv a)
{   /* sin t = cos(t - pi/2), on the outward-rounded shifted interval */
    kc_iv whole = { -1.0, 1.0 };
    if (!kc_finite(a)) return whole;
    arb_t x, h, s; arb_init(x); arb_init(h); arb_init(s);
    arb_mul_2exp_si(h, ARB_PI, -1);
    kc_iv sh;
    arb_set_d(x, a.lo); arb_sub(s, x, h, KC_PREC); sh.lo = kc_down(s);
    arb_set_d(x, a.hi); arb_sub(s, x, h, KC_PREC); sh.hi = kc_up(s);
    arb_clear(x); arb_clear(h); arb_clear(s);
    return kc_cos(sh);
}
