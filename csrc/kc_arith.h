/* kc_arith.h: outward-rounded interval arithmetic for the C checker kernel (containment tested against Arb by csrc/test_arith.c).
 *
 * Status: trusted.  Part of the trusted base (listed in docs/code-tour.md section 0; code standards in docs/decisions.md, 2026-10-16; spec v0.3 section 1.2).
 * Implements: the arithmetic contract of docs/code-tour.md section 0.1, as a bit-for-bit reimplementation of
 * kisscrit/core/interval.py (the Python reference): an interval is a pair of doubles [lo, hi]; the IEEE
 * operations + - * / sqrt are correctly rounded and stepped one ulp outward; sin, cos, atan, asin, acos and pi
 * are evaluated as Arb balls at KC_PREC bits and rounded outward to doubles.  Every function returns an
 * interval that contains the exact real result for every point of its argument intervals.
 *
 * Infinite endpoints always mean "unbounded range of finite reals".  Two conventions are shared with the
 * reference (and with the search, which is untrusted): a product with an exactly-zero factor [0, 0] is
 * exactly [0, 0], and division by an interval that may contain 0 gives (-inf, inf) unless the numerator is
 * exactly [0, 0].  Exact zeros arise only from structural zeros of gradients, never from rounded operations.
 */
#ifndef KC_ARITH_H
#define KC_ARITH_H
#include <flint/arb.h>

#define KC_PREC 64            /* working precision of the Arb balls, as ctx.prec = 64 in interval.py */

typedef struct { double lo, hi; } kc_iv;

/* constants: outward-rounded enclosures of pi, pi/2, 2 pi (interval.py PI, HALF_PI, TWO_PI) */
extern kc_iv KC_PI, KC_HALF_PI, KC_TWO_PI;

/* kc_init must be called once before any other function: it computes the constants and checks the
 * floating-point environment (round-to-nearest, no FMA contraction); returns 0 on success, else a code. */
int kc_init(void);

/* rounding of an Arb ball to doubles: the largest double <= every point of the ball, and the smallest >= */
double kc_down(const arb_t x);
double kc_up(const arb_t x);

kc_iv kc_point(double x);                 /* the exact interval [x, x] */
int   kc_is_zero(kc_iv a);                /* a is exactly [0, 0] */
int   kc_finite(kc_iv a);
kc_iv kc_add(kc_iv a, kc_iv b);
kc_iv kc_sub(kc_iv a, kc_iv b);
kc_iv kc_neg(kc_iv a);
kc_iv kc_mul(kc_iv a, kc_iv b);
kc_iv kc_div(kc_iv a, kc_iv b);
kc_iv kc_square(kc_iv a);
kc_iv kc_sqrt(kc_iv a);
kc_iv kc_clamp(kc_iv a, double lo, double hi);
kc_iv kc_hull(kc_iv a, kc_iv b);
kc_iv kc_cos(kc_iv a);
kc_iv kc_sin(kc_iv a);
kc_iv kc_atan(kc_iv a);
kc_iv kc_asin(kc_iv a);
kc_iv kc_acos(kc_iv a);

#endif
