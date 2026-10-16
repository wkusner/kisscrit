#define _POSIX_C_SOURCE 200809L
/* test_arith.c: containment tests of kc_arith against Arb at 256 bits (spec v0.3 section 4.3).
 *
 * Status: untrusted (a test).  For each function, N random intervals over the ranges the face model uses,
 * plus the special endpoints 0, +-pi/2, +-pi, +-1; for every interval, the exact function value at the two
 * endpoints and at 8 interior points, as a 256-bit Arb ball, must lie inside the kc enclosure.  The basic
 * operations are tested on random endpoint pairs the same way.  Exits nonzero on the first failure.
 * Usage: test_arith [N]   (default 1000000)
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846   /* the float seed only; never decides a verdict */
#endif
#include <time.h>
#include <flint/arb.h>
#include "kc_arith.h"

#define HI_PREC 256
static unsigned long long rng_state = 0x9E3779B97F4A7C15ULL;
static double urand(void)                      /* uniform in [0, 1), xorshift64* */
{
    rng_state ^= rng_state >> 12; rng_state ^= rng_state << 25; rng_state ^= rng_state >> 27;
    return (double)((rng_state * 0x2545F4914F6CDD1DULL) >> 11) / 9007199254740992.0;
}
static double uni(double a, double b) { return a + (b - a) * urand(); }

static const double SPECIAL[] = { 0.0, M_PI / 2, -M_PI / 2, M_PI, -M_PI, 1.0, -1.0, 2 * M_PI, 0.5, -0.5 };
#define NSPECIAL (sizeof(SPECIAL) / sizeof(SPECIAL[0]))

/* does the 256-bit ball v lie inside [lo, hi]?  (the exact value is in v, so this implies containment) */
static int inside(const arb_t v, double lo, double hi)
{
    arb_t L, H; arb_init(L); arb_init(H); arb_set_d(L, lo); arb_set_d(H, hi);
    int ok = (isinf(lo) && lo < 0 ? 1 : arb_ge(v, L)) && (isinf(hi) && hi > 0 ? 1 : arb_le(v, H));
    arb_clear(L); arb_clear(H);
    return ok;
}

typedef kc_iv (*kc_fun)(kc_iv);
typedef void (*arb_fun)(arb_t, const arb_t, slong);

static long test_unary(const char *name, kc_fun f, arb_fun g, double lo_range, double hi_range, long n, double *ns)
{
    arb_t x, v; arb_init(x); arb_init(v);
    long fails = 0; double t_sum = 0; struct timespec t0, t1;
    for (long i = 0; i < n; i++) {
        double a, b;
        if (i < (long)(NSPECIAL * NSPECIAL)) { a = SPECIAL[i % NSPECIAL]; b = SPECIAL[i / NSPECIAL]; }
        else { a = uni(lo_range, hi_range); b = a + fabs(uni(-1e-9, 1.0)) * (urand() < 0.5 ? 1 : 0.01); }
        if (b < a) { double t = a; a = b; b = t; }
        if (g == arb_asin || g == arb_acos) {           /* these are defined on [-1, 1] only: the kc functions clamp */
            a = a < -1 ? -1 : (a > 1 ? 1 : a);
            b = b < -1 ? -1 : (b > 1 ? 1 : b);
        }
        kc_iv A = { a, b };
        clock_gettime(CLOCK_MONOTONIC, &t0);
        kc_iv R = f(A);
        clock_gettime(CLOCK_MONOTONIC, &t1);
        t_sum += (t1.tv_sec - t0.tv_sec) * 1e9 + (t1.tv_nsec - t0.tv_nsec);
        for (int j = 0; j < 10; j++) {
            double t = (j == 0) ? a : (j == 1) ? b : uni(a, b);
            arb_set_d(x, t); g(v, x, HI_PREC);
            if (!inside(v, R.lo, R.hi)) {
                if (fails < 5) printf("  FAIL %s on [%.17g, %.17g] at %.17g: enclosure [%.17g, %.17g]\n", name, a, b, t, R.lo, R.hi);
                fails++;
            }
        }
    }
    arb_clear(x); arb_clear(v);
    *ns = t_sum / n;
    return fails;
}

static long test_binary(const char *name, int op, long n)
{
    arb_t x, y, v; arb_init(x); arb_init(y); arb_init(v);
    long fails = 0;
    for (long i = 0; i < n; i++) {
        double a = uni(-10, 10), b = a + fabs(uni(0, 2)), c = uni(-10, 10), d = c + fabs(uni(0, 2));
        if (urand() < 0.1) { a = uni(-1e-6, 1e-6); b = a + 1e-7; }
        kc_iv A = { a, b }, B = { c, d }, R;
        if (op == 0) R = kc_add(A, B); else if (op == 1) R = kc_sub(A, B); else if (op == 2) R = kc_mul(A, B); else R = kc_div(A, B);
        for (int j = 0; j < 8; j++) {
            double s = uni(a, b), t = uni(c, d);
            arb_set_d(x, s); arb_set_d(y, t);
            if (op == 0) arb_add(v, x, y, HI_PREC); else if (op == 1) arb_sub(v, x, y, HI_PREC);
            else if (op == 2) arb_mul(v, x, y, HI_PREC); else { if (c <= 0 && 0 <= d) continue; arb_div(v, x, y, HI_PREC); }
            if (!inside(v, R.lo, R.hi)) {
                if (fails < 5) printf("  FAIL %s on [%g,%g] op [%g,%g] at (%g,%g): [%g,%g]\n", name, a, b, c, d, s, t, R.lo, R.hi);
                fails++;
            }
        }
    }
    arb_clear(x); arb_clear(y); arb_clear(v);
    return fails;
}

int main(int argc, char **argv)
{
    long n = argc > 1 ? atol(argv[1]) : 1000000;
    int rc = kc_init();
    if (rc) { printf("kc_init failed: code %d (1 rounding mode, 2 FMA contraction, 3 double format)\n", rc); return 1; }
    printf("environment ok: round-to-nearest, no FMA contraction; KC_PI = [%.17g, %.17g]\n", KC_PI.lo, KC_PI.hi);
    long total = 0; double ns;
    struct { const char *name; kc_fun f; arb_fun g; double lo, hi; } U[] = {
        { "cos",  kc_cos,  arb_cos,  -8, 8 }, { "sin",  kc_sin,  arb_sin,  -8, 8 },
        { "atan", kc_atan, arb_atan, -50, 50 }, { "asin", kc_asin, arb_asin, -1, 1 }, { "acos", kc_acos, arb_acos, -1, 1 },
    };
    for (unsigned i = 0; i < sizeof(U) / sizeof(U[0]); i++) {
        long f = test_unary(U[i].name, U[i].f, U[i].g, U[i].lo, U[i].hi, n, &ns);
        printf("%-5s %ld intervals x 10 points: %ld failures, %.0f ns per call\n", U[i].name, n, f, ns);
        total += f;
    }
    const char *B[] = { "add", "sub", "mul", "div" };
    for (int op = 0; op < 4; op++) {
        long f = test_binary(B[op], op, n / 4);
        printf("%-5s %ld pairs x 8 points: %ld failures\n", B[op], n / 4, f);
        total += f;
    }
    printf(total ? "FAILED: %ld containment failures\n" : "all containment tests passed\n", total);
    return total ? 1 : 0;
}
