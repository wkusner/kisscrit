/* kc_leaf.h: the leaf tests of the C checker kernel, and the per-chart driver that applies them to every leaf
 * and every contraction slab of a chart's proof tree.
 *
 * Status: trusted.  Part of the trusted base (transfer prompt v1.2, code standards; spec v0.3 sections 4.1-4.2).
 * Implements: checklist item F4 (the kernel's leaf tests): the exclusion tests of docs/lemmas.md L15 (vertex
 * sums equal 2 pi at a realization) and L16 (mean-value inclusion with given multipliers), the model's validity
 * flag (L3, L6-L8, L10'), and bisection; paper Appendix D (certification).  Each test is a bit-for-bit
 * reimplementation of the matching function of kisscrit/core/checker.py (the Python reference):
 *   kc_leaf_invalid_vsum   Checker.test_invalid_vsum on G.corners(lo, hi, tighten=0 or 1);
 *   kc_leaf_context        Checker.linear_context (the centre m, F(m), J(B));
 *   kc_leaf_hint           Checker.test_linear_mu and checker.combination;
 *   kc_leaf_exclude        Checker.exclude with its linear-programming tests left out (see below);
 *   kc_check_chart         the loop of checker.check_certificate over one chart, with the same left-outs.
 *
 * What the kernel leaves out, and why.  The reference also tries tests whose multipliers mu or vectors y are
 * chosen by floating-point linear programs (test_linear, test_linear_combo and test_linear_bound through scipy
 * HiGHS or the search's simplex; test_gordan and test_motzkin through scipy HiGHS).  The kernel does not run
 * them: certificate format v2 carries no Gordan or Motzkin witness (that is format v3, checklist F5), and a
 * linear program is outside the kernel's arithmetic.  A box that the kernel's own tests cannot exclude is
 * returned as KC_V_UNCLOSED, and the front end (kisscrit/ckernel.py, untrusted) hands it to the unchanged
 * Python reference.  So the kernel never rejects a box for good: it either proves the exclusion, or defers.
 *
 * Order of the tests (spec v0.3 section 4.2).  For format v2 the only witness is the multiplier hint, so step 1
 * ("the witness test, without tightening") is the hint test, run exactly where the reference runs it: on a
 * linear leaf with its own box as B, and on a slab with the contraction's box as B.  Then step 2 (invalid and
 * vertex sums, without tightening), step 4 (the same with monotonicity tightening, when the chart has a
 * pentagon), step 6 (bisection with the reference's depth and piece limits), step 7 (unclosed).  Steps 3 and 5
 * need a linear program in format v2 and are left to the fallback.
 *
 * Soundness does not depend on the order or on the limits: every closed box is closed by an interval proof
 * (an enclosure that excludes 0, or a validity flag that is 0 only when no valid point exists).  The order and
 * the limits only make the kernel close a subset of what the reference closes: on any box, the reference runs
 * the kernel's tests in the same order with other tests interleaved, and its bisection tree is the kernel's
 * tree with some nodes closed earlier, so it never needs more pieces than the kernel.
 */
#ifndef KC_LEAF_H
#define KC_LEAF_H
#include <stdint.h>
#include "kc_model.h"
#include "kc_replay.h"

/* checker.py MAX_DEPTH and MAX_PIECES: bisection depth and pieces allowed below one leaf or slab */
#define KC_MAX_DEPTH  40
#define KC_MAX_PIECES 400

/* Verdict of one unit (a leaf or a contraction slab).  The closed verdicts name the test that closed it, by
 * the reference's counter names (Checker.bump). */
#define KC_V_UNCLOSED     0   /* the kernel's tests did not exclude it: hand it to the Python fallback */
#define KC_V_INVALID      1   /* 'invalid':        no valid point in the box (validity flag 0, untightened) */
#define KC_V_VSUM         2   /* 'vsum':           a vertex sum minus 2 pi excludes 0 (untightened) */
#define KC_V_HINT         3   /* 'linear-hint':    the mean-value inclusion with the certificate's hint excludes 0 */
#define KC_V_INVALID_T    4   /* 'invalid-tight':  as 1, with monotonicity tightening */
#define KC_V_VSUM_T       5   /* 'vsum-tight':     as 2, with monotonicity tightening */
#define KC_V_BISECT       6   /* bisection, every piece closed by one of 1, 2, 4, 5 */
#define KC_V_TERMINAL     7   /* a terminal leaf no wider than the tolerance (recorded, not excluded) */
#define KC_V_TERM_WIDE    8   /* a terminal leaf wider than the tolerance: a failure, as in the reference */
#define KC_V_OPEN         9   /* an open leaf: a failure, as in the reference */

/* Counters per chart, indexed by these constants (the reference's Checker.stats keys) */
#define KC_C_INVALID      0   /* 'invalid' */
#define KC_C_VSUM         1   /* 'vsum' */
#define KC_C_HINT         2   /* 'linear-hint' */
#define KC_C_INVALID_T    3   /* 'invalid-tight' */
#define KC_C_VSUM_T       4   /* 'vsum-tight' */
#define KC_C_BISECT       5   /* 'bisect' */
#define KC_NCOUNT         6

/* Errors of the leaf tests (negative; the replay's KC_RERR_* codes are passed through unchanged) */
#define KC_LERR_MODEL   -20   /* kc_corners returned an error (the reference raises ArithmeticError or KeyError);
                               * the model's code is in the driver's *model_err */
#define KC_LERR_HINTROW -21   /* a hint names a vertex >= n (the reference raises IndexError in combination) */
#define KC_LERR_CAP     -22   /* more units than the caller's verdict array holds */
#define KC_LERR_SIZE    -23   /* the chart's box dimension differs from the graph's parameter count */
#define KC_LERR_NOMEM   -24   /* allocation failed */
#define KC_LERR_ARG     -25   /* depth or piece limit above the reference's, or negative */

/* Workspace for one graph: corner arrays, vertex sums and their gradients, the linear context, the bisection
 * boxes.  Allocated once per chart by kc_leaf_ws_new; not shared between threads. */
typedef struct kc_leaf_ws kc_leaf_ws;

kc_leaf_ws *kc_leaf_ws_new(const kc_graph *g);
void kc_leaf_ws_free(kc_leaf_ws *w);

/* Checker.test_invalid_vsum on the corners of [lo, hi] with or without tightening (G.corners(lo, hi,
 * tighten=...), clamp bound from lo).  *rule receives 0 (not excluded), 1 ('invalid') or 2 ('vsum').
 * Returns KC_OK or a kc_corners error code (positive, kc_model.h). */
int kc_leaf_invalid_vsum(kc_leaf_ws *w, const double *lo, const double *hi, int tighten, int *rule);

/* Checker.linear_context for the box B = [Blo, Bhi]: CB = G.corners(Blo, Bhi, tighten=False), the centre
 * m_j = 0.5 (Blo_j + Bhi_j), F(m) = the vertex sums of G.corners(m, m, Blo, Bhi), and J(B) = the vertex-sum
 * gradients of CB.  Stored in the workspace for kc_leaf_hint.  Returns KC_OK or a kc_corners error code. */
int kc_leaf_context(kc_leaf_ws *w, const double *Blo, const double *Bhi);

/* Checker.test_linear_mu on [lo, hi] (a sub-box of the context's B) with the hint rows[i] -> vals[i],
 * i < count (built as a Python dict, so a repeated row keeps its last value, then taken in increasing row
 * order as sorted(mu.items())).  Requires kc_leaf_context first.  *closed receives 1 iff the interval
 * s = mu.F(m) + sum_j (mu^T J(B))_j Y_j excludes 0; *s (may be NULL) receives s.  count = 0 gives *closed = 0
 * (the reference's `if not mu: return False`).  Returns KC_OK or KC_LERR_HINTROW. */
int kc_leaf_hint(kc_leaf_ws *w, const double *lo, const double *hi, const uint16_t *rows, const float *vals,
                 int64_t count, int *closed, kc_iv *s);

/* Checker.exclude(lo, hi) without its linear-programming tests: untightened invalid/vsum, tightened
 * invalid/vsum (when the graph has a pentagon), then bisection of the widest coordinate at its midpoint, the
 * lower half first, down to max_depth and max_pieces (each at most the reference's).  *verdict receives
 * KC_V_INVALID, KC_V_VSUM, KC_V_INVALID_T, KC_V_VSUM_T, KC_V_BISECT or KC_V_UNCLOSED; counts (KC_NCOUNT
 * entries, may be NULL) are incremented as the reference's Checker.bump would be, for a closed box only.
 * Returns KC_OK, a kc_corners error code, or KC_LERR_ARG. */
int kc_leaf_exclude(kc_leaf_ws *w, const double *lo, const double *hi, int max_depth, int max_pieces,
                    int *verdict, int64_t *counts);

/* The whole chart.  Replays the proof tree (kc_replay) from the decoded streams and the starting box, and for
 * every unit in the reference's order (for each node in preorder: a contraction's slabs in order, or the leaf)
 * writes its verdict to verdict[0 .. *nunits) and, if unit_kind is not NULL, its kind to unit_kind (0 for a
 * slab, the leaf event 1..8 for a leaf).  Per contraction node the context of its box is computed first
 * (as check_certificate does, whether or not a slab follows).  Leaves 1-6: the hint test for a linear leaf
 * with a nonempty hint (B = its own box), then kc_leaf_exclude.  Slabs: the hint test with the contraction's
 * box as B, then kc_leaf_exclude.  Terminal leaves: KC_V_TERM_WIDE iff max_j (hi_j - lo_j) > tol, the maximum
 * taken in coordinate order with Python's max semantics; else KC_V_TERMINAL.  Open leaves: KC_V_OPEN.
 * counts[KC_NCOUNT] receives the chart's counters (closed units only).  unit_ns, if not NULL, receives the
 * wall time of each unit in nanoseconds (CLOCK_MONOTONIC; measurement only, never part of a verdict).
 * Returns 0 when the replay completes, else the first error: a KC_RERR_* replay code (after the units before
 * it), KC_LERR_MODEL (with *model_err set), KC_LERR_HINTROW, KC_LERR_CAP, KC_LERR_SIZE, KC_LERR_NOMEM or
 * KC_LERR_ARG.  On an error, the units written so far keep their verdicts and nothing after them is decided. */
int kc_check_chart(const kc_graph *g, int n, const double *lo0, const double *hi0,
                   const uint8_t *ev, int64_t nev, const uint16_t *q, int64_t nq,
                   const float *mu, int64_t nmu, int has_mu, double tol, int max_depth, int max_pieces,
                   int8_t *verdict, int8_t *unit_kind, int64_t cap, int64_t *nunits, int64_t *counts, int64_t *unit_ns,
                   int *model_err);

#endif
