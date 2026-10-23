/* kc_replay.h: the replay of one chart's proof tree, for the C checker kernel.
 *
 * Status: trusted.  Part of the trusted base (listed in docs/code-tour.md section 0; code standards in docs/decisions.md, 2026-10-16; spec v0.3 sections 4.1-4.2).
 * Implements: checklist item F3b, docs/certificate-format.md section 3 (the proof tree: contraction rounds with
 * grid values and multiplier hints, the removed slabs, the leaf events, the midpoint split and the order of the
 * children), as a bit-for-bit reimplementation of kisscrit/core/certificate.py:replay (the Python reference).
 *
 * Contract.  Given one chart's decoded streams (events uint8, q uint16, mu float32, as certificate.unpack returns
 * them) and its starting box [lo0, hi0] of dimension n, kc_replay_next returns, one call per node, exactly the
 * sequence of nodes that certificate.replay yields, in the same preorder: the same node kinds, the same boxes as
 * IEEE bit patterns, the same contraction slabs (slab box, the box it was cut from, coordinate k, side) in the
 * same order, and the same multiplier hints (row indices and float32 values).  Where the reference raises an
 * exception, kc_replay_next returns the matching error code at the same point of the sequence (after the same
 * nodes).  tests/test_creplay.py compares the two on the `cond` certificates and on malformed streams.
 *
 * The replay is a pure function of its inputs: the state is the caller's kc_replay struct and the caller's
 * workspace (size from kc_replay_work_doubles); there are no globals and no allocation.  It uses only IEEE
 * double + - * and comparisons, as the reference does; no elementary function, no FLINT.
 *
 * Validity of a node: the pointers in a kc_node (boxes, slabs) point into the workspace and stay valid only
 * until the next call of kc_replay_next; hint rows and values point into the caller's q and mu arrays.
 */
#ifndef KC_REPLAY_H
#define KC_REPLAY_H
#include <stdint.h>

#define KC_REPLAY_MAXN 64          /* box dimension limit (measured on `cond`: n <= 14) */

/* event codes of docs/certificate-format.md section 3 (certificate.py EV_ROUND, EV_SPLIT, EV_CONTRACT, LEAF) */
#define KC_EV_ROUND     15
#define KC_EV_SPLIT     16         /* 16 + k, k < 48 */
#define KC_EV_CONTRACT  64         /* 64 + k */
#define KC_LEAF_LINEAR   3         /* the leaf that carries a hint; leaves are the events 1..8 */
#define KC_QMAX     65536

/* node kinds: a leaf is its event code 1..8 (1 invalid, 2 vsum, 3 linear, 4 gordan-inherited, 5 gordan,
 * 6 motzkin, 7 terminal, 8 open); the other two kinds take the codes of the events that start them */
#define KC_NODE_CONTRACT 15
#define KC_NODE_SPLIT    16

/* return values of kc_replay_next */
#define KC_REPLAY_NODE        1    /* a node was written to *out */
#define KC_REPLAY_DONE        0    /* the tree is complete and every stream fully consumed */
/* errors, each where the reference raises: ValueError with the message shown, or IndexError */
#define KC_RERR_OUTSIDE      -1    /* 'contraction outside the box' */
#define KC_RERR_EMPTY_ROUND  -2    /* 'empty contraction round' */
#define KC_RERR_END_INSIDE   -3    /* 'event stream ended inside the tree' */
#define KC_RERR_UNKNOWN      -4    /* 'unknown event %d' */
#define KC_RERR_TRAILING     -5    /* 'trailing events' */
#define KC_RERR_INDEX        -6    /* IndexError: a q entry read past the end of q, or a coordinate k >= n */
/* errors of this implementation only (the reference has no such limit); never hit with the documented sizes */
#define KC_RERR_LIMIT        -7    /* n > KC_REPLAY_MAXN, n < 0, or stack_cap < 1 (from kc_replay_init) */
#define KC_RERR_STACK        -8    /* the node stack is full (cannot happen with stack_cap >= kc_replay_stack_need) */

/* A multiplier hint (certificate.py read_mu).  present = 0 is the reference's None; present = 1 is a dict,
 * possibly empty: the map rows[j] -> vals[j] for j < count, built in order, so for a repeated row the last entry
 * gives the value (Python dict semantics; tests/test_creplay.py builds the dict the same way).  count is
 * min(c, rows available in q, values available in mu) for the stored count c, as the reference's slicing and
 * zip give.  rows and vals are NULL when count = 0. */
typedef struct {
    int present;
    int64_t count;
    const uint16_t *rows;
    const float *vals;
} kc_hint;

/* A slab removed by a contraction round: the box [lo, hi], the box [cut_lo, cut_hi] it was cut from,
 * coordinate k, side -1 (below the new lower bound) or +1 (above the new upper bound), and its hint. */
typedef struct {
    const double *lo, *hi, *cut_lo, *cut_hi;
    int k, side;
    kc_hint hint;
} kc_slab;

/* One node.  kind: a leaf code 1..8, KC_NODE_CONTRACT or KC_NODE_SPLIT.
 *   leaf:     [lo, hi] is the leaf box; hint is the linear leaf's hint (present only for kind 3 with a mu stream).
 *   split:    [lo, hi] is the box, k the coordinate.
 *   contract: [lo, hi] is the box at the start of the round, [new_lo, new_hi] the box after it, and
 *             slabs[0 .. nslabs) the removed region in the reference's order. */
typedef struct {
    int kind, k;
    const double *lo, *hi;
    const double *new_lo, *new_hi;
    int nslabs;
    const kc_slab *slabs;
    kc_hint hint;
} kc_node;

/* Replay state.  Treat as opaque: set by kc_replay_init, advanced by kc_replay_next. */
typedef struct {
    int n, has_mu, phase, err;
    const uint8_t *ev; int64_t nev, i;
    const uint16_t *q; int64_t nq, qpos;
    const float *mu; int64_t nmu, mpos;
    double *stack; int64_t stack_cap, stack_len;      /* stack entry j: lo at stack[2nj], hi at stack[2nj + n] */
    double *lo, *hi, *base_lo, *base_hi, *cur_lo, *cur_hi, *slab_store;
    kc_hint hint_lower[KC_REPLAY_MAXN], hint_upper[KC_REPLAY_MAXN];   /* this round's hints per coordinate */
    kc_slab slabs[2 * KC_REPLAY_MAXN];
} kc_replay;

/* A stack capacity that always suffices: 1 + the number of split events in ev (each split adds one entry). */
int64_t kc_replay_stack_need(const uint8_t *ev, int64_t nev);

/* The number of doubles of workspace kc_replay_init needs: 2 n stack_cap (stack) + 6 n (box, round box, cut box)
 * + 8 n^2 (up to 2 n slabs of 4 n doubles each). */
int64_t kc_replay_work_doubles(int n, int64_t stack_cap);

/* Start the replay of a chart.  has_mu = 0 means the chart has no "mu" entry (the reference then reads no hint
 * at a linear leaf, and hint values are empty); mu may then be NULL with nmu = 0.  work has
 * kc_replay_work_doubles(n, stack_cap) doubles.  Returns 0, or KC_RERR_LIMIT. */
int kc_replay_init(kc_replay *r, int n, const double *lo0, const double *hi0,
                   const uint8_t *ev, int64_t nev, const uint16_t *q, int64_t nq,
                   const float *mu, int64_t nmu, int has_mu, double *work, int64_t stack_cap);

/* The next node: KC_REPLAY_NODE (written to *out), KC_REPLAY_DONE, or an error code.  After DONE or an error,
 * every further call returns the same value. */
int kc_replay_next(kc_replay *r, kc_node *out);

/* sizeof(kc_replay), for the ctypes front end */
int64_t kc_replay_sizeof(void);

#endif
