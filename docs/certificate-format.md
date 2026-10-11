# Certificate format, version 1

This document specifies the proof tree that the completeness search writes for each candidate stressed graph. It is everything an independent checker needs to re-verify the search.

2026-10-11. Written by `kisscrit/completeness/runner.py` (`run(..., certificate=True)`) and decoded by `kisscrit/core/certificate.py`.

## 1. What a certificate claims

For a candidate stressed graph G (a plane graph with 3-, 4- and 5-gon faces) and a range θ ∈ [θ_lo, θ_hi] of contact angles, the certificate partitions the face-parameter space of G, region by region, into boxes. Each box is one of:

- **excluded:** no realization of G in that box carries a strictly positive balanced stress;
- **terminal:** a small box left for identification with a known critical configuration;
- **open:** not examined, because the box budget ran out.

The certificate is a proof for G exactly when every chart is complete and there are no open leaves. Terminal leaves then still need identification (checklist item C4).

## 2. File layout

Certificates live one per line in `results/completeness/certificates/certificates-<route>-<k>-of-<K>.jsonl`:

```
{"code": <canonical code of G>, "certificate": {
   "version": 1,
   "theta": [θ_lo, θ_hi],                  radians
   "tol": <terminal width>,
   "faces": [[v0, v1, ...], ...],          oriented face cycles (counterclockwise seen from outside)
   "complete": true | false,
   "charts": [ {
       "modes": {"<face index>": "a" | "b" | "c", ...},   pentagon regions (model.py docstring)
       "lo0": [...], "hi0": [...],                       starting box (floats, exact JSON round trip)
       "events": <base64 of zlib of uint8 array>,
       "q":      <base64 of zlib of little-endian uint16 array>,
       "n_events": <int>,
       "complete": true | false
     } or {"modes": {...}, "complete": false, "not_run": true}, ... ] } }
```

Parameter vector p:
- p[0] = θ;
- one parameter per quadrilateral: the corner at boundary positions 0 and 2;
- two per pentagon: the corners at formula positions 1 and 2 of its chart. Region 'a' uses boundary positions 1 and 2, region 'b' uses 3 and 4, region 'c' uses 4 and 0.

The parameter index of each face follows the face order. A pentagon whose boundary position 4 is a vertex of degree 2 has that straight corner fixed at π in the starting box, so lo0 = hi0 = π for that coordinate.

## 3. The proof tree

`events` lists the nodes of the tree in depth-first preorder. Every node is a box. Processing a node:

1. **Contractions (zero or more).** Event `64 + k` sets coordinate k to the grid interval
   [g(qlo), g(qhi)], with `qlo = q[i]` and `qhi = 65536 − q[i+1]`. The next two entries of the `q` array are consumed. The grid is relative to the box at the start of this node (before any contraction of this node):
   - g(0) = lo;
   - g(65536) = hi;
   - otherwise g(q) = lo + (hi − lo) · (q / 65536), evaluated in IEEE double precision in exactly this order.

   The region removed by the node's contractions must itself be excluded. It is the original box minus the contracted box, which the checker splits into at most 2n slabs.
2. **Then exactly one of:**

   | event | meaning |
   |---|---|
   | 1 | excluded: **invalid** — no point of the box describes convex equilateral faces in this chart (or violates a lemma used as a box condition; see `docs/lemmas.md`) |
   | 2 | excluded: **vertex sum** — some vertex angle sum cannot equal 2π on the box |
   | 3 | excluded: **linear** — the linearized vertex equations with interval remainder (mean-value form at the box midpoint) have no solution in the box |
   | 4 | excluded: **Gordan**, with a vector inherited from an ancestor |
   | 5 | excluded: **Gordan** — some y with y_u·t_u(e) + y_w·t_w(e) > 0 for every edge on the whole box |
   | 6 | excluded: **Motzkin** — every edge term is either exactly 0 (both end vectors are zero) or > 0 on the box, and at least one is > 0 |
   | 7 | **terminal** |
   | 8 | **open** |
   | 16 + k | **split** coordinate k at m = 0.5·(lo_k + hi_k), in double precision. The subtree of the upper child [m, hi_k] follows, then that of the lower child [lo_k, m]. |

Evaluated exactly as stated, the splits and contractions leave no gaps or overlaps (up to shared faces). Coverage therefore follows from the tree structure alone. `kisscrit/core/certificate.py:replay` reconstructs every node box exactly.

## 4. What the checker must verify

For every leaf and every removed slab, the checker re-derives the stated exclusion from its own interval evaluation of the model. It recomputes multipliers and Gordan or Motzkin vectors itself. The search's multipliers are not stored: they come from floating-point linear programs and carry no rigour.

The checker may replace an exclusion with a different one, or subdivide a box and exclude the pieces. It may not accept a box without an interval proof. Contraction slabs carry no reason code. Any exclusion is acceptable for them, typically the linear test on the parent box, which the search used.

## 5. Assumptions the tree does not record

These are fixed by the code version, which is the tag of the run:

- the candidate list and its completeness — candidate generation, and the lemmas in `docs/lemmas.md`;
- the pentagon region cover — Lemma P, `scripts/completeness/lemma_pentagon.py`;
- the face-angle model itself (`kisscrit/completeness/model.py`, and the compiled equivalent in `fast.py`).

## 6. Size

Measured on pilot candidates: 1.5–3.7 bytes per box after compression, typically 3–15 kB per candidate, and about 100 kB for the hardest known critical graph.
