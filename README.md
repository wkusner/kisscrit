# kisscrit

Critical configurations of the min-distance (kissing / Tammes) function for N points on the sphere S^{d-1}.

This is the code that supports W. Kusner, *Critical radii of kissing configurations: algebraicity, and twelve spheres above r = 1* (in preparation). It is research code, written with Claude (Anthropic) as an assistant, and it is numerical unless a function says otherwise.

## The problem

For X = (x_1, …, x_N) in (S^{d-1})^N, let c(X) = max_{i<j} ⟨x_i, x_j⟩. The kissing radius r(X) = s/(1−s), with s = sin(θ/2), increases exactly when c decreases. r is a min-type function (Baryshnikov–Bubenik–Kahle). Its critical configurations are the Clarke-critical points of c.

Clarke criticality means the contact graph A (the pairs attaining c) carries a balanced stress σ ≥ 0, σ ≠ 0. At every vertex, Σ_j σ_ij x_j is parallel to x_i.

For a fixed contact graph A, the critical configurations are the solutions of a square polynomial system (`kisscrit.kkt`):

    |x_i|² = 1,   ⟨x_i, x_j⟩ = c  (ij ∈ A),   Σ_j σ_ij x_j = μ_i x_i,   Σ σ_ij = 1,

with σ > 0 and every non-contact cosine strictly below c. The system has N·d + 1 + |A| + N unknowns. It is square after removing the d(d−1)/2 torque identities and the O(d) gauge. So isolated critical radii are algebraic, and their degrees are bounded by Bézout or BKK on this system.

## What is here

| module | contents |
|---|---|
| `geometry` | conventions, tangent bases on S^{d-1}, great-circle moves, r ↔ c |
| `stress` | balanced-stress LP (Clarke criticality), stress-space dimension, rigidity kernel modulo rotations, Lagrangian Hessian, co-index (`analyze`) |
| `kkt` | the square system and a Levenberg–Marquardt solver with analytic Jacobian (any d) |
| `ascent` | monotone ascent of r by LP directions (the flow used for basins) |
| `search` | `Census` of types; random ascent; contact-graph surgery; enumeration over all graphs on n ≤ 7 vertices (`graphs_on`, `solve_graph`); rattler strata; Euler-characteristic check against M_{0,N} |
| `symmetry` | the symmetry group in O(d) and the labeled count N!/\|G ∩ SO(d)\| |
| `exact` | arbitrary-precision Newton on the KKT system (`refine_kkt`); LLL minimal polynomials with a height test (`minpoly`); conversion from c to r |
| `certify` | interval (Krawczyk) certification of a nondegenerate critical configuration on S²: existence and uniqueness in a box, exact contact graph, positive stress, co-index via the bordered Hessian |
| `hull` | (d = 3) enumeration of the triangulations of S² on n vertices by flips; θ*(T), the largest minimum angle of an inscribed convex realization of T (local optimization, so a lower bound) |
| `packing` | the same toolkit for equal balls in containers (unit square/cube, flat tori, disk/ball): ascent, square KKT systems with gauge, co-index, Morse–Bott family and rattler-stratum counting, Euler checks, high-precision refinement |

## Quick start

```bash
pip install -e .            # numpy, scipy, networkx, mpmath; python-flint for minpoly
python -m kisscrit census --N 6 --d 3 --starts 30
python -m kisscrit analyze data/n12_above1.json
pytest -q
```

```python
import numpy as np, kisscrit as kc
from kisscrit.search import graphs_on, solve_graph
rng = np.random.default_rng(0)
C = kc.Census(6, 3)
for G in graphs_on(6):
    for X in solve_graph(G, 3, rng, starts=30):
        C.add(X)
for t in C.sorted():
    print(t['r'], t['coindex'], t['labeled'])
```

## Status of results (numerical unless stated)

- **N = 3 … 6, d = 3.** Every graph with minimum degree ≥ 2 was solved from random starts. Counted with rattler strata, the signed counts give χ(M_{0,N}) for N = 4 and 6. N = 5 has a degenerate two-dimensional maximum set at r = 1 + √2: twenty triangles glued at thirty square pyramids, with χ = −10. With it, N = 5 checks as well. See `docs/small-n.md`.
- **N = 12, d = 3, above r = 1.** There are twelve critical types. Contact graphs and data are in `data/`. The local maximum M at r = 1.0010121… is certified exactly elsewhere (paper, §6).
- **N = 12 certification.** `scripts/certify_n12.py` certifies the ten saddles above r = 1 in about 10 s each (results/certified_n12.json).
- **Containers.** `scripts/pack_explore.py`, `scripts/pack_maxima.py`, `scripts/pack_euler.py`; see `docs/packing.md`.
- **Hull triangulations for N = 12.** All 7595 types were enumerated. The critical configurations above r = 1 have hulls T0 (DOD only), T1, T5 and T7. See `docs/hull.md`.

## Caveats

Numerical KKT solutions from random starts can miss solutions. Treat the counts as lower bounds unless an Euler-characteristic check closes them. Co-indices come from eigenvalue signs in double precision. Degenerate critical sets (`nzero > 0`), and rattler regions that can hold more than one rattler, are reported and not counted.

## License

MIT. See `LICENSE`.
