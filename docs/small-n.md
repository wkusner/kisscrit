# Critical configurations for small N (d = 3)

*2026-10-10. Numerical, in double precision. Reproduce with `scripts/graph_census.py`, `scripts/surgery_census.py` and `scripts/euler.py`. The data are in `results/census_n*.json`.*

## Method

For N ≤ 7 we took every graph on N vertices with minimum degree ≥ 2, up to isomorphism (the networkx atlas). For each graph we solved the square KKT system from 20–30 starts, half random and half from a penalty seed (`search.graph_start`). A solution was kept only if every edge carries positive stress and every non-edge is strictly longer. Contact-graph surgery was then run on the result.

Each type is then given a rattler stratum for every M > N. Its number of components is the number of components of the region where a rattler fits, which is numerical.

The signed count is

  Σ (−1)^{co-index} · #labeled copies · #rattler components.

It is compared with χ(M_{0,N}) = ∏_{k=2}^{N−2}(1−k). For N ≥ 3, BConf(N) at small r is homotopy equivalent to M_{0,N}.

## Results

| N | r | contacts | co-index | labeled copies in BConf(N) | description |
|---|---|---|---|---|---|
| 3 | 3 + 2√3 = 6.464102 | 3 | 0 | 1 | equilateral triangle on a great circle |
| 4 | 2 + √6 = 4.449490 | 6 | 0 | 2 | regular tetrahedron |
| 4 | 1 + √2 = 2.414214 | 4 | 1 | 3 | square on a great circle |
| 5 | 1 + √2 | 6–8 | — | — | degenerate maximum set (see below) |
| 5 | 1.425920 | 5 | 2 | 12 | regular pentagon on a great circle |
| 6 | 1 + √2 | 12 | 0 | 30 | octahedron |
| 6 | (3 + √21)/4 = 1.895644 | 9 | 1 | 120 | trigonal prism |
| 6 | (pentagon + rattler) | 5 | 2 | 12 · 6 · 2 = 144 | Morse–Bott: the rattler sits in either polar cap |
| 6 | 1 | 6 | 3 | 60 | regular hexagon on a great circle |
| 7 | 1.691321 | 12 | 0 | 1680 | the maximum |
| 7 | 1.665320 | 11 | 1 | 2520 | |
| 7 | 1.348026 | 10 | 2 | 2520 | found by surgery, not by the graph solve |
| 7 | 1.323194 | 9 | 3 | 1680 | |
| 7 | (hexagon + rattler) | 6 | 3 | 60 · 7 · 2 = 840 | |
| 7 | (pentagon + 2 rattlers) | 5 | 2 | 12 · 42 = 504 | one rattler in each polar cap |
| 7 | 0.766422 | 7 | 4 | 360 | regular heptagon on a great circle |

Signed counts:

| N | signed count | χ(M_{0,N}) |
|---|---|---|
| 4 | 2 − 3 = −1 | −1 |
| 5 | −10 + 12 = 2 | 2 |
| 6 | 30 − 120 + 144 − 60 = −6 | −6 |
| 7 | 1680 − 2520 + 2520 − 1680 − 840 + 504 + 360 = 24 | 24 |

**N = 5.** The set where r = 1 + √2 (minimum angle 90°) is not a point. By Rankin's argument it contains an antipodal pair: two poles and three equatorial points with gaps ≥ 90°. For each labeled choice this is a triangle of configurations, and there are 20 such triangles in BConf(5). Their corners are the 30 square pyramids, and each corner lies on exactly two triangles. So the maximum set is a 2-complex of 20 triangles glued at 30 vertices, with χ = 20 − 30 = −10. With the 12 pentagons, the count is 2 = χ(M_{0,5}).

The soft-min picture in `research-morse-complex-small-n` (30 maxima, 60 edge saddles, 20 bipyramids) is a perturbation of this degenerate set. The min-type function has no isolated saddles there.

## What is not done

- **N = 8.** The atlas stops at 7 vertices. Ascent plus single-edge surgery found only the maximum (the square antiprism, r = 1.549598). Saddles need either multi-edge surgery or a graph generator (nauty `geng`). There are about 11,000 graphs with minimum degree ≥ 2 on 8 vertices, roughly 4–8 CPU-hours at the current speed.
- **Counts are not certified.** They are lower bounds closed by the Euler characteristic: a missing pair of types with opposite signs would not show up in the check.
- **Labeled counts assume trivial stabilizers.** N!/|G ∩ SO(3)| is correct only if no rattler region is fixed by a symmetry. The polar caps of the pentagon and hexagon are handled by hand above.
