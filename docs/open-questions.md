# Open questions and future projects for kisscrit

Purpose: a running list of open mathematical questions and planned projects that use this toolkit. Each entry says what is known, where the code and data live, and what the next step would be.

Version 1, 2026-10-11. Add entries at the end, and mark resolved ones rather than deleting them.

## 1. Completeness of the N = 12 census above r = 1

**Status:** open. The certifier and the sharded runner exist: `kisscrit/completeness/`, `docs/completeness-runbook.md`. In pilots it excludes about 98% of candidates. It does not yet resolve:
- low-edge graphs (E = 17–19), whose realization spaces have dimension 22 − E = 3–5;
- the range 1 < r < 1.0001, where many graphs are nearly realizable because of the kissing configurations at r = 1.

**Next:**
- run the gen phase of the conditional route to get the candidate count;
- add first-order (Taylor-model) Gordan certificates;
- find a perturbation argument for the range just above r = 1.

## 2. Higher attaching maps at N = 12

**Status:** the joins of the co-index-1 saddles are certified (`scripts/certify_joins.py`). The cells of dimension 2 and 3, and so b1, b2 and b3 just above r = 1, are not determined. For 1 < r < 1.0002111 the Morse inequalities give only b2 ≥ 12!/60 − 2.

## 3. Exactness of the remaining critical values

**Status:** several of the critical values have no integer relation of degree ≤ 90 at 1600 digits. Settling them needs either an exact elimination on the square systems, or about 10^4 digits.

## 4. Edwards entropy versus basin weights (future project)

**Question:** do the jammed states reached by a compression protocol occur with Edwards' flat measure (each labelled state equally likely), with some other universal measure, or with one that is protocol-specific? Can the Morse structure of the packing function predict basin weights?

**First look:** free planar clusters modulo similarity (2026-10-11).
- *Setup:* Kendall shape space CP^{N−2}, with ρ = min distance / radius of gyration. Jammed states are local maxima of ρ. Starting points are Gaussian, which is uniform on shape space.
- *Code and data:* `scripts/edwards/shape.py`, `run.py`, `analyze.py`; results in `results/edwards/`.
- *Flat measure rejected:* for N = 6–10 under both protocols tested (Morse ascent and a local optimizer). At N = 7 the hexagonal flower takes 99.7% of the Morse runs where Edwards predicts 14%.
- *Compactivity weighting not enough:* the one-parameter weighting p ∝ e^{βρ}/|Aut| fits N = 8 but is rejected at N = 9 and 10.
- *A degenerate case:* at N = 4 the maximum is a one-parameter family of rhombi, i.e. Morse–Bott.
- *Caveat:* everything above is numerical. The test for a local maximum is a perturbation test, and the samples are 300–2000 runs per N.

**Next steps:**
1. Census the saddles for N = 6–9, using the Euler check χ(CP^{N−2}) = N − 1.
2. Test whether a basin's weight is predicted by the saddles on its boundary.
3. Add protocols: sticky growth, soft-potential energy minimization, and a tapping Markov chain.
4. Compare other normalizations: by enclosing circle (the disk container) and by diameter.

**Further settings, all feasible with this toolkit:**
- *Containers:* square, disk and tori. Lists of maxima for N ≤ 10 are already in `docs/packing.md`.
- *Quasi-1D hard disks in a channel:* below width 1 + √3/2 diameters, jammed states can be counted exactly by transfer matrices. Between 1 + √3/2 and 2, next-nearest neighbours interact and the transfer matrix is larger. At width 2 the disks can pass each other and the problem is no longer quasi-1D.
- *Spheres in a narrow cylinder:* columnar packings.
- *Points on the sphere:* N ≤ 8, using the complete census in `results/census_n*.json`.

**References:**
- Martiniani, Schrenk, Ramola, Chakraborty, Frenkel, *Numerical test of the Edwards conjecture shows that all packings are equally probable at jamming*, Nature Physics 2017 (arXiv:1610.06328).
- Bowles and Ashwin on quasi-1D channels (arXiv:1102.0352).
- *Marginally jammed states of hard disks in a one-dimensional channel* (arXiv:2005.02901).

## 5. Sticky growth in flat geometry

**Observation:** a free cluster in the plane can always grow by homothety, v_i = x_i / r. So σ·b = 0 for every self-stress σ, and sticky growth never breaks down inside a free flat cluster. In flat containers and on tori, every breakdown involves the boundary or a cluster that wraps around.

**Next:** check this numerically with the container code, counting the release events that touch the boundary.
