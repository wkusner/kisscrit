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

**Status:** the joins of the co-index-1 saddles are checked in double precision with an assumed rounding margin (1e-13) by `scripts/certify_joins.py` (untrusted; not interval arithmetic). The cells of dimension 2 and 3, and so b1, b2 and b3 just above r = 1, are not determined. For 1 < r < 1.0002111 the Morse inequalities give only b2 ≥ 12!/60 − 2.

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
- *Points on the sphere:* N ≤ 8, using the numerical census (lower bounds, see `docs/small-n.md`) in `results/census_n*.json`.

**References:**
- Martiniani, Schrenk, Ramola, Chakraborty, Frenkel, *Numerical test of the Edwards conjecture shows that all packings are equally probable at jamming*, Nature Physics (DOI 10.1038/nphys4168; arXiv:1610.06328). Title and authors checked against the arXiv page 2026-10-23; volume and year not confirmed.
- R. K. Bowles and S. S. Ashwin, *Edwards entropy and compactivity in a model of granular matter* (arXiv:1102.0352; DOI 10.1103/PhysRevE.83.031302). Title and authors checked against the arXiv page 2026-10-23. The abstract describes hard discs between hard walls and does not mention a channel; the earlier description "on quasi-1D channels" was not supported.
- Y. Zhang, M. J. Godfrey and M. A. Moore, *Marginally jammed states of hard disks in a one-dimensional channel* (arXiv:2005.02901; DOI 10.1103/PhysRevE.102.042614). Title and authors checked against the arXiv page 2026-10-23.

## 5. Sticky growth in flat geometry

**Observation:** a free cluster in the plane can always grow by homothety, v_i = x_i / r. So σ·b = 0 for every self-stress σ, and sticky growth never breaks down inside a free flat cluster. In flat containers and on tori, every breakdown involves the boundary or a cluster that wraps around.

**Next:** check this numerically with the container code, counting the release events that touch the boundary.

## 6. Stored symmetry of `data/n12_above1.json`, entry 6 (r = 1.0003062547) (added 2026-10-22, CLEAN-A)

**Question for Wöden:** the stored symmetry of entry 6 looks wrong. Is it, and should the data file be corrected (the file was not edited)?

**Evidence (computed by the DOM-1 session, `docs/session-log.md`, networkx 3.7, numpy 2.5.3; first noticed by the EXACT-1 session, open question 2; neither reproduced in the CLEAN-A session):**
- The entry stores `sym: 2, sym_proper: 1`, but its stored 20-contact graph has |Aut| = 1.
- `kisscrit.symmetry.group_orders(X)` at its default tolerance 1e-6 gives (1, 1). The stored X satisfies its contacts to 2.8e-16 and has only the trivial isometry at Gram tolerances up to 3e-4.
- The stored `sym: 2, sym_proper: 1` appears only at tolerance 1e-3 or more: a reflection (0 5)(2 11)(4 6)(8 10) with maximal point residual 5.3e-4. It swaps the contact (9, 11) with the non-contact (2, 9), whose Gram gap is 6.1e-4.
- The other 11 entries: stored `sym` equals `group_orders` at 1e-6 and equals |Aut| of the stored contact graph.

**Reading (inferred):** a loose tolerance produced the stored value; contacts and coordinates are consistent. Not decided: whether the true critical configuration at this r has a mirror symmetry (the two nearly equal 19- and 21-contact neighbours do).

**AUD-2 (2026-10-23, sync document S-07 and review R6; scratch computation, floating point with stated tolerances, networkx 3.7):** reproduced DOM-1's numbers (|Aut| = 1; `group_orders` (1, 1) at 1e-6, 1e-4 and 3e-4; (2, 1) only at tolerance ≥ 1e-3; smallest non-contact gap 6.07e-4 at pair (2, 9)) and adds: by the Krawczyk enclosure of `results/certified_n12.json` (radius 1e-40; untrusted script `scripts/certify_n12.py`, interval arithmetic) no configuration in the box has a reflection symmetry, since a reflection would move points by about 5e-4. So |Sym| = 1 for the true configuration as well (computed with an untrusted interval script; not proof).

**Resolved 2026-10-23 (decision D-033, Wöden accepted the AUD-2 fixes):** `data/n12_above1.json` entry 6 now stores `sym: 1, sym_proper: 1` (commit DOCS-1-9). Still to follow: the paper's table row (|Sym| = 1) and the generated text of `data/exact/critical-configs-n12.md` and `scripts/exact/table_md.py`, which still say the data file gives order 2.

**Why it matters:** the paper's table of the twelve configurations, `data/n12_above1.json`, and any symmetric-subspace method applied to this value.
**Next:** recompute `sym` for this entry at tolerance 1e-6 and decide; a code session if the file changes (Seed in the DOM-1 log: recompute `sym` -> twelve).

## 7. `scripts/hull12`: inputs that a clean clone does not have (added 2026-10-23, deferred by the maintainer)

**Status: deferred.** The maintainer decided on 2026-10-23 to leave this as it is until a deeper review of `scripts/hull12`; no fix is planned before that review.

**What was seen (read from the code and the tracked file list in the CLEAN-C session, 2026-10-23; the scripts were not run):**
- `flipwalk.py`, `sample_hulls.py`, `seeded.py`, `survey.py` and `thetastar.py` open `tri12.json` by a path relative to the working directory. No file of that name is tracked. `scripts/hull12/README.md` says it is `data/triangulations_n12.json`; `cmp` shows that file is byte-identical to `data/completeness/triangulations12.json` (computed here, 2026-10-23).
- `flipwalk.py` also opens `seeded.json` (tracked as `results/seeded.json`) and `ts_a.json` (survey output, not tracked); `survey.py` reads `ts_a.json` and `ts_b.json`.
- These scripts are marked untrusted (exploration, no part of the proof). As they stand they do not run from a clean clone; `sample_hulls.py` would after `tri12.json` is placed in the working directory (not tested).
- CLEAN-C-1 only changed how the two scripts find the repository root; it did not touch these file names.

**To decide after the review:** whether the scripts take the paths as arguments, read the tracked data file directly, or are moved to `attic/`.
