"""Direct tests of kisscrit/completeness/candidates.py (CLEAN-9, survey SURVEY-1; trusted code, tests are not).

What is checked, on triangulation T0 (index 0 of data/completeness/triangulations12.json) and on the pilot's
gen file results/completeness/gen/gen-cond-0-of-4.jsonl.gz, which holds exactly the candidates of T0
(shard 0 of 4 of the `cond` route handles triangulation 0 only):

1. the gen file has the counts it records: one end marker {"done_triangulation": 0, "new": 9689}, 9,689 code
   lines, 9,689 distinct codes, all from triangulation 0;
2. `candidates_of_triangulation` is deterministic and reproduces the gen file: its first K candidates, codes and
   faces, equal the first K code lines of the file in order (the file lists each code at its first appearance);
3. a seeded sample of gen records recomputes to the stored code (canonical code of the rotation system) and is
   `lp_feasible`;
4. `lp_feasible` on one known-feasible case (the icosahedron: 5 triangles at every vertex, alpha = 72 degrees,
   the top of the allowed range) and one known-infeasible case (a triangulation with a degree-6 vertex: six
   triangle corners of size alpha sum to 360 degrees only for alpha = 60 degrees < 70.5288 degrees).

The LP discard does nothing on this prefix (computed at e4ae530: none of the first 175,308 tilings of T0 is
rejected by `lp_feasible`), so case 2 does not test the LP; case 4 does.  Case 2 does test `tilings`, `_vmask`
(mutants of VMASK are detected) and the canonical-code deduplication.  Regenerating all of T0 takes minutes and
is not done here (K below is small); it is not in the suite.  The test reads the files by absolute path, so it
does not depend on the working directory.
These tests skip, and fail under KISSCRIT_REQUIRE_ALL=1 (tests/conftest.py), if the gen file is absent.
"""
import gzip, itertools, json, pathlib, random

import pytest

from kisscrit.completeness import candidates as C
from kisscrit.completeness.canon import canonical_code, rotation_of, code_string
from kisscrit.completeness.structure import Structure, orient_triangulation, merge_faces

ROOT = pathlib.Path(__file__).resolve().parents[1]
GEN = ROOT / 'results' / 'completeness' / 'gen' / 'gen-cond-0-of-4.jsonl.gz'
TRI = ROOT / 'data' / 'completeness' / 'triangulations12.json'

SEED = 20261022      # fixed seed of the sample
K_PREFIX = 500       # candidates of T0 compared in order with the gen file
N_SAMPLE = 150       # gen records recomputed
N_T0_CODES = 9689    # the number of candidates of T0 recorded in the gen file (computed, zcat | wc)


@pytest.fixture(scope='module')
def triangulations():
    return json.load(open(TRI))


@pytest.fixture(scope='module')
def gen():
    if not GEN.exists():
        pytest.skip('no gen file results/completeness/gen/gen-cond-0-of-4.jsonl.gz in the repository')
    recs = [json.loads(line) for line in gzip.open(GEN, 'rt')]
    return [r for r in recs if 'code' in r], [r for r in recs if 'code' not in r]


def test_gen_file_counts(gen):
    codes, markers = gen
    assert markers == [{'done_triangulation': 0, 'new': N_T0_CODES}]
    assert len(codes) == N_T0_CODES
    assert len({r['code'] for r in codes}) == N_T0_CODES
    assert {r['tri'] for r in codes} == {0}


def test_t0_prefix_reproduces_gen_file_in_order(gen, triangulations):
    codes, _ = gen
    got = list(itertools.islice(C.candidates_of_triangulation(triangulations[0]), K_PREFIX))
    assert len(got) == K_PREFIX
    assert [c for c, _ in got] == [r['code'] for r in codes[:K_PREFIX]]
    assert [f for _, f in got] == [r['faces'] for r in codes[:K_PREFIX]]


def test_seeded_sample_recomputes_and_is_lp_feasible(gen):
    codes, _ = gen
    rng = random.Random(SEED)
    sample = rng.sample(codes, N_SAMPLE)
    for r in sample:
        S = Structure(r['faces'])
        assert code_string(canonical_code(rotation_of(S))) == r['code']
        assert all(3 <= k <= 5 for k in S.sizes)
        assert C.lp_feasible(S)


def _all_triangle_structure(F):
    OT = orient_triangulation([tuple(sorted(t)) for t in F])
    return Structure(merge_faces(OT, [(i,) for i in range(len(OT))]))


def _degrees(F):
    d = [0] * 12
    for t in F:
        for v in t: d[v] += 1
    return d


def test_lp_feasible_known_feasible_and_infeasible(triangulations):
    ico = [F for F in triangulations if _degrees(F) == [5] * 12]
    assert len(ico) == 1                                    # the icosahedron is the only 5-regular triangulation
    assert C.lp_feasible(_all_triangle_structure(ico[0]))   # 5 alpha = 360 at alpha = 72 degrees, the top of the range
    deg6 = next(F for F in triangulations if max(_degrees(F)) == 6)
    S = _all_triangle_structure(deg6)
    assert any(sum(1 for (_, fi, _) in S.rot[v] if S.sizes[fi] == 3) == 6 for v in range(S.n))
    assert not C.lp_feasible(S)                             # that vertex forces alpha = 60 degrees < 70.5288 degrees
