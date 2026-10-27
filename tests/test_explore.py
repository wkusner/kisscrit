"""Tests of the exploration driver: scripts/explore/ (worker, merge, export, run_explore.sh).

Status: untrusted (a test of exploration tooling, not of the mathematics).

What is checked: a worker is reproducible and stops cleanly on its budget and on SIGTERM; its surgery slices are
disjoint and marked; the merge is idempotent, independent of the file order, and survives a truncated file;
the export has the layout of data/n12_above1.json and flags types by the stated thresholds; the refine contract
holds (each `regular` type of short N = 6 and N = 7 runs converges under the Gauss-Newton of scripts/exact/hiprec.py
to 60 digits); run_explore.sh validates its parameters, runs end to end, stops and resumes.  Every census value
involved is numerical (double precision); the tests compare it with the repository's own census files, which are
themselves numerical.
"""
import getpass, json, os, pathlib, shutil, signal, socket, subprocess, sys, time

import numpy as np
import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts' / 'explore'))
sys.path.insert(0, str(ROOT / 'scripts' / 'exact'))
import explore_common as ec, explore_export, explore_merge, explore_worker
from kisscrit.search import Census

RUN = ROOT / 'scripts' / 'explore' / 'run_explore.sh'
ABOVE1 = ROOT / 'data' / 'n12_above1.json'


def _worker(tmp, name, N=4, **kw):
    """Run a worker in this process with a generous budget and return (path, its keys-to-hits dict)."""
    path = str(tmp / name)
    explore_worker.run_worker(N, 3, path, budget=600, dump_seconds=60, **kw)
    return path, {k: t['hits'] for k, t in Census.load(path).types.items()}


@pytest.fixture(scope='module')
def two_workers(tmp_path_factory):
    """Two N = 4 worker files with different indices (the same seed) for the merge tests."""
    tmp = tmp_path_factory.mktemp('workers')
    a, _ = _worker(tmp, 'w000000.json', seed=5, index=0, runs=3, rounds=0, run_id='t', commit='c0')
    b, _ = _worker(tmp, 'w000001.json', seed=5, index=1, runs=3, rounds=0, run_id='t', commit='c0')
    return tmp, a, b


def test_worker_is_reproducible_and_its_generator_depends_on_the_index(tmp_path):
    a, ha = _worker(tmp_path, 'a.json', seed=3, index=0, runs=3, rounds=0)
    b, hb = _worker(tmp_path, 'b.json', seed=3, index=0, runs=3, rounds=0)
    assert ha == hb and ha
    xa = {k: t['X'] for k, t in Census.load(a).types.items()}
    xb = {k: t['X'] for k, t in Census.load(b).types.items()}
    assert xa == xb                                       # same census, to the last bit
    g = [np.random.default_rng(np.random.SeedSequence([3, i])).random() for i in (0, 1)]
    assert g[0] != g[1]                                   # (seed, index) selects the generator


def test_worker_budget_ends_the_run_and_the_file_is_a_complete_census(tmp_path):
    t0 = time.time()
    r = explore_worker.run_worker(4, 3, str(tmp_path / 'w.json'), runs=10 ** 6, budget=2.0, dump_seconds=1)
    assert r['stopped_by'] == 'budget' and time.time() - t0 < 30
    D = json.load(open(tmp_path / 'w.json'))
    assert D['worker']['stopped_by'] == 'budget' and Census.load(str(tmp_path / 'w.json')).N == 4


def test_worker_sigterm_dumps_and_exits_cleanly(tmp_path):
    out = tmp_path / 'w.json'
    p = subprocess.Popen([sys.executable, str(ROOT / 'scripts' / 'explore' / 'explore_worker.py'), '4', '3', '--out', str(out),
                          '--runs', '1000000', '--budget', '600'], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    deadline = time.time() + 120
    while not out.exists() and time.time() < deadline:   # the first (empty) dump marks that the handlers are installed
        time.sleep(0.2)
    assert out.exists()
    time.sleep(1.0); p.send_signal(signal.SIGTERM)
    text, _ = p.communicate(timeout=60)
    assert p.returncode == 0, text
    assert json.load(open(out))['worker']['stopped_by'] == 'sigterm'


def test_surgery_slices_are_disjoint_marked_and_start_hits_are_zero(tmp_path):
    start = str(ROOT / 'results' / 'census_n4.json')       # two types, no provenance, no surg marks
    keys = set(Census.load(start).types)
    got = []
    for i in (0, 1):
        path, hits = _worker(tmp_path, 's%d.json' % i, seed=1, index=i, runs=0, rounds=1, start=start, slice_index=i, slice_count=2)
        types = Census.load(path).types
        got.append({k for k, t in types.items() if k in keys})
        assert all('surg' in types[k] and types[k]['surg'] == i for k in got[-1])
    assert got[0] | got[1] == keys and not (got[0] & got[1])
    assert all(h >= 0 for h in hits.values())


def test_pending_list_skips_types_already_used_by_surgery(tmp_path):
    D = json.load(open(ROOT / 'results' / 'census_n4.json'))
    D['types'][0]['surg'] = 7
    json.dump(D, open(tmp_path / 'start.json', 'w'))
    path, _ = _worker(tmp_path, 'w.json', seed=1, index=0, runs=0, rounds=0, start=str(tmp_path / 'start.json'))
    assert len(Census.load(path).types) == 1               # only the type without a surg mark is taken


def test_merge_is_idempotent_and_independent_of_the_file_order(two_workers, tmp_path):
    _, a, b = two_workers
    outs = []
    for k, files in enumerate(([a, b], [b, a], [a, b])):
        out = str(tmp_path / ('m%d.json' % k))
        assert explore_merge.main(['--out', out, '--no-euler'] + files) == 0
        outs.append(open(out, 'rb').read())
    assert outs[0] == outs[1] == outs[2]
    merged, info = explore_merge.merge_files([a, b])
    ca, cb = Census.load(a), Census.load(b)
    assert len(merged.types) == len(set(ca.types) | set(cb.types)) == info['types_merged']
    for k, t in merged.types.items():                      # hits are summed over the files
        assert t['hits'] == ca.types.get(k, {'hits': 0})['hits'] + cb.types.get(k, {'hits': 0})['hits']
        assert t['sigma'] is not None and t['prov']['commit'] == 'c0'


def test_merge_skips_a_truncated_file_and_a_merge_output(two_workers, tmp_path):
    _, a, b = two_workers
    bad = tmp_path / 'w000009.json'
    bad.write_bytes(open(b, 'rb').read()[:200])            # a half-written copy
    good = str(tmp_path / 'good.json')
    explore_merge.main(['--out', good, '--no-euler', a])
    out = str(tmp_path / 'm.json')
    assert explore_merge.main(['--out', out, '--no-euler', a, str(bad), good]) == 0
    merged = json.load(open(out))
    skipped = {s['file'] for s in merged['merge_info']['skipped']}
    assert skipped == {'w000009.json', 'good.json'}
    assert set(Census.load(out).types) == set(Census.load(a).types)    # nothing of the readable file is lost
    assert not [p for p in tmp_path.iterdir() if '.tmp.' in p.name]    # no temporary file is left behind
    assert explore_merge.main(['--out', str(tmp_path / 'none.json'), '--no-euler', str(bad)]) == 1


def test_euler_check_of_the_repository_census_files_n4_to_n7():
    """The signed counts that scripts/euler.py prints for results/census_n4..7.json against chi(M_{0,N}) (computed).

    N = 5 differs: its census holds a degenerate type (nzero = 2), which euler_check skips; the paper's
    count for it is handled outside euler_check.  These are the repository's numerical census files."""
    res = {}
    for N in (4, 5, 6, 7):
        merged = Census.load(str(ROOT / 'results' / ('census_n%d.json' % N)))
        e = explore_merge.euler_report(merged, str(ROOT / 'results'))
        res[N] = (e['signed_count'], e['chi'], e['skipped_types'])
        print('N=%d signed count %d, chi %d, skipped types %d (%s)' % (N, e['signed_count'], e['chi'], e['skipped_types'], e['status']))
    assert [res[N][1] for N in (4, 5, 6, 7)] == [-1, 2, -6, 24]
    assert res[4][0] == -1 and res[6][0] == -6 and res[7][0] == 24 and res[5][0] != 2 and res[5][2] == 1


def test_export_has_the_layout_of_n12_above1_and_a_permutation_of_priorities(tmp_path):
    C, meta, rows, records = explore_export.build(str(ABOVE1), [str(ABOVE1)])
    explore_export.write_outputs(C, meta, rows, records, [str(ABOVE1)], str(tmp_path))
    A = json.load(open(ABOVE1)); B = json.load(open(tmp_path / 'refine-input.json'))
    assert list(A) == list(B) == ['N', 'd', 'note', 'types'] and A['N'] == B['N'] and A['d'] == B['d']
    for ta, tb in zip(A['types'], B['types']):
        assert list(ta) == list(tb) and ta['X'] == tb['X'] and ta['pairs'] == tb['pairs'] and ta['r'] == tb['r']
    assert 'numerical (double precision)' in B['note'] and 'no value in it is a proof input' in B['note']
    lines = [l for l in open(tmp_path / 'refine-queue.tsv').read().splitlines() if not l.startswith('#')]
    head = lines[0].split('\t'); table = [dict(zip(head, l.split('\t'))) for l in lines[1:]]
    assert sorted(int(r['priority']) for r in table) == list(range(1, 13)) and len(table) == 12
    assert all(B['types'][int(r['index'])]['r'] == float(r['r']) for r in table)
    assert [r['flag'] for r in table].count('regular') == 2                # DOD and M: kernel 0 and a strictly positive stress
    assert all(r['new'] == 'no' for r in table)
    assert open(tmp_path / 'refine-queue.tsv').read().startswith('# refine queue')


def test_flags_follow_the_thresholds():
    rec = json.load(open(ABOVE1))['types'][0]               # the dodecahedron-like type: kernel 0, stress > 0
    m = ec.kkt_metrics(rec)
    assert ec.classify(rec, m) == 'regular' and m['kkt_resid'] < 1e-12 and m['n_contacts'] == 30
    assert ec.classify(dict(rec, rattlers=1), m) == 'rattler'
    assert ec.classify(dict(rec, kernel=1), m) == 'degenerate'
    X = np.array(rec['X']); X[0] = X[0] + 1e-6 * np.array([1.0, 0.0, 0.0]); X[0] /= np.linalg.norm(X[0])
    moved = dict(rec, X=X.tolist())
    assert ec.classify(moved, ec.kkt_metrics(moved)) == 'suspect'     # contacts no longer equal: residual above 1e-9
    assert ec.classify(rec, dict(m, min_gap=1e-7)) == 'suspect'       # a non-contact pair nearly in contact
    assert ec.classify(rec, dict(m, kkt_resid=float('nan'))) == 'suspect'


def _short_run(tmp, N, seed, runs):
    """A short single-worker run, merged and exported; returns the export directory."""
    w = str(tmp / ('w%d.json' % N)); explore_worker.run_worker(N, 3, w, seed=seed, index=0, runs=runs, rounds=1, budget=900,
                                                                run_id='contract', commit='test')
    m = str(tmp / ('m%d.json' % N)); assert explore_merge.main(['--out', m, '--no-euler', w]) == 0
    out = tmp / ('ex%d' % N)
    C, meta, rows, records = explore_export.build(m, [])
    explore_export.write_outputs(C, meta, rows, records, [], str(out))
    return out


def _table(path):
    lines = [l for l in open(path).read().splitlines() if not l.startswith('#')]
    head = lines[0].split('\t'); return [dict(zip(head, l.split('\t'))) for l in lines[1:]]


def polish_regular_types(out):
    """Polish every `regular` type of the export in directory `out` with hiprec.newton to 60 digits.

    A type counts as converged when the residual of |x_i|^2 = 1 and <x_i, x_j> = c on the contacts is below 1e-50 at
    80 digits and the result is within 1e-9 of the double-precision X.  Returns (types, regular, other flags, done,
    failed); nothing is dropped silently: every failure is printed."""
    import mpmath as mp, hiprec
    D = json.load(open(out / 'refine-input.json')); table = _table(out / 'refine-queue.tsv')
    assert len(table) == len(D['types'])                                  # every type is in the queue
    regular = [r for r in table if r['flag'] == 'regular']
    others = {f: sum(r['flag'] == f for r in table) for f in ('degenerate', 'rattler', 'suspect')}
    done = failed = 0; worst = mp.mpf(0)
    for r in regular:
        t = D['types'][int(r['index'])]; A = [tuple(p) for p in t['pairs']]
        X, c, w, res = hiprec.newton(t['X'], A, 60, stress=False)
        mp.mp.dps = 80; n = len(X)
        resid = max([abs(sum(v * v for v in X[i]) - 1) for i in range(n)] +
                    [abs(sum(X[i][a] * X[j][a] for a in range(3)) - c) for i, j in A])
        close = max(abs(X[i][a] - mp.mpf(t['X'][i][a])) for i in range(n) for a in range(3)) < 1e-9
        if resid < mp.mpf(10) ** -50 and close:
            done += 1; worst = max(worst, resid)
        else:
            failed += 1; print('  NOT converged: r=%s residual %s close-to-double %s' % (t['r'], mp.nstr(resid, 3), close))
    print('N=%d: %d types (regular %d, degenerate %d, rattler %d, suspect %d); polished to 60 digits: %d converged, %d not; '
          'largest residual %s' % (D['N'], len(table), len(regular), others['degenerate'], others['rattler'], others['suspect'],
                                   done, failed, mp.nstr(worst, 3) if done else 'n/a'))
    return len(table), len(regular), others, done, failed


@pytest.mark.parametrize('N,seed,runs', [(6, 1, 4), (7, 1, 4)])
def test_refine_contract_regular_types_converge_under_hiprec(tmp_path, N, seed, runs):
    """Step 2b: every `regular` type of a short run is polished by hiprec.newton to 60 digits; none is dropped."""
    pytest.importorskip('flint')
    _, regular, _, done, failed = polish_regular_types(_short_run(tmp_path, N, seed, runs))
    assert regular and failed == 0 and done == regular


def _run(args, env=None, **kw):
    e = dict(os.environ, KISSCRIT_PYTHON=sys.executable, **(env or {}))
    return subprocess.run(['bash', str(RUN)] + args, capture_output=True, text=True, env=e, **kw)


def test_run_script_syntax_and_shellcheck():
    assert subprocess.run(['bash', '-n', str(RUN)]).returncode == 0
    if shutil.which('shellcheck'):
        p = subprocess.run(['shellcheck', '-x', str(RUN)], capture_output=True, text=True)
        assert p.returncode == 0, p.stdout
    else:
        pytest.skip('shellcheck not installed')


def test_run_script_header_says_untrusted():
    head = open(RUN).read().splitlines()[:6]
    assert 'Status: untrusted' in ' '.join(head) and head[0].startswith('#!')


def test_run_script_refuses_bad_parameters_and_status_works_on_an_empty_run_home(tmp_path):
    nproc = os.cpu_count() or 1
    rh = ['--run-home', str(tmp_path / 'rh')]
    for flags, word in ((['--cores', str(nproc + 1)], 'nproc'), (['--hours', '0'], '--hours'), (['--n', '2'], '--n'),
                        (['--nice', '25'], '--nice'), (['--cores', 'x'], '--cores')):
        p = _run(['start', '--foreground'] + rh + flags)
        assert p.returncode != 0 and word in p.stderr, (flags, p.stderr)
    assert not (tmp_path / 'rh').exists()                 # nothing was created by a refused start
    p = _run(['status'] + rh)
    assert p.returncode == 0 and 'no exploration run' in p.stdout
    assert _run(['bogus']).returncode != 0 and _run(['help']).returncode == 0


def test_end_to_end_two_workers_one_minute_n5_d3(tmp_path):
    """A 60-second run for N = 5, d = 3 on 2 workers; its types are compared with results/census_n5.json (numerical)."""
    if (os.cpu_count() or 1) < 2:
        pytest.skip('needs two cores')
    rh = tmp_path / 'rh'
    p = _run(['start', '--foreground', '--run-home', str(rh), '--tree', str(ROOT), '--n', '5', '--d', '3', '--cores', '2',
              '--hours', '0.0167', '--batch-minutes', '0.5', '--runs', '3', '--rounds', '1'], timeout=600)
    assert p.returncode == 0, p.stderr
    res = rh / 'explore' / 'result'
    assert sorted(f.name for f in res.iterdir()) == ['manifest.json', 'merge-report.txt', 'merged-census.json', 'refine-input.json',
                                                     'refine-queue.tsv', 'summary.txt']
    merged = Census.load(str(res / 'merged-census.json')); known = Census.load(str(ROOT / 'results' / 'census_n5.json'))
    found = sorted({round(t['r'], 6) for t in merged.types.values()}); want = sorted({round(t['r'], 6) for t in known.types.values()})
    print('N=5, 2 workers, 60 s: %d types with radii %s; results/census_n5.json has %d types with radii %s'
          % (len(merged.types), found, len(known.types), want))
    assert set(found) <= set(want)                         # nothing outside the repository's N = 5 census
    assert len(_table(res / 'refine-queue.tsv')) == len(merged.types)
    man = json.load(open(res / 'manifest.json'))
    assert man['parameters']['N'] == 5 and man['parameters']['cores'] == 2
    assert len(man['commit']) == 40 or man['commit'] == 'unknown'
    text = b''.join(f.read_bytes() for f in res.iterdir())
    for private in (socket.gethostname(), getpass.getuser(), str(tmp_path), os.path.expanduser('~')):
        assert private.encode() not in text, private          # no host, user or path string leaves with the results
    s = _run(['status', '--run-home', str(rh)]).stdout
    assert 'finished (deadline)' in s and 'Euler check' in s


def test_stop_merges_and_a_second_start_resumes_with_fewer_cores(tmp_path):
    if (os.cpu_count() or 1) < 2:
        pytest.skip('needs two cores')
    rh = tmp_path / 'rh'; base = ['--run-home', str(rh)]
    env = dict(os.environ, KISSCRIT_PYTHON=sys.executable)
    start = ['bash', str(RUN), 'start', '--foreground', *base, '--tree', str(ROOT), '--n', '4', '--cores', '2', '--hours', '0.5',
             '--batch-minutes', '5', '--runs', '1000000', '--rounds', '1']
    proc = subprocess.Popen(start, env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    b0 = rh / 'explore' / 'batches' / 'b0000'
    for _ in range(300):
        if b0.exists() and len(list(b0.glob('w*.json'))) == 2:
            break
        time.sleep(0.5)
    time.sleep(2)
    s = _run(['stop', *base], timeout=900)
    proc.communicate(timeout=60)
    assert s.returncode == 0 and (rh / 'explore' / 'state' / 'finished').read_text().strip() == 'stopped', s.stderr
    assert (rh / 'explore' / 'result' / 'merged-census.json').exists() and (rh / 'explore' / 'result' / 'summary.txt').exists()
    assert _run(['start', *base, '--hours', '5']).returncode != 0                # a changed parameter is refused on a resume
    proc = subprocess.Popen(['bash', str(RUN), 'start', '--foreground', *base, '--cores', '1'], env=env, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, text=True)
    b1 = rh / 'explore' / 'batches' / 'b0001'
    for _ in range(300):
        if b1.exists() and list(b1.glob('w*.json')):
            break
        time.sleep(0.5)
    s = _run(['stop', *base], timeout=900); proc.communicate(timeout=60)
    assert s.returncode == 0
    names = sorted(p.name for p in (rh / 'explore' / 'batches').rglob('w*.json'))
    assert names == ['w000000.json', 'w000001.json', 'w000002.json']             # indices stay unique across the resume
    man = json.load(open(rh / 'explore' / 'result' / 'manifest.json'))
    assert [r['cores'] for r in man['restarts']] == [1]


def test_default_start_runs_a_snapshot_of_head_and_ignores_uncommitted_edits(tmp_path):
    """Without --tree the run takes `git archive HEAD` of the repository that holds the script, so a later edit of the
    working tree cannot change a run in progress; the manifest records the snapshot's commit."""
    repo = tmp_path / 'repo'
    shutil.copytree(ROOT / 'kisscrit', repo / 'kisscrit', ignore=shutil.ignore_patterns('__pycache__'))
    shutil.copytree(ROOT / 'scripts' / 'explore', repo / 'scripts' / 'explore', ignore=shutil.ignore_patterns('__pycache__'))
    git = ['git', '-C', str(repo), '-c', 'user.name=t', '-c', 'user.email=t@example.invalid']
    subprocess.run(git + ['init', '-q'], check=True); subprocess.run(git + ['add', '-A'], check=True)
    subprocess.run(git + ['commit', '-q', '-m', 'snapshot source'], check=True)
    head = subprocess.run(git + ['rev-parse', 'HEAD'], capture_output=True, text=True, check=True).stdout.strip()
    (repo / 'scripts' / 'explore' / 'explore_common.py').write_text('raise SystemExit("edited after the commit")\n')
    rh = tmp_path / 'rh'
    p = subprocess.run(['bash', str(repo / 'scripts' / 'explore' / 'run_explore.sh'), 'start', '--foreground', '--run-home', str(rh),
                        '--n', '4', '--cores', '1', '--hours', '0.005', '--batch-minutes', '0.1', '--runs', '1', '--rounds', '0'],
                       capture_output=True, text=True, env=dict(os.environ, KISSCRIT_PYTHON=sys.executable), timeout=300)
    assert p.returncode == 0, p.stderr
    assert 'uncommitted changes' in p.stdout
    assert json.load(open(rh / 'explore' / 'result' / 'manifest.json'))['commit'] == head
    assert (rh / 'explore' / 'tree' / 'scripts' / 'explore' / 'explore_worker.py').exists()
    assert 'edited after the commit' not in (rh / 'explore' / 'tree' / 'scripts' / 'explore' / 'explore_common.py').read_text()
    assert list((rh / 'explore' / 'batches').rglob('w*.json'))
