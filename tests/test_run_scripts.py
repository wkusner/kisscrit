"""Tests of the run scripts scripts/full/run.sh, run_state.py and run_repair.py (checklist F8, milestone M1).

Status: untrusted (tests of run tooling; no part of the proof).

Everything runs in a temporary HOME with a local bare git repository standing in for GitHub, the system Python
standing in for the conda env (KISSCRIT_PYTHON), dummy token files and fake crontab/systemctl programs.  Nothing
is pushed, no real token exists, and no route is run in the repository tree.  What is really run here: setup
(git clone mode), a self-test with a pytest subset, a route-test run killed and resumed, tmux start/stop.  What
is only dry-run or faked: Miniforge and conda, gh login, systemd (fake systemctl; the unit file is checked with
systemd-analyze verify when it is installed), cron (fake crontab).
"""
import json, os, pathlib, shutil, signal, subprocess, sys, time, gzip

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]
RUN = ROOT / 'scripts' / 'full' / 'run.sh'
HELPER = ROOT / 'scripts' / 'full' / 'run_state.py'
sys.path.insert(0, str(HELPER.parent))
import run_state, run_repair  # noqa: E402

pytestmark = pytest.mark.skipif(sys.platform != 'linux' or not shutil.which('git') or not shutil.which('flock'),
                                reason='run.sh targets Linux and needs git and flock')


def sh(args, env, check=False, timeout=900, **kw):
    """Run bash run.sh ARGS (the copy given as args[0] if it is a path) and return the CompletedProcess."""
    script = str(RUN)
    if args and args[0].endswith('run.sh'): script, args = args[0], args[1:]
    p = subprocess.run(['bash', script] + list(args), env=env, capture_output=True, text=True, timeout=timeout, **kw)
    if check: assert p.returncode == 0, p.stdout + p.stderr
    return p


def snapshot_repo(base):
    """A git repository holding the current working tree (tracked and new files) as one commit, and a bare clone
    of it standing in for GitHub.  Returns (work repo, bare url, commit)."""
    src = base / 'src'; src.mkdir()
    files = subprocess.run(['git', '-C', str(ROOT), 'ls-files', '-z', '--cached', '--others', '--exclude-standard'],
                           capture_output=True, check=True).stdout
    subprocess.run(['tar', '--null', '-T', '-', '-cf', str(base / 'src.tar')], input=files, cwd=ROOT, check=True)
    subprocess.run(['tar', '-xf', str(base / 'src.tar'), '-C', str(src)], check=True)
    g = lambda *a: subprocess.run(['git', '-C', str(src)] + list(a), check=True, capture_output=True, text=True).stdout.strip()
    g('init', '-q', '-b', 'main'); g('add', '-A'); g('-c', 'user.name=t', '-c', 'user.email=t@t', 'commit', '-q', '-m', 'snapshot')
    subprocess.run(['git', 'clone', '-q', '--bare', str(src), str(base / 'remote.git')], check=True)
    return src, 'file://%s' % (base / 'remote.git'), g('rev-parse', 'HEAD')


def base_env(home, **extra):
    """The environment of a run.sh call in the temporary HOME (system Python, fake crontab and systemctl)."""
    env = {k: v for k, v in os.environ.items() if not k.startswith('KISSCRIT_')}
    env.update(HOME=str(home), KISSCRIT_PYTHON=sys.executable, KISSCRIT_TMUX_SOCKET='kc-test-%d' % os.getpid(),
               KISSCRIT_CRONTAB=str(home / 'fakebin' / 'crontab'), KISSCRIT_SYSTEMCTL=str(home / 'fakebin' / 'systemctl'),
               KISSCRIT_LOGINCTL=str(home / 'fakebin' / 'loginctl'), GIT_CONFIG_NOSYSTEM='1')
    env.update(extra)
    return env


def fake_bins(home):
    """crontab keeps its table in ~/.fake-crontab; systemctl logs its arguments; loginctl says Linger=no."""
    fb = home / 'fakebin'; fb.mkdir(parents=True)
    (fb / 'crontab').write_text('#!/bin/bash\nF=$HOME/.fake-crontab\ncase "$1" in -l) [ -f "$F" ] && cat "$F" || '
                                '{ echo "no crontab for $USER" >&2; exit 1; } ;; -) cat > "$F" ;; esac\n')
    (fb / 'systemctl').write_text('#!/bin/bash\necho "$@" >> $HOME/.fake-systemctl.log\n'
                                  'case "$*" in *is-active*) echo inactive; exit 3 ;; esac\nexit 0\n')
    (fb / 'loginctl').write_text('#!/bin/bash\necho no\n')
    for f in fb.iterdir(): f.chmod(0o755)


@pytest.fixture(scope='module')
def world(tmp_path_factory):
    """A temporary HOME after `run.sh setup` (git clone mode) of a snapshot of this working tree, with a passed
    self-test record for its commit written by hand (the real self-test has its own test below)."""
    base = tmp_path_factory.mktemp('runsh')
    src, url, commit = snapshot_repo(base)
    home = base / 'home'; home.mkdir(); fake_bins(home)
    env = base_env(home)
    p = sh(['setup', '--repo-url', url, '--skip-preflight'], env, check=True)
    assert 'cloning' in p.stdout
    rh = home / 'kisscrit-run'
    rec = rh / 'selftest' / ('%s.json' % commit); rec.parent.mkdir(parents=True)
    rec.write_text(json.dumps(dict(commit=commit, passed=True, subset=True, peak_rss_kib=200000,
                                   note='written by tests/test_run_scripts.py, not by run.sh selftest')))
    return dict(base=base, src=src, url=url, commit=commit, home=home, env=env, rh=rh)


RUN_ARGS = ['--shards', '20', '--budget', '20000', '--need-disk-gb', '1', '--poll', '1', '--allow-partial-selftest']


# --------------------------------------------------------------------------------------------- static checks
def test_bash_syntax_and_shellcheck():
    subprocess.run(['bash', '-n', str(RUN)], check=True)
    if shutil.which('shellcheck'):
        p = subprocess.run(['shellcheck', '-x', str(RUN)], capture_output=True, text=True)
        assert p.returncode == 0, p.stdout
    else:
        pytest.skip('shellcheck not installed')


def test_help_lists_every_subcommand():
    out = subprocess.run(['bash', str(RUN), 'help'], capture_output=True, text=True).stdout
    for s in ('setup', 'preflight', 'selftest', 'start', 'status', 'stop', 'check-tokens', 'zenodo-sandbox-test'):
        assert s in out


# --------------------------------------------------------------------------------------------- dry runs
def test_dry_run_setup_installs_nothing_and_uses_wget_without_curl(tmp_path):
    home = tmp_path / 'h'; home.mkdir()
    env = base_env(home); env.pop('KISSCRIT_PYTHON')
    env['KISSCRIT_MINIFORGE'] = str(tmp_path / 'nomf'); env['KISSCRIT_HIDE'] = 'curl'
    p = sh(['setup', '--dry-run', '--skip-preflight', '--repo-url', 'file:///nonexistent'], env)
    assert p.returncode == 0, p.stderr
    assert 'wget -q --tries=3 -O' in p.stdout and 'Miniforge3-Linux' in p.stdout and 'curl -fsSL' not in p.stdout
    assert 'conda create -y -q -n kisscrit-run -c conda-forge --override-channels python=3.12' in p.stdout
    for pkg in ('libflint', 'c-compiler', 'python-flint', 'numba', 'gh', 'git', 'make', 'pytest', 'sympy'):
        assert ' %s' % pkg in p.stdout
    env.pop('KISSCRIT_HIDE')
    p = sh(['setup', '--dry-run', '--skip-preflight', '--repo-url', 'file:///nonexistent'], env)
    assert 'curl -fsSL --retry 3 -o' in p.stdout
    assert not (tmp_path / 'nomf').exists() and not (home / 'kisscrit-run').exists(), 'a dry run changed something'


def test_dry_run_every_subcommand(world, tmp_path):
    env = dict(world['env']); rh = world['rh']; before = sorted(str(p) for p in rh.rglob('*'))
    stub = tmp_path / 'pf.sh'; stub.write_text('#!/bin/bash\necho stub\n')
    env['KISSCRIT_PREFLIGHT'] = str(stub)
    cases = {
        ('preflight', '--dry-run'): 'pf.sh',
        ('selftest', '--commit', world['commit'], '--dry-run'): 'make -C',
        ('start', '--route', 'test', '--commit', world['commit'], '--dry-run', '--survival', 'tmux', '--run-id', 'dry') + tuple(RUN_ARGS): 'worktree add',
        ('start', '--route', 'test', '--commit', world['commit'], '--dry-run', '--survival', 'systemd', '--run-id', 'dry2') + tuple(RUN_ARGS): 'worktree add',
        ('status', '--dry-run'): 'no run yet',
        ('stop', '--dry-run'): None,
        ('check-tokens', '--dry-run'): 'Zenodo tokens',
        ('_run', '--dry-run'): 'not dry-run',
    }
    for args, want in cases.items():
        p = sh(list(args), env)
        if args[0] == 'stop':
            assert 'no current run' in p.stderr; continue
        assert want in p.stdout, (args, p.stdout, p.stderr)
        if args[0] in ('selftest', 'start', 'preflight'):
            assert '[dry-run]' in p.stdout, args
    p = sh(['zenodo-sandbox-test', '--dry-run'], env)
    assert p.returncode != 0 and 'M3' in p.stderr
    after = sorted(str(p) for p in rh.rglob('*'))
    assert before == after, 'a dry run changed files under ~/kisscrit-run: %s' % (set(after) ^ set(before))


def test_systemd_unit_in_dry_run(world):
    p = sh(['start', '--route', 'test', '--commit', world['commit'], '--dry-run', '--survival', 'systemd', '--run-id', 'u'] + RUN_ARGS, world['env'])
    out = p.stdout
    for line in ('ExecStart=/bin/bash', 'scripts/full/run.sh _run', 'Restart=on-failure', 'WantedBy=default.target',
                 'KillMode=control-group', 'StandardOutput=append:'):
        assert line in out, line
    assert '--user daemon-reload' in out and '--user enable kisscrit-run.service' in out and '--user restart kisscrit-run.service' in out


# --------------------------------------------------------------------------------------------- refusals
def test_preflight_failure_stops_with_a_clear_message(world, tmp_path):
    env = dict(world['env']); rep = world['home'] / 'kisscrit-preflight-report.txt'
    stub = tmp_path / 'pf.sh'
    stub.write_text('#!/bin/bash\nprintf "  [ ok ] something\\n  [FAIL] Zenodo (https://zenodo.org/api/records?size=1): no response (exit 6)\\n'
                    '  [FAIL] only 3 GB free where the results go (x; want 50+)\\n" > "$HOME/kisscrit-preflight-report.txt"\n')
    env['KISSCRIT_PREFLIGHT'] = str(stub)
    p = sh(['preflight'], env)
    assert p.returncode == 3
    assert 'cannot be reached' in p.stdout and 'disk is too small' in p.stdout and 'zenodo.org' in p.stdout
    stub.write_text('#!/bin/bash\nprintf "  [ ok ] fine\\n  [warn] tmux missing\\n" > "$HOME/kisscrit-preflight-report.txt"\n')
    p = sh(['preflight'], env)
    assert p.returncode == 0 and 'no failures (1 warnings' in p.stdout and 'open files: soft' in p.stdout
    rep.unlink()


def test_token_permissions_refused_and_never_printed(tmp_path):
    home = tmp_path / 'h'; (home / '.config' / 'kisscrit').mkdir(parents=True)
    env = base_env(home)
    dummy = 'DUMMY-NOT-A-TOKEN-1234567890'
    t = home / '.config' / 'kisscrit' / 'zenodo-sandbox.token'; t.write_text(dummy + '\n'); t.chmod(0o644)
    p = sh(['check-tokens'], env)
    assert p.returncode == 1 and 'REFUSED' in p.stdout and 'mode 644' in p.stdout and 'zenodo.token: missing' in p.stdout
    t.chmod(0o600)
    p = sh(['check-tokens'], env)
    assert p.returncode == 0 and 'mode 600 (ok)' in p.stdout
    t.chmod(0o640)
    assert sh(['check-tokens'], env).returncode == 1
    for q in (p,):
        assert dummy not in q.stdout + q.stderr


def test_commit_must_be_a_tag_or_full_hash(world):
    for ref in ('main', world['commit'][:12], 'HEAD'):
        p = sh(['start', '--route', 'test', '--commit', ref, '--run-id', 'x'] + RUN_ARGS, world['env'])
        assert p.returncode != 0 and 'tag or a full 40-character commit hash' in p.stderr, ref


def test_route_full_refused_until_F6(world):
    p = sh(['start', '--route', 'full', '--commit', world['commit'], '--run-id', 'f'] + RUN_ARGS, world['env'])
    assert p.returncode != 0 and 'F6-READY' in p.stderr


def test_dirty_clone_refused(world):
    f = world['rh'] / 'repo' / 'kisscrit' / '__init__.py'; old = f.read_text()
    try:
        f.write_text(old + '\n# local edit\n')
        p = sh(['start', '--route', 'test', '--commit', world['commit'], '--run-id', 'd'] + RUN_ARGS, world['env'])
        assert p.returncode != 0 and 'uncommitted code changes' in p.stderr
    finally:
        f.write_text(old)


def test_existing_outputs_need_fresh(world):
    """The snapshot holds the cond results; start refuses them without --fresh and leaves nothing behind;
    with --fresh (and --prepare-only) it removes them on the run branch in one commit."""
    env = world['env']; rh = world['rh']
    p = sh(['start', '--route', 'cond', '--commit', world['commit'], '--run-id', 'c1', '--survival', 'tmux'] + RUN_ARGS, env)
    assert p.returncode == 2 and '--fresh' in p.stdout
    assert not (rh / 'runs' / 'c1').exists()
    br = subprocess.run(['git', '-C', str(rh / 'repo'), 'branch', '--list', 'run/*'], capture_output=True, text=True).stdout
    assert 'run/cond-c1' not in br
    p = sh(['start', '--route', 'cond', '--commit', world['commit'], '--run-id', 'c2', '--survival', 'tmux', '--fresh', '--prepare-only'] + RUN_ARGS, env)
    assert p.returncode == 0, p.stdout + p.stderr
    tree = rh / 'runs' / 'c2' / 'tree'
    log = subprocess.run(['git', '-C', str(tree), 'log', '--format=%s', world['commit'] + '..HEAD'], capture_output=True, text=True).stdout
    assert '--fresh: remove' in log and 'run record' in log
    assert not list((tree / 'results' / 'completeness' / 'certificates').glob('certificates-cond-*'))
    names = subprocess.run(['git', '-C', str(tree), 'diff', '--name-only', world['commit'] + '..HEAD'], capture_output=True, text=True).stdout.split()
    assert names and all(n.startswith('results/') for n in names)
    assert (rh / 'runs' / 'c2' / 'logs' / 'fresh-removed.txt').read_text().count('cond') > 100
    rec = json.loads((tree / 'results' / 'full' / 'run-record.json').read_text())
    assert rec['commit'] == world['commit'] and rec['verdict_label'].endswith('; not certified')
    assert rec['checker']['script'] == 'scripts/completeness/check_certificates.py' and rec['checker']['default'] is True
    assert rec['checker']['name'] == 'checker-2-reference' and len(rec['package_list_sha256']) == 64
    # a second new run needs --new while this one is unfinished
    p = sh(['start', '--route', 'test', '--commit', world['commit'], '--run-id', 'c3'] + RUN_ARGS, env)
    assert p.returncode != 0 and '--new' in p.stderr


def test_partial_selftest_needs_the_flag(world):
    args = [a for a in RUN_ARGS if a != '--allow-partial-selftest']
    p = sh(['start', '--new', '--route', 'test', '--commit', world['commit'], '--run-id', 'p1', '--prepare-only'] + args, world['env'])
    assert p.returncode != 0 and 'pytest subset only' in p.stderr


# --------------------------------------------------------------------------------------------- helper units
def test_repair_cuts_partial_lines_and_orphans(tmp_path):
    R = tmp_path / 'results'; c = R / 'completeness'
    for d in ('cert', 'certificates', 'checks', 'gen'): (c / d).mkdir(parents=True)
    line = lambda code, **k: (json.dumps(dict(code=code, **k)) + '\n').encode()
    (c / 'cert' / 'cert-test-0-of-2.jsonl').write_bytes(line('a', status='excluded') + line('b', status='excluded') + b'{"code": "c", "sta')
    (c / 'certificates' / 'certificates-test-0-of-2.jsonl').write_bytes(line('a') + line('b') + line('c'))   # c: no result
    (c / 'cert' / 'cert-test-1-of-2.jsonl').write_bytes(line('d', status='excluded'))                         # d: no certificate
    (c / 'checks' / 'checker-2-reference-certificates-test-0-of-2.jsonl').write_bytes(line('a', verified=True) + b'{"co')
    good = gzip.compress(line('x', faces=[]) + (json.dumps(dict(done_triangulation=3000, new=1)) + '\n').encode())
    torn = gzip.compress(line('y', faces=[]) + line('z', faces=[]))[:-9]
    (c / 'gen' / 'gen-test-0-of-1.jsonl.gz').write_bytes(good + torn)
    log = run_repair.repair(str(R), 'test', str(tmp_path / 'q'))
    reasons = sorted(e['reason'] for e in log)
    assert len(log) == 5, reasons
    codes = lambda p: [json.loads(l)['code'] for l in open(p)]
    assert codes(c / 'cert' / 'cert-test-0-of-2.jsonl') == ['a', 'b']
    assert codes(c / 'certificates' / 'certificates-test-0-of-2.jsonl') == ['a', 'b']
    assert codes(c / 'cert' / 'cert-test-1-of-2.jsonl') == []
    assert codes(c / 'checks' / 'checker-2-reference-certificates-test-0-of-2.jsonl') == ['a']
    with gzip.open(c / 'gen' / 'gen-test-0-of-1.jsonl.gz', 'rt') as f:
        g = [json.loads(l) for l in f]
    assert g[0]['code'] == 'x' and g[1]['done_triangulation'] == 3000 and len(g) in (2, 3, 4)
    assert all(json.loads(l) for l in gzip.open(c / 'gen' / 'gen-test-0-of-1.jsonl.gz', 'rt'))
    q = list((tmp_path / 'q').iterdir()); assert len(q) == 5
    assert run_repair.repair(str(R), 'test', str(tmp_path / 'q')) == [], 'repair is not idempotent'


def test_workers_cap():
    GiB = 1 << 20
    assert run_state.workers(64, 256 * GiB, 500 * 1024, 2, 0, 40) == (64, 21, 2, 209)
    assert run_state.workers(64, 16 * GiB, 1 * GiB, 2, 0, 40)[0] == 6            # memory caps the workers
    assert run_state.workers(64, 256 * GiB, 500 * 1024, 2, 8, 40)[:3] == (8, 2, 2)
    assert run_state.workers(2, 8 * GiB, 200 * 1024, 2, 0, 2)[:3] == (2, 1, 2)
    assert run_state.workers(64, 256 * GiB, 500 * 1024, 2, 0, 3)[1] == 3          # never more shards than exist


def test_jset_jget_roundtrip(tmp_path):
    f = tmp_path / 'x.json'
    run_state.dump(str(f), {})
    subprocess.run([sys.executable, str(HELPER), 'jset', str(f), 'a.b', '3'], check=True)
    subprocess.run([sys.executable, str(HELPER), 'jappend', str(f), 'l', '{"r": "x y"}'], check=True)
    out = subprocess.run([sys.executable, str(HELPER), 'jget', str(f), 'a.b'], capture_output=True, text=True).stdout
    assert out.strip() == '3' and json.loads(f.read_text())['l'] == [{'r': 'x y'}]


# --------------------------------------------------------------------------------------------- real runs
def _wait(pred, timeout, what):
    t0 = time.time()
    while not pred():
        assert time.time() - t0 < timeout, 'timed out waiting for ' + what
        time.sleep(0.5)


def test_kill_and_resume_route_test(world):
    """A real route-test run (shards 0 and 1 of 20), killed with SIGKILL to its whole process group in the
    middle of certify, with a torn line appended as a kill during a write would leave it; `start` resumes it to
    the end.  Each shard then has every candidate once, in the result file, the certificate file and the check
    file, and the run branch touches only results/ (never the certificates)."""
    env = dict(world['env']); rh = world['rh']; C = world['commit']
    p = sh(['start', '--new', '--route', 'test', '--commit', C, '--run-id', 'k1', '--only-shards', '0,1',
            '--survival', 'foreground', '--prepare-only'] + RUN_ARGS, env)
    assert p.returncode == 0, p.stdout + p.stderr
    tree = rh / 'runs' / 'k1' / 'tree'; cdir = tree / 'results' / 'completeness'
    res0 = cdir / 'cert' / 'cert-test-0-of-20.jsonl'
    proc = subprocess.Popen(['bash', str(RUN), 'start'], env=env, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                            start_new_session=True)
    try:
        _wait(lambda: res0.exists() and len(res0.read_bytes().splitlines()) >= 4, 600, 'four result lines in shard 0')
    finally:
        os.killpg(proc.pid, signal.SIGKILL); proc.wait()
    def live():   # members of the killed group that are not zombies (the container's init reaps them late)
        out = subprocess.run(['ps', '-o', 'stat=', '-g', str(proc.pid)], capture_output=True, text=True).stdout
        return [s for s in out.split() if not s.startswith('Z')]
    _wait(lambda: not live(), 30, 'the killed processes to exit')
    st = sh(['status'], env).stdout
    assert 'processes: none' in st and 'lock: free' in st and 'phase: certify' in st
    with open(res0, 'ab') as f: f.write(b'{"code": "torn by the test", "stat')       # a write cut by the kill
    p = sh(['start'], env, timeout=1800)
    assert p.returncode == 0, p.stdout[-3000:] + p.stderr[-3000:]
    assert (rh / 'runs' / 'k1' / 'state' / 'done').exists()
    rep = (rh / 'runs' / 'k1' / 'logs' / 'repair.jsonl').read_text()
    assert 'partial or malformed tail' in rep
    for k in (0, 1):
        codes = lambda p: [json.loads(l)['code'] for l in open(p)]
        r = codes(cdir / 'cert' / ('cert-test-%d-of-20.jsonl' % k))
        c = codes(cdir / 'certificates' / ('certificates-test-%d-of-20.jsonl' % k))
        ch = [json.loads(l) for l in open(cdir / 'checks' / ('checker-2-reference-certificates-test-%d-of-20.jsonl' % k))]
        assert len(r) == len(set(r)) > 0 and sorted(r) == sorted(c) == sorted(x['code'] for x in ch)
        assert all(x['verified'] for x in ch)
    man = json.loads((tree / 'results' / 'full' / 'manifest.json').read_text())
    assert man['cert']['duplicates'] == man['certificates']['duplicates'] == man['checks']['duplicates'] == 0
    assert man['complete'] is False                       # 2 of 20 shards: never complete
    rec = json.loads((tree / 'results' / 'full' / 'run-record.json').read_text())
    assert len(rec['resumes']) == 2 and rec['phases']['check']['wall_seconds'] >= 0
    g = lambda *a: subprocess.run(['git', '-C', str(tree)] + list(a), capture_output=True, text=True).stdout
    names = g('diff', '--name-only', C + '..HEAD').split()
    assert names and all(n.startswith('results/') for n in names)
    assert not [n for n in names if n.startswith('results/completeness/certificates')]
    subjects = g('log', '--format=%s', C + '..HEAD')
    assert 'certify test shard 0/20' in subjects and 'certify test shard 1/20' in subjects and 'check test shard 1/20' in subjects
    assert not g('status', '--porcelain', '--', '.', ':(exclude)results')
    assert "not certified" in (tree / 'results' / 'full' / 'report.txt').read_text()
    assert (tree / 'results' / 'full' / 'report-interim.txt').exists()
    p = sh(['status', '--final'], env)
    assert 'manifest: complete=false' in p.stdout and 'not certified' in p.stdout


@pytest.mark.skipif(not shutil.which('tmux'), reason='tmux not installed')
def test_tmux_start_stop_and_reboot_line(world):
    env = dict(world['env']); rh = world['rh']; C = world['commit']
    p = sh(['start', '--new', '--route', 'test', '--commit', C, '--run-id', 't1', '--only-shards', '2',
            '--survival', 'tmux'] + RUN_ARGS, env)
    assert p.returncode == 0, p.stdout + p.stderr
    cron = (world['home'] / '.fake-crontab').read_text()
    assert cron.count('kisscrit-run @reboot line') == 1 and 'runs/t1/tree/scripts/full/run.sh start --boot' in cron
    _wait(lambda: (rh / 'runs' / 't1' / 'state' / 'pid').exists(), 60, 'the pid file')
    st = sh(['status'], env).stdout
    assert 'tmux session kisscrit-run' in st and ': running' in st and '@reboot crontab line: present' in st
    p = sh(['stop', '--reason', 'test stop'], env)
    assert p.returncode == 0 and 'test stop' in p.stdout
    st = sh(['status'], env).stdout
    assert 'processes: none' in st and 'stopped:' in st and 'test stop' in st and ': absent' in st
    rec = json.loads((rh / 'runs' / 't1' / 'tree' / 'results' / 'full' / 'run-record.json').read_text())
    assert rec['stops'][-1]['reason'] == 'test stop'
    # after a reboot the @reboot line runs `start --boot`: a stopped run stays stopped
    p = sh([str(rh / 'runs' / 't1' / 'tree' / 'scripts' / 'full' / 'run.sh'), 'start', '--boot'], env)
    assert p.returncode == 0 and 'not restarted' in p.stdout
    assert subprocess.run(['tmux', '-L', env['KISSCRIT_TMUX_SOCKET'], 'has-session', '-t', 'kisscrit-run'],
                          capture_output=True).returncode != 0
    # a second start in tmux replaces the tagged crontab line instead of adding one
    p = sh(['start', '--new', '--route', 'test', '--commit', C, '--run-id', 't2', '--only-shards', '3',
            '--survival', 'tmux'] + RUN_ARGS, env)
    assert p.returncode == 0
    cron = (world['home'] / '.fake-crontab').read_text()
    assert cron.count('kisscrit-run @reboot line') == 1 and 'runs/t2/' in cron
    # one run at a time: while it holds the lock, start says so and does not start another copy
    _wait(lambda: 'lock: held' in sh(['status'], env).stdout, 60, 'the lock')
    p = sh(['start'], env)
    assert 'already running' in p.stdout
    assert sh(['stop', '--reason', 'end of test'], env).returncode == 0


def test_systemd_launch_with_fake_systemctl(world):
    env = dict(world['env']); C = world['commit']
    log = world['home'] / '.fake-systemctl.log'
    if log.exists(): log.unlink()
    p = sh(['start', '--new', '--route', 'test', '--commit', C, '--run-id', 's1', '--only-shards', '4',
            '--survival', 'systemd'] + RUN_ARGS, env)
    assert p.returncode == 0, p.stdout + p.stderr
    unit = world['home'] / '.config' / 'systemd' / 'user' / 'kisscrit-run.service'
    text = unit.read_text()
    assert 'runs/s1/tree/scripts/full/run.sh _run' in text and 'Restart=on-failure' in text
    calls = log.read_text().splitlines()
    assert calls[:3] == ['--user daemon-reload', '--user enable kisscrit-run.service', '--user restart kisscrit-run.service']
    if shutil.which('systemd-analyze'):
        v = subprocess.run(['systemd-analyze', 'verify', '--man=no', str(unit)], capture_output=True, text=True)
        bad = [l for l in (v.stdout + v.stderr).splitlines() if 'kisscrit-run.service' in l and 'Failed to' not in l
               and 'bus' not in l.lower()]
        assert not bad, v.stdout + v.stderr
    st = sh(['status'], env).stdout
    assert 'systemd user service kisscrit-run.service' in st
    p = sh(['stop', '--reason', 'systemd test'], env)
    assert p.returncode == 0
    assert '--user stop kisscrit-run.service' in log.read_text() and '--user disable kisscrit-run.service' in log.read_text()


def test_real_selftest_with_a_pytest_subset(world):
    """The self-test for real: C build and make test, a pytest subset, selftest.py, peak RSS (about 2 minutes)."""
    if not (pathlib.Path('/usr/include/flint/arb.h').exists() or os.environ.get('CPATH')):
        pytest.skip('FLINT headers not found')
    env = dict(world['env']); C = world['commit']
    p = sh(['selftest', '--commit', C, '--pytest-subset', 'tests/test_code_standards.py'], env, timeout=1800)
    assert p.returncode == 0, p.stdout + p.stderr
    rec = json.loads((world['rh'] / 'selftest' / ('%s.json' % C)).read_text())
    assert rec['passed'] is True and rec['subset'] is True and rec['peak_rss_kib'] > 10000
    assert rec['parts']['c_kernel'] == 'passed' and rec['parts']['selftest_py'] == 'passed'
    assert rec['parts']['pytest']['failed'] == 0 and rec['parts']['pytest']['passed'] > 0


def test_a_lone_copy_of_run_sh_works(world, tmp_path):
    """docs/setup-guide.md B1 copies run.sh alone to the server; status and check-tokens work from that copy, and
    status hands over to the run tree's own copy."""
    lone = tmp_path / 'run.sh'; shutil.copy(RUN, lone)
    p = sh([str(lone), 'status'], world['env'])
    assert p.returncode == 0 and ('run ' in p.stdout or 'no run yet' in p.stdout), p.stdout + p.stderr
    assert sh([str(lone), 'check-tokens'], world['env']).returncode == 0
