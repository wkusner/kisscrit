#!/usr/bin/env python3
"""Helpers for scripts/full/run.sh: run record, repair before resume, progress, worker cap, peak RSS, hashes.

Status: untrusted (run tooling, checklist F8; no part of the proof).  Nothing here decides whether a candidate is
excluded: the search writes the certificates and the checker named in the run record decides.

run.sh calls this file with one subcommand per call:
  env CONFIG                  print the run configuration as shell assignments (for eval)
  jget FILE KEY               print one value of a JSON file (dotted key; '' if absent)
  jset FILE KEY JSON          set a value (dotted key) in a JSON file, written atomically
  jappend FILE KEY JSON       append a value to a list (dotted key) in a JSON file
  sha256 FILE ...             SHA-256 of files, with Python hashlib (Ubuntu 26.04's sha256sum is a slow rewrite)
  outputs TREE ROUTE          list result files of a route that already exist (start refuses them without --fresh)
  workers ...                 the number of worker processes and of concurrent shards (see workers())
  rss -- CMD ...              run CMD and print the peak resident set size of its largest descendant, in KiB
  progress ...                write results/full/progress.json (rates only, no time estimate)
  versions                    print the versions of Python and the packages the run imports, as JSON
  sizes TREE ROUTE            bytes of the route's result files by kind

The repair before a resume is scripts/full/run_repair.py.  All writes go through atomic_write (temporary file in the same directory, fsync, rename), so a kill leaves either
the old file or the new one.
"""
import argparse, collections, glob, hashlib, json, os, resource, shlex, subprocess, sys, time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from run_repair import atomic_write, route_files, _split_jsonl  # noqa: E402  (the same directory)


# ------------------------------------------------------------------------------------------------ JSON files
def load(path):
    """The JSON object in path, or {} if the file does not exist."""
    if not os.path.exists(path): return {}
    with open(path) as f: return json.load(f)


def dump(path, obj):
    """Write obj as indented JSON, atomically."""
    atomic_write(path, json.dumps(obj, indent=1, sort_keys=True) + '\n')


def _walk(obj, key, create):
    """The parent dict and last key of a dotted key."""
    parts = key.split('.')
    for p in parts[:-1]:
        if p not in obj or not isinstance(obj[p], dict):
            if not create: return None, None
            obj[p] = {}
        obj = obj[p]
    return obj, parts[-1]


def cmd_jget(a):
    """Print a value: strings as they are, other values as JSON, nothing if absent."""
    parent, k = _walk(load(a.file), a.key, False)
    if parent is None or k not in parent: return
    v = parent[k]
    print(v if isinstance(v, str) else json.dumps(v))


def cmd_jset(a):
    """Set a dotted key to a JSON value (a bare word that is not JSON is stored as a string)."""
    obj = load(a.file); parent, k = _walk(obj, a.key, True)
    try: parent[k] = json.loads(a.value)
    except ValueError: parent[k] = a.value
    dump(a.file, obj)


def cmd_jappend(a):
    """Append a JSON value to the list at a dotted key (created if absent)."""
    obj = load(a.file); parent, k = _walk(obj, a.key, True)
    try: v = json.loads(a.value)
    except ValueError: v = a.value
    parent.setdefault(k, []).append(v)
    dump(a.file, obj)


ENV_KEYS = ['run_id', 'route', 'commit', 'ref', 'source', 'branch', 'tree', 'shards', 'gen_shards', 'budget',
            'workers', 'wps', 'checker_cmd', 'checker_default', 'only_shards', 'survival', 'python', 'poll',
            'need_disk_gb', 'progress_every', 'push']


def cmd_env(a):
    """Shell assignments KC_<KEY>=<quoted value> for the run configuration."""
    c = load(a.config)
    for k in ENV_KEYS:
        v = c.get(k, '')
        if isinstance(v, bool): v = int(v)
        if isinstance(v, list): v = ' '.join(map(str, v))
        print('KC_%s=%s' % (k.upper(), shlex.quote(str(v))))


def sha256_file(path, block=1 << 20):
    """SHA-256 of a file, read in 1 MiB blocks."""
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(block), b''): h.update(b)
    return h.hexdigest()


def cmd_sha256(a):
    """Print '<hex>  <path>' per file, like sha256sum."""
    for p in a.files: print('%s  %s' % (sha256_file(p), p))


# ------------------------------------------------------------------------------------------------ result files
def cmd_outputs(a):
    """Print every existing result file of the route, one per line, relative to the tree."""
    for kind, files in route_files(os.path.join(a.tree, 'results'), a.route).items():
        for f in files: print(os.path.relpath(f, a.tree))


# ------------------------------------------------------------------------------------------------ resources
def workers(nproc, mem_avail_kib, rss_kib, wps, requested, shards):
    """(total workers W, concurrent shards P, workers per shard w).

    W = min(requested or nproc, floor(0.8 * MemAvailable / per-process memory)); per-process memory is twice
    the measured peak RSS of one certify worker (safety factor 2, at least 256 MiB), because the self-test
    measures it on route test, whose candidates are lighter than the heavy ones of cond and full.  Each shard
    runs one driver process plus w pool workers, so P = max(1, floor(W / (w + 1))), capped by the shard count."""
    per = max(2 * rss_kib, 256 * 1024)
    by_mem = max(1, int(0.8 * mem_avail_kib // per))
    W = max(1, min(requested or nproc, by_mem))
    w = max(1, min(wps, W))
    P = max(1, min(shards, W // (w + 1) if W >= w + 1 else 1))
    return W, P, w, by_mem


def cmd_workers(a):
    """Print 'W P w by_mem' for run.sh."""
    print(*workers(a.nproc, a.mem_avail_kib, a.rss_kib, a.wps, a.requested, a.shards))


def cmd_rss(a):
    """Run the command, then print the peak RSS (KiB) of its largest waited-for descendant (Linux: ru_maxrss of
    RUSAGE_CHILDREN is in KiB) and the command's exit status."""
    t0 = time.time()
    rc = subprocess.call(a.cmd[1:] if a.cmd and a.cmd[0] == '--' else a.cmd)
    kib = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    if sys.platform == 'darwin': kib //= 1024           # bytes on macOS
    print(json.dumps(dict(peak_rss_kib=kib, exit=rc, seconds=round(time.time() - t0, 1))))
    sys.exit(rc)


def vmhwm_kib(pattern, cwd):
    """The largest VmHWM (peak RSS, KiB) over the running processes whose command line contains pattern and
    whose working directory is cwd (run.sh starts run_shard.py from the run tree with a relative path)."""
    best = 0
    for d in glob.glob('/proc/[0-9]*'):
        try:
            cmd = open(os.path.join(d, 'cmdline'), 'rb').read().replace(b'\0', b' ').decode(errors='replace')
            if pattern not in cmd or os.path.realpath(os.readlink(os.path.join(d, 'cwd'))) != os.path.realpath(cwd): continue
            for line in open(os.path.join(d, 'status')):
                if line.startswith('VmHWM:'): best = max(best, int(line.split()[1]))
        except OSError:
            continue
    return best


# ------------------------------------------------------------------------------------------------ progress
def cmd_progress(a):
    """Write results/full/progress.json: phase, shards, candidates by status, boxes, CPU seconds and measured
    rates.  No time estimate (decision of 2026-10-15)."""
    cfg = load(a.config); st = load(a.state)
    results = os.path.join(cfg['tree'], 'results'); route = cfg['route']; K = int(cfg['shards'])
    shards = [int(s) for s in cfg.get('only_shards') or []] or list(range(K))
    cdir = os.path.join(results, 'completeness')
    by_status = collections.Counter(); boxes = 0; cpu = 0.0; n = 0
    for k in shards:
        p = os.path.join(cdir, 'cert', 'cert-%s-%d-of-%d.jsonl' % (route, k, K))
        if not os.path.exists(p): continue
        good, _ = _split_jsonl(p)
        for l in good:
            d = json.loads(l); by_status[d['status']] += 1; n += 1
            boxes += d.get('stats', {}).get('boxes', 0); cpu += d.get('seconds', 0.0)
    done = [k for k in shards if os.path.exists(os.path.join(cdir, 'cert', 'done-%s-%d-of-%d.txt' % (route, k, K)))]
    checked = verified = 0
    for f in glob.glob(os.path.join(cdir, 'checks', '*-certificates-%s-*-of-%d.jsonl' % (route, K))):
        good, _ = _split_jsonl(f)
        for l in good:
            checked += 1; verified += bool(json.loads(l).get('verified'))
    now = time.time()
    t_cert = st.get('phases', {}).get('certify', {}).get('first_start_epoch')
    wall_h = (now - t_cert) / 3600 if t_cert else None
    out = dict(run_id=cfg['run_id'], route=route, commit=cfg['commit'], branch=cfg.get('branch'),
               updated_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime(now)),
               phase=st.get('phase', 'not started'), shards_total=len(shards), shards_done=len(done),
               candidates_certified=n, by_status=dict(by_status), boxes=boxes, search_cpu_seconds=round(cpu, 1),
               checked=checked, verified_by_checker=verified,
               rates=dict(candidates_per_wall_hour=round(n / wall_h, 1) if wall_h and wall_h > 1e-3 else None,
                          search_cpu_seconds_per_candidate=round(cpu / n, 3) if n else None,
                          boxes_per_cpu_second=round(boxes / cpu, 1) if cpu else None),
               peak_rss_kib_seen=max(st.get('peak_rss_kib_seen', 0), vmhwm_kib('scripts/completeness/run_shard.py', cfg['tree'])),
               workers=cfg.get('workers'), concurrent_shards=cfg.get('concurrent_shards'), wps=cfg.get('wps'),
               verdict=cfg.get('verdict_label'),
               note='rates only; no time estimate (docs/decisions.md, 2026-10-15)')
    dump(a.out, out)
    if out['peak_rss_kib_seen'] > st.get('peak_rss_kib_seen', 0):
        st['peak_rss_kib_seen'] = out['peak_rss_kib_seen']; dump(a.state, st)
    print('progress: phase %s, %d of %d shards, %d candidates %s, %d checked (%d verified by the checker)' % (
        out['phase'], len(done), len(shards), n, dict(by_status), checked, verified))


def cmd_sizes(a):
    """Bytes of the route's result files by kind, as one line 'gen 1.5 MB; cert ...' (only this route's files)."""
    parts = []
    for kind, files in route_files(os.path.join(a.tree, 'results'), a.route).items():
        if kind in ('done', 'claims', 'summary', 'run'): continue
        n = sum(os.path.getsize(f) for f in files)
        parts.append('%s %.1f MB in %d files' % (kind, n / 1e6, len(files)))
    print('; '.join(parts))


def cmd_versions(a):
    """Versions of Python and the packages the run imports."""
    v = dict(python=sys.version.split()[0], executable=sys.executable)
    for mod, name in (('numpy', 'numpy'), ('scipy', 'scipy'), ('numba', 'numba'), ('flint', 'python-flint'),
                      ('networkx', 'networkx'), ('mpmath', 'mpmath'), ('sympy', 'sympy'), ('pytest', 'pytest')):
        try: v[name] = __import__(mod).__version__
        except Exception as e: v[name] = 'import failed (%s)' % type(e).__name__
    print(json.dumps(v))


def main():
    """Command line entry point (one subcommand per call)."""
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0]); sp = ap.add_subparsers(dest='cmd', required=True)
    p = sp.add_parser('env'); p.add_argument('config'); p.set_defaults(f=cmd_env)
    p = sp.add_parser('jget'); p.add_argument('file'); p.add_argument('key'); p.set_defaults(f=cmd_jget)
    p = sp.add_parser('jset'); p.add_argument('file'); p.add_argument('key'); p.add_argument('value'); p.set_defaults(f=cmd_jset)
    p = sp.add_parser('jappend'); p.add_argument('file'); p.add_argument('key'); p.add_argument('value'); p.set_defaults(f=cmd_jappend)
    p = sp.add_parser('sha256'); p.add_argument('files', nargs='+'); p.set_defaults(f=cmd_sha256)
    p = sp.add_parser('outputs'); p.add_argument('tree'); p.add_argument('route'); p.set_defaults(f=cmd_outputs)
    p = sp.add_parser('workers')
    for k in ('nproc', 'mem_avail_kib', 'rss_kib', 'wps', 'requested', 'shards'):
        p.add_argument('--' + k.replace('_', '-'), dest=k, type=int, required=True)
    p.set_defaults(f=cmd_workers)
    p = sp.add_parser('rss'); p.add_argument('cmd', nargs=argparse.REMAINDER); p.set_defaults(f=cmd_rss)
    p = sp.add_parser('progress'); p.add_argument('--config', required=True); p.add_argument('--state', required=True)
    p.add_argument('--out', required=True); p.set_defaults(f=cmd_progress)
    p = sp.add_parser('versions'); p.set_defaults(f=cmd_versions)
    p = sp.add_parser('sizes'); p.add_argument('tree'); p.add_argument('route'); p.set_defaults(f=cmd_sizes)
    a = ap.parse_args(); a.f(a)


if __name__ == '__main__':
    main()
