#!/usr/bin/env python3
"""Repair before resume (scripts/full/run.sh): make the append-only result files of a route safe to resume
after a kill, moving everything removed to a quarantine directory.

Status: untrusted (run tooling, checklist F8; no part of the proof).  It never adds or changes a result: it only
removes damaged or unpaired lines, which the next run of run_shard.py or check_certificates.py then redoes.

  python3 scripts/full/run_repair.py RESULTS_DIR ROUTE QUARANTINE_DIR [--log FILE]

The underlying cause, writes that a kill can cut, is in scripts/completeness/run_shard.py and
check_certificates.py (read-only for F8); a version of run_shard.py with atomic writes is proposed, not written.
"""
import argparse, glob, gzip, json, os, time, zlib


def atomic_write(path, data, mode='w'):
    """Write data to path through a temporary file and a rename (POSIX rename is atomic in one filesystem)."""
    d = os.path.dirname(os.path.abspath(path)); os.makedirs(d, exist_ok=True)
    tmp = os.path.join(d, '.%s.tmp.%d' % (os.path.basename(path), os.getpid()))
    with open(tmp, mode) as f:
        f.write(data); f.flush(); os.fsync(f.fileno())
    os.replace(tmp, path)


def route_files(results, route):
    """The result files of a route under results/completeness, by kind (the names run_shard.py and
    check_certificates.py write), and the run files under results/full that a new run would overwrite."""
    R = os.path.join(results, 'completeness')
    g = lambda base, pat: sorted(glob.glob(os.path.join(base, pat)))
    F = os.path.join(results, 'full')
    return dict(gen=g(R, 'gen/gen-%s-*.jsonl.gz' % route), cert=g(R, 'cert/cert-%s-*.jsonl' % route),
                done=g(R, 'cert/done-%s-*.txt' % route), certificates=g(R, 'certificates/certificates-%s-*.jsonl' % route),
                checks=g(R, 'checks/*-certificates-%s-*.jsonl' % route), claims=g(R, 'claims/*-%s-*.txt' % route),
                summary=g(R, 'summary-%s.json' % route),
                run=[p for p in g(F, '*') if os.path.basename(p) in (
                    'progress.json', 'run-record.json', 'manifest.json', 'check-summary.json', 'report.txt',
                    'report-interim.txt', 'conda-list.txt', 'zenodo-state.json', 'zenodo-manifest.json')])


def _split_jsonl(path):
    """(complete lines, bad tail) of a JSON-lines file: a line counts only if it ends in a newline and parses;
    everything from the first bad line on is the tail (a kill can only damage the end of an append-only file)."""
    with open(path, 'rb') as f: data = f.read()
    good, pos = [], 0
    while pos < len(data):
        nl = data.find(b'\n', pos)
        if nl < 0: break
        line = data[pos:nl + 1]
        try: json.loads(line)
        except ValueError: break
        good.append(line); pos = nl + 1
    return good, data[pos:]


def _read_gz_tolerant(path):
    """(complete lines, damaged) of an appended gzip file: members are read one by one and reading stops at the
    first truncated or corrupt member; a line without its newline at the end is dropped."""
    with open(path, 'rb') as f: data = f.read()
    out, damaged = b'', False
    while data:
        d = zlib.decompressobj(zlib.MAX_WBITS | 16)
        try:
            chunk = d.decompress(data)
        except zlib.error:
            damaged = True; break
        if not d.eof:                       # member cut short by a kill: keep what it decoded so far
            out += chunk; damaged = True; break
        out += chunk; data = d.unused_data
    lines = out.split(b'\n')
    tail = lines.pop()                      # text after the last newline: empty for an intact file
    good = []
    for i, l in enumerate(lines):
        try: json.loads(l); good.append(l + b'\n')
        except ValueError: damaged = True; break
    if tail: damaged = True
    return good, damaged


def _quarantine(path, qdir, payload, why, log):
    """Save payload (bytes) under qdir with the reason, and log it."""
    os.makedirs(qdir, exist_ok=True)
    q = os.path.join(qdir, '%s.%s.%d' % (os.path.basename(path), time.strftime('%Y%m%dT%H%M%S'), os.getpid()))
    with open(q, 'wb') as f: f.write(payload)
    log.append(dict(file=path, quarantined=q, bytes=len(payload), reason=why))


def repair(results, route, qdir):
    """Make the append-only result files of a route consistent before a resume, moving anything removed to qdir.

    (1) A partial last line (no newline, or not JSON) of a .jsonl file is cut off: run_shard.py and
        check_certificates.py read their own output back with json.loads and would stop on it.
    (2) A damaged gzip member at the end of a gen file is dropped and the complete lines are rewritten as one
        member; the triangulation whose `done_triangulation` line is missing is redone by run_shard.py.
    (3) run_shard.py writes the certificate line first and the result line second.  A kill between the two
        leaves a certificate with no result; on resume the candidate is redone and a second certificate is
        appended, a duplicate in the manifest (V5).  So a certificate line whose code has no result line in the
        matching cert file is removed, and a result line without a certificate is removed too (redone).
    Returns the list of actions (empty when nothing was damaged)."""
    log = []
    F = route_files(results, route)
    for p in F['cert'] + F['certificates'] + F['checks']:
        good, tail = _split_jsonl(p)
        if tail:
            _quarantine(p, qdir, tail, 'partial or malformed tail of a JSON-lines file', log)
            atomic_write(p, b''.join(good), 'wb')
    for p in F['gen']:
        good, damaged = _read_gz_tolerant(p)
        if damaged:
            with open(p, 'rb') as f: _quarantine(p, qdir, f.read(), 'damaged gzip member (original file)', log)
            atomic_write(p, gzip.compress(b''.join(good)), 'wb')
    cdir = os.path.join(results, 'completeness')
    tags = sorted({os.path.basename(f)[len('certificates-'):] for f in F['certificates']} |
                  {os.path.basename(f)[len('cert-'):] for f in F['cert']})        # <route>-<k>-of-<K>.jsonl
    for tag in tags:
        cf = os.path.join(cdir, 'certificates', 'certificates-' + tag)
        rf = os.path.join(cdir, 'cert', 'cert-' + tag)
        res = [l for l in open(rf, 'rb')] if os.path.exists(rf) else []
        cer = [l for l in open(cf, 'rb')] if os.path.exists(cf) else []
        rcodes = {json.loads(l)['code'] for l in res}
        ccodes = {json.loads(l)['code'] for l in cer}
        orphan_c = [l for l in cer if json.loads(l)['code'] not in rcodes]
        orphan_r = [l for l in res if json.loads(l)['code'] not in ccodes]
        if orphan_c:
            _quarantine(cf, qdir, b''.join(orphan_c), 'certificate without a result line (kill between the writes)', log)
            atomic_write(cf, b''.join(l for l in cer if json.loads(l)['code'] in rcodes), 'wb')
        if orphan_r:
            _quarantine(rf, qdir, b''.join(orphan_r), 'result line without a certificate', log)
            atomic_write(rf, b''.join(l for l in res if json.loads(l)['code'] in ccodes), 'wb')
    return log


def main():
    """Command line: repair one route and print one line per action."""
    ap = argparse.ArgumentParser(description='repair result files before a resume')
    ap.add_argument('results'); ap.add_argument('route'); ap.add_argument('quarantine'); ap.add_argument('--log')
    a = ap.parse_args()
    log = repair(a.results, a.route, a.quarantine)
    for e in log:
        print('repair: %s: %s (%d bytes kept in %s)' % (e['file'], e['reason'], e['bytes'], e['quarantined']))
    if a.log and log:
        with open(a.log, 'a') as f:
            for e in log: f.write(json.dumps(dict(e, time=time.strftime('%Y-%m-%dT%H:%M:%S%z'))) + '\n')
    if not log: print('repair: nothing to repair')


if __name__ == '__main__':
    main()
