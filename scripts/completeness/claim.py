#!/usr/bin/env python3
"""Claim the next free shard through git, so several machines can share the work without a server.

Status: untrusted (exploration or tooling; no part of the proof).

A claim is a small file results/completeness/claims/<phase>-<route>-<k>-of-<K>.txt naming the machine.  Claiming
commits and pushes that file; if the push is rejected (someone else pushed first), we pull, and if the shard
was taken in the meantime we try the next one.

  python scripts/completeness/claim.py --phase certify --route cond --of 40 --name <machine-name>
prints the claimed shard number (and nothing else on stdout); exits with status 1 if none are free, 3 if git
fails.
"""
import argparse, os, subprocess, sys, datetime

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
CL = os.path.join(ROOT, 'results', 'completeness', 'claims')


def git(*a, check=True):
    return subprocess.run(['git', '-C', ROOT] + list(a), check=check, capture_output=True, text=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--phase', required=True, choices=['gen', 'certify'])
    ap.add_argument('--route', default='cond')
    ap.add_argument('--of', type=int, required=True)
    ap.add_argument('--name', required=True)
    a = ap.parse_args()
    os.makedirs(CL, exist_ok=True)
    for attempt in range(20):
        r = git('pull', '--rebase', '--quiet', check=False)
        if r.returncode != 0:
            print('git pull failed (network? uncommitted changes?):\n' + r.stderr, file=sys.stderr); sys.exit(3)
        taken = set(os.listdir(CL))
        free = [k for k in range(a.of) if '%s-%s-%d-of-%d.txt' % (a.phase, a.route, k, a.of) not in taken]
        if not free:
            print('no free shards', file=sys.stderr); sys.exit(1)
        k = free[0]; fn = os.path.join(CL, '%s-%s-%d-of-%d.txt' % (a.phase, a.route, k, a.of))
        open(fn, 'w').write('%s %s\n' % (a.name, datetime.datetime.now().isoformat(timespec='seconds')))
        git('add', fn); git('commit', '--quiet', '-m', 'claim %s %s shard %d/%d (%s)' % (a.phase, a.route, k, a.of, a.name))
        if git('push', '--quiet', check=False).returncode == 0:
            print(k); return
        git('reset', '--hard', '--quiet', 'HEAD~1')     # lost the race: drop our claim commit and retry
    print('could not claim a shard (push kept failing)', file=sys.stderr); sys.exit(1)


if __name__ == '__main__':
    main()
