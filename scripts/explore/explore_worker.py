"""One worker of the overnight exploration: random ascent and contact-graph surgery, then a census file.

Status: untrusted (exploration or tooling; no part of the proof).

usage: python scripts/explore/explore_worker.py N d --out FILE [--seed S] [--index I] [--runs R] [--rounds K]
           [--budget SECONDS] [--near M] [--dump-seconds T] [--start MERGED --slice-index i --slice-count W
           --surgery-types M] [--commit HASH] [--run-id ID] [--batch B]

The worker calls kisscrit.search.random_ascent and kisscrit.search.surgery unchanged and writes a census file in
the format of Census.dump (N, d, types) plus a "worker" entry that Census.load ignores.  Everything it finds is
numerical (double precision) and unverified.

* Reproducible: the random generator is numpy.random.default_rng(SeedSequence([seed, index])), so a worker with
  the same seed, index, runs, rounds and start file and enough time finds the same census (the ascents are
  deterministic given the generator; numerics are single-threaded).  With a short budget the census depends on
  the machine's speed, because the budget ends the work wherever it is.
* Order of work: without --start, R ascents and then K surgery rounds on the worker's own finds (the order of
  scripts/ascent_census.py).  With --start (a merged census), the worker first runs surgery on its slice of the
  types that no surgery has used yet (position i of every W in the order of decreasing r, at most M of them),
  then spends the rest of the budget on ascents.  A type whose surgery call completed carries "surg" = the worker
  index; a call cut short by the budget marks nothing, so those types come back in a later batch.
* Output: dumped at least every --dump-seconds, on SIGTERM and when the budget ends; each dump is written to a
  temporary name and renamed, so the file is always a complete census.  The "worker" entry says why the worker
  stopped: finished, budget (a time limit cut some phase) or sigterm.  "hits" counts the worker's own finds
  (types taken from --start begin at 0).  Each new type gets "prov": commit, seed, worker index, batch, run id
  and the UTC time at which the worker first saw it (to the resolution of one ascent or one surgery call).
"""
import os
# One thread per process: the driver runs one worker per core, and a threaded BLAS or OpenMP would oversubscribe.
# These must be set before numpy is imported, which explore_common does.
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ[_v] = '1'
import argparse, copy, signal, sys, time
import numpy as np
from explore_common import atomic_write_json, load_census_file, utc_now
from kisscrit.search import Census, random_ascent, surgery


class _Stop(Exception):
    """Raised from a signal handler to leave the current phase; the census is dumped by the caller."""


class _Worker:
    """State of one run: the census, the clock, and what has been done so far."""

    def __init__(self, N, d, out, seed, index, dump_seconds, commit, run_id, batch, slice_info):
        self.C = Census(N, d); self.out = out; self.seed, self.index = seed, index
        self.dump_seconds, self.commit, self.run_id, self.batch = dump_seconds, commit, run_id, batch
        self.slice_info = slice_info
        self.t0 = time.time(); self.started_utc = utc_now(); self.last_dump = self.t0
        self.limit = self.t0; self.terminated = False; self.dumping = False
        self.ascents = 0; self.surgery_done = False; self.stopped_by = 'finished'

    def stamp_new(self):
        """Give every type that has no provenance yet the worker's provenance and the current UTC time."""
        for t in self.C.types.values():
            if 'prov' not in t:
                t['prov'] = dict(commit=self.commit, seed=self.seed, worker=self.index, batch=self.batch,
                                 run_id=self.run_id, first_found_utc=utc_now())

    def dump(self):
        """Write the census atomically.  A dump already in progress (the alarm fired inside the main loop's dump)
        is not interrupted by a second one."""
        if self.dumping:
            return
        self.dumping = True
        try:
            self.stamp_new()
            meta = dict(index=self.index, seed=self.seed, batch=self.batch, run_id=self.run_id, commit=self.commit,
                        slice=self.slice_info, started_utc=self.started_utc, written_utc=utc_now(),
                        ascents=self.ascents, surgery_done=self.surgery_done, stopped_by=self.stopped_by)
            atomic_write_json(dict(N=self.C.N, d=self.C.d, types=self.C.sorted(), worker=meta), self.out)
            self.last_dump = time.time()
        finally:
            self.dumping = False

    def arm(self):
        """Schedule the next alarm: at the end of the current phase or after dump_seconds, whichever is first."""
        signal.setitimer(signal.ITIMER_REAL, max(0.01, min(self.dump_seconds, self.limit - time.time())))

    def on_alarm(self, signum, frame):
        """Alarm: leave the phase if its time is up, otherwise dump if the dump interval has passed."""
        if time.time() >= self.limit:
            raise _Stop()
        if time.time() - self.last_dump >= self.dump_seconds:
            self.dump()
        self.arm()

    def on_term(self, signum, frame):
        """SIGTERM: leave the current phase and finish with a final dump."""
        self.terminated = True; self.stopped_by = 'sigterm'
        raise _Stop()

    def phase(self, limit, work):
        """Run work() until it returns or the phase limit (absolute time) passes.  Returns True if it returned."""
        self.limit = limit; self.arm()
        try:
            work(); return True
        except _Stop:
            if not self.terminated:
                self.stopped_by = 'budget'
            return False
        finally:
            signal.setitimer(signal.ITIMER_REAL, 0)


def _start_slice(path, N, d, slice_index, slice_count, surgery_types):
    """Copies of this worker's slice of the types of `path` that no surgery has used yet, with hits reset to 0."""
    S, _ = load_census_file(path)
    if (S.N, S.d) != (N, d):
        raise ValueError('start census is for N=%d d=%d, not N=%d d=%d' % (S.N, S.d, N, d))
    pending = [kv for kv in sorted(S.types.items(), key=lambda kv: (-kv[1]['r'], kv[0])) if 'surg' not in kv[1]]
    mine = pending[slice_index::slice_count]
    if surgery_types:
        mine = mine[:surgery_types]
    out = {}
    for key, t in mine:
        t = copy.deepcopy(t); t['hits'] = 0; out[key] = t
    return out


def run_worker(N, d, out, seed=0, index=0, runs=50, rounds=1, budget=600.0, near=8, dump_seconds=180.0, start=None,
               slice_index=0, slice_count=1, surgery_types=0, surgery_share=0.7, commit='unknown', run_id='', batch=0):
    """Run one worker and write its census to `out`.  Returns the "worker" metadata of the final dump."""
    if not (N >= 3 and d >= 2 and runs >= 0 and rounds >= 0 and budget >= 0 and near >= 1 and dump_seconds > 0
            and slice_count >= 1 and 0 <= slice_index < slice_count and index >= 0 and 0 < surgery_share <= 1):
        raise ValueError('invalid worker parameter')
    w = _Worker(N, d, out, seed, index, dump_seconds, commit, run_id, batch,
                [slice_index, slice_count] if start else None)
    rng = np.random.default_rng(np.random.SeedSequence([seed, index]))
    deadline = w.t0 + budget

    def do_ascents():
        # One ascent per call, so a stop loses at most the ascent in progress and w.ascents counts finished ones.
        while w.ascents < runs:
            random_ascent(w.C, 1, rng)
            w.ascents += 1
            w.stamp_new()

    def do_surgery():
        # Types in the census when the call starts are its queue; they count as used once the call returns.
        queue = list(w.C.types)
        if queue and rounds:
            surgery(w.C, max_rounds=rounds, near=near)
            for key in queue:
                w.C.types[key]['surg'] = index
        w.surgery_done = True

    old_alarm = signal.signal(signal.SIGALRM, w.on_alarm); old_term = signal.signal(signal.SIGTERM, w.on_term)
    try:
        w.dump()   # an empty census at once: the file exists from the start, and a signal after this point is handled
        if start:
            w.C.types.update(_start_slice(start, N, d, slice_index, slice_count, surgery_types))
            w.phase(min(deadline, w.t0 + surgery_share * budget), do_surgery)
            if not w.terminated:
                w.phase(deadline, do_ascents)
        else:
            w.phase(deadline, do_ascents)
            if not w.terminated and w.ascents >= runs:
                w.phase(deadline, do_surgery)
    except _Stop:
        pass   # a signal arrived between two phases; on_term has recorded why
    finally:
        signal.signal(signal.SIGTERM, signal.SIG_IGN)   # a second SIGTERM must not interrupt the final dump
        signal.setitimer(signal.ITIMER_REAL, 0)
        w.dump()
        signal.signal(signal.SIGALRM, old_alarm); signal.signal(signal.SIGTERM, old_term)
    return dict(ascents=w.ascents, surgery_done=w.surgery_done, stopped_by=w.stopped_by, types=len(w.C.types),
                seconds=round(time.time() - w.t0, 1))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('N', type=int, help='number of points')
    ap.add_argument('d', type=int, help='ambient dimension (3: points on S^2)')
    ap.add_argument('--out', required=True, help='census file to write')
    ap.add_argument('--seed', type=int, default=0, help='base seed')
    ap.add_argument('--index', type=int, default=0, help='worker index; the generator is seeded from (seed, index)')
    ap.add_argument('--runs', type=int, default=50, help='random ascents')
    ap.add_argument('--rounds', type=int, default=1, help='surgery rounds per surgery call')
    ap.add_argument('--budget', type=float, default=600.0, help='wall-clock budget in seconds')
    ap.add_argument('--near', type=int, default=8, help='nearest non-contacts that surgery tries (as in scripts/ascent_census.py)')
    ap.add_argument('--dump-seconds', type=float, default=180.0, help='longest time between two dumps')
    ap.add_argument('--start', help='merged census whose unused types give this worker its surgery slice')
    ap.add_argument('--slice-index', type=int, default=0, help='this worker takes types i, i+W, i+2W, ... of the pending list')
    ap.add_argument('--slice-count', type=int, default=1, help='W, the number of workers sharing the pending list')
    ap.add_argument('--surgery-types', type=int, default=0, help='at most this many slice types (0: all of the slice)')
    ap.add_argument('--surgery-share', type=float, default=0.7, help='fraction of the budget that surgery may use when --start is given')
    ap.add_argument('--commit', default='unknown', help='kisscrit commit recorded in the provenance')
    ap.add_argument('--run-id', default='', help='run identifier recorded in the provenance')
    ap.add_argument('--batch', type=int, default=0, help='batch number recorded in the provenance')
    a = ap.parse_args(argv)
    try:
        r = run_worker(a.N, a.d, a.out, a.seed, a.index, a.runs, a.rounds, a.budget, a.near, a.dump_seconds, a.start,
                       a.slice_index, a.slice_count, a.surgery_types, a.surgery_share, a.commit, a.run_id, a.batch)
    except ValueError as e:
        ap.error(str(e))
    print('worker %d: %d types, %d ascents, surgery %s, stopped by %s, %.0f s' % (
        a.index, r['types'], r['ascents'], 'done' if r['surgery_done'] else 'not finished', r['stopped_by'], r['seconds']), flush=True)


if __name__ == '__main__':
    main()
