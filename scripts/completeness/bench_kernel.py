#!/usr/bin/env python3
"""Timing of the C checker kernel (checklist F4, spec v0.3 section 4.4 "Performance gate").

Status: untrusted (a measurement script; it reads certificates and prints rates, it decides nothing).

Measured on one core, sequentially, on the certificates given:
  * kernel time per unit (a unit is a leaf or a contraction slab), as the wall time of kc_check_chart over all
    run charts divided by the number of units, repeated --repeats times (median, min, max over the repeats);
    measured without the per-unit clock;
  * a breakdown per unit kind and per verdict (closed units by rule; unclosed units separately), from one more
    pass with the kernel's per-unit clock (kc_check_chart's unit_ns; CLOCK_MONOTONIC, about 2 clock reads per
    unit of overhead);
  * the fallback: the Python reference's time on the units the kernel leaves unclosed (ckernel.check_certificate
    'fallback_seconds', which includes the Python replay of those charts), once;
  * checker CPU time = kernel time (median) + fallback time, against the search's CPU time on the same
    candidates when --cert-records gives the search's result lines (fields 'seconds' and stats['boxes']), and
    per search box.  The search's 0.87 ms per box in checklist F4 is quoted, not re-measured, unless the search
    ran on this machine (then its records are a measurement).
Spec v0.3 section 4.4: "Target: checker CPU time <= search CPU time on the same candidates."  This script
prints the ratio; it does not decide the gate.

  python3 scripts/completeness/bench_kernel.py --files <certificates.jsonl> [--cert-records <cert.jsonl>] --repeats 7
"""
import argparse, collections, json, os, platform, statistics, subprocess, sys, time

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)


def load(files, limit=None):
    """Certificate lines (dicts) from the given files."""
    out = []
    for p in files:
        with open(p) as f:
            for line in f:
                out.append(json.loads(line))
                if limit and len(out) >= limit: return out
    return out


def prepared(certs):
    """Per run chart: the C graph, the decoded streams, the chart and its tolerance (built once, outside the clock)."""
    from kisscrit import ckernel as kk
    out = []
    for d in certs:
        cert = d['certificate']
        for ch in cert['charts']:
            if ch.get('not_run'): continue
            cg = kk.CGraph(cert['faces'], ch['modes'])
            out.append((d['code'], cg, kk.streams(ch), ch, float(cert['tol'])))
    return out


def one_pass(charts, timing=False):
    """Kernel wall time (s) over all charts, units, and per-unit records (kind, verdict, ns) if timing."""
    from kisscrit import ckernel as kk
    total = 0.0; units = 0; recs = []
    for code, cg, dec, ch, tol in charts:
        r = kk.check_chart(cg, ch, tol, timing=timing, decoded=dec)
        total += r['seconds']; units += len(r['verdict'])
        if r['rc'] != 0: recs.append(('error', r['error'], 0))
        if timing:
            for k, v, ns in zip(r['kind'].tolist(), r['verdict'].tolist(), r['unit_ns'].tolist()):
                recs.append((k, v, ns))
    return total, units, recs


def machine():
    """A one-line description of the machine and compiler."""
    cpu = ''
    try:
        for line in open('/proc/cpuinfo'):
            if line.startswith('model name'): cpu = line.split(':', 1)[1].strip(); break
    except OSError:
        cpu = platform.processor()
    cc = subprocess.run(['cc', '--version'], capture_output=True, text=True).stdout.splitlines()[:1]
    return dict(cpu=cpu, cores=os.cpu_count(), platform=platform.platform(), python=platform.python_version(),
                compiler=cc[0] if cc else '?', flags='-O2 -std=c11 -fno-fast-math -ffp-contract=off -fPIC (csrc/Makefile)')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--files', nargs='+', required=True, help='certificate files (one JSON line per candidate)')
    ap.add_argument('--cert-records', nargs='*', default=[], help="the search's result lines for the same candidates")
    ap.add_argument('--repeats', type=int, default=7)
    ap.add_argument('--limit', type=int, default=None)
    ap.add_argument('--out', default=None, help='write the summary JSON here (outside the repository results/)')
    a = ap.parse_args()
    from kisscrit import ckernel as kk
    certs = load(a.files, a.limit)
    charts = prepared(certs)
    one_pass(charts[:5])                                          # warm-up (library load, page faults)
    times = []; units = None
    for _ in range(a.repeats):
        t, u, _ = one_pass(charts); times.append(t); units = u
    _, _, recs = one_pass(charts, timing=True)
    by = collections.defaultdict(list)
    for k, v, ns in recs:
        if k == 'error': continue
        kind = 'slab' if k == 0 else 'leaf %d' % k
        by[(kind, kk.VERDICT[v])].append(ns)
    breakdown = [dict(kind=k, verdict=v, n=len(x), median_ns=int(statistics.median(x)), mean_ns=int(statistics.mean(x)))
                 for (k, v), x in sorted(by.items())]
    closed = [ns for k, v, ns in recs if k != 'error' and v in kk.CLOSED]
    unclosed = [ns for k, v, ns in recs if k != 'error' and v == 0]
    # the fallback, once
    fb_s = 0.0; fb_units = 0; fb_closed = 0; c_s = 0.0; verified = 0
    for d in certs:
        w = kk.check_certificate(d['certificate'], structural=False)
        fb_s += w['kernel']['fallback_seconds']; fb_units += w['kernel']['fallback']
        fb_closed += w['kernel']['fallback_closed']; verified += bool(w['verified'])
    med = statistics.median(times)
    out = dict(machine=machine(), certificates=len(certs), charts=len(charts), units=units, repeats=a.repeats,
               kernel_seconds=dict(median=med, min=min(times), max=max(times), all=times),
               kernel_ns_per_unit=dict(median=1e9 * med / units, min=1e9 * min(times) / units,
                                       max=1e9 * max(times) / units),
               closed_units=len(closed), closed_ns_per_unit_mean=(sum(closed) / len(closed)) if closed else None,
               unclosed_units=len(unclosed), unclosed_ns_per_unit_mean=(sum(unclosed) / len(unclosed)) if unclosed else None,
               fallback=dict(units=fb_units, closed=fb_closed, seconds=fb_s,
                             note='Python reference on unclosed units, incl. the Python replay of their charts'),
               front_end_verified_certificates=verified, breakdown=breakdown)
    if a.cert_records:
        recs = {}
        for p in a.cert_records:
            for line in open(p):
                r = json.loads(line); recs[r['code']] = r
        codes = [d['code'] for d in certs]
        have = [c for c in codes if c in recs]
        s_sec = sum(recs[c]['seconds'] for c in have); s_box = sum(recs[c]['stats']['boxes'] for c in have)
        checker = med + fb_s
        out['search'] = dict(candidates=len(have), seconds=s_sec, boxes=s_box,
                             ms_per_box=1e3 * s_sec / s_box if s_box else None,
                             hosts=sorted({recs[c].get('host') for c in have}))
        out['checker_cpu_seconds'] = checker
        out['ratio_checker_to_search'] = checker / s_sec if s_sec else None
        out['kernel_ns_per_search_box'] = 1e9 * med / s_box if s_box else None
    out['checklist_search_ms_per_box'] = dict(value=0.87, source='docs/checklist.md F4 (quoted, not re-measured)')
    out['kernel_ns_per_unit_vs_0.87ms'] = out['kernel_ns_per_unit']['median'] / 0.87e6
    text = json.dumps(out, indent=1)
    print(text)
    if a.out:
        if os.path.abspath(a.out).startswith(os.path.join(ROOT, 'results')): sys.exit('refusing to write under results/')
        open(a.out, 'w').write(text)


if __name__ == '__main__':
    main()
