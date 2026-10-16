"""Build the LaTeX tables of Appendix B from results/maxima_*.json and results/pack_euler.json.

Status: untrusted (exploration or tooling; no part of the proof).

"""
import json, os, sys
R = os.path.join(os.path.dirname(__file__), '..', 'results')
NAMES = [('sq', 'square'), ('disk', 'disk'), ('tor', 'square torus'), ('hex', 'hexagonal torus'),
         ('cube', 'cube'), ('ball', 'ball'), ('tor3', 'cubic $3$-torus')]


def poly_tex(co):
    if not co:
        return '---'
    co = list(co)
    if co[-1] < 0:
        co = [-c for c in co]
    terms = []
    for k in range(len(co) - 1, -1, -1):
        c = co[k]
        if c == 0: continue
        mon = '' if k == 0 else ('\\rho' if k == 1 else '\\rho^{%d}' % k)
        a = abs(c); cs = ('%d' % a) if (a != 1 or k == 0) else ''
        sign = '-' if c < 0 else '+'
        terms.append((sign, cs + mon))
    s = ('-' if terms[0][0] == '-' else '') + terms[0][1]
    for sg, t in terms[1:]:
        s += sg + t
    return '$' + s + '$'


def maxima_table():
    rows = []
    for w, name in NAMES:
        p = os.path.join(R, 'maxima_%s.json' % w)
        if not os.path.exists(p): continue
        D = json.load(open(p)); first = True
        for n in sorted(D, key=int):
            v = D[n]
            note = []
            if v['rattlers']: note.append('%d rattler%s' % (len(v['rattlers']), 's' if len(v['rattlers']) > 1 else ''))
            if v['nzero']: note.append('degenerate')
            rows.append('%s & %s & $%.10f$ & %s & %d & %d%s\\\\' % (
                name if first else '', n, v['best_rho'], poly_tex(v['poly']), v['k'], v['n_maxima'],
                (' (' + ', '.join(note) + ')') if note else ''))
            first = False
        rows.append('\\hline')
    return '\n'.join(rows)


def euler_table():
    p = os.path.join(R, 'pack_euler.json')
    if not os.path.exists(p): return ''
    D = json.load(open(p)); rows = []
    for w, name in NAMES:
        ks = sorted([k for k in D if k.rsplit('_', 1)[0] == w], key=lambda k: int(k.rsplit('_', 1)[1]))
        for k in ks:
            v = D[k]; n = k.rsplit('_', 1)[1]
            tgt = '---' if v['target'] is None else str(v['target'])
            ok = 'yes' if v['target'] is not None and v['signed'] == v['target'] and not v['skipped'] else 'no'
            if v['target'] is None: ok = '---'
            rows.append('%s & %s & %d & %d & %s & %s\\\\' % (name, n, v['types'], v['signed'], tgt, ok))
    return '\n'.join(rows)


if __name__ == '__main__':
    print(maxima_table()); print(); print(euler_table())
