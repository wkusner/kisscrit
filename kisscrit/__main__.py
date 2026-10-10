"""Command line: python -m kisscrit {census,analyze,triangulations}"""
import argparse, json, sys
import numpy as np


def main(argv=None):
    p = argparse.ArgumentParser(prog='kisscrit')
    sub = p.add_subparsers(dest='cmd', required=True)
    c = sub.add_parser('census', help='critical types of N points on S^{d-1}')
    c.add_argument('--N', type=int, required=True); c.add_argument('--d', type=int, default=3)
    c.add_argument('--starts', type=int, default=30, help='random starts per graph (graph mode)')
    c.add_argument('--mode', choices=['graphs', 'ascent'], default='graphs',
                   help='graphs: solve the KKT system for every graph on N vertices (N <= 7); '
                        'ascent: random monotone ascent + contact-graph surgery (any N)')
    c.add_argument('--runs', type=int, default=200, help='ascent runs (ascent mode)')
    c.add_argument('--seed', type=int, default=0); c.add_argument('--out', default=None)
    a = sub.add_parser('analyze', help='analyze configurations in a census JSON (or a list of X)')
    a.add_argument('path')
    t = sub.add_parser('triangulations', help='enumerate triangulations of S^2 on n vertices')
    t.add_argument('--n', type=int, required=True); t.add_argument('--out', default=None)
    args = p.parse_args(argv)

    import kisscrit as kc
    if args.cmd == 'census':
        rng = np.random.default_rng(args.seed); C = kc.Census(args.N, args.d)
        if args.mode == 'graphs':
            from kisscrit.search import graphs_on, solve_graph
            for G in graphs_on(args.N):
                for X in solve_graph(G, args.d, rng, args.starts):
                    C.add(X, 'graph')
        else:
            kc.random_ascent(C, args.runs, rng); kc.surgery(C)
        for rec in C.sorted():
            print('r=%.12f  k=%2d  coindex=%d  nzero=%d  rattlers=%d  |G|=%s  labeled=%s' % (
                rec['r'], rec['k'], rec['coindex'], rec['nzero'], rec['rattlers'], rec.get('sym'), rec.get('labeled')))
        if args.out:
            C.dump(args.out)
    elif args.cmd == 'analyze':
        D = json.load(open(args.path))
        items = D['types'] if isinstance(D, dict) and 'types' in D else D
        for it in items:
            X = np.array(it['X'] if isinstance(it, dict) else it)
            r = kc.analyze(X)
            print('r=%.12f  k=%2d  critical=%s  coindex=%s  kernel=%s' % (
                kc.r_of_c(r['c']), r['k'], r['critical'], r.get('coindex'), r.get('kernel')))
    elif args.cmd == 'triangulations':
        from kisscrit.hull import enumerate_triangulations
        T = enumerate_triangulations(args.n)
        print(len(T), 'triangulations')
        if args.out:
            json.dump(T, open(args.out, 'w'))


if __name__ == '__main__':
    main()
