r"""Targeted UV scan (one configuration, p+ = 0.37 k+): the pairs of points that can diverge in Group I Row V and
Group III Rows III, V, VI, per term and per equation.   python3 uvtarget.py <filter> -> uvtarget_<filter>.json"""
import json
import os
import re
import sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import uvscan as UV  # noqa: E402
import uvscan2 as U2  # noqa: E402
import surfscan as SS  # noqa: E402

PAIRS = {
    'Row V': [('x', 'z'), ('x', 'y'), ('y', 'z'), ('x', 'yp'), ('xp', 'wp')],
    'Row III': [('yp', 'z'), ('x', 'z'), ('y', 'w'), ('xp', 'yp'), ('z', 'w'), ('x', 'w'), ('xp', 'z')],
    'Row VI': [('yp', 'z'), ('x', 'z'), ('y', 'w'), ('xp', 'yp'), ('xp', 'z')],
    'Row I': [('z', 'w'), ('z', 'wp'), ('x', 'z'), ('x', 'w'), ('xp', 'wp')],
}


def main(only):
    out = {}
    fn = os.path.join(HERE, 'uvtarget_%s.json' % re.sub(r'\W+', '_', only))
    for key, texs in U2.groups():
        if only not in key:
            continue
        try:
            ts = [UV.Term(t) for t in texs]
        except Exception as ex:
            out[key] = {'skipped': str(ex)[:200]}
            print(key, 'skipped', str(ex)[:80], flush=True)
            continue
        U2.ACTIVE[:] = [f for f, _, pats in SS.SURF.values() if any(SS.has(t, pats) for t in texs)]
        rowname = 'Row III' if ('Group III' in key and 'Row V' in key and 'Row VI' not in key) else 'Row VI' if 'Row VI' in key else ('Row V' if 'Row V' in key else ('Row I' if 'Row I /' in key else 'Row III'))
        names = sorted(set().union(*[set(t.pts) for t in ts]), key=UV.NAMES.index)
        pairs = []
        for P, Q in PAIRS[rowname]:
            if P in names and Q in names:
                pairs.append(tuple(sorted((P, Q), key=UV.NAMES.index)))
        field = UV.Field(0)
        pts, kv = UV.base_points(0)
        res = {'sum': U2.scan_group(ts, field, pts, kv, 0.37, pairs=pairs)}
        if len(ts) > 1:
            res['terms'] = [U2.scan_group([t], field, pts, kv, 0.37, pairs=pairs) for t in ts]
        out[key] = res
        print(key, '| sum:', {k: v['s'] for k, v in res['sum'].items()},
              '| terms:', [{k: v['s'] for k, v in r.items() if v['s'] < -1.7} for r in res.get('terms', [])], flush=True)
        json.dump(out, open(fn, 'w'), indent=1)


if __name__ == '__main__':
    main(sys.argv[1])
