"""Row comparison on the collected singlet forms (fast): your rows vs the reference rows, U = U^T and general U.
Output: compare_rows_fast.json, rows_collected.pkl"""
import json
import os
import pickle
from collections import OrderedDict
import env  # noqa: F401
from core import Config
from render_singlet import collect, numeric_value
from compare_rows import user_rows, ref_rows

HERE = os.path.dirname(os.path.abspath(__file__))


def merge(*colls):
    out = OrderedDict()
    for c in colls:
        for kk, d in c.items():
            dd = out.setdefault(kk, OrderedDict())
            for ws, v in d.items():
                dd[ws] = dd.get(ws, 0) + v
    return out


def rel(a, b, sym, seeds=(3, 7, 11)):
    r = []
    for s in seeds:
        cfg = Config(s, sym)
        va, vb = numeric_value(a, cfg), numeric_value(b, cfg)
        r.append(abs(va - vb) / max(abs(va), abs(vb), 1e-12))
    return max(r)


if __name__ == '__main__':
    U, F, four_user = user_rows()
    Rf = ref_rows()
    UC = OrderedDict((k, collect(v)) for k, v in U.items())
    RC = OrderedDict((k, collect(v)) for k, v in Rf.items())
    pickle.dump(dict(user=UC, ref=RC), open(os.path.join(HERE, 'rows_collected.pkl'), 'wb'))
    out = OrderedDict()
    for k in RC:
        a = UC.get(k, OrderedDict())
        out[k] = dict(sym=rel(a, RC[k], True), general=rel(a, RC[k], False))
    g1u = merge(*[UC[k] for k in ('G1-I', 'G1-II', 'G1-III', 'G1-IV', 'G1-V', 'G2-IV')])
    g1r = merge(*[RC[k] for k in ('G1-I', 'G1-II', 'G1-III', 'G1-IV', 'G1-V', "G1-I'", 'G2-IV')])
    out['Group I + G2-IV'] = dict(sym=rel(g1u, g1r, True), general=rel(g1u, g1r, False))
    allu = merge(*UC.values())
    allr = merge(*RC.values())
    out['all rows'] = dict(sym=rel(allu, allr, True), general=rel(allu, allr, False))
    for k, v in out.items():
        print('%-16s U=U^T: %.1e   general U: %.1e' % (k, v['sym'], v['general']))
    json.dump(out, open(os.path.join(HERE, 'compare_rows_fast.json'), 'w'), indent=1)
