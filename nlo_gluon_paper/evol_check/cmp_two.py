r"""Evolution_two_rho.tex, row by row, against the log(k+/Lambda) two-rho terms (J) of your Group rows (pointwise,
same labels; log(k+/Lambda) = 1, epsilon = 0.37 finite)."""
import os, sys, re, json
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'full_check'))
import ll  # noqa
import three_check as T3
import texeval as TE
from rows_user import rows

MAP = {43: ['G1 II(1) #0'], 55: ['G1 V(2a) #2', 'G1 V(2a) #3'], 58: ['G1 V(2b UV) #0'],
       63: ['G1 V(2b fin) #3', 'G1 V(2b fin) #4', 'G1 V(2b fin) #5'], 71: ['G2 I(2b) #0'], 82: ['G2 IV #1'],
       98: ['G3 III(2) #1']}
J = {lab: t for lab, t, cc in T3.row_terms('J', True)}


def prep(t):
    t = TE.unwrap_colours(t) if hasattr(TE, 'unwrap_colours') else t
    t = t.replace(r'\Bar', r'\bar')
    return t


def val(t, c):
    return TE.term_value(t, T3.env(c))


def run(user_text=None, seeds=(1, 2, 3)):
    U = {r[3]: prep(r[4]) for r in rows('Evolution_two_rho.tex')}
    if user_text:
        U.update(user_text)
    out = {}
    for ln, labs in MAP.items():
        res = []
        for s in seeds:
            c = T3.configs(s)
            a = val(U[ln], c)
            b = sum(val(J[l], c) for l in labs)
            res.append((a, b))
        rel = max(abs(a - b) / max(abs(b), 1e-300) for a, b in res)
        ratio = [complex(a / b) for a, b in res]
        out[ln] = dict(ours=labs, rel=float(rel), ratio=[str(np.round(r, 6)) for r in ratio])
    return out


if __name__ == '__main__':
    o = run()
    tag = '_corrected' if os.environ.get('EVOL_DIR') else ''
    json.dump(o, open(os.path.join(HERE, 'cmp_two%s.json' % tag), 'w'), indent=1)
    for k, v in o.items():
        print(k, v)
