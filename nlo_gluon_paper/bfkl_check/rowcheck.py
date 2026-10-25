"""Row by row: your log(vee/k+) terms (user_lvk.py) vs the exact coefficient of log(vee/k+) of the same row of
Finite_terms_1.tex (C_row of ../finite_v2/split.py, validated there: row - C_row L_v has a limit for vee -> infinity).
Also: the exact coefficient obtained directly from the file (../finite/ft.py) by d/dlog(vee) at vee/k+ = 1e7..1e9."""
import json
import numpy as np
import env  # noqa: F401
import split as SP
import simp_exact as SE
import user_lvk as UL
import ft

# the file's rows as functions of vee (Lambda -> 0): closed forms of simp_exact (checked against ft there)
FILE = {k: v[0] for k, v in SP.ROWS.items()}
FILE['G1-II(2a)'] = lambda c, kp, V: sum(ft.G1_II_2a(c, kp, V, 0.0))
CROW = {k: v[1] for k, v in SP.ROWS.items()}
CROW['G1-II(2a)'] = lambda c, kp: sum(ft.G1_II_2a(c, kp, np.e * kp, 0.0))      # 2a is proportional to L_v


def dlog(f, c, kp, V1=1e7, V2=1e9):
    return (f(c, kp, V2) - f(c, kp, V1)) / np.log(V2 / V1)


def run(seeds=(1, 2, 3)):
    out = {}
    for name, uf in UL.ROWS.items():
        res = []
        for s in seeds:
            c = SP.cfg(s)
            kp = 1.0
            u = complex(uf(c, kp))
            C = complex(CROW[name](c, kp))
            d = complex(dlog(FILE[name], c, kp)) if name in FILE else None
            res.append(dict(user=abs(u), exact=abs(C), rel_user_vs_exact=abs(u - C) / max(abs(C), 1e-300),
                            exact_vs_dlog=(abs(C - d) / max(abs(C), 1e-300)) if d is not None else None))
        out[name] = res
        r = res[0]
        print('%-13s |user| %.3e  |exact| %.3e  |user-exact|/|exact| %s   exact vs d/dlog(vee) of the file: %s'
              % (name, r['user'], r['exact'], ' '.join('%.2e' % x['rel_user_vs_exact'] for x in res),
                 ' '.join('%.1e' % x['exact_vs_dlog'] for x in res if x['exact_vs_dlog'] is not None)))
    return out


if __name__ == '__main__':
    json.dump(run(), open('rowcheck.json', 'w'), indent=1)
