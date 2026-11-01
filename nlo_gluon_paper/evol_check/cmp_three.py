r"""Evolution_three_rho.tex, every displayed row (first form and rewritten form), against the log(k+/Lambda) three-rho
terms (J) of your Group rows.  Pointwise after symmetrization over the 6 placements of the three charges (and with the
c.c. for Groups I and III), log(k+/Lambda) = 1.  Output: cmp_three.json"""
import os, sys, json, re
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'full_check'))
import ll  # noqa
import three_check as T3
import texeval as TE
from rows_user import rows

G1III = ['G1 III N1 #1', 'G1 III N1 #4', 'G1 III N1 #6', 'G1 III N2 #4', 'G1 III N2 #5', 'G1 III N3a #1',
         'G1 III N3a #2', 'G1 III N3a #3']
MAP = {(44,): ['G1 II(2a) #0', 'G1 II(2a) #2'],
       (54, 57): G1III, (64, 67): G1III,
       (74,): ['G1 IV #2'], (79,): ['G1 IV #2'],
       (87,): ['G1 V(1) #2', 'G1 V(1) #3'], (91, 94): ['G1 V(1) #2', 'G1 V(1) #3'],
       (102,): ['G2 I(1) A #2', 'G2 I(1) B #2'], (106,): ['G2 I(1) A #2', 'G2 I(1) B #2'],
       (121,): ['G3 I #2'], (125,): ['G3 I #2'], (132,): ['G3 II #2'], (136,): ['G3 II #2'],
       (145,): ['G3 III(1) #1'], (149,): ['G3 III(1) #1'], (160,): ['G3 V #1'], (164,): ['G3 V #1'],
       (171,): ['G3 VI #1'], (175,): ['G3 VI #1']}
JF = {lab: t for lab, t, cc in T3.row_terms('J', True)}      # with D (F3) for G1 III
JW = {lab: t for lab, t, cc in T3.row_terms('J', False)}     # as written in Group_I.tex


KK = re.compile(r'\\mathcal\{K\}\s*\(([^()]*)\)\s*\\cdot\s*\\mathcal\{K\}\s*\(([^()]*)\)')


def kdot(t):
    """K(A).K(B) -> (A).(B)/((A)^2 (B)^2): a scalar before any colour factor"""
    return KK.sub(lambda m: r'\frac{(%s)\cdot(%s)}{(%s)^{2}\,(%s)^{2}}' % (m.group(1), m.group(2), m.group(1), m.group(2)), t)


def prep(t):
    t = TE.unwrap_colours(t).replace(r'\Bar', r'\bar')
    return kdot(T3.z_form(t))


def sv(t, c, cc):
    return T3.sym_value(t, c, cc)


def run(U=None, seeds=(1, 2)):
    if U is None:
        U = {r[3]: prep(r[4]) for r in rows('Evolution_three_rho.tex', stop=r'\section{All three')}
    out = {}
    for lines, labs in MAP.items():
        cc = not labs[0].startswith('G2')
        key = '+'.join(map(str, lines))
        try:
            res = {}
            for s in seeds:
                c = T3.configs(s)
                u = sum(sv(U[l], c, cc) for l in lines)
                w = sum(sv(JW[l], c, cc) for l in labs)
                f = sum(sv(JF[l], c, cc) for l in labs)
                res[s] = (u, w, f)
            rw = max(abs(u - w) / abs(w) for u, w, f in res.values())
            rf = max(abs(u - f) / abs(f) for u, w, f in res.values())
            out[key] = dict(ours=labs, rel_vs_written=float(rw), rel_vs_D=float(rf),
                            ratio=[str(np.round(complex(u / w), 6)) for u, w, f in res.values()])
        except Exception as e:  # noqa
            out[key] = dict(ours=labs, error=str(e)[:200])
        print(key, out[key])
    return out


if __name__ == '__main__':
    tag = '_corrected' if os.environ.get('EVOL_DIR') else ''
    json.dump(run(), open(os.path.join(HERE, 'cmp_three%s.json' % tag), 'w'), indent=1)
