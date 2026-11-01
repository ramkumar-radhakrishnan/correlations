r"""Evolution_three_rho.tex, Sec. 'All three rho term with c.c.' (consolidated list) and Sec. 'Cancellation of all
three rho's' (alpha, alpha', beta, beta' blocks).
 (1) each genuine row of the list vs the log(k+/Lambda) three-rho terms of your Group rows (with D, fix F3);
     the 'From four rho's' part vs the three-rho lists of 4rho.tex;
 (2) the whole list cancels (symmetrized over the charges, with c.c.);
 (3) each alpha/alpha'/beta/beta' block, written '= 0', vanishes.   Output: cmp_list.json"""
import os, sys, json, re
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import cmp_three as C
T3, TE = C.T3, C.TE
from rows_user import rows

ROW = {'Group I, Row II': ['G1 II(2a) #0', 'G1 II(2a) #2'], 'Group I, Row III': C.G1III,
       'Group I, Row IV': ['G1 IV #2'], 'Group I, Row V': ['G1 V(1) #2', 'G1 V(1) #3'],
       'Group II, Row I': ['G2 I(1) A #2', 'G2 I(1) B #2'], 'Group III, Row I': ['G3 I #2'],
       'Group III, Row II': ['G3 II #2'], 'Group III, Row III': ['G3 III(1) #1'], 'Group III, Row V': ['G3 V #1'],
       'Group III, Row VI': ['G3 VI #1']}
UNBAL = []
_tr = TE.translate


def _balanced(t):
    code = _tr(t)
    d = code.count('(') - code.count(')')
    if d > 0:                       # a \Bigg[ that is never closed in your text
        UNBAL.append(t[:80])
        code += ')' * d
    return code


TE.translate = _balanced
LOGB = r"\left[\log{\left[\frac{k^{+}}{\Lambda}\right]}+\log{\left[\frac{\vee}{k^{+}}\right]}\right]"


def prep(t):
    t = TE.unwrap_colours(t, cols=('blue', 'red', 'green', 'orange', 'black', 'violet', 'gray', 'pink', 'purple', 'yellow'))
    for s in (r"\alpha'", r'\alpha', r"\beta'", r'\beta'):
        t = t.replace(s, ' ')
    t = re.sub(r'&\s*\\-', '& -', t).replace('U ^', 'U^')
    t = t.replace(LOGB, T3.LKL)
    t = re.sub(r'=\s*0\s*$', '', t.strip())
    return C.prep(t)


def main(seeds=(1, 2)):
    R = rows('Evolution_three_rho.tex')
    out = {'rows': {}, 'blocks': {}}
    lst = [(e[2] or e[1], prep(e[4])) for e in R if e[0].startswith('All three')]
    four = [t for lab, t in T3.four_rho_three('J')]
    for s in seeds:
        c = T3.configs(s)
        unit = abs(TE.fp.G**4 / ((2 * np.pi)**3 * 16 * np.pi**5))
        # (1) row by row
        groups = {}
        for lab, t in lst:
            groups.setdefault(lab, []).append(t)
        for lab, ts in groups.items():
            if lab.startswith('From four'):
                u = sum(T3.sym_value(t, c, False) for t in ts)
                o = sum(T3.sym_value(t, c, False) for t in four)
                uc = sum(T3.sym_value(t, c, True) for t in ts)
                oc = sum(T3.sym_value(t, c, True) for t in four)
            else:
                labs = ROW[lab]
                cc = not lab.startswith('Group II,')
                u = uc = sum(T3.sym_value(t, c, False) for t in ts)       # the list is written with its c.c.
                o = oc = sum(T3.sym_value(C.JF[l], c, cc) for l in labs)
            out['rows'].setdefault(lab, {})[s] = dict(rel=float(abs(u - o) / abs(o)), rel_cc=float(abs(uc - oc) / abs(oc)),
                                                      ratio=str(np.round(complex(u / o), 6)))
        # (2) the whole list
        vals = [T3.sym_value(t, c, False) for lab, t in lst]
        tot, big = sum(vals), max(abs(v) for v in vals)
        out.setdefault('total', {})[s] = dict(residual=float(abs(tot) / unit), largest_term=float(big / unit))
        # (3) cancellation blocks
        for e in R:
            if e[0].startswith('Cancellation'):
                key = '%s / %s (line %d)' % (e[1], e[2], e[3])
                t = prep(e[4])
                try:
                    v0, v1 = T3.sym_value(t, c, False), T3.sym_value(t, c, True)
                    out['blocks'].setdefault(key, {})[s] = dict(direct=float(abs(v0) / unit), with_cc=float(abs(v1) / unit))
                except Exception as ex:  # noqa
                    out['blocks'].setdefault(key, {})[s] = dict(error=str(ex)[:150])
    return out


if __name__ == '__main__':
    o = main()
    o['unbalanced'] = sorted(set(UNBAL))
    print('unbalanced:', o['unbalanced'])
    for k in ('rows', 'total', 'blocks'):
        for a, b in o[k].items():
            print(k, a, b)
    tag = '_corrected' if os.environ.get('EVOL_DIR') else ''
    json.dump(o, open(os.path.join(HERE, 'cmp_list%s.json' % tag), 'w'), indent=1)
