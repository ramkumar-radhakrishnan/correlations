"""Structure-by-structure table of the two-rho terms of the corrected list against B - J (seed 1 values, residuals for 3 seeds).
Output: table2.json  (key, LaTeX of the kernel, contributions per source group, B - J, residual)"""
import json
from collections import OrderedDict, defaultdict
import numpy as np
import env  # noqa: F401
import two_rho_engine as TE
import yz_test as YZ

GROUP = OrderedDict([
    ('G1 II', ['G1 II 2a', 'G1 II(1)']), ('G1 III', ['G1 III']), ('G1 IV', ['G1 IV']),
    ('G1 V', ['G1 V(1)', 'G1 V(2a)', 'G1 V(2b UV)', 'G1 V(2b fin)']), ('G2 I(1)', ['G2 I(1)']), ('G2 I(2)', ['G2 I(2)']),
    ('G2 IV', ['G2 IV']), ('G3 I', ['G3 I']), ('G3 II', ['G3 II']),
    ('G3 III', ['G3 III(1)', 'G3 III(2)']), ('G3 III yz', ['G3 III(1) yz', 'G3 III(2) yz']),
    ('G3 V', ['G3 V']), ('G3 V yz', ['G3 V yz']), ('G3 VI', ['G3 VI']), ('G3 VI yz', ['G3 VI yz'])])
INV = {s: g for g, ss in GROUP.items() for s in ss}
PT = {'p1': 'x', 'p2': "x'"}


def v(a, b):
    a, b = PT.get(a, a), PT.get(b, b)
    return '%s-%s' % (a, b)


def tex_factor(k):
    if k[0] == 'KK':
        (a, b), (c, d) = k[1], k[2]
        return r'\mathcal K(%s)\cdot\mathcal K(%s)' % (v(a, b), v(c, d))
    if k[0] == 'inv2':
        a, b = k[1]
        return r'\frac{1}{(%s)^2}' % v(a, b)
    raise ValueError(k)


def tex_key(key):
    kind = key[0]
    if kind == 'pt':
        return ''.join(tex_factor(k) for k in key[1:])
    if kind == 'P':
        return r'\pi P\;' + ''.join(tex_factor(k) for k in key[1:])
    if kind == 'LOG':
        a, b = key[1], key[2]
        return r'-\pi\log(%s)^2\;' % v(a, b) + ''.join(tex_factor(k) for k in key[3:])
    if kind == 'ell0':
        return r'\int_zK(%s,%s,z)=0\;' % (PT.get(key[1], key[1]), PT.get(key[1], key[1])) + ''.join(tex_factor(k) for k in key[2:])
    raise ValueError(key)


if __name__ == '__main__':
    J, B4, Bref = TE.targets()
    FT = YZ.all2(-1)
    rows = []
    accs = []
    for seed in (1, 2, 3):
        groups = dict(FT)
        groups['B - J'] = TE.neg(Bref) + J[-1]['terms']          # stored with minus sign: residual = file - (B - J)
        accs.append(TE.keyed(groups, seed))
    keys = sorted(accs[0].keys(), key=lambda k: (['pt', 'LOG', 'P', 'ell0'].index(k[0]), str(k)))
    for key in keys:
        d = accs[0][key]
        con = defaultdict(complex)
        for lab, val in d.items():
            if lab == 'B - J':
                continue
            con[INV[lab]] += val
        target = -d.get('B - J', 0.0)
        res = [abs(sum(a[key].values())) / max(1.0, sum(abs(x) for x in a[key].values())) for a in accs]
        rows.append(dict(key=str(key), tex=tex_key(key),
                         groups={g: [float(np.real(c)), float(np.imag(c))] for g, c in con.items() if abs(c) > 1e-9},
                         target=[float(np.real(target)), float(np.imag(target))], residual=max(res)))
    print(len(rows), 'structures; max residual', max(r['residual'] for r in rows if not r['key'].startswith("('ell0'")))
    json.dump(rows, open('table2.json', 'w'), indent=1)
    for r in rows[:80]:
        print('%-90s %s  B-J=%.3f  res=%.1e' % (r['tex'][:90], ','.join(r['groups']), r['target'][0], r['residual']))
