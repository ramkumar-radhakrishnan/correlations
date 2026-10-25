"""BFKL_terms_check.tex from static.tex and the JSON outputs"""
import json
import re

J = lambda f: json.load(open(f))


def e(v):
    if v == 0:
        return '$0$'
    m, x = ('%.1e' % v).split('e')
    return r'$%s\times10^{%d}$' % (m, int(x))


R = J('rowcheck.json')
ORDER = [
    ('G1 Row II part 1', None, r'coefficient $1$ instead of $2$ (analytic)', r'\no'),
    ('G1 Row II 2a', 'G1-II(2a)', None, None),
    ('G1 Row III', 'G1-III', None, None),
    ('G1 Row IV', 'G1-IV', None, None),
    ('G1 Row V part 1', 'G1-V(1)', None, None),
    ('G1 Row V 2a', 'G1-V(2a)', None, None),
    ('G1 Row V 2b UV', None, r'same as the file (analytic)', r'\ok'),
    ('G1 Row V 2b finite', 'G1-V(2b fin)', None, None),
    ('G2 Row I part 1', 'G2-I(1)', None, None),
    ('G2 Row I part 2', 'G2-I(2)', None, None),
    ('G2 Row IV', None, r'NA; the file has $-\frac{\alpha_s}{2\pi}2N_c\Lv{\rm LO}[W]+\Lv H_4$', r'\miss'),
    ('G3 Row I', 'G3-I', None, None),
    ('G3 Row II', 'G3-II', None, None),
    ('G3 Row III part 1', 'G3-III(1)', None, None),
    ('G3 Row III part 2', 'G3-III(2)', None, None),
    ('G3 Row V', 'G3-V', None, None),
    ('G3 Row VI', 'G3-VI', None, None),
]
rows = []
for name, key, txt, verdict in ORDER:
    if key is None:
        rows.append(r'%s & %s & %s & --- \\' % (name, txt, verdict))
        continue
    r = R[key]
    rel = [x['rel_user_vs_exact'] for x in r]
    if all(x['user'] == 0 for x in r):
        a, v = r'NA (yours $=0$)', r'\miss'
    elif max(rel) < 1e-12:
        a, v = ', '.join(e(x) for x in rel), r'\ok'
    else:
        a, v = ', '.join('%.2g' % x for x in rel), r'\no'
    d = [x['exact_vs_dlog'] for x in r if x['exact_vs_dlog'] is not None]
    rows.append(r'%s & %s & %s & %s\\' % (name, a, v, (r'$\le$' + e(max(d))) if d else '---'))
C = J('cw_test.json')
F = J('ft_asym.json')
rep = {
    'ROWTAB': '\n'.join(rows),
    'CWREL': r'$\le$' + e(max(c['rel'] for c in C)),
    'CWGROW': '; '.join(', '.join('%.2f' % g for g in c['your_kernel_with_cutoff']) for c in C),
    'FTASYM': ', '.join('%.4f' % v['q2I_over_minus4pi'] for v in F['gauss']) + ' (Gaussian $h$), '
              + ', '.join('%.4f' % v['q2I_over_minus4pi'] for v in F['lorentz']) + r' ($h=(1+\Delta^2)^{-2}$)',
}
tex = open('static.tex').read()
for k, v in rep.items():
    tex = tex.replace('@@%s@@' % k, v)
left = re.findall(r'@@[A-Z0-9]+@@', tex)
assert not left, left
open('BFKL_terms_check.tex', 'w').write(tex)
print('written')
