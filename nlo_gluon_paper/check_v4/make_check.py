"""Writes Evolution_terms_check.tex from check_static.tex and the results of check_all.py, check_singlet.py,
check_m1.py, m1_singlet.py, make_corrected.py."""
import json
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
J = lambda f: json.load(open(os.path.join(HERE, f)))


def e(v):
    if v == 0:
        return '$0$'
    m, x = ('%.1e' % v).split('e')
    return r'$%s\times10^{%d}$' % (m, int(x))


def verdict(r, kind):
    good = {'target': r'\ok{}', 'same': r'\textcolor{ForestGreen}{agrees}', 'ident': r'\textcolor{ForestGreen}{identical}'}[kind]
    if r['sym'] < 1e-10 and r['general'] < 1e-10:
        return good + ('' if kind == 'ident' else ', any $U$')
    if r['sym'] < 1e-10:
        return good + r' for $U^{bc}=U^{cb}$'
    return r'\textcolor{red}{differs}'


def row(name, r, extra=None, kind='target'):
    cells = [name] + ([extra] if extra else []) + [e(r['sym']), e(r['general']), verdict(r, kind)]
    return ' & '.join(cells) + r'\\'


GEN = [('(1.1) Group I, Row II', 'Group I, Row II'),
       ('(1.2) Group I, Row V, $d^3N_{2a}$', 'Group I, Row V, 2a'),
       ('(1.3) Group I, Row V, $d^3N_{2b\\,\\rm UV}$', 'Group I, Row V, 2b UV'),
       ('(1.4) Group I, Row V, $d^3N_{2b\\,\\rm finite}$', 'Group I, Row V, 2b finite (first)'),
       ('(1.5) Group I, Row V, second ``$d^3N_{2b\\,\\rm finite}$\'\' vs counterterm', 'Group I, Row V, 2b finite (second) = counterterm'),
       ('(1.3) + (1.5), relative to (1.3) (vanishes)', 'Group I, Row V, 2b UV + second "2b finite" (should vanish)'),
       ('(1.2)--(1.5) together vs Row V', 'Group I, Row V, all four'),
       ('(1.6) Group II, Row I', 'Group II, Row I'),
       ('(1.7) Group II, Row IV', 'Group II, Row IV'),
       ('(1.8) Group III, Row III, second part (as written)', 'Group III, Row III (second part), as written'),
       ('(1.8) with its c.c.\\ added once more', 'Group III, Row III (second part), written + c.c.')]


def sec2_name(k, n):
    sub, subsub, part = [x.strip() for x in k.split(' / ')]
    sub = sub.replace('\\,', '').replace(' (with c.c.)', '')
    sub = 'From four $\\rho$\'s' if sub.startswith('From') else sub
    s = '(2.%d) %s' % (n, sub)
    if subsub:
        s += ', ' + subsub
    if part:
        s += ', $N_%s$' % part[1]
    return s


if __name__ == '__main__':
    A = J('check_all.json')
    S = J('check_singlet.json')
    M = J('check_m1.json')
    MS = J('m1_singlet.json')
    C = J('corrected.json')
    rep = {}
    g = A['genuine']
    rep['GENTAB'] = '\n'.join(row(n, g[k], kind='same') for n, k in GEN) + '\n\\midrule\n' + row('all of Sec.~1', A['genuine, all'], kind='same')
    lines = []
    for n, (k, r) in enumerate(A['sec2_vs_old'].items(), 1):
        lines.append(row(sec2_name(k, n), r, '%d / %d' % (r['n_new'], r['n_old']), kind='ident'))
    rep['SEC2TAB'] = '\n'.join(lines)
    s3, s4 = A['sec3'], A['sec4']
    L = [row('Secs.~1+2 vs JIMWLK$\\otimes$LO', A['Secs 1+2 vs JIMWLKxLO']),
         row('Secs.~1+2 + \\texttt{4rho.tex} lists (P1) vs JIMWLK$\\otimes$LO', A['Secs 1+2 + 4rho.tex lists (M1) vs JIMWLKxLO']),
         r'\midrule',
         row('(3.1) + (3.2) vs (3.3)', s3['dN1 + dN2 vs collected form'], kind='same'),
         row('(3.3) vs JIMWLK$\\otimes$LO', s3['collected form vs JIMWLKxLO']),
         row('(3.3) vs Secs.~1+2 as written', s3['collected form vs Secs 1+2']),
         r'\midrule',
         row('(4.1) vs (4.2)', s4['first form vs collected form'], kind='same'),
         row('(4.2) vs BFKL$\\otimes$LO', s4['collected form vs BFKLxLO']),
         row('(4.1) vs BFKL$\\otimes$LO', s4['first form vs BFKLxLO'])]
    rep['S34TAB'] = '\n'.join(L)
    sj = lambda k: dict(sym=S['U=U^T'][k], general=S['general U'][k])
    rep['SINGTAB'] = '\n'.join([
        row('(3.4) vs singlet projection of (3.3)', sj('J_singlet_vs_proj_sec3'), kind='same'),
        row('(3.4) vs Eq.~(3.1) of \\texttt{Evolution\\_singlet.pdf}', sj('J_singlet_vs_note')),
        row('(4.3) vs singlet projection of (4.2)', sj('B_singlet_vs_proj_sec4'), kind='same'),
        row('(4.3) vs Eq.~(3.2) of \\texttt{Evolution\\_singlet.pdf}', sj('B_singlet_vs_note'))])
    rep['J3'] = e(s3['collected form vs JIMWLKxLO']['general'])
    rep['B4'] = e(s4['collected form vs BFKLxLO']['general'])
    rep['SJP'] = e(S['general U']['J_singlet_vs_proj_sec3'])
    rep['SBP'] = e(S['general U']['B_singlet_vs_proj_sec4'])
    rep['SJN'] = e(S['general U']['J_singlet_vs_note'])
    rep['SBN'] = e(S['general U']['B_singlet_vs_note'])
    pt = M['|Secs1+2 - J| / |J|, U=U^T, pointwise']
    rep['M1PT'] = '%d--%d\\%%' % (round(100 * min(pt)), round(100 * max(pt)))
    ms = MS['U=U^T']
    rep['M1SING'] = '%d--%d\\%%' % (round(100 * min(ms)), round(100 * max(ms)))
    rep['M1FIX'] = e(A['Secs 1+2 + 4rho.tex lists (M1) vs JIMWLKxLO']['sym'])
    rep['FOURCHK'] = e(M['explicit list vs checked 4rho.tex lists']['sym'])
    rep['UVCT'] = e(g['Group I, Row V, 2b UV + second "2b finite" (should vanish)']['sym'])
    rep['G3CC'] = e(g['Group III, Row III (second part), written + c.c.']['sym'])
    rep['CORR'] = e(max(C['S12_vs_J']['sym'], C['S12_vs_sec3']['sym']))
    rep['NT'] = '624'
    body = open(os.path.join(HERE, 'fourlists_body.tex')).read()
    body = re.sub(r'\\color\{\w+\}', r'\\color{blue}', body).strip()
    rep['FOUR'] = body
    tex = open(os.path.join(HERE, 'check_static.tex')).read()
    for k, v in rep.items():
        tex = tex.replace('@@%s@@' % k, v)
    left = re.findall(r'@@[A-Z0-9]+@@', tex)
    assert not left, left
    open(os.path.join(HERE, 'Evolution_terms_check.tex'), 'w').write(tex)
    print('written')
