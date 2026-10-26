"""BFKL_final.tex from static.tex, the json outputs and your Evolution_terms.tex (Secs. 3, 4 and their singlet forms)."""
import json
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
J = lambda f: json.load(open(os.path.join(HERE, f)))
EVOL = os.path.normpath(os.path.join(HERE, '..', 'bfkl_split', 'input', 'Evolution_terms.tex'))


def e(v):
    if v == 0:
        return '$0$'
    m, x = ('%.1e' % v).split('e')
    return r'$%s\times10^{%d}$' % (m, int(x))


def elist(vs):
    return ', '.join(e(v) for v in vs)


def evol_blocks():
    L = open(EVOL).read().splitlines()
    blocks = []
    i = 0
    while i < len(L):
        if L[i].startswith(r'\begin{align}'):
            j = i
            while not L[j].startswith(r'\end{align}'):
                j += 1
            blocks.append((i + 1, '\n'.join(L[i:j + 1])))
            i = j
        i += 1
    return blocks


def body(block):
    return block[block.index(r'\Bigg['):block.rindex(r'\Bigg]') + len(r'\Bigg]')]


def main():
    s = open(os.path.join(HERE, 'static.tex')).read()
    rep = {}
    t3 = J('three_rho_engine.json')
    yz = J('yz_test.json')
    t2 = J('two_rho_engine.json')
    rep['T3_W'] = r'%s--%s' % (e(min(t3['as written']['rel_residual'])), e(max(t3['as written']['rel_residual'])))
    rep['T3_A'] = r'$\le$' + e(max(t3['with F3']['rel_residual']))
    rep['T3_AB'] = r'$\le$' + e(max(yz['three-rho |sum|/sum|terms|']))
    rep['T2_W'] = '%d of %d' % (t2['as written'][0]['n_bad'], t2['as written'][0]['n_keys'])
    rep['T2_A'] = '%d of %d' % (t2['with F3'][0]['n_bad'], t2['with F3'][0]['n_keys'])
    rep['T2_AB'] = '%d of %d' % (max(yz['two-rho vs B-J, III(2) yz sign -1 (bad keys)']), t2['with F3'][0]['n_keys'])
    rep['T3_AB_LIST'] = elist(yz['three-rho |sum|/sum|terms|']) + ' (three configurations)'
    rep['T3_W_LIST'] = elist(t3['as written']['rel_residual']) + ' (five configurations)'
    ch = J('checks.json')
    rep['G1III_ID'] = e(ch['G1 III: N2 + N3a(2,3) = p8 K_h'])
    rep['FZ'] = e(ch["X'.F = (1/pi) int_z F_z"])
    rep['YZ3'] = e(ch["y'->z terms: symmetric three-rho part / sum|terms|"])
    rep['CT'] = e(ch['2b finite term 3 = K(z,z,x) - K(y,z,x)'])
    rep['LEMMA'] = ', '.join('%.6f' % r['ratio'][1] for r in ch['Fourier lemma (1+zeta) FT / limit'])
    ph = J('phi_table.json')
    for n, k in enumerate(['G1 II 2a', 'G1 III (F3)', 'K_IV family (G1 IV, G2 I(1) A, G3 I, G3 II)', 'G1 V(1)',
                           'K_C family c.c. (G3 III(1), V, VI)'], 1):
        rep['PHI_%d' % n] = e(ph[k])
    # row table
    co = J(os.path.join('..', 'bfkl_check', 'corrected_check.json'))
    va = J('validate.json')
    ROWS = [('Group I Row II 2a', 'G1-II(2a)', 'G1 II 2a'), ('Group I Row III', 'G1-III', 'G1 III (as written)'),
            ('Group I Row IV', 'G1-IV', 'G1 IV'), ('Group I Row V part 1', 'G1-V(1)', 'G1 V(1)'),
            ('Group I Row V 2a', 'G1-V(2a)', 'G1 V(2a)'), ('Group I Row V 2b finite', 'G1-V(2b fin)', 'G1 V(2b fin)'),
            ('Group II Row I part 1', 'G2-I(1)', 'G2 I(1)'), ('Group II Row I part 2', 'G2-I(2)', 'G2 I(2)'),
            ('Group III Row I', 'G3-I', 'G3 I'), ('Group III Row II', 'G3-II', 'G3 II'),
            ('Group III Row III part 1', 'G3-III(1)', 'G3 III(1)'), ('Group III Row III part 2', 'G3-III(2)', 'G3 III(2)'),
            ('Group III Row V', 'G3-V', 'G3 V'), ('Group III Row VI', 'G3-VI', 'G3 VI')]
    rows = []
    for name, kc, kv in ROWS:
        rows.append(r'%s & %s & %s\\' % (name, e(co[kc]['vs_exact']), e(va[kv])))
    rows.append(r'Group I Row III with (A) & --- & %s\\' % e(va['G1 III (F3)']))
    rep['ROWTAB'] = '\n'.join(rows)
    # two-rho table
    T = J('table2.json')
    ABBR = {'G1 II': 'I.II', 'G1 III': 'I.III', 'G1 IV': 'I.IV', 'G1 V': 'I.V', 'G2 I(1)': 'II.I(1)', 'G2 I(2)': 'II.I(2)',
            'G2 IV': 'II.IV', 'G3 I': 'III.I', 'G3 II': 'III.II', 'G3 III': 'III.III', 'G3 III yz': r'III.III$^{*}$',
            'G3 V': 'III.V', 'G3 V yz': r'III.V$^{*}$', 'G3 VI': 'III.VI', 'G3 VI yz': r'III.VI$^{*}$'}
    sections = [('pt', r'\multicolumn{4}{@{}l}{\emph{Pointwise in $z$}}'),
                ('LOG', r'\multicolumn{4}{@{}l}{\emph{$z$-independent colour, integrated over $z$: logarithms}}'),
                ('P', r'\multicolumn{4}{@{}l}{\emph{$z$-independent colour, integrated over $z$: pole terms}}'),
                ('ell0', r'\multicolumn{4}{@{}l}{\emph{$z$-independent colour, scaleless ($\int_zK(a,a,z)=0$), compared anyway}}')]
    out = [r'''{\small
\begin{longtable}{@{}>{$}l<{$}p{6.2cm}rr@{}}
\caption{Two-$\rho$ terms of the corrected list vs.\ $\BB-\JJ$, kernel structure by kernel structure (the charges at $x,x'$, both assignments).}\label{tab:two}\\
\toprule \text{structure} & rows contributing & $\BB-\JJ$ & diff.\\\midrule\endfirsthead
\toprule \text{structure} & rows contributing & $\BB-\JJ$ & diff.\\\midrule\endhead''']
    for kind, head in sections:
        out.append(head + r'\\')
        for r in T:
            if not r['key'].startswith("('%s'" % kind):
                continue
            tex = r['tex']
            tex = re.sub(r"\\mathcal K\((\w'?)-z\)\\cdot\\mathcal K\((\w'?)-z\)", r'K(\1,\2,z)', tex)
            tex = re.sub(r"\\frac\{1\}\{\((\w'?)-z\)\^2\}", r'K(\1,\1,z)', tex)
            tex = tex.replace(r'\mathcal K', r'\KK').replace(r'\;', r'\,')
            if kind == 'ell0':
                tex = re.sub(r'\\int_zK\((\w\'?),\1,z\)=0\\,', r'\\textstyle\\int_zK(\1,\1,z)\\,', tex)
            grp = ', '.join(ABBR[g] for g in r['groups'])
            if not r['groups'] and abs(r['target'][0]) < 1e-9:
                continue
            tv = r['target'][0]
            tvs = '0' if abs(tv) < 5e-4 else '%.3f' % tv
            out.append(r'%s & %s & %s & %s\\' % (tex, grp if grp else '---', tvs, e(r['residual'])))
    out.append(r'\bottomrule' + '\n' + r'\end{longtable}}')
    rep['TWOTAB'] = '\n'.join(out)
    nzero = sum(1 for r in T if not r['groups'] and abs(r['target'][0]) < 1e-9)
    rep['NKEYS'] = str(len(T))
    rep['MAXRES'] = e(max(r['residual'] for r in T)) + ('; %d structures are zero on both sides and are not listed' % nzero if nzero else '')
    # your equations
    B = dict(evol_blocks())
    LV, LK = r'\new{\log\Big[\frac{\vee}{k^{+}}\Big]}', r'\new{\log\Big[\frac{k^{+}}{\Lambda}\Big]}'
    b_op = B[802]
    b_op = b_op.replace(r'\frac{d^{3}N_{\mathrm{NLO}}}{d^{3}k}', r'\frac{d^{3}N_{\mathrm{BFKL}}}{d^{3}k}', 1)
    assert b_op.count(r'\log\Big[\frac{\vee}{\Lambda}\Big]') == 1
    long_ = r"-\,f^{abc}\,f^{dce}\,U^{ha}(w')\,U^{hd}(w)\,\rho^{b}(x')\,\rho^{e}(x) +\,f^{abc}"
    assert b_op.count(long_) == 1
    b_op = b_op.replace(long_, r"-\,f^{abc}\,f^{dce}\,U^{ha}(w')\,U^{hd}(w)\,\rho^{b}(x')\,\rho^{e}(x)\notag\\" + "\n" + r"& +\,f^{abc}")
    rep['BFKL_OP'] = b_op.replace(r'\log\Big[\frac{\vee}{\Lambda}\Big]', LV)
    INT_OLD = r'\int_{x,,z,w,w'
    INT_NEW = r'\int_{x,\new{y},z,w,w'
    b_s = B[831].replace(r'\frac{d^{3}N_{\mathrm{NLO}}}{d^{3}k}', r'\frac{d^{3}N_{\mathrm{BFKL}}}{d^{3}k}', 1)
    assert b_s.count(r'\log\Big[\frac{\vee}{\Lambda}\Big]') == 1 and b_s.count(INT_OLD) == 1
    b_s = b_s.replace(r'\end{align}', r'\label{eq:BS}' + '\n' + r'\end{align}')
    rep['BFKL_S'] = b_s.replace(r'\log\Big[\frac{\vee}{\Lambda}\Big]', LV).replace(INT_OLD, INT_NEW)
    j_s = B[737].replace(r'\frac{d^{3}N_{\mathrm{NLO}}}{d^{3}k}', r'\frac{d^{3}N_{\mathrm{JIMWLK}}}{d^{3}k}', 1)
    assert j_s.count(r'\log\Big[\frac{\vee}{\Lambda}\Big]') == 1 and j_s.count(INT_OLD) == 1
    rep['JIMWLK_S'] = j_s.replace(r'\log\Big[\frac{\vee}{\Lambda}\Big]', LK).replace(INT_OLD, INT_NEW)
    for k, v in rep.items():
        assert '@@%s@@' % k in s, k
        s = s.replace('@@%s@@' % k, v)
    left = re.findall(r'@@\w+@@', s)
    assert not left, left
    open(os.path.join(HERE, 'BFKL_final.tex'), 'w').write(s)
    print('written BFKL_final.tex')


if __name__ == '__main__':
    main()
