"""Finite_terms_simplified.tex from static.tex + final_checks.json + compact_simpl.pkl"""
import json
import os
import pickle
import re
import env  # noqa: F401
import render

HERE = os.path.dirname(os.path.abspath(__file__))
J = lambda f: json.load(open(os.path.join(HERE, f)))


def e(v):
    if v == 0:
        return '$0$'
    m, x = ('%.1e' % v).split('e')
    return r'$%s\times10^{%d}$' % (m, int(x))


ROWK = [
    ('G1-II (2a, 2b)', r'Group I, Row II (2a, 2b), colour $\OO_{\rm II}$ ($+$c.c.)', r'\KK_{\rm II}'),
    ('G1-III', r'Group I, Row III, colour $\OO_{\rm III}$ ($+$c.c.)', r'\KK_{\rm III}'),
    ('G1-IV', r'Group I, Row IV, colour $\OO_{\rm IV}$ ($+$c.c.)', r'p_8\KK_{\rm IV}'),
    ('G1-V(1)', r'Group I, Row V part 1, colour $\OO_{\rm V1}$ ($+$c.c.); $\KK_{\rm V1}$ is the whole of Eq.~\eqref{eq:V1} before $\OO_{\rm V1}$, '
                r'with its prefactor $\frac{ig^4}{32\pi^5k^+}$ and its two phases', r'\KK_{\rm V1}'),
    ('G2-I(1) A', r"Group II, Row I part 1, block A (no c.c.; block B is $[\dots]^*_{w\leftrightarrow w'}$ of this)", r'(-p_{16})\KK_{\rm IV}'),
    ('G3-I', r'Group III, Row I, colour $\OO_{\rm G3I}$ ($+$c.c.)', r'(-p_8)\KK_{\rm IV}'),
    ('G3-II', r'Group III, Row II, colour $\OO_{\rm G3II}$ ($+$c.c.)', r'(-p_8)\KK_{\rm IV}'),
    ('G3-III(1)', r'Group III, Row III part 1, colour $\OO_{\rm III(1)}$ ($+$c.c.)', r'\tfrac{c_{\rm III(1)}}{k^+}\KK_C'),
    ('G3-V', r'Group III, Row V, colour $\OO_{\rm V}$ of Eq.~\eqref{eq:Crow} ($+$c.c.)', r'\tfrac{c_{\rm V}}{k^+}\KK_C'),
    ('G3-VI', r'Group III, Row VI, colour $\OO_{\rm VI}$ ($+$c.c.)', r'\tfrac{c_{\rm VI}}{k^+}\KK_C'),
]


def remnants():
    S = pickle.load(open(os.path.join(HERE, 'compact_simpl.pkl'), 'rb'))
    out = []
    for key, title, kname in ROWK:
        g = S[key]['groups']
        out.append(r'\paragraph{%s.}' % title)
        out.append('\\begin{align}\n\\OO^{(2)}\\text{-part}={}' + render.row_block(g, kname) + '\n\\end{align}')
    return '\n'.join(out)


if __name__ == '__main__':
    F = J('final_checks.json')
    sp_ = F['sympy']
    split = F['split']
    rep = {}
    cc = split['cc']
    rows = []
    for k, v in cc.items():
        name = k.replace("<->", r"\leftrightarrow ")
        name = re.sub(r'\[(.*)\]', lambda m: r'[$' + m.group(1) + '$]', name)
        verdict = 'self-conjugate' if v < 1e-10 else 'not self-conjugate'
        rows.append(r'%s & %s (%s)\\' % (name, e(v), verdict))
    rows.append("G2-IV LO structure [$y\\leftrightarrow y'$, $x\\leftrightarrow x'$] & %s (self-conjugate)\\\\" % e(split['g2iv_cc']))
    rep['CCTAB'] = '\n'.join(rows)
    rep['CCSUM'] = '%s / $\\ge1$' % e(max(v for k, v in cc.items() if k.startswith('G2')))
    rep['B0'] = '$%s$' % sp_['beta0_rows_II_V']
    rep['POLE'] = '$%s$' % sp_['pole_units_as2pi_m2eps_Nc'].replace('log', r'\log').replace('*', '')
    rep['FIN'] = '$%s$' % sp_['finite_minus_expected']
    rep['LVCOEF'] = "$\\frac{\\alpha_sN_c}{2\\pi}(2\\ell_\\Delta-\\ell_X-\\ell_{X'})$"
    lim = split['limit']
    rep['SPLIT'] = r'$\le$' + e(max(r['rel'] for v in lim.values() for r in v))
    rep['G2I2'] = e(F['G2_I2_collected'])
    cc_ = F['compact_colour']
    rep['FAMIV'] = e(cc_['K_IV family'])
    rep['G1G3'] = e(cc_['G1IV-G3I'])
    rep['G2G3'] = e(cc_['G2A+G3II'])
    rep['FAMC'] = e(cc_['K_C family'])
    rep['SYMCHK'] = r'$\le$' + e(max(cc_.values()))
    cs = J('compact_simpl.json')
    rep['REMCHK'] = r'$\le$' + e(max(v for r in cs.values() for v in r['check'].values()))
    cf = F['closed_forms_vs_file']
    vals = [max(max(x) for x in v) for k, v in cf.items() if isinstance(v, list) and v and isinstance(v[0], list)]
    rep['CF'] = r'$\le$' + e(max(vals)) + r' $=O(\Lambda/k^+)$'
    rep['MAST'] = e(cf['masters'])
    rep['REMNANTS'] = remnants()
    tex = open(os.path.join(HERE, 'static.tex')).read()
    for k, v in rep.items():
        tex = tex.replace('@@%s@@' % k, v)
    left = re.findall(r'@@[A-Z0-9]+@@', tex)
    assert not left, left
    open(os.path.join(HERE, 'Finite_terms_simplified.tex'), 'w').write(tex)
    print('written')
