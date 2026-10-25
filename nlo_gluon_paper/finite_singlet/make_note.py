"""Finite_terms_singlet.tex from static.tex + project.pkl + project.json + log2_origin.json"""
import json
import os
import pickle
import re
import numpy as np
import env  # noqa: F401
import srender as R
import project as PR

HERE = os.path.dirname(os.path.abspath(__file__))
J = lambda f: json.load(open(os.path.join(HERE, f)))


def e(v):
    if v == 0:
        return '$0$'
    m, x = ('%.1e' % v).split('e')
    return r'$%s\times10^{%d}$' % (m, int(x))


LAB = {"x'": "x'", "y'": "y'", 'x': 'x', 'y': 'y'}
ROWS3 = [
    ('II', r"Group I Row II", r"\KK_{\rm II}", r"FS~(7.1)", r"\phw", 'c.c', ["x'", 'x', 'y'],
     r"$\OO_{\rm II}$"),
    ('III', r"Group I Row III", r"\KK_{\rm III}", r"FS~(7.2)", r"\phw", 'c.c', ["x'", 'x', 'y'],
     r"$\OO_{\rm III}$"),
    ('V1', r"Group I Row V part 1", r"\KK_{\rm V1}", r"the whole of FS~(7.4) before $\OO_{\rm V1}$, with its prefactor $\frac{ig^4}{32\pi^5k^+}$ and both phases", '', 'c.c',
     ["x'", 'y', "y'"], r"$\OO_{\rm V1}$"),
    ('famIV', r"The $\KK_{\rm IV}$ family (Group I Row IV, Group III Rows I, II, Group II Row I part 1)", r"p_8\KK_{\rm IV}",
     r"FS~(7.3), (8.1)", r"\phw", 'c.c. (for the Group II term: block B)', ["x'", 'x', 'y'],
     r"$\OO_{\rm IV}-\OO_{\rm G3I}-\OO_{\rm G3II}-\frac12\OO_A$"),
    ('famC', r"The $\KK_C$ family (Group III Rows III part 1, V, VI)", r"\frac{ig^4}{16\pi^5k^+}\KK_C", r"FS~(7.11), (8.3)",
     r"e^{-ik\cdot(y'-w)}", r"c.c. ($w\leftrightarrow y'$)", ["x'", 'x', 'y'],
     r"$\OO_{\rm III(1)}+2\OO_{\rm V}+2\OO_{\rm VI}$"),
]


def three_tex(S):
    out = []
    for key, title, kname, kref, phase, cc, labs, colname in ROWS3:
        groups = S['three'][key]
        out.append(r'\paragraph{%s.} Kernel $%s$: %s. Colour %s with the general-$U$ placement; charges $(%s)$; $+$%s.'
                   % (title, kname, kref, colname, ','.join(LAB[l] for l in labs), cc))
        lines = []
        for n, (pair, v) in enumerate(groups.items()):
            keep, gone = pair
            other = [l for l in labs if l not in pair][0]
            head = ('' if n == 0 else '+') + r'%s\big|_{%s\to %s}\,\PP(%s,%s)\Big[' % (kname, LAB[gone], LAB[keep], LAB[keep], LAB[other])
            body = R.bracket(v['comb'], indent=r'&\qquad\qquad')
            lines.append(head + body + r'\Big]')
        blk = r'\frac{d^3N^{(2)}}{d^3k}={}&\pre\int\,%s\Big\{' % phase + '\\notag\\\\\n&\\quad'.join(lines) + r'\Big\}+{\rm c.c.}'
        out.append('\\begin{align}\n' + blk + '\n\\end{align}')
    return '\n'.join(out)


def lo_check():
    """Tr{(U^dag(x)-U^dag(y))(U(x')-U(y'))} vs the projected LO colour (comb), general U"""
    S = pickle.load(open(os.path.join(HERE, 'project.pkl'), 'rb'))
    d = S['two']['LO']['comb']
    rel = []
    for s in (1, 2, 3):
        U = PR.randU(s)
        a = np.trace((U['x'] - U['y']).T @ (U["x'"] - U["y'"]))
        b = PR.eval_dict(d, U)
        rel.append(abs(a - b) / abs(a))
    return float(max(rel))


if __name__ == '__main__':
    S = pickle.load(open(os.path.join(HERE, 'project.pkl'), 'rb'))
    P = J('project.json')
    L2 = J('log2_origin.json')
    rep = {}
    rep['L2FIN'] = r'$\le$' + e(max(r['finite'] for r in L2['psiA_plus_psiB_vs_rowII']))
    rep['L2POL'] = r'$\le$' + e(max(r['pole'] for r in L2['psiA_plus_psiB_vs_rowII']))
    rep['L2UNR'] = r'$\le$' + e(max(max(r['finite_unrestricted'], r['pole_unrestricted']) for r in L2['psiA_plus_psiB_vs_rowII']))
    rep['L2PCS'] = r'$\le$' + e(max(L2['pieces'].values()))
    rep['SMALLA'] = ', '.join('%.4f' % r['value_over_a_log'] for r in L2['compact'])
    gs = P['general_sources']
    rep['GSMAX'] = r'$\le$' + e(max(v['general_vs_fixed23'] for v in gs.values()))
    rep['GSSYM'] = r'$\le$' + e(max(v['general_vs_file_for_symmetric_U'] for v in gs.values()))
    rep['G1MAX'] = r'$\le$' + e(max(P['group1_as_written_vs_first_principles'].values()))
    rep['FPTWO'] = r'$\le$' + e(max(v['colour_vs_first_principles'] for v in P['two'].values()))
    rep['FPMAX'] = r'$\le$' + e(max(max(v['general_vs_fixed23'] for v in gs.values()),
                                    max(P['group1_as_written_vs_first_principles'].values()),
                                    max(v['colour_vs_first_principles'] for v in P['two'].values())))
    rep['ENG'] = r'$\le$' + e(max(v['pair_rule_vs_engine'] for v in P['three'].values()))
    tr = [v['traces_vs_direct'] for v in P['two'].values()] + [g['traces_vs_direct'] for v in P['three'].values() for g in v['groups'].values()]
    tc = [v['combined_vs_direct'] for v in P['two'].values()] + [g['combined_vs_direct'] for v in P['three'].values() for g in v['groups'].values()]
    rep['TRMAX'] = r'$\le$' + e(max(tr))
    rep['TRCMB'] = r'$\le$' + e(max(tc))
    rep['LOCHK'] = 'agreement ' + e(lo_check())
    rep['LOTR'] = '$' + R.bracket(S['two']['LO']['comb']).replace(r'\,{\rm Tr}', r'{\rm Tr}', 1) + '$'
    for k in ('V2a', 'V2b1', 'V2b2', 'G2I2', 'III2'):
        rep[k] = R.bracket(S['two'][k]['comb'])
    rep['THREE'] = three_tex(S)
    tex = open(os.path.join(HERE, 'static.tex')).read()
    for k, v in rep.items():
        tex = tex.replace('@@%s@@' % k, v)
    left = re.findall(r'@@[A-Z0-9]+@@', tex)
    assert not left, left
    open(os.path.join(HERE, 'Finite_terms_singlet.tex'), 'w').write(tex)
    print('written')
