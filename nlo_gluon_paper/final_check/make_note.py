"""Writes Final_evolution_check.tex from final_static.tex, the books of render.py (books.pkl), the user's Sec. 2 lists
(coloured verbatim) and the json results of the checks."""
import json
import pickle
import re
import sys
import sympy as sp
import env  # noqa
from engine import tex_colour
sys.path.insert(0, env.V2)
from two_from_final import coef_tex, kern_tex  # noqa: E402
import make_doc as MD  # noqa: E402
from parse_user import parse_chunk  # noqa: E402
from render import klass_key, klass_term, COL  # noqa: E402

D = pickle.load(open('books.pkl', 'rb'))
R = json.load(open('render.json'))
T = json.load(open('total_check.json'))
C2 = json.load(open('compare_sec2.json'))
TH = json.load(open('three_cancel.json'))
TG = json.load(open('targets_check.json'))


def e(x, d=0):
    if x == 0:
        return '0'
    m, p = ('%.*e' % (d, x)).split('e')
    return r'%s\times10^{%d}' % (m, int(p))


def entries(bk, group, sign=1):
    out = []
    for key in MD.ordered(bk):
        for net in bk[key]:
            c = sp.expand(net['groups'].get(group, 0))
            if c != 0:
                out.append((key, sign * c, net['col']))
    return out


def eqn(head, ents, label, norm=r'\NN'):
    L = [r'\begin{align}', r"%s={}&%s\int_{x,x',z,w,w'}e^{-ik\cdot(w'-w)}\Bigg[\notag\\" % (head, norm)]
    body = []
    for key, c, col in ents:
        k = klass_key(key)
        body.append(r'&\quad{\color{%s}%s\;%s\;%s}\notag\\' % (COL[k], coef_tex(c), kern_tex(key), tex_colour(col)))
    if not body:
        body = [r'&\quad 0\notag\\']
    body[-1] = body[-1][:-len(r'\notag\\')] + r'\Bigg]\label{%s}' % label
    return '\n'.join(L + body + [r'\end{align}'])


# ------------------------------------------------------------------ genuine terms (expanded, z-integrand form)
GEN_TITLES = {
    'G1-II': ('Group I, Row II', r'G1\text{-}II'),
    'G1-V 2a': ('Group I, Row V, part 2a', r'G1\text{-}V\,2a'),
    'G1-V 2b endpoints': (r'Group I, Row V, part 2b finite: the two endpoint terms ($B_2$, $B_4$)', r'G1\text{-}V\,2b\,{\rm end}'),
    'G1-V 2b UV': (r'Group I, Row V, part 2b UV (bracket written as $\frac1\pi\int_x[K(z,z,x)-K(y,z,x)]$)', r'G1\text{-}V\,2b\,{\rm UV}'),
    'G1-V 2b CT': (r'Group I, Row V, part 2b finite: the counterterm', r'G1\text{-}V\,2b\,{\rm CT}'),
    'G2-I': ('Group II, Row I', r'G2\text{-}I'),
    'G2-IV': ('Group II, Row IV', r'G2\text{-}IV'),
    'G3-III(2)': ('Group III, Row III, second part', r'G3\text{-}III(2)'),
}
gen = []
for i, g in enumerate(D['gen']):
    title, short = GEN_TITLES[g]
    gen.append(r'\paragraph{%s.}' % title)
    gen.append(eqn(r'\frac{d^3N}{d^3k}\bigg|_{%s}' % short, entries(D['bkA'], g), 'eq:gen%d' % i))
GEN = '\n'.join(gen)

FOUR = eqn(r'\frac{d^3N}{d^3k}\bigg|_{4\rho\text{-lists}}', entries(D['bkA'], '4rho'), 'eq:fourlists')
JSUM = eqn(r'\frac{d^3N}{d^3k}\bigg|_{{\rm JIMWLK}\otimes{\rm LO}}', entries(D['bkA'], '-JIMWLK', -1), 'eq:Jsum')
BSUM = eqn(r'\frac{d^3N}{d^3k}\bigg|_{{\rm BFKL}\otimes{\rm LO}}', entries(D['bkP'], '-BFKL', -1), 'eq:Bsum', norm=r'\NNv')

P_TITLES = {'3:I.II': ('Group I, Row II (three $\\rho$; the $\\alpha$ part, fix F4)', r'I.II'),
            '3:I.III': ('Group I, Row III (three $\\rho$; with the sign of $p^+>k^+$, fix P1)', r'I.III'),
            '3:I.IV': ('Group I, Row IV (three $\\rho$)', r'I.IV'),
            '3:II.I': ('Group II, Row I (three $\\rho$)', r'II.I'),
            '3:III.I': ('Group III, Row I (three $\\rho$)', r'III.I'),
            '3:III.II': ('Group III, Row II (three $\\rho$)', r'III.II'),
            '3:III.III': ('Group III, Row III (three $\\rho$)', r'III.III'),
            '3:III.V': ('Group III, Row V (three $\\rho$)', r'III.V'),
            '3:III.VI': ('Group III, Row VI (three $\\rho$)', r'III.VI'),
            'g:II.I': ('Group II, Row I (genuine two $\\rho$)', r'II.I^{\rm gen}'),
            'g:II.IV': ('Group II, Row IV (genuine two $\\rho$)', r'II.IV^{\rm gen}'),
            'g:III.III': ('Group III, Row III, second part (genuine two $\\rho$)', r'III.III^{\rm gen}'),
            'A,B': ('From four $\\rho$\'s: remainders $d^3N_a$, $d^3N_b$', r'4\rho,\,a\text{-}b'),
            '4rho': ('From four $\\rho$\'s: the two-$\\rho$ lists of \\texttt{4rho.tex}', r'4\rho\text{-lists}')}
plist = []
for i, k in enumerate(D['P']):
    title, short = P_TITLES[k]
    ents = entries(D['bkP'], k)
    plist.append(r'\paragraph{%s.}' % title)
    plist.append(eqn(r'\frac{d^3N}{d^3k}\bigg|^{\log\vee}_{%s}' % short, ents, 'eq:P%d' % i, norm=r'\NNv'))
PLIST = '\n'.join(plist)


# ------------------------------------------------------------------ the user's Sec. 2, coloured verbatim
def colour_sec2(path):
    tex = open(path).read()
    tex = tex[tex.index(r'\section{Symmetrizing three and four'):tex.index(r'\end{document}')]
    tex = tex.replace("\\section{Symmetrizing three and four $\\rho$'s}", '')
    tex = re.sub(r'\\subsection\{([^}]*)\}', r'\\subsubsection*{\1}', tex)
    tex = re.sub(r'\\subsubsection\{([^}]*)\}', r'\\paragraph{\1}', tex)
    tex = tex.replace("\\textbf{From three $\\rho$'s}", "\\subsubsection*{From three $\\rho$'s}")

    def col_align(m):
        body = m.group(1)
        parts = re.split(r'(\\notag\s*\\\\|&)', body)
        out = []
        for p in parts:
            if r'\rho^' in p and r'\mathcal{K}' in p:
                t = parse_chunk(p, 'x')
                k = klass_term(t)
                core = p.rstrip()
                tail = ''
                if core.endswith(r'\Bigg]'):
                    core, tail = core[:-len(r'\Bigg]')], r'\Bigg]'
                out.append(r'{\color{%s}%s}%s ' % (COL[k], core, tail))
            else:
                out.append(p)
        return r'\begin{align}' + ''.join(out) + r'\end{align}'
    tex = re.sub(r'\\begin\{align\}(.*?)\\end\{align\}', col_align, tex, flags=re.S)
    return tex


SEC2 = colour_sec2('input/Evolution_terms.tex')


# ------------------------------------------------------------------ tables
def struct_table(per):
    rows = []
    for key, k, nb, bad in per:
        key = eval(key) if isinstance(key, str) else key
        em, soft = MD.split_key(key)
        name = r'\;'.join(em + soft).replace(r'\KK', r'\mathcal K')
        rows.append(r'$%s$ & {\color{%s}%s} & %d & %d\\' % (name, COL[k], {'J': 'JIMWLK only', 'C': 'charge--charge', 'B': 'BFKL only'}[k], nb, bad))
    return '\n'.join(rows)


TABA = struct_table(R['perA'])
TABP = struct_table(R['perP'])
nofour_bad = [r for r in R['perNoFour'] if r[3]]
NOFOUR = ', '.join(r'$%s$' % r'\;'.join(MD.split_key(eval(r[0]))[0] + MD.split_key(eval(r[0]))[1]).replace(r'\KK', r'\mathcal K') for r in nofour_bad)

cnt = R['counts']
sub = {
    'SEC2DIFF': ', '.join('%s: %d' % (k, v['differing_kernel_keys']) for k, v in C2.items()),
    'SEC2MAX': e(max(v['max_rel_pointwise'] for v in C2.values())),
    'N2': str(sum(v['n_user'] for v in C2.values())),
    'TOT_NO4': e(T['sec2+genuine (as in Evolution_terms.tex)']['max_rel_pointwise'], 1),
    'TOT_4': e(T['sec2+genuine+4rho.tex lists']['max_rel_pointwise']),
    'DIFF_NO4': str(T['sec2+genuine (as in Evolution_terms.tex)']['differing']),
    'NKEY': str(T['sec2+genuine (as in Evolution_terms.tex)']['kernel_monomials']),
    'A_STRUCT': str(R['A']['structures']), 'A_BASIS': str(R['A']['basis']), 'A_BAD': str(R['A']['nonzero']),
    'A4_BASIS': str(R['A_without_4rho']['basis']), 'A4_BAD': str(R['A_without_4rho']['nonzero']),
    'P_STRUCT': str(R['P']['structures']), 'P_BASIS': str(R['P']['basis']), 'P_BAD': str(R['P']['nonzero']),
    'NOFOUR': NOFOUR,
    'TH_N': str(TH['n_terms']), 'TH_RES': e(max(r['max_residual'] for r in TH['cancellation_U_symmetric'])),
    'TH_NZ': str(sum(r['nonzero'] for r in TH['cancellation_U_symmetric'])),
    'APPF_SENT': '%.2f--%.2f' % tuple(sorted(abs(x['AppF_as_sent'] - x['JIMWLKxLO']) / abs(x['JIMWLKxLO']) for x in TG)[::len(TG) - 1]),
    'APPF_CORR': e(max(abs(x['AppF_corrected'] - x['JIMWLKxLO']) / abs(x['JIMWLKxLO']) for x in TG)),
    'APPE': e(max(abs(x['AppE'] - x['BFKLxLO']) / abs(x['BFKLxLO']) for x in TG)),
    'CNT_USER': '%d / %d' % (cnt['user_sec2'].get('J', 0), cnt['user_sec2'].get('C', 0)),
    'CNT_GEN': '%d / %d' % (cnt['gen'].get('J', 0), cnt['gen'].get('C', 0)),
    'CNT_FOUR': '%d / %d' % (cnt['four'].get('J', 0), cnt['four'].get('C', 0)),
    'CNT_J': '%d / %d' % (cnt['J'].get('J', 0), cnt['J'].get('C', 0)),
    'CNT_P': '%d / %d' % (cnt['P'].get('B', 0), cnt['P'].get('C', 0)),
    'CNT_B': '%d / %d' % (cnt['B'].get('B', 0), cnt['B'].get('C', 0)),
    'GEN': GEN, 'FOUR': FOUR, 'SEC2': SEC2, 'JSUM': JSUM, 'PLIST': PLIST, 'BSUM': BSUM, 'TABA': TABA, 'TABP': TABP,
}
s = open('final_static.tex').read()
for k, v in sub.items():
    s = s.replace('@@%s@@' % k, v)
left = re.findall(r'@@[A-Z0-9_]+@@', s)
assert not left, left
open('Final_evolution_check.tex', 'w').write(s)
print('wrote Final_evolution_check.tex')
