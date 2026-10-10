"""Generate LogVeeKplus_BFKL.tex from pk_static.tex, rows.tex and the computed tables.
Inputs: reduced2.pkl, reduced3.pkl (reduceP.py), checks.json (pk_checks.py), states3.json (states3.py),
bfkl_singlet.json (bfkl_singlet.py)."""
import json
import pickle
import sympy as sp
import env  # noqa
from engine import tex_colour
import userP as UP
import psec

GL = ['I.II', 'I.III', 'I.IV', 'I.V', 'II.I', 'II.IV', 'III.I', 'III.II', 'III.III', 'III.V', 'III.VI']


def glabel(g):
    row, kind = g.split()
    return r'%s$_{%s}$' % (row, 'P' if kind == 'P' else r'\Lambda')


def gorder(g):
    row, kind = g.split()
    return (GL.index(row), kind)


def ctex(c, plus=True):
    c = sp.nsimplify(sp.expand(c))
    s = sp.latex(c).replace('N_{c}', 'N_c')
    s = s.replace(r'\frac{N_c}{2}', r'\tfrac{N_c}{2}').replace(r'\frac{1}{2}', r'\tfrac12').replace(r'\frac{1}{4}', r'\tfrac14')
    if plus and not s.startswith('-'):
        s = '+' + s
    return s


def cpre(c, first=False):
    """coefficient in front of a colour structure: +1 -> '+', -1 -> '-'"""
    c = sp.nsimplify(sp.expand(c))
    if c == 1:
        return '' if first else '+'
    if c == -1:
        return '-'
    return ctex(c, plus=not first) + r'\,'


def kname(key):
    out = []
    for k in key:
        if k[0] == 'inv2':
            a, b = k[1]
            out.append(r'K(%s,%s,z)' % (a, a) if b == 'z' else r'\frac{1}{(%s-%s)^2}' % (a, b))
        else:
            (a, b), (c, d) = k[1], k[2]
            if b == 'z' and d == 'z':
                out.append(r'K(%s,%s,z)' % (a, c))
            else:
                out.append(r'\KK(%s-%s)\cdot\KK(%s-%s)' % (a, b, c, d))
    out.sort(key=lambda s: s.startswith('K('))
    return r'\;'.join(out)


def colour(col):
    return tex_colour(col).replace(r'\,', r'\,\allowbreak ')


def tables(R, target):
    out = []
    for n, (key, d) in enumerate(R.items()):
        out.append(r'\subsection*{Structure %d: $%s$}' % (n + 1, kname(key)))
        out.append(r'\addcontentsline{toc}{subsection}{Structure %d}' % (n + 1))
        if target:
            hasB = any('-BFKL' in r['contrib'] for r in d['table'])
            hasJ = any('+JIMWLK' in r['contrib'] for r in d['table'])
            if not hasB and not hasJ:
                out.append(r'Neither BFKL$\otimes$LO nor JIMWLK$\otimes$LO has a term with this kernel: the $\log(\vee/k^+)$ terms cancel among themselves.')
            out.append(r'\begin{longtable}{@{}>{\raggedright\arraybackslash$}p{0.33\textwidth}<{$}>{\raggedright\arraybackslash}p{0.38\textwidth}>{$}c<{$}>{$}c<{$}>{$}c<{$}@{}}')
            out.append(r'\text{colour structure} & contributions (row: coefficient) & \text{sum} & \text{BFKL} & \text{JIMWLK}\\ \hline\endhead')
        else:
            out.append(r'\begin{longtable}{@{}>{\raggedright\arraybackslash$}p{0.42\textwidth}<{$}>{\raggedright\arraybackslash}p{0.44\textwidth}>{$}c<{$}@{}}')
            out.append(r'\text{colour structure} & contributions (row: coefficient) & \text{sum}\\ \hline\endhead')
        for row in d['table']:
            if not row['contrib']:
                continue
            notes = [(g, v) for g, v in row['contrib'].items() if g not in ('-BFKL', '+JIMWLK')]
            notes.sort(key=lambda gv: gorder(gv[0]))
            s = sp.nsimplify(sp.expand(sum(v for _, v in notes)))
            items = r';\ '.join(r'%s:\,$%s$' % (glabel(g), ctex(v, plus=False)) for g, v in notes) or r'--'
            if target:
                b = -row['contrib'].get('-BFKL', 0)
                j = row['contrib'].get('+JIMWLK', 0)
                ok = r'\ \text{\ok}' if sp.simplify(s - (b - j)) == 0 else r'\ \text{\bad}'
                out.append(r'%s & %s & %s & %s & %s%s\tabularnewline[2pt]' % (colour(row['col']), items, ctex(s, plus=False) if notes else '0',
                                                                              ctex(b, plus=False), ctex(j, plus=False), ok))
            else:
                ok = r'\ \text{\ok}' if s == 0 else r'\ \text{\bad}'
                out.append(r'%s & %s & %s%s\tabularnewline[2pt]' % (colour(row['col']), items, ctex(s, plus=False), ok))
        out.append(r'\end{longtable}\addtocounter{table}{-1}')
    return '\n'.join(out)


def bdecomp(R):
    out = [r'\begin{align}', r"\text{BFKL}\otimes\text{LO}={}&\NP\int e^{-ik\cdot(w'-w)}\int_{x,x',z}\Big\{\notag\\"]
    for key, d in R.items():
        terms = [(row['col'], -row['contrib']['-BFKL']) for row in d['table'] if '-BFKL' in row['contrib']]
        if not terms:
            continue
        out.append(r'&\quad +%s\,\Big[\notag\\' % kname(key))
        for j, (col, c) in enumerate(terms):
            end = r'\Big]' if j == len(terms) - 1 else ''
            out.append(r'&\qquad\qquad %s%s%s\notag\\' % (cpre(c), tex_colour(col), end))
    out[-1] = out[-1].rstrip('\\').rstrip().replace(r'\notag', r'\Big\}\label{eq:Bdecomp}')
    out.append(r'\end{align}')
    return '\n'.join(out)


def relations(R):
    out = []
    for n, (key, d) in enumerate(R.items()):
        if not d['relations']:
            continue
        out.append(r'\paragraph{Structure %d, $%s$.}' % (n + 1, kname(key)))
        out.append(r'\begin{align*}')
        lines = []
        for k, combo in d['relations']:
            rhs = ''.join(r'%s%s' % (cpre(c, first=(j == 0)), tex_colour(col)) for j, (col, c) in enumerate(combo))
            lines.append(r'&%s\\&\qquad=%s' % (tex_colour(d['cols'][k]), rhs))
        out.append(r'\\'.join(lines))
        out.append(r'\end{align*}')
    return '\n'.join(out)


def sci(x):
    if x == 0:
        return '0'
    m, e = ('%.1e' % x).split('e')
    e = int(e)
    return m if e == 0 else r'%s\times10^{%d}' % (m, e)


def psections():
    s = UP.PSRC['GI.III.P']
    UP.PSRC['GI.III.P'] = UP.mk('GI.III.P', 'GI.III', UP.I, r'K^{P}_{\rm I.III}', [(-1, s.kexpl[0][1]), (2, s.kexpl[1][1])])
    out = []
    for key in UP.PSRC:
        t = psec.section(key)
        t = t.replace('+- ', '-').replace('+-', '-')
        t = t.replace(r'\Big[-1\,', r'\Big[-').replace(r'\Big[ - 1\,', r'\Big[-').replace(r'+ - 1\,', '-')
        t = t.replace(r'\frac{d^3N}{d^3k}\bigg|_{\rm', r'P\bigg|_{\rm').replace(r'\frac{d^3N^{(2\rho)}}{d^3k}\bigg|_{\rm', r'P^{(2\rho)}\bigg|_{\rm')
        t = t.replace(r'\NN', r'\NP')
        out.append(t)
    return '\n'.join(out)


def main():
    R2 = pickle.load(open('reduced2.pkl', 'rb'))
    R3 = pickle.load(open('reduced3.pkl', 'rb'))
    ch = json.load(open('checks.json'))
    st3 = json.load(open('states3.json'))
    sing = json.load(open('bfkl_singlet.json'))
    T = open('pk_static.tex').read()
    T = T.replace('@@ROWS@@', open('rows.tex').read())
    rows = []
    for r in ch['rowsP']:
        w, c = r['written'], r['corrected']
        rows.append(r'%s & %d & %d / %d & %d / %d (max $%s$)\\' % (r['row'], r['rho'], w[1], w[0], c[1], c[0], sci(c[2])))
    st2 = ch['states2']
    names = ['as written', '+ P1 (Group I Row III: sign of $K_h$)', r'+ F3 carried over (Row III part 2 $\to\mathcal D$)',
             '+ F5 carried over (Group III Row III part 2: sign)']
    states = []
    for nm, (n2, v2), (n3, v3) in zip(names, st2, st3):
        states.append(r'%s & %d / %d & %d / %d\\' % (nm, v2[1], v2[0], v3[1], v3[0]))
    seeds = 'At two other random configurations (seeds %s) the corrected two-$\\rho$ total again agrees in all structures (largest difference $%s$).' % (
        ' and '.join(str(s) for s, _ in ch['seeds']), sci(max(v[2] for _, v in ch['seeds'])))
    gu = ch['generalU']
    genu = r'%d of %d structures differ' % (gu[1], gu[0])
    sg = []
    for r in sing:
        sg.append(r'%s & %d & $%.8f$ & $%.8f$\\' % ('symmetric' if r['sym'] else 'general', r['seed'], r['op'], r['eq']))
    n4 = 12
    rep = {
        '@@PSECTIONS@@': psections(),
        '@@T3TABLES@@': tables(R3, False),
        '@@T2TABLES@@': tables(R2, True),
        '@@BDECOMP@@': bdecomp(R2),
        '@@RELATIONS@@': relations(R2) + relations(R3),
        '@@ROWSP@@': '\n'.join(rows),
        '@@STATES@@': '\n'.join(states),
        '@@SEEDS@@': seeds,
        '@@GENU@@': genu,
        '@@SINGLET@@': '\n'.join(sg),
        '@@S3W@@': str(st3[0][1][1]), '@@S3N@@': str(st3[0][1][0]),
        '@@S2W@@': str(st2[0][1][1]), '@@S2N@@': str(st2[0][1][0]),
        '@@S2NC@@': str(len(R2)), '@@N3NET@@': str(sum(len(d['nets']) for d in R3.values())),
        '@@N4@@': str(n4), '@@S3NC@@': str(len(R3)),
    }
    for k, v in rep.items():
        assert k in T, k
        T = T.replace(k, v)
    assert '@@' not in T, T[T.index('@@') - 50:T.index('@@') + 50]
    open('LogVeeKplus_BFKL.tex', 'w').write(T)
    print('written LogVeeKplus_BFKL.tex')


if __name__ == '__main__':
    main()
