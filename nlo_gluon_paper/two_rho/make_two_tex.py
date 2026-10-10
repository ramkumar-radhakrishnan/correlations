"""Generate TwoRho_evolution.tex: all two-rho terms proportional to log(vee/Lambda) vs JIMWLK x LO.
Inputs: reduced.pkl (reduce.py), diag_rows.json (diag_rows.py), states.json (states.py),
appF.json (appF.py), variants.json (variants.py), two_static.tex."""
import json
import pickle
import os
import sympy as sp
import core  # noqa: F401  (paths)
from engine import tex_colour, Nc
import user4

GNAME = {'G I, Row II': r'I.II', 'G I, Row III': r'I.III', 'G I, Row IV': r'I.IV', 'G I, Row V': r'I.V',
         'G II, Row I': r'II.I', 'G II, Row IV': r'II.IV', 'G III, Row I': r'III.I', 'G III, Row II': r'III.II',
         'G III, Row III': r'III.III', 'G III, Row V': r'III.V', 'G III, Row VI': r'III.VI',
         '4rho.tex 2rho': r'4$\rho$', '4rho remainders A,B': r'A,B', 'D': r'$\mathcal D$'}
GORDER = list(GNAME)


def ctex(c, plus=True):
    c = sp.nsimplify(sp.expand(c))
    s = sp.latex(c).replace('N_{c}', 'N_c')
    s = s.replace(r'\frac{N_c}{2}', r'\tfrac{N_c}{2}').replace(r'\frac{1}{2}', r'\tfrac12').replace(r'\frac{3 N_c}{2}', r'\tfrac{3N_c}{2}')
    if plus and not s.startswith('-'):
        s = '+' + s
    return s


def kname(key):
    """canonical kernel key -> 'KK(x-w).KK(x'-w') K(a,b,z)'"""
    lo, soft = None, None
    for k in key:
        if k[0] == 'inv2':
            soft = (k[1][0], k[1][0])
        else:
            (a, b), (c, d) = k[1], k[2]
            if b == 'z' and d == 'z':
                soft = (a, c)
            else:
                lo = r'\KK(%s-%s)\cdot\KK(%s-%s)' % (a, b, c, d)
    return r'%s\;K(%s,%s,z)' % (lo, soft[0], soft[1])


def colour(col):
    return tex_colour(col).replace(r'\,', r'\,\allowbreak ')


def match_tables(R):
    out = []
    for n, (key, d) in enumerate(R.items()):
        has_j = any('-JIMWLK' in row['contrib'] for row in d['table'])
        out.append(r'\subsection*{Structure %d: $%s$}' % (n + 1, kname(key)))
        out.append(r'\addcontentsline{toc}{subsection}{Structure %d}' % (n + 1))
        if not has_j:
            out.append(r"JIMWLK$\otimes$LO has no term with this kernel: the terms of the notes must cancel among themselves.")
        out.append(r'\begin{longtable}{@{}>{\raggedright\arraybackslash$}p{0.36\textwidth}<{$}>{\raggedright\arraybackslash}p{0.40\textwidth}>{$}c<{$}>{$}c<{$}@{}}')
        out.append(r'\text{colour structure } B & contributions of the notes (row: coefficient) & \text{sum} & \text{JIMWLK}\\ \hline\endhead')
        for row in d['table']:
            if not row['contrib']:
                continue
            notes = [(g, v) for g, v in row['contrib'].items() if g != '-JIMWLK']
            notes.sort(key=lambda gv: GORDER.index(gv[0]))
            jim = -row['contrib'].get('-JIMWLK', 0)
            s = sp.nsimplify(sp.expand(sum(v for _, v in notes)))
            items = r';\ '.join(r'%s:\,$%s$' % (GNAME[g], ctex(v, plus=False)) for g, v in notes) or r'--'
            ok = r'\ \text{\ok}' if sp.simplify(s - jim) == 0 else r'\ \text{\bad}'
            out.append(r'%s & %s & %s & %s%s\tabularnewline[2pt]' % (colour(row['col']), items, ctex(s, plus=False) if notes else '0',
                                                                       ctex(jim, plus=False), ok))
        out.append(r'\end{longtable}\addtocounter{table}{-1}')
    return '\n'.join(out)


def jimwlk_decomp(R):
    out = [r'\begin{align}', r"\frac{d^3N}{d^3k}\bigg|_{\rm JIMWLK\otimes LO}={}&\NN\int e^{-ik\cdot(w'-w)}\int_{x,x',z}\Big\{\notag\\"]
    parts = []
    for key, d in R.items():
        terms = [(row['col'], -row['contrib']['-JIMWLK']) for row in d['table'] if '-JIMWLK' in row['contrib']]
        if not terms:
            continue
        inner = r'\Big[' + ''.join(r'%s\,%s' % (ctex(c), tex_colour(col)) if k else r'%s\,%s' % (ctex(c, plus=False), tex_colour(col))
                                   for k, (col, c) in enumerate(terms)) + r'\Big]'
        parts.append((kname(key), terms))
    for k, (kn, terms) in enumerate(parts):
        out.append(r'&\quad +%s\,\Big[\notag\\' % kn)
        for j, (col, c) in enumerate(terms):
            end = r'\Big]' if j == len(terms) - 1 else ''
            out.append(r'&\qquad\qquad %s\,%s%s\notag\\' % (ctex(c), tex_colour(col), end))
    out[-1] = out[-1].rstrip('\\').rstrip()
    out[-1] = out[-1].replace(r'\notag', r'\Big\}\label{eq:Jdecomp}')
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
            rhs = ''.join(r'%s\,%s' % (ctex(c, plus=(j > 0)), tex_colour(col)) for j, (col, c) in enumerate(combo))
            lines.append(r'&%s\\&\qquad=%s' % (tex_colour(d['cols'][k]), rhs))
        out.append(r'\\'.join(lines))
        out.append(r'\end{align*}')
    return '\n'.join(out)


def user4_list():
    """the lists of 4rho.tex in the common notation, with the two corrections marked"""
    def kern_tex(kern):
        out = []
        for f in kern:
            a, b, c, d = f[1:]
            if (b, d) == ('z', 'z'):
                out.append(r'K(%s,%s,z)' % (a, c))
            else:
                out.append(r'\KK(%s-%s)\cdot\KK(%s-%s)' % (a, b, c, d))
        return r'\,'.join(out)
    L = [r'\begin{longtable}{@{}>{\raggedright\arraybackslash}p{0.11\textwidth} >{$}r<{$} >{\raggedright\arraybackslash$}p{0.30\textwidth}<{$} >{\raggedright\arraybackslash$}p{0.44\textwidth}<{$}@{}}',
         r'term & \text{coeff.} & \text{kernel} & \text{colour}\\ \hline\endhead']
    for name, lst, mult in (('first part', user4.F1, r'N_c'), ('$dN^{(1)}$', user4.F2, ''), ('$dN^{(2)}$', user4.F3, '')):
        for n, (c, kern, fs, Us, rs) in enumerate(lst):
            lab = '%s, %d' % (name, n + 1)
            facs = [('f',) + tuple(f) for f in fs] + [('U', i, j, p) for i, j, p in Us] + [('rho', i, p) for i, p in rs]
            col = tex_colour(facs, keep_order=True).replace(r'\,', r'\,\allowbreak ')
            cc = ctex(sp.nsimplify(c)) + (r'\,' + mult if mult else '')
            if name == 'first part' and n == 5:
                L.append(r'\textcolor{red}{%s (drop)} & \textcolor{red}{%s} & %s & %s\tabularnewline' % (lab, cc, kern_tex(kern), col))
            elif name == '$dN^{(1)}$' and n == 1:
                c2, kern2, fs2, Us2, rs2 = user4.F2_2_FIXED
                facs2 = [('f',) + tuple(f) for f in fs2] + [('U', i, j, p) for i, j, p in Us2] + [('rho', i, p) for i, p in rs2]
                L.append(r'\textcolor{red}{%s (as printed)} & %s & %s & %s\tabularnewline' % (lab, cc, kern_tex(kern), col))
                L.append(r'\textcolor{blue}{%s (corrected)} & %s & %s & %s\tabularnewline' % (lab, cc, kern_tex(kern2), tex_colour(facs2, keep_order=True)))
            else:
                L.append(r'%s & %s & %s & %s\tabularnewline' % (lab, cc, kern_tex(kern), col))
    L.append(r'\end{longtable}\addtocounter{table}{-1}')
    return '\n'.join(L)


def d_section():
    """step-by-step reordering of D, with the generator of the previous note"""
    import make_tex as MT
    import cancel_data as CD
    from engine import Rho, pieces_for, make_piece, simplify, canonical_rename
    src = CD.D_SOURCE
    src.kname = r'K_{\mathcal D}'
    sres = {'src': src, 'subs': []}
    for n, (sign, colour, word) in enumerate(src.subterms):
        rhos = [Rho(i, p) for i, p in word]
        pcs = []
        for coeff, X, Y, Z, tag in pieces_for(src.form, rhos):
            c0, facs0, (old, new) = make_piece(colour, sign * coeff, X, Y, Z)
            c1, facs1, rules = simplify(c0, facs0)
            facs1c = canonical_rename(facs1) if facs1 else facs1
            pcs.append(dict(tag=tag, X=(X.idx, X.pos), Y=(Y.idx, Y.pos), Z=(Z.idx, Z.pos), coeff_before=c0, facs_before=facs0,
                            subst=(old, new), coeff_after=c1, facs_after=facs1c, rules=rules))
        sres['subs'].append(dict(sign=sign, colour=colour, word=word, pieces=pcs))
    txt, res = MT.gen_source(src, sres)
    txt = txt.replace(r'\subsection{Replacement $\mathcal D$ of Row III, part 2}', r'\subsubsection*{The two-$\rho$ terms of $\mathcal D$}')
    txt = txt.replace(r'Eq.~\eqref{eq:plain}', r'Eq.~(3.4) of the previous note')
    return txt


def sci(x):
    if x == 0:
        return '0'
    m, e = ('%.1e' % x).split('e')
    e = int(e)
    return m if e == 0 else r'%s\times10^{%d}' % (m, e)


def row_table():
    rows = json.load(open('diag_rows.json'))
    out = []
    for r in rows:
        w, c = r['written'], r['corrected']
        fw, fc = sci(w[0]), sci(c[0])
        verdict = 'agree' if c[0] < 1e-10 else r'agree after $\int_z$ (see text)'
        if r['row'] in ('Group II, Row IV', 'Group I, all rows'):
            verdict = r'agree after $\int_z$; the sum agrees pointwise'
        out.append(r'%s & $%s$ & $%s$ & %.1f & %s\\' % (r['row'], fw, fc, w[1], verdict))
    return '\n'.join(out)


def state_table():
    st = json.load(open('states.json'))
    out = []
    for s in st:
        name = (s['state'].replace('->', r'$\to$').replace('F_alpha', r'$\mathfrak F_\alpha$').replace('U(w)A3', r'$U(w)\mathbb A^{(3)}$')
                .replace('-8Nc', r'$-8N_c$').replace('4rho.tex', r'\texttt{4rho.tex}'))
        name = name.replace(' D', r' $\mathcal D$')
        out.append(r'%s & %d & %d \\' % (name, s['bad'], s['keys']))
    return '\n'.join(out)


def appf_table():
    a = json.load(open('appF.json'))
    out = []
    for r in a:
        out.append(r'%s & %d & $%.6f$ & $%.6f$ & $%.6f$ & $%.4f$\\' % ('symmetric' if r['sym'] else 'general', r['seed'], r['target'], r['corrected'],
                                                                      r['written'], r['bfkl']))
    return '\n'.join(out)


def variants_text():
    v = json.load(open('variants.json'))
    def f(k):
        return r'%d of %d structures differ, largest difference $%s$' % (v[k][1], v[k][0], sci(v[k][2]))
    out = [r'\begin{tabular}{@{}l l@{}}', r'variant & total minus JIMWLK$\otimes$LO\\ \hline']
    out.append(r'correction of Row II as $\mathfrak F\to\mathfrak F_\beta$ in $U(w)\mathbb A^{(3)}$ & %s\\' % f('simple'))
    out.append(r'correction of Row II in commutator form (Eq.~(7.3) of the three-$\rho$ note) & %s\\' % f('commutator'))
    for k in v:
        if k.startswith('D '):
            out.append(r'$\mathcal D$ with charge order $\rho(%s)\rho(%s)\rho(%s)$ & %s\\' % (tuple(k[2:].split()) + (f(k),)))
    for k in v:
        if k.startswith('seed'):
            out.append(r'other random configuration (%s) & %s\\' % (k, f(k)))
    out.append(r'general (non-symmetric) adjoint $U$ & %s\\' % f('general U'))
    out.append(r'\end{tabular}')
    return '\n'.join(out)


def main():
    R = pickle.load(open('reduced.pkl', 'rb'))
    T = open('two_static.tex').read()
    rep = {
        '@@MATCHTABLES@@': match_tables(R),
        '@@JDECOMP@@': jimwlk_decomp(R),
        '@@RELATIONS@@': relations(R),
        '@@USER4@@': user4_list(),
        '@@DSECTION@@': d_section(),
        '@@ROWTABLE@@': row_table(),
        '@@STATETABLE@@': state_table(),
        '@@APPFTABLE@@': appf_table(),
        '@@VARIANTS@@': variants_text(),
        '@@NSTRUCT@@': str(len(R)),
        '@@NBASIS@@': str(sum(len(d['table']) for d in R.values())),
    }
    for k, v in rep.items():
        assert k in T, k
        T = T.replace(k, v)
    assert '@@' not in T
    open('TwoRho_evolution.tex', 'w').write(T)
    print('written TwoRho_evolution.tex')


if __name__ == '__main__':
    main()
