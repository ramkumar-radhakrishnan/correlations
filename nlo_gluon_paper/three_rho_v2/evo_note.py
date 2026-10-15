"""Write Evolution_from_final.tex: how the two-rho terms give JIMWLK x LO (coefficient of log(vee/Lambda)) and
BFKL x LO (coefficient of log vee).  Inputs: evo.pkl (evo_build.py), ThreeRho_final.tex."""
import pickle
import re
import sys
import sympy as sp
import env
from book import Book
from engine import tex_colour
import two_from_final as TF
sys.path.insert(0, env.ALL2)
import make_doc as MD
import sets as S2

GCOL = {'3rho': r'\text{3$\rho$ eq.}', 'gen': r'\text{genuine}', '4rho': r'\text{\texttt{4rho.tex}}'}
GNAME = {'3rho': r'from the three-$\rho$ equation (all rows)', 'gen': r'genuine two-$\rho$ terms', '4rho': r'two-$\rho$ lists of \texttt{4rho.tex}'}
PNAME = {'3rho': r'three-$\rho$ rows ($\vee$ end)', 'gen': r'genuine rows ($\vee$ end)', '4rho': r'two-$\rho$ lists of \texttt{4rho.tex} ($P=A$)'}


def entries_block(rows):
    """rows: list of (title, terms) -> align body lines, one entry per collected term"""
    body, n = [], 0
    for title, terms in rows:
        B = Book()
        for t in terms:
            B.add(t, 'g')
        B.finalize()
        bk = B.dump()
        ents = [(key, sp.expand(net['groups'].get('g', 0)), net['col']) for key in MD.ordered(bk) for net in bk[key]]
        ents = [e for e in ents if e[1] != 0]
        body.append(r'&\quad\underline{\text{%s}}\notag\\' % title)
        if not ents:
            body.append(r'&\quad 0\notag\\')
        for key, c, col in ents:
            n += 1
            body.append(r'&\quad %s\;%s\;%s\notag\\' % (TF.coef_tex(c), TF.kern_tex(key), tex_colour(col)))
    return body, n


def align(lhs, body, label):
    out = [r'\begin{align}', r"%s={}&\int e^{-ik(w'-w)}\,\Bigg[\notag\\" % lhs] + body
    out[-1] = out[-1].replace(r'\notag\\', r'\Bigg]\label{%s}' % label)
    out.append(r'\end{align}')
    return '\n'.join(out)


def target_eq(R, keys, tname, lhs, label, tag):
    body = []
    for n, key in enumerate(keys):
        terms = [(row['col'], -row['contrib'][tname]) for row in R[key]['table'] if tname in row['contrib']]
        if not terms:
            continue
        body.append(r'&\quad +\underbrace{%s}_{%s_{%d}}\Big[\notag\\' % (TF.kern_tex(key), tag, n + 1))
        for j, (col, c) in enumerate(terms):
            end = r'\Big]' if j == len(terms) - 1 else ''
            body.append(r'&\qquad\qquad %s\;%s%s\notag\\' % (TF.coef_tex(c), tex_colour(col), end))
    return align(lhs, body, label)


def sum_tables(R, keys, tname, tcol, tag):
    out = []
    for n, key in enumerate(keys):
        rows = [r for r in R[key]['table'] if r['contrib']]
        out.append(r'\subsubsection*{Structure $%s_{%d}$: $%s$}' % (tag, n + 1, TF.kern_tex(key)))
        if not rows:
            out.append(r'Every colour structure has coefficient 0 in each of the three sets and in %s.' % tcol)
            continue
        if not any(tname in r['contrib'] for r in rows):
            out.append(r'%s has no term with this kernel; the terms of the three sets cancel.' % tcol)
        out.append(r'\begin{longtable}{@{}>{\raggedright\arraybackslash$}p{0.42\textwidth}<{$}>{$}c<{$}>{$}c<{$}>{$}c<{$}>{$}c<{$}>{$}c<{$}@{}}')
        out.append(r'\text{colour structure} & %s & %s & %s & \text{sum} & \text{%s}\\ \hline\endhead' % (GCOL['3rho'], GCOL['gen'], GCOL['4rho'], tcol))
        for r in rows:
            c = r['contrib']
            a, g, f = (c.get(k, 0) for k in ('3rho', 'gen', '4rho'))
            j = -c.get(tname, 0)
            s = sp.nsimplify(sp.expand(a + g + f))
            ok = r'\ \checkmark' if sp.simplify(s - j) == 0 else r'\ \times'
            out.append(r'%s & %s & %s & %s & %s & %s%s\tabularnewline[2pt]' % (
                MD.colour(r['col']), MD.ctex(a, plus=False), MD.ctex(g, plus=False), MD.ctex(f, plus=False),
                MD.ctex(s, plus=False), MD.ctex(j, plus=False), ok))
        out.append(r'\end{longtable}\addtocounter{table}{-1}')
    return '\n'.join(out)


def relations(R, keys, tag):
    out = []
    for n, key in enumerate(keys):
        d = R[key]
        if not d['relations']:
            continue
        out.append(r'\paragraph{Structure $%s_{%d}$.}' % (tag, n + 1))
        out.append(r'\begin{align*}')
        lines = []
        for k, combo in d['relations']:
            rhs = ''.join(r'%s%s' % (MD.cpre(c, first=(j == 0)), tex_colour(col)) for j, (col, c) in enumerate(combo))
            lines.append(r'&%s\\&\qquad=%s' % (tex_colour(d['cols'][k]), rhs))
        out.append(r'\\'.join(lines))
        out.append(r'\end{align*}')
    return '\n'.join(out)


def worked(R, key, tname, tcol, names, tag, n):
    """one structure in full: each set before the reduction, the identities, the sum, the target"""
    d = R[key]
    nets = d['nets']
    cols = d['cols']
    rel = {k: combo for k, combo in d['relations']}
    kt = TF.kern_tex(key)
    out = [r'The kernel is $%s$ (structure $%s_{%d}$). Every term below multiplies it.' % (kt, tag, n)]
    for g in ('3rho', 'gen', '4rho'):
        terms = [(cols[k], net['groups'][g]) for k, net in enumerate(nets) if net['groups'].get(g, 0) != 0]
        out.append(r'\paragraph{%s.}' % names[g][0].upper().__add__(names[g][1:]))
        if not terms:
            out.append(r'No terms with this kernel.')
            continue
        out.append(r'\begin{align*}')
        out.append(r'\\'.join(r'&%s\;%s' % (TF.coef_tex(c), tex_colour(col)) for col, c in terms))
        out.append(r'\end{align*}')
    used = [k for k, net in enumerate(nets) if k in rel and any(v != 0 for gg, v in net['groups'].items() if gg != tname)]
    if used:
        out.append(r'\paragraph{Colour identities.} The colour structures are not independent. In this structure we use')
        out.append(r'\begin{align*}')
        out.append(r'\\'.join(r'&%s\\&\qquad=%s' % (tex_colour(cols[k]), ''.join(r'%s%s' % (MD.cpre(c, first=(j == 0)), tex_colour(col)) for j, (col, c) in enumerate(rel[k])))
                              for k in used))
        out.append(r'\end{align*}')
    out.append(r'\paragraph{Sum of the three sets, in independent colour structures, and %s.}' % tcol)
    out.append(r'\begin{center}\small\begin{tabular}{@{}l c c c c c@{}}')
    out.append(r'colour structure & 3$\rho$ eq. & genuine & \texttt{4rho.tex} & sum & %s\\ \hline' % tcol)
    for r in d['table']:
        c = r['contrib']
        if not c:
            continue
        a, g, f = (c.get(k, 0) for k in ('3rho', 'gen', '4rho'))
        j = -c.get(tname, 0)
        s = sp.nsimplify(sp.expand(a + g + f))
        out.append(r'$%s$ & $%s$ & $%s$ & $%s$ & $%s$ & $%s$\\' % (tex_colour(r['col']), MD.ctex(a, plus=False), MD.ctex(g, plus=False),
                                                              MD.ctex(f, plus=False), MD.ctex(s, plus=False), MD.ctex(j, plus=False)))
    out.append(r'\end{tabular}\end{center}')
    return '\n'.join(out)


def find_key(keys, tex):
    for n, k in enumerate(keys):
        if TF.kern_tex(k) == tex:
            return n + 1, k
    raise KeyError(tex)


def main():
    E = pickle.load(open('evo.pkl', 'rb'))
    RA, RP = E['A']['reduced'], E['P']['reduced']
    kA, kP = MD.ordered(RA), MD.ordered(RP)
    A = S2.A_sets()
    P = E['P_sets']
    gen_rows = [('Group I, Row II (fix F6)', A['g:I.II']), ('Group I, Row V', A['g:I.V']), ('Group II, Row I', A['g:II.I']),
                ('Group II, Row IV', A['g:II.IV']), ('Group III, Row III, part 2 (fix F5)', A['g:III.III'])]
    genA, nga = entries_block(gen_rows)
    fourA, nfa = entries_block([(r'\texttt{4rho.tex} (fixes F1, F2)', A['4rho'])])
    prow = [('Group I, Row II (fix F4: the $\\mathfrak F_\\alpha$ part, with the sign reversed)', P['3:I.II']),
            ('Group I, Row III (with fix P1)', P['3:I.III']), ('Group I, Row IV', P['3:I.IV']), ('Group II, Row I', P['3:II.I']),
            ('Group III, Row I', P['3:III.I']), ('Group III, Row II', P['3:III.II']), ('Group III, Row III, part 1', P['3:III.III']),
            ('Group III, Row V', P['3:III.V']), ('Group III, Row VI', P['3:III.VI']),
            ("Four-$\\rho$ remainders $d^3N_a$, $d^3N_b$ ($P=A$)", P['A,B']),
            ('Genuine: Group II, Row I', P['g:II.I']), ('Genuine: Group II, Row IV', P['g:II.IV']),
            ('Genuine: Group III, Row III, part 2 (fix F5)', P['g:III.III'])]
    plist, npl = entries_block(prow)
    nA, keyA = find_key(kA, r"\mathcal{K}(x-w)\cdot\mathcal{K}(x'-w')\,\mathcal{K}(x-z)\cdot\mathcal{K}(x'-z)")
    nA2, keyA2 = find_key(kA, r"\mathcal{K}(x-w)\cdot\mathcal{K}(x'-w')\,\mathcal{K}(x-z)\cdot\mathcal{K}(w-z)")
    nP, keyP = find_key(kP, r"\mathcal{K}(w-z)\cdot\mathcal{K}(w'-z)\,\mathcal{K}(x-z)\cdot\mathcal{K}(x'-z)")
    nP2, keyP2 = find_key(kP, r"\mathcal{K}(x-w)\cdot\mathcal{K}(x'-w')\,\mathcal{K}(x-z)\cdot\mathcal{K}(x'-z)")
    struct = '\n'.join(r'$S_{%d}$ & $%s$ & %s\\' % (n + 1, TF.kern_tex(k), r'\checkmark' if any('-JIMWLK' in r['contrib'] for r in RA[k]['table']) else '--')
                       for n, k in enumerate(kA))
    structP = '\n'.join(r"$S'_{%d}$ & $%s$ & %s\\" % (n + 1, TF.kern_tex(k), r'\checkmark' if any('-BFKL' in r['contrib'] for r in RP[k]['table']) else '--')
                        for n, k in enumerate(kP))
    nb = lambda R: sum(1 for k in R for r in R[k]['table'] if r['contrib'])
    rep = {
        '@@GENA@@': align(r'\frac{d^3N^{2\rho}_{\rm gen}}{d^3k}\Big|_{\log(\vee/\Lambda)}\Big/\NO\LL', genA, 'eq:gen'),
        '@@FOURA@@': align(r'\frac{d^3N^{2\rho}_{\rm 4\rho}}{d^3k}\Big|_{\log(\vee/\Lambda)}\Big/\NO\LL', fourA, 'eq:four'),
        '@@PLIST@@': align(r'P\big/\NO', plist, 'eq:Plist'),
        '@@JEQ@@': target_eq(RA, kA, '-JIMWLK', r'\J', 'eq:Jexp', 'S'),
        '@@BEQ@@': target_eq(RP, kP, '-BFKL', r'\B', 'eq:Bexp', "S'"),
        '@@JTAB@@': sum_tables(RA, kA, '-JIMWLK', 'JIMWLK', 'S'),
        '@@BTAB@@': sum_tables(RP, kP, '-BFKL', 'BFKL', "S'"),
        '@@RELA@@': relations(RA, kA, 'S'), '@@RELP@@': relations(RP, kP, "S'"),
        '@@WA@@': worked(RA, keyA, '-JIMWLK', r'JIMWLK$\otimes$LO', GNAME_CAP(GNAME), 'S', nA),
        '@@WA2@@': worked(RA, keyA2, '-JIMWLK', r'JIMWLK$\otimes$LO', GNAME_CAP(GNAME), 'S', nA2),
        '@@WP@@': worked(RP, keyP, '-BFKL', r'BFKL$\otimes$LO', GNAME_CAP(PNAME), "S'", nP),
        '@@WP2@@': worked(RP, keyP2, '-BFKL', r'BFKL$\otimes$LO', GNAME_CAP(PNAME), "S'", nP2),

        '@@STRUCTA@@': struct, '@@STRUCTP@@': structP,
        '@@NBA@@': str(nb(RA)), '@@NBP@@': str(nb(RP)), '@@NSA@@': str(len(kA)), '@@NSP@@': str(len(kP)),
        '@@NGEN@@': str(nga), '@@NFOUR@@': str(nfa), '@@NPL@@': str(npl),
    }
    T = open('evo_static.tex').read()
    for k, v in rep.items():
        assert k in T, k
        T = T.replace(k, v)
    left = re.findall(r'@@[^@]*@@', T)
    assert not left, left
    open('Evolution_from_final.tex', 'w').write(T)
    print('written Evolution_from_final.tex', dict(gen=nga, four=nfa, P=npl, basisA=nb(RA), basisP=nb(RP)))


def GNAME_CAP(d):
    return {k: v for k, v in d.items()}


if __name__ == '__main__':
    main()
