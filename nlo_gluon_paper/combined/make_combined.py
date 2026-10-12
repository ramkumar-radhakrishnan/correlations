"""Write RowII_FT_and_evolution.tex: one note with
   Part I   the Fourier transform of Group I, Row II (ft_part.tex, new);
   Part II  the cancellation of the three-rho terms (../three_rho/ThreeRho_cancellation.tex);
   Part III the two-rho terms vs JIMWLK x LO (../two_rho/TwoRho_evolution.tex);
   Part IV  the coefficient of log(vee) vs BFKL x LO (../bfkl/pk_static.tex, rows.tex, new P tables);
   appendices: the reorderings (../three_rho/ThreeRho_to_TwoRho.tex and the region-P sources), colour relations, checks.
Labels of the included notes get the prefixes R:, C:, T:, B:.
Inputs: ft_checks.json (ft/ft_checks.py), allstates.json (check_all.py), reducedP2.pkl, reducedP3.pkl (build_P.py + reduceP.py),
and the json files of ../bfkl."""
import json
import os
import pickle
import re
import sympy as sp
import env  # noqa
import make_pk_tex as MPK
from engine import tex_colour

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = env.ROOT
PATHS = {'R': os.path.join(ROOT, 'three_rho', 'ThreeRho_to_TwoRho.tex'),
         'C': os.path.join(ROOT, 'three_rho', 'ThreeRho_cancellation.tex'),
         'T': os.path.join(ROOT, 'two_rho', 'TwoRho_evolution.tex')}
BF = os.path.join(ROOT, 'bfkl')

# labels of the included notes that are defined in sections we do not include -> labels of this note
MAP = {
    'R:sec:conv': 'sec:conv', 'R:eq:NN': 'eq:NN', 'R:eq:comm': 'eq:comm', 'R:eq:U1': 'eq:U1', 'R:eq:U2': 'eq:U2',
    'R:eq:Usym': 'eq:Usym', 'R:app:F': 'sec:ft-steps', 'R:sec:intro': 'part:three',
    'C:sec:conv': 'sec:conv', 'C:eq:cc': 'eq:cc', 'C:app:scripts': 'app:scripts',
    'T:sec:conv': 'sec:conv', 'T:app:scripts': 'app:scripts', 'T:app:rel': 'app:rel',
    'B:sec:conv': 'sec:conv', 'B:eq:comm': 'eq:comm', 'B:eq:plain': 'R:eq:plain', 'B:eq:Aanti': 'R:eq:Aanti',
    'B:eq:antiA': 'R:eq:antiA', 'B:sec:reorder': 'app:reorderP', 'B:app:rel': 'app:rel', 'B:sec:answer': 'sec:bfkl-overview',
}
LABRE = re.compile(r'\\(label|ref|eqref|pageref)\{([^}]*)\}')


def prefix(text, p):
    def f(m):
        lab = m.group(2)
        if lab.startswith('!'):
            return '\\%s{%s}' % (m.group(1), lab[1:])
        new = '%s:%s' % (p, lab)
        return '\\%s{%s}' % (m.group(1), MAP.get(new, new))
    return LABRE.sub(f, text)


def body(path):
    t = open(path).read()
    return t[t.index(r'\tableofcontents') + len(r'\tableofcontents'):t.rindex(r'\end{document}')]


def sections(text):
    """{label: text} for every \\section, in order; appendix sections included."""
    text = text.replace(r'\appendix', '')
    parts = re.split(r'(?m)^(?=\\section\{)', text)
    out = {}
    for s in parts:
        if not s.startswith(r'\section{'):
            continue
        m = re.search(r'\\label\{([^}]*)\}', s)
        out[m.group(1)] = s
    return out


def rep(text, pairs, src):
    for a, b in pairs:
        if a not in text:
            raise SystemExit('replacement not found in %s: %s' % (src, a[:80]))
        text = text.replace(a, b)
    return text


def retitle(sec, title):
    return re.sub(r'^\\section\{.*?\}\\label', lambda m: r'\section{%s}\label' % title, sec, count=1, flags=re.S)


# ------------------------------------------------------------------------------------------ Part I
def ft_part():
    ch = json.load(open(os.path.join(HERE, 'ft', 'ft_checks.json')))
    assert ch['partial_fractions'] and ch['E_plus_D'] and ch['A_plus_B_full'] == ['0', '0']
    xi = '\n'.join(r'$10^{%d}$ & $%.6f$ & $%.6f$ & $%.6f$\\' % (round(__import__('math').log10(r['vee'])), r['main'], r['full'], r['msg'])
                   for r in ch['xi'])
    lv = '\n'.join(r'%d & $10^{%d}$ & $%s$ & $%.2f$\\' % (r['trial'] + 1, round(__import__('math').log10(r['vee'])), MPK.sci(r['diff']), r['size'])
                   for r in ch['large_vee'])
    t = open(os.path.join(HERE, 'ft_part.tex')).read()
    t = t.replace('@@XITABLE@@', xi).replace('@@LARGEV@@', lv).replace('@@COMM@@', '$' + MPK.sci(ch['commutator']['identity']) + '$')
    assert '@@' not in t
    return t


# ------------------------------------------------------------------------------------------ Part II
def part_C():
    S = sections(body(PATHS['C']))
    keep = ['sec:answer', 'sec:method', 'sec:entries', 'sec:networks', 'sec:diagnosis', 'sec:corrections', 'sec:after',
            'sec:pointwise', 'sec:Usym']
    t = '\n'.join(S[k] for k in keep)
    t = retitle(t, 'Overview')
    t = rep(t, [
        (r'The previous note (\texttt{ThreeRho\_to\_TwoRho}) reordered every', r'Apps.~\ref{!R:sec:identity}--\ref{!R:sec:summary} reorder every'),
        (r'(previous note, App.~A)', r'(Part~\ref{!part:ft}, Eq.~\eqref{!eq:Fab})'),
        (r'Before computing the two-$\rho$ terms of Row II, check which form the derivation of $\mathbb A^{(3)}$ gives.',
         r'Part~\ref{!part:jimwlk} found that both forms give the same two-$\rho$ terms (Sec.~\ref{!T:sec:num}). '
         r'Sec.~\ref{!sec:alpha} shows what correction~1 means: the $\mathfrak F_\alpha$ part is the commutator $[\mathbb N^{(2)},\mathbb A^{(1)}]$, '
         r'and correction~1 is the $p^+$ ordering of the virtual factor.'),
        (r'The previous note derived its two-$\rho$ terms from the terms as written.',
         r'Apps.~\ref{!R:sec:fourrho}--\ref{!R:sec:summary} derive the two-$\rho$ terms from the terms as written.'),
    ], 'C')
    return prefix(t, 'C')


# ------------------------------------------------------------------------------------------ Part III
def part_T():
    S = sections(body(PATHS['T']))
    keep = ['sec:answer', 'sec:which', 'sec:rowV', 'sec:all', 'sec:target', 'sec:match', 'sec:num']
    t = '\n'.join(S[k] for k in keep)
    t = retitle(t, 'Overview')
    t = rep(t, [
        (r'(note \texttt{ThreeRho\_to\_TwoRho}, including the remainders A1--B2 of the four-$\rho$ terms)',
         r'(Apps.~\ref{!R:sec:fourrho}--\ref{!R:sec:summary}, including the remainders A1--B2 of the four-$\rho$ terms)'),
        (r'is replaced by $\mathcal D$ (three-$\rho$ note).', r'is replaced by $\mathcal D$ (Part~\ref{!part:three}).'),
        (r'in the $U(w)\mathbb A^{(3)}$ term (three-$\rho$ note).', r'in the $U(w)\mathbb A^{(3)}$ term (Part~\ref{!part:three}).'),
        (r'so either form can be used.', r'so either form can be used. F4 is the $p^+$ ordering of the virtual factor (Sec.~\ref{!sec:alpha}).'),
        (r'Secs.~4--8 of \texttt{ThreeRho\_to\_TwoRho.pdf}', r'Apps.~\ref{!R:sec:fourrho}--\ref{!R:sec:summary}'),
        (r'The three-$\rho$ note changed two Group I rows', r'Part~\ref{!part:three} changed two Group I rows'),
        (r'Eq.~(3.4) of the previous note', r'Eq.~\eqref{!R:eq:plain}'),
        (r'from the three-$\rho$ note]', r'from App.~\ref{!R:sec:fourrho}]'),
        (r'This supersedes the remark in \texttt{ThreeRho\_cancellation} (Sec.~on the first correction)',
         r'This supersedes the remark at the end of correction~1 in Sec.~\ref{!C:sec:corrections}'),
        (r'Eq.~(7.3) of the three-$\rho$ note', r'Eq.~\eqref{!C:eq:comm}'),
    ], 'T')
    return prefix(t, 'T')


def rel_T():
    S = sections(body(PATHS['T']))
    t = S['app:rel']
    t = t[t.index('\n', t.index(r'\label{app:rel}')) + 1:]
    return prefix(t, 'T')


# ------------------------------------------------------------------------------------------ Part IV
GLP = ['A,B', '4rho', 'I.II', 'I.III', 'I.IV', 'I.V', 'II.I', 'II.IV', 'III.I', 'III.II', 'III.III', 'III.V', 'III.VI']


def glabelP(g):
    return {'A,B': 'A1--B2', '4rho': r'\texttt{4rho}'}.get(g, g)


def tablesP(R, target):
    out = []
    for n, (key, d) in enumerate(R.items()):
        out.append(r'\subsection*{Structure %d: $%s$}' % (n + 1, MPK.kname(key)))
        out.append(r'\addcontentsline{toc}{subsection}{Structure %d}' % (n + 1))
        if target:
            if not any('-BFKL' in r['contrib'] for r in d['table']):
                out.append(r'BFKL$\otimes$LO has no term with this kernel: the terms cancel among themselves.')
            out.append(r'\begin{longtable}{@{}>{\raggedright\arraybackslash$}p{0.38\textwidth}<{$}>{\raggedright\arraybackslash}p{0.40\textwidth}>{$}c<{$}>{$}c<{$}@{}}')
            out.append(r'\text{colour structure} & contributions (row: coefficient) & \text{sum} & \text{BFKL}\\ \hline\endhead')
        else:
            out.append(r'\begin{longtable}{@{}>{\raggedright\arraybackslash$}p{0.42\textwidth}<{$}>{\raggedright\arraybackslash}p{0.44\textwidth}>{$}c<{$}@{}}')
            out.append(r'\text{colour structure} & contributions (row: coefficient) & \text{sum}\\ \hline\endhead')
        for row in d['table']:
            if not row['contrib']:
                continue
            notes = [(g, v) for g, v in row['contrib'].items() if g != '-BFKL']
            notes.sort(key=lambda gv: GLP.index(gv[0]))
            s = sp.nsimplify(sp.expand(sum(v for _, v in notes)))
            items = r';\ '.join(r'%s:\,$%s$' % (glabelP(g), MPK.ctex(v, plus=False)) for g, v in notes) or r'--'
            if target:
                b = -row['contrib'].get('-BFKL', 0)
                ok = r'\ \text{\ok}' if sp.simplify(s - b) == 0 else r'\ \text{\bad}'
                out.append(r'%s & %s & %s & %s%s\tabularnewline[2pt]' % (MPK.colour(row['col']), items, MPK.ctex(s, plus=False) if notes else '0',
                                                                      MPK.ctex(b, plus=False), ok))
            else:
                ok = r'\ \text{\ok}' if s == 0 else r'\ \text{\bad}'
                out.append(r'%s & %s & %s%s\tabularnewline[2pt]' % (MPK.colour(row['col']), items, MPK.ctex(s, plus=False), ok))
        out.append(r'\end{longtable}\addtocounter{table}{-1}')
    return '\n'.join(out)


P3_TEXT = r"""\section{The symmetrized three-\texorpdfstring{$\rho$}{rho} terms of \texorpdfstring{$\log\vee$}{log vee} cancel}\label{sec:three}
All symmetrized three-$\rho$ terms of the coefficient of $\log\vee$, after the fixes, sorted by kernel structure and colour network. Every entry is renamed so that the charge tied to $w$ is $x$, the one tied to $w'$ is $x'$, and the third is $y$; if several namings are equivalent, the entry is averaged over them. Contributions are labelled by row: A1--B2 are the three-$\rho$ remainders of the four-$\rho$ terms (for them $P=A$), the other labels are the rows of Sec.~\ref{sec:rows} (their region-$P$ sources, App.~\ref{!app:reorderP}; for I.II the term $-\mathfrak F_\alpha$ of fix F4; for I.III the corrected kernel of P1). Coefficients are in units of $\frac{1}{(2\pi)^3}\frac{g^4}{16\pi^5k^+}\log\vee$ and include the c.c.\ where the row has one. In every table the sum of each row is zero: there are @@P3NC@@ kernel structures and @@P3NET@@ colour networks.

@@P3TABLES@@
"""

P2_TEXT = r"""\section{The two-\texorpdfstring{$\rho$}{rho} terms of \texorpdfstring{$\log\vee$}{log vee}: comparison with BFKL\texorpdfstring{$\otimes$}{x}LO}\label{sec:two}
All two-$\rho$ terms of the coefficient of $\log\vee$, after the fixes:
\begin{itemize}
\item the two-$\rho$ terms of the four-$\rho$ terms: those of the remainders A1--B2 (App.~\ref{!R:sec:fourrho}) and the lists of \texttt{4rho.tex} with F1, F2 (label \texttt{4rho}); for all of them $P=A$;
\item the two-$\rho$ terms of the region-$P$ three-$\rho$ sources (App.~\ref{!app:reorderP}), with P1 in Row III;
\item Group I, Row II: the two-$\rho$ terms of $-\mathfrak F_\alpha$ in the $U(w)\mathbb A^{(3)}$ sub-terms (fix F4, Sec.~\ref{!sec:alpha});
\item the genuine two-$\rho$ terms with $P\neq0$: Group II Rows I and IV, and Group III Row III part 2 (with F5) (Sec.~\ref{sec:rows}).
\end{itemize}
They are grouped by kernel structure and colour network. For each structure the colour networks are reduced to an independent basis, with coefficients that are exact polynomials in $N_c$ (checked with SU(3) and SU(4); the relations are in App.~\ref{!app:rel}). The last column is BFKL$\otimes$LO, Eq.~\eqref{eq:Bop}, in the same basis, and a check mark means that the sum equals it. Units as in Sec.~\ref{sec:three}. The genuine Group II Row IV term is written as a $z$ integrand, $K(w,w',z)-\frac12K(w,w,z)-\frac12K(w',w',z)$.

@@P2TABLES@@
"""


def part_B():
    T = open(os.path.join(BF, 'pk_static.tex')).read()
    S = sections(T[T.index(r'\begin{document}'):T.rindex(r'\end{document}')])
    rows = open(os.path.join(BF, 'rows.tex')).read()
    rows = rep(rows, [
        (r'As in the previous note (fix F6)', r'As in Part~\ref{!part:jimwlk} (fix F6)'),
        (r'the position-space kernel of the three-$\rho$ note (App.~A there).', r'the coordinate-space kernel of Part~\ref{!part:ft}, Eq.~\eqref{!eq:Fab}.'),
        (r'In the previous note the $\Phi$ terms', r'In Part~\ref{!part:three} the $\Phi$ terms'),
        (r'minus Eq.~(4.6) of the previous note', r'minus Eq.~\eqref{!T:eq:rowV}'),
    ], 'rows')
    i0 = rows.index(r'\textbf{Fix F4 does not change $B$.}')
    i1 = rows.index('\n', i0)
    rows = rows[:i0] + (r'\textbf{Fix F4 does not change $B$.} F4 removes $\mathfrak F_\alpha$ from the $\log(\vee/\Lambda)$ terms of the '
                        r'$U(w)\mathbb A^{(3)}$ sub-terms ($A\to A-\alpha$) and leaves the $\log(\vee/k^+)$ terms as they are, so '
                        r'$P=0\to-\alpha$: with the $p^+$-ordered virtual factor, $\mathfrak F_\alpha$ multiplies $-\log(\vee/k^+)$ '
                        r'(Sec.~\ref{!sec:alpha}, Eq.~\eqref{!eq:F4meaning}). The region-$P$ comparison (Sec.~\ref{sec:checks}) confirms '
                        r'$P_{\rm I.II}=-\alpha$ independently. The $\log(\vee/k^+)$ terms of Row II are correct as written, provided the bracket '
                        r'contains $-\bk^i\log(\vee/k^+)$ (Eq.~\eqref{!eq:Bfull}).') + rows[i1:]
    ch = json.load(open(os.path.join(BF, 'checks.json')))
    st3 = json.load(open(os.path.join(BF, 'states3.json')))
    sing = json.load(open(os.path.join(BF, 'bfkl_singlet.json')))
    rowsP = []
    for r in ch['rowsP']:
        w, c = r['written'], r['corrected']
        rowsP.append(r'%s & %d & %d / %d & %d / %d (max $%s$)\\' % (r['row'], r['rho'], w[1], w[0], c[1], c[0], MPK.sci(c[2])))
    st2 = ch['states2']
    names = ['as written', '+ P1 (Group I Row III: sign of $K_h$)', r'+ F3 carried over (Row III part 2 $\to\mathcal D$)',
             '+ F5 carried over (Group III Row III part 2: sign)']
    states = [r'%s & %d / %d & %d / %d\\' % (nm, v2[1], v2[0], v3[1], v3[0]) for nm, (n2, v2), (n3, v3) in zip(names, st2, st3)]
    seeds = 'At two other random configurations (seeds %s) the corrected two-$\\rho$ total again agrees in all structures (largest difference $%s$).' % (
        ' and '.join(str(s) for s, _ in ch['seeds']), MPK.sci(max(v[2] for _, v in ch['seeds'])))
    gu = ch['generalU']
    sg = [r'%s & %d & $%.8f$ & $%.8f$\\' % ('symmetric' if r['sym'] else 'general', r['seed'], r['op'], r['eq']) for r in sing]
    RP2 = pickle.load(open(os.path.join(HERE, 'reducedP2.pkl'), 'rb'))
    RP3 = pickle.load(open(os.path.join(HERE, 'reducedP3.pkl'), 'rb'))
    for R in (RP2, RP3):
        assert all(r['total'] == 0 for d in R.values() for r in d['table']), 'P tables do not close'
    p3 = P3_TEXT.replace('@@P3TABLES@@', tablesP(RP3, False)).replace('@@P3NC@@', str(len(RP3))).replace(
        '@@P3NET@@', str(sum(len(d['nets']) for d in RP3.values())))
    p2 = P2_TEXT.replace('@@P2TABLES@@', tablesP(RP2, True))
    bf = S['sec:bfkl']
    bf = rep(bf, [
        (r'In the basis of Sec.~\ref{sec:two}:', r'In the basis of Sec.~\ref{sec:two} (the same structures and colour basis as the comparison there):'),
        (r'(unlike Appendix F; see the previous note)', r'(unlike Appendix F; see Part~\ref{!part:jimwlk})'),
        (r'Adding the coefficients of $\log(\vee/\Lambda)$ (JIMWLK$\otimes$LO, previous note) and of $\log(\vee/k^+)$ (this note) gives BFKL$\otimes$LO, Eq.~\eqref{eq:Bop}: the projectile evolves with $\log\vee$, and the target with $\log(k^+/\Lambda)$.',
         r'The coefficient of $\log(\vee/\Lambda)$ is JIMWLK$\otimes$LO (Part~\ref{!part:jimwlk}); adding that of $\log(\vee/k^+)$ gives BFKL$\otimes$LO, Eq.~\eqref{eq:Bop} (Sec.~\ref{sec:two}): the projectile evolves with $\log\vee$, and the target with $\log(k^+/\Lambda)$.'),
    ], 'B:bfkl')
    bf = bf.replace('@@BDECOMP@@', MPK.bdecomp(RP2)).replace('@@SINGLET@@', '\n'.join(sg))
    chk = S['sec:checks']
    chk = rep(chk, [(r'i.e.\ fix F4 of the previous note.', r'i.e.\ fix F4 (Sec.~\ref{!sec:alpha}).')], 'B:checks')
    chk = chk.replace('@@ROWSP@@', '\n'.join(rowsP)).replace('@@STATES@@', '\n'.join(states)).replace('@@SEEDS@@', seeds).replace(
        '@@GENU@@', r'%d of %d structures differ' % (gu[1], gu[0]))
    which = rep(S['sec:which'], [(r'The second equation is the result of \texttt{TwoRho\_evolution}.',
                                 r'The second equation is the result of Part~\ref{!part:jimwlk}.')], 'B:which')
    t = '\n'.join([which, S['sec:rows'].replace('@@ROWS@@', rows), S['sec:four'].replace('@@N4@@', '12'), p3, p2, bf, chk])
    assert '@@' not in t, t[t.index('@@') - 80:t.index('@@') + 80]
    reorder = rep(S['sec:reorder'], [
        (r'are worked out step by step in \texttt{ThreeRho\_to\_TwoRho}, Secs.~4--7 (with the corrections of \texttt{ThreeRho\_cancellation}: $\mathcal D$, and the $\alpha$ part of Row II).',
         r'are worked out step by step in Apps.~\ref{!R:sec:fourrho}--\ref{!R:sec:G3} (with the corrections of Part~\ref{!part:three}: $\mathcal D$, and the $\alpha$ part of Row II).'),
        (r'Group I Rows II and V have $P=0$, and their $B=-A$.',
         r'As written, Group I Rows II and V have $P=0$, and their $B=-A$; with F4, Row II has $P=-\alpha$ (Sec.~\ref{!sec:alpha}).'),
    ], 'B:reorder')
    reorder = re.sub(r'^\\section\{.*?\}\\label\{sec:reorder\}', r'\\section{Region $p^+\\gg k^+$: reordering the three-$\\rho$ terms}\\label{sec:reorder}',
                     reorder, count=1, flags=re.S)
    reorder = reorder.replace('@@PSECTIONS@@', MPK.psections())
    rel = MPK.relations(RP2) + '\n' + MPK.relations(RP3)
    return prefix(t, 'B'), prefix(reorder, 'B'), rel


# ------------------------------------------------------------------------------------------ appendices from ThreeRho_to_TwoRho
def apps_R():
    S = sections(body(PATHS['R']))
    keep = ['sec:identity', 'sec:fourrho', 'sec:G1', 'sec:G2', 'sec:G3', 'sec:summary']
    t = '\n'.join(S[k] for k in keep)
    t = rep(t, [(r'We do not list them again, and we do not check their cancellation in this note.',
                 r'We do not list them again; their cancellation is checked in Part~\ref{!part:three}.')], 'R')
    t = t.replace(r'\label{sec:identity}', r'\label{sec:identity}\label{!app:reorder}', 1)
    chk = retitle(S['app:checks'], 'Checks of the reordering')
    return prefix(t, 'R'), prefix(chk, 'R')


SCRIPTS = r"""\section{Scripts}\label{app:scripts}
The numbers and tables come from the scripts of the folders \texttt{three\_rho}, \texttt{two\_rho}, \texttt{bfkl} (see their README files) and of this folder:
\begin{itemize}
\item \texttt{ft/ft\_checks.py}: Part~\ref{part:ft} (partial fractions, the split of the bracket, the $\xi$ integral of \texttt{main.tex}, the commutator identity with explicit charges);
\item \texttt{allsets.py}: all terms by row, for the coefficients of $\log(\vee/\Lambda)$ ($A$) and of $\log\vee$ ($P$), as written and corrected;
\item \texttt{check\_all.py}: pointwise totals of $A$ vs JIMWLK$\otimes$LO and of $P$ vs BFKL$\otimes$LO (two-$\rho$; symmetrized three-$\rho$);
\item \texttt{build\_P.py} and \texttt{reduce\_run.py}: the exact tables of Secs.~\ref{B:sec:three}, \ref{B:sec:two};
\item \texttt{make\_combined.py}: writes this file from \texttt{combined\_static.tex}, \texttt{ft\_part.tex} and the notes of the other folders.
\end{itemize}
"""


def main():
    T = open(os.path.join(HERE, 'combined_static.tex')).read()
    st = json.load(open(os.path.join(HERE, 'allstates.json')))
    partB, reorderP, relP = part_B()
    rapps, rchk = apps_R()
    RP2 = pickle.load(open(os.path.join(HERE, 'reducedP2.pkl'), 'rb'))
    T = T.replace('@@P2W@@', str(st['P2_written'][0])).replace('@@P2N@@', str(st['P2_written'][1]))
    T = T.replace('@@P3W@@', str(st['P3_written'][0])).replace('@@P3N@@', str(st['P3_written'][1]))
    T = T.replace('@@P2NC@@', str(len(RP2)))
    rowsA = [('$\\log(\\vee/\\Lambda)$', 'two-$\\rho$ vs JIMWLK$\\otimes$LO', 'A2'), ('$\\log(\\vee/\\Lambda)$', 'symmetrized three-$\\rho$ (sum $=0$)', 'A3'),
             ('$\\log\\vee$', 'two-$\\rho$ vs BFKL$\\otimes$LO', 'P2'), ('$\\log\\vee$', 'symmetrized three-$\\rho$ (sum $=0$)', 'P3')]
    T = T.replace('@@STATESALL@@', '\n'.join(r'%s & %s & %d / %d & %d / %d (max $%s$)\\' % (
        a, b, st[k + '_written'][0], st[k + '_written'][1], st[k + '_corr'][0], st[k + '_corr'][1], MPK.sci(st[k + '_corr'][2])) for a, b, k in rowsA))
    rel = (r'\section{Relations between colour networks}\label{app:rel}' + '\n'
           r'For the kernel structures in which the simplified colour structures are linearly dependent, each dependent one equals the '
           r'following combination of basis structures (exact; checked for SU(3) and SU(4)).' + '\n'
           r'\subsection{Part~\ref{part:jimwlk}: coefficient of $\log(\vee/\Lambda)$, two-$\rho$}' + '\n' + rel_T() + '\n'
           r'\subsection{Part~\ref{part:bfkl}: coefficient of $\log\vee$}' + '\n' + relP)
    rep_ = {
        '%%@@FTPART@@': ft_part(),
        '%%@@PART2@@': part_C(),
        '%%@@PART3@@': part_T(),
        '%%@@PART4@@': partB,
        '%%@@APPS@@': '\n'.join([rapps, reorderP, rel, rchk, SCRIPTS]),
    }
    for k, v in rep_.items():
        assert k in T, k
        T = T.replace(k, v)
    assert '@@' not in T, T[T.index('@@') - 80:T.index('@@') + 80]
    for bad in ('previous note', 'this note', r'three-$\rho$ note', r'\texttt{ThreeRho\_to\_TwoRho.pdf}'):
        for m in re.finditer(re.escape(bad), T):
            print('check wording:', T[max(0, m.start() - 100):m.end() + 60].replace('\n', ' '))
    open(os.path.join(HERE, 'RowII_FT_and_evolution.tex'), 'w').write(T)
    print('written RowII_FT_and_evolution.tex', len(T))


if __name__ == '__main__':
    main()
