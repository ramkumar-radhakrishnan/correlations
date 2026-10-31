r"""Writes, from your BFKL.tex and Other_terms.tex (input/):

  BFKL_rows.tex       all terms proportional to log(vee/k+), row by row: yours (black), the ones still inside
                      Other_terms.tex and the ones that come out of line 1 by exact identities (blue), the line that is not
                      in the rows (red); the identities; three-rho cancellation; two-rho sum = BFKL x LO; poles; your files.
  Remaining_rows.tex  Other_terms.tex without them, Lambda -> 0 (blue: changed or added, red: removed).

Only vee > k+, p+ > Lambda are used; Lambda -> 0 where possible, vee exact.  Numbers from check_rows.json,
exact_check.json, three_bfkl.json, pole_bfkl.json, compare_4rho.json, compare_final.json."""
import json
import os
import re
import rows as R
import check_rows as CR
import exact_check as EC
import assemble as AS

HERE = R.HERE
UB, UO, EX = CR.UB, CR.UO, CR.EX
HB = CR.HB
J = lambda f: json.load(open(os.path.join(HERE, f)))
CRJ, T3J, PBJ, C4J, CFJ, ECJ = (J('check_rows.json'), J('three_bfkl.json'), J('pole_bfkl.json'), J('compare_4rho.json'),
                                J('compare_final.json'), J('exact_check.json'))
FEC = open(os.path.join(HERE, '..', 'final_check', 'Final_evolution_check.tex')).read()


def tex_num(x):
    if x == 0:
        return '0'
    m, e = ('%.1e' % x).split('e')
    return r'%s\times10^{%d}' % (m, int(e))


PREAMBLE = r'''\documentclass[10pt]{article}
\usepackage[a4paper,margin=0.3in]{geometry}
\usepackage{etoolbox}
\AtBeginEnvironment{align}{\small}
\usepackage{amsmath,amssymb,amsfonts,mathtools}
\usepackage[dvipsnames]{xcolor}
\usepackage{array,booktabs,float,longtable}
\usepackage[colorlinks=true,linkcolor=blue]{hyperref}
\allowdisplaybreaks
\numberwithin{equation}{section}
\setlength{\jot}{3pt}
\setlength{\emergencystretch}{3em}
\definecolor{JIM}{RGB}{0,130,60}
\definecolor{CC}{RGB}{140,0,170}
\definecolor{BFK}{RGB}{215,95,0}
\newcommand{\new}[1]{{\color{blue}#1}}
\newcommand{\old}[1]{{\color{red}#1}}
\newcommand{\nc}{\color{blue}}
\newcommand{\rc}{\color{red}}
\newcommand{\NNv}{\mathcal N\log\frac{\vee}{k^+}}
\def\mf#1{\mathbb{#1}}
'''

GROUPS = [
    ('Group I (with c.c.)', 'gI', [
        ('Row I', [], 'four'),
        ('Row II', ['G1 II(1)', 'G1 II(2a)', 'G1 II(2b)'], None),
        ('Row III', ['G1 III N1', 'G1 III N2', 'G1 III N3a', 'G1 III N3b'], None),
        ('Row IV', ['G1 IV'], None),
        ('Row V', ['G1 V(1)', 'G1 V(2a)', 'G1 V(2b UV)', 'G1 V(2b fin)'], None)]),
    ('Group II (no c.c.)', 'gII', [
        ('Row I', ['G2 I(1) A', 'G2 I(1) B', 'G2 I(2a)', 'G2 I(2b)', 'G2 I(2c)', 'G2 I(2d)', 'G2 I(2e)', 'G2 I(2f)'], None),
        ('Row II', [], 'four'), ('Row III', [], 'four'),
        ('Row IV', ['G2 IV'], None)]),
    ('Group III (with c.c.)', 'gIII', [
        ('Row I', ['G3 I'], None), ('Row II', ['G3 II'], None),
        ('Row III', ['G3 III(1)', 'G3 III(2)'], None),
        ('Row IV', [], 'four'),
        ('Row V', ['G3 V'], None), ('Row VI', ['G3 VI'], None)])]

SUB = {'G1 II(1)': 'first part', 'G1 II(2a)': r'second part, $d^3N_{2a}$', 'G1 II(2b)': r'second part, $d^3N_{2b}$',
       'G1 III N1': r'$d^3N_1$', 'G1 III N2': r'$d^3N_2$', 'G1 III N3a': r'$d^3N_{3a}$', 'G1 III N3b': r'$d^3N_{3b}$',
       'G1 IV': '', 'G1 V(1)': 'first part', 'G1 V(2a)': r'second part, $d^3N_{2a}$',
       'G1 V(2b UV)': r'second part, $d^3N_{2b\,\rm UV}$', 'G1 V(2b fin)': r'second part, $d^3N_{2b\,\rm finite}$',
       'G2 I(1) A': 'first part, first line', 'G2 I(1) B': 'first part, second line',
       'G2 I(2a)': r'second part, $d^3N_{2a}$', 'G2 I(2b)': r'second part, $d^3N_{2b}$', 'G2 I(2c)': r'second part, $d^3N_{2c}$',
       'G2 I(2d)': r'second part, $d^3N_{2d}$', 'G2 I(2e)': r'second part, $d^3N_{2e}$', 'G2 I(2f)': r'second part, $d^3N_{2f}$',
       'G2 IV': '', 'G3 I': '', 'G3 II': '', 'G3 III(1)': 'first part', 'G3 III(2)': 'second part', 'G3 V': '', 'G3 VI': ''}


def c5(name, t):
    if name != 'G3 III(2)':
        return t
    for m in (r"\int_{w,w',y',z,y,x}", r"\int_{w,y',z,y,x}"):
        t = t.replace(m, r"\old{%s}\,\new{\int_{w,x',y',z,x}}" % m)
    return t


def lab(kind, name):
    return 'eq:%s:%s' % (kind, re.sub(r'[^A-Za-z0-9]', '', name))


def align(B, items, label):
    if not items:
        return ''
    s = AS.assemble(B, items)
    s = re.sub(r'\n\s*\n', '\n', s)
    return s.replace('\n\\end{align}', r' \label{%s}' % label + '\n\\end{align}')


def your_items(name):
    if name not in UB:
        return []
    return [(t, 'red' if (name, n) in CR.REMOVED else None) for n, t in enumerate(UB[name].terms)]


def extracted_items(name):
    return [(c5(name, s['B']), 'blue') for s in EX.get(name, []) if s.get('B') is not None]


def line1_items(name):
    """Group III Rows III, V, VI: from line 1 by the exact identities (I2) and (I3)"""
    if name not in CR.G3_LINE1:
        return []
    return [(c5(name, EC.negate(EX[name][1]['B'])), 'blue'), (HB[name], 'blue')]


def remaining_items(name):
    items = []
    for n, s in enumerate(EX[name]):
        if n == 0 and name in CR.G3_LINE1:
            items += [(c5(name, t), st) for t, st in EC.g3_display(name)]
        elif n == 0 and name == 'G2 IV':
            items += EC.g2iv_display()
        elif s.get('F') is None:
            items.append((UO[name].terms[n], 'red'))
        elif (name, n) == ('G1 III N3b', 3):
            items.append((s['F'], 'red'))
        else:
            items.append((c5(name, s['F']), None))
    return items


# ------------------------------------------------------------------------------------------------ notes
NB = {
    'G1 II(2a)': r'''As you have it.''',
    'G1 III N1': r'''Your line is correct. \new{Added}: lines 1, 2 and 4 of $d^3N_1$ in \texttt{Other\_terms.tex} still contain $\log\frac{\vee}{k^+}$: $\log(\frac{\vee}{k^+}-1)=\log\frac{\vee}{k^+}+\log(1-\frac{k^+}{\vee})$, $\log\frac{\vee(x-z)^2-k^+(x-w)^2}{(\Lambda+k^+)(x-z)^2-k^+(x-w)^2}=\log\frac{\vee}{k^+}+\dots$, $\log\frac{k^+(x-w)^2-\vee(x-z)^2}{k^+(x-w)^2-\Lambda(x-z)^2}=\log\frac{\vee}{k^+}+\dots$.''',
    'G1 III N2': r'''Your line is correct. \new{Added}: the $\log\frac{\vee}{k^+}$ parts of all eight lines of $d^3N_2$ in \texttt{Other\_terms.tex}.''',
    'G1 IV': r'''Your line is correct. \new{Added}: the $\log\frac{\vee}{k^+}$ parts of $\log\frac{(\vee-k^+)(x-z)^2+k^+(x-w)^2}{\Lambda(x-z)^2+k^+(x-w)^2}$ and $\log(\frac{\vee}{k^+}-1)$ (both lines in \texttt{Other\_terms.tex}).''',
    'G1 V(2b UV)': r'''\textbf{Remove (\old{red}).} Row V has $0<p^+<k^+$: no $\vee$ and no $\log\frac{\vee}{k^+}$. Before ``Setting $\Lambda\to0$'' your bracket is $2\log\frac{\vee}{\Lambda}-2\log\frac{\vee}{k^+}+2\log(1-\frac{\Lambda}{k^+})-\frac{11}{6}=2\log\frac{k^+}{\Lambda}+2\log(1-\frac{\Lambda}{k^+})-\frac{11}{6}$. The $2\log\frac{k^+}{\Lambda}$ is in your JIMWLK terms, and $2\log(1-\frac{\Lambda}{k^+})-\frac{11}{6}$ is in \texttt{Other\_terms.tex}. This line comes from the bracket $\big[\log\frac{\vee}{k^+}+\log\frac{k^+}{\Lambda}\big]-2\log\frac{\vee}{k^+}$ after ``Setting $\Lambda\to0$'', which lacks a factor 2 on the first bracket (E1).''',
    'G2 I(1) A': r'''Not in your \texttt{BFKL.tex}. \new{Added}: the $\log\frac{\vee}{k^+}$ parts of all three lines (\texttt{Other\_terms.tex}).''',
    'G2 I(1) B': r'''\new{Added}: the $\log\frac{\vee}{k^+}$ parts of all three lines.''',
    'G2 I(2b)': r'''Not in your \texttt{BFKL.tex}; all three lines are \new{added}.''',
    'G2 I(2c)': r'''Not in your \texttt{BFKL.tex}; both lines are \new{added}.''',
    'G2 I(2d)': r'''Your line is correct. \new{Added}: the two lines of \texttt{Other\_terms.tex}.''',
    'G2 I(2e)': r'''Your line is correct. \new{Added}: the two lines of \texttt{Other\_terms.tex}.''',
    'G2 I(2f)': r'''Your line is correct. \new{Added}: lines 1 and 3 of \texttt{Other\_terms.tex}.''',
    'G2 IV': r'''\textbf{Both of your lines are correct.} They come out of line 1 of \texttt{Other\_terms.tex} by exact identities (Sec.~\ref{sec:exact}).
\begin{itemize}
\item Line 1 is the $\frac{d\xi}{\xi}$ part of the subtraction in $\frac{\xi}{(1-\xi)_+}$, $-g(1)\log\frac{\vee}{k^+}$, Eq.~\eqref{eq:I2xi}.
\item Line 2, $H_4$, is the $\log(y'-y)^2$ term of $\frac{1-\xi}{\xi}$ at $y'=y$, Eq.~\eqref{eq:I4}.
\end{itemize}
\new{Added}: the $\log\frac{\vee}{k^+}$ part of $\log(\frac{\vee}{k^+}-1)$ in line 2 of \texttt{Other\_terms.tex}. It is exactly minus your line 1, poles included, so the BFKL term of this row is $H_4$. Line 1 of \texttt{Other\_terms.tex} is rewritten to match (\texttt{Remaining\_rows.pdf}).''',
    'G3 I': r'''Your line is correct. \new{Added}: the $\log\frac{\vee}{k^+}$ parts of both lines of \texttt{Other\_terms.tex}.''',
    'G3 II': r'''Your line is correct. \new{Added}: the $\log\frac{\vee}{k^+}$ parts of both lines of \texttt{Other\_terms.tex}.''',
    'G3 III(1)': r'''Not in your \texttt{BFKL.tex}. \new{Added}:
\begin{itemize}
\item[(i)] the $\log\frac{\vee}{k^+}$ part of $\log(\frac{\vee}{k^+}-1)=\log\frac{\vee}{k^+}+\log(1-\frac{k^+}{\vee})$ in line 2;
\item[(ii)] from line 1, the subtraction of $\left[\frac{1}{p^+}\right]_+$: $-f(0)\log\frac{\vee}{k^+}$, Eq.~\eqref{eq:I2}. It is exactly minus (i);
\item[(iii)] from line 1, the term (B): the $-\frac{1}{k^+}\delta_{im}\delta_{jk'}$ term with $G(y')\to G(z)$, Eq.~\eqref{eq:I3}.
\end{itemize}
(i) and (ii) cancel, so the BFKL term of this part is (B).''',
    'G3 III(2)': r'''\new{Added}: (i) and (ii) as in the first part (they cancel), and (iii) (B) from the $-\frac{p^+}{k^{+2}}$ term, Eqs.~\eqref{eq:I3}, \eqref{eq:I5}.''',
    'G3 V': r'''\new{Added}: (i), (ii), (iii) as in Row III: (i) and (ii) cancel, and (B) comes from the $-\frac{1}{k^+}\delta_{im}\delta_{jk'}$ term.''',
    'G3 VI': r'''\new{Added}: (i), (ii), (iii) as in Row III: (i) and (ii) cancel, and (B) comes from the $-\frac{1}{k^+}\delta_{im}\delta_{jk'}$ term.''',
}
NONE_B = {'G1 II(1)': r'''No BFKL term. The only $\log\frac{\vee}{k^+}$ is in $\log(1-\frac{k^+}{\vee})\log\frac{\vee}{k^+}$, whose coefficient depends on $\vee$. For every $\vee>k^+$ it is bounded, $0\le-\log(1-\frac{k^+}{\vee})\log\frac{\vee}{k^+}\le\log^22$, with the maximum at $\vee=2k^+$. It stays in the remaining terms, as in your \texttt{Other\_terms.tex}.''',
          'G1 II(2b)': r'No BFKL term ($R^i$ depends on $\vee$ only through $\frac{k^+}{\vee}$).',
          'G1 III N3a': 'No BFKL term.', 'G1 III N3b': 'No BFKL term.',
          'G1 V(1)': 'No BFKL term.', 'G1 V(2a)': 'No BFKL term.', 'G1 V(2b fin)': 'No BFKL term.',
          'G2 I(2a)': r'No BFKL term: $\vee$ cancels between numerator and denominator of the first logarithm.'}

NR = {
    'G1 II(1)': r'''Unchanged. $\log(1-\frac{k^+}{\vee})\log\frac{\vee}{k^+}$ stays here: its coefficient depends on $\vee$, and it is $\le\log^22$ for all $\vee>k^+$. The pole bracket times $\frac{11}{3}$ is the $\beta_0$ term.''',
    'G1 II(2b)': r'''Unchanged; $R^i$ as in your file.''',
    'G1 III N1': r'''\new{Blue}: what is left of the logarithms after the $\log\frac{\vee}{k^+}$ parts are taken out, with $\Lambda\to0$.''',
    'G1 III N2': r'''\new{Blue}: as in $d^3N_1$.''',
    'G1 III N3a': r'''Only $\Lambda\to0$ (\new{blue}); $\mathcal A_0$ is $\mathcal A$ at $\Lambda=0$, Eq.~\eqref{eq:A0}.''',
    'G1 III N3b': r'''$\Lambda\to0$. \textbf{C4}: the last line (\old{red}) is, for $\Lambda\to0$, $\log\frac{(z-w)^2}{(x-w)^2}$ times the edge kernel $\frac{(x'-w')\cdot(y-z)}{(x'-w')^2(y-z)^2(z-w)^2}$: not integrable at $z\to w$. Its log-partner was replaced by $\mathcal D$; its finite part belongs to the $\mathcal D$ region (open).''',
    'G1 IV': r'''\new{Blue}: the remainders of the two logarithms.''',
    'G1 V(1)': 'Unchanged.', 'G1 V(2a)': 'Unchanged.',
    'G1 V(2b UV)': r'''$2\log(1-\frac{\Lambda}{k^+})\to0$; $-\frac{11}{6}$ is the $\beta_0$ term that cancels the $\frac{11}{3}$ of Row II part 1.''',
    'G1 V(2b fin)': r'''$\Lambda\to0$ in the limits of the plus-prescribed integrals.''',
    'G2 I(1) A': r'''\new{Blue}: remainders. Line 2 (\old{red}) is entirely BFKL: $\log\frac{\vee}{\Lambda+k^+}=\log\frac{\vee}{k^+}-\log(1+\frac{\Lambda}{k^+})\to\log\frac{\vee}{k^+}$.''',
    'G2 I(1) B': r'''As the first line.''',
    'G2 I(2a)': r'''$\Lambda\to0$ in the second logarithm.''',
    'G2 I(2b)': r'''\new{Blue}: remainders.''', 'G2 I(2c)': r'''\new{Blue}: remainders.''',
    'G2 I(2d)': r'''\new{Blue}: remainders.''', 'G2 I(2e)': r'''\new{Blue}: remainders.''',
    'G2 I(2f)': r'''\new{Blue}: remainders; $\frac{1}{k^++\Lambda}\to\frac{1}{k^+}$.''',
    'G2 IV': r'''Line 1 is rewritten exactly (\texttt{BFKL\_rows.pdf}, Sec.~2):
\begin{itemize}
\item $\frac{\xi}{(1-\xi)_+}\to\xi+\frac{\xi^2}{(1-\xi)_+}$. The difference, $-g(1)\log\frac{\vee}{k^+}$, is your \texttt{BFKL.tex} line 1.
\item The first \new{blue} line subtracts the $\log(y'-y)^2$ term of $\frac{1-\xi}{\xi}$ at $y'=y$. After the $y'$ integral, that term is $\big[\log\frac{\vee}{k^+}-(1-\frac{k^+}{\vee})\big]H_4\big|_{\log\to1}$. Its $\log\frac{\vee}{k^+}$ part is your $H_4$.
\item The second \new{blue} line is the rest, $-(1-\frac{k^+}{\vee})H_4\big|_{\log\to1}$.
\end{itemize}
Line 2: $\log(\frac{\vee}{k^+}-1)\to\log(1-\frac{k^+}{\vee})$. In line 3, $\frac{k^+}{k^++\Lambda}\to1$.''',
    'G3 I': r'''\new{Blue}: remainders.''', 'G3 II': r'''\new{Blue}: remainders.''',
    'G3 III(1)': r'''Line 1 is rewritten exactly (\texttt{BFKL\_rows.pdf}, Sec.~2):
\begin{itemize}
\item $\left[\frac{1}{p^+}\right]_+\to\left[\frac{k^+}{p^+(p^++k^+)}\right]_++\frac{1}{p^++k^+}$. The difference, $-f(0)\log\frac{\vee}{k^+}$, is in the BFKL list.
\item In the $-\frac{1}{k^+}\delta_{im}\delta_{jk'}$ term, $G(y')\to G(y')-G(z)$ (\new{blue} line). The $G(z)$ part is (B).
\end{itemize}
Line 2: $\log(\frac{\vee}{k^+}-1)\to\log(1-\frac{k^+}{\vee})$.''',
    'G3 III(2)': r'''As the first part, with the $-\frac{p^+}{k^{+2}}$ term in place of $-\frac{1}{k^+}\delta_{im}\delta_{jk'}$. The last \new{blue} line is the part of its $G(z)$ term that is not $\log\frac{\vee}{k^+}$, Eq.~(I5). \textbf{C5}: the measure is $\int_{w,x',y',z,x}$.''',
    'G3 V': r'''As Row III.''', 'G3 VI': r'''As Row III.''',
}


def rows_tex(kind):
    out = []
    for gname, gl, rows in GROUPS:
        out.append(r'\section{%s}\label{sec:%s%s}' % (gname, kind, gl))
        for rname, names, special in rows:
            out.append(r'\subsection{%s}' % rname)
            if special == 'four':
                out.append(r'Only four $\rho$\textquotesingle s (Sec.~\ref{sec:four}).' if kind == 'B'
                           else r'Only four $\rho$\textquotesingle s: no remaining terms.')
                continue
            for name in names:
                sub = SUB[name]
                out.append(r'\paragraph{%s%s.}' % (rname, (', ' + sub) if sub else ''))
                if kind == 'B':
                    yi, xi, li = your_items(name), extracted_items(name), line1_items(name)
                    if not yi and not xi:
                        out.append(NONE_B.get(name, 'No BFKL term.'))
                        continue
                    if NB.get(name):
                        out.append(NB[name])
                    if yi:
                        out.append(r'\noindent Your \texttt{BFKL.tex}:')
                        out.append(align(UB[name], yi, lab('Byour', name)))
                    if xi:
                        out.append(r'\noindent \new{From \texttt{Other\_terms.tex}}%s:' % (' (i)' if li else ''))
                        out.append(align(UO[name], xi, lab('Bnew', name)))
                    if li:
                        out.append(r'\noindent \new{From line 1, by the exact identities} (ii), (iii):')
                        out.append(align(UO[name], li, lab('Bline1', name)))
                else:
                    if name not in UO:
                        out.append(r'Nothing is left: this part has only $\log\frac{k^+}{\Lambda}$ and $\log\frac{\vee}{k^+}$ terms.')
                        continue
                    if NR.get(name):
                        out.append(NR[name])
                    out.append(align(UO[name], remaining_items(name), lab('F', name)))
                    if name == 'G1 II(2b)':
                        out.append('where\n\\begin{align}' + R.align_text('Other_terms', 50, 55) + '\\end{align}')
    return '\n'.join(out)


# ------------------------------------------------------------------------------------------------ pieces
def user_align(name):
    a, b = R.B_RANGES[name]
    return '\\begin{align}' + R.align_text('BFKL', a, b) + '\\end{align}'


def comment_section():
    L = open(os.path.join(HERE, 'input', 'BFKL.tex')).read()
    i = L.index(r'\section{BFKL contribution}') + len(r'\section{BFKL contribution}')
    j = L.index(r'\end{comment}')
    return L[i:j]


def three_table(key):
    d = T3J[key]['seed11']
    rows, tot = [], 0.0
    for labl in T3J[key]['terms']:
        v = d[labl]
        tot += v
        rows.append(r'%s & $%+.4f$\\' % (labl.replace('#', r'line~').replace('_', r'\_'), v))
    return rows, tot


def sec(t, a, b):
    i = t.index(a)
    return t[i:t.index(b, i)]


BF = sec(FEC, r'\section{BFKL: the coefficient', r'\section{Every term that depends')
BF = BF[BF.index(r'\paragraph{Group I, Row II'):]
BF = BF.replace(r'\bigg|^{\log\vee}', r'\bigg|^{\log\frac{\vee}{k^+}}')
BF = BF.replace(r' the $\alpha$ part, fix F4', '').replace(r'; with the sign of $p^+>k^+$, fix P1', '')
BF = BF.replace(r'\paragraph{The sum.}', r'\subsection{The sum}')
BF = BF.replace(r'\eqref{eq:BFKL}', r'\eqref{eq:BFKLapp}')
BF = BF.replace(r'Eq.~\eqref{eq:BFKLapp} before the charge average', r'your App.~E formula, Eq.~\eqref{eq:BFKLapp}, before the charge average')

pb = [v for k, v in PBJ.items() if k != 'pole_items']
WF, WG = T3J['final']['worst'], T3J['as given']['worst']
NF, NG = len(T3J['final']['terms']), len(T3J['as given']['terms'])
A_ok = max(v['A'] for v in CRJ.values())
A_given = CRJ['G1 V(2b UV)']['A_given']
B_ok = max(v['B'] for v in CRJ.values())
C_ok = max(v['C'] for v in CRJ.values())
C4 = max(v['relative'] for v in C4J.values())
CF = max(max(v) for v in CFJ.values())
EC_SEEDS = [ECJ[k] for k in ECJ if k != 'fourier']
EC_ID = max(r[n]['identity'] for r in EC_SEEDS for n in r)
EC_B = max(r[n]['GzFT_vs_B'] for r in EC_SEEDS for n in r if n != 'G2 IV')
EC_T0 = max(r['G2 IV']['term0_plus_line2'] for r in EC_SEEDS)
FOU = ECJ['fourier']
PL = max(p['list_pole_relative'] for p in pb)
PLF = max(p['list_finite_relative'] for p in pb)


def fourier_line():
    j1 = [v for k, v in FOU.items() if k.startswith('int J1')]
    lg = [v for k, v in FOU.items() if k.startswith('FT log')]
    return (r'$\int_0^\infty J_1(qr)\,dr\cdot q$ with regulator $e^{-ar}$, $a=0.1,0.03,0.01$: %s; '
            r'$2\pi\int_0^\infty rJ_0(qr)\log r^2dr\big/(-\frac{4\pi}{q^2})$ with $e^{-ar^2}$: %s (both $\to1$ as $a\to0$; $q=0.7,1.9$).'
            % ('; '.join(', '.join('%.4f' % x for x in v) for v in j1), '; '.join(', '.join('%.4f' % x for x in v) for v in lg)))


BFKL = PREAMBLE + r'''
\title{BFKL evolution of the projectile:\\ all terms $\propto\log\frac{\vee}{k^+}$ from your \texttt{BFKL.tex} and \texttt{Other\_terms.tex},\\ with $\vee$ exact (exact identities only)}
\date{}
\begin{document}
\maketitle
\noindent Notation as in your files; $\mathcal N=\frac{1}{(2\pi)^3}\frac{g^4}{16\pi^5}\frac{1}{k^+}$, $K(a,b,z)=\mathcal K(a-z)\cdot\mathcal K(b-z)$, $U^{bc}=U^{cb}$. Groups I and III are taken with their c.c., Group II without, everywhere. In the rows: black = your \texttt{BFKL.tex}, \new{blue} = added, \old{red} = to be removed. The remaining terms are in \texttt{Remaining\_rows.pdf}. The same analysis from your Group I, II, III and \texttt{4rho.tex} directly: \texttt{full\_check/Checks\_JIMWLK\_BFKL.pdf} and \texttt{full\_check/Remaining\_RC\_DGLAP\_finite.pdf}.

\section{Answer}\label{sec:answer}
\paragraph{BFKL is there.} Your cross sections give
\begin{equation}
\frac{d^3N_{\rm BFKL}}{d^3k}=\mathcal N\log\frac{\vee}{k^+}\;{\rm BFKL}\otimes{\rm LO},
\end{equation}
your App.~E formula, Eq.~\eqref{eq:BFKLfinal}. Only $\vee>k^+$ and $p^+>\Lambda$ are used, with $\Lambda\to0$ wherever possible and $\vee$ exact. Every step is an exact identity, Sec.~\ref{sec:exact}; nothing is expanded or dropped, and what is not $\propto\log\frac{\vee}{k^+}$ stays, exactly, in the remaining terms.

\paragraph{Your Group I, II, III files} are correct (Sec.~\ref{sec:groups}). The one change that matters is E1, in Group I Row V 2b UV. The others are typos and measures.

\paragraph{Where the BFKL terms are.}
\begin{itemize}
\item In explicit logarithms of \texttt{Other\_terms.tex}, e.g.\ $\log(\frac{\vee}{k^+}-1)=\log\frac{\vee}{k^+}+\log(1-\frac{k^+}{\vee})$.
\item In line 1 of Group III Rows III, V, VI and Group II Row IV, the $p^+$ (resp.\ $\xi$) integral. Three exact identities bring them out:
\begin{itemize}
\item the subtraction of the plus prescription gives $-f(0)\log\frac{\vee}{k^+}$, which cancels the $\log\frac{\vee}{k^+}$ of line 2 exactly (I2);
\item the 2-D Fourier transform of the $y'$ integral at $y'=z$ gives the terms (B), with $\int_0^{\vee-k^+}\frac{dp^+}{p^++k^+}=\log\frac{\vee}{k^+}$ (I3);
\item the Fourier transform of $\log(y'-y)^2$ gives $H_4$ (I4).
\end{itemize}
\end{itemize}
The test that the split is complete: the remaining terms stay finite for $\vee\to\infty$. This fixes the coefficient of $\log\frac{\vee}{k^+}$ uniquely.

\paragraph{Changes to your files.}
\begin{itemize}
\item \texttt{BFKL.tex}: remove only Group I Row V 2b UV (\old{red}, E1). Your two lines of Group II Row IV are correct: they are (I2) and (I4).
\item Add (\new{blue}): the $\log\frac{\vee}{k^+}$ parts still in \texttt{Other\_terms.tex}; and, for Group III Rows III, V, VI, the terms $-f(0)\log\frac{\vee}{k^+}$ and (B) from line 1.
\item \texttt{Other\_terms.tex}: line 1 of Group III Rows III, V, VI and Group II Row IV in its exact rewritten form, and $\log(1-\frac{k^+}{\vee})$ in line 2 (\texttt{Remaining\_rows.pdf}).
\end{itemize}

\paragraph{Result.}
\begin{itemize}
\item Row by row, BFKL list $+$ remaining terms $=$ the row minus its JIMWLK part. The pointwise terms agree to $@@AOK@@$ ($O(\Lambda)$ at $\Lambda=10^{-9}$). The rewriting of line 1 is checked numerically, with the $p^+$, $\xi$ integrals done, to $@@ECID@@$ (Sec.~\ref{sec:exact}).
\item The @@NF@@ three-$\rho$ terms \textbf{cancel}, to $@@WF@@$ (Sec.~\ref{sec:three}). Your \texttt{BFKL.tex} alone leaves $@@WG@@$.
\item The list equals, row by row, the list whose two-$\rho$ terms add up to \textbf{BFKL$\otimes$LO} exactly (to $@@BOK@@$; 9 kernel structures, 109 colour structures, 0 differences; Sec.~\ref{sec:two}).
\item The list has no pole. Your Group II Row IV line 1 and the $\log\frac{\vee}{k^+}$ part of line 2 cancel, poles included ($@@PL@@$).
\item The result is your commented ``BFKL contribution'' to $@@CF@@$. In colour-singlet form it is your App.~E formula (Sec.~\ref{sec:final}).
\end{itemize}

\paragraph{Correction of my previous version.} There I split only the explicit logarithms and kept every $p^+$, $\xi$ integral untouched. Then the $\log\frac{\vee}{k^+}$ inside line 1 stays in the remaining terms, which grow with $\vee$. That split is not complete, and that is why BFKL$\otimes$LO did not come out.

(B) and $H_4$ are not approximations. They follow from the exact identities (I2)--(I5). Earlier I described them through the region of large $p^+$; that was a description, not a step of the calculation.

\section{The exact identities}\label{sec:exact}
Nothing is expanded: each line below is an identity, and every remainder is kept in \texttt{Remaining\_rows.pdf}.

\paragraph{(I1) Explicit logarithms.} $\log(\frac{\vee}{k^+}-1)=\log\frac{\vee}{k^+}+\log(1-\frac{k^+}{\vee})$, $\log\frac{(\vee-k^+)A+k^+B}{\Lambda A+k^+B}=\log\frac{\vee}{k^+}+\log\frac{(1-\frac{k^+}{\vee})A+\frac{k^+}{\vee}B}{B}+O(\Lambda)$, \dots\ (the table is in \texttt{LL\_terms.pdf}, Sec.~2).

\paragraph{(I2) The subtraction of the plus prescription.} With your $\left[\frac{1}{p^+}\right]_+f=\frac{f(p^+)-f(0)}{p^+}$, define $\left[\frac{k^+}{p^+(p^++k^+)}\right]_+f\equiv\frac{k^+[f(p^+)-f(0)]}{p^+(p^++k^+)}$. Since $\frac{1}{p^+}=\frac{1}{p^++k^+}+\frac{k^+}{p^+(p^++k^+)}$ and $\int_0^{\vee-k^+}\frac{dp^+}{p^++k^+}=\log\frac{\vee}{k^+}$,
\begin{equation}
\int_0^{\vee-k^+}dp^+\left[\frac{1}{p^+}\right]_+f(p^+)=\int_0^{\vee-k^+}dp^+\left(\left[\frac{k^+}{p^+(p^++k^+)}\right]_++\frac{1}{p^++k^+}\right)f(p^+)-f(0)\log\frac{\vee}{k^+}.
\label{eq:I2}
\end{equation}
For $\xi$, with your $\frac{d\xi}{\xi^2}\frac{\xi}{(1-\xi)_+}g=\frac{d\xi}{\xi(1-\xi)}[g(\xi)-g(1)]$, define $\frac{d\xi}{\xi^2}\frac{\xi^2}{(1-\xi)_+}g\equiv\frac{d\xi}{1-\xi}[g(\xi)-g(1)]$. Then $\frac{1}{\xi(1-\xi)}=\frac{1}{\xi}+\frac{1}{1-\xi}$ gives
\begin{equation}
\int_{\frac{k^+}{\vee}}^1\frac{d\xi}{\xi^2}\frac{\xi}{(1-\xi)_+}g=\int_{\frac{k^+}{\vee}}^1\frac{d\xi}{\xi^2}\Big(\xi+\frac{\xi^2}{(1-\xi)_+}\Big)g-g(1)\log\frac{\vee}{k^+}.
\label{eq:I2xi}
\end{equation}
Your line 2 is $f(0)\big[\log\frac{k^+}{\Lambda}+\log(\frac{\vee}{k^+}-1)\big]$ (resp.\ $g(1)[\dots]$). It came from this subtraction, so its $\log\frac{\vee}{k^+}$ and the $-f(0)\log\frac{\vee}{k^+}$ above cancel exactly, and $\log(1-\frac{k^+}{\vee})$ is left. In Group II Row IV, $-g(1)\log\frac{\vee}{k^+}$ is your \texttt{BFKL.tex} line 1.

\paragraph{(I3) The $y'$ integral at $y'=z$.} In line 1 of Group III Rows III(1), V, VI, take the $-\frac{1}{k^+}\delta_{im}\delta_{jk'}$ term; in Row III(2), the $-\frac{p^+}{k^{+2}}$ term. Its coefficient tends to a constant at large $p^+$. Write
\begin{equation*}
G^{k'}(y')=\frac{(x'-y')^{k'}}{(x'-y')^2}\times\big[\text{colour bracket with }U(y')\big]=G^{k'}(z)+\big[G^{k'}(y')-G^{k'}(z)\big].
\end{equation*}
With $\int_u e^{-iq\cdot u}\frac{u^m}{u^2}=-2\pi i\frac{q^m}{q^2}$ and $q=k\frac{p^++k^+}{k^+}$, exactly,
\begin{equation}
\int_{y'}e^{-ik\cdot(y'-w)}e^{-ik\cdot\frac{p^+}{k^+}(y'-z)}\frac{(y'-z)^m}{(y'-z)^2}\,G^{k'}(z)=\frac{k^+}{p^++k^+}\int_{w'}e^{-ik\cdot(w'-w)}\frac{(w'-z)^m}{(w'-z)^2}\,G^{k'}(z).
\label{eq:I3}
\end{equation}
For Rows III(1), V, VI: $\int_0^{\vee-k^+}\frac{dp^+}{2\pi}\big(-\frac{1}{k^+}\big)\frac{k^+}{p^++k^+}=-\frac{1}{2\pi}\log\frac{\vee}{k^+}$, exactly. This is (B), with nothing left over.

\paragraph{(I5) Row III(2).} With $A=(x-z)^2$, $B=(x-w)^2$,
\begin{equation}
-\frac{p^+}{k^{+2}}\,\frac{k^+}{p^++k^+}\,\frac{1}{p^+A+k^+B}=-\frac{1}{k^+A\,(p^++k^+)}+\frac{B}{A\,(p^++k^+)(p^+A+k^+B)}.
\label{eq:I5}
\end{equation}
The first term gives (B). The second falls like $\frac{1}{p^{+2}}$ and stays in the remaining terms.

\paragraph{(I4) Group II Row IV.} In the $\frac{1-\xi}{\xi}$ term of line 1, write $H(y,y')=H(y,y)+[H(y,y')-H(y,y)]$, with $H(y,y')=\frac{(x-y)\cdot(x'-y')}{(x-y)^2(x'-y')^2}[U^{e'c}(y')-U^{e'c}(x')][U^{ec}(y)-U^{ce}(x)]\rho^{e'}(x')\rho^e(x)$. The constant part of the pole bracket gives $\delta^{(2)}(k/\xi)=0$ for $k\neq0$, and $\int_{y'}e^{-ik\cdot(y'-y)\frac{1}{\xi}}\log\frac{(y'-y)^2\mu^2}{4}=-\frac{4\pi\xi^2}{k^2}$. Also $\int_{\frac{k^+}{\vee}}^1d\xi\,\frac{1-\xi}{\xi}=\log\frac{\vee}{k^+}-1+\frac{k^+}{\vee}$ and $\int_{w,w'}e^{-ik\cdot(w'-w)}\mathcal K(w-z)\cdot\mathcal K(w'-z)=\frac{4\pi^2}{k^2}$. So the $H(y,y)$ part is, exactly,
\begin{equation}
\frac{g^2}{4\pi^2}\frac{N_c}{(2\pi)^3}\frac{g^2}{2\pi^2}\frac{1}{k^+}\,\frac{d_\perp}{2}\,\frac{4\pi}{k^2}\Big[\log\frac{\vee}{k^+}-1+\frac{k^+}{\vee}\Big]\int_{x,x',z}H(z,z)
=\Big[\log\frac{\vee}{k^+}-1+\frac{k^+}{\vee}\Big]\,H_4\Big|_{\log\frac{\vee}{k^+}\to1},
\label{eq:I4}
\end{equation}
at $d_\perp=2$ (no pole multiplies it). Its $\log\frac{\vee}{k^+}$ term is your $H_4$, Group II Row IV line 2.

\paragraph{The split is complete.} After these steps every remaining $p^+$ ($\xi$) integrand falls at least like $\frac{1}{p^{+2}}$ after the transverse integrals:
\begin{itemize}
\item $G(y')-G(z)$ and $H(y,y')-H(y,y)$ vanish at the coincident point;
\item the $\frac{1}{p^++k^+}$ and $\left[\frac{k^+}{p^+(p^++k^+)}\right]_+$ terms carry the $\frac{k^+}{p^++k^+}$ of (I3);
\item the remaining terms are convergent oscillating integrals.
\end{itemize}
So the remaining terms stay finite for $\vee\to\infty$, and the coefficient of $\log\frac{\vee}{k^+}$ is unique. $\vee$ stays exact in the remaining terms.

\paragraph{Numerical checks} (\texttt{exact\_check.py}). Two configurations, $\vee=20k^+$, the $p^+$ and $\xi$ integrals done numerically, pointwise in the transverse points:
\begin{itemize}
\item line 1 $=$ (rewritten line 1) $+$ ($G(y')-G(z)$ term) $+$ ($-f(0)\log\frac{\vee}{k^+}$) $+$ ($G(z)$ term), for Group III Rows III(1), III(2), V, VI; and line 1 $=$ (rewritten line 1) $+$ (your line 1) for Group II Row IV: $\le@@ECID@@$;
\item the $G(z)$ term in the form \eqref{eq:I3}, $p^+$ integral done, $=$ (B) ($+$ the second term of \eqref{eq:I5} in Row III(2)): $\le@@ECB@@$;
\item your Group II Row IV line 1 $+$ the $\log\frac{\vee}{k^+}$ part of line 2 $=0$: $@@ECT0@@$;
\item the Fourier transforms: @@FOU@@
\end{itemize}

@@ROWS@@

\section{Four \texorpdfstring{$\rho$}{rho}'s}\label{sec:four}
Your \texttt{BFKL.tex}, unchanged: the three-$\rho$ remainders $d^3N_a$, $d^3N_b$ (they enter Sec.~\ref{sec:three}) and the two-$\rho$ lists $d^3N_1$, $d^3N_2$. Numerically, $d^3N_1+d^3N_2$ equals the four-$\rho$ two-$\rho$ terms that entered the checked BFKL sum (Eqs.~\eqref{eq:P12}, \eqref{eq:P13} below) to $@@C4@@$ (charges symmetrized, with c.c.).
@@FOURA@@
@@FOURB@@
@@FOUR1@@
@@FOUR2@@

\section{The three-\texorpdfstring{$\rho$}{rho} terms cancel}\label{sec:three}
Same test as in \texttt{LL\_terms.pdf}, Sec.~7. Every three-$\rho$ term of the list (your lines, the blue lines, (B), and your $d^3N_a$, $d^3N_b$) is
\begin{itemize}
\item renamed to the frame $w,w',z$, charges $x,y,x'$;
\item summed over the 6 placements of the charges and, for Groups I and III, its c.c.;
\item evaluated with random points, random SU(3) adjoint Wilson lines ($U^{bc}=U^{cb}$), random charges, and $\log\frac{\vee}{k^+}=1$.
\end{itemize}
Your Group I Row II line 2 has the $z$ integral done; it is written back as a $z$ integral. Values at one configuration, in units of $i\,\mathcal N\,e^{-ik\cdot(w'-w)}$:
\begin{longtable}{@{}lr@{}}\toprule
term & value\\\midrule
@@TT@@
\midrule
total & $@@TTOT@@$\\\bottomrule
\end{longtable}
\noindent Largest relative residual over three configurations: $@@WF@@$. Your \texttt{BFKL.tex} alone (without the blue terms): $@@WG@@$.

\section{The two-\texorpdfstring{$\rho$}{rho} terms: BFKL}\label{sec:two}
Row by row, the list above equals the $\log\frac{\vee}{k^+}$ list of \texttt{LL\_terms.pdf} ($@@BOK@@$; check B). In each row of Group III Rows III, V, VI and Group II Row IV, the pair (i), (ii) cancels; what is left there is (B), resp.\ $H_4$. The two-$\rho$ terms of that list were added up and compared with BFKL$\otimes$LO kernel by kernel. Row by row they are the genuine two-$\rho$ terms and the commutator terms of the three-$\rho$ terms, $[\rho^a(u),\rho^b(v)]=if^{abc}\rho^c(u)\delta(u-v)$ (orange / purple: 168 / 162):
@@BF@@

\section{The result}\label{sec:final}
Your commented ``BFKL contribution'' (both forms) is the sum, Eq.~\eqref{eq:Bsum}, to $@@CF@@$ (charges symmetrized, with c.c.):
@@COMMENT@@
This is your App.~E formula,
\begin{align}
\frac{dN}{d^2w\,d^2w'\,dy}\bigg|_{\rm BFKL}={}&-\frac{\alpha_s^2N_c\,\delta Y}{8\pi^6(N_c^2-1)}\int_{x,y,z}\langle\rho^f(x)\rho^f(y)\rangle\Big[
\mathcal K(x-w)\cdot\mathcal K(y-w')K_D(x,y,z){\rm Tr}\{(U^\dagger(x)-U^\dagger(w))(U(y)-U(w'))\}\notag\\
&-2\,\mathcal K(z-w)\cdot\mathcal K(z-w')K(x,y,z){\rm Tr}\{(U^\dagger(z)-U^\dagger(w))(U(z)-U(w'))\}\notag\\
&+\mathcal K(z-w)\cdot\mathcal K(y-w')\big(2K(x,y,z)-K(x,x,z)\big){\rm Tr}\{(U^\dagger(z)-U^\dagger(w))(U(y)-U(w'))\}\notag\\
&+\mathcal K(x-w)\cdot\mathcal K(z-w')\big(2K(x,y,z)-K(y,y,z)\big){\rm Tr}\{(U^\dagger(x)-U^\dagger(w))(U(z)-U(w'))\}\Big],
\label{eq:BFKLapp}
\end{align}
$K_D=K(x,x,z)-2K(x,y,z)+K(y,y,z)$, $\delta Y=\log\frac{\vee}{k^+}$, $\mathcal N\delta Y\leftrightarrow\frac{\alpha_s^2\delta Y}{8\pi^6(N_c^2-1)}$. Term by term, your singlet form is this, with $\mathcal K(z-w)=-\mathcal K(w-z)$ and $x\leftrightarrow y$ in the last two lines. It agrees with BFKL$\otimes$LO to $8\times10^{-15}$ (\texttt{Final\_evolution\_check.pdf}). Hence
\begin{equation}
\frac{d^3N_{\rm BFKL}}{d^3k}=\frac{dN}{d^2w\,d^2w'\,dy}\bigg|_{\rm BFKL}\otimes\int_{w,w'}e^{-ik\cdot(w'-w)},\qquad \delta Y=\log\frac{\vee}{k^+}.
\label{eq:BFKLfinal}
\end{equation}

\section{Poles}\label{sec:pole}
With $E(r)=\frac{2}{\epsilon}+2\gamma_E+\log\frac{r^2\mu^2}{4}$ and $\mathcal S$ the LO structure (renaming as in \texttt{Final\_evolution\_check.pdf}, Sec.~6.1):
\begin{itemize}
\item The list has two terms with a pole: your Group II Row IV line 1, $+\frac{g^4N_c}{8\pi^4}E(w-w')\mathcal S\,\frac{1}{(2\pi)^3k^+}\log\frac{\vee}{k^+}$, and the $\log\frac{\vee}{k^+}$ part of line 2 of \texttt{Other\_terms.tex}, which is exactly minus it. Pointwise, their poles cancel (relative $@@PL@@$) and so do their finite parts ($@@PLF@@$).
\item The removed Group I Row V 2b UV line, $-\frac{g^4N_c}{16\pi^4}\big[E(x-w)+E(x'-w')\big]\mathcal S$ (with c.c.) times $\frac{1}{(2\pi)^3k^+}\log\frac{\vee}{k^+}$, has a pole of its own. With it the list would not be pole free.
\item The poles $\times\log\frac{k^+}{\Lambda}$ of Group I Row II, Group I Row V 2b UV (with E1) and Group II Row IV are in the ratio $1:-2:1$ and cancel (\texttt{LL\_terms.pdf}).
\end{itemize}

\section{Your Group I, II, III files}\label{sec:groups}
Every step of your \texttt{Group\_I.tex}, \texttt{Group\_II.tex}, \texttt{Group\_III.tex} was checked numerically (\texttt{rows\_check/Group\_check.pdf}, with corrected copies). Groups I and III are taken with their complex conjugate, Group II without, in every test here (three $\rho$, two $\rho$, poles).
\begin{itemize}
\item \textbf{M1 (changes the result), Group I Row V 2b UV}, Eqs.~(1.26), (1.32): the bracket is $\big[2\log\frac{k^+}{\Lambda}+2\log(1-\frac{\Lambda}{k^+})-\frac{11}{6}\big]$, with no $\log\frac{\vee}{k^+}$ (Row V has $0<p^+<k^+$). Your exact Eq.~(1.25) has $2\log\frac{\vee}{\Lambda}-2\log\frac{\vee}{k^+}=2\log\frac{k^+}{\Lambda}$. This is the red line here (E1).
\item M2, Group II Row IV, Eq.~(1.16): the phase of the second term is $e^{-ik\cdot(y'-y)\frac{p^++k^+}{k^+}}$ (a typo; Eq.~(1.17) is right).
\item M3, Group III Row III, first line: the overall factor is $-\frac{4}{(2\pi)^3}$, as in Rows V and VI.
\item M4, Group III Row III part 2: the measures are $\int_{w,w',x',y',z,x}$ and, after the $w'$ integral, $\int_{w,x',y',z,x}$ (C5).
\item M5, Group I Row II, Eq.~(1.6): the measure is $\int_{w,w',x'}$.
\item M6: typos ($k\{+\}$ in Group II Eq.~(1.13); brace mismatches in Group I; two commas; a repeated sentence in Group II).
\item Every other step agrees to $10^{-16}$--$5\times10^{-12}$; the plus-prescription steps agree up to $O(\Lambda)$.
\end{itemize}
The rows used here are those of your Group files: \texttt{Other\_terms.tex} $+$ \texttt{BFKL.tex} reproduce them (check A). The latest upload differs from the files of the split only in colour marks and in the Group I Row V 2b UV line (M1).

\section{Checks}\label{sec:checks}
Scripts in \texttt{bfkl\_rows/}:
\begin{center}\small
\begin{tabular}{@{}p{0.62\textwidth}p{0.32\textwidth}@{}}\toprule
check & result\\\midrule
A. \texttt{BFKL.tex} $+$ \texttt{Other\_terms.tex} $=$ row $-$ JIMWLK, row by row, pointwise (\texttt{check\_rows.py}) & $\le @@AOK@@$ ($=O(\Lambda)$, $\Lambda=10^{-9}$)\\
\quad with the red line (G1 V 2b UV) & $@@AGIVEN@@$\\
\quad terms with an unevaluated integral: identical text & all\\
line 1 of G3 III, V, VI, G2 IV $=$ rewritten line 1 $+$ the BFKL terms (\texttt{exact\_check.py}) & $\le@@ECID@@$; (B): $\le@@ECB@@$\\
B. BFKL list (black $+$ blue) $=$ the list of \texttt{LL\_terms.pdf} (two-$\rho$ sum BFKL$\otimes$LO) & $\le @@BOK@@$\\
C. remaining terms $=$ the rest of the rows, $\Lambda\to0$ & $\le @@COK@@$\\
three-$\rho$ terms cancel: list / your \texttt{BFKL.tex} alone (\texttt{three\_bfkl.py}) & $@@WF@@$ / $@@WG@@$\\
four-$\rho$ two-$\rho$ lists vs the checked ones (\texttt{compare\_4rho.py}) & $@@C4@@$\\
your ``BFKL contribution'' vs the checked sum (\texttt{compare\_final.py}) & $@@CF@@$\\
poles of the list cancel / finite parts (\texttt{pole\_bfkl.py}) & $@@PL@@$ / $@@PLF@@$\\
two-$\rho$ sum vs BFKL$\otimes$LO (\texttt{../final\_check/render.py}) & 0 of 109 differ\\\bottomrule
\end{tabular}
\end{center}
\end{document}
'''

rows_t3, tot_t3 = three_table('final')
REP = {'@@AOK@@': tex_num(A_ok), '@@AGIVEN@@': '%.2f' % A_given, '@@BOK@@': tex_num(B_ok), '@@COK@@': tex_num(C_ok),
       '@@ECID@@': tex_num(EC_ID), '@@ECB@@': tex_num(EC_B), '@@ECT0@@': tex_num(EC_T0), '@@FOU@@': fourier_line(),
       '@@PL@@': tex_num(PL), '@@PLF@@': tex_num(PLF),
       '@@NF@@': str(NF), '@@NG@@': str(NG), '@@WF@@': tex_num(WF), '@@WG@@': '%.2f' % WG,
       '@@C4@@': tex_num(C4), '@@CF@@': tex_num(CF), '@@TTOT@@': '%+.1e' % tot_t3}
out = BFKL
for k, v in REP.items():
    out = out.replace(k, v)
out = out.replace('@@ROWS@@', rows_tex('B')).replace('@@TT@@', '\n'.join(rows_t3)).replace('@@BF@@', BF)
out = out.replace('@@COMMENT@@', comment_section())
for k, n in (('@@FOURA@@', '4rho a'), ('@@FOURB@@', '4rho b'), ('@@FOUR1@@', '4rho 1'), ('@@FOUR2@@', '4rho 2')):
    out = out.replace(k, user_align(n))
assert '@@' not in out

# ------------------------------------------------------------------------------------------------ remaining terms
A0 = (r"\mathcal{A}_{0} & = \frac{1}{\sqrt{(x-w)^{2}\,(x-z)^{2}-\left[(x-w)\cdot (x-z)\right]^{2}}}\, \Bigg\{"
      r"\mathrm{arc\,tan}\,\left[\frac{(x-z)^{2}-(x-w)\cdot (x-z)}{\sqrt{(x-w)^{2}\,(x-z)^{2}-\left[(x-w)\cdot (x-z)\right]^{2}}}\right]"
      r"\notag \\ & + \mathrm{arc\,tan}\,\left[\frac{(x-w)\cdot (x-z)}{\sqrt{(x-w)^{2}\,(x-z)^{2}-\left[(x-w)\cdot (x-z)\right]^{2}}}\right] \Bigg\}")
DEFS = r'''\paragraph{Notation.} Your plus prescriptions, and two more:
\begin{align*}
\int dp^+\left[\frac{1}{p^+}\right]_+f&=\int dp^+\,\frac{f(p^+)-f(0)}{p^+},&
\int dp^+\left[\frac{k^+}{p^+(p^++k^+)}\right]_+f&=\int dp^+\,\frac{k^+\,[f(p^+)-f(0)]}{p^+(p^++k^+)},\\
\int\frac{d\xi}{\xi^2}\,\frac{\xi}{(1-\xi)_+}\,g&=\int\frac{d\xi}{\xi(1-\xi)}\,[g(\xi)-g(1)],&
\int\frac{d\xi}{\xi^2}\,\frac{\xi^2}{(1-\xi)_+}\,g&=\int\frac{d\xi}{1-\xi}\,[g(\xi)-g(1)].
\end{align*}'''
REM = PREAMBLE + r'''
\title{The remaining terms: your \texttt{Other\_terms.tex} without the BFKL terms, $\Lambda\to0$}
\date{}
\begin{document}
\maketitle
\noindent Companion of \texttt{BFKL\_rows.pdf}. Your \texttt{Other\_terms.tex}, row by row: \new{blue} = changed or added, \old{red} = removed. Only $\vee>k^+$ and $p^+>\Lambda$ are used; $\Lambda\to0$ wherever possible, $\vee$ exact; every change is an exact identity (\texttt{BFKL\_rows.pdf}, Sec.~2).

\section{What changed}\label{sec:what}
\begin{itemize}
\item Every explicit logarithm that contains $\log\frac{\vee}{k^+}$ is replaced by its remainder, with $\Lambda\to0$ (its $\log\frac{\vee}{k^+}$ part is in \texttt{BFKL\_rows.pdf}):
\begin{equation*}
\log\Big(\frac{\vee}{k^+}-1\Big)\to\log\Big(1-\frac{k^+}{\vee}\Big),\qquad
\log\frac{(\vee-k^+)A+k^+B}{\Lambda A+k^+B}\to\log\frac{(1-\frac{k^+}{\vee})A+\frac{k^+}{\vee}B}{B},\qquad \log\frac{\vee}{\Lambda+k^+}\to0,\ \dots
\end{equation*}
\item Line 1 of Group III Rows III, V, VI and Group II Row IV is rewritten exactly.
\begin{itemize}
\item $\left[\frac{1}{p^+}\right]_+\to\left[\frac{k^+}{p^+(p^++k^+)}\right]_++\frac{1}{p^++k^+}$, and $\frac{\xi}{(1-\xi)_+}\to\xi+\frac{\xi^2}{(1-\xi)_+}$. The difference, $-f(0)\log\frac{\vee}{k^+}$ (resp.\ your Group II Row IV line 1), is in the BFKL list; it cancels the $\log\frac{\vee}{k^+}$ of line 2.
\item In the term whose coefficient tends to a constant at large $p^+$, $G(y')\to G(y')-G(z)$. The $G(z)$ part is (B).
\item In Group II Row IV, the $\log(y'-y)^2$ term of $\frac{1-\xi}{\xi}$ at $y'=y$ is subtracted, and its non-logarithmic part $-(1-\frac{k^+}{\vee})H_4|_{\log\to1}$ is added. Its $\log\frac{\vee}{k^+}$ part is your $H_4$.
\end{itemize}
\item $\Lambda\to0$ everywhere else:
\begin{itemize}
\item logarithms without a $\log\frac{\vee}{k^+}$ part;
\item $\mathcal A\to\mathcal A_0$;
\item the limits of the plus-prescribed integrals;
\item $\frac{1}{k^++\Lambda}\to\frac{1}{k^+}$, $\frac{k^+}{k^++\Lambda}\to1$, $2\log(1-\frac{\Lambda}{k^+})\to0$.
\end{itemize}
\item Removed (\old{red}):
\begin{itemize}
\item Group II Row I line 2 (both blocks), which is entirely BFKL;
\item (C4) the edge-kernel line of Group I Row III $d^3N_{3b}$.
\end{itemize}
\item C5: the measure of Group III Row III part 2.
\end{itemize}
Row by row, \texttt{BFKL\_rows.pdf} $+$ this file $=$ every row minus its JIMWLK part. The pointwise terms agree to $@@COK@@$; the rewriting of line 1 is checked in \texttt{exact\_check.py}. These terms stay finite for $\vee\to\infty$ and $\Lambda\to0$, apart from the open finite part of the $\mathcal D$ region (C4).
@@DEFS@@
\begin{align}
@@A0@@
\label{eq:A0}
\end{align}

@@ROWS@@
\end{document}
'''
REM = REM.replace('@@COK@@', tex_num(C_ok)).replace('@@A0@@', A0).replace('@@DEFS@@', DEFS).replace('@@ROWS@@', rows_tex('F'))
if __name__ == '__main__':
    open(os.path.join(HERE, 'BFKL_rows.tex'), 'w').write(out)
    open(os.path.join(HERE, 'Remaining_rows.tex'), 'w').write(REM)
    print('written BFKL_rows.tex, Remaining_rows.tex')
