r"""Writes, from your BFKL.tex and Other_terms.tex (input/):

  BFKL_rows.tex       all terms proportional to log(vee/k+): yours (black), the ones still inside Other_terms.tex
                      (blue), the two that are not BFKL terms (red); three-rho cancellation; two-rho sum; final form.
  Remaining_rows.tex  Other_terms.tex without them, Lambda -> 0 (blue: changed, red: removed).

Numbers from check_rows.json, three_bfkl.json, pole_bfkl.json, compare_4rho.json, compare_final.json."""
import json
import os
import re
import rows as R
import extract as X
import check_rows as CR
import assemble as AS
import gen as GN

HERE = R.HERE
UB, UO, EX = CR.UB, CR.UO, CR.EX
HB = CR.HB
J = lambda f: json.load(open(os.path.join(HERE, f)))
CRJ, T3J, PBJ, C4J, CFJ = J('check_rows.json'), J('three_bfkl.json'), J('pole_bfkl.json'), J('compare_4rho.json'), J('compare_final.json')
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
    ('Group II', 'gII', [
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


def neg(t):
    t = t.lstrip()
    assert t.startswith('-'), t[:30]
    return '+' + t[1:]


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
    return [(s['B'], 'blue') for s in EX.get(name, []) if s.get('B') is not None]


def remaining_items(name):
    items = []
    for n, s in enumerate(EX[name]):
        if s.get('F') is None:
            items.append((UO[name].terms[n], 'red'))
        elif (name, n) == ('G1 III N3b', 3):
            items.append((s['F'], 'red'))
        else:
            items.append((c5(name, s['F']), None))
    if name in HB:
        items.append((neg(HB[name]), 'blue'))
    if name == 'G2 IV':
        items.append(('-' + re.sub(r'^\s*&?\s*\+?\s*', '', UB['G2 IV'].terms[1]), 'blue'))
    return items


# ------------------------------------------------------------------------------------------------ notes
NB = {
    'G1 II(2a)': r'''As you have it.''',
    'G1 III N1': r'''Your line is correct. \new{Added}: lines 1, 2 and 4 of $d^3N_1$ in \texttt{Other\_terms.tex} still contain $\log\frac{\vee}{k^+}$: $\log(\frac{\vee}{k^+}-1)=\log\frac{\vee}{k^+}+\log(1-\frac{k^+}{\vee})$, $\log\frac{\vee(x-z)^2-k^+(x-w)^2}{(\Lambda+k^+)(x-z)^2-k^+(x-w)^2}=\log\frac{\vee}{k^+}+\dots$, $\log\frac{k^+(x-w)^2-\vee(x-z)^2}{k^+(x-w)^2-\Lambda(x-z)^2}=\log\frac{\vee}{k^+}+\dots$.''',
    'G1 III N2': r'''Your line is correct. \new{Added}: the $\log\frac{\vee}{k^+}$ parts of all eight lines of $d^3N_2$ in \texttt{Other\_terms.tex}.''',
    'G1 IV': r'''Your line is correct. \new{Added}: the $\log\frac{\vee}{k^+}$ parts of $\log\frac{(\vee-k^+)(x-z)^2+k^+(x-w)^2}{\Lambda(x-z)^2+k^+(x-w)^2}$ and $\log(\frac{\vee}{k^+}-1)$ (both lines in \texttt{Other\_terms.tex}).''',
    'G1 V(2b UV)': r'''\textbf{Remove (\old{red}).} Row V has no $\log\frac{\vee}{k^+}$. Before ``Setting $\Lambda\to0$'' your bracket is $2\log\frac{\vee}{\Lambda}-2\log\frac{\vee}{k^+}+2\log(1-\frac{\Lambda}{k^+})-\frac{11}{6}=2\log\frac{k^+}{\Lambda}+2\log(1-\frac{\Lambda}{k^+})-\frac{11}{6}$. The $2\log\frac{k^+}{\Lambda}$ is in your JIMWLK terms, and $2\log(1-\frac{\Lambda}{k^+})-\frac{11}{6}$ is in \texttt{Other\_terms.tex}. Nothing is left for BFKL; this line comes from the bracket $\big[\log\frac{\vee}{k^+}+\log\frac{k^+}{\Lambda}\big]-2\log\frac{\vee}{k^+}$ after ``Setting $\Lambda\to0$'', which lacks a factor 2 on the first bracket (E1). Its pole cancels the pole of the red line of Group II Row IV, Sec.~\ref{sec:pole}.''',
    'G2 I(1) A': r'''Not in your \texttt{BFKL.tex}. \new{Added}: the $\log\frac{\vee}{k^+}$ parts of all three lines (\texttt{Other\_terms.tex}).''',
    'G2 I(1) B': r'''\new{Added}: the $\log\frac{\vee}{k^+}$ parts of all three lines.''',
    'G2 I(2b)': r'''Not in your \texttt{BFKL.tex}; all three lines are \new{added}.''',
    'G2 I(2c)': r'''Not in your \texttt{BFKL.tex}; both lines are \new{added}.''',
    'G2 I(2d)': r'''Your line is correct. \new{Added}: the two lines of \texttt{Other\_terms.tex}.''',
    'G2 I(2e)': r'''Your line is correct. \new{Added}: the two lines of \texttt{Other\_terms.tex}.''',
    'G2 I(2f)': r'''Your line is correct. \new{Added}: lines 1 and 3 of \texttt{Other\_terms.tex}.''',
    'G2 IV': r'''\textbf{Line 2 ($H_4$) is correct.} It is the BFKL term of this row, contained in the $\xi$ integral of line 1 of \texttt{Other\_terms.tex}; there it must be subtracted (Remaining terms, C1).
\textbf{Line 1: remove (\old{red}).} It is the $-\log\frac{\vee}{k^+}$ hidden in the subtraction of $\frac{\xi}{(1-\xi)_+}$ in line 1 of \texttt{Other\_terms.tex}. That one is cancelled exactly by the $+\log\frac{\vee}{k^+}$ inside $\log(\frac{\vee}{k^+}-1)$ of line 2 of \texttt{Other\_terms.tex}. Line 1 of \texttt{Other\_terms.tex} still contains it, so writing it here as well counts it twice. Together with the red line of Group I Row V it is pole free, but the two leave a term that is not in BFKL$\otimes$LO (Sec.~\ref{sec:pole}).''',
    'G3 I': r'''Your line is correct. \new{Added}: the $\log\frac{\vee}{k^+}$ parts of both lines of \texttt{Other\_terms.tex}.''',
    'G3 II': r'''Your line is correct. \new{Added}: the $\log\frac{\vee}{k^+}$ parts of both lines of \texttt{Other\_terms.tex}.''',
    'G3 III(1)': r'''Not in your \texttt{BFKL.tex}. \new{Added}: the term (B), the $\log\frac{\vee}{k^+}$ hidden in the $-\frac{1}{k^+}$ term of line 1 (region $y'\to z$; \texttt{Finite\_terms\_rows.pdf}, Sec.~2). The $\log(\frac{\vee}{k^+}-1)$ of line 2 is \emph{not} a BFKL term: it cancels the $\log\frac{\vee}{k^+}$ of the subtraction in $\left[\frac{1}{p^+}\right]_+$.''',
    'G3 III(2)': r'''\new{Added}: (B) of this part (from the $-\frac{p^+}{k^{+2}}$ term of line 1).''',
    'G3 V': r'''\new{Added}: (B), as in Row III.''',
    'G3 VI': r'''\new{Added}: (B), as in Row III.''',
}
NONE_B = {'G1 II(1)': r'No BFKL term: $\log(1-\frac{k^+}{\vee})\log\frac{\vee}{k^+}\to0$ for $\vee\to\infty$.',
          'G1 II(2b)': r'No BFKL term ($R^i$ depends on $\vee$ only through $\frac{k^+}{\vee}$).',
          'G1 III N3a': 'No BFKL term.', 'G1 III N3b': 'No BFKL term.',
          'G1 V(1)': 'No BFKL term.', 'G1 V(2a)': 'No BFKL term.', 'G1 V(2b fin)': 'No BFKL term.',
          'G2 I(2a)': r'No BFKL term: $\vee$ cancels between numerator and denominator of the first logarithm.'}

NR = {
    'G1 II(1)': r'''Unchanged. ($\log(1-\frac{k^+}{\vee})\log\frac{\vee}{k^+}\to0$ for $\vee\to\infty$; the pole bracket times $\frac{11}{3}$ is the $\beta_0$ term.)''',
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
    'G2 IV': r'''\textbf{C1}: the \new{blue} line $-\log\frac{\vee}{k^+}H_4$ takes out the BFKL term contained in the $\xi$ integral of line 1 (it is in the BFKL list). Line 2 keeps $\log(\frac{\vee}{k^+}-1)$: it cancels the $-\log\frac{\vee}{k^+}$ of the subtraction in $\frac{\xi}{(1-\xi)_+}$ of line 1. In line 3, $\frac{k^+}{k^++\Lambda}\to1$. Lines 1, 2 and the blue line are the DGLAP term.''',
    'G3 I': r'''\new{Blue}: remainders.''', 'G3 II': r'''\new{Blue}: remainders.''',
    'G3 III(1)': r'''\textbf{C3}: the \new{blue} line is $-$(B), the BFKL term contained in the $-\frac{1}{k^+}$ term of line 1 (it is in the BFKL list). Line 2 keeps $\log(\frac{\vee}{k^+}-1)$. Lines 1, 2 and the blue line together are finite for $\vee\to\infty$ (\texttt{Finite\_terms\_rows.pdf}, Sec.~2).''',
    'G3 III(2)': r'''C3 as in the first part. \textbf{C5}: the measure is $\int_{w,x',y',z,x}$.''',
    'G3 V': r'''C3 as in Row III.''', 'G3 VI': r'''C3 as in Row III.''',
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
                    yi, xi = your_items(name), extracted_items(name)
                    if name in HB:
                        xi = xi + [(HB[name], 'blue')]
                    if not yi and not xi:
                        out.append(NONE_B.get(name, 'No BFKL term.'))
                        continue
                    if NB.get(name):
                        out.append(NB[name])
                    if yi:
                        out.append(r'\noindent Your \texttt{BFKL.tex}:')
                        out.append(align(UB[name], yi, lab('Byour', name)))
                    if xi:
                        out.append(r'\noindent \new{From \texttt{Other\_terms.tex}}:')
                        out.append(align(UO[name], xi, lab('Bnew', name)))
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
BF = BF.replace(r'\label{eq:Bsum}', r'\label{eq:Bsum}').replace(r'\paragraph{The sum.}', r'\subsection{The sum}')
BF = BF.replace(r'\eqref{eq:BFKL}', r'\eqref{eq:BFKLapp}')
BF = BF.replace(r'Eq.~\eqref{eq:BFKLapp} before the charge average', r'your App.~E formula, Eq.~\eqref{eq:BFKLapp}, before the charge average')

pb = list(PBJ.values())
WF, WG = T3J['final']['worst'], T3J['as given']['worst']
NF, NG = len(T3J['final']['terms']), len(T3J['as given']['terms'])
A_ok = max(v['A'] for k, v in CRJ.items() if k not in ('G1 V(2b UV)', 'G2 IV'))
B_ok = max(v['B'] for v in CRJ.values())
C_ok = max(v['C'] for v in CRJ.values())
C4 = max(v['relative'] for v in C4J.values())
CF = max(max(v) for v in CFJ.values())

BFKL = PREAMBLE + r'''
\title{BFKL evolution of the projectile:\\ all terms $\propto\log\frac{\vee}{k^+}$ from your \texttt{BFKL.tex} and \texttt{Other\_terms.tex}}
\date{}
\begin{document}
\maketitle
\noindent Notation as in your files; $\mathcal N=\frac{1}{(2\pi)^3}\frac{g^4}{16\pi^5}\frac{1}{k^+}$, $K(a,b,z)=\mathcal K(a-z)\cdot\mathcal K(b-z)$, $U^{bc}=U^{cb}$. In the rows: black = your \texttt{BFKL.tex}, \new{blue} = added (the $\log\frac{\vee}{k^+}$ still inside \texttt{Other\_terms.tex}), \old{red} = to be removed. The remaining terms are in \texttt{Remaining\_rows.pdf}.

\section{Answer}\label{sec:answer}
\begin{itemize}
\item \textbf{Your two files are complete.} Row by row, \texttt{BFKL.tex} $+$ \texttt{Other\_terms.tex} $=$ the row minus its $\log\frac{k^+}{\Lambda}$ (JIMWLK) part, numerically (the differences, $\le@@AOK@@$ at $\Lambda=10^{-9}$, are the $O(\Lambda)$ remainders of the $\log\frac{k^+}{\Lambda}$ lines); the terms with an unevaluated $p^+$, $\xi$ or momentum integral are identical to those of your Group files. The only exceptions are the two red lines below.
\item \textbf{Many BFKL terms are still in \texttt{Other\_terms.tex}}: every logarithm there that grows with $\vee$ contains $\log\frac{\vee}{k^+}$,
\begin{equation*}
\log\Big(\frac{\vee}{k^+}-1\Big)=\log\frac{\vee}{k^+}+\log\Big(1-\frac{k^+}{\vee}\Big),\qquad
\log\frac{(\vee-k^+)A+k^+B}{\Lambda A+k^+B}=\log\frac{\vee}{k^+}+\log\frac{(1-\frac{k^+}{\vee})A+\frac{k^+}{\vee}B}{B}+O(\Lambda),\qquad\dots
\end{equation*}
($A=(x-z)^2$, $B=(x-w)^2$; the full table is Sec.~2 of \texttt{LL\_terms.pdf}). They are added in \new{blue}: Group I Rows III and IV, Group II Row I (both parts), Group III Rows I and II. The exceptions are the hidden logarithms (Group III Rows III, V, VI and Group II Row IV), where $\log(\frac{\vee}{k^+}-1)$ cancels against $\left[\frac{1}{p^+}\right]_+$ (resp.\ $\frac{\xi}{(1-\xi)_+}$); there the BFKL term is (B) (resp.\ $H_4$), added in blue for Group III.
\item \textbf{Remove two lines} (\old{red}):
\begin{itemize}
\item Group I Row V, 2b UV: Row V has no $\log\frac{\vee}{k^+}$; this line comes from the missing factor 2 after ``Setting $\Lambda\to0$'' (E1).
\item Group II Row IV, line 1, $-\log\frac{\vee}{k^+}\times$(pole bracket): this is the $\log\frac{\vee}{k^+}$ hidden in $\frac{\xi}{(1-\xi)_+}$; it cancels against the one in $\log(\frac{\vee}{k^+}-1)$ of line 2 of \texttt{Other\_terms.tex}, so it is not a BFKL term.
\end{itemize}
Their poles cancel each other (to $@@PSUM@@$), but their sum is not zero ($@@PFIN@@$ relative): it would add $\log\frac{\vee}{k^+}\log\frac{(w-w')^4}{(x-w)^2(x'-w')^2}$ times the LO structure, which is not in BFKL$\otimes$LO (Sec.~\ref{sec:pole}).
\item \textbf{With these changes:}
\begin{itemize}
\item every row of the BFKL list equals the $\log\frac{\vee}{k^+}$ list of \texttt{LL\_terms.pdf} (to $@@BOK@@$);
\item the @@NF@@ three-$\rho$ terms $\propto\log\frac{\vee}{k^+}$ (rows $+$ your $d^3N_a$, $d^3N_b$) \textbf{cancel}, to $@@WF@@$ (Sec.~\ref{sec:three}). Without the blue terms (your \texttt{BFKL.tex} alone, @@NG@@ terms) they leave $@@WG@@$;
\item your four-$\rho$ two-$\rho$ lists $d^3N_1$, $d^3N_2$ are equal to the ones in the checked sum ($@@C4@@$);
\item the two-$\rho$ terms add up to \textbf{BFKL$\otimes$LO} exactly (9 kernel structures, 109 colour structures, 0 differences; Sec.~\ref{sec:two});
\item the result is your (commented) ``BFKL contribution'', both forms, to $@@CF@@$; in colour-singlet form it is your App.~E formula (Sec.~\ref{sec:final}).
\end{itemize}
\end{itemize}
So
\begin{equation}
\frac{d^3N_{\rm BFKL}}{d^3k}=\mathcal N\log\frac{\vee}{k^+}\;{\rm BFKL}\otimes{\rm LO},
\end{equation}
Eq.~\eqref{eq:BFKLfinal}. The remaining terms (\texttt{Remaining\_rows.pdf}) are your \texttt{Other\_terms.tex} minus the blue terms here, with $\Lambda\to0$, the subtractions $-$(B) and $-\log\frac{\vee}{k^+}H_4$, and C4, C5.

@@ROWS@@

\section{Four \texorpdfstring{$\rho$}{rho}'s}\label{sec:four}
Your \texttt{BFKL.tex}, unchanged: the three-$\rho$ remainders $d^3N_a$, $d^3N_b$ (they enter Sec.~\ref{sec:three}) and the two-$\rho$ lists $d^3N_1$, $d^3N_2$. Numerically, $d^3N_1+d^3N_2$ equals the four-$\rho$ two-$\rho$ terms that entered the checked BFKL sum (Eqs.~\eqref{eq:P12}, \eqref{eq:P13} below) to $@@C4@@$ (charges symmetrized, with c.c.).
@@FOURA@@
@@FOURB@@
@@FOUR1@@
@@FOUR2@@

\section{The three-\texorpdfstring{$\rho$}{rho} terms cancel}\label{sec:three}
Same test as in \texttt{LL\_terms.pdf}, Sec.~7: every three-$\rho$ term of the BFKL list (your lines, the blue lines, (B), and your $d^3N_a$, $d^3N_b$) is renamed to the frame $w,w',z$, charges $x,y,x'$; summed over the 6 placements of the charges and, for Groups I and III, its c.c.; random points, random SU(3) adjoint Wilson lines ($U^{bc}=U^{cb}$), random charges; $\log\frac{\vee}{k^+}=1$. Your Group I Row II line 2 has the $z$ integral done; it is written back as a $z$ integral. Values at one configuration, in units of $i\,\mathcal N\,e^{-ik\cdot(w'-w)}$:
\begin{longtable}{@{}lr@{}}\toprule
term & value\\\midrule
@@TT@@
\midrule
total & $@@TTOT@@$\\\bottomrule
\end{longtable}
\noindent Largest relative residual over three configurations: $@@WF@@$. Your \texttt{BFKL.tex} alone (without the blue terms): $@@WG@@$.

\section{The two-\texorpdfstring{$\rho$}{rho} terms: BFKL}\label{sec:two}
Every row of the BFKL list above equals the $\log\frac{\vee}{k^+}$ list of \texttt{LL\_terms.pdf} ($@@BOK@@$), whose two-$\rho$ terms were added up and compared with BFKL$\otimes$LO kernel by kernel. They are, row by row (genuine two-$\rho$ terms and the commutator terms of the three-$\rho$ terms, $[\rho^a(u),\rho^b(v)]=if^{abc}\rho^c(u)\delta(u-v)$; orange / purple: 168 / 162):
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
$K_D=K(x,x,z)-2K(x,y,z)+K(y,y,z)$, $\delta Y=\log\frac{\vee}{k^+}$, $\mathcal N\delta Y\leftrightarrow\frac{\alpha_s^2\delta Y}{8\pi^6(N_c^2-1)}$: term by term, your singlet form is this with $\mathcal K(z-w)=-\mathcal K(w-z)$ and $x\leftrightarrow y$ in the last two lines. It agrees with BFKL$\otimes$LO to $8\times10^{-15}$ (\texttt{Final\_evolution\_check.pdf}). Hence
\begin{equation}
\frac{d^3N_{\rm BFKL}}{d^3k}=\frac{dN}{d^2w\,d^2w'\,dy}\bigg|_{\rm BFKL}\otimes\int_{w,w'}e^{-ik\cdot(w'-w)},\qquad \delta Y=\log\frac{\vee}{k^+}.
\label{eq:BFKLfinal}
\end{equation}

\section{The two red lines}\label{sec:pole}
Your Group I Row V 2b UV line and Group II Row IV line 1 both carry $\log\frac{\vee}{k^+}$ times a pole bracket. With $E(r)=\frac{2}{\epsilon}+2\gamma_E+\log\frac{r^2\mu^2}{4}$ and $\mathcal S$ the LO structure (renaming as in \texttt{Final\_evolution\_check.pdf}, Sec.~6.1), they are $-\frac{g^4N_c}{16\pi^4}\big[E(x-w)+E(x'-w')\big]\mathcal S$ (with c.c.) and $+\frac{g^4N_c}{8\pi^4}E(w-w')\mathcal S$, times $\frac{1}{(2\pi)^3k^+}\log\frac{\vee}{k^+}$. Pointwise: the poles cancel (relative $@@PSUM@@$), the sum does not ($@@PFIN@@$): it is $\frac{g^4N_c}{16\pi^4}\log\frac{(w-w')^4}{(x-w)^2(x'-w')^2}\mathcal S$. Neither line is in your rows (\texttt{BFKL.tex} $+$ \texttt{Other\_terms.tex} exceeds the row by exactly these lines, check A of Sec.~\ref{sec:checks}), and BFKL$\otimes$LO is reproduced without them.

\section{Checks}\label{sec:checks}
Scripts in \texttt{bfkl\_rows/}:
\begin{center}\small
\begin{tabular}{@{}p{0.62\textwidth}p{0.32\textwidth}@{}}\toprule
check & result\\\midrule
A. \texttt{BFKL.tex} $+$ \texttt{Other\_terms.tex} $=$ row $-$ JIMWLK, row by row (\texttt{check\_rows.py}) & $\le @@AOK@@$ ($=O(\Lambda)$, $\Lambda=10^{-9}$) except G1 V 2b UV, G2 IV (the red lines)\\
\quad terms with an unevaluated integral: identical text & all\\
B. BFKL list (black $+$ blue, without red) $=$ $\log\frac{\vee}{k^+}$ list of \texttt{LL\_terms.pdf} & $\le @@BOK@@$\\
C. remaining terms $=$ \texttt{Finite\_terms\_rows.pdf} & $\le @@COK@@$\\
three-$\rho$ terms cancel: final / your \texttt{BFKL.tex} alone (\texttt{three\_bfkl.py}) & $@@WF@@$ / $@@WG@@$\\
four-$\rho$ two-$\rho$ lists vs the checked ones (\texttt{compare\_4rho.py}) & $@@C4@@$\\
your ``BFKL contribution'' vs the checked sum (\texttt{compare\_final.py}) & $@@CF@@$\\
red lines: poles / sum (\texttt{pole\_bfkl.py}) & $@@PSUM@@$ / $@@PFIN@@$\\
two-$\rho$ sum vs BFKL$\otimes$LO (\texttt{../final\_check/render.py}) & 0 of 109 differ\\\bottomrule
\end{tabular}
\end{center}
\end{document}
'''

rows_t3, tot_t3 = three_table('final')
REP = {'@@AOK@@': tex_num(A_ok), '@@BOK@@': tex_num(B_ok), '@@COK@@': tex_num(C_ok),
       '@@PSUM@@': tex_num(max(p['pole_sum_relative'] for p in pb)), '@@PFIN@@': '%.2f' % min(p['finite_sum_relative'] for p in pb),
       '@@NF@@': str(NF), '@@NG@@': str(NG), '@@WF@@': tex_num(WF), '@@WG@@': '%.2f' % WG,
       '@@C4@@': tex_num(C4), '@@CF@@': tex_num(CF), '@@TTOT@@': '%+.1e' % tot_t3}
out = BFKL
for k, v in REP.items():
    out = out.replace(k, v)
out = out.replace('@@ROWS@@', rows_tex('B')).replace('@@TT@@', '\n'.join(rows_t3)).replace('@@BF@@', BF)
out = out.replace('@@COMMENT@@', comment_section())
for k, n in (('@@FOURA@@', '4rho a'), ('@@FOURB@@', '4rho b'), ('@@FOUR1@@', '4rho 1'), ('@@FOUR2@@', '4rho 2')):
    out = out.replace(k, user_align(n))
open(os.path.join(HERE, 'BFKL_rows.tex'), 'w').write(out)

# ------------------------------------------------------------------------------------------------ remaining terms
A0 = (r"\mathcal{A}_{0} & = \frac{1}{\sqrt{(x-w)^{2}\,(x-z)^{2}-\left[(x-w)\cdot (x-z)\right]^{2}}}\, \Bigg\{"
      r"\mathrm{arc\,tan}\,\left[\frac{(x-z)^{2}-(x-w)\cdot (x-z)}{\sqrt{(x-w)^{2}\,(x-z)^{2}-\left[(x-w)\cdot (x-z)\right]^{2}}}\right]"
      r"\notag \\ & + \mathrm{arc\,tan}\,\left[\frac{(x-w)\cdot (x-z)}{\sqrt{(x-w)^{2}\,(x-z)^{2}-\left[(x-w)\cdot (x-z)\right]^{2}}}\right] \Bigg\}")
REM = PREAMBLE + r'''
\title{The remaining terms: your \texttt{Other\_terms.tex} without the BFKL terms, $\Lambda\to0$}
\date{}
\begin{document}
\maketitle
\noindent Companion of \texttt{BFKL\_rows.pdf}. Your \texttt{Other\_terms.tex}, row by row: \new{blue} = changed or added, \old{red} = removed.

\section{What changed}\label{sec:what}
\begin{itemize}
\item Every logarithm that contains $\log\frac{\vee}{k^+}$ is replaced by its remainder (its $\log\frac{\vee}{k^+}$ part is in \texttt{BFKL\_rows.pdf}), with $\Lambda\to0$; $\vee$ is kept exact:
\begin{equation*}
\log\Big(\frac{\vee}{k^+}-1\Big)\to\log\Big(1-\frac{k^+}{\vee}\Big),\qquad
\log\frac{(\vee-k^+)A+k^+B}{\Lambda A+k^+B}\to\log\frac{(1-\frac{k^+}{\vee})A+\frac{k^+}{\vee}B}{B},\qquad \log\frac{\vee}{\Lambda+k^+}\to0,\ \dots
\end{equation*}
\item Exception, the hidden logarithms: in Group III Rows III, V, VI and Group II Row IV, line 2 keeps $\log(\frac{\vee}{k^+}-1)$, because it cancels the $\log\frac{\vee}{k^+}$ of the subtraction in $\left[\frac{1}{p^+}\right]_+$ (resp.\ $\frac{\xi}{(1-\xi)_+}$). The BFKL term inside line 1 is subtracted instead (\new{blue} lines, C3 and C1).
\item $\Lambda\to0$ everywhere else: logarithms without a large part, $\mathcal A\to\mathcal A_0$, the limits of the plus-prescribed integrals, $\frac{1}{k^++\Lambda}\to\frac{1}{k^+}$, $\frac{k^+}{k^++\Lambda}\to1$, $2\log(1-\frac{\Lambda}{k^+})\to0$.
\item Removed (\old{red}): Group II Row I line 2 (both blocks), which is entirely BFKL; and (C4) the edge-kernel line of Group I Row III $d^3N_{3b}$.
\item C5: the measure of Group III Row III part 2.
\end{itemize}
Row by row this is \texttt{Finite\_terms\_rows.pdf} ($@@COK@@$), and \texttt{BFKL\_rows.pdf} $+$ this file $=$ every row minus its JIMWLK part. The terms are finite for $\epsilon\to0$, $\vee\to\infty$, $\Lambda\to0$, apart from the open finite part of the $\mathcal D$ region (C4). Running coupling and DGLAP are in Group I Row II part 1, Group I Row V 2b UV and Group II Row IV, as before (\texttt{Finite\_terms\_final.pdf}, Secs.~4, 5).
\begin{align}
@@A0@@
\label{eq:A0}
\end{align}

@@ROWS@@
\end{document}
'''
REM = REM.replace('@@COK@@', tex_num(C_ok)).replace('@@A0@@', A0).replace('@@ROWS@@', rows_tex('F'))
open(os.path.join(HERE, 'Remaining_rows.tex'), 'w').write(REM)
print('written BFKL_rows.tex, Remaining_rows.tex')
