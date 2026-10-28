r"""Writes, from your BFKL.tex and Other_terms.tex (input/):

  BFKL_rows.tex       all terms proportional to an explicit log(vee/k+): yours (black), the ones still inside
                      Other_terms.tex (blue), the three lines that are not in the rows (red); three-rho test; what the
                      two-rho terms give; poles.
  Remaining_rows.tex  Other_terms.tex without them, Lambda -> 0 (blue: changed, red: removed).

Only vee > k+, p+ > Lambda are used; Lambda -> 0 where possible, vee exact; the p+ and xi integrals are kept as written.
Numbers from check_rows.json, three_bfkl.json, delta_two_rho.json, pole_bfkl.json, compare_4rho.json, compare_final.json."""
import ast
import json
import os
import re
import rows as R
import extract as X
import check_rows as CR
import delta_two_rho as DT
import assemble as AS
import gen as GN

HERE = R.HERE
UB, UO, EX = CR.UB, CR.UO, CR.EX
J = lambda f: json.load(open(os.path.join(HERE, f)))
CRJ, T3J, PBJ, C4J, CFJ = J('check_rows.json'), J('three_bfkl.json'), J('pole_bfkl.json'), J('compare_4rho.json'), J('compare_final.json')
DTJ = J('delta_two_rho.json')


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
    return [(c5(name, s['B']), 'blue') for s in EX.get(name, []) if s.get('B') is not None]


def remaining_items(name):
    items = []
    for n, s in enumerate(EX[name]):
        if s.get('F') is None:
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
    'G2 IV': r'''\textbf{Both lines: remove (\old{red}).} Neither is a logarithm of your row: with them, \texttt{BFKL.tex} $+$ \texttt{Other\_terms.tex} exceeds the row by exactly these two lines (check A, Sec.~\ref{sec:checks}).
Line 1, $-\log\frac{\vee}{k^+}\times$(line 2 of \texttt{Other\_terms.tex}), is the $\frac{d\xi}{\xi}$ part of the subtraction in $\frac{\xi}{(1-\xi)_+}$ of line 1 of \texttt{Other\_terms.tex}, $-e^{-ik\cdot(y'-y)}\int_{k^+/\vee}^{1}\frac{d\xi}{\xi}$. It is not taken out there: the $\xi$ integral is kept as you wrote it, so it would be counted twice.
Line 2 ($H_4$) is what the region $y'\to y$ of that $\xi$ integral gives for $\vee\gg k^+$; with $\vee$ exact it is not there.
\new{Added}: the $\log\frac{\vee}{k^+}$ part of $\log(\frac{\vee}{k^+}-1)$ in line 2 of \texttt{Other\_terms.tex}. It is the BFKL term of this row, and it carries the pole bracket (Sec.~\ref{sec:pole}).''',
    'G3 I': r'''Your line is correct. \new{Added}: the $\log\frac{\vee}{k^+}$ parts of both lines of \texttt{Other\_terms.tex}.''',
    'G3 II': r'''Your line is correct. \new{Added}: the $\log\frac{\vee}{k^+}$ parts of both lines of \texttt{Other\_terms.tex}.''',
    'G3 III(1)': r'''Not in your \texttt{BFKL.tex}. \new{Added}: the $\log\frac{\vee}{k^+}$ part of $\log(\frac{\vee}{k^+}-1)=\log\frac{\vee}{k^+}+\log(1-\frac{k^+}{\vee})$ in line 2 of \texttt{Other\_terms.tex}. Line 1, the $p^+$ integral up to $\vee-k^+$ with $\left[\frac{1}{p^+}\right]_+$, is kept as you wrote it. It contains no explicit logarithm and gives no BFKL term.''',
    'G3 III(2)': r'''\new{Added}: the $\log\frac{\vee}{k^+}$ part of line 2, as in the first part.''',
    'G3 V': r'''\new{Added}: the $\log\frac{\vee}{k^+}$ part of line 2, as in Row III.''',
    'G3 VI': r'''\new{Added}: the $\log\frac{\vee}{k^+}$ part of line 2, as in Row III.''',
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
    'G2 IV': r'''Line 2: $\log(\frac{\vee}{k^+}-1)\to$ \new{$\log(1-\frac{k^+}{\vee})$}. Its $\log\frac{\vee}{k^+}$ part is in the BFKL list. Line 1 is unchanged: the $\xi$ integral from $\frac{k^+}{\vee}$ with $\frac{\xi}{(1-\xi)_+}$, as you wrote it. In line 3, $\frac{k^+}{k^++\Lambda}\to1$.''',
    'G3 I': r'''\new{Blue}: remainders.''', 'G3 II': r'''\new{Blue}: remainders.''',
    'G3 III(1)': r'''Line 2: $\log(\frac{\vee}{k^+}-1)\to$ \new{$\log(1-\frac{k^+}{\vee})$}. Its $\log\frac{\vee}{k^+}$ part is in the BFKL list. Line 1 is unchanged: the $p^+$ integral up to $\vee-k^+$ with $\left[\frac{1}{p^+}\right]_+$, as you wrote it.''',
    'G3 III(2)': r'''As the first part. \textbf{C5}: the measure is $\int_{w,x',y',z,x}$.''',
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
                    yi, xi = your_items(name), extracted_items(name)
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


def three_table(key):
    d = T3J[key]['seed11']
    rows, tot = [], 0.0
    for labl in T3J[key]['terms']:
        v = d[labl]
        tot += v
        rows.append(r'%s & $%+.4f$\\' % (labl.replace('#', r'line~').replace('_', r'\_'), v))
    return rows, tot


# ------------------------------------------------------------------------------------------------ two-rho difference
PTN = {'p1': 'x', 'p2': "x'"}


def _v(a, b):
    return '%s-%s' % (PTN.get(a, a), PTN.get(b, b))


def _factor(k):
    if k[0] == 'KK':
        (a, b), (c, d) = k[1], k[2]
        return r'\mathcal K(%s)\cdot\mathcal K(%s)' % (_v(a, b), _v(c, d))
    if k[0] == 'inv2':
        return r'\frac{1}{(%s)^2}' % _v(*k[1])
    raise ValueError(k)


def tex_key(key):
    if key[0] == 'pt':
        return ''.join(_factor(k) for k in key[1:])
    if key[0] == 'P':
        return r'\pi P\;' + ''.join(_factor(k) for k in key[1:])
    if key[0] == 'LOG':
        return r'-\pi\log(%s)^2\;' % _v(key[1], key[2]) + ''.join(_factor(k) for k in key[3:])
    raise ValueError(key)


SRC = {'+line 2 G3 III(1)': 'III(1)', '+line 2 G3 III(2)': 'III(2)', '+line 2 G3 V': 'V', '+line 2 G3 VI': 'VI',
       '+line 2 G2 IV': r'II.IV', '-(B) G3 III(1)': r'\old{B}III(1)', '-(B) G3 III(2)': r'\old{B}III(2)',
       '-(B) G3 V': r'\old{B}V', '-(B) G3 VI': r'\old{B}VI', '-H4': r'\old{$H_4$}'}


def delta_table():
    rows = []
    for e in DTJ['Delta'][0]['nonzero']:
        key = ast.literal_eval(e['key'])
        src = ', '.join(SRC[s] for s, v in e['sources'].items() if abs(complex(*v)) > 1e-9)
        rows.append(r'$%s$ & $%+.3f$ & {\footnotesize\leavevmode %s}\\' % (tex_key(key), e['total'][0], src))
    return rows


def withdrawn_B():
    HB = DT.hidden_B_text()
    out = [r'\begin{align}']
    names = ['G3 III(1)', 'G3 III(2)', 'G3 V', 'G3 VI']
    for k, n in enumerate(names):
        t = HB[n].replace("y'", "w'")
        t = t.replace('&', r'&\rc ')          # the colour holds only inside one alignment cell
        out.append(r'&\text{%s:}\quad\rc %s%s' % (n.replace('G3 ', 'Row '), t, r'\notag\\' if k < 3 else r'\label{eq:Bold}'))
    out.append(r'\end{align}')
    return '\n'.join(out)


pb = [v for k, v in PBJ.items() if k != 'pole_items']
WF, WG, WL = T3J['final']['worst'], T3J['as given']['worst'], T3J['line 2']['worst']
WR = T3J['without the line-2 terms']['worst']
NF, NG = len(T3J['final']['terms']), len(T3J['as given']['terms'])
A_ok = max(v['A'] for v in CRJ.values())
A_given = min(CRJ[k]['A_given'] for k in ('G1 V(2b UV)', 'G2 IV'))
B_ok = max(v['B'] for v in CRJ.values())
C_ok = max(v['C'] for v in CRJ.values())
C4 = max(v['relative'] for v in C4J.values())
CF = max(max(v) for v in CFJ.values())
VAL = max(DTJ['validation'].values())
ND = DTJ['Delta'][0]['n_nonzero']
assert all(r['n_nonzero'] == r['n_keys'] == ND for r in DTJ['Delta'])
NL = DTJ['line 2 terms only'][0]['n_nonzero']
NBH = DTJ['(B) and H4 only'][0]['n_nonzero']
PR = max(abs(p['list_pole_over_red_G1V'] - 1) for p in pb)
assert all(list(p['list_poles']) == ['extracted G2 IV #2'] for p in pb)
CANCEL = WF < 1e-12
THREE_SENT = (r'they \textbf{cancel}, to $%s$' % tex_num(WF)) if CANCEL else \
    (r'they \textbf{do not cancel}: the residual is up to $%s$ times the largest term (three configurations). It is exactly the sum of the three blue line-2 terms of Group III Rows III(1), V, VI; all other terms cancel among themselves, to $%s$' % ('%.2f' % WF if WF > 0.01 else tex_num(WF), tex_num(WR)))

BFKL = PREAMBLE + r'''
\title{BFKL evolution of the projectile:\\ all terms $\propto\log\frac{\vee}{k^+}$ from your \texttt{BFKL.tex} and \texttt{Other\_terms.tex},\\ with $\vee$ exact and no large-$p^+$ step}
\date{}
\begin{document}
\maketitle
\noindent Notation as in your files; $\mathcal N=\frac{1}{(2\pi)^3}\frac{g^4}{16\pi^5}\frac{1}{k^+}$, $K(a,b,z)=\mathcal K(a-z)\cdot\mathcal K(b-z)$, $U^{bc}=U^{cb}$. In the rows: black = your \texttt{BFKL.tex}, \new{blue} = added (the $\log\frac{\vee}{k^+}$ still inside \texttt{Other\_terms.tex}), \old{red} = to be removed. The remaining terms are in \texttt{Remaining\_rows.pdf} and \texttt{Rest\_terms.pdf}.

\section{Answer}\label{sec:answer}
\paragraph{Rule.} The only constraints are $\vee>k^+$ and $p^+>\Lambda$. $\Lambda\to0$ wherever possible; $\vee$ is exact. Only explicit logarithms are split, and exactly:
\begin{equation*}
\log\Big(\frac{\vee}{k^+}-1\Big)=\log\frac{\vee}{k^+}+\log\Big(1-\frac{k^+}{\vee}\Big),\qquad
\log\frac{(\vee-k^+)A+k^+B}{\Lambda A+k^+B}=\log\frac{\vee}{k^+}+\log\frac{(1-\frac{k^+}{\vee})A+\frac{k^+}{\vee}B}{B}+O(\Lambda),\qquad\dots
\end{equation*}
($A=(x-z)^2$, $B=(x-w)^2$; the full table is Sec.~2 of \texttt{LL\_terms.pdf}). Every $p^+$ and $\xi$ integral is kept as you wrote it: $\left[\frac{1}{p^+}\right]_+$, $\frac{\xi}{(1-\xi)_+}$, upper limit $\vee-k^+$, lower limit $\frac{k^+}{\vee}$. There is no large-$p^+$ (small-$\xi$) step, so these integrals give no BFKL term.

\paragraph{Withdrawn from the previous version.}
\begin{itemize}
\item The terms (B) of Group III Rows III (both parts), V, VI, Eq.~\eqref{eq:Bold}, and $H_4$ of Group II Row IV. Both came from the region of large $p^+$ (resp.\ $\xi\to\frac{k^+}{\vee}$).
\item The subtractions $-$(B) and $-\log\frac{\vee}{k^+}H_4$ in the remaining terms.
\item The rule that $\log(\frac{\vee}{k^+}-1)$ in line 2 of these rows is not a BFKL term. It relied on the same region.
\item \texttt{Why\_log\_vee.pdf}.
\end{itemize}

\paragraph{Result.}
\begin{itemize}
\item \textbf{Your two files are complete.} Row by row, \texttt{BFKL.tex} (without the red lines) $+$ \texttt{Other\_terms.tex} $=$ the row minus its $\log\frac{k^+}{\Lambda}$ (JIMWLK) part, numerically. The differences, $\le@@AOK@@$ at $\Lambda=10^{-9}$, are the $O(\Lambda)$ remainders. The terms with an unevaluated $p^+$, $\xi$ or momentum integral are identical to those of your Group files.
\item \textbf{Remove three lines} (\old{red}):
\begin{itemize}
\item Group I Row V 2b UV: Row V has no $\vee$; this line comes from the missing factor 2 after ``Setting $\Lambda\to0$'' (E1).
\item Group II Row IV, both lines: $-\log\frac{\vee}{k^+}\times$(line 2) and $\log\frac{\vee}{k^+}H_4$. They are statements about the $\xi$ integral of line 1, not logarithms of the row.
\end{itemize}
With these lines, \texttt{BFKL.tex} $+$ \texttt{Other\_terms.tex} exceeds the row (by $\ge@@AGIVEN@@$ relative).
\item \textbf{Add} (\new{blue}) every $\log\frac{\vee}{k^+}$ that is still in \texttt{Other\_terms.tex}:
\begin{itemize}
\item Group I Rows III and IV;
\item Group II Row I (both parts) and Row IV (line 2);
\item Group III Rows I, II, III (both parts), V, VI (line 2).
\end{itemize}
Then the list is exactly the $\log\frac{\vee}{k^+}$ part of your rows (to $@@BOK@@$).
\item \textbf{What this list gives.}
\begin{itemize}
\item \textbf{Three $\rho$'s} (@@NF@@ terms, Sec.~\ref{sec:three}): @@THREE@@. Your \texttt{BFKL.tex} alone leaves $@@WG@@$, and the three line-2 terms of Group III Rows III, V, VI alone leave $@@WL@@$.
\item \textbf{Two $\rho$'s: not BFKL$\otimes$LO.} The previous list had (B) and $H_4$ instead of the line-2 terms, and its two-$\rho$ terms added up to BFKL$\otimes$LO. The change is $+$(the line-2 terms of Group III Rows III, V, VI and Group II Row IV) $-$(B) $-H_4$. Its two-$\rho$ terms do not cancel: they are non-zero in all @@ND@@ kernel structures in which they appear (@@NL@@ from the line-2 terms, @@NBH@@ from (B) and $H_4$; Sec.~\ref{sec:two}).
\item \textbf{Pole.} The blue term of Group II Row IV is the only term of the list with $\frac{1}{\epsilon}$. Its pole equals that of the red Group I Row V line (to $@@PR@@$), so nothing in the list cancels it (Sec.~\ref{sec:pole}).
\item Your commented ``BFKL contribution'' (your App.~E formula, BFKL$\otimes$LO) is the sum of the previous list, not of this one.
\end{itemize}
\end{itemize}
So, with $\vee$ exact and the $p^+$, $\xi$ integrals kept as written, the coefficient of the explicit $\log\frac{\vee}{k^+}$ is \emph{not} BFKL$\otimes$LO.

The rows themselves are unchanged (check A). Only the bookkeeping changes. The rest of their $\vee$ dependence is in the remaining terms:
\begin{itemize}
\item the $p^+$ integrals up to $\vee-k^+$ (line 1 of Group III Rows III, V, VI);
\item the $\xi$ integrals from $\frac{k^+}{\vee}$ (lines 1 and 3 of Group II Row IV);
\item $\log(1-\frac{k^+}{\vee})$ and $\log(1-\frac{k^+}{\vee})\log\frac{\vee}{k^+}$;
\item $R^i$ of Group I Row II.
\end{itemize}

\paragraph{A second exact bookkeeping (for your decision).} The split above is not the only exact one. In the subtraction of $\left[\frac{1}{p^+}\right]_+$ one may use $\frac{1}{p^+}=\frac{1}{p^++k^+}+\frac{k^+}{p^+(p^++k^+)}$. With $\int_0^{\vee-k^+}\frac{dp^+}{p^++k^+}=\log\frac{\vee}{k^+}$ this gives, exactly,
\begin{equation*}
\int_0^{\vee-k^+}dp^+\left[\frac{1}{p^+}\right]_+f(p^+)=\int_0^{\vee-k^+}dp^+\,\frac{k^+\,[f(p^+)-f(0)]}{p^+(p^++k^+)}+\int_0^{\vee-k^+}dp^+\,\frac{f(p^+)}{p^++k^+}-f(0)\log\frac{\vee}{k^+}.
\end{equation*}
For $\frac{\xi}{(1-\xi)_+}$ the analogue is $\frac{1}{\xi(1-\xi)}=\frac{1}{\xi}+\frac{1}{1-\xi}$, which gives $-g(1)\log\frac{\vee}{k^+}$.
\begin{itemize}
\item Your line~2 is $f(0)\big[\log\frac{k^+}{\Lambda}+\log(\frac{\vee}{k^+}-1)\big]$ (resp.\ with $g(1)$), so the explicit $-f(0)\log\frac{\vee}{k^+}$ cancels its $\log\frac{\vee}{k^+}$ exactly. Group III Rows III, V, VI and Group II Row IV then have no BFKL term at all.
\item Your red line~1 of Group II Row IV is exactly this $-g(1)\log\frac{\vee}{k^+}$. It belongs to this bookkeeping, with line~1 of \texttt{Other\_terms.tex} rewritten accordingly.
\end{itemize}
Both bookkeepings are exact; they differ only in which $\vee$ dependence is called an explicit logarithm:
\begin{center}\small
\begin{tabular}{@{}p{0.34\textwidth}p{0.27\textwidth}p{0.3\textwidth}@{}}\toprule
 & this file (explicit logarithms only) & second bookkeeping (subtraction split)\\\midrule
BFKL term of G3 III, V, VI, G2 IV & $\log\frac{\vee}{k^+}$ part of line 2 & none\\
three-$\rho$ terms & do not cancel ($@@WF@@$) & cancel ($@@WR@@$)\\
two-$\rho$ terms vs BFKL$\otimes$LO & differ in @@ND@@ structures & differ in @@NBH@@ structures ($-$(B), $-H_4$)\\
pole $\times\log\frac{\vee}{k^+}$ in the list & yes (G2 IV) & none\\\bottomrule
\end{tabular}
\end{center}
In neither is the two-$\rho$ sum BFKL$\otimes$LO. That needs (B) and $H_4$, i.e.\ the large-$p^+$ region of the integrals that both keep unevaluated.

@@ROWS@@

\section{Four \texorpdfstring{$\rho$}{rho}'s}\label{sec:four}
Your \texttt{BFKL.tex}, unchanged: the three-$\rho$ remainders $d^3N_a$, $d^3N_b$ (they enter Sec.~\ref{sec:three}) and the two-$\rho$ lists $d^3N_1$, $d^3N_2$. Numerically, $d^3N_1+d^3N_2$ equals the four-$\rho$ two-$\rho$ terms of the checked sum (\texttt{Final\_evolution\_check.pdf}, Eqs.~P12, P13) to $@@C4@@$ (charges symmetrized, with c.c.).
@@FOURA@@
@@FOURB@@
@@FOUR1@@
@@FOUR2@@

\section{The three-\texorpdfstring{$\rho$}{rho} terms}\label{sec:three}
Same test as in \texttt{LL\_terms.pdf}, Sec.~7. Every three-$\rho$ term of the list (your lines, the blue lines, your $d^3N_a$, $d^3N_b$) is
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
\noindent Largest relative residual over three configurations: $@@WF@@$.
\begin{itemize}
\item Your \texttt{BFKL.tex} alone (without the blue terms): $@@WG@@$.
\item The three blue line-2 terms of Group III Rows III(1), V, VI alone: $@@WL@@$.
\end{itemize}

\section{The two-\texorpdfstring{$\rho$}{rho} terms}\label{sec:two}
Every row of the list above equals the $\log\frac{\vee}{k^+}$ part of your rows (check B). The previous list was the same except in four places:
\begin{itemize}
\item it had (B), Eq.~\eqref{eq:Bold}, and $H_4$ (your red line 2 of Group II Row IV);
\item it did not have the blue line-2 terms of Group III Rows III(1), III(2), V, VI and Group II Row IV.
\end{itemize}
Its two-$\rho$ terms (genuine two-$\rho$ terms and the commutator terms of the three-$\rho$ terms, $[\rho^a(u),\rho^b(v)]=if^{abc}\rho^c(u)\delta(u-v)$, with all c.c.'s) were added up kernel by kernel and are BFKL$\otimes$LO (\texttt{Final\_evolution\_check.pdf}, \texttt{BFKL\_final.pdf}). The withdrawn terms (B) were:
@@BOLD@@
Hence
\begin{equation}
\frac{d^3N}{d^3k}\bigg|^{\rm two\ \rho}_{\log\frac{\vee}{k^+}}=\mathcal N\log\frac{\vee}{k^+}\Big[{\rm BFKL}\otimes{\rm LO}+\Delta_2\Big],
\qquad \Delta_2=\Big[\new{\text{line 2 of G3 III(1), III(2), V, VI, G2 IV}}-\old{({\rm B})}-\old{H_4}\Big]_{\rm two\ \rho}.
\end{equation}
$\Delta_2$ is computed with the reordering engine of \texttt{../two\_rho} (the transcriptions of \texttt{../bfkl\_final/rows.py}). Each transcription equals its LaTeX term here pointwise to $@@VAL@@$ (\texttt{delta\_two\_rho.py}). Each kernel monomial is evaluated at random points, the colour numerically; colour factors without $U(z)$ are integrated over $z$ in dimensional regularization, $\int_zK(a,b,z)=\pi[P-\log(a-b)^2]$ (MS-bar, with $P$ containing $\frac{2}{\epsilon}$). $\Delta_2$ is non-zero in every one of its @@ND@@ structures, for three configurations. Values at one configuration (units of $\mathcal N\log\frac{\vee}{k^+}$, $x$ and $x'$ the two charges):
\begin{longtable}{@{}p{0.5\textwidth}rp{0.38\textwidth}@{}}\toprule
structure & $\Delta_2$ & from (blue line 2 of III(1), \dots, II.IV; \old{B} = $-$(B); \old{$H_4$} = $-H_4$)\\\midrule
@@DT@@
\bottomrule
\end{longtable}
So BFKL$\otimes$LO is reproduced only with (B) and $H_4$ in place of the line-2 terms, i.e.\ only with the large-$p^+$ step.

\section{Poles}\label{sec:pole}
With $E(r)=\frac{2}{\epsilon}+2\gamma_E+\log\frac{r^2\mu^2}{4}$ and $\mathcal S$ the LO structure (renaming as in \texttt{Final\_evolution\_check.pdf}, Sec.~6.1):
\begin{itemize}
\item \textbf{Your two red lines} with a pole are Group I Row V 2b UV and Group II Row IV line 1:
\begin{equation*}
-\frac{g^4N_c}{16\pi^4}\big[E(x-w)+E(x'-w')\big]\mathcal S\ \text{(with c.c.)}\qquad\text{and}\qquad+\frac{g^4N_c}{8\pi^4}E(w-w')\mathcal S,
\end{equation*}
times $\frac{1}{(2\pi)^3k^+}\log\frac{\vee}{k^+}$. Pointwise, their poles cancel (relative $@@PSUM@@$), but their sum does not ($@@PFIN@@$): it is $\frac{g^4N_c}{16\pi^4}\log\frac{(w-w')^4}{(x-w)^2(x'-w')^2}\mathcal S$. Neither line is in your rows (check A).
\item \textbf{The list} has one term with a pole: the blue term of Group II Row IV, $-\frac{g^4N_c}{8\pi^4}E(w-w')\mathcal S\times\frac{1}{(2\pi)^3k^+}\log\frac{\vee}{k^+}$. It is minus your red line 1, so its pole is that of the red Group I Row V line (to $@@PR@@$, three configurations), and nothing in the list cancels it.
\item In the rows, the other pole terms that depend on $\vee$ are all in the remaining terms:
\begin{itemize}
\item Group I Row II part 1: the pole bracket $\times\log(1-\frac{k^+}{\vee})$;
\item Group II Row IV: the pole bracket $\times$ the $\xi$ integral of line 1 (from $\frac{k^+}{\vee}$), and $\times\log(1-\frac{k^+}{\vee})$ of line 2.
\end{itemize}
\end{itemize}

\section{Checks}\label{sec:checks}
Scripts in \texttt{bfkl\_rows/}:
\begin{center}\small
\begin{tabular}{@{}p{0.62\textwidth}p{0.32\textwidth}@{}}\toprule
check & result\\\midrule
A. \texttt{BFKL.tex} (without red) $+$ \texttt{Other\_terms.tex} $=$ row $-$ JIMWLK, row by row (\texttt{check\_rows.py}) & $\le @@AOK@@$ ($=O(\Lambda)$, $\Lambda=10^{-9}$)\\
\quad the same with the red lines & G1 V 2b UV, G2 IV: $\ge@@AGIVEN@@$\\
\quad terms with an unevaluated integral: identical text & all\\
B. BFKL list (black $+$ blue, without red) $=$ $\log\frac{\vee}{k^+}$ part of the rows (exact split, \texttt{../llsplit/gen.py}) & $\le @@BOK@@$\\
C. remaining terms $=$ the rest of the rows, $\Lambda\to0$ & $\le @@COK@@$\\
three-$\rho$ terms: list / your \texttt{BFKL.tex} alone / line-2 terms alone / second bookkeeping (\texttt{three\_bfkl.py}) & $@@WF@@$ / $@@WG@@$ / $@@WL@@$ / $@@WR@@$\\
transcriptions of the line-2 terms, (B), $H_4$ vs the LaTeX (\texttt{delta\_two\_rho.py}) & $\le@@VAL@@$\\
two-$\rho$ change $\Delta_2$: non-zero structures (\texttt{delta\_two\_rho.py}) & @@ND@@ of @@ND@@\\
four-$\rho$ two-$\rho$ lists vs the checked ones (\texttt{compare\_4rho.py}) & $@@C4@@$\\
your ``BFKL contribution'' vs the checked sum of the previous list (\texttt{compare\_final.py}) & $@@CF@@$\\
red lines: poles / sum; pole of the list vs red G1 V (\texttt{pole\_bfkl.py}) & $@@PSUM@@$ / $@@PFIN@@$; $@@PR@@$\\\bottomrule
\end{tabular}
\end{center}
\end{document}
'''

rows_t3, tot_t3 = three_table('final')
REP = {'@@AOK@@': tex_num(A_ok), '@@AGIVEN@@': '%.2f' % A_given, '@@BOK@@': tex_num(B_ok), '@@COK@@': tex_num(C_ok),
       '@@PSUM@@': tex_num(max(p['pole_sum_relative'] for p in pb)), '@@PFIN@@': '%.2f' % min(p['finite_sum_relative'] for p in pb),
       '@@PR@@': tex_num(PR), '@@VAL@@': tex_num(VAL),
       '@@NF@@': str(NF), '@@NG@@': str(NG), '@@WF@@': ('%.2f' % WF) if WF > 0.01 else tex_num(WF), '@@WG@@': '%.2f' % WG,
       '@@WL@@': ('%.2f' % WL) if WL > 0.01 else tex_num(WL), '@@THREE@@': THREE_SENT, '@@WR@@': tex_num(WR),
       '@@ND@@': str(ND), '@@NL@@': str(NL), '@@NBH@@': str(NBH),
       '@@C4@@': tex_num(C4), '@@CF@@': tex_num(CF), '@@TTOT@@': '%+.1e' % tot_t3}
out = BFKL
for k, v in REP.items():
    out = out.replace(k, v)
out = out.replace('@@ROWS@@', rows_tex('B')).replace('@@TT@@', '\n'.join(rows_t3))
out = out.replace('@@DT@@', '\n'.join(delta_table())).replace('@@BOLD@@', withdrawn_B())
for k, n in (('@@FOURA@@', '4rho a'), ('@@FOURB@@', '4rho b'), ('@@FOUR1@@', '4rho 1'), ('@@FOUR2@@', '4rho 2')):
    out = out.replace(k, user_align(n))
assert '@@' not in out

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
The only constraints are $\vee>k^+$ and $p^+>\Lambda$; $\Lambda\to0$ wherever possible, $\vee$ exact, no large-$p^+$ step.
\begin{itemize}
\item Every explicit logarithm that contains $\log\frac{\vee}{k^+}$ is replaced by its remainder, with $\Lambda\to0$ (its $\log\frac{\vee}{k^+}$ part is in \texttt{BFKL\_rows.pdf}):
\begin{equation*}
\log\Big(\frac{\vee}{k^+}-1\Big)\to\log\Big(1-\frac{k^+}{\vee}\Big),\qquad
\log\frac{(\vee-k^+)A+k^+B}{\Lambda A+k^+B}\to\log\frac{(1-\frac{k^+}{\vee})A+\frac{k^+}{\vee}B}{B},\qquad \log\frac{\vee}{\Lambda+k^+}\to0,\ \dots
\end{equation*}
This includes line 2 of Group III Rows III, V, VI and Group II Row IV.
\item Every $p^+$ and $\xi$ integral is kept as you wrote it. This covers $\left[\frac{1}{p^+}\right]_+$ up to $\vee-k^+$ in line 1 of Group III Rows III, V, VI, and $\frac{\xi}{(1-\xi)_+}$ from $\frac{k^+}{\vee}$ in line 1 of Group II Row IV. Nothing is subtracted from them. The subtractions $-$(B) and $-\log\frac{\vee}{k^+}H_4$ of the previous version are withdrawn.
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
Row by row, \texttt{BFKL\_rows.pdf} $+$ this file $=$ every row minus its JIMWLK part ($@@COK@@$). The $\vee$ dependence that is not an explicit $\log\frac{\vee}{k^+}$ is here:
\begin{itemize}
\item the $p^+$ and $\xi$ integrals above;
\item $\log(1-\frac{k^+}{\vee})$;
\item $\log(1-\frac{k^+}{\vee})\log\frac{\vee}{k^+}$ (Group I Row II part 1; $\le\log^22$ for all $\vee>k^+$);
\item $R^i$ (Group I Row II part 2).
\end{itemize}
Open: the finite part of the $\mathcal D$ region (C4).
\begin{align}
@@A0@@
\label{eq:A0}
\end{align}

@@ROWS@@
\end{document}
'''
REM = REM.replace('@@COK@@', tex_num(C_ok)).replace('@@A0@@', A0).replace('@@ROWS@@', rows_tex('F'))
if __name__ == '__main__':
    open(os.path.join(HERE, 'BFKL_rows.tex'), 'w').write(out)
    open(os.path.join(HERE, 'Remaining_rows.tex'), 'w').write(REM)
    print('written BFKL_rows.tex, Remaining_rows.tex')
