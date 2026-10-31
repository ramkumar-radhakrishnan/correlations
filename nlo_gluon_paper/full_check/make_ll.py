r"""Writes Checks_JIMWLK_BFKL.tex from your new Group_I/II/III.tex and 4rho.tex (input/):
the checks of every step of your files, every term proportional to log(k+/Lambda) (JIMWLK) and to log(vee/k+) (BFKL)
row by row, the exact identities that bring the BFKL terms out of line 1 of Group III Rows III, V, VI and Group II
Row IV, the four- and three-rho cancellations, the two-rho sums and the comparison with your two formulas.

The split of each term is made by ../llsplit/gen.py on your new files (ll.py), checked in ll_checks.py; the two-rho
sections and the targets are those of ../final_check/Final_evolution_check.tex, with log(k+/Lambda), log(vee/k+).
The remaining terms, running coupling, DGLAP and the finite terms: make_rest.py."""
import json
import sys
import os
import re
import ll  # noqa: F401  (your new files; must come first)
import assemble as AS
import gen as GN
import logsplit as LS
import three_check as T3

HERE = os.path.dirname(os.path.abspath(__file__))
import bterms as BT  # noqa: E402
BLK, SPL = GN.all_splits()
LKL, LVK, LVK1 = LS.LKL, LS.LVK, LS.LVK1
FEC = open(os.path.join(HERE, '..', 'final_check', 'Final_evolution_check.tex')).read()
J = lambda f: json.load(open(os.path.join(HERE, f)))
CS, PC, FF, SL, RV = J('ll_split.json'), J('ll_poles.json'), J('four_file.json'), J('same_lists.json'), J('rows_vs_files.json')
TC = json.load(open(os.environ.get('TC_JSON', os.path.join(HERE, 'll_three.json'))))   # TC_JSON: only for a trial run
STEPS = {g: J('step_%s.json' % g) for g in ('Group_I', 'Group_II', 'Group_III')}
ECJ = json.load(open(os.path.join(HERE, '..', 'bfkl_rows', 'exact_check.json')))
TGJ = json.load(open(os.path.join(HERE, '..', 'final_check', 'targets_check.json')))


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
\newcommand{\NN}{\mathcal N\log\frac{k^+}{\Lambda}}
\newcommand{\NNv}{\mathcal N\log\frac{\vee}{k^+}}
\newcommand{\Lv}{\log{\left(\frac{\vee}{k^{+}}\right)}}
\def\mf#1{\mathbb{#1}}
'''

# ------------------------------------------------------------------------------------------------ rows
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

HIDDEN = {'G3 III(1)': 'III(1)', 'G3 III(2)': 'III(2)', 'G3 V': 'V', 'G3 VI': 'VI'}


HB = BT.hidden_B_text()     # the terms (B), Eq. (I3): keyed by row, with w'
LINE1_ROWS = list(HIDDEN) + ['G2 IV']


def negate(t):
    t = t.strip()
    if t.startswith('&'):
        t = t[1:].strip()
    if t.startswith('+'):
        return '-' + t[1:]
    if t.startswith('-'):
        return '+' + t[1:]
    return '-' + t


def explicit_B(name):
    """(i): the log(vee/k+) part of log(vee/k+ - 1) in line 2 (I1)"""
    return GN.split_term(BLK[name].terms[1], hidden=False)['B']


def remainder_line2(name):
    """line 2 of a row with BFKL terms in line 1: [log(k+/Lambda) + log(vee/k+ - 1)] -> log(1 - k+/vee)"""
    return GN.split_term(BLK[name].terms[1], hidden=False)['F']
H4 = (r"\frac{1}{(2\pi)^{3}}\,\frac{g^{4}\,N_{c}}{8\pi^{5}}\,\frac{1}{k^{+}}\,%s\,\int_{x,x',z,w,w'}\,e^{-ik\cdot (w'-w)}\,"
      r"\frac{(w-z)\cdot(w'-z)}{(w-z)^{2}\,(w'-z)^{2}}\,\frac{(x-z)\cdot(x'-z)}{(x-z)^{2}\,(x'-z)^{2}}\notag \\ & "
      r"\left[U^{e'c}(z) - U^{e'c}(x')\right] \times\, \left[U^{ec}(z) - U^{ce}(x)\right]\, \rho^{e'}(x') \rho^{e}(x)")


def c5(name, t):
    if name != 'G3 III(2)':
        return t
    for m in (r'\int_{w,w\',y\',z,y,x}', r'\int_{w,y\',z,y,x}'):
        t = t.replace(m.replace("\\'", "'"), r"\old{%s}\,\new{\int_{w,x',y',z,x}}" % m.replace("\\'", "'"))
    return t


def neg(t):
    t = t.lstrip()
    assert t.startswith('-'), t[:30]
    return '+' + t[1:]


def j_items(name):
    items = []
    for n, s in enumerate(SPL[name]):
        if s.get('J') is None:
            continue
        t = c5(name, s['J'])
        if name == 'G1 V(2b UV)':
            t = t.replace(r'2\,' + LKL, r'\new{2}\,' + LKL)
        if (name, n) == ('G1 III N3a', 1):
            items += [(t, 'red'), (T3.D_TERM % LKL, 'blue')]
            continue
        items.append((t, None))
    return items


def b_items(name):
    items = [(c5(name, s['B']), None) for s in SPL[name] if s.get('B') is not None]
    if name in LINE1_ROWS:
        b = c5(name, explicit_B(name))
        items += [(b, None), (negate(b), 'blue')]                         # (i), (ii)
        items.append((HB[name], 'blue') if name in HIDDEN else (H4 % LVK, 'blue'))   # (iii)
    return items


A0 = (r"\mathcal{A}_{0} & = \frac{1}{\sqrt{(x-w)^{2}\,(x-z)^{2}-\left[(x-w)\cdot (x-z)\right]^{2}}}\, \Bigg\{"
      r"\mathrm{arc\,tan}\,\left[\frac{(x-z)^{2}-(x-w)\cdot (x-z)}{\sqrt{(x-w)^{2}\,(x-z)^{2}-\left[(x-w)\cdot (x-z)\right]^{2}}}\right]"
      r"\notag \\ & + \mathrm{arc\,tan}\,\left[\frac{(x-w)\cdot (x-z)}{\sqrt{(x-w)^{2}\,(x-z)^{2}-\left[(x-w)\cdot (x-z)\right]^{2}}}\right] \Bigg\}")
LAM0 = [(r'\int_{\Lambda}^{k^{+}-\Lambda}', r'\int_{0}^{k^{+}}'), (r'\int_{k^{+}-\Lambda}^{\Lambda}', r'\int_{k^{+}}^{0}'),
        (r'\frac{1}{k^{+} + \Lambda}', r'\frac{1}{k^{+}}'),
        (r'\int_{\frac{k^{+}}{\vee}}^{\frac{k^{+}}{k^{+}+\Lambda}}', r'\int_{\frac{k^{+}}{\vee}}^{1}'),
        (r'\mathcal{A}', r'\mathcal{A}_{0}')]


def lam0(t):
    for a, b in LAM0:
        t = t.replace(a, b)
    assert r'\Lambda' not in t, t[:200]
    return t


def lab(kind, name):
    return 'eq:%s:%s' % (kind, re.sub(r'[^A-Za-z0-9]', '', name))


def align(name, items, kind):
    if not items:
        return ''
    s = AS.assemble(BLK[name], items)
    s = s.replace('\\begin{align}\n', '\\begin{align}\n', 1)
    s = re.sub(r'\n\s*\n', '\n', s)
    return s.replace('\n\\end{align}', r' \label{%s}' % lab(kind, name) + '\n\\end{align}')


# notes in the LL document
NOTE_LL = {
    'G1 II(1)': r'''Your bracket $\big[\frac{11}{3}-2\log\frac{k^+}{\Lambda}+\log(1-\frac{k^+}{\vee})\big]$ multiplies the pole bracket; only its $-2\log\frac{k^+}{\Lambda}$ is a large logarithm. The term $\log(1-\frac{k^+}{\vee})\log\frac{\vee}{k^+}$ is not a BFKL term: its coefficient depends on $\vee$, it is bounded for every $\vee>k^+$ ($0\le-\log(1-\frac{k^+}{\vee})\log\frac{\vee}{k^+}\le\log^22$, maximum at $\vee=2k^+$) and vanishes for $\vee\to\infty$; it stays with the remaining terms.''',
    'G1 II(2a)': r'''Line 1 has $\big[\log\frac{\vee}{k^+}+\log\frac{k^+}{\Lambda}\big]$ and enters both lists. In lines 2 and 3 the $z$ integral is already done (hence $\frac{1}{32\pi^4}$); the measure should read $\int_{x,y,x'}$.''',
    'G1 II(2b)': r'''No large logarithm ($R^i$ depends on $\vee$ only through $\frac{k^+}{\vee}$).''',
    'G1 III N1': r'''Lines 1 and 4 (your green $\frac{1}{\Lambda}$ brackets) are not logarithms. Together with line 1 of $d^3N_2$ and line 1 of $d^3N_{3a}$ they cancel exactly (\texttt{Remaining\_RC\_DGLAP\_finite.pdf}).''',
    'G1 III N2': r'''Line 1 is a $\frac{1}{\Lambda}$ term (see $d^3N_1$). Your two red lines cancel each other and are left out.''',
    'G1 III N3a': r'''\old{Red}: the first logarithmic line, with the edge kernel $\frac{(x'-w')\cdot(y-z)}{(x'-w')^2(y-z)^2(z-w)^2}$, is not integrable at $z\to w$. \new{Blue}: it is replaced by $\mathcal D$ (fix F3, \texttt{BFKL\_final.pdf} Sec.~2.1): $U^{bc}(z)\to U^{bc}(w)$ and the kernel $\mathcal K(x'-w')\cdot\mathcal K(y-w)\,\mathcal K(x-z)\cdot\mathcal K(w-z)$. Without this change the three-$\rho$ terms do not cancel (Sec.~\ref{sec:three}).''',
    'G1 III N3b': r'''No large logarithm. ($\mathcal A$ and the logarithms of $d^3N_{3a}$, $d^3N_{3b}$ have finite limits for $\Lambda\to0$; \texttt{Remaining\_RC\_DGLAP\_finite.pdf}.)''',
    'G1 V(2b UV)': r'''Your bracket $\big[2\log\frac{k^+}{\Lambda}+2\log(1-\frac{\Lambda}{k^+})-\frac{11}{6}\big]$ is correct (E1 is now in your file; Row V has $0<p^+<k^+$, so no $\vee$). The poles of the $\log\frac{k^+}{\Lambda}$ terms cancel: Group I Row II : Row V 2b UV : Group II Row IV $=@@P1@@:@@P2@@:1$, sum $@@PSUM@@$ (relative; \texttt{ll\_checks.py poles}).''',
    'G1 V(2b fin)': r'''Lines 4--6 carry $\log\frac{k^+-\Lambda}{\Lambda}=\log\frac{k^+}{\Lambda}+\log(1-\frac{\Lambda}{k^+})$.''',
    'G2 I(2a)': r'''No large logarithm: $\vee$ cancels between numerator and denominator of the first logarithm, and the second has the finite limit $\log\frac{(x-w)^2}{(x'-w')^2}$.''',
    'G2 IV': r'''\textbf{BFKL terms in line 1.} Line 1 is a $\xi$ integral down to $\xi=\frac{k^+}{\vee}$; it contains $\log\frac{\vee}{k^+}$ terms that are brought out by exact identities (Sec.~\ref{sec:exact}):
\begin{itemize}
\item[(i)] line 2: $\log(\frac{\vee}{k^+}-1)=\log\frac{\vee}{k^+}+\log(1-\frac{k^+}{\vee})$; its $\log\frac{\vee}{k^+}$ part is the first line below;
\item[(ii)] line 1, the subtraction in $\frac{\xi}{(1-\xi)_+}$: $-g(1)\log\frac{\vee}{k^+}$, Eq.~\eqref{eq:I2xi}, exactly minus (i), poles included (second line, \new{blue});
\item[(iii)] line 1, the $\log(y'-y)^2$ term of $\frac{1-\xi}{\xi}$ at $y'=y$: $\log\frac{\vee}{k^+}H_4$, Eq.~\eqref{eq:I4} (third line, \new{blue}).
\end{itemize}
(i) and (ii) cancel: the BFKL term of this row is $\log\frac{\vee}{k^+}H_4$. Its $\log\frac{k^+}{\Lambda}$ term is as you have it.''',
    'G3 III(1)': r'''\textbf{BFKL terms in line 1.} Line 1 is a $p^+$ integral up to $\vee-k^+$; exact identities bring out (Sec.~\ref{sec:exact}):
\begin{itemize}
\item[(i)] line 2: the $\log\frac{\vee}{k^+}$ part of $\log(\frac{\vee}{k^+}-1)$ (first line below);
\item[(ii)] line 1, the subtraction of $\left[\frac{1}{p^+}\right]_+$: $-f(0)\log\frac{\vee}{k^+}$, Eq.~\eqref{eq:I2}, exactly minus (i) (\new{blue});
\item[(iii)] line 1, the $-\frac{1}{k^+}\delta_{im}\delta_{jk'}$ term with $G(y')\to G(z)$: the term (B), Eq.~\eqref{eq:I3} (\new{blue}).
\end{itemize}
(i) and (ii) cancel, so the BFKL term of this part is (B). Its $\log\frac{k^+}{\Lambda}$ term is as you have it.''',
    'G3 III(2)': r'''(i), (ii), (iii) as in the first part; (B) comes from the $-\frac{p^+}{k^{+2}}$ term, Eqs.~\eqref{eq:I3}, \eqref{eq:I5}. \textbf{C5}: the measure is $\int_{w,x',y',z,x}$ (the row has $\rho(x')\rho(x)$; there is no $y$ and no $w'$).''',
    'G3 V': r'''(i), (ii), (iii) as in Row III: (i) and (ii) cancel, and (B) comes from the $-\frac{1}{k^+}\delta_{im}\delta_{jk'}$ term.''',
    'G3 VI': r'''(i), (ii), (iii) as in Row III: (i) and (ii) cancel, and (B) comes from the $-\frac{1}{k^+}\delta_{im}\delta_{jk'}$ term.''',
}
pc = list(PC.values())
NOTE_LL['G1 V(2b UV)'] = (NOTE_LL['G1 V(2b UV)'].replace('@@P1@@', '%.0f' % pc[0]['poles_relative_to_G2_IV']['G1 II(1)'])
                          .replace('@@P2@@', '%.0f' % pc[0]['poles_relative_to_G2_IV']['G1 V(2b UV)'])
                          .replace('@@PSUM@@', tex_num(max(p['sum_relative'] for p in pc))))



def rows_tex(kind):
    out = []
    for gname, gl, rows in GROUPS:
        out.append(r'\section{%s}\label{sec:%s%s}' % (gname, kind, gl))
        for rname, names, special in rows:
            out.append(r'\subsection{%s}' % rname)
            if special == 'four':
                out.append(r'Only four $\rho$\textquotesingle s: Sec.~\ref{sec:four} (with the same coefficient for both logarithms).')
                continue
            for name in names:
                sub = SUB[name]
                out.append(r'\paragraph{%s%s.}' % (rname, (', ' + sub) if sub else ''))
                note = NOTE_LL.get(name)
                if note:
                    out.append(note)
                if kind == 'LL':
                    J, B = j_items(name), b_items(name)
                    out.append(r'\noindent Terms $\propto\log\frac{k^+}{\Lambda}$:' + ('' if J else ' none.'))
                    out.append(align(name, J, 'J'))
                    out.append(r'\noindent Terms $\propto\log\frac{\vee}{k^+}$:' + ('' if B else ' none.'))
                    out.append(align(name, B, 'B'))
    return '\n'.join(out)


# ------------------------------------------------------------------------------------------------ log table
def log_table():
    seen, lines = set(), []
    for s, cL, cV, R, Rx in LS._raw:
        k = LS.norm(s)
        if k in seen or (R is None and Rx is None and not (cL and cV)):
            continue
        seen.add(k)
        rhs = []
        if cL:
            rhs.append({1: '', -1: '-'}.get(cL, '%d' % cL) + r'\log\frac{k^+}{\Lambda}')
        if cV:
            rhs.append({1: '', -1: '-'}.get(cV, '%d' % cV) + r'\log\frac{\vee}{k^+}')
        rem = Rx if Rx is not None else R
        if rem is not None:
            rhs.append(rem)
        r = ' + '.join(rhs).replace('+ -', '- ').replace(r'+ \left(-\log', r'- \log') or '0'
        if r'- \log{\left(1' in r and r.endswith(r'\right)}\right)'):
            r = r[:-len(r'\right)')]
        lim = ''
        if Rx is not None and R != Rx:
            lim = r'\notag\\&\qquad\xrightarrow{\Lambda\to0}\ ' + (' + '.join(rhs[:-1] + ([R] if R else [])) or '0')
        lines.append(s + ' &= ' + r + lim)
    return '\\begin{align}\n' + ' \\notag\\\\\n'.join(lines) + '\n\\end{align}\n'


# ------------------------------------------------------------------------------------------------ four rho
def four_rho():
    L = open(os.path.join(HERE, 'input', '4rho.tex')).read()
    out = []
    for head in (r'We rewrite them (by some simplification)', r'Rewriting with some simplification, we get'):
        i = L.index(head)
        t = L[L.index(r'\begin{align}', i):L.index(r'\end{align}', i) + len(r'\end{align}')]
        t = t.replace(r'\log{\left(\frac{\vee}{\Lambda}\right)}',
                      r'\new{\left[\log{\left(\frac{k^{+}}{\Lambda}\right)}+\log{\left(\frac{\vee}{k^{+}}\right)}\right]}')
        out.append(t)
    return out


# ------------------------------------------------------------------------------------------------ three rho
def three_table(part):
    d = TC[part]['seed11_values_imag']
    rows = []
    tot = 0.0
    for labl in TC[part]['terms']:
        v = d[labl]
        cc = v['cc_added']
        tv = v['direct'] + (v['cc'] if cc else 0.0)
        tot += tv
        if labl.startswith('(B)'):
            ref = r'Eq.~\eqref{%s}' % lab('B', [k for k, h in HIDDEN.items() if h == labl.split()[-1]][0])
            name = r'(B), Group III Row %s' % labl.split()[-1]
        elif labl.startswith('4rho'):
            ref = r'Eq.~\eqref{eq:four%s}' % ('a' if 'rewrite' in labl else 'b')
            name = r'\texttt{4rho.tex}, %s Wilson lines' % ('two' if 'rewrite' in labl else 'three')
        else:
            nm, n = labl.rsplit(' #', 1)
            ref = r'Eq.~\eqref{%s}' % lab(part, nm)
            g = {'G1': 'I', 'G2': 'II', 'G3': 'III'}[nm[:2]]
            name = 'Group %s %s, line %d' % (g, nm[3:].replace('N1', r'$d^3N_1$').replace('N2', r'$d^3N_2$')
                                               .replace('N3a', r'$d^3N_{3a}$').replace('N3b', r'$d^3N_{3b}$'), int(n) + 1)
            if labl == 'G1 III N3a #1' and part == 'J':
                name += r' ($\mathcal D$)'
        ccs = '%+.4f' % v['cc'] if cc else (r'(in the term)' if cc is None else r'(self-conj.)')
        rows.append(r'%s & %s & $%+.4f$ & %s & $%+.4f$\\' % (name, ref, v['direct'], ('$%s$' % ccs) if cc else ccs, tv))
    return rows, tot


def three_section():
    out = [r'''\section{The three-\texorpdfstring{$\rho$}{rho} terms cancel}\label{sec:three}
All terms with three $\rho$\textquotesingle s of Secs.~\ref{sec:LLgI}--\ref{sec:four}, for each logarithm separately. Which terms have three $\rho$\textquotesingle s is decided numerically ($\rho\to2\rho$ multiplies them by 8).
\paragraph{Method.}
\begin{itemize}
\item The charges are integration variables, so only the totally symmetric part of $\rho\rho\rho$ can cancel against other terms. The rest (commutators, $[\rho^a(u),\rho^b(v)]=if^{abc}\rho^c(u)\delta(u-v)$) are two-$\rho$ terms; they are in Secs.~\ref{sec:two} and \ref{sec:bfkl}. The symmetric part is evaluated with classical charges.
\item In every term the integration variables are renamed so that the phase is $e^{-ik\cdot(w'-w)}$, the soft gluon is $z$ and the charges are $x,y,x'$: in Group I Row V part 1 the measured gluon of the amplitude is at $z$ (line 3) or $x$ (line 4); in Group III Rows III, V, VI and in (B) the conjugate one is at $y'$.
\item Each term is summed over the 6 placements of its three charges and, for Groups I and III, its c.c.\ is added (conjugate, $w\leftrightarrow w'$). Group II Row I is self-conjugate as written: its second line is the c.c.\ of the first. The \texttt{4rho.tex} lists are self-conjugate as written.
\item In Group I Row II $d^3N_{2a}$, lines 2 and 3, the $z$ integral is done. For a pointwise test it is written back: $\log\frac{(w-y)^2}{(x-y)^2}=-\frac1\pi\int_z\big[K(w,y,z)-K(x,y,z)\big]$, $K(a,b,z)=\mathcal K(a-z)\cdot\mathcal K(b-z)$.
\item Random points, random SU(3) adjoint Wilson lines with $U^{bc}=U^{cb}$, random charges; three configurations.
\end{itemize}
Values at one configuration, in units of $i\,\mathcal N\,e^{-ik\cdot(w'-w)}$, with the large logarithm set to 1 (all values are real in these units):''']
    for part, title in (('J', r'$\log\frac{k^+}{\Lambda}$'), ('B', r'$\log\frac{\vee}{k^+}$')):
        rows, tot = three_table(part)
        out.append(r'\paragraph{Three $\rho$\textquotesingle s $\propto$ %s (%d terms).}' % (title, len(rows)))
        out.append(r'''\begin{longtable}{@{}llrrr@{}}\toprule
term & & direct & c.c. & sum\\\midrule
''' + '\n'.join(rows) + r'''
\midrule
total & & & & $%+.1e$\\\bottomrule
\end{longtable}''' % tot)
        V = TC[part]['variants']
        out.append(r'\noindent Largest relative residual over three configurations: ' + '; '.join(
            r'%s: $%s$' % (k, tex_num(v['relative'])) for k, v in V.items()) + '.')
    out.append(r'''\paragraph{So:}
\begin{itemize}
\item \textbf{The three-$\rho$ terms $\propto\log\frac{k^+}{\Lambda}$ cancel}, with $\mathcal D$ in place of the edge-kernel term of Group I Row III $d^3N_{3a}$. With the edge kernel they do not.
\item \textbf{The three-$\rho$ terms $\propto\log\frac{\vee}{k^+}$ cancel.} The (B) terms of Group III Rows III(1) and V cancel each other, that of Row VI vanishes by itself after symmetrization; so the three-$\rho$ check does not test (B). (B) is tested by the two-$\rho$ sum, Sec.~\ref{sec:bfkl}.
\item Adding a c.c.\ to the \texttt{4rho.tex} lists would double them: they already contain it.
\item Group I Row IV and Group III Row I cancel each other line by line, and so do several pairs of lines in Group I Row III.
\end{itemize}''')
    return '\n'.join(out)


# ------------------------------------------------------------------------------------------------ reused sections
def sec(t, a, b):
    i = t.index(a)
    j = t.index(b, i)
    return t[i:j]


def once(t, a, b):
    assert t.count(a) == 1, (a, t.count(a))
    return t.replace(a, b)



def rel_sent():
    r = [abs(t['AppF_as_sent'] - t['JIMWLKxLO']) / abs(t['JIMWLKxLO']) for t in TGJ]
    c = [abs(t['AppF_corrected'] - t['JIMWLKxLO']) / abs(t['JIMWLKxLO']) for t in TGJ]
    b = [abs(t['AppE'] - t['BFKLxLO']) / abs(t['BFKLxLO']) for t in TGJ]
    return min(r), max(r), max(c), max(b)


RS0, RS1, RC, RB = rel_sent()
JIM_SENT = (r'''The JIMWLK formula of your last message is the one checked before (your App.~F): $M=(U^\dagger(x)-U^\dagger(t))(U(y)-U(t))$ in every term of $dN^{(1)}$, $dN^{(3)}$, $M=(U^\dagger(y)-U^\dagger(t))(U(x)-U(t))$ in every term of $dN^{(2)}$, $dN^{(4)}$, and $dN^{(5)}$ with an overall minus sign. (In your message the signs in front of lines 2--4 of $dN^{(5)}$ came through as bullets; they are read as in your App.~F, $-,-,+$, as below.) As sent, it differs from JIMWLK$\otimes$LO, i.e.\ from the two-$\rho$ sum of Sec.~\ref{sec:two}, by @@RS0@@--@@RS1@@ (relative; six configurations, symmetric and general $U$). Corrected as below, it agrees to $@@RC@@$.'''
            .replace('@@RS0@@', '%.2f' % RS0).replace('@@RS1@@', '%.2f' % RS1).replace('@@RC@@', tex_num(RC)))


TARGETS = sec(FEC, r'\section{The targets}', r'\section{All two-$\rho$')
TARGETS = once(TARGETS, r'\subsection{JIMWLK$\otimes$LO: your formula, corrected (M2)}',
               r'\subsection{JIMWLK$\otimes$LO: your formula, corrected (M2)}' + '\n' + JIM_SENT)
TARGETS = once(TARGETS, r'\section{The targets}\label{sec:targets}', r'''\section{The targets}\label{sec:targets}
In both formulas $\delta Y$ is the evolution interval: $\delta Y=\log\frac{k^+}{\Lambda}$ for JIMWLK and $\delta Y=\log\frac{\vee}{k^+}$ for BFKL. With $\langle\rho^a\rho^b\rangle=\delta^{ab}\langle\rho^f\rho^f\rangle/(N_c^2-1)$ and $\int_{w,w'}e^{-ik\cdot(w'-w)}$, $\mathcal N\,\delta Y\leftrightarrow\frac{\alpha_s^2\,\delta Y}{8\pi^6(N_c^2-1)}$.''')

TWO = sec(FEC, r'\section{All two-$\rho$', r'\section{BFKL: the coefficient')
TWO = once(TWO, r'\section{All two-$\rho$ terms proportional to $\log\frac{\vee}{\Lambda}$}\label{sec:two}',
           r'\section{The two-$\rho$ terms $\propto\log\frac{k^+}{\Lambda}$: JIMWLK}\label{sec:two}')
TWO = once(TWO, r'your Sec.~2: 208 / 127', r'from reordering three and four $\rho$\textquotesingle s: 208 / 127')
TWO = once(TWO, r'\subsection{Genuine terms (Sec.~1 of \texttt{Evolution\_terms.tex})}',
           r'\subsection{Genuine two-$\rho$ terms of the rows}')
TWO = once(TWO, r'To compare term by term, every',
           r'These are the two-$\rho$ terms $\propto\log\frac{k^+}{\Lambda}$ of Secs.~\ref{sec:LLgI}--\ref{sec:LLgIII} (Group I Rows II and V, Group II Rows I and IV, Group III Row III part 2), written in a common notation. To compare term by term, every')
TWO = once(TWO, r'\item Your genuine terms are equal to the corrected ones of the earlier check: Group III Row III with the sign F5, Group I Row II with the bracket of F6 (numerically, rel.\ $\le10^{-15}$).',
           r'\item Row V 2b UV is taken with the coefficient $2\log\frac{k^+}{\Lambda}$ (E1).')
TWO = once(TWO, r'So the $\log\frac{\vee}{\Lambda}$ part of Row V 2b', r'So the $\log\frac{k^+}{\Lambda}$ part of Row V 2b')
TWO = once(TWO, r'(missing in \texttt{Evolution\_terms.tex}, M1)', '')
TWO = once(TWO, r'with the two fixes found earlier:', r'with the two fixes found earlier (both still needed in this \texttt{4rho.tex}):')
TWO = once(TWO, r'\subsection{Your Sec.~2 (from four and three $\rho$'"'"'s), coloured}',
           r'\subsection{From reordering the three and four $\rho$\textquotesingle s}')
TWO = once(TWO, r'Verbatim from \texttt{Evolution\_terms.tex}; only the colours are added. Every row is correct (Sec.~\ref{sec:checks}).',
           r'''The commutator terms of the three-$\rho$ terms of Sec.~\ref{sec:three} and of the four-$\rho$ terms. They are your lists of \texttt{Evolution\_terms.tex}, Sec.~2, which were checked row by row against these rows (0 differing kernel monomials in every row, \texttt{Final\_evolution\_check.pdf}), here with $\log\frac{\vee}{\Lambda}\to\new{\log\frac{k^+}{\Lambda}}$ and coloured.''')
TWO = TWO.replace(r'\log\Big[\frac{\vee}{\Lambda}\Big]', r'\new{\log\Big[\frac{k^{+}}{\Lambda}\Big]}')
TWO = once(TWO, r'This is $\mathcal N\log\frac{\vee}{\Lambda}\times\mathcal J$', r'This is $\mathcal N\log\frac{k^+}{\Lambda}\times\mathcal J$')

BF = sec(FEC, r'\section{BFKL: the coefficient', r'\section{Every term that depends')
i = BF.index(r'\paragraph{Group I, Row II')
BF = r'''\section{The two-$\rho$ terms $\propto\log\frac{\vee}{k^+}$: BFKL}\label{sec:bfkl}
The coefficient of $\log\frac{\vee}{k^+}$ of every row is its $\log\frac{\vee}{k^+}$ list of Secs.~\ref{sec:LLgI}--\ref{sec:LLgIII}.
\begin{itemize}
\item It is zero for Group I Rows II (part 1) and V, and for $\mathcal D$.
\item For the four-$\rho$ terms, their three-$\rho$ remainders and the \texttt{4rho.tex} lists it equals their $\log\frac{k^+}{\Lambda}$ terms.
\item Group III Rows III, V, VI: the terms (B); Group II Row IV: $H_4$.
\end{itemize}
The two-$\rho$ terms (genuine and from reordering the three- and four-$\rho$ terms), row by row (orange / purple: 168 / 162):
''' + BF[i:]
BF = BF.replace(r'\bigg|^{\log\vee}', r'\bigg|^{\log\frac{\vee}{k^+}}')
BF = BF.replace(r' the $\alpha$ part, fix F4', '').replace(r'; with the sign of $p^+>k^+$, fix P1', '')
BF = BF.replace(r'(three $\rho$; )', r'(three $\rho$)')
assert r'\log\vee' not in BF.replace(r'\bigg|^{\log\frac', ''), [m.start() for m in re.finditer(r'\\log\\vee', BF)][:3]

# ------------------------------------------------------------------------------------------------ your files, step by step
sys.path.insert(0, os.path.normpath(os.path.join(HERE, '..', 'rows_check')))
import pe  # noqa: E402
pe.LINES = {g: open(os.path.join(HERE, 'input', '%s.tex' % g)).read().split('\n') for g in ('Group_I', 'Group_II', 'Group_III')}
import check_all as CA  # noqa: E402
import eqnum as EN  # noqa: E402


def step_rows():
    rows = []
    for g, gname in (('Group_I', 'I'), ('Group_II', 'II'), ('Group_III', 'III')):
        for name, kind, e1, e2, opt in CA.STEPS[g]:
            res = STEPS[g][name]
            if kind == 'plus':
                val = r'$%s$ / $%s$' % (tex_num(max(r['rel_L1e-06'] for r in res)), tex_num(max(r['rel_L1e-07'] for r in res)))
            elif kind == 'int':
                k = [x for x in res[0] if x.startswith('rel')][0]
                val = '$%s$' % tex_num(max(r[k] for r in res))
            else:
                val = '$%s$' % tex_num(max(r['rel'] for r in res))
            mark = r' \textbf{(M2)}' if 'after the z integral' in name else ''
            disp = (name.replace('Eq. (111)', 'Eq. (1.16)').replace('(q+ = k+ - p+)', '($q^+=k^+-p^+$)')
                    .replace('P_gg', r'$\mathcal P_{gg}$').replace('f f U U = N_c U', '$ffUU=N_cU$')
                    .replace('p+ ', '$p^+$ ').replace('xi', r'$\xi$').replace("w'", "$w'$").replace('->', r'$\to$'))
            kname = {'int': 'integrals', 'plus': r'plus, $\Lambda\to0$', 'pt': r'fixed $p^+$'}[kind]
            eqs = lambda rr: ', '.join(EN.eq(g, a, b) for a, b in rr)
            rows.append(r'%s & %s%s & %s $\to$ %s & %s & %s\\' % (gname, disp, mark, eqs(e1), eqs(e2), kname, val))
    return '\n'.join(rows)


UVJ = STEPS['Group_I']['Row V 2b: UV term'][0]


def four_numbers():
    """equation numbers of 4rho.tex as it compiles: {first line of the align: '(s.n)--(s.m)'}"""
    L = open(os.path.join(HERE, 'input', '4rho.tex')).read().split('\n')
    sec_, n, out, i = 0, 0, {}, 0
    while i < len(L):
        s_ = L[i]
        if re.match(r'\s*\\section\{', s_):
            sec_ += 1
            n = 0
        if r'\begin{align}' in s_ and not s_.lstrip().startswith('%'):
            j = i
            while r'\end{align}' not in L[j]:
                j += 1
            body = re.sub(r'(?<!\\)%[^\n]*', '', '\n'.join(L[i:j + 1]))
            body = body[body.index(r'\begin{align}') + 13:body.rindex(r'\end{align}')]
            nums = []
            for r in re.split(r'\\\\', body):
                if r.strip() and not re.search(r'\\(notag|nonumber)', r):
                    n += 1
                    nums.append('%d.%d' % (sec_, n))
            out[i + 1] = ('(%s)' % nums[0]) if len(nums) == 1 else ('(%s)--(%s)' % (nums[0], nums[-1]) if nums else '(none)')
            i = j + 1
        else:
            i += 1
    return out


FOURNUM = four_numbers()
FOURWHAT = {66: r'four $\rho$\textquotesingle s, Eq.~(4rho)', 73: r'its symmetrized sum $=0$', 77: r'four Wilson lines, rewritten',
            83: r'four Wilson lines, $=0$', 90: r'three $\rho$: first form', 95: r'three $\rho$', 100: r'three $\rho$',
            107: r'three $\rho$', 114: r'three $\rho$', 121: r'three $\rho$', 128: r'three $\rho$, $=0$',
            134: r'two $\rho$, first part', 143: r'two $\rho$, second part (1)', 166: r'two $\rho$, second part (2)'}


def four_rows():
    rows = []
    for e in FF['equations']:
        ln = e['line_new']
        what = FOURWHAT.get(ln, '')
        if 'text_identical' in e:
            res = r'identical text (without the log bracket); logs %d/%d' % (e['logs_new'], e['logs_old']) if e['text_identical'] \
                else r'\textbf{text differs}'
        elif e.get('both_zero'):
            res = r'both $0$'
        elif e['rel'] < 1e-12:
            res = r'$%s$' % tex_num(e['rel'])
        else:
            q = e['ratio'][0][0]
            res = r'\textbf{ratio %.4f} $=1/\log\frac{\vee}{\Lambda}$: the log bracket is missing' % q
            assert abs(q * e['log_vee_over_Lambda'] - 1) < 1e-6
        rows.append(r'{}%s & %s & %s\\' % (FOURNUM.get(ln, '?'), what, res))      # {}: a row may not start with '(' after \midrule
    return '\n'.join(rows)


FC4 = max(abs(r['residual']) / r['scale'] for r in FF['four_rho_cancel'])
TL = FF['two_rho_lists']


# ------------------------------------------------------------------------------------------------ exact identities
EC_SEEDS = [ECJ[k] for k in ECJ if k != 'fourier']
EC_ID = max(r[n]['identity'] for r in EC_SEEDS for n in r)
EC_B = max(r[n]['GzFT_vs_B'] for r in EC_SEEDS for n in r if n != 'G2 IV')
EC_T0 = max(r['G2 IV']['term0_plus_line2'] for r in EC_SEEDS)
FOU = ECJ['fourier']


def fourier_line():
    j1 = [(k.split('q=')[1].split()[0], v) for k, v in FOU.items() if k.startswith('int J1')]
    lg = [(k.split('q=')[1].split()[0], v) for k, v in FOU.items() if k.startswith('FT log')]
    a = '; '.join(r'$q=%s$: %s' % (q, ', '.join('%.4f' % x for x in v)) for q, v in j1)
    b = '; '.join(r'$q=%s$: %s $\to$ %.5f' % (q, ', '.join('%.4f' % x for x in v[:-1]), v[-1]) for q, v in lg)
    return (r'$q\int_0^\infty J_1(qr)\,e^{-ar}dr$ at $a=10^{-2},3\times10^{-3},10^{-3}$ (must tend to 1): ' + a +
            r'. $2\pi\int_0^\infty rJ_0(qr)\log r^2\,e^{-ar^2}dr\big/(-\frac{4\pi}{q^2})$ at $a=4,2,1\times10^{-3}$, extrapolated to $a=0$ (must be 1): ' + b + '.')


# ------------------------------------------------------------------------------------------------ the document
N_J = len(TC['J']['terms'])
N_B = len(TC['B']['terms'])
WJ = TC['J']['variants']['final']['relative']
WB = TC['B']['variants']['final']['relative']
WE = TC['J']['variants']['edge kernel (as written) instead of D']['relative']
WNB = TC['B']['variants']['without (B)']['relative']
PSUM = max(p['sum_relative'] for p in pc)
RVA = max(v['A'] for k, v in RV.items() if k != 'G1 V(2b UV)')
RVB = max(v['B'] for v in RV.values())
RVC = max(v['C'] for v in RV.values())
SLW = max(max(v.values()) for v in SL.values())

DOC = PREAMBLE + r'''
\title{Your Group I, II, III and \texttt{4rho.tex}: every step checked;\\ the terms $\propto\log\frac{k^+}{\Lambda}$ (JIMWLK) and $\propto\log\frac{\vee}{k^+}$ (BFKL), row by row;\\ three- and four-$\rho$ cancellation; the two-$\rho$ sums and your two formulas}
\date{}
\begin{document}
\maketitle
\noindent Notation as in your files: $\mathcal K(X)=X/X^2$, $K(a,b,z)=\mathcal K(a-z)\cdot\mathcal K(b-z)$, $U^{bc}=U^{cb}$,
\begin{equation*}
\mathcal N=\frac{1}{(2\pi)^3}\frac{g^4}{16\pi^5}\frac{1}{k^+}.
\end{equation*}
Groups I and III are taken with their complex conjugate, Group II without, in every test (three $\rho$, two $\rho$, poles). \new{Blue}: new or corrected; \old{red}: removed. Only $\vee>k^+$ and $p^+>\Lambda$ are used; $\Lambda\to0$ wherever possible, $\vee$ exact, and every step is an exact identity. The remaining terms, the running coupling, DGLAP and the finite terms are in the companion file \texttt{Remaining\_RC\_DGLAP\_finite.pdf}.

\section{Answer}\label{sec:answer}
\paragraph{Your files.} Every step of \texttt{Group\_I.tex}, \texttt{Group\_II.tex}, \texttt{Group\_III.tex} and every equation of \texttt{4rho.tex} was checked numerically (Sec.~\ref{sec:files}). E1 (Group I Row V 2b UV) is now correct in your file. What is left to change:
\begin{description}
\item[M2] Group II Row IV, Eq.~(1.16): phase of the second term (typo; the next equation is right).
\item[4$\rho$] \texttt{4rho.tex}, Eq.~@@EQ134@@ (two $\rho$, first part): the factor $\new{\big[\log\frac{k^+}{\Lambda}+\log\frac{\vee}{k^+}\big]}$ is missing (the earlier file had $\log\frac{\vee}{\Lambda}$); F1 and F2 are still needed in the two-$\rho$ lists.
\item[F3] Group I Row III $d^3N_{3a}$: the logarithmic line with the edge kernel is replaced by $\mathcal D$. Without it the three-$\rho$ terms $\propto\log\frac{k^+}{\Lambda}$ do not cancel ($@@WE@@$).
\item[C5] Group III Row III part 2: the measure is $\int_{w,x',y',z,x}$ (M4).
\item[M3, M5, M6] a prefactor sign in an unnumbered first line, a measure, typos (Sec.~\ref{sec:files}).
\item[M2'] Your JIMWLK formula: $M$ at the wrong points and the sign of $dN^{(5)}$ (Sec.~\ref{sec:targets}). Your BFKL formula is correct.
\end{description}
\paragraph{The split.} Every row is split exactly, term by term,
\begin{equation}
\frac{d^3N_{\rm row}}{d^3k}=\log\frac{k^+}{\Lambda}\;\mathfrak T_{\rm row}+\log\frac{\vee}{k^+}\;\mathfrak P_{\rm row}+R_{\rm row},
\label{eq:split}
\end{equation}
with $\log\frac{\vee}{\Lambda}=\log\frac{k^+}{\Lambda}+\log\frac{\vee}{k^+}$ and the identities of Sec.~\ref{sec:logs}. $R_{\rm row}$ (the remaining terms) is finite for $\vee\to\infty$ and $\Lambda\to0$, which fixes $\mathfrak T$ and $\mathfrak P$ uniquely. Each of the @@NCHK@@ pointwise-evaluable terms equals its three parts to $@@ID@@$.
\paragraph{Results.}
\begin{itemize}
\item \textbf{Three $\rho$\textquotesingle s $\propto\log\frac{k^+}{\Lambda}$} (@@NJ@@ terms): \textbf{zero} after symmetrization, $@@WJ@@$ (Sec.~\ref{sec:three}).
\item \textbf{Three $\rho$\textquotesingle s $\propto\log\frac{\vee}{k^+}$} (@@NB@@ terms): \textbf{zero}, $@@WB@@$.
\item \textbf{Four $\rho$\textquotesingle s} (\texttt{4rho.tex}, with $\big[\log\frac{k^+}{\Lambda}+\log\frac{\vee}{k^+}\big]$): \textbf{zero} after symmetrization, $@@FC4@@$ (Sec.~\ref{sec:four}). Their three- and two-$\rho$ remainders enter both lists.
\item \textbf{Two $\rho$\textquotesingle s $\propto\log\frac{k^+}{\Lambda}$ $=$ JIMWLK$\otimes$LO}, exactly: 13 kernel structures, 154 independent colour structures, 0 differences (Sec.~\ref{sec:two}). It is your JIMWLK formula $\sum_{i=1}^5dN^{(i)}$ after the correction M2' (Sec.~\ref{sec:targets}).
\item \textbf{Two $\rho$\textquotesingle s $\propto\log\frac{\vee}{k^+}$ $=$ BFKL$\otimes$LO}, exactly: 9 kernel structures, 109 colour structures, 0 differences (Sec.~\ref{sec:bfkl}). It is your BFKL formula, Eq.~(LOBFKL), as you wrote it.
\end{itemize}
Therefore
\begin{equation}
\frac{d^3N_{\rm LL}}{d^3k}=\mathcal N\Big[\log\frac{k^+}{\Lambda}\;{\rm JIMWLK}\otimes{\rm LO}+\log\frac{\vee}{k^+}\;{\rm BFKL}\otimes{\rm LO}\Big]
=-\mathcal N\log\frac{\Lambda}{k^+}\;{\rm JIMWLK}\otimes{\rm LO}+\mathcal N\log\frac{\vee}{k^+}\;{\rm BFKL}\otimes{\rm LO}.
\label{eq:LL}
\end{equation}
In your formulas $\delta Y=\log\frac{k^+}{\Lambda}$ for JIMWLK (the target evolves from $y_T=\log\Lambda$ to the gluon) and $\delta Y=\log\frac{\vee}{k^+}$ for BFKL (the projectile evolves from the gluon to $y_P=\log\vee$). Written with $\log\frac{\Lambda}{k^+}$, the JIMWLK term has the coefficient $-{\rm JIMWLK}\otimes{\rm LO}$.

\paragraph{Where the BFKL terms are.} In the explicit logarithms of the rows, and, in Group III Rows III, V, VI and Group II Row IV, inside line 1 (the $p^+$, resp.\ $\xi$, integral). There three exact identities bring them out: the subtraction of the plus prescription gives $-f(0)\log\frac{\vee}{k^+}$, which cancels the $\log\frac{\vee}{k^+}$ of line 2 exactly (I2); the Fourier transform of the $y'$ integral at $y'=z$ gives the terms (B) (I3); the Fourier transform of $\log(y'-y)^2$ gives $H_4$ (I4). Nothing is expanded or dropped.

\section{Your files, step by step}\label{sec:files}
\subsection{Group I, II, III}
Every step that turns one of your equations into the next is checked numerically with your equations as written (read from your \LaTeX): the $p^+$, $q^+$, $\xi$ integrals, the $w$, $w'$ integrals over $\delta^{(2)}$, the Wilson-line reductions, the splits (Group I Row V into 2a, 2b and 2b into UV and finite; Group II Row I part 2 into 2a--2f) and every plus prescription. Pointwise in all transverse coordinates, random points, random adjoint SU(3) Wilson lines ($U^{bc}=U^{cb}$), random charges, two configurations; remaining longitudinal integrals numerically. Equation numbers as your files compile. ``plus'': the plus-prescription steps hold up to $O(\Lambda)$; the two numbers are at $\Lambda=10^{-6}$ and $10^{-7}$.
{\small
\begin{longtable}{@{}lp{0.36\textwidth}p{0.22\textwidth}lr@{}}\toprule
G & step & Eqs. & kind & rel.\ difference\\\midrule
@@STEPS@@
\bottomrule
\end{longtable}}
\noindent Group I Row V 2b, UV term (your E1 bracket): coefficient of the pole bracket vs.\ the exact $p^+$ integral, relative difference $@@UV@@$.
\paragraph{Mistakes left in the Group files.}
\begin{enumerate}
\item[\textbf{M2}] \textbf{Group II, Row IV, Eq.~(1.16):} in the second term the phase must be $e^{-ik\cdot(y'-y)\frac{p^++k^+}{k^+}}$, not $e^{-ik\cdot(y'-y)\frac{p^+}{(p^++k^+)}}$. A typo: Eq.~(1.14) before and Eq.~(1.17) after have the right phase. The step fails by $@@M2A@@$ as written and agrees to $@@M2B@@$ with the phase corrected.
\item[\textbf{M3}] \textbf{Group III, Row III, first (unnumbered) line:} the prefactor is $\new{-}\frac{4}{(2\pi)^3}$, as in Rows V and VI ($-\frac{2}{(2\pi)^3}$ per term). Eqs.~(1.5)--(1.10) below it are right.
\item[\textbf{M4}] \textbf{Group III, Row III, second part, Eqs.~(1.8)--(1.10) (C5):} the measure is $\int_{w,w',x',y',z,x}$ in Eq.~(1.8) and $\int_{w,x',y',z,x}$ after the $w'$ integral (Eqs.~(1.9), (1.10)).
\item[\textbf{M5}] \textbf{Group I, Row II, Eq.~(1.6):} $\int_{w,w'}\to\int_{w,w',x'}$.
\item[\textbf{M6}] \textbf{Typos that do not change anything:} Group~II Eq.~(1.13): $k\{+\}\to k^+$ in the phase. Group~I Eq.~(1.9): the opening \texttt{\textbackslash Bigg\textbackslash\{} is never closed; Eqs.~(1.19), (1.20), (1.22), (1.24): a closing \texttt{\textbackslash Bigg\textbackslash\}} closes nothing. Group~I Eqs.~(1.16), (1.17), (1.30): $\int_{x,z,y,x',y'w'}\to\int_{x,z,y,x',y',w'}$; Eqs.~(1.1), (1.2): a stray comma. Group~II Row~I: ``We write it in two parts'' twice.
\end{enumerate}
\noindent The final forms of all rows (the ones used below) are the same as in the files of the earlier checks, apart from your corrected E1 bracket: their $\log\frac{k^+}{\Lambda}$ parts, $\log\frac{\vee}{k^+}$ parts and remainders agree row by row and pointwise ($@@SLW@@$, \texttt{same\_lists.py}). Your earlier \texttt{BFKL.tex} $+$ \texttt{Other\_terms.tex} are these rows minus their JIMWLK parts ($\le@@RVA@@$, i.e.\ $O(\Lambda)$ at $\Lambda=10^{-9}$; \texttt{rows\_vs\_files.py}).

\subsection{\texttt{4rho.tex}}
Each equation against the earlier, checked \texttt{4rho.tex}, pointwise (random charges, points and Wilson lines; $\log\frac{\vee}{\Lambda}$ there, $\big[\log\frac{k^+}{\Lambda}+\log\frac{\vee}{k^+}\big]$ here). The two equations written as ``prefactor $\times\rho\rho\rho\rho[A]+[B]$'' are compared as text.
{\small
\begin{longtable}{@{}llp{0.55\textwidth}@{}}\toprule
Eq. & & new vs.\ earlier file\\\midrule
@@FOUR@@
\bottomrule
\end{longtable}}
\begin{itemize}
\item \textbf{The four-$\rho$ terms cancel} after symmetrization over the four charges: residual / size $@@FC4@@$ (three configurations; \texttt{three\_rho/four\_check.py}, whose transcription equals your Eq.~(4rho)).
\item \textbf{Two $\rho$, first part, Eq.~@@EQ134@@:} the prefactor is $\frac{g^4}{8\pi^3}\frac{1}{2\pi^2}$; it must be $\frac{g^4}{8\pi^3}\frac{\new{\left[\log\frac{k^+}{\Lambda}+\log\frac{\vee}{k^+}\right]}}{2\pi^2}$, as in the second part. As written it is a log-free (finite) term, and the evolution would miss it.
\item \textbf{The two-$\rho$ lists} vs.\ the checked lists (Final\_evolution\_check.pdf, Eq.~(fourlists)), charges symmetrized, with c.c., three configurations: with the log restored and F1, F2: $@@TLOK@@$; with the log restored and F2 only: @@TLF2@@; without the log (F2 only): @@TLNO@@. F2 is always needed: as written, the index $b'$ occurs four times in one term.
\item[F1] first part: drop \old{$-8N_c\,\mathcal K(x'-w')\cdot\mathcal K(x'-w)\,K(y,x',z)\,\rho\rho\,UU$} (the term starting with $-8$);
\item[F2] second part, $dN^{(1)}$, term 2: $f^{abc}f^{de\new{g}}\rho^{b'}(x')\rho^{c}(y)U^{a\new{g}}(y)U^{db'}(w')U^{eb}(z)$.
\end{itemize}

\section{How the logarithms are split: exact identities}\label{sec:logs}
\paragraph{(I1) Explicit logarithms.} A term carries at most one logarithmic factor. Your brackets:
\begin{align}
\left[\log{\left(\frac{k^{+}}{\Lambda}\right)}+\log{\left(\frac{\vee}{k^{+}}-1\right)}\right]&=\log\frac{k^+}{\Lambda}+\log\frac{\vee}{k^+}+\log\Big(1-\frac{k^+}{\vee}\Big),\notag\\
\left[\log{\left[\frac{\vee}{k^{+}}\right]}+\log{\left[\frac{k^{+}}{\Lambda}\right]}\right]&=\log\frac{k^+}{\Lambda}+\log\frac{\vee}{k^+}.\notag
\end{align}
All other logarithms ($A=(x-z)^2$ etc.; the remainder exactly and, where it depends on $\Lambda$, its limit):
@@LOGTABLE@@
Logarithms of $\epsilon$, $\mu$, $\log(1-\frac{k^+}{\vee})$ and $\log(1-\frac{\Lambda}{k^+})$ are finite and are not split. Three special cases:
\begin{itemize}
\item Group I Row II part 1: $\big[\frac{11}{3}-2\log\frac{k^+}{\Lambda}+\log(1-\frac{k^+}{\vee})\big]$ times the pole bracket; only $-2\log\frac{k^+}{\Lambda}$ is large.
\item Group I Row V 2b UV: $\big[2\log\frac{k^+}{\Lambda}+2\log(1-\frac{\Lambda}{k^+})-\frac{11}{6}\big]$.
\item Group I Row V 2b finite, line 6: $\big[-\frac{11}{6}+2\log\frac{k^+-\Lambda}{\Lambda}\big]$.
\end{itemize}

\paragraph{(I2) The subtraction of the plus prescription.}\label{sec:exact} With your $\left[\frac{1}{p^+}\right]_+f=\frac{f(p^+)-f(0)}{p^+}$, define $\left[\frac{k^+}{p^+(p^++k^+)}\right]_+f\equiv\frac{k^+[f(p^+)-f(0)]}{p^+(p^++k^+)}$. Since $\frac{1}{p^+}=\frac{1}{p^++k^+}+\frac{k^+}{p^+(p^++k^+)}$ and $\int_0^{\vee-k^+}\frac{dp^+}{p^++k^+}=\log\frac{\vee}{k^+}$,
\begin{equation}
\int_0^{\vee-k^+}dp^+\left[\frac{1}{p^+}\right]_+f(p^+)=\int_0^{\vee-k^+}dp^+\left(\left[\frac{k^+}{p^+(p^++k^+)}\right]_++\frac{1}{p^++k^+}\right)f(p^+)-f(0)\log\frac{\vee}{k^+}.
\label{eq:I2}
\end{equation}
For $\xi$, with your $\frac{d\xi}{\xi^2}\frac{\xi}{(1-\xi)_+}g=\frac{d\xi}{\xi(1-\xi)}[g(\xi)-g(1)]$, define $\frac{d\xi}{\xi^2}\frac{\xi^2}{(1-\xi)_+}g\equiv\frac{d\xi}{1-\xi}[g(\xi)-g(1)]$. Then $\frac{1}{\xi(1-\xi)}=\frac{1}{\xi}+\frac{1}{1-\xi}$ gives
\begin{equation}
\int_{\frac{k^+}{\vee}}^1\frac{d\xi}{\xi^2}\frac{\xi}{(1-\xi)_+}g=\int_{\frac{k^+}{\vee}}^1\frac{d\xi}{\xi^2}\Big(\xi+\frac{\xi^2}{(1-\xi)_+}\Big)g-g(1)\log\frac{\vee}{k^+}.
\label{eq:I2xi}
\end{equation}
Your line 2 is $f(0)\big[\log\frac{k^+}{\Lambda}+\log(\frac{\vee}{k^+}-1)\big]$ (resp.\ $g(1)[\dots]$): it came from this subtraction. So its $\log\frac{\vee}{k^+}$, (i), and $-f(0)\log\frac{\vee}{k^+}$, (ii), cancel exactly, poles included; $\log\frac{k^+}{\Lambda}$ (JIMWLK) and $\log(1-\frac{k^+}{\vee})$ (remaining) are left.

\paragraph{(I3) The $y'$ integral at $y'=z$.} In line 1 of Group III Rows III(1), V, VI take the $-\frac{1}{k^+}\delta_{im}\delta_{jk'}$ term; in Row III(2) the $-\frac{p^+}{k^{+2}}$ term. Its coefficient tends to a constant at large $p^+$. Write
\begin{equation*}
G^{k'}(y')=\frac{(x'-y')^{k'}}{(x'-y')^2}\times\big[\text{colour bracket with }U(y')\big]=G^{k'}(z)+\big[G^{k'}(y')-G^{k'}(z)\big].
\end{equation*}
With $\int_u e^{-iq\cdot u}\frac{u^m}{u^2}=-2\pi i\frac{q^m}{q^2}$ and $q=k\frac{p^++k^+}{k^+}$, exactly,
\begin{equation}
\int_{y'}e^{-ik\cdot(y'-w)}e^{-ik\cdot\frac{p^+}{k^+}(y'-z)}\frac{(y'-z)^m}{(y'-z)^2}\,G^{k'}(z)=\frac{k^+}{p^++k^+}\int_{w'}e^{-ik\cdot(w'-w)}\frac{(w'-z)^m}{(w'-z)^2}\,G^{k'}(z).
\label{eq:I3}
\end{equation}
For Rows III(1), V, VI: $\int_0^{\vee-k^+}\frac{dp^+}{2\pi}\big(-\frac{1}{k^+}\big)\frac{k^+}{p^++k^+}=-\frac{1}{2\pi}\log\frac{\vee}{k^+}$, exactly. This is (B), (iii), with nothing left over. The $G(y')-G(z)$ part stays in the remaining terms.

\paragraph{(I5) Row III(2).} With $A=(x-z)^2$, $B=(x-w)^2$,
\begin{equation}
-\frac{p^+}{k^{+2}}\,\frac{k^+}{p^++k^+}\,\frac{1}{p^+A+k^+B}=-\frac{1}{k^+A\,(p^++k^+)}+\frac{B}{A\,(p^++k^+)(p^+A+k^+B)}.
\label{eq:I5}
\end{equation}
The first term gives (B). The second falls like $\frac{1}{p^{+2}}$ and stays in the remaining terms.

\paragraph{(I4) Group II Row IV.} In the $\frac{1-\xi}{\xi}$ term of line 1 write $H(y,y')=H(y,y)+[H(y,y')-H(y,y)]$, $H(y,y')=\frac{(x-y)\cdot(x'-y')}{(x-y)^2(x'-y')^2}[U^{e'c}(y')-U^{e'c}(x')][U^{ec}(y)-U^{ce}(x)]\rho^{e'}(x')\rho^e(x)$. The constant part of the pole bracket gives $\delta^{(2)}(k/\xi)=0$ for $k\neq0$, and $\int_{y'}e^{-ik\cdot(y'-y)\frac{1}{\xi}}\log\frac{(y'-y)^2\mu^2}{4}=-\frac{4\pi\xi^2}{k^2}$. Also $\int_{\frac{k^+}{\vee}}^1d\xi\,\frac{1-\xi}{\xi}=\log\frac{\vee}{k^+}-1+\frac{k^+}{\vee}$ and $\int_{w,w'}e^{-ik\cdot(w'-w)}\mathcal K(w-z)\cdot\mathcal K(w'-z)=\frac{4\pi^2}{k^2}$. So the $H(y,y)$ part is, exactly,
\begin{equation}
\frac{g^2}{4\pi^2}\frac{N_c}{(2\pi)^3}\frac{g^2}{2\pi^2}\frac{1}{k^+}\,\frac{d_\perp}{2}\,\frac{4\pi}{k^2}\Big[\log\frac{\vee}{k^+}-1+\frac{k^+}{\vee}\Big]\int_{x,x',z}H(z,z)
=\Big[\log\frac{\vee}{k^+}-1+\frac{k^+}{\vee}\Big]\,H_4\Big|_{\log\frac{\vee}{k^+}\to1},
\label{eq:I4}
\end{equation}
at $d_\perp=2$ (no pole multiplies it). Its $\log\frac{\vee}{k^+}$ term is $\log\frac{\vee}{k^+}H_4$, (iii) of Group II Row IV.

\paragraph{The split is complete.} After these steps every remaining $p^+$ ($\xi$) integrand falls at least like $\frac{1}{p^{+2}}$ after the transverse integrals ($G(y')-G(z)$ and $H(y,y')-H(y,y)$ vanish at the coincident point; the $\frac{1}{p^++k^+}$ and $\left[\frac{k^+}{p^+(p^++k^+)}\right]_+$ terms carry the $\frac{k^+}{p^++k^+}$ of (I3)). So the remaining terms stay finite for $\vee\to\infty$, and the coefficient of $\log\frac{\vee}{k^+}$ is unique; $\vee$ stays exact in the remaining terms.

\paragraph{Numerical checks} (\texttt{bfkl\_rows/exact\_check.py}; two configurations, $\vee=20k^+$, the $p^+$ and $\xi$ integrals done numerically, pointwise in the transverse points; line 1 of these rows is the same in your new files):
\begin{itemize}
\item line 1 $=$ (rewritten line 1) $+$ ($G(y')-G(z)$ term) $+$ ($-f(0)\log\frac{\vee}{k^+}$) $+$ ($G(z)$ term), Group III Rows III(1), III(2), V, VI; line 1 $=$ (rewritten line 1) $+$ ($-g(1)\log\frac{\vee}{k^+}$), Group II Row IV: $\le@@ECID@@$;
\item the $G(z)$ term in the form \eqref{eq:I3}, $p^+$ integral done, $=$ (B) ($+$ the second term of \eqref{eq:I5} in Row III(2)): $\le@@ECB@@$;
\item Group II Row IV: (i) $+$ (ii) $=0$: $@@ECT0@@$;
\item the Fourier transforms: @@FOU@@
\end{itemize}

@@ROWS@@

\section{Four \texorpdfstring{$\rho$}{rho}'s}\label{sec:four}
Group I Row I, Group II Rows II and III and Group III Row IV have only four $\rho$\textquotesingle s: your \texttt{4rho.tex}. Every four-$\rho$ term carries $\big[\log\frac{k^+}{\Lambda}+\log\frac{\vee}{k^+}\big]$: the same coefficient for both logarithms.
\begin{itemize}
\item The four-$\rho$ terms cancel after symmetrization (Sec.~\ref{sec:files}: $@@FC4@@$).
\item Their three-$\rho$ remainders, your final forms, enter both three-$\rho$ lists of Sec.~\ref{sec:three}:
\end{itemize}
@@FOUR_A@@
@@FOUR_B@@
\begin{itemize}
\item Their two-$\rho$ remainders (your ``Two $\rho$ contribution'', with the log bracket restored in the first part and F1, F2) enter both two-$\rho$ sums, Secs.~\ref{sec:two} and \ref{sec:bfkl}.
\end{itemize}

@@THREE@@

@@TWO@@

@@BF@@

@@TARGETS@@
\paragraph{Your BFKL formula} (Eq.~(LOBFKL) of your message) is correct as written: with $K_D=K(x,x,z)-2K(x,y,z)+K(y,y,z)$, $\delta Y=\log\frac{\vee}{k^+}$ and the singlet projection above, it equals BFKL$\otimes$LO, i.e.\ the two-$\rho$ sum of Sec.~\ref{sec:bfkl}, to $@@RB@@$ (six configurations).

\section{Poles}\label{sec:poles}
With $E(r)=\frac{2}{\epsilon}+2\gamma_E+\log\frac{r^2\mu^2}{4}$:
\begin{itemize}
\item The poles $\times\log\frac{k^+}{\Lambda}$ of Group I Row II, Group I Row V 2b UV and Group II Row IV are in the ratio $@@P1@@:@@P2@@:1$ and cancel ($@@PSUM@@$, three configurations; \texttt{ll\_checks.py poles}).
\item The $\log\frac{\vee}{k^+}$ list has two terms with a pole, (i) and (ii) of Group II Row IV; they cancel exactly. So the BFKL terms are pole free.
\end{itemize}

\section{Checks}\label{sec:checks}
Scripts in \texttt{full\_check/} unless stated:
\begin{center}\small
\begin{tabular}{@{}p{0.62\textwidth}p{0.32\textwidth}@{}}\toprule
check & result\\\midrule
every step of your Group files (\texttt{step\_check.py}) & Sec.~\ref{sec:files}\\
every equation of \texttt{4rho.tex} vs.\ the checked file (\texttt{four\_file\_check.py}) & Sec.~\ref{sec:files}\\
every evaluable term $=$ its $\log\frac{k^+}{\Lambda}$, $\log\frac{\vee}{k^+}$ and remaining parts, exactly in $\Lambda$, $\vee$ (\texttt{ll\_checks.py split}; @@NCHK@@ terms, 3 configurations, 2 sets of $k^+,\vee,\Lambda$) & $@@ID@@$\\
remaining parts: exact remainder $\to$ its $\Lambda\to0$ form & $@@L0@@$ at $\Lambda=10^{-9}$\\
$\frac{1}{\Lambda}$ lines of Group I Row III add up to zero; your two red lines & $@@POW@@$; $0$\\
the three parts equal those of the earlier files, row by row (\texttt{same\_lists.py}) & $@@SLW@@$\\
earlier \texttt{BFKL.tex} $+$ \texttt{Other\_terms.tex} $=$ new rows $-$ JIMWLK; their BFKL list $=$ this one; remaining parts (\texttt{rows\_vs\_files.py}) & $@@RVA@@$; $@@RVB@@$; $@@RVC@@$\\
exact identities of line 1 (\texttt{bfkl\_rows/exact\_check.py}) & $\le@@ECID@@$; (B): $\le@@ECB@@$\\
poles $\times\log\frac{k^+}{\Lambda}$ cancel (\texttt{ll\_checks.py poles}) & $@@PSUM@@$\\
three $\rho$\textquotesingle s $\propto\log\frac{k^+}{\Lambda}$ / $\propto\log\frac{\vee}{k^+}$ cancel (\texttt{ll\_checks.py three}) & $@@WJ@@$ / $@@WB@@$\\
the same with the edge kernel instead of $\mathcal D$ / without (B) & $@@WE@@$ / $@@WNB@@$\\
four $\rho$\textquotesingle s cancel (\texttt{four\_file\_check.py}) & $@@FC4@@$\\
two-$\rho$ lists of \texttt{4rho.tex} (log restored, F1, F2) vs.\ the checked ones & $@@TLOK@@$\\
two $\rho$\textquotesingle s $\propto\log\frac{k^+}{\Lambda}$ vs JIMWLK$\otimes$LO; $\propto\log\frac{\vee}{k^+}$ vs BFKL$\otimes$LO, by colour structure (\texttt{final\_check/render.py}) & 0 of 154; 0 of 109 differ\\
your JIMWLK formula: as sent / corrected; your BFKL formula (\texttt{final\_check/targets\_check.py}) & @@RS0@@--@@RS1@@ / $@@RC@@$; $@@RB@@$\\
\bottomrule
\end{tabular}
\end{center}
\end{document}
'''
M2S = STEPS['Group_II']['Row IV: after the z integral -> simplifying']
M2C = STEPS['Group_II']['Row IV: the same, phase of the 2nd term of Eq. (111) corrected']
REP = {'@@NCHK@@': str(CS['n_checked']), '@@ID@@': tex_num(CS['worst_identity']), '@@L0@@': tex_num(CS['worst_Lambda_to_0']),
       '@@POW@@': tex_num(CS['power_terms_sum_rel']),
       '@@NJ@@': str(N_J), '@@NB@@': str(N_B), '@@WJ@@': tex_num(WJ), '@@WB@@': tex_num(WB), '@@WE@@': '%.2f' % WE,
       '@@WNB@@': tex_num(WNB),
       '@@PSUM@@': tex_num(PSUM), '@@P1@@': '%.0f' % pc[0]['poles_relative_to_G2_IV']['G1 II(1)'],
       '@@P2@@': '%.0f' % pc[0]['poles_relative_to_G2_IV']['G1 V(2b UV)'],
       '@@FC4@@': tex_num(FC4), '@@TLOK@@': tex_num(max(TL['log restored, F1, F2'])),
       '@@TLF2@@': '%.2f--%.2f' % (min(TL['log restored, F2']), max(TL['log restored, F2'])),
       '@@TLNO@@': '%.2f--%.2f' % (min(TL['F2']), max(TL['F2'])),
       '@@EQ134@@': FOURNUM[134], '@@SLW@@': tex_num(SLW), '@@RVA@@': tex_num(RVA), '@@RVB@@': tex_num(RVB),
       '@@RVC@@': tex_num(RVC),
       '@@ECID@@': tex_num(EC_ID), '@@ECB@@': tex_num(EC_B), '@@ECT0@@': tex_num(EC_T0),
       '@@M2A@@': '%.2f' % max(max(r.values()) for r in M2S), '@@M2B@@': tex_num(max(max(r.values()) for r in M2C)),
       '@@UV@@': tex_num(max(max(v for v in r.values() if isinstance(v, float)) for r in STEPS['Group_I']['Row V 2b: UV term'])),
       '@@RS0@@': '%.2f' % RS0, '@@RS1@@': '%.2f' % RS1, '@@RC@@': tex_num(RC), '@@RB@@': tex_num(RB)}
def main():
    out = DOC
    for k, v in REP.items():
        out = out.replace(k, v)
    FA, FB = four_rho()
    FA = FA.replace(r'\end{align}', r'\label{eq:foura}' + '\n' + r'\end{align}', 1)
    FB = FB.replace(r'\end{align}', r'\label{eq:fourb}' + '\n' + r'\end{align}', 1)
    for k, v in (('@@LOGTABLE@@', log_table()), ('@@ROWS@@', rows_tex('LL')), ('@@FOUR_A@@', FA), ('@@FOUR_B@@', FB),
                 ('@@THREE@@', three_section()), ('@@TWO@@', TWO), ('@@BF@@', BF), ('@@TARGETS@@', TARGETS),
                 ('@@STEPS@@', step_rows()), ('@@FOUR@@', four_rows()), ('@@FOU@@', fourier_line())):
        out = out.replace(k, v)
    assert '@@' not in out, out[out.index('@@') - 100:out.index('@@') + 50]
    open(os.path.join(HERE, 'Checks_JIMWLK_BFKL.tex'), 'w').write(out)
    print('written Checks_JIMWLK_BFKL.tex')


if __name__ == '__main__':
    main()
