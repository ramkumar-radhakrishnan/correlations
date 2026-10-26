r"""Writes the two documents of this directory from your Group_I/II/III.tex and 4rho.tex (input/):

  LL_terms.tex           every term proportional to log(k+/Lambda) (JIMWLK) and to log(vee/k+) (BFKL), row by row;
                         the three-rho terms and their cancellation; the two-rho terms and their sums.
  Finite_terms_rows.tex  everything else, row by row, with Lambda -> 0 wherever possible.

The split of each term is made by gen.py (exact, checked term by term in check_split.py).  The two-rho sections and
the targets are those of ../final_check/Final_evolution_check.tex, rewritten with log(k+/Lambda) and log(vee/k+)."""
import json
import os
import re
import assemble as AS
import gen as GN
import logsplit as LS
import three_check as T3

HERE = os.path.dirname(os.path.abspath(__file__))
BLK, SPL = GN.all_splits()
LKL, LVK, LVK1 = LS.LKL, LS.LVK, LS.LVK1
FEC = open(os.path.join(HERE, '..', 'final_check', 'Final_evolution_check.tex')).read()
CS = json.load(open(os.path.join(HERE, 'check_split.json')))
TC = json.load(open(os.environ.get('TC_JSON', os.path.join(HERE, 'three_check.json'))))
PC = json.load(open(os.path.join(HERE, 'pole_check.json')))


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

HIDDEN = {'G3 III(1)': 'III(1)', 'G3 III(2)': 'III(2)', 'G3 V': 'V', 'G3 VI': 'VI'}


def hidden_terms():
    """the terms (B) of BFKL_final.pdf, as written there (with y')"""
    S = open(os.path.join(HERE, '..', 'bfkl_final', 'static.tex')).read()
    i = S.index(r'\paragraph{The new terms}')
    blk = S[S.index(r'\begin{align}', i) + len(r'\begin{align}'):S.index(r'\end{align}', i)]
    parts = re.split(r'\\text\{Row ([IV()12]+):\}\\quad&', blk)[1:]
    out = {}
    for k in range(0, len(parts), 2):
        t = parts[k + 1].replace(r'\Lv', LVK).replace(r'\nc', '')
        t = re.sub(r'[,.]\s*(\\\\)?\s*$', '', t.strip())
        out[parts[k]] = t
    return out


HB = hidden_terms()
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
    if name in HIDDEN:
        items.append((HB[HIDDEN[name]], 'blue'))
    if name == 'G2 IV':
        items.append((H4 % LVK, 'blue'))
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


def f_items(name):
    items = []
    for n, s in enumerate(SPL[name]):
        if s.get('power') or s.get('red'):
            items.append((BLK[name].terms[n], 'red'))
        elif (name, n) == ('G1 III N3b', 3):
            items.append((lam0(s['F']), 'red'))
        elif s.get('F') is not None:
            items.append((lam0(c5(name, s['F'])), None))
    if name in HIDDEN:
        items.append((neg(HB[HIDDEN[name]]), 'blue'))
    if name == 'G2 IV':
        items.append(('-' + H4 % LVK, 'blue'))
    return items


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
    'G1 II(1)': r'''Your bracket $\big[\frac{11}{3}-2\log\frac{k^+}{\Lambda}+\log(1-\frac{k^+}{\vee})\big]$ multiplies the pole bracket; only its $-2\log\frac{k^+}{\Lambda}$ is a large logarithm. The term $\log(1-\frac{k^+}{\vee})\log\frac{\vee}{k^+}$ vanishes for $\vee\to\infty$ and is not a BFKL term; it stays with the finite terms.''',
    'G1 II(2a)': r'''Line 1 has $\big[\log\frac{\vee}{k^+}+\log\frac{k^+}{\Lambda}\big]$ and enters both lists. In lines 2 and 3 the $z$ integral is already done (hence $\frac{1}{32\pi^4}$); the measure should read $\int_{x,y,x'}$.''',
    'G1 II(2b)': r'''No large logarithm ($R^i$ depends on $\vee$ only through $\frac{k^+}{\vee}$).''',
    'G1 III N1': r'''Lines 1 and 4 (your green $\frac{1}{\Lambda}$ brackets) are not logarithms. Together with line 1 of $d^3N_2$ and line 1 of $d^3N_{3a}$ they cancel exactly (other file).''',
    'G1 III N2': r'''Line 1 is a $\frac{1}{\Lambda}$ term (see $d^3N_1$). Your two red lines cancel each other and are left out.''',
    'G1 III N3a': r'''\old{Red}: the first logarithmic line, with the edge kernel $\frac{(x'-w')\cdot(y-z)}{(x'-w')^2(y-z)^2(z-w)^2}$, is not integrable at $z\to w$. \new{Blue}: it is replaced by $\mathcal D$ (fix F3, \texttt{BFKL\_final.pdf} Sec.~2.1): $U^{bc}(z)\to U^{bc}(w)$ and the kernel $\mathcal K(x'-w')\cdot\mathcal K(y-w)\,\mathcal K(x-z)\cdot\mathcal K(w-z)$. Without this change the three-$\rho$ terms do not cancel (Sec.~\ref{sec:three}).''',
    'G1 III N3b': r'''No large logarithm. ($\mathcal A$ and the logarithms of $d^3N_{3a}$, $d^3N_{3b}$ have finite limits for $\Lambda\to0$; other file.)''',
    'G1 V(2b UV)': r'''\textbf{E1.} Your bracket $\big[\log\frac{k^+}{\Lambda}+2\log(1-\frac{\Lambda}{k^+})-\frac{11}{6}\big]$ must be $\big[\new{2}\log\frac{k^+}{\Lambda}+2\log(1-\frac{\Lambda}{k^+})-\frac{11}{6}\big]$. The error enters in the equation after ``Setting $\Lambda\to0$'': there your blue bracket $\big[\log\frac{\vee}{k^+}+\log\frac{k^+}{\Lambda}\big]$ replaces $2\log\frac{\vee}{\Lambda}$ of the equation before, so it needs a factor \new{2}; then $2\log\frac{\vee}{\Lambda}-2\log\frac{\vee}{k^+}=2\log\frac{k^+}{\Lambda}$. (The counterterm, line 6 of $d^3N_{2b\,\rm finite}$, has $2\log\frac{k^+-\Lambda}{\Lambda}$, as it should.) With the 2 the poles of the $\log\frac{k^+}{\Lambda}$ terms cancel: Group I Row II : Row V 2b UV : Group II Row IV $=@@P1@@:@@P2@@:1$, sum $@@PSUM@@$ (relative); with your coefficient the sum is $@@PW@@$ (\texttt{pole\_check.py}).''',
    'G1 V(2b fin)': r'''Lines 4--6 carry $\log\frac{k^+-\Lambda}{\Lambda}=\log\frac{k^+}{\Lambda}+\log(1-\frac{\Lambda}{k^+})$.''',
    'G2 I(2a)': r'''No large logarithm: $\vee$ cancels between numerator and denominator of the first logarithm, and the second has the finite limit $\log\frac{(x-w)^2}{(x'-w')^2}$.''',
    'G2 IV': r'''\textbf{Hidden logarithm.} In line 2, $\big[\log\frac{k^+}{\Lambda}+\log(\frac{\vee}{k^+}-1)\big]$, the $\log\frac{\vee}{k^+}$ is cancelled by the one in the $\frac{\xi}{(1-\xi)_+}$ integral of line 1 (lower limit $\frac{k^+}{\vee}$). It is not a BFKL term. The BFKL term of this row is \new{$\log\frac{\vee}{k^+}H_4$} (blue): at $\xi\to\frac{k^+}{\vee}$ the phase $e^{-ik\cdot(y'-y)\frac1\xi}$ probes $|y'-y|\to0$, and the large-momentum tail of $\log(y'-y)^2$ gives $y,y'\to z$ (\texttt{Finite\_terms\_final.pdf}, C1).''',
    'G3 III(1)': r'''\textbf{Hidden logarithm.} The $\log\frac{\vee}{k^+}$ in $\big[\log\frac{k^+}{\Lambda}+\log(\frac{\vee}{k^+}-1)\big]$ is cancelled by the one in the $\left[\frac{1}{p^+}\right]_+$ integral of line 1 (its subtraction $-\frac{f(0)}{p^+}$ up to $\vee-k^+$). It is not a BFKL term. The BFKL term of the row is the \new{blue} term (B) of \texttt{BFKL\_final.pdf}: the terms of the $p^+$ bracket that tend to a constant at large $p^+$ give $\log\frac{\vee}{k^+}$ after $\int_{y'}$, with $y'\to z$.''',
    'G3 III(2)': r'''Hidden logarithm as in the first part; the \new{blue} term is (B). \textbf{C5}: the measure is $\int_{w,x',y',z,x}$ (the row has $\rho(x')\rho(x)$; there is no $y$ and no $w'$).''',
    'G3 V': r'''Hidden logarithm as in Row III; the \new{blue} term is (B).''',
    'G3 VI': r'''Hidden logarithm as in Row III; the \new{blue} term is (B).''',
}
pc = list(PC.values())
NOTE_LL['G1 V(2b UV)'] = (NOTE_LL['G1 V(2b UV)'].replace('@@P1@@', '%.0f' % pc[0]['poles_relative_to_G2_IV']['G1 II(1)'])
                          .replace('@@P2@@', '%.0f' % pc[0]['poles_relative_to_G2_IV']['G1 V(2b UV)'])
                          .replace('@@PSUM@@', tex_num(max(p['sum_corrected_relative'] for p in pc)))
                          .replace('@@PW@@', '%.1f' % max(p['sum_as_written_relative'] for p in pc)))

NOTE_F = {
    'G1 II(1)': r'''The pole bracket times $\frac{11}{3}$ is the $\beta_0$ term; it cancels against the $-\frac{11}{6}$ of Row V 2b UV (\texttt{Finite\_terms\_final.pdf}, Sec.~4). $\log(1-\frac{k^+}{\vee})\to0$ for $\vee\to\infty$.''',
    'G1 II(2b)': r'''Unchanged: no large logarithm. $R^i$ as in your file.''',
    'G1 III N1': r'''\old{Red}: the $\frac{1}{\Lambda}$ lines (your green brackets). Lines 1 and 4 here, line 1 of $d^3N_2$ and line 1 of $d^3N_{3a}$ add up to zero exactly, pointwise (\texttt{check\_split.py}: $@@POW@@$).''',
    'G1 III N2': r'''\old{Red}: line 1 is a $\frac{1}{\Lambda}$ line (see $d^3N_1$); your two red lines cancel each other.''',
    'G1 III N3a': r'''\old{Red}: line 1 is a $\frac{1}{\Lambda}$ line (see $d^3N_1$). The edge-kernel line has no finite part: its $\log\frac{k^+-\Lambda}{\Lambda}=\log\frac{k^+}{\Lambda}+\log(1-\frac{\Lambda}{k^+})$ and $\log(1-\frac{\Lambda}{k^+})\to0$. $\mathcal A_0$ is $\mathcal A$ at $\Lambda=0$, Eq.~\eqref{eq:A0}.''',
    'G1 III N3b': r'''\textbf{C4.} The last line (\old{red}) is, for $\Lambda\to0$, $\log\frac{(z-w)^2}{(x-w)^2}$ times the edge kernel $\frac{(x'-w')\cdot(y-z)}{(x'-w')^2(y-z)^2(z-w)^2}$: not integrable at $z\to w$. The edge kernel is replaced by $\mathcal D$ (LL file), and its finite part belongs to the region of $\mathcal D$, which is not yet computed (Q3, open).''',
    'G1 V(2b UV)': r'''Only $-\frac{11}{6}$ is left (E1, $2\log(1-\frac{\Lambda}{k^+})\to0$). It is the $\beta_0$ term that cancels the $\frac{11}{3}$ of Row II part 1.''',
    'G1 V(2b fin)': r'''$\Lambda\to0$ in the limits of the plus-prescribed integrals (regular at the end points). In line 6, $2\log\frac{k^+-\Lambda}{\Lambda}$ went to the LL file; $-\frac{11}{6}$ is the counterterm constant.''',
    'G2 I(2f)': r'''$\frac{1}{k^++\Lambda}\to\frac{1}{k^+}$.''',
    'G2 IV': r'''\textbf{C1.} The \new{blue} line, $-\log\frac{\vee}{k^+}H_4$, is the BFKL term contained in the $\xi$ integral of line 1; it is now in the BFKL list. Line 2 keeps $\log(\frac{\vee}{k^+}-1)$, which cancels the $\log\frac{\vee}{k^+}$ of the $\frac{\xi}{(1-\xi)_+}$ integral. Lines 1, 2 and the blue line are the DGLAP term (\texttt{Finite\_terms\_final.pdf}, Sec.~5). In line 3, $\frac{k^+}{k^++\Lambda}\to1$.''',
    'G3 III(1)': r'''\textbf{C3.} The \new{blue} line is $-\log\frac{\vee}{k^+}\times$(B): the BFKL term of the $\left[\frac{1}{p^+}\right]_+$ integral (region $y'\to z$), now in the BFKL list. Line 2 keeps $\log(\frac{\vee}{k^+}-1)$, which cancels the $-\log\frac{\vee}{k^+}$ of the subtraction in $\left[\frac{1}{p^+}\right]_+$. Lines 1, 2 and the blue line together are finite for $\vee\to\infty$.''',
    'G3 III(2)': r'''C3 as in the first part. \textbf{C5}: measure.''',
    'G3 V': r'''C3 as in Row III.''',
    'G3 VI': r'''C3 as in Row III.''',
}
NOTE_F['G1 III N1'] = NOTE_F['G1 III N1'].replace('@@POW@@', tex_num(CS['power_terms_sum_rel']))


def rows_tex(kind):
    out = []
    for gname, gl, rows in GROUPS:
        out.append(r'\section{%s}\label{sec:%s%s}' % (gname, kind, gl))
        for rname, names, special in rows:
            out.append(r'\subsection{%s}' % rname)
            if special == 'four':
                out.append(r'Only four $\rho$\textquotesingle s: Sec.~\ref{sec:four} (with the same coefficient for both logarithms).'
                           if kind == 'LL' else r'Only four $\rho$\textquotesingle s (\texttt{LL\_terms.pdf}, Sec.~6); no finite terms in your files.')
                continue
            for name in names:
                sub = SUB[name]
                out.append(r'\paragraph{%s%s.}' % (rname, (', ' + sub) if sub else ''))
                note = (NOTE_LL if kind == 'LL' else NOTE_F).get(name)
                if note:
                    out.append(note)
                if kind == 'LL':
                    J, B = j_items(name), b_items(name)
                    out.append(r'\noindent Terms $\propto\log\frac{k^+}{\Lambda}$:' + ('' if J else ' none.'))
                    out.append(align(name, J, 'J'))
                    out.append(r'\noindent Terms $\propto\log\frac{\vee}{k^+}$:' + ('' if B else ' none.'))
                    out.append(align(name, B, 'B'))
                else:
                    F = f_items(name)
                    out.append(align(name, F, 'F') if F else 'Nothing is left.')
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


TARGETS = sec(FEC, r'\section{The targets}', r'\section{All two-$\rho$')
TARGETS = once(TARGETS, r'\subsection{JIMWLK$\otimes$LO: your formula, corrected (M2)}',
               r'\subsection{JIMWLK$\otimes$LO: your formula (App.~F), corrected (M2)}')
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

# ------------------------------------------------------------------------------------------------ LL document
N_J = len(TC['J']['terms'])
N_B = len(TC['B']['terms'])
WJ = TC['J']['variants']['as in the note']['relative']
WB = TC['B']['variants']['as in the note']['relative']
WE = TC['J']['variants']['edge kernel instead of D']['relative']

LL = PREAMBLE + r'''
\title{The terms proportional to $\log\frac{k^+}{\Lambda}$ (JIMWLK) and $\log\frac{\vee}{k^+}$ (BFKL),\\ row by row from \texttt{Group\_I/II/III.tex} and \texttt{4rho.tex}}
\date{}
\begin{document}
\maketitle
\noindent Notation as in your files: $\mathcal K(X)=X/X^2$, $K(a,b,z)=\mathcal K(a-z)\cdot\mathcal K(b-z)$, $U^{bc}=U^{cb}$,
\begin{equation*}
\mathcal N=\frac{1}{(2\pi)^3}\frac{g^4}{16\pi^5}\frac{1}{k^+}.
\end{equation*}
\new{Blue}: new or corrected; \old{red}: removed. The finite terms (everything that is not $\propto\log\frac{k^+}{\Lambda}$ or $\log\frac{\vee}{k^+}$) are in the separate file \texttt{Finite\_terms\_rows.pdf}.

\section{Answer}\label{sec:answer}
Every row is split exactly, term by term,
\begin{equation}
\frac{d^3N_{\rm row}}{d^3k}=\log\frac{k^+}{\Lambda}\;\mathfrak T_{\rm row}+\log\frac{\vee}{k^+}\;\mathfrak P_{\rm row}+R_{\rm row},
\label{eq:split}
\end{equation}
with $\log\frac{\vee}{\Lambda}=\log\frac{k^+}{\Lambda}+\log\frac{\vee}{k^+}$ and the splits of Sec.~\ref{sec:logs}; $\vee>p^+$, $k^+>\Lambda$; $\Lambda\to0$ in $R_{\rm row}$ wherever possible, $\vee$ exact. Numerically, for every one of the @@NCHK@@ terms that can be evaluated pointwise, your term $=$ its three parts to $@@ID@@$, and the remainder has its $\Lambda\to0$ limit (\texttt{check\_split.py}).
\begin{itemize}
\item \textbf{Three $\rho$\textquotesingle s $\propto\log\frac{k^+}{\Lambda}$} (@@NJ@@ terms): \textbf{cancel}, relative residual $@@WJ@@$ (Sec.~\ref{sec:three}).
\item \textbf{Three $\rho$\textquotesingle s $\propto\log\frac{\vee}{k^+}$} (@@NB@@ terms): \textbf{cancel}, $@@WB@@$.
\item \textbf{Four $\rho$\textquotesingle s} (your \texttt{4rho.tex}): they have $\log\frac{\vee}{\Lambda}$, i.e.\ the same coefficient for both logarithms. They cancel after symmetrization (Sec.~\ref{sec:four}); their three- and two-$\rho$ remainders enter both lists.
\item \textbf{Two $\rho$\textquotesingle s $\propto\log\frac{k^+}{\Lambda}$} $=$ \textbf{JIMWLK$\otimes$LO}, exactly: 13 kernel structures, 154 independent colour structures, 0 differences (Sec.~\ref{sec:two}).
\item \textbf{Two $\rho$\textquotesingle s $\propto\log\frac{\vee}{k^+}$} $=$ \textbf{BFKL$\otimes$LO}, exactly: 9 kernel structures, 109 colour structures, 0 differences (Sec.~\ref{sec:bfkl}).
\end{itemize}
Therefore
\begin{equation}
\frac{d^3N_{\rm LL}}{d^3k}=\mathcal N\Big[\log\frac{k^+}{\Lambda}\;{\rm JIMWLK}\otimes{\rm LO}+\log\frac{\vee}{k^+}\;{\rm BFKL}\otimes{\rm LO}\Big],
\label{eq:LL}
\end{equation}
with the targets of Sec.~\ref{sec:targets}. The target evolves (JIMWLK) between $y_T=\log\Lambda$ and the gluon, the projectile (BFKL) between the gluon and $y_P=\log\vee$.

\paragraph{What has to change in your files} (shown in red/blue in the rows):
\begin{itemize}
\item[E1] Group I Row V 2b UV: the bracket must have $\new{2}\log\frac{k^+}{\Lambda}$ (from $2\log\frac{\vee}{\Lambda}-2\log\frac{\vee}{k^+}$; the factor 2 is lost after ``Setting $\Lambda\to0$''). With your coefficient 1 the poles $\times\log\frac{k^+}{\Lambda}$ do not cancel (Sec.~\ref{sec:LLgI}, Row V).
\item[F3] Group I Row III $d^3N_{3a}$: the logarithmic line with the edge kernel is replaced by $\mathcal D$. Without it the three-$\rho$ terms $\propto\log\frac{k^+}{\Lambda}$ do not cancel ($@@WE@@$).
\item[H] Hidden logarithms. In Group III Rows III (both parts), V, VI and Group II Row IV the $\log\frac{\vee}{k^+}$ inside $\big[\log\frac{k^+}{\Lambda}+\log(\frac{\vee}{k^+}-1)\big]$ is \emph{not} a BFKL term: it cancels against the one inside $\left[\frac{1}{p^+}\right]_+$ (resp.\ $\frac{\xi}{(1-\xi)_+}$) of the same row. The BFKL terms of these rows are the terms (B) of \texttt{BFKL\_final.pdf} and $H_4$ (blue). Their $\log\frac{k^+}{\Lambda}$ parts are as you have them.
\item[C5] Group III Row III part 2: the measure is $\int_{w,x',y',z,x}$.
\item[F1, F2] The two-$\rho$ lists of \texttt{4rho.tex} still need the two fixes of \texttt{Final\_evolution\_check.pdf} (Sec.~\ref{sec:two}).
\item[M2] Your JIMWLK formula (App.~F) needs the correction of Sec.~\ref{sec:targets}; your BFKL formula is correct.
\end{itemize}
Everything else is as you wrote it.

\section{How the logarithms are split}\label{sec:logs}
A term carries at most one logarithmic factor. Your brackets:
\begin{align}
\left[\log{\left(\frac{k^{+}}{\Lambda}\right)}+\log{\left(\frac{\vee}{k^{+}}-1\right)}\right]&=\log\frac{k^+}{\Lambda}+\log\frac{\vee}{k^+}+\log\Big(1-\frac{k^+}{\vee}\Big),\notag\\
\left[\log{\left[\frac{\vee}{k^{+}}\right]}+\log{\left[\frac{k^{+}}{\Lambda}\right]}\right]&=\log\frac{k^+}{\Lambda}+\log\frac{\vee}{k^+},\notag
\end{align}
except in the rows with a hidden logarithm (H above), where the first bracket contributes only $\log\frac{k^+}{\Lambda}$ and $\log(\frac{\vee}{k^+}-1)$ stays in the finite terms. All other logarithms ($A=(x-z)^2$ etc.; the remainder is given exactly and, where it depends on $\Lambda$, its limit):
@@LOGTABLE@@
Logarithms of $\epsilon$, $\mu$, $\log(1-\frac{k^+}{\vee})$ and $\log(1-\frac{\Lambda}{k^+})$ are finite and are not split. Three special cases:
\begin{itemize}
\item Group I Row II part 1: $\big[\frac{11}{3}-2\log\frac{k^+}{\Lambda}+\log(1-\frac{k^+}{\vee})\big]$ times the pole bracket; only $-2\log\frac{k^+}{\Lambda}$ is large.
\item Group I Row V 2b UV: $\big[\new{2}\log\frac{k^+}{\Lambda}+2\log(1-\frac{\Lambda}{k^+})-\frac{11}{6}\big]$ (E1).
\item Group I Row V 2b finite, line 6: $\big[-\frac{11}{6}+2\log\frac{k^+-\Lambda}{\Lambda}\big]$.
\end{itemize}
Every line of every row below is in your form; only the logarithm is replaced by its large part.

@@ROWS@@

\section{Four \texorpdfstring{$\rho$}{rho}'s}\label{sec:four}
Group I Row I, Group II Rows II and III and Group III Row IV have only four $\rho$\textquotesingle s. They are your \texttt{4rho.tex}, Eq.~(4rho), with
\begin{equation*}
\log\frac{\vee}{\Lambda}=\new{\log\frac{k^+}{\Lambda}+\log\frac{\vee}{k^+}}:
\end{equation*}
every four-$\rho$ term has the same coefficient for both logarithms.
\begin{itemize}
\item The four-$\rho$ terms cancel after symmetrization: your demonstration in \texttt{4rho.tex} (two and four Wilson lines; three Wilson lines by $x\leftrightarrow y$); checked numerically before ($4.7\times10^{-15}$, \texttt{Final\_evolution\_check.pdf}).
\item Their three-$\rho$ remainders, your final forms, enter both three-$\rho$ lists of Sec.~\ref{sec:three}:
\end{itemize}
@@FOUR_A@@
@@FOUR_B@@
\begin{itemize}
\item Their two-$\rho$ remainders (your ``Two $\rho$ contribution'', with F1, F2) enter both two-$\rho$ sums, Secs.~\ref{sec:two} and \ref{sec:bfkl}.
\end{itemize}

@@THREE@@

@@TWO@@

@@BF@@

@@TARGETS@@

\section{Checks}\label{sec:checks}
Scripts in \texttt{llsplit/}:
\begin{center}\small
\begin{tabular}{@{}p{0.62\textwidth}p{0.32\textwidth}@{}}\toprule
check & result\\\midrule
every evaluable term $=$ its $\log\frac{k^+}{\Lambda}$, $\log\frac{\vee}{k^+}$ and finite parts, exactly in $\Lambda$, $\vee$ (\texttt{check\_split.py}; @@NCHK@@ terms, 3 configurations, 2 sets of $k^+,\vee,\Lambda$) & $@@ID@@$\\
finite parts: exact remainder $\to$ its $\Lambda\to0$ form & $@@L0@@$ at $\Lambda=10^{-9}$\\
$\frac{1}{\Lambda}$ lines of Group I Row III add up to zero; your two red lines & $@@POW@@$; $0$\\
$\log\frac{\vee}{k^+}-\log\frac{k^+}{\Lambda}$ coefficient of each three-$\rho$ row vs the exact row coefficient of the earlier check (\texttt{cross\_check.py}) & $\le@@CC@@$\\
poles $\times\log\frac{k^+}{\Lambda}$ cancel, with E1 / as written (\texttt{pole\_check.py}) & $@@PSUM@@$ / $@@PW@@$\\
three $\rho$\textquotesingle s $\propto\log\frac{k^+}{\Lambda}$ / $\propto\log\frac{\vee}{k^+}$ cancel (\texttt{three\_check.py}) & $@@WJ@@$ / $@@WB@@$\\
the same with the edge kernel instead of $\mathcal D$ & $@@WE@@$\\
two $\rho$\textquotesingle s $\propto\log\frac{k^+}{\Lambda}$ vs JIMWLK$\otimes$LO; $\propto\log\frac{\vee}{k^+}$ vs BFKL$\otimes$LO, exact by colour structure (\texttt{../final\_check/render.py}) & 0 of 154; 0 of 109 differ\\
four $\rho$\textquotesingle s cancel (\texttt{Final\_evolution\_check.pdf}) & $4.7\times10^{-15}$\\
\bottomrule
\end{tabular}
\end{center}
\end{document}
'''
CC = json.load(open(os.path.join(HERE, 'cross_check.json')))
REP = {'@@NCHK@@': str(CS['n_checked']), '@@ID@@': tex_num(CS['worst_identity']), '@@L0@@': tex_num(CS['worst_Lambda_to_0']),
       '@@POW@@': tex_num(CS['power_terms_sum_rel']), '@@CC@@': tex_num(max(CC.values())),
       '@@NJ@@': str(N_J), '@@NB@@': str(N_B), '@@WJ@@': tex_num(WJ), '@@WB@@': tex_num(WB), '@@WE@@': tex_num(WE),
       '@@PSUM@@': tex_num(max(p['sum_corrected_relative'] for p in pc)),
       '@@PW@@': '%.1f' % max(p['sum_as_written_relative'] for p in pc)}
for k, v in REP.items():
    LL = LL.replace(k, v)
FA, FB = four_rho()
FA = FA.replace(r'\end{align}', r'\label{eq:foura}' + '\n' + r'\end{align}', 1)
FB = FB.replace(r'\end{align}', r'\label{eq:fourb}' + '\n' + r'\end{align}', 1)
for k, v in (('@@LOGTABLE@@', log_table()), ('@@ROWS@@', rows_tex('LL')), ('@@FOUR_A@@', FA), ('@@FOUR_B@@', FB),
             ('@@THREE@@', three_section()), ('@@TWO@@', TWO), ('@@BF@@', BF), ('@@TARGETS@@', TARGETS)):
    LL = LL.replace(k, v)
open(os.path.join(HERE, 'LL_terms.tex'), 'w').write(LL)

# ------------------------------------------------------------------------------------------------ finite document
FIN = PREAMBLE + r'''
\title{The finite terms of \texttt{Group\_I/II/III.tex}:\\ everything that is not $\propto\log\frac{k^+}{\Lambda}$ or $\log\frac{\vee}{k^+}$, with $\Lambda\to0$}
\date{}
\begin{document}
\maketitle
\noindent Companion of \texttt{LL\_terms.pdf}. Each row is $\log\frac{k^+}{\Lambda}\,\mathfrak T_{\rm row}+\log\frac{\vee}{k^+}\,\mathfrak P_{\rm row}+R_{\rm row}$ exactly; here is $R_{\rm row}$, row by row, in your form. \new{Blue}: new; \old{red}: removed or cancelling.

\section{Rules and summary}\label{sec:rules}
\begin{itemize}
\item Each logarithm is replaced by its remainder (\texttt{LL\_terms.pdf}, Sec.~2): $\big[\log\frac{k^+}{\Lambda}+\log(\frac{\vee}{k^+}-1)\big]\to\log(1-\frac{k^+}{\vee})$, $\log\frac{\vee}{\Lambda+k^+}\to0$, $\log\frac{(\vee-k^+)A+k^+B}{\Lambda A+k^+B}\to\log\frac{(1-\frac{k^+}{\vee})A+\frac{k^+}{\vee}B}{B}$, etc.
\item $\Lambda\to0$ wherever possible (the exact remainders differ by terms that vanish for $\Lambda\to0$: $@@L0@@$ at $\Lambda=10^{-9}$); $\vee$ is kept exact.
\begin{itemize}
\item In the limits of the plus-prescribed integrals of Group I Row V 2b and Group II Row IV: $\int_\Lambda^{k^+-\Lambda}\to\int_0^{k^+}$, $\int_{k^+-\Lambda}^{\Lambda}\to\int_{k^+}^0$, $\frac{k^+}{k^++\Lambda}\to1$.
\item $\mathcal A\to\mathcal A_0$,
\begin{align}
@@A0@@
\label{eq:A0}
\end{align}
\end{itemize}
\item \textbf{Hidden logarithms} (Group II Row IV, Group III Rows III, V, VI): $\big[\log\frac{k^+}{\Lambda}+\log(\frac{\vee}{k^+}-1)\big]\to\log(\frac{\vee}{k^+}-1)$, which cancels the $\log\frac{\vee}{k^+}$ of the $\left[\frac{1}{p^+}\right]_+$ (resp.\ $\frac{\xi}{(1-\xi)_+}$) integral; $\left[\frac{1}{p^+}\right]_+$ is kept unevaluated. The BFKL term contained in that integral is subtracted (\new{blue}, C1, C3).
\item \textbf{Cancelling lines} (\old{red}): the $\frac{1}{\Lambda}$ lines of Group I Row III (sum zero, $@@POW@@$), your two red lines of $d^3N_2$, and the edge-kernel line of $d^3N_{3b}$ (C4, Q3).
\item \textbf{C5}: the measure of Group III Row III part 2.
\item \textbf{Running coupling and DGLAP} are in Group I Row II part 1 ($\frac{11}{3}$), Group I Row V 2b UV ($-\frac{11}{6}$) and Group II Row IV. Their analysis is unchanged (\texttt{Finite\_terms\_final.pdf}, Secs.~4, 5): the $\beta_0$ terms cancel, and Group II Row IV with the blue line is the DGLAP term.
\end{itemize}

@@ROWS@@

\section{Checks}\label{sec:checks}
\texttt{check\_split.py}: for each of the @@NCHK@@ terms without a remaining $p^+$, $\xi$ or momentum integral, your term $=$ (its $\log\frac{k^+}{\Lambda}$ part) $+$ (its $\log\frac{\vee}{k^+}$ part) $+$ (the exact remainder) to $@@ID@@$, at 3 configurations and $(k^+,\vee,\Lambda)=(1,37,0.013),(0.6,900,3\times10^{-4})$; the remainder written here is its $\Lambda\to0$ limit ($@@L0@@$ at $\Lambda=10^{-9}$). The terms with a remaining integral (Group I Row II 2b, Row V; Group II Row IV lines 1, 3; Group III line 1 of Rows III, V, VI) have no $\Lambda$ or large logarithm outside the limits treated above, and are copied.
\end{document}
'''
for k, v in REP.items():
    FIN = FIN.replace(k, v)
FIN = FIN.replace('@@A0@@', A0).replace('@@ROWS@@', rows_tex('F'))
open(os.path.join(HERE, 'Finite_terms_rows.tex'), 'w').write(FIN)
print('written LL_terms.tex, Finite_terms_rows.tex')
