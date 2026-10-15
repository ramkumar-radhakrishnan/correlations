"""Writes ThreeRho_cancellation_v3.tex: the three-rho terms of input/three_rho_all_v3.tex, organised by kernel class
(alpha, alpha', beta, beta') and number of Wilson lines (2, 3, 4), with every colour network in its own colour."""
import json
import os
import re
import numpy as np
import parse3 as P
import networks3 as N

HERE = os.path.dirname(os.path.abspath(__file__))
COLS = ['RoyalBlue', 'BrickRed', 'ForestGreen', 'DarkOrchid', 'BurntOrange', 'TealBlue', 'Mahogany', 'Magenta',
        'OliveGreen', 'NavyBlue', 'RedViolet', 'Sepia']
KTEX = {
    'alpha': r"\mathcal{K}(x'-w')\cdot\mathcal{K}(x-w)\;\mathcal{K}(x-z)\cdot\mathcal{K}(y-z)",
    "alpha'": r"\mathcal{K}(x'-w')\cdot\mathcal{K}(x-w)\;\mathcal{K}(x'-z)\cdot\mathcal{K}(y-z)",
    'beta': r"\mathcal{K}(x-w)\cdot\mathcal{K}(x'-w')\;\mathcal{K}(w-z)\cdot\mathcal{K}(y-z)",
    "beta'": r"\mathcal{K}(x-w)\cdot\mathcal{K}(x'-w')\;\mathcal{K}(w'-z)\cdot\mathcal{K}(y-z)",
}
NAME = {'alpha': r'$\alpha$', "alpha'": r"$\alpha'$", 'beta': r'$\beta$', "beta'": r"$\beta'$"}
WORD = {2: 'Two', 3: 'Three', 4: 'Four'}


def leg_tex(l, idx):
    chain, q = l
    if not chain:
        return r'\rho^{%s}(%s)' % (idx, q)
    return r'\big[%s\rho(%s)\big]^{%s}' % (''.join('U(%s)' % p for p in chain), q, idx)


def net_tex(legs):
    return r'f^{abc}\,' + r'\,'.join(leg_tex(l, i) for l, i in zip(legs, 'abc'))


def colour_written(e):
    t = e['text']
    i = t.index('f^{')
    s = t[i:].strip()
    s = re.sub(r'(\\,|\s)+$', '', s)
    return s


def num(c):
    c = complex(c)
    assert abs(c.imag) < 1e-12
    v = c.real
    if abs(v - round(v)) < 1e-12:
        return '%+d' % round(v)
    return '%+g' % v


def build():
    E = P.parse()
    R = N.analyse(E)
    chk = json.load(open(os.path.join(HERE, 'cancel_check.json')))
    out = []
    w = out.append
    w(r'''\documentclass[10pt]{article}
\usepackage[a4paper,margin=0.6in,landscape]{geometry}
\usepackage{amsmath,amssymb}
\usepackage[dvipsnames]{xcolor}
\usepackage{booktabs,longtable,array}
\usepackage[colorlinks=true,linkcolor=blue!50!black]{hyperref}
\allowdisplaybreaks
\setlength{\LTpre}{2pt}\setlength{\LTpost}{6pt}
\renewcommand{\arraystretch}{1.25}
\title{Three-$\rho$ terms $\propto\log(\vee/\Lambda)$: they cancel exactly\\[2pt]\large by kernel class ($\alpha,\alpha',\beta,\beta'$) and number of Wilson lines}
\date{}
\begin{document}
\maketitle
\sloppy
''')
    nterms = len(E)
    nnet = sum(len(v) for v in R.values())
    w(r'''\section*{Result}
\begin{itemize}
\item \textbf{All three-$\rho$ terms cancel exactly.} The list has %d terms. After symmetrisation of the charges they form %d distinct colour networks, and the coefficients of every network add up to zero.
\item \textbf{The cancellation is term against term.} In every network the cancelling terms have literally the same colour structure: the same charges, rotated by the same chain of Wilson lines, up to the order of the indices of $f^{abc}$ (a sign). No identity is needed: neither $U^TU=1$ nor $f^{abc}U^{aa'}U^{bb'}U^{cc'}=f^{a'b'c'}$.
\item \textbf{Each class and each number of Wilson lines cancels separately.} The four kernel classes $\alpha,\alpha',\beta,\beta'$ do not mix. Within a class every network contains terms with the same number of Wilson lines, so the two-, three- and four-Wilson-line terms each cancel among themselves (Secs.~\ref{sec:alpha}--\ref{sec:betap}).
\item \textbf{Numerical check.} With random SU(3) Wilson lines and charges, the residual of every kernel monomial is below $3\times10^{-15}$, while single terms are of order one (Sec.~\ref{sec:num}).
\item \textbf{Same as the corrected list.} Group by group, the list equals the corrected list of \texttt{ThreeRho\_final.pdf} (T1--T3, F3, F4 applied).
\item \textbf{Convention $U^{bc}=U^{cb}$.} The list is written with your convention. Read literally with a general (non-symmetric) adjoint $U$, the sums do not vanish. They do vanish for a general $U$ if every chain is read as $U(p_1)\cdots U(p_n)\rho$, with each $U$ oriented from $f$ towards the charge. So the cancellation does not rely on $U^2=1$; only the index order of the written placements assumes $U=U^T$.
\end{itemize}
''' % (nterms, nnet))
    w(r'''\section*{Conventions}
\begin{itemize}
\item All terms carry the common factor ${\cal N}\equiv\frac{1}{(2\pi)^3}\frac{ig^4}{16\pi^5}\frac{1}{k^+}\log\frac{\vee}{\Lambda}\int_{x,y,z,x',w,w'}e^{-ik(w'-w)}$. The four-$\rho$ block $d^3N_a$ has $-{\cal N}$ (its prefactor is $-\frac{i}{(2\pi)^3}\frac{g^4}{8\pi^3}\frac{\log}{2\pi^2}$), $d^3N_b$ has $+{\cal N}$.
\item The coefficient $c$ of a term (in units of ${\cal N}$) is its sign $\times\,2$ (for ``$2\,$'') $\times\,\frac12$ (for ``$\frac12$'') $\times\,2$ (for an anticommutator). Only the fully symmetrised part of the charges is compared, so $\frac12\rho\{\rho,\rho\}$ counts as $\rho\rho\rho$ and $\{\rho,\rho\}\rho$ as $2\rho\rho\rho$.
\item Within a class the kernel is common to all terms, and the charges sit at $x$, $x'$, $y$ with fixed roles. The terms of a class cancel if and only if their colour factors, as functions of $x,x',y,w,w',z$, add up to zero.
\item A colour network is written $f^{abc}L_1^aL_2^bL_3^c$. Each leg $L$ is a charge rotated by a chain of Wilson lines, e.g.\ $\big[U(w)U(w')\rho(x')\big]^c=U^{cd}(w)U^{de}(w')\rho^e(x')$ (with $U^{bc}=U^{cb}$). A term written in another index order is mapped to its network by reordering $f^{abc}$, which gives the sign $s=\pm1$; its contribution to the network is $c\,s$.
\item Labels: $4\rho_a$, $4\rho_b$ = from four $\rho$'s ($d^3N_a$, $d^3N_b$); I.II.1 = Group I Row II part 1, \dots; III.VI = Group III Row VI. The number in brackets is the position of the term in its block, as sent.
\item In each subsection, every network has its own colour. Its terms are listed together, and the last line of each network is the sum of their contributions.
\end{itemize}
''')
    # summary table
    w(r'\section*{Summary}' + '\n' + r'\begin{center}\begin{tabular}{lcccc}\toprule' + '\n'
      + r'class & Wilson lines & terms & networks & sum of every network\\\midrule' + '\n')
    for cls, nets in R.items():
        for nU in (2, 3, 4):
            ns = [n for n in nets if n['terms'][0][0]['nU'] == nU]
            nt = sum(len(n['terms']) for n in ns)
            ok = all(abs(sum(c for _, c, _ in n['terms'])) < 1e-12 for n in ns)
            w(r'%s & %d & %d & %d & %s\\' % (NAME[cls], nU, nt, len(ns), r'$0$ (all)' if ok else r'\textbf{nonzero}') + '\n')
        w(r'\midrule' + '\n')
    w(r'total & & %d & %d & \\\bottomrule\end{tabular}\end{center}' % (nterms, nnet) + '\n')

    labels = {'alpha': 'sec:alpha', "alpha'": 'sec:alphap', 'beta': 'sec:beta', "beta'": 'sec:betap'}
    for cls, nets in R.items():
        w(r'\section{%s terms}\label{%s}' % (NAME[cls], labels[cls]) + '\n')
        w(r'Kernel: $%s$. %d terms, %d networks.' % (KTEX[cls], sum(len(n['terms']) for n in nets), len(nets)) + '\n')
        for nU in (2, 3, 4):
            ns = [n for n in nets if n['terms'][0][0]['nU'] == nU]
            assert all(e['nU'] == nU for n in ns for e, _, _ in n['terms'])
            w(r'\subsection{%s Wilson lines}' % WORD[nU] + '\n')
            w(r'\begin{longtable}{@{}>{\raggedright}p{0.06\linewidth}p{0.075\linewidth}>{\centering}p{0.04\linewidth}p{0.62\linewidth}>{\centering}p{0.03\linewidth}>{\centering\arraybackslash}p{0.06\linewidth}@{}}' + '\n')
            w(r'\toprule net & source & $c$ & colour factor as written & $s$ & $c\,s$\\\midrule\endhead' + '\n')
            for k, n in enumerate(ns):
                col = COLS[k % len(COLS)]
                tag = r'$\mathsf N_{%d}$' % (k + 1)
                w(r'\multicolumn{6}{@{}l}{\color{%s}%s $=\;%s$}\\' % (col, tag, net_tex(n['legs'])) + '\n')
                parts = []
                for e, cs, legs0 in n['terms']:
                    c = (e['pref'] / 1j) * e['sign'] * e['mult']
                    s = int(round((cs / c).real))
                    src = r'%s\,(%d)' % (P.SHORT[e['block']], e['n'])
                    w(r'\color{%s}%s & \color{%s}%s & \color{%s}$%s$ & \color{%s}$%s$ & \color{%s}$%s$ & \color{%s}$%s$\\' % (
                        col, tag, col, src, col, num(c), col, colour_written(e), col, '+' if s > 0 else '-', col, num(cs)) + '\n')
                    parts.append(num(cs))
                tot = sum(cs for _, cs, _ in n['terms'])
                w(r'\multicolumn{6}{@{}r}{\color{%s}\textbf{sum for %s:} $%s=%s$}\\[2pt]' % (
                    col, tag, ''.join(parts).lstrip('+'), '0' if abs(tot) < 1e-12 else num(tot)) + '\n')
            w(r'\bottomrule\end{longtable}' + '\n')

    # numerical check
    w(r'\section{Numerical check}\label{sec:num}' + '\n')
    w(r'''Each term is reduced to its symmetrised charge part. The three charges are placed on three points $p_1,p_2,p_3$ in all $3!$ ways, the terms are collected by kernel monomial, and the colour factors are evaluated with random SU(3) adjoint Wilson lines and random charges (three configurations). Script: \texttt{check\_cancel.py}.
\begin{center}\begin{tabular}{llccc}\toprule
Wilson lines & set & monomials & non-cancelling (3 configurations) & largest residual / largest term\\\midrule
''')
    for name, d in chk.items():
        lab = r'$U=U^T$ (your convention)' if 'symmetric' in name else r'general $U$, as written'
        for k, v in d.items():
            kk = {'all terms': 'all', 'alpha': r'$\alpha$', "alpha'": r"$\alpha'$", 'beta': r'$\beta$', "beta'": r"$\beta'$"}[k]
            w(r'%s & %s & %d & %s & $%.1e$ / $%.1e$\\' % (lab, kk, v[0]['monomials'], ', '.join(str(x['nonzero']) for x in v),
                                                        max(x['max_residual'] for x in v), v[0]['max_single']) + '\n')
        w(r'\midrule' + '\n')
    w(r'''\bottomrule\end{tabular}\end{center}
With a general $U$ and every chain read as $U(p_1)\cdots U(p_n)\rho$ (oriented from $f$ to the charge), the sum over the terms of each class vanishes to $10^{-14}$ (\texttt{general\_u.py}).

\paragraph{Reproducing.} \texttt{python3 parse3.py} (reads the list; checks that every index occurs twice and every kernel matches its label), \texttt{python3 check\_cancel.py}, \texttt{python3 networks3.py}, \texttt{python3 general\_u.py}, \texttt{python3 cmp\_v2.py} (comparison with the corrected list of \texttt{three\_rho\_v2}), \texttt{python3 make\_note3.py}; then \texttt{pdflatex ThreeRho\_cancellation\_v3.tex} twice.
\end{document}
''')
    open(os.path.join(HERE, 'ThreeRho_cancellation_v3.tex'), 'w').write(''.join(out))


if __name__ == '__main__':
    build()
    print('wrote ThreeRho_cancellation_v3.tex')
