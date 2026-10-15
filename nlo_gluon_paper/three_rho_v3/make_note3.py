"""Writes ThreeRho_cancellation_v3.tex: the three-rho terms of input/three_rho_all_v3.tex, organised by kernel class
(alpha, alpha', beta, beta') and number of Wilson lines (2, 3, 4).  For every class and number of Wilson lines there is
one continuous equation: the common prefactor, f^{abc} and the kernel in front, then every term as written (one per
line), coloured by the set of terms it cancels with, and finally the sum of every coloured set, which is zero."""
import json
import os
import re
import parse3 as P
import networks3 as N

HERE = os.path.dirname(os.path.abspath(__file__))
COLS = ['RoyalBlue', 'BrickRed', 'ForestGreen', 'DarkOrchid', 'BurntOrange', 'TealBlue', 'Mahogany', 'Magenta',
        'OliveGreen', 'NavyBlue', 'RedViolet', 'Sepia']
KTEX = {
    'alpha': r"\mathcal{K}(x'-w')\cdot\mathcal{K}(x-w)\,\mathcal{K}(x-z)\cdot\mathcal{K}(y-z)",
    "alpha'": r"\mathcal{K}(x'-w')\cdot\mathcal{K}(x-w)\,\mathcal{K}(x'-z)\cdot\mathcal{K}(y-z)",
    'beta': r"\mathcal{K}(x-w)\cdot\mathcal{K}(x'-w')\,\mathcal{K}(w-z)\cdot\mathcal{K}(y-z)",
    "beta'": r"\mathcal{K}(x-w)\cdot\mathcal{K}(x'-w')\,\mathcal{K}(w'-z)\cdot\mathcal{K}(y-z)",
}
NAME = {'alpha': r'\alpha', "alpha'": r"\alpha'", 'beta': r'\beta', "beta'": r"\beta'"}
WORD = {2: 'Two', 3: 'Three', 4: 'Four'}
FABC = re.compile(r'f\^\{\s*a\s*b\s*c\s*\}')
PREF = (r"\frac{i\,f^{abc}}{(2\pi)^{3}}\,\frac{g^{4}}{16\pi^{5}}\,\frac{1}{k^{+}}\,\log\Big(\frac{\vee}{\Lambda}\Big)"
        r"\int_{x,y,z,x',w,w'}e^{-ik(w'-w)}\,")


def written_parts(e):
    """(sign, magnitude tex, charge-and-Wilson-line tex) of a term as written; the sign includes the overall minus
    of the four-rho block d^3N_a, which is absorbed so that every term has the common prefactor i f^{abc} ..."""
    t = e['text']
    t2 = re.sub(r'^(\\,)+', '', t.lstrip('+-').strip()).strip()
    mag = ''
    if re.match(r'^2\s*\\,', t2):
        mag = r'2\,'
    elif t2.startswith(r'\frac{1}{2}'):
        mag = r'\tfrac{1}{2}\,'
    m = FABC.search(t)
    body = t[m.end():].strip()
    body = re.sub(r'^(\\,|\s)+', '', body)
    body = re.sub(r'(\\,|\s)+$', '', body)
    sign = e['sign'] * (1 if e['pref'] == 1j else -1)
    return sign, mag, body


def leg_tex(l, idx):
    chain, q = l
    if not chain:
        return r'\rho^{%s}(%s)' % (idx, q)
    return r'\big[%s\rho(%s)\big]^{%s}' % (''.join('U(%s)' % p for p in chain), q, idx)


def net_tex(legs):
    return r'\,'.join(leg_tex(l, i) for l, i in zip(legs, 'abc'))


def num(v):
    v = complex(v).real
    return '%+d' % round(v) if abs(v - round(v)) < 1e-12 else '%+g' % v


def src(e):
    return r'\color{gray}\scriptsize %s\,(%d)' % (P.SHORT[e['block']], e['n'])


def build():
    E = P.parse()
    R = N.analyse(E)
    chk = json.load(open(os.path.join(HERE, 'cancel_check.json')))
    out = []
    w = out.append
    w(r'''\documentclass[10pt]{article}
\usepackage[a4paper,margin=0.7in]{geometry}
\usepackage{amsmath,amssymb}
\usepackage[dvipsnames]{xcolor}
\usepackage{booktabs}
\usepackage[colorlinks=true,linkcolor=blue!50!black]{hyperref}
\allowdisplaybreaks
\setlength{\jot}{3pt}
\title{Cancellation of all three-$\rho$ terms $\propto\log(\vee/\Lambda)$}
\date{}
\begin{document}
\maketitle
\sloppy
''')
    nnet = sum(len(v) for v in R.values())
    w(r'''\section*{Result}
\begin{itemize}
\item \textbf{All %d terms cancel exactly.} After symmetrisation of the charges they form %d sets of identical colour structures, and every set adds up to zero.
\item \textbf{The cancellation is term against term.} The terms of a set have the same charges, rotated by the same chain of Wilson lines; they differ only in the order of the indices of $f^{abc}$, i.e.\ in a sign. No identity is used.
\item \textbf{Each class and each number of Wilson lines cancels separately.} The four classes $\alpha,\alpha',\beta,\beta'$ do not mix, and each set contains terms with a single number of Wilson lines. So the two-, three- and four-Wilson-line terms cancel among themselves.
\item \textbf{Numerical check.} With random SU(3) Wilson lines and charges, every kernel monomial has a residual below $3\times10^{-15}$ (Sec.~\ref{sec:num}).
\item \textbf{Same as the corrected list.} The list equals the corrected list of \texttt{ThreeRho\_final.pdf}. It uses your convention $U^{bc}=U^{cb}$; for a general $U$ the cancellation holds when every chain is read as $U(p_1)\cdots U(p_n)\rho$ (oriented from $f$ to the charge).
\end{itemize}

\paragraph{How to read the equations.}
\begin{itemize}
\item Every term of the list carries $f^{abc}$ and the common prefactor, which are written in front. The kernel is common to all terms of a class. The overall minus sign of the four-$\rho$ block $d^3N_a$ is absorbed into the signs of its terms.
\item Each term is written as in the list, one per line, with its source and its position in that block in grey at the right. $4\rho_a$, $4\rho_b$: from four $\rho$'s ($d^3N_a$, $d^3N_b$); I.II.1: Group I Row II part 1; and so on.
\item Terms of the same colour cancel each other. The last lines of each equation give, for every colour, the sum of the coefficients times the common colour structure; it is zero. There $(\cdots)$ stands for the prefactor and the kernel of the first line.
\item The charges are compared in their fully symmetrised form, so $\tfrac12\rho\{\rho,\rho\}$ counts as $\rho\rho\rho$ and $\{\rho,\rho\}\rho$ as $2\rho\rho\rho$. A term with its $f$ indices in another order picks up the sign of the reordering.
\item In the sums, $\big[U(w)U(w')\rho(x')\big]^c=U^{cd}(w)U^{de}(w')\rho^e(x')$.
\end{itemize}

\begin{center}\begin{tabular}{lccc}\toprule
class & Wilson lines & terms & sets (each sums to $0$)\\\midrule
''' % (len(E), nnet))
    for cls, nets in R.items():
        for nU in (2, 3, 4):
            ns = [n for n in nets if n['terms'][0][0]['nU'] == nU]
            w(r'$%s$ & %d & %d & %d\\' % (NAME[cls], nU, sum(len(n['terms']) for n in ns), len(ns)) + '\n')
        w(r'\midrule' + '\n')
    w(r'total & & %d & %d\\\bottomrule\end{tabular}\end{center}' % (len(E), nnet) + '\n')

    w(r'\section{Cancellation of all three $\rho$'"'"'s}' + '\n')
    for cls, nets in R.items():
        w(r'\subsection{$%s$ terms}' % NAME[cls] + '\n')
        for nU in (2, 3, 4):
            ns = [n for n in nets if n['terms'][0][0]['nU'] == nU]
            assert all(e['nU'] == nU for n in ns for e, _, _ in n['terms'])
            w(r'\subsubsection{%s Wilson lines}' % WORD[nU] + '\n')
            w(r'%d terms in %d colours.' % (sum(len(n['terms']) for n in ns), len(ns)) + '\n')
            w(r'\begin{align}' + '\n')
            w(r'\frac{d^{3}N^{(%s)}_{%d\,U}}{d^{3}k}={}&%s%s\notag\\' % (NAME[cls], nU, PREF, KTEX[cls]) + '\n')
            rows = []
            for k, n in enumerate(ns):
                col = COLS[k % len(COLS)]
                for e, cs, legs0 in n['terms']:
                    sign, mag, body = written_parts(e)
                    rows.append((col, '%s%s%s' % ('+' if sign > 0 else '-', mag, body), src(e)))
            for r_, (col, txt, lab) in enumerate(rows):
                lead = r'&\times\Bigg[' if r_ == 0 else r'&\qquad'
                close = r'\;\Bigg]' if r_ == len(rows) - 1 else ''
                w(r'%s{\color{%s}%s}%s\tag*{%s}\\' % (lead, col, txt, close, lab) + '\n')
            # the sums
            for k, n in enumerate(ns):
                col = COLS[k % len(COLS)]
                parts = ''.join(num(cs) for _, cs, _ in n['terms'])
                tot = sum(cs for _, cs, _ in n['terms'])
                lead = r'={}&\big(\cdots\big)\Bigg[' if k == 0 else r'&\qquad'
                end = r'\Bigg]=0.' if k == len(ns) - 1 else r'\notag\\'
                w(r'%s{\color{%s}\big(%s\big)\,%s}%s' % (lead, col, parts.lstrip('+'), net_tex(n['legs']), end) + '\n')
                assert abs(tot) < 1e-12
            w(r'\end{align}' + '\n')

    w(r'\section{Numerical check}\label{sec:num}' + '\n')
    w(r'''Each term is reduced to its symmetrised charge part. The three charges are placed on three points in all $3!$ ways, the terms are collected by kernel monomial, and the colour factors are evaluated with random SU(3) adjoint Wilson lines and random charges, in three configurations (\texttt{check\_cancel.py}).
\begin{center}\begin{tabular}{llcc}\toprule
Wilson lines & set & kernel monomials & largest residual (largest single term)\\\midrule
''')
    for name, d in chk.items():
        lab = r'$U=U^T$ (your convention)' if 'symmetric' in name else r'general $U$, as written'
        for k, v in d.items():
            kk = {'all terms': 'all', 'alpha': r'$\alpha$', "alpha'": r"$\alpha'$", 'beta': r'$\beta$', "beta'": r"$\beta'$"}[k]
            w(r'%s & %s & %d & $%.1e$ ($%.1e$)\\' % (lab, kk, v[0]['monomials'], max(x['max_residual'] for x in v), v[0]['max_single']) + '\n')
        w(r'\midrule' + '\n')
    w(r'''\bottomrule\end{tabular}\end{center}
For a general $U$, with every chain read as $U(p_1)\cdots U(p_n)\rho$ (oriented from $f$ to the charge), the sum over each class vanishes to $10^{-14}$ (\texttt{general\_u.py}).
\end{document}
''')
    open(os.path.join(HERE, 'ThreeRho_cancellation_v3.tex'), 'w').write(''.join(out))


if __name__ == '__main__':
    build()
    print('wrote ThreeRho_cancellation_v3.tex')
