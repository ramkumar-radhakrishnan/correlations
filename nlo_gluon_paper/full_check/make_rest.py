r"""Writes Remaining_RC_DGLAP_finite.tex from your new Group_I/II/III.tex (input/):

  1. the remaining terms: every row minus its log(k+/Lambda) (JIMWLK) and log(vee/k+) (BFKL) terms
     (Checks_JIMWLK_BFKL.pdf), Lambda -> 0 wherever possible, vee exact, line 1 of Group III Rows III, V, VI and
     Group II Row IV in its exact rewritten form (../bfkl_rows/exact_check.py);
  2. the running coupling;  3. DGLAP  (the bookkeeping of ../dglap_rc/checks.py, checks.json);
  4. the finite terms, step by step: what is left of each row after 2 and 3;  5. the result.
The remainders are those of ll_checks.py split (exact, checked term by term)."""
import json
import os
import re
import make_ll as ML          # imports ll first (your new files), and the split
import gen as GN
import exact_check as EC       # ../bfkl_rows: the exact rewriting of line 1 (its line 1 = yours, rows_vs_files.py)

HERE = ML.HERE
BLK, SPL = ML.BLK, ML.SPL
tex_num = ML.tex_num
DG = os.path.join(HERE, '..', 'dglap_rc')
CHK = json.load(open(os.path.join(DG, 'checks.json')))
DGTEX = open(os.path.join(DG, 'Finite_terms_DGLAP_RC.tex')).read()
CS = ML.CS


def lam0(t):
    return ML.lam0(t)


def remaining_items(name):
    """(text, style) of the remaining terms of one row"""
    items = []
    for n, s in enumerate(SPL[name]):
        if name in ML.LINE1_ROWS and n == 0:
            items += [(ML.c5(name, t), st) for t, st in (EC.g3_display(name) if name != 'G2 IV' else EC.g2iv_display())]
        elif name in ML.LINE1_ROWS and n == 1:
            items.append((lam0(ML.c5(name, ML.remainder_line2(name))), None))
        elif s.get('power') or s.get('red'):
            items.append((BLK[name].terms[n], 'red'))
        elif (name, n) == ('G1 III N3b', 3):
            items.append((lam0(s['F']), 'red'))
        elif s.get('F') is not None:
            items.append((lam0(ML.c5(name, s['F'])), None))
    return items


NR = {
    'G1 II(1)': r'''Your bracket without its $-2\log\frac{k^+}{\Lambda}$ (JIMWLK). $\log(1-\frac{k^+}{\vee})\log\frac{\vee}{k^+}$ stays here: its coefficient depends on $\vee$, and $0\le-\log(1-\frac{k^+}{\vee})\log\frac{\vee}{k^+}\le\log^22$ for all $\vee>k^+$. The pole bracket times $\frac{11}{3}$ is the $\beta_0$ term (Sec.~\ref{sec:rc}).''',
    'G1 II(2b)': r'''Unchanged: no large logarithm; $R^i$ as in your file.''',
    'G1 III N1': r'''\new{Blue}: what is left of the logarithms after their $\log\frac{k^+}{\Lambda}$, $\log\frac{\vee}{k^+}$ parts are taken out, with $\Lambda\to0$. \old{Red}: the $\frac{1}{\Lambda}$ lines (your green brackets). Lines 1 and 4 here, line 1 of $d^3N_2$ and line 1 of $d^3N_{3a}$ add up to zero exactly, pointwise ($@@POW@@$).''',
    'G1 III N2': r'''As $d^3N_1$. \old{Red}: line 1 is a $\frac{1}{\Lambda}$ line; your two red lines cancel each other.''',
    'G1 III N3a': r'''\old{Red}: line 1 is a $\frac{1}{\Lambda}$ line (see $d^3N_1$). The edge-kernel line has no remaining part: $\log\frac{k^+-\Lambda}{\Lambda}=\log\frac{k^+}{\Lambda}+\log(1-\frac{\Lambda}{k^+})$. $\mathcal A_0$ is $\mathcal A$ at $\Lambda=0$, Eq.~\eqref{eq:A0}.''',
    'G1 III N3b': r'''$\Lambda\to0$. \textbf{C4}: the last line (\old{red}) is, for $\Lambda\to0$, $\log\frac{(z-w)^2}{(x-w)^2}$ times the edge kernel $\frac{(x'-w')\cdot(y-z)}{(x'-w')^2(y-z)^2(z-w)^2}$: not integrable at $z\to w$. Its logarithmic partner was replaced by $\mathcal D$ (F3); its finite part belongs to the $\mathcal D$ region (open, Q3).''',
    'G1 IV': r'''\new{Blue}: the remainders of the two logarithms.''',
    'G1 V(1)': 'Unchanged.', 'G1 V(2a)': 'Unchanged.',
    'G1 V(2b UV)': r'''$2\log(1-\frac{\Lambda}{k^+})\to0$; $-\frac{11}{6}$ is the $\beta_0$ term that cancels the $\frac{11}{3}$ of Row II part 1 (Sec.~\ref{sec:rc}).''',
    'G1 V(2b fin)': r'''$\Lambda\to0$ in the limits of the plus-prescribed integrals (regular at the end points). In line 6, $2\log\frac{k^+-\Lambda}{\Lambda}\to2\log\frac{k^+}{\Lambda}$ (JIMWLK) $+2\log(1-\frac{\Lambda}{k^+})\to0$.''',
    'G2 I(1) A': r'''\new{Blue}: remainders. Line 2 is entirely BFKL ($\log\frac{\vee}{\Lambda+k^+}=\log\frac{\vee}{k^+}-\log(1+\frac{\Lambda}{k^+})\to\log\frac{\vee}{k^+}$): nothing of it is left.''',
    'G2 I(1) B': r'''As the first line.''',
    'G2 I(2a)': r'''$\Lambda\to0$ in the second logarithm.''',
    'G2 I(2b)': r'''\new{Blue}: remainders.''', 'G2 I(2c)': r'''\new{Blue}: remainders.''',
    'G2 I(2d)': r'''\new{Blue}: remainders.''', 'G2 I(2e)': r'''\new{Blue}: remainders.''',
    'G2 I(2f)': r'''\new{Blue}: remainders; $\frac{1}{k^++\Lambda}\to\frac{1}{k^+}$.''',
    'G2 IV': r'''Line 1 is rewritten exactly (\texttt{Checks\_JIMWLK\_BFKL.pdf}, Sec.~3):
\begin{itemize}
\item $\frac{\xi}{(1-\xi)_+}\to\xi+\frac{\xi^2}{(1-\xi)_+}$, Eq.~(I2$\xi$). The difference, $-g(1)\log\frac{\vee}{k^+}$, is (ii) of the BFKL list; it cancels (i), the $\log\frac{\vee}{k^+}$ of line 2.
\item The first \new{blue} line subtracts the $\log(y'-y)^2$ term of $\frac{1-\xi}{\xi}$ at $y'=y$. After the $y'$ integral that term is $\big[\log\frac{\vee}{k^+}-(1-\frac{k^+}{\vee})\big]H_4\big|_{\log\to1}$, Eq.~(I4); its $\log\frac{\vee}{k^+}$ part is the BFKL term $\log\frac{\vee}{k^+}H_4$.
\item The second \new{blue} line is the rest, $+(1-\frac{k^+}{\vee})H_4\big|_{\log\to1}$, which stays here.
\end{itemize}
Line 2: $\big[\log\frac{k^+}{\Lambda}+\log(\frac{\vee}{k^+}-1)\big]\to\log(1-\frac{k^+}{\vee})$. In line 3, $\frac{k^+}{k^++\Lambda}\to1$. Lines 1, 2 with the two blue lines are the DGLAP term (Sec.~\ref{sec:dglap}).''',
    'G3 I': r'''\new{Blue}: remainders.''', 'G3 II': r'''\new{Blue}: remainders.''',
    'G3 III(1)': r'''Line 1 is rewritten exactly (\texttt{Checks\_JIMWLK\_BFKL.pdf}, Sec.~3):
\begin{itemize}
\item $\left[\frac{1}{p^+}\right]_+\to\left[\frac{k^+}{p^+(p^++k^+)}\right]_++\frac{1}{p^++k^+}$, Eq.~(I2). The difference, $-f(0)\log\frac{\vee}{k^+}$, is (ii) of the BFKL list; it cancels (i), the $\log\frac{\vee}{k^+}$ of line 2.
\item In the $-\frac{1}{k^+}\delta_{im}\delta_{jk'}$ term, $G(y')\to G(y')-G(z)$ (\new{blue} line). The $G(z)$ part is (B), Eq.~(I3), in the BFKL list.
\end{itemize}
Line 2: $\big[\log\frac{k^+}{\Lambda}+\log(\frac{\vee}{k^+}-1)\big]\to\log(1-\frac{k^+}{\vee})$.''',
    'G3 III(2)': r'''As the first part, with the $-\frac{p^+}{k^{+2}}$ term in place of $-\frac{1}{k^+}\delta_{im}\delta_{jk'}$. The last \new{blue} line is the part of its $G(z)$ term that is not $\log\frac{\vee}{k^+}$, Eq.~(I5). \textbf{C5}: the measure is $\int_{w,x',y',z,x}$.''',
    'G3 V': r'''As Row III.''', 'G3 VI': r'''As Row III.''',
}
NR['G1 III N1'] = NR['G1 III N1'].replace('@@POW@@', tex_num(CS['power_terms_sum_rel']))


def rows_tex():
    out = []
    for gname, gl, rows in ML.GROUPS:
        out.append(r'\section{Remaining terms: %s}\label{sec:F%s}' % (gname, gl))
        for rname, names, special in rows:
            out.append(r'\subsection{%s}' % rname)
            if special == 'four':
                out.append(r'Only four $\rho$\textquotesingle s: no remaining terms (\texttt{4rho.tex} has only $\big[\log\frac{k^+}{\Lambda}+\log\frac{\vee}{k^+}\big]$ terms).')
                continue
            for name in names:
                sub = ML.SUB[name]
                out.append(r'\paragraph{%s%s.}' % (rname, (', ' + sub) if sub else ''))
                items = remaining_items(name)
                if not items:
                    out.append(r'Nothing is left: this part has only $\log\frac{k^+}{\Lambda}$ and $\log\frac{\vee}{k^+}$ terms.')
                    continue
                if NR.get(name):
                    out.append(NR[name])
                out.append(ML.align(name, items, 'F'))
    return '\n'.join(out)


# ------------------------------------------------------------------------------------------------ RC and DGLAP
def sec(t, a, b):
    i = t.index(a)
    return t[i:t.index(b, i)]


def once(t, a, b):
    assert t.count(a) == 1, (a[:70], t.count(a))
    return t.replace(a, b)


RCDG = sec(DGTEX, r'\section{Running coupling}', r'\section{The finite terms}')
RCDG = once(RCDG, r'The poles of the remaining terms, in units of LO (\texttt{checks.json}):',
            r'The poles of the remaining terms (Secs.~\ref{sec:FgI}--\ref{sec:FgIII}), in units of LO (\texttt{dglap\_rc/checks.json}):')
RCDG = once(RCDG, r'G2 Row IV, lines 1, 2 &', r'G2 Row IV, lines 1, 2 (with the blue lines) &')
RCDG = once(RCDG, r'Group II Row IV, lines 1 and 2 of the remaining terms, with the pole removed',
            r'Group II Row IV, lines 1 and 2 of the remaining terms with the two blue lines (Sec.~\ref{sec:FgII}, Row IV), with the pole removed')
RCDG = once(RCDG, r'together with the $\log\frac{\vee}{k^+}H_4$ subtraction:',
            r'i.e.\ $\big[$your lines 1 and 2 without $\log\frac{k^+}{\Lambda}\big]-\log\frac{\vee}{k^+}H_4$ (the blue lines and $\xi+\frac{\xi^2}{(1-\xi)_+}$ add up to this exactly):')
RCDG = once(RCDG, r'\item The last two lines ($-\log\frac{\vee}{k^+}H_4$) remove the $\log\frac{\vee}{k^+}$ that the $\xi$ integral produces at $\xi\to\frac{k^+}{\vee}$ (it is in the BFKL terms); with them the DGLAP term is finite for $\vee\to\infty$.',
            r'\item The last two lines ($-\log\frac{\vee}{k^+}H_4$) remove the $\log\frac{\vee}{k^+}$ that the $\xi$ integral contains, Eq.~(I4) (it is in the BFKL terms); with them the DGLAP term is finite for $\vee\to\infty$.')
RCDG = once(RCDG, r'Line 3 of the row (linear polarization) is finite and not DGLAP; it is in Sec.~\ref{sec:fin}.',
            r'Line 3 of the row (linear polarization) is finite and not DGLAP; it is a finite term (Sec.~\ref{sec:fin}).')
RCDG = RCDG.replace(r'(Sec.~\ref{sec:aabkl})', r'(arXiv:1411.2869)')
assert 'Finite_terms_split' not in RCDG.replace(r'\_', '_'), 'old file name in the reused text'
assert 'sec:aabkl' not in RCDG

PART1 = r'''\paragraph{Notation.} $a\equiv\frac{k^+}{\vee}$, $\beta_0=\frac{11}{3}N_c$ ($n_f=0$), $d_\perp=2-\epsilon$, $c_0\equiv2e^{-\gamma_E}$, $\ell_r\equiv\log\frac{r^2\mu^2}{c_0^2}=2\gamma_E+\log\frac{r^2\mu^2}{4}$, $\hat E_r\equiv\frac2\epsilon+\ell_r$; $X=x-y$, $X'=x'-y'$, $\Delta=y'-y$. LO in your notation, and LO with a weight $f$:
\begin{align*}
{\rm LO}[f]&\equiv\frac{1}{(2\pi)^{3}}\frac{g^{2}}{2\pi^{2}}\frac{1}{k^{+}}\int_{x,y,x',y'}e^{-ik\cdot(y'-y)}\,f\,\frac{(x-y)\cdot(x'-y')}{(x-y)^{2}(x'-y')^{2}}\left[U^{e'c}(y')-U^{e'c}(x')\right]\left[U^{ec}(y)-U^{ce}(x)\right]\rho^{e'}(x')\rho^{e}(x),
\end{align*}
${\rm LO}\equiv{\rm LO}[1]$; ${\rm LO}_\xi[f]$ is the same with $e^{-ik\cdot(y'-y)}\to e^{-ik\cdot(y'-y)\frac1\xi}$. $P^{(r)}_{gg,+}(\xi)=2N_c\big[\frac{\xi}{(1-\xi)_+}+\frac{1-\xi}{\xi}+\xi(1-\xi)\big]$ with the plus prescription of $[0,1]$ restricted to $[a,1]$ ($\int_a^1\frac{d\xi\,g(\xi)}{(1-\xi)_+}=\int_a^1\frac{g(\xi)-g(1)}{1-\xi}d\xi+g(1)\log(1-a)$), $P_{gg}=P^{(r)}_{gg,+}+\frac{\beta_0}{2}\delta(1-\xi)$. Only three rows have poles or $\beta_0$: Group I Row II part 1, Group I Row V 2b UV, Group II Row IV (lines 1, 2); all other remaining terms are finite. Their conversion to units of LO was checked against your rows before (\texttt{finite\_v2/final\_checks.json}); the bookkeeping is \texttt{dglap\_rc/checks.py}: $\beta_0$ sum $@@B0@@$, poles @@POLES@@.'''

FIN = r'''\section{The finite terms, step by step}\label{sec:fin}
What is left of every row after the evolution (\texttt{Checks\_JIMWLK\_BFKL.pdf}), the running coupling (Sec.~\ref{sec:rc}) and DGLAP (Sec.~\ref{sec:dglap}):
\begin{enumerate}
\item \textbf{Group I Row II part 1}, Eq.~\eqref{eq:F:G1II1}. Its pole bracket times $\frac{11}{3}$: the $\beta_0$ term, cancelled by Row V 2b UV, Eq.~(\ref{eq:b1})$+$(\ref{eq:b2})$=0$. Its constant $-\frac{67}{9}+\frac{\pi^2}{3}$: Eq.~\eqref{eq:cmw}. Its pole bracket times $\log(1-\frac{k^+}{\vee})$: Q1, Eq.~\eqref{eq:q1}. Left as a finite term:
\begin{equation}
\frac{d^3N_{V_\vee}}{d^3k}=-\frac{\alpha_s}{2\pi}N_c\log\Big(1-\frac{k^+}{\vee}\Big)\log\Big(\frac{\vee}{k^+}-1\Big)\,{\rm LO}\qquad(\text{with c.c.}),
\label{eq:Vvee}
\end{equation}
the sum of its $-\log^2(1-\frac{k^+}{\vee})$ and $-\log(1-\frac{k^+}{\vee})\log\frac{\vee}{k^+}$ terms. It is bounded for every $\vee>k^+$ and vanishes for $\vee\to\infty$.
\item \textbf{Group I Row II part 2}: $d^3N_{2a}$ has only $\log\frac{k^+}{\Lambda}$ and $\log\frac{\vee}{k^+}$ terms, nothing is left; $d^3N_{2b}$, Eq.~\eqref{eq:F:G1II2b}, is finite and unchanged.
\item \textbf{Group I Row III}, Eqs.~\eqref{eq:F:G1IIIN1}--\eqref{eq:F:G1IIIN3b}: finite (red lines removed: the $\frac{1}{\Lambda}$ lines add up to zero, your red lines cancel, C4).
\item \textbf{Group I Row IV}, Eq.~\eqref{eq:F:G1IV}; \textbf{Row V} parts 1, 2a, 2b finite, Eqs.~\eqref{eq:F:G1V1}, \eqref{eq:F:G1V2a}, \eqref{eq:F:G1V2bfin}: finite.
\item \textbf{Group I Row V 2b UV}, Eq.~\eqref{eq:F:G1V2bUV}: nothing is left. $-\frac{11}{6}$ is the $\beta_0$ term (Sec.~\ref{sec:rc}); $2\log(1-\frac{\Lambda}{k^+})\to0$.
\item \textbf{Group II Row I}, Eqs.~\eqref{eq:F:G2I1A}--\eqref{eq:F:G2I2f}: finite.
\item \textbf{Group II Row IV}, Eq.~\eqref{eq:F:G2IV}: lines 1, 2 and the two blue lines are the DGLAP term, Eq.~\eqref{eq:dglap}. Line 3 (linear polarization, $\frac{k^ik^j}{k^2}$) is left as a finite term.
\item \textbf{Group III Rows I, II, III, V, VI}, Eqs.~\eqref{eq:F:G3I}--\eqref{eq:F:G3VI}: finite. In Rows III, V, VI the $p^+$ integrands of line 1 fall like $\frac{1}{p^{+2}}$; the rows are finite for $\vee\to\infty$.
\end{enumerate}
All these terms are finite for $\epsilon\to0$, $\Lambda\to0$ and $\vee\to\infty$, apart from the open finite part of the $\mathcal D$ region (Q3). The terms with three $\rho$\textquotesingle s among them need $\langle\rho\rho\rho\rangle$ (zero for a Gaussian projectile) or their commutator (two-$\rho$) parts.

\section{Result}\label{sec:res}
\begin{equation}
\boxed{\begin{aligned}
\frac{d^3N_h}{d^3k}={}&\int_{k^+/\vee}^1\frac{d\zeta}{\zeta^2}\,D_{h/g}(\zeta,\mu_F)\,\frac{d^3N_g}{d^3k}\Big(\frac{k^+}{\zeta},\frac{\mathbf k}{\zeta}\Big),\\
\frac{d^3N_g}{d^3k}={}&{\rm LO}\big|_{\alpha_s(\mu_R)}+\frac{\alpha_s\beta_0}{4\pi}\log\frac{\mu_R^2}{\mu_F^2}\,{\rm LO}+\frac{\alpha_s}{2\pi}N_c\Big[\frac{67}{9}-\frac{\pi^2}{3}\Big]{\rm LO}\qquad\text{(running coupling, Sec.~\ref{sec:rc})}\\
&+\mathcal N\Big[\log\frac{k^+}{\Lambda}\,{\rm JIMWLK}\otimes{\rm LO}+\log\frac{\vee}{k^+}\,{\rm BFKL}\otimes{\rm LO}\Big]\qquad\text{(evolution, \texttt{Checks\_JIMWLK\_BFKL.pdf})}\\
&+\frac{d^3N_{\rm DGLAP}}{d^3k}\qquad\text{(Eq.~\eqref{eq:dglap})}\\
&+\frac{\alpha_sN_c}{4\pi}\log\Big(1-\frac{k^+}{\vee}\Big)\,{\rm LO}\big[\hat E_X+\hat E_{X'}\big]\qquad\text{(Q1, Eq.~\eqref{eq:q1})}\\
&+\frac{d^3N_{V_\vee}}{d^3k}+\text{the finite terms of Sec.~\ref{sec:fin}}+O(\alpha_s^3).
\end{aligned}}
\label{eq:result}
\end{equation}
\begin{itemize}
\item $\mathcal N\big[\dots\big]$ is your JIMWLK formula (corrected, M2'; $\delta Y=\log\frac{k^+}{\Lambda}$) plus your BFKL formula ($\delta Y=\log\frac{\vee}{k^+}$), charge-averaged and with $\int_{w,w'}e^{-ik\cdot(w'-w)}$.
\item Q1 and $V_\vee$ vanish for $\vee\to\infty$ at fixed $k^+$. Both come from the limit $\xi\ge\frac{k^+}{\vee}$ in the virtual term of Group I Row II; if the virtual term is taken with $\int_0^1d\xi$ (as in AABKL, arXiv:1411.2869), both are absent.
\end{itemize}

\section{Open points}\label{sec:open}
\begin{itemize}
\item[Q1] The pole $\propto\log(1-\frac{k^+}{\vee})$, Eq.~\eqref{eq:q1}: neither DGLAP nor running coupling (see above).
\item[Q2] A finite four-$\rho$ term is not in your files (\texttt{4rho.tex} has only the $\big[\log\frac{k^+}{\Lambda}+\log\frac{\vee}{k^+}\big]$ terms).
\item[Q3] The finite part of the $\mathcal D$ region (the edge-kernel lines of Group I Row III, C4 and F3).
\end{itemize}
'''

DEFS = r'''\paragraph{Plus prescriptions.} Yours, and two more (Eqs.~(I2), (I2$\xi$)):
\begin{align*}
\int dp^+\left[\frac{1}{p^+}\right]_+f&=\int dp^+\,\frac{f(p^+)-f(0)}{p^+},&
\int dp^+\left[\frac{k^+}{p^+(p^++k^+)}\right]_+f&=\int dp^+\,\frac{k^+\,[f(p^+)-f(0)]}{p^+(p^++k^+)},\\
\int\frac{d\xi}{\xi^2}\,\frac{\xi}{(1-\xi)_+}\,g&=\int\frac{d\xi}{\xi(1-\xi)}\,[g(\xi)-g(1)],&
\int\frac{d\xi}{\xi^2}\,\frac{\xi^2}{(1-\xi)_+}\,g&=\int\frac{d\xi}{1-\xi}\,[g(\xi)-g(1)].
\end{align*}'''

DOC = ML.PREAMBLE + r'''
\title{After JIMWLK and BFKL: the remaining terms, the running coupling, DGLAP\\ and the finite terms of your Group I, II, III}
\date{}
\begin{document}
\maketitle
\noindent Companion of \texttt{Checks\_JIMWLK\_BFKL.pdf} (whose equations (I1)--(I5) are cited here). Your rows, without their $\log\frac{k^+}{\Lambda}$ (JIMWLK) and $\log\frac{\vee}{k^+}$ (BFKL) terms: \new{blue} = changed or added, \old{red} = removed or cancelling. Only $\vee>k^+$ and $p^+>\Lambda$ are used; $\Lambda\to0$ wherever possible, $\vee$ exact; every change is an exact identity. Groups I and III with c.c., Group II without.

\section{What is done}\label{sec:what}
\begin{enumerate}
\item \textbf{Remaining terms} (Secs.~\ref{sec:FgI}--\ref{sec:FgIII}). Each row $=\log\frac{k^+}{\Lambda}\mathfrak T+\log\frac{\vee}{k^+}\mathfrak P+R$ exactly; here is $R$.
\begin{itemize}
\item Every logarithm is replaced by its remainder (Sec.~3 of \texttt{Checks\_JIMWLK\_BFKL.pdf}), with $\Lambda\to0$: $\big[\log\frac{k^+}{\Lambda}+\log(\frac{\vee}{k^+}-1)\big]\to\log(1-\frac{k^+}{\vee})$, $\log\frac{\vee}{\Lambda+k^+}\to0$, $\log\frac{(\vee-k^+)A+k^+B}{\Lambda A+k^+B}\to\log\frac{(1-\frac{k^+}{\vee})A+\frac{k^+}{\vee}B}{B}$, \dots
\item Line 1 of Group III Rows III, V, VI and Group II Row IV is rewritten exactly, (I2)--(I5): what it contains of $\log\frac{\vee}{k^+}$ is in the BFKL list.
\item $\Lambda\to0$ everywhere else: $\mathcal A\to\mathcal A_0$; the limits of the plus-prescribed integrals ($\int_\Lambda^{k^+-\Lambda}\to\int_0^{k^+}$, $\frac{k^+}{k^++\Lambda}\to1$); $\frac{1}{k^++\Lambda}\to\frac{1}{k^+}$; $2\log(1-\frac{\Lambda}{k^+})\to0$.
\item Removed (\old{red}): the $\frac{1}{\Lambda}$ lines of Group I Row III (they add up to zero), your two red lines of $d^3N_2$ (they cancel), and the edge-kernel line of $d^3N_{3b}$ (C4).
\end{itemize}
Each of the @@NCHK@@ pointwise-evaluable terms equals its three parts to $@@ID@@$, and the remainder has the $\Lambda\to0$ limit written here ($@@L0@@$ at $\Lambda=10^{-9}$; \texttt{ll\_checks.py split}). The rewriting of line 1 is checked in \texttt{bfkl\_rows/exact\_check.py} ($\le@@ECID@@$). These terms stay finite for $\vee\to\infty$ and $\Lambda\to0$, apart from the open finite part of the $\mathcal D$ region (Q3).
\item \textbf{Running coupling} (Sec.~\ref{sec:rc}) and \textbf{DGLAP} (Sec.~\ref{sec:dglap}) are taken out of Group I Row II part 1, Group I Row V 2b UV and Group II Row IV, with the coupling counterterm of LO.
\item \textbf{The finite terms} are what is left (Sec.~\ref{sec:fin}); Sec.~\ref{sec:res} is the result.
\end{enumerate}
@@DEFS@@
\begin{align}
@@A0@@
\label{eq:A0}
\end{align}

@@ROWS@@

@@RCDG@@

@@FIN@@
\end{document}
'''


def main():
    pk = [k for k in CHK if k.startswith('pole -')][0]
    fk = [k for k in CHK if k.startswith('finite -')][0]
    part1 = (PART1.replace('@@B0@@', '$%s$' % CHK['beta0: G1 II part 1 + G1 V 2b UV'])
             .replace('@@POLES@@', r'$-\frac{\alpha_s}{2\pi}\big(-\frac{2}{\epsilon}\big)\big[P_{gg}\otimes{\rm LO}-N_c\log(1-a){\rm LO}\big]=%s$, finite parts $-$(CMW $+V_\vee+$ Q1 $+$ DGLAP) $=%s$' % (CHK[pk], CHK[fk])))
    out = DOC.replace('@@DEFS@@', DEFS).replace('@@A0@@', ML.A0).replace('@@ROWS@@', rows_tex())
    rc = RCDG.replace(r'\section{Running coupling}\label{sec:rc}', r'\section{Running coupling}\label{sec:rc}' + '\n' + part1, 1)
    out = out.replace('@@RCDG@@', rc).replace('@@FIN@@', FIN)
    for k, v in {'@@NCHK@@': str(CS['n_checked']), '@@ID@@': tex_num(CS['worst_identity']),
                 '@@L0@@': tex_num(CS['worst_Lambda_to_0']), '@@ECID@@': tex_num(ML.EC_ID)}.items():
        out = out.replace(k, v)
    assert '@@' not in out, out[out.index('@@') - 100:out.index('@@') + 50]
    open(os.path.join(HERE, 'Remaining_RC_DGLAP_finite.tex'), 'w').write(out)
    print('written Remaining_RC_DGLAP_finite.tex')


if __name__ == '__main__':
    main()
