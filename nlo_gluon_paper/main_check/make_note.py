r"""Writes Calculation_main_check.tex from pending.json and eq_diff.json."""
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
PEND = json.load(open(os.path.join(HERE, 'pending.json')))
EQ = json.load(open(os.path.join(HERE, 'eq_diff.json')))
STEPS = {g: json.load(open(os.path.join(HERE, '..', 'full_check', 'step_%s.json' % g))) for g in ('Group_I', 'Group_II', 'Group_III')}


def eq_rows(g):
    out = []
    for r in EQ[g]:
        st = r['status']
        if st == 'changed':
            st = r'rewritten; checked step by step (item 4 below)'
        else:
            st = st.replace(' and the measure', ' up to the measure').replace('same', 'identical').replace('up to the overall sign', r'up to the overall sign (\textbf{M1} of the Group II/III check, the fix you made)')
        out.append(r'%s & %d--%d & %s\\' % (r['eq'], r['lines'][0], r['lines'][1], st))
    return '\n'.join(out)


DOC = r'''\documentclass[10pt]{article}
\usepackage[a4paper,margin=0.6in]{geometry}
\usepackage{amsmath,amssymb,mathtools,booktabs,longtable,array}
\usepackage[dvipsnames]{xcolor}
\usepackage[colorlinks=true,linkcolor=blue]{hyperref}
\newcommand{\new}[1]{{\color{blue}#1}}
\newcommand{\old}[1]{{\color{red}#1}}
\def\mf#1{\mathbb{#1}}
\def\v#1{\mathbf{#1}}
\allowdisplaybreaks
\title{Your new \texttt{Calculation\_main.tex}: check, and consistency with Group I, II, III}
\date{}
\begin{document}
\maketitle

\section{Answer}
\begin{itemize}
\item \textbf{Your Group I, II, III cross sections agree with \texttt{Calculation\_main.tex}}, once the corrections listed earlier (\texttt{groups23/main/Main\_corrections.pdf}) are made. Your Group files already use the correct forms everywhere they differ from what is written in \texttt{Calculation\_main.tex} (Sec.~\ref{sec:chain}).
\item \textbf{In the new file, two of the earlier corrections are done:}
\begin{itemize}
\item E9, \texttt{stseve} [line 355]: the energy denominators are now $\big(\frac{\v p^2}{2p^+}\big)\big(\frac{\v r^2}{2r^+}+\frac{\v q^2}{2q^+}\big)\big(\frac{\v k^2}{2k^+}\big)$. Correct.
\item M8 [line 398]: ``Setting $\Lambda=0$, wherever possible'' $\to$ ``Simplifying''. Correct: nothing is set to zero in the next line.
\end{itemize}
Apart from these two lines the file is identical to the one checked before.
\item \textbf{All other corrections are still needed} (Sec.~\ref{sec:pending}). Seven of them change results if \texttt{Calculation\_main.tex} is used as written (E1--E7). Your Group files are not affected by them, because the Group files were built from the correct forms.
\item \textbf{Remove the text after \texttt{\textbackslash end\{document\}}} (lines 1965--2031): it is a copy of your message to me. \LaTeX\ ignores it, but it should not be in the file.
\end{itemize}

\section{What is still to be corrected in \texttt{Calculation\_main.tex}}\label{sec:pending}
Line numbers are those of your new file (the same as before). The marked-up copy is \texttt{Calculation\_main\_corrected.pdf} (red: as written, blue: corrected).
{\small
\begin{longtable}{@{}p{0.05\textwidth}p{0.20\textwidth}p{0.38\textwidth}p{0.30\textwidth}@{}}\toprule
& where & what & your Group files\\\midrule
\multicolumn{4}{@{}l}{\textbf{Errors that change results}}\\
E1 & \texttt{onegrho12} [683] & overall sign: $-\frac{ig^2f^{abc}}{8\pi^3}\to+\frac{ig^2f^{abc}}{8\pi^3}$; its $\xi$-form below is right & Group I Row III uses the right sign\\
E2 & \texttt{Eq:onegrhoX1} [701] & overall sign: $+\frac{ig^2f^{abc}}{4\pi^3}\to-\frac{ig^2f^{abc}}{4\pi^3}$ & Group I Row III uses the right sign\\
E3 & $\mf B_3$ [1276--1300] & the non-instantaneous $\Theta(p^+-k^+-\Lambda)$ entry needs a minus sign; then every $1/\Lambda$ term cancels & right in Group I Row III: its $\frac1\Lambda$ lines add up to zero ($3\times10^{-16}$)\\
E4 & $\mf A^{(3)}_2$, momentum space [1252] & counts $\frac{\pi^2}{3}-\frac{67}{9}$ twice and misses the pole of $\Psi^B$; use the bracket of \texttt{grhoposition} & Group I Row II part 1 uses \texttt{grhoposition} (right)\\
E5 & unitarity, $\mf A^{(3)\dagger}_1$ [1087--1101] & $\mf N^{(2)}\mf A^{(1)}=\mf A^{(1)}\mf N^{(2)}$ is false; keep $-[\mf N^{(2)},\mf A^{(1)}_2]$ & its commutator is a two-$\rho$ $\log$ term of Group I Row II; with it the two-$\rho$ sum is JIMWLK$\otimes$LO\\
E6 & $\mf B_2$, $(k^+,p^+)$ forms [1260, 1268] & $\Theta(k^+-\Lambda-p^+)\to\Theta(p^+-\Lambda)\Theta(k^+-\Lambda)\Theta(\vee-k^+-p^+)$ & Groups II, III use $\Lambda<p^+<\vee-k^+$ (right)\\
E7 & Appendix F [1857--1922] & JIMWLK formula: $M$ at the points of $\mathcal K(u,v,t)$; $dN^{(5)}$ with $+$ & M2$'$ of \texttt{Checks\_JIMWLK\_BFKL.pdf}\\
\multicolumn{4}{@{}l}{\textbf{Typos (the next line is right)}}\\
E8 & \texttt{part2} [270] & $\delta^{ni}\to\delta^{ji}$ & not used\\
E10 & \texttt{onegrho22} [731] & ``$+$'' $\to$ product & not used\\
E11 & \texttt{A1coeff} [1247] & $\mf A^a_{2i}(k^+,\v z)\to\mf A^a_{2i}(k^+,-\v z)$ & right in the Groups\\
\multicolumn{4}{@{}l}{\textbf{Exact forms} (replace every $\Lambda\to0$, $\vee\to\infty$ step; only $\vee>p^+,k^+>\Lambda$)}\\
X1--X4 & \texttt{gtwoloop1}, \texttt{instgg}, \texttt{onerho2}, \texttt{tworho2}, \texttt{combinedgrhorho} & exact $\xi$ integrals, $\ln\frac{\Lambda}{k^+}\to\ln\frac{\Lambda}{k^+-\Lambda}$ & the Groups keep $\Lambda$, $\vee$ exact\\
X5 & \texttt{threerhoafter} [411], after \texttt{zetainteg} [432] & ``Using the argument when $\vee$ is very large or $\Lambda$ is very small'': not needed, the $\frac1\Lambda-\frac1\vee$ terms cancel exactly (your M8 fix is the first part of X5) & three-$\rho$ terms\\
X6, X7 & \texttt{grhoposition}, \texttt{combinedgrhorho(pos)} [496--540] & exact $\xi$ limits $[\frac{k^+}{\vee},\frac{k^+}{k^++\Lambda}]$; keep $\Lambda$ & Group I Row II part 1 has $\log(1-\frac{k^+}{\vee})$: exact\\
\multicolumn{4}{@{}l}{\textbf{Only for a general adjoint $U$} (with $U^{bc}=U^{cb}$ written $=$ corrected)}\\
G1--G4 & $\mathbb X_1$, \texttt{eq: A1a}, $\bar{\mf C}\mf A$, \texttt{eq:CS2} & index placements & as in your convention\\
\multicolumn{4}{@{}l}{\textbf{Minor}}\\
M1--M7 & & notation, typography, order-$g^4$ normalisation ($\frac38$, not $-\frac18$; not used) & --\\
\bottomrule
\end{longtable}}
\noindent Programmatically (\texttt{pending.py}): the ``as written'' text of every item above is still in your file; only E9 and the M8 sentence are no longer found.

\section{\texttt{Calculation\_main.tex} $\to$ your Group I, II, III}\label{sec:chain}
The chain of checks:
\begin{enumerate}
\item \textbf{Wave functions and coefficients} of \texttt{Calculation\_main.tex} ($\Psi$'s, $\mf A^{(1)}$, $\mf A^{(3)}$, $\mf B_2$, $\mf B_3$, $\mf C_1$, $\mathbb X_{1,2}$, Eq.~(CS2)) were rebuilt from the matrix elements m1--m5 and compared point by point (\texttt{groups23/main}). The results are the corrected forms of Sec.~\ref{sec:pending}.
\item \textbf{The first lines of every row} of your Group files (the products of these coefficients, before any $p^+$, $w$, $w'$ integral) were compared numerically with that first-principles evaluation: Group I in \texttt{group1} (A3, B3, B2, C1 and the cross-section formula), Groups II and III in \texttt{groups23} (random points, random adjoint Wilson lines, random charges; exact $p^+$ integrals).
\item \textbf{Your current Group files have the same first lines.} Equation by equation (\texttt{eq\_diff.py}; colours, spaces and the lists of points under $\int$ ignored):
\end{enumerate}
{\small
\begin{longtable}{@{}lll@{}}\toprule
Eq. & lines & against the version checked against \texttt{Calculation\_main.tex}\\\midrule
\multicolumn{3}{@{}l}{\textbf{Group I} (checked: \texttt{group1/input/Group\_I\_v2.tex})}\\
@@G1@@
\multicolumn{3}{@{}l}{\textbf{Group II} (checked: \texttt{groups23/input/Group\_II.tex})}\\
@@G2@@
\multicolumn{3}{@{}l}{\textbf{Group III} (checked: \texttt{groups23/input/Group\_III.tex})}\\
@@G3@@
\bottomrule
\end{longtable}}
\begin{enumerate}\setcounter{enumi}{3}
\item \textbf{Every later step} of your current files (the rewritten equations above included) was checked numerically, from one equation to the next (\texttt{full\_check/step\_check.py}, \texttt{Checks\_JIMWLK\_BFKL.pdf}, Sec.~2): all agree to $10^{-16}$--$6\times10^{-12}$ (plus prescriptions up to $O(\Lambda)$), apart from the typo M2 (Group II Eq.~(1.16)). The changed equations of Group II Row IV (1.15)--(1.21) are the $d_\perp$-general form (\texttt{groups23/rowIV\_d}); those of Group I are the $\log$ rewritings and E1.
\end{enumerate}
So the Group files follow from \texttt{Calculation\_main.tex} with its corrections. Group III Row III's first unnumbered line still has the old sign (M3 of \texttt{Checks\_JIMWLK\_BFKL.pdf}); the explicit lines below it are right.

\section{Files}
\begin{itemize}
\item \texttt{Calculation\_main\_corrected.pdf}: your new file with every remaining correction marked in place (red/blue); its source \texttt{Calculation\_main\_corrected.tex} keeps your preamble (revtex4-2, tikzit); the PDF is built from a copy with the article class and the diagrams as boxes.
\item \texttt{pending.py} $\to$ \texttt{pending.json}; \texttt{eq\_diff.py} $\to$ \texttt{eq\_diff.json}.
\item The explanations of E1--E11, X1--X7, G1--G4: \texttt{groups23/main/Main\_corrections.pdf}.
\end{itemize}
\end{document}
'''
out = DOC.replace('@@G1@@', eq_rows('Group_I')).replace('@@G2@@', eq_rows('Group_II')).replace('@@G3@@', eq_rows('Group_III'))
open(os.path.join(HERE, 'Calculation_main_check.tex'), 'w').write(out)
print('written Calculation_main_check.tex')
