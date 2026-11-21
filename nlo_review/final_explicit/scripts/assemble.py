# Assemble Final_terms_explicit.tex from the author's file + generated/hand-written parts.
import sys, os, re
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from texlib import LINES, phase, L
OUT = '/home/user/correlations/nlo_review/final_explicit/Final_terms_explicit.tex'
pre = L(1, 33)   # author's preamble up to \date{}
add = r"""
% ---- added for Final_terms_explicit.tex (colours of the cancellations; red = changes, notes) ----
\definecolor{cI4}{rgb}{0.0,0.50,0.0}      % exact: Group I Row IV vs Group III Row I
\definecolor{cUV}{rgb}{0.0,0.50,0.50}     % exact: 11/3 UV terms, Group I Rows II and V
\definecolor{pA}{rgb}{0.85,0.45,0.0}      % pair A (tail only)
\definecolor{pB}{rgb}{0.75,0.0,0.55}      % pair B (tail only)
\definecolor{pC}{rgb}{0.45,0.2,0.75}      % pair C (tail only)
\definecolor{pD}{rgb}{0.55,0.30,0.10}     % pair D (tail only)
\definecolor{cPr}{rgb}{0.0,0.45,0.75}     % problematic term: (l-E) + 24-term tail (after the w integration)
\definecolor{cJt}{rgb}{0.45,0.45,0.0}     % joint large-distance tail (requires joint integration)
\newcommand{\ph}[1]{{\color{red}#1}}
\newcommand{\NOTE}{{\color{red}\textbf{NOTE}}}
\allowdisplaybreaks
"""
head = pre+add+r"""
\begin{document}
\maketitle
\noindent{\sloppy\color{red}[\textbf{Phase changed to $e^{-ik\cdot(w-w')}$ on the author's instruction (2026-11-20); the author's file had $e^{-ik\cdot(w'-w)}$.} The author's file \texttt{Final\_terms\_author\_2026-11-20.tex} used $e^{-ik\cdot(\text{primed}-\text{unprimed})}$ uniformly (all 57 phases, including the recoil phases; inventory by script, Sec.~\ref{sec:F}); here every phase is $e^{-ik\cdot(\text{unprimed}-\text{primed})}$, the convention of \texttt{fixed\_files/Final\_terms\_fixed.tex} and of the paper draft, and every changed phase is printed in red. This is $k\to-k$ in every term: the parts even in $k$ are unchanged, the parts odd in $k$ (symmetrised three-$\rho$ parts) change sign; the terms below are, in this convention, the verified terms.
\textbf{This file:} Secs.~1--4 are the author's BFKL, JIMWLK, DGLAP and running-coupling sections, unchanged except for the phases. Sec.~5 is the author's statement on the finite terms, followed by: Sec.~\ref{sec:conv} conventions and plus prescriptions; Sec.~\ref{sec:A} check of every finite term of the author's file against the verified files (with the \NOTE{}s); Sec.~\ref{sec:B} every finite term written explicitly; Sec.~\ref{sec:C} the cancellations, in colours (legend at the start of that section); Sec.~\ref{sec:D} what remains, with every three-$\rho$ product symmetrised exactly, Group by Group and Row by Row; Sec.~\ref{sec:E} the same list averaged over the projectile, in colour-singlet form; Sec.~\ref{sec:F} numerical checks. \textbf{Differences from the verified terms found in the author's finite terms} (details in Sec.~\ref{sec:A}): (1) the last term of $N_{3b}$ in Group I Row III is not equivalent to the verified recoil form (verified form in red, Eq.~\eqref{eq:B-I3}); (2) the remnant term of Group II Row IV is missing (added in red, Eq.~\eqref{eq:B-II4rem}). Everything else is equivalent to the verified terms (how verified: Sec.~\ref{sec:A}). \textbf{Plus prescriptions} (author's request of 2026-11-21, ``Rewrite all the $\mathcal{E}(s)$ and other things with plus prescription''): every longitudinal integral of the finite terms is written as a $p^{+}$ integral with the author's plus prescriptions ($\big[\frac{1}{p^{+}}\big]_{+}$, $\big[\frac{k^{+}}{p^{+}}\big]_{+}$, $\big[\frac{k^{+}}{k^{+}-p^{+}}\big]_{+}$, $\big[\frac{k^{+}}{p^{+}(p^{+}+k^{+})}\big]_{+}$) and the recoil phases $e^{-ik\cdot(\dots)p^{+}/k^{+}}$, in place of the closed-form functions of the previous version; each replacement is an identity (checked numerically, Sec.~\ref{sec:F}). \textbf{Two prescriptions} (author's instruction of 2026-11-20): Group I Row III is given with principal values (standing-wave dressing, as in the author's file) and with $i\epsilon$ (out-state dressing, $\eta=+1$), which adds the finite term $\Delta N$, Eq.~\eqref{eq:B-DN}; pure glue ($N_{f}=0$).]\par}
"""
body, n1 = phase(L(34, 151).replace(r'\begin{document}', '').replace(r'\maketitle', ''), 7)
fin_intro = L(152, 153)
parts = [open(os.path.join(HERE, p)).read() for p in ['conv.tex', 'partA.tex', 'partB.tex', 'partC.tex', 'partD.tex', 'partE.tex', 'partF.tex'] if os.path.exists(os.path.join(HERE, p))]
tail = '\n'.join(LINES[335:])   # \end{document} and the author's text after it (not compiled)
op = os.path.join(HERE, 's4_operator.log')
optxt = open(op).read().strip() if os.path.exists(op) else 'pending'
import re as _re
m = _re.search(r'= (\d+) cases\): .*?norm < 1\) ([0-9.e+-]+); .*?: ([0-9.e+-]+); control .*?: rel dev ([0-9.]+)', optxt)
if m:
    fmt = lambda g: _re.sub(r'([0-9.]+)e-0*([0-9]+)', lambda mm: mm.group(1)+r'\times10^{-'+mm.group(2)+'}', g)
    optxt = (r'SU(3) fundamental charges as explicit $3\times3$ generators on three sites ($27$-dimensional space), $\delta$-functions as Kronecker deltas of the sites, random Wilson lines (equal at coinciding sites, and independent): for all $27$ assignments of the three charge points of each of the thirteen cubic colour structures to the sites ($%s$ cases), the colour operator of the row as written in Sec.~\ref{sec:B} equals the tensor $C^{abc}$ of Sec.~\ref{sec:D} contracted with the ordered charges, and equals (i)+(ii)+(iii): maximal deviation $%s$ (relative; absolute for operators of norm below one, e.g.~Row I.II at $x=y$, where the operator vanishes). Group III Row VI: author\'s simplified colour operator equals the earlier form with $U(z)U(z)$: $%s$. Control: omitting the contact terms at $x\'=x$ gives a relative deviation $%s$.' % tuple(fmt(g) for g in m.groups())).replace("\\'", "'")
parts = [p.replace('OPERATORCHECK', optxt) for p in parts]
doc = head+body+'\n'+fin_intro+'\n'+'\n'.join(parts)+'\n'+tail
open(OUT, 'w').write(doc)
print('written', OUT, 'phases replaced in Secs. 1-4:', n1)
