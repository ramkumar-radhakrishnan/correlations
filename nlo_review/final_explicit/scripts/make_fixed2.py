# make_fixed2.py -- builds Final_terms_explicit_2_fixed.tex from the author's file Final_terms_explicit_2_author.tex
# (2026-11-23): (A) corrections of Secs. 1-7 in red with short red notes; (B) new Sec. 8 "Final result with factorised
# Wilson-line structures": Sec. 6 (rho-explicit) and Sec. 7 (singlet) transformed row by row: colour structures written
# with differences of Wilson lines (s16_sec8data.py, s16c_singlet.json; checked in s16c, s16e) and the longitudinal
# integrals written with xi (checked in s16d).  Every replacement asserts its number of occurrences.
# Run: python3 -I scripts/make_fixed2.py  (after s16c_sec8_colour.py)
import os, re, json, sys
HERE = os.path.dirname(os.path.abspath(__file__)); TOP = os.path.dirname(HERE); sys.path.insert(0, HERE)
from s16_sec8data import CUBIC8, QUAD8, TEXPT
from s16_lib import bracket_after
SRC = open(os.path.join(TOP, 'Final_terms_explicit_2_author.tex')).read()
SING = json.load(open(os.path.join(HERE, 's16c_singlet.json')))
def rep(s, old, new, n=1):
    c = s.count(old)
    assert c == n, ('count %d != %d for %r' % (c, n, old[:120]))
    return s.replace(old, new)
def rrep(s, pat, new, n):
    r, c = re.subn(pat, new, s)
    assert c == n, ('regex count %d != %d for %r' % (c, n, pat[:120]))
    return r
RED = lambda t: r'{\color{red}' + t + '}'
s = SRC
# ======================================================================================== (A) fixes of Secs. 1-7
s = rep(s, r'\allowdisplaybreaks', r'''\allowdisplaybreaks
\newcommand{\fixnote}[1]{\par\smallskip\noindent\fcolorbox{red}{white}{\parbox{\dimexpr\linewidth-2\fboxsep-2\fboxrule\relax}{\color{red}\small\textbf{NOTE (2026-11-23).} #1}}\par\smallskip}''')
SUMMARY = r'''
{\color{red}
\section*{Corrections of 2026-11-23 (in red) and the new Sec.~\ref{sec:F}}
\noindent This is the author's file of 2026-11-23 with corrections in red (each with a short red note) and a new Sec.~\ref{sec:F}, ``Final result with factorised Wilson-line structures'' (author's request: ``write the Wilson lines $[U-U]\times[U-U]$ \dots\ and also write longitudinal momentum fraction in terms of $\xi$''). Secs.~1--4 are unchanged except one symbol in Sec.~2.

\noindent\textbf{Items of our previous check (\texttt{Final\_terms\_explicit\_1\_checked.pdf}):}
\begin{enumerate}
\item Recoil factor $\frac{p^{+}}{k^{+}}$ outside the exponent in the last term of Eq.~\eqref{eq:B-I3}: \textbf{fixed by the author}.
\item $i\epsilon$ term $\Delta N$ of Group I Row III: superseded by the author's decision (principal value only, 2026-11-22); the averaged $\Delta N$ rows have been \textbf{deleted by the author}; the principal-value prescription is now stated in the conventions (Sec.~\ref{sec:conv2}).
\item $\xi$ undefined in Group I Row V (second part), Eq.~\eqref{eq:B-I5two}: \textbf{not fixed}; the definition $\xi\equiv p^{+}/k^{+}$ is added in red.
\item Eq.~\eqref{eq:D-ord}: (a) the second factor $\frac{i}{2}C^{abc}$ in the third contact term: \textbf{not fixed}, corrected in red; (b) the subscript ${}_{\rm (iii)}$: \textbf{fixed by the author}; (c) ``three Wilson lines'': \textbf{not fixed}, corrected to ``three colour charges $\rho$''.
\item The one-$\rho$ (double commutator) terms are not written in Sec.~6: the sentence saying so is \textbf{missing}; added in red.
\item Conventions (phase, plus prescriptions, distributions, limits and their order, prescription of Row III, pure glue): \textbf{missing}; re-inserted (Sec.~\ref{sec:conv2}, red).
\item Terms that cancel only in their large-distance tails must be integrated together (pairs A--D; the $w$ integral before the $z$ integral in Row I.III; the joint tail of Rows II.I, II.IV, III.III): \textbf{missing}; added (Sec.~\ref{sec:conv2} and the start of Sec.~6).
\item Sec.~7 was commented out: \textbf{fixed by the author} (all rows compiled); the opening $\langle$ of the averaging formula: \textbf{not fixed}, corrected (also in the new definition of $\sigma$); the averaging rule is added in red.
\item UV: no change required.
\end{enumerate}
\noindent\textbf{New findings in this check (corrected in red):} (a) In Sec.~6 the outer brackets around the colour tensor are missing in Eqs.~\eqref{eq:D-II.I-3}, \eqref{eq:D-II.I-2}, \eqref{eq:D-III.II-3}, \eqref{eq:D-III.II-2}: as typeset the factor $2$ (resp.\ $i$, $\frac{i}{2}$) and the kernel in front multiply only the first of the two terms of the tensor and the charges only the second (for Group II Row I the colour tensor read as typeset differs by $32\%$, \texttt{s16e}). The same was already in the previous version of the author's file; our previous TeX comparison ignored brackets and did not find it. (b) Sec.~7 defines $\sigma(x,y,z)$ but the equations use $o(\dots)$: replaced by $\sigma$ (red). (c) In the contact (two-$\rho$) equations of Sec.~7 the terms of the bracket follow each other without a sign between them; a red $+$ is inserted (the signs are in the coefficients $(\pm i)$, $(\pm\frac{i}{2})$; the same convention was used in our verified file). (d) Sec.~2: ``where $\mathcal{K}(x,y,z)$'' should be $K(x,y,z)$ (red).

\noindent\textbf{Checks of the author's Secs.~5--7.} Every labelled equation was compared with the previous version and with our verified file (\texttt{s16a\_texcompare.py}): apart from the items above the only differences are rewritings already checked before. All colour factors of Sec.~7, as printed (TeX parsed), were compared with direct contraction of the verified colour tensors of Sec.~6 (\texttt{s16b\_author2\_sec7.py}, random $SU(3)$ and $SU(4)$ Wilson lines, coefficient and bracket together): 61 items, all equal (largest relative deviation $8\times10^{-15}$).

\noindent\textbf{New Sec.~\ref{sec:F}:} both forms, the $\rho$-explicit one (Sec.~6 form, Sec.~\ref{sec:F1}) and the projectile-averaged colour-singlet one (Sec.~7 form, Sec.~\ref{sec:F2}), Group by Group and Row by Row, with every colour structure written as sums of products of two differences of Wilson lines (times $f$, $d$ and one further Wilson line where unavoidable), all longitudinal integrals in $\xi$, checks in Sec.~\ref{sec:F3}.
}
'''
s = rep(s, r'\maketitle', r'\maketitle' + SUMMARY)
# Sec. 2: K vs mathcal K
s = rep(s, r'where $\mathcal{K}(x,y,z) = \frac{(x-z)\cdot (y-z)}{(x-z)^{2}\,(y-z)^{2}}$',
        r'where $' + RED(r'K(x,y,z)') + r' = \frac{(x-z)\cdot (y-z)}{(x-z)^{2}\,(y-z)^{2}}$ ' + RED(r'(\NOTE: $K$, as in Sec.~1; $\mathcal{K}$ is the vector $\mathcal{K}(x)=x/x^{2}$.)'))
# Conventions (Sec. 5)
CONV = r'''
{\color{red}
\subsection*{Conventions of the finite terms (re-inserted 2026-11-23)}\label{sec:conv2}
\begin{itemize}
\item \textbf{Phase} $e^{-ik\cdot(\text{unprimed}-\text{primed})}$ in all terms (amplitude coordinates unprimed, conjugate amplitude primed), with the recoil phases $e^{-ik\cdot(z-w)p^{+}/k^{+}}$ (Group I Row III), $e^{-ik\cdot(x-z)p^{+}/k^{+}}$ (Group I Row V), $e^{-ik\cdot(z-y')p^{+}/k^{+}}$ (Group III Rows III, V, VI), $e^{-ik\cdot(y-y')(1+p^{+}/k^{+})}$ (Group II Row IV). \textbf{Complex conjugate} implied in Groups I and III, not in Group II. The charges $\rho$ are operators in the order written. $\int_{x}\equiv\int d^{2}x$, $\mathcal{K}(x)\equiv\frac{x}{x^{2}}$, $K(x,y,z)=\mathcal{K}(x-z)\cdot\mathcal{K}(y-z)$; adjoint Wilson lines $U^{ab}$ are real, $U^{\dagger}=U^{T}$; $(T^{e})_{bc}=-if^{ebc}$, $(D^{e})_{bc}=d^{ebc}$, ${\rm Tr}\,T^{a}T^{b}=N_{c}\delta^{ab}$.
\item \textbf{Limits and their order.} $\Lambda\to0$ and $\vee\to\infty$ after the logarithmic (BFKL, JIMWLK) terms have been taken out; $\vee\to\infty$ before $\epsilon\to0$; terms $O(\Lambda/k^{+})$, $O(k^{+}/\vee)$ dropped. Upper limits of the $p^{+}$ integrals: $k^{+}$ in Group I, $\infty$ in Group III ($\vee-k^{+}$ with $\vee\to\infty$); Group II Row IV is written with its upper limit $\vee-k^{+}$.
\item \textbf{Prescription of the vanishing energy denominators} (those of $\mathbb{B}_{3}$, Group I Row III and its c.c.): \textbf{principal value only} (author, 2026-11-22); there is no $i\epsilon$ term $\Delta N$ in the result.
\item \textbf{Pure glue}: $N_{f}=0$, $\beta_{0}=\frac{11N_{c}}{3}$.
\item \textbf{Plus prescriptions} act on the whole $p^{+}$-dependent function $f$ of the term (recoil phase and energy denominators):
\begin{align}\label{eq:plus2}
\int_{0}^{k^{+}}dp^{+}\Big[\frac{1}{p^{+}}\Big]_{+}f(p^{+})&=\int_{0}^{k^{+}}\frac{dp^{+}}{p^{+}}\big[f(p^{+})-f(0)\big]\,,\qquad\Big[\frac{k^{+}}{p^{+}}\Big]_{+}=k^{+}\Big[\frac{1}{p^{+}}\Big]_{+}\,,\notag\\
\int_{0}^{k^{+}}dp^{+}\Big[\frac{k^{+}}{k^{+}-p^{+}}\Big]_{+}f(p^{+})&=\int_{0}^{k^{+}}dp^{+}\,\frac{k^{+}}{k^{+}-p^{+}}\big[f(p^{+})-f(k^{+})\big]\,,\notag\\
\int_{0}^{\infty}dp^{+}\Big[\frac{k^{+}}{p^{+}\,(p^{+}+k^{+})}\Big]_{+}f(p^{+})&=\int_{0}^{\infty}dp^{+}\,\frac{k^{+}}{p^{+}\,(p^{+}+k^{+})}\big[f(p^{+})-f(0)\big]\,.
\end{align}
\item \textbf{Distributions (instantaneous terms of Group III).} $\int_{0}^{\infty}dp^{+}\,e^{-ik\cdot(z-y')p^{+}/k^{+}}=k^{+}\big[\pi\delta(s)-i\,\mathcal{P}\frac{1}{s}\big]$, $s=k\cdot(z-y')$, defined by integration against the transverse kernels, with the two terms of the instantaneous bracket under one $y'$ integral, done before the $p^{+}$ integral.
\item \textbf{Terms that must be integrated together} (their large-distance tails, each of which alone makes the integral divergent, cancel only in these combinations; none of them is dropped): pair A: Group I Row III (term with $\big[U(w)-U(y)\big]$) with Group III Row V; pair B: Group I Row V first part ($\big[\frac{k^{+}}{p^{+}}\big]_{+}$ terms) with Group III Row III first part ($\big[\frac{k^{+}}{p^{+}(p^{+}+k^{+})}\big]_{+}$ term); pair C: Group I Row V second part ($\big[\frac{k^{+}}{p^{+}}\big]_{+}$, $\big[\frac{k^{+}}{k^{+}-p^{+}}\big]_{+}$ terms) with Group III Row III second part ($\big[\frac{k^{+}}{p^{+}(p^{+}+k^{+})}\big]_{+}$ term); pair D: Group I Row III (the $\int_{0}^{k^{+}}dp^{+}\big[\frac{1}{p^{+}}\big]_{+}e^{-ik\cdot(z-w)p^{+}/k^{+}}$ part of the term with $\big[U(z)-U(y)\big]$) with Group III Row VI. Each pair under one $z$ integral, the Group III member with $y'$ renamed $w'$ (pair D also $x\leftrightarrow y$). In Group I Row III the 24-term kernel and the part $-\frac{1}{2}\log\frac{(z-w)^{2}}{(x-w)^{2}}-\int_{0}^{k^{+}}dp^{+}\big[\frac{1}{p^{+}}\big]_{+}e^{-ik\cdot(z-w)p^{+}/k^{+}}$ are integrable only with the $w$ integral done before the $z$ integral. Group II Row IV, Group II Row I (second part) and Group III Row III (second part) have to be integrated together over the position of the gluon pair (joint tail, weights $1:1:-2$).
\item \textbf{Collinear point of Group I Row V, second part:} at $z\to x$ this term is finite only together with its $\frac{11}{6}$ line and with a rotation-invariant (angular-average-first) treatment of the point $z=x$ (our UV check of 2026-11-22).
\end{itemize}
}
'''
s = rep(s, r'and then we take $\epsilon \to 0$. ', r'and then we take $\epsilon \to 0$. ' + CONV)
# Sec. 5, Row I.II: the antisymmetrised kernel
s = rep(s, r'''\label{eq:B-I2b}\\''', r'''\label{eq:B-I2b}\\''')
s = rep(s, r'''\log\left[\frac{\mathbf{q}^{2}}{(\mathbf{q}-\mathbf{p})^{2}}\right]
\end{align}
\subsubsection*{Row III}''', r'''\log\left[\frac{\mathbf{q}^{2}}{(\mathbf{q}-\mathbf{p})^{2}}\right]
\end{align}
\fixnote{The kernel of $N_{2b}$ is the antisymmetrised one, $\frac{1}{2}\big[\frac{R^{i}(\mathbf{q},\mathbf{p})}{\mathbf{p}^{2}}-\frac{R^{i}(\mathbf{q},\mathbf{q}-\mathbf{p})}{(\mathbf{q}-\mathbf{p})^{2}}\big]$ (correct, as required): each of the two parts alone behaves as $\log(1/r)/r^{2}$ at $\mathbf{p}\to\mathbf{q}$ ($r=|\mathbf{p}-\mathbf{q}|$) and at $\mathbf{p}\to0$ ($r=|\mathbf{p}|$), which is not integrable over $d^{2}p$; in the antisymmetrised combination these terms cancel and the kernel is $O(\log(1/r)/r)$, integrable (\texttt{scripts/s15\_N2b\_p\_to\_q.py}). The same kernel is used in Secs.~6--8.}
\subsubsection*{Row III}''')
# Sec. 5, Row I.V second part: xi
s = rep(s, r'''\rho^{b'}(x')\,\rho^{f}(y)\,,\label{eq:B-I5two}
\end{align}''', r'''\rho^{b'}(x')\,\rho^{f}(y)\,,\label{eq:B-I5two}
\end{align}
\fixnote{Here $\xi\equiv\frac{p^{+}}{k^{+}}$ (undefined before). With it the equation is identical to the $p^{+}$ form of Eq.~\eqref{eq:D-I.V-2o} (previous check s11, (3): $4\times10^{-16}$).}''')
# Sec. 6: eq:D-ord third contact term
s = rep(s, r'''&+\frac{i}{2}C^{abc}f^{acr}\delta^{(2)}(p_{1}-p_{3})''', r'''&''' + RED(r'+f^{acr}') + r'''\delta^{(2)}(p_{1}-p_{3})''')
s = rep(s, r'''It stands for any $c$-number color tensor multiplying an ordered product of three Wilson lines. ''',
        r'''It stands for any $c$-number color tensor multiplying an ordered product of three ''' + RED(r'colour charges $\rho$ (the Wilson lines are inside $C^{abc}$)') + r'''.
\fixnote{(i) In the third contact term of Eq.~\eqref{eq:D-ord} the factor $\frac{i}{2}C^{abc}$ appeared a second time inside the bracket; removed (the identity as corrected is checked at the operator level, \texttt{s16e}: $4\times10^{-16}$; all rows below use the correct form). (ii) ``three Wilson lines'' $\to$ ``three colour charges $\rho$''. (iii) \textbf{One-$\rho$ terms:} the last lines of Eqs.~\eqref{eq:D-ord}--\eqref{eq:D-r21} (double commutators, one charge) are not written in the rows below. They vanish after the average over a colour-invariant projectile ($\langle\rho^{s}\rangle=0$), so Sec.~7 is complete; as operator identities the rows below equal those of Sec.~5 only together with them. (iv) \textbf{What remains:} Group I Row IV and Group III Row I cancel exactly, and the $\frac{11}{3}$ terms of Group I Rows II and V cancel exactly; all other terms of Sec.~5 remain. The terms whose large-distance tails cancel only in combination (pairs A--D, the Row I.III kernel, the joint tail) remain in their rows and must be integrated together as stated in Sec.~\ref{sec:conv2}. (v) In Eqs.~\eqref{eq:D-II.I-3}, \eqref{eq:D-II.I-2}, \eqref{eq:D-III.II-3}, \eqref{eq:D-III.II-2} the colour tensor is a difference of two terms; the red braces enclose it (without them the coefficient multiplies only the first term).}''')
# Sec. 6: missing outer brackets (red braces)
LB, RB = RED(r'\Big\{'), RED(r'\Big\}')
II1a = (r"U^{ga}(w')\notag \\ & \Big[f^{kbc}\,U^{gk}(w)-f^{gkl}\,U^{kb}(z)\,U^{lc}(x)\Big]-U^{ga}(x')\,U^{mh}(z)\,U^{mb}(y)\Big[f^{khc}\,U^{gk}(w)-f^{gkl}\,U^{kh}(z)\,U^{lc}(x)\Big]",
        r"U^{ga}(w')\Big[f^{kbc}\,U^{gk}(w)-f^{gkl}\,U^{kb}(z)\,U^{lc}(x)\Big]-U^{ga}(x')\,U^{mh}(z)\,U^{mb}(y)\Big[f^{khc}\,U^{gk}(w)-f^{gkl}\,U^{kh}(z)\,U^{lc}(x)\Big]")
II1b = (r"U^{gb}(w)\notag \\ & \Big[f^{kca}\,U^{gk}(w')-f^{gkl}\,U^{kc}(z)\,U^{la}(x)\Big]-U^{gb}(x')\,U^{mh}(z)\,U^{mc}(y)\Big[f^{kha}\,U^{gk}(w')-f^{gkl}\,U^{kh}(z)\,U^{la}(x)\Big]",
        r"U^{gb}(w)\Big[f^{kca}\,U^{gk}(w')-f^{gkl}\,U^{kc}(z)\,U^{la}(x)\Big]-U^{gb}(x')\,U^{mh}(z)\,U^{mc}(y)\Big[f^{kha}\,U^{gk}(w')-f^{gkl}\,U^{kh}(z)\,U^{la}(x)\Big]")
III2 = (r"U^{ma}(y)\,U^{gb}(x')\notag \\ & U^{mh}(z)\,\Big[f^{khc}\,U^{gk}(w)-f^{gkl}\,U^{kh}(z)\,U^{lc}(x)\Big]-U^{gb}(x')\Big[f^{kac}\,U^{gk}(w)-f^{gkl}\,U^{ka}(z)\,U^{lc}(x)\Big]",
        r"U^{ma}(y)\,U^{gb}(x')\,U^{mh}(z)\Big[f^{khc}\,U^{gk}(w)-f^{gkl}\,U^{kh}(z)\,U^{lc}(x)\Big]-U^{gb}(x')\Big[f^{kac}\,U^{gk}(w)-f^{gkl}\,U^{ka}(z)\,U^{lc}(x)\Big]")
AUTHOR_C = {'II.I 1a': II1a, 'II.I 1b': II1b, 'III.II': III2}
sec6a, sec6b = s.index(r'\section{Finite terms: After cancellation'), s.index(r'\section{Average over the projectile}')
S6 = s[sec6a:sec6b]
for k, (c3, c2) in AUTHOR_C.items():
    S6 = rep(S6, c3, LB + c3 + RB)
    S6 = rep(S6, c2, LB + c2 + RB)
s = s[:sec6a] + S6 + s[sec6b:]
# Sec. 7: averaging formula, sigma, o -> sigma, '+' between contact groups
sec7a = s.index(r'\section{Average over the projectile}'); S7 = s[sec7a:s.index(r'\end{document}')]
S7o = S7
S7 = rep(S7, r'''(N_{c}^{2}-1)}\,d^{efg}\{\rho^{e}(x)\rho^{f}(y)\rho^{g}(z)\}_{\rm sym}\big\rangle''', r'''(N_{c}^{2}-1)}\,''' + RED(r'\big\langle') + r'''d^{efg}\{\rho^{e}(x)\rho^{f}(y)\rho^{g}(z)\}_{\rm sym}\big\rangle''')
S7 = rep(S7, r'''$$\sigma(x,y,z) = d^{efg}''', r'''$$\sigma(x,y,z) = ''' + RED(r'\big\langle') + r'''d^{efg}''')
S7 = rep(S7, r'''\subsection{Group I}''', r'''\fixnote{(i) The opening $\langle$ was missing in the second formula and in the definition of $\sigma$. (ii) $n(x,y)=n(y,x)$ is the two-point function at $x\neq y$ (the contact parts are in the two-$\rho$ rows); $\sigma$ is real and totally symmetric; $\langle\rho^{a}(x)\rangle=0$; the $f^{abc}$ part of the three-point function vanishes; for $SU(2)$, $d^{abc}=0$ and the $\sigma$ terms are absent. (iii) With these, for the colour tensor $C^{abc}$ of a row (one-$\rho$ terms average to zero):
$\langle C^{abc}\rho^{a}(p_{1})\rho^{b}(p_{2})\rho^{c}(p_{3})\rangle=\frac{N_{c}}{(N_{c}^{2}-4)(N_{c}^{2}-1)}C^{abc}d^{abc}\,\sigma(p_{1},p_{2},p_{3})+\frac{i}{2(N_{c}^{2}-1)}\big[\Phi_{12}\delta^{(2)}(p_{1}-p_{2})n(p_{1},p_{3})+\Phi_{23}\delta^{(2)}(p_{2}-p_{3})n(p_{1},p_{2})-\Phi_{13}\delta^{(2)}(p_{1}-p_{3})n(p_{1},p_{2})\big]$, $\Phi_{ij}=C^{abc}f^{abc}$ with the Wilson lines at $p_{i}=p_{j}$; for the anticommutator types the factor $2$ on the first term and $i[\Phi_{12}\delta_{12}n-\Phi_{13}\delta_{13}n]$, resp.\ $i[\Phi_{23}\delta_{23}n-\Phi_{13}\delta_{13}n]$. (iv) The equations below wrote $o(\dots)$ for $\sigma(\dots)$; replaced (red). (v) In the contact brackets a red $+$ is inserted between the terms. (vi) All colour factors of this section were checked numerically (summary at the top).}
\subsection{Group I}''')
S7 = rep(S7, r'\,o(', r'\,' + RED(r'\sigma') + '(', 13)
S7 = rrep(S7, r'\\notag\\\\\n& \\Big\((?=[+-](?:i|\\frac\{1\}\{2\}i)\\Big\)\\delta)', lambda m: '\\notag\\\\\n& ' + RED('+') + '\\Big(', 20)
s = s[:sec7a] + S7 + s[sec7a + len(S7o):]
# ======================================================================================== (B) new Sec. 8
def xi_groupI(t):
    t = re.sub(r'\\cdot\{(\([^{}]*\))\}', r'\\cdot\1', t)
    t = t.replace(r'\int_{0}^{k^{+}}dp^{+}\Big[\frac{1}{p^{+}}\Big]_{+}e^{-ik\cdot(z-w)\frac{p^{+}}{k^{+}}}', r'\int_{0}^{1}d\xi\Big[\frac{1}{\xi}\Big]_{+}e^{-ik\cdot(z-w)\,\xi}')
    t = t.replace(r'\frac{1}{k^{+\,2}}\int_{0}^{k^{+}}\frac{dp^{+}}{2\pi}', r'\frac{1}{k^{+}}\int_{0}^{1}\frac{d\xi}{2\pi}')
    t = t.replace(r'\frac{k^{+}\,e^{-ik\cdot(x-z)\frac{p^{+}}{k^{+}}}}{p^{+}\,(y-x)^{2}+(k^{+}-p^{+})\,(y-z)^{2}}', r'\frac{e^{-ik\cdot(x-z)\,\xi}}{\xi\,(y-x)^{2}+(1-\xi)\,(y-z)^{2}}')
    t = t.replace(r'e^{-ik\cdot(x-z)\frac{p^{+}}{k^{+}}}', r'e^{-ik\cdot(x-z)\,\xi}').replace(r'e^{-ik\cdot(z-x)\frac{p^{+}}{k^{+}}}', r'e^{-ik\cdot(z-x)\,\xi}')
    t = t.replace(r'\Big[\frac{k^{+}}{p^{+}}\Big]_{+}', r'\Big[\frac{1}{\xi}\Big]_{+}').replace(r'\Big[\frac{k^{+}}{k^{+}-p^{+}}\Big]_{+}', r'\Big[\frac{1}{1-\xi}\Big]_{+}')
    t = t.replace(r'\frac{p^{+}(k^{+}-p^{+})}{k^{+\,2}}', r'\xi(1-\xi)').replace(r'+\frac{k^{+}-p^{+}}{k^{+}}\,\mathcal{K}', r'+(1-\xi)\,\mathcal{K}').replace(r'+\frac{p^{+}}{k^{+}}\,\mathcal{K}', r'+\xi\,\mathcal{K}')
    return t
def xi_groupII4(t):
    t = re.sub(r'\\cdot\{(\([^{}]*\))\}', r'\\cdot\1', t)
    t = rep(t, r'\int_{0}^{\vee-k^{+}}\frac{dp^{+}}{2\pi}\,\frac{p^{+}}{(p^{+}+k^{+})^{2}}', r'\int_{k^{+}/\vee}^{1}\frac{d\xi}{2\pi}\,\frac{1-\xi}{\xi}')
    t = rep(t, r"e^{-ik\cdot(y-y')}\,e^{-ik\cdot(y-y')\frac{p^{+}}{k^{+}}}", r"e^{-ik\cdot(y-y')\frac{1}{\xi}}")
    return t
def xi_groupIII(t):
    t = re.sub(r'\\cdot\{(\([^{}]*\))\}', r'\\cdot\1', t)
    # Row III second part, pieces 1-3 (specific forms first)
    t = t.replace(r'\frac{1}{(2\pi)^{3}}\,\frac{g^{4}}{4\pi^{4}}\int_{0}^{\infty}\frac{dp^{+}}{2\pi}', r'\frac{1}{(2\pi)^{3}}\,\frac{g^{4}}{4\pi^{4}}\,\frac{1}{k^{+}}\int_{0}^{\infty}\frac{d\xi}{2\pi}')
    t = t.replace(r'\frac{1}{k^{+\,2}}\int_{0}^{\infty}\frac{dp^{+}}{2\pi}', r'\frac{1}{k^{+}}\int_{0}^{\infty}\frac{d\xi}{2\pi}')
    t = t.replace(r'\frac{p^{+}}{p^{+}\,(x-z)^{2}+k^{+}\,(x-w)^{2}}', r'\frac{\xi}{\xi\,(x-z)^{2}+(x-w)^{2}}')
    t = t.replace(r'\frac{1}{p^{+}\,(x-z)^{2}+k^{+}\,(x-w)^{2}}', r'\frac{1}{\xi\,(x-z)^{2}+(x-w)^{2}}')
    t = t.replace(r'\frac{p^{+}\big[(x-z)^{2}-(x-w)^{2}\big]}{(z-w)^{2}\,(p^{+}+k^{+})^{2}}', r'\frac{\xi\big[(x-z)^{2}-(x-w)^{2}\big]}{(z-w)^{2}\,(1+\xi)^{2}}')
    t = t.replace(r'\frac{p^{+}}{k^{+}(p^{+}+k^{+})}', r'\frac{\xi}{1+\xi}')
    t = t.replace(r'\frac{(x-z)^{2}-(x-w)^{2}}{(p^{+}+k^{+})\,(x-w)^{2}}', r'\frac{(x-z)^{2}-(x-w)^{2}}{(1+\xi)\,(x-w)^{2}}')
    t = t.replace(r'-\frac{1}{k^{+}}\Big[(x-z)', r'-\Big[(x-z)')
    # generic Group III
    t = t.replace(r'\frac{1}{k^{+}}\int_{0}^{\infty}\frac{dp^{+}}{2\pi}', r'\frac{1}{k^{+}}\int_{0}^{\infty}\frac{d\xi}{2\pi}')
    t = t.replace(r"e^{-ik\cdot(z-y')\frac{p^{+}}{k^{+}}}", r"e^{-ik\cdot(z-y')\,\xi}")
    t = t.replace(r'\frac{1}{p^{+}+k^{+}}', r'\frac{1}{1+\xi}')
    t = t.replace(r'\Big[\frac{k^{+}}{p^{+}\,(p^{+}+k^{+})}\Big]_{+}', r'\Big[\frac{1}{\xi\,(1+\xi)}\Big]_{+}')
    t = t.replace(r'-\frac{1}{k^{+}}\,\mathcal{K}', r'-\mathcal{K}').replace(r'+\frac{1}{k^{+}}\,\mathcal{K}', r'+\mathcal{K}')
    return t
def wrapC(t, lab):
    """line breaks in the long factorised tensors"""
    t = t.replace(r"\bigg]+U^{ga}(x')\bigg[", r"\bigg]\notag\\&\qquad+U^{ga}(x')\bigg[").replace(r"\bigg]+U^{gb}(x')\bigg[", r"\bigg]\notag\\&\qquad+U^{gb}(x')\bigg[")
    t = t.replace(r"f^{khc}+U^{gb}(x')", r"f^{khc}\notag\\&\qquad+U^{gb}(x')")
    t = t.replace(r"f^{khc}-f^{gml}", r"f^{khc}\notag\\&\qquad\quad-f^{gml}").replace(r"f^{kha}-f^{gml}", r"f^{kha}\notag\\&\qquad\quad-f^{gml}")
    t = t.replace(r"\Big[U^{kb}(w)-U^{kb}(y)\Big]\bigg]", r"\Big[U^{kb}(w)-U^{kb}(y)\Big]\bigg]\notag\\&\qquad")
    t = t.replace(r"\bigg]\bigg[\Big[U^{ac}(w)", r"\bigg]\notag\\&\qquad\times\bigg[\Big[U^{ac}(w)")
    return t
def C8(lab): return wrapC(CUBIC8[lab]['tex'], lab)
def wrap_terms(tex, per=2, first=1):
    parts = re.split(r'(?<=\\Big\\\})(?=[+-])', tex)
    o = ''
    for n, p in enumerate(parts):
        if n and (n == first or (n > first and (n - first) % per == 0)): o += r'\notag\\&\qquad'
        o += p
    return o
# ---------- 8.1: transformed copy of Sec. 6
s6 = s[s.index(r'\section{Finite terms: After cancellation'):s.index(r'\section{Average over the projectile}')]
F1 = s6[s6.index(r'\subsection{Group I}'):]
F1 = F1.replace(LB, '').replace(RB, '')   # the red braces of the fixes are not needed (the tensor is replaced)
F1 = re.sub(r'\\fixnote\{.*?\}\n', '', F1, flags=re.S)
def Bg(t): return r'\Bigg\{' + t + r'\Bigg\}'
# Group I
F1 = rep(F1, r"2\,\Big[U^{ga}(x')-U^{ga}(w')\Big]\times \Big[f^{hbc}\,U^{gh}(w)-f^{gde}\,U^{db}(x)\,U^{ec}(y)\Big]", r"2\," + Bg(C8('I.II 2b')))
F1 = rep(F1, r"i\,\Big[U^{ga}(x')-U^{ga}(w')\Big]\Big[f^{hbc}\,U^{gh}(w)-f^{gde}\,U^{db}(x)\,U^{ec}(y)\Big]", r"i\," + Bg(C8('I.II 2b')))
for u, lab in (('z', 'I.III C_D (u=z)'), ('w', 'I.III C_A (u=w)')):
    F1 = rep(F1, r"\Big[U^{ga}(x')-U^{ga}(w')\Big]f^{ghd}\,U^{db}(x)\Big[U^{hc}(%s)-U^{hc}(y)\Big]" % u, Bg(C8(lab)), 2)
F1 = rep(F1, r"2\,\Big[U^{ga}(x')-U^{ga}(w')\Big]f^{gde}\Big[U^{dc}(x)\,U^{eb}(z)-U^{eb}(y)\,U^{dc}(y')\Big]", r"2\," + Bg(C8('I.V first')))
F1 = rep(F1, r"\frac{ig^{4}\,f^{gde}}{16\pi^{4}}", r"\frac{ig^{4}}{16\pi^{4}}")
F1 = rep(F1, r"i\,\Big[U^{ga}(x')-U^{ga}(w')\Big]\Big[U^{dc}(x)\,U^{eb}(z)-U^{eb}(y)\,U^{dc}(y')\Big]", r"i\," + Bg(C8('I.V first')))
F1 = rep(F1, r"\Big[U^{ab'}(x')-U^{ab'}(w')\Big]\Big[-f^{adc}\,f^{fbe}\,U^{db}(x)\,U^{ce}(z)+N_{c}\,U^{af}(y)\Big]", Bg(QUAD8['I.V N2']['tex']))
# Group II
F1 = rep(F1, r"2\,U^{ga}(w')\notag \\ & \Big[f^{kbc}\,U^{gk}(w)-f^{gkl}\,U^{kb}(z)\,U^{lc}(x)\Big]-U^{ga}(x')\,U^{mh}(z)\,U^{mb}(y)\Big[f^{khc}\,U^{gk}(w)-f^{gkl}\,U^{kh}(z)\,U^{lc}(x)\Big]", r"2\," + Bg(C8('II.I 1a')))
F1 = rep(F1, r"2\,U^{gb}(w)\notag \\ & \Big[f^{kca}\,U^{gk}(w')-f^{gkl}\,U^{kc}(z)\,U^{la}(x)\Big]-U^{gb}(x')\,U^{mh}(z)\,U^{mc}(y)\Big[f^{kha}\,U^{gk}(w')-f^{gkl}\,U^{kh}(z)\,U^{la}(x)\Big]", r"2\," + Bg(C8('II.I 1b')))
F1 = rep(F1, r"i\," + II1a[1], r"i\," + Bg(C8('II.I 1a')))
F1 = rep(F1, r"i\," + II1b[1], r"i\," + Bg(C8('II.I 1b')))
QII = (r"\Bigg[f^{c^{\prime}bd^{\prime}}\,f^{cbd}\, U^{ac^{\prime}}(w^{\prime})\,U^{ac}(w)\,\rho^{d^{\prime}}(x^{\prime})\rho^{d}(x)- f^{c^{\prime}bd^{\prime}}\,f^{acd}\,U^{ac^{\prime}}(w^{\prime})\,U^{cb}(z)\,U^{de}(x)\,\rho^{d^{\prime}}(x^{\prime})\rho^{e}(x) \notag \\ & -f^{ac^{\prime}d^{\prime}}\,f^{cbd}\,U^{ac}(w)\,U^{c^{\prime}b}(z)\,U^{d^{\prime}e^{\prime}}(x^{\prime})\,\rho^{e^{\prime}}(x^{\prime})\,\rho^{d}(x) +N_{c}\,U^{de^{\prime}}(x^{\prime})\,U^{de}(x)\,\rho^{e^{\prime}}(x^{\prime})\,\rho^{e}(x)\Bigg]")
F1 = rep(F1, QII, wrapC(QUAD8['II.I second']['tex'], '') + r"\rho^{d'}(x')\,\rho^{d}(x)", 5)
# Group III
F1 = rep(F1, r"\frac{i}{2}\," + III2[1], r"\frac{i}{2}\," + Bg(C8('III.II')))
F1 = rep(F1, III2[0], Bg(C8('III.II')))
for v, tag in (("y'", "(y')"), ('z', '(z)')):
    F1 = rep(F1, r"f^{mnk}\Big[U^{ma}(%s)-U^{ma}(x')\Big]\Big[U^{nc}(z)\,U^{kb}(w)-U^{kb}(y)\,U^{nc}(x)\Big]" % v, Bg(C8('III.III L' + tag)), 2)
    F1 = rep(F1, r"\Big[U^{ma}(%s)-U^{ma}(x')\Big]U^{kb}(x)\Big[U^{lc}(y)-U^{lc}(w)\Big]" % v, Bg(C8('III.V L' + tag)), 2)
    F1 = rep(F1, r"\Big[U^{ma}(%s)-U^{ma}(x')\Big]U^{kb}(y)\Big[U^{nc}(x)-U^{nc}(z)\Big]" % v, Bg(C8('III.VI L' + tag)), 2)
F1 = rep(F1, r"\frac{ig^{4}\,f^{mkl}}{4\pi^{4}}", r"\frac{ig^{4}}{4\pi^{4}}", 2)
F1 = rep(F1, r"\frac{ig^{4}\,f^{mnk}}{4\pi^{4}}", r"\frac{ig^{4}}{4\pi^{4}}", 2)
# Group III Row III second part (genuine two-rho): L(v) and R factors
LFAC = lambda v: r"\,f^{c'd'a}\,U^{d'b}(z)\Big[U^{c'e'}(%s)-U^{c'e'}(x')\Big]\rho^{e'}(x')" % v
RFAC = r"\bigg[\Big[U^{ac}(w)-U^{ac}(z)\Big]f^{cbd}-f^{ace}\,U^{cb}(z)\Big[U^{ed}(x)-U^{ed}(z)\Big]\bigg]\rho^{d}(x)"
F1 = rep(F1, r"\Big[f^{c'd'a}\,U^{c'e'}(y')\,U^{d'b}(z)\,\rho^{e'}(x')-f^{c'd'a}\,U^{d'b}(z)\,U^{c'e'}(x')\,\rho^{e'}(x')\Big]", LFAC("y'"), 3)
F1 = rep(F1, r"\Big[f^{c'd'a}\,U^{c'e'}(z)\,U^{d'b}(z)\,\rho^{e'}(x')-f^{c'd'a}\,U^{d'b}(z)\,U^{c'e'}(x')\,\rho^{e'}(x')\Big]", LFAC('z'), 1)
F1 = rep(F1, r"\left[f^{c'd'a}\,U^{c'e'}(z)\,U^{d'b}(z)\,\rho^{e'}(x')-f^{c'd'a}\, U^{d'b}(z)\,U^{c'e'}(x')\,\rho^{e'}(x')\right]", LFAC('z'), 1)
F1 = rep(F1, r"\Big[f^{cbd}\,U^{ac}(w)\,\rho^{d}(x)-f^{acd}\,U^{cb}(z)\,U^{de}(x)\,\rho^{e}(x)\Big]", RFAC, 3)
F1 = rep(F1, r"\left[f^{cbd}\,U^{ac}(w)\,\rho^{d}(x) -f^{acd}\,U^{cb}(z)\,U^{de}(x)\,\rho^{e}(x)\right]", RFAC, 1)
# xi in the kernels
g1, g2, g3 = F1.index(r'\subsection{Group I}'), F1.index(r'\subsection{Group II}'), F1.index(r'\subsection{Group III}')
F1 = xi_groupI(F1[g1:g2]) + xi_groupII4(F1[g2:g3]) + xi_groupIII(F1[g3:])
F1 = re.sub(r"(\\bigg\\\})(\\times\s*)?((?:2|i|\\frac\{i\}\{2\})\\,)?(\\Bigg\\\{)", lambda m: m.group(1) + r"\notag\\&\quad\times " + (m.group(3) or '') + m.group(4), F1)
F1 = re.sub(r"\\Bigg\\\}(\\times)?\\,\\Big\[f\^\{abr\}", lambda m: r"\Bigg\}\notag\\&\quad\times\Big[f^{abr}", F1)
F1 = F1.replace(r'\label{eq:D-', r'\label{eq:H-')
assert 'p^{+}' not in F1, F1[F1.index('p^{+}') - 300:F1.index('p^{+}') + 100]
# ---------- 8.2: transformed copy of Sec. 7
s7 = s[s.index(r'\section{Average over the projectile}'):s.index(r'\end{document}')]
F2 = s7[s7.index(r'\subsection{Group I}'):]
PLAN = [  # same order as in s16b
    ('eq:E-I.II-3', ['I.II 2b']), ('eq:E-I.II-2', ['I.II 2b'] * 2), ('eq:E-I.II-2o', ['I.II N1']),
    ('eq:E-I.III-3', ['I.III C_D (u=z)', 'I.III C_A (u=w)']), ('eq:E-I.III-2', ['I.III C_D (u=z)'] * 3 + ['I.III C_A (u=w)'] * 3),
    ('eq:E-I.V-3', ['I.V first']), ('eq:E-I.V-2', ['I.V first'] * 2), ('eq:E-I.V-2o', ['I.V N2', 'I.V 11/6']),
    ('eq:E-II.I-3', ['II.I 1a', 'II.I 1b']), ('eq:E-II.I-2', ['II.I 1a'] * 2 + ['II.I 1b'] * 2), ('II.I genuine', ['II.I second'] * 5),
    ('eq:E-II.IVrem', ['II.IV', 'II.IV rem']), ('eq:E-III.II-3', ['III.II']), ('eq:E-III.II-2', ['III.II'] * 3),
    ('eq:E-III.III-3', ["III.III L(y')", 'III.III L(z)']), ('eq:E-III.III-2', ["III.III L(y')"] * 2 + ['III.III L(z)'] * 2),
    ('eq:E-III.III-2o', ["III.III second L(y')"] * 3 + ['III.III second L(z)'] * 2),
    ('eq:E-III.V-3', ["III.V L(y')", 'III.V L(z)']), ('eq:E-III.V-2', ["III.V L(y')"] * 3 + ['III.V L(z)'] * 3),
    ('eq:E-III.VI-3', ["III.VI L(y')", 'III.VI L(z)']), ('eq:E-III.VI-2', ["III.VI L(y')"] * 3 + ['III.VI L(z)'] * 3)]
PT = {"x'": 'xp', "w'": 'wp', "y'": 'yp', 'x': 'x', 'y': 'y', 'z': 'z', 'w': 'w'}
PTS = {'I.II 2b': ('xp', 'x', 'y'), 'I.III C_D (u=z)': ('xp', 'x', 'y'), 'I.III C_A (u=w)': ('xp', 'x', 'y'), 'I.V first': ('xp', 'y', 'yp'),
       'II.I 1a': ('xp', 'y', 'x'), 'II.I 1b': ('x', 'xp', 'y'), 'III.II': ('y', 'xp', 'x')}
for v in ("(y')", '(z)'):
    for r_ in ('III.III L', 'III.V L', 'III.VI L'): PTS[r_ + v] = {'III.III L': ('xp', 'y', 'x'), 'III.V L': ('xp', 'x', 'y'), 'III.VI L': ('xp', 'y', 'x')}[r_]
pat = re.compile(r"\\Big\((?P<c3>[^()]*?)\\Big\)\\,\{\\color\{red\}\\sigma\}\((?P<oa>[^()]*)\)"
                 r"|\\Big\((?P<cc>[^()]*?)\\Big\)\\delta\^\{\(2\)\}\((?P<pq>[^()]*)\)\\,n\((?P<na>[^()]*)\)"
                 r"|\\frac\{n\((?P<qa>[^()]*)\)\}\{N_\{c\}\^\{2\}-1\}")
def block_span(T, lab):
    if lab == 'II.I genuine':
        a = T.index(r'\label{eq:E-II.I-2}'); return a, T.index(r'\subsubsection*{Row II}', a)
    i = T.index(r'\label{%s}' % lab); return T.rfind(r'\begin{align}', 0, i), i
nrep = 0
for lab, structs in reversed(PLAN):
    a, b = block_span(F2, lab); blk = F2[a:b]
    ms = list(pat.finditer(blk)); assert len(ms) == len(structs), (lab, len(ms))
    edits = []
    for m, st in zip(ms, structs):
        j = m.end()
        while blk[j] in ' \\,' and not blk.startswith(r'\Big', j) and not blk.startswith(r'{\rm', j): j += 1
        bare = blk.startswith(r'{\rm Tr}', j)
        if bare:
            b0 = j + len(r'{\rm Tr}'); cont, end = bracket_after(blk, b0); start = j
        else:
            cont, end = bracket_after(blk, j); start = j
        if m.group('c3') is not None:
            d = SING['cubic'][st]['3rho']
            new_head = r'\Big(' + d['coef'] + r'\Big)\,\sigma(' + m.group('oa') + ')'
            edits.append((m.start(), end, new_head + r'\Big[' + wrap_terms(d['tex']) + r'\Big]'))
        elif m.group('cc') is not None:
            p, q = [PT[x.strip()] for x in m.group('pq').split('-')]
            pts = PTS[st]; pair = ''.join(sorted('%d%d' % (pts.index(p) + 1, pts.index(q) + 1)))
            d = SING['cubic'][st]['contacts'][pair]
            edits.append((start, end, r'\Big[' + wrap_terms(d['tex']) + r'\Big]'))
        else:
            Q = QUAD8[st]['sing']
            edits.append((start, end, Q if bare else r'\Big[' + wrap_terms(Q) + r'\Big]'))
        nrep += 1
    for st_, en_, new in reversed(edits): blk = blk[:st_] + new + blk[en_:]
    F2 = F2[:a] + blk + F2[b:]
assert nrep == 61, nrep
assert r'\color{red}\sigma' not in F2 and 'o(' not in F2.replace(r'\log(', '')
g1, g2, g3 = F2.index(r'\subsection{Group I}'), F2.index(r'\subsection{Group II}'), F2.index(r'\subsection{Group III}')
F2 = xi_groupI(F2[g1:g2]) + xi_groupII4(F2[g2:g3]) + xi_groupIII(F2[g3:])
F2 = F2.replace(r'\label{eq:E-', r'\label{eq:S-').replace(RED('+'), '+')
F2 = re.sub(r'\\fixnote\{.*?\}\n', '', F2, flags=re.S)
assert 'p^{+}' not in F2, F2[F2.index('p^{+}') - 300:F2.index('p^{+}') + 100]
# headings: subsection -> subsubsection etc. inside Sec. 8
def demote(T, tag):
    T = T.replace(r'\subsubsection*{', r'\paragraph*{')
    T = T.replace(r'\subsection{Group I}', r'\subsubsection{Group I%s}' % tag).replace(r'\subsection{Group II}', r'\subsubsection{Group II%s}' % tag).replace(r'\subsection{Group III}', r'\subsubsection{Group III%s}' % tag)
    return T
F1 = demote(F1, '')
F2 = demote(F2, '')
INTRO = open(os.path.join(HERE, 'sec8_intro.tex')).read()
CHECKS = open(os.path.join(HERE, 'sec8_checks.tex')).read()
SEC8 = (r'\section{Final result with factorised Wilson-line structures}\label{sec:F}' + '\n' + INTRO
        + '\n' + r'\subsection{$\rho$-explicit form (after cancellation and symmetrisation)}\label{sec:F1}' + '\n'
        + r'\noindent Prefactors, transverse kernels and phases as in Sec.~6; the colour tensors are written with differences of Wilson lines and the longitudinal integrals with $\xi$. Complex conjugate implied in Groups I and III (written as $+{\rm c.c.}$).' + '\n'
        + F1 + '\n' + r'\subsection{Average over the projectile: colour-singlet form}\label{sec:F2}' + '\n'
        + r'\noindent The same rows averaged over the projectile (Sec.~7 form; $n(x,y)$, $\sigma(x,y,z)$ as defined in Sec.~7), the traces written with differences of Wilson lines, the longitudinal integrals with $\xi$.' + '\n'
        + F2 + '\n' + CHECKS + '\n')
s = s.replace(r'\end{document}', SEC8 + r'\end{document}')
open(os.path.join(TOP, 'Final_terms_explicit_2_fixed.tex'), 'w').write(s)
print('written; Sec. 8: %d chars, %d singlet items replaced' % (len(SEC8), nrep))
