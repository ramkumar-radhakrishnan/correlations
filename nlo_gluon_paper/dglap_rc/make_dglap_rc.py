r"""Finite_terms_DGLAP_RC.tex: from the remaining terms of ../bfkl_split/Finite_terms_split.tex
  Sec. 2  running coupling (beta_0 terms, coupling counterterm, RG term, constant 67/9 - pi^2/3)
  Sec. 3  DGLAP (all poles, the DGLAP finite term, the pole that is left over)
  Sec. 4  the finite terms, row by row, in the format of Finite_terms_1.tex
  Sec. 5  the result
  Sec. 6  comparison with AABKL (arXiv:1411.2869)
Bookkeeping checked in checks.py (sympy)."""
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.normpath(os.path.join(HERE, '..', 'bfkl_split')))
sys.path.insert(0, os.path.normpath(os.path.join(HERE, '..', 'bfkl_check')))
import make_split as MS  # noqa: E402
import make_corrected as MC  # noqa: E402

OUT = os.path.join(HERE, 'Finite_terms_DGLAP_RC.tex')
LV = r"\log{\left[\frac{\vee}{k^{+}}\right]}"
L1A = r"\log{\left[1-\frac{k^{+}}{\vee}\right]}"
COL_II = (r"\left[U^{ab'}(x')\,\rho^{b'}(x') - U^{ab'}(w')\,\rho^{b'}(x')\right]\times\,\left[U^{ab}(w)\,\rho^{b}(x) - "
          r"U^{ab}(x)\,\rho^{b}(x)\right]")
COL_UV = r"\Big[U^{ab'}(x')-U^{ab'}(w')\Big]\Big[U^{af}(y)-U^{af}(z)\Big]\rho^{b'}(x')\,\rho^{f}(y)"
K_II = r"\frac{(x'-w')\cdot (x-w)}{(x'-w')^{2}\, (x-w)^{2}}"

INTRO = r"""
\section{What is done here}
\begin{itemize}
\item Starting point: the remaining terms of \texttt{Finite\_terms\_split.pdf} (Secs.~7--9), i.e.\ every row of \texttt{Finite\_terms\_1.tex} without its BFKL ($\log\frac{\vee}{k^+}$) terms, with $\left[\frac{1}{p^+}\right]_+$ kept as $\left[\frac{1}{p^+}\right]_+^{(k^+)}$.
\item Only three of them contain $\frac1\epsilon$ poles or $\beta_0$: Group I Row II part 1, Group I Row V 2b UV, Group II Row IV. All other rows are finite.
\item From these three rows, and the coupling counterterm of LO, the running coupling (Sec.~\ref{sec:rc}) and the DGLAP term (Sec.~\ref{sec:dglap}) are taken out. What is left, with all other rows, are the finite terms (Sec.~\ref{sec:fin}). Sec.~\ref{sec:res} is the result, Sec.~\ref{sec:aabkl} the comparison with the paper.
\end{itemize}
\paragraph{Notation.} $a\equiv\frac{k^+}{\vee}$, $\beta_0=\frac{11}{3}N_c$ ($n_f=0$), $d_\perp=2-\epsilon$, $c_0\equiv2e^{-\gamma_E}$, $\ell_r\equiv\log\frac{r^2\mu^2}{c_0^2}=2\gamma_E+\log\frac{r^2\mu^2}{4}$, $\hat E_r\equiv\frac2\epsilon+\ell_r$; $X=x-y$, $X'=x'-y'$, $\Delta=y'-y$. LO in your notation, and LO with a weight $f$:
\begin{align*}
{\rm LO}[f]&\equiv\frac{1}{(2\pi)^{3}}\frac{g^{2}}{2\pi^{2}}\frac{1}{k^{+}}\int_{x,y,x',y'}e^{-ik\cdot(y'-y)}\,f\,\frac{(x-y)\cdot(x'-y')}{(x-y)^{2}(x'-y')^{2}}\left[U^{e'c}(y')-U^{e'c}(x')\right]\left[U^{ec}(y)-U^{ce}(x)\right]\rho^{e'}(x')\rho^{e}(x),
\end{align*}
${\rm LO}\equiv{\rm LO}[1]$, and ${\rm LO}_\xi[f]$ = the same with $e^{-ik\cdot(y'-y)}\to e^{-ik\cdot(y'-y)\frac1\xi}$. $P^{(r)}_{gg,+}(\xi)=2N_c\big[\frac{\xi}{(1-\xi)_+}+\frac{1-\xi}{\xi}+\xi(1-\xi)\big]$ with the plus prescription on $[a,1]$, $P_{gg}=P^{(r)}_{gg,+}+\frac{\beta_0}{2}\delta(1-\xi)$. The conversion of the three rows to units of LO was checked against your file before (\texttt{../finite\_v2/final\_checks.json}); the bookkeeping below is checked in \texttt{checks.py}.
"""

RC = r"""
\section{Running coupling}\label{sec:rc}
\paragraph{The $\beta_0$ terms.} Two rows have $\beta_0$: the $\frac{11}{3}$ of Group I Row II part 1 and the $-\frac{11}{6}$ of Group I Row V 2b UV:
\begin{align}
\frac{d^{3}N^{\beta_0}_{1}}{d^{3}k}&=\frac{1}{(2\pi)^{3}}\,\frac{g^{4}\,N_{c}}{32\pi^{4}}\,\frac{1}{k^{+}}\,\frac{11}{3}\int_{w,w',x,x'}e^{-ik(w'-w)}\,@@KII@@\left[-\frac{2}{\epsilon}-2\gamma-\log{\left[\frac{(x-w)^{2}\,\mu^{2}_{\overline{MS}}}{4}\right]}\right]\notag\\
&\quad\times@@COLII@@\ \ +{\rm c.c.}\;=\;\frac{\alpha_s\beta_0}{4\pi}\,{\rm LO}\big[\hat E_X+\hat E_{X'}\big],\label{eq:b1}\\
\frac{d^{3}N^{\beta_0}_{2b\,\rm UV}}{d^{3}k}&=-\frac{1}{(2\pi)^{3}}\,\frac{g^{4}N_{c}}{16\pi^{4}}\,\frac{1}{k^{+}}\,\frac{11}{6}\int_{y,z,x',w'}e^{-ik\cdot(w'-z)}\left[\frac{2}{\epsilon}+2\gamma_{E}+\log\frac{(y-z)^{2}\mu^{2}_{\overline{MS}}}{4}\right]\frac{(x'-w')\cdot(y-z)}{(x'-w')^{2}\,(y-z)^{2}}\notag\\
&\quad\times@@COLUV@@\ \ +{\rm c.c.}\;=\;-\frac{\alpha_s\beta_0}{4\pi}\,{\rm LO}\big[\hat E_X+\hat E_{X'}\big].\label{eq:b2}
\end{align}
With $y\to x$, $z\to w$ in \eqref{eq:b2}, the measure, the phase, the kernel and the bracket are those of \eqref{eq:b1}, the colour factor changes sign, and $\frac{g^4N_c}{16\pi^4}\frac{11}{6}=\frac{g^4N_c}{32\pi^4}\frac{11}{3}$. So
\begin{equation}
\boxed{\frac{d^{3}N^{\beta_0}_{1}}{d^{3}k}+\frac{d^{3}N^{\beta_0}_{2b\,\rm UV}}{d^{3}k}=0}
\end{equation}
for every $a$: pole, $\gamma_E$, $\mu$ and $\log(x-w)^2$ all cancel (\texttt{checks.json}). The virtual correction (Row II) and the real emission at short distance (Row V 2b UV) carry opposite $\beta_0$ terms, and no $\beta_0\log(x-w)^2$, i.e.\ no running of the coupling at the emission distance, is left.

\paragraph{The coupling of LO.} LO is $\propto g^2$. In $\overline{MS}$, $g^2\to g^2(\mu_R)\big[1-\frac{\alpha_s\beta_0}{4\pi}\frac{2}{\epsilon}\big]$ ($d_\perp=2-\epsilon$), which adds
\begin{equation}
\frac{d^3N}{d^3k}\bigg|_{\rm CT}=-\frac{\alpha_s\beta_0}{4\pi}\,\frac{2}{\epsilon}\,{\rm LO}=\frac{\alpha_s}{2\pi}\Big(-\frac{2}{\epsilon}\Big)\frac{\beta_0}{2}\,{\rm LO}.
\label{eq:ct}
\end{equation}
This pole is the $\frac{\beta_0}{2}\delta(1-\xi)$ part of the DGLAP pole (Sec.~\ref{sec:dglap}). After the coupling is renormalized at $\mu_R$ and the fragmentation function at $\mu_F$ (both $\overline{MS}$), the only $\beta_0$ term of the result is the renormalization-group term
\begin{equation}
\boxed{\frac{d^3N}{d^3k}\bigg|_{\beta_0}=\frac{\alpha_s\beta_0}{4\pi}\,\log\frac{\mu_R^2}{\mu_F^2}\;{\rm LO}\Big|_{\alpha_s(\mu_R)}}
\label{eq:rg}
\end{equation}
It makes the result independent of $\mu_R$: $\frac{d}{d\log\mu_R^2}{\rm LO}|_{\alpha_s(\mu_R)}=-\frac{\alpha_s\beta_0}{4\pi}{\rm LO}$. It vanishes for $\mu_R=\mu_F$; then LO is taken with $\alpha_s(\mu_F)$.

\paragraph{The constant.} The non-logarithmic part of the bracket of Row II part 1, $-\frac{67}{9}+\frac{\pi^2}{3}$, gives
\begin{align}
&\frac{1}{(2\pi)^{3}}\frac{g^{4}N_{c}}{32\pi^{4}}\frac{1}{k^{+}}\int_{w,w',x,x'}e^{-ik(w'-w)}@@KII@@\Big[-\frac{67}{9}+\frac{\pi^{2}}{3}\Big]\notag\\
&\qquad\times@@COLII@@+{\rm c.c.}\;=\;\boxed{\frac{\alpha_s}{2\pi}\,N_c\Big[\frac{67}{9}-\frac{\pi^2}{3}\Big]{\rm LO}}
\label{eq:cmw}
\end{align}
$N_c(\frac{67}{9}-\frac{\pi^2}{3})=2K$ with $K=N_c(\frac{67}{18}-\frac{\pi^2}{6})$, the constant of the two-loop cusp anomalous dimension that the CMW scheme absorbs into the coupling. It is the finite companion of $\beta_0$ in this row.
"""

DGLAP = r"""
\section{DGLAP}\label{sec:dglap}
\paragraph{All poles.} The poles of the remaining terms, in units of LO (\texttt{checks.json}):
\begin{center}\small
\begin{tabular}{@{}ll@{}}\hline
G1 Row II part 1 $+$ c.c. & $\frac{\alpha_sN_c}{4\pi}\frac{4}{\epsilon}\big[\frac{11}{3}+\log(1-a)\big]{\rm LO}$\\
G1 Row V 2b UV $+$ c.c. & $\frac{\alpha_sN_c}{2\pi}\frac{4}{\epsilon}\big[-\frac{11}{6}+2\log(1-\frac{\Lambda}{k^+})\big]{\rm LO}$ \quad($\log(1-\frac{\Lambda}{k^+})\to0$)\\
G2 Row IV, lines 1, 2 & $\frac{\alpha_s}{2\pi}\big(-\frac{2}{\epsilon}\big)\int_a^1\frac{d\xi}{\xi^2}P^{(r)}_{gg,+}(\xi)\,{\rm LO}_\xi$\\
coupling counterterm \eqref{eq:ct} & $\frac{\alpha_s}{2\pi}\big(-\frac{2}{\epsilon}\big)\frac{\beta_0}{2}{\rm LO}$\\\hline
sum & $\frac{\alpha_s}{2\pi}\big(-\frac{2}{\epsilon}\big)\Big[\int_a^1\frac{d\xi}{\xi^2}P_{gg}(\xi)\,{\rm LO}_\xi-N_c\log(1-a)\,{\rm LO}\Big]$\\\hline
\end{tabular}
\end{center}
The first term is the DGLAP pole: it is absorbed in $D_{h/g}(\zeta,\mu_F)$ ($\overline{MS}$), with the convolution $\int_{k^+/\vee}^1\frac{d\zeta}{\zeta^2}D_{h/g}(\zeta,\mu_F)\frac{d^3N_g}{d^3k}(\frac{k^+}{\zeta},\frac{\mathbf k}{\zeta})$. The second is not (end of this section).

\paragraph{The DGLAP term.} Group II Row IV, lines 1 and 2 of the remaining terms, with the pole removed (your bracket $\big[-\frac{d_\perp}{\epsilon}-d_\perp\gamma_E-\frac{d_\perp}{2}\log\frac{(y'-y)^2\mu^2}{4}-1\big]=-\hat E_\Delta+O(\epsilon)$ becomes $-\ell_\Delta$; the $-1$ cancels against $-\frac{d_\perp}{\epsilon}$), together with the $\log\frac{\vee}{k^+}H_4$ subtraction:
\begin{align}
\frac{d^{3}N_{\rm DGLAP}}{d^{3}k}&=\frac{g^{2}}{4\pi^{2}}\,\frac{N_{c}}{(2\pi)^{3}}\,\frac{g^{2}}{2\pi^{2}}\,\frac{1}{k^{+}}\int_{x,y,x',y'}\frac{(x-y)\cdot(x'-y')}{(x-y)^{2}\,(x'-y')^{2}}\left[-2\gamma_{E}-\log\left(\frac{(y'-y)^{2}\,\mu_{F}^{2}}{4}\right)\right]\notag\\
&\quad\times\Bigg[\int_{\frac{k^{+}}{\vee}}^{1}\frac{d\xi}{\xi^{2}}\,e^{-ik\cdot(y'-y)\frac{1}{\xi}}\left[\xi(1-\xi)+\frac{1-\xi}{\xi}\right]+\int_{\frac{k^{+}}{\vee}}^{1}d\xi\,\frac{\frac{1}{\xi}e^{-ik\cdot(y'-y)\frac{1}{\xi}}-e^{-ik\cdot(y'-y)}}{1-\xi}\notag\\
&\qquad+\log\left(1-\frac{k^{+}}{\vee}\right)e^{-ik\cdot(y'-y)}\Bigg]\notag\\
&\quad\times\left[U^{e'c}(y')-U^{e'c}(x')\right]\left[U^{ec}(y)-U^{ce}(x)\right]\rho^{e'}(x')\,\rho^{e}(x)\notag\\
&\quad-\frac{1}{(2\pi)^{3}}\,\frac{g^{4}N_{c}}{8\pi^{5}}\,\frac{1}{k^{+}}\,@@LV@@\int_{x,x',z,w,w'}e^{-ik\cdot(w'-w)}\,\frac{(w-z)\cdot(w'-z)}{(w-z)^{2}(w'-z)^{2}}\,\frac{(x-z)\cdot(x'-z)}{(x-z)^{2}(x'-z)^{2}}\notag\\
&\quad\times\left[U^{e'c}(z)-U^{e'c}(x')\right]\left[U^{ec}(z)-U^{ce}(x)\right]\rho^{e'}(x')\,\rho^{e}(x)\notag\\
&=\boxed{-\frac{\alpha_s}{2\pi}\int_{k^+/\vee}^{1}\frac{d\xi}{\xi^{2}}\,P^{(r)}_{gg,+}(\xi)\,{\rm LO}_\xi\big[\ell_\Delta(\mu_F)\big]-\log\frac{\vee}{k^+}\,H_4,}\qquad \ell_\Delta(\mu_F)=\log\frac{(y'-y)^2\mu_F^2}{c_0^2}.
\label{eq:dglap}
\end{align}
\begin{itemize}
\item The bracket in $\xi$ is the plus distribution $\frac{1}{2N_c}P^{(r)}_{gg,+}$ on $[\frac{k^+}{\vee},1]$; its last term is your $\log(1-\frac{k^+}{\vee})$ line.
\item $\mu_F^2\frac{\partial}{\partial\mu_F^2}D_{h/g}=\frac{\alpha_s}{2\pi}P_{gg}\otimes D_{h/g}$ cancels the $\mu_F$ dependence of $\ell_\Delta(\mu_F)$ and of \eqref{eq:rg}.
\item The last two lines ($-\log\frac{\vee}{k^+}H_4$) remove the $\log\frac{\vee}{k^+}$ that the $\xi$ integral produces at $\xi\to\frac{k^+}{\vee}$ (it is in the BFKL terms); with them the DGLAP term is finite for $\vee\to\infty$.
\item Line 3 of the row (linear polarization) is finite and not DGLAP; it is in Sec.~\ref{sec:fin}.
\end{itemize}

\paragraph{The pole that is left (as before, Q1).} The second term of the pole, with its finite partner from Row II part 1:
\begin{equation}
\frac{d^3N}{d^3k}\bigg|_{\rm Q1}=\frac{\alpha_sN_c}{4\pi}\log\Big(1-\frac{k^+}{\vee}\Big)\,{\rm LO}\big[\hat E_X+\hat E_{X'}\big].
\label{eq:q1}
\end{equation}
It is not DGLAP (it is not $P_{gg}\otimes$LO) and not running coupling. It comes from the limit $\xi\ge\frac{k^+}{\vee}$ in the virtual term of Row II part 1 and vanishes for $k^+\ll\vee$, i.e.\ at mid rapidity in the high-energy limit. In AABKL the virtual term has $\int_0^1d\xi$ and there is no such term (Sec.~\ref{sec:aabkl}).
"""

RESULT = r"""
\section{Result}\label{sec:res}
\begin{equation}
\boxed{\begin{aligned}
\frac{d^3N_h}{d^3k}={}&\int_{k^+/\vee}^1\frac{d\zeta}{\zeta^2}\,D_{h/g}(\zeta,\mu_F)\,\frac{d^3N_g}{d^3k}\Big(\frac{k^+}{\zeta},\frac{\mathbf k}{\zeta}\Big),\\
\frac{d^3N_g}{d^3k}={}&{\rm LO}\big|_{\alpha_s(\mu_R)}+\frac{\alpha_s\beta_0}{4\pi}\log\frac{\mu_R^2}{\mu_F^2}\,{\rm LO}+\frac{\alpha_s}{2\pi}N_c\Big[\frac{67}{9}-\frac{\pi^2}{3}\Big]{\rm LO}\qquad\text{(running coupling, Sec.~\ref{sec:rc})}\\
&+\frac{1}{(2\pi)^{3}}\frac{g^{4}}{16\pi^{5}}\frac{1}{k^{+}}\int_{x,x',z,w,w'}e^{-ik(w'-w)}\Big[\log\frac{k^+}{\Lambda}\,\mathcal J+\log\frac{\vee}{k^+}\,\mathcal B\Big]\qquad\text{(evolution: JIMWLK, BFKL)}\\
&+\frac{d^3N_{\rm DGLAP}}{d^3k}\qquad\text{(Eq.~\eqref{eq:dglap})}\\
&+\frac{\alpha_sN_c}{4\pi}\log\Big(1-\frac{k^+}{\vee}\Big)\,{\rm LO}\big[\hat E_X+\hat E_{X'}\big]\qquad\text{(Q1, Eq.~\eqref{eq:q1})}\\
&+\text{finite terms (Sec.~\ref{sec:fin})}+\text{four-}\rho\text{ term (Q2)}+O(\alpha_s^3).
\end{aligned}}
\end{equation}
\begin{itemize}
\item $\mathcal J$, $\mathcal B$: the brackets of Secs.~3 and 4 of \texttt{Evolution\_terms.tex}; $\log\frac{\vee}{\Lambda}\mathcal J$ (\texttt{Evolution\_terms.tex}) plus the BFKL terms $\log\frac{\vee}{k^+}(\mathcal B-\mathcal J)$ (\texttt{Finite\_terms\_split.pdf}).
\item The finite terms contain three-$\rho$ rows. Their fully symmetric part needs $\langle\rho\rho\rho\rangle$ (zero for a Gaussian projectile); their two-$\rho$ parts are written out in \texttt{Finite\_terms\_simplified.pdf}, Sec.~9.
\item The four-$\rho$ finite term is not in \texttt{Finite\_terms\_1.tex} (\texttt{Finite\_terms\_simplified.pdf}, Sec.~10).
\end{itemize}
"""

AABKL = r"""
\section{Comparison with AABKL (arXiv:1411.2869)}\label{sec:aabkl}
T.~Altinoluk, N.~Armesto, G.~Beuf, A.~Kovner, M.~Lublinsky: forward single-inclusive production in the hybrid formalism. The projectile is a collinear parton with PDF $f(x_p)$, LO is $O(\alpha_s^0)$, NLO is one splitting. The relevant channel is $g\to g$, their Eqs.~(4.20), (4.22) and App.~C. Dictionary: $\xi_{\rm AABKL}=1-\xi$ (they use the fraction of the unobserved gluon), $x_p\leftrightarrow\frac{k^+}{\vee}$, $A^i(y-z)=-\frac{1}{2\pi}\frac{(y-z)^i}{(y-z)^2}\leftrightarrow-\frac{1}{2\pi}\mathcal K(y-z)$, $s_A[y,\bar y]\leftrightarrow$ the LO colour factor, $(y,\bar y)\leftrightarrow(y,y')$.

\paragraph{1. DGLAP term: the same.}
\begin{itemize}
\item Splitting function: their $\frac{1-\xi}{[\xi]_+}+\frac{\xi}{1-\xi}+\xi(1-\xi)$ at $\xi\to1-\xi$ is $\frac{1}{2N_c}P^{(r)}_{gg,+}$, with the plus prescription at the soft end in both (\texttt{checks.json}: difference $0$).
\item Kinematic limit: their real term has $\int_0^{1-x_p}d\xi$, yours $\int_{k^+/\vee}^1d\xi$; both produce $\log(1-x_p)$, resp.\ $\log(1-\frac{k^+}{\vee})$, through the plus distribution.
\item The logarithm: their after-target term $2C_A\big[A^i(y-z)A^i(\bar y-z)-(A^i(z)A^i(z))_{\mu^2}\big](1-\xi)^2s_A[(1-\xi)y,(1-\xi)\bar y]$ integrates over $z$ to $\log\big((y-\bar y)^2\mu^2\big)$ times the rescaled dipole; yours is $\ell_\Delta=\log\frac{(y'-y)^2\mu_F^2}{c_0^2}$ times ${\rm LO}_\xi$ (the rescaling is your $\frac1\xi$ in the phase). Their App.~C, Eq.~(7.13): $-\frac{g^2}{2(2\pi)^2}\log\big(\mu^2_{CO}(y-\bar y)^2\big)P\otimes\dots$, the same structure as Eq.~\eqref{eq:dglap}.
\item Scheme: their Eq.~(7.22), $R=2e^{\Psi(1)}\simeq1.1229$ relates their scheme to $\overline{MS}$. This is exactly your $c_0=2e^{-\gamma_E}$ (\texttt{checks.json}: $R-c_0=0$), so their $\log(r^2\mu^2)$ and your $\ell_r=\log\frac{r^2\mu^2}{c_0^2}$ are the same $\overline{MS}$ logarithm.
\item Hadron convolution: $\int_{x_F}^1\frac{d\zeta}{\zeta^2}D(\zeta)\frac{d\bar\sigma}{d^2p}(\frac{p}{\zeta},\frac{x_F}{\zeta})$ vs.\ $\int_{k^+/\vee}^1\frac{d\zeta}{\zeta^2}D_{h/g}(\zeta)\frac{d^3N_g}{d^3k}(\frac{k^+}{\zeta},\frac{\mathbf k}{\zeta})$: the same.
\end{itemize}

\paragraph{2. Running coupling: different, because LO is different.}
\begin{itemize}
\item AABKL: LO is $O(\alpha_s^0)$, so there is no coupling to renormalize and no $\beta_0\log\mu_R^2$ term. $\beta_0$ enters only as the $\frac{\beta_0}{2}\delta(1-\xi)$ part of $P_{gg}$ in the DGLAP evolution of $f^g$ and $D^g$. It comes from their virtual term: $\int_0^1d\xi\big[\frac{1-\xi}{[\xi]_+}+\frac{\xi}{[1-\xi]_+}+\xi(1-\xi)\big]=-\frac{11}{6}$, and $2C_A\cdot\frac{11}{6}=\beta_0$ (\texttt{checks.json}).
\item Yours: LO is $O(\alpha_s)$. The $\beta_0$ of the virtual term (Row II) cancels against the short-distance real emission (Row V 2b UV), and the $\frac{\beta_0}{2}\delta(1-\xi)$ of $P_{gg}$ comes from the coupling counterterm of LO. In the end the same $P_{gg}$ is absorbed in $D_{h/g}$.
\item The RG term \eqref{eq:rg} and the constant $\frac{\alpha_s}{2\pi}N_c(\frac{67}{9}-\frac{\pi^2}{3}){\rm LO}$ have no counterpart in AABKL.
\end{itemize}

\paragraph{3. The finite terms: same colour structures and topologies, different kernels.}
{\small\renewcommand{\arraystretch}{1.25}
\begin{longtable}{@{}>{\raggedright\arraybackslash}p{0.34\textwidth}>{\raggedright\arraybackslash}p{0.27\textwidth}>{\raggedright\arraybackslash}p{0.33\textwidth}@{}}\hline
AABKL, Eq.~(4.22) & yours, Sec.~\ref{sec:fin} & comparison\\\hline\endhead
before the target: $2C_A[A\!\cdot\!A-(AA)_{\mu^2}]\,s_A[y,\bar y]$ & G2 Row I part 2, colour $f\!f\,U(w')U(w)$ & same topology and colour; theirs collinear divergent (PDF subtraction), yours finite: the parent is emitted at $|x-w|$ and the energy denominators are $p^+(x-z)^2+k^+(x-w)^2$\\
interference: $\frac{N_c^3}{N_c^2-1}A\!\cdot\!A\{sss+\dots-\frac{1}{N_c^2}\overline X\}$ & G3 Row III part 2; mixed colour terms of G2 Row I part 2 & same topology (parent at $(1-\xi)y+\xi z$ vs.\ your $y'$), colour the same after projection, kernels differ\\
virtual: $-x_pf\int_0^1d\xi(\dots)\{C_A[AA+AA-2(AA)_{\mu^2}]s_A+\dots\}$ & G1 Row V 2a, 2b finite: colour $-f\!f\,U(x)U(z)+N_cU(x)$ & colour the same; theirs $\int_0^1d\xi$, yours $\int_0^{k^+}dp^+$; your Row II has $\xi\ge\frac{k^+}{\vee}$ (Q1, $V_\vee$)\\
last term: $\int\frac{d\xi}{\xi}[A_\xi-A_\xi]^2\{\dots\}$ (Ioffe time) & evolution, $\log\frac{k^+}{\Lambda}\mathcal J+\log\frac{\vee}{k^+}\mathcal B$ & same soft real$-$virtual (BK/JIMWLK) structure; theirs regulated by the Ioffe time, yours by the $p^+$ limits\\
--- & linear polarization (G2 Row IV, line 3) & absent: their incoming gluon is unpolarized\\
--- & parent emitted after the target by the rotated charge (G2 Row I part 2, G3 Rows V, VI) & absent: their parent exists before the target\\
--- & G1 Rows III, IV, G2 Row I part 1, G3 Rows I, II (emission by the charges around the splitting) & absent: no colour charges in their projectile\\
--- & three-$\rho$, four-$\rho$ terms & absent: one parton\\
--- & $\log\frac{\vee}{k^+}\mathcal B$ (BFKL) & absent: their projectile evolution sits in the PDF\\
--- & $V_\vee$, Q1 & absent: their virtual term has no $x_p$ limit\\\hline
\end{longtable}}

\paragraph{In short.}
\begin{itemize}
\item \textbf{DGLAP}: the same as AABKL, in the splitting function, the plus prescription, the kinematic limit, the logarithm and the scheme constant ($R=c_0$).
\item \textbf{Running coupling}: none in AABKL (LO is $\alpha_s^0$). Yours has the RG term and the constant $\frac{67}{9}-\frac{\pi^2}{3}$; the $\beta_0$ of $P_{gg}$ is the same in both.
\item \textbf{Finite terms}: they resemble AABKL in colour structure and topology (before the target, interference, virtual), not in the kernels. The difference is the projectile: colour charges that emit the gluon at LO, so your kernels contain the emission distance and there is no collinear divergence to subtract into a PDF. Linear polarization, emission after the target, the terms with emission by the charges around the splitting, three- and four-$\rho$ terms, $V_\vee$ and Q1 have no counterpart.
\end{itemize}
"""


def finite_sections():
    t = MS.remaining_sections()
    t = t.replace(r'\section{Remaining terms: ', r'\section{Finite terms: ')
    aligns = list(re.finditer(r'\\begin\{align\}.*?\\end\{align\}', t, flags=re.S))

    def find(marker):
        m = [x for x in aligns if marker in x.group(0)]
        assert len(m) == 1, (marker, len(m))
        return m[0]

    # Group I Row II part 1: only V_vee is left
    a = find(r"\frac{11}{3}+\log{\left[1-\frac{k^{+}}{\vee}\right]}")
    new_ii = (r"\begin{align}" + "\n" + r"\frac{d^{3}N_{1\,\rm NLO}^{\rm}}{d^{3}k} & = \frac{1}{(2\pi)^{3}}\, \frac{g^{4}\,N_{c}}{32\pi^{4}}\,\frac{1}{k^{+}}\,"
              r"\int_{w, w',x,x'} e^{-ik(w^\prime-w)}\, " + K_II + r"\Bigg[\log^{2}{\left[1-\frac{k^{+}}{\vee}\right]}+" + L1A + r"\," + LV
              + r"\Bigg]\notag \\ & \times" + COL_II + "\n" + r"\end{align}" + "\n"
              + r"\new{Part 1: what is left is $V_\vee=-\frac{\alpha_s}{2\pi}N_c\log(1-\frac{k^+}{\vee})\log(\frac{\vee}{k^+}-1)\,{\rm LO}$ (with c.c.); it vanishes for $k^+\ll\vee$. The $\frac{11}{3}$ is in Sec.~\ref{sec:rc} (cancels), the constant $-\frac{67}{9}+\frac{\pi^2}{3}$ in Eq.~\eqref{eq:cmw}, the pole bracket times $\log(1-\frac{k^+}{\vee})$ in Eq.~\eqref{eq:q1}.}" + "\n")
    t = t[:a.start()] + new_ii + t[a.end():]
    old_note = r"\new{Part 1: the remaining $\log\big[1-\frac{k^+}{\vee}\big]\log\big[\frac{\vee}{k^+}\big]$ vanishes for $\vee\to\infty$; it is not a $\log\frac{\vee}{k^+}$ term. Part 2b: no change.}"
    assert t.count(old_note) == 1
    t = t.replace(old_note, r"\new{Part 2b: no change.}")
    aligns = list(re.finditer(r'\\begin\{align\}.*?\\end\{align\}', t, flags=re.S))
    # Group I Row V 2b UV: nothing is left
    a = find(r"2b\,\rm UV")
    t = t[:a.start()] + (r"\new{$d^3N_{2b\,\rm UV}$: nothing is left. Its $-\frac{11}{6}$ is in Sec.~\ref{sec:rc} (cancels against Row II part 1); "
                         r"its $2\log(1-\frac{\Lambda}{k^+})$ term vanishes for $\Lambda\to0$.}" + "\n") + t[a.end():]
    aligns = list(re.finditer(r'\\begin\{align\}.*?\\end\{align\}', t, flags=re.S))
    # Group II Row IV: only line 3 (linear polarization)
    a = find(r"\frac{\xi}{(1-\xi)_{+}}")
    L = a.group(0)[len(r'\begin{align}'):a.group(0).rindex(r'\end{align}')].split(r'\notag \\')
    lin = [l for l in L if r"\int_{\frac{k^{+}}{\vee}}^{\frac{k^{+}}{k^{+}+\Lambda}}" in l]
    assert len(lin) == 1
    k = L.index(lin[0])
    head = r"\frac{d^{3}N_{\rm NLO}^{\rm}}{d^{3}k}  & = " + re.sub(r'^\+\s*\\,\s*', '', lin[0].split('&', 1)[1].strip())
    new_iv = r"\begin{align}" + "\n" + head + r"\notag \\" + L[k + 1].rstrip() + "\n" + r"\end{align}" + "\n"
    note_iv = (r"\new{\textbf{Row IV.} Only line 3 (linear polarization) is left here; it is finite. Lines 1, 2 and the $\log\frac{\vee}{k^+}H_4$ subtraction "
               r"are the DGLAP term, Eq.~\eqref{eq:dglap}; their pole is in Sec.~\ref{sec:dglap}.}" + "\n")
    rest = t[a.end():]
    assert rest.lstrip().startswith(MS.G2IV_NOTE.strip()[:40]), rest[:200]
    i = rest.index(MS.G2IV_NOTE.strip()) + len(MS.G2IV_NOTE.strip())
    t = t[:a.start()] + new_iv + note_iv + rest[i:]
    for a_, b_ in ((r"definition in Sec.~6", r"definition in \texttt{Finite\_terms\_split.pdf}, Sec.~6"),
                   (r"BFKL term of this row (Sec.~4)", r"BFKL term of this row (\texttt{Finite\_terms\_split.pdf}, Sec.~4)"),
                   (r"all of it is in Sec.~2", r"all of it is in the BFKL terms, \texttt{Finite\_terms\_split.pdf}, Sec.~2")):
        t = t.replace(a_, b_)
    return t


def assemble():
    src = open(MS.SRC).read()
    pre = src[:src.index(r'\begin{document}')]
    pre = pre.replace(r"\title{Finite terms}", r"\title{Finite terms: DGLAP, running coupling and the finite terms}")
    pre = pre.replace(r"\usepackage{comment}", r"\usepackage{comment}" + "\n" + r"\usepackage{longtable}\usepackage{array}" + "\n"
                      + r"\newcommand{\new}[1]{{\color{blue}#1}}" + "\n" + r"\newcommand{\old}[1]{{\color{red}#1}}" + "\n"
                      + r"\newcommand{\nc}{\color{blue}}" + "\n" + r"\newcommand{\rc}{\color{red}}" + "\n" + r"\allowdisplaybreaks")
    rc = RC.replace('@@KII@@', K_II).replace('@@COLII@@', COL_II).replace('@@COLUV@@', COL_UV)
    dg = DGLAP.replace('@@LV@@', LV)
    fin = finite_sections()
    fin_intro = r"""
\section{The finite terms}\label{sec:fin}
The remaining terms of \texttt{Finite\_terms\_split.pdf}, Secs.~7--9, without the running coupling and the DGLAP term. Group I Row II part 1, Row V 2b UV and Group II Row IV are reduced as described under each; every other row is as in \texttt{Finite\_terms\_split.pdf} (blue: changed with respect to \texttt{Finite\_terms\_1.tex}). All of them are finite for $\epsilon\to0$, $\vee\to\infty$ and $\Lambda\to0$; the terms with $\log(1\pm\frac{\Lambda}{k^+})$ vanish for $\Lambda\to0$.
"""
    doc = (pre + r"\begin{document}" + "\n" + r"\maketitle" + "\n" + r"\sloppy" + "\n" + INTRO + rc + dg + fin_intro + fin + RESULT
           + AABKL + r"\end{document}" + "\n")
    open(OUT, 'w').write(doc)
    print('written', OUT)


if __name__ == '__main__':
    assemble()
