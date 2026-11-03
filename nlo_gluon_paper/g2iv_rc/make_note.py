r"""Group II Row IV (input/G2IV.tex, your equation): DGLAP, running coupling and finite terms, step by step.
  python3 make_note.py && pdflatex G2IV_DGLAP_RC.tex (twice)"""
import os
import re
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
blk = open(os.path.join(HERE, 'input', 'G2IV.tex')).read()
body = blk[blk.index(r'\begin{align}') + len(r'\begin{align}'):blk.index(r'\end{align}')]
lines = [x.strip() for x in body.split('\\notag \\\\\n')]
assert len(lines) == 6, len(lines)
L1, L2, L3a, L3b, L4, L5 = lines
L1 = L1.replace(r'\frac{d^{3}N_{\rm NLO}^{\rm}}{d^{3}k}  & =', r'(L1)\quad & ')
assert L1.startswith('(L1)')

# ---- the bookkeeping, symbolically, with d_perp = 2 - eps kept (h = d_perp/2); units alpha_s/2pi
eps, lmu, Pr, b0, LO, h = sp.symbols('epsilon ell_Delta_mu P_r beta0 LO h')
mF, mR, lDF = sp.symbols('m_F m_R ellD_F')            # m_X = log(mu^2/mu_X^2), ellD_F = ell_Delta(mu_F)
Pgg = Pr + b0 / 2 * LO                                 # P_gg x LO = P^(r) x LO + beta0/2 LO
row = Pr * (-2 / eps - h * lmu)                        # (L1)+(L4): your bracket is exactly -2/eps - (d_perp/2) ell_Delta(mu)
ct_F = Pgg * (2 / eps + h * mF)                        # fragmentation counterterm at mu_F
ct_R = -(b0 / 2) * (2 / eps + h * mR) * LO             # coupling counterterm at mu_R
tot = sp.expand((row + ct_F + ct_R).subs(lmu, lDF + mF))
claim = sp.expand(-Pgg * h * lDF + (b0 / 2) * LO * h * (lDF + mF - mR))   # ellD_R = ellD_F + m_F - m_R
assert sp.simplify(tot - claim) == 0
# exact bracket: -d/eps - d*gE - (d/2) log - 1 with d = 2 - eps  ==  -2/eps - (d/2)(2 gE + log)
gE, Lg = sp.symbols('gamma_E L')
d = 2 - eps
assert sp.simplify((-d / eps - d * gE - d / 2 * Lg - 1) - (-2 / eps - d / 2 * (2 * gE + Lg))) == 0
CHECK = 'exact bracket identity and pole/mu bookkeeping with d_perp = 2 - eps verified (sympy)'

T = r"""\documentclass[10pt]{article}
\usepackage[a4paper,margin=0.35in]{geometry}
\usepackage{amsmath,amssymb,amsfonts,mathtools}
\usepackage[dvipsnames]{xcolor}
\usepackage{booktabs}
\usepackage[colorlinks=true,linkcolor=blue]{hyperref}
\usepackage{etoolbox}
\AtBeginEnvironment{align}{\small}
\allowdisplaybreaks
\numberwithin{equation}{section}
\setlength{\emergencystretch}{3em}
\newcommand{\new}[1]{{\color{blue}#1}}
\newcommand{\old}[1]{{\color{red}#1}}
\newcommand{\LO}{{\rm LO}}
\title{Group II Row IV: DGLAP, running coupling and finite terms, step by step}
\date{}
\begin{document}
\maketitle
\noindent Your equation (Group II Row IV, no c.c.). Conventions:
\begin{itemize}
\item $a\equiv\frac{k^+}{\vee}$, $\alpha_s=\frac{g^2}{4\pi}$, $\beta_0=\frac{11}{3}N_c$ ($n_f=0$), $c_0=2e^{-\gamma_E}$, $\Delta\equiv y'-y$.
\item $\ell_\Delta(\mu)\equiv2\gamma_E+\log\frac{(y'-y)^2\mu^2}{4}$.
\item \textbf{$d_\perp=2-\epsilon$ everywhere; no step sets $d_\perp=2$.}
\end{itemize}
\new{Blue} = new, \old{red} = removed (with respect to the previous version of this note).

\section{Answer}
\begin{center}\fbox{\begin{minipage}{0.97\textwidth}
\begin{itemize}
\item \textbf{No $\log\frac{\vee}{k^+}$.} You are right: none is left. The previous version showed \old{$-\log\frac{\vee}{k^+}H_4$} because I had done the $y'$ integral of your subtraction line (L2) on its own.
\begin{itemize}
\item Your line (L1) contains $+\log\frac{\vee}{k^+}$ inside its $\xi$ integral, from $\xi\to\frac{k^+}{\vee}$.
\item Your line (L2) is exactly the subtraction that removes it.
\item Separately the two are $\pm\log\frac{\vee}{k^+}$; together they are finite for $\vee\to\infty$.
\end{itemize}
Here (L2) is kept as you wrote it, inside the DGLAP term, and no $\log\frac{\vee}{k^+}$ appears (Step 7). The $+\log\frac{\vee}{k^+}H_4$ of your BFKL list is the counterpart of this subtraction, so nothing is counted twice.
\item \textbf{DGLAP.} Your bracket is \emph{exactly} $-\frac2\epsilon-\frac{d_\perp}{2}\ell_\Delta(\mu)$: the $+1$ of $-\frac{d_\perp}{\epsilon}$ cancels your $-1$. The pole of the row is $\frac{\alpha_s}{2\pi}(-\frac2\epsilon)P^{(r)}_{gg,+}\otimes\LO$, with $P^{(r)}_{gg,+}=P_{gg}-\frac{\beta_0}{2}\delta(1-\xi)$.
\begin{itemize}
\item The $P_{gg}$ part goes into $D_{h/g}(\zeta,\mu_F)$, which obeys $\mu_F^2\partial_{\mu_F^2}D=\frac{\alpha_s}{2\pi}P_{gg}\otimes D$.
\end{itemize}
\item \textbf{Running coupling.} The $-\frac{\beta_0}{2}\delta(1-\xi)$ part, together with the $\overline{MS}$ counterterm of the coupling of LO, leaves
\begin{equation*}
\frac{\alpha_s\beta_0}{4\pi}\,\LO\Big[\frac{d_\perp}{2}\,\ell_\Delta(\mu_R)\Big]=\frac{\alpha_s\beta_0}{4\pi}\,\LO\Big[\frac{d_\perp}{2}\log\frac{(y'-y)^2\mu_R^2}{c_0^2}\Big].
\end{equation*}
This is the running of the coupling of LO with the transverse distance $|y'-y|$, conjugate to $\mathbf k$. At $\epsilon\to0$ it turns $\alpha_s(\mu_R)$ into $\alpha_s(c_0/|y'-y|)$.
\item \textbf{Finite terms.} Eq.~\eqref{eq:final2} gives all of them in your notation:
\begin{itemize}
\item the DGLAP term with your subtraction (L2);
\item your (L3);
\item the running-coupling term;
\item your (L5).
\end{itemize}
In practice: your row with $-\frac{d_\perp}{\epsilon}$ and $-1$ removed (they add up to exactly $-\frac{2}{\epsilon}$, which goes into $D_{h/g}$ and the coupling) and $\mu_{\overline{MS}}\to\mu_F$, plus the two \new{$\frac{11}{12}$} terms.
\end{itemize}
\textbf{Other $\beta_0$ terms.} Group I Row II part 1 ($\frac{11}{3}$) and Group I Row V $d^3N_{2b\,\rm UV}$ ($-\frac{11}{6}$) have $\beta_0$ terms of their own. These cancel each other exactly. Keep both or drop both: in your latest \texttt{Other\_terms.tex}, Row II part 1 is commented out but Row V $2b$ UV is not.
\end{minipage}}\end{center}

\section{Notation}
\begin{align}
\LO[f]&\equiv\frac{1}{(2\pi)^{3}}\frac{g^{2}}{2\pi^{2}}\frac{1}{k^{+}}\int_{x,y,x',y'}e^{-ik\cdot(y'-y)}\,f\,\frac{(x-y)\cdot(x'-y')}{(x-y)^{2}(x'-y')^{2}}\left[U^{e'c}(y')-U^{e'c}(x')\right]\left[U^{ec}(y)-U^{ce}(x)\right]\rho^{e'}(x')\rho^{e}(x),\label{eq:LO}
\end{align}
$\LO\equiv\LO[1]$. $\LO_\xi[f]$: the same with $e^{-ik\cdot(y'-y)}\to e^{-ik\cdot(y'-y)\frac1\xi}$, so $\LO_\xi=\frac1\xi\frac{d^3N_{\rm LO}}{d^3k}(\frac{k^+}{\xi},\frac{\mathbf k}{\xi})$ ($d^3k=dk^+d^2\mathbf k$). $\LO^{(y)}_\xi[f]$: the same as $\LO_\xi[f]$ with $\frac{(x'-y')}{(x'-y')^2}\to\frac{(x'-y)}{(x'-y)^2}$ and $U^{e'c}(y')\to U^{e'c}(y)$, i.e.\ the $y'=y$ value used in your line (L2).
\begin{equation}
P^{(r)}_{gg,+}(\xi)=2N_c\Big[\frac{\xi}{(1-\xi)_+}+\frac{1-\xi}{\xi}+\xi(1-\xi)\Big],\qquad P_{gg}(\xi)=P^{(r)}_{gg,+}(\xi)+\frac{\beta_0}{2}\,\delta(1-\xi),\label{eq:P}
\end{equation}
with the plus distribution of $[0,1]$: $\int_a^1d\xi\frac{g(\xi)}{(1-\xi)_+}=\int_a^1d\xi\frac{g(\xi)-g(1)}{1-\xi}+g(1)\log(1-a)$.

\section{Step by step}
\paragraph{Step 0: your row, line by line.}
\begin{align}
""" + L1 + r"""\notag\\
(L2)\quad""" + L2 + r"""\notag\\
(L3)\quad""" + L3a + r"""\notag\\
""" + L3b + r"""\notag\\
(L4)\quad""" + L4 + r"""\notag\\
(L5)\quad""" + L5 + r"""
\end{align}

\paragraph{Step 1: units of LO.} $\frac{g^2}{4\pi^2}=\frac{\alpha_s}{\pi}$. The $x,y,x',y'$ integrands of (L1), (L4) are those of \eqref{eq:LO}, (L1) with $e^{-ik\cdot\Delta}\to e^{-ik\cdot\Delta/\xi}$. With $\mathcal B_\Delta$ your bracket:
\begin{equation}
({\rm L1})=\frac{\alpha_sN_c}{\pi}\int_a^1\frac{d\xi}{\xi^2}\Big[\xi(1-\xi)+\frac{1-\xi}{\xi}+\xi+\frac{\xi^2}{(1-\xi)_+}\Big]_{\rm yours}\LO_\xi[\mathcal B_\Delta],\qquad
({\rm L4})=\frac{\alpha_sN_c}{\pi}\log(1-a)\,\LO[\mathcal B_\Delta].
\end{equation}

\paragraph{Step 2: the bracket, exactly.} With $d_\perp=2-\epsilon$, $-\frac{d_\perp}{\epsilon}=-\frac{2}{\epsilon}+1$. So, with no expansion in $\epsilon$:
\begin{equation}
\mathcal B_\Delta=-\frac{d_\perp}{\epsilon}-d_\perp\gamma_E-\frac{d_\perp}{2}\log\frac{\Delta^2\mu^2}{4}-1=-\frac{2}{\epsilon}-\frac{d_\perp}{2}\Big[2\gamma_E+\log\frac{\Delta^2\mu^2}{4}\Big]=-\frac2\epsilon-\frac{d_\perp}{2}\,\ell_\Delta(\mu).\label{eq:B}
\end{equation}
The $+1$ of $-\frac{d_\perp}{\epsilon}$ and your $-1$ cancel exactly (\texttt{make\_note.py}, sympy).

\paragraph{Step 3: (L4) completes the plus distribution.} Your prescription, $\int\frac{d\xi}{\xi^2}\frac{\xi^2}{(1-\xi)_+}G=\int\frac{d\xi}{1-\xi}[G(\xi)-G(1)]$, has no $\log(1-a)$. On $[a,1]$:
\begin{equation*}
\int_a^1\frac{d\xi}{\xi^2}\Big[\xi+\frac{\xi^2}{(1-\xi)_+}\Big]_{\rm yours}G=\int_a^1\frac{d\xi}{\xi^2}\frac{\xi}{(1-\xi)_+}G-\log(1-a)\,G(1),
\end{equation*}
with $\LO_{\xi=1}=\LO$ (\texttt{other\_check/checks.py}, item 3). (L4) adds the missing term:
\begin{equation}
({\rm L1})+({\rm L4})=\frac{\alpha_s}{2\pi}\int_a^1\frac{d\xi}{\xi^2}\,P^{(r)}_{gg,+}(\xi)\,\LO_\xi\Big[-\frac2\epsilon-\frac{d_\perp}{2}\ell_\Delta(\mu)\Big].\label{eq:L14}
\end{equation}

\paragraph{Step 4: split $P^{(r)}_{gg,+}$.} With $P^{(r)}_{gg,+}=P_{gg}-\frac{\beta_0}{2}\delta(1-\xi)$:
\begin{equation}
({\rm L1})+({\rm L4})=\underbrace{\frac{\alpha_s}{2\pi}\int_a^1\frac{d\xi}{\xi^2}\,P_{gg}(\xi)\,\LO_\xi\Big[-\frac2\epsilon-\frac{d_\perp}{2}\ell_\Delta(\mu)\Big]}_{\rm DGLAP}\;+\;\underbrace{\frac{\alpha_s\beta_0}{4\pi}\,\LO\Big[\frac2\epsilon+\frac{d_\perp}{2}\ell_\Delta(\mu)\Big]}_{\rm running\ coupling}.\label{eq:split}
\end{equation}

\paragraph{Step 5: DGLAP.} $\int_a^1\frac{d\xi}{\xi^2}P\,\LO_\xi=\int_{k^+/\vee}^1\frac{d\xi}{\xi^3}P(\xi)\frac{d^3N_{\rm LO}}{d^3k}(\frac{k^+}{\xi},\frac{\mathbf k}{\xi})$. This is a gluon of momentum $\frac{k}{\xi}$ that turns into one with the fraction $\xi$. For a hadron, $\frac{d^3N_h}{d^3k}=\int_{k^+/\vee}^1\frac{d\zeta}{\zeta^3}D^{\rm bare}_{h/g}(\zeta)\frac{d^3N_g}{d^3k}(\frac{k^+}{\zeta},\frac{\mathbf k}{\zeta})$.
\begin{itemize}
\item With $z=\zeta\xi$ the pole becomes $\frac{\alpha_s}{2\pi}(-\frac2\epsilon)\int_a^1\frac{dz}{z^3}\frac{d^3N_{\rm LO}}{d^3k}(\frac kz)[P_{gg}\otimes D](z)$, with $[P\otimes D](z)=\int_z^1\frac{d\zeta}{\zeta}P(\frac z\zeta)D(\zeta)$. This needs the plus distribution of Step 3 (\texttt{other\_check/checks.py}, item 5).
\item The $\overline{MS}$ counterterm of $D$ at $\mu_F$ is
\begin{equation}
D^{\rm bare}_{h/g}=D_{h/g}(\mu_F)+\frac{\alpha_s}{2\pi}\Big[\frac2\epsilon+\frac{d_\perp}{2}\log\frac{\mu^2}{\mu_F^2}\Big]P_{gg}\otimes D_{h/g}(\mu_F),\qquad\mu_F^2\frac{\partial D_{h/g}}{\partial\mu_F^2}=\frac{\alpha_s}{2\pi}P_{gg}\otimes D_{h/g}.\label{eq:D}
\end{equation}
For $\mu_F=\mu$ it is the pure pole. The log of $\mu/\mu_F$ carries the same $\frac{d_\perp}{2}$ as your bracket.
\end{itemize}
With $\ell_\Delta(\mu)=\ell_\Delta(\mu_F)+\log\frac{\mu^2}{\mu_F^2}$, the DGLAP term of \eqref{eq:split} plus \eqref{eq:D} leaves, exactly,
\begin{equation}
-\frac{\alpha_s}{2\pi}\int_{k^+/\vee}^1\frac{d\xi}{\xi^2}\,P_{gg}(\xi)\,\LO_\xi\Big[\frac{d_\perp}{2}\ell_\Delta(\mu_F)\Big].\label{eq:dglap0}
\end{equation}

\paragraph{Step 6: running coupling.} The $\overline{MS}$ counterterm of the coupling of LO at $\mu_R$ is $-\frac{\alpha_s\beta_0}{4\pi}\big[\frac2\epsilon+\frac{d_\perp}{2}\log\frac{\mu^2}{\mu_R^2}\big]\LO$. Added to the second term of \eqref{eq:split}:
\begin{equation}
\boxed{\frac{\alpha_s\beta_0}{4\pi}\,\LO\Big[\frac{d_\perp}{2}\log\frac{(y'-y)^2\mu_R^2}{c_0^2}\Big].}\label{eq:rc}
\end{equation}
The pole and $\mu$ cancel exactly. Since $\alpha_s(Q)=\alpha_s(\mu_R)[1-\frac{\alpha_s\beta_0}{4\pi}\log\frac{Q^2}{\mu_R^2}]$, at $\epsilon\to0$ this is $\LO|_{\alpha_s(\mu_R)}\to\LO[\alpha_s(c_0/|y'-y|)]$: the coupling of LO runs with $|y'-y|$, and $\mu_R$ drops out.
\begin{itemize}
\item Split at $\mu_F$:\vspace{-1ex}\begin{equation*}\eqref{eq:rc}=\frac{\alpha_s\beta_0}{4\pi}\frac{d_\perp}{2}\log\frac{\mu_R^2}{\mu_F^2}\LO+\frac{\alpha_s\beta_0}{4\pi}\LO\Big[\frac{d_\perp}{2}\ell_\Delta(\mu_F)\Big].\end{equation*}The first term is the renormalization-group term of \texttt{Other\_terms\_analysis.pdf}.
\item The bookkeeping of Steps 5--6 is checked symbolically with $d_\perp=2-\epsilon$ kept (\texttt{make\_note.py}).
\end{itemize}

\paragraph{Step 7: (L2) and the absence of $\log\frac{\vee}{k^+}$.} Only the $\frac{1-\xi}{\xi}$ term of $P_{gg}$ matters at small $\xi$.
\begin{itemize}
\item \textbf{Where the hidden log is.} In \eqref{eq:dglap0} this term gives $\int_a\frac{d\xi}{\xi^2}\frac{1-\xi}{\xi}\LO_\xi[\dots]$. For $\xi\to0$ the phase $e^{-ik\cdot\Delta/\xi}$ restricts $|y'-y|\lesssim\xi/|\mathbf k|$, so $\LO_\xi\propto\xi^2$. The integral then behaves as $\int_a\frac{d\xi}{\xi}=\log\frac{\vee}{k^+}$.
\item \textbf{How (L2) removes it.} (L2) is the same term at $y'=y$, i.e.\ $\LO^{(y)}_\xi$.
\begin{itemize}
\item Its $\mu$ and $\gamma_E$ parts multiply $\int_{y'}e^{-ik\cdot(y'-y)/\xi}\propto\delta^{(2)}(\mathbf k)$ and vanish for $\mathbf k\neq0$. So exactly
\begin{equation}
({\rm L2})=+\frac{\alpha_s}{2\pi}\,2N_c\int_a^1\frac{d\xi}{\xi^2}\frac{1-\xi}{\xi}\,\LO^{(y)}_\xi\Big[\frac{d_\perp}{2}\ell_\Delta(\mu_F)\Big].
\end{equation}
\item Combined with \eqref{eq:dglap0}, the small-$\xi$ integrand becomes $\frac{1-\xi}{\xi}\big(\LO_\xi-\LO^{(y)}_\xi\big)$.
\item This difference contains $U(y')-U(y)$ and $\frac{x'-y'}{(x'-y')^2}-\frac{x'-y}{(x'-y)^2}$, both $O(y'-y)=O(\xi)$. The $\xi$ integral therefore converges at $\xi\to0$, and the sum has a finite limit for $\vee\to\infty$.
\end{itemize}
\end{itemize}
The DGLAP term of the row is therefore
\begin{equation}
\boxed{\frac{d^3N}{d^3k}\bigg|_{\rm DGLAP}=-\frac{\alpha_s}{2\pi}\int_{k^+/\vee}^1\frac{d\xi}{\xi^2}\Big\{P_{gg}(\xi)\,\LO_\xi\Big[\frac{d_\perp}{2}\ell_\Delta(\mu_F)\Big]-2N_c\frac{1-\xi}{\xi}\,\LO^{(y)}_\xi\Big[\frac{d_\perp}{2}\ell_\Delta(\mu_F)\Big]\Big\}}\label{eq:dglap}
\end{equation}
with no $\log\frac{\vee}{k^+}$. (Doing the $y'$ integral of (L2) on its own, with the two-dimensional Fourier transform, gives (L2) $=-\frac{d_\perp}{2}\big[\log\frac{\vee}{k^+}-(1-\frac{k^+}{\vee})\big]H_4$, as in the previous version. That is the same subtraction written out: its $\log\frac{\vee}{k^+}$ cancels the hidden one of the first term, so it must not be listed separately.)

\paragraph{Step 8: (L3) and (L5).} Both have no pole and are kept as you wrote them, with $d_\perp=2-\epsilon$.
\begin{itemize}
\item (L3) $=-(1-\frac{k^+}{\vee})H_4$ is bounded for every $\vee$.
\item (L5): the $\Delta$ angular average of its bracket vanishes (not DGLAP), and its measure is $d\xi$, so it is finite for $\vee\to\infty$.
\end{itemize}

\section{Result in your notation}
After $D_{h/g}(\zeta,\mu_F)$ takes the $P_{gg}$ pole and the coupling counterterm takes the $\beta_0$ pole, your row is, exactly ($d_\perp=2-\epsilon$):
\begin{align}
\frac{d^{3}N^{\rm G2\,IV}}{d^{3}k}={}&\frac{g^{2}}{4\pi^{2}}\frac{N_{c}}{(2\pi)^{3}}\frac{g^{2}}{2\pi^{2}}\frac{1}{k^{+}}\int_{x,y,x',y'}\frac{(x-y)\cdot(x'-y')}{(x-y)^{2}(x'-y')^{2}}\left[\old{-\frac{d_{\perp}}{\epsilon}}-d_{\perp}\gamma_{E}-\frac{d_{\perp}}{2}\log\frac{(y'-y)^{2}\new{\mu_{F}^{2}}}{4}\old{-1}\right]\notag\\
&\quad\times\Bigg[\int_{\frac{k^{+}}{\vee}}^{1}\frac{d\xi}{\xi^{2}}e^{-ik\cdot(y'-y)\frac{1}{\xi}}\left[\xi(1-\xi)+\frac{(1-\xi)}{\xi}+\xi+\frac{\xi^{2}}{(1-\xi)_{+}}\right]+\Big(\log\Big(1-\frac{k^{+}}{\vee}\Big)\new{+\frac{11}{12}}\Big)e^{-ik\cdot(y'-y)}\Bigg]\notag\\
&\quad\times\left[U^{e'c}(y')-U^{e'c}(x')\right]\left[U^{ec}(y)-U^{ce}(x)\right]\rho^{e'}(x')\rho^{e}(x)\notag\\
&-\frac{g^{2}}{4\pi^{2}}\frac{N_{c}}{(2\pi)^{3}}\frac{g^{2}}{2\pi^{2}}\frac{1}{k^{+}}\int_{x,y,x',y'}\frac{(x-y)\cdot(x'-y)}{(x-y)^{2}(x'-y)^{2}}\left[-\frac{d_{\perp}}{2}\log\frac{(y'-y)^{2}\mu^{2}_{\overline{MS}}}{4}\right]\Bigg[\int_{\frac{k^{+}}{\vee}}^{1}\frac{d\xi}{\xi^{2}}e^{-ik\cdot(y'-y)\frac{1}{\xi}}\frac{(1-\xi)}{\xi}\Bigg]\notag\\
&\quad\times\left[U^{e'c}(y)-U^{e'c}(x')\right]\left[U^{ec}(y)-U^{ce}(x)\right]\rho^{e'}(x')\rho^{e}(x)\notag\\
&-\frac{1}{(2\pi)^{3}}\frac{g^{4}N_{c}}{8\pi^{5}}\frac{1}{k^{+}}\Big(1-\frac{k^{+}}{\vee}\Big)\int_{x,x',z,w,w'}e^{-ik\cdot(w'-w)}\frac{(w-z)\cdot(w'-z)}{(w-z)^{2}(w'-z)^{2}}\frac{(x-z)\cdot(x'-z)}{(x-z)^{2}(x'-z)^{2}}\notag\\
&\quad\times\left[U^{e'c}(z)-U^{e'c}(x')\right]\left[U^{ec}(z)-U^{ce}(x)\right]\rho^{e'}(x')\rho^{e}(x)\notag\\
&\new{+\frac{g^{2}}{4\pi^{2}}\frac{N_{c}}{(2\pi)^{3}}\frac{g^{2}}{2\pi^{2}}\frac{1}{k^{+}}\,\frac{11}{12}\int_{x,y,x',y'}e^{-ik\cdot(y'-y)}\frac{(x-y)\cdot(x'-y')}{(x-y)^{2}(x'-y')^{2}}\left[d_{\perp}\gamma_{E}+\frac{d_{\perp}}{2}\log\frac{(y'-y)^{2}\mu_{R}^{2}}{4}\right]}\notag\\
&\quad\new{\times\left[U^{e'c}(y')-U^{e'c}(x')\right]\left[U^{ec}(y)-U^{ce}(x)\right]\rho^{e'}(x')\rho^{e}(x)}\notag\\
&+\frac{g^{2}}{4\pi^{2}}\frac{1}{(2\pi)^{3}}\frac{g^{2}}{2\pi^{2}}\frac{N_{c}}{k^{+}}\int_{\frac{k^{+}}{\vee}}^{1}d\xi\frac{(1-\xi)}{\xi}\int_{x,y,x',y'}e^{-ik\cdot(y'-y)\frac{1}{\xi}}\Bigg[\frac{(x-y)\cdot(x'-y')}{(x-y)^{2}(x'-y')^{2}}-\frac{(x-y)\cdot(y'-y)\,(x'-y')\cdot(y'-y)}{(x-y)^{2}(x'-y')^{2}(y'-y)^{2}}\,d_{\perp}\Bigg]\notag\\
&\quad\times\left[U^{e'c}(y')-U^{e'c}(x')\right]\left[U^{ec}(y)-U^{ce}(x)\right]\rho^{e'}(x')\rho^{e}(x).\label{eq:final2}
\end{align}
\begin{itemize}
\item Lines 1--5: the DGLAP term \eqref{eq:dglap}. Lines 1--3 are your (L1)$+$(L4) without the pole, with $\mu\to\mu_F$ and \new{$+\frac{11}{12}$}, where $\frac{11}{12}=\frac{\beta_0}{2}\cdot\frac{1}{2N_c}$ is the $\delta(1-\xi)$ of $P_{gg}$. Lines 4--5 are your (L2), unchanged; its $\mu_{\overline{MS}}$ drops out for $\mathbf k\neq0$.
\item Lines 6--7: your (L3), unchanged.
\item Lines 8--9: the running coupling \eqref{eq:rc}, $\frac{\alpha_s\beta_0}{4\pi}=\frac{g^2N_c}{4\pi^2}\cdot\frac{11}{12}$.
\item Lines 10--11: your (L5), unchanged, with $d_\perp$.
\end{itemize}
The two \new{$\frac{11}{12}$} terms add up to $\frac{\alpha_s\beta_0}{4\pi}\frac{d_\perp}{2}\log\frac{\mu_R^2}{\mu_F^2}\LO$ and cancel for $\mu_R=\mu_F$. No term contains $\log\frac{\vee}{k^+}$, and the sum of lines 1--5 is finite for $\vee\to\infty$ (Step 7).
\end{document}
"""
open(os.path.join(HERE, 'G2IV_DGLAP_RC.tex'), 'w').write(T)
print('written;', CHECK)
