r"""Builds Other_terms_analysis.tex from your Other_terms.tex (input/) and the text below.
  python3 make_note.py && pdflatex Other_terms_analysis.tex (twice)
Your equations are copied verbatim; the only edits (Lambda -> 0 in the limits of d^3N_2b finite) are asserted and
shown red (removed) / blue (new)."""
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
src = open(os.path.join(HERE, 'input', 'Other_terms.tex')).read()
B = [m.group(1) for m in re.finditer(r'\\begin\{align\}(.*?)\\end\{align\}', src, re.S)]
assert len(B) == 28


def blk(i, edits=()):
    s = B[i]
    for a, b in edits:
        assert s.count(a) == 1, (i, a)
        s = s.replace(a, b)
    return '\\begin{align}' + s + '\\end{align}\n'


R, Bl = r'{\color{red}', r'{\color{blue}'
e2bfin = [(r'\int_{\Lambda}^{k^{+}-\Lambda}', R + r'\int_{\Lambda}^{k^{+}-\Lambda}}' + Bl + r'\int_{0}^{k^{+}}}'),
          (r'\int_{k^{+}-\Lambda}^{\Lambda}', R + r'\int_{k^{+}-\Lambda}^{\Lambda}}' + Bl + r'\int_{k^{+}}^{0}}')]

# Group II Row IV, your five lines
g2iv = B[21].split('\\notag \\\\\n')
assert len(g2iv) == 6
L1, L2, L3a, L3b, L4, L5 = [x.strip() for x in g2iv]
L1 = L1.replace(r'\frac{d^{3}N_{\rm NLO}^{\rm}}{d^{3}k}  & =', r'(L1)\quad & ').strip()
assert L1.strip().startswith('(L1)')

CLO = r"\left[U^{e'c}(y') - U^{e'c}(x')\right] \left[U^{ec}(y) - U^{ce}(x)\right] \rho^{e'}(x') \rho^{e}(x)"
H4 = (r"\frac{1}{(2\pi)^{3}}\,\frac{g^{4}\,N_{c}}{8\pi^{5}}\,\frac{1}{k^{+}}\int_{x,x',z,w,w'}e^{-ik\cdot (w'-w)}\,"
      r"\frac{(w-z)\cdot(w'-z)}{(w-z)^{2}\,(w'-z)^{2}}\,\frac{(x-z)\cdot(x'-z)}{(x-z)^{2}\,(x'-z)^{2}}\\&\qquad\times"
      r"\left[U^{e'c}(z) - U^{e'c}(x')\right]\left[U^{ec}(z) - U^{ce}(x)\right] \rho^{e'}(x') \rho^{e}(x)")
PREF = r"\frac{g^{2}}{4\pi^{2}}\,\frac{N_{c}}{(2\pi)^{3}}\,\frac{g^{2}}{2\pi^{2}}\,\frac{1}{k^{+}}"
KK = r"\frac{(x-y)\cdot (x'-y')}{(x-y)^{2}\,(x'-y')^{2}}"
BR = r"\left[\xi\,(1-\xi)+\frac{(1-\xi)}{\xi}+\xi+\frac{\xi^{2}}{(1-\xi)_{+}}\right]"

T = r"""\documentclass[10pt]{article}
\usepackage[a4paper,margin=0.25in]{geometry}
\usepackage{amsmath,amssymb,amsfonts,mathtools}
\usepackage[dvipsnames]{xcolor}
\usepackage{array,booktabs,float}
\usepackage[colorlinks=true,linkcolor=blue]{hyperref}
\usepackage{etoolbox}
\AtBeginEnvironment{align}{\small}
\allowdisplaybreaks
\numberwithin{equation}{section}
\setlength{\jot}{3pt}
\setlength{\emergencystretch}{3em}
\def\mf#1{\mathbb{#1}}
\def\wt#1{\widetilde{#1}}
\newcommand{\new}[1]{{\color{blue}#1}}
\newcommand{\old}[1]{{\color{red}#1}}
\newcommand{\LO}{{\rm LO}}
\newcommand{\Ehat}{\hat E}
\title{Your \texttt{Other\_terms.tex}: the $\mu$ coefficients, DGLAP from Group II Row IV,\\ the finite terms, and the running coupling}
\date{}
\begin{document}
\maketitle
\noindent Your file lists, row by row, what is left after the $\log\frac{k^+}{\Lambda}$ (JIMWLK) and $\log\frac{\vee}{k^+}$ (BFKL) terms are taken out. Groups I and III are with c.c., Group II is without. \new{Blue} = new, \old{red} = removed. Only $\vee>k^+$ and $p^+>\Lambda$ are used. $\vee$ is kept exact, and $\Lambda\to0$ is taken only where the limit exists term by term.

\paragraph{Your file is the checked list.} After removing the markup, your file is identical, equation by equation, to Secs.~2--4 (the remaining terms) of \texttt{Remaining\_RC\_DGLAP\_finite.pdf} (\texttt{diff\_check.py}). The differences are only:
\begin{itemize}
\item Group I Row V $d^3N_{2b\,\rm UV}$ keeps $2\log(1-\frac{\Lambda}{k^+})$, and $d^3N_{2b\,\rm finite}$ keeps $\Lambda$ in two limits. Both go to $0$ for $\Lambda\to0$ (Sec.~\ref{sec:fin}).
\item $d^3N_{2b\,\rm finite}$ writes the $\frac{11}{6}$ line with $+$, where the earlier file had $-[-\frac{11}{6}]$; the two are equal.
\end{itemize}
So everything below uses your equations as they are.

\section{Notation}
\begin{itemize}
\item $a\equiv\frac{k^+}{\vee}$, $\beta_0=\frac{11}{3}N_c$ ($n_f=0$), $\alpha_s=\frac{g^2}{4\pi}$, $d_\perp=2-\epsilon$.
\item $c_0\equiv2e^{-\gamma_E}$, $\ell_r(\mu)\equiv\log\frac{r^2\mu^2}{c_0^2}=2\gamma_E+\log\frac{r^2\mu^2}{4}$, $\hat E_r\equiv\frac{2}{\epsilon}+\ell_r(\mu)$.
\item $X=x-y$, $X'=x'-y'$, $\Delta=y'-y$; $\mu\equiv\mu_{\overline{MS}}$.
\end{itemize}
LO in your notation, and LO with a weight $f$ in the integrand:
\begin{align}
\LO[f]&\equiv\frac{1}{(2\pi)^{3}}\,\frac{g^{2}}{2\pi^{2}}\,\frac{1}{k^{+}}\int_{x,y,x',y'}e^{-ik\cdot(y'-y)}\,f\,"""+KK+r"""\,"""+CLO+r""",\label{eq:LO}
\end{align}
with $\LO\equiv\LO[1]$. ${\rm LO}_\xi[f]$ is the same with $e^{-ik\cdot(y'-y)}\to e^{-ik\cdot(y'-y)\frac1\xi}$. Since $\LO$ at momentum $(p^+,\mathbf p)$ carries $\frac{1}{p^+}e^{-i\mathbf p\cdot(y'-y)}$,
\begin{equation}
{\rm LO}_\xi=\frac{1}{\xi}\,\frac{d^3N_{\rm LO}}{d^3k}\Big(\frac{k^+}{\xi},\frac{\mathbf k}{\xi}\Big),\qquad d^3k\equiv dk^+d^2\mathbf k.\label{eq:LOxi}
\end{equation}
Your $H_4$ (the line that you write with $(1-\frac{k^+}{\vee})$):
\begin{equation}\begin{split}
H_4\equiv{}&"""+H4+r""".\end{split}\label{eq:H4}
\end{equation}
\paragraph{Splitting functions.}
\begin{equation}
P^{(r)}_{gg,+}(\xi)=2N_c\Big[\frac{\xi}{(1-\xi)_+}+\frac{1-\xi}{\xi}+\xi(1-\xi)\Big],\qquad P_{gg}(\xi)=P^{(r)}_{gg,+}(\xi)+\frac{\beta_0}{2}\,\delta(1-\xi).
\end{equation}
Here $\frac{1}{(1-\xi)_+}$ is the plus distribution of $[0,1]$. Cut at $\xi\ge a$ it acts as
\begin{equation}
\int_a^1d\xi\,\frac{g(\xi)}{(1-\xi)_+}=\int_a^1d\xi\,\frac{g(\xi)-g(1)}{1-\xi}+g(1)\log(1-a).\label{eq:plusstd}
\end{equation}
Your prescription has no $\log(1-a)$: $\int\frac{d\xi}{\xi^2}\frac{\xi^2}{(1-\xi)_+}G=\int\frac{d\xi}{1-\xi}[G(\xi)-G(1)]$. On $[a,1]$ the two are related exactly by (\texttt{checks.py}, item 3)
\begin{equation}
\int_a^1\frac{d\xi}{\xi^2}\Big[\xi+\frac{\xi^2}{(1-\xi)_+}\Big]_{\rm yours}G(\xi)=\int_a^1\frac{d\xi}{\xi^2}\,\frac{\xi}{(1-\xi)_+}\,G(\xi)\;-\;\log(1-a)\,G(1).\label{eq:plusrel}
\end{equation}

%=====================================================================
\section{(1) All coefficients of $\mu$}\label{sec:mu}
\paragraph{Where $\mu$ appears.} In your file $\mu$ appears in five places: the two orange brackets (Group I Row II part 1; Group I Row V $d^3N_{2b\,\rm UV}$) and lines (L1), (L2), (L4) of Group II Row IV. Nowhere else: the $d_\perp$ of Group I Row V 2a, 2b finite, Group II Row I $d^3N_{2f}$ and Group II Row IV (L5) multiplies finite integrands, so $d_\perp\to2$ there exactly for $\epsilon\to0$. Every $\mu$ comes in the form $\hat E_r=\frac2\epsilon+\ell_r(\mu)$, i.e.\ from $(\mu^2)^{\epsilon/2}/\epsilon$. So \textbf{the coefficient of $\log\mu^2$ equals the coefficient of $\frac{2}{\epsilon}$}, source by source.

\paragraph{Step 1: Group I Row II, part 1.} The colour factor with $w\to y$, $w'\to y'$ and $a\to c$, $b'\to e'$, $b\to e$ ($U^{bc}=U^{cb}$):
\begin{equation*}
\left[U^{ab'}(x')\rho^{b'}(x')-U^{ab'}(y')\rho^{b'}(x')\right]\left[U^{ab}(y)\rho^{b}(x)-U^{ab}(x)\rho^{b}(x)\right]=-\left[U^{e'c}(y')-U^{e'c}(x')\right]\left[U^{ec}(y)-U^{ce}(x)\right]\rho^{e'}(x')\rho^{e}(x).
\end{equation*}
The prefactor is $\frac{g^4N_c}{32\pi^4}=\frac{g^2}{2\pi^2}\cdot\frac{\alpha_sN_c}{4\pi}$; the measure, the phase and the kernel become those of \eqref{eq:LO}; your orange bracket is $-\hat E_{x-w}\to-\hat E_X$. The c.c.\ is the same with $X\to X'$. So
\begin{align}
\frac{d^3N_{1}}{d^3k}\bigg|_{\rm G1\,II}+{\rm c.c.}={}&\frac{\alpha_sN_c}{4\pi}\Big[\frac{11}{3}+\log\Big(1-\frac{k^+}{\vee}\Big)\Big]\LO\big[\hat E_X+\hat E_{X'}\big]+\frac{\alpha_s}{2\pi}N_c\Big[\frac{67}{9}-\frac{\pi^2}{3}\Big]\LO\notag\\
&-\frac{\alpha_s}{2\pi}N_c\log\Big(1-\frac{k^+}{\vee}\Big)\log\Big(\frac{\vee}{k^+}-1\Big)\LO.\label{eq:G1II}
\end{align}
The last term is $\log^2(1-\frac{k^+}{\vee})+\log(1-\frac{k^+}{\vee})\log\frac{\vee}{k^+}=\log(1-\frac{k^+}{\vee})\log(\frac{\vee}{k^+}-1)$.

\paragraph{Step 2: Group I Row V, $d^3N_{2b\,\rm UV}$.} With $y\to x$, $z\to y$, $w'\to y'$ the colour factor is $+$ that of LO. The prefactor is $\frac{g^4N_c}{16\pi^4}=\frac{g^2}{2\pi^2}\cdot\frac{\alpha_sN_c}{2\pi}$, and the orange bracket is $+\hat E_{y-z}\to\hat E_X$:
\begin{equation}
\frac{d^3N_{2b\,\rm UV}}{d^3k}+{\rm c.c.}=\frac{\alpha_sN_c}{2\pi}\Big[2\log\Big(1-\frac{\Lambda}{k^+}\Big)-\frac{11}{6}\Big]\LO\big[\hat E_X+\hat E_{X'}\big]\;\xrightarrow{\Lambda\to0}\;-\frac{\alpha_s\beta_0}{4\pi}\LO\big[\hat E_X+\hat E_{X'}\big].\label{eq:2bUV}
\end{equation}

\paragraph{Step 3: Group II Row IV, (L1) and (L4).} $\frac{g^2}{4\pi^2}=\frac{\alpha_s}{\pi}$, and your bracket is
\begin{equation}
-\frac{d_\perp}{\epsilon}-d_\perp\gamma_E-\frac{d_\perp}{2}\log\frac{\Delta^2\mu^2}{4}-1=-\frac{2}{\epsilon}+1-2\gamma_E-\log\frac{\Delta^2\mu^2}{4}-1+O(\epsilon)=-\hat E_\Delta+O(\epsilon).\label{eq:bracket}
\end{equation}
The $-1$ cancels the $+1$ of $-\frac{d_\perp}{\epsilon}$. With \eqref{eq:plusrel}, your (L4) is exactly the $\log(1-a)$ term of the standard plus distribution \eqref{eq:plusstd}:
\begin{equation}
({\rm L1})+({\rm L4})=\frac{\alpha_sN_c}{\pi}\int_a^1\frac{d\xi}{\xi^2}\Big[\xi(1-\xi)+\frac{1-\xi}{\xi}+\frac{\xi}{(1-\xi)_+}\Big]{\rm LO}_\xi\big[-\hat E_\Delta\big]=\frac{\alpha_s}{2\pi}\int_a^1\frac{d\xi}{\xi^2}\,P^{(r)}_{gg,+}(\xi)\,{\rm LO}_\xi\big[-\hat E_\Delta\big].\label{eq:L14}
\end{equation}

\paragraph{Step 4: Group II Row IV, (L2).} (L2) has $U(y)$ in place of $U(y')$ and the kernel $\frac{(x-y)\cdot(x'-y)}{(x-y)^2(x'-y)^2}$. So $y'$ appears only in the phase and in the log:
\begin{equation}
\int_{y'}e^{-ik\cdot(y'-y)\frac1\xi}\log\frac{(y'-y)^2\mu^2}{4}=-\frac{4\pi\xi^2}{\mathbf k^2}+\log\frac{\mu^2}{4}\,(2\pi)^2\xi^2\,\delta^{(2)}(\mathbf k)\qquad(\text{\texttt{checks.py}, item 4}).
\end{equation}
For $\mathbf k\neq0$ the $\mu$ term is zero: \textbf{(L2) has no $\mu$}. (L2) has no pole either, so $d_\perp\to2$ in it.

\paragraph{Step 5: the counterterms.} In $\overline{MS}$, the coupling of LO renormalized at $\mu_R$ and the fragmentation function at $\mu_F$ add (Sec.~\ref{sec:dglap}, Step 6)
\begin{equation}
-\frac{\alpha_s\beta_0}{4\pi}\Big[\frac{2}{\epsilon}+\log\frac{\mu^2}{\mu_R^2}\Big]\LO\quad\text{(coupling)},\qquad
+\frac{\alpha_s}{2\pi}\Big[\frac{2}{\epsilon}+\log\frac{\mu^2}{\mu_F^2}\Big]\int_a^1\frac{d\xi}{\xi^2}P_{gg}(\xi)\,{\rm LO}_\xi\quad\text{(fragmentation)}.\label{eq:cts}
\end{equation}

\paragraph{The table.} Coefficient of $\log\mu^2$, which is also the coefficient of $\frac{2}{\epsilon}$, in units of $\frac{\alpha_s}{2\pi}$. Here $P\otimes\LO\equiv\int_a^1\frac{d\xi}{\xi^2}P^{(r)}_{gg,+}(\xi)\,{\rm LO}_\xi$ (\texttt{checks.py}, item 1):
\begin{center}
\begin{tabular}{@{}lll@{}}\toprule
source & coefficient of $\log\mu^2$ (and of $\frac2\epsilon$) & what it is\\\midrule
G1 Row II part 1 $+$ c.c., \eqref{eq:G1II} & $N_c\big[\frac{11}{3}+\log(1-a)\big]\LO=\big[\beta_0+N_c\log(1-a)\big]\LO$ & $\beta_0$; Q1\\
G1 Row V $2b$ UV $+$ c.c., \eqref{eq:2bUV} & $2N_c\big[2\log(1-\frac{\Lambda}{k^+})-\frac{11}{6}\big]\LO\to-\beta_0\,\LO$ & $\beta_0$\\
G2 Row IV (L1)$+$(L4), \eqref{eq:L14} & $-P\otimes\LO$ & DGLAP\\
G2 Row IV (L2) & $0$ & ---\\
G2 Row IV (L3), (L5); all other rows & $0$ & ---\\\midrule
sum of your file & $N_c\log(1-a)\,\LO-P\otimes\LO$ \;\;($+4N_c\log(1-\frac{\Lambda}{k^+})\LO\to0$) & \\\midrule
coupling counterterm & $-\frac{\beta_0}{2}\LO$ & \\
fragmentation counterterm & $+P\otimes\LO+\frac{\beta_0}{2}\LO$ & \\\midrule
total & $N_c\log(1-a)\,\LO$ & Q1\\\bottomrule
\end{tabular}
\end{center}
\begin{itemize}
\item The two $\beta_0$ entries cancel exactly: pole, $\gamma_E$, $\mu$, and also $\log(x-w)^2$ against $\log(y-z)^2$, since \eqref{eq:G1II} and \eqref{eq:2bUV} have the same $\LO[\hat E_X+\hat E_{X'}]$.
\item The DGLAP entry $-P\otimes\LO$, with the $\frac{\beta_0}{2}$ of the coupling counterterm, is $-P_{gg}\otimes\LO$, and it is removed by the fragmentation counterterm.
\item \textbf{What is left is $\frac{\alpha_s}{2\pi}N_c\log(1-\frac{k^+}{\vee})\log\mu^2\,\LO$, together with its pole} (Q1, Eq.~\eqref{eq:q1}). It comes from the cut $\xi\ge\frac{k^+}{\vee}$ in the virtual term of Group I Row II. It vanishes for $\vee\to\infty$ at fixed $k^+$, but it is not zero at finite $\vee$ (Sec.~\ref{sec:open}).
\item After the counterterms, the scales that remain are $\mu_R$ and $\mu_F$: $\frac{\alpha_s\beta_0}{4\pi}\log\frac{\mu_R^2}{\mu_F^2}\LO$ (Sec.~\ref{sec:rc}) and $-\frac{\alpha_s}{2\pi}P\otimes\LO[\ell_\Delta(\mu_F)]$ (Sec.~\ref{sec:dglap}).
\end{itemize}

%=====================================================================
\section{(2) DGLAP from Group II Row IV, step by step}\label{sec:dglap}
\paragraph{Step 0: your row.} Your five lines are (L1); (L2), the subtraction at $y'=y$; (L3), the $H_4$ line; (L4), the $\log(1-\frac{k^+}{\vee})$ line; and (L5), the linear-polarization line:
\begin{align}
"""+L1+r"""\notag\\
(L2)\quad"""+L2.replace('&', '&', 1)+r"""\notag\\
(L3)\quad"""+L3a+r"""\notag\\
"""+L3b+r"""\notag\\
(L4)\quad"""+L4+r"""\notag\\
(L5)\quad"""+L5+r"""
\end{align}

\paragraph{Step 1: (L1) in units of LO.} $\frac{g^2}{4\pi^2}=\frac{\alpha_s}{\pi}$, and the $x,y,x',y'$ integrand of (L1) is that of \eqref{eq:LO} with $e^{-ik\cdot(y'-y)}\to e^{-ik\cdot(y'-y)\frac1\xi}$. So
\begin{equation}
({\rm L1})=\frac{\alpha_sN_c}{\pi}\int_a^1\frac{d\xi}{\xi^2}\Big[\xi(1-\xi)+\frac{1-\xi}{\xi}+\xi+\frac{\xi^2}{(1-\xi)_+}\Big]_{\rm yours}{\rm LO}_\xi\big[\mathcal B_\Delta\big],\qquad({\rm L4})=\frac{\alpha_sN_c}{\pi}\log(1-a)\,\LO\big[\mathcal B_\Delta\big],
\end{equation}
where $\mathcal B_\Delta$ is your bracket \eqref{eq:bracket}.

\paragraph{Step 2: (L4) completes the plus distribution.} By \eqref{eq:plusrel}, the bracket of (L1) is the standard one minus $\log(1-a)\,\delta(1-\xi)$, and ${\rm LO}_{\xi=1}=\LO$. (L4) adds this term back:
\begin{equation}
({\rm L1})+({\rm L4})=\frac{\alpha_s}{2\pi}\int_a^1\frac{d\xi}{\xi^2}\,P^{(r)}_{gg,+}(\xi)\,{\rm LO}_\xi\big[\mathcal B_\Delta\big].
\end{equation}
This is needed in Step 5. Without (L4) the pole does not factorize (\texttt{checks.py}, item 5: $30\%$ and $110\%$ mismatch).

\paragraph{Step 3: the pole.} By \eqref{eq:bracket}, $\mathcal B_\Delta=-\frac2\epsilon-\ell_\Delta(\mu)+O(\epsilon)$. The $O(\epsilon)$ multiplies a finite integral and drops out. Write $\ell_\Delta(\mu)=\ell_\Delta(\mu_F)+\log\frac{\mu^2}{\mu_F^2}$:
\begin{equation}
({\rm L1})+({\rm L4})=\frac{\alpha_s}{2\pi}\Big[-\frac2\epsilon-\log\frac{\mu^2}{\mu_F^2}\Big]\int_a^1\frac{d\xi}{\xi^2}P^{(r)}_{gg,+}(\xi)\,{\rm LO}_\xi\;-\;\frac{\alpha_s}{2\pi}\int_a^1\frac{d\xi}{\xi^2}P^{(r)}_{gg,+}(\xi)\,{\rm LO}_\xi\big[\ell_\Delta(\mu_F)\big].\label{eq:step3}
\end{equation}

\paragraph{Step 4: the convolution.} By \eqref{eq:LOxi},
\begin{equation}
\int_a^1\frac{d\xi}{\xi^2}\,P(\xi)\,{\rm LO}_\xi=\int_{k^+/\vee}^1\frac{d\xi}{\xi^3}\,P(\xi)\,\frac{d^3N_{\rm LO}}{d^3k}\Big(\frac{k^+}{\xi},\frac{\mathbf k}{\xi}\Big).
\end{equation}
This describes a gluon of momentum $(\frac{k^+}{\xi},\frac{\mathbf k}{\xi})$, produced at LO, that turns into a gluon with the fraction $\xi$ of its momentum. The lower limit says that the parent has $\frac{k^+}{\xi}\le\vee$. This is the structure of a fragmentation function.

\paragraph{Step 5: fragmentation.} For a hadron $h$ ($d^3k=dk^+d^2\mathbf k$, all components scaled by $\zeta$):
\begin{equation}
\frac{d^3N_h}{d^3k}=\int_{k^+/\vee}^1\frac{d\zeta}{\zeta^3}\,D^{\rm bare}_{h/g}(\zeta)\,\frac{d^3N_g}{d^3k}\Big(\frac{k^+}{\zeta},\frac{\mathbf k}{\zeta}\Big).\label{eq:frag}
\end{equation}
For the gluon of momentum $\frac{k}{\zeta}$ the cut in Step 4 is $\xi\ge\frac{a}{\zeta}$. Insert \eqref{eq:step3} and change variables $z=\zeta\xi$. The region $a\le\zeta\le1$, $\frac{a}{\zeta}\le\xi\le1$ is $a\le z\le\zeta\le1$, and $\frac{d\zeta}{\zeta^3}\frac{d\xi}{\xi^3}=\frac{dz}{z^3}\frac{d\zeta}{\zeta}$. Then
\begin{equation}
\int_a^1\frac{d\zeta}{\zeta^3}D(\zeta)\int_{a/\zeta}^1\frac{d\xi}{\xi^3}P(\xi)\frac{d^3N_{\rm LO}}{d^3k}\Big(\frac{k}{\zeta\xi}\Big)=\int_a^1\frac{dz}{z^3}\,\frac{d^3N_{\rm LO}}{d^3k}\Big(\frac{k}{z}\Big)\,\big[P\otimes D\big](z),\qquad\big[P\otimes D\big](z)=\int_z^1\frac{d\zeta}{\zeta}\,P\Big(\frac{z}{\zeta}\Big)D(\zeta).\label{eq:conv}
\end{equation}
The cut $\xi\ge\frac{a}{\zeta}$ becomes $z\ge a$, independent of $\zeta$. So $[P\otimes D]$ is the ordinary DGLAP convolution, with the plus distribution of $[0,1]$. This holds only if the cut plus distribution is \eqref{eq:plusstd}, which is why (L4) belongs here (\texttt{checks.py}, item 5: $2\times10^{-12}$).

\paragraph{Step 6: $\overline{MS}$ fragmentation function.} Add the coupling counterterm \eqref{eq:cts}. Its pole is $\frac{\alpha_s}{2\pi}(-\frac2\epsilon)\frac{\beta_0}{2}\LO$, the $\delta(1-\xi)$ part of $P_{gg}$; the $\beta_0$ of your Group I rows cancels by itself (Sec.~\ref{sec:rc}). The pole of \eqref{eq:step3} becomes $\frac{\alpha_s}{2\pi}(-\frac2\epsilon-\log\frac{\mu^2}{\mu_F^2})P_{gg}\otimes\LO$. In \eqref{eq:frag} it is removed by
\begin{equation}
D^{\rm bare}_{h/g}(z)=D_{h/g}(z,\mu_F)+\frac{\alpha_s}{2\pi}\Big[\frac2\epsilon+\log\frac{\mu^2}{\mu_F^2}\Big]\big[P_{gg}\otimes D_{h/g}\big](z,\mu_F).\label{eq:Dbare}
\end{equation}
\begin{itemize}
\item $D^{\rm bare}$ does not depend on $\mu_F$, so
\begin{equation}
\mu_F^2\frac{\partial D_{h/g}(z,\mu_F)}{\partial\mu_F^2}=\frac{\alpha_s}{2\pi}\big[P_{gg}\otimes D_{h/g}\big](z,\mu_F).
\end{equation}
This is the DGLAP equation.
\item Its $\frac{\beta_0}{2}$ part, with the coupling counterterm, leaves the renormalization-group term $\frac{\alpha_s\beta_0}{4\pi}\log\frac{\mu_R^2}{\mu_F^2}\LO$ (Sec.~\ref{sec:rc}).
\end{itemize}

\paragraph{Step 7: the DGLAP part of the row.} The DGLAP part of Group II Row IV is
\begin{equation}
\boxed{\frac{d^3N}{d^3k}\bigg|_{\rm DGLAP\ pole}=\frac{\alpha_s}{2\pi}\Big(-\frac{2}{\epsilon}-\log\frac{\mu^2}{\mu_F^2}\Big)\int_{k^+/\vee}^1\frac{d\xi}{\xi^3}\,P^{(r)}_{gg,+}(\xi)\,\frac{d^3N_{\rm LO}}{d^3k}\Big(\frac{k^+}{\xi},\frac{\mathbf k}{\xi}\Big).}
\end{equation}
With the $\frac{\beta_0}{2}$ of the coupling counterterm it is absorbed into $D_{h/g}(\zeta,\mu_F)$ through \eqref{eq:Dbare}. Its finite partner $-\frac{\alpha_s}{2\pi}P^{(r)}_{gg,+}\otimes\LO[\ell_\Delta(\mu_F)]$ stays in the gluon cross section. Its $\mu_F$ dependence is cancelled by the DGLAP evolution of $D_{h/g}$.

\paragraph{Step 8: (L2)$+$(L3).} The subtraction (L2) after the $y'$ integral (Step 4 of Sec.~\ref{sec:mu}; $\frac{d_\perp}{2}\to1$):
\begin{itemize}
\item $\int_{y'}e^{-ik\cdot(y'-y)\frac1\xi}\big[-\log\frac{(y'-y)^2\mu^2}{4}\big]=\frac{4\pi\xi^2}{\mathbf k^2}$;
\item $\int_a^1\frac{d\xi}{\xi^2}\frac{1-\xi}{\xi}\,\xi^2=\log\frac{\vee}{k^+}-\big(1-\frac{k^+}{\vee}\big)$;
\item $\int_{w,w'}e^{-ik\cdot(w'-w)}\frac{(w-z)\cdot(w'-z)}{(w-z)^2(w'-z)^2}=\frac{(2\pi)^2}{\mathbf k^2}$.
\end{itemize}
With the prefactor ratio $1$ (\texttt{checks.py}, item 2) and $y\to z$,
\begin{equation}
({\rm L2})=-\Big[\log\frac{\vee}{k^+}-\Big(1-\frac{k^+}{\vee}\Big)\Big]H_4,\qquad({\rm L3})=-\Big(1-\frac{k^+}{\vee}\Big)H_4,\qquad\boxed{({\rm L2})+({\rm L3})=-\log\frac{\vee}{k^+}\,H_4.}
\end{equation}
The BFKL list has $+\log\frac{\vee}{k^+}H_4$ for this row (\texttt{BFKL\_check.pdf}, C2). (L2)$+$(L3) removes from (L1) the $\log\frac{\vee}{k^+}$ of its small-$\xi$ end, where $e^{-ik\cdot(y'-y)/\xi}$ forces $y'\to y$. So the finite part of (L1)$+$(L4) together with (L2)$+$(L3) stays finite for $\vee\to\infty$.

\paragraph{The remaining terms of Group II Row IV.} After the DGLAP pole is absorbed in $D_{h/g}(\zeta,\mu_F)$ (Step 7), and with $\epsilon\to0$ exact in what is finite, your row is
\begin{align}
\frac{d^{3}N^{\rm G2\,IV}_{\rm rem}}{d^{3}k}  & = """+PREF+r"""\int_{x,y,x',y'}"""+KK+r"""\,\new{\left[-2\gamma_{E}-\log\left(\frac{(y'-y)^{2}\,\mu^{2}_{F}}{4}\right)\right]}\notag\\
&\quad\times\Bigg[\int_{\frac{k^{+}}{\vee}}^{1}\frac{d\xi}{\xi^{2}}\,e^{-ik\cdot (y'-y)\frac{1}{\xi}}"""+BR+r"""+\log\Big(1-\frac{k^{+}}{\vee}\Big)e^{-ik\cdot(y'-y)}\Bigg]\notag\\
&\qquad\times"""+CLO+r"""\notag\\
&\quad-\,"""+PREF+r"""\int_{x,y,x',y'}\frac{(x-y)\cdot (x'-y)}{(x-y)^{2}\,(x'-y)^{2}}\left[-\log\left(\frac{(y'-y)^{2}\,\mu^{2}}{4}\right)\right]\notag\\
&\qquad\times\int_{\frac{k^{+}}{\vee}}^{1}\frac{d\xi}{\xi^{2}}\,e^{-ik\cdot (y'-y)\frac{1}{\xi}}\,\frac{(1-\xi)}{\xi}\left[U^{e'c}(y) - U^{e'c}(x')\right] \left[U^{ec}(y) - U^{ce}(x)\right] \rho^{e'}(x') \rho^{e}(x)\notag\\
&\quad-\frac{1}{(2\pi)^{3}}\,\frac{g^{4}\,N_{c}}{8\pi^{5}}\,\frac{1}{k^{+}}\Big(1-\frac{k^{+}}{\vee}\Big)\int_{x,x',z,w,w'}e^{-ik\cdot (w'-w)}\,\frac{(w-z)\cdot(w'-z)}{(w-z)^{2}\,(w'-z)^{2}}\,\frac{(x-z)\cdot(x'-z)}{(x-z)^{2}\,(x'-z)^{2}}\notag\\
&\qquad\times\left[U^{e'c}(z) - U^{e'c}(x')\right]\left[U^{ec}(z) - U^{ce}(x)\right] \rho^{e'}(x') \rho^{e}(x)\notag\\
&\quad+"""+PREF+r"""\int_{\frac{k^{+}}{\vee}}^{1}d\xi\,\frac{(1-\xi)}{\xi}\int_{x,y,x',y'}e^{-ik\cdot (y'-y)\frac{1}{\xi}}\notag\\
&\qquad\times\Bigg["""+KK+r"""-\new{2}\,\frac{(x-y)\cdot (y'-y)\,(x'-y')\cdot (y'-y)}{(x-y)^{2}\,(x'-y')^{2}\,(y'-y)^{2}}\Bigg]\notag\\
&\qquad\times"""+CLO+r""".\label{eq:G2IVrem}
\end{align}
\begin{itemize}
\item Lines 1--2 are (L1)$+$(L4) without the pole, at $\mu_F$: $-\frac{\alpha_s}{2\pi}\int_{k^+/\vee}^1\frac{d\xi}{\xi^2}P^{(r)}_{gg,+}(\xi)\,{\rm LO}_\xi[\ell_\Delta(\mu_F)]$. Your bracket $[\xi+\frac{\xi^2}{(1-\xi)_+}]$ and the $\log(1-\frac{k^+}{\vee})$ term are kept as you wrote them.
\item Lines 3--6 are (L2) and (L3), $=-\log\frac{\vee}{k^+}H_4$. The $\mu$ in (L2) drops out (Step 4 of Sec.~\ref{sec:mu}).
\item Lines 7--8 are (L5) with $d_\perp=2$. They are finite and not DGLAP: the $\Delta$ angular average of the bracket is zero, and the Fourier transform gives the linearly polarized structure $\frac{2k^ik^j}{\mathbf k^2}-\delta^{ij}$. They are also finite for $\vee\to\infty$: the measure is $d\xi$, not $\frac{d\xi}{\xi^2}$, so the small-$\xi$ end behaves as $\int d\xi\,\xi(1-\xi)$.
\end{itemize}

%=====================================================================
\section{(3) The other finite terms of the cross section}\label{sec:fin}
\textbf{Approximations: none.} $\vee$ is exact. The only limits are $\Lambda\to0$ in the two limits of $d^3N_{2b\,\rm finite}$ and in $2\log(1-\frac{\Lambda}{k^+})$, and $\epsilon\to0$ in integrands without poles ($d_\perp\to2$). Both limits exist term by term: the integrands of $d^3N_{2b\,\rm finite}$ are regular at $p^+,q^+\to0,k^+$, the first through the plus prescription and the second because $[(k^+-q^+)(y-x)^2+q^+(y-z)^2]^{-1}$ is finite there. Your equations are copied unchanged.

\subsection{From Group I Row II part 1 (with c.c.)}
After the $\beta_0$ term (cancelled, Sec.~\ref{sec:rc}) and Q1 (Eq.~\eqref{eq:q1}), \eqref{eq:G1II} leaves two finite terms:
\begin{align}
\frac{d^3N_{\rm const}}{d^3k}&=\frac{1}{(2\pi)^{3}}\,\frac{g^{4}N_{c}}{32\pi^{4}}\,\frac{1}{k^{+}}\int_{w,w',x,x'}e^{-ik(w'-w)}\,\frac{(x'-w')\cdot (x-w)}{(x'-w')^{2}\,(x-w)^{2}}\Big[-\frac{67}{9}+\frac{\pi^{2}}{3}\Big]\notag\\
&\qquad\times\left[U^{ab'}(x')\,\rho^{b'}(x') - U^{ab'}(w')\,\rho^{b'}(x')\right]\left[U^{ab}(w)\,\rho^{b}(x) - U^{ab}(x)\,\rho^{b}(x)\right]+{\rm c.c.}=\frac{\alpha_s}{2\pi}\,N_c\Big[\frac{67}{9}-\frac{\pi^2}{3}\Big]\LO,\label{eq:const}\\
\frac{d^3N_{V_\vee}}{d^3k}&=\frac{1}{(2\pi)^{3}}\,\frac{g^{4}N_{c}}{32\pi^{4}}\,\frac{1}{k^{+}}\int_{w,w',x,x'}e^{-ik(w'-w)}\,\frac{(x'-w')\cdot (x-w)}{(x'-w')^{2}\,(x-w)^{2}}\,\log\Big(1-\frac{k^{+}}{\vee}\Big)\log\Big(\frac{\vee}{k^{+}}-1\Big)\notag\\
&\qquad\times\left[U^{ab'}(x')\,\rho^{b'}(x') - U^{ab'}(w')\,\rho^{b'}(x')\right]\left[U^{ab}(w)\,\rho^{b}(x) - U^{ab}(x)\,\rho^{b}(x)\right]+{\rm c.c.}\notag\\
&=-\frac{\alpha_s}{2\pi}\,N_c\log\Big(1-\frac{k^+}{\vee}\Big)\log\Big(\frac{\vee}{k^+}-1\Big)\LO.\label{eq:Vvee}
\end{align}
$V_\vee$ is bounded for every $\vee>k^+$ and vanishes for $\vee\to\infty$.

\subsection{Group I Row II part 2, $d^3N_{2b}$ (with c.c.)}\label{sub:first}
Unchanged (three $\rho$'s, no large logarithm):
"""+blk(1)+r"""where
"""+blk(2)+r"""
\subsection{Group I Row III (with c.c.)}
"""+blk(3)+blk(4)+blk(5)+blk(6)+r"""where
"""+blk(7)+r"""
\subsection{Group I Row IV (with c.c.)}
"""+blk(8)+r"""
\subsection{Group I Row V (with c.c.)}
\paragraph{First part.}
"""+blk(9)+r"""
\paragraph{Second part, $d^3N_{2a}$.}
"""+blk(10)+r"""
\paragraph{Second part, $d^3N_{2b\,\rm UV}$.} Nothing is left: $-\frac{11}{6}$ is the $\beta_0$ term, which cancels that of Row II part 1 (Sec.~\ref{sec:rc}), and $2\log(1-\frac{\Lambda}{k^+})\to0$.
\paragraph{Second part, $d^3N_{2b\,\rm finite}$} ($\Lambda\to0$ in the limits):
"""+blk(12, e2bfin)+r"""
The last line (with $\frac{11}{6}$) is the subtraction of the UV counterterm, written with the kernel $\frac{(y-x)\cdot(y-z)}{(y-x)^2(x-z)^2}$. Alone it diverges logarithmically at $x\to z$. Only the sum of all lines of $d^3N_{2b\,\rm finite}$ is finite, because the $\frac{1}{(x-z)^2}$ terms cancel for every $p^+$ (\texttt{rowV/Part2b\_MSbar.pdf}). Keep the lines together.

\subsection{Group II Row I}\label{sub:G2I}
"""+blk(13)+blk(14)+blk(15)+blk(16)+blk(17)+blk(18)+blk(19)+blk(20)+r"""
\subsection{Group II Row IV}
Equation~\eqref{eq:G2IVrem}.

\subsection{Group III (with c.c.)}\label{sub:G3}
\paragraph{Row I.}
"""+blk(22)+r"""
\paragraph{Row II.}
"""+blk(23)+r"""
\paragraph{Row III, first part.}
"""+blk(24)+r"""
\paragraph{Row III, second part.}
"""+blk(25)+r"""
\paragraph{Row V.}
"""+blk(26)+r"""
\paragraph{Row VI.}
"""+blk(27)+r"""
In Rows III, V, VI the $p^+$ integrands fall like $\frac{1}{p^{+2}}$, so these rows are finite for $\vee\to\infty$.
\paragraph{Three $\rho$'s.} Group I Row II $d^3N_{2b}$, Rows III, IV, V (first part), Group II Row I (first part) and Group III Rows I, II, III (first part), V, VI have three $\rho$'s. In the cross section they need $\langle\rho\rho\rho\rangle$, which is zero for a Gaussian projectile, or their two-$\rho$ commutator parts ($[\rho^a(x),\rho^b(y)]=if^{abc}\rho^c(x)\delta^{(2)}(x-y)$).

\subsection{The cross section}
\begin{equation}
\boxed{\begin{aligned}
\frac{d^3N_h}{d^3k}={}&\int_{k^+/\vee}^1\frac{d\zeta}{\zeta^3}\,D_{h/g}(\zeta,\mu_F)\,\frac{d^3N_g}{d^3k}\Big(\frac{k^+}{\zeta},\frac{\mathbf k}{\zeta}\Big)\qquad(d^3k=dk^+d^2\mathbf k),\\
\frac{d^3N_g}{d^3k}={}&\LO\big|_{\alpha_s(\mu_R)}+\frac{\alpha_s\beta_0}{4\pi}\log\frac{\mu_R^2}{\mu_F^2}\,\LO+\frac{\alpha_s}{2\pi}N_c\Big[\frac{67}{9}-\frac{\pi^2}{3}\Big]\LO\\
&+\mathcal N\Big[\log\frac{k^+}{\Lambda}\,{\rm JIMWLK}\otimes\LO+\log\frac{\vee}{k^+}\,{\rm BFKL}\otimes\LO\Big]\\
&+\frac{d^3N^{\rm G2\,IV}_{\rm rem}}{d^3k}\ \text{\eqref{eq:G2IVrem}}+\frac{d^3N_{V_\vee}}{d^3k}\ \text{\eqref{eq:Vvee}}+\text{the terms of Secs.~\ref{sub:first}--\ref{sub:G2I}, \ref{sub:G3}}\\
&+\frac{\alpha_sN_c}{4\pi}\log\Big(1-\frac{k^+}{\vee}\Big)\LO\big[\hat E_X+\hat E_{X'}\big]\quad\text{(Q1)}\;+\;O(\alpha_s^3).
\end{aligned}}\label{eq:result}
\end{equation}
\begin{itemize}
\item The Jacobian is $\frac{1}{\zeta^3}$ for $d^3k=dk^+d^2\mathbf k$; by \eqref{eq:LOxi} this is $\frac{1}{\zeta^2}$ in the ${\rm LO}_\zeta$ form. \old{The $\frac{d\zeta}{\zeta^2}$} in \texttt{Remaining\_RC\_DGLAP\_finite.pdf} (DGLAP and Result sections) and in the earlier finite-term notes (\texttt{Finite\_terms\_DGLAP\_RC}, \texttt{\_final}, \texttt{\_singlet}, \texttt{\_simplified}, \texttt{\_exact}, \texttt{Comparison\_AABKL}, \texttt{BFKL\_origin\_and\_AABKL}) is \new{$\frac{d\zeta}{\zeta^3}$}; corrected in all of them. AABKL's $\frac{d\zeta}{\zeta^2}$ is right for $\frac{d\sigma}{d^2p\,d\eta}$.
\item $\frac{d^3N^{\rm G2\,IV}_{\rm rem}}{d^3k}$ contains $-\log\frac{\vee}{k^+}H_4$, and the BFKL term of this row contains $+\log\frac{\vee}{k^+}H_4$.
\end{itemize}

%=====================================================================
\section{(4) Running coupling}\label{sec:rc}
\begin{enumerate}
\item \textbf{Where $\beta_0$ appears in your file.} In two places only, both orange: the $\frac{11}{3}$ of Group I Row II part 1 (virtual) and the $-\frac{11}{6}$ of Group I Row V $d^3N_{2b\,\rm UV}$ (real, the two gluons of the loop at short distance $x\to z$). By \eqref{eq:G1II} and \eqref{eq:2bUV},
\begin{equation}
\frac{\alpha_s\beta_0}{4\pi}\LO\big[\hat E_X+\hat E_{X'}\big]-\frac{\alpha_s\beta_0}{4\pi}\LO\big[\hat E_X+\hat E_{X'}\big]=0
\end{equation}
exactly, for every $\vee$: the pole, $\gamma_E$, $\mu$ and the logarithm of the transverse size all cancel. With $y\to x$, $z\to w$ the two rows have the same measure, phase, kernel and colour bracket up to a sign, and $\frac{g^4N_c}{16\pi^4}\frac{11}{6}=\frac{g^4N_c}{32\pi^4}\frac{11}{3}$. So these rows produce \textbf{no running coupling}: no $\beta_0\log(r^2\mu^2)$ term with $r=|x-w|$ survives, so nothing in them fixes the scale of $\alpha_s$ at a transverse size of the emission.
\item \textbf{The coupling of LO.} LO is $\propto g^2$. In $\overline{MS}$, $g^2\to g^2(\mu_R)\big[1-\frac{\alpha_s\beta_0}{4\pi}(\frac2\epsilon+\log\frac{\mu^2}{\mu_R^2})\big]$ gives the counterterm \eqref{eq:cts}. Its pole is the $\frac{\beta_0}{2}\delta(1-\xi)$ part of $P_{gg}$ in the DGLAP pole (Sec.~\ref{sec:dglap}, Step 6). After $D_{h/g}$ is renormalized at $\mu_F$, the only $\beta_0$ term left is the renormalization-group term
\begin{equation}
\boxed{\frac{d^3N}{d^3k}\bigg|_{\beta_0}=\frac{\alpha_s\beta_0}{4\pi}\log\frac{\mu_R^2}{\mu_F^2}\,\LO\big|_{\alpha_s(\mu_R)}.}
\end{equation}
Since $\frac{d}{d\log\mu_R^2}\LO|_{\alpha_s(\mu_R)}=-\frac{\alpha_s\beta_0}{4\pi}\LO$, the cross section does not depend on $\mu_R$ at this order. For $\mu_R=\mu_F$ the term is zero, and LO is taken with $\alpha_s(\mu_F)$.
\item \textbf{The constant.} The non-logarithmic part of Row II part 1 is \eqref{eq:const}, $\frac{\alpha_s}{2\pi}N_c(\frac{67}{9}-\frac{\pi^2}{3})\LO$. Here $N_c(\frac{67}{9}-\frac{\pi^2}{3})=2K$, where $K=N_c(\frac{67}{18}-\frac{\pi^2}{6})$ is the two-loop cusp constant ($n_f=0$).
\begin{itemize}
\item The constant depends on the scheme of $\alpha_s$. Under $\alpha_s\to\alpha_s(1+c\,\frac{\alpha_s}{2\pi})$ it shifts by $-c$.
\item It is absorbed completely by $\alpha_s\to\alpha_s(1+2K\frac{\alpha_s}{2\pi})$. That is twice the CMW shift $K\frac{\alpha_s}{2\pi}$ of soft-gluon emission, so here it is not the CMW scheme by itself. In $\overline{MS}$ it stays as a finite term.
\end{itemize}
\item \textbf{The $\frac{11}{6}$ in $d^3N_{2b\,\rm finite}$} has no pole and no $\mu$. It is the UV counterterm subtracted from the rest of $2b$ (Sec.~\ref{sec:fin}), not a running-coupling term. It must stay with the other lines of $d^3N_{2b\,\rm finite}$.
\item \textbf{Running coupling in the evolution.} The scale of $\alpha_s$ inside JIMWLK and BFKL (e.g.\ Balitsky's prescription) comes from the $\beta_0$ terms of the NLO kernels. They are of order $\alpha_s^2\delta Y$, beyond this calculation ($n_f=0$, one emission). At this order, running coupling enters only through $\alpha_s(\mu_R)$ in LO and the renormalization-group term above.
\end{enumerate}

\section{Open point}\label{sec:open}
\paragraph{Q1.} The pole and $\mu$ that the counterterms do not remove:
\begin{equation}
\frac{d^3N}{d^3k}\bigg|_{\rm Q1}=\frac{\alpha_sN_c}{4\pi}\log\Big(1-\frac{k^+}{\vee}\Big)\LO\big[\hat E_X+\hat E_{X'}\big]\label{eq:q1}
\end{equation}
\begin{itemize}
\item It is neither DGLAP (not $P_{gg}\otimes\LO$) nor running coupling (no $\beta_0$).
\item It comes from the cut $\xi\ge\frac{k^+}{\vee}$ in the virtual term of Group I Row II part 1, and vanishes for $\vee\to\infty$ at fixed $k^+$.
\item With the virtual term integrated over $\int_0^1d\xi$, as in AABKL (arXiv:1411.2869), both Q1 and $V_\vee$ are absent.
\end{itemize}

\section*{Reproducing}
In \texttt{nlo\_gluon\_paper/other\_check/}: \texttt{python3 diff\_check.py} (your file against the checked list); \texttt{python3 checks.py} ($\mu$ table, prefactors, plus prescriptions, Fourier identity, convolution $\to$ \texttt{checks.json}); \texttt{python3 make\_note.py} (this file).
\end{document}
"""
open(os.path.join(HERE, 'Other_terms_analysis.tex'), 'w').write(T)
print('written')
