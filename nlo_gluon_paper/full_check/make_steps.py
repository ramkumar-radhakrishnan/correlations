r"""Writes Steps_line1_and_mu.tex:
  Part A: from your own equations, step by step, where the extra log(vee/k+) terms of Group III Rows III, V, VI and
          Group II Row IV come from (no approximation; only exact identities).
  Part B: the coefficient in front of log mu^2 of the five sources (A3 and the shockwave through the loop, each with
          its c.c., and the real term), the coupling of LO, and what is running coupling and what is DGLAP.
Your equations are copied verbatim from input/; the numbers come from ../bfkl_rows/exact_check.json, mu_check.json,
ll_poles.json and ../dglap_rc/checks.json."""
import json
import os
import re
import sys
import make_ll as ML          # your new files (ll first)
import exact_check as EC

HERE = ML.HERE
tex_num = ML.tex_num
MU = json.load(open(os.path.join(HERE, 'mu_check.json')))
LINES = {g: open(os.path.join(HERE, 'input', '%s.tex' % g)).read().split('\n') for g in ('Group_I', 'Group_II', 'Group_III')}


def user_eq(g, a, b, num):
    """your align environment (lines a..b, 1-based, \\begin and \\end included), with your number shown"""
    t = '\n'.join(LINES[g][a - 1:b])
    t = t[t.index(r'\begin{align}') + len(r'\begin{align}'):t.rindex(r'\end{align}')]
    t = re.sub(r'\\label\{[^}]*\}', '', t)
    t = ML.GN.strip_user_colours(t)
    return '\\begin{align}' + t.rstrip() + r'\tag{your %s}' % num + '\n\\end{align}\n'


def disp(t, tag=None):
    t = ML.GN.strip_user_colours(t).strip()
    if t.startswith('&'):
        t = t[1:]
    t = re.sub(r'(\\notag|\\nonumber)?\s*\\\\\s*$', '', t)
    return '\\begin{align}\n&' + t + (r'\tag{%s}' % tag if tag else '') + '\n\\end{align}\n'


ECS = [ML.ECJ[k] for k in ML.ECJ if k != 'fourier']
G3V = EC.g3_display('G3 V')
G32 = EC.g3_display('G3 III(2)')
H4 = ML.H4 % ML.LVK


def mu_rows():
    c = MU['coefficient of log mu^2, units alpha_s/pi x LO']
    names = [('A3 (Row II part 1) + c.c.', r'A3, Group I Row II part 1 $+$ c.c., Eq.~\eqref{eq:muA3}'),
             ('shockwave through the loop (Row V 2b UV) + c.c.', r'shockwave through the loop, Group I Row V 2b UV $+$ c.c., Eq.~\eqref{eq:muSW}'),
             ('real, Group II Row IV line 1', r'real, Group II Row IV line 1, Eq.~\eqref{eq:muR1}'),
             ('real, Group II Row IV line 2', r'real, Group II Row IV line 2, Eq.~\eqref{eq:muR2}'),
             ('coupling counterterm of LO (mu_R)', r'coupling of LO, Eq.~\eqref{eq:ct}')]
    f = lambda s: '$' + s.replace('N_c', 'N_c').replace('*', '') + '$'
    out = []
    for k, lab in names:
        r = c[k]
        out.append(r'%s & %s & %s & %s & %s & %s\\' % (lab, f(r['beta0 (11/3)']), f(r['x log(k+/Lambda)']),
                                                      f(r['x log(vee/k+)']), f(r['x log(1-k+/vee)']), f(r['x P conv'])))
    return '\n'.join(out)


DOC = ML.PREAMBLE + r'''
\newcommand{\Kc}{\mathcal K}
\title{Step by step from your equations:\\ where the extra $\log\frac{\vee}{k^+}$ terms of Group III Rows III, V, VI and Group II Row IV come from;\\ and the $\mu$ bookkeeping (running coupling, DGLAP)}
\date{}
\begin{document}
\maketitle
\noindent Your notation throughout. Your plus prescription:
\begin{equation}
\int dp^+\left[\frac{1}{p^+}\right]_+f(p^+)\equiv\int dp^+\,\frac{f(p^+)-f(0)}{p^+},\qquad
\int\frac{d\xi}{\xi^2}\,\frac{\xi}{(1-\xi)_+}\,g(\xi)\equiv\int\frac{d\xi}{\xi(1-\xi)}\,\big[g(\xi)-g(1)\big].
\label{eq:plus}
\end{equation}
Your \texttt{BFKL.tex} calls the Group III rows ``III, IV, V''; in \texttt{Group\_III.tex} they are Rows III, V, VI (Row IV has only four $\rho$'s). Equation numbers ``your (1.x)'' are those of your files as they compile.

\paragraph{Your ``first part of the calculation''.} The \texttt{.tex} file you mention did not reach me: the last files I received are \texttt{Group\_I/II/III.tex} and \texttt{4rho.tex} (2.30.24 AM). Please upload it again and I will check it against the three Groups.

\section*{Part A. The extra terms in line 1}
\section{The answer in one paragraph}
Nothing is added to your cross section. Every ``extra'' term is a piece of your line 1, the $p^+$ (resp.\ $\xi$) integral that runs up to $\vee-k^+$ (resp.\ down to $\frac{k^+}{\vee}$). It is written out because it grows like $\log\frac{\vee}{k^+}$. There are three of them in each row:
\begin{description}
\item[(i)] The $\log\frac{\vee}{k^+}$ in your line 2, from $\log(\frac{\vee}{k^+}-1)=\log\frac{\vee}{k^+}+\log(1-\frac{k^+}{\vee})$.
\item[(ii)] Minus the same, hidden in the subtraction $-f(0)$ of $\left[\frac{1}{p^+}\right]_+$ in line 1. (i) and (ii) cancel exactly. So the $\log\frac{\vee}{k^+}$ you see in line 2 is \emph{not} a $\vee$ dependence of the row: before your plus prescription, Eq.~(your 1.12), there is no such logarithm (Step 3).
\item[(iii)] The term (B), hidden in the $-\frac{1}{k^+}\delta_{im}\delta_{jk'}$ term of line 1. This term has no $\frac{1}{p^+}$, but after the $y'$ integral it becomes \emph{exactly} $\frac{k^+}{p^++k^+}\times$ (a $p^+$-independent integral) $+$ (a part that falls faster). Then $\int_0^{\vee-k^+}\frac{dp^+}{p^++k^+}=\log\frac{\vee}{k^+}$ (Step 4). This is the genuine BFKL term of the row.
\end{description}
In Group II Row IV, (i) and (ii) are the same with $\xi$. (iii) is $H_4$, from the $\log(y'-y)^2$ of the pole bracket in the $\frac{1-\xi}{\xi}$ term (Sec.~\ref{sec:G2IV}).

\section{Group III Row V, step by step}\label{sec:V}
\paragraph{Step 1: your equation after the $w'$ integral.}
@@EQ112@@
After the Wilson lines are simplified (your Eq.~(1.14)), the colour factor is
\begin{equation}
\mathcal C(y')=\Big[U^{e'c'}(y')\,\rho^{e'}(x')-U^{c'e'}(x')\,\rho^{e'}(x')\Big]\times\Big[U^{ce}(x)\,U^{ad}(y)\,\rho^{e}(x)\rho^{d}(y)-U^{ad}(w)\,U^{ce}(x)\,\rho^{e}(x)\rho^{d}(y)\Big],
\label{eq:C}
\end{equation}
with $f^{c'ca}$ and $C\equiv\frac{1}{(2\pi)^{3}}\frac{ig^{4}f^{c'ca}}{4\pi^{4}}\frac{1}{k^{+}}$ in front. Only $U^{e'c'}(y')$ depends on $y'$.

\paragraph{Step 2: names.} Contract the Kronecker deltas of the three terms of the bracket. For each $p^+$,
\begin{align}
F_a(p^+)&=\int_{w,y',z,y,x,x'}e^{-ik\cdot (y'-w)}\, e^{-ik\cdot \frac{p^{+}}{k^{+}}(y'-z)}\,\frac{(y'-z)\cdot (x'-y')\,(y-w)\cdot (x-z)}{(y'-z)^{2}(x'-y')^{2}(y-w)^{2}(x-z)^{2}}\,\mathcal C(y') &&(\delta_{k'm}\delta_{ij}),\notag\\
F_b(p^+)&=\int_{w,y',z,y,x,x'}e^{-ik\cdot (y'-w)}\, e^{-ik\cdot \frac{p^{+}}{k^{+}}(y'-z)}\,\frac{(y'-z)\cdot (x-z)\,(x'-y')\cdot (y-w)}{(y'-z)^{2}(x'-y')^{2}(y-w)^{2}(x-z)^{2}}\,\mathcal C(y') &&(\delta_{jm}\delta_{ik'}),\notag\\
F_c(p^+)&=\int_{w,y',z,y,x,x'}e^{-ik\cdot (y'-w)}\, e^{-ik\cdot \frac{p^{+}}{k^{+}}(y'-z)}\,\frac{(y'-z)\cdot (y-w)\,(x'-y')\cdot (x-z)}{(y'-z)^{2}(x'-y')^{2}(y-w)^{2}(x-z)^{2}}\,\mathcal C(y') &&(\delta_{im}\delta_{jk'}).
\label{eq:F}
\end{align}
Your Eq.~(1.12), with the Wilson lines simplified, is then
\begin{equation}
\frac{d^3N_{\rm V}}{d^3k}=C\int_{\Lambda}^{\vee-k^{+}}\frac{dp^{+}}{2\pi}\bigg[\frac{F_a(p^+)}{p^++k^+}-\frac{F_b(p^+)}{p^+}-\frac{F_c(p^+)}{k^+}\bigg].
\label{eq:start}
\end{equation}

\paragraph{Step 3: your plus prescription, and what it does to $\vee$.} In your notation, exactly,
\begin{equation}
\int_\Lambda^{\vee-k^+}dp^+\,\frac{F_b(p^+)}{p^+}=\int_\Lambda^{\vee-k^+}dp^+\left[\frac{1}{p^+}\right]_+F_b(p^+)+F_b(0)\int_\Lambda^{\vee-k^+}\frac{dp^+}{p^+},\qquad
\int_\Lambda^{\vee-k^+}\frac{dp^+}{p^+}=\log\frac{k^+}{\Lambda}+\log\Big(\frac{\vee}{k^+}-1\Big).
\label{eq:step3}
\end{equation}
With $\int_\Lambda\to\int_0$ in the plus-prescribed term ($O(\Lambda)$), Eq.~\eqref{eq:start} is your Eq.~(1.14):
\begin{align}
\text{line 1}&=C\int_{0}^{\vee-k^{+}}\frac{dp^{+}}{2\pi}\bigg[\frac{F_a(p^+)}{p^++k^+}-\frac{F_b(p^+)-F_b(0)}{p^+}-\frac{F_c(p^+)}{k^+}\bigg],\notag\\
\text{line 2}&=-\frac{C}{2\pi}\,F_b(0)\,\Big[\log\frac{k^+}{\Lambda}+\log\Big(\frac{\vee}{k^+}-1\Big)\Big].
\label{eq:lines}
\end{align}
($\frac{C}{2\pi}=\frac{1}{(2\pi)^3}\frac{ig^4f^{c'ca}}{8\pi^5}\frac{1}{k^+}$, and $F_b(0)$, without the $p^+$ phase, is the integral of your line 2.)

\textbf{The point.} The left side of Eq.~\eqref{eq:step3} does \emph{not} grow like $\log\vee$: its integrand $\frac{F_b(p^+)}{p^+}$ falls at least like $\frac{1}{p^{+2}}$ (Step 4). The right side has $+F_b(0)\log\frac{\vee}{k^+}$ explicitly. So the plus-prescribed integral must contain $-F_b(0)\log\frac{\vee}{k^+}$: it subtracts $\frac{F_b(0)}{p^+}$ all the way up to $p^+=\vee-k^+$. Exactly: since $\frac{1}{p^+}=\frac{1}{p^++k^+}+\frac{k^+}{p^+(p^++k^+)}$ and $\int_0^{\vee-k^+}\frac{dp^+}{p^++k^+}=\log\frac{\vee}{k^+}$,
\begin{equation}
\int_0^{\vee-k^+}dp^+\,\frac{F_b(p^+)-F_b(0)}{p^+}=\int_0^{\vee-k^+}dp^+\bigg[\frac{k^+\big(F_b(p^+)-F_b(0)\big)}{p^+(p^++k^+)}+\frac{F_b(p^+)}{p^++k^+}\bigg]\;-\;F_b(0)\,\log\frac{\vee}{k^+}.
\label{eq:I2}
\end{equation}
This is identity (I2). Both integrands on the right fall like $\frac{1}{p^{+2}}$, so they are finite for $\vee\to\infty$. Put Eq.~\eqref{eq:I2} into line 1, Eq.~\eqref{eq:lines}:
\begin{align}
\text{(i)}\;&=-\frac{C}{2\pi}\,F_b(0)\log\frac{\vee}{k^+}\qquad\text{(the }\log\tfrac{\vee}{k^+}\text{ of line 2)},\notag\\
\text{(ii)}\;&=+\frac{C}{2\pi}\,F_b(0)\log\frac{\vee}{k^+}\qquad\text{(from }-\left[\tfrac{1}{p^+}\right]_+\text{ in line 1)},\qquad \text{(i)}+\text{(ii)}=0.
\label{eq:i_ii}
\end{align}
Written out, (i) is the first and (ii) the second $\log\frac{\vee}{k^+}$ line of Group III Row V in \texttt{Checks\_JIMWLK\_BFKL.pdf}.

\paragraph{Step 4: the $y'$ integral, exactly.} The two phases combine to one:
\begin{equation}
e^{-ik\cdot (y'-w)}\,e^{-ik\cdot\frac{p^+}{k^+}(y'-z)}=e^{-ik\cdot(z-w)}\,e^{-iq\cdot(y'-z)},\qquad q=t\,k,\quad t=\frac{p^++k^+}{k^+}.
\end{equation}
In each $F_n$ the integrand is $\frac{(y'-z)^m}{(y'-z)^2}G_n^m(y')$, where $G_n(y')$ collects $\frac{(x'-y')^{k'}}{(x'-y')^2}$ and $\mathcal C(y')$. Write, identically,
\begin{equation}
G_n(y')=G_n(z)+\big[G_n(y')-G_n(z)\big].
\label{eq:split}
\end{equation}
For the $G_n(z)$ part the $y'$ integral is a pure 2-D Fourier transform. With $\int d^2r\,e^{-iq\cdot r}\frac{r^m}{r^2}=-2\pi i\frac{q^m}{q^2}$ (from $\frac{r^m}{r^2}=\frac12\partial^m\log r^2$ and $\int d^2r\,e^{-iq\cdot r}\log r^2=-\frac{4\pi}{q^2}$):
\begin{equation}
\int_{y'}e^{-iq\cdot(y'-z)}\frac{(y'-z)^m}{(y'-z)^2}=-2\pi i\,\frac{1}{t}\,\frac{k^m}{k^2}=\frac{k^+}{p^++k^+}\int_{w'}e^{-ik\cdot(w'-z)}\frac{(w'-z)^m}{(w'-z)^2}.
\label{eq:I3}
\end{equation}
This is identity (I3). It is exact; no large-$p^+$ limit is taken. Hence, for $n=a,b,c$,
\begin{equation}
F_n(p^+)=\frac{k^+}{p^++k^+}\,F_n^{z}+F_n^{\Delta}(p^+),
\label{eq:Fsplit}
\end{equation}
with
\begin{align}
F_c^{z}&=\int_{w,w',z,y,x,x'}e^{-ik\cdot (w'-w)}\,\frac{(w'-z)\cdot (y-w)\,(x'-z)\cdot (x-z)}{(w'-z)^{2}(y-w)^{2}(x'-z)^{2}(x-z)^{2}}\,\mathcal C(z),\notag\\
F_c^{\Delta}(p^+)&=\int_{w,y',z,y,x,x'}e^{-ik\cdot (y'-w)}e^{-ik\cdot\frac{p^+}{k^+}(y'-z)}\frac{(y'-z)^m(y-w)^m}{(y'-z)^2(y-w)^2}\,\frac{(x-z)^{k'}}{(x-z)^2}\bigg[\frac{(x'-y')^{k'}}{(x'-y')^2}\mathcal C(y')-\frac{(x'-z)^{k'}}{(x'-z)^2}\mathcal C(z)\bigg],
\label{eq:Fc}
\end{align}
where $\mathcal C(z)$ is $\mathcal C(y')$ with $U^{e'c'}(y')\to U^{e'c'}(z)$, and similarly for $a$, $b$. $F_n^{\Delta}(p^+)$ falls faster than $\frac{1}{p^+}$:
\begin{itemize}
\item $G_n(y')-G_n(z)$ vanishes at $y'=z$, so the product with $\frac{(y'-z)^m}{(y'-z)^2}$ is bounded there, and its Fourier transform at $q=tk$ falls like $\frac{1}{q^2}\propto\frac{k^{+2}}{(p^++k^+)^2}$;
\item at $y'=x'$ the factor $\frac{(x'-y')^{k'}}{(x'-y')^2}$ gives a $\frac{1}{q}$ term, but with the phase $e^{-iq\cdot(x'-z)}$, which oscillates in $p^+$. Its $p^+$ integral converges.
\end{itemize}
So $\int^{\vee}dp^+\,\frac{F_n^\Delta}{k^+}$, $\int^\vee dp^+\frac{F_n(p^+)}{p^++k^+}$ and $\int^\vee dp^+\frac{F_n(p^+)}{p^+}$ are all finite for $\vee\to\infty$. The only $p^+$ integrand of the row that falls like $\frac{1}{p^+}$ is the $G(z)$ part of the $-\frac{1}{k^+}$ term:
\begin{equation}
-\frac{C}{k^+}\int_0^{\vee-k^+}\frac{dp^+}{2\pi}\,\frac{k^+}{p^++k^+}\,F_c^z=-\frac{C}{2\pi}\,F_c^z\,\log\frac{\vee}{k^+}\qquad\text{exactly}.
\label{eq:iii}
\end{equation}

\paragraph{Step 5: (iii), the term (B), written out.} Eq.~\eqref{eq:iii} with Eq.~\eqref{eq:Fc}:
@@BV@@
This is the blue (iii) line of Group III Row V in \texttt{Checks\_JIMWLK\_BFKL.pdf}.

\paragraph{Step 6: the row.} Steps 3--5 give, exactly,
\begin{equation}
\frac{d^3N_{\rm V}}{d^3k}=\underbrace{-\frac{C}{2\pi}F_b(0)\log\frac{k^+}{\Lambda}}_{\text{JIMWLK}}\;+\;\underbrace{\text{(i)}+\text{(ii)}}_{=0}\;+\;\underbrace{\text{(B)}_{\rm V}}_{\text{BFKL}}\;+\;R_{\rm V},
\end{equation}
where $R_{\rm V}$ is line 1 with $\left[\frac{1}{p^+}\right]_+$ replaced as in Eq.~\eqref{eq:I2}, and with $G(y')\to G(y')-G(z)$ in the $-\frac{1}{k^+}$ term, plus line 2 with $\log(1-\frac{k^+}{\vee})$:
@@RV@@
Every $p^+$ integrand in $R_{\rm V}$ falls at least like $\frac{1}{p^{+2}}$ or oscillates, so $R_{\rm V}$ is finite for $\vee\to\infty$. ($\left[\frac{k^+}{p^+(p^++k^+)}\right]_+F\equiv\frac{k^+[F(p^+)-F(0)]}{p^+(p^++k^+)}$.) This is the uniqueness argument: had we put a different amount of $\log\frac{\vee}{k^+}$ into the BFKL list, $R$ would grow like $\log\vee$.

\paragraph{Numerical check} (\texttt{bfkl\_rows/exact\_check.py}; line 1 of your new file is identical to the one checked there). Two configurations, $\vee=20k^+$, the $p^+$ integral done numerically, pointwise in the transverse points:
\begin{itemize}
\item line 1 $=$ (rewritten line 1) $+$ ($G(y')-G(z)$ term) $+$ (ii) $+$ ($G(z)$ term), Rows III(1), III(2), V, VI: $\le@@ECID@@$;
\item the $G(z)$ term after the Fourier identity \eqref{eq:I3} and the $p^+$ integral $=$ (B): $\le@@ECB@@$;
\item the Fourier transforms themselves: @@FOU@@
\end{itemize}

\section{Group III Rows VI and III (first part)}
Your Eqs.~(1.15)--(1.17) (Row VI) and (1.5)--(1.7) (Row III part 1) have the same three terms in the bracket. The steps are those of Sec.~\ref{sec:V} word for word. Only the colour factor $\mathcal C(y')$ changes; it is $\big[U^{e'c'}(y')U^{d'b}(z)\rho^{e'}(x')-U^{bd'}(z)U^{c'e'}(x')\rho^{e'}(x')\big]\times[\dots]$ (and $f^{c'd'a}$). The results:
@@BVI@@
@@BIII1@@
\noindent (i) and (ii) cancel as in Eq.~\eqref{eq:i_ii}.

\section{Group III Row III, second part}
Your Eq.~(1.10), line 1:
@@EQ110@@
The bracket has six terms. After the $y'$ integral, Eq.~\eqref{eq:Fsplit}, each $G(z)$ part gets $\frac{k^+}{p^++k^+}$; the common factor is $\frac{1}{p^+(x-z)^2+k^+(x-w)^2}\sim\frac{1}{p^+}$. The large-$p^+$ power of each term, exactly:
\begin{center}\renewcommand{\arraystretch}{1.5}
\begin{tabular}{lll}\toprule
term & coefficient $\times\frac{1}{p^+A+k^+B}\times\frac{k^+}{p^++k^+}$ & $\int^{\vee}dp^+$\\\midrule
$\frac{\delta_{k'm}d_\perp p^+}{2(z-w)^2(p^++k^+)^2}[\dots]$ & $\sim p^{+-3}$ & finite\\
$\frac{p^+\delta_{k'm}}{k^+(p^++k^+)}[\dots]$ & $\sim p^{+-2}$ & finite\\
$\frac{\delta_{k'm}}{p^++k^+}[\dots]$ & $\sim p^{+-3}$ & finite\\
$-\left[\frac{1}{p^+}\right]_+[\dots]$ & $F(p^+)$ part $\sim p^{+-3}$; $-F(0)$ part: (ii), as in Eq.~\eqref{eq:I2} & (ii) $=-$(i)\\
$-\frac{p^+}{k^{+2}}[\dots]$ & $\sim p^{+-1}$ & $\log\frac{\vee}{k^+}$: (iii)\\
$-\frac{1}{k^+}[\dots]$ & $\sim p^{+-2}$ & finite\\\bottomrule
\end{tabular}
\end{center}
with $A=(x-z)^2$, $B=(x-w)^2$. For the $-\frac{p^+}{k^{+2}}$ term, the partial fractions (I5) are exact:
\begin{equation}
-\frac{p^+}{k^{+2}}\,\frac{k^+}{p^++k^+}\,\frac{1}{p^+A+k^+B}=-\frac{1}{k^+A}\,\frac{1}{p^++k^+}+\frac{B}{A\,(p^++k^+)(p^+A+k^+B)}.
\label{eq:I5}
\end{equation}
The first term gives $\int_0^{\vee-k^+}\frac{dp^+}{p^++k^+}=\log\frac{\vee}{k^+}$, i.e.\ (B); the second falls like $\frac{1}{p^{+2}}$ and stays in the remaining terms. Written out:
@@BIII2@@
In the remaining terms (\texttt{Remaining\_RC\_DGLAP\_finite.pdf}, Group III Row III part 2) the $-\frac{p^+}{k^{+2}}$ term becomes its $G(y')-G(z)$ part, and the second term of Eq.~\eqref{eq:I5} is added:
@@R32@@

\section{Group II Row IV}\label{sec:G2IV}
\paragraph{Step 1: your equations.} With $\xi=\frac{k^+}{p^++k^+}$, your Eq.~(1.19), and after the plus prescription, your Eq.~(1.21):
@@EQ121@@
Call $g(\xi)$ the $\xi$-dependent part of the integrand of line 1 at fixed transverse points: $g(\xi)\propto e^{-ik\cdot(y'-y)\frac1\xi}$. Line 2 is $g(1)$ times $\int_{\frac{k^+}{\vee}}^{\frac{k^+}{k^++\Lambda}}\frac{d\xi}{\xi(1-\xi)}=\log\frac{k^+}{\Lambda}+\log(\frac{\vee}{k^+}-1)$: your plus prescription, Eq.~\eqref{eq:plus}, in the $\xi$ form of Eq.~\eqref{eq:step3}.

\paragraph{Step 2: (i) and (ii).} $\frac{1}{\xi(1-\xi)}=\frac{1}{\xi}+\frac{1}{1-\xi}$ and $\int_{\frac{k^+}{\vee}}^1\frac{d\xi}{\xi}=\log\frac{\vee}{k^+}$ give, exactly,
\begin{equation}
\int_{\frac{k^+}{\vee}}^1\frac{d\xi}{\xi^2}\frac{\xi}{(1-\xi)_+}g=\int_{\frac{k^+}{\vee}}^1\frac{d\xi}{\xi^2}\Big(\xi+\frac{\xi^2}{(1-\xi)_+}\Big)g\;-\;g(1)\log\frac{\vee}{k^+},\qquad \int\frac{d\xi}{\xi^2}\frac{\xi^2}{(1-\xi)_+}g\equiv\int\frac{d\xi}{1-\xi}\big[g(\xi)-g(1)\big].
\label{eq:I2xi}
\end{equation}
The last term, times the prefactor of line 1, is (ii). It is your \texttt{BFKL.tex} line 1 for this row; it is exactly minus (i), the $\log\frac{\vee}{k^+}$ of line 2, poles included ($@@ECT0@@$).

\paragraph{Step 3: (iii), $H_4$.} Which $\xi$ integrand of line 1 falls like $\frac{1}{\xi}$ for $\xi\to\frac{k^+}{\vee}\to0$? The $y'$ integral is a Fourier transform at $q=\frac{k}{\xi}$, large for small $\xi$. Write $H(y,y')=\frac{(x-y)\cdot(x'-y')}{(x-y)^2(x'-y')^2}[U^{e'c}(y')-U^{e'c}(x')][U^{ec}(y)-U^{ce}(x)]\rho^{e'}(x')\rho^e(x)=H(y,y)+[H(y,y')-H(y,y)]$.
\begin{itemize}
\item The constants of the pole bracket times $H(y,y)$: $\int_{y'}e^{-ik\cdot(y'-y)\frac1\xi}=(2\pi)^2\xi^2\delta^{(2)}(k)=0$ for $k\neq0$.
\item $-\frac{d_\perp}{2}\log\frac{(y'-y)^2\mu^2}{4}$ times $H(y,y)$: $\int_{y'}e^{-ik\cdot(y'-y)\frac1\xi}\log\frac{(y'-y)^2\mu^2}{4}=-\frac{4\pi\xi^2}{k^2}$ ($\mu$ drops out for $k\neq0$). With the measure $\frac{d\xi}{\xi^2}$, the three terms of the bracket give $\int d\xi\,\xi(1-\xi)$, $\int d\xi\,\xi$ (the $\xi$ of Eq.~\eqref{eq:I2xi}), $\int d\xi\frac{\xi}{1-\xi}\cdot$(subtracted): finite; and
\begin{equation}
\int_{\frac{k^+}{\vee}}^1d\xi\,\frac{1-\xi}{\xi}=\log\frac{\vee}{k^+}-1+\frac{k^+}{\vee}.
\end{equation}
\item $H(y,y')-H(y,y)$ vanishes at $y'=y$: its transform falls faster, finite.
\end{itemize}
So the only $\log\frac{\vee}{k^+}$ of line 1 besides (ii) is the $\frac{1-\xi}{\xi}$ term with $\log(y'-y)^2$ at $y'=y$. With $\frac{4\pi}{k^2}=\frac{1}{\pi}\int_{w,w'}e^{-ik\cdot(w'-w)}\Kc(w-z)\cdot\Kc(w'-z)$ and $y\to z$, its $\log\frac{\vee}{k^+}$ term is, exactly (identity (I4)),
@@H4@@
your \texttt{BFKL.tex} line 2 of this row. Its rest, $\big(-1+\frac{k^+}{\vee}\big)H_4\big|_{\log\to1}$, stays in the remaining terms. ($d_\perp=2$ here: no pole multiplies it.)

\paragraph{Step 4: the row.} Group II Row IV $=\log\frac{k^+}{\Lambda}$ (your line 2, JIMWLK) $+$ (i) $+$ (ii) $+\log\frac{\vee}{k^+}H_4$ $+$ a rest that is finite for $\vee\to\infty$; (i) $+$ (ii) $=0$.

\section*{Part B. The coefficient in front of $\log\mu^2$: running coupling and DGLAP}
\section{The five sources of \texorpdfstring{$\mu$}{mu}}\label{sec:mu}
As your advisor says, $\mu$ enters through A3 and the shockwave through the loop, each with its c.c., and through the real term: five sources. Each comes in the combination $\hat E_r=\frac{2}{\epsilon}+\log\frac{r^2\mu^2}{c_0^2}$ ($c_0=2e^{-\gamma_E}$), so the coefficient of $\log\mu^2$ is the coefficient of $\frac{2}{\epsilon}$. In units of
\begin{equation*}
{\rm LO}\equiv\frac{1}{(2\pi)^{3}}\frac{g^{2}}{2\pi^{2}}\frac{1}{k^{+}}\int_{x,y,x',y'}e^{-ik\cdot(y'-y)}\frac{(x-y)\cdot(x'-y')}{(x-y)^{2}(x'-y')^{2}}\left[U^{e'c}(y')-U^{e'c}(x')\right]\left[U^{ec}(y)-U^{ce}(x)\right]\rho^{e'}(x')\rho^{e}(x)
\end{equation*}
($\alpha_s=\frac{g^2}{4\pi}$; $L\equiv\log\frac{k^+}{\Lambda}$, $l_1\equiv\log(1-\frac{k^+}{\vee})$, $X=x-y$, $\Delta=y'-y$):
\paragraph{A3, your Group I Eq.~(1.2).} Prefactor $\frac{1}{(2\pi)^3}\frac{g^4N_c}{32\pi^4}\frac{1}{k^+}=\frac{\alpha_sN_c}{4\pi}\times$(LO prefactor). With $w'\to y'$, $w\to y$ the colour factor is $-$(that of LO) and the bracket is $-\hat E_X$:
\begin{equation}
\text{A3}+\text{c.c.}=\frac{\alpha_sN_c}{4\pi}\,{\rm LO}\Big[(\hat E_X+\hat E_{X'})\Big(\frac{11}{3}-2L+l_1\Big)\Big]+\dots\quad\Rightarrow\quad\log\mu^2:\ \frac{\alpha_s}{\pi}\,\frac{N_c}{2}\Big(\frac{11}{3}-2L+l_1\Big){\rm LO}.
\label{eq:muA3}
\end{equation}
\paragraph{Shockwave through the loop, your Group I Eq.~(1.32).} Prefactor $\frac{1}{(2\pi)^3}\frac{g^4N_c}{16\pi^4}\frac{1}{k^+}=\frac{\alpha_sN_c}{2\pi}\times$(LO prefactor). With $y\to x$, $z\to w$ (dummy variables) the phase, kernel and colour factor are those of A3 with the opposite sign of the colour factor, and the bracket is $+\hat E$:
\begin{equation}
\text{SW}+\text{c.c.}=\frac{\alpha_sN_c}{2\pi}\,{\rm LO}\Big[(\hat E_X+\hat E_{X'})\Big(-\frac{11}{6}+2L+2\log(1-\tfrac{\Lambda}{k^+})\Big)\Big]\quad\Rightarrow\quad\log\mu^2:\ \frac{\alpha_s}{\pi}N_c\Big(-\frac{11}{6}+2L\Big){\rm LO}.
\label{eq:muSW}
\end{equation}
\paragraph{Real, your Group II Eq.~(1.21).} Prefactor $\frac{g^2}{4\pi^2}N_c\times$(LO prefactor)$=\frac{\alpha_s}{\pi}N_c$. The bracket $\big[-\frac{d_\perp}{\epsilon}-d_\perp\gamma_E-\frac{d_\perp}{2}\log\frac{(y'-y)^2\mu^2}{4}-1\big]=-\hat E_\Delta+O(\epsilon)$ ($d_\perp=2-\epsilon$). Line 1 with your $\frac{\xi}{(1-\xi)_+}$ is, by Eq.~\eqref{eq:I2xi} and the $[0,1]$ plus distribution restricted to $[\frac{k^+}{\vee},1]$ ($\int_a^1\frac{d\xi\,h}{(1-\xi)_+}=\int_a^1\frac{h(\xi)-h(1)}{1-\xi}d\xi+h(1)\log(1-a)$, the $\overline{MS}$ one), $P^{(r)}_{gg,+}\otimes{\rm LO}-2N_c\big(l_1+\log\frac{\vee}{k^+}\big){\rm LO}$:
\begin{align}
\text{line 1}&\Rightarrow\ \log\mu^2:\ -\frac{\alpha_s}{2\pi}\Big[\int_{k^+/\vee}^1\frac{d\xi}{\xi^2}P^{(r)}_{gg,+}(\xi)\,{\rm LO}_\xi-2N_c\Big(l_1+\log\frac{\vee}{k^+}\Big){\rm LO}\Big],\label{eq:muR1}\\
\text{line 2}&\Rightarrow\ \log\mu^2:\ -\frac{\alpha_s}{2\pi}\,2N_c\Big(L+\log\frac{\vee}{k^+}+l_1\Big){\rm LO}.\label{eq:muR2}
\end{align}
$P^{(r)}_{gg,+}=2N_c\big[\frac{\xi}{(1-\xi)_+}+\frac{1-\xi}{\xi}+\xi(1-\xi)\big]$ is the real part of $P_{gg}$; it has no $\delta(1-\xi)$.
\paragraph{The coupling of LO.} LO is $\propto g^2$. In $\overline{MS}$, $g^2=g^2(\mu_R)\big[1-\frac{\alpha_s\beta_0}{4\pi}\hat E\big]$, $\beta_0=\frac{11}{3}N_c$:
\begin{equation}
\text{CT}=-\frac{\alpha_s\beta_0}{4\pi}\,\hat E\;{\rm LO}\quad\Rightarrow\quad\log\mu^2:\ -\frac{\alpha_s}{\pi}\,\frac{\beta_0}{4}\,{\rm LO}.
\label{eq:ct}
\end{equation}

\paragraph{Sum them up.} Coefficient of $\log\mu^2$ in units of $\frac{\alpha_s}{\pi}{\rm LO}$ (\texttt{mu\_check.py}; $P^*\equiv\int\frac{d\xi}{\xi^2}P^{(r)}_{gg,+}{\rm LO}_\xi/{\rm LO}$):
{\small
\begin{longtable}{@{}p{0.36\textwidth}cccccc@{}}\toprule
source & $\beta_0$ ($\frac{11}{3}$) & $\times L$ & $\times\log\frac{\vee}{k^+}$ & $\times l_1$ & $\times P^*$\\\midrule
@@MUTAB@@
\midrule
sum of the five sources & $0$ & $0$ & $0$ & $\frac{N_c}{2}$ & $-\frac12$\\
$+$ coupling of LO & $-\frac{11N_c}{12}=-\frac{\beta_0}{4}$ & $0$ & $0$ & $\frac{N_c}{2}$ & $-\frac12$\\\bottomrule
\end{longtable}}
\noindent So:
\begin{enumerate}
\item \textbf{$\beta_0$: the four sources of the single-gluon final state cancel exactly}: $\frac{N_c}{2}\cdot\frac{11}{3}-N_c\cdot\frac{11}{6}=0$. Not only $\log\mu^2$: the pole, $\gamma_E$ and $\log(x-w)^2$ cancel too, pointwise after the renaming (\texttt{dglap\_rc/checks.json}: $@@B0@@$). This is what your advisor expected.
\item \textbf{$\times\log\frac{k^+}{\Lambda}$ (evolution): three sources}, A3 : shockwave : real $=@@EV@@$, sum $0$. They cancel only with the real term. Pointwise the poles are in the ratio $1:-2:1$ and cancel to $@@PSUM@@$ (\texttt{ll\_poles.json}). These are the three sources of the evolution.
\item \textbf{$\times\log\frac{\vee}{k^+}$}: $+N_c$ (line 1, (ii)) and $-N_c$ (line 2, (i)) cancel: the BFKL terms have no $\mu$.
\item What is left: $-\frac{\alpha_s}{2\pi}\log\mu^2\big[\int P^{(r)}_{gg,+}\otimes{\rm LO}+\frac{\beta_0}{2}{\rm LO}\big]=-\frac{\alpha_s}{2\pi}\log\mu^2\,P_{gg}\otimes{\rm LO}$ (DGLAP, with $P_{gg}=P^{(r)}_{gg,+}+\frac{\beta_0}{2}\delta(1-\xi)$), and $+\frac{\alpha_s}{2\pi}N_c\,l_1\log\mu^2\,{\rm LO}$ (Q1, from A3 only; $l_1\to0$ for $\vee\to\infty$).
\end{enumerate}

\section{Why \texorpdfstring{$\beta_0$}{beta0} is not in the NLO terms, and whether Group II Row IV is running coupling}\label{sec:why}
\begin{itemize}
\item \textbf{The NLO diagrams carry no $\beta_0\log\mu^2$}, because of item 1: A3 ($+\frac{11}{3}$) and the shockwave through the loop ($2\times(-\frac{11}{6})$) cancel, and the real term has no $\beta_0$. So at this order nothing in your rows fixes a scale for the coupling; in particular no $\beta_0\log(x-w)^2$ survives.
\item \textbf{$\beta_0$ is still in the result: through the coupling of LO.} LO is $O(g^2)$, so $g^2$ must be renormalized; Eq.~\eqref{eq:ct} is the only $\beta_0$. Its pole is the $\frac{\beta_0}{2}\delta(1-\xi)$ part of $P_{gg}$, which the real term lacks. After $\overline{MS}$ renormalization of $g^2$ at $\mu_R$ and of $D_{h/g}$ at $\mu_F$, the $\beta_0$ term of the result is
\begin{equation}
\frac{\alpha_s\beta_0}{4\pi}\,\log\frac{\mu_R^2}{\mu_F^2}\;{\rm LO}\Big|_{\alpha_s(\mu_R)},
\end{equation}
and the result is $\mu_R$-independent: $\frac{d}{d\log\mu_R^2}{\rm LO}|_{\alpha_s(\mu_R)}=-\frac{\alpha_s\beta_0}{4\pi}{\rm LO}$. This is the running coupling: the coupling of LO runs, the NLO diagrams do not add to it.
\item \textbf{Group II Row IV is not running coupling.} Its $\log\mu^2$ (Eqs.~\eqref{eq:muR1}, \eqref{eq:muR2}) multiplies $\int\frac{d\xi}{\xi^2}P^{(r)}_{gg,+}(\xi){\rm LO}_\xi$, a convolution in $\xi$ with the \emph{real} splitting function. It has no $\delta(1-\xi)$ part and no $\beta_0$: the plus distribution integrates the $\frac{\xi}{1-\xi}$ singularity against $g(\xi)-g(1)$, and its $\delta(1-\xi)$ coefficient is zero in the $\overline{MS}$ convention, Eq.~\eqref{eq:muR1}. A running-coupling term would be $\beta_0\log\mu^2\times{\rm LO}$ (local in $\xi$). What Group II Row IV gives is the $\mu_F$ dependence of the fragmentation function (DGLAP), and $-2N_cL\log\mu^2$, which cancels against Group I (item 2).
\item With both: $-\frac{\alpha_s}{2\pi}\log\mu^2\,P_{gg}\otimes{\rm LO}$ is exactly what $\mu_F^2\frac{\partial}{\partial\mu_F^2}D_{h/g}=\frac{\alpha_s}{2\pi}P_{gg}\otimes D_{h/g}$ cancels. The coupling counterterm completes $P^{(r)}_{gg,+}$ to the full $P_{gg}$.
\end{itemize}
\end{document}
'''


def main():
    FOU = ML.fourier_line()
    import make_rest as MR
    rv = ML.AS.assemble(ML.BLK['G3 V'], MR.remaining_items('G3 V'), lhs=r'R_{\rm V}&=')
    r32 = disp(ML.AS.mark(G32[-1][0].strip(), r'\nc'))
    pc = list(ML.PC.values())
    ev = MU['evolution part, x L: A3 : shockwave : real']
    rep = {
        '@@EQ112@@': user_eq('Group_III', 104, 106, '(1.12)'),
        '@@EQ110@@': disp(ML.BLK['G3 III(2)'].terms[0], 'your (1.10), line 1'),
        '@@EQ121@@': user_eq('Group_II', 131, 133, '(1.21)'),
        '@@BV@@': disp(ML.HB['G3 V'], r'(B)$_{\rm V}$'), '@@BVI@@': disp(ML.HB['G3 VI'], r'(B)$_{\rm VI}$'),
        '@@BIII1@@': disp(ML.HB['G3 III(1)'], r'(B)$_{\rm III(1)}$'), '@@BIII2@@': disp(ML.HB['G3 III(2)'], r'(B)$_{\rm III(2)}$'),
        '@@RV@@': rv, '@@R32@@': r32, '@@H4@@': disp(H4, r'$\log\frac{\vee}{k^+}H_4$'),
        '@@ECID@@': tex_num(ML.EC_ID), '@@ECB@@': tex_num(ML.EC_B), '@@ECT0@@': tex_num(ML.EC_T0), '@@FOU@@': FOU,
        '@@MUTAB@@': mu_rows(), '@@EV@@': ':'.join(e.replace('*', '') for e in ev),
        '@@PSUM@@': tex_num(max(p['sum_relative'] for p in pc)),
        '@@B0@@': str(json.load(open(os.path.join(HERE, '..', 'dglap_rc', 'checks.json')))['beta0: G1 II part 1 + G1 V 2b UV'])}
    out = DOC
    for k, v in rep.items():
        out = out.replace(k, v)
    assert '@@' not in out, out[out.index('@@') - 80:out.index('@@') + 40]
    open(os.path.join(HERE, 'Steps_line1_and_mu.tex'), 'w').write(out)
    print('written Steps_line1_and_mu.tex')


if __name__ == '__main__':
    main()
