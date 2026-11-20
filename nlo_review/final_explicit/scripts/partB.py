# Part B: every finite term written explicitly (paper phase convention), from the author's rows.
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from texlib import *
o = []
A = o.append
NPH = [0]
def ph(s, n):
    s2, k = phase(s, n); NPH[0] += k; return s2
nl = r"\notag\\"+"\n"
A(r"""\section{Finite terms written explicitly}\label{sec:B}
\noindent Each term of the author's ``Three $\rho$ terms'' and ``Two $\rho$ terms'' lists is written out below in the convention of Sec.~\ref{sec:conv} (phase $e^{-ik\cdot(\text{unprimed}-\text{primed})}$, in red where it differs from the author's file). Where the author's term contains a plus prescription or a $p^{+}$ integral that can be done in closed form, the integral is done here and every piece is displayed: the terms with $\mathcal{E}$ and $\Psi$ are the \emph{soft remainders} of the plus prescriptions (they are the members of the pairs A--D of Sec.~\ref{sec:C}), the terms with $\widehat{\mathcal{E}}$, $\mathcal{J}$ and the $\xi$ integrals are the \emph{regular} terms, and the terms with $\pi\delta(s)-i\mathcal{P}\frac1s$ are the instantaneous terms. For every row it is stated whether it is the author's term unchanged (only the phase), the author's term with the integral done (equivalent, verified as stated in Sec.~\ref{sec:A}), or a corrected term (red). The complex conjugate is implied in Groups I and III.

\subsection*{Three $\rho$ terms}
\subsection{Group I (with c.c.)}
\subsubsection*{Row I} \textsc{NA}
\subsubsection*{Row II}
The author's term, with the kernel written in its antisymmetrised (integrable) form; both colour structures are antisymmetric under $x\leftrightarrow y$, i.e.\ under $\mathbf{p}\to\mathbf{q}-\mathbf{p}$, so this equals the author's form (Sec.~\ref{sec:A}, item A3):
\begin{align}
\frac{d^{3}N^{\rm I.II}_{2b}}{d^{3}k} & = -\frac{1}{(2\pi)^{3}}\,\frac{g^{4}}{64\pi^{6}}\,\frac{1}{k^{+}}\int_{w,w'}"""+PH_WW+r"""\int_{x,y,x'}"""+nl+r"""&\quad\times"""+R_I2b+nl+r"""&\quad\times\Big[U^{ab'}(x')-U^{ab'}(w')\Big]\rho^{b'}(x')\Big[f^{hbc}\,U^{ah}(w)\left\{\rho^{b}(x),\rho^{c}(y)\right\}-f^{abc}\,U^{bd}(x)\,U^{ce}(y)\left\{\rho^{d}(x),\rho^{e}(y)\right\}\Big]\,,\label{eq:B-I2b}\\
R^{i}(\mathbf{q},\mathbf{p}) & = -2q^{i}\log\left(\frac{\mathbf{p}^{2}}{\mathbf{q}^{2}}\right)+\frac{(2\mathbf{p}^{2}-2\,\mathbf{q}\cdot\mathbf{p}+\mathbf{q}^{2})\,q^{i}+\mathbf{q}^{2}p^{i}}{(\mathbf{q}-\mathbf{p})^{2}}\,\log\left[\frac{\mathbf{q}^{2}}{(\mathbf{q}-\mathbf{p})^{2}}\right]\quad\text{(the author's $R^{i}$)}\,.\notag
\end{align}
""")
# ---- Row III
A(r"""\subsubsection*{Row III}
The author's combined kernel (lines 1--4 of the next equation, unchanged up to the phase) is verified to be the sum of the separate terms $N_{1}$, $N_{2}$, $N_{3a}$ and of the first three terms of $N_{3b}$ of the earlier version (\texttt{Final\_terms\_checked}); without the third term of $N_{3b}$ it is the combined kernel $\mathcal{M}^{ij}$ of the verified files, called the 24-term kernel below (Sec.~\ref{sec:A}, item A1). Its last line contains the \emph{third term of $N_{3b}$}, $-\frac12\log\frac{(z-w)^{2}}{(x-w)^{2}}\,\mathcal{K}(y-z)\cdot\mathcal{K}(w-z)\,\mathcal{K}(x'-w')\cdot\mathcal{K}(x-w)$ (it is the term $\frac{1}{2(z-w)^{2}}\log\frac{(z-w)^{2}}{(x-w)^{2}}\mathcal{K}(x'-w')\cdot\mathcal{K}(x-w)\,\mathcal{K}(y-z)\cdot(z-w)$ of the author's last line). The \emph{last term of $N_{3b}$} is given in the verified form (red); the author's form of it is shown in orange directly after it and is \textbf{not} equivalent (\textbf{NOTE} A2). The $w$ integral has to be done before the $z$ integral (Sec.~\ref{sec:C4}).
\begin{align}
""")
blk, n = phase(L(172, 176), 1); NPH[0] += n
blk = blk.replace(r"\frac{d^{3}N_{\rm NLO}}{d^{3}k} &=", r"\frac{d^{3}N^{\rm I.III}}{d^{3}k} &=")
blk = blk.replace(r"\mathrm{arc\,cos}", r"\arccos")
A(blk+"\n")
A(r"""&{\color{red}-\frac{1}{(2\pi)^{3}}\,\frac{ig^{4}\,f^{abd}}{8\pi^{5}}\,\frac{1}{k^{+}}\int_{x,y,z,w,x',w'}"""+PH_WW+r"""\,\mathcal{E}\big(k\cdot\ph{(z-w)}\big)\,\mathcal{K}(x-z)\cdot\mathcal{K}(w-z)\,\mathcal{K}(x'-w')\cdot\mathcal{K}(y-w)}\notag\\
&{\color{red}\quad\times\Big[U^{ab'}(x')-U^{ab'}(w')\Big]\,U^{de}(x)\Big[U^{bc}(w)-U^{bc}(y)\Big]\rho^{b'}(x')\rho^{e}(x)\rho^{c}(y)}\label{eq:B-I3}
\end{align}
{\color{orange}Author's last term of $N_{3b}$ (phase changed), not used: $-\frac{1}{(2\pi)^{3}}\frac{ig^{4}f^{abd}}{16\pi^{5}}\frac{1}{k^{+}}\int_{x,y,z,w,x',w'}e^{-ik\cdot(w-w')}\log\frac{(z-w)^{2}}{(x-w)^{2}}\,\mathcal{K}(x'-w')\cdot\mathcal{K}(y-w)\,\mathcal{K}(x-z)\cdot\mathcal{K}(w-z)\times\big[U^{ab'}(x')-U^{ab'}(w')\big]U^{de}(x)\big[U^{bc}(w)-U^{bc}(y)\big]\rho^{b'}(x')\rho^{e}(x)\rho^{c}(y)$.}
Here $\mathcal{E}(s)=\int_{0}^{1}\frac{dt}{t}(e^{-ist}-1)$, Eq.~\eqref{eq:Ecal}. The red term is the recoil form of the last term of $N_{3b}$: with $\lambda=\Lambda/k^{+}$, the last term of $N_{3b}$ together with its $\log\frac{k^{+}}{\Lambda}$ partner of $N_{3a}$ is, exactly in $\lambda$ (substitution of the soft-gluon position for the hard-gluon position under the Fourier integral),
\begin{align}
-\int_{w,z}e^{-ik\cdot w}\,\mathcal{K}(x-z)\cdot\mathcal{K}(w-z)\,\Phi(w)\int_{\lambda}^{1-\lambda}\frac{dt}{t}\,e^{-it\,k\cdot(z-w)}=-\int_{w,z}e^{-ik\cdot w}\,\mathcal{K}(x-z)\cdot\mathcal{K}(w-z)\,\Phi(w)\Big[\log\frac{k^{+}}{\Lambda}+\mathcal{E}\big(k\cdot(z-w)\big)+O(\lambda)\Big],\notag
\end{align}
$\Phi(w)=\mathcal{K}(x'-w')\cdot\mathcal{K}(y-w)\,[U(x')-U(w')]^{ab'}\rho^{b'}(x')\,f^{abd}U^{de}(x)[U^{bc}(w)-U^{bc}(y)]\rho^{e}(x)\rho^{c}(y)$ (units $\frac{1}{(2\pi)^{3}}\frac{ig^{4}}{8\pi^{5}}\frac{1}{k^{+}}$); the $\log\frac{k^{+}}{\Lambda}$ term is in the JIMWLK list, the $\mathcal{E}$ term is the finite term (paper draft, Eq.~(N3bexact)).

\paragraph{Row III in the two prescriptions (author's instruction of 2026-11-20: ``Do it in two ways with $i\epsilon$ and PV. We are in pure glue'').} The vanishing energy denominators $p^{+}(x-z)^{2}-k^{+}(x-w)^{2}$ of $\mathbb{B}_{3}$ (the only vanishing energy denominators of the calculation; they enter only this row and its complex conjugate) are taken (i) as principal values (standing-wave dressing of the observed gluon): Eq.~\eqref{eq:B-I3} as written; (ii) with $+i\eta0$, $\eta=+1$ (out-state dressing; $\eta=-1$: in-state dressing of the final state; PV is the average). With (ii) the $i\pi\theta$ pieces of the separate terms $N_{1},N_{2},N_{3a}$ combine exactly (paper notes, ref3\_prescriptions Eq.~(pr:shift)) into
\begin{align}
\log\frac{(x-z)^{2}}{(x-w)^{2}}\;\longrightarrow\;\log\frac{(x-z)^{2}}{(x-w)^{2}}-i\eta\pi\qquad\text{in the first two lines of Eq.~\eqref{eq:B-I3} only}\notag
\end{align}
(the $\arccos$ term, $\log\frac{(z-w)^{2}}{(x-w)^{2}}$ and the last term are unchanged; the pole on $(x-z)^{2}=(x-w)^{2}$ cancels for every $\eta$). The bracket multiplying that logarithm is, exactly,
\begin{align}
&\bigg[\frac{\mathcal{K}(x'-w')\cdot\mathcal{K}(y-z)}{2\,(z-w)^{2}}+\mathcal{K}(x'-w')\cdot\mathcal{K}(z-w)\,\mathcal{K}(y-z)\cdot\mathcal{K}(x-z)\notag\\
&\qquad-\mathcal{K}(x'-w')\cdot\mathcal{K}(x-w)\,\mathcal{K}(y-z)\cdot\mathcal{K}(z-w)-\mathcal{K}(x'-w')\cdot\mathcal{K}(x-w)\,\mathcal{K}(y-z)\cdot\mathcal{K}(x-z)\bigg]\notag\\
&\quad=\frac{1}{2\,(z-w)^{2}}\Big[\cos2\phi\;\mathcal{K}(x'-w')\cdot\mathcal{K}(y-z)+\sin2\phi\;\mathcal{K}(x'-w')\times\mathcal{K}(y-z)\Big]\,,\label{eq:B-Bij}\\
&\cos2\phi=\frac{[(x-w)\cdot(x-z)]^{2}-[(x-w)\times(x-z)]^{2}}{(x-w)^{2}(x-z)^{2}}\,,\qquad\sin2\phi=\frac{2\,[(x-w)\cdot(x-z)]\,[(x-w)\times(x-z)]}{(x-w)^{2}(x-z)^{2}}\,,\qquad a\times b\equiv a^{1}b^{2}-a^{2}b^{1}\,,\notag
\end{align}
$\phi$ the angle from $x-w$ to $x-z$. Row III with the $i\epsilon$ prescription is therefore Eq.~\eqref{eq:B-I3} plus
\begin{align}
\frac{d^{3}\Delta N^{\rm I.III}}{d^{3}k}&="""+PF_DN+r"""\int_{w,w'}"""+PH_WW+r"""\int_{x,y,x'}\Bigg\{\int_{z}\bigg[\frac{\mathcal{K}(x'-w')\cdot\mathcal{K}(y-z)}{2\,(z-w)^{2}}+\mathcal{K}(x'-w')\cdot\mathcal{K}(z-w)\,\mathcal{K}(y-z)\cdot\mathcal{K}(x-z)\notag\\
&\qquad-\mathcal{K}(x'-w')\cdot\mathcal{K}(x-w)\,\mathcal{K}(y-z)\cdot\mathcal{K}(z-w)-\mathcal{K}(x'-w')\cdot\mathcal{K}(x-w)\,\mathcal{K}(y-z)\cdot\mathcal{K}(x-z)\bigg]\notag\\
&\qquad\times f^{abd}\Big[U^{ab'}(x')-U^{ab'}(w')\Big]\rho^{b'}(x')\,U^{de}(x)\Big[U^{bc}(z)-U^{bc}(y)\Big]\rho^{e}(x)\rho^{c}(y)\notag\\
&\qquad-\int_{z}\theta\big(\rho_{0}-|z-w|\big)\,\frac{\mathcal{K}(x'-w')\cdot\mathcal{K}(y-w)}{2\,(z-w)^{2}}\notag\\
&\qquad\times f^{abd}\Big[U^{ab'}(x')-U^{ab'}(w')\Big]\rho^{b'}(x')\,U^{de}(x)\Big[U^{bc}(w)-U^{bc}(y)\Big]\rho^{e}(x)\rho^{c}(y)\notag\\
&\qquad+\pi\Big[\log\frac{\rho_{0}}{|x-w|}{\color{gray}\;+\;\log\frac{k^{+}}{\Lambda}}\Big]\mathcal{K}(x'-w')\cdot\mathcal{K}(y-w)\notag\\
&\qquad\times f^{abd}\Big[U^{ab'}(x')-U^{ab'}(w')\Big]\rho^{b'}(x')\,U^{de}(x)\Big[U^{bc}(w)-U^{bc}(y)\Big]\rho^{e}(x)\rho^{c}(y)\Bigg\}+{\rm c.c.}\label{eq:B-DN}
\end{align}
($\eta=+1$; $\rho_{0}>0$ arbitrary, the result does not depend on it; prefactor $\frac{ig^{4}}{8\pi^{5}}\times(-i\eta\pi)=\frac{\eta g^{4}}{8\pi^{4}}$). The first term alone has a non-integrable $1/(z-w)^{2}$ at $z\to w$ (where $\phi\to0$ and the bracket tends to $\frac{\mathcal{K}(x'-w')\cdot\mathcal{K}(y-w)}{2(z-w)^{2}}$); it is cut off by the window $|p^{+}-k^{+}|>\Lambda$ of $\mathbb{B}_{3}$, which gives $\int_{|z-w|<\rho_{0}}\frac{d^{2}z}{(z-w)^{2}}\big|_{\rm window}=2\pi\log\frac{\rho_{0}k^{+}}{\Lambda|x-w|}$; this is written as the subtraction in the second line plus the third line. {\color{gray}The term with $\log\frac{k^{+}}{\Lambda}$ (gray) is a forward Coulomb phase and cancels exactly against its complex conjugate (Sec.~\ref{sec:C7})}; what remains is finite and integrable. Properties (ref3\_prescriptions, a.2, a.7): $\int d^{3}k\,\frac{d^{3}\Delta N}{d^{3}k}=0$ (the on-shell part of $\mathbb{B}_{3}$ is anti-Hermitian; the total number of gluons is unchanged); $\Delta N$ is (real kernel)$\,+\,$c.c., while the principal-value terms of the row are $i\times$(real kernel)$\,+\,$c.c., so $\Delta N$ does not cancel against its complex conjugate and its parity in $k$ is the opposite of that of the principal-value terms: its symmetrised three-$\rho$ part is \emph{even} in $k$, its commutator (two-$\rho$) part \emph{odd} (Secs.~\ref{sec:D}, \ref{sec:E}). $\Delta N$ is linear in $\eta$. No other row has a vanishing energy denominator, so nothing else cancels $\Delta N$. Its size for a realistic configuration has not been evaluated.

\subsubsection*{Row IV}
The author's term (unchanged up to the phase); it cancels Group III Row I exactly (Sec.~\ref{sec:C1}):
\begin{align}
""")
blk, n = phase(L(182, 182), 1); NPH[0] += n
blk = blk.replace(r"\frac{d^{3}N_{\rm NLO}}{d^{3}k} & =", r"\frac{d^{3}N^{\rm I.IV}}{d^{3}k} & =")
A(blk+r"\label{eq:B-I4}"+"\n"+r"\end{align}"+"\n")
# ---- Row V first part
A(r"""\subsubsection*{Row V}
The author's term with the $p^{+}$ integrals done: with $t=p^{+}/k^{+}$, $\int_{0}^{k^{+}}\frac{dp^{+}}{2\pi}\frac{1}{k^{+2}}e^{-its}=\frac{1}{2\pi k^{+}}\mathcal{J}(s)$, $\int_{0}^{k^{+}}\frac{dp^{+}}{2\pi}\frac{1}{k^{+2}}\big[\frac{k^{+}}{p^{+}}\big]_{+}e^{-its}=\frac{1}{2\pi k^{+}}\mathcal{E}(s)$ and, for the second term, $e^{-ik\cdot(x-w')}=e^{-ik\cdot(z-w')}e^{-is}$, $\int_{0}^{k^{+}}\frac{dq^{+}}{2\pi}\big[\frac{1}{q^{+}}\big]_{+}e^{-ik\cdot(z-x)q^{+}/k^{+}}=\frac{1}{2\pi}\mathcal{E}(-s)$:
\begin{align}
\frac{d^{3}N^{\rm I.V}_{1}}{d^{3}k} & =\frac{1}{(2\pi)^{3}}\,\frac{ig^{4}}{32\pi^{5}}\,\frac{1}{k^{+}}\int_{x,z,y,y',x',w'}"""+PH_ZW+nl+r"""&\quad\times\Big["""+KV1_J+nl+r"""&\qquad"""+KV1_E+r"""\Big]"""+nl+r"""&\quad\times\Big[U^{ab'}(x')-U^{ab'}(w')\Big]\rho^{b'}(x')\bigg[f^{adc}\,U^{db}(x)\,U^{ce}(z)\,\{\rho^{e}(y),\rho^{b}(y')\}-f^{abc}\,U^{cd}(y)\,U^{be}(y')\,\{\rho^{d}(y),\rho^{e}(y')\}\bigg]\,,\qquad """+S_G1V+r""".\label{eq:B-I5}
\end{align}
The first term is regular, the second and third are the soft remainders of $\big[\frac{k^{+}}{p^{+}}\big]_{+}$ and $\big[\frac{1}{q^{+}}\big]_{+}$ (soft gluon at $x$ and at $z$, respectively; pair B).
""")
# ---- Group II Row I first part
A(r"""\subsection{Group II}
\subsubsection*{Row I (first part)}
The author's terms (unchanged up to the phase); Group II has no complex conjugate:
\begin{align}
""")
for a, b, lab in [(196, 197, 'eq:B-II1a'), (200, 201, 'eq:B-II1b')]:
    blk, n = phase(L(a, b), 1); NPH[0] += n
    blk = blk.replace(r"N_{1a\,\rm NLO}", r"N^{\rm II.I}_{1a}").replace(r"N_{1b\,\rm NLO}", r"N^{\rm II.I}_{1b}")
    A(blk+r"\label{%s}" % lab+("\\\\\n" if lab.endswith('1a') else "\n"))
A(r"""\end{align}
\subsubsection*{Rows II--IV} \textsc{NA}
\subsection{Group III (with c.c.)}
\subsubsection*{Row I}
The author's term (unchanged up to the phase); it cancels Group I Row IV exactly (Sec.~\ref{sec:C1}):
\begin{align}
""")
blk, n = phase(L(212, 214), 1); NPH[0] += n
blk = blk.replace(r"\frac{d^{3}N_{\rm NLO}}{d^{3}k} & =", r"\frac{d^{3}N^{\rm III.I}}{d^{3}k} & =")
A(blk+r"\label{eq:B-III1}"+"\n"+r"\end{align}"+"\n")
A(r"""\subsubsection*{Row II}
The author's term (unchanged up to the phase):
\begin{align}
""")
blk, n = phase(L(218, 219), 1); NPH[0] += n
blk = blk.replace(r"\frac{d^{3}N_{\rm NLO}}{d^{3}k} & =", r"\frac{d^{3}N^{\rm III.II}}{d^{3}k} & =")
A(blk+r"\label{eq:B-III2}"+"\n"+r"\end{align}"+"\n")
# ---- Group III Row III first part
III_intro = r"""With $t=p^{+}/k^{+}$ and $\vee\to\infty$ (Sec.~\ref{sec:A}, item A4): $\int_{0}^{\infty}\frac{dp^{+}}{2\pi}\frac{1}{p^{+}+k^{+}}e^{-its}=\frac{1}{2\pi}\widehat{\mathcal{E}}(s)$, $\int_{0}^{\infty}\frac{dp^{+}}{2\pi}\big[\frac{k^{+}}{p^{+}(p^{+}+k^{+})}\big]_{+}e^{-its}=\frac{1}{2\pi}\Psi(s)$, $\int_{0}^{\infty}\frac{dp^{+}}{2\pi}\frac{1}{k^{+2}}e^{-its}=\frac{1}{2\pi k^{+}}\big[\pi\delta(s)-i\mathcal{P}\frac{1}{s}\big]$, Eqs.~\eqref{eq:Ehat}--\eqref{eq:distr}."""
def g3row(name, pref, colL, colLz, colR, lab, intro, labname):
    return (r"""\subsubsection*{%s}
%s
\begin{align}
\frac{d^{3}N%s}{d^{3}k} & =\frac{1}{(2\pi)^{3}}\,%s\,\frac{1}{k^{+}}\int_{w,x',y',z,y,x}""" % (name, intro, lab, pref)+PH_WY+nl+
            r"""&\quad\times\bigg\{\widehat{\mathcal{E}}(s)\Big[\mathcal{K}(y'-z)\cdot\mathcal{K}(x'-y')\,\mathcal{K}(y-w)\cdot\mathcal{K}(x-z)-\mathcal{K}(y'-z)\cdot\mathcal{K}(x-z)\,\mathcal{K}(x'-y')\cdot\mathcal{K}(y-w)\Big]"""+nl+
            r"""&\qquad\qquad\times"""+colL+nl+
            r"""&\qquad-\Psi(s)\,\mathcal{K}(y'-z)\cdot\mathcal{K}(x-z)\,\mathcal{K}(x'-y')\cdot\mathcal{K}(y-w)"""+colL+nl+
            r"""&\qquad-"""+DIST+r"""\mathcal{K}(y'-z)\cdot\mathcal{K}(y-w)\Big[\mathcal{K}(x-z)\cdot\mathcal{K}(x'-y')"""+colL+nl+
            r"""&\qquad\qquad-\mathcal{K}(x-z)\cdot\mathcal{K}(x'-z)"""+colLz+r"""\Big]\bigg\}"""+nl+
            r"""&\quad\times """+colR+r""",\qquad """+S_G3+r""".\label{%s}
\end{align}
""" % labname)
A(g3row('Row III (first part)', r"\frac{ig^{4}\,f^{c'd'a}}{16\pi^{5}}",
        r"\,U^{d'b}(z)\Big[U^{c'e'}(y')-U^{c'e'}(x')\Big]\rho^{e'}(x')", r"\,U^{d'b}(z)\Big[U^{c'e'}(z)-U^{c'e'}(x')\Big]\rho^{e'}(x')",
        r"\Big[U^{ac}(w)\,\{\rho^{c}(y),\rho^{b}(x)\}-U^{cb}(z)\,U^{ad}(y)\,U^{ce}(x)\,\{\rho^{d}(y),\rho^{e}(x)\}\Big]", r"^{\rm III.III}_{1}",
        III_intro+r""" The first line is regular, the second is the soft remainder of $\big[\frac{k^{+}}{p^{+}(p^{+}+k^{+})}\big]_{+}$ (pair B), the third is the instantaneous term $-\frac{1}{k^{+}}\delta_{im}\delta_{jk'}$ with its subtraction at $y'=z$; the two terms of the bracket of the third line belong to one $y'$ integral, done before the $p^{+}$ integral.""", 'eq:B-III3one'))
A(r"""\subsubsection*{Row IV} \textsc{NA}
""")
A(g3row('Row V', r"\frac{ig^{4}\,f^{c'ca}}{8\pi^{5}}",
        r"\,\Big[U^{c'e'}(y')-U^{c'e'}(x')\Big]\rho^{e'}(x')", r"\,\Big[U^{c'e'}(z)-U^{c'e'}(x')\Big]\rho^{e'}(x')",
        r"U^{ce}(x)\Big[U^{ad}(y)-U^{ad}(w)\Big]\rho^{e}(x)\rho^{d}(y)", r"^{\rm III.V}",
        r"""As Row III (first part); the second line is the soft remainder of pair A.""", 'eq:B-III5'))
A(g3row('Row VI', r"\frac{ig^{4}\,f^{c'd'a}}{8\pi^{5}}",
        r"\,\Big[U^{c'e'}(y')-U^{c'e'}(x')\Big]\rho^{e'}(x')", r"\,\Big[U^{c'e'}(z)-U^{c'e'}(x')\Big]\rho^{e'}(x')",
        r"U^{ad}(y)\Big[U^{d'e}(x)-U^{d'e}(z)\Big]\rho^{d}(y)\rho^{e}(x)", r"^{\rm III.VI}",
        r"""As Row III (first part); the second line is the soft remainder of pair D. The colour structure is the author's (which equals the earlier form $[U^{c'e'}(y')U^{d'b}(z)-U^{d'b}(z)U^{c'e'}(x')]\rho^{e'}(x')[U^{cb}(z)U^{ad}(y)U^{ce}(x)\rho^{d}(y)\rho^{e}(x)-U^{ad}(y)\rho^{d}(y)\rho^{b}(x)]$ by $U^{d'b}(z)U^{cb}(z)=\delta^{d'c}$; checked at the operator level, Sec.~\ref{sec:F}).""", 'eq:B-III6'))
# ================= Two rho =================
A(r"""\subsection*{Two $\rho$ terms}
\subsection{Group I (with c.c.)}
\subsubsection*{Row I} \textsc{NA}
\subsubsection*{Row II}
The author's term (unchanged up to the phase). The virtual $\frac{11}{3}$ terms of this row are not in the author's file because they cancel Row V exactly; they are displayed in Sec.~\ref{sec:C2}.
\begin{align}
""")
blk, n = phase(L(246, 247), 1); NPH[0] += n
blk = blk.replace(r"\frac{d^{3}N_{1\,\rm NLO}}{d^{3}k} & =", r"\frac{d^{3}N^{\rm I.II}_{1}}{d^{3}k} & =")
A(blk+r"\label{eq:B-I2a}"+"\n"+r"\end{align}"+"\n")
A(r"""\subsubsection*{Rows III, IV} \textsc{NA}
\subsubsection*{Row V}
The author's two parts $N_{2a}$ and $N_{2b\,\rm rest}$ have the same kernel; their colour factors add up to $-f^{adc}f^{fbe}U^{db}(x)U^{ce}(z)+N_{c}U^{af}(y)$. With $\xi=p^{+}/k^{+}$, the plus prescriptions acting on the whole function $e^{-i\xi s}/[\xi(y-x)^{2}+(1-\xi)(y-z)^{2}]$ (Eq.~\eqref{eq:plus}), and $\int_{0}^{1}d\xi\,[\xi(1-\xi)-2]=-\frac{11}{6}$ for the last term of $N_{2b\,\rm rest}$:
\begin{align}
\frac{d^{3}N^{\rm I.V}_{2}}{d^{3}k} & =\frac{1}{(2\pi)^{3}}\,\frac{g^{4}}{16\pi^{5}}\,\frac{1}{k^{+}}\int_{x,y,z,x',w'}"""+PH_ZW+nl+r"""&\quad\times\Bigg\{\int_{0}^{1}d\xi\,\frac{e^{-i\xi s}}{\xi\,(y-x)^{2}+(1-\xi)\,(y-z)^{2}}\bigg[\xi(1-\xi)\,\frac{(y-x)^{2}-(y-z)^{2}}{(x-z)^{2}}\,\mathcal{K}(x'-w')\cdot\mathcal{K}(z-x)"""+nl+r"""&\qquad+(1-\xi)\,\mathcal{K}(x'-w')\cdot\mathcal{K}(z-x)\Big[(y-z)\cdot\mathcal{K}(x-z)+\frac{1}{2}(y-z)\cdot\mathcal{K}(y-x)-\frac{(y-x)^{2}-(y-z)^{2}}{2\,(x-z)^{2}}\Big]"""+nl+r"""&\qquad+\xi\,\mathcal{K}(x'-w')\cdot\mathcal{K}(z-x)\Big[(y-x)\cdot\mathcal{K}(x-z)-\frac{1}{2}(y-x)\cdot\mathcal{K}(y-z)-\frac{(y-x)^{2}-(y-z)^{2}}{2\,(x-z)^{2}}\Big]"""+nl+r"""&\qquad+\frac{(y-x)^{2}}{(y-z)^{2}}\Big[-\frac{(y-z)\cdot\mathcal{K}(x'-w')}{(x-z)^{2}}+\frac{1}{2}\,(y-z)\cdot\mathcal{K}(x'-w')\,\mathcal{K}(z-x)\cdot\mathcal{K}(y-x)\Big]"""+nl+r"""&\qquad+\frac{(y-z)^{2}}{(y-x)^{2}}\Big[-\frac{(y-x)\cdot\mathcal{K}(x'-w')}{(x-z)^{2}}-\frac{1}{2}\,(y-x)\cdot\mathcal{K}(x'-w')\,\mathcal{K}(z-x)\cdot\mathcal{K}(y-z)\Big]"""+nl+r"""&\qquad-(y-x)\cdot\mathcal{K}(z-x)\,\mathcal{K}(x-z)\cdot\mathcal{K}(x'-w')+\frac{1}{2}\,(y-x)\cdot\mathcal{K}(z-x)\,\mathcal{K}(y-z)\cdot\mathcal{K}(x'-w')"""+nl+r"""&\qquad-(y-z)\cdot\mathcal{K}(z-x)\,\mathcal{K}(x-z)\cdot\mathcal{K}(x'-w')-\frac{1}{2}\,\mathcal{K}(y-x)\cdot\mathcal{K}(x'-w')\,(y-z)\cdot\mathcal{K}(z-x)\bigg]"""+nl+r"""&\qquad-\frac{\mathcal{E}(s)}{(y-z)^{2}}\Big[-\frac{(y-z)\cdot\mathcal{K}(x'-w')}{(x-z)^{2}}+\frac{1}{2}\,(y-z)\cdot\mathcal{K}(x'-w')\,\mathcal{K}(z-x)\cdot\mathcal{K}(y-x)\Big]"""+nl+r"""&\qquad-\frac{e^{-is}\,\mathcal{E}(-s)}{(y-x)^{2}}\Big[-\frac{(y-x)\cdot\mathcal{K}(x'-w')}{(x-z)^{2}}-\frac{1}{2}\,(y-x)\cdot\mathcal{K}(x'-w')\,\mathcal{K}(z-x)\cdot\mathcal{K}(y-z)\Big]\Bigg\}"""+nl+r"""&\qquad\times\Big[U^{ab'}(x')-U^{ab'}(w')\Big]\Big[-f^{adc}\,f^{fbe}\,U^{db}(x)\,U^{ce}(z)+N_{c}\,U^{af}(y)\Big]\rho^{b'}(x')\,\rho^{f}(y)"""+nl+r"""&+\frac{1}{(2\pi)^{3}}\,\frac{g^{4}N_{c}}{16\pi^{5}}\,\frac{1}{k^{+}}\int_{x,y,z,x',w'}"""+PH_ZW+r"""\,e^{-is}\,\frac{11}{6}\,\mathcal{K}(x'-w')\cdot\mathcal{K}(y-x)\,\frac{(y-x)\cdot\mathcal{K}(y-z)}{(x-z)^{2}}"""+nl+r"""&\qquad\times\Big[U^{ab'}(x')-U^{ab'}(w')\Big]\Big[U^{af}(y)-U^{af}(x)\Big]\rho^{b'}(x')\,\rho^{f}(y)\,,\label{eq:B-I5two}
\end{align}
$"""+S_G1V+r"""$. The $\xi$ integral (lines 1--7) is the regular term; lines 8 and 9 are the soft remainders of $\big[\frac{k^{+}}{p^{+}}\big]_{+}$ (soft gluon at $x$) and of $\big[\frac{k^{+}}{k^{+}-p^{+}}\big]_{+}$ (soft gluon at $z$), pair C; the last line is the collinear subtraction of $N_{2b\,\rm rest}$ ($e^{-ik\cdot(z-w')}e^{-is}=e^{-ik\cdot(x-w')}$). In the notation of the paper draft, lines 1--7 are $\mathcal{K}^{i}(x'-w')\mathcal{K}^{m}(z-x)\mathcal{F}^{im}_{\rm reg}$ and lines 8, 9 are $-\mathcal{E}(s)\mathcal{K}^{i}\mathcal{K}^{m}B_{2}^{im}/(y-z)^{2}-e^{-is}\mathcal{E}(-s)\mathcal{K}^{i}\mathcal{K}^{m}B_{1}^{im}/(y-x)^{2}$, written out here.
""")
A(r"""\subsection{Group II}
\subsubsection*{Row I (second part)}
The author's terms (unchanged up to the phase). They have to be kept together: on the surface $(x'-w')^{2}(x-z)^{2}=(x-w)^{2}(x'-z)^{2}$ each of $N_{2a}$, $N_{2b}$, $N_{2c}$, $N_{2d}$ has a simple pole; the poles cancel in $N_{2a}+N_{2b}$ and in $N_{2c}+N_{2d}$ (Sec.~\ref{sec:A}, item A6).
\begin{align}
""")
for a, b, nm, lab in [(281, 282, '2a', 'eq:B-II2a'), (285, 287, '2b', 'eq:B-II2b'), (290, 294, '2c', 'eq:B-II2c'), (297, 303, '2d', 'eq:B-II2d'), (306, 307, '2e', 'eq:B-II2e')]:
    blk, n = phase(L(a, b), 1); NPH[0] += n
    blk = blk.replace(r"N_{%s\,\rm NLO}" % nm, r"N^{\rm II.I}_{%s}" % nm)
    blk = blk.replace(r"+\frac{(x-w)\cdot\mathcal{K}(x-z)\,\big[(x^{\prime}-z)^{2}-(x^{\prime}-w^{\prime})^{2}\big]}{4\,(z-w^{\prime})^{2}}\notag", r"\notag\\ &\qquad+\frac{(x-w)\cdot\mathcal{K}(x-z)\,\big[(x^{\prime}-z)^{2}-(x^{\prime}-w^{\prime})^{2}\big]}{4\,(z-w^{\prime})^{2}}\notag", 1)
    blk = blk.replace(r"+\frac{(x^{\prime}-w^{\prime})\cdot(z-w^{\prime})\,\big[(x-z)^{2}-(x-w)^{2}\big]}{2\,(z-w)^{2}\,(z-w^{\prime})^{2}}\notag", r"\notag\\ &\qquad+\frac{(x^{\prime}-w^{\prime})\cdot(z-w^{\prime})\,\big[(x-z)^{2}-(x-w)^{2}\big]}{2\,(z-w)^{2}\,(z-w^{\prime})^{2}}\notag", 1)
    A(blk+r"\label{%s}" % lab+("\\\\\n" if nm != '2e' else "\n"))
A(r"""\end{align}
\subsubsection*{Rows II, III} \textsc{NA}
\subsubsection*{Row IV}
The author's term (unchanged up to the phase; equal to the $\xi$ form $\frac{1}{(2\pi)^{3}}\frac{g^{4}N_{c}}{8\pi^{4}}\frac{1}{k^{+}}\int_{k^{+}/\vee}^{1}d\xi\,\frac{1-\xi}{\xi}\int e^{-ik\cdot(y-y')/\xi}[\dots]$ of the verified files by $\xi=k^{+}/(p^{+}+k^{+})$), and the remnant term of the verified files, which is missing in the author's file (red, \textbf{NOTE} A5):
\begin{align}
""")
blk, n = phase(L(315, 315), 2); NPH[0] += n
blk = blk.replace(r"\frac{d^{3}N}{d^{3}k}= &", r"\frac{d^{3}N^{\rm II.IV}_{\rm fin}}{d^{3}k}= &")
A(blk+r"\label{eq:B-II4}\\"+"\n")
A(r"""{\color{red}\frac{d^{3}N^{\rm II.IV}_{\rm rem}}{d^{3}k}}=&{\color{red}-\frac{1}{(2\pi)^{3}}\,\frac{g^{4}\,N_{c}}{2\pi^{3}}\,\frac{1}{k^{2}}\,\frac{1}{k^{+}}\int_{x,x',z}\mathcal{K}(x-z)\cdot\mathcal{K}(x'-z)\Big[U^{ce'}(z)-U^{ce'}(x')\Big]\Big[U^{ce}(z)-U^{ce}(x)\Big]\rho^{e'}(x')\rho^{e}(x)}\label{eq:B-II4rem}
\end{align}
The DGLAP and running-coupling terms of this row are those of Secs.~3 and 4 (the author's, phase changed). The remnant goes with the author's convention for the subtraction in the DGLAP term ($[\frac{1-\xi}{\xi}]_{+}$ subtracting $\frac{1-\xi}{\xi}\times$ the value at $\xi\to0$); it is the $-1$ of $\int_{k^{+}/\vee}^{1}d\xi\,\frac{1-\xi}{\xi}=\log\frac{\vee}{k^{+}}-1$ (the logarithm is in the BFKL list).
\subsection{Group III (with c.c.)}
\subsubsection*{Rows I, II} \textsc{NA}
\subsubsection*{Row III (second part)}
The author's term with the plus prescription evaluated (the subtraction includes the phase and the energy denominator, Eq.~\eqref{eq:plus}); the coefficient of the last bracket of the first line group is $-\big[\frac{k^{+}}{p^{+}(p^{+}+k^{+})}\big(1-\frac{p^{+}(x-z)^{2}+k^{+}(x-w)^{2}}{k^{+}(x-w)^{2}}\big)+\frac{1}{p^{+}+k^{+}}\big]=\frac{(x-z)^{2}-(x-w)^{2}}{(p^{+}+k^{+})\,(x-w)^{2}}$, and the subtraction at $p^{+}=0$ gives the $\Psi$ term (pair C). The last term is the author's (phase changed).
\begin{align}
\frac{d^{3}N^{\rm III.III}_{2}}{d^{3}k} & =\frac{1}{(2\pi)^{3}}\,\frac{g^{4}}{4\pi^{4}}\int_{0}^{\infty}\frac{dp^{+}}{2\pi}\int_{w,x',y',z,x}"""+PH_WY+r"""\,e^{-i\frac{p^{+}}{k^{+}}s}\,\frac{1}{p^{+}\,(x-z)^{2}+k^{+}\,(x-w)^{2}}"""+nl+r"""&\quad\times\Bigg\{\mathcal{K}(y'-z)\cdot\mathcal{K}(x'-y')\bigg[\frac{p^{+}\big[(x-z)^{2}-(x-w)^{2}\big]}{(z-w)^{2}\,(p^{+}+k^{+})^{2}}"""+nl+r"""&\qquad+\frac{p^{+}}{k^{+}(p^{+}+k^{+})}\Big((x-z)\cdot\mathcal{K}(z-w)-\frac{1}{2}(x-z)\cdot\mathcal{K}(x-w)-\frac{(x-z)^{2}-(x-w)^{2}}{2\,(z-w)^{2}}\Big)"""+nl+r"""&\qquad+\frac{1}{p^{+}+k^{+}}\Big((x-w)\cdot\mathcal{K}(z-w)+\frac{1}{2}(x-w)\cdot\mathcal{K}(x-z)-\frac{(x-z)^{2}-(x-w)^{2}}{2\,(z-w)^{2}}\Big)\bigg]"""+nl+r"""&\qquad+\frac{(x-z)^{2}-(x-w)^{2}}{(p^{+}+k^{+})\,(x-w)^{2}}\Big[(x-w)\cdot\mathcal{K}(x'-y')\,\mathcal{K}(z-w)\cdot\mathcal{K}(y'-z)+\frac{1}{2}\,\mathcal{K}(x-z)\cdot\mathcal{K}(y'-z)\,(x-w)\cdot\mathcal{K}(x'-y')\Big]"""+nl+r"""&\qquad-\frac{1}{k^{+}}\Big[(x-z)\cdot\mathcal{K}(y'-z)\,\mathcal{K}(z-w)\cdot\mathcal{K}(x'-y')-\frac{1}{2}(x-z)\cdot\mathcal{K}(y'-z)\,\mathcal{K}(x-w)\cdot\mathcal{K}(x'-y')"""+nl+r"""&\qquad\qquad+(x-w)\cdot\mathcal{K}(y'-z)\,\mathcal{K}(z-w)\cdot\mathcal{K}(x'-y')"""+nl+r"""&\qquad\qquad+\frac{1}{2}\,\mathcal{K}(x-z)\cdot\mathcal{K}(x'-y')\,(x-w)\cdot\mathcal{K}(y'-z)\Big]\Bigg\}"""+nl+r"""&\qquad\times\Big[f^{c'd'a}\,U^{c'e'}(y')\,U^{d'b}(z)\,\rho^{e'}(x')-f^{c'd'a}\,U^{d'b}(z)\,U^{c'e'}(x')\,\rho^{e'}(x')\Big]"""+nl+r"""&\qquad\times\Big[f^{cbd}\,U^{ac}(w)\,\rho^{d}(x)-f^{acd}\,U^{cb}(z)\,U^{de}(x)\,\rho^{e}(x)\Big]"""+nl+r"""&-\frac{1}{(2\pi)^{3}}\,\frac{g^{4}}{8\pi^{5}}\,\frac{1}{k^{+}}\int_{w,x',y',z,x}"""+PH_WY+r"""\,\frac{\Psi(s)}{(x-w)^{2}}"""+nl+r"""&\qquad\times\Big[(x-w)\cdot\mathcal{K}(x'-y')\,\mathcal{K}(z-w)\cdot\mathcal{K}(y'-z)+\frac{1}{2}\,\mathcal{K}(x-z)\cdot\mathcal{K}(y'-z)\,(x-w)\cdot\mathcal{K}(x'-y')\Big]"""+nl+r"""&\qquad\times\Big[f^{c'd'a}\,U^{c'e'}(y')\,U^{d'b}(z)\,\rho^{e'}(x')-f^{c'd'a}\,U^{d'b}(z)\,U^{c'e'}(x')\,\rho^{e'}(x')\Big]"""+nl+r"""&\qquad\times\Big[f^{cbd}\,U^{ac}(w)\,\rho^{d}(x)-f^{acd}\,U^{cb}(z)\,U^{de}(x)\,\rho^{e}(x)\Big]"""+nl+r"""&-\frac{1}{(2\pi)^{3}}\,\frac{g^{4}}{4\pi^{4}}\,\frac{1}{k^{+\,2}}\int_{0}^{\infty}\frac{dp^{+}}{2\pi}\int_{w,x',y',z,x}"""+PH_WY+r"""\,e^{-i\frac{p^{+}}{k^{+}}s}\,\frac{p^{+}}{p^{+}\,(x-z)^{2}+k^{+}\,(x-w)^{2}}"""+nl+r"""&\qquad\times\bigg\{\Big[(x-z)\cdot\mathcal{K}(x'-y')\,\mathcal{K}(z-w)\cdot\mathcal{K}(y'-z)-\frac{1}{2}(x-z)\cdot\mathcal{K}(x'-y')\,\mathcal{K}(x-w)\cdot\mathcal{K}(y'-z)\Big]"""+nl+r"""&\qquad\times\Big[f^{c'd'a}\,U^{c'e'}(y')\,U^{d'b}(z)\,\rho^{e'}(x')-f^{c'd'a}\,U^{d'b}(z)\,U^{c'e'}(x')\,\rho^{e'}(x')\Big]"""+nl+r"""&\qquad-\Big[(x-z)\cdot\mathcal{K}(x'-z)\,\mathcal{K}(z-w)\cdot\mathcal{K}(y'-z)-\frac{1}{2}(x-z)\cdot\mathcal{K}(x'-z)\,\mathcal{K}(x-w)\cdot\mathcal{K}(y'-z)\Big]"""+nl+r"""&\qquad\times\Big[f^{c'd'a}\,U^{c'e'}(z)\,U^{d'b}(z)\,\rho^{e'}(x')-f^{c'd'a}\,U^{d'b}(z)\,U^{c'e'}(x')\,\rho^{e'}(x')\Big]\bigg\}"""+nl+r"""&\qquad\times\Big[f^{cbd}\,U^{ac}(w)\,\rho^{d}(x)-f^{acd}\,U^{cb}(z)\,U^{de}(x)\,\rho^{e}(x)\Big]"""+nl)
blk = L(325, 325)
i0 = blk.index(r"& +\frac{1}{(2\pi)^{3}}\,\frac{g^{4}}{8\pi^{5}}")
blk, n = phase(blk[i0:], 1); NPH[0] += n
blk = blk.replace(r"\right]\times \left[f^{cbd}", r"\right]\notag\\ &\qquad\times\left[f^{cbd}")
A(blk+r"""\,,\label{eq:B-III3two}
\end{align}
$"""+S_G3+r"""$; in the third term the $y'$ integral is done before the $p^{+}$ integral, and the two brackets stay together. The last term is what is left of the term added back for the subtraction in the third term after its $\frac{1}{(x-z)^{2}}\log\frac{\vee}{k^{+}}$ part has been assigned to the BFKL list ($\int_{0}^{\infty}dp^{+}\frac{(x-w)^{2}}{(x-z)^{2}(p^{+}+k^{+})[p^{+}(x-z)^{2}+k^{+}(x-w)^{2}]}=\frac{1}{k^{+}}\frac{(x-w)^{2}}{(x-z)^{2}[(x-z)^{2}-(x-w)^{2}]}\log\frac{(x-z)^{2}}{(x-w)^{2}}$, the form of the author's last term).
\subsubsection*{Rows IV--VI} \textsc{NA}
""")
open('partB.tex', 'w').write(''.join(o))
print('phases replaced in author lines:', NPH[0])
