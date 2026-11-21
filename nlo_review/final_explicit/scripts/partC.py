# Part C: the cancellations, in colours.
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from texlib import *
o = []; A = o.append
nl = r"\notag\\"+"\n"
def col(c, s): return r"{\color{%s}%s}" % (c, s)
A(r"""\section{Cancellations}\label{sec:C}
\noindent\textbf{Colour legend.} Each cancelling set has its own colour; all members of a set are printed in the same colour.
\begin{center}
\renewcommand{\arraystretch}{1.25}\small
\begin{tabular}{|l|p{8.3cm}|p{5.2cm}|l|}\hline
colour & members & type & section\\\hline
"""+col('cI4', r'\rule{2em}{1ex} green')+r""" & Group I Row IV $\leftrightarrow$ Group III Row I & \textbf{exact} (operator identity) & \ref{sec:C1}\\
"""+col('cUV', r'\rule{2em}{1ex} teal')+r""" & Group I Row II $\frac{11}{3}$ term $\leftrightarrow$ Group I Row V $N_{2b\,\rm UV}$ ($-\frac{11}{6}$) & \textbf{exact} (pole, $2\gamma_{E}$, log) & \ref{sec:C2}\\
"""+col('pA', r'\rule{2em}{1ex} orange')+r""" & pair A: Group I Row III (last term of $N_{3b}$) $\leftrightarrow$ Group III Row V ($\big[\frac{k^{+}}{p^{+}(p^{+}+k^{+})}\big]_{+}$ term) & \textbf{tail only} & \ref{sec:C3}\\
"""+col('pB', r'\rule{2em}{1ex} magenta')+r""" & pair B: Group I Row V first part ($\big[\frac{k^{+}}{p^{+}}\big]_{+}$ terms) $\leftrightarrow$ Group III Row III first part ($\big[\frac{k^{+}}{p^{+}(p^{+}+k^{+})}\big]_{+}$ term) & \textbf{tail only} & \ref{sec:C3}\\
"""+col('pC', r'\rule{2em}{1ex} violet')+r""" & pair C: Group I Row V second part ($\big[\frac{k^{+}}{p^{+}}\big]_{+}$, $\big[\frac{k^{+}}{k^{+}-p^{+}}\big]_{+}$ soft terms) $\leftrightarrow$ Group III Row III second part ($\big[\frac{k^{+}}{p^{+}(p^{+}+k^{+})}\big]_{+}$ term) & \textbf{tail only} & \ref{sec:C3}\\
"""+col('pD', r'\rule{2em}{1ex} brown')+r""" & pair D: Group I Row III (third term of $N_{3b}$, part $\int_{0}^{k^{+}}dp^{+}\big[\frac{1}{p^{+}}\big]_{+}e^{-ik\cdot(z-w)\frac{p^{+}}{k^{+}}}$) $\leftrightarrow$ Group III Row VI ($\big[\frac{k^{+}}{p^{+}(p^{+}+k^{+})}\big]_{+}$ term) & \textbf{tail only} & \ref{sec:C3}\\
"""+col('cPr', r'\rule{2em}{1ex} blue')+r""" & Group I Row III: third term of $N_{3b}$, part $\ell-\int_{0}^{k^{+}}dp^{+}\big[\frac{1}{p^{+}}\big]_{+}e^{-ik\cdot(z-w)\frac{p^{+}}{k^{+}}}$ $\leftrightarrow$ the 24-term kernel of the same row & \textbf{tail only, requires the $w$ integration first} & \ref{sec:C4}\\
"""+col('cJt', r'\rule{2em}{1ex} olive')+r""" & Group II Row IV $\leftrightarrow$ Group II Row I (second part) $\leftrightarrow$ Group III Row III (second part), $1:1:-2$ & \textbf{tail only, requires joint integration} & \ref{sec:C6}\\
"""+col('gray', r'\rule{2em}{1ex} gray')+r""" & $i\epsilon$ prescription only: the $\log\frac{k^{+}}{\Lambda}$ (Coulomb) term of $\Delta N$ (Group I Row III) $\leftrightarrow$ its complex conjugate & \textbf{exact} (operator identity) & \ref{sec:C7}\\\hline
\end{tabular}
\end{center}
``Exact'' means that the members are equal and opposite as integrands (after relabelling of integration variables), so that they can be removed. ``Tail only'' means that the members are different functions whose $1/z^{2}$ tails at large distance (each of which alone makes the integral divergent) are equal and opposite; the members are \emph{not} removed: what is left is an integrable remainder, which is part of the finite result and is written out below. Units: $\frac{1}{(2\pi)^{3}}\frac{g^{4}}{16\pi^{5}}\frac{1}{k^{+}}$ appears in all prefactors (as $\frac{1}{(2\pi)^{3}}\frac{g^{4}}{8\pi^{4}}\frac{1}{k^{+}}\int\frac{dp^{+}}{2\pi}$ in the terms with a $p^{+}$ integral).

""")
# ---------------- C1
A(r"""\subsection{Group I Row IV and Group III Row I: exact}\label{sec:C1}
\begingroup\color{cI4}
\begin{align}
""")
b1, _ = phase(L(182, 182), 1)
b1 = b1.replace(r"\frac{d^{3}N_{\rm NLO}}{d^{3}k} & =", r"\frac{d^{3}N^{\rm I.IV}}{d^{3}k} & =")
A(b1+r"\label{eq:C-I4}\\"+"\n")
b2, _ = phase(L(212, 214), 1)
b2 = b2.replace(r"\frac{d^{3}N_{\rm NLO}}{d^{3}k} & =", r"\frac{d^{3}N^{\rm III.I}}{d^{3}k} & =")
A(b2+r"\label{eq:C-III1}"+"\n"+r"""\end{align}
\endgroup
\noindent\textbf{Cancels exactly} (both with c.c., same phase, same integration variables). The kernels are identical (the four terms in the brackets are the same, in a different order), the prefactors are $\pm\frac{ig^{4}}{8\pi^{5}}$, and the colour operators are equal: in Group III Row I, $[U^{ab'}(x')U^{c'd'}(y)U^{c'b}(z)-U^{d'b}(z)U^{d'e'}(y)U^{ac'}(w')\dots]=\big[U^{ab'}(x')-U^{ab'}(w')\big]U^{gb}(z)U^{gd'}(y)\rho^{b'}(x')\rho^{d'}(y)$, so that
\begin{align}
&{\color{cI4}\big[U^{ab'}(x')-U^{ab'}(w')\big]U^{gb}(z)U^{gd'}(y)\,\rho^{b'}(x')\rho^{d'}(y)\Big[f^{cbd}\,U^{ac}(w)\,\rho^{d}(x)-f^{acd}\,U^{cb}(z)\,U^{de}(x)\,\rho^{e}(x)\Big]}\notag\\
&\quad{\color{cI4}=\big[U^{ab'}(x')-U^{ab'}(w')\big]\rho^{b'}(x')\Big[f^{cdc'}\,U^{ac}(w)\,U^{bd}(z)\,U^{be}(y)\,\rho^{e}(y)\rho^{c'}(x)-f^{abc}\,U^{be}(y)\,U^{cd}(x)\,\rho^{e}(y)\rho^{d}(x)\Big]}\,,\notag
\end{align}
the colour operator of Group I Row IV (first term: relabel $b\to d$, $g\to b$, $d'\to e$, $d\to c'$; second term: $U^{gb}(z)U^{cb}(z)=\delta^{gc}$). This uses only relabelling of colour indices and $U^{T}U=1$; the order of the charges ($x'$, $y$, $x$) is the same, so it is an operator identity, and it holds also in colour-singlet form and for the symmetrised parts. (Numerically: colour tensors equal to $1.6\times10^{-15}$, $SU(3)$ and $SU(4)$, Sec.~\ref{sec:F}.) Neither row appears in Secs.~\ref{sec:D}, \ref{sec:E}.

""")
# ---------------- C2
A(r"""\subsection{The $\frac{11}{3}$ ultraviolet terms of Group I Rows II and V: exact}\label{sec:C2}
These terms are not in the author's finite terms (they have been cancelled there); they are shown here explicitly (Divergent\_terms\_fixed, Eqs.~(dv-II), (dv-V); fix\_poles\_dglap Steps 1, 2, 7). Row II (one-$\rho$ part of $\mathbb{A}^{(3)}$, its $\frac{11}{3}$ part) and Row V ($N_{2b\,\rm UV}$, its $-\frac{11}{6}$ part), both with c.c.:
\begingroup\color{cUV}
\begin{align}
\frac{d^{3}N^{\rm I.II}_{11/3}}{d^{3}k}&=\frac{1}{(2\pi)^{3}}\,\frac{g^{4}N_{c}}{32\pi^{4}}\,\frac{1}{k^{+}}\,\frac{11}{3}\int_{w,w',x,x'}"""+PH_WW+r"""\,\mathcal{K}(x'-w')\cdot\mathcal{K}(x-w)\Big[-\frac{2}{\epsilon}-2\gamma_{E}-\log\frac{(x-w)^{2}\mu^{2}_{\overline{MS}}}{4}\Big]"""+nl+r"""&\qquad\times\Big[U^{ab'}(x')-U^{ab'}(w')\Big]\rho^{b'}(x')\Big[U^{ab}(w)-U^{ab}(x)\Big]\rho^{b}(x)\,,\label{eq:C-UV2}\\
\frac{d^{3}N^{\rm I.V}_{2b\,\rm UV}}{d^{3}k}\bigg|_{-\frac{11}{6}}&=\frac{1}{(2\pi)^{3}}\,\frac{g^{4}N_{c}}{16\pi^{4}}\,\frac{1}{k^{+}}\Big(-\frac{11}{6}\Big)\int_{x,y,x',w'}e^{-ik\cdot\ph{(x-w')}}\Big[\frac{2}{\epsilon}+2\gamma_{E}+\log\frac{(y-x)^{2}\mu^{2}_{\overline{MS}}}{4}\Big]\mathcal{K}(x'-w')\cdot\mathcal{K}(y-x)"""+nl+r"""&\qquad\times\Big[U^{ab'}(x')-U^{ab'}(w')\Big]\Big[U^{af}(y)-U^{af}(x)\Big]\rho^{b'}(x')\,\rho^{f}(y)\label{eq:C-UV5}\\
&=\frac{1}{(2\pi)^{3}}\,\frac{g^{4}N_{c}}{32\pi^{4}}\,\frac{1}{k^{+}}\,\frac{11}{3}\int_{w,w',x,x'}"""+PH_WW+r"""\,\mathcal{K}(x'-w')\cdot\mathcal{K}(x-w)\Big[+\frac{2}{\epsilon}+2\gamma_{E}+\log\frac{(x-w)^{2}\mu^{2}_{\overline{MS}}}{4}\Big]"""+nl+r"""&\qquad\times\Big[U^{ab'}(x')-U^{ab'}(w')\Big]\rho^{b'}(x')\Big[U^{ab}(w)-U^{ab}(x)\Big]\rho^{b}(x)=-\frac{d^{3}N^{\rm I.II}_{11/3}}{d^{3}k}\,.\notag
\end{align}
\endgroup
\noindent The second form follows by renaming $x\to w$, $y\to x$: $\frac{g^{4}N_{c}}{16\pi^{4}}(-\frac{11}{6})=-\frac{g^{4}N_{c}}{32\pi^{4}}\frac{11}{3}$ and $[U(y)-U(x)]\to[U(x)-U(w)]=-[U(w)-U(x)]$, with the order $\rho(x')\rho(x)$ unchanged. \textbf{Cancels exactly}: pole, $2\gamma_{E}$ and logarithm, provided the same $\mu_{\overline{MS}}$ is used in both rows. Both poles are ultraviolet (UV against UV). Consequently there is no finite term $\frac{\alpha_{s}\beta_{0}}{4\pi}[\log(x-w)^{2}+\log(x'-w')^{2}-\log(w-w')^{2}]\times$LO. The remaining parts of the two rows are Eq.~\eqref{eq:B-I2a} ($-\frac{67}{9}+\frac{\pi^{2}}{3}$) and the $\frac{11}{6}$ collinear subtraction in Eq.~\eqref{eq:B-I5two} (last line), which are finite.

""")
# ---------------- C3 pairs
PREF_I3 = r"\frac{1}{(2\pi)^{3}}\,\frac{ig^{4}f^{abd}}{8\pi^{5}}\,\frac{1}{k^{+}}"
WW = r"\int_{w,w'}"+PH_WW+r"\int_{x,y,z,x'}"
CA = r"\Big[U^{ab'}(x')-U^{ab'}(w')\Big]\rho^{b'}(x')\,U^{de}(x)\Big[U^{bc}(w)-U^{bc}(y)\Big]\rho^{e}(x)\rho^{c}(y)"
CD = r"\Big[U^{ab'}(x')-U^{ab'}(w')\Big]\rho^{b'}(x')\,U^{de}(x)\Big[U^{bc}(z)-U^{bc}(y)\Big]\rho^{e}(x)\rho^{c}(y)"
CB = r"\Big[U^{ab'}(x')-U^{ab'}(w')\Big]\rho^{b'}(x')\Big[f^{adc}\,U^{db}(z)\,U^{ce}(w)\,\{\rho^{e}(y),\rho^{b}(x)\}-f^{abc}\,U^{cd}(y)\,U^{be}(x)\,\{\rho^{d}(y),\rho^{e}(x)\}\Big]"
CC = r"\Big[U^{ab'}(x')-U^{ab'}(w')\Big]\Big[N_{c}\,U^{ah}(x)-f^{adc}f^{hbe}\,U^{db}(z)\,U^{ce}(w)\Big]\rho^{b'}(x')\,\rho^{h}(x)"
EW = SE('z-w', True)    # Group I soft integral (dp+/(2pi); prefactor x 2pi)
PW = SP("z-w'", True)   # Group III soft integral
PREF2_I3 = r"\frac{1}{(2\pi)^{3}}\,\frac{ig^{4}f^{abd}}{4\pi^{4}}\,\frac{1}{k^{+}}"
PB2 = r"\frac{1}{(2\pi)^{3}}\,\frac{ig^{4}}{8\pi^{4}}\,\frac{1}{k^{+}}"
PC2 = r"\frac{1}{(2\pi)^{3}}\,\frac{g^{4}}{8\pi^{4}}\,\frac{1}{k^{+}}"
G1S = DPK+PL1   # int_0^k+ dp+/(2pi) [1/p+]_+
G3S = DPI+PL3   # int_0^inf dp+/(2pi) [k+/(p+(p++k+))]_+

def disp(color, rows, label=None):
    body = (r"\notag\\"+"\n").join(rows)
    if label: body += r"\label{%s}" % label
    else: body += r"\notag"
    if color: return r"\begingroup\color{%s}"%color+"\n"+r"\begin{align}"+"\n"+body+"\n"+r"\end{align}"+"\n"+r"\endgroup"+"\n"
    return r"\begin{align}"+"\n"+body+"\n"+r"\end{align}"+"\n"
def pairsum(color, pref, measure, V, kw, kwp, colour, label, sign='-'):
    c = lambda t: r"{\color{%s}%s}" % (color, t)
    return disp(None, [
        r"\frac{d^{3}N^{\rm %s}}{d^{3}k}&=" % label[-1]+c(sign+pref+measure),
        "&\\qquad"+c(r"\times"+V+r"\bigg\{\Big["+kw+"-"+kwp+r"\Big]"),
        "&\\qquad\\qquad"+c(r"\times"+EW),
        "&\\qquad"+c(r"+"+kwp+r"\Big["+EW),
        "&\\qquad\\qquad"+c("-"+PW+r"\Big]\bigg\}"),
        "&\\qquad"+c(r"\times"+colour)], label)
PB_I = r"\frac{1}{(2\pi)^{3}}\,\frac{ig^{4}}{16\pi^{5}}\,\frac{1}{k^{+}}"
PC_I = r"\frac{1}{(2\pi)^{3}}\,\frac{g^{4}}{16\pi^{5}}\,\frac{1}{k^{+}}"
PAIRS = r"\paragraph{Pair A}"+"\n"+disp('pA', [
 r"&\text{Group I Row III, last term of $N_{3b}$, Eq.~\eqref{eq:B-I3}:}",
 r"&-"+PREF2_I3+G1S+WW+EP('z-w'), r"&\qquad\times\mathcal{K}(x-z)\cdot\mathcal{K}(w-z)\,\mathcal{K}(x'-w')\cdot\mathcal{K}(y-w)",
 r"&\qquad\times"+CA,
 r"&\text{Group III Row V, Eq.~\eqref{eq:B-III5}, soft remainder:}",
 r"&-\frac{1}{(2\pi)^{3}}\,\frac{ig^{4}f^{c'ca}}{4\pi^{4}}\,\frac{1}{k^{+}}"+G3S+r"\int_{w,y',z,y,x,x'}"+PH_WY+PH_G3, r"&\qquad\times\mathcal{K}(y'-z)\cdot\mathcal{K}(x-z)\,\mathcal{K}(x'-y')\cdot\mathcal{K}(y-w)",
 r"&\qquad\times\Big[U^{c'e'}(y')-U^{c'e'}(x')\Big]\rho^{e'}(x')\,U^{ce}(x)\Big[U^{ad}(y)-U^{ad}(w)\Big]\rho^{e}(x)\rho^{d}(y)",
 r"&\overset{y'\to w'}{=}+"+PREF2_I3+G3S+WW+EP("z-w'"), r"&\qquad\times\mathcal{K}(x-z)\cdot\mathcal{K}(w'-z)\,\mathcal{K}(x'-w')\cdot\mathcal{K}(y-w)",
 r"&\qquad\times"+CA+"."])+ \
 r"Sum (the integrable remainder, Eq.~\eqref{eq:C-pairform} with $\kappa(u)=\mathcal{K}(x-z)\cdot\mathcal{K}(u-z)$):"+"\n"+ \
 pairsum('pA', PREF2_I3, WW, r"\mathcal{K}(x'-w')\cdot\mathcal{K}(y-w)", r"\mathcal{K}(x-z)\cdot\mathcal{K}(w-z)", r"\mathcal{K}(x-z)\cdot\mathcal{K}(w'-z)", CA, 'eq:C-A')+ \
 r"\paragraph{Pair B}"+"\n"+disp('pB', [
 r"&\text{Group I Row V, first part, soft remainders, Eq.~\eqref{eq:B-I5}:}",
 r"&-"+PF_I5+DPK+r"\int_{x,z,y,y',x',w'}"+PH_ZW+r"\Big[",
 r"&\qquad"+PLK+EP('x-z')+r"\,\mathcal{K}(x'-w')\cdot\mathcal{K}(y-z)\,\mathcal{K}(z-x)\cdot\mathcal{K}(y'-x)",
 r"&\qquad+"+PLK+r"e^{-ik\cdot\ph{(x-z)}}"+EP('z-x')+r"\,\mathcal{K}(x'-w')\cdot\mathcal{K}(y'-x)\,\mathcal{K}(z-x)\cdot\mathcal{K}(y-z)\Big]",
 r"&\qquad\times\Big[U^{ab'}(x')-U^{ab'}(w')\Big]\rho^{b'}(x')\Big[f^{adc}\,U^{db}(x)\,U^{ce}(z)\,\{\rho^{e}(y),\rho^{b}(y')\}-f^{abc}\,U^{cd}(y)\,U^{be}(y')\,\{\rho^{d}(y),\rho^{e}(y')\}\Big]",
 r"&=-"+PB2+G1S+WW+EP('z-w'), r"&\qquad\times\mathcal{K}(x-z)\cdot\mathcal{K}(w-z)\,\mathcal{K}(x'-w')\cdot\mathcal{K}(y-w)",
 r"&\qquad\times"+CB,
 r"&\text{Group III Row III, first part, soft remainder, Eq.~\eqref{eq:B-III3one}:}",
 r"&-\frac{1}{(2\pi)^{3}}\,\frac{ig^{4}f^{c'd'a}}{8\pi^{4}}\,\frac{1}{k^{+}}"+G3S+r"\int_{w,x',y',z,y,x}"+PH_WY+PH_G3, r"&\qquad\times\mathcal{K}(y'-z)\cdot\mathcal{K}(x-z)\,\mathcal{K}(x'-y')\cdot\mathcal{K}(y-w)",
 r"&\qquad\times U^{d'b}(z)\Big[U^{c'e'}(y')-U^{c'e'}(x')\Big]\rho^{e'}(x')",
 r"&\qquad\times\Big[U^{ac}(w)\,\{\rho^{c}(y),\rho^{b}(x)\}-U^{cb}(z)\,U^{ad}(y)\,U^{ce}(x)\,\{\rho^{d}(y),\rho^{e}(x)\}\Big]",
 r"&\overset{y'\to w'}{=}+"+PB2+G3S+WW+EP("z-w'"), r"&\qquad\times\mathcal{K}(x-z)\cdot\mathcal{K}(w'-z)\,\mathcal{K}(x'-w')\cdot\mathcal{K}(y-w)",
 r"&\qquad\times"+CB+"."])+ \
 r"The Group I member is written in common variables by renaming, in its first term, $x\to z$ (soft), $z\to w$ (hard), $y'\to x$, and in its second term $x\to w$ (hard), $y\leftrightarrow y'$, then $y'\to x$ (the recoil phase becomes $e^{-ik\cdot(z-w)p^{+}/k^{+}}$ in both, $e^{-ik\cdot(z-w')}e^{-ik\cdot(x-z)}=e^{-ik\cdot(x-w')}$, and $\frac{1}{k^{+\,2}}\big[\frac{k^{+}}{p^{+}}\big]_{+}=\frac{1}{k^{+}}\big[\frac{1}{p^{+}}\big]_{+}$); the two terms are then equal (in the second both the kernel and the colour factor change sign). Sum:"+"\n"+ \
 pairsum('pB', PB2, WW, r"\mathcal{K}(x'-w')\cdot\mathcal{K}(y-w)", r"\mathcal{K}(x-z)\cdot\mathcal{K}(w-z)", r"\mathcal{K}(x-z)\cdot\mathcal{K}(w'-z)", CB, 'eq:C-B')+ \
 r"\paragraph{Pair C}"+"\n"+disp('pC', [
 r"&\text{Group I Row V, second part, soft remainders, Eq.~\eqref{eq:B-I5two}:}",
 r"&-\frac{1}{(2\pi)^{3}}\,\frac{g^{4}}{8\pi^{4}}\,\frac{1}{k^{+\,2}}"+DPK+r"\int_{x,y,z,x',w'}"+PH_ZW+r"\bigg\{",
 r"&\qquad"+PLK+r"\frac{"+EP('x-z')+r"}{(y-z)^{2}}\Big[-\frac{(y-z)\cdot\mathcal{K}(x'-w')}{(x-z)^{2}}+\frac{1}{2}\,(y-z)\cdot\mathcal{K}(x'-w')\,\mathcal{K}(z-x)\cdot\mathcal{K}(y-x)\Big]",
 r"&\qquad+"+PLKQ+r"\frac{"+EP('x-z')+r"}{(y-x)^{2}}\Big[-\frac{(y-x)\cdot\mathcal{K}(x'-w')}{(x-z)^{2}}-\frac{1}{2}\,(y-x)\cdot\mathcal{K}(x'-w')\,\mathcal{K}(z-x)\cdot\mathcal{K}(y-z)\Big]\bigg\}",
 r"&\qquad\times\Big[U^{ab'}(x')-U^{ab'}(w')\Big]\Big[-f^{adc}\,f^{fbe}\,U^{db}(x)\,U^{ce}(z)+N_{c}\,U^{af}(y)\Big]\rho^{b'}(x')\,\rho^{f}(y)",
 r"&=-"+PC2+G1S+r"\int_{w,w'}"+PH_WW+r"\int_{x,z,x'}",
 r"&\qquad\times\,"+EP('z-w')+r"\,\mathcal{K}(x'-w')\cdot\mathcal{K}(x-w)\Big[\mathcal{K}(x-z)\cdot\mathcal{K}(w-z)-\frac{2}{(z-w)^{2}}\Big]",
 r"&\qquad\times"+CC,
 r"&\text{Group III Row III, second part, soft remainder, Eq.~\eqref{eq:B-III3two}:}",
 r"&-\frac{1}{(2\pi)^{3}}\,\frac{g^{4}}{4\pi^{4}}\,\frac{1}{k^{+}}"+G3S+r"\int_{w,x',y',z,x}"+PH_WY+PH_G3,
 r"&\qquad\times\frac{1}{(x-w)^{2}}\Big[(x-w)\cdot\mathcal{K}(x'-y')\,\mathcal{K}(z-w)\cdot\mathcal{K}(y'-z)+\frac{1}{2}\,\mathcal{K}(x-z)\cdot\mathcal{K}(y'-z)\,(x-w)\cdot\mathcal{K}(x'-y')\Big]",
 r"&\qquad\times\Big[f^{c'd'a}\,U^{c'e'}(y')\,U^{d'b}(z)\,\rho^{e'}(x')-f^{c'd'a}\,U^{d'b}(z)\,U^{c'e'}(x')\,\rho^{e'}(x')\Big]\Big[f^{cbd}\,U^{ac}(w)\,\rho^{d}(x)-f^{acd}\,U^{cb}(z)\,U^{de}(x)\,\rho^{e}(x)\Big]",
 r"&\overset{y'\to w'}{=}+"+PC2+G3S+r"\int_{w,w'}"+PH_WW+r"\int_{x,z,x'}",
 r"&\qquad\times\,"+EP("z-w'")+r"\,\mathcal{K}(x'-w')\cdot\mathcal{K}(x-w)\Big[\mathcal{K}(x-z)\cdot\mathcal{K}(w'-z)-2\,\mathcal{K}(w-z)\cdot\mathcal{K}(w'-z)\Big]",
 r"&\qquad\times"+CC+"."])+ \
 r"The Group I member in common variables: first term $x\to z$, $z\to w$, $y\to x$; second term $p^{+}\to k^{+}-p^{+}$ (which turns $\big[\frac{k^{+}}{k^{+}-p^{+}}\big]_{+}e^{-ik\cdot(x-z)p^{+}/k^{+}}$, subtraction at $p^{+}=k^{+}$, into $\big[\frac{k^{+}}{p^{+}}\big]_{+}e^{-ik\cdot(x-z)}e^{-ik\cdot(z-x)p^{+}/k^{+}}$, subtraction at $p^{+}=0$), then $x\to w$, $y\to x$; both give $\frac{1}{k^{+}}\big[\frac{1}{p^{+}}\big]_{+}e^{-ik\cdot(w-w')}e^{-ik\cdot(z-w)p^{+}/k^{+}}\,\mathcal{K}(x'-w')\cdot\mathcal{K}(x-w)\big[\frac12\mathcal{K}(x-z)\cdot\mathcal{K}(w-z)-\frac{1}{(z-w)^{2}}\big]$ under $\frac{1}{(2\pi)^{3}}\frac{g^{4}}{8\pi^{4}}\int_{0}^{k^{+}}\frac{dp^{+}}{2\pi}$, with the colour factor $[U(x')-U(w')]^{ab'}[N_{c}U^{ah}(x)-f^{adc}f^{hbe}U^{db}(z)U^{ce}(w)]$ (using $f^{adc}f^{hbe}U^{db}(w)U^{ce}(z)=f^{adc}f^{hbe}U^{db}(z)U^{ce}(w)$). Sum:"+"\n"+ \
 pairsum('pC', PC2, r"\int_{w,w'}"+PH_WW+r"\int_{x,z,x'}", r"\mathcal{K}(x'-w')\cdot\mathcal{K}(x-w)", r"\Big(\mathcal{K}(x-z)\cdot\mathcal{K}(w-z)-\frac{2}{(z-w)^{2}}\Big)", r"\Big(\mathcal{K}(x-z)\cdot\mathcal{K}(w'-z)-2\,\mathcal{K}(w-z)\cdot\mathcal{K}(w'-z)\Big)", CC, 'eq:C-C')+ \
 r"\paragraph{Pair D}"+"\n"+disp('pD', [
 r"&\text{Group I Row III, soft part of the third term of $N_{3b}$ (Sec.~\ref{sec:C4}):}",
 r"&+"+PREF2_I3+G1S+WW+EP('z-w'), r"&\qquad\times\mathcal{K}(y-z)\cdot\mathcal{K}(w-z)\,\mathcal{K}(x'-w')\cdot\mathcal{K}(x-w)",
 r"&\qquad\times"+CD,
 r"&\text{Group III Row VI, Eq.~\eqref{eq:B-III6}, soft remainder:}",
 r"&-\frac{1}{(2\pi)^{3}}\,\frac{ig^{4}f^{c'd'a}}{4\pi^{4}}\,\frac{1}{k^{+}}"+G3S+r"\int_{w,y',z,y,x,x'}"+PH_WY+PH_G3, r"&\qquad\times\mathcal{K}(y'-z)\cdot\mathcal{K}(x-z)\,\mathcal{K}(x'-y')\cdot\mathcal{K}(y-w)",
 r"&\qquad\times\Big[U^{c'e'}(y')-U^{c'e'}(x')\Big]\rho^{e'}(x')\,U^{ad}(y)\Big[U^{d'e}(x)-U^{d'e}(z)\Big]\rho^{d}(y)\rho^{e}(x)",
 r"&\overset{y'\to w',\,x\leftrightarrow y}{=}-"+PREF2_I3+G3S+WW+EP("z-w'"), r"&\qquad\times\mathcal{K}(y-z)\cdot\mathcal{K}(w'-z)\,\mathcal{K}(x'-w')\cdot\mathcal{K}(x-w)",
 r"&\qquad\times"+CD+"."])+ \
 "Sum:\n"+pairsum('pD', PREF2_I3, WW, r"\mathcal{K}(x'-w')\cdot\mathcal{K}(x-w)", r"\mathcal{K}(y-z)\cdot\mathcal{K}(w-z)", r"\mathcal{K}(y-z)\cdot\mathcal{K}(w'-z)", CD, 'eq:C-D', sign='+')
A(r"""\subsection{The soft remainders: pairs A--D (tail only)}\label{sec:C3}
Each Group I soft remainder (soft gluon absorbed by the observed gluon in the amplitude: virtual, recoil on $w$, $\int_{0}^{k^{+}}\frac{dp^{+}}{2\pi}\big[\frac{1}{p^{+}}\big]_{+}e^{-ik\cdot(z-w)p^{+}/k^{+}}$, upper limit $k^{+}$) has one Group III partner (soft gluon emitted by the observed gluon in the conjugate amplitude: real, recoil on $w'$, $\int_{0}^{\infty}\frac{dp^{+}}{2\pi}\big[\frac{k^{+}}{p^{+}(p^{+}+k^{+})}\big]_{+}e^{-ik\cdot(z-w')p^{+}/k^{+}}$, upper limit $\infty$) with the same coefficient, the same colour structure and the same kernel, except that $w$ is replaced by $w'$ in the kernel (fix\_soft Steps 2--4; paper Table~1). To see this, the Group III members are written with the integration variable $y'$ renamed $w'$, and the Group I Row V members with the row variables renamed (soft gluon at $z$, hard gluon at $w$). The colour identities used (operator identities: relabelling of colour indices, $U^{d'b}(z)U^{cb}(z)=\delta^{d'c}$, $f^{c'd'a}f^{dd'a}=N_{c}\delta^{c'd}$, no reordering of charges; checked numerically, Sec.~\ref{sec:F}) are
\begin{align}
&f^{c'ca}\Big[U^{c'e'}(w')-U^{c'e'}(x')\Big]\rho^{e'}(x')\,U^{ce}(x)\Big[U^{ad}(y)-U^{ad}(w)\Big]\rho^{e}(x)\rho^{d}(y)\notag\\
&\qquad=-f^{abd}"""+CA+r""",\notag\\
&f^{c'd'a}U^{d'b}(z)\Big[U^{c'e'}(w')-U^{c'e'}(x')\Big]\rho^{e'}(x')\Big[U^{ac}(w)\{\rho^{c}(y),\rho^{b}(x)\}-U^{cb}(z)U^{ad}(y)U^{ce}(x)\{\rho^{d}(y),\rho^{e}(x)\}\Big]\notag\\
&\qquad=-"""+CB+r""",\notag\\
&f^{c'd'a}U^{d'b}(z)\Big[U^{c'e'}(w')-U^{c'e'}(x')\Big]\rho^{e'}(x')\Big[f^{cbd}U^{ac}(w)\rho^{d}(x)-f^{acd}U^{cb}(z)U^{de}(x)\rho^{e}(x)\Big]\notag\\
&\qquad=-"""+CC+r""",\notag\\
&f^{c'd'a}\Big[U^{c'e'}(w')-U^{c'e'}(x')\Big]\rho^{e'}(x')\,U^{ad}(x)\Big[U^{d'e}(y)-U^{d'e}(z)\Big]\rho^{d}(x)\rho^{e}(y)\notag\\
&\qquad=+f^{abd}"""+CD+r""".\label{eq:C-colid}
\end{align}
(the last one is the colour structure of Group III Row VI with $y'\to w'$ and $x\leftrightarrow y$). Each pair has the form (the Group I soft integral, the Group III soft integral; $dp^{+}=2\pi\frac{dp^{+}}{2\pi}$)
\begin{align}
\kappa(w)\int_{0}^{k^{+}}dp^{+}\Big[\frac{1}{p^{+}}\Big]_{+}e^{-ik\cdot\ph{(z-w)}\frac{p^{+}}{k^{+}}}-\kappa(w')\int_{0}^{\infty}dp^{+}\Big[\frac{k^{+}}{p^{+}\,(p^{+}+k^{+})}\Big]_{+}e^{-ik\cdot\ph{(z-w')}\frac{p^{+}}{k^{+}}}\notag
\end{align}
times a common factor, which is rewritten exactly as
\begin{align}
&\kappa(w)\int_{0}^{k^{+}}dp^{+}\Big[\frac{1}{p^{+}}\Big]_{+}e^{-ik\cdot\ph{(z-w)}\frac{p^{+}}{k^{+}}}\notag\\
&\qquad-\kappa(w')\int_{0}^{\infty}dp^{+}\Big[\frac{k^{+}}{p^{+}\,(p^{+}+k^{+})}\Big]_{+}e^{-ik\cdot\ph{(z-w')}\frac{p^{+}}{k^{+}}}\notag\\
&\quad=\big[\kappa(w)-\kappa(w')\big]\int_{0}^{k^{+}}dp^{+}\Big[\frac{1}{p^{+}}\Big]_{+}e^{-ik\cdot\ph{(z-w)}\frac{p^{+}}{k^{+}}}\notag\\
&\qquad+\kappa(w')\bigg[\int_{0}^{k^{+}}dp^{+}\Big[\frac{1}{p^{+}}\Big]_{+}e^{-ik\cdot\ph{(z-w)}\frac{p^{+}}{k^{+}}}-\int_{0}^{\infty}dp^{+}\Big[\frac{k^{+}}{p^{+}\,(p^{+}+k^{+})}\Big]_{+}e^{-ik\cdot\ph{(z-w')}\frac{p^{+}}{k^{+}}}\bigg]\,.\label{eq:C-pairform}
\end{align}
At large $z$ each member behaves like $\pm\frac{T}{z^{2}}\big[-\log|k\cdot z|-\gamma_{E}-\frac{i\pi}{2}{\rm sgn}(k\cdot z)\big]$ (both soft integrals tend to this, Eq.~\eqref{eq:asym}; not integrable); in the sum, $\kappa(w)-\kappa(w')=O(|w-w'|/|z|^{3})$ while the Group I soft integral is bounded by $\log|z|+{\rm const}$, and the difference of the two soft integrals (last bracket) is $\log\frac{|k\cdot(z-w')|}{|k\cdot(z-w)|}+O(1/|k\cdot z|)$ (the odd parts $-\frac{i\pi}{2}{\rm sgn}$ cancel), $O(1/|z|)$ outside the strip $|k\cdot z|\lesssim|k||w-w'|$ and bounded by logarithms inside it: \emph{the tails cancel, the remainder~\eqref{eq:C-pairform} is integrable}, for any $U(z)$ (the colour factor, with all its dependence on $U(z)$, is common to both terms). Numerically (fix\_soft Step 3), $z^{2}\times$(sum of the two members) falls like $1/|z|$ in every direction, for all four pairs, $SU(2)$ and $SU(3)$, two- and three-$\rho$ sectors. The complex-conjugate halves pair in the same way.

"""+PAIRS+r"""The kernels of the members in common variables were checked pointwise against the row forms (Sec.~\ref{sec:F}).

""")
# ---------------- C4 problematic term
A(r"""\subsection{The problematic third term of $N_{3b}$ (Group I Row III): tail only, after the $w$ integration}\label{sec:C4}
The third term of $N_{3b}$, contained in the last line of Eq.~\eqref{eq:B-I3}, is
\begin{align}
&"""+PREF_I3+WW+r"""\notag\\
&\quad\times\Big(-\frac12\log\frac{(z-w)^{2}}{(x-w)^{2}}\Big)\,\mathcal{K}(y-z)\cdot\mathcal{K}(w-z)\,\mathcal{K}(x'-w')\cdot\mathcal{K}(x-w)\notag\\
&\quad\times"""+CD+r""".\notag
\end{align}
With $a_{t}=(x-w)-t(x-z)$, $-\frac12\log\frac{(z-w)^{2}}{(x-w)^{2}}=-\int_{0}^{1}\frac{dt}{t}\big[1-\frac{(x-w)\cdot a_{t}}{a_{t}^{2}}\big]$: the LO kernel minus the LO kernel with the source shifted by the recoil; the soft integral is cut off at $t\simeq|x-w|/|z-w|$ without a phase, and nothing of this term is in the JIMWLK list. It is integrable at $z\to w$ ($\sim\log|z-w|/|z-w|$), but at large $z$ it behaves like $-\frac{1}{z^{2}}\log\frac{|z|}{|x-w|}$ times $\mathcal{K}(x'-w')\cdot\mathcal{K}(x-w)\times$colour: isotropic, without the odd part $-\frac{i\pi}{2}{\rm sgn}$ of the soft integral. It is split exactly as
\begin{align}
-\frac12\log\frac{(z-w)^{2}}{(x-w)^{2}}={\color{pD}\int_{0}^{k^{+}}dp^{+}\Big[\frac{1}{p^{+}}\Big]_{+}e^{-ik\cdot\ph{(z-w)}\frac{p^{+}}{k^{+}}}}+{\color{cPr}\Big[-\frac12\log\frac{(z-w)^{2}}{(x-w)^{2}}-\int_{0}^{k^{+}}dp^{+}\Big[\frac{1}{p^{+}}\Big]_{+}e^{-ik\cdot\ph{(z-w)}\frac{p^{+}}{k^{+}}}\Big]}\,.\notag
\end{align}
{\color{pD}The first part is the Group I member of pair D} (its $\log|z|$ tail, for the $U(z)$-independent and the $U(z)$-linear colour parts alike, is cancelled by Group III Row VI, Eq.~\eqref{eq:C-D}). {\color{cPr}The second part is combined with the 24-term kernel of the same row (lines 1--4 of Eq.~\eqref{eq:B-I3} without the third term of $N_{3b}$):}
\begingroup\color{cPr}
\begin{align}
\widehat{\mathcal{M}}&\equiv"""+KZ_24[0]+nl+r"&\quad"+KZ_24[1]+nl+r"&\quad"+KZ_24[2]+nl+r"&\quad"+KZ_24[3]+nl+r"&\quad"+KZ_24[4]+nl+r"&\quad"+KZ_24[5]+nl+r"""&\quad+"""+KZ_LE+r"""\,.\label{eq:C-Mhat}
\end{align}
At fixed $x,y,x',w'$ and $|z|\to\infty$, with $\varphi\in(-\pi,\pi)$ the angle from $x-w$ to $z$, $e_{\perp}$ the unit vector $(x-w)/|x-w|$ rotated by $+\frac{\pi}{2}$, $\theta$ the angle between $k$ and $z$:
\begin{align}
\text{(24 terms)}&\;\longrightarrow\;\frac{\varphi}{z^{2}}\,\frac{e_{\perp}\cdot\mathcal{K}(x'-w')}{|x-w|}+O\Big(\frac{\log|z|}{|z|^{3}}\Big)\quad\text{(from the arccos terms; zero angular mean, not zero pointwise)},\notag\\
\text{(last line)}&\;\longrightarrow\;\frac{\mathcal{K}(x'-w')\cdot\mathcal{K}(x-w)}{z^{2}}\Big[\gamma_{E}+\log\big(|k|\,|x-w|\big)+\log|\cos\theta|+\frac{i\pi}{2}{\rm sgn}(\cos\theta)\Big]+\big(\text{terms}\propto e^{-ik\cdot(z-w)}\big)\,.\notag
\end{align}
Pointwise in $w$ the sum still has a $1/z^{2}$ tail. After the Fourier integration over $w$ ($r=x-w$, $e^{-ik\cdot w}=e^{-ik\cdot x}e^{ik\cdot r}$, convergence factor $e^{-\epsilon r}$), with $\chi(\theta)\equiv\log|2\cos\theta|+\frac{i\pi}{2}{\rm sgn}(\cos\theta)$,
\begin{align}
\int d^{2}r\,e^{ik\cdot r}\,\frac{r^{i}}{r^{2}}\Big[\gamma_{E}+\log(|k||r|)+\log|\cos\theta|+\frac{i\pi}{2}{\rm sgn}(\cos\theta)\Big]&=+\frac{2\pi i\,k^{i}}{|k|^{2}}\,\chi(\theta)\,,\qquad
\int d^{2}r\,e^{ik\cdot r}\,\varphi\,\frac{e_{\perp}^{i}}{|r|}=-\frac{2\pi i\,k^{i}}{|k|^{2}}\,\chi(\theta)\,:\notag
\end{align}
the two cancel, and
\begin{align}
\int d^{2}w\,e^{-ik\cdot w}\,\widehat{\mathcal{M}}=O\Big(\frac{\log|z|}{|z|^{3}}\Big)+\big(\text{terms}\propto e^{-ik\cdot z}\big)\qquad(|z|\to\infty)\,,\notag
\end{align}
in every direction and for any $U(z)$ (the colour factor does not depend on $w$).
\endgroup
\noindent\textbf{Type: tail only, requires the $w$ integration before the $z$ integral} (paper Eqs.~(saw), (sawid), (Mhatclaim); fix\_soft Step 5: at $|z|=2560$ the two parts of $z^{2}\int_{w}e^{-ik\cdot w}\widehat{\mathcal{M}}$ are $6.933-2.690i$ and $-6.906+2.690i$, prediction $\pm(6.909-2.687i)$; $z^{2}\times$sum falls like $1/|z|$). The large-$z$ form of the 24-term kernel is established numerically; the Fourier identity is analytic. \emph{Answer to ``what cancels the problematic term'':} its $\log|z|$ part is cancelled by the soft remainder of Group III Row VI (pair D), its remaining $1/z^{2}$ part by the arccos (``odd-angle'') tail of the other 24 terms of its own row, after the $w$ integration.

""")
# ---------------- C5 within-row
A(r"""\subsection{Cancellations inside single rows (no colour)}\label{sec:C5}
(i) Group I Row III: the poles on $(x-z)^{2}=(x-w)^{2}$ of the separate terms $N_{1},N_{2},N_{3a}$ (principal values, $\log|\cdot|$) cancel in the combined kernel (author's form; kernel value $0.19948358$ and $0.19948353$ at $\pm10^{-7}$ from the line, paper Sec.~5.3). (ii) Group II Row I, second part: poles on $(x'-w')^{2}(x-z)^{2}=(x-w)^{2}(x'-z)^{2}$ cancel in $N_{2a}+N_{2b}$ and $N_{2c}+N_{2d}$ (Sec.~\ref{sec:A}, A6). (iii) Group III, instantaneous terms: the terms with $\mathcal{K}(x'-y')$ and $\mathcal{K}(x'-z)$ each grow like $\log\vee$ (BFKL term) and only their difference is finite. (iv) Group I Row II, $N_{2b}$: the non-integrable points $\mathbf{p}\to0$ and $\mathbf{q}\to\mathbf{p}$ of the author's kernel cancel in the antisymmetrised kernel. (v) Group I Row V, second part: the collinear point $z\to x$ of the regular term is integrable (angular integral first) by virtue of the $\frac{11}{6}$ term.

""")
# ---------------- C6 joint tail
A(r"""\subsection{The joint large-distance tail of Group II Row IV: tail only, requires joint integration}\label{sec:C6}
In the region in which the observed gluon, the unobserved gluon and their parent are all far from the projectile together ($w,w',z,y'\to\infty$ at fixed $w-z$, $w'-z$, $y'-z$; $z=Z\,n$, $n=(\cos\phi,\sin\phi)$, $\phi$ the angle between $k$ and $n$), three rows have $1/Z^{2}$ tails. With the far Wilson lines equal, the colour factors of the three rows are $+N_{c}\Gamma$ (Group II Row IV), $+N_{c}\Gamma$ (Group II Row I, all of the second part), $-N_{c}\Gamma$ (Group III Row III, second part), $\Gamma=[U(Y)-U(x')]^{ce'}[U(Y)-U(x)]^{ce}\rho^{e'}(x')\rho^{e}(x)$, and the three rows have the same tail density per $dt$ ($t=p^{+}/k^{+}$), $\varrho(t)=\frac{(t^{2}+t+1)^{2}+t^{2}\cos2\phi}{t\,(1+t)^{4}}$, with weights $+1,+1,-2$ (Group III with c.c.). In units
\begin{align}
\mathsf{u}\equiv2N_{c}\,\frac{1}{(2\pi)^{3}}\,\frac{g^{4}}{16\pi^{4}}\,\frac{1}{k^{+}}\,\frac{4\pi}{k^{2}}\int\frac{d^{2}Y}{Y^{2}}\int_{x,x'}\Big[U^{ce'}(Y)-U^{ce'}(x')\Big]\Big[U^{ce}(Y)-U^{ce}(x)\Big]\rho^{e'}(x')\rho^{e}(x)\notag
\end{align}
the tails of the finite terms as written are (fix\_tail, Steps 2--6):
\begin{center}
\renewcommand{\arraystretch}{1.3}\small
\begin{tabular}{|l|c|}\hline
{\color{cJt}Group II Row IV: DGLAP $+\frac{1}{12}$, RC $-\frac{11}{12}$, remnant Eq.~\eqref{eq:B-II4rem} $-1$, tensor term Eq.~\eqref{eq:B-II4} $+\frac16\cos2\phi$} & {\color{cJt}$-\frac{11}{6}+\frac16\cos2\phi$}\\
{\color{cJt}Group II Row I, second part, Eqs.~\eqref{eq:B-II2a}--\eqref{eq:B-II2e}} & {\color{cJt}$-\frac{11}{6}+\frac16\cos2\phi$}\\
{\color{cJt}Group III Row III, second part, Eq.~\eqref{eq:B-III3two}, with c.c.} & {\color{cJt}$+\frac{11}{3}-\frac13\cos2\phi$}\\
Group II Row I, first part; Group III Row III, first part & $0$\\\hline
sum & $0$\\\hline
\end{tabular}
\end{center}
and the BFKL and JIMWLK coefficients of the three rows carry $1+1-2=0$ times $\mathsf{u}$ each. \textbf{Type: tail only, requires joint integration}: each row alone is logarithmically divergent at large pair position $Y$; the three rows have to be integrated together under one integral over the position of the pair (with $Y=\frac{y+y'}{2}$ of Row IV identified with $z$ of the other rows). This holds for any target uniform at large distance, and for a colour-neutral projectile for any target; for a coloured projectile and a target not uniform at large distance a remainder proportional to the total charge is left (the non-singlet term of the JIMWLK formula). Note that the remnant of Group II Row IV (missing in the author's file, \NOTE{} A5) is needed for this cancellation.
""")
A(r"""\subsection{The Coulomb logarithm of $\Delta N$ ($i\epsilon$ prescription): exact, against the complex conjugate}\label{sec:C7}
In the $i\epsilon$ prescription the term $\Delta N$, Eq.~\eqref{eq:B-DN}, contains
\begingroup\color{gray}
\begin{align}
T_{\Lambda}&="""+PF_DN+r"""\,\pi\log\frac{k^{+}}{\Lambda}\int_{w,w'}"""+PH_WW+r"""\int_{x,y,x'}\mathcal{K}(x'-w')\cdot\mathcal{K}(y-w)\,f^{abd}\Big[U^{ab'}(x')-U^{ab'}(w')\Big]\rho^{b'}(x')\,U^{de}(x)\rho^{e}(x)\Big[U^{bc}(w)-U^{bc}(y)\Big]\rho^{c}(y)\notag\\
&=-"""+PF_DN+r"""\,\pi\log\frac{k^{+}}{\Lambda}\int_{w,w'}"""+PH_WW+r"""\,f^{abd}\,\mathcal{A}^{a}(w')\,\bar{Q}^{d}\,\mathcal{A}^{b}(w)\,,\label{eq:C-Coul}
\end{align}
\endgroup
with $\mathcal{A}^{b}(w)\equiv\int_{y}\mathcal{K}(y-w)\big[U^{bc}(w)-U^{bc}(y)\big]\rho^{c}(y)$ (the LO amplitude, a vector) and $\bar{Q}^{d}\equiv\int_{x}U^{de}(x)\rho^{e}(x)$ (the colour charge of the projectile after the collision): a forward (Coulomb) scattering phase of the observed gluon. Its complex conjugate is $\propto\int e^{+ik\cdot(w-w')}f^{abd}\mathcal{A}^{b}(w)\bar{Q}^{d}\mathcal{A}^{a}(w')$ (hermitian charges, real $U$); renaming $w\leftrightarrow w'$, $a\leftrightarrow b$ (which keeps the order of the charges: $\mathcal{A}(w')$, $\bar{Q}$, $\mathcal{A}(w)$) and $f^{bad}=-f^{abd}$ give $-T_{\Lambda}$. \textbf{Exact cancellation} of {\color{gray}$T_{\Lambda}+{\rm c.c.}=0$} (operator identity; checked with $SU(3)$ charges on three sites: $|T(w,w')^{\dagger}+T(w',w)|/|T|=9\times10^{-16}$, Sec.~\ref{sec:F}). It exists only in the $i\epsilon$ prescription; in the principal-value prescription there is no such term.
""")
open('partC.tex', 'w').write(''.join(o))
print('ok')
