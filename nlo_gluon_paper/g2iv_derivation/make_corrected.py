r"""Your derivation of Group II Row IV (input/derivation_body.tex) with the corrections, red = removed, blue = new.
  python3 make_corrected.py -> G2IV_derivation_corrected.tex     (every edit asserted)"""
import os

HERE = os.path.dirname(os.path.abspath(__file__))
t = open(os.path.join(HERE, 'input', 'derivation_body.tex')).read()


def rep(old, new, n=1):
    global t
    assert t.count(old) == n, (old[:70], t.count(old))
    t = t.replace(old, new)


# C1 + C4: the log(vee/k+) H4 line has no sign; it is +, and with d_perp/2
rep(r"&   \frac{1}{(2\pi)^{3}}\,\frac{g^{4}\,N_{c}}{8\pi^{5}}\,\frac{1}{k^{+}}\,\log{\left(\frac{\vee}{k^{+}}\right)}",
    r"&   \new{+\,\frac{d_{\perp}}{2}}\,\frac{1}{(2\pi)^{3}}\,\frac{g^{4}\,N_{c}}{8\pi^{5}}\,\frac{1}{k^{+}}\,\log{\left(\frac{\vee}{k^{+}}\right)}")
# C4: the (1 - k+/vee) H4 line with d_perp/2
rep(r"&   -\frac{1}{(2\pi)^{3}}\,\frac{g^{4}\,N_{c}}{8\pi^{5}}\,\frac{1}{k^{+}}\,\left(1-\frac{k^{+}}{\vee}\right)",
    r"&   -\new{\frac{d_{\perp}}{2}}\,\frac{1}{(2\pi)^{3}}\,\frac{g^{4}\,N_{c}}{8\pi^{5}}\,\frac{1}{k^{+}}\,\left(1-\frac{k^{+}}{\vee}\right)")
# C2: wording
rep(r"Since, we have handled the UV divergence, now we can safely use plus prescription to integrate wrt $\xi$.",
    r"Since \old{we have handled the UV divergence}\new{we have isolated the collinear (infrared) divergence: the $z$ integral "
    r"diverges at large $|z-y|$, not at $z\to y$, where its integrand goes only like $1/|z-y|$}, now we can safely use "
    r"plus prescription to integrate wrt $\xi$. \new{Here $\int\frac{d\xi}{\xi^2}\frac{\xi}{(1-\xi)_+}G\equiv\int\frac{d\xi}"
    r"{\xi(1-\xi)}[G(\xi)-G(1)]$ (your prescription). With it, $\int_{k^+/\vee}^{k^+/(k^++\Lambda)}\frac{d\xi}{\xi(1-\xi)}=\log"
    r"\frac{k^+}{\Lambda}+\log(\frac{\vee}{k^+}-1)$ exactly, which is the blue line below.}")
# C3: Lambda -> 0 in the upper limits (annotate after the equation with the plus prescription)
rep(r"\end{align} We can split this as :",
    r"\end{align}" "\n"
    r"\new{\noindent In the first line the upper limit $\frac{k^+}{k^++\Lambda}$ of the plus-prescribed integral is set to $1$. "
    r"This drops $\int_{k^+/(k^++\Lambda)}^{1}\frac{d\xi}{\xi(1-\xi)}[G(\xi)-G(1)]=O(\frac{\Lambda}{k^+})$, i.e.\ it is the "
    r"limit $\Lambda\to0$. The same holds for the upper limit of the linear-polarization line, $\frac{k^+}{k^++\Lambda}\to1$, in the "
    r"next equation.}" "\n\n" r"We can split this as :")
# C3: mark the LP upper limit in the split
old_lp = (r"& +\,\frac{g^{2}}{4\pi^{2}}\, \frac{1}{(2\pi)^{3}}\,\frac{g^{2}}{2\pi^{2}}\, \frac{N_{c}}{k^{+}}\,\int_{\frac{k^{+}}{\vee}}^{1}\,d\xi")
rep(old_lp, old_lp.replace(r"\int_{\frac{k^{+}}{\vee}}^{1}",
                           r"\int_{\frac{k^{+}}{\vee}}^{\old{\frac{k^{+}}{k^{+}+\Lambda}}\to\new{1}}"))
# closing remark after the split
assert t.rstrip().endswith(r"\end{align}")
t = t.rstrip() + "\n" + (
    r"\new{\noindent\textbf{Check of the split.} The two lines with $\pm\log\frac{\vee}{k^+}$ times the bracket cancel. "
    r"Your $y'=y$ subtraction line gives, after the $y'$ integral, $-\frac{d_\perp}{2}\big[\log\frac{\vee}{k^+}"
    r"-(1-\frac{k^+}{\vee})\big]H_4$. Together with the two $H_4$ lines, $-\frac{d_\perp}{2}(1-\frac{k^+}{\vee})H_4$ and "
    r"$+\frac{d_\perp}{2}\log\frac{\vee}{k^+}H_4$, it adds up to zero. Hence the last equation equals the one before it, "
    r"exactly at $d_\perp=2-\epsilon$. Here $H_4$ is the $\int_{x,x',z,w,w'}$ structure of those lines.}" "\n")

head = r"""\documentclass[10pt]{article}
\usepackage[a4paper,margin=0.3in]{geometry}
\usepackage{amsmath,amssymb,amsfonts,mathtools}
\usepackage[dvipsnames]{xcolor}
\usepackage[colorlinks=true,linkcolor=blue]{hyperref}
\usepackage{etoolbox}
\AtBeginEnvironment{align}{\small}
\allowdisplaybreaks
\setlength{\emergencystretch}{3em}
\def\mf#1{\mathbb{#1}}
\newcommand{\new}[1]{{\color{blue}#1}}
\newcommand{\old}[1]{{\color{red}#1}}
\title{Group II Row IV from $\mathfrak C_1$: your derivation, corrected}
\date{}
\begin{document}
\maketitle
\noindent\new{\textbf{Corrections} (blue = new, red = removed). Every other step is correct as written; the checks are in \texttt{checks.py}.}
\begin{enumerate}
\item[\new{C1}] \new{Last equation, the $\log\frac{\vee}{k^+}\,H_4$ line: the sign at the start of the line was missing. It is $+$; without it the split does not add up, and LaTeX prints the line as a product with the previous one.}
\item[\new{C2}] \new{The $\frac1\epsilon$ is the collinear (infrared) divergence of the $z$ integral (large $|z-y|$), not an ultraviolet one. Only the wording changes.}
\item[\new{C3}] \new{The upper limit $\frac{k^+}{k^++\Lambda}\to1$, in the plus-prescribed integral and in the linear-polarization line, is the limit $\Lambda\to0$: it drops $O(\frac{\Lambda}{k^+})$.}
\item[\new{C4}] \new{With $d_\perp=2-\epsilon$ kept: the two $H_4$ lines carry $\frac{d_\perp}{2}$, like the $y'=y$ subtraction line. Then the last equation equals the one before it exactly.}
\end{enumerate}
\new{\noindent\textbf{Checked and correct:}
\begin{itemize}
\item the prefactor $4\cdot\frac{g^2}{2\pi^2P}\cdot\frac{g^2\,p^+k^+}{8\pi^2P^3}=\frac{g^4p^+k^+}{4\pi^4P^4}$ ($P=p^++k^+$);
\item the $w,w'$ integrals: Jacobian $(P/k^+)^2$ each, $w-z=\frac{P}{k^+}(y-z)$, phase $e^{-ik(y'-y)P/k^+}$;
\item the $i,j$ contraction and the bracket with $\frac{1}{k^+}$ (symbolic, and index by index in $d=2,3$);
\item the Wilson-line reduction to $N_c[U(y')-U(x')][U(y)-U(x)]$;
\item the $z$ integral: trace $-\pi[\frac2\epsilon+2\gamma_E+\log\frac{(y'-y)^2\mu^2}{4}]$, traceless part $-\pi\frac{(y'-y)^m(y'-y)^{m'}}{(y'-y)^2}$;
\item the drop-out of the $\frac{2}{k^+}$ terms;
\item the change of variables $\xi=\frac{k^+}{p^++k^+}$, with its limits and sign;
\item the blue $\log\frac{k^+}{\Lambda}+\log(\frac{\vee}{k^+}-1)$, which is exact with your plus prescription.
\end{itemize}}

"""
tail = "\n\\end{document}\n"
open(os.path.join(HERE, 'G2IV_derivation_corrected.tex'), 'w').write(head + t + tail)
print('written')
