"""Your Finite_terms.tex with the changes: red = remove, blue = add.  Output: Finite_terms_corrected.tex
Every replacement is an exact substring of your file and must occur the stated number of times."""
import os

HERE = os.path.dirname(os.path.abspath(__file__))
T = open(os.path.join(HERE, 'input', 'Finite_terms.tex')).read()
L = T.split('\n')


def rep(i, old, new, count=1):
    """replace in line i (1-based) of your file"""
    n = L[i - 1].count(old)
    assert n == count, (i, old[:60], n)
    L[i - 1] = L[i - 1].replace(old, new)


def add_after(i, text):
    L[i - 1] = L[i - 1] + '\n' + text


PP = r'\left[\frac{1}{p^{+}}\right]_{+}'
PPK = r'\new{\left[\frac{1}{p^{+}}\right]_{+}^{(k^{+})}}'
LV = r'\log{\left[\frac{\vee}{k^{+}}\right]}'
NOTE_PP = (r'\new{\textbf{Change C2.} $\left[\frac{1}{p^+}\right]_+\to\left[\frac{1}{p^+}\right]_+^{(k^+)}$ (kept '
           r'unevaluated). $\int_0^{\vee-k^+}dp^+\left[\frac{1}{p^+}\right]_+F=\int_0^{\vee-k^+}dp^+\left[\frac{1}{p^+}'
           r'\right]_+^{(k^+)}F-F(0)\log\frac{\vee}{k^+}$, and $-F(0)\log\frac{\vee}{k^+}$ is in the BFKL list. '
           r'With $\left[\frac{1}{p^+}\right]_+$ it is counted twice.')
NOTE_B = (r" \textbf{Change C3.} The last (blue) line takes out the term (B) of this row (the region $y'\to z$), which is in "
          r"the BFKL list (\texttt{BFKL\_final.pdf}, terms (B)). After $\int_{y'}$, at large $p^+$, the "
          r"$-\frac{1}{k^+}\delta_{im}\delta_{jk'}$ term equals minus the $p^+$ integrand of the blue line, up to terms that "
          r"fall like $\frac{1}{p^{+2}}$ (Fourier lemma). Each alone grows like $\log\frac{\vee}{k^+}$; their sum is finite for "
          r"$\vee\to\infty$. The $p^+$ integral of the blue line is $-$(term (B)).}")
NOTE_B2 = NOTE_B.replace(r"$-\frac{1}{k^+}\delta_{im}\delta_{jk'}$ term",
                         r"$-\frac{p^+}{k^{+2}}\big[\dots\big]$ term ($\to-\frac{1}{k^{+2}(x-z)^2}\big[\dots\big]$)")
assert NOTE_B2 != NOTE_B

# ---------------------------------------------------------------- preamble, key
rep(17, r'\usepackage{comment}',
    r'\usepackage{comment}' + '\n' + r'\newcommand{\new}[1]{{\color{blue}#1}}' + '\n'
    + r'\newcommand{\old}[1]{{\color{red}#1}}' + '\n' + r'\newcommand{\nc}{\color{blue}}' + '\n'
    + r'\newcommand{\rc}{\color{red}}' + '\n' + r'\allowdisplaybreaks')
rep(30, r'\title{Finite terms}', r'\title{Finite terms \new{(checked: red = remove, blue = add)}}')
rep(35, r'\maketitle', r'\maketitle' + '\n' + r'\sloppy' + '\n' + r'''\new{\noindent\textbf{Key.} Red: remove. Blue: add. The checks and the derivations are in \texttt{Finite\_terms\_final.pdf}. Every row not marked is correct as written: it equals the row of \texttt{Finite\_terms\_1.tex} minus its BFKL ($\log\frac{\vee}{k^+}$) and JIMWLK ($\log\frac{\vee}{\Lambda}$) terms, term by term. The changes are:
\begin{itemize}
\item[C1] Group II Row IV: the line $-\log\frac{\vee}{k^+}H_4$ is missing. $H_4$ is in the BFKL list, and the $\xi$ integral of line 1 also contains it.
\item[C2] Group III Rows III (both parts), V, VI: $\left[\frac{1}{p^+}\right]_+\to\left[\frac{1}{p^+}\right]_+^{(k^+)}$, where
$\int_0^{\vee-k^+}dp^+\left[\frac{1}{p^+}\right]_+^{(k^+)}F(p^+)\equiv\int_0^{\vee-k^+}dp^+\,\frac{F(p^+)-F(0)\frac{k^+}{p^++k^+}}{p^+}$.
\item[C3] Group III Rows III (both parts), V, VI: subtract the terms (B) (the region $y'\to z$). They are in the BFKL list.
\item[C4] Group I Row III: the first term of $d^3N_{3a}$ and the last term of $d^3N_{3b}$ (edge kernel $\frac{(x'-w')\cdot(y-z)}{(x'-w')^2(y-z)^2(z-w)^2}$) are removed. The edge kernel is replaced by $\mathcal D$ (fix (A)). The last term of $d^3N_{3b}$ is not integrable at $z\to w$.
\item[C5] Group III Row III part 2: the measure is $\int_{w,x',y',z,x}$.
\end{itemize}
The DGLAP and running-coupling terms sit in Group I Row II part 1, Group I Row V 2b UV and Group II Row IV (notes under each).}''')

# ---------------------------------------------------------------- Group I Row II part 1
add_after(42, r'''\new{\textbf{No change.} In the final result: the $\frac{11}{3}$ term is running coupling and cancels against $-\frac{11}{6}$ of Row V 2b UV. The constant $-\frac{67}{9}+\frac{\pi^2}{3}$ gives $\frac{\alpha_s}{2\pi}N_c\big(\frac{67}{9}-\frac{\pi^2}{3}\big){\rm LO}$. The pole bracket times $\log(1-\frac{k^+}{\vee})$ is Q1. $\log^2+\log\log$ is $V_\vee$, a finite term. (\texttt{Finite\_terms\_final.pdf}, Secs.~3, 4.)}''')
add_after(55, r'\new{\textbf{No change.}}')

# ---------------------------------------------------------------- Group I Row III: C4
rep(64, r'\Bigg\{\notag \\ &  -\frac{ig^{4}f^{abd}}{8\pi^{5}}', r'\Bigg\{\notag \\ & \rc -\frac{ig^{4}f^{abd}}{8\pi^{5}}')
rep(65, r'(z-w)^{2}}\notag \\ & \left[U^{ab', r'(z-w)^{2}}\notag \\ & \rc \left[U^{ab')
rep(88, r'\notag \\ & -\frac{ig^{4}f^{abd}}{16\pi^{5}}', r'\notag \\ & \rc -\frac{ig^{4}f^{abd}}{16\pi^{5}}')
rep(89, r'\,(z-w)^{2}}\notag \\ & \left[U^{ab', r'\,(z-w)^{2}}\notag \\ & \rc \left[U^{ab')
rep(89, r'\right]\Bigg\}', r'\right]\color{black}\Bigg\}')
add_after(94, r'''\new{\textbf{Change C4.} The two red terms have the edge kernel $\frac{(x'-w')\cdot(y-z)}{(x'-w')^{2}(y-z)^{2}(z-w)^{2}}$. Their $\log\frac{\vee}{k^+}$ and $\log\frac{\vee}{\Lambda}$ partners are replaced by $\mathcal D$ (fix (A), \texttt{BFKL\_final.pdf}; fix F3 of \texttt{Evolution\_terms.tex}). The edge kernel therefore does not belong to the row, and its finite terms go with it.
\begin{itemize}
\item The first term of $d^3N_{3a}$ vanishes for $\Lambda\to0$.
\item The last term of $d^3N_{3b}$ is, for $\Lambda\to0$, $\log\frac{(z-w)^2}{(x-w)^2}$ times the edge kernel. Near $z\to w$ it behaves like $\frac{\log(1/|z-w|)}{(z-w)^2}$, so $\int_z$ diverges. Numerically: the angular average of the row behaves like $\frac{\log(1/r)}{r^2}$, $r=|z-w|$; without this term it is integrable (\texttt{edge\_check.json}).
\item In its place belongs the finite part of the $\mathcal D$ region. It needs the $p^+$ integral of this region with the kernel of $\mathcal D$. That integral is not in your notes or in mine (open, Q3).
\end{itemize}
All other terms of Row III are finite. On the circle $(x-z)^2=(x-w)^2$ the $\frac{1}{[(x-z)^2-(x-w)^2]^n}$ cancel among them (\texttt{edge\_check.json}).}''')

# ---------------------------------------------------------------- Group I Row V
add_after(113, r'''\new{\textbf{No change.} In the final result: $-\frac{11}{6}$ cancels the $\frac{11}{3}$ of Row II part 1 (running coupling). $2\log(1-\frac{\Lambda}{k^+})\to0$. Nothing of this term is left among the finite terms.}''')

# ---------------------------------------------------------------- Group II Row IV: C1
add_after(153, r'''\notag \\ & \nc -\,\frac{1}{(2\pi)^{3}}\,\frac{g^{4}\,N_{c}}{8\pi^{5}}\,\frac{1}{k^{+}}\,\log{\left[\frac{\vee}{k^{+}}\right]}\,\int_{x,x',z,w,w'}\,e^{-ik\cdot (w'-w)}\,\frac{(w-z)\cdot(w'-z)}{(w-z)^{2}\,(w'-z)^{2}}\,\frac{(x-z)\cdot(x'-z)}{(x-z)^{2}\,(x'-z)^{2}}\notag \\ & \nc \left[U^{e'c}(z) - U^{e'c}(x')\right] \times\, \left[U^{ec}(z) - U^{ce}(x)\right]\, \rho^{e'}(x') \rho^{e}(x)''')
add_after(154, r'''\new{\textbf{Change C1.} At $\xi\to\frac{k^+}{\vee}$ the phase $e^{-ik\cdot(y'-y)\frac1\xi}$ of line 1 probes $|y'-y|\sim\frac{k^+}{\vee|k|}$. The large-momentum tail of $\log(y'-y)^2$ then gives $\log\frac{\vee}{k^+}H_4$, with $H_4$ the bracket of the blue line ($y,y'\to z$: the measured gluon and its conjugate are emitted by the soft gluon $z$). $\log\frac{\vee}{k^+}H_4$ is one of the BFKL terms. The blue line takes it out here, so that it is not counted twice.
Lines 1 and 2 with the blue line make up the DGLAP term; line 3 (linear polarization) is finite (\texttt{Finite\_terms\_final.pdf}, Sec.~4).}''')

# ---------------------------------------------------------------- Group III: C2, C3, C5
CB_III1 = (r"\left[U^{e'c'}(z)\,U^{d'b}(z)\,\rho^{e'}(x')-U^{bd'}(z)\,U^{c'e'}(x')\,\rho^{e'}(x')\right]")
KB = (r"\frac{(y'-z)\cdot (y-w)\,(x'-z)\cdot (x-z)}{(y'-z)^{2}\,(y-w)^{2}\,(x'-z)^{2}\,(x-z)^{2}}")
rep(167, PP, PPK)
add_after(167, r'''\notag \\ & \nc +\,\frac{1}{(2\pi)^{3}}\, \frac{ig^{4}\,f^{c'd'a}}{8\pi^{4}}\, \frac{1}{k^{+}}\,\int_{0}^{\vee-k^{+}}\, \frac{dp^{+}}{2\pi}\,\frac{1}{p^{+}+k^{+}}\, \int_{w,x',y',z,y,x}\, e^{-ik\cdot (y'-w)}\,''' + KB + r'''\notag \\ & \nc ''' + CB_III1 + r'''\times\left[U^{ac}(w)\, \{\rho^{c}(y),\rho^{b}(x)\}-U^{bc}(z)\,U^{ad}(y)\,U^{ce}(x)\,\{\rho^{d}(y),\rho^{e}(x)\}\right]''')
add_after(168, NOTE_PP + NOTE_B)

rep(171, PP, PPK)
rep(171, r"\int_{w,w',y',z,y,x}", r"\int_{\old{w,w',y',z,y,x}\new{w,x',y',z,x}}")
rep(172, r"\int_{w,y',z,y,x}", r"\int_{\old{w,y',z,y,x}\new{w,x',y',z,x}}")
add_after(173, r'''\notag \\ & \nc +\,\frac{1}{(2\pi)^{3}}\,\frac{g^{4}}{4\pi^{4}}\,\int_{0}^{\vee-k^{+}}\, \frac{dp^{+}}{2\pi}\,\frac{1}{k^{+}\,(p^{+}+k^{+})}\,\int_{w,x',y',z,x}\, e^{-ik\cdot (y'-w)}\,\frac{(x'-z)\cdot (x-z)}{(x'-z)^{2}\,(x-z)^{2}}\notag \\ & \nc \bigg[\frac{(y'-z)\cdot (z-w)}{(y'-z)^{2}\,(z-w)^{2}}-\frac{(y'-z)\cdot (x-w)}{2\,(y'-z)^{2}\,(x-w)^{2}}\bigg]\left[f^{c'd'a}\,U^{e'c'}(z)\,U^{d'b}(z)\,\rho^{e'}(x')-f^{c'd'a}\, U^{bd'}(z)\,U^{c'e'}(x')\,\rho^{e'}(x')\right]\notag \\ & \nc \times\left[f^{cbd}\,U^{ac}(w)\,\rho^{d}(x) -f^{acd}\,U^{bc}(z)\,U^{de}(x)\,\rho^{e}(x)\right]''')
add_after(174, NOTE_PP.replace(r'counted twice.', r"counted twice. Here $F(p^+)=\frac{e^{-ik\cdot(y'-z)p^+/k^+}}{p^+(x-z)^2+k^+(x-w)^2}$.")
          + NOTE_B2[:-1] + r" \textbf{Change C5.} Measure: the row has $\rho(x')\rho(x)$ and the kernels depend on $w,x',y',z,x$; there is no $y$ and no $w'$.}")

rep(179, PP, PPK)
add_after(179, r'''\notag \\ & \nc +\,\frac{1}{(2\pi)^{3}}\, \frac{ig^{4}\,f^{c'ca}}{4\pi^{4}}\, \frac{1}{k^{+}}\,\int_{0}^{\vee-k^{+}}\, \frac{dp^{+}}{2\pi}\,\frac{1}{p^{+}+k^{+}}\, \int_{w,y',z,y,x,x'}\, e^{-ik\cdot (y'-w)}\,''' + KB + r'''\notag \\ & \nc \left[U^{e'c'}(z)\,\rho^{e'}(x')-U^{c'e'}(x')\,\rho^{e'}(x')\right]\times\left[U^{ce}(x)\,U^{ad}(y)\,\rho^{e}(x)\rho^{d}(y)-U^{ad}(w)\,U^{ce}(x)\,\rho^{e}(x)\rho^{d}(y)\right]''')
add_after(180, NOTE_PP + NOTE_B)

rep(183, PP, PPK)
add_after(183, r'''\notag \\ & \nc +\,\frac{1}{(2\pi)^{3}}\, \frac{ig^{4}\,f^{c'd'a}}{4\pi^{4}}\, \frac{1}{k^{+}}\,\int_{0}^{\vee-k^{+}}\, \frac{dp^{+}}{2\pi}\,\frac{1}{p^{+}+k^{+}}\, \int_{w,y',z,y,x,x'}\, e^{-ik\cdot (y'-w)}\,''' + KB + r'''\notag \\ & \nc ''' + CB_III1 + r'''\times\left[U^{bc}(z)\,U^{ad}(y)\,U^{ce}(x)\,\rho^{d}(y)\,\rho^{e}(x)-U^{ad}(y)\,\rho^{d}(y)\,\rho^{b}(x)\right]''')
add_after(184, NOTE_PP + NOTE_B)

if __name__ == '__main__':
    out = '\n'.join(L)
    assert out.count(r'\left[\frac{1}{p^{+}}\right]_{+}\,') == 0
    open(os.path.join(HERE, 'Finite_terms_corrected.tex'), 'w').write(out)
    print('written Finite_terms_corrected.tex')
