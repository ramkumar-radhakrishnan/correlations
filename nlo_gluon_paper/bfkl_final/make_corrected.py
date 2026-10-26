"""bfkl_corrected.tex: your bfkl.tex with the two changes, red = remove, blue = add.
 (A) Group I Row III, d^3N_3a: the first term (kernel X'.Y/(z-w)^2) -> the D term (fix F3, as in Evolution_terms.tex).
 (B) Group III Rows III (both parts), V, VI: the log(vee/k+) terms from the region y' -> z.
 (C) the one-line notes under the rows."""
import os
HERE = os.path.dirname(os.path.abspath(__file__))
s = open(os.path.join(HERE, 'input', 'bfkl.tex')).read()


def rep(old, new, count=1):
    global s
    assert s.count(old) == count, (s.count(old), old[:80])
    s = s.replace(old, new)


rep(r'\numberwithin{equation}{section}',
    r'''\numberwithin{equation}{section}
\newcommand{\new}[1]{{\color{blue}#1}}
\newcommand{\old}[1]{{\color{red}#1}}
\newcommand{\nc}{\color{blue}}
\newcommand{\rc}{\color{red}}
\allowdisplaybreaks''')
rep(r'\title{BFKL evolution of the projectile}', r'\title{BFKL evolution of the projectile\\[2pt] \large \new{(corrected: red = remove, blue = add)}}')

# (A) Group I Row III, d^3N_3a, first term
CL = r"\left[U^{ab'}(x')\,\rho^{b'}(x')-U^{ab'}(w')\, \rho^{b'}(x')\right]"
old_A = (r'''\Bigg\{\notag \\ &  -\frac{ig^{4}f^{abd}}{8\pi^{5}}\,\frac{1}{k^{+}}\, \log{\left[\frac{k^{+}}{\vee}\right]}\,\int_{x,y,z,x'}
\frac{(x'-w')\cdot(y-z)}{(x'-w')^{2}(y-z)^{2}(z-w)^{2}}\notag \\ & ''' + CL +
         r'''\times \left[U^{de}(x)\,U^{bc}(z)\, \rho^{e}(x)\,\rho^{c}(y)-U^{de}(x)\,U^{bc}(y)\,\rho^{e}(x)\,\rho^{c}(y)\right]\notag \\ & -\frac{ig^{4}f^{abd}}{8\pi^{5}}''')
new_A = (r'''\Bigg\{\notag \\ & \rc -\frac{ig^{4}f^{abd}}{8\pi^{5}}\,\frac{1}{k^{+}}\, \log{\left[\frac{k^{+}}{\vee}\right]}\,\int_{x,y,z,x'}
\frac{(x'-w')\cdot(y-z)}{(x'-w')^{2}(y-z)^{2}(z-w)^{2}}\notag \\ & \rc ''' + CL +
         r'''\times \left[U^{de}(x)\,U^{bc}(z)\, \rho^{e}(x)\,\rho^{c}(y)-U^{de}(x)\,U^{bc}(y)\,\rho^{e}(x)\,\rho^{c}(y)\right]\notag \\
& \nc -\frac{ig^{4}f^{abd}}{8\pi^{5}}\,\frac{1}{k^{+}}\, \log{\left[\frac{k^{+}}{\vee}\right]}\,\int_{x,y,z,x'}
\frac{(x'-w')\cdot(y-w)\,(x-z)\cdot(w-z)}{(x'-w')^{2}\,(y-w)^{2}\,(x-z)^{2}\,(w-z)^{2}}\notag \\ & \nc ''' + CL +
         r'''\times \left[U^{de}(x)\,U^{bc}(w)\, \rho^{e}(x)\,\rho^{c}(y)-U^{de}(x)\,U^{bc}(y)\,\rho^{e}(x)\,\rho^{c}(y)\right]\notag \\ & -\frac{ig^{4}f^{abd}}{8\pi^{5}}''')
rep(old_A, new_A)

NOTE_A = r'''
\new{\textbf{Row III.} Red: the term with $\frac{(x'-w')\cdot(y-z)}{(x'-w')^2(y-z)^2(z-w)^2}$ is what is left of the edge kernel $\frac{(x-z)^2+(x-w)^2}{[(x-z)^2-(x-w)^2]^2}$ of part~2. It is not integrable at $z\to w$, and with it the three-$\rho$ terms do not cancel. Blue: the $\mathcal D$ term that replaces the edge kernel (fix F3). Your \texttt{Evolution\_terms.tex} already has $\mathcal D$ in Group~I Row~III (the $\mathcal K(x-z)\cdot\mathcal K(w-z)$ terms of $d^3N_{2}$); this is its $\log\frac{\vee}{k^+}$ partner, from $\log\frac{k^+}{\Lambda}=\log\frac{\vee}{\Lambda}-\log\frac{\vee}{k^+}$. The other terms of $d^3N_2$, $d^3N_{3a}$ add up to $+\frac{ig^4f^{abd}}{8\pi^5}\frac{1}{k^+}\log\frac{\vee}{k^+}\frac{(x'-w')\cdot(z-w)\,(y-z)\cdot(x-z)}{(x'-w')^2(z-w)^2(y-z)^2(x-z)^2}$ (unchanged).}
'''
rep(r'''\subsection{Row IV}
\begin{align}
    \frac{d^{3}N_{\rm NLO}^{\rm}}{d^{3}k}  & = -\, \frac{1}{(2\pi)^{3}}\,\frac{ig^{4}}{16\pi^{5}}''',
    NOTE_A + r'''\subsection{Row IV}
\begin{align}
    \frac{d^{3}N_{\rm NLO}^{\rm}}{d^{3}k}  & = -\, \frac{1}{(2\pi)^{3}}\,\frac{ig^{4}}{16\pi^{5}}''')

# (B) Group III: y' -> z terms
KPC = r'''\frac{(y'-z)\cdot (y-w)\,(x'-z)\cdot (x-z)}{(y'-z)^{2}\,(y-w)^{2}\,(x'-z)^{2}\,(x-z)^{2}}'''
LC = r'''\left[U^{e'c'}(z)\,U^{d'b}(z)\,\rho^{e'}(x')-U^{bd'}(z)\,U^{c'e'}(x')\,\rho^{e'}(x')\right]'''

III1_old = r'''\left[U^{e'c'}(y')\,U^{d'b}(z)\,\rho^{e'}(x')-U^{bd'}(z)\,U^{c'e'}(x')\,\rho^{e'}(x')\right]\times \left[U^{ac}(w)\, \{\rho^{c}(y),\rho^{b}(x)\}-U^{bc}(z)\,U^{ad}(y)\,U^{ce}(x)\,\{\rho^{d}(y),\rho^{e}(x)\}\right]
\end{align}'''
III1_new = r'''\left[U^{e'c'}(y')\,U^{d'b}(z)\,\rho^{e'}(x')-U^{bd'}(z)\,U^{c'e'}(x')\,\rho^{e'}(x')\right]\times \left[U^{ac}(w)\, \{\rho^{c}(y),\rho^{b}(x)\}-U^{bc}(z)\,U^{ad}(y)\,U^{ce}(x)\,\{\rho^{d}(y),\rho^{e}(x)\}\right]\notag \\
& \nc -\frac{1}{(2\pi)^{3}}\, \frac{ig^{4}\,f^{c'd'a}}{16\pi^{5}}\, \frac{1}{k^{+}}\,\log{\left[\frac{\vee}{k^{+}}\right]}\,\int_{w,x',y',z,y,x}\, e^{-ik\cdot (y'-w)}\,''' + KPC + r'''\notag \\
& \nc ''' + LC + r'''\times \left[U^{ac}(w)\, \{\rho^{c}(y),\rho^{b}(x)\}-U^{bc}(z)\,U^{ad}(y)\,U^{ce}(x)\,\{\rho^{d}(y),\rho^{e}(x)\}\right]
\end{align}'''
rep(III1_old, III1_new)

III2_old = r'''& \left[f^{cbd}\,U^{ac}(w)\,\rho^{d}(x) -f^{acd}\,U^{bc}(z)\,U^{de}(x)\,\rho^{e}(x)\right]
\end{align}'''
III2_new = r'''& \left[f^{cbd}\,U^{ac}(w)\,\rho^{d}(x) -f^{acd}\,U^{bc}(z)\,U^{de}(x)\,\rho^{e}(x)\right]\notag \\
& \nc -\frac{1}{(2\pi)^{3}}\,\frac{g^{4}}{8\pi^{5}}\,\frac{1}{k^{+}}\,\log{\left[\frac{\vee}{k^{+}}\right]}\,\int_{w,x',y',z,x}\, e^{-ik\cdot (y'-w)}\,\bigg[\frac{(x'-z)\cdot (x-z)\,(y'-z)\cdot (z-w)}{(x'-z)^{2}\,(x-z)^{2}\,(y'-z)^{2}\,(z-w)^{2}}\notag \\
& \nc -\frac{(x'-z)\cdot (x-z)\,(y'-z)\cdot (x-w)}{2\,(x'-z)^{2}\,(x-z)^{2}\,(y'-z)^{2}\,(x-w)^{2}}\bigg]\, \left[f^{c'd'a}\,U^{e'c'}(z)\,U^{d'b}(z)\,\rho^{e'}(x')-f^{c'd'a}\, U^{bd'}(z)\,U^{c'e'}(x')\,\rho^{e'}(x')\right]\times\notag \\
& \nc \left[f^{cbd}\,U^{ac}(w)\,\rho^{d}(x) -f^{acd}\,U^{bc}(z)\,U^{de}(x)\,\rho^{e}(x)\right]
\end{align}'''
rep(III2_old, III2_new)

V_old = r'''\left[U^{ce}(x)\,U^{ad}(y)\,\rho^{e}(x)\rho^{d}(y)-U^{ad}(w)\,U^{ce}(x)\,\rho^{e}(x)\rho^{d}(y)\right]
\end{align}'''
V_new = r'''\left[U^{ce}(x)\,U^{ad}(y)\,\rho^{e}(x)\rho^{d}(y)-U^{ad}(w)\,U^{ce}(x)\,\rho^{e}(x)\rho^{d}(y)\right]\notag \\
& \nc -\frac{1}{(2\pi)^{3}}\, \frac{ig^{4}\,f^{c'ca}}{8\pi^{5}}\, \frac{1}{k^{+}}\,\log{\left[\frac{\vee}{k^{+}}\right]}\,\int_{w,y',z,y,x,x'}\, e^{-ik\cdot (y'-w)}\,''' + KPC + r'''\notag \\
& \nc \left[U^{e'c'}(z)\,\rho^{e'}(x')-U^{c'e'}(x')\,\rho^{e'}(x')\right]\times \left[U^{ce}(x)\,U^{ad}(y)\,\rho^{e}(x)\rho^{d}(y)-U^{ad}(w)\,U^{ce}(x)\,\rho^{e}(x)\rho^{d}(y)\right]
\end{align}'''
rep(V_old, V_new)

VI_old = r'''\left[U^{bc}(z)\,U^{ad}(y)\,U^{ce}(x)\,\rho^{d}(y)\,\rho^{e}(x)-U^{ad}(y)\,\rho^{d}(y)\,\rho^{b}(x)\right]
\end{align}'''
VI_new = r'''\left[U^{bc}(z)\,U^{ad}(y)\,U^{ce}(x)\,\rho^{d}(y)\,\rho^{e}(x)-U^{ad}(y)\,\rho^{d}(y)\,\rho^{b}(x)\right]\notag \\
& \nc -\frac{1}{(2\pi)^{3}}\, \frac{ig^{4}\,f^{c'd'a}}{8\pi^{5}}\, \frac{1}{k^{+}}\,\log{\left[\frac{\vee}{k^{+}}\right]}\,\int_{w,y',z,y,x,x'}\, e^{-ik\cdot (y'-w)}\,''' + KPC + r'''\notag \\
& \nc ''' + LC + r'''\times \left[U^{bc}(z)\,U^{ad}(y)\,U^{ce}(x)\,\rho^{d}(y)\,\rho^{e}(x)-U^{ad}(y)\,\rho^{d}(y)\,\rho^{b}(x)\right]
\end{align}'''
rep(VI_old, VI_new)

NOTE_B = r'''
\new{\textbf{Group III Rows III, V, VI (blue).} In \texttt{Finite\_terms\_1.tex} these rows have the phase $e^{-ik\cdot(y'-z)\frac{p^+}{k^+}}$ under $\int_0^{\vee-k^+}dp^+$. The terms whose $p^+$ integrand tends to a constant at large $p^+$, the $-\frac{1}{k^+}\delta_{im}\delta_{jk'}$ term of Rows III(1), V, VI and the $-\frac{p^+}{k^{+2}}\big[\dots\big]$ term of Row III part 2, give $\log\frac{\vee}{k^+}$ only after the $y'$ integral, from $y'\to z$. For large $p^+$,
\begin{equation*}
\int_{y'}e^{-ik\cdot(y'-w)}e^{-ik\cdot(y'-z)\frac{p^+}{k^+}}\,\frac{(y'-z)^m}{(y'-z)^2}\,G_m(y')\;\to\;\frac{k^+}{p^++k^+}\int_{y'}e^{-ik\cdot(y'-w)}\frac{(y'-z)^m}{(y'-z)^2}\,G_m(z),
\end{equation*}
so $\int_0^{\vee-k^+}\frac{dp^+}{k^+}$ gives $\log\frac{\vee}{k^+}$ times the blue terms: $y'\to z$ in everything except $\frac{(y'-z)^m}{(y'-z)^2}$. The measured gluon of the conjugate amplitude ($y'$) is emitted by the soft gluon $z$; these are the $\mathcal K(x-w)\cdot\mathcal K(w'-z)$ terms of your BFKL contribution. The same mechanism gave $H_4$ in Group~II Row~IV.}
'''
rep(r'''\subsection{Row VI}''', r'''\subsection{Row VI}''')
rep(r'''\end{document}''', NOTE_B + r'''
\new{With these changes the three-$\rho$ terms cancel and the two-$\rho$ terms add up to $\frac{1}{(2\pi)^3}\frac{g^4}{16\pi^5}\frac{1}{k^+}\log\big[\frac{\vee}{k^+}\big]\int_{x,x',z,w,w'}e^{-ik\cdot(w'-w)}\big[\mathcal B-\mathcal J\big]$ ($\mathcal B$: your BFKL contribution; $\mathcal J$: your JIMWLK contribution). See \texttt{BFKL\_final.pdf}.}
\end{document}''')
open(os.path.join(HERE, 'bfkl_corrected.tex'), 'w').write(s)
print('written bfkl_corrected.tex')
