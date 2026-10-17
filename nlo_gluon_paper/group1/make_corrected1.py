"""Writes Group_I_corrected.tex: input/Group_I.tex with every correction of Group_I_check.pdf marked in place
(removed text in red, new text in blue).  Every edit is asserted to match the expected number of times."""
import os

HERE = os.path.dirname(os.path.abspath(__file__))
s = open(os.path.join(HERE, 'input', 'Group_I.tex')).read()


def rep(old, new, count=1):
    global s
    n = s.count(old)
    assert n == count, (n, count, old[:80])
    s = s.replace(old, new)


R = lambda t: r'{\color{red}' + t + '}'
B = lambda t: r'{\color{blue}' + t + '}'

# ---------------------------------------------------------------- header
rep(r'\maketitle', r'''\maketitle
\begingroup\color{blue}\noindent Corrected copy. Removed text is in red, new text in blue. The reasons and the numerical checks are in \texttt{Group\_I\_check.pdf}; G1--G6 refer to its Sec.~1.\par\endgroup''')

# ---------------------------------------------------------------- Row II, part 1 (G6, G5a)
rep(r'\log{\left[1-\frac{k^{+}}{\vee}\right]^{2}}',
    R(r'\log{\left[1-\frac{k^{+}}{\vee}\right]^{2}}') + r'\,' + B(r'\log^{2}{\left[1-\frac{k^{+}}{\vee}\right]}'), count=2)
anchor = r'''\begin{align}
   \frac{d^{3}N_{2\,\rm NLO}^{\rm}}{d^{3}k} &  = -\,'''
rep(anchor, r'''\begingroup\color{blue}\noindent\textbf{Exact bracket (G5a; no limit $\Lambda\to0$).} For $\Lambda<k^{+}<\vee-\Lambda$, with $\delta=\Lambda/k^{+}$ and $P_{\epsilon}=-\frac{2}{\epsilon}-2\gamma-\log\frac{(x-w)^{2}\mu^{2}}{4}$, the big bracket above is exactly
\begin{align*}
&P_{\epsilon}\Big[\frac{11}{3}+\log\frac{\Lambda^{2}}{k^{+2}-\Lambda^{2}}+\log\Big(1-\frac{k^{+}}{\vee}\Big)-8\delta+2\delta^{2}-\frac{4}{3}\delta^{3}\Big]-\frac{67}{9}+\frac{\pi^{2}}{3}+\log^{2}\Big(1-\frac{k^{+}}{\vee}\Big)+\log\Big(1-\frac{k^{+}}{\vee}\Big)\log\frac{\vee}{k^{+}}+R_{\delta},\\
&R_{\delta}=\frac{46}{3}\delta-\frac{4}{3}\delta^{2}+\frac{8}{9}\delta^{3}+\Big(\frac{22}{3}-8\delta+2\delta^{2}-\frac{4}{3}\delta^{3}\Big)\log(1-\delta)-\Big(8\delta-2\delta^{2}+\frac{4}{3}\delta^{3}\Big)\log\delta-2\log\delta\log(1-\delta)\\
&\qquad-4\,\mathrm{Li}_{2}(\delta)-\log^{2}(1+\delta)+\log(1+\delta)\log\delta .
\end{align*}
The bracket above is its $\delta\to0$ limit. $P_{\epsilon}=\frac{1}{\pi}\int_{z}\big[K(x,w,z)-K(w,w,z)\big]$ exactly in $d_{\perp}=2-\epsilon$, $\overline{\rm MS}$: the same scheme as Row V, part 2b.\par\endgroup
''' + anchor)

# ---------------------------------------------------------------- Row II, part 2 (G2, G6)
rep(r'\frac{g^{4}\,f^{abc}}{64\pi^{6}}', r'\frac{g^{4}\,' + R(r'f^{abc}') + r'}{64\pi^{6}}', count=2)
rep(r'\int_{x,y}\, e^{ik\cdot (w-y)}', r'\int_{x,y' + B(r',x^{\prime}') + r'}\, e^{ik\cdot (w-y)}')
rep(r'\int_{x,y} e^{ik\cdot(w-y)}', r'\int_{x,y' + B(r',x^{\prime}') + r'} e^{ik\cdot(w-y)}')
rep(r'\log{\left[1-\frac{k^{+}}{\vee}\right]}\right]\Bigg]',
    r'\log{\left[1-\frac{k^{+}}{\vee}\right]}\right]' + B(r'-k^{i}\,\log{\left[\frac{\vee}{k^{+}}\right]}') + r'\Bigg]')
rep(r'U^{ab}(w)\,\left\{ \rho^{b}(\mathbf{x}),\,\rho^{c}(\mathbf{y})\right\}-U^{bd}(x)\, U^{ce}(y)',
    R(r'U^{ab}(w)') + B(r'U^{ah}(w)\,f^{hbc}') + r'\,\left\{ \rho^{b}(\mathbf{x}),\,\rho^{c}(\mathbf{y})\right\}-'
    + B(r'f^{abc}') + r'\,U^{bd}(x)\, U^{ce}(y)')
rep(r'U^{ab}(w)\left\{\rho^{b}(x),\rho^{c}(y)\right\}-U^{bd}(x)\,U^{ce}(y)',
    R(r'U^{ab}(w)') + B(r'U^{ah}(w)\,f^{hbc}') + r'\left\{\rho^{b}(x),\rho^{c}(y)\right\}-' + B(r'f^{abc}') + r'\,U^{bd}(x)\,U^{ce}(y)')

# ---------------------------------------------------------------- Row II, part 2a (G3)
anchor = r'''\begin{align}
\frac{d^{3}N_{2\,b}}{d^{3}k} & ='''
rep(anchor, r'''\begingroup\color{blue}\noindent\textbf{Operator ordering (G3).} \texttt{Calculation\_main.tex} drops $-[\mathbb{N}^{(2)},\mathbb{A}^{(1)}_{2}]$ in the unitarity relation (its correction E5). With the virtual factor ordered in $p^{+}$, part 2a is
\begin{align*}
\frac{d^{3}N_{2\,a}}{d^{3}k}={}&-\frac{1}{(2\pi)^{3}}\frac{i\,g^{4}}{32\pi^{5}}\frac{1}{k^{+}}\int_{w,w'}e^{-ik\cdot(w'-w)}\int_{x,y,x'}\frac{(x'-w')^{i}}{(x'-w')^{2}}\Big[U^{ab'}(x')-U^{ab'}(w')\Big]\rho^{b'}(x')\\
&\times\Big\{\Big[\log\Big(\frac{\vee}{\Lambda}\Big)\mathfrak{F}^{i}_{\beta}-\log\Big(\frac{\vee}{k^{+}}\Big)\mathfrak{F}^{i}\Big]U^{ah}(w)f^{hbc}\{\rho^{b}(x),\rho^{c}(y)\}-\Big[\log\Big(\frac{\vee}{\Lambda}\Big)-\log\Big(\frac{\vee}{k^{+}}\Big)\Big]\mathfrak{F}^{i}f^{abc}U^{bd}(x)U^{ce}(y)\{\rho^{d}(x),\rho^{e}(y)\}\Big\},
\end{align*}
where $\mathfrak{F}^{i}_{\beta}$ is the first line of $\mathfrak{F}^{i}$ (the two terms with $\int_{z}\frac{(w-z)\cdot(\cdots)}{\cdots}$). Only the three-$\rho$ term $\propto\log(\vee/\Lambda)$ of the $U^{ah}(w)$ part changes.\par\endgroup
''' + anchor)

# ---------------------------------------------------------------- Row II, part 2b: exact R (G5b)
anchor = r'\subsection{Row III'
rep(anchor, r'''\begingroup\color{blue}\noindent\textbf{Exact $R^{i}$ (G5b).} With $E^{i}=\frac{k^{2}p^{i}-k\cdot p\,k^{i}}{(k-p)^{2}}$, $C^{i}=\frac{(2p^{2}-2k\cdot p+k^{2})k^{i}+k^{2}p^{i}}{(k-p)^{2}}$, $D^{i}=\frac{p^{2}k^{i}+k^{2}p^{i}}{(k-p)^{2}}$ and $\Lambda<k^{+}<\vee-\Lambda$:
\begin{align*}
R^{i}_{\rm exact}={}&-2k^{i}\log\frac{p^{2}+\frac{\Lambda}{k^{+}}k^{2}}{k^{2}+\frac{k^{+}}{\vee}p^{2}}+C^{i}\log\Big[\frac{k^{2}}{(k-p)^{2}}+\frac{k^{+}}{\vee}\Big(1-\frac{k^{2}}{(k-p)^{2}}\Big)\Big]-D^{i}\log\Big(1-\frac{k^{+}}{\vee}\Big)\\
&+E^{i}\log\frac{k^{+}}{k^{+}-\Lambda}-C^{i}\log\Big[1+\frac{\Lambda}{k^{+}}\frac{k^{2}}{(k-p)^{2}}\Big]+k^{i}\log\Big(1+\frac{\Lambda}{k^{+}}\Big).
\end{align*}
Part 2a does not change.\par\endgroup
''' + anchor)

# ---------------------------------------------------------------- Row III (G1, G4)
rep(r'& +\frac{ig^{4}\, f^{abd}}{4\pi^{4}}\, \frac{1}{k^{+}}\,\int_{\Lambda + k^{+}}^{\vee}',
    r'& ' + R('+') + B('-') + r'\frac{ig^{4}\, f^{abd}}{4\pi^{4}}\, \frac{1}{k^{+}}\,\int_{\Lambda + k^{+}}^{\vee}')
rep(r'Integrating wrt $p^{+}$ and we write it in three parts:',
    r'''\begingroup\color{blue}\noindent\textbf{Sign (G1).} The non-instantaneous entry of $\mathbb{B}_{3}$ for $p^{+}>k^{+}+\Lambda$ needs a minus sign (correction E3 of \texttt{Calculation\_main.tex}). Hence the sign of the fourth line above, and of every line of $d^{3}N_{2\,\rm NLO}$ below.\par\endgroup
Integrating wrt $p^{+}$ and we write it in three parts:''')
rep(r'\frac{d^{3}N_{2\,\rm NLO}^{\rm}}{d^{3}k}  & = \frac{1}{(2\pi)^{3}}  \int_{w, w^\prime} e^{-ik(w^\prime-w)}\, \Bigg\{',
    r'\frac{d^{3}N_{2\,\rm NLO}^{\rm}}{d^{3}k}  & = ' + B('-') + r'\frac{1}{(2\pi)^{3}}  \int_{w, w^\prime} e^{-ik(w^\prime-w)}\, \Bigg\{')
rep(r'\subsection{Row IV', r'''\begingroup\color{blue}\noindent\textbf{Principal values (G4).} The integrand has a pole at $p^{+}_{*}=k^{+}(x-w)^{2}/(x-z)^{2}$. When it lies inside the range, the $p^{+}$ integral is a principal value, and every logarithm of a ratio of $k^{+}(x-w)^{2}-c\,(x-z)^{2}$ factors in $d^{3}N_{1}$, $d^{3}N_{2}$, $d^{3}N_{3a}$ must be read as $\log|\cdots|$.

\noindent\textbf{Green terms.} With the sign of $d^{3}N_{2}$ corrected, the four green terms cancel exactly:
$\big[\frac{1}{\Lambda}+\frac{1}{k^{+}-\vee}\big]+\big[\frac{1}{\Lambda}-\frac{1}{k^{+}-\Lambda}\big]-\big[\frac{1}{\Lambda}-\frac{1}{\vee-k^{+}}\big]-\big[\frac{1}{\Lambda}-\frac{1}{k^{+}-\Lambda}\big]=0$, all multiplying $\frac{ig^{4}f^{abd}}{8\pi^{5}}\frac{(x'-w')\cdot(y-z)}{(x'-w')^{2}(y-z)^{2}[(x-z)^{2}-(x-w)^{2}]}$ and the same colour factor. As written they leave $2\big[\frac{1}{\Lambda}-\frac{1}{\vee-k^{+}}\big]$.\par\endgroup
\subsection{Row IV''')

# ---------------------------------------------------------------- Row V (G5c, G5d, G6)
rep("x',y'w'}", "x',y'" + B(',') + "w'}", count=s.count("x',y'w'}"))
rep(r'N_{c}\,S^{ad}(z)', r'N_{c}\,' + R(r'S^{ad}(z)') + r'\,' + B(r'U^{ad}(z)'))
rep(r'\int_{0}^{k^{+}}\, \frac{dp^{+}}{2\pi}', R(r'\int_{0}^{k^{+}}') + B(r'\int_{\Lambda}^{k^{+}-\Lambda}') + r'\, \frac{dp^{+}}{2\pi}', count=2)
rep(r'\int_{0}^{k^{+}}\, \frac{dq^{+}}{2\pi}', R(r'\int_{0}^{k^{+}}') + B(r'\int_{\Lambda}^{k^{+}-\Lambda}') + r'\, \frac{dq^{+}}{2\pi}')
rep(r'\int^{0}_{k^{+}}', R(r'\int^{0}_{k^{+}}') + B(r'\int^{\Lambda}_{k^{+}-\Lambda}'))
note = (r'\begingroup\color{blue}\noindent\textbf{Limits (G5c).} The plus prescription is taken on $\Lambda<p^{+}<k^{+}-\Lambda$: '
        r'$\int_{\Lambda}^{k^{+}-\Lambda}dp^{+}\,[\frac{k^{+}}{p^{+}}]_{+}f\equiv\int_{\Lambda}^{k^{+}-\Lambda}dp^{+}\frac{k^{+}}{p^{+}}[f(p^{+})-f(0)]$. '
        r'Then the endpoint terms, with $\log\frac{k^{+}-\Lambda}{\Lambda}=\log\frac{\vee}{\Lambda}+\log\frac{k^{+}-\Lambda}{\vee}$, are exact. '
        r'With $\int_{0}^{k^{+}}$ the result is off by $O(\Lambda/k^{+})$.\par\endgroup' + '\n')
rep(r'Using $+$ prescription, we get', note + r'Using $+$ prescription, we get')
rep(r'We will handle the $p^{+}$, using $+$ prescription for the first part. We get the following:',
    note + r'We will handle the $p^{+}$, using $+$ prescription for the first part. We get the following:')
rep(r'Setting $\Lambda\rightarrow 0$, we get',
    R(r'Setting $\Lambda\rightarrow 0$, we get') + ' ' + B(r'(G5d) Keep the exact form above: the form below drops $\frac{4\Lambda}{k^{+}}-\frac{\Lambda^{2}}{k^{+2}}+\frac{2\Lambda^{3}}{3k^{+3}}$ but keeps $2\log(1-\frac{\Lambda}{k^{+}})$, which is of the same order.'))

open(os.path.join(HERE, 'Group_I_corrected.tex'), 'w').write(s)
print('wrote Group_I_corrected.tex')
