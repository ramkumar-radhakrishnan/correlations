"""Write Group_II_corrected.tex."""
L = open('input/Group_II.tex').read().split('\n')
B = r"{\color{blue}%s}"


def rep(n, old, new):
    s = L[n - 1]
    assert old in s, (n, old)
    L[n - 1] = s.replace(old, new, 1)


for n in (54, 58, 62):
    rep(n, r"\int_{z,w,w^{\prime}}", r"\int_{%s z,w,w^{\prime}}" % (B % "x,x',"))
rep(92, r"\int_{x,y,y^{\prime}}", r"\int_{x,%s,y,y^{\prime}}" % (B % "x'"))
rep(98, r"\int_{z,w,w^{\prime}}\,e^{-ik\,\left(y'-y\right)\, \frac{(k^{+}+p^{+})}{k{+}}}\,\int_{x,y,y^{\prime}}",
    r"\int_{%s}\,e^{-ik\,\left(y'-y\right)\, \frac{(k^{+}+p^{+})}{%s}}\,\int_{x,%s,y,y^{\prime}}" % (B % "z", B % "k^{+}", B % "x'"))
rep(104, r"\int_{z,w,w^{\prime}}", r"\int_{%s}" % (B % "z"))
rep(104, r"\int_{x,y,y^{\prime}}", r"\int_{x,%s,y,y^{\prime}}" % (B % "x'"))
# lines 106-116: the z integration and the fragmentation function -> red (removed), then the corrected version
start = [i for i, l in enumerate(L) if l.startswith('Integrating wrt $z$')][0]
end = [i for i, l in enumerate(L) if 'is the Fragmentation function' in l][0]
old = L[start:end + 1]
red = []
for l in old:
    if l.strip().startswith(r'\begin{align}'):
        red.append(r'\begingroup\color{red}' + l)
    elif l.strip().startswith(r'\end{align}'):
        red.append(l + r'\endgroup')
    else:
        red.append(l)
red[0] = r'{\color{red}\textbf{Removed (see the corrected version below).}} \begingroup\color{red}' + red[0] + r'\endgroup'
red[-1] = r'\begingroup\color{red}' + red[-1] + r'\endgroup'
NEW = r"""
{\color{blue}
\paragraph{Corrected: the $z$ integration, with $d_\perp$ kept general.} Every transverse trace is $\delta^{mm}=\delta^{ii}=d_\perp$; $d_\perp$ is not set to 2. The derivation, step by step, is in \texttt{rowIV\_d/RowIV\_dperp.pdf}. The $z$ integral is ultraviolet finite and infrared (collinear) divergent. In $d=2-\epsilon$ ($\overline{\rm MS}$ factor $(\mu^2e^{\gamma_E}/4\pi)^{\epsilon/2}$, the scheme of Group I, Row II), with the $d$-dimensional kernels $\mathcal K_d^m(r)=\pi^{1-\frac d2}\Gamma(\frac d2)\,r^m/|r|^d$, exactly,
\begin{align}
\int_{z}\mathcal K_d^{m}(y-z)\,\mathcal K_d^{m'}(y'-z)
&=\frac{\pi^{2-\frac d2}}{2}\big[(y'-y)^2\big]^{1-\frac d2}\Big(\frac{\mu^2e^{\gamma_E}}{4\pi}\Big)^{\frac\epsilon2}\Big[\Gamma\Big(\frac d2-1\Big)\delta^{mm'}-2\Gamma\Big(\frac d2\Big)\frac{(y'-y)^{m}(y'-y)^{m'}}{(y'-y)^{2}}\Big]\notag\\
&=-\frac{\pi}{2}\,E\,\delta^{mm'}-\pi\,\frac{(y'-y)^{m}(y'-y)^{m'}}{(y'-y)^{2}}+O(\epsilon),
\end{align}
with $E\equiv\frac{2}{\epsilon}+2\gamma_{E}+\log\frac{(y'-y)^{2}\mu^{2}_{\overline{MS}}}{4}=-\big[\frac{1}{\bar\epsilon}-\log\frac{(y'-y)^{2}\mu^{2}e^{2\gamma_E}}{4}\big]\equiv-\bar E$ ($\epsilon<0$: infrared). Its trace is $-\frac{\pi}{2}d_\perp E-\pi$. The $w,w'$ integrations give $(\frac{p^++k^+}{k^+})^2$ in any $d_\perp$ (the Jacobian $(\frac{P}{k^+})^d$ times the scaling $(\frac{P}{k^+})^{1-d}$ of $\mathcal K_d$, per side). The integral is symmetric in $m\leftrightarrow m'$, so the term $\frac{2}{k^+}(\delta_{mk'}\delta_{m'k}-\delta_{m'k'}\delta_{km})$ drops out for any $d_\perp$. With $\xi=\frac{k^{+}}{p^{+}+k^{+}}$, $\hat\Phi(\xi)=\xi(1-\xi)+\frac{\xi}{1-\xi}+\frac{1-\xi}{\xi}$ and exact limits ($\Lambda<p^{+}<\vee-k^{+}$ is $\frac{k^{+}}{\vee}<\xi<\frac{k^{+}}{k^{+}+\Lambda}$):
\begin{align}
\frac{d^{3}N_{\rm NLO}}{d^{3}k} ={}& \frac{1}{(2\pi)^{3}}\,\frac{g^{4}N_{c}}{16\pi^{4}}\,\frac{1}{k^{+}}\int_{k^{+}/\vee}^{k^{+}/(k^{+}+\Lambda)}\frac{d\xi}{\xi^{2}}\int_{x,x',y,y'}e^{-ik\cdot(y'-y)/\xi}\Bigg\{\frac{(x-y)\cdot(x'-y')}{(x-y)^{2}(x'-y')^{2}}\Big[d_\perp\,\hat\Phi(\xi)\,\bar E-2\Big(\frac{\xi}{1-\xi}+\frac{1-\xi}{\xi}\Big)\Big]\notag\\
&-2\,d_\perp\,\xi(1-\xi)\,\frac{(x-y)\cdot(y'-y)\,(x'-y')\cdot(y'-y)}{(x-y)^{2}(x'-y')^{2}(y'-y)^{2}}\Bigg\}\Big[U^{ce'}(y')-U^{ce'}(x')\Big]\rho^{e'}(x')\,\Big[U^{ce}(y)-U^{ce}(x)\Big]\rho^{e}(x).
\end{align}
The term $-2(\frac{\xi}{1-\xi}+\frac{1-\xi}{\xi})$ comes from the traceless part of the $z$ integral contracted with $\delta_{mm'}$ and does not depend on $d_\perp$. For $d_\perp=2-\epsilon$ the braces become $2\{X\cdot X'[\hat\Phi\bar E+\xi(1-\xi)]-2\xi(1-\xi)(X\cdot\hat\Delta)(X'\cdot\hat\Delta)\}$ (your ``$+1$'' in the second line of (1.15) is $d_\perp\bar E$ at $d_\perp=2-\epsilon$; note $\int_{1}^{k^+/\vee}=-\int_{k^+/\vee}^{1}$); for $d_\perp=2$ they become $2\{X\cdot X'[\hat\Phi\bar E-\frac{\xi}{1-\xi}-\frac{1-\xi}{\xi}]-2\xi(1-\xi)(X\cdot\hat\Delta)(X'\cdot\hat\Delta)\}$. Use the same $d_\perp$ in all rows.
\paragraph{Corrected: the endpoint $\xi\to1$ (soft unobserved gluon).} With $\frac{\xi}{1-\xi}=\frac{1}{1-\xi}-1$ and, exactly, $\int_{\xi_0}^{\xi_1}d\xi\,\frac{h(\xi)}{1-\xi}=\int_{\xi_0}^{\xi_1}d\xi\,\frac{h(\xi)-h(1)}{1-\xi}+h(1)\log\frac{1-\xi_0}{1-\xi_1}$, where
\begin{align}
\log\frac{1-\xi_{0}}{1-\xi_{1}}=\log\frac{(\vee-k^{+})(k^{+}+\Lambda)}{\vee\,\Lambda}=\log\Big[\frac{\vee}{\Lambda}\Big]-\log\Big[\frac{\vee}{k^{+}}\Big]+\log\Big[1-\frac{k^{+}}{\vee}\Big]+\log\Big[1+\frac{\Lambda}{k^{+}}\Big],
\end{align}
the soft endpoint is (its coefficient is $d_\perp\bar E-2$)
\begin{align}
\frac{d^{3}N_{\rm NLO}}{d^{3}k}\bigg|_{\rm endpoint} ={}& \frac{1}{(2\pi)^{3}}\,\frac{g^{4}N_{c}}{16\pi^{4}}\,\frac{1}{k^{+}}\left[\log\Big[\frac{\vee}{\Lambda}\Big]-\log\Big[\frac{\vee}{k^{+}}\Big]+\log\Big[1-\frac{k^{+}}{\vee}\Big]+\log\Big[1+\frac{\Lambda}{k^{+}}\Big]\right]\notag\\
&\times\int_{x,x',y,y'}e^{-ik\cdot(y'-y)}\,\frac{(x-y)\cdot(x'-y')}{(x-y)^{2}(x'-y')^{2}}\big(d_\perp\bar E-2\big)\Big[U^{ce'}(y')-U^{ce'}(x')\Big]\rho^{e'}(x')\,\Big[U^{ce}(y)-U^{ce}(x)\Big]\rho^{e}(x),
\end{align}
and the rest is the $\xi$ integral above with $\frac{1}{\xi^{2}(1-\xi)}e^{-ik\cdot(y'-y)/\xi}\to\frac{1}{1-\xi}\big[\frac{1}{\xi^{2}}e^{-ik\cdot(y'-y)/\xi}-e^{-ik\cdot(y'-y)}\big]$ in the $\frac{1}{1-\xi}$ terms, over the same limits. Your $\log\frac{\vee}{\Lambda}+\log\frac{k^{+}}{\vee}$ misses $\log(1-\frac{k^{+}}{\vee})+\log(1+\frac{\Lambda}{k^{+}})$, and your $\int_{1}^{k^{+}/\vee}$ uses $\xi_1=1$ instead of $\frac{k^{+}}{k^{+}+\Lambda}$.
\paragraph{Corrected: the fragmentation form.} With the LO cross section (emitter at $y$)
\begin{align}
\frac{d^{3}N_{\rm LO}}{d^{3}k}(k^{+},\mathbf k)=\frac{1}{(2\pi)^{3}}\frac{g^{2}}{2\pi^{2}k^{+}}\int_{x,x',y,y'}e^{-ik\cdot(y'-y)}\frac{(x-y)\cdot(x'-y')}{(x-y)^{2}(x'-y')^{2}}\Big[U^{ce'}(y')-U^{ce'}(x')\Big]\rho^{e'}(x')\Big[U^{ce}(y)-U^{ce}(x)\Big]\rho^{e}(x),
\end{align}
the pole part of this row is
\begin{align}
\frac{d^{3}N_{\rm NLO}}{d^{3}k}\bigg|_{\rm pole}=\frac{\alpha_{s}}{2\pi}\,\frac{1}{\bar\epsilon}\int_{k^{+}/\vee}^{k^{+}/(k^{+}+\Lambda)}\frac{d\xi}{\xi^{3}}\,\frac{d_\perp}{2}\,P^{(r)}_{gg}(\xi)\,\frac{d^{3}N_{\rm LO}}{d^{3}k}\Big(\frac{k^{+}}{\xi},\frac{\mathbf k}{\xi}\Big),\qquad P^{(r)}_{gg}(\xi)=2N_{c}\,\hat\Phi(\xi).
\end{align}
The measure is $d\xi/\xi^{3}$, not $d\xi/\xi^{2}$: $\frac{d^{3}N_{\rm LO}}{d^{3}k}(k^{+}/\xi,\mathbf k/\xi)$ carries $\xi/k^{+}$. This is the real ($\xi<1$) part of the collinear counterterm of the gluon fragmentation function, with the polarization trace $\frac{d_\perp}{2}$. The $\delta(1-\xi)$ part of $P_{gg}$ ($\frac{\beta_0}{2}$ and the soft endpoint) is not produced by this row; it comes with the virtual terms of Group I. So $D_{g/g}\otimes{\rm LO}$ with $D_{g/g}=\delta(1-\xi)+O(g^{2})$ cannot be written for this row alone: its $\delta(1-\xi)$ term is the LO cross section itself. The constant in the logarithm is $e^{2\gamma_{E}}/4$, not $4\pi^{2}$.
}
"""
L[start:end + 1] = red + NEW.split('\n')
HEAD = r"""{\color{blue}\noindent\textbf{Corrections (blue; checked against Calculation\_main.tex, see Groups\_II\_III\_check.pdf):}
\begin{itemize}
\item Row I (parts 1 and 2) and Row IV up to the $z$ integration: correct. Only the integration measures are completed.
\item Row IV, $z$ integration and fragmentation form: replaced by the exact result (no limits on $\vee,k^+,\Lambda$).
\item Not changed here: index placements of the adjoint Wilson lines. They are correct with your convention $U^{bc}=U^{cb}$; for a general adjoint $U$ see Sec.~5 of the check note.
\end{itemize}}"""
i = L.index(r'\maketitle')
L.insert(i + 1, HEAD)
T = '\n'.join(L)
T = T.replace(r'\title{Group - II crossection derivation}', r'\title{Group - II crossection derivation {\color{blue}(corrected)}}')
open('Group_II_corrected.tex', 'w').write(T)
print('ok')
