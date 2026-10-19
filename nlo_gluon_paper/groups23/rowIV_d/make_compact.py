"""RowIV_compact.tex: the final Row IV result in compact form, checked against the term-by-term form (make_final.py)."""
import numpy as np
from make_final import terms, q


def compact_value(c, kp, Lam, vee):
    """(d Eb - 2) XX [L ph(1) + int dxi (ph/xi^2 - ph(1))/(1-xi) + int dxi (1-2xi)/xi^3 ph] + d int dxi (1-xi)/xi ph [Eb XX - 2 XD]"""
    x0, x1 = kp / vee, kp / (kp + Lam)
    ph = lambda x: np.exp(-1j * c['a'] / x)
    L = np.log(vee / Lam) - np.log(vee / kp) + np.log(1 - kp / vee) + np.log(1 + Lam / kp)
    A = (c['dp'] * c['Eb'] - 2) * c['sXX'] * (L * ph(1.0) + q(lambda x: (ph(x) / x**2 - ph(1.0)) / (1 - x), x0, x1)
                                             + q(lambda x: (1 - 2 * x) / x**3 * ph(x), x0, x1))
    B = c['dp'] * q(lambda x: (1 - x) / x * ph(x), x0, x1) * (c['Eb'] * c['sXX'] - 2 * c['sXD'])
    return A + B


def terms_value(c, kp, Lam, vee):
    x0, x1 = kp / vee, kp / (kp + Lam)
    ph = lambda x: np.exp(-1j * c['a'] / x)
    L = np.log(vee / Lam) - np.log(vee / kp) + np.log(1 - kp / vee) + np.log(1 + Lam / kp)
    tot = 0
    for s, _, kind, f in terms:
        if kind == 'end':
            tot += s * L * ph(1.0) * f(c)
        elif kind == 'plus':
            tot += s * f(c) * q(lambda x: (ph(x) / x**2 - ph(1.0)) / (1 - x), x0, x1)
        else:
            tot += s * q(lambda x, f=f: ph(x) / x**2 * f(c, x), x0, x1)
    return tot


rng = np.random.default_rng(8)
worst = 0
for t in range(6):
    kp, Lam = rng.uniform(0.5, 2), rng.uniform(0.01, 0.3)
    vee = kp + rng.uniform(1, 30)
    c = dict(dp=rng.uniform(1.5, 3.5), Eb=rng.normal() * 3, sXX=rng.normal(), sXD=rng.normal(), a=rng.uniform(0.1, 2))
    a, b = compact_value(c, kp, Lam, vee), terms_value(c, kp, Lam, vee)
    worst = max(worst, abs(a - b) / abs(b))

doc = r'''\documentclass[11pt]{article}
\usepackage[a4paper,margin=0.7in]{geometry}
\usepackage{amsmath,amssymb}
\begin{document}
\section*{Group II, Row IV: final result, compact form ($\delta^{mm}=d_\perp$, exact limits)}
\begin{align}
\frac{d^{3}N_{\rm IV}}{d^{3}k}={}&\frac{1}{(2\pi)^{3}}\,\frac{g^{4}N_{c}}{16\pi^{4}}\,\frac{1}{k^{+}}\int_{x,x',y,y'}\mathcal{O}\,\Bigg\{\big(d_{\perp}\bar{E}-2\big)\,\frac{(x-y)\cdot(x'-y')}{(x-y)^{2}(x'-y')^{2}}\notag\\
&\times\bigg[L\,e^{-ik\cdot(y'-y)}+\int_{\xi_0}^{\xi_1}\!d\xi\,\frac{\xi^{-2}e^{-ik\cdot(y'-y)/\xi}-e^{-ik\cdot(y'-y)}}{1-\xi}+\int_{\xi_0}^{\xi_1}\!d\xi\,\frac{1-2\xi}{\xi^{3}}\,e^{-ik\cdot(y'-y)/\xi}\bigg]\notag\\
&+d_{\perp}\int_{\xi_0}^{\xi_1}\!d\xi\,\frac{1-\xi}{\xi}\,e^{-ik\cdot(y'-y)/\xi}\bigg[\bar{E}\,\frac{(x-y)\cdot(x'-y')}{(x-y)^{2}(x'-y')^{2}}-2\,\frac{(x-y)\cdot(y'-y)\,(x'-y')\cdot(y'-y)}{(x-y)^{2}(x'-y')^{2}(y'-y)^{2}}\bigg]\Bigg\},
\end{align}
where
\begin{align*}
&\bar{E}=\frac{1}{\bar\epsilon}-\log\frac{(y'-y)^{2}\mu^{2}_{\overline{MS}}\,e^{2\gamma_{E}}}{4},\qquad
L=\log\frac{\vee}{\Lambda}-\log\frac{\vee}{k^{+}}+\log\Big(1-\frac{k^{+}}{\vee}\Big)+\log\Big(1+\frac{\Lambda}{k^{+}}\Big),\\
&\xi_0=\frac{k^{+}}{\vee},\qquad \xi_1=\frac{k^{+}}{k^{+}+\Lambda},\qquad
\mathcal{O}=\Big[U^{ce'}(y')-U^{ce'}(x')\Big]\rho^{e'}(x')\,\Big[U^{ce}(y)-U^{ce}(x)\Big]\rho^{e}(x).
\end{align*}
\end{document}
'''
open('RowIV_compact.tex', 'w').write(doc)
print('wrote RowIV_compact.tex; compact vs ten terms, max rel. diff %.1e' % worst)
