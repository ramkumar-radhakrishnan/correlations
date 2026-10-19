"""RowIV_final.tex: the final Row IV result written as a sum of separate terms.
Each term is stored with its LaTeX and its numerical integrand; their sum is checked against the compact form
of RowIV_dperp.pdf, Eq. (1.5), integrated over the exact xi range (scalar test values for the transverse factors)."""
import numpy as np
from scipy import integrate

PH = r"e^{-ik\cdot(y'-y)}"
PHX = r"e^{-ik\cdot(y'-y)/\xi}"
XX = r"\frac{(x-y)\cdot(x'-y')}{(x-y)^{2}\,(x'-y')^{2}}"
XD = r"\frac{(x-y)\cdot(y'-y)\,(x'-y')\cdot(y'-y)}{(x-y)^{2}\,(x'-y')^{2}\,(y'-y)^{2}}"
EB = r"\bigg[\frac{1}{\bar\epsilon}-\log\frac{(y'-y)^{2}\mu^{2}_{\overline{MS}}\,e^{2\gamma_{E}}}{4}\bigg]"
COL = r"\Big[U^{ce'}(y')-U^{ce'}(x')\Big]\rho^{e'}(x')\,\Big[U^{ce}(y)-U^{ce}(x)\Big]\rho^{e}(x)"
PRE = r"\frac{1}{(2\pi)^{3}}\,\frac{g^{4}N_{c}}{16\pi^{4}}\,\frac{1}{k^{+}}"
LOGS = (r"\bigg[\log\Big(\frac{\vee}{\Lambda}\Big)-\log\Big(\frac{\vee}{k^{+}}\Big)+\log\Big(1-\frac{k^{+}}{\vee}\Big)"
        r"+\log\Big(1+\frac{\Lambda}{k^{+}}\Big)\bigg]")
XI = r"\int_{\frac{k^{+}}{\vee}}^{\frac{k^{+}}{k^{+}+\Lambda}}"
INT = r"\int_{x,x',y,y'}"
PLUS = r"\frac{1}{1-\xi}\bigg[\frac{1}{\xi^{2}}\,e^{-ik\cdot(y'-y)/\xi}-e^{-ik\cdot(y'-y)}\bigg]"

# (sign, latex body after the prefactor, kind, numeric integrand of xi or endpoint factor)
# numeric model: phase(xi) = exp(-i a / xi), XX -> sXX, XD -> sXD, Ebar -> Eb, d_perp -> dp
terms = [
    (+1, r"d_{\perp}\," + LOGS + "|" + INT + PH + XX + EB, 'end', lambda c: c['dp'] * c['Eb'] * c['sXX']),
    (-1, r"2\," + LOGS + "|" + INT + PH + XX, 'end', lambda c: 2 * c['sXX']),
    (+1, r"d_{\perp}" + XI + r"d\xi" + INT + PLUS + "|" + XX + EB, 'plus', lambda c: c['dp'] * c['Eb'] * c['sXX']),
    (+1, r"d_{\perp}" + XI + r"\frac{d\xi}{\xi^{2}}" + INT + PHX + r"\,\xi(1-\xi)" + "|" + XX + EB, 'reg',
     lambda c, x: c['dp'] * c['Eb'] * c['sXX'] * x * (1 - x)),
    (+1, r"d_{\perp}" + XI + r"\frac{d\xi}{\xi^{2}}" + INT + PHX + r"\,\frac{1-\xi}{\xi}" + "|" + XX + EB, 'reg',
     lambda c, x: c['dp'] * c['Eb'] * c['sXX'] * (1 - x) / x),
    (-1, r"d_{\perp}" + XI + r"\frac{d\xi}{\xi^{2}}" + INT + PHX + "|" + XX + EB, 'reg',
     lambda c, x: c['dp'] * c['Eb'] * c['sXX']),
    (-1, r"2" + XI + r"d\xi" + INT + PLUS + "|" + XX, 'plus', lambda c: 2 * c['sXX']),
    (+1, r"2" + XI + r"\frac{d\xi}{\xi^{2}}" + INT + PHX + "|" + XX, 'reg', lambda c, x: 2 * c['sXX']),
    (-1, r"2" + XI + r"\frac{d\xi}{\xi^{2}}" + INT + PHX + r"\,\frac{1-\xi}{\xi}" + "|" + XX, 'reg',
     lambda c, x: 2 * c['sXX'] * (1 - x) / x),
    (-1, r"2\,d_{\perp}" + XI + r"\frac{d\xi}{\xi^{2}}" + INT + PHX + r"\,\xi(1-\xi)" + "|" + XD, 'reg',
     lambda c, x: 2 * c['dp'] * c['sXD'] * x * (1 - x)),
]


def q(f, a, b):
    return (integrate.quad(lambda x: f(x).real, a, b, limit=400, epsrel=1e-13)[0]
            + 1j * integrate.quad(lambda x: f(x).imag, a, b, limit=400, epsrel=1e-13)[0])


def check():
    worst = 0.0
    rng = np.random.default_rng(3)
    for t in range(5):
        kp, Lam = rng.uniform(0.5, 2), rng.uniform(0.01, 0.3)
        vee = kp + rng.uniform(1, 30)
        c = dict(dp=rng.uniform(1.5, 3.5), Eb=rng.normal() * 3, sXX=rng.normal(), sXD=rng.normal(), a=rng.uniform(0.1, 2))
        x0, x1 = kp / vee, kp / (kp + Lam)
        ph = lambda x: np.exp(-1j * c['a'] / x)
        Phi = lambda x: x * (1 - x) + x / (1 - x) + (1 - x) / x
        compact = q(lambda x: ph(x) / x**2 * ((c['dp'] * Phi(x) * c['Eb'] - 2 * (x / (1 - x) + (1 - x) / x)) * c['sXX']
                                                - 2 * c['dp'] * x * (1 - x) * c['sXD']), x0, x1)
        L = np.log(vee / Lam) - np.log(vee / kp) + np.log(1 - kp / vee) + np.log(1 + Lam / kp)
        tot = 0
        for s, _, kind, f in terms:
            if kind == 'end':
                tot += s * L * ph(1.0) * f(c)
            elif kind == 'plus':
                tot += s * f(c) * q(lambda x: (ph(x) / x**2 - ph(1.0)) / (1 - x), x0, x1)
            else:
                tot += s * q(lambda x, f=f: ph(x) / x**2 * f(c, x), x0, x1)
        worst = max(worst, abs(tot - compact) / abs(compact))
    return worst


rel = check()
lines = []
for n, (s, body, kind, f) in enumerate(terms):
    sign = '' if (n == 0 and s > 0) else ('+' if s > 0 else '-')
    lead = r'\frac{d^{3}N_{\rm IV}}{d^{3}k}={}&' if n == 0 else '&'
    head, mid = body.split('|')
    lines.append(lead + sign + PRE + r'\,' + head + r'\notag\\')
    lines.append(r'&\qquad\times' + mid + r'\notag\\')
    lines.append(r'&\qquad\times' + COL + r'\notag\\[8pt]')
lines[-1] = lines[-1].replace(r'\notag\\[8pt]', '')
doc = r'''\documentclass[10pt]{article}
\usepackage[a4paper,margin=0.5in]{geometry}
\usepackage{amsmath,amssymb}
\allowdisplaybreaks
\begin{document}
\section*{Group II, Row IV: final result ($\delta^{mm}=d_\perp$, exact limits)}
\begin{align}
''' + '\n'.join(lines) + r'''
\end{align}
\end{document}
'''
open('RowIV_final.tex', 'w').write(doc)
print('wrote RowIV_final.tex; sum of terms vs compact form, max rel. diff %.1e' % rel)
