# `nlo/` — longitudinal (`p^+`) integration of the NLO gluon-production term

Mathematica code for the `p^+` integral in

$$
\frac{d^{3}N_{\rm NLO}}{d^{3}k}=-\frac{1}{(2\pi)^{3}}\frac{ig^{4}f^{c'd'a}}{4\pi^{4}}\frac{1}{k^{+}}
\int_{\Lambda}^{\vee-k^{+}}\frac{dp^{+}}{2\pi}\int_{w,y',z,y,x,x'}
e^{-ik\cdot(y'-w)}e^{-ik\cdot\frac{p^{+}}{k^{+}}(y'-z)}\,
\frac{(y'-z)^{m}(x'-y')^{k'}(y-w)^{i}(x-z)^{j}}{(y'-z)^{2}(x'-y')^{2}(y-w)^{2}(x-z)^{2}}
\Big[\tfrac{\delta_{k'm}\delta_{ij}}{p^{+}+k^{+}}-\tfrac{\delta_{jm}\delta_{ik'}}{p^{+}}-\tfrac{\delta_{im}\delta_{jk'}}{k^{+}}\Big]\times(\text{color}),
$$

in the window $\vee \gg k^{+}\gg\Lambda>0$.

* `NLO_pplus_integral.nb` — **the notebook**. Self-contained (it does not load the
  package); open it and *Evaluation ▸ Evaluate Notebook*. Sections: notation, the
  integral by brute force, closed form, the $\vee\gg k^+\gg\Lambda$ limits,
  transverse structure, color structure, assembly, numerical audit.
* `NLOPplusIntegral.wl` — the same content as a loadable package.
* `nlo_pplus_demo.wls` — script driver for the package
  (`wolframscript -file nlo/nlo_pplus_demo.wls`).
* `pplus_scalar_integral.nb` / `.wls` — the scalar form of the same integral,
  `Integrate[Exp[-I a p/k] (A/(p+k) - B/p - 1/k), {p, Lam, V-k}]`, with its
  closed form, its $\vee\gg k\gg\Lambda$ limit and an exact-vs-asymptotic plot.
* `make_notebook.py`, `make_scalar_notebook.py` — regenerate the two `.nb` files
  from a cell list; edit and re-run rather than hand-editing notebook source.

## What actually has to be integrated

The only `p^+` dependence is the phase and the three poles, so the `p^+` integral
factorizes out of the transverse integrals *and* out of the color operators —
those pass through untouched. Writing

$$\omega \equiv k\cdot(y'-z),\qquad b\equiv\omega/k^{+},$$

the whole longitudinal job is the three integrals `PPlusKernel` returns
(each already carrying the $1/2\pi$ of the measure):

| pole | $2\pi\int_{\Lambda}^{\vee-k^{+}}\!dp^{+}\,e^{-ibp^{+}}\times$ pole |
|---|---|
| $1/(p^{+}+k^{+})$ | $e^{ibk^{+}}\big[\mathrm{Ei}(-ib\vee)-\mathrm{Ei}(-ib(k^{+}+\Lambda))\big]$ |
| $1/p^{+}$ | $\mathrm{Ei}\big(-ib(\vee-k^{+})\big)-\mathrm{Ei}(-ib\Lambda)$ |
| $1/k^{+}$ | $\big[e^{-ib\Lambda}-e^{-ib(\vee-k^{+})}\big]\big/(i b k^{+})$ |

`ExpIntegralEi` is evaluated just off the imaginary axis, so no branch cut is
crossed for either sign of $\omega$; the $-i0$ that puts it there is the usual
$p^{+}\to p^{+}-i0$ prescription.

The delta structures collapse the transverse numerator to three scalar pairings
(this is what `TransverseData` returns as `N1`, `N2`, `N3`):

| structure | contraction |
|---|---|
| $\delta_{k'm}\delta_{ij}$ | $[(y'-z)\cdot(x'-y')]\,[(y-w)\cdot(x-z)]$ |
| $\delta_{jm}\delta_{ik'}$ | $[(y'-z)\cdot(x-z)]\,[(y-w)\cdot(x'-y')]$ |
| $\delta_{im}\delta_{jk'}$ | $[(y'-z)\cdot(y-w)]\,[(x-z)\cdot(x'-y')]$ |

## The strongly ordered limit

Using $\mathrm{Ei}(-ix)=\gamma_{E}+\ln x-i\pi/2+O(x)$ and
$\mathrm{Ei}(-iX)\to-i\pi\,\mathrm{sgn}(X)$, `PPlusKernelAsymptotic` gives, for
$|\omega|\Lambda/k^{+}\ll1\ll|\omega|\vee/k^{+}$,

$$2\pi I_{1}=-e^{i\omega}\big[\mathrm{Ei}(-i\omega)+i\pi\,\mathrm{sgn}\,\omega\big],\qquad
2\pi I_{2}=\ln\frac{k^{+}}{|\omega|\Lambda}-\gamma_{E}-\frac{i\pi}{2}\mathrm{sgn}\,\omega,\qquad
2\pi I_{3}=\frac{1-e^{-i\omega(\vee-k^{+})/k^{+}}}{i\omega}.$$

Three things worth noting, all of them checked numerically in the demo:

1. **$\vee$ drops out of $I_1$ and $I_2$.** The phase, not the upper cutoff,
   provides the regulator once $|\omega|\vee\gg k^{+}$. The one surviving cutoff
   dependence is the single rapidity log $\ln(1/\Lambda)$ in $I_2$ — i.e. it rides
   entirely on the $-\delta_{jm}\delta_{ik'}$ structure (`RapidityLog`), which is
   the piece a JIMWLK-type evolution has to absorb.
2. **$I_3$ has no $\vee\to\infty$ limit pointwise** — its boundary term keeps
   oscillating with unit modulus. It averages to zero only against a smooth
   transverse profile, so keep it exact rather than taking a limit under the
   $\int_{y',z}$.
3. **Do not switch the phase off in $I_3$.** At $\omega\to0$
   (`PPlusKernelCollinear`) the $\vee$ log returns in $I_1,I_2$ and
   $I_3\to(\vee-k^{+}-\Lambda)/k^{+}$ diverges *linearly* — that term is finite
   only because of the phase.

## Verification

`CheckPPlusKernel[]` compares the quoted closed form against `Integrate` and
against direct `NIntegrate` over `p^+` for several $(\omega,k^{+},\Lambda,\vee)$,
and the asymptotics against the exact result. Agreement is at machine precision
for the exact forms ($\sim10^{-16}$) and at the expected $O(\Lambda,k^{+}/\omega\vee)$
level ($\sim10^{-9}$) for the asymptotics; the reassembled transverse factor
matches a direct one-dimensional `NIntegrate` of the whole bracket to $5\times10^{-17}$.

## One note on the measure

As written, the measure is $\int_{w,y',z,y,x}$ but the integrand contains $x'$
(in $(x'-y')$ and in $\rho^{e'}(x')$). The code integrates over
$\{w,y',z,y,x,x'\}$; if $x'$ is meant to be tied to one of the other points,
substitute for it before calling `NLOResult`.
