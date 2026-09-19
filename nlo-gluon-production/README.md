# $\beta_0$ logs in NLO single-inclusive gluon production: DGLAP or running coupling?

A transposition of Kovner–Lublinsky–Skokov–Zhao, *Not all that is $\beta_0$ is
$\beta$-function* ([arXiv:2308.15545](https://arxiv.org/abs/2308.15545)), from the
NLO JIMWLK Hamiltonian to single-inclusive gluon production at $O(g^4)$.

**Read [`NLO_gluon_production_DGLAP.pdf`](NLO_gluon_production_DGLAP.pdf).**

## What it establishes

The two $\rho\rho$ terms of the $O(g^4)$ production cross section map one-to-one onto
the two $\beta_0$-carrying NLO JIMWLK kernels, and the mapping is exact rather than
analogical:

- $\mathfrak{C}_1$ *is* the light-cone $g\to gg$ splitting vertex — its square is
  $P_{gg}(\zeta)/2N_c$, identically KLSZ's combination, and its $\delta^{(2)}$ places
  the parent at the momentum-weighted centroid of the daughters.
- $\mathfrak{B}_2$, in the merging limit, factorizes *exactly* into
  $\mathfrak{A}^{(1)}\otimes\mathfrak{C}_1$ — prefactors, $4\pi^2$ and all.

So the coefficient of the transverse logarithm is forced to be $\beta_0^g = 11N_c/3$,
half of it running coupling and half of it projectile DGLAP.

Two results are specific to production and not in KLSZ:

- the measurement phase carries **no** dependence on the pair separation $Z$, so $k_T$
  does not regulate the collinear integral — it enters only through the
  emitter–emission distance;
- consequently the residual DGLAP-like log is $\ln(Q_s^2/k_T^2)$, large in the
  **saturation region $k_T \lesssim Q_s$** and switching off above it. The resummed
  correction is a power $(Q_s^2/k_T^2)^{\alpha_s\beta_0^g/4\pi}$ on the Wilson line —
  a factor of $2$–$4$ near $k_T\sim Q_s/3$.

Section 10 of the note lists what is *not* settled, in the order worth attacking.

## Checks

```bash
python3 verify.py     # 14 checks: sympy + numpy + scipy
```

Covers the splitting-function algebra, the $-11/6$ double-plus integral, the merging
factorization (residue and prefactors), the $SU(3)$ colour collapse, and the
$\overline{\rm MS}$ transverse integral.

## Rebuilding the PDF

`pdf/document.html` is the typeset source; KaTeX renders the math and headless
Chromium prints it (no LaTeX toolchain). The `<!--FIGURE-->` placeholder is filled at
build time from `pdf/figure_factor.svg`.

```bash
cd pdf
npm install katex          # provides node_modules/katex/dist/*
pip install playwright
python3 make_figure.py     # regenerates figure_factor.svg
python3 build_pdf.py       # -> ../NLO_gluon_production_DGLAP.pdf
```

`build_pdf.py` fails loudly on any KaTeX parse error or unrendered `$...$` fragment.
