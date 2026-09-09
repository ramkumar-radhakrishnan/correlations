# `trijet/` — extracting the rapidity divergence from LO trijet production at small x

Working notes on how to slice the divergent part out of the leading-order **trijet**
(`γ* + A → q q̄ g + X`) cross-section in the CGC, so that what remains is a finite resolved
trijet plus a finite virtual dijet, with the single logarithmic rapidity divergence — and only
that — absorbed into B-JIMWLK / BK evolution of the target.

Two documents:

1. **[`trijet_rapidity_divergence.pdf`](trijet_rapidity_divergence.pdf)** (18 pp) — the method:
   divergence inventory, the slicing construction, lifetime ordering, algorithm, pitfalls.
2. **[`trijet_step_by_step.pdf`](trijet_step_by_step.pdf)** (13 pp) — the calculation:
   every integral done explicitly, every cancellation shown, ending in two tables that list
   **all divergent terms** and **all finite terms**. Verification scripts in
   [`checks/`](checks/).

Contents:

| § | Topic |
|---|---|
| 0 | The master decomposition in one page |
| 1 | Frame, kinematics, the three distinct rapidity scales, colour correlators |
| 2 | Complete divergence inventory of the trijet: rapidity (D1), soft (D2), collinear (D3), UV (D4) |
| 3 | Two-cutoff phase-space slicing (Harris–Owens) translated to the CGC, with a full dictionary |
| 4 | The slice region by region: **S**, **HC**, **H̄C̄**, and the two cancellations |
| 5 | Lifetime ordering: why a naive `z_g < z_f` cut over-subtracts, the Lund-plane wedge, the constrained evolution equations |
| 6 | Master formula, step-by-step algorithm, and the checks to run |
| 7 | Reference formulas for every finite piece |
| 8 | Pitfalls |
| A | Full evaluation of the over-subtraction integral → `−(α_s N_c/2π) ln²(Q_f² r_bb'²)` |
| B | Glossary |

Sources: Caucal–Salazar–Schenke–Venugopalan, arXiv:2208.13872v2 (CGC conventions, NLO
decomposition, lifetime ordering); Harris–Owens, Phys. Rev. D **65** (2002) 094032 (slicing).

## Building the PDF

`pdf/document.html` is the typeset source; KaTeX renders the math and headless Chromium
prints it to PDF (no LaTeX toolchain required). Figures are inline SVG.

```bash
npm install katex          # provides node_modules/katex/dist/* at the repo root
pip install playwright
python3 pdf/build_pdf.py   # -> trijet_rapidity_divergence.pdf
```

`build_pdf.py` locates the pre-installed Chromium and fails loudly on any KaTeX parse error
or unrendered `$...$` fragment.


## The step-by-step calculation

`trijet_step_by_step.pdf` starts from the LO trijet cross-section and works through:

| Step | What is done |
|---|---|
| 1 | The trijet written out; gluon phase space; the three interference channels |
| 2 | Channel A (`R₂×R₂`, `C_F Ξ_LO`): eikonal split, **Integrals 1–3** → CSSV (B.7), (B.10), the `1/ε` coefficient |
| 3 | Channel B (`R₂×R₂'`, `Ξ_NLO,3`): **Integrals 4–5** → CSSV (B.20), (B.22) |
| 4 | Channel C (gluon crosses the shockwave): rapidity divergence only |
| 5 | Virtual inputs (quoted): (B.23), (B.30) |
| 6 | **`1/ε` cancellation** — KLN, coefficient = `∫dξ P_gq(ξ)`; `μ` dependence cancels too |
| 7 | **`ln²z₀` cancellation** — in-cone vs out-of-cone, and real vs virtual |
| 8 | Full assembly → CSSV (B.25); residual `z₀`-dependence is a *single* log |
| 9 | Rapidity split at `z_f` → divergent (B.26),(B.32) + finite (B.27),(B.33) |
| 10 | **Integral 6**: the kernel identities → `H_LL ⊗ dσ_LO`; BK at large `N_c` |
| 11 | **Integral 7**: the lifetime triangle → `−(α_s N_c/2π) ln²(Q_f² r_bb'²)` |
| 12–13 | The two answer tables |

Result: of the trijet's nine divergent structures, six cancel and three combine into
`ln(z_f/z_0) H_LL ⊗ dσ_LO` — the only thing that enters B-JIMWLK/BK.

```bash
python3 pdf/build_calc_pdf.py   # -> trijet_step_by_step.pdf
```
