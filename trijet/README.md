# `trijet/` — extracting the rapidity divergence from LO trijet production at small x

Working notes on how to slice the divergent part out of the leading-order **trijet**
(`γ* + A → q q̄ g + X`) cross-section in the CGC, so that what remains is a finite resolved
trijet plus a finite virtual dijet, with the single logarithmic rapidity divergence — and only
that — absorbed into B-JIMWLK / BK evolution of the target.

**Read [`trijet_rapidity_divergence.pdf`](trijet_rapidity_divergence.pdf)** (18 pp).

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
