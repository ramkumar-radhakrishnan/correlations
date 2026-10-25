# Why BFKL comes with log(∨/k⁺) and JIMWLK with log(∨/Λ)

[`Rapidity_logs.pdf`](Rapidity_logs.pdf) (source `Rapidity_logs.tex`) is a passage ready to send to a collaborator, with a figure and references.

The figure is also available on its own:
- `Rapidity_figure.pdf`, `.png`, `.svg`;
- the TikZ source is `figure_body.tex`, wrapped by `Rapidity_figure.tex`.

The argument:
- The measured gluon k⁺ splits the soft window: log(∨/Λ) = log(∨/k⁺) + log(k⁺/Λ).
- Plus momentum is conserved, so the BFKL cascade ρ → p⁺ → k⁺ needs p⁺ > k⁺. BFKL therefore only gets log(∨/k⁺).
- The rows give log(∨/Λ) 𝒥 + log(∨/k⁺)(ℬ − 𝒥) = log(k⁺/Λ) 𝒥 + log(∨/k⁺) ℬ.

```
pdflatex Rapidity_figure.tex
pdflatex Rapidity_logs.tex   # twice
```

## Where the BFKL terms come from, and forward (AABKL) vs mid rapidity

[`BFKL_origin_and_AABKL.pdf`](BFKL_origin_and_AABKL.pdf) (source `.tex`; figure `Origin_figure.tex`, `origin_figure_body.tex`).

- The input terms ∝ log(∨/Λ) add up to JIMWLK⊗LO only. In every input term the measured gluon is emitted by a valence charge; in 23 of the 35 BFKL⊗LO terms it is emitted by the soft gluon.
- Each row is a p⁺ integral. The energy denominator splits 1/p⁺ into a p⁺ ≪ k⁺ part (log(k⁺/Λ)) and a p⁺ ≫ k⁺ part (log(∨/k⁺)). The BFKL terms live only in the second part.
- Cutoff independence fixes both equations, as in AABKL's derivation of BK:
  - −∂/∂lnΛ gives JIMWLK, which acts on U;
  - ∂/∂ln∨ at fixed k⁺ gives BFKL, which acts on ⟨ρρ⟩.
- The coefficient of log(∨/Λ) alone depends on the bookkeeping; the coefficients of log(k⁺/Λ) and log(∨/k⁺) do not.
- AABKL evolve the projectile over ln(1/x_p) = log(∨/k⁺), which is O(1) at forward rapidity. At mid rapidity it is ½ log(∨/Λ) and gives BFKL.
- Compared with AABKL:
  - same: target evolution, final state, colour structures;
  - different: LO and the whole projectile side;
  - scheme-dependent: the finite constants (Ioffe time vs sharp cutoffs).
