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
