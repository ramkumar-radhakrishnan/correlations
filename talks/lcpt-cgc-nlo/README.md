# Theory seminar: LCPT and the soft gluon wave function in the CGC at NLO

Beamer slides on light-cone perturbation theory and
R. Radhakrishnan, *Soft Gluon Wave Function and Evolution Operator in the CGC at
Next-to-Leading Order*, arXiv:2607.18373.

**Audience assumption:** QCD experts with no small-$x$ / CGC background.
**Length:** 30–40 min (26 main slides + 5 backup).

## Build

```
pdflatex talk.tex && pdflatex talk.tex
```
Needs `beamer`, `beamertheme-metropolis`, `tikz`.
Put your name on the title page: replace `\author{\textit{[speaker]}}`.

## Suggested timing (35 min target)

| Slides | Section | Time |
|---|---|---|
| 1–2 | Title, roadmap | 2 min |
| 3–9 | **Part I** — light-front dynamics, LC gauge, LCPT rules | 12 min |
| 10–14 | **Part II** — saturation, CGC, why the wave function | 7 min |
| 15–26 | **Part III** — constructing $\Omega$, results, diagonalisation | 13 min |
| 30–31 | Summary / take-home | 2 min |

If you are short on time, the compressible slides are 4 (Dirac's forms),
13 (eikonal Hamiltonian details), 20 (explicit unitarity constraints) and
25 (structure of the $\mathcal{O}(g^2)$ coefficients) — state the result and move on.

## Backup slides (appendix, unnumbered)

- LCPT matrix elements used in the calculation
- The $\Theta$-function kinematic window for $a^\dagger a$ transitions
- Rayleigh–Schrödinger PT with non-commuting matrix elements
- Equal-time vs. light-front dictionary

The dictionary slide (last) is the most likely thing to be asked for by an
equal-time-minded audience.
