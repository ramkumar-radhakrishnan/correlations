# Theory seminar: *Chiral soliton lattice in inhomogeneous magnetic fields*

A 30–40 minute talk for a QCD-theory audience, built as a story in four acts:
**a wall → a lattice → a bent field → what we learned.**

```
pdflatex csl_seminar      # twice (metropolis needs the second pass for the TOC/progress bar)
```
Output: `csl_seminar.pdf` — 25 content slides + 6 backup slides. Requires
`beamer`, `beamertheme-metropolis`, `appendixnumberbeamer` (all in TeX Live).

Edit the speaker name on line 49 of `csl_seminar.tex`.

## Speaker notes

Every frame carries a `\note{}` with a target duration and the one sentence to
land. To get a notes version, add after `\begin{document}`:

```latex
\setbeameroption{show notes}                 % notes on separate pages
% \setbeameroption{show notes on second screen=right}   % for a dual-screen setup
```

## Timing plan (~39 min + questions)

| # | Act | Slides | Target |
|---|-----|--------|--------|
| 0 | Hook + roadmap | The setting; The story in four steps | 3 min |
| I | The anomaly builds a wall | Why only π⁰; anomalous term; one wall → a stack | 6 min |
| II | The lattice, solved exactly | pendulum; gap equation; genuine phase; phase diagram | 8 min |
| III | **What happens when the field bends** | three reasons; finite box; chiral limit; two theorems; separation of variables; μ(x); massive pions; numerics; 2D results; domain wall; energy costs; the Matsubara gem | 18 min |
| IV | Where this leaves us | summary; outlook; thank you | 4 min |

**To cut to 30 min**, drop (in this order):
1. *A worked family: separation of variables* — mention that closed forms exist and move on.
2. *Bonus: a non-uniform chemical potential* — one sentence on the previous slide.
3. *How much does bending cost?* — the qualitative statement is repeated in the summary.

Do **not** cut *The chiral limit* → *Two theorems*: that pair is the core of the
new work.

## Figures (`figs/`)

| file | source |
|---|---|
| `csl_profile.pdf`, `csl_kofB.pdf`, `csl_phase.pdf` | recomputed from the Brauner–Yamamoto formulae (`make_figures.py`) |
| `csl_finitevol.pdf` | solutions of the EoM satisfying the natural BC, computed here |
| `csl_dwenergy.pdf` | Eq. (A2)/(A5) of arXiv:2608.27319, computed here |
| `csl_geometries.pdf` | schematic of the three field geometries |
| `paper_fig1..6.png` | Figs. 1–6 of arXiv:2608.27319 |

Regenerate the computed ones with

```
python3 make_figures.py      # needs numpy, scipy, matplotlib
```

Sanity checks built into these curves (all pass):
`E(k)/k → 1` at `k → 1`; the finite-volume BC solution at `L̄ = 5, H̄ = 2` gives
`k = 0.804` and a peak gradient `2.487`, matching Fig. 2 of the paper;
`⟨E_DW⟩/E₀ = 0.4998` at `α = 1` (exactly 1/2) and `→ 1` as `α → 0`;
`B_BEC(chiral) = 0.14 GeV²` at `μ = 900 MeV`.

## Sources

- D. T. Son, M. A. Stephanov, PRD **77** (2008) 014021 [arXiv:0710.1084]
- T. Brauner, N. Yamamoto, JHEP **04** (2017) 132 [arXiv:1609.05213]
- T. Brauner, R. Radhakrishnan, arXiv:2608.27319
