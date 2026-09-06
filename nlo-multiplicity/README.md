# Gluon multiplicity moments and fluctuations from the NLO CGC evolution operator

Research analysis of the program announced at the end of the Introduction of
[arXiv:2607.18373](https://arxiv.org/abs/2607.18373), *Soft Gluon Wave Function and Evolution
Operator in the CGC at Next-to-Leading Order*: **moments of the gluon multiplicity operator and
particle-number fluctuations at order `g⁴`**.

**Start with [`docs/00-scope-and-caveats.md`](docs/00-scope-and-caveats.md)** — Parts I–IX were
first drafted without access to the paper and then reconciled against v2; that file records
exactly which conclusions survived and which coefficient identifications changed.

The whole analysis is also typeset as a single 40-page document,
[`NLO_CGC_multiplicity_analysis.pdf`](NLO_CGC_multiplicity_analysis.pdf):

```
pip install reportlab
python3 pdf/build_pdf.py
```

## The three results that carry the project

1. **`Ω|0⟩` is a truncated displaced squeezed state.** `Ω` itself is *not* Gaussian — Eq. (4.1)
   carries an `O(g)` trilinear `C` and `O(g²)` quartic/sextic structures. But it is written fully
   normal-ordered, so those all annihilate the soft vacuum:
   ```
   Ω|0⟩ = N|0⟩ + ∫A₂ a†|0⟩ + ∫∫B₂ a†a†|0⟩          exact through O(g²)
   ```
   Every moment then follows from a Fredholm determinant:
   ```
   G(z) = det[1 − (z−1)K]^(−1/2) · exp{ ½(z−1) w̄ [1−(z−1)K]^(−1) w }
   C_p  = (p−1)!/2 · Tr(Kᵖ)  +  p!/2 · w̄ K^(p−1) w        ⇒   C_p = O(g^{2p})
   ```
   verified against four quantum-optics limits (coherent, thermal, squeezed vacuum, displaced
   thermal).

2. **The leading non-Poissonian effect is displacement–squeezing interference, not `|B₂|²`.**
   Because `m = ⟨aa⟩_c = B₂` is `O(g²)` while `n = ⟨a†a⟩_c = B₂†B₂` is `O(g⁴)`:
   ```
   C₂ = 2 Re ∫ A₂*(k) B₂(k,q) A₂*(q)  +  ∫ |B₂(k,q)|²   +  O(g⁶)      [complete at O(g⁴)]
   ```
   The first term is an interference and is **not sign-definite** — a sub-Poissonian region would
   be genuinely non-classical gluon statistics.

3. **Angular structure separates new from known physics.** Bose enhancement peaks near-side
   (`k ≈ +q`, `O(g⁶)`); pair emission from `B₂` peaks **away-side** (`k ≈ −q`, `O(g⁴)`). At NLO
   the away-side term should dominate — the sharpest open prediction here.

## The two limits stated honestly

* `C₂` is complete at `O(g⁴)`. **`C₃` and `C₄` are not** — the three-gluon component of `Ω|0⟩`
  first appears at `O(g³)` and contributes at the same order. `κ₃` and `κ₄` are complete only
  through `O(g⁴)` (`⟨N⟩ + 3C₂`, `⟨N⟩ + 7C₂`).
* Classical colour-charge fluctuations and the new quantum fluctuation are **the same order in
  `g` and carry the same `1/N_dof` suppression**. They cannot be separated by power counting —
  only by angular and `k_⊥` structure.

## Assessment

**One PRD paper**, ~70% confidence, conditional on a two-week go/no-go test. Not JHEP as scoped.
More than a note. Main risks: scooping ([arXiv:2605.03038](https://arxiv.org/abs/2605.03038),
[arXiv:2412.05085](https://arxiv.org/abs/2412.05085)), the effect being numerically negligible,
and the "wave-function multiplicity is not an observable" referee report.

Note also that the multiplicity program needs **only Eqs. (5.1)–(5.3)** of the paper; `C`, `D`,
`E` come alive only once the target Wilson lines act (`a → U a`), which is Project E.

See [`docs/08-publication.md`](docs/08-publication.md) and [`docs/09-roadmap.md`](docs/09-roadmap.md).

## Contents

| File | Part |
|---|---|
| [`docs/00-scope-and-caveats.md`](docs/00-scope-and-caveats.md) | Notation, verified coefficient dictionary, what changed after reading the paper |
| [`docs/01-formalism.md`](docs/01-formalism.md) | I — structure of `Ω`, Bogoliubov reduction, law of total cumulance |
| [`docs/02-moments.md`](docs/02-moments.md) | II — `⟨N⟩ … ⟨N⁴⟩`, generating function, order counting |
| [`docs/03-cumulants.md`](docs/03-cumulants.md) | III — `κ₂, κ₃, κ₄`, factorial cumulants, physical meaning |
| [`docs/04-correlations.md`](docs/04-correlations.md) | IV — `n`-gluon correlations, `c₂ = 1/k_NBD`, KNO |
| [`docs/05-models.md`](docs/05-models.md) | V — seven model calculations, integrals, complexity |
| [`docs/06-numerics.md`](docs/06-numerics.md) | VI — twelve proposed plots, prioritized |
| [`docs/07-literature.md`](docs/07-literature.md) | VII — survey: known vs open |
| [`docs/08-publication.md`](docs/08-publication.md) | VIII — ranked, critical publication assessment |
| [`docs/09-roadmap.md`](docs/09-roadmap.md) | IX — phased research plan with times and risks |
| [`docs/10-jhep-assessment.md`](docs/10-jhep-assessment.md) | Addendum — will the paper itself make it into JHEP? (~65–75%) |

A separate, equation-by-equation error review of the draft manuscript lives in
[`../review/draft-review-omega-g2.md`](../review/draft-review-omega-g2.md).
