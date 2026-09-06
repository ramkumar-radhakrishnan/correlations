# Scope, notation, and honest caveats

This directory contains a research analysis of the program announced at the end of the
Introduction of

> R. Radhakrishnan, *Soft Gluon Wave Function and Evolution Operator in the CGC at
> Next-to-Leading Order*, arXiv:2607.18373 [hep-ph] (2026),

namely: **moments of the gluon multiplicity operator and particle-number fluctuations at
order `g^4`**.

## Status — reconciled against the paper (v2, 57 pp.)

Parts I–IX were first drafted without the paper (session egress blocked arxiv.org and all
mirrors). The paper has since been read. **The structural conclusions survived; the coefficient
dictionary did not.** This section records both, because the corrections matter.

### What was wrong

I guessed that the five labels `A…E` ran over the five Fock structures available to an `O(g²)`
normal-ordered operator (`a†`, `a`, `a†a`, `a†a†`, `aa`). **They do not.** The paper's Eq. (4.1)
organizes by *operator length*, not by order in `g`, and one of the coefficients sits at `O(g)`
with three operators — a possibility I did not anticipate:

| Paper (Eq. 4.1) | Fock structure | Order | Content |
|---|---|---|---|
| `N` (blackboard `N`) | `1` | `N − 1 = O(g²)` | normalization / virtual correction |
| `A₁`, `A₂` | `a`, `a†` | `O(g)` | soft emission/absorption off the valence charge, `∝ ρ` |
| `B₁`, `B₂`, `B₃` | `aa`, `a†a†`, `a†a` | `O(g²)` | pair emission, pair absorption, rescattering |
| `C₁`, `C₂` | `a†a†a`, `a†aa` | **`O(g)`** | **soft→soft splitting from the three-gluon vertex; carries no `ρ`** |
| `D₁,₂,₃` | four operators | `O(g²)` | |
| `E₁,₂,₃` | six operators | `O(g²)` | built from `C·C` |

The `O(g)` trilinear `C` is the substantive surprise. It means **`Ω` is not a Gaussian
(quadratic-exponent) unitary at any order** — my §1.5 claim, as a statement about `Ω`, is false.

### What survived, and why

The Gaussian machinery is nevertheless *correct for the observable*, for a reason that is now
verifiable rather than assumed. `Ω` in Eq. (4.1) is written **fully normal-ordered**, so every
`C`, `D`, `E` structure carries at least one annihilation operator on the right and annihilates
the soft vacuum. Acting on `|0⟩`, only three terms survive:

```
Ω|0⟩  =  N|0⟩  +  ∫ A₂ a†|0⟩  +  ∫∫ B₂ a†a†|0⟩          (exact through O(g²))
```

which is precisely a truncated **displaced squeezed state**. So the identifications used
throughout Parts II–IV hold, with this dictionary:

```
   this document          paper (Eq. 4.1, 5.1–5.3)      order
   ─────────────────────────────────────────────────────────────
   α   (displacement)  →  A₂  [Eq. (5.2)]               O(g),  ∝ ρ
   m   (pair/anomalous)→  B₂  [Eq. (5.3)]               O(g²), ∝ ρρ ⊕ ρ
   n   (occupation)    →  B₂†B₂                          O(g⁴)
   Z   (normalization) →  N   [Eq. (5.1)]               1 + O(g²)
   my "A"              →  A₂           my "D"  →  B₂
   my "B"              →  A₁           my "E"  →  B₁
   my "C" (a†a)        →  B₃
```

Eq. (5.2) is `A₂ᵇⱼ(p) = −√2 g pʲ ρᵇ(−p)/(√p⁺ p²)` — the Weizsäcker–Williams field, as assumed.
Eq. (5.3) for `B₂` contains exactly the two pieces predicted in §1.2: a `ρρ` two-source term and
an `f^{acb}ρ` term from the three-gluon vertex, plus an instantaneous `H_gg-inst` contribution.

**Consequences for the research program:**

* The headline formula is unchanged, now with verified labels:
  `C₂ = 2Re ∫ A₂* B₂ A₂* + ∫ |B₂|²`, complete at `O(g⁴)`.
* `C₃` incompleteness is **confirmed**: `Ω|0⟩` has no three-gluon component at `O(g²)`, so the
  three-gluon amplitude first appears at `O(g³)` and is not in this operator.
* **The wave-function multiplicity program needs only Eqs. (5.1)–(5.3)** — three of the paper's
  coefficients. `C`, `D`, `E` are irrelevant to `Ω|0⟩`. They become essential the moment the
  target Wilson lines act (`a → U a`), because then the annihilation operators no longer kill the
  state. That is Project E in Part VIII, and it is the natural home for the rest of the operator.

### What is still approximate here

Numerical prefactors, phases, and the light-cone normalization of `a, a†` are not reproduced in
Parts I–IV. Formulas are correct in structure, scaling and parametric order; assembling them into
publishable expressions still requires one pass against Eqs. (5.1)–(5.3).

## Caveat 2 — what I am *not* doing

I am not summarizing the paper. Parts I–IV build the observable-side formalism the paper does
not contain; Parts V–IX are a research assessment. Part VIII is deliberately unflattering where
the physics warrants it.

## Conventions used throughout

| Symbol | Meaning |
|---|---|
| `ρ^a(x)` | valence color charge density operator, transverse position `x`, adjoint index `a` |
| `[ρ^a(x), ρ^b(y)] = i g f^{abc} δ²(x−y) ρ^c(x)` | valence charge algebra (`su(N_c)` current algebra) |
| `a^{a†}_i(k)` | soft gluon creation operator: transverse momentum `k`, color `a`, transverse polarization `i = 1,2` |
| `[a^a_i(k), a^{b†}_j(q)] = (2π)² δ² (k−q) δ^{ab} δ_{ij}` | soft Fock algebra (light-cone normalization absorbed) |
| `b^a_i(k) = i g (k_i / k²) ρ^a(k)` | Weizsäcker–Williams field (momentum space), `O(g)` |
| `N = ∫_k a^{a†}_i(k) a^a_i(k)` | soft gluon number operator |
| `⟨·⟩_ρ` | expectation in the soft Fock space at **fixed** valence configuration |
| `⟨⟨·⟩⟩` | additionally averaged over the CGC weight functional `W[ρ]` |
| `S_⊥` | transverse area; `Q_s` saturation scale; `μ²` MV color charge density |
| dilute counting | `ρ = O(1)`, `g` explicit |
| dense counting | `g ρ = O(1)`, i.e. `Q_s >> Λ_QCD` |

Symbols introduced later, listed here so they are not mistaken for the paper's own notation:
`w`, `w̄` = Nambu-space displacement vector `(α, α*)` and its conjugate; `K` = Nambu covariance;
`Γ(k,q)` = connected two-gluon correlation; `Θ` = the `ρ`-stripped splitting kernel;
`Ξ` = the anti-hermitian generator in `Ω = exp(Ξ)`; `Z[ρ]` = the normalization of `Ω`.

Powers of `g` are quoted in **dilute counting** unless stated. Where the distinction matters —
and it matters a great deal for the publication case — both countings are given.

## File map

| File | Content |
|---|---|
| `01-formalism.md` | Part I — structure of `Ω`, Bogoliubov reduction, what is needed for observables |
| `02-moments.md` | Part II — `⟨N⟩ … ⟨N⁴⟩` |
| `03-cumulants.md` | Part III — `κ₂, κ₃, κ₄`, factorial cumulants, what is genuinely new |
| `04-correlations.md` | Part IV — `n`-gluon distributions, generating function |
| `05-models.md` | Part V — model calculations and their integrals |
| `06-numerics.md` | Part VI — proposed plots |
| `07-literature.md` | Part VII — literature survey, what is already done |
| `08-publication.md` | Part VIII — ranked, critical publication assessment |
| `09-roadmap.md` | Part IX — concrete research plan |
| `10-jhep-assessment.md` | Addendum — will the paper itself make it into JHEP? |
