# Part III — Connected correlations and cumulants

## 3.1 Conversion table

With `κ_n = Σ_p S(n,p) C_p` (Stirling numbers of the second kind):

```
κ₁  =  C₁
κ₂  =  C₁ +  C₂
κ₃  =  C₁ + 3C₂ +  C₃
κ₄  =  C₁ + 7C₂ + 6C₃ + C₄
κ₅  =  C₁ + 15C₂ + 25C₃ + 10C₄ + C₅
```

Read these as: **a Poisson distribution has `C_p = 0` for `p ≥ 2`, hence `κ_n = ⟨N⟩` for all `n`.**
Every deviation of `κ_n` from `⟨N⟩` is a factorial cumulant, i.e. a genuine correlation.

This is why the factorial cumulants `C_p` — not `κ_n` — are the right observables here. They are
the *background-free* quantities: they are exactly zero for the leading-order coherent state, so
any nonzero value is unambiguously an NLO or a colour-fluctuation effect. Reporting `κ₄` instead
means reporting a number whose leading contribution is `⟨N⟩` and whose physics is a per-mille
correction on top.

**Recommendation: make `C₂` the headline observable of the follow-up paper, and quote
`κ₂ = ⟨N⟩ + C₂` only for contact with the heavy-ion cumulant literature.**

## 3.2 The cumulants, order by order

### `κ₂` — variance

```
κ₂  =  ⟨N⟩  +  C₂
     =  ⟨N⟩  +  [ 2Re ∫ A₂* B₂ A₂*  +  ∫|B₂|² ]_ρ-averaged  +  Var_ρ[ μ(ρ) ]
        └─LO─┘   └──── quantum, O(g⁴), NEW ────┘            └─ classical, O(g⁴), known ─┘
```

* **`⟨N⟩` (Poisson part).** Present at LO. Trivially the coherent-state result. Not new.
* **`Var_ρ[μ]`.** Colour-charge fluctuations of the classical WW field. This is the CGC origin of
  the negative binomial distribution (Gelis–Lappi–McLerran "Glittering Glasmas"). Known.
* **`2Re∫A₂*B₂A₂* + ∫|B₂|²`.** The genuinely new `O(g⁴)` piece. Note the first term is an
  **interference** and therefore **not sign-definite** — the fixed-source distribution can be
  driven sub-Poissonian. That would be a striking and quotable result, because every existing
  CGC multiplicity prediction is super-Poissonian (NBD). If a region of `k_⊥` or `Q_s` gives
  `C₂ < 0`, that is a headline.

> **Physical meaning of `κ₂`:** event-by-event fluctuation of the soft gluon number at fixed
> valence configuration plus fluctuation of the valence configuration itself. `κ₂/⟨N⟩ = 1`
> exactly for a coherent state; `> 1` for Bose-enhanced/thermal-like emission; `< 1` for
> destructive displacement–squeezing interference (a genuinely non-classical, "sub-Poissonian
> light" situation with no analogue in the classical CGC).

### `κ₃` — skewness

```
κ₃  =  ⟨N⟩ + 3C₂ + C₃
```

* `C₃` is `O(g⁶)` and **incomplete** with the `O(g²)` operator (Part II §2.6): the missing
  three-gluon amplitude contributes at the same order.
* The `3C₂` term *is* complete at `O(g⁴)` and is not negligible — in fact for `N_dof >> 1`,
  `3C₂ >> C₃`, since `C_p/⟨N⟩ᵖ ~ N_dof^{1−p}`. So **`κ₃` is numerically dominated by pieces you
  can compute completely**, even though `C₃` itself is not.

  This is a legitimate and defensible publication strategy: *quote `κ₃` complete through
  `O(g⁴)`, i.e. `κ₃ = ⟨N⟩ + 3C₂ + O(g⁶)`, and state that the `O(g⁶)` remainder is suppressed by
  `1/N_dof`.* That is honest and it is a real result.
* The classical part `κ₃^ρ[μ]` (third cumulant of the WW field over the colour ensemble) is
  complete and computable; in MV it is the `⟨ρ⁶⟩` connected topology.

> **Physical meaning of `κ₃`:** asymmetry of `P(N)`. Positive skew is the hallmark of a
> fluctuating-source (compound Poisson / NBD) picture: rare large-`Q_s` configurations produce a
> long high-multiplicity tail. A measured/computed skew *smaller* than the NBD value would
> indicate that quantum pair correlations partially cancel the source fluctuations.

### `κ₄` — kurtosis

```
κ₄  =  ⟨N⟩ + 7C₂ + 6C₃ + C₄
```

Same structure: complete through `O(g⁴)` (`⟨N⟩ + 7C₂`), incomplete at `O(g⁶)` and beyond.

> **Physical meaning of `κ₄`:** tail weight. In heavy-ion phenomenology `κ₄` and the ratios
> `κ₄/κ₂`, `κ₃/κ₂` are the standard tools for identifying non-trivial correlation structure
> (they are the same observables used for baryon-number fluctuations and critical-point
> searches). A CGC-based prediction for the *initial-state* `κ₄/κ₂` of gluons is a genuinely
> useful input to that discussion, because it fixes the non-thermal baseline.

### Higher cumulants

`C_p` for `p ≥ 5` are all `O(g^{2p})` and increasingly incomplete. They are **not** worth
computing from this operator. The only exception: the purely classical tower `κ_p^ρ[μ]` (colour
fluctuations of the LO WW field) can be computed to all `p` with no extra effort once you sample
`ρ` on a lattice — you get the whole distribution for free. Report that as a distribution, not as
a list of cumulants.

## 3.3 Which cumulants are nonzero at leading order

| Object | LO (coherent state, fixed `ρ`) | LO + colour average | NLO (`O(g⁴)`) |
|---|---|---|---|
| `⟨N⟩` | `∫\|A₂\|²` ≠ 0 | ≠ 0 | correction |
| `C₂` | **0** | `Var_ρ[μ]` ≠ 0 | **first genuine quantum contribution** |
| `C₃` | **0** | `κ₃^ρ[μ]` ≠ 0 | partial |
| `C₄` | **0** | `κ₄^ρ[μ]` ≠ 0 | partial |
| `κ₂` | `= ⟨N⟩` (Poisson) | `> ⟨N⟩` (NBD) | `≷` |
| `κ₃, κ₄` | `= ⟨N⟩` | `> ⟨N⟩` | `≷` |

**The single sharpest statement in the whole program:**

> At leading order and at fixed valence colour charge, the soft gluon multiplicity distribution is
> *exactly Poissonian*: all factorial cumulants vanish. The first departure from Poisson
> statistics arises at `O(g⁴)` and is controlled entirely by the pair-emission (squeezing)
> coefficient `B₂` of the NLO evolution operator, interfering with the classical
> Weizsäcker–Williams displacement `A₂`.

That is a clean, falsifiable, one-sentence abstract for a paper. It is also exactly what the
Introduction of arXiv:2607.18373 promises to deliver.

## 3.4 What contains information beyond the coherent-state approximation

Rank-ordered by how much is genuinely beyond LO:

1. **The anomalous correlator `m = ⟨aa⟩_c = B₂`.** Zero for any coherent state. It is the *only*
   `O(g²)` object in the theory that is invisible at LO. Everything new in `C₂` traces back to it.
2. **The back-to-back structure of `m_{k,q}`** (peaked at `q ≈ −k`). No coherent state and no
   Bose-enhancement mechanism produces this. See Part IV §4.3 — this is the discriminator.
3. **The sign of `C₂`.** Coherent → 0; classical fluctuating source → strictly `> 0`; quantum
   interference → can be `< 0`. A negative region is unreachable by any classical CGC calculation.
4. **The `⟨ρρρ⟩` (odderon) contribution to `⟨N⟩` at `O(g³)`.** Vanishes in MV; nonzero with a
   C-odd weight. This is orthogonal to everything above and is essentially virgin territory.
5. **`n = B₂†B₂`** — the true "pair occupancy". Beyond LO but `O(g⁴)` in the *kernel*, hence
   `O(g⁸)` in `C₂`. Numerically irrelevant. Do not build a project on it.

## 3.5 The one thing that will get the paper into trouble

**Gluon number in the light-cone wave function is not an observable.** It is gauge-dependent
(defined in `A⁺=0`), scheme-dependent (depends on the `k⁺` cutoff separating soft from valence —
the paper's `Λ` and `∨`), and IR-divergent. A referee will ask, correctly:

> "You have computed the fluctuations of a number operator in a particular gauge and a particular
> factorization scheme. What is measured?"

There are three defensible answers, in increasing order of work and of persuasiveness:

* **(a) It is an intermediate quantity.** Present `C₂` as a property of the CGC wave function,
  analogous to how the unintegrated gluon distribution is presented. Acceptable, but weak, and
  will draw a "limited physics impact" report.
* **(b) Show that the cutoff dependence is JIMWLK.** Demonstrate that `∂⟨N⟩/∂Y` and
  `∂C₂/∂Y` are generated by the JIMWLK Hamiltonian acting on the lower-order results. This makes
  the scheme dependence *physics* rather than an artefact. **Strongly recommended; cheap.**
* **(c) Compute produced-gluon cumulants.** Scatter the dressed state off a target: insert
  Wilson lines, `a → U a`, and compute the multiplicity of *produced* gluons. This converts the
  whole calculation into an observable. It roughly doubles the work (you need `⟨U U U U⟩`-type
  target correlators, i.e. quadrupoles) but it is what turns a technical note into a physics
  paper — and it is the natural consumer of the paper's `C`, `D`, `E` coefficients, which play no
  role in `Ω|0⟩`.

**My assessment: (b) is mandatory, (c) is what separates a PRD from a JHEP.** See Part VIII.

---

## Summary of Part III

* `C₂` is the right headline observable: exactly zero at LO, complete at `O(g⁴)`, sign-indefinite.
* `κ₃, κ₄` are complete only through `O(g⁴)` (`⟨N⟩ + 3C₂`, `⟨N⟩ + 7C₂`); their `O(g⁶)` pieces
  need the `O(g³)` operator. Say so explicitly rather than hoping nobody checks.
* Everything genuinely new traces to the anomalous correlator `m = B₂`, and its distinguishing
  signature is back-to-back momentum structure.
* The wave-function-vs-observable problem is the paper's main vulnerability; fix it with the
  JIMWLK consistency check at minimum.
