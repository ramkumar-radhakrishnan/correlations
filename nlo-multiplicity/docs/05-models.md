# Part V — Model calculations before the full CGC average

Ordered from "do this in an afternoon" to "this is the actual paper". For each: what to derive,
what the integrals are, the cost, and what gets plotted.

Throughout, `A₂` and `B₂` are the paper's Eqs. (5.2) and (5.3).

---

## Model 1 — Fixed classical source (single configuration)

**Setup.** `ρ^a(x)` a fixed c-number field; no averaging at all. `[ρ,ρ] → 0`.

**Derive.**
```
α_k = A₂,k[ρ] = i g (k_i/k²) ρ^a(k) × (norm)
m_{kq} = B₂,{kq}[ρ]
C₂ = 2Re ∫_{k,q} A₂*_k B₂,{kq} A₂*_q + ∫_{k,q} |B₂,{kq}|²
```
Everything is a deterministic functional. `P(N)` from the Fredholm determinant.

**Integrals.** Two 2-D transverse integrals ⇒ **4-D**. With colour and polarization sums: a 4-D
quadrature times `O((N_c²−1)²·4)` index terms.

**Cost.** Seconds. Adaptive quadrature or a `256²` grid.

**Plot.** `P(N)` versus `N` for a single "event", overlaid with the Poisson of the same mean. The
message: *how far from Poisson does one configuration get?* This is the cleanest possible display
of the new physics because there is no colour averaging to hide behind.

**Value.** This is the single most important sanity plot in the whole program and it is nearly
free. **Do it first.** If `C₂` comes out numerically negligible for a realistic single
configuration, the project's headline claim is in trouble and you want to know that on day one.

---

## Model 2 — Toy one-dimensional source

**Setup.** Replace `d²k → dk`, `1/k² → 1/k`, one colour. An Abelian scalar toy.

**Derive.** Everything analytically. `A₂(k) = g ρ(k)/k`, choose `B₂(k,q) = g² ρ(k+q) f(k,q)/(kq)`
for a model splitting kernel `f`.

**Integrals.** 2-D. Often closed form.

**Cost.** Instant; much of it symbolic.

**Plot.** `C₂` versus the IR cutoff; `C₂` versus source width. Demonstrates the sign structure of
the interference term without any colour algebra.

**Value.** High as a *debugging and intuition* tool, zero as a publishable result. Use it to
verify the Fredholm-determinant code against exact answers before turning on 2-D.

---

## Model 3 — Finite number of point colour charges

**Setup.** `ρ^a(x) = Σ_{n=1}^{N_s} g T^a_n δ²(x − x_n)`, with `N_s` charges at fixed or random
positions and random colour orientations.

**Derive.**
```
A₂,i^a(k) = i g Σ_n (k_i/k²) T^a_n e^{−ik·x_n}
```
The `1/N_s` suppression of connected correlations becomes explicit: `c₂ ∝ 1/N_s` up to
interference phases.

**Integrals.** Same 4-D, but the source integrals are replaced by finite sums ⇒ **the transverse
integrals become analytic** for point sources, leaving only sums over `n`. `O(N_s⁴)` colour terms
for `C₂`.

**Cost.** `N_s ≤ 100`: milliseconds. `N_s = 10⁴`: still trivial with FFT.

**Plot.** `c₂ · N_s` versus `N_s` — should approach a constant, verifying the `1/N_dof` scaling
law. Also `c₂` vs `N_s` at fixed total charge, separating "more sources" from "stronger sources".

**Value.** This is the cleanest demonstration that the effect scales as advertised, and it
pre-empts a referee asking "is your effect just `1/N`?". **Include it as a figure.**

---

## Model 4 — MV / Gaussian colour charge distribution

**Setup.** `W[ρ] ∝ exp{ −∫d²x ρ^a(x)ρ^a(x) / (2 g²μ²) }`, so
`⟨ρ^a(x)ρ^b(y)⟩ = g²μ² δ^{ab} δ²(x−y)`.

**Derive.** Two equivalent routes:

**(4a) Analytic Wick reduction.** Every correlator factorizes into pairings. For `C₂` you need
`⟨ρ⁴⟩` = 3 pairings. Each pairing produces a distinct transverse convolution. Schematically:

```
⟨ 2Re ∫ A₂*_k B₂,{kq} A₂*_q ⟩  →  g⁶ μ⁴ S_⊥ ∫ d²k d²q  (k_i q_j / k²q²) Θ_{ij}(k,q) × [colour factor]
```
with `Θ` the `ρ`-stripped splitting kernel. Colour factors: `f^{abc}f^{abc} = N_c(N_c²−1)`,
`Tr(T^aT^b) = …` — a finite table, ~10 entries.

**Integrals.** After the delta functions from `⟨ρρ⟩` are used, `C₂` reduces to **one or two 2-D
integrals ⇒ 2-D to 4-D**. `C₃` would be 4-D to 8-D; `C₄` up to 12-D.

**Cost.** `C₂`: adaptive 4-D quadrature, seconds to minutes; or Vegas with `10⁶` points. `C₃`,
`C₄`: Vegas with `10⁷–10⁸` points, hours — and with severe cancellation, so error bars will be
the limiting factor. **This is why I recommend route (4b) for anything above `C₂`.**

**(4b) Direct lattice sampling — the recommended production method.**
```
for each of N_conf configurations:
    sample ρ^a(x) as Gaussian white noise on an L×L transverse lattice
    FFT  →  ρ^a(k)
    build α(k) = A₂[ρ], and the kernel m(k,q) = B₂[ρ]
    compute μ = Σ|α|², C₂[ρ] from the O(g⁴) formula
accumulate the ensemble cumulants of the full P(N) via the law of total cumulance
```
**Cost.** `O(N_conf · L² log L)` for the displacement; the kernel `m(k,q)` is the bottleneck at
`O(L⁴)` if stored densely. Mitigations: (i) `m` is a convolution in `k+q` for the `ρρ` piece, so
it is `O(L² log L)` with FFT; (ii) evaluate `Tr(Kᵖ)` stochastically with Hutchinson trace
estimators (random probe vectors) rather than forming `K` — this is standard in lattice QCD and
reduces the cost to `O(N_probe · L² log L)`.

With `L = 512`, `N_conf = 10⁵`, `N_probe = 20`: **a few CPU-hours, or minutes on a GPU.** This is
entirely feasible for a single author. Higher cumulants need `N_conf ~ 10⁶–10⁷`, still feasible.

**Plot.** `c₂` vs `Q_s` (=`g²μ√S_⊥`-ish); `P(N)` vs Poisson and vs NBD; KNO function.

**Value.** **This is the paper.** MV is the standard, referees expect it, and it is where the
`1/k_NBD` connection is made.

---

## Model 5 — Translationally invariant source

**Setup.** Take `S_⊥ → ∞` with `μ²` fixed; all correlators diagonal in momentum,
`⟨ρ(k)ρ(q)⟩ = (2π)²δ²(k+q) g²μ² S_⊥`.

**Effect.** Each momentum delta function removes 2 dimensions. `C₂` collapses from 4-D to **2-D**;
`C₃` from 8-D to 4-D; `C₄` from 12-D to 6-D.

**Cost.** `C₂`: milliseconds. **`C₄`: 6-D Vegas, minutes.** This is the *only* setting in which the
fourth cumulant is comfortably computable.

**Caveat.** Translational invariance kills all impact-parameter/geometry fluctuations, which are a
large part of the physical `Var_ρ[μ]`. So this model gives you the *quantum* part cleanly and
throws away the classical part. **That is actually a feature**: it is the ideal setting for
isolating the new `O(g⁴)` effect. Present it as "the quantum contribution in the
translationally-invariant limit", then add geometry in Model 4.

**Plot.** `c₂, c₃, c₄` versus `k_⊥` and versus `Q_s`, with the classical contribution switched off
by construction.

---

## Model 6 — Large `N_c`

**Setup.** Keep only leading `N_c` in every colour contraction.

**Effect.** The colour-factor table collapses to a handful of terms; `f^{abc}f^{abd} = N_c δ^{cd}`
dominates and quadrupole-type contractions factorize into dipoles.

**Cost.** Reduces the colour bookkeeping by roughly an order of magnitude; no change in integral
dimensionality.

**Value.** Essential as a **cross-check**, not as the main result. The `1/N_c²` corrections to
`c₂` are interesting in their own right (they are `O(10%)` for `N_c=3`) and a plot of
`c₂(N_c)/c₂(∞)` is a cheap extra figure. **Do not present large-`N_c` as the primary result** — for
multiplicity fluctuations the subleading-`N_c` terms carry the Bose-enhancement structure, and a
referee will ask about exactly that.

---

## Model 7 — Gaussian (non-local) colour charge distribution / MV with `Q_s`-fluctuations

**Setup.** `⟨ρρ⟩ = G(x−y)` with a nontrivial correlator, or MV with `μ²` itself fluctuating
event-by-event (`μ² → μ²(1 + δ)` with `δ` Gaussian or log-normal).

**Why.** `μ²`-fluctuations are the dominant known source of `κ₃, κ₄` in phenomenology and they are
*not* in the standard MV weight. Including them is what makes contact with the "fluctuating `Q_s`"
literature, and it lets you show explicitly that the quantum `C₂` survives on top of them.

**Integrals.** One extra 1-D integral over `δ` per observable.

**Cost.** Negligible.

**Plot.** `c₂^total` decomposed into quantum and classical parts as a function of the width of the
`Q_s` fluctuation. **This is your defence against the "your effect is swamped" referee report** —
it quantifies exactly when it is swamped and when it is not.

---

## Complexity summary

| Model | `C₂` dim | `C₃` dim | `C₄` dim | Cost for `C₂` | Cost for `C₄` |
|---|---|---|---|---|---|
| 1. Fixed source | 4 | 6 | 8 | seconds | minutes |
| 2. 1-D toy | 2 | 3 | 4 | instant | instant |
| 3. `N_s` point charges | sums | sums | sums | ms | seconds |
| 4. MV, analytic Wick | 4 | 8 | 12 | minutes | hours–days ✗ |
| 4b. MV, lattice sampling | — | — | — | CPU-hours | CPU-hours ✓ |
| 5. Translationally inv. | 2 | 4 | 6 | ms | minutes ✓ |
| 6. Large `N_c` | as above | as above | as above | ÷10 colour | ÷10 |
| 7. `Q_s`-fluctuating MV | 4+1 | 8+1 | 12+1 | minutes | via 4b |

**Recommended sequence: 2 → 1 → 3 → 5 → 4b → 7 → 6.** That order front-loads all the debugging and
puts the expensive production run (4b) after every analytic cross-check exists.

---

## Two technical warnings

1. **The `k`-integrals are IR divergent** at LO (`∫d²k/k²`). Every quantity must be defined with an
   explicit IR regulator, and you must show which combinations are IR-finite. `c₂ = C₂/⟨N⟩²`
   involves a ratio of divergent quantities — **check whether the divergences cancel in the
   ratio.** If they do, that is a result worth stating. If they do not, `c₂` is cutoff-dependent
   and must be quoted at fixed `k_⊥` rather than integrated. My expectation: the ratio is *not*
   IR-finite, because `C₂` and `⟨N⟩²` have different IR weights. **Plan for `k_⊥`-differential
   observables from the start.**
2. **Cancellations in `C_p`.** `C₄` is a difference of quantities of order `⟨N⟩⁴`; naive Monte
   Carlo loses ~4 digits. Always compute `C_p` from the *cumulant* estimator directly (or from the
   Fredholm determinant), never as a difference of sampled moments.
