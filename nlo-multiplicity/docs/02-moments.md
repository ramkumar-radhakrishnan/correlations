# Part II — Moments of the gluon multiplicity operator

## 2.0 Strategy

Do **not** compute `⟨N⟩, ⟨N²⟩, ⟨N³⟩, ⟨N⁴⟩` by brute-force Wick contraction of `Ω† N^p Ω`. The
combinatorics grow like `(2p)!` and the `N^p` products generate a forest of contact terms from
`[a, a†]` that all cancel in the end.

Instead:

1. Compute the **factorial moments** `F_p = ⟨:N^p:⟩` — normal-ordered, so no contact terms.
2. Get them from the **generating function**, which at fixed `ρ` is a closed-form Fredholm
   determinant (§2.2). This is exact, not perturbative.
3. Convert to ordinary moments by Stirling numbers (§2.5).

This reduces "compute the fourth moment" from a diagrammatic nightmare to evaluating
`Tr(K⁴)` and `w̄K³w` for two known kernels.

---

## 2.1 The three objects

Everything reduces to (fixed `ρ`; index `k` collects `{k_⊥, a, i}`; `∫_k ≡ ∫ d²k/(2π)²` with
colour/polarization sums):

```
α_k      =  ⟨Ψ| a_k |Ψ⟩                            =  A₂,k         + O(g³)      [O(g)]
m_{kl}   =  ⟨Ψ| a_k a_l |Ψ⟩ − α_k α_l              =  B₂,{kl}      + O(g⁴)      [O(g²)]
n_{kl}   =  ⟨Ψ| a†_k a_l |Ψ⟩ − α*_k α_l            =  (B₂†B₂)_{kl} + O(g⁶)      [O(g⁴)]
```

`A₂` and `B₂` are the paper's coefficients, Eqs. (5.2) and (5.3); see `00-scope §Status` for the
dictionary. `A₂ᵇⱼ(p) = −√2 g pʲρᵇ(−p)/(√p⁺ p²)` is the Weizsäcker–Williams field.

Derivation: these are read straight off `Ω|0⟩ = N|0⟩ + ∫A₂a†|0⟩ + ∫∫B₂a†a†|0⟩` (Part I §1.5),
which is exact through `O(g²)` because the fully normal-ordered `C`, `D`, `E` structures all
annihilate the soft vacuum. Equivalently, in Bogoliubov language
`Ω†a_kΩ = ∫_q[u_{kq}a_q + v_{kq}a†_q] + α_k` with `α = A₂ + O(g³)`, `v = B₂ + O(g⁴)`,
`u = 1 + O(g²)`; unitarity (`A₁ = −A₂*`, `B₁ = −B₂*`, Eqs. 4.4–4.5) guarantees `uu† − vv† = 1`.

**Order hierarchy — the central fact.**
```
m = O(g²)     is ONE power of the pair amplitude
n = O(g⁴)     is TWO powers
```
Consequently the leading fluctuation effect is the **interference of the classical
Weizsäcker–Williams displacement with the quantum pair-emission amplitude**, not the modulus
squared of the two-gluon amplitude. This determines which diagrams matter and is the reason the
promised `O(g^4)` calculation is feasible with an `O(g²)` operator.

## 2.2 Exact generating function at fixed `ρ`

Assemble the Nambu-space displacement and covariance

```
w  =  ( α , α* )ᵀ ,      w̄ = ( α* , α ) ,      K  =  [ [ n , m ] , [ m* , nᵀ ] ]
```

Then for the displaced squeezed (Gaussian) soft state,

```
                                                        ⎧ (z−1)                     ⎫
G(z) ≡ ⟨ z^N ⟩_ρ  =  det[ 1 − (z−1) K ]^(−1/2)  ·  exp ⎨ ─────  w̄ [1 − (z−1)K]^(−1) w ⎬
                                                        ⎩   2                       ⎭
```

This is **exact for Gaussian `Ω|0⟩`**, i.e. exact at fixed classical `ρ` through `O(g²)`, to all
orders in the kernels. Checks:

| Limit | `G(z)` | Distribution |
|---|---|---|
| `K = 0` (coherent) | `exp[(z−1)|α|²]` | Poisson ✓ |
| `α = 0`, `m = 0` (thermal) | `det[1 + (1−z)n]^(−1)` | Bose–Einstein ✓ |
| single-mode squeezed vacuum | `[cosh²r − z² sinh²r]^(−1/2)` | even-`n` only ✓ |
| displaced thermal | `(1−tn)^(−1) exp[t|α|²/(1−tn)]` | Laguerre ✓ |

All four limits reproduce the textbook quantum-optics results, so the formula is verified.

**Immediate corollary (factorial cumulants, all orders):** writing `t = z−1`,

```
log G(1+t)  =  −½ Tr log(1 − tK)  +  (t/2) w̄ (1 − tK)^(−1) w
            =  ½ Σ_{p≥1} tᵖ Tr(Kᵖ)/p   +   ½ Σ_{p≥1} tᵖ w̄ K^(p−1) w
```

hence

```
┌──────────────────────────────────────────────────────────────┐
│   C_p  =  (p−1)!/2 · Tr(Kᵖ)   +   p!/2 · w̄ K^(p−1) w        │
└──────────────────────────────────────────────────────────────┘
```

Since `K = O(g²)` and `w = O(g)`, **`C_p = O(g^{2p})`** for every `p`. This is the clean
statement of the order counting the paper's Introduction is pointing at.

## 2.3 First moment `⟨N⟩`

```
⟨N⟩_ρ  =  C_1  =  ½ Tr K  +  ½ w̄w  =  Tr n  +  |α|²
```

Expanded in `g`:

```
⟨N⟩_ρ  =  ∫_k |A₂^{(1)}_k|²                                       O(g²)   ⟨ρρ⟩
        +  2 Re ∫_k A₂^{(1)*}_k A₂^{(2)}_k                        O(g³)   ⟨ρρρ⟩
        +  ∫_k |A₂^{(2)}_k|²  +  2Re∫_k A₂^{(1)*}_k A₂^{(3)}_k
           +  ∫_{k,q} |B₂,{kq}|²                                  O(g⁴)   ⟨ρ⁴⟩
```

**Reading of the three orders.**

* `O(g²)` — the classical Weizsäcker–Williams gluon number,
  `⟨N⟩_LO = ∫ d²k/(2π)² (g²/k²) ⟨ρ^a(k)ρ^a(−k)⟩ · (N_c²−1)`-type. Log-divergent in the IR
  (regulated by `Q_s` once saturation is resummed) and in the UV (regulated by the `k⁺` cutoff /
  matching to the hard sector). **Not new.**
* `O(g³)` — requires `⟨ρρρ⟩`. **Vanishes identically for any C-even weight functional, in
  particular MV and Gaussian JIMWLK.** Nonzero only in the presence of an odderon-type odd
  correlator `⟨ρ^aρ^bρ^c⟩ ∝ d^{abc}`. This is a small, clean, genuinely unexploited calculation.
* `O(g⁴)` — three physically distinct pieces:
  1. `|A₂^{(2)}|²`: two-source emission squared — classical, enhanced by `ρ⁴`;
  2. `2Re A₂^{(1)*}A₂^{(3)}`: virtual/self-energy correction interfering with the classical field —
     **this is where the `ln(1/x)` rapidity divergence lives**, and it must exponentiate into
     JIMWLK acting on the LO result. Verifying that is a mandatory consistency check;
  3. `∫|B₂|²`: real pair emission — the genuinely new quantum piece.

**Consistency requirement (do this first, before anything else).** The rapidity-divergent part of
(2) plus the real emission (3) must combine into `∂_Y ⟨N⟩_LO = H_JIMWLK ⊳ ⟨N⟩_LO`. If it does not,
either the kernel extraction or the cutoff scheme is wrong. This check is cheap and it is the
single most convincing technical validation you can put in a paper.

## 2.4 Higher factorial moments

From `F_p` in terms of `C_p` (exponential/Bell relations):

```
F_1  =  C_1
F_2  =  C_2 + C_1²
F_3  =  C_3 + 3 C_1 C_2 + C_1³
F_4  =  C_4 + 4 C_1 C_3 + 3 C_2² + 6 C_1² C_2 + C_1⁴
```

with, from §2.2, at leading nontrivial order in `g`:

```
C_2  =  Tr(n²) + Tr(m m†) + 2 α†nα + 2 Re (α* m α*)

      →  ∫_{k,q} |B₂,{kq}|²  +  2 Re ∫_{k,q} A₂*_k B₂,{kq} A₂*_q       [ O(g⁴), complete ]

C_3  =  Tr(K³) + 3 w̄ K² w                                             [ O(g⁶), INCOMPLETE — §2.6 ]

C_4  =  (3!/2) Tr(K⁴) + 12 w̄ K³ w                                     [ O(g⁸), INCOMPLETE — §2.6 ]
```

Written out, the complete `O(g⁴)` result is

```
┌─────────────────────────────────────────────────────────────────────────────────┐
│  C₂ = 2 Re ∫_{k,q} A₂^{a*}_i(k) B₂^{ab}_{ij}(k,q) A₂^{b*}_j(q)                  │
│         +   ∫_{k,q} | B₂^{ab}_{ij}(k,q) |²        +  O(g⁶)                      │
└─────────────────────────────────────────────────────────────────────────────────┘
```

This is the whole `O(g^4)` fluctuation calculation, at fixed `ρ`, in one line. Everything after
this is colour averaging and numerics.

**Which colour correlators appear.** With `A₂ ∝ gρ` and `B₂ ∝ g²ρρ ⊕ g²ρ`:

| Term | `ρ` content | MV status |
|---|---|---|
| `Re A₂*B₂A₂*` with `B₂ ∝ ρρ` | `⟨ρρρρ⟩` | 3 Wick pairings |
| `Re A₂*B₂A₂*` with `B₂ ∝ ρ` (`g→gg`) | `⟨ρρρ⟩` | **vanishes** in MV |
| `∫|B₂|²` with `B₂ ∝ ρρ` | `⟨ρρρρ⟩` | 3 pairings |
| `∫|B₂|²` with `B₂ ∝ ρ` | `⟨ρρ⟩` | trivial |
| `ρ`-commutator reorderings | `g⟨ρρρ⟩ → ⟨ρρ⟩` | trivial, but **must be kept** |

Note the last row: the current-algebra commutator converts a vanishing MV correlator into a
surviving one. Dropping `[ρ,ρ]` is therefore **not** a harmless simplification at this order.
(The paper's Eq. (4.28)–(4.29) does exactly this decomposition.)

## 2.5 Ordinary moments

```
⟨N⟩    =  F₁
⟨N²⟩   =  F₁ + F₂
⟨N³⟩   =  F₁ + 3F₂ + F₃
⟨N⁴⟩   =  F₁ + 7F₂ + 6F₃ + F₄
```
(Stirling numbers of the second kind: `⟨Nᵖ⟩ = Σ_j S(p,j) F_j`.)

**Do not quote ordinary moments as results.** They are dominated by trivial powers of `⟨N⟩` —
`⟨N⁴⟩ ≈ ⟨N⟩⁴` to a part in `10⁴` — so all the physics sits in a cancellation of four leading
digits. Publish `C_p`, not `⟨N^p⟩`. This is also the numerically stable choice (Part VI).

## 2.6 What is complete and what is not — read this before proposing the paper

| Quantity | Leading order | Complete with the `O(g²)` `Ω`? |
|---|---|---|
| `⟨N⟩` | `O(g²)` | ✔ complete through `O(g⁴)` |
| `C₂` (variance beyond Poisson) | `O(g⁴)` | ✔ **complete** |
| `C₃` | `O(g⁶)` | ✘ **incomplete** |
| `C₄` | `O(g⁸)` | ✘ **incomplete** |

**Why `C₃` is incomplete.** `Ω|0⟩` has no three-gluon component at `O(g²)` (Part I §1.5), so the
three-gluon amplitude first arises at `O(g³)` — and it contributes to `C₃` at
`O(α*³ T) = O(g³·g³) = O(g⁶)`, i.e. *the same order* as the Gaussian contribution `w̄K²w`. The
`O(g²)` operator simply does not contain it. By contrast `C₂` is safe: a cubic amplitude first
corrects `m` at `O(g⁴)`, which is subleading to `m = B₂ = O(g²)`.

**Consequences for the project.** You can honestly publish:

* `⟨N⟩` to `O(g⁴)`, with the JIMWLK-exponentiation check;
* `C₂` (equivalently the variance) **complete** at `O(g⁴)` — this is the headline;
* `C₃, C₄` **only** in their classical colour-fluctuation part (§2.7), which is complete at its
  own leading order because it needs nothing but the LO displacement;
* the *Gaussian part* of `C₃, C₄` as a partial result, clearly labelled as such.

You cannot honestly publish "the third and fourth cumulants at NLO". A referee who knows LCPT
will spot the missing three-gluon amplitude immediately. **Either restrict the claim to `C₂`, or
extend `Ω` to `O(g³)` first** — see Part VIII, Project B.

## 2.7 The colour average, and why it does not rescue the higher cumulants (but is still needed)

Apply the law of total cumulance (§1.8). The `κ_p^ρ[μ]` terms — cumulants of the *fixed-`ρ` mean*
`μ(ρ) = ⟨N⟩_ρ` over the CGC ensemble — need only the LO displacement, so they are complete at
their leading order for every `p`. In MV,

```
μ(ρ)  =  ∫_k (g²/k²) ρ^a(k) ρ^a(−k) / (2π)²   →   κ_p^ρ[μ]  ∝  g^{2p} · (contractions of ρ^{2p})
```

and each connected cumulant is suppressed by one power of the number of degrees of freedom
`N_dof ≃ (N_c²−1) S_⊥ Q_s²` per order:

```
κ_p^ρ[μ] / μᵖ   ~   1 / N_dof^{p−1}
```

**Crucial and uncomfortable fact:** the classical term `κ_p^ρ[μ]` and the quantum term
`E_ρ[C_p^Fock]` are **the same order in `g`** and carry **the same `1/N_dof^{p−1}` suppression**.
They cannot be separated by power counting. Anyone claiming "we computed the new quantum
fluctuation" must show a *structural* discriminator. There is one, and it is good:

> **Momentum-space signature.** Bose enhancement (the classical/`n`-type correlation) peaks at
> `k ≃ +q` — a near-side, same-momentum enhancement. Pair emission from the squeezing kernel
> (`m`-type) peaks at `k ≃ −q` — **back-to-back**, because the two gluons are emitted from the
> same vertex and share the source recoil. Plotting the two-gluon correlation as a function of
> the angle between `k` and `q` separates them cleanly.

That observation is, I think, the strongest single piece of new physics available in this
program, and it survives the "your effect is swamped" objection because the two effects live at
different angles rather than at different magnitudes.

---

## Summary of Part II

* Use `G(z) = det[1−(z−1)K]^{−1/2} exp{…}` — exact at fixed `ρ`, verified against four
  quantum-optics limits.
* `C_p = (p−1)!/2 · Tr(Kᵖ) + p!/2 · w̄K^{p−1}w`, so `C_p = O(g^{2p})`.
* `⟨N⟩` complete to `O(g⁴)`; `C₂` complete at `O(g⁴)`; `C₃, C₄` **not** complete without the
  `O(g³)` coefficients of `Ω`.
* The `O(g³)` term in `⟨N⟩` needs `⟨ρρρ⟩` and vanishes in MV — an odderon-sensitive observable.
* Classical and quantum fluctuations are degenerate in `g` and in `1/N_dof`; separate them by the
  near-side (`k≃+q`) vs back-to-back (`k≃−q`) momentum structure.
