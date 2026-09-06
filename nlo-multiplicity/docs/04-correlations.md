# Part IV — Correlation functions

## 4.1 Single-gluon distribution

```
dN/d²k  ≡  Σ_{a,i} ⟨Ψ| a^{a†}_i(k) a^a_i(k) |Ψ⟩
        =  Σ_{a,i} [ |α^a_i(k)|²  +  n^{aa}_{ii}(k,k) ]
```

At LO this is the Weizsäcker–Williams / unintegrated gluon distribution:

```
dN/d²k |_LO  =  (g²/k²) ⟨ ρ^a(k) ρ^a(−k) ⟩ / (2π)²   ~   (N_c²−1) g²μ² S_⊥ / k²
```

with the familiar `1/k²` tail and IR log divergence cut off by saturation at `k ~ Q_s`. Nothing
new; it is the normalization against which everything else is measured.

The `O(g⁴)` correction `n(k,k) = (B₂†B₂)(k,k)` is the **pair-production contribution to the gluon
spectrum**. It is the real-emission part of the NLO unintegrated gluon distribution and it is
worth plotting on its own: it has a *different* `k`-dependence from the LO term (harder, because
pair emission from the `g→gg` vertex populates larger `k`), so the ratio
`n(k,k) / |α(k)|²` rises with `k_⊥`.

## 4.2 Two-gluon correlation

```
dN₂/d²k d²q  =  ⟨ a†_k a†_q a_k a_q ⟩
             =  (dN/d²k)(dN/d²q)  +  Γ(k,q)
```

The connected piece, from the Gaussian Wick expansion of §2.2:

```
Γ(k,q)  =  |n(k,q)|²  +  |m(k,q)|²
          +  2 Re[ α*(k) n(k,q) α(q) ]                 ← Bose enhancement / HBT, peaks at k ≈ +q
          +  2 Re[ α*(k) α*(q) m(k,q) ]                ← pair emission,          peaks at k ≈ −q
```

with the orders

```
|m(k,q)|²          O(g⁴)
2Re[α* α* m]       O(g⁴)      ← LEADING, sign-indefinite
2Re[α* n α]        O(g⁶)
|n(k,q)|²          O(g⁸)
```

**This is the central result of Part IV.** At `O(g⁴)` the connected two-gluon correlation of the
NLO wave function is

```
┌──────────────────────────────────────────────────────────────────────┐
│  Γ(k,q) |_{O(g⁴)}  =  |B₂(k,q)|²  +  2 Re[ A₂*(k) A₂*(q) B₂(k,q) ]   │
└──────────────────────────────────────────────────────────────────────┘
```

### 4.3 Why the angular structure is the discriminator

Both known CGC mechanisms and the new one produce a connected two-gluon correlation of the same
parametric size (`~ ⟨N⟩²/N_dof`), so magnitude alone proves nothing. The angular structure does:

| Mechanism | Correlator | Peaks at | Origin |
|---|---|---|---|
| Bose enhancement in projectile WF | `α* n α` | `k ≈ +q` (near side, `Δφ ≈ 0`) | identical-boson symmetrization; the Altinoluk–Armesto–Beuf–Kovner–Lublinsky ridge mechanism |
| HBT | `α* n α` | `k ≈ +q` | same object, different kinematic limit |
| **Pair emission / squeezing** | `α*α* m` | `k ≈ −q` (**away side, `Δφ ≈ π`**) | both gluons from one vertex; transverse momentum conservation against the source |
| Colour-charge (classical) fluctuations | enters via `Var_ρ` | broad, `Δφ`-independent | fluctuating `Q_s` |

So: **plot `Γ(k,q)` versus `Δφ` at fixed `|k| = |q|`.** The `O(g⁴)` NLO term contributes an
away-side ridge that the glasma-graph/Bose-enhancement calculations do not produce. That is a
concrete, checkable, publishable prediction, and it connects directly to the ridge literature.

**Honest caveat:** an away-side correlation from momentum conservation is *not* by itself
surprising, and referees will point out that momentum-conservation-induced away-side correlations
are ubiquitous and often treated as background. The claim must therefore be quantitative — the
*normalization and `Q_s`-dependence* of the away-side term, not its existence.

## 4.4 Three-gluon correlation

```
dN₃/d²k₁d²k₂d²k₃ = ⟨ a†₁a†₂a†₃ a₁a₂a₃ ⟩
```

Connected part, at leading order in the kernels:

```
Γ₃(k₁,k₂,k₃) |_leading  =  2 Re [ α*₁ α*₂ m₁₂ · (α*₃ … ) ]-type  +  Re[ m₁₂ n₂₃ … ]
                        ~  O(g⁶)
```

**Incomplete** for the same reason `C₃` is: an `O(g³)` three-gluon amplitude in `Ω|0⟩` contributes
a genuine three-gluon emission vertex at `O(g⁶)`. Physically that is the `1→3` splitting, and it
is *precisely* the object one would expect to dominate a three-gluon correlation.

**Recommendation: do not build a project around the three-gluon correlation from this operator.**
The one honest three-gluon statement available is the *Gaussian prediction*: with `Ω|0⟩`
truncated at two-gluon components, `Γ₃` is completely determined by `Γ₂` (no independent
three-body correlation). Stating and testing that — "the NLO wave function predicts Gaussian-state
three-body factorization up to `O(g⁶)`" — is a legitimate, if modest, result, and it is a
genuinely useful null prediction for the ridge community.

## 4.5 Factorial moments and normalized cumulants

```
F_p  =  ⟨N(N−1)…(N−p+1)⟩              factorial moments
C_p                                     factorial cumulants  (§2.2)

f_p  =  F_p / ⟨N⟩ᵖ                     normalized factorial moments   →  1  for Poisson
c_p  =  C_p / ⟨N⟩ᵖ                     normalized factorial cumulants →  0  for Poisson
```

`c_p` is the correct dimensionless observable. Expected scaling:

```
c_p   ~   1 / N_dof^{p−1}  ,        N_dof ≃ (N_c²−1) S_⊥ Q_s² / (2π)²
```

so `c₂ ~ 1/[(N_c²−1)S_⊥Q_s²]`. **This is the same `1/k` that appears in the negative binomial fits
to `pp` and `AA` multiplicity data**, and connecting the NLO calculation to it is the most direct
route to phenomenological relevance:

> `k_NBD = 1/c₂` — computing `c₂` from the NLO wave function is computing the NBD parameter `k`
> from first principles, including its `O(g⁴)` quantum correction.

That framing is much more likely to be read and cited than "we computed a factorial cumulant".

## 4.6 Generating function

Fixed `ρ` (exact through `O(g²)`, from §2.2):

```
G_ρ(z)  =  det[1 − (z−1)K]^{−1/2} exp{ ½(z−1) w̄[1−(z−1)K]^{−1}w }
```

Full distribution:

```
G(z)  =  ∫ Dρ  W_Y[ρ]  G_ρ(z)          P(N) = (1/N!) d^N G/dz^N |_{z=0}
```

Three tractable evaluations:

1. **`K = 0`** (LO): `G_ρ(z) = e^{(z−1)μ(ρ)}`, so `G(z) = ∫Dρ W[ρ] e^{(z−1)μ(ρ)}` — a **compound
   Poisson** distribution. If `μ(ρ)` were Gamma-distributed this gives exactly the negative
   binomial. Deriving the actual `μ`-distribution in MV and showing how close it is to Gamma is a
   clean, self-contained calculation and directly connects to the "KNO scaling from a nearly
   Gaussian action" line of work.
2. **`K ≠ 0`, leading order:** expand `log G_ρ` to `O(g⁴)`; `P(N)` is a Poisson deformed by a
   single factorial-cumulant correction — analytic, and gives the shape distortion in closed form.
3. **Numerical:** sample `ρ`, build `K` and `w`, evaluate the Fredholm determinant on a
   momentum grid, invert the `z`-transform by FFT. Gives full `P(N)` including tails. See Part VI.

**KNO scaling** is directly accessible: `Ψ(N/⟨N⟩) = ⟨N⟩ P(N)`. Whether KNO scaling survives the
`O(g⁴)` correction — i.e. whether `c_p` are `⟨N⟩`-independent — is a well-posed, interesting, and
to my knowledge unanswered question. **KNO holds iff `c_p` depends on `Y` and `Q_s` only through
`⟨N⟩`.** Since `c₂ ~ 1/(S_⊥Q_s²)` and `⟨N⟩ ~ S_⊥Q_s²/α_s`, we get `c₂ ~ α_s/⟨N⟩` — so KNO is
**violated by a running-`α_s` factor** at LO, and the `O(g⁴)` term modifies this. That is a sharp,
quotable prediction and it is testable against the existing KNO literature.

## 4.7 Which observables are worth the effort

| Observable | Theory interest | Experimental contact | Verdict |
|---|---|---|---|
| `dN/d²k` at `O(g⁴)` | medium (NLO ugd) | indirect | do it, it is a byproduct |
| `c₂` = `1/k_NBD` | **high** | **direct** (NBD fits, `pp`/`pA`/`AA`) | **primary target** |
| `Γ(k,q)` vs `Δφ` | **high** | direct (ridge, away-side) | **primary target** |
| KNO violation from `c_p` | high | direct (H1, ALICE, CMS) | strong secondary |
| `Γ₃` | low (incomplete) | weak | null prediction only |
| `⟨N⟩` at `O(g³)`, odderon | medium-high (novel) | weak (needs C-odd observable) | interesting side project |
| `κ₄`, `κ₄/κ₂` | medium | contact with fluctuation programme | secondary |

---

## Summary of Part IV

* Connected two-gluon correlation at `O(g⁴)`: `Γ = |B₂|² + 2Re[A₂*A₂*B₂]`.
* Its angular signature is **away-side** (`k ≈ −q`), distinguishing it from Bose enhancement and
  HBT, which are near-side. This is the cleanest new prediction available.
* Normalized factorial cumulant `c₂ = 1/k_NBD`: a first-principles NLO computation of the negative
  binomial parameter. Best route to phenomenological impact.
* The generating function is a Fredholm determinant at fixed `ρ`, compound-Poisson after colour
  averaging; KNO scaling is violated at the level of `α_s`, computably.
* Three-gluon correlations are not honestly computable here — publish the Gaussian-factorization
  null prediction instead.
