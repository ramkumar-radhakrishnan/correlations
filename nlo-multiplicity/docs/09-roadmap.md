# Part IX — Research roadmap

Three projects, in the order they should be attempted. Times assume one person working on this as
their main project, with the `Ω` derivation already in hand.

---

## Phase 0 — Two weeks. Decide whether to proceed at all.

This phase costs almost nothing and can kill the project before you invest in it. **Do not skip
it.**

| # | Task | Output |
|---|---|---|
| 0.1 | Read [arXiv:2605.03038](https://arxiv.org/abs/2605.03038) (squeezed-state shockwave radiation) and [arXiv:2412.05085](https://arxiv.org/abs/2412.05085) (Born–Oppenheimer RG) in full | Scooping assessment. If either contains the displaced-squeezed identification *and* multiplicity statistics, reframe immediately |
| 0.2 | Read [arXiv:1503.07126](https://arxiv.org/abs/1503.07126) and [arXiv:1804.02910](https://arxiv.org/abs/1804.02910), looking specifically for an `⟨aa⟩`-type (anomalous) correlator | Determines whether the away-side result (Project C) is new |
| 0.3 | Read [Lublinsky–Mulian 1610.03453](https://arxiv.org/abs/1610.03453) (= Ref. [19]), looking for the two-gluon Fock amplitude | Determines how much of `B₂` is already public |
| 0.4 | Read [hep-ph/0609227](https://arxiv.org/abs/hep-ph/0609227) (*One gluon, two gluon*, = Ref. [63]) | Prepare the one-paragraph "how is this different" answer |
| 0.5 | Implement Model 2 (1-D toy) + the Fredholm-determinant generating function; verify against the four quantum-optics limits in §2.2 | Working, validated code. ~2 days |
| 0.6 | Run **Model 5 (translationally invariant)** with the actual `B₂` from Eq. (5.3); evaluate `C₂` and compare to `Var_ρ[μ]` | **THE GO/NO-GO NUMBER.** If quantum/classical `< 10⁻²` everywhere, reframe around angular structure |

**Decision point.** Proceed to Phase 1 only if 0.1–0.4 leave the ground clear and 0.6 is not
catastrophic.

---

## Project A(+C+D) — "Non-Poissonian gluon statistics from the NLO CGC wave function"

### 1. Required derivations

| # | Derivation | Difficulty | Time |
|---|---|---|---|
| A1 | Read `α, m, n` off `Ω\|0⟩ = N\|0⟩ + ∫A₂a†\|0⟩ + ∫∫B₂a†a†\|0⟩`; confirm `C, D, E` drop out | Easy | 1 wk |
| A2 | Verify unitarity relations (`A₁=−A₂*`, `B₁=−B₂*`, `uu†−vv†=1`) against Eqs. (4.4)–(4.10) | Easy — but it is your main internal cross-check | 1 wk |
| A3 | `ρ`-ordering: fix a prescription, compute the commutator-induced `O(g)` reshuffles (cf. the paper's Eqs. 4.28–4.29), prove observable-level ordering independence for `⟨N⟩` at `O(g⁴)` | **Hard, error-prone** | 4–6 wk |
| A4 | `⟨N⟩` to `O(g⁴)`: assemble real (`∫\|B₂\|²`), virtual (`\|N\|²`), and interference terms | Medium | 3 wk |
| A5 | **JIMWLK consistency**: show the `ln(1/x)` divergence of A4 = `H_JIMWLK ⊳ ⟨N⟩_LO` | Medium — high payoff | 3 wk |
| A6 | `C₂ = 2Re∫A₂*B₂A₂* + ∫\|B₂\|²` with full colour algebra; MV Wick reduction of `⟨ρ⁴⟩` | Medium | 4 wk |
| A7 | IR analysis: which of `⟨N⟩`, `C₂`, `c₂`, `C₂(k_⊥)` are IR-finite | Medium — determines the claims | 2 wk |
| A8 | `Γ(k,q)` decomposition into near-side (`α*nα`) and away-side (`α*α*m`) | Easy given A6 | 2 wk |
| A9 | `O(g³)` odderon term (Project D): `2Re∫A₂^{(1)*}A₂^{(2)}` with `⟨ρρρ⟩ ∝ d^{abc}` | Easy | 3 wk |
| A10 | Gaussian-state null prediction for `Γ₃` (three-body factorization up to `O(g⁶)`) | Easy | 1 wk |

**Total derivation time: ~6 months**, dominated by A3.

### 2. Required numerics

| # | Computation | Model | Cost |
|---|---|---|---|
| N1 | Fredholm-determinant `P(N)` code, validated | 2, 1 | done in Phase 0 |
| N2 | `c₂` vs `Q_s`, quantum/classical decomposition (**P3, P2**) | 5, then 4b | CPU-hours |
| N3 | `c₂` vs `Λ_IR` and vs `k_⊥` (**P5, P7**) | 5 | CPU-hours |
| N4 | `Γ(Δφ)` near/away decomposition (**P8**) | 5, 4b | CPU-hours |
| N5 | Sign map in `(k_⊥,Q_s)` (**P12**) | 5 | CPU-hours |
| N6 | Full MV lattice sampling, `P(N)`, KNO (**P1, P10**) | 4b, `L=512`, `N_conf=10⁵` | ~1 GPU-day |
| N7 | `N_s`-scaling validation (**P9**), large-`N_c` check | 3, 6 | minutes |
| N8 | `Q_s`-fluctuation robustness (**P7-model**) | 7 | minutes |

**Total numerics time: ~2 months** including code development, run mostly in parallel with the
derivations.

### 3. Estimated difficulty
**Moderate-to-hard.** No single step is beyond a competent practitioner; the risk is accumulated
algebra error in A3/A6. Mitigate with the 1-D toy (Model 2) and the `N_s`-scaling check (N7) as
independent validations, and with the JIMWLK check (A5) as the global validation.

### 4. Estimated completion time
**6–9 months** to a submittable draft. Add 2–3 months if A3 goes badly (it might).

### 5. Publication target
**Physical Review D**, regular article, ~30 pages with appendices.
Fallback: EPJC. Stretch: JHEP only if Project E material is folded in — do not plan on it.

### 6. Risks

| Risk | Probability | Impact | Mitigation |
|---|---|---|---|
| Scooped by the KLM group or the 2605.03038 authors | Medium | Fatal | Phase 0.1; move fast; post to arXiv early |
| Quantum effect numerically negligible | Medium | High | Phase 0.6 go/no-go; reframe around angular structure |
| `ρ`-ordering algebra error | Medium-high | High | A2 unitarity check; 1-D toy; ordering-independence proof |
| `c₂` IR-divergent, no clean observable | Medium | Medium | A7 early; go `k_⊥`-differential |
| Away-side result turns out to be known | Medium | Medium | Phase 0.2 |
| "Not an observable" referee report | **High** | Medium | A5 (JIMWLK check) + explicit framing; accept one round of revision |
| "Uses only 3 of 6 coefficient families" | Medium | Low | State it; point at Project E |

### 7. Expected scientific impact
**Modest but real.** Realistic citation expectation: 15–40 over five years. It will be read by the
small-`x` formalism community (~50 people) and, *if* the squeezed-state and away-side framing
lands, by the correlations/ridge community (~200 people). The `k_NBD` at NLO result has the
longest legs, because it is the one number a phenomenologist can pick up and use.

It will **not** be a high-impact paper. It is a solid, careful, first-of-its-kind calculation that
establishes a formalism others will use. That is a good outcome for a follow-up to a formalism
paper, and it is the correct thing to write next.

---

## Project E — "Produced-gluon multiplicity cumulants at NLO" (second paper)

1. **Derivations.** Insert target Wilson lines: `a → U a`; compute `⟨U U U U⟩` (quadrupole) in the
   Gaussian/MV approximation; handle the rapidity subtraction and the real–virtual cancellation in
   the production amplitude; assemble `C₂^produced`. This is where the paper's `C`, `D`, `E`
   coefficients finally do work. **Hard.**
2. **Numerics.** Quadrupole correlators on a lattice (standard technology, e.g. existing JIMWLK
   codes); `C₂^produced` vs `k_⊥`, `Y`, target `Q_s`.
3. **Difficulty.** High. This is the contentious part of the NLO CGC literature.
4. **Time.** 18–24 months.
5. **Target.** JHEP or PRD.
6. **Risks.** Negativity/scheme problems endemic to NLO `pA` production; the possibility that the
   result is dominated by known glasma-graph physics; a crowded and critical referee pool.
   Also: coordinate with Ref. [53] (single inclusive production at NLO, in preparation).
7. **Impact.** Substantially higher than Project A *if* it works — this would be a genuinely
   observable NLO fluctuation prediction. Also substantially more likely to stall.

---

## Project B — "Third factorial cumulant and the `O(g³)` evolution operator" (third paper, or never)

1. **Derivations.** Extend `Ω` to `O(g³)`: the three-gluon component of `Ω|0⟩`, plus `O(g³)`
   corrections to `A₂` and `B₂`. This is LCPT at the order Kovner–Lublinsky–Mulian worked at, but
   in normal-ordered Fock form. Then `C₃` complete at `O(g⁶)`.
2. **Numerics.** `C₃` in the translationally invariant model (6-D) and MV (lattice).
3. **Difficulty.** High — the operator extraction is the whole job.
4. **Time.** 12–24 months.
5. **Target.** JHEP (the operator extension is formalism, which JHEP takes) or PRD.
6. **Risks.** "This is KLM in a different basis" is a real and possibly fatal objection. Also: by
   the time this is done, `C₃` may be of little interest if `C₂` turned out small.
7. **Impact.** Formal. Would be cited by the same ~50 people.

---

## Recommended sequence and honest expectation

```
Phase 0 (2 weeks)  →  Project A+C+D  →  PRD, ~9 months from now
                          ↓
                   reassess: was C₂ big? was the away-side result new?
                          ↓
              yes → Project E (JHEP, +2 years)
              no  → Project B (formalism, PRD/JHEP, +18 months)
                    or move to a different problem
```

**One paper is near-certain. Two is plausible. Three is optimistic and I would not plan on it.**

The most valuable thing you can do in the next two weeks is Phase 0.6 — the go/no-go number. Every
strategic decision in this document hinges on whether the `O(g⁴)` quantum fluctuation is a
percent-level or a per-mille-level effect, and that is currently unknown.
