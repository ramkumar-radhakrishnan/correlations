# Prompt: independent check of the final NLO single inclusive gluon cross section

*Copy everything below the line into a new conversation (or give it to a human referee), together with the files listed in "Material provided".*

---

You are an expert in small-x QCD, the Colour Glass Condensate, light-cone perturbation theory (LCPT) and NLO calculations. Your task is to **check, independently and critically, the final result for single inclusive gluon production at mid rapidity at next-to-leading order in the CGC**, and to write a referee-style review of it. Do not trust any statement in the material because it is written there; re-derive or re-compute what you check. Do not assume anything that is not determined by the material: when something is ambiguous, say so and state what you would need.

## 1. Physics set-up (as used in the material)

- Dilute projectile described by colour-charge density operators ρ^a(x) that do **not** commute, [ρ^a(x), ρ^b(y)] = i f^{abc} ρ^c(x) δ²(x−y); dense target described by adjoint eikonal Wilson lines U^{ab}(x) (real, orthogonal). Pure Yang–Mills (N_f = 0).
- Projectile gluon modes split by longitudinal momentum: soft modes Λ < p⁺ < ∨ (∨ is the upper cut-off, Λ the lower); modes above ∨ are valence, below Λ belong to the target.
- The projectile state is dressed by an LCPT dressing operator Ω (coefficients 𝔑, 𝔄₁, 𝔄₂ incl. the order-g³ coefficient 𝔄₂⁽³⁾, 𝔅₁, 𝔅₂, 𝔅₃, ℭ₁, ℭ₂). The observable is
  dN/d³k = (2π)⁻³ ∫ e^{−ik·(w−w′)} ⟨ Ω†S†Ω a†(w′) a(w) Ω†SΩ ⟩,
  with phase convention e^{−ik·(unprimed − primed)} everywhere.
- At NLO the result is organised in three groups (Group I: 2Re[LO*×NLO] one-gluon final state, with c.c.; Group II: diagonal two-gluon terms, no c.c.; Group III: interference terms, with c.c.), each split into rows.
- Final result = LO + running-coupling (RC) term + BFKL term (∝ log(∨/k⁺), projectile evolution) + JIMWLK term (∝ log(k⁺/Λ), target evolution) + DGLAP term (collinear, MS-bar number-density fragmentation function, CDR, counterterm with d-dimensional Jacobian ξ^{−d}) + finite terms (limits Λ → 0 and ∨ → ∞, with ∨ → ∞ taken before ε → 0).
- Finite terms are given (a) with ordered charges, (b) after the exact cancellations, with every three-ρ product split exactly into the symmetrised product plus commutator (two-ρ contact) terms plus double-commutator (one-ρ) terms, (c) averaged over the projectile with ⟨ρ^a(x)ρ^b(y)⟩ = δ^{ab} n(x,y)/(N_c²−1) and ⟨{ρ^aρ^bρ^c}_sym⟩ = d^{abc} N_c/((N_c²−4)(N_c²−1)) σ(x,y,z), σ ≡ ⟨d^{efg}{ρ^eρ^fρ^g}_sym⟩ (no f-type part, no model for σ, σ-terms absent for SU(2)), in colour-singlet form with traces of Wilson lines, and (d) with colour structures factorised as products of Wilson-line differences [U−U]×[U−U]. Longitudinal integrals are written with plus prescriptions, ∫dp⁺[1/p⁺]₊ f(p⁺) = ∫dp⁺ (f(p⁺)−f(0))/p⁺, and with ξ = p⁺/k⁺, except in Group II Row IV and the DGLAP term, where ξ = k⁺/(p⁺+k⁺) is the observed gluon's momentum fraction of its parent. The energy denominators of 𝔅₃ (Group I Row III) are principal values.

## 2. Material provided

1. **Final result:** `Final_terms_explicit_2_fixed.tex` (+ PDF) — all sections (BFKL, JIMWLK, DGLAP, RC, finite terms in forms (a)–(d)).
2. **Derivation input:** `Calculation_main_1_checked.tex` (wave functions, dressing operator Ω and its coefficients, the general expression for the spectrum), `Group_I_checked.tex`, `Group_II_checked.tex`, `Group_II_Row_IV_checked.tex`, `Group_III_checked.tex`, `4rho_checked.tex`, the evolution files (`BFKL_checked.tex`, `Evolution_JIMWLK_terms_checked.tex`, `Evolution_two_rho_checked.tex`, `Evolution_three_rho_checked.tex`), `Divergent_terms_fixed.tex`, `Limiting_cases_final_fixed.tex`.
3. **Paper draft** (`main.pdf`) and the notes on the checks already made (`fix_*.pdf`, `ref3_*.pdf`) — use these only to know what has been claimed; check the claims yourself.
4. Optionally: the referee reports already received.

## 3. What to check (do each item; mark each as VERIFIED (you re-derived/re-computed it), PLAUSIBLE (consistent but not re-derived), WRONG (with the correct form), or OPEN)

### A. Inputs of the derivation
1. LO: light-cone wave function, normalisation 𝔑₂, amplitude 𝔄⁽¹⁾, LO spectrum and its prefactor α_s/(4π⁴k⁺); its k_T-factorised limit; single-quark limit (Gunion–Bertsch α_sC_F/(π²k²)).
2. The order-g³ wave-function coefficient 𝔄₂⁽³⁾: operator ordering of its three-ρ part (𝔄⁽¹⁾𝔑₂, not 𝔑₂𝔄), the one-loop constant 67/18 and the cancellation of the cut-off double logs; Fourier transform to coordinate space.
3. 𝔅₂, 𝔅₃, ℭ₁, ℭ₂: matrix elements, supports in p⁺ (kinematics, not coefficients, decide the limits), signs; unitarity relations of Ω.
4. Assembly of each Group/Row from Ω†S†Ω a†a Ω†SΩ: which terms come with c.c., charge ordering (conjugate amplitude on the left), and the four-ρ terms (must contribute only to the evolution terms).

### B. Logarithmic terms
5. BFKL term: coefficient of log(∨/k⁺) equals the one-step projectile (KLWMIJ/BFKL) evolution of LO, for non-commuting ρ, including normalisation; colour-singlet form.
6. JIMWLK term: coefficient of log(k⁺/Λ) equals −H_JIMWLK acting on LO; colour-singlet form.
7. Three-ρ parts of the evolution terms: their symmetrised (incl. d^{abc}) parts cancel; reason (C-parity plus factorised form).

### C. Divergences, DGLAP, running coupling
8. Pole bookkeeping: every 1/ε pole, UV or IR, and its cancellation (Row I.II vs Row I.V 11/3 terms; collinear pole vs fragmentation-function counterterm; coupling counterterm); the role of the order of limits ∨ → ∞ before ε → 0.
9. DGLAP term: splitting function P̃_gg incl. (β₀/2)δ(1−ξ), β₀ = 11N_c/3; ξ → 0 subtraction; scheme (MS-bar number density, CDR) and the shift Δ to the Chirilli–Xiao–Yuan convention; reproduce at least one standard NLO coefficient function (e⁺e⁻ → gX or H → gX) with the same conventions.
10. Running coupling: LO + RC = LO with g²(μ); separation of μ_R and μ_F; which "11/6" numbers are β₀/(2N_c) and which are splitting-function integrals.

### D. Finite terms (the main part of the review)
11. Equation by equation, check that each finite term of each Group/Row follows from the derivation input in the limits Λ → 0, ∨ → ∞ (prefactors, signs, factors of 2π, plus prescriptions, integration limits, phases incl. recoil phases e^{−ik·(…)p⁺/k⁺}, colour structures, charge ordering, c.c.).
12. Exact cancellations: Group I Row IV vs Group III Row I; the 11/3 terms of Rows I.II and I.V. Verify they are exact.
13. Tail-only cancellations: the soft remainders of the plus prescriptions (pairs between Group I Rows III, V and Group III Rows III, V, VI), the "problematic" third term of N₃b in Group I Row III (needs the w integral before z), and the large-distance tail of Group II Row IV (1:1:−2 with Group II Row I second part and Group III Row III). Verify that the terms are integrable only in these combinations and that nothing has been dropped that only cancels in the tail.
14. Integrability (finiteness) of every remaining term: power-count every coincidence of integration points (UV/short distance) and the large-distance behaviour, including the vanishing of colour factors at coincidences; the Group I Row II N₂b kernel in momentum space (p → q, p → 0: only the antisymmetrised kernel ½[R(q,p)/p² − R(q,q−p)/(q−p)²] is integrable — check that this is what is used); the conditionally convergent collinear point of Group I Row V. State whether the result is absolutely convergent once grouped as prescribed, and whether integration-order or angular prescriptions can shift it.
15. No residual dependence on μ, Λ or ∨ in the finite terms.
16. Symmetrisation identities (ordered, ρ{ρ,ρ}, {ρ,ρ}ρ) and their application row by row; the contact terms; the one-ρ terms.
17. Projectile average and colour-singlet forms: re-derive at least three rows completely (one from each group); check the d^{abc} normalisation; check every trace identity used; check the [U−U]×[U−U] factorised forms are exact identities.
18. Parity in k: which parts are even/odd; the claim that the symmetrised cubic part is odd in k configuration by configuration and C-odd in projectile and target.
19. ξ variables: check every rewriting from p⁺ to ξ (Jacobians, plus prescriptions, limits), in particular Group II Row IV, where ξ = k⁺/(p⁺+k⁺) as in the DGLAP term (Jacobian dp⁺ = k⁺dξ/ξ²).

### E. Limits and comparison with known results
20. Dilute–dilute limit: real emission vs the Bartels–Sabio Vera–Schwennsen quasi-multi-Regge gg matrix element; virtual part vs the one-loop Lipatov vertex and the Regge trajectories (Fadin–Lipatov, Del Duca–Schmidt); the 67/9 − π²/3 constant (it should not appear as a net LO-type constant of Group I).
21. Single quark / colour-singlet dipole projectile computed exactly in colour space: uniqueness of the split into symmetrised and contact parts.
22. Comparison with the classical first saturation correction (McLerran–Skokov; Kovchegov–Skokov; Li–Skokov), including the known open factor-i phase difference in the nonlinear Weizsäcker–Williams structures.

## 4. How to work
- Re-derive by hand/computer algebra where possible; check numerically where not (random SU(2)/SU(3) Wilson lines, random points, explicit charge-space matrices for non-commuting ρ, quadrature for longitudinal integrals). Quote the numbers you obtain.
- For every discrepancy: give the equation, what is written, what you obtain, how you obtained it, and its effect on the result (sign, magnitude, which parity sector).
- Do not "fix" anything by assumption; if a convention is not stated, list it as an ambiguity.

## 5. Output
1. **Verdict** (one paragraph): is the final cross section correct, correct up to listed corrections, or not established.
2. **Table** of all items 1–22 with status VERIFIED / PLAUSIBLE / WRONG / OPEN and one line each.
3. **Errors** found, each with the corrected expression.
4. **Major concerns** (physics, scheme, prescriptions, convergence) and **minor points** (typos, notation, presentation).
5. **What must be done** before the result can be published, ordered by importance.
