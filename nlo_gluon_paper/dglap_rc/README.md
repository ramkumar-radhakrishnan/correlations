# DGLAP, running coupling and the finite terms (from Finite_terms_split, remaining terms)

[`Finite_terms_DGLAP_RC.pdf`](Finite_terms_DGLAP_RC.pdf) (source `Finite_terms_DGLAP_RC.tex`).

- **Running coupling (Sec. 2).**
  - The β₀ terms of G1 Row II part 1 (11/3) and G1 Row V 2b UV (−11/6) cancel exactly: pole, γ_E, μ and log(x−w)² all go.
  - β₀ enters through the coupling counterterm of LO. Its pole is the (β₀/2)δ(1−ξ) part of P_gg.
  - What remains is the renormalization-group term (α_s β₀/4π) log(μ_R²/μ_F²) LO, plus the constant (α_s/2π) N_c (67/9 − π²/3) LO.
- **DGLAP (Sec. 3).**
  - All poles add up to (α_s/2π)(−2/ε)[P_gg ⊗ LO − N_c log(1−k⁺/∨) LO]. The first term is absorbed in D_{h/g}(μ_F).
  - The DGLAP finite term is −(α_s/2π) ∫_{k⁺/∨}^1 dξ/ξ² P^{(r)}_{gg,+} LO_ξ[ℓ_Δ(μ_F)] − log(∨/k⁺) H₄.
  - The leftover pole (Q1) comes with its finite partner.
- **Finite terms (Secs. 4–7):** row by row, in the format of Finite_terms_1.tex.
- **Result (Sec. 8)** and **comparison with AABKL, arXiv:1411.2869 (Sec. 9).**
  - DGLAP is the same as in AABKL: splitting function, plus prescription, kinematic limit, logarithm, and scheme constant (their R = 2e^{Ψ(1)} = our c₀).
  - AABKL have no running-coupling term, because their LO is α_s⁰.
  - The finite terms match AABKL in colour structure and topology, not in kernels.

`checks.py` → `checks.json`: the sympy bookkeeping (β₀ cancellation, poles, finite parts) and the AABKL constants (−11/6, splitting function, R = c₀).

```
python3 checks.py; python3 make_dglap_rc.py; pdflatex Finite_terms_DGLAP_RC.tex   # three times
```
