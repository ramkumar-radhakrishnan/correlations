# The three-ρ terms ∝ log(∨/Λ): check of the list, fixes, two-ρ terms

[`ThreeRho_check.pdf`](ThreeRho_check.pdf) (25 pages) checks the list of all three-ρ terms ∝ log(∨/Λ) with c.c. The list is [`input/three_rho_all.tex`](input/three_rho_all.tex), verbatim as sent: 160 terms in common variables, kernel classes α, α′, β, β′, γ, γ′.

## Result

- **As sent, the terms do not cancel.**
  - After symmetrization, all 36 kernel monomials are left over (pointwise).
  - With exact coefficients, 24 of 50 colour networks do not cancel.
- **Five changes make them cancel exactly.** After them, all 40 networks cancel and the pointwise residual is 0.
  - **T1** Group III Row I, term 15 (β′): ρ^h(x′) → ρ^c(x′). As written, the index h occurs three times.
  - **T2** Group III Row II, term 13 (β′): U^{hl}(y) → U^{hl}(x).
  - **T3** Group I Row II, part 1, the α terms: kernel 𝒦(x′−w)·𝒦(y−z) → 𝒦(x−z)·𝒦(y−z) (the α kernel).
  - **F4** Group I Row II: remove the α part of the U(w)𝔸⁽³⁾ sub-terms: part 1, terms 1, 2, 5, 6, and their c.c. This is the commutator [ℕ⁽²⁾, 𝔸⁽¹⁾]; with the p⁺ ordering it multiplies −log(∨/k⁺).
  - **F3** Group I Row III: replace the γ, γ′ terms (edge kernel φ) by 𝒟. The 8 terms are in [`input/D_terms.tex`](input/D_terms.tex), in the same notation.
  - Apart from T1–T3, the list equals the earlier transcription as written. F3 and F4 are the fixes of the earlier three-ρ note.
- **The two-ρ terms need the operator order of the charges.** The list gives the order of the symmetrized check, which differs in two ways:
  - every c.c. term must have its charges reversed (ABC → CBA flips the sign of the two-ρ terms);
  - the anticommutator rows (Group I Row V, Group II Row I, Group III Row III) must be ½·A{B,C}, not 2·ABC.

  Reordering the list literally gives wrong two-ρ terms in 11 of 12 sources. The note gives the list again in operator order, with the user's indices. Its two-ρ terms equal those of the original expressions, row by row.
- **Two-ρ terms:** they are listed per source. Together with the genuine two-ρ terms and the `4rho.tex` lists, they give JIMWLK⊗LO exactly. BFKL⊗LO is the coefficient of log∨ (see [`../two_rho_all`](../two_rho_all)).

## Reproducing

```
python3 parse.py          # parses the list; reports kernels that differ from their class label
python3 check3.py         # row-by-row comparison with the earlier transcription (symmetrized)
python3 cmp_entries.py    # term-by-term comparison for Group III Rows I, II (finds T1, T2)
python3 check3b.py        # cancellation as sent and after each fix -> states3.json
python3 networks.py       # exact colour-network tables -> networks.pkl
python3 two_cmp.py        # two-rho terms: literal order vs operator order -> two_cmp.json
python3 oporder.py        # operator order of every term; check against the original expressions
python3 make_note.py      # writes ThreeRho_check.tex from note_static.tex
pdflatex ThreeRho_check.tex   # three times
```

The scripts use the modules of [`../two_rho`](../two_rho), [`../three_rho`](../three_rho), [`../bfkl`](../bfkl), [`../combined`](../combined) and [`../two_rho_all`](../two_rho_all). They require `numpy` and `sympy`.
