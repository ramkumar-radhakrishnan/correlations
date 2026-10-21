# Final check of the evolution terms

[`Final_evolution_check.pdf`](Final_evolution_check.pdf) (32 pages). It checks the files of this upload, `input/`:
- `Evolution_terms.tex`
- `Evolution_three_rho.tex`
- `Group_I.tex`, `Group_II.tex`, `Group_III.tex`

It also uses `4rho.tex` from an earlier upload.

## Result

- **Four ρ's** cancel after symmetrization (4.7·10⁻¹⁵).
- **Three ρ's:** the list "All three ρ term with c.c." (152 terms) is identical to the v3 list and cancels exactly.
- **Two ρ's ∝ log(∨/Λ):**
  - Sec. 2 of `Evolution_terms.tex` (335 terms) agrees row by row with the checked lists: 0 differing kernel monomials in every row.
  - The genuine terms of Sec. 1 are correct. In Row V 2b the UV term and the counterterm cancel each other pointwise.
  - **M1:** the two-ρ lists of `4rho.tex` (with F1, F2) are missing. Without them, 18 of 152 independent colour structures differ (all charge–charge structures with the LO kernel K(x−w)·K(x′−w′)).
  - With them, the sum is **JIMWLK⊗LO exactly**: 13 kernel structures, 154 independent colour structures, 0 differences.
- **M2:** the JIMWLK formula as sent (App. F) has M at the wrong points in dN⁽¹⁾–dN⁽⁴⁾ and the wrong sign of dN⁽⁵⁾. The corrected formula is in Eq. (3.1).
- **BFKL** is not in the log(∨/Λ) terms. It is the coefficient of log ∨ (the log(∨/Λ) and log(∨/k⁺) terms together), and it equals BFKL⊗LO exactly: 9 kernel structures, 109 colour structures. App. E is correct as sent.
- **Final:** d³N_LL = 𝒩[log(k⁺/Λ)·JIMWLK⊗LO + log(∨/k⁺)·BFKL⊗LO].
- **DGLAP:**
  - The three pole terms ∝ log(∨/Λ) (Group I Row II, Group I Row V 2b UV, Group II Row IV) cancel in the pole, γ_E and μ. They leave πN_c log[(w−w′)⁴/((x−w)²(x′−w′)²)] × LO structure.
  - Group II Row IV with Group I and the charge renormalization gives P_gg and the DGLAP equation.
- **M3:** `Group_I.tex` (1.5) has a sign typo, (w−z)·(z−x) → (w−z)·(x−z); as written the 𝔉_β term vanishes. The f^{abc}UU term has the wrong overall sign.
- **M4:** the row-by-row display at the top of `Evolution_three_rho.tex` still lacks F4 (Row II) and 𝒟 (Row III). The list that is used has both.

Colours in the note:
- green: JIMWLK-only kernels (the unobserved gluon attached to the measured gluon);
- purple: charge–charge kernels (in both);
- orange: BFKL-only kernels (the measured gluon emitted by the gluon at z).

## Reproducing

```
python3 compare_sec2.py    # Sec. 2 of Evolution_terms.tex vs the checked lists -> compare_sec2.json
python3 genuine_new.py     # genuine terms of Sec. 1, encoded from the new text, vs the corrected ones
python3 total_check.py     # sum vs JIMWLK x LO, with and without the 4rho.tex lists -> total_check.json
python3 three_cancel.py    # the three-rho list: identical to v3? cancels? -> three_cancel.json
python3 targets_check.py   # App. F (as sent, corrected) and App. E vs the targets -> targets_check.json
python3 render.py          # exact tables (basis reduction), log(vee/Lambda) and log(vee) -> render.json, books.pkl
python3 make_note.py       # writes Final_evolution_check.tex from final_static.tex
pdflatex Final_evolution_check.tex   # three times
```

The scripts use `../two_rho`, `../three_rho`, `../two_rho_all`, `../three_rho_v2`, `../three_rho_v3` and `../bfkl`.
