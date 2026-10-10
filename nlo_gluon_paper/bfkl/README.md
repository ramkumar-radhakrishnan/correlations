# Terms ∝ log(∨/k⁺) and the BFKL evolution

[`LogVeeKplus_BFKL.pdf`](LogVeeKplus_BFKL.pdf) does for the terms proportional to log(∨/k⁺) what [`../three_rho`](../three_rho) and [`../two_rho`](../two_rho) did for log(∨/Λ):
- it collects the four-, three- and two-ρ terms row by row;
- it reorders the three-ρ terms and checks that their symmetrized parts cancel;
- it compares the two-ρ terms with BFKL;
- it writes out the BFKL terms explicitly.

## Result

- **What multiplies log(∨/k⁺).** With the rules log(Λ/k⁺) = −log(∨/Λ) + log(∨/k⁺) and log((∨−k⁺)/Λ) = log(∨/Λ) + log(1−k⁺/∨), the coefficient of log(∨/k⁺) is **BFKL⊗LO − JIMWLK⊗LO**.
  - ∨ appears only as the upper limit of the p⁺ integrals. So the coefficient of log∨ (log(∨/Λ) plus log(∨/k⁺)) is the large-p⁺ limit, which is BFKL⊗LO.
  - Equivalently: d³N_LL = log(k⁺/Λ)·JIMWLK⊗LO + log(∨/k⁺)·BFKL⊗LO.
- **Row by row.** For every row, the log(∨/k⁺) terms are B = P − A:
  - A is the row's log(∨/Λ) coefficient (as in the Evolution files);
  - P is its coefficient of log∨, the large-p⁺ limit of the integrand.
  - Many are explicit through log(Λ/k⁺), log(∨/(Λ+k⁺)), log(∨/k⁺−1).
  - Others are hidden in logarithms of transverse distances or in the plus-prescription integrals (Group II Row IV; Group III Rows III, V, VI).
- **Four-ρ:** there are no four-ρ terms ∝ log(∨/k⁺).
- **As written:** after symmetrization, 30 of 42 three-ρ structures are left over, and 35 of 49 two-ρ structures differ from BFKL⊗LO − JIMWLK⊗LO.
- **Groups II and III are correct** row by row in region P. Group III Rows III, V, VI must have the overall sign of `Evolution_three_rho.tex`.
- **Fixes:**
  - **P1 (new), Group I Row III:** the large-p⁺ limit of the 𝔹̄₃† term integrated over k⁺+Λ < p⁺ < ∨ must be +2K_h, where K_h = 𝒦(x′−w′)·𝒦(z−w) K(x,y,z). The notes give −2K_h. The p⁺ integration itself is correct (checked numerically); the sign error is in the integrand for p⁺ > k⁺.
  - **F3, F5, F6 carry over** from the log(∨/Λ) note, because the log(∨/k⁺) terms of a row contain minus its log(∨/Λ) terms:
    - F3: Row III part 2 → 𝒟;
    - F5: the sign of Group III Row III part 2;
    - F6: the Row II bracket.
  - **F4 does not change the log(∨/k⁺) terms.** The α part of Row II multiplies log(∨/Λ) only; the region-P comparison confirms F4 independently.
- **After the fixes:**
  - the symmetrized three-ρ terms cancel exactly (6 kernel structures, 68 colour networks);
  - the two-ρ terms equal BFKL⊗LO − JIMWLK⊗LO exactly. This holds pointwise in z and with exact coefficients, in all 18 kernel structures.
- **BFKL terms:**
  - Sec. 9 of the note gives BFKL⊗LO in operator form and fully written out.
  - In the colour-singlet projection it equals Eq. (LOBFKL) of Appendix E to machine precision, for symmetric and general U. **Appendix E is correct as written.**

## Reproducing

The scripts use [`../two_rho`](../two_rho) (term format, the notes' log(∨/Λ) terms, JIMWLK target) and [`../three_rho`](../three_rho) (engine, sources, reordering generator). The region-P reference rows and the BFKL target are in `reference/`.

```
python3 check_refP.py     # region-P reference rows = BFKL x LO (two-rho, pointwise)
python3 bfkl_singlet.py   # BFKL x LO (operator form) = Eq. (LOBFKL) in the singlet projection -> bfkl_singlet.json
python3 pk_checks.py      # row by row vs region P, cumulative fixes, seeds, general U -> checks.json (~5 min)
python3 states3.py        # three-rho states -> states3.json
python3 build_books.py 3  # symmetrized three-rho terms by kernel structure / colour network -> book3.pkl
python3 build_books.py 2  # two-rho terms + BFKL - JIMWLK -> book2.pkl
python3 reduceP.py book2.pkl reduced2.pkl   # exact coefficients, basis reduction
python3 reduceP.py book3.pkl reduced3.pkl
python3 make_pk_tex.py    # writes LogVeeKplus_BFKL.tex from pk_static.tex and rows.tex
pdflatex LogVeeKplus_BFKL.tex   # twice
```

| File | Role |
|---|---|
| `userP.py` | The coefficients of log∨ (P) of the rows of the notes: three-ρ sources and genuine two-ρ terms |
| `lsets.py` | The log(∨/k⁺) terms B = P − A, by row, as written and corrected |
| `refP.py`, `reference/rowsP_sym.py`, `reference/bsym.py` | Region-P reference rows and the BFKL target |
| `acc3.py`, `acc4.py` | Pointwise fingerprints for symmetrized three- and four-ρ terms |
| `book.py`, `build_books.py`, `reduceP.py` | Exact tables by kernel structure and colour network |
| `pk_checks.py`, `states3.py`, `check4.py`, `check_refP.py`, `bfkl_singlet.py` | The checks quoted in the note |
| `g1P_fit.py` | The fit that identifies the Group I fixes in region P (K_h sign, α part of Row II) |
| `row3_num.py` | Numerical p⁺ integration of Group I Row III (the integration is correct; the ∨ dependence equals the large-p⁺ limit) |
| `psec.py`, `make_pk_tex.py`, `pk_static.tex`, `rows.tex` | LaTeX generation |

Requires `numpy`, `sympy` and `scipy` (`row3_num.py` only).
