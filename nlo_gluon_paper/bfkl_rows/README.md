# BFKL terms and remaining terms from BFKL.tex and Other_terms.tex

The input is your `input/BFKL.tex` (the terms ∝ log(∨/k⁺), with three and four ρ's, and the three-ρ remainders of the four-ρ terms) and `input/Other_terms.tex` (the rest of the cross section after JIMWLK). The outputs:
- [`BFKL_rows.pdf`](BFKL_rows.pdf) (source `BFKL_rows.tex`) has the BFKL list row by row: your lines, the log(∨/k⁺) parts taken out of `Other_terms.tex` (blue), and the two lines that must go (red). It then shows the three-ρ cancellation, the two-ρ terms and the result, BFKL⊗LO.
- [`Remaining_rows.pdf`](Remaining_rows.pdf) (source `Remaining_rows.tex`) has your `Other_terms.tex` without the BFKL terms, row by row, with Λ → 0 (blue = changed or added, red = removed).
- [`Rest_terms.pdf`](Rest_terms.pdf) (source `Rest_terms.tex`) has the rest of the terms (everything apart from JIMWLK and BFKL), all black, in your notation. It is `Remaining_rows.pdf` without the colour marks, and the two subtractions of a hidden BFKL term are written as integrals: ∫₀^{∨−k⁺}dp⁺/(p⁺+k⁺) for Group III Rows III, V, VI and ∫_{k⁺/∨}^1 dξ/ξ for Group II Row IV.
- [`Why_log_vee.pdf`](Why_log_vee.pdf) (source `Why_log_vee.tex`) explains, step by step with your Group III Row VI, where the log(∨/k⁺) of Group III comes from.
  - The log(∨/k⁺ − 1) of line 2 is cancelled by the subtraction of [1/p⁺]₊ in line 1.
  - The BFKL term (B) sits in the −1/k⁺ term of line 1, where the large-p⁺ phase forces y′ → z and F(p⁺) ≈ k⁺F^z/(p⁺+k⁺).
  - A toy model with numbers is included (`why_toy.py`).

## Result

- **BFKL.tex + Other_terms.tex = row − JIMWLK**, row by row, to O(Λ). The two exceptions are the red lines below.
- **Two red lines in BFKL.tex.** These are Group I Row V 2b UV (from E1) and Group II Row IV line 1 (the hidden log(∨/k⁺) of that line, cancelled by the log(∨/k⁺ − 1) of line 2).
  - Their poles cancel each other (1.5e-15).
  - Their sum is not zero (0.23–0.51 relative). It is log(∨/k⁺)·(g⁴N_c/16π⁴) log[(w−w′)⁴/((x−w)²(x′−w′)²)]·𝒮, which is not part of BFKL⊗LO.
- **Hidden log(∨/k⁺) in Other_terms.tex.** Several logarithms in `Other_terms.tex` still contain log(∨/k⁺), e.g. log(∨/k⁺ − 1) and log[((∨−k⁺)A + k⁺B)/(ΛA + k⁺B)]. Their log(∨/k⁺) parts move to the BFKL list. This applies to:
  - Group I: Row III d³N₁ and d³N₂, Row IV.
  - Group II: Row I (1), (2b)–(2f).
  - Group III: Rows I, II.
- **(B), Group III Rows III, V, VI.** The BFKL term of these rows comes from the −1/k⁺ term (region y′ → z) and is added to the BFKL list. Line 2 keeps log(∨/k⁺ − 1), which cancels the subtraction in [1/p⁺]₊.
- **C1, Group II Row IV.** Your H₄ line in `BFKL.tex` is correct. `Other_terms.tex` must subtract −log(∨/k⁺)H₄ (the BFKL term inside the ξ integral of line 1).
- **Three ρ's** (35 terms with these changes): they cancel, with residual 3.5e-16 (`three_bfkl.json`). Your `BFKL.tex` alone (9 terms) leaves 0.39. The three-ρ parts of the (B) terms cancel among themselves, so the test passes without them too. (B) matters for the two-ρ terms: without it, 9 of 67 structures differ from BFKL⊗LO (`../bfkl_final/`).
- **Four ρ's:** your d³N₁ + d³N₂ equal the checked four-ρ two-ρ lists to 3.2e-15 (`compare_4rho.json`).
- **Two ρ's:** your commented "BFKL contribution" (explicit, collected and colour-singlet forms) is the sum, to 2.1e-16 (`compare_final.json`). It equals App. E, BFKL⊗LO.
- **Remaining terms = Finite_terms_rows.pdf** of `../llsplit/`, to 3.4e-16.

## Scripts

```
python3 extract.py        # split every term of Other_terms.tex into its log(vee/k+) part and the Lambda -> 0 remainder
python3 check_rows.py     # A: BFKL.tex + Other_terms.tex = row - JIMWLK; B: BFKL list = LL_terms; C: remaining = Finite_terms_rows
python3 three_bfkl.py     # three-rho terms of the BFKL list cancel (pointwise, symmetrized, with c.c.)
python3 compare_4rho.py   # your four-rho two-rho lists vs the checked ones
python3 compare_final.py  # your "BFKL contribution" vs the checked sum
python3 pole_bfkl.py      # the two red lines: poles cancel, finite sum does not
python3 make_bfkl.py      # writes BFKL_rows.tex and Remaining_rows.tex
python3 make_rest.py      # writes Rest_terms.tex (Remaining_rows.tex without colour marks; subtractions as integrals)
python3 why_toy.py        # toy model: p+ integrals of the three terms of Group III Row VI line 1 vs vee
python3 make_why.py       # writes Why_log_vee.tex
```
