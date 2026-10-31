# BFKL terms and remaining terms from BFKL.tex and Other_terms.tex (∨ exact, exact identities only)

The inputs are your `input/BFKL.tex` and `input/Other_terms.tex`:
- `BFKL.tex` has the terms ∝ log(∨/k⁺), with three and four ρ's, and the three-ρ remainders of the four-ρ terms.
- `Other_terms.tex` has the rest of the cross section after JIMWLK.

The same analysis, made directly from your new `Group_I/II/III.tex` and `4rho.tex`, is in [`../full_check/`](../full_check/README.md).

**Rule.** The only constraints are ∨ > k⁺ and p⁺ > Λ. Λ → 0 wherever possible, and ∨ is exact. Every step is an exact identity; nothing is expanded or dropped.

## Outputs

- [`BFKL_rows.pdf`](BFKL_rows.pdf) (source `BFKL_rows.tex`) has the BFKL list row by row:
  - your lines (black);
  - the log(∨/k⁺) parts still in `Other_terms.tex`, and the terms that come out of line 1 by exact identities (blue);
  - the one line that is not in the rows (red).

  It also gives the identities, the three-ρ test, the two-ρ sum, the result and the poles.
- [`Remaining_rows.pdf`](Remaining_rows.pdf) (source `Remaining_rows.tex`) has `Other_terms.tex` without the BFKL terms, row by row, with Λ → 0 (blue = changed, red = removed).

## Result

- **BFKL is there.** d³N_BFKL = 𝒩 log(∨/k⁺) BFKL⊗LO, your App. E formula.
- **Where the BFKL terms are.**
  - In explicit logarithms of `Other_terms.tex`, e.g. log(∨/k⁺ − 1) = log(∨/k⁺) + log(1 − k⁺/∨).
  - In line 1 (the p⁺ or ξ integral) of Group III Rows III, V, VI and Group II Row IV. Exact identities bring them out:
    - **(I2)** [1/p⁺]₊ = [k⁺/(p⁺(p⁺+k⁺))]₊ + 1/(p⁺+k⁺) − f(0)/(p⁺+k⁺) under the integral, and ∫₀^{∨−k⁺} dp⁺/(p⁺+k⁺) = log(∨/k⁺). The −f(0) log(∨/k⁺) cancels the log(∨/k⁺) of line 2 exactly. For ξ: 1/(ξ(1−ξ)) = 1/ξ + 1/(1−ξ).
    - **(I3)** G(y′) = G(z) + [G(y′) − G(z)]. With ∫_u e^{−iq·u} u^m/u² = −2πi q^m/q² (q = k(p⁺+k⁺)/k⁺), the G(z) part is k⁺/(p⁺+k⁺) times the w′ form. Its p⁺ integral gives the terms (B) exactly.
    - **(I5)** the partial fractions of Row III part 2.
    - **(I4)** the Fourier transform of log(y′−y)² gives H₄, Group II Row IV.
  - The remaining terms then stay finite for ∨ → ∞, which fixes the coefficient of log(∨/k⁺) uniquely.
- **Checks** (`check_rows.json`, `exact_check.json`, `three_bfkl.json`, `pole_bfkl.json`):
  - BFKL.tex (without the red line) + Other_terms.tex = row − JIMWLK, row by row, to O(Λ) (≤ 2.5e-9 at Λ = 1e-9).
  - The rewriting of line 1 holds to ≤ 3.7e-15, with the p⁺ and ξ integrals done numerically. The G(z) term equals (B) to ≤ 3.4e-15. Your two Group II Row IV lines are exactly (I2) and (I4).
  - The BFKL list equals the list whose two-ρ sum is BFKL⊗LO (≤ 8.8e-17).
  - Its 41 three-ρ terms cancel (3.5e-16). Your BFKL.tex alone leaves 0.39.
  - Its poles cancel (2.9e-16).
  - Your commented "BFKL contribution" equals the sum (2.1e-16). Your four-ρ two-ρ lists equal the checked ones (3.2e-15).
- **Red:** only Group I Row V 2b UV (E1: Row V has no ∨).

**Correction of the previous version (commit 57f5646).** That version split only the explicit logarithms and kept the p⁺, ξ integrals untouched, so the log(∨/k⁺) inside line 1 stayed in the remaining terms, which then grow with ∨. That split was incomplete; this is why BFKL⊗LO did not come out. (B) and H₄ are not approximations: they follow from the exact identities (I2)–(I5). Its `Rest_terms` and `delta_two_rho` files are removed; they are in the git history.

## Scripts

```
python3 extract.py        # split every explicit log of Other_terms.tex into its log(vee/k+) part and the Lambda -> 0 remainder
python3 exact_check.py    # the identities (I2)-(I5) on line 1, p+ and xi integrals numerically (--fourier: the 2-D Fourier transforms)
python3 check_rows.py     # A: BFKL.tex + Other_terms.tex = row - JIMWLK; B: BFKL list = checked list; C: remaining = rest
python3 three_bfkl.py     # three-rho terms of the BFKL list (pointwise, symmetrized, with c.c.)
python3 pole_bfkl.py      # poles of the list
python3 compare_4rho.py   # your four-rho two-rho lists vs the checked ones
python3 compare_final.py  # your "BFKL contribution" vs the checked sum
python3 make_bfkl.py      # writes BFKL_rows.tex and Remaining_rows.tex
```
