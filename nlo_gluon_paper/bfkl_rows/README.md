# BFKL terms and remaining terms from BFKL.tex and Other_terms.tex (∨ exact, no large-p⁺ step)

The inputs are your `input/BFKL.tex` and `input/Other_terms.tex`:
- `BFKL.tex` has the terms ∝ log(∨/k⁺), with three and four ρ's, and the three-ρ remainders of the four-ρ terms.
- `Other_terms.tex` has the rest of the cross section after JIMWLK.

**Rule.** The only constraints are ∨ > k⁺ and p⁺ > Λ. Λ → 0 wherever possible, and ∨ is exact.
- Only explicit logarithms are split, and exactly, e.g. log(∨/k⁺ − 1) = log(∨/k⁺) + log(1 − k⁺/∨).
- Every p⁺ and ξ integral is kept as written: [1/p⁺]₊, ξ/(1−ξ)₊, upper limit ∨ − k⁺, lower limit k⁺/∨.
- There is no large-p⁺ (small-ξ) step, so these integrals give no BFKL term.

**Withdrawn from the previous version:**
- the terms (B) of Group III Rows III, V, VI;
- H₄ of Group II Row IV;
- the subtractions −(B) and −log(∨/k⁺)H₄;
- the "hidden logarithm" rule, which kept log(∨/k⁺ − 1) whole in line 2 of those rows;
- `Why_log_vee.pdf`. It is still in the git history, at commit 5710f0b.

All of these came from the large-p⁺ region.

## Outputs

- [`BFKL_rows.pdf`](BFKL_rows.pdf) (source `BFKL_rows.tex`) has the BFKL list row by row:
  - your lines;
  - the log(∨/k⁺) parts still in `Other_terms.tex` (blue);
  - the three lines that are not in the rows (red).

  It then gives the three-ρ test, what the two-ρ terms give, and the poles.
- [`Remaining_rows.pdf`](Remaining_rows.pdf) (source `Remaining_rows.tex`) has `Other_terms.tex` without the BFKL terms, row by row, with Λ → 0 (blue = changed, red = removed).
- [`Rest_terms.pdf`](Rest_terms.pdf) (source `Rest_terms.tex`) is the same list, all black, in your notation. It holds everything apart from JIMWLK and the explicit log(∨/k⁺).

## Result

- **BFKL.tex (without the red lines) + Other_terms.tex = row − JIMWLK**, row by row, to O(Λ) (≤ 2.5e-9 at Λ = 1e-9; `check_rows.json`).
- **Three red lines in BFKL.tex:**
  - Group I Row V 2b UV (E1: Row V has no ∨).
  - Both lines of Group II Row IV: −log(∨/k⁺)×(line 2), and log(∨/k⁺)H₄. They describe the ξ integral of line 1 for ∨ ≫ k⁺; they are not logarithms of the row.
  - With them, BFKL.tex + Other_terms.tex exceeds the row (0.66 and 1.0 relative).
- **Blue: every explicit log(∨/k⁺) still in Other_terms.tex:**
  - Group I Rows III, IV;
  - Group II Row I, and Row IV line 2;
  - Group III Rows I, II, and line 2 of Rows III (both parts), V, VI.

  Then the list is exactly the log(∨/k⁺) part of the rows (`check_rows.json`, B).
- **Group I Row II part 1:** log(1 − k⁺/∨) log(∨/k⁺) stays in the remaining terms. Its coefficient depends on ∨, and it is ≤ log²2 for all ∨ > k⁺ (maximum at ∨ = 2k⁺).
- **Three ρ's** (35 terms): they do not cancel: the residual is up to 1.01 times the largest term (three configurations). It is exactly the sum of the three line-2 terms of Group III Rows III(1), V, VI; all other terms cancel among themselves, to 3.6e-16 (`three_bfkl.json`).
  - Your BFKL.tex alone leaves 0.39.
  - The line-2 terms of Group III Rows III(1), V, VI alone leave 1.11.
- **Two ρ's: not BFKL⊗LO.** The previous list (with (B), H₄, without the line-2 terms) had two-ρ sum BFKL⊗LO. The new list differs from it by +(line-2 terms) − (B) − H₄. Its two-ρ terms are non-zero in all 41 kernel structures where they appear (32 from the line-2 terms, 9 from (B), H₄; `delta_two_rho.json`).
  - Every engine transcription used there equals its LaTeX term to ≤ 5e-14.
- **Poles.** The blue line-2 term of Group II Row IV is the only term of the list with 1/ε. Its pole equals that of the red Group I Row V line (ratio 1 to 2e-15), so nothing in the list cancels it (`pole_bfkl.json`).
- **Not changed:**
  - your four-ρ two-ρ lists equal the checked ones (3.2e-15, `compare_4rho.json`);
  - your commented "BFKL contribution" equals the checked sum of the previous list, i.e. BFKL⊗LO (2.1e-16, `compare_final.json`).

So, with ∨ exact and the integrals kept, the coefficient of the explicit log(∨/k⁺) is not BFKL⊗LO. The rest of the ∨ dependence of the rows is in the remaining terms:
- the p⁺ integrals up to ∨ − k⁺;
- the ξ integrals from k⁺/∨;
- log(1 − k⁺/∨);
- Rⁱ.

**A second exact bookkeeping (your decision).** Split the subtraction of [1/p⁺]₊ with 1/p⁺ = 1/(p⁺+k⁺) + k⁺/(p⁺(p⁺+k⁺)), or for ξ/(1−ξ)₊ with 1/(ξ(1−ξ)) = 1/ξ + 1/(1−ξ).
- This gives −f(0) log(∨/k⁺) exactly, which cancels the line-2 log(∨/k⁺). Group III Rows III, V, VI and Group II Row IV then have no BFKL term.
- Your red line 1 of Group II Row IV is exactly this −g(1) log(∨/k⁺).
- Then the three-ρ terms cancel (3.6e-16), and there is no pole in the list.
- The two-ρ terms still differ from BFKL⊗LO, in the 9 structures of −(B) − H₄.

In neither bookkeeping is the two-ρ sum BFKL⊗LO without (B) and H₄.

## Scripts

```
python3 extract.py        # split every explicit log of Other_terms.tex into its log(vee/k+) part and the Lambda -> 0 remainder
python3 check_rows.py     # A: BFKL.tex + Other_terms.tex = row - JIMWLK; B: BFKL list = log(vee/k+) part; C: remaining = rest
python3 three_bfkl.py     # three-rho terms of the BFKL list (pointwise, symmetrized, with c.c.)
python3 delta_two_rho.py  # two-rho terms of the change: + line-2 terms - (B) - H4 (engine of ../two_rho, validated vs the LaTeX)
python3 pole_bfkl.py      # red lines: poles cancel, sum does not; the pole of the list
python3 compare_4rho.py   # your four-rho two-rho lists vs the checked ones
python3 compare_final.py  # your "BFKL contribution" vs the checked sum of the previous list
python3 make_bfkl.py      # writes BFKL_rows.tex and Remaining_rows.tex
python3 make_rest.py      # writes Rest_terms.tex (Remaining_rows.tex without colour marks)
```
