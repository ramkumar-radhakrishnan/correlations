# log(k⁺/Λ) (JIMWLK) and log(∨/k⁺) (BFKL) terms, and the finite terms, from Group_I/II/III.tex and 4rho.tex

The input is your `input/Group_I.tex`, `Group_II.tex`, `Group_III.tex` and `4rho.tex`. There are two outputs:
- [`LL_terms.pdf`](LL_terms.pdf) (source `LL_terms.tex`) has every term ∝ log(k⁺/Λ) and ∝ log(∨/k⁺), row by row in your form, then the three-ρ cancellation and the two-ρ sums.
- [`Finite_terms_rows.pdf`](Finite_terms_rows.pdf) (source `Finite_terms_rows.tex`) has everything else, row by row, with Λ → 0 wherever possible.

## Result

- **The split.** log(∨/Λ) = log(k⁺/Λ) + log(∨/k⁺). Every row is split exactly into log(k⁺/Λ)𝔗 + log(∨/k⁺)𝔓 + R.
  - We take ∨ > p⁺ and k⁺ > Λ, set Λ → 0 in R wherever possible, and keep ∨ exact.
  - Each of the 79 pointwise-evaluable terms equals its three parts to 1.5e-14 (`check_split.json`).
- **Three ρ's ∝ log(k⁺/Λ)** (22 terms, rows + `4rho.tex`): they cancel, with residual 7.7e-16 (`three_check.json`).
- **Three ρ's ∝ log(∨/k⁺)** (35 terms): they cancel, with residual 1.5e-16.
- **Four ρ's** (`4rho.tex`): they carry log(∨/Λ), so both logarithms get the same coefficient. They cancel, and their three- and two-ρ remainders enter both lists.
- **Two ρ's ∝ log(k⁺/Λ) = JIMWLK⊗LO**, and **two ρ's ∝ log(∨/k⁺) = BFKL⊗LO**, both exactly (`../final_check/render.py`).
- **Final:** d³N_LL = 𝒩[log(k⁺/Λ)·JIMWLK⊗LO + log(∨/k⁺)·BFKL⊗LO].

### Changes to your files
- **E1, Group I Row V 2b UV:** the coefficient must be 2 log(k⁺/Λ). The factor 2 is lost after "Setting Λ → 0".
  - With the 2, the poles × log(k⁺/Λ) of G1 II : G1 V 2b UV : G2 IV are 1 : −2 : 1, and they cancel (8e-16).
  - With your coefficient they leave 0.5 (`pole_check.json`).
- **F3, Group I Row III d³N₃ₐ:** replace the edge-kernel log line by 𝒟. Without it the three-ρ terms ∝ log(k⁺/Λ) leave 0.59.
- **H (hidden logarithms), Group III Rows III, V, VI and Group II Row IV:**
  - The log(∨/k⁺) in [log(k⁺/Λ) + log(∨/k⁺ − 1)] cancels against [1/p⁺]₊ (resp. ξ/(1−ξ)₊), so it is not a BFKL term.
  - The BFKL terms of these rows are (B) (region y′ → z) and H₄.
  - In the finite file they are subtracted (C3, C1). Sec. 2 of `Finite_terms_rows.pdf` (`blue_section.tex`) derives where both logarithms come from, with your [1/p⁺]₊f = (f(p⁺) − f(0))/p⁺.
- **C5, Group III Row III part 2:** the measure is ∫_{w,x′,y′,z,x}.
- **F1 and F2** (the `4rho.tex` two-ρ lists) and **M2** (the JIMWLK formula, App. F) are still needed.
- **Finite file:** the 1/Λ lines of G1 III cancel exactly, and so do your two red lines. The last line of d³N₃b (edge kernel) is not integrable (C4); its finite replacement, the 𝒟 region, is still open (Q3).

## Scripts

```
python3 check_split.py   # every term = J + B + exact remainder; remainder -> Lambda -> 0 form; 1/Lambda and red lines cancel
python3 cross_check.py   # (B - J) per three-rho row vs the exact row coefficient of ../finite_v2/split.py
python3 pole_check.py    # poles x log(k+/Lambda): cancel with E1
python3 three_check.py   # three-rho terms of both lists cancel (pointwise, symmetrized, with c.c.)
python3 make_ll.py       # writes LL_terms.tex and Finite_terms_rows.tex
```

- `blocks.py` cuts your final row equations into terms.
- `logsplit.py` is the table of logarithm splits, and `gen.py` does the term-by-term split.
- `texeval.py` evaluates your LaTeX numerically, and `assemble.py` rebuilds the align environments.
