# Your Group I, II, III and 4rho.tex: checks, JIMWLK, BFKL, remaining terms, running coupling, DGLAP, finite terms

The inputs are your latest uploads, in `input/`: `Group_I.tex`, `Group_II.tex`, `Group_III.tex` and `4rho.tex`.
- Groups I and III are taken with their complex conjugate, and Group II without, in every test.
- Only ∨ > k⁺ and p⁺ > Λ are used. Λ → 0 wherever possible, ∨ is exact, and every step is an exact identity.

## Outputs

- [`Checks_JIMWLK_BFKL.pdf`](Checks_JIMWLK_BFKL.pdf) (source `Checks_JIMWLK_BFKL.tex`) contains:
  - every step of your files, checked;
  - the exact identities;
  - every term ∝ log(k⁺/Λ) (JIMWLK) and ∝ log(∨/k⁺) (BFKL), row by row;
  - the four- and three-ρ cancellation;
  - the two-ρ sums, compared with your JIMWLK and BFKL formulas.
- [`Remaining_RC_DGLAP_finite.pdf`](Remaining_RC_DGLAP_finite.pdf) (source `Remaining_RC_DGLAP_finite.tex`) contains:
  - the remaining terms, row by row;
  - the running coupling;
  - DGLAP;
  - the finite terms, step by step;
  - the result.

- [`Comparison_AABKL.pdf`](Comparison_AABKL.pdf) (source `Comparison_AABKL.tex`): the final result against AABKL, arXiv:1411.2869 (forward rapidity): what agrees, the forward limit, what differs at mid rapidity.
- [`Steps_line1_and_mu.pdf`](Steps_line1_and_mu.pdf) (source `Steps_line1_and_mu.tex`) contains:
  - starting from your own Eqs. (1.12)–(1.14), step by step, where the extra log(∨/k⁺) terms of Group III Rows III, V, VI and Group II Row IV come from: (i), (ii) and (iii);
  - the coefficient in front of log μ² from the five sources (A3 and the shockwave through the loop, each with its c.c., and the real term), plus the coupling of LO;
  - why β₀ is not in the NLO terms, and why Group II Row IV is DGLAP and not running coupling.

## Result

- **Your files, step by step** (`step_Group_*.json`, `four_file.json`):
  - Every step passes: 1e-16 to 6e-12, and O(Λ) for the plus prescriptions. E1 is now correct in your Group I.
  - The exceptions are the typo M2 (Group II Eq. 1.16, 0.35 as written, 2e-16 corrected) and the measures and typos M3–M6.
  - `4rho.tex`: every equation equals the checked earlier file, except one. The first two-ρ part, Eq. (1.12), has lost its factor [log(k⁺/Λ) + log(∨/k⁺)]. Its value is exactly 1/log(∨/Λ) times the earlier one.
  - The two-ρ lists still need F1 and F2. With the log restored and F1, F2, they equal the checked lists to 6e-15.
- **The split.** Every row = log(k⁺/Λ)𝔗 + log(∨/k⁺)𝔓 + R, exactly. Each of the 79 pointwise-evaluable terms equals its parts to 1.5e-14 (`ll_split.json`).
  - Inside line 1 of Group III Rows III, V, VI and Group II Row IV, the BFKL terms come out by exact identities. (I2) is the subtraction of the plus prescription, which cancels the line-2 log(∨/k⁺). (I3) is the y′ = z Fourier identity, which gives (B). (I4) is the Fourier transform of log(y′−y)², which gives H₄.
  - These identities are checked in `../bfkl_rows/exact_check.py` (≤ 3.7e-15).
- **Four ρ's cancel** after symmetrization (1.3e-16).
- **Three ρ's ∝ log(k⁺/Λ)** (22 terms) **cancel** (7.7e-16), with 𝒟 (F3). With your edge kernel they leave 0.59.
- **Three ρ's ∝ log(∨/k⁺)** (35 terms) **cancel** (1.5e-16) (`ll_three.json`).
- **Two ρ's ∝ log(k⁺/Λ) = JIMWLK⊗LO** (0 of 154 colour structures differ).
  - Your JIMWLK formula, as sent, differs by 0.22–0.47.
  - Corrected (M at the right points, + sign of dN⁽⁵⁾), it agrees to 2e-15.
- **Two ρ's ∝ log(∨/k⁺) = BFKL⊗LO** (0 of 109 differ). It is your BFKL formula, as written (8e-15).
- **d³N_LL = 𝒩[log(k⁺/Λ) JIMWLK⊗LO + log(∨/k⁺) BFKL⊗LO].** In terms of log(Λ/k⁺), the JIMWLK coefficient is −JIMWLK⊗LO.
- **Poles × log(k⁺/Λ)** of Group I Row II, Group I Row V 2b UV and Group II Row IV are in the ratio 1 : −2 : 1 and cancel (8e-16). The BFKL terms are pole free.
- **The new rows equal the earlier ones.** Their J/B/remaining parts equal the earlier files row by row (0, `same_lists.json`). The earlier BFKL.tex + Other_terms.tex equal the new rows minus JIMWLK (O(Λ), `rows_vs_files.json`).
- **Running coupling.**
  - β₀ (G1 II part 1) + β₀ (G1 V 2b UV) = 0.
  - What is left is the RG term (α_sβ₀/4π) log(μ_R²/μ_F²) LO and (α_s/2π) N_c(67/9 − π²/3) LO.
- **DGLAP.** −(α_s/2π) ∫_{k⁺/∨}¹ dξ/ξ² P_gg LO_ξ[ℓ_Δ(μ_F)] − log(∨/k⁺) H₄.
- **Finite terms.** Every remaining row, without RC and DGLAP. From Group I Row II part 1 only V_∨ = −(α_s/2π) N_c log(1−k⁺/∨) log(∨/k⁺−1) LO is left.
- **Open:**
  - Q1 is the pole ∝ log(1−k⁺/∨).
  - Q2 is a finite four-ρ term, which is not in your files.
  - Q3 is the finite part of the 𝒟 region.

## Scripts

```
python3 step_check.py Group_I     # every step of your file (same for Group_II, Group_III) -> step_<group>.json
python3 four_file_check.py        # 4rho.tex vs the checked file; four-rho cancellation; two-rho lists -> four_file.json
python3 ll_checks.py split poles three   # the split, poles x log(k+/Lambda), three-rho cancellation -> ll_*.json
python3 same_lists.py             # J/B/remaining parts vs the earlier files -> same_lists.json
python3 rows_vs_files.py          # earlier BFKL.tex + Other_terms.tex vs the new rows -> rows_vs_files.json
python3 make_ll.py                # writes Checks_JIMWLK_BFKL.tex
python3 make_rest.py              # writes Remaining_RC_DGLAP_finite.tex
python3 mu_check.py               # coefficient of log mu^2 of every source -> mu_check.json
python3 make_steps.py             # writes Steps_line1_and_mu.tex
```
`ll.py` points the split machinery of `../llsplit` at your new files.

## Your Group files, corrected (`groups/`)

`python3 groups/make_groups.py` writes the `_marked` copies, and `python3 groups/make_groups.py clean` writes the clean ones. Each version is in `.tex` and `.pdf`.

- **`Group_I`, `Group_II`, `Group_III`** are your current files with every remaining correction applied. Your own colours are kept and there are no correction marks.
- **`Group_*_marked`** shows the same corrections in place: red = remove, blue = new. Your own colours are remapped: blue → teal, green → olive, red → violet, orange → brown.
- **Group I:**
  - Eq. (1.6) has the measure ∫_{w,w′,x′}.
  - A stray comma is removed in Eqs. (1.1), (1.2).
  - The commas in ∫_{x,z,y,x′,y′,w′} are fixed.
  - Braces: Eq. (1.9) is closed, and the extra \Bigg\} is removed in Eqs. (1.19), (1.20), (1.22), (1.24).
  - E1 was already correct.
- **Group II:**
  - The phase in Eq. (1.16) (M2) is corrected.
  - k{+} → k⁺ in Eq. (1.13).
  - The repeated sentence is removed.
- **Group III:**
  - The first line of Row III has sign − (M3).
  - The measures of Row III part 2 are corrected (C5).
