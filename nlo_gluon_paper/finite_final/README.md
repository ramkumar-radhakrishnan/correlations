# Finite terms: check, DGLAP, running coupling, and the form of AABKL

> **Restored (exact identities):** the terms (B) (y′ → z) of Group III Rows III, V, VI, H₄ of Group II Row IV, their subtractions and the cancellation of the line-2 log(∨/k⁺) follow from exact identities (I2)–(I5) with ∨ exact; no large-p⁺ step is needed (`../bfkl_rows/exact_check.py`, `../full_check/Checks_JIMWLK_BFKL.pdf`, Sec. 3). The results here stand. The latest files are analysed in [`../full_check/`](../full_check/README.md).

The input is your `input/Finite_terms.tex`: every row of `Finite_terms_1.tex` minus its JIMWLK and BFKL terms. There are two outputs:
- [`Finite_terms_final.pdf`](Finite_terms_final.pdf) (source `Finite_terms_final.tex`, 20 pages) is the note.
- [`Finite_terms_corrected.pdf`](Finite_terms_corrected.pdf) (source `Finite_terms_corrected.tex`) is your file with the changes: red = remove, blue = add.

## Result

- **20 of 27 equations are correct as written.** Each is identical, term by term, to the remaining terms of `../bfkl_split/Finite_terms_split.tex` (`compare_file.json`).
- **Changes:**
  - **C1, Group II Row IV:** add the line −log(∨/k⁺)H₄. H₄ is in the BFKL list, and the ξ integral of line 1 also contains it.
  - **C2, Group III Rows III (both parts), V, VI:** [1/p⁺]₊ → [1/p⁺]₊^{(k⁺)}, still unevaluated. Otherwise −F(0)log(∨/k⁺) is counted twice.
  - **C3, the same rows:** subtract the terms (B) of `../bfkl_final` (region y′ → z), one line per row. With that line the row is finite for ∨ → ∞.
  - **C4, Group I Row III:** remove the two edge-kernel terms (first term of d³N₃ₐ, last of d³N₃b); the edge kernel is replaced by 𝒟 (fix (A)).
    - The last term of d³N₃b is not integrable at z → w (`edge_check.json`).
    - Its replacement, the finite part of the 𝒟 region, is open (Q3).
  - **C5, Group III Row III part 2:** the measure is ∫_{w,x′,y′,z,x}.
- **Running coupling and DGLAP** (`checks_rc.json`, all zero): the same as in `../dglap_rc`.
  - Running coupling: β₀ of G1 II part 1 + G1 V 2b UV = 0. What is left is (α_sβ₀/4π)log(μ_R²/μ_F²)LO + (α_s/2π)N_c(67/9 − π²/3)LO.
  - DGLAP: −(α_s/2π)∫ dξ/ξ² P_gg LO_ξ[ℓ_Δ(μ_F)] − log(∨/k⁺)H₄.
- **log(k⁺/Λ) and log(∨/k⁺).** Your two lists give log(∨/Λ)𝒥 + log(∨/k⁺)(ℬ − 𝒥), which equals log(k⁺/Λ)𝒥 + log(∨/k⁺)ℬ: JIMWLK (target) multiplies log(k⁺/Λ), BFKL (projectile) multiplies log(∨/k⁺).
- **AABKL form** (Sec. 11 of the note, `checks_xi.json`). Put ξ = p⁺/(p⁺+k⁺) and y′ = (1−ξ)u + ξz. Then Group III becomes AABKL's real term:
  - ∫₀^{1−k⁺/∨}dξ, with the parent at their (1−ξ)ȳ + ξz;
  - C2 together with your log(1−k⁺/∨) lines is their 1/[ξ]₊;
  - C3 is a plus prescription at ξ → 1, and what it subtracts is BFKL;
  - Group I Row V is their virtual ∫₀¹dξ;
  - Group II Row IV is their DGLAP term;
  - JIMWLK is their soft (∫dξ/ξ) term.

## Scripts

| File | Role |
|---|---|
| `compare_file.py` | your file vs the remaining terms, equation by equation → `compare_file.json` |
| `edge_check.py` | C4: Group I Row III near z → w and on the circle (x−z)² = (x−w)² → `edge_check.json` |
| `checks_rc.py` | running coupling, poles, DGLAP, AABKL constants (sympy) → `checks_rc.json` |
| `checks_xi.py` | C2, C3 and the ξ form: plus prescriptions, pointwise maps, ξ → 1 limit → `checks_xi.json` |
| `make_corrected.py` | `Finite_terms_corrected.tex` |
| `make_note.py`, `static.tex` | `Finite_terms_final.tex` (Secs. 4, 5 from `../dglap_rc`) |

```
python3 compare_file.py; python3 edge_check.py; python3 checks_rc.py; python3 checks_xi.py
python3 make_corrected.py; pdflatex Finite_terms_corrected.tex   # twice
python3 make_note.py; pdflatex Finite_terms_final.tex             # three times
```
