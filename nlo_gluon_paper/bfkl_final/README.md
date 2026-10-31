# BFKL terms of the projectile: three-ρ cancellation, two-ρ sum, colour-singlet form

> **Restored (exact identities):** the terms (B) (y′ → z) of Group III Rows III, V, VI, H₄ of Group II Row IV, their subtractions and the cancellation of the line-2 log(∨/k⁺) follow from exact identities (I2)–(I5) with ∨ exact; no large-p⁺ step is needed (`../bfkl_rows/exact_check.py`, `../full_check/Checks_JIMWLK_BFKL.pdf`, Sec. 3). The results here stand. The latest files are analysed in [`../full_check/`](../full_check/README.md).

The input is your file `input/bfkl.tex` (the log(∨/k⁺) list). There are two outputs:
- [`BFKL_final.pdf`](BFKL_final.pdf) (source `BFKL_final.tex`, 11 pages) is the note.
- [`bfkl_corrected.pdf`](bfkl_corrected.pdf) (source `bfkl_corrected.tex`) is your file with the changes: red = remove, blue = add.

## Result

- **Your file** is the list from `../bfkl_split/Finite_terms_split.pdf` (Secs. 2–4), term by term. Every row equals the log(∨/k⁺) coefficient of the same row of `Finite_terms_1.tex`, pointwise in the transverse coordinates.
- **As it stands, it does not give BFKL.**
  - The symmetrized three-ρ terms leave 10⁻³–4×10⁻² of their size.
  - 33 of 75 two-ρ kernel structures differ from BFKL⊗LO − JIMWLK⊗LO.
- **(A) Group I Row III.** The first term of d³N₃ₐ, with (x′−w′)·(y−z)/((x′−w′)²(y−z)²(z−w)²), is the edge-kernel remnant.
  - It is not integrable at z → w.
  - It is replaced by the 𝒟 term (fix F3), which your `Evolution_terms.tex` already uses.
- **(B) Group III Rows III (both parts), V, VI.** The p⁺ terms with the phase e^{−ik·(y′−z)p⁺/k⁺} that tend to a constant at large p⁺ give log(∨/k⁺), from y′ → z, after the y′ integral.
  - These terms are −δ_im δ_jk′/k⁺ in Rows III(1), V, VI and −(p⁺/k⁺²)[…] in Row III(2).
  - This is the same mechanism as H₄ in Group II Row IV.
  - The new terms are: the measured gluon of the conjugate amplitude is emitted by the soft gluon. They have the kernel 𝒦(y′−z)·𝒦(y−w) K(x′,x,z) and y′ → z in the colour.
- **With (A) and (B):**
  - The symmetrized three-ρ terms cancel.
    - Analytically: a table in Φ[p,q] = f^{abc} C_Lᵃ (U(p)ρ(x))ᵇ (U(q)ρ(y))ᶜ, with every row checked separately.
    - Numerically: to 10⁻¹⁷, with all c.c.'s and operator orders.
  - The two-ρ terms give log(∨/k⁺)(ℬ − 𝒥) in all 67 kernel structures (largest difference 3×10⁻¹⁵).
    - ℬ is your Sec. 4 (commented out; it equals App. E, LOBFKL).
    - 𝒥 is your Sec. 3.
  - Adding the log(∨/Λ)𝒥 of `Evolution_terms.tex` gives log(k⁺/Λ)𝒥 + log(∨/k⁺)ℬ.
  - The BFKL term is 𝒩 log(∨/k⁺) ∫ℬ. Its singlet form is your Sec. 4 singlet form, which equals EC (4.3) and LOBFKL.
- **Corrections to earlier notes** (Sec. 8 of the note):
  - `../bfkl_check`: the three-ρ and two-ρ claims hold only with (A) and (B).
  - `../bfkl_split`: Sec. 4 statement; the remaining terms of Group III Rows III, V, VI and Group I Row III.
  - `../finite_v2`: Sec. 5 statement on the e^{−iκ(1−a)/a} term.
  - `../dglap_rc`: the finite terms of the same rows. The DGLAP and running-coupling terms are unaffected.

## Scripts

| File | Role |
|---|---|
| `rows.py` | every row of your file (and A, B) in the reordering engine of `../two_rho/core.py`, your charge order, c.c. |
| `validate.py` | engine vs. the numerical transcriptions (`../bfkl_check/corrected_lvk.py`, `../finite_v2/split.py`) → `validate.json` |
| `three_rho_test.py`, `three_rho_engine.py` | symmetrized three-ρ terms (as written / with A) → `three_rho_engine.json` |
| `phi_table.py` | the Φ table, row by row → `phi_table.json` |
| `two_rho_engine.py` | two-ρ terms vs ℬ − 𝒥, key by key; targets checked against JIMWLK⊗LO and App. E → `two_rho_engine.json` |
| `rowcmp.py` | row by row vs `../bfkl` (region-P analysis) → `rowcmp.json` |
| `yz_test.py`, `table2.py` | with A and B: three-ρ and two-ρ, structure table → `yz_test.json`, `table2.json` |
| `checks.py` | identities: Row III simplification, F_z, (B) three-ρ, the Fourier lemma, the 2b counterterm → `checks.json` |
| `make_corrected.py`, `make_note.py`, `static.tex` | the two documents |

```
python3 validate.py; python3 three_rho_test.py; python3 three_rho_engine.py; python3 phi_table.py
python3 two_rho_engine.py; python3 rowcmp.py; python3 yz_test.py; python3 table2.py; python3 checks.py
python3 make_corrected.py; pdflatex bfkl_corrected.tex   # twice
python3 make_note.py; pdflatex BFKL_final.tex             # three times
```
