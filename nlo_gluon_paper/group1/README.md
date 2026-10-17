# Group I cross section (one gluon in the final state): check against `Calculation_main.tex`

[`Group_I_check.pdf`](Group_I_check.pdf) (8 pages) lists what is correct, the mistakes and their fixes. [`Group_I_corrected.pdf`](Group_I_corrected.pdf) is your file with every correction marked in place: changes in blue, removed text in red (source `Group_I_corrected.tex`). The file as received is in `input/`.

## How it was checked

As for Groups II and III. Every displayed equation of `Group_I.tex` was transcribed into code (`user1.py`) and compared with a first-principles evaluation (`fp1.py`). The reference is Eq. (CS1) of `Calculation_main.tex` with its coefficients 𝔸⁽¹⁾, 𝔸⁽³⁾, 𝔹₃, 𝔹₂ and ℭ₁. The comparison uses:
- random transverse points, random charges and random SU(3) adjoint Wilson lines (symmetric, and general);
- exact (principal-value) p⁺ integrals, with no limit in ∨, k⁺ or Λ.

The instantaneous and two-ρ entries of 𝔹₃ were rebuilt from the matrix elements m1 and m5 (`wf3.py`). They equal (2π)³ times the light-cone wave functions.

## Result

- **Correct as written** (also for a general U):
  - Row II: part 1 and the rewriting; the prefactor of part 2; the split into 2a and 2b, with 𝔉 and Rⁱ.
  - Row III: the explicit form (with 𝔹₃ as written) and every p⁺ integral (N1, N2, N3a, N3b, including the arctan 𝒜), to 10⁻¹⁴.
  - Row IV: everything, exactly.
  - Row V: part 1 (new) up to the plus-prescription limits; part 2 as checked before.
- **G1. Row III:** the p⁺ > k⁺ non-instantaneous line has the wrong sign. This is correction E3 of 𝔹₃.
  - So every line of N2 changes sign.
  - With the fix, all green terms (1/Λ, 1/(∨−k⁺), 1/(k⁺−Λ)) cancel exactly. As written, 2·[1/Λ − 1/(∨−k⁺)] is left over.
- **G2. Row II, part 2:**
  - The bracket of (1.3) lacks −kⁱ log(∨/k⁺). The split 2a + 2b already assumes it.
  - In (1.3) and (1.8), the colour factor should be U^{ah}(w) f^{hbc}{ρ,ρ} − f^{abc} U U {ρ,ρ}.
- **G3. Row II:** the commutator [ℕ₂, 𝔸⁽¹⁾] is missing (correction E5, fix F4). In 2a, the U^{ah}(w) term needs log(∨/Λ) 𝔉 → log(∨/Λ) 𝔉_β.
- **G4. Row III:** the p⁺ integrand has a pole at p⁺ = k⁺(x−w)²/(x−z)². The logarithms must be read as log|·| (principal value).
- **G5. Exact forms (no limits):**
  - Row II brackets: exact forms X6, X4 and X7, and the exact Rⁱ.
  - Row V plus prescription (part 1 and 2a): keep Λ < p⁺ < k⁺ − Λ. Off by O(Λ/k⁺) otherwise.
  - (1.28): keep the Λ terms of (1.27).
- **G6. Notation:**
  - log[1 − k⁺/∨]² is log².
  - Add ∫_{x′}.
  - S^{ad} → U^{ad}.
  - Put the phase inside ∫dp⁺.
  - Add a missing comma.
- **The 1/ε bracket of Row II** equals (1/π)∫_z[K(x,w,z) − K(w,w,z)] exactly in d⊥ = 2 − ε, MS-bar. This is the scheme of Row V part 2b. Fix F6 of the two-ρ note concerned a different convention (d = 2 − 2ε).
- **Four-ρ parts of Rows III and IV** are not in `Group_I.tex`. They are given in Sec. 10 of the note, checked against first principles.

## Reproducing

```
python3 checks1.py          # all checks -> checks1.json (~15 s)
python3 make_note1.py       # Group_I_check.tex from check_static.tex
pdflatex Group_I_check.tex  # three times
python3 make_corrected1.py  # Group_I_corrected.tex from input/Group_I.tex (every edit asserted)
pdflatex Group_I_corrected.tex   # twice
```

| File | Contents |
|---|---|
| `fp1.py` | first principles: conjugate LO amplitude, 𝔹₃ (coordinate form, sign E3 as an option), rows II–V |
| `user1.py` | the equations of `Group_I.tex` as written (explicit forms and after the integrations) |
| `wf3.py` | 𝔹₃: instantaneous and two-ρ entries against the wave functions from m1, m5 |
| `integ.py` | principal-value p⁺ integration |
| `checks1.py` | all comparisons |

The scripts use `../groups23/fp.py` and `../three_rho/numerics.py`. They require `numpy`, `scipy` and `sympy`.
