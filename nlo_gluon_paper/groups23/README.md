# Group II and Group III cross sections: check against `Calculation_main.tex`

[`Groups_II_III_check.pdf`](Groups_II_III_check.pdf) lists the mistakes and their fixes. The corrected notes are [`Group_II_corrected.pdf`](Group_II_corrected.pdf) and [`Group_III_corrected.pdf`](Group_III_corrected.pdf) (sources `*.tex`; changes in blue, removed text in red). The notes as received are in `input/`.

Every displayed equation was transcribed into code (`user23*.py`) and compared with a first-principles evaluation (`fp.py`). The reference is built from the coefficients of `Calculation_main.tex` (𝔸⁽¹⁾, 𝔹₂, ℂ₁), the two-gluon amplitude and Eq. (CS2). The comparison uses:
- random transverse points;
- random SU(3) adjoint Wilson lines;
- random charges;
- exact p⁺ integrals over Λ < p⁺ < ∨ − k⁺, with no limits taken in ∨, k⁺ or Λ.

## Result

- **Correct:**
  - Group II Row I: parts 1 and 2, including all six integrated parts 2a–2f.
  - Group II Row IV up to the reduction of the Wilson lines.
  - Group III Rows I and II.
  - Group III Rows III, V and VI apart from their overall sign: the w′ integration, the log(∨/Λ) coefficients and the Row V simplification are right.
- **M1:** Group III Rows III, V and VI have the wrong overall sign in every line. In Eq. (CS2) they enter as −4 ℭ′†𝔅 and −2 ℭ′†𝕏′₁,₂, because ℭ = −ℭ′. This is the same error as fix F5.
- **M2:** the plus prescription extends the p⁺ integral down to 0. With lower limit Λ it is exact; from 0 it is off by 0.3–10 % in the tests.
- **M3:** Group II Row IV, from the z integration on.
  - The constant in the log should be e^{2γ_E}/4, not 4π².
  - The ξ limits should be k⁺/∨ < ξ < k⁺/(k⁺+Λ); the endpoint log is log(∨/Λ) − log(∨/k⁺) + log(1−k⁺/∨) + log(1+Λ/k⁺).
  - The fragmentation measure should be dξ/ξ³.
  - The δ(1−ξ) term of D_{g/g} is the LO cross section itself, which this row does not contain.
- **M4:** index placements of the adjoint Wilson lines. They are correct only for U = Uᵀ; a table gives the general-U forms.
- **M5:** typos in the integration measures.
- **M6:** in `Calculation_main.tex`:
  - 𝕏₁ needs U^{db}(v), not U^{†db}(v);
  - the Θ of 𝔹₂ should read p⁺ < ∨ − k⁺;
  - two smaller slips.

## Reproducing

```
python3 checks23.py            # -> checks23.json (~2 min)
python3 make_note.py           # Groups_II_III_check.tex from check_static.tex
python3 correct2.py            # Group_II_corrected.tex from input/Group_II.tex
python3 correct3.py            # Group_III_corrected.tex from input/Group_III.tex
pdflatex <file>.tex            # twice (three times for the check note)
```

| File | Contents |
|---|---|
| `fp.py` | first principles: coefficients, amplitude pieces 𝔅, 𝕏′₁, 𝕏′₂, ℭ′, products |
| `user23.py`, `user23b.py`, `user23c.py` | the equations as written (before and after the integrations) |
| `fixed23.py` | the same with the index placements for a general U |
| `run1.py` … `run5.py` | the individual comparisons (also collected in `checks23.py`) |

`fp.py` uses `../three_rho/numerics.py`. Requires `numpy`, `scipy`, `sympy`, and `pymupdf` only for previews.
