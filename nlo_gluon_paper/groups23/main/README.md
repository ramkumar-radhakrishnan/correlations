# Corrections to `Calculation_main.tex`

[`Main_corrections.pdf`](Main_corrections.pdf) (9 pages) lists every correction with the corrected equation: red is the text as written, blue is the correction. [`Calculation_main_corrected.tex`](Calculation_main_corrected.tex) is the file itself with all corrections marked in place; [`Calculation_main_corrected.pdf`](Calculation_main_corrected.pdf) is that file compiled. The file as received is [`../input/Calculation_main.tex`](../input/Calculation_main.tex).

## How it was checked

- **Wave functions.** Every wave function was rebuilt from the matrix elements m1–m5 (Rayleigh–Schrödinger, ½ per identical two-gluon state) and compared with the notes point by point before any integration (SU(3), random momenta).
- **Loop integrals.** The transverse loop integrals were redone in dim reg (d = 2 − ε, MS-bar) with Feynman parameters.
- **Longitudinal integrals.** Every p⁺/ξ integral was redone exactly and compared with quadrature. Only ∨ > p⁺, k⁺ > Λ is used; no Λ → 0 or ∨ → ∞.
- **Ω coefficients.** The coefficients of Ω were compared with the verified wave functions.

## Result

**Errors that change results:**
- **E1**, onegrho12: overall sign; its ξ-form is right.
- **E2**, Eq:onegrhoX1: overall sign (− instead of +).
- **E3**, 𝔹₃: the non-instantaneous Θ(p⁺−k⁺−Λ) entry needs a minus, in momentum space and in both coordinate forms. With this sign every 1/Λ term cancels.
  - As written, 𝔹₃/ψ = −8π³ for this entry and +8π³ for the other one.
- **E4**, 𝔸⁽³⁾₂ in momentum space: π²/3 − 67/9 is counted twice and the pole of Ψ^B is missing. Replace with the bracket of grhoposition.
- **E5**, unitarity: N⁽²⁾𝔸⁽¹⁾ = 𝔸⁽¹⁾N⁽²⁾ is false; keep −[N⁽²⁾, 𝔸⁽¹⁾₂].
- **E6**, 𝔹₂ in the (k⁺, p⁺) forms: Θ(k⁺−Λ−p⁺) → Λ < p⁺, Λ < k⁺, k⁺+p⁺ < ∨. Otherwise 𝔹₂ = 8π³ × wave function, term by term.
- **E7**, Appendix F:
  - in dN⁽¹⁾–dN⁽⁴⁾ the Wilson-line matrix must sit at the points of 𝒦(u,v,t);
  - dN⁽⁵⁾ has sign +.

**Typos (the next line is right):**
- **E8**, part2: δ^{ni} → δ^{ji}.
- **E9**, stseve: the energy denominators.
- **E10**, onegrho22: "+" → product.
- **E11**, A1coeff: argument −z.

**Exact forms replacing every Λ → 0 / ∨ → ∞ step (X1–X7):**
- **gtwoloop1, onerho2:** closed forms with Li₂(Λ/k⁺).
- **instgg:** k⁺/Λ − k⁺/(k⁺−Λ) + ln((k⁺−Λ)/Λ).
- **grho3 + tworho2:** exactly the notes' result with ln(Λ/k⁺) → ln(Λ/(k⁺−Λ)). All 1/Λ terms cancel identically.
- **Three-ρ terms:** the 1/Λ − 1/∨ terms cancel exactly, so no approximation is needed.
- **combinedgrho / grhoposition and combinedgrhorho:** exact ξ limits k⁺/∨ < ξ < k⁺/(k⁺+Λ), with k⁺ < ∨ − Λ for the Ψ^B and R parts.

**Only for a general adjoint U (G1–G4):** 𝕏₁, eq: A1a, the item C̄A, and seven placements in CS2. With the convention U^{bc} = U^{cb} the written and corrected forms are equal.

**Minor (M1–M8):** notation in the dim-reg steps, typos, and the order-g⁴ normalisation coefficient (3/8, not −1/8; not used).

**Verified as written:** the whole vacuum two-gluon section, Ψ³–Ψ⁵ and the 4+5 combination including 4+5onerho (pole, ln μ and finite parts), the three-ρ chain, the final sums, 𝔹₂, N, and the cross-section pieces listed in the note.

## Reproducing

```
python3 run_checks.py           # all checks -> main_checks.json, exact_forms.json (~2 min)
python3 make_note.py            # Main_corrections.tex from main_static.tex
pdflatex Main_corrections.tex   # three times
python3 make_corrected_main.py  # Calculation_main_corrected.tex (+ _build.tex for this machine)
pdflatex Calculation_main_corrected_build.tex   # twice; article class, diagrams replaced by boxes
```

| File | Contents |
|---|---|
| `lcpt1.py` | one gluon → one gluon, non-instantaneous, from m2, m3 (E1, E10) |
| `b3_check.py` | non-instantaneous entries of 𝔹₃ against the wave functions (E3) |
| `lcpt_vac.py` | vacuum one-gluon terms at order g³ from m1–m5 (instgg, gg-uu, 3grhorho, psif, psiff, 4+5tworho) |
| `b2mom_check.py` | momentum-space 𝔹₂ against the vacuum two-gluon wave functions |
| `onerho_dimreg.py` | commutator parts of 3grhorho and psif+psiff in dim reg (X3, 4+5onerho) |
| `tworho2_check.py` | tworho2 against the exact ξ integral |
| `exact_forms.py` | all exact forms X1–X7 against quadrature |
| `make_corrected_main.py` | the marked-up copy (every edit asserted to match) |

The scripts use `../../three_rho/numerics.py`. They require `numpy`, `sympy` and `mpmath`, plus `pymupdf` for previews only. `Calculation_main_corrected.tex` keeps the original preamble (revtex4-2, tikzit), which is not installed here. The PDF is built from `Calculation_main_corrected_build.tex`, which differs only in the preamble.
