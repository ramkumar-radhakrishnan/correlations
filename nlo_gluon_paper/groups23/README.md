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

## Corrections to `Calculation_main.tex`

The full list of corrections to `Calculation_main.tex` (all sections, not only those used by Groups II and III), with a marked-up copy of the file, is in [`main/`](main/README.md).

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

## Group II, Row IV with d⊥ kept general

[`rowIV_d/RowIV_dperp.pdf`](rowIV_d/RowIV_dperp.pdf) (5 pages; source `rowIV_d/RowIV_dperp.tex`) redoes Row IV from the polarization sum on, with d⊥ never set to 2. Every transverse trace is δ^{mm} = δ^{ii} = d⊥. The corrected block of `Group_II_corrected.tex` (Eqs. 1.17–1.22) now uses this form.

- **Polarization sum.** Σ_ij V V = [d⊥ − 2P/p⁺ − 2P/k⁺] δ_km δ_k′m′ + [(P/p⁺)² + (P/k⁺)²] δ_mm′ δ_kk′ + (2P²/p⁺k⁺) δ_mk′ δ_km′. Checked for d = 2…5.
- **w, w′ integrations.** The factor (P/k⁺)² holds in any d. The Jacobian (P/k⁺)^d is offset by the scaling (P/k⁺)^{1−d} of the d-dimensional kernel 𝒦_d.
- **Tensor z integral (exact in d).** ∫_z 𝒦_d^m(y−z) 𝒦_d^{m′}(y′−z) = (π^{2−d/2}/2) (Δ²)^{1−d/2} [Γ(d/2−1) δ^{mm′} − 2Γ(d/2) Δ̂Δ̂] = −(π/2) E δ^{mm′} − π Δ̂Δ̂ + O(ε). Checked by quadrature in d = 2 (traceless part) and d = 3.
- **Result.** With Ē = 1/ε̄ − log(Δ²μ²e^{2γ_E}/4), the braces are {X·X′[d⊥ Φ̂(ξ) Ē − 2(ξ/(1−ξ) + (1−ξ)/ξ)] − 2 d⊥ ξ(1−ξ) (X·Δ̂)(X′·Δ̂)}, with prefactor (1/(2π)³)(g⁴N_c/16π⁴)(1/k⁺) ∫_{k⁺/∨}^{k⁺/(k⁺+Λ)} dξ/ξ².
  - The endpoint coefficient is d⊥Ē − 2.
  - The pole is (α_s/2π)(1/ε̄) ∫ dξ/ξ³ (d⊥/2) P_gg^{(r)} ⊗ LO.
  - At d⊥ = 2 − ε this reduces to the previous result.

- **Final result as a sum of terms:** [`rowIV_d/RowIV_final.pdf`](rowIV_d/RowIV_final.pdf) (source `RowIV_final.tex`, written by `make_final.py`). It has ten terms: two endpoint terms and eight ξ integrals. Their sum equals the compact form to 3×10⁻¹⁵.

```
python3 rowIV_d/check_d.py      # -> rowIV_d/check_d.json (~2 min)
python3 rowIV_d/make_rowIV.py   # RowIV_dperp.tex from rowIV_static.tex (run inside rowIV_d)
```

### Fragmentation function and DGLAP from Row IV

[`rowIV_d/RowIV_FF_DGLAP.pdf`](rowIV_d/RowIV_FF_DGLAP.pdf) (6 pages; source `RowIV_FF_DGLAP.tex`, written by `make_frag.py` from `frag_static.tex` and `frag_checks.json`).

- **Your ξ form** is correct. The only change is the measure ∫_{x,x′,y,y′}, in both terms.
- **Convolution form (exact, any d⊥).** d³N_IV = (α_s/2π) ∫_{ξ₀}^{ξ₁} dξ/ξ³ {P_gg^{(r)}(ξ) 𝓛(k/ξ) + 2N_c ξ(1−ξ) 𝓛_D(k/ξ)}.
  - 𝓛 is LO with the weight (d⊥Ē − 2)/2. In CDR this weight is 1/ε̄ − ℓ, so the residue is exactly P_gg^{(r)}.
  - 𝓛_D is LO with the traceless tensor δ − d⊥Δ̂Δ̂. It has no pole.
- **Plus prescription with exact limits.**
  - Row IV gives P^{(r)}_{gg,+} ⊗ LO plus a δ(1−ξ) soft term 2N_c log((k⁺+Λ)/Λ).
  - The regular remainder is R_Λ = O(Λ/k⁺).
- **MS-bar counterterm.** Γ_gg = δ(1−ξ) − (α_s/2π)(1/ε̄) P_gg.
- **DGLAP.** μ² ∂D/∂μ² = (α_s/2π) ∫_ξ¹ dz/z P_gg(z) D(ξ/z, μ). It is derived two ways: from the μ independence of the cross section, and from that of the bare D.
- **δ(1−ξ) tally with Group I** (Row II and Row V 2b UV, both with c.c., plus the charge renormalization of LO):
  - δY cancels with the weights 1:2:2.
  - The 11/3 and the powers of δ cancel between Row II and 2b.
  - β₀/2 = 11N_c/6 is the charge renormalization. It also follows from the momentum sum rule.
  - Exact remainder: N_c[3 log((1+δ)/(1−δ)) − log(1−k⁺/∨)] LO − R_Λ, which is O(Λ/k⁺, k⁺/∨). This is open.

```
python3 rowIV_d/frag_checks.py  # -> rowIV_d/frag_checks.json
python3 rowIV_d/make_frag.py    # RowIV_FF_DGLAP.tex
```
