# Two-ρ terms ∝ log(∨/Λ) and the evolution equations

[`TwoRho_evolution.pdf`](TwoRho_evolution.pdf) collects every two-ρ term proportional to log(∨/Λ) and compares the sum with the BFKL and JIMWLK evolution of the LO cross section (Appendices E and F of the main calculation).

The terms come from three places:
1. reordering the three-ρ terms (the note in [`../three_rho`](../three_rho)), including the remainders A1–B2 of the four-ρ terms;
2. reordering the four-ρ terms ([`input/4rho.tex`](input/4rho.tex), "Two ρ contribution");
3. the genuine two-ρ terms of [`input/Evolution_two_rho.tex`](input/Evolution_two_rho.tex). This includes the log(∨/Λ) part of Group I, Row V (part 2b), which is extracted in Sec. 4 of the note.

## Result

- **The coefficient of log(∨/Λ) is JIMWLK⊗LO, not BFKL⊗LO.**
  - With log(Λ/k⁺) = −log(∨/Λ) + log(∨/k⁺), the coefficient of log(∨/Λ) is the soft-gluon limit p⁺ ≪ k⁺. That limit evolves the target, which is JIMWLK.
  - BFKL⊗LO multiplies log(∨/k⁺).
  - You can also see this in the terms: in every log(∨/Λ) term the measured gluon is emitted by a charge. In BFKL⊗LO it is also emitted by the new gluon at z.
- **As written:** 32 of the 38 kernel structures differ from JIMWLK⊗LO.
- **After six fixes the sum equals JIMWLK⊗LO exactly.** The match is pointwise in z, holds in all 13 kernel structures, and the coefficients are exact polynomials in N_c. No two-ρ term is left over.
  - **F1** `4rho.tex`, first part: drop the −8N_c 𝒦(x′−w′)·𝒦(x′−w) K(y,x′,z) term.
  - **F2** `4rho.tex`, second part, dN⁽¹⁾ term 2: one index is repeated four times. The correct term is 3 f^{abc} f^{deg} ρ^{b′}(x′) ρ^c(y) U^{ag}(y) U^{db′}(w′) U^{eb}(z).
  - **F3** Group I, Row III: part 2 is replaced by 𝒟 (three-ρ note). Its two-ρ terms need the charge order ρ(x′)ρ(x)ρ(y).
  - **F4** Group I, Row II: 𝔉 → 𝔉_β in the U(w)𝔸⁽³⁾ term (three-ρ note). The commutator form gives the same two-ρ terms.
  - **F5** Group III, Row III, part 2: reverse the overall sign.
  - **F6** Group I, Row II (genuine two-ρ):
    - The bracket −2/ε − 2γ − log((x−w)²μ²/4) must be (1/π)∫_z[K(x,w,z) − K(w,w,z)].
    - That equals 1/ε̄ − log(4π²μ²(x−w)²).
    - The (x−w) dependence was already right; only the pole and the constant change.
- **Row by row:**
  - Group II and Group III agree with an independent leading-log computation, after F5.
  - Group I agrees after F3, F4 and F6.
- **Appendix F as pasted is still wrong in two places:** M sits at the wrong points in dN⁽¹⁾–dN⁽⁴⁾, and dN⁽⁵⁾ has the wrong sign. The corrected form equals −δY H_JIMWLK on LO exactly.
- **Convention:** U^{bc} = U^{cb} throughout, as in the notes.

## Reproducing

The scripts use the modules of [`../three_rho`](../three_rho) (`engine.py`, `sources.py`, `numerics.py`, …) and the reference rows in `../three_rho/reference`. The JIMWLK Hamiltonian acting on LO is in `reference/hsym.py`.

```
python3 build_tables.py   # all two-rho terms after the fixes, grouped by kernel and colour network -> book.pkl
python3 reduce.py         # exact coefficients (SU(3) and SU(4)), basis reduction -> reduced.pkl
python3 states.py         # number of differing structures, as written and after each fix -> states.json
python3 diag_rows.py      # row-by-row comparison with the reference -> diag_rows.json
python3 appF.py           # Appendix F as written / corrected vs JIMWLK x LO, and BFKL x LO -> appF.json
python3 variants.py       # other orders of D, other forms of fix F4, other seeds, non-symmetric U -> variants.json
python3 make_two_tex.py   # writes TwoRho_evolution.tex from two_static.tex
pdflatex TwoRho_evolution.tex   # twice
```

Each step takes from a few seconds to about ten minutes (`variants.py`).

| File | Role |
|---|---|
| `core.py` | Term format, renaming, Weyl symmetrization, pointwise colour fingerprints, kernel keys |
| `users.py` | Two-ρ terms from reordering the three- and four-ρ terms; genuine two-ρ terms of `Evolution_two_rho.tex`; 𝒟 |
| `user4.py` | Transcription of the two-ρ lists of `4rho.tex` (and fixes F1, F2) |
| `target.py` | JIMWLK⊗LO |
| `reference.py` | Region-T reference rows (independent leading-log computation) |
| `tables.py`, `build_tables.py`, `reduce.py` | Exact tables by kernel structure and colour network |
| `compare.py`, `rowcmp.py`, `rowcmp2.py`, `fullpt.py` | Pointwise comparison helpers |
| `states.py`, `diag_rows.py`, `appF.py`, `variants.py` | The checks quoted in the note |
| `check_ref.py` | Reference rows = JIMWLK⊗LO |
| `check_user4.py`, `check_user4b.py`, `fix_f22.py` | Checks of the `4rho.tex` lists against the exact Weyl reordering; search for fix F2 |
| `check_prev.py` | The two-ρ terms of the three-ρ note vs the exact Weyl reordering (needs `../three_rho/results.pkl` from `process.py`) |
| `g3test.py`, `g1N.py`, `fullN.py` | Checks of fixes F5 and F6 |
| `make_two_tex.py`, `two_static.tex` | LaTeX generation |

Requires `numpy` and `sympy`.
