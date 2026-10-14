# All two-ρ terms ∝ log(∨/Λ), and how JIMWLK and BFKL come out

[`TwoRho_all.pdf`](TwoRho_all.pdf) (43 pages) lists every two-ρ term proportional to log(∨/Λ), source by source and term by term. The terms come from three places:
- reordering the three-ρ terms of `Evolution_three_rho.tex`, including 𝒟;
- the genuine two-ρ terms of `Evolution_two_rho.tex`;
- reordering the four-ρ terms: the lists of `4rho.tex`, plus the remainders A1–B2.

The note adds them up and compares the sums with JIMWLK⊗LO and BFKL⊗LO. No limits are put on ∨, p⁺, k⁺, Λ; only ∨ > p⁺, k⁺ > Λ is used.

## Result

- **The terms.** There are 2386 ordered two-ρ terms before collecting. After the colour algebra, they collect into 440 entries (kernel × colour structure × exact coefficient) from 18 sources.
- **Sum of the log(∨/Λ) terms = JIMWLK⊗LO exactly.**
  - This holds in all 13 kernel structures and all 154 independent colour structures, with exact N_c coefficients, pointwise in z.
  - The fixes F1–F6 are needed.
- **Why log(∨/Λ) carries JIMWLK and not BFKL, without any limit.**
  - Each row is an exact p⁺ integral.
  - Λ enters only through those limits where a gluon carries momentum Λ, i.e. is softer than the measured gluon.
  - ∨ enters only through the upper limit, where the unobserved gluon is harder than the measured one.
  - The exact identity ∫ = 𝔗·log(∨/Λ) + (𝔓 − 𝔗)·log(∨/k⁺) + R holds with unique coefficients:
    - 𝔗 is the integrand at the Λ end (JIMWLK);
    - 𝔓 is the integrand at the ∨ end (BFKL);
    - R contains neither logarithm.
  - Your Group II, Row IV shows both ends explicitly.
- **BFKL is the coefficient of log∨,** i.e. coefficient of log(∨/Λ) + coefficient of log(∨/k⁺).
  - Row by row, this is P = A for the four-ρ terms, P = 0 for Group I Row V, 𝒟 and the genuine Group I terms, and the ∨-end integrand otherwise.
  - The sum is BFKL⊗LO exactly: 9 kernel structures, 114 independent colour structures. Fix P1 is also needed.
- **Result:**
  - d³N_LL = 𝒩₀[log(∨/Λ)·𝒥 + log(∨/k⁺)·(ℬ − 𝒥)] = 𝒩₀[log(k⁺/Λ)·𝒥 + log(∨/k⁺)·ℬ].
  - Equivalently: 𝒩₀[½log(∨/Λ)(𝒥 + ℬ) + ½log(k⁺²/(∨Λ))(𝒥 − ℬ)]. The second term vanishes only for k⁺ = √(∨Λ).
- **Kernel structures:**
  - those with the unobserved gluon attached to the measured gluon occur only among the log(∨/Λ) terms;
  - those with the measured gluon emitted by the gluon at z occur only among the log∨ terms;
  - charge–charge structures occur in both, with different coefficients (`common.py`).

## Reproducing

The scripts use the term definitions of [`../two_rho`](../two_rho), [`../three_rho`](../three_rho) and [`../bfkl`](../bfkl).

```
python3 checks.py      # pointwise: sum A = JIMWLK x LO, sum P = BFKL x LO, row values -> checks.json (~1 min)
python3 build.py A     # exact tables, log(vee/Lambda) terms vs JIMWLK x LO -> bookA.pkl, reducedA.pkl
python3 build.py P     # exact tables, log(vee) terms vs BFKL x LO -> bookP.pkl, reducedP.pkl
python3 common.py      # charge-only parts of JIMWLK x LO and BFKL x LO differ -> common.json
python3 make_doc.py    # writes TwoRho_all.tex from static.tex
pdflatex TwoRho_all.tex   # three times
```

| File | Role |
|---|---|
| `sets.py` | All two-ρ terms with one label per source: coefficient of log(∨/Λ) (`A_sets`) and of log∨ (`P_sets`), with the fixes |
| `checks.py` | Pointwise tests (Tables 3 and 5 of the note) |
| `build.py` | Exact tables by kernel structure and colour structure; basis reduction |
| `common.py` | Comparison of the charge-only parts of JIMWLK⊗LO and BFKL⊗LO |
| `make_doc.py`, `static.tex` | LaTeX generation |

Requires `numpy`, `sympy` and `pymupdf` (only for previews).
