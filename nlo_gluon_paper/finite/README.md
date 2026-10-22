# Finite terms

Check, segregation (running coupling, DGLAP), simplification and final result of the finite terms
(`input/Finite_terms.tex`), with the comparison to arXiv:1411.2869 (`input/1411.2869v2.pdf`).

| Note | Content |
|---|---|
| [`Finite_terms_exact.pdf`](Finite_terms_exact.pdf) | **Current version.** Exact in ∨ (a = k⁺/∨ kept), only Λ → 0. Final result, Eq. (11.1); AABKL's final result next to it, Secs. 12–13; questions Q1–Q11. |
| [`Finite_terms_final.pdf`](Finite_terms_final.pdf) | Earlier version, with k⁺/∨ → 0 (superseded). |

## Scripts

| File | Role |
|---|---|
| `ft.py` | `Finite_terms.tex` transcribed term by term |
| `verify.py` | every row against the exact p⁺ integral of its first-principles integrand → `verify.json` |
| `simp_exact.py` | exact-∨ forms of all rows (closed forms with E₁) against the file at ∨/k⁺ = 1.3–300 → `simp_exact.json` |
| `checks_exact.py` | pole tally with log(1−k⁺/∨) (sympy), four-ρ lines with range ∨−k⁺ → `checks_exact.json` |
| `checks.py` | three-ρ relations and the earlier checks → `checks.json` |
| `make_exact.py` | `exact_static.tex` + JSON → `Finite_terms_exact.tex` |
| `simp.py`, `crows_w.py`, `make_final.py` | the earlier (∨ → ∞) version |

```
python3 verify.py; python3 simp.py; python3 checks.py      # file check, three-rho relations
python3 simp_exact.py; python3 checks_exact.py              # exact-vee forms and tallies
python3 make_exact.py; pdflatex Finite_terms_exact.tex      # twice
```
