# Check of the log(∨/k⁺) terms (`Finite_terms_1_1.tex`) against BFKL

[`BFKL_terms_check.pdf`](BFKL_terms_check.pdf) (6 pages; source `BFKL_terms_check.tex`). The input as received is `input/Finite_terms_1_1.tex`.

## Result

As written, the list does not give BFKL:
- the symmetrized three-ρ terms do not cancel. The z-independent structure C_w = f^{abc} C_Lᵃ (U(w)ρ(x))ᵇ (U(w)ρ(y))ᶜ comes with the kernel −X'·Y/(2(z−w)²), whose z integral diverges at z → w;
- the two-ρ terms leave a 1/ε pole and lack the kernel 𝒦(z−w)·𝒦(z−w')K(x,y,z) of the second line of Eq. (LOBFKL).

Row by row, against the exact coefficient of log(∨/k⁺) in the same row of `Finite_terms_1.tex`:
- 5 rows are exactly right: G1 II 2a; G1 V parts 1, 2a, 2b UV, 2b finite;
- 6 rows are wrong: G1 III, IV; G2 I parts 1, 2; G3 I, II;
- 5 rows are marked NA but contain log(∨/k⁺): G2 IV; G3 III(1), III(2), V, VI;
- G1 II part 1 needs the coefficient 2, not 1.

The cause is log(∨/k⁺) hidden inside other logarithms, e.g. log[((∨−k⁺)A + k⁺B)/(ΛA + k⁺B)]. It also sits inside log(∨/k⁺ − 1) and at the lower end ξ → k⁺/∨ of ξ integrals.

With the exact coefficients:
- the three-ρ terms cancel;
- the two-ρ terms add up to log(∨/k⁺)(ℬ − 𝒥), not ℬ alone;
- Eq. (LOBFKL) itself is ℬ and is correct.

Together with the log(∨/Λ)𝒥 of the evolution list, the total is log(k⁺/Λ)𝒥 + log(∨/k⁺)ℬ.

## Correction to `../finite_v2` and `../finite_singlet`

The lower end of the DGLAP integral (Group II Row IV) contains log(∨/k⁺)H₄, a BFKL term (Sec. 8 of the note). FS (11.1) and `Finite_terms_singlet.pdf` (7.1) had counted ½log(∨/Λ)H₄ twice. Both are corrected; the changes are in blue.

## Scripts

| File | Contents |
|---|---|
| `user_lvk.py` | your list, term by term (coefficient of log(∨/k⁺) per row) |
| `rowcheck.py` | Table 1: your terms vs the exact C_row, plus d/dlog∨ of the file → `rowcheck.json` |
| `cw_test.py` | Sec. 5: the C_w structure, exact vs your kernel → `cw_test.json` |
| `ft_asym.py` | Sec. 8: the large-q tail of the Fourier transform of log Δ² → `ft_asym.json` |
| `zw_sing.py` | auxiliary: the three-ρ rows near z → w → `zw_sing.json` (not used in the note) |
| `make_note.py`, `static.tex` | the note |

```
python3 rowcheck.py; python3 cw_test.py; python3 ft_asym.py
python3 make_note.py; pdflatex BFKL_terms_check.tex   # three times
```
