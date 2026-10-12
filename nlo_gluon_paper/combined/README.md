# Group I Row II in coordinate space, the three-ρ cancellation, and the JIMWLK / BFKL evolution

[`RowII_FT_and_evolution.pdf`](RowII_FT_and_evolution.pdf) (107 pages) puts everything about the large logarithms in one note:

- **Part I (new):** the Fourier transform of the three-ρ term of Group I, Row II, with every step written out.
- **Part II:** the cancellation of the symmetrized three-ρ terms ∝ log(∨/Λ), from [`../three_rho`](../three_rho).
- **Part III:** the two-ρ terms ∝ log(∨/Λ) compared with JIMWLK⊗LO, from [`../two_rho`](../two_rho).
- **Part IV:** the coefficient of log∨ compared with BFKL⊗LO. It reuses [`../bfkl`](../bfkl); the comparison tables are new.
- **Part V:** the leading-log cross section.
- **Appendices:** the step-by-step reorderings, colour relations and checks.

## Result

- **Fourier transform.** The bracket in the message lacks the term −k^i log(∨/k⁺) that `Group_I.tex` and `main.tex` have. The ξ-integral form of `main.tex` confirms that the term is needed.
  - With the full bracket, 𝔅 = A·log(k⁺/Λ) + R, where A = −E − 2k − D.
  - Only the part odd under x↔y contributes. The transform of A/(k²p²) is 2πi(𝔉_β + 𝔉_α):
    - 𝔉_β = 𝒦(x−w)∫K(w,y,z) − 𝒦(y−w)∫K(w,x,z);
    - 𝔉_α = −[𝒦(x−w) − 𝒦(y−w)]∫K(x,y,z).
  - The prefactor is −(i/2)𝒩 for log(∨/Λ) and +(i/2)𝒩_P for log(∨/k⁺).
- **Nothing was missing in the earlier input.** `Evolution_three_rho.tex` has exactly A·log(∨/Λ) for this row. The other pieces of the full bracket contain no Λ.
- **Three-ρ terms ∝ log(∨/Λ).**
  - As written, 20 of 48 colour networks do not cancel. Every left-over comes from Group I (Row II α part; Row III part 2 / 𝒟).
  - After fixes F3 and F4, all networks cancel.
- **What F4 means.** The −2k^i term of the bracket is exactly the commutator [ℕ^(2), 𝔸^(1)]. `main.tex` sets this commutator to zero.
  - F4 is the p⁺ ordering of the virtual factor: soft part (log(k⁺/Λ)) to the left of 𝔸^(1), hard part (log(∨/k⁺)) to the right.
  - With this ordering, the gluon exchanged between the two charges (𝔉_α) comes with log(∨/k⁺), as part of the projectile evolution (BFKL).
- **Two-ρ terms ∝ log(∨/Λ)** = JIMWLK⊗LO exactly after F1–F6. This is the target evolution.
- **Coefficient of log∨** (log(∨/Λ) plus log(∨/k⁺)), after P1, F4, F5 and F1, F2:
  - two-ρ terms = BFKL⊗LO exactly (9 kernel structures, exact coefficients);
  - symmetrized three-ρ terms cancel (4 structures, 44 networks).
- **Together:** d³N_LL = log(k⁺/Λ)·JIMWLK⊗LO + log(∨/k⁺)·BFKL⊗LO.

## Reproducing

The scripts use the modules of `../two_rho`, `../three_rho` and `../bfkl`.

```
python3 ft/ft_checks.py        # Part I: partial fractions, split of the bracket, xi integral, commutator identity -> ft/ft_checks.json
python3 check_all.py           # pointwise totals: log(vee/Lambda) vs JIMWLK x LO, log(vee) vs BFKL x LO -> allstates.json
python3 build_P.py 3           # coefficient of log(vee), symmetrized three-rho terms -> bookP3.pkl
python3 build_P.py 2           # coefficient of log(vee), two-rho terms and BFKL x LO -> bookP2.pkl
python3 reduce_run.py bookP2.pkl reducedP2.pkl   # exact coefficients, basis reduction
python3 reduce_run.py bookP3.pkl reducedP3.pkl
python3 make_combined.py       # writes RowII_FT_and_evolution.tex
pdflatex RowII_FT_and_evolution.tex   # three times
```

`make_combined.py` builds the document from:
- the hand-written parts in `combined_static.tex` and `ft_part.tex`;
- the generated `.tex` files of `../three_rho`, `../two_rho` and the templates of `../bfkl`.

Labels from those notes get the prefixes `R:`, `C:`, `T:`, `B:`.

| File | Role |
|---|---|
| `ft_part.tex` | Part I (Fourier transform, the α commutator) |
| `combined_static.tex` | Summary, conventions, part introductions, Part V |
| `allsets.py` | All terms by row: coefficient of log(∨/Λ) (A) and of log∨ (P), as written and corrected |
| `check_all.py`, `build_P.py`, `reduce_run.py` | The checks and exact tables quoted in Parts IV and V |
| `ft/ft_checks.py` | The checks of Part I |
| `make_combined.py` | LaTeX assembly |

Requires `numpy`, `sympy` and `scipy`.
