# Check of your Evolution_two_rho.tex, Evolution_three_rho.tex and Evolution_terms.tex

The inputs are your uploads, in `input/`. The note is [`Evolution_check.pdf`](Evolution_check.pdf) (source `Evolution_check.tex`). The corrected copies are in [`corrected/`](corrected):
- `*_marked.pdf`: red = remove, blue = new; your own colours are remapped (blue → teal, green → olive, red → violet, orange → brown);
- `*.pdf`: clean.

The reference is the log(k⁺/Λ) part (J list) of every row of your current Group I/II/III files, from the exact split of `../full_check`. Every equation is evaluated pointwise:
- random points and classical charges, U = Uᵀ, log(k⁺/Λ) = 1, ε = 0.37 kept finite in the pole brackets;
- three-ρ terms symmetrized over their charges;
- Groups I and III with their c.c.

## Result

- **Evolution_two_rho.tex is correct.**
  - All 7 rows equal the J two-ρ terms of the Group rows (≤ 2e-16). This includes Row V 2b UV/finite, G2 IV with d⊥, and the sign of G3 III part 2.
  - T1: stray comma in ∫_{w,w′,x,x′,}.
  - T2: G3 III part 2 measure ∫_{w,y′,z,y,x} → ∫_{w,x′,y′,z,x}.
- **Evolution_three_rho.tex:**
  - Every row agrees (≤ 3e-13), except Group I Row III d³N₂:
    - E1 (first form): the first two terms have the wrong sign. With the right sign they cancel the 4th and 5th terms (same kernels), and the edge term is replaced by 𝒟 (F3).
    - E2 (rewritten form): as written it is the [(x−w)²+(x−z)²] term, which cancels inside Row III. It must be 𝒟: 1/8π⁵, kernel 𝒦(x′−w′)·𝒦(y−w) 𝒦(x−z)·𝒦(w−z), and U^{bc}(z) → U^{bc}(w).
  - The list "All three ρ term with c.c." agrees row by row, already with 𝒟, and cancels: 1.3e-13 against terms of size 120–260.
  - All 12 α/α′/β/β′ blocks vanish.
  - Typos:
    - T3: G1 II log term, ∫_{x,y,x′,z} → ∫_{x,y,x′}.
    - T4: comma in ∫_{…,y′w′}.
    - T5: G2 I rewritten has 𝒦·𝒦·𝒦·𝒦 → 𝒦·𝒦 𝒦·𝒦.
    - T6: "From four ρ's" should carry log(k⁺/Λ) only.
    - T7: missing \Bigg].
    - T8: `\-` typed for a minus sign.
- **Evolution_terms.tex is correct.**
  - It is identical to the checked file apart from the title and the log.
  - Sec. 1 + Sec. 2 = Sec. 3 = JIMWLK⊗LO (8e-17).
  - T9: singlet-form measure ∫_{x,,z,w,w′} → ∫_{x,y,z,w,w′}.
  - T10: define K(a,b,z).
  - The commented-out BFKL section should carry log(∨/k⁺).
- **The corrected copies were re-checked** with the same programs (`*_corrected.json`): every row agrees, and the list still cancels.

## Scripts

```
python3 cmp_two.py        # two-rho rows vs Group rows (J)              -> cmp_two.json
python3 cmp_three.py      # three-rho rows (both forms) vs Group rows   -> cmp_three.json
python3 cmp_list.py       # "All three rho" list, total, alpha/beta blocks -> cmp_list.json
python3 terms_check.py    # Evolution_terms.tex: Sec.1 + Sec.2 = Sec.3 = JIMWLK x LO -> terms_check.json
python3 measures.py       # measures vs integrands
python3 make_corrected.py; python3 make_corrected.py marked     # corrected/
EVOL_DIR=$PWD/corrected python3 cmp_two.py    # also cmp_three.py, cmp_list.py: the corrected copies -> *_corrected.json
```
`rows_user.py` parses your equations row by row. The scripts use `../full_check/ll.py`, `../llsplit` and `../check_v4`.
