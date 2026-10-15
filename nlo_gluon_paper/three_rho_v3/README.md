# Three-ρ terms ∝ log(∨/Λ): exact cancellation by class and number of Wilson lines

[`ThreeRho_cancellation_v3.pdf`](ThreeRho_cancellation_v3.pdf) (source [`ThreeRho_cancellation_v3.tex`](ThreeRho_cancellation_v3.tex)) checks the list of all three-ρ terms with c.c. as sent ([`input/three_rho_all_v3.tex`](input/three_rho_all_v3.tex)) and shows how they cancel.

The terms are organised by kernel class (α, α′, β, β′) and by the number of Wilson lines (2, 3, 4). Each class and number of Wilson lines is one continuous equation: the common prefactor, f^{abc} and the kernel are in front, and every term follows as written, one per line, with its source at the right. Terms that cancel each other have the same colour, and the last lines give each colour's sum, which is zero.

## Result

- **All terms cancel exactly.** There are 152 terms. After symmetrisation of the charges they form 40 colour networks: α 8, α′ 8, β 12, β′ 12. The coefficients of every network add up to zero.
- **Term against term.** In each network the cancelling terms have literally the same colour structure: the same charges with the same chain of Wilson lines, up to the order of the indices of f^{abc}. No identity is needed, neither UᵀU = 1 nor f-invariance.
- **Separate cancellation.** The classes do not mix, and within a class each network contains terms with a single number of Wilson lines. So the 2-, 3- and 4-Wilson-line terms each cancel among themselves.

  | class | 2 U | 3 U | 4 U |
  |---|---|---|---|
  | α, α′ | 4 terms / 2 networks | 18 / 4 | 6 / 2 |
  | β, β′ | 6 / 2 | 36 / 8 | 6 / 2 |
- **Numerical check.** With random SU(3) Wilson lines and charges, every kernel monomial has residual ≤ 2.4·10⁻¹⁵, while single terms are O(1).
- **Same as the corrected list.** Group by group, the list equals the corrected list of [`../three_rho_v2/ThreeRho_final.pdf`](../three_rho_v2/ThreeRho_final.pdf) (T1–T3, F3, F4).
- **Convention.** The list uses U^{bc} = U^{cb}. Read literally with a general adjoint U, it does not cancel. It does cancel for a general U if every chain is read as U(p₁)⋯U(pₙ)ρ, oriented from f towards the charge. So the cancellation does not use U² = 1.

## Parsing conventions

- **Coefficients.** The coefficient of a term is its sign × 2 (for "2\,") × ½ (for "½") × 2 (for an anticommutator). Only the symmetrised charge part is compared: ½ρ{ρ,ρ} counts as ρρρ, and {ρ,ρ}ρ as 2ρρρ.
- **Typed forms.** `\-` at the start of a term is read as a minus sign; leading `\quad` is ignored.
- **Commented lines.** Lines starting with `%` are dropped.

## Reproducing

```
python3 parse3.py        # reads the list; index and kernel consistency
python3 check_cancel.py  # pointwise cancellation, U = Uᵀ and general U -> cancel_check.json
python3 networks3.py     # colour networks and their sums -> networks3.json
python3 general_u.py     # general U, chains read as U(p1)...U(pn) rho
python3 cmp_v2.py        # group-by-group comparison with ../three_rho_v2 (corrected list)
python3 make_note3.py    # writes ThreeRho_cancellation_v3.tex
pdflatex ThreeRho_cancellation_v3.tex   # twice
```

The scripts use the modules of `../two_rho`, `../three_rho`, `../bfkl` (through `env.py`) and `../three_rho_v2` (for `cmp_v2.py`). They require `numpy` and `sympy`.
