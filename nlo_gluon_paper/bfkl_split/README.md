# Evolution_terms.tex check, and Finite_terms_1.tex split into BFKL terms + remaining terms

[`Finite_terms_split.pdf`](Finite_terms_split.pdf) (source `Finite_terms_split.tex`). Your files as received are in `input/`.

## 1. Evolution_terms.tex: is Sec. 1 + Sec. 2 = Sec. 3?

Yes: exactly, with no left-over, for U^{bc} = U^{cb}. The relative difference is ≤ 1×10⁻¹⁶ at 5 random configurations (`evol_check.json`).
- Sec. 3 = JIMWLK⊗LO.
- The new list d³N₂ of "From 4ρ's" equals the two-ρ lists of 4rho.tex that were missing before. Without it the difference is 2.4 %.

## 2. Finite_terms_1.tex: BFKL terms and remaining terms, row by row

For every row: row = BFKL terms + remaining terms, exactly, for any ∨ and Λ.

- **BFKL terms** (Secs. 2–4): the log(∨/k⁺) terms of every row. Each was checked against ∂/∂log∨ of the same row (`../bfkl_check/corrected_check.json`).
- **Remaining terms** (Secs. 7–9): the rows in the form of your file, with every ∨-dependent log replaced by its remainder. The identities are checked in `split_logs.py` / `split_logs.json`.
  - Group II Row IV: its ξ integral is rewritten by hand.
  - Group III Rows III, V, VI: the [1/p⁺]₊ terms are kept unevaluated; only the subtraction moves from F(0) to F(0)·k⁺/(p⁺+k⁺). Exactly, ∫₀^{∨−k⁺} dp⁺ [1/p⁺]₊ F = ∫₀^{∨−k⁺} dp⁺ [1/p⁺]₊^{(k⁺)} F − F(0)·log(∨/k⁺); the last term is the BFKL term of the row (`pplus_check.py`).

## 3. The sum (Sec. 5)

The BFKL terms add up to log(∨/k⁺)(ℬ − 𝒥). Together with Evolution_terms.tex this gives:

  log(∨/Λ)𝒥 + log(∨/k⁺)(ℬ − 𝒥) = log(k⁺/Λ)𝒥 + log(∨/k⁺)ℬ

## Scripts

| File | Contents |
|---|---|
| `evol_check.py` | Sec. 1 + Sec. 2 vs Sec. 3 vs JIMWLK⊗LO (uses `../check_v4`) → `evol_check.json` |
| `split_logs.py` | the log identities (evaluated from the LaTeX strings), the p⁺ integrals, the ξ identity → `split_logs.json` |
| `pplus_check.py` | the [1/p⁺]₊ and [1/p⁺]₊^{(k⁺)} p⁺ integrals done numerically for growing ∨, and the identity between them → `pplus_check.json` |
| `make_split.py` | builds `Finite_terms_split.tex` (BFKL part from `../bfkl_check/make_corrected.py`) |

```
python3 evol_check.py; python3 split_logs.py; python3 pplus_check.py; python3 make_split.py
pdflatex Finite_terms_split.tex   # twice
```

## Correction (see `../bfkl_final`)

The log(∨/k⁺) (BFKL) terms used here miss two contributions:
- **(A) Group I Row III:** the edge-kernel term (x′−w′)·(y−z)/((x′−w′)²(y−z)²(z−w)²) must be replaced by 𝒟 (fix F3, as in `Evolution_terms.tex`).
- **(B) Group III Rows III, V, VI:** a log(∨/k⁺) from the region y′ → z appears after the y′ integral. It comes from the −δ_im δ_jk′/k⁺ term, and in Row III part 2 from the −(p⁺/k⁺²)[…] term.

Consequences:
- Only with (A) and (B) do the three-ρ terms cancel and the two-ρ terms give log(∨/k⁺)(ℬ − 𝒥).
- The remaining/finite terms of these rows still contain these logs.
- See `../bfkl_final/BFKL_final.pdf`, Sec. 8.
