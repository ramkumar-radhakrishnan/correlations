# Check of your BFKL.tex

Your file is in `input/`. The note is [`BFKL_check.pdf`](BFKL_check.pdf) (source `BFKL_check.tex`). The corrected copies are [`corrected/BFKL_marked.pdf`](corrected/BFKL_marked.pdf) (red = remove, blue = new) and [`corrected/BFKL.pdf`](corrected/BFKL.pdf) (clean).

## Result

- **The BFKL result is correct.**
  - The sum of all two-ρ lists is BFKL⊗LO exactly (1.2e-15, U = Uᵀ).
  - Both forms of "BFKL contribution" are BFKL⊗LO.
  - The colour-singlet form is the checked one (App. E).
- **Row by row:**
  - Every first form equals the log(∨/k⁺) list of its row (`../full_check/Checks_JIMWLK_BFKL.pdf`, eq. (B:…)), except Group I Row II.
  - Every "This can be written as" list is exactly the two-ρ (commutator) part of its first form (≤ 5e-14).
  - The genuine two-ρ rows equal their first forms.
  - The "From 4ρ's" two-ρ lists are identical to those of Evolution_terms.tex.
- **Corrections:**
  - **C1:** Group I Row II, first form: it is a copy of Row III d³N₁. It is replaced by eq. (B:G1II2a), with ∫_{x,y,x′} in the log term. The written list was already right.
  - **C2:** Group II Row IV, third line (H₄): the + sign is missing. It was also missing in my Checks note; fixed there.
  - **C3:** colour-singlet measure: ∫_{x,,z,w,w′} → ∫_{x,y,z,w,w′}.

## Method

- `reorder.py` computes the two-ρ part of a three-ρ term pointwise:
  - every pair of charges gives ½[earlier, later], with [ρᵃ(x), ρᵇ(y)] = i f^{abc} ρᶜ(x) δ(x−y);
  - anticommutators are kept ordered, and the c.c. reverses the order.
- `rule_vs_engine.py` shows that this rule equals the reordering engine of `../two_rho/core.py`.
- `test_tool*.py` reproduce the checked J lists of Evolution_terms.tex.
- `texfix.py` does the index sums (δ^{ij}, (x−w)^i …) term by term. texeval summed the whole expression over every index, so a term without that index was counted twice. That gave spurious differences in four rows (Group I Row IV, Group II Row I, Group III Rows I, II); with the fix they agree exactly.

## Scripts

```
python3 rerun_first.py      # first forms vs the log(vee/k+) lists     -> cmp_first.json
python3 cmp_two_rows.py     # genuine two-rho rows                     -> cmp_two_rows.json
python3 run_reorder.py      # three-rho rows: written vs two-rho part  -> cmp_reorder.json
python3 rerun_rows.py       # the four rows with mixed index terms     -> cmp_reorder_fixed.json
python3 two_total.py        # all two-rho lists vs BFKL x LO           -> two_total.json
python3 rule_vs_engine.py   # ordering rule vs the engine              -> rule_vs_engine.json
python3 make_corrected.py; python3 make_corrected.py marked
```

The other scripts (`cmp_written*.py`, `cmp_keyed*.py`, `dpairs.py`, `zdep*.py`, `dherm.py`, `dwsep.py`, `ccfit.py`, `old_vs_new.py`, `engine_tex*.py`, `g3ii_*.py`, `cc_test.py`, `total_tex.py`) are the diagnostics that located the evaluator problem. They ran before `texfix.py`, so their D ≠ 0 results are superseded.
