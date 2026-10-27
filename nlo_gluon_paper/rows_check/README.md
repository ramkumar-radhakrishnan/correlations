# Check of Group_I.tex, Group_II.tex, Group_III.tex (every step)

The inputs are your files, copied to `input/`. There are two outputs:
- [`Group_check.pdf`](Group_check.pdf) (source `Group_check.tex`) lists the mistakes with their fixes, then a table with every step checked numerically.
- [`Group_I_corrected.pdf`](Group_I_corrected.pdf), [`Group_II_corrected.pdf`](Group_II_corrected.pdf) and [`Group_III_corrected.pdf`](Group_III_corrected.pdf) are your files with each correction marked in place: blue = new, red = remove. Your own colours are remapped to teal, olive, violet and brown.

## Result

- **M1 (changes the result): Group I Row V 2b UV, Eqs. (1.26) and (1.32).**
  - The bracket must be [2 log(k⁺/Λ) + 2 log(1 − Λ/k⁺) − 11/6]. Your file has log(k⁺/Λ) − log(∨/k⁺).
  - Row V has 0 < p⁺ < k⁺, so ∨ cannot appear. Your exact Eq. (1.25) is 2 log(∨/Λ) − 2 log(∨/k⁺) = 2 log(k⁺/Λ). The step "Setting Λ → 0" turns 2 log(∨/Λ) into log(∨/k⁺) + log(k⁺/Λ).
  - Your 2b finite subtracts the same counterterm with 2 log((k⁺−Λ)/Λ) − 11/6. Numerically, 2b = 2b finite + counterterm to 4e-16.
- **M2: Group II Row IV, Eq. (1.16).** The phase of the second term should be e^{−ik(y′−y)(p⁺+k⁺)/k⁺}. As written the next step fails by 0.27–0.35; with this phase it agrees to 2e-16. This is a typo only, since Eq. (1.17) is right.
- **M3: Group III Row III, first line.** It needs −4/(2π)³, as you already have in Rows V and VI.
- **M4: Group III Row III part 2.** The measures in Eqs. (1.8)–(1.10) should be ∫_{w,w′,x′,y′,z,x} before the w′ integral and ∫_{w,x′,y′,z,x} after it (C5).
- **M5: Group I Row II, Eq. (1.6).** The measure should be ∫_{w,w′,x′}.
- **M6: typos.**
  - Group II Eq. (1.13) has k{+}.
  - Brace mismatches (\Bigg\{ … \Bigg\}) in Group I Eqs. (1.9), (1.19), (1.20), (1.22) and (1.24).
  - Commas: y′w′ in the measure, and a stray comma in ∫_{w,w′,x,x′,}.
  - Group II repeats one sentence.
- **Every other step is correct**: Group I Rows II (part 1), III (exactly in Λ, principal value), IV and V; Group II Rows I and IV; Group III Rows I, II, III, V and VI. These agree to 1e-16–5e-12; the plus-prescription steps agree up to O(Λ).
- **The log(∨/k⁺) (BFKL) terms of Group III Rows III, V, VI and Group II Row IV are not missing.** They sit inside line 1 of those rows, in the p⁺ (or ξ) integral, and appear only after the y′ integral. The corrected copies write each of these equations, unchanged, as (B) + [equation − (B)], in blue (Group II: H₄).
- **Still open from before:** log|·| in Group I Row III; the edge kernel of Row III (F3 / C4).

## Scripts

```
python3 check_all.py Group_I     # every step of the file, two configurations -> check_Group_I.json (same for Group_II, Group_III)
python3 measures.py              # points that appear vs points integrated -> measures.json
python3 make_corrected.py        # Group_*_corrected.tex (every edit asserted)
python3 make_note.py             # Group_check.tex from the json files
python3 outline.py Group_I       # line ranges of the equations; eqnum.py: their printed numbers
```

`pe.py` evaluates your equations as written. It builds on `../llsplit/texeval.py`, adds the longitudinal variables p⁺, q⁺ and ξ as arrays, your plus prescriptions and the δ⁽²⁾ substitutions, and contracts transverse indices term by term. `quick.py` runs single steps.
