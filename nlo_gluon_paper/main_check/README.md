# Your new Calculation_main.tex: check, and consistency with Group I, II, III

Input: `input/Calculation_main.tex` (the latest upload).

- [`Calculation_main_check.pdf`](Calculation_main_check.pdf) (source `Calculation_main_check.tex`) is the note.
- [`Calculation_main_corrected.pdf`](Calculation_main_corrected.pdf) (source `Calculation_main_corrected.tex`) is your file with every remaining correction marked in place: red as written, blue corrected.

## Result

- Compared with the file checked before (`../groups23/input/Calculation_main.tex`), only two lines changed. Both are earlier corrections, now done correctly:
  - **E9** (stseve, line 355): the energy denominators.
  - **M8** (line 398): "Setting Λ = 0" became "Simplifying".
- Your message to me is pasted after `\end{document}` (lines 1965–2031). It should be removed.
- All other corrections of `../groups23/main/Main_corrections.pdf` are still needed (`pending.json`):
  - errors that change results: E1–E7;
  - typos: E8, E10, E11;
  - exact forms: X1–X7;
  - general-U placements: G1–G4.
- Your Group I, II, III files use the correct forms everywhere these differ, so they agree with Calculation_main.tex once it is corrected.
  - Their first lines are identical to the versions that were checked numerically against the first-principles coefficients. The only differences are the M1 sign fix in Group III and the measures (`eq_diff.json`).
  - Every later step is checked in `../full_check`.

## Scripts

```
python3 pending.py              # which corrections are still needed -> pending.json
python3 eq_diff.py              # Group files vs the versions checked against Calculation_main.tex -> eq_diff.json
python3 make_corrected_main.py  # Calculation_main_corrected.tex (+ _build.tex for this machine)
python3 make_note.py            # Calculation_main_check.tex
```
