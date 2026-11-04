# Nov 4 uploads: JIMWLK and BFKL redone, left-over terms, Other_terms.tex

Inputs (`input/`): your Group_I.tex, Group_II.tex, Group_III.tex, 4rho.tex, Other_terms.tex of Nov 4.
Group_I, Group_III, 4rho are byte-identical to `../full_check/input` (Oct 31); Group_II differs only in
Eqs. (1.13), (1.16) (the typo fixes M2). Other_terms: 27/30 equations identical to `../final_terms/input` (Nov 3).

* `Summary_Nov4.pdf/.tex`: cover page: what matches, what still has to change (4rho Eq. 1.12 + F1, F2; F3; M3-M6;
  JIMWLK formula M2'), left-over terms, and the Other_terms.tex check (Group II Row IV: the log must sit inside the
  y' integral of d^3N^y_{xi LO}; equivalently +8 pi N_c xi(1-xi)/k^2 d^3N^y_LO).
* `Checks_JIMWLK_BFKL.pdf/.tex`: every step of your files, the split, 3- and 4-rho cancellation, two-rho = JIMWLK x LO
  and BFKL x LO, your formulas.
* `Remaining_RC_DGLAP_finite.pdf/.tex`: remaining terms row by row, running coupling, DGLAP, finite terms.

Rerun: `python3 step_check.py Group_I` (and Group_II) were rerun on the new files (`step_Group_*.json`); Group I is
identical to full_check, Group II differs only in the now-passing Eq. (1.16). The split / three-rho / four-rho /
same-list results (`ll_*.json`, `four_file.json`, `same_lists.json`, `rows_vs_files.json`, `step_Group_III.json`,
`mu_check.json`) are copied from `../full_check`: their inputs are identical (the split reads Group II lines
48-50, 65-82, 131-133 only). Then `python3 make_ll.py; python3 make_rest.py` and pdflatex (3x).
