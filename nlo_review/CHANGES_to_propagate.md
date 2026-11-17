# Changes to propagate (red = `\chg{}`, teal = `\rnote{}` in the .tex files)

## Final_terms.tex
1. Preamble: add `\newcommand{\chg}[1]{{\color{red}#1}}` and `\newcommand{\rnote}[1]{{\color{teal}[\textit{#1}]}}`.
2. After `\maketitle`: new "Conventions" block (definitions of K^i(x), K(x,y,z); phase convention; the three plus prescriptions; scheme; principal value).
3. JIMWLK colour-singlet formula: `\int_{x,,z,w,w'}` -> `\int_{x,y,z,w,w'}`; below it `\mathcal{K}(x,y,z)` -> `K(x,y,z)`.
4. DGLAP: after your P~_gg formula add beta_0 = 11 N_c/3 and the statement that both plus prescriptions act on
   G(xi) = (1/xi^2) int e^{-ik.(y'-y)/xi}[...], with G(0) = -(4 pi/k^2) int_{x,x',z} K(x-z).K(x'-z)[U(z)-U(x')][U(z)-U(x)] rho rho.
5. Running coupling: note that this term is minus the (beta_0/2) delta(1-xi) part of DGLAP.
6. "Finite terms are obtained by ...": add "after removing the poles in MS-bar".
7. Three-rho heading: note on symmetrising rho rho rho.
8. Group I Row II (3 rho): note on the p -> 0 region of -2 q^i log(p^2/q^2)/p^2.
9. Group I Row III (3 rho):
   - third block (prefactor -g^4 f^{abd}/(8 pi^4), no log) DELETED: principal value;
   - second block unchanged, with a note that it is not yet the consistent finite remainder.

## Limiting_GIII_RI.tex
1. Preamble macros as above.
2. Last term of N_3b: `U^{bc}(z)` -> `U^{bc}(w)`.
3. After the `log(-A) = log A + i pi` lines: note that for principal value every i pi is absent,
   so the last block of the final sum (prefactor g^4 f^{abd}/8 pi^4) is dropped.
4. End of file: status note on the last term of N_3b / second block of the sum.

## Group_III_RI.tex
1. Preamble macros as above.
2. Last term of N_3b: `U^{bc}(z)` -> `U^{bc}(w)`.
3. End of file: status note; logs of negative arguments are log|.| for principal value.

## Soft sector: Group I Row III and Group III Row V (folder `rowIII_fix/`)

Corrected copies of your latest `Other_terms.tex`, `Limiting_Cases.tex`, `Final_terms.tex` are in `rowIII_fix/` (changes in red, notes in teal).

1. **Group I, Row III, last term.** The term `log[(z-w)^2/(x-w)^2] K(x'-w').K(y-w) K(x-z).K(w-z)` is replaced by
   `(1/(z-w)^2) log[(z-w)^2/(x-w)^2] K^i(x'-w') { K^i(y-z)[U(z)-U(y)] - (x-w)^2/((x-w)^2+(z-w)^2) K^i(y-w)[U(w)-U(y)] }`
   plus the contact term `(pi^3/6) K(x'-w').K(y-w)[U(w)-U(y)]`, same prefactor `-(ig^4 f^{abd}/16 pi^5 k+)`.
2. **Group III, Row V.** The `-( [k+/(p+(p+ + k+))]_+ + 1/(p+ + k+) )` term is taken out of the p+ integral and replaced by the closed form
   `-(ig^4 f^{c'ca}/8 pi^4 k+) ln^2[2 e^{-gamma_E}/(|k||x-w'|)] K(x'-w').K(y-w) [U(w')-U(x')] U(x) [U(y)-U(w)]`.
   In `Other_terms.tex` the `log(1-k+/vee)` term of this row is removed with it.
3. **`Limiting_Cases.tex`, Row III summed form.** The `i pi` block (prefactor `g^4 f^{abd}/8 pi^4`, no `i`) is removed: principal value.
4. **Evolution list.** No change to the cancellation; two notes added in `Evolution_three_rho.tex`.
5. **Not changed, open:** Group III Rows III and VI (same `-1/p+` structure, colour factor contains `U(z)`).

Numerical checks: `rowIII_fix/checks/` (`v1_rowIII.py`, `v2_rowV.py`, `v3_j0.py`).
