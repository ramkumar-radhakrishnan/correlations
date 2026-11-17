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
