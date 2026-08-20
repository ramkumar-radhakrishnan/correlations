# Omega at O(g^2): gluons + quarks

`omega_quark.tex` / `omega_quark.pdf` — construction of the soft light-cone wave function and of
the unitary evolution operator `Omega` of a fast-moving hadron in the CGC **including quark degrees
of freedom**, through `O(g^2)`. Companion to the pure Yang-Mills draft ("Paper I").

Contents:

| Section | Content |
|---|---|
| 2 | Conventions, quark sector of the LC Hamiltonian, Born-Oppenheimer + eikonal, `rho = rho_g + rho_q` |
| 3 | LCPT rules and the complete list of matrix elements (quark vertex, instantaneous quark/gluon exchange) |
| 4 | Complete normal-ordered ansatz for `Omega_q`; power counting; **all unitarity relations** at `O(g)` and `O(g^2)` |
| 5 | Action of `Omega` on Fock states; universal `(2 pi)` counting rule `K = (2 pi)^{3(n+m)/2} Psi` |
| 6 | LCWFs from LCPT: vacuum, 1-, 2-, 3-gluon and `q qbar` incoming states; normalization |
| 7 | **All coefficients** of `Omega` through `O(g^2)` |
| 8 | Diagonalization of the LC Hamiltonian at `O(g)` and `O(g^2)` |
| 9 | Hermitian generator `G`, `Omega = exp(i G)` |
| App. A-C | Spinor identities, instantaneous vertices, and a list of corrections/additions relative to the preliminary notes |

Build: `pdflatex omega_quark.tex` (twice; needs `empheq`, `slashed`, `mathtools`, `tikz`).
