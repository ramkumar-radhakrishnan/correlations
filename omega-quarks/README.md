# Omega at O(g^2): gluons + quarks

## Current manuscript (JHEP style, Paper I format)

`omega_quarks_JHEP.tex` / `omega_quarks_JHEP.pdf` — **the paper**, written in the structure and
style of the updated Paper I: *Soft gluon wave function and evolution operator with quark
interactions in the CGC at next-to-leading order* (51 pp). Single self-contained LaTeX file; build
with `pdflatex omega_quarks_JHEP.tex` (run twice). The same source split by section is in `jhep/`
(`main.tex` + `sec*.tex`, `appendices.tex`, `refs.tex`).

Sections: 1 Introduction · 2 Background and setup (LC QCD Hamiltonian with quarks, `rho = rho_g + rho_q`)
· 3 LCWFs and LCPT (i-epsilon prescription, reduction formula, unitarity identity) · 4 Ansatz for
`Omega_q`, unitarity relations, action on Fock states, and all LCWFs with diagrams (vacuum, 1/2/3
gluons, single quark/antiquark, q-qbar pair, quark-gluon) · 5 All coefficients and the generator `G`
· 6 Diagonalization at `O(g)` and `O(g^2)` · 7 Conclusions · App. A-E.

Completions relative to the earlier drafts: the coefficients `F`, `G` (gluon emission by soft
quark/antiquark, `O(g)`), `P`, `Pbar` (quark rescattering on `rho`, `O(g^2)`), `Q1,2,3`, `K1,2`;
uniform `d b` ordering for pair annihilation; t/u-channel quark exchange in `B1`; s-channel
instantaneous term in `D3`; quark loop in the gluon norm (`11 N_c - 4 T_F N_f`).

## Earlier working drafts


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
