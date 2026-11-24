# Classical entanglement and entropy of two coupled quantum rotors

This paper repeats the analysis of Duan, Kovner and Skokov (arXiv:2301.05735)
for two rotors coupled by lambda*cos(theta1 - theta2). It follows their sections
one by one:

| DKS                                              | here                                           |
|--------------------------------------------------|------------------------------------------------|
| Sec. II: classical entanglement entropy          | Sec. III: arcsine distribution, S_CE = ln[2 pi^2 lambda/(Delta |omega1-omega2|)] |
| Sec. III: quantum entropy from WKB               | Sec. IV: WKB plus the Jacobi-Anger expansion; eigenvalues J_j(gamma)^2/2 |
| Sec. IV: ground state, small-C excited state, KL | Sec. V: ground state and small-gamma excited state, -(1-2F)ln(1-2F) - 2F ln F |

The perturbation theory and the density matrix from the authors' notes are in
Appendix A, with the corrections listed below. Appendix B evaluates the quantum
entropy at large gamma (the phase average that gives 1 - ln 2). Appendix C lists
the numerical checks.

Files:

* `rotor_classical_entanglement.tex` / `.pdf` - the paper (REVTeX preprint).
  Compile with pdflatex three times.
* `checks.py` - every numerical check quoted in the paper; its output is in `checks.out`.
* `figs.py` - Figs. 1 and 2.
* `rotor_classical_entanglement.zip` - the tex file, both figures, the PDF and the scripts.

## Corrections to the notes (Entanglement_entropy.pdf)

1. **Eq. (33), last line.** The normalization term in the density matrix has
   coefficient lambda^2 I^2/(4 hbar^4)[1/(n+1)^2 + 1/(n-1)^2], not .../(8 hbar^4).
   The density matrix contains N^2 = 1 - eps, not N = 1 - eps/2. With 1/8,
   Tr rho = 1 + eps/2, which is wrong.
2. **Eqs. (17)-(18).** The remainders should read O(lambda^3), not O(lambda^2).
3. **Second-order state, Eq. (30).** It also requires n != +-2, not only n != +-1.
4. **Typos.**
   * Eq. (2): the derivatives should be written d^2/d(theta)^2.
   * Eq. (13): "[H.J]" should be [H,J].
   * Under Eq. (1): "m in Z" (the notes have "m \ni Z").
