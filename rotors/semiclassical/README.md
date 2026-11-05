# Semiclassical entanglement of two coupled quantum rotors

Standalone paper covering the semiclassical regime (deep rotation, generic
energies and the separatrix, deep libration), plus a section on Wigner and
Husimi distributions. Target journal: Annals of Physics.

* `rotor_semiclassical.tex` / `.pdf` - the paper (elsarticle, compile with pdflatex three times).
* `wigner_husimi.py` - numerics for Sec. 7: the relative Wigner function and its
  checks, Wehrl entropies of the quantum and classical momentum distributions,
  the fringe-model prediction D(sigma), and the deep-rotation (Bessel) table.
  Writes `wigner_husimi.npz`; the printed output is in `wigner_husimi.out`.
* `figs_wigner.py` - Figs. 5 and 6 from `wigner_husimi.npz`.
* `rotorlib.py` - shared exact-diagonalisation and classical-quadrature routines
  (same as `../expanded/rotorlib.py`).
* Figs. 1-4 are copies of figures generated in `../expanded/`:
  `fig1_rotation.pdf` = `fig2_rotation.pdf` (`figs_rot.py`),
  `fig2_separatrix.pdf` = `fig3_separatrix.pdf` (`figs.py`),
  `fig3_fringes.pdf` = `fig4_fringes.pdf` and `fig4_libration.pdf` = `fig5_libration.pdf` (`figs_rest.py`).

All numerics use hbar = I = 1 and need Python 3 with numpy, scipy and matplotlib
(`python3.13` in this environment).

## Main new results (Sec. 7)

* The Wigner function of the reduced state is the marginal of the two-rotor
  Wigner function and equals p_k / 2pi, so its entropy is the von Neumann
  entropy identically: using the Wigner function in place of the classical
  distribution gives back the quantum answer.
* The constant 1 - ln 2 is the entropy of the chord (interference) terms of the
  Wigner function, which survive the integration over the relative angle.
* Wehrl entropy in deep rotation: S_W(sigma) - S_q = (1 - ln 2) - D(sigma),
  D(sigma) = (2/pi) int_0^{pi/2} f(exp(-2 sigma^2 b^2)) db,
  f(A) = ln[(1 + sqrt(1 - A^2))/2] + 1 - sqrt(1 - A^2), D ~ 0.0802/sigma.
  This does not depend on gamma; it is checked against Bessel weights up to gamma = 1e4.
* Quantum and classical distributions smoothed by the same coherent state agree
  to O(1/sigma), but both entropies depend on sigma. The coherent-state width
  plays the role of the arbitrary cell in the oscillator problem.
