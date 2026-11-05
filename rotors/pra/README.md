# Phys. Rev. A version

`rotor_pra.tex` / `rotor_pra.pdf`: REVTeX 4.2 (pra, two-column), 11 pages.
Compile with pdflatex three times.

This is the condensed PRA version of the material in `../expanded/` and
`../semiclassical/`. The main text states the results; the derivations are in
Appendices A-E. It contains:

* the exact reduction (Mathieu sectors, a diagonal reduced density matrix, the ln 2
  coming from exchange symmetry);
* the weak-coupling result, including the O(lambda^2) normalization term;
* the classical reduced distribution, with the cell fixed at h;
* deep rotation (Bessel weights, 1 - ln 2), generic energies and the separatrix,
  and deep libration;
* Wigner and Husimi distributions (Sec. VIII);
* the motivation from the bosonic Josephson junction and from stacked planar
  dipolar rotors.

Figures are copies:

| here | source |
|------|--------|
| `fig1_weak.pdf` | `../expanded/fig1_weak.pdf` (`../expanded/figs.py`) |
| `fig2_rotation.pdf` - `fig5_libration.pdf` | `../semiclassical/fig1-4_*.pdf` (generated in `../expanded/`) |
| `fig6_wigner.pdf`, `fig7_husimi.pdf` | `../semiclassical/fig5_wigner.pdf`, `fig6_husimi.pdf` (`../semiclassical/figs_wigner.py`) |

Still to complete before submission:

* author list, affiliations and acknowledgments;
* the full references for arXiv:2410.16338 and arXiv:2606.02656;
* the published version of the DKS paper, if it has appeared.
