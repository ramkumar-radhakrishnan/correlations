# `gluon_nlo/` — review of single-inclusive gluon production at NLO (mid-rapidity, dilute–dense CGC)

- **[`review/NLO_gluon_review.pdf`](review/NLO_gluon_review.pdf)** (8 pp, LaTeX source alongside): independent checks
  of the final cross section (`Final_terms_original.tex`, as supplied) and comparison with the literature.
- **[`checks/`](checks/)**: every number in the review. `./checks/run_all.sh` (about 1 min) regenerates
  `checks/RESULTS.txt`. Requires `numpy`, `scipy`, `opt_einsum`.
  - `texeval.py` parses the supplied TeX directly into tensor contractions, so nothing is retyped.
  - `check_bfkl.py`: (1.1)=(1.2); singlet (1.3); (1.3) = exact one-gluon (dipole BFKL) evolution.
  - `check_jimwlk.py`: (2.1)+(2.2)=(2.3)=(2.4)= -H_JIMWLK acting on the LO cross section.
  - `kernels_g2.py`, `probe_g2.py`, `sym_g2.py`: Group II (7.3)-(7.8) singular surfaces cancel; exchange symmetry.
  - `misc_checks.py`: U -> 1, U -> constant, colour-factor symmetry.
  - `check_row4.py`: the Group II Row IV (DGLAP/RC) derivation in `Group_II_Row_IV_original.tex`.
  - `kernels_g1r3.py`, `probe_g1r3*.py`, `slopes_g1r3.py`: Group I row III (6.4)-(6.8) singularities and integrability.
  - `probe_71.py`, `probe_g2_limits.py`: integrability of (7.1) and (7.3)-(7.8).
  - `pconv.py`, `terms_r5g3.py`, `pnum.py`, `colour.py`: Group I row V and Group III (p+ integrals, plus
    prescriptions, colour factors on smooth random fields); probes `vee_g3.py`, `parity_check.py`, `diag_611.py`,
    `probe_r5g3*.py`, `probe_far_*.py` (outputs `*.out`, about 30 min on 4 cores).
