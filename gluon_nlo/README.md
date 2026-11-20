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
