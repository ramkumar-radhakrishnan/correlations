# correlations

## `csl/` — digitized chiral-soliton-lattice Hamiltonian

Generalization of the one-dimensional CSL lattice Hamiltonian derived in
`CSL_quantumcomp.pdf` (one qubit per site, two sites) to **arbitrary qubits per
site** ($n_q$, i.e. $Q=2^{n_q}$ field values) and **arbitrary lattice size** $N$,
with worked-out results for $(n_q,N)=(1,4)$, $(2,2)$ and $(2,4)$.

Start with **[`csl/DERIVATION.md`](csl/DERIVATION.md)**, or the typeset PDF
`csl/CSL_lattice_Hamiltonian_general.pdf` (build source in `csl/pdf/`). Raw generated output is
in `csl/RESULTS.txt`.

```
cd csl
python3 verify.py    # 57 checks, incl. exact reproduction of eq. (1.21)
python3 run_all.py   # regenerates RESULTS.txt
```

Requires `numpy` only.

## `nlo-gluon-production/` — $\beta_0$ logs in NLO gluon production

Transposition of Kovner–Lublinsky–Skokov–Zhao ([arXiv:2308.15545](https://arxiv.org/abs/2308.15545))
to single-inclusive gluon production at $O(g^4)$: which of the $\beta_0$ transverse
logarithms are running-coupling and which are DGLAP.

Read **[`nlo-gluon-production/NLO_gluon_production_DGLAP.pdf`](nlo-gluon-production/NLO_gluon_production_DGLAP.pdf)**
(build source in `nlo-gluon-production/pdf/`).

```
cd nlo-gluon-production
python3 verify.py    # 14 checks of the splitting-function and factorization algebra
```

Requires `sympy`, `numpy`, `scipy`.
