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

## `trijet/` — rapidity divergence of LO trijet production at small x

How to phase-space-slice the divergences out of the leading-order `γ*+A → q q̄ g + X`
cross-section in the CGC: finite resolved trijet + finite virtual dijet, with the rapidity
(slow-gluon) logarithm absorbed into B-JIMWLK/BK evolution — including the lifetime-ordering
constraint that prevents over-subtracting it.

- **[`trijet/trijet_rapidity_divergence.pdf`](trijet/trijet_rapidity_divergence.pdf)** (18 pp) — the method.
- **[`trijet/trijet_step_by_step.pdf`](trijet/trijet_step_by_step.pdf)** (13 pp) — the calculation done
  explicitly, ending in complete tables of the divergent and the finite terms.
- **[`trijet/trijet_divergence_separation.pdf`](trijet/trijet_divergence_separation.pdf)** (10 pp) —
  the separation done directly on the user's own LO trijet expressions by plain subtraction of the
  `1/ξ` slow-gluon pole: no slicing, no cutoffs, full analytic divergent and finite parts.
- **[`trijet/checks/`](trijet/checks/)** — numerical and symbolic verification of every integral,
  every soft limit, and the full algebraic assembly.

Sources in `trijet/pdf/`; index in [`trijet/README.md`](trijet/README.md).
