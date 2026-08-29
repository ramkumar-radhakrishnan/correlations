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

## `nlo/` — longitudinal `p^+` integration for NLO gluon production

Mathematica code performing the $\int_\Lambda^{\vee-k^+} dp^+$ integral of the
NLO gluon-production integrand in closed form (exponential integrals), plus its
$\vee\gg k^+\gg\Lambda$ asymptotics and the reassembled integrand.

See **[`nlo/README.md`](nlo/README.md)**.

```
wolframscript -file nlo/nlo_pplus_demo.wls
```
