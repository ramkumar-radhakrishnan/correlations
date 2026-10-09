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

## `nlo_gluon_paper/` — single inclusive gluon production at NLO at mid rapidity

Draft paper on the NLO single inclusive gluon spectrum in dilute–dense scattering, in the light-cone
wave function approach. It covers:

- the NLO wave functions and the coefficients of the evolution operator $\Omega$;
- the cross section at order $g^4$;
- the leading-log factorization $\log(\vee/k^+)\,[\mathrm{BFKL}\otimes\mathrm{LO}]+\log(k^+/\Lambda)\,[\mathrm{JIMWLK}\otimes\mathrm{LO}]$;
- the finite part, row by row;
- running coupling and DGLAP.

Start with **[`nlo_gluon_paper/README.md`](nlo_gluon_paper/README.md)**. The PDF is
`nlo_gluon_paper/NLO_gluon_midrapidity.pdf`.

```
cd nlo_gluon_paper
pdflatex NLO_gluon_midrapidity.tex   # three times
```
