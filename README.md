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

## `nlo-multiplicity/` — gluon multiplicity moments from the NLO CGC evolution operator

Research analysis of the programme announced at the end of the Introduction of
[arXiv:2607.18373](https://arxiv.org/abs/2607.18373): moments of the gluon multiplicity operator
and particle-number fluctuations at order $g^4$. Formalism (Parts I–IV), model calculations and
proposed numerics (Parts V–VI), literature survey (VII), and a critical publication assessment
and roadmap (VIII–IX), plus an addendum on the source paper's JHEP prospects.

Start with **[`nlo-multiplicity/README.md`](nlo-multiplicity/README.md)**, or the typeset PDF
`nlo-multiplicity/NLO_CGC_multiplicity_analysis.pdf` (build source in `nlo-multiplicity/pdf/`).

```
pip install reportlab
python3 nlo-multiplicity/pdf/build_pdf.py
```

## `review/` — error review of the draft manuscript

[`review/draft-review-omega-g2.md`](review/draft-review-omega-g2.md) — equation-by-equation
review of the `Omega at g^2` draft: errors verified by re-derivation, internal inconsistencies,
conceptual gaps, and typography.
