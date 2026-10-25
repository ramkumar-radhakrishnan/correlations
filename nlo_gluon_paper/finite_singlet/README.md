# The log² term, and the final finite result in colour-singlet form

[`Finite_terms_singlet.pdf`](Finite_terms_singlet.pdf) (8 pages; source `Finite_terms_singlet.tex`).

## The log² term

- It appears once in `Finite_terms_1.tex`, in Group I Row II part 1: `log²(1−k⁺/∨) + log(1−k⁺/∨) log(∨/k⁺)`. Nothing cancels it.
- It is the finite part of |Ψᴬ_gρ⟩ + |Ψᴮ_gρ⟩. It comes from the lower limit ξ ≥ k⁺/∨ of the ξ integral of |Ψᴮ⟩ (the parent gluon below ∨). This was re-derived analytically and checked numerically (`log2_origin.py`).
- The same limit gives the log(1−a) of the pole, i.e. the uncancelled pole Q1. The two come and go together:
  - if the regular part of |Ψᴮ⟩ is integrated from ξ = 0 (AABKL), both vanish and only 67/9 − π²/3 is left.
- Written compactly: −log²(1−a) − L_∨ log(1−a) = −log(1−k⁺/∨) log(∨/k⁺ − 1). This is ≈ a log(1/a), i.e. small for k⁺ ≪ ∨.
- It is not running coupling. `../finite_v2/` (FS (4.3), (11.1)) now lists it separately as V_∨.

## Singlet projection

The projection is ⟨ρᵃ(u)ρᵇ(v)⟩ = δᵃᵇ⟨ρᶠ(u)ρᶠ(v)⟩/(N_c²−1). The resulting adjoint traces use `../singlet/trace.py`. Differences of Wilson lines are combined after the projection.

The following pieces are projected:
- LO and the LO-like terms;
- the two-ρ rows (G1 V 2a, 2b; G2 I(2); G3 III(2));
- the two-ρ parts of every three-ρ row, merge by merge (G1 II, III, V(1); the K_IV and K_C families).

The Wilson lines use the index placement for a general U:
- Group I: as written.
- Groups II and III: the M4 placements of `../groups23/fixed23.py`.

| File | Contents |
|---|---|
| `log2_origin.py` | Ψᴬ + Ψᴮ vs the Row II part 1 bracket; the restricted and unrestricted ξ integrals → `log2_origin.json` |
| `gsrc.py` | general-U colour operators of the Group II/III three-ρ rows; check vs `fixed23.py` → `gsrc.json` |
| `project.py` | projection, trace forms, differences, and all colour checks → `project.pkl`, `project.json` |
| `srender.py`, `make_note.py`, `static.tex` | the note |

```
python3 log2_origin.py
python3 gsrc.py
python3 project.py
python3 make_note.py
pdflatex Finite_terms_singlet.tex   # three times
```
