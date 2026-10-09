# Three-ρ terms ∝ log(∨/Λ): two-ρ terms and cancellation

This folder has two notes:
1. [`ThreeRho_to_TwoRho.pdf`](ThreeRho_to_TwoRho.pdf): the two-ρ terms obtained by reordering the three-ρ terms.
2. [`ThreeRho_cancellation.pdf`](ThreeRho_cancellation.pdf): the check that the symmetrized three-ρ terms cancel among each other (see below).

## Note 1: two-ρ terms

[`ThreeRho_to_TwoRho.pdf`](ThreeRho_to_TwoRho.pdf) is a 36-page note. It takes every three-ρ term proportional to log(∨/Λ) and writes it as its fully symmetrized part plus commutators. The commutators produce the two-ρ terms, which are listed with every step.

The three-ρ terms come from two places:
- the remainders of reordering the four-ρ terms (two and three Wilson lines);
- the genuine three-ρ terms of Groups I, II and III, from [`input/Evolution_three_rho.tex`](input/Evolution_three_rho.tex).

## What the note contains

- **Secs. 2–3:** conventions, and the reordering identity derived step by step:
  ABC = S(ABC) + ¼({[A,B],C} + {[A,C],B} + {[B,C],A}) + (one-ρ terms),
  plus its forms A{B,C} and {B,C}A.
- **Secs. 4–7:** for every source and every sub-term:
  - the commutators and δ-functions;
  - each two-ρ piece before and after the colour algebra;
  - the collected result.
- **Sec. 8:** all two-ρ terms in one place, with the kernels written out. The complex conjugate is included as (P ↔ P′) for Groups I and III.
- **App. A:** the coordinate-space kernel of Row II of Group I.
- **App. B:** numerical checks.

## Reproducing

```
python3 build_doc.py          # regenerates ThreeRho_to_TwoRho.tex from sources.py
pdflatex ThreeRho_to_TwoRho.tex   # twice
python3 process.py            # colour-algebra check of all 144 pieces (takes several minutes)
python3 check_operator.py     # whole reordering with explicit charge operators
python3 check_identity.py     # the reordering identity with random matrices
```

| File | Role |
|---|---|
| `sources.py` | Transcription of all three-ρ terms |
| `engine.py` | Commutators, δ-functions, colour identities |
| `make_tex.py`, `static_tex.py`, `build_doc.py` | LaTeX generation |

Requires `numpy` and `sympy`.

## Note 2: do the symmetrized three-ρ terms cancel?

[`ThreeRho_cancellation.pdf`](ThreeRho_cancellation.pdf) is a 19-page note. It adds up the symmetrized parts S(ρρρ) of all three-ρ terms:
- the four-ρ remainders;
- Group I and Group III, each with its c.c.;
- Group II.

How the check works:
- Every term is written in common integration variables.
- The terms fall into six kernel classes.
- Within a class, the colour structures are grouped into networks, and the coefficients of each network are summed exactly.

Result:
- **As written:** 48 networks; 28 cancel and 20 do not. Every left-over traces back to Group I.
  - Row II: the α part of the two-charge kernel in the U(w)𝔸⁽³⁾ term has no partner.
  - Row III, part 2: the edge-kernel Φ terms have no partner.
  - A term 𝒟 is missing: a soft gluon exchanged between the measured gluon and a rotated charge.
- **Correct as written:** the four-ρ remainders, Group II and Group III. They agree row by row with an independent leading-log computation (`reference/`).
- **After two changes in Group I,** all networks cancel exactly:
  1. Row II keeps only 𝔉_β in the U(w)𝔸⁽³⁾ term;
  2. Row III part 2 is replaced by 𝒟.

  A pointwise numerical test confirms this: the relative residual is 10⁻¹⁶ after the changes, and 5–40% as written.
- **Convention:** everything uses U^{bc} = U^{cb}, as the notes do (Sec. 10 of the note).

```
python3 crosscheck.py         # independent checks vs the reference -> crosscheck.json (~20 s)
python3 pointwise_check.py    # pointwise numerical test -> pointwise.json
python3 make_cancel_tex.py    # writes ThreeRho_cancellation.tex from cancel_static.tex
pdflatex ThreeRho_cancellation.tex   # twice
```

| File | Role |
|---|---|
| `cancel_data.py` | Entries, renaming to common variables, kernel classes, colour networks, exact sums, the term 𝒟 |
| `cancel3.py`, `classes3.py`, `classes4.py`, `classes6.py` | Term lists and class bookkeeping (numerical format) |
| `corrections.py` | The two corrections (and the equivalent commutator form of the first) |
| `compare_ref.py`, `fourrho.py`, `four_check.py`, `g1split.py`, `comm.py`, `crosscheck.py` | Comparison with the reference rows in `reference/` |
| `pointwise.py`, `pointwise_check.py` | Pointwise test |
| `check_pbw4.py` | Check of the three-ρ part of a product of four charges |
| `make_cancel_tex.py`, `cancel_static.tex` | LaTeX generation |
