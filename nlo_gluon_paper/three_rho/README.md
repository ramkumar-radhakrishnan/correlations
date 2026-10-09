# Two-ρ terms from the three-ρ terms ∝ log(∨/Λ)

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
