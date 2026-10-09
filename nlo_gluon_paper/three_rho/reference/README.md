# Reference rows (region T, leading log)

Independent construction of the log(∨/Λ) rows used as a cross-check in
`ThreeRho_cancellation.pdf` (Sec. 6.2 there).

The rows are built directly from the leading-log light-cone wave function, with the operator order kept:
- two-gluon amplitude pieces for the real rows (Groups II, III);
- one-gluon pieces for the virtual rows (Group I).

They were checked earlier against an exact numerical evaluation of the operator products.

| File | Contents |
|---|---|
| `rowsT_sym.py` | the rows: `real_rows()`, `virt_rows()` |
| `sym.py` | symbolic term class (Wilson lines, f's, ordered charges, kernel) |
| `weyl.py` | Weyl decomposition of ordered charge products |

The `__main__` block of `rowsT_sym.py` compares with the numerical model `rowsT.py`, which is not included here.
