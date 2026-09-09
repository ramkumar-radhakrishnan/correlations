# Verification of the divergence-extraction steps

Every integral and every algebraic assembly used in
`../trijet_step_by_step.pdf` is checked here. Requires `numpy`, `scipy`, `sympy`.

```
python3 checks1.py   # eikonal transverse integral, xi integral, DGLAP integral
python3 checks2.py   # the lifetime-ordering "triangle" integral
python3 checks3.py   # symbolic assembly  (B.7)+(B.23)+(B.24) = (B.25)
python3 checks4.py   # the CGC dipole-kernel identities (3.27) and K_dip
```

Expected results:

| check | statement | status |
|---|---|---|
| 1a | `∫_a^∞ J₀(x)dx/x = ln(2/a) − γ_E + O(a²)` | agrees to `O(a²)` |
| 1b | `∫d²C/(2π)² e^{−iC·Δ}/C² Θ(C²>Λ²) = (1/2π)ln(c₀/ΛΔ)` — CSSV (B.9) | agrees |
| 1c | `∫_ε^1 dξ/ξ (A+2lnξ) = A ln(1/ε) − ln²ε` | exact |
| 1d | `∫_ε^1 dξ (1+(1−ξ)²)/ξ = 2ln(1/ε) − 3/2` — fixes the `1/ε` coefficient | exact |
| 2 | `∫d²z (r_bb'²/r_zb²r_zb'²) ln(min·Q_f²)Θ(min·Q_f²−1) = π ln²(r_bb'²Q_f²)` — CSSV (4.51) | agrees, **no residual single log** |
| 3 | `V + C + O − (B.25) = 0` symbolically; `1/ε`, `ln²z₀`, `ln z₀ ln z_i`, `ln μ²` all cancel; `ln R²` coefficient `−3/2`; residual `z₀`-dependence is exactly `ln z₀ · ln(r_xx'²r_yy'²/(r_xy²r_x'y'²))` | **exactly 0** |
| 4 | identity CSSV (3.27), and `(1/π)∫d²z K_dip = ln(r_xy²r_x'y'²/(r_xx'²r_yy'²))` | agrees |


## Verification for `trijet_divergence_separation.pdf`

These operate directly on the user-supplied `Trijet_LO.pdf` expressions (1.1)-(1.24).

```
python3 checkW.py    # soft (xi -> 0) limits of the Phi.Phi* spinor structures
python3 checkA.py    # exact xi-dependence and the subtracted combinations
python3 checkW2.py   # the four-term W algebra, and UV finiteness at z -> y
python3 checkK.py    # the eikonal-kernel factorization identity
```

| check | statement | status |
|---|---|---|
| W | all `eps.eps`, Reg×Inst, Inst×Inst structures vanish as `xi -> 0` (as `xi`, `xi`, `xi²`) | confirmed |
| A | `A_qbar(z,0) = A_q(z,0) = A_qqbar(z,0) = 16 z² zbar² [z²+zbar²]` — one common weight | exact |
| A | `A(z,xi) - A(z,0) = O(xi)` in all three channels, so the subtraction is finite | exact |
| W2 | every `S` and every `+1` cancels in the four-term `W` combinations | exact |
| W2 | `W_qbar -> 0` at `z -> y` and at `zbar -> ybar`; `W_qqbar -> 0` at `z -> x` — UV finiteness | identically 0 |
| K | `K_qq + K_qbqb - K_qqb - K_qbq = A·Abar` (no coincidence limit needed) | exact |
| K | coincidence limit `= (x-y)²/[(z-x)²(z-y)²]` — the BK dipole kernel | exact |
