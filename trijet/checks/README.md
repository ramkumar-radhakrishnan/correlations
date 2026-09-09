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
