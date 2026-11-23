# Verification for `DIS_NLO.pdf` and for `../../paper/`

Checks of the massive trijet note `DIS_NLO.pdf` (Secs. 1-4) and of every formula used in the JHEP-format
paper in [`../../paper/`](../../paper/). Sections 1-3 of the note are covered by the scripts in
[`../massive/`](../massive/) (run from that directory); the two scripts here cover Sec. 4 and the limits.

```
python3 checkJ_momentum_space.py   # Sec. 4 (momentum space, elastic part) + collinear and soft weights (~5 s)
python3 checkK_limits.py           # numerical limits of the radial functions                          (~10 s)
```

They need `sympy`, `numpy` and `scipy`.

| script | statement | result |
|---|---|---|
| `../massive/checkF` | full Feynman amplitude with massive Dirac spinors (both emitters, regular + instantaneous) `= q^+ sqrt(zeta eta) x` the two-component expressions (2.3), (2.24), (2.32), (2.46) | exact, one common constant |
| `../massive/checkG` | radial functions (2.11)-(2.14) vs. direct Fourier transform; (2.25) `G_E = E_f K1(E_f D)/D` | agree to `1e-10` |
| `../massive/checkH` | helicity sums (3.3)-(3.5), (3.24)-(3.26), (3.49)-(3.51), symbolic in `zeta, xi, m` | every coefficient of the note reproduced |
| checkJ (1) | changes of variables (4.3)-(4.7) and shifted momenta (4.14) | exact |
| checkJ (2) | `f1 - f2` of (4.10); identity (1.9) for `E_f^2` | exact |
| checkJ (3) | momentum-space sums (4.19), (4.21), (4.23) from the amplitudes (4.17), random complex `T`, three `(zeta, xi, m)` sets incl. `m = 0` | exact |
| checkJ (4) | after x after term: `sum |V_gamma V_g|^2 = 8 [(zeta^2+zb^2) k^2 + m^2] zb^2 [(1+(1-z)^2) p^2 + z^4 m^2]`, i.e. LO factor times the massive quasi-collinear splitting function `(1+(1-z)^2)/z - 2z(1-z)m^2/(p^2+z^2m^2)` | exact |
| checkJ (5) | soft limit: the three channels carry the same weight for `[zeta^2+zb^2] eps^2 K1 K1` and for `m^2 K0 K0` | `[1,1,1]`, `[1,1,1]` |
| checkK (1) | `m -> 0`: `G11 -> Q K1(QD)/(D Y^2)` quadratically in `m`; flip terms vanish | confirmed |
| checkK (2) | `xi -> 0`: `Y^2 G11 -> eps K1(eps R)/(zeta zb R)`, `Y^2 G01 -> K0(eps R)/(zeta zb)` | confirmed |
| checkK (3) | `Q -> 0` (photoproduction) finite for `m > 0` | confirmed |
| checkK (4) | dead cone in coordinate space: `om K1(om Y)/Y -> 1/Y^2` only for `om Y << 1` | confirmed |
