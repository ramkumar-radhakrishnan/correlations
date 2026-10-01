# Verification for `../../massive/Trijet_simplified_massive.tex`

These scripts check the massive rewrite of `Trijet_simplified.tex` (LO trijet, massive quarks).
Run them from this directory, because `checkE` and `checkF` load `lb.py` from here. They need `sympy`, `numpy` and `scipy`.

```
python3 checkE_vertices.py        # massive two-component vertices fitted to Dirac spinors   (~1 min)
python3 checkF_full_amplitude.py  # full amplitude, both emitters, reg + inst, vs Sec. 2     (~3 min)
python3 checkG_fourier.py         # the new radial functions G, A, G_E                       (~1 min)
python3 checkH_helicity_sums.py   # every helicity/polarisation sum of Sec. 3                (~2 min)
python3 checkI_identities.py      # the Reg x Inst contractions quoted in Sec. 3
```

| script | statement | result |
|---|---|---|
| checkE | `sum_lam' [ubar eps v][vbar eps v] = c sum_lam (phi k + m mu)(tau p + xi omega mu)`, with `mu^i_{ab} = delta_{b,-a}(2a delta^{i1} - i delta^{i2})`; the quark line is the same with `xi^2/(zeta+xi)` | unique fit; `m^2` terms then match with no free parameter left; same constant `c` in both channels |
| checkF | Full Feynman amplitude (both emitters, including the instantaneous `gamma^+` pieces) `= q^+ sqrt(zeta eta) x` the two-component expressions of Sec. 2 | exact: 16 helicity/polarisation combinations, 4 mass/fraction sets (incl. `m = 0`), 2 momenta each |
| checkG | Feynman-parameter form of `G^(ab)` vs. exact p-integral plus Hankel transform | agree to `1e-13` |
| checkG | `G^(11)(m=0) = Q K1(QD)/(D Y^2)`; instantaneous FT `= E K1(E D)/D`, `E^2 = Q^2 + m^2(1/zeta + 1/eta)` | exact to `1e-15` |
| checkG | `Y -> 0`: `Y^2 G^(11)`, `Y^2 G^(01)` tend to the after-shock-wave limits; `G^(10) - A^(10)`, `G^(00) - A^(00)` finite | confirmed |
| checkH | all `H`, `H~`, `H^x` of Sec. 3, symbolic in `zeta, xi, m`; massless Reg x Reg of the original note reproduced | exact |
| checkI | `sum phi tau (delta + 2i eps lam)^* = 4 zeta eta (delta - 2i eps lam)` and the quark and interference analogues | exact |

Among other things, checkH shows that the massless Regular x Instantaneous and Instantaneous x Instantaneous
expressions of the original note (L289, L293, L385, L389, L474) had the wrong tensor structure or the wrong
coefficient. The corrected forms are in the rewritten note and listed in its Appendix B.
