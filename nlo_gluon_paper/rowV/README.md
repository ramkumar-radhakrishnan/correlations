# Group I, Row V: derivation check, UV divergences, part 2b in dimensional regularization

[`RowV_UV.pdf`](RowV_UV.pdf) (9 pages) checks `Group_V.tex` and treats the UV divergence of part 2b in $\overline{\rm MS}$.

## Result

- **The derivation is correct.** These steps were checked:
  - the prefactor;
  - the contraction of the ℂ₁ vertex with the 𝔹₂ bracket (all nine products);
  - the w integration;
  - the colour algebra;
  - the split into 2a and 2b;
  - the plus prescription of 2a.

  Small notational points:
  - S^{ad} → U^{ad} in the identity f f U U = N_c U (it holds for any adjoint U);
  - after the w integration, the phase e^{-ik·(z−x)p⁺/k⁺} belongs inside ∫dp⁺.
- **Part 2a has no UV divergence.** Its colour operator vanishes linearly at x = z, in the plus-prescribed integral and in the endpoint terms alike.
- **Part 2b is logarithmically UV divergent at x → z.**
  - The coefficient is the splitting function times LO: Φ̂(ξ) = ξ(1−ξ) + ξ/(1−ξ) + (1−ξ)/ξ.
  - Φ̂ does not depend on d_⊥ when the angular average is done in d_⊥ dimensions.
  - After the plus prescription, the endpoint (JIMWLK) terms carry 2 log((k⁺−Λ)/Λ) of this divergence.
  - The remainder carries ∫₀¹[ξ(1−ξ) − 2] = −11/6: the β₀ part.
- **MS-bar result** (per ordering, with the factor 2 of the row):
  - d³N_2b|_UV = { −(α_s β₀/4π)[1/ε + log(μ̄²/k²)] − (α_s N_c/4π)(67/9) } × LO;
  - the rest has the plus-prescription form of 2a (Eq. 5.9 of the note).
- **Running coupling.** Per ordering, Group I Row II (part 1) has the opposite β₀ term.
  - The sum has no pole, no μ dependence and no 67/9: (α_s β₀/4π) log((x−w)² k² e^{2γ_E}/4) × LO.
  - This is equivalent to α_s → α_s · α_s(c₀/(x−w)²)/α_s(k²).
  - The coupling renormalization of LO itself must then come from the self-energy of the measured gluon after the target (diagonal 𝔹̄₃, "G1-I′"). This is open.

## Reproducing

```
python3 rowV_checks.py     # all numerical and symbolic checks -> rowV_checks.json (~30 s)
python3 make_rowV_tex.py   # writes RowV_UV.tex from RowV_static.tex
pdflatex RowV_UV.tex       # twice
```

`rowV_checks.py` uses `../three_rho/numerics.py` (SU(3) structure constants, random adjoint Wilson lines). Requires `numpy`, `scipy`, `sympy`.

## Part 2b in dimensional regularization, without plus prescription

[`Part2b_MSbar.pdf`](Part2b_MSbar.pdf) (7 pages; source `Part2b_MSbar.tex`) splits part 2b, as written, into a UV term and a UV-finite rest.

- **Counterterm.** The short-distance limit of the integrand, [p⁺(k⁺−p⁺)/k⁺² + p⁺/(k⁺−p⁺) + (k⁺−p⁺)/p⁺] 𝒦(x'−w')·𝒦(y−z) N_c[U(y)−U(z)], times K(z,z,x) − K(y,z,x) = (y−x)·(y−z)/((y−x)²(x−z)²).
  - Its angular average is exactly θ(|y−z|−|x−z|)/(x−z)².
- **UV term** (d_⊥ = 2−ε, MS-bar).
  - The x integral gives [2/ε + 2γ_E + log((y−z)²μ²/4)] = (1/π)∫_x[K(z,z,x) − K(y,z,x)], the same form as G1-II after F6.
  - The p⁺ integral is exact: 2log((k⁺−Λ)/Λ) − 11/6 + 4Λ/k⁺ − Λ²/k⁺² + 2Λ³/(3k⁺³).
  - The prefactor is (α_s N_c/2π) × LO.
  - The β₀ part cancels the β₀ part of G1-II (part 1) per ordering, pole included.
- **Rest.** The original integrand minus the counterterm. It is UV finite at every p⁺ and is evaluated with d_⊥ = 2.

```
python3 part2b_checks.py     # -> part2b_checks.json (~1 min)
python3 make_part2b_tex.py   # writes Part2b_MSbar.tex from Part2b_static.tex
pdflatex Part2b_MSbar.tex    # twice
```
