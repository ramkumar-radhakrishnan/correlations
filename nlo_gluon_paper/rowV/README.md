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

## Plus prescription of the UV-finite part of 2b

[`PlusPrescription_2b.pdf`](PlusPrescription_2b.pdf) (5 pages; source `PlusPrescription_2b.tex`) applies the plus prescription to d³N_{2b finite NLO}/d³k (part 2b minus its counterterm), as written in [`input_2b_finite.tex`](input_2b_finite.tex).

- **Exact, on Λ < p⁺ < k⁺−Λ.** No limit Λ → 0 is taken.
  - ∫ (k⁺/p⁺) f = ∫ [k⁺/p⁺]₊ f + k⁺ f(0) ln((k⁺−Λ)/Λ).
  - ∫ (k⁺/(k⁺−p⁺)) f = ∫ [k⁺/(k⁺−p⁺)]₊ f + k⁺ f(k⁺) ln((k⁺−Λ)/Λ).
  - ln((k⁺−Λ)/Λ) = log(∨/Λ) − log(∨/k⁺) + log(1−Λ/k⁺).
- **Singular coefficients.**
  - −(k⁺−p⁺)/p⁺ = 1 − k⁺/p⁺.
  - −p⁺/(k⁺−p⁺) = 1 − k⁺/(k⁺−p⁺).
  - Counterterm: p⁺(k⁺−p⁺)/k⁺² + p⁺/(k⁺−p⁺) + (k⁺−p⁺)/p⁺ = p⁺(k⁺−p⁺)/k⁺² − 2 + k⁺/p⁺ + k⁺/(k⁺−p⁺). Its 1/p⁺ and 1/(k⁺−p⁺) parts go entirely into the endpoint terms, so only p⁺(k⁺−p⁺)/k⁺² − 2 stays under the p⁺ integral.
- **Endpoint terms**, with prefactor (1/(2π)³)(g⁴/16π⁵)(1/k⁺) ln((k⁺−Λ)/Λ), in the notation K(a,b,c) = (a−c)·(b−c)/((a−c)²(b−c)²):
  - p⁺ → 0: e^{−ik(w′−z)} 𝒦(x′−w′)·𝒦(y−z) [(K(z,z,x) − ½K(y,z,x)) 𝒪_x − (K(z,z,x) − K(y,z,x)) 𝒪_z].
  - p⁺ → k⁺: e^{−ik(w′−x)} 𝒦(x′−w′)·𝒦(y−x) (K(x,x,z) − ½K(y,x,z)) 𝒪_x − e^{−ik(w′−z)} 𝒦(x′−w′)·𝒦(y−z) (K(z,z,x) − K(y,z,x)) 𝒪_z.
  - Combined, after x ↔ z in the second line: e^{−ik(w′−z)} 𝒦(x′−w′)·𝒦(y−z) N_c[U(x′)−U(w′)] {(1/(x−z)²)[U(z)−U(x)] + K(y,z,x)[U(y) − (3/2)U(z) + ½U(x)]} ρρ.
- **UV finite piece by piece.** The plus-prescribed integral and both endpoint terms are each finite at x → z.
- **Large distances.** The separated pieces fall like 1/(2r²) at large |x−z|. Each is logarithmically divergent there, and the divergences cancel in the sum. So transverse integrals must be done on the sum, or with the same large-distance cutoff in both pieces.
- **Checks.** Plus + endpoint reproduces the original p⁺ integral to relative 10⁻¹⁵ at random transverse points. The endpoint coefficients agree with the limits of the integrand.

```
python3 plus2b_checks.py     # exactness, endpoint limits, UV -> plus2b_checks.json
python3 plus2b_ir.py         # large-distance behaviour -> plus2b_ir.json
python3 make_plus2b_tex.py   # writes PlusPrescription_2b.tex from plus2b_static.tex
pdflatex PlusPrescription_2b.tex    # twice
```


## Part 2b: UV term + finite term, d⊥ kept general

[`Part2b_UV_finite.pdf`](Part2b_UV_finite.pdf) (4 pages; `make_uvfin.py` from `uvfin_static.tex`, `part2b_dperp_checks.json` and `part2b_checks.json`).

- **Short-distance coefficient.** At x → z the coefficient is (d⊥/d)ξ(1−ξ) + ξ/(1−ξ) + (1−ξ)/ξ, where d is the dimension of the x integral. Checked for d⊥ ≠ d.
- **UV term.**
  [2log(∨/Λ) − 2log(∨/k⁺) + 2log(1−Λ/k⁺) − 2 + 4Λ/k⁺][2/ε + 2γ_E + log((y−z)²μ²/4)] + (d⊥/2)[1/6 − Λ²/k⁺² + 2Λ³/3k⁺³][2/ε + 2γ_E + log((y−z)²μ²/4) + 1].
  With d⊥ = 2−ε it reduces to the −11/6 form of `Part2b_MSbar.pdf`.
- **Finite term.** It is 2b minus the counterterm. The counterterm has (d⊥/2)ξ(1−ξ) + ξ/(1−ξ) + (1−ξ)/ξ, kernel (y−x)·(y−z)/((y−x)²(x−z)²), colour N_c[U(y) − U(z)], and no phase.

```
python3 part2b_dperp_checks.py   # -> part2b_dperp_checks.json
python3 make_uvfin.py            # Part2b_UV_finite.tex
```


## Part 2b: UV term with the phase kept

[`Part2b_UV_phase.pdf`](Part2b_UV_phase.pdf) (3 pages; `make_uvphase.py` from `uvphase_static.tex` and `part2b_phase_checks.json`).

- **Subtracted term.** The UV term keeps e^{−ik·(z−x)p⁺/k⁺}. The subtracted term is the leading 1/(x−z)² of all six brackets, with denominator and colour at x = z. The kernel is plain 1/(x−z)², because the phase makes large distances converge.
- **x integral.** π[2/ε − log(p⁺² k²/(k⁺² μ²))] (MS-bar, d = 2−ε). The pole is the same as before.
- **After the exact p⁺ integral.** [2log((1−λ)/λ) − 2 + 4λ][2/ε − log(k²/μ²)] − 2A(λ) + (d⊥/2)([1/6 − λ² + 2λ³/3][2/ε − log(k²/μ²) + 1] − 2B(λ)). A and B are given in closed form; A contains log²(Λ/k⁺).
- **Relation to the split without the phase.** The difference is ∫d²x[e^{iq·r}/r² − K-kernel] = −π log(q²(y−z)²e^{2γ_E}/4), which is finite. The total is unchanged.

```
python3 part2b_phase_checks.py   # -> part2b_phase_checks.json
python3 make_uvphase.py          # Part2b_UV_phase.tex
```


## Part 2 (2a, 2b UV, 2b finite with plus prescription): check of your derivation

[`Part2_check.pdf`](Part2_check.pdf) (4 pages; `make_part2check.py` from `part2check_static.tex` and `part2_user_check.json`).

- **Correct:**
  - the colour identity and the 2a/2b split;
  - 2a with the plus prescription, up to O(Λ/k⁺);
  - 2b UV;
  - 2b finite before "Simplifying".
- **E1.** The B₂ pieces of 2b finite (the [k⁺/q⁺]₊ term, the −B₂ term and the q⁺ → 0 endpoint) must have colour N_c[U(y) − U(x)], not U(z).
- **E2.** The p⁺ → 0 endpoint of B₄ in 2b finite has phase e^{−ik·(w′−z)}, not e^{−ik·(w′−x)}, and colour U(x).
- **E3.** The p⁺-integrated counterterm has the bracket [2L − 11/6] (+4λ − λ² + 2λ³/3), not [11/6 − L/3]. Here L = log((k⁺−Λ)/Λ).
- **Minor:** Λ → 0 is applied inconsistently in 2b UV; the plus definition needs the factor k⁺; the phase belongs inside the integrals; and the integrated counterterm is UV finite only in the sum. The recommended grouping avoids the last point.
- **Numerical check.** At random transverse points, the user's 2b finite is off by O(1). With E1–E3 fixed it agrees to O(Λ/k⁺). The corrected form with exact limits agrees to 10⁻¹³.

```
python3 part2_user_check.py   # -> part2_user_check.json
python3 make_part2check.py    # Part2_check.tex
```
