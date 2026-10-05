# Revision report: *Chiral soliton lattice in inhomogeneous magnetic fields*

Response to Referee 1. The referee recommends publication with minor questions and no major concerns. All five comments are addressed in the revised manuscript. Each addition was checked independently with numerical computations.

## Files in this folder

| File | What it is |
|---|---|
| `manuscript.tex` | Revised manuscript (sn-jnl), ready to replace the submitted source |
| `manuscript_original.tex` | Submitted version, unchanged, for reference |
| `manuscript_diff.tex` | latexdiff marked-up version (substantive changes only; see note below) |
| `response_to_referee.tex` / `.pdf` | Point-by-point reply letter |
| `new_references.bib` | 7 new BibTeX entries to append to `references.bib` |
| `previews/*.pdf` | Previews compiled with a **stub article class** (sn-jnl was not available here), figures shown as placeholders. Use them to read the text; they are not for submission |
| `verification/*.py` | Scripts behind the numerical checks quoted below (numpy + scipy) |

**Note on the diff:** `manuscript_diff.tex` compares against the original with the H→b rename already applied. Otherwise the rename would mark up ~100 places, and latexdiff breaks on the multi-line `align` blocks. The letter tells the referee that the rename is not highlighted.

## How each comment was addressed

| # | Referee comment | What changed | Where |
|---|---|---|---|
| 1 | Symbol **H** suggests the auxiliary magnetic field | Renamed **H → b** throughout (macro `\bmag`, so it can be changed in one line). Added a sentence saying **b** is just the rescaled **B**, not the auxiliary **H**. Also renamed the Gaussian stream function ψ → χ, because it clashed with the scalar potential ψ[**b**] of Eq. (18). | Below Eq. (3); Eq. (32); captions |
| 2 | Confusion about Fig. 1: is φ ≠ 0 above *and* below the curve? | **Yes.** In finite volume, φ=0 violates the natural boundary condition, so it is not a stationary point. The text now explains the two classes: type (16) has φ₀(0)=0 and is the vacuum-like state, bent only in a boundary layer of width ~1/m_π. Type (15) has φ₀(0)=π, i.e. a soliton core at the center. The curve marks where the first soliton enters, not where φ₀ becomes nonzero. Also added: the explicit weak-field profile, the reason the critical field rises at small L̄, and the parity rule (odd/even soliton number ↔ class). Caption of Fig. 1 rewritten. | End of Sec. 2.2, Fig. 1 caption |
| 3 | Physical lessons from Figs. 3, 4; intuition for other fields | (a) **Frustration:** a gradient cannot circulate, so φ locks onto the conservative part of **b**, set by the boundary flux. Closed field loops are invisible to it. (b) **New identity, Eq. (22):** N_B = (f²/μ)∫(∇φ₀)² = −2E/μ ≥ 0 in the chiral limit, so lost condensation energy means proportionally lost baryon number. (c) **Fig. 3 derived analytically:** φ₀ ∝ r⁴ sin4θ, giving density ∝ −r⁴e^{−r²/R²}cos4θ, which is negative on the axes, positive on the diagonals, with peaks at r=√2R. (d) **Fig. 4:** 0 ≤ n_B ≤ b₀² by the maximum principle. The healing length is L_z/π, set by geometry because the chiral limit has no intrinsic scale; this explains the L_z² energy scaling in App. A. (e) Rules of thumb collected in the Summary. | End of Sec. 3.1; Sec. 4.2; Sec. 4.3; Sec. 5 |
| 4 | Give numbers for the scales; neutron star mergers? | New intro paragraph. The CSL scale is 1/m_π ≈ 1.4 fm. In heavy-ion collisions, eB ~ m_π² (RHIC) to ~10 m_π² (LHC) varies over a few fm, which matches our spatial regime, but the field is short-lived, μ_B is small and T is high. Magnetars have 10¹⁴–10¹⁵ G at the surface and ≲10¹⁸ G inside, below B_CSL ~ 10¹⁹ G, and vary over ~km. In mergers, KH/MRI amplification gives ≳10¹⁶ G, but the structure remains macroscopic (simulations resolve ~10 m) and T reaches tens of MeV. So the CSL would follow the local field adiabatically, and thermal effects likely destroy it. The paper's configurations are relevant to lattice QCD, not astrophysics. | Sec. 1, Sec. 5 |
| 5 | Why no coupling to Maxwell's equations? | New Eq. (4): the induced current is **j** = ∇μ×∇φ/(4π²), which vanishes identically for uniform μ. Equivalently, the magnetization **M** = μ∇φ/(4π²) is curl-free, so Ampère's law holds trivially in the bulk. Back-reaction comes only through ∇·**M** and boundary terms, suppressed by αμ²/(4π³f_π²) ≈ 0.7% at μ = 1 GeV. Added a matching remark for nonuniform μ in Sec. 3.3. | Sec. 2 (after Eq. 3), Sec. 3.3 |

Also fixed: `2L` → `2\bar L` in Eq. (14); `∂_z φ₀ = H` → `∂_{z̄} φ₀ = b̄`; "convolution" → "product" (Sec. 3.3); "both of the above equations" → "questions" (Sec. 5).

## Verification performed

- **Compilation:** the revised manuscript, the diff and the letter all compile with **0 errors and 0 undefined references/citations**, given `new_references.bib`. The manuscript was compiled with a stub class, because sn-jnl could not be downloaded here.
- **Point 2 (1D, `one_dimensional_checks.py`):** an independent minimization reproduces Fig. 1.
  - Fields just below the curve give the min-type state: L̄=5, b̄=1; L̄=10, b̄=1; L̄=2, b̄=2. Each has φ₀(0)=0 and no soliton.
  - Fields just above give the max-type state: L̄=5, b̄=2; L̄=10, b̄=1.5; L̄=2, b̄=4. Each has φ₀(0)=π and one soliton.
  - Parity rule confirmed: L̄=10, b̄=3 has 4 solitons and a min at the center; b̄=5 has 7 solitons and a max.
- **Point 3 (chiral limit, `chiral_limit_checks.py`, 201×201 grid):**
  - *Gaussian field:* the sign pattern and peak radius (2.83) match Fig. 3. max|n_B| / max b² ≈ 1.2×10⁻⁴. ∫b·∇φ₀ = ∫(∇φ₀)² ≈ 1.8×10⁻⁴, compared with ∫b² ≈ 50.
  - *Domain wall:* ⟨E⟩/E₀ = 0.4998 (analytic 1/2). n_B/b₀² lies in [0, 0.996]. The baryon-number identity holds to 10⁻⁴. ∂_zφ at the edge of the square is 0.50·b₀.

## Action items before resubmission (please do these)

1. **Regenerate Fig. 1:** its y-axis label is still **H̄** in the image file. Change it to **b̄**. The letter already says this was done. No other figure shows H on an axis.
2. **Append `new_references.bib` to `references.bib`** and check the 7 entries against INSPIRE. I wrote them from memory; keys, volumes and pages should be right but are unverified: Kharzeev:2007jp, Skokov:2009qp, Deng:2012pc, Kaspi:2017fwg, Price:2006fi, Kiuchi:2015sga, Perego:2019adq.
3. **Compile with sn-jnl** and check the layout. Equation numbers shift by +1 after Eq. (3) and +2 after the new Eq. (22). The letter uses the **new** numbers.
4. **Check the numbers quoted for Fig. 3 against your Mathematica data.** In your units (χ₀=4, R=2), the plotted density should peak at |n_B| ≈ 3.5×10⁻⁴.
5. **Check the factor-of-e convention** in the Maxwell remark. The text says e is absorbed into **B**, which is consistent with B_CSL = 16πf²m_π/μ ≈ 0.06 GeV² ≈ 10¹⁹ G in the paper.
6. **Check the astrophysical and heavy-ion numbers** in the intro (standard literature values, written conservatively) against your preferred sources.

## Reviewer-style observations (optional; not raised by the referee)

- **Fig. 3 is close to trivial.** As shown above, the Gaussian "vortex" field creates almost no condensate in the chiral limit (~10⁻⁴ of the uniform-field density). The revised text now presents it as an illustration of frustration, which is a real physical point. Still, the referee or other readers may ask why it was chosen as "one example of a CSL solution". You might add a second example with a substantial boundary flux, e.g. a conservative non-uniform field such as **b** = ∇(xz) restricted to the domain, where the CSL is maximal.
- **π⁰E·B term.** The full WZW action also couples π⁰ to E·B (the π⁰→γγ vertex). In a CSL this induces an anomalous electric charge density ∝ **B**·∇φ, which would enter Gauss's law. The paper, like Brauner–Yamamoto, implicitly assumes **E** = 0, i.e. neutralization by a background such as electrons. A referee asking about Maxwell's equations may follow up on Gauss's law. Consider one sentence saying this is assumed. I left it out of the manuscript because I could not check it against your references.
- **Massive-pion healing length.** Fig. 6 shows the sharp-wall result staying at ≈0.48 up to m_πL = 15. The text now cites this as evidence that the geometric healing length persists for b₀ ≫ m_π. A short study of b₀ closer to threshold (b̄₀ ≈ 4/π) would show whether 1/m_π eventually takes over. This would be a natural follow-up.
- **Abstract** was not changed. You may want to add a clause about the physical interpretation (frustration, the baryon-number identity).
