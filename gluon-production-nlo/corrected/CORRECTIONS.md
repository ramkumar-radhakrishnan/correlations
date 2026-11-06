# Corrected files — what was fixed (all changes marked in red in the tex/PDFs)

Every correction below is marked `{\color{red}...}` in the tex files and was
re-verified with the operator-level engine (`../verification-scripts/`).

## The decisive new check (this pass)

The final JIMWLK expression (`Evolution_JIMWLK__terms`, N1+N2) was compared with
the JIMWLK Hamiltonian acting on the LO cross-section,
H = λ ∫_z Σ_{p,q} K(p,q,z) [J_L J_L + J_R J_R − 2 J_L U(z) J_R]:
**exact agreement**, one fitted constant λ identical for SU(2) and SU(3),
relative residual ~7·10⁻¹⁶ (the swapped J_L↔J_R convention does *not* fit).
This certifies the collected two-ρ log(k⁺/Λ) terms as the JIMWLK evolution of
the LO cross-section, and (with the exact three-ρ cancellation) fixes which
version of each disputed block is correct. Script: `hjimwlk.py`.

## Group_I_corrected.tex

1. **Row II, N2a** — the logarithmic kernel applies only to the
   f^{abc}U^{bd}(x)U^{ce}(y) channel; for the f^{hbc}U^{ah}(w) channel the
   correct kernel is the ∫_z form (the structures 𝒦_i(x−w)𝒦(x−z)·𝒦(y−z) and
   𝒦_i(y−w)𝒦(y−z)·𝒦(x−z) must be absent). A corrected red equation was added.
   *Proof:* corrected N2a (+h.c.) ≡ collected Row-II blocks, 0 residual structures.
2. **Row III, N2** — the two user-red terms ∓(log[∨/Λ]+log[1−k⁺/∨])(…) are equal
   and opposite; removed (red note).
3. **Row III, N3a, 2nd term** — the log[(k⁺−Λ)/Λ] term with kernel
   (x′−w′)·(y−z)/((x′−w′)²(y−z)²(z−w)²) replaced (red) by the term with kernel
   𝒦(x′−w′)·𝒦(y−w) 𝒦(x−z)·𝒦(w−z) and first Wilson line U^{bc}(w).
   The 1/(z−w)² structure cannot be produced/cancelled by any other 3ρ source.
   *Proof:* K2-term + replacement (+h.c.) ≡ collected Row-III block, 0 residual.
   Red note: the five-piece p⁺ master (the B3 transcription from Ω) must be
   amended upstream to match.
4. Missing commas `x',y'w'` → `x',y',w'` (9 places).

## Evolution_three_rho_corrected.tex

1. **Row III five-term list** — signs of the first two
   [(x−w)²−(x−z)²]⁻² structures corrected to −,+ (they come from N2, which
   enters with an overall minus); the four such structures now cancel pairwise
   against terms four and five.
2. **Row III rewritten N2** — U^{bc}(z) → U^{bc}(w) (only this version matches
   the term it descends from and the collected block; engine-checked).
3. Red verification note added to the cancellation section: the complete
   collected 3ρ set cancels identically for arbitrary N_c (the incomplete β,β′
   sub-displays are covered by the global check).
4. Comma typo as above (1 place).
*Re-check after edits:* collection still cancels, 0 surviving structures.

## 4rho_corrected.tex

- **Master equation, second group**: 8 Wilson lines transposed (red) so the
  group follows the third group's orientation convention (U(z): cross-contracted
  index second; U(w′): external index first). *Proof:* a scan over all 2¹²
  orientation assignments shows the fully symmetrised 4ρ sector cancels for this
  assignment and (up to transposing both members of a contracted pair) no other.
- Red caveat: the 3ρ/2ρ reordering remainders quoted below the master equation
  are certified independently but do not yet follow mechanically from the master
  as printed; its transcription from Ω should be revisited. They, not the master,
  feed the evolution files.

## Group_II_corrected.tex
- Duplicated sentence "We write it in two parts:" removed (Row I).

## Group_III_corrected.tex
- Row II garbled factor Ā^{†(1)a}_i(k⁺,−U(w′)^{c′b}(z) → Ā^{†(1)a}_i(k⁺,−w′)U^{c′b}(z) (red).

## Evolution_two_rho_corrected.tex / BFKL_corrected.tex / Evolution_JIMWLK__terms_corrected.tex
- **No corrections needed**; red verification notes added (two-ρ file ≡ L_J parts
  of the Group results; BFKL and JIMWLK adjoint finals ≡ their colour-singlet
  *adjoint-trace* forms exactly; JIMWLK final ≡ −H·N_LO exactly).

## Other_terms_corrected.tex
- Red simplification added: closed form 𝒜₀ = θ/|(x−w)∧(x−z)| with
  θ = angle at w between (x−w) and (z−w) (exact).
- N2f note strengthened: the 1/(z−w)² coefficient cancellation is algebraically
  exact (identically, not just numerically).

## Calculation_main.tex
- Unchanged (no errors identified). The uploaded source referenced TikZiT figure
  files that were not provided; a minimal `tikzit.sty` is included so the PDF
  compiles (missing figures render as red placeholder boxes).

## Consistency statement
With these corrections the chain is consistent end-to-end:
Group files → three-ρ collection (cancels exactly, any N_c) and two-ρ L_J set
→ JIMWLK final = −log(k⁺/Λ)·H_JIMWLK·(LO), and BFKL/JIMWLK singlet trace forms
are exact rewritings. The only item left open upstream is the re-transcription
of the 4ρ master equation and of the Row-III p⁺ master formula from Ω, flagged
in red in the respective files.
