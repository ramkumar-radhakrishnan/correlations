# Single inclusive gluon production at NLO at mid rapidity

Draft paper: [`NLO_gluon_midrapidity.pdf`](NLO_gluon_midrapidity.pdf) (52 pages).

It covers the scattering of a dilute projectile on a dense target, at order $g^4$ relative to the
order-$g^2$ leading result, in the light-cone wave function approach. The soft window is
$\Lambda<k^+<\vee$ and the rapidity interval is $\delta Y=\log(\vee/\Lambda)$.

| Section | Content |
|---|---|
| II | Light-cone Hamiltonian; evolution operator $\Omega$ and its unitarity relations; cross section $\langle 0\vert\Omega^\dagger S^\dagger\Omega\, a^\dagger a\,\Omega^\dagger S\Omega\vert 0\rangle$; LO |
| III | NLO light-cone wave functions of the valence vacuum and of one soft gluon |
| IV | Coefficients $\mathbb N,\ \mathbb A^{(1)},\ \mathbb A^{(3)},\ \mathbb B_{2},\ \mathbb B_3,\ \mathbb C_1$ in coordinate space |
| V | Cross section at $g^4$: Group I (virtual), Group II (real), Group III (interference) |
| VI | Leading logs: $\log(\vee/k^+)\,[\mathrm{BFKL}\otimes\mathrm{LO}]+\log(k^+/\Lambda)\,[\mathrm{JIMWLK}\otimes\mathrm{LO}]$; three- and four-charge terms cancel |
| VII | Finite part, row by row, after exact subtraction of $\delta Y_T=\log(k^+/\Lambda)$ and $\delta Y_P=\log(\vee/k^+)$ |
| VIII | Running coupling ($\beta_0$), CMW constant, DGLAP and fragmentation |
| IX | Summary and list of open items |
| A–E | Conventions and Fourier transforms; longitudinal functions; explicit BFKL/JIMWLK action; leading-log content of every term; numerical checks |

## Build

```
pdflatex NLO_gluon_midrapidity.tex   # run three times (table of contents, references)
```

The document uses REVTeX 4.2 (`aps,prd`) when `revtex4-2.cls` is installed and falls back to
the standard `article` class otherwise. The PDF here was built with the fallback.

Open items are printed as red `[Note: ...]` remarks and marked `%TODO` in the sources. To hide
them, set `\shownotesfalse` in `NLO_gluon_midrapidity.tex`. Author names and affiliations are
placeholders.

`python3 flatten.py` writes `NLO_gluon_midrapidity_flat.tex`, a single-file version (e.g. for Overleaf or
arXiv).

## Corrections relative to the draft notes

- **Sign of the $p^+>k^+$ non-instantaneous one-gluon amplitude** (Eq. III.12, and $\mathbb B_3$ in Eq. IV.7):
  - The sign is fixed by requiring that the double pole $1/(p^+-k^+)^2$ cancels on both sides of $p^+=k^+$.
  - With this sign, every term $\propto 1/\Lambda$ in G1-III cancels.
- **Unitarity relation for $\mathbb A^{(3)\dagger}_1$:**
  - It contains the commutator $-[\mathbb N_2,\mathbb A^{(1)}_2]$, which is not zero.
  - The commutator contributes only at leading log.
- **Explicit JIMWLK $\otimes$ LO (App. C):**
  - The matrix $M$ is evaluated at the points where the derivatives act.
  - The term $dN^{(5)}$ enters with a positive sign.

## Open items

Section IX lists six open items:

1. The finite part of the diagonal term G1-I′.
2. G1-V part 2b in dimensional regularization, including its weight in the $\beta_0$ balance.
3. The adjoint index placement in G2-I(2), G2-IV and G3-III(2).
4. The bookkeeping of the commutator $[\mathbb N_2,\bar{\mathbb A}^{(1)}]$ in Group I.
5. The three-charge LL terms of G3-III part 1.
6. A direct re-derivation of the sign in Eq. (III.12).

## Follow-up notes

- [`three_rho/`](three_rho/README.md): the three-ρ terms ∝ log(∨/Λ) reordered into symmetric parts plus two-ρ terms, with every step written out and checked numerically.
- [`two_rho/`](two_rho/README.md): all two-ρ terms ∝ log(∨/Λ) compared with JIMWLK⊗LO and BFKL⊗LO. They give JIMWLK⊗LO exactly after six fixes; BFKL⊗LO multiplies log(∨/k⁺), not log(∨/Λ).
- [`bfkl/`](bfkl/README.md): all terms ∝ log(∨/k⁺) (four-, three- and two-ρ), reordered and compared with BFKL. Their coefficient is BFKL⊗LO − JIMWLK⊗LO exactly after one new fix (Group I Row III) and three carried-over fixes. The BFKL terms are written out, and Appendix E is confirmed.
- [`combined/`](combined/README.md): one note with everything above. It adds the Fourier transform of Group I Row II, step by step, and shows that the α part of that row is the commutator [ℕ⁽²⁾, 𝔸⁽¹⁾]: fix F4 is the p⁺ ordering of the virtual factor. It compares the coefficient of log∨ directly with BFKL⊗LO. Result: d³N_LL = log(k⁺/Λ)·JIMWLK⊗LO + log(∨/k⁺)·BFKL⊗LO.
- [`two_rho_all/`](two_rho_all/README.md): every two-ρ term ∝ log(∨/Λ) in one file: the genuine terms, those from reordering three ρ's and those from reordering four ρ's. There are 440 entries from 18 sources.
  - Their sum is JIMWLK⊗LO exactly.
  - Adding the log(∨/k⁺) coefficients, row by row, gives the coefficient of log∨, which is BFKL⊗LO exactly.
  - The note shows, with no limits on ∨, k⁺, Λ, why log(∨/Λ) carries JIMWLK: Λ enters only where a gluon is softer than the measured one. Group II Row IV is worked out as an example.
- [`three_rho_v2/`](three_rho_v2/README.md): check of the full list of three-ρ terms ∝ log(∨/Λ) (with c.c.).
  - As sent it does not cancel. Three typos (T1–T3) and the earlier fixes F3 and F4 make it cancel exactly.
  - The two-ρ terms need the operator order: c.c. terms reversed, and anticommutators as ½·A{B,C}.
  - The list is given in that order; its two-ρ terms, together with the genuine and `4rho.tex` terms, give JIMWLK⊗LO exactly.
- [`rowV/`](rowV/README.md): Group I, Row V. The derivation of `Group_V.tex` is checked. Part 2a is UV finite. Part 2b is UV divergent with the splitting function as coefficient; in MS-bar its β₀ part is −(α_sβ₀/4π)[1/ε + log(μ̄²/k²)] − (α_sN_c/4π)(67/9), times LO, per ordering, and the remainder takes the plus prescription. Combined with Group I Row II, the poles cancel, leaving the running-coupling ratio α_s(c₀/(x−w)²)/α_s(k²). [`PlusPrescription_2b.pdf`](rowV/PlusPrescription_2b.pdf) does the plus prescription of the UV-finite rest exactly on Λ < p⁺ < k⁺−Λ. The endpoint terms multiply ln((k⁺−Λ)/Λ), and every piece is UV finite.
- [`groups23/`](groups23/README.md): every equation of `Group_II.tex` and `Group_III.tex` checked against `Calculation_main.tex`, pointwise and after the integrations, with no limits in ∨, k⁺ or Λ.
- [`group1/`](group1/README.md): every equation of `Group_I.tex` (one gluon in the final state) checked against `Calculation_main.tex`, as for Groups II and III. Rows IV and V (up to the plus-prescription limits) are correct. In Row III the p⁺ > k⁺ non-instantaneous line has the wrong sign (𝔹₃, E3); with the fix every 1/Λ term cancels. Row II part 2 lacks −kⁱ log(∨/k⁺) and needs the commutator [ℕ₂, 𝔸⁽¹⁾] (F4). The logs of Row III are principal values (log|·|).
- [`groups23/main/`](groups23/main/README.md): all corrections to `Calculation_main.tex` (sign errors in onegrho12, onegrhoX1 and 𝔹₃, the double counting in 𝔸⁽³⁾₂, the unitarity commutator, ranges, Appendix F, typos), exact forms for every Λ → 0 step, and a marked-up corrected copy of the file.
- [`three_rho_v3/`](three_rho_v3/README.md): the full three-ρ list with c.c. (152 terms) cancels exactly: 40 colour networks, each summing to zero, organised by kernel class (α, α′, β, β′) and number of Wilson lines (2, 3, 4).
  - Most lines are correct.
  - Group III Rows III, V and VI have the wrong overall sign (same error as F5).
  - The plus prescription must keep the lower limit Λ.
  - Group II Row IV from the z integration on is corrected (constant in the log, exact ξ limits and endpoint log, dξ/ξ³).
  - Index placements are given for a general adjoint U.
  - Corrected copies of both notes are included.
- [`final_check/`](final_check/README.md): final check of the evolution terms of the last upload (`Evolution_terms.tex`, `Evolution_three_rho.tex`, Groups I–III).
  - Three and four ρ's cancel.
  - The two-ρ log(∨/Λ) terms give JIMWLK⊗LO exactly once the two-ρ lists of `4rho.tex` are added (they are missing from `Evolution_terms.tex`).
  - The JIMWLK formula as sent needs two corrections.
  - BFKL⊗LO is the coefficient of log ∨.
  - The pole terms in the log(∨/Λ) sector cancel; DGLAP comes from Group II Row IV.
  - Every two-ρ term is coloured by kernel class (JIMWLK only, charge–charge, BFKL only).
- [`bfkl_final/`](bfkl_final/README.md): your `bfkl.tex` (the log(∨/k⁺) list), checked to the end.
  - Three-ρ cancellation, written out with Φ[p,q] structures and numerically.
  - The two-ρ sum, structure by structure, against ℬ − 𝒥.
  - The BFKL term and its colour-singlet form.
  - Two contributions were missing: (A) 𝒟 in place of the edge kernel in Group I Row III; (B) the y′ → z logarithms of Group III Rows III, V, VI.
  - With them, the three-ρ terms cancel to 10⁻¹⁷ and all 67 two-ρ structures equal log(∨/k⁺)(ℬ − 𝒥). Total: log(k⁺/Λ)𝒥 + log(∨/k⁺)ℬ.
- [`finite_final/`](finite_final/README.md): your `Finite_terms.tex` (the terms apart from BFKL and JIMWLK), checked.
  - 20 of 27 equations are correct as written. Changes C1–C5: the −log(∨/k⁺)H₄ line, [1/p⁺]₊^{(k⁺)}, the y′ → z subtraction, the edge-kernel terms of Group I Row III, one measure.
  - Running coupling and DGLAP re-checked; the final finite terms are listed.
  - In ξ = p⁺/(p⁺+k⁺) the rows take the form of AABKL Eq. (4.22): real ∫₀^{1−x_p}dξ with [ξ]₊, virtual ∫₀¹dξ, DGLAP, soft term = JIMWLK.
- [`llsplit/`](llsplit/README.md): from your new `Group_I/II/III.tex` and `4rho.tex`, every term ∝ log(k⁺/Λ) (JIMWLK) and ∝ log(∨/k⁺) (BFKL), row by row (`LL_terms.pdf`), and all other terms with Λ → 0 (`Finite_terms_rows.pdf`).
  - The three-ρ terms of both lists cancel. The two-ρ terms are JIMWLK⊗LO and BFKL⊗LO.
  - Changes needed: E1 (2 log(k⁺/Λ) in Group I Row V 2b UV), 𝒟 in Group I Row III, and the hidden logarithms of Group III Rows III, V, VI and Group II Row IV.
- [`bfkl_rows/`](bfkl_rows/README.md): your `BFKL.tex` and `Other_terms.tex`. The BFKL terms are row by row (`BFKL_rows.pdf`) and the remaining terms are row by row with Λ → 0 (`Remaining_rows.pdf`).
  - Two lines of `BFKL.tex` must go: Group I Row V 2b UV and Group II Row IV line 1.
  - The hidden log(∨/k⁺) parts of `Other_terms.tex` and (B) of Group III Rows III, V, VI join the BFKL list, and `Other_terms.tex` subtracts H₄ (C1).
  - The three-ρ terms then cancel (3.5e-16), and the two-ρ terms are your "BFKL contribution" = BFKL⊗LO.
  - `Rest_terms.pdf`: everything apart from JIMWLK and BFKL, all black. `Why_log_vee.pdf`: step by step (your Group III Row VI) where the log(∨/k⁺) of Group III comes from.
