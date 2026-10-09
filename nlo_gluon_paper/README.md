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
