# NLO single-inclusive gluon production — verification report (pass 1)

**Files checked:** `Calculation_main.tex`, `Group_I_corrected.tex`, `Group_II_corrected.tex`,
`Group_III_corrected.tex`, `4rho_corrected.tex` (uploads of 2026-11-06).

**Conventions used throughout** (as stated by the author): adjoint Wilson lines are real
orthogonal, $U^{\dagger\,bc}=U^{cb}$; $d_\perp = 2-\epsilon$; no limits taken
(in particular $\log(\vee/k^+-1)$ is kept distinct from $\log(\vee/k^+)$);
$\rho$'s are non-commuting and symmetrized via
$\rho^a(x)\rho^b(y)=\tfrac12\{\rho,\rho\}+\tfrac12[\rho,\rho]$.

Shorthand below: $L_J \equiv \log(k^+/\Lambda)$ (JIMWLK-class log),
$L_B \equiv \log(\vee/k^+)$ (BFKL-class log),
$A=(x-z)^2$, $B=(x-w)^2$, $C=(z-w)^2$ and primed analogues.

---

## 1. What was verified and how

**Scope of this pass.** I verified the *internal derivation steps* inside the five files:
every "integrating w.r.t. $p^+$" step, the $\delta^{(2)}$ (i.e. $w$, $w'$) integrations and
their Jacobians/phases, the plus-prescription extractions, the UV ($1/\epsilon$)
extractions, the Wilson-line reductions, and the $4\rho$ cancellation claims.
I did **not** yet re-derive each row's *starting* expression from the evolution operator
$\Omega$ (amplitude assembly), nor the Fourier-transform appendix identities (one
exception noted below), nor the final $2\rho$ assembly into BFKL/JIMWLK — see §5.

**Method.** The nontrivial $p^+$ integrations were verified *numerically*: random 2D
transverse configurations $(x,y,z,w,w',x',y')$ are drawn, the pre-integration integrand
(exactly as written in the file, including the full tensor contractions) is integrated
over $p^+$ with multi-precision quadrature (principal value where the denominators
vanish inside the range), and compared to the claimed post-integration expression at
the same configuration. The common colour/$\rho$ brackets factor out and drop. Scripts
are in `verification-scripts/` (`g2_row1_part2.py`, `g1_row3.py`, `misc_checks.py`).

---

## 2. Verdicts, file by file

### Group III (two gluons in final state, with ℂ₁ pieces)

| Step | Verdict |
|---|---|
| Row I, eq. (1.2) → (1.3): $\int dp^+$ of $\frac{1}{p^+A+k^+B}\big[\frac{\delta^{ij}(A-B)}{2C(p^++k^+)}+\frac{1}{p^+}(\cdots)+\frac{1}{k^+}(\cdots)\big]$ | **Correct** (verified analytically term by term; all three output structures incl. the $\log[\vee/(k^++\Lambda)]$ $\delta^{ij}$ piece and the blue $[L_J+\log(\vee/k^+-1)]$ piece have the right coefficients and signs) |
| Row II: same longitudinal bracket, same result | **Correct** (identical integral) |
| Row III part 1: $w'$ integration (Jacobian, $(w'-z)\to\frac{p^++k^+}{k^+}(y'-z)$, phase $e^{-ik(y'-w)}e^{-ik\frac{p^+}{k^+}(y'-z)}$) and the $[1/p^+]_+$ extraction with coefficient $-\frac{ig^4f}{16\pi^5}\frac{1}{k^+}[L_J+\log(\vee/k^+-1)]$ and the correct $\delta_{jm}\delta_{ik'}$ contraction $(y'-z)\!\cdot\!(x-z)\,(x'-y')\!\cdot\!(y-w)$ | **Correct** |
| Row III part 2: plus-prescription extraction (the $-\frac{g^4}{8\pi^5}\frac1{k^+}[{\rm blue}]$ term, $F(0)=1/(k^+B)$, both tensor structures) | **Correct** |
| Row V: $w'$ integration; $[1/p^+]_+$ extraction $\log[(\vee-k^+)/\Lambda]$ (later rewritten as the blue combination — same thing); Wilson contraction $U^{d'b}(z)U^{cb}(z)=\delta^{d'c}$, $f^{c'd'a}\to f^{c'ca}$ | **Correct** |
| Row VI: mirrors Row V | **Correct** |

### Group II

| Step | Verdict |
|---|---|
| Row I part 1: $\int dp^+$ (both the unprimed and the $w\leftrightarrow w'$-swapped halves, mirror signs) | **Correct** (same integral as Group III Row I; prefactors $\mp ig^4/16\pi^5$, $\pm ig^4/32\pi^5$ consistent) |
| Row I part 2 (the $\mathbb{B}_2^\dagger\mathbb{B}_2$ / $N_2$ piece): eq. (56) → results **2a–2f** | **Correct — verified numerically to $10^{-22}$ relative**, *exactly*, at finite $\Lambda$ (no dropped terms). This includes all bilinear tensor brackets, all five log structures and the power term $[1/\vee - 1/(k^++\Lambda)]$ in 2f. |
| Row IV: $w,w'$ integrations of the double $\delta^{(2)}$ (Jacobians $(p^++k^+)^4/k^{+4}$, kernel factors $(k^+/(p^++k^+))^2$, phase $e^{-ik(y'-y)(p^++k^+)/k^+}$); the tensor-product bracket $[\delta\delta\,d_\perp p^+/(k^++p^+)^2 + \frac{2}{k^+}(\delta\delta-\delta\delta)+ (p^+/k^{+2}+1/p^+)\delta\delta]$ | **Correct** (expanded independently; incl. the subtle fact that the $\frac{2}{k^+}$ term drops in both the $\delta^{mm'}$ and $r^mr^{m'}$ projections — your eq. after "Simplifying" handles this correctly) |
| Row IV: Wilson reduction $f^{c'd'a}f^{cda}U^{d'b}(z)U^{db}(z)=N_c\delta^{c'c}$ | **Correct** |
| Row IV: the $z$-integral identity $\int_z \frac{(y-z)^m(y'-z)^{m'}}{(y-z)^2(y'-z)^2} = -\frac{\pi}{2}[\frac{2}{\epsilon}+2\gamma_E+\log\frac{(y'-y)^2\mu^2}{4}]\delta^{mm'}-\pi\frac{(y'-y)^m(y'-y)^{m'}}{(y'-y)^2}$ | **Coefficients confirmed numerically** (traceless part $=-\pi\,r^mr^{m'}/r^2$ to 0.3%; trace $\log$ coefficient $=-\pi$ per the identity). The exact scheme constants ($2\gamma_E$, $\log 4$) are not distinguishable numerically from the $1/\epsilon$ mapping; they are standard $\overline{\rm MS}$ bookkeeping. |
| Row IV: $\xi$ substitution ($p^+/(k^++p^+)^2\to\xi(1-\xi)/k^+$ etc., $dp^+=-k^+d\xi/\xi^2$, reversed limits), plus prescription in $\xi$, and the red endpoint log $\log\frac{(\vee-k^+)(k^++\Lambda)}{\vee\Lambda}$ from $\int_{k^+/\vee}^{k^+/(k^++\Lambda)}\frac{d\xi}{1-\xi}$ | **Correct** (the red log is exact) |
| Row IV: $\mathcal{P}_{gg}$ assembly and pole isolation, incl. the "$-1$" in $[-d_\perp/\epsilon-d_\perp\gamma_E-\frac{d_\perp}{2}\log(\cdot)-1]$ | **Correct** (coefficients $g^4N_c/8\pi^4$ consistent) |

### Group I (one gluon in final state)

| Step | Verdict |
|---|---|
| Row II $N_1$: the rewrite $\frac{11}{3}+\log\frac{\Lambda}{k^+}-\log\frac{\vee}{\Lambda}+\log(\frac{\vee}{k^+}-1) = \frac{11}{3}-2L_J+\log(1-\frac{k^+}{\vee})$ | **Correct** (exact identity; also consistent with $\Psi_{g\rho}$, eq. (grhoposition) of the main file) |
| Row II $N_2$: split into $2a$ ($\propto L_J$) and $2b$ ($R^i$, no cutoff logs): I verified **symbolically** that [eq. (58) bracket] $- R^i$ is *exactly* $\propto L_J$ with **zero** residual $\vee$-dependence, and that its coefficient equals $-\big[\frac{k^2p^i-(k\cdot p)k^i}{(k-p)^2}+2k^i+\frac{p^2k^i+k^2p^i}{(k-p)^2}\big]$ | **Correct** as a momentum-space split. The *position-space* form of 2a (the $\log\frac{(w-y)^2}{(x-y)^2}$ kernel) relies on your FT appendix and was **not independently re-derived** in this pass. |
| Row III: eq. (85) [all five $p^+$ pieces, incl. the vector-square denominator $[k^+(x-w)-p^+(x-z)]^2$] → $N_1$ + $N_2$ + $N_{3a}$ + $N_{3b}$ with the arctan $\mathcal{A}$ | **Correct — verified numerically to $10^{-13}$**, with principal-value integration and logs read as $\log\lvert\cdot\rvert$ (see §4). Piece-by-piece: pieces 1+2 ↔ $N_1$ items 1–6; piece 3 ↔ $N_1$ items 7–8; piece 4 ↔ $N_2$; piece 5 ↔ $N_{3a}+N_{3b}$. The red pair inside $N_2$ cancels identically (it is written as $-X+X$) and can simply be deleted. |
| Row IV: $\int dp^+$ | **Correct** (it is the exact sign-mirror of Group III Row I, which was verified; all claimed signs are flipped consistently) |
| Row V part 1: $w$ integration (phase $e^{-ik(w'-z)}e^{-ik(z-x)p^+/k^+}$), the $q^+=k^+-p^+$ change of variables for the $\frac{k^+}{k^+-p^+}\delta_{ij}\delta_{km}$ term (phase $\to e^{-ik(w'-x)}e^{-ik(x-z)q^+/k^+}$, contraction $(x'-w')\cdot(y'-x)(z-x)\cdot(y-z)$), both $\log[(k^+-\Lambda)/\Lambda]$ extractions with coefficient $-\frac{ig^4}{32\pi^5 k^+}$ and the correct phase on each | **Correct** (by inspection of all kernels/phases/signs) |
| Row V part 2: eq. (158)→(162) (every longitudinal factor rescales by $p^+(k^+-p^+)/k^{+2}$ correctly); the identity $f^{abc}f^{def}U^{be}U^{cf}=N_cU^{ad}$ usage; the $2a$/$2b$ split; the $q^+$ integrals with reversed limits $\int_{k^+-\Lambda}^{\Lambda}$ and the $[k^+/q^+]_+$ bookkeeping of eq. (→233) incl. which phase ($e^{-ik(w'-z)}$ vs $e^{-ik(w'-x)}$) accompanies each extraction; the counterterm added back with coefficient $[-\frac{11}{6}+2L_J]$ | **Correct** (by inspection) |
| Row V 2b-UV: the longitudinal factor $2\log\frac{\vee}{\Lambda}-2\log\frac{\vee}{k^+}+2\log(1-\frac{\Lambda}{k^+})-\frac{11}{6}+\frac{4\Lambda}{k^+}-\frac{\Lambda^2}{k^{+2}}+\frac{2\Lambda^3}{3k^{+3}}$ | **Correct and exact**: it equals $\int_\Lambda^{k^+-\Lambda}\frac{dp^+}{k^+}\big[\frac{p^+(k^+-p^+)}{k^{+2}}+\frac{p^+}{k^+-p^+}+\frac{k^+-p^+}{p^+}\big]$ identically (verified analytically, including the power-suppressed terms). The overall normalization of the accompanying $[\frac{2}{\epsilon}+2\gamma_E+\log\frac{(y-z)^2\mu^2}{4}]$ factor rests on the FT appendix constant (same caveat as above). |

### 4ρ file (`4rho_corrected.tex`)

| Claim | Verdict |
|---|---|
| $4\rho$, two-Wilson-line bracket (12 terms) $=0$ | **Correct** — the sum vanishes identically once one uses that the $U$ matrix *elements* commute (they are c-number functions of the target field); the terms pair up as $(-AB-BA+2AB)$ etc. |
| $4\rho$, four-Wilson-line bracket (6 surviving terms) $=0$ | **Correct** — verified symbolically: the kernel is symmetric under $x\leftrightarrow y$, and under $x\leftrightarrow y$ together with the dummy relabel $d\leftrightarrow e$ the six terms map pairwise into minus each other. (Caveat: this uses the $x\!\leftrightarrow\!y$ relabel *through* the $\rho^d(x)\rho^c(u)\rho^e(y)$ ordering; that is legitimate at the level you use it because these are the already-symmetrized/ordering-insensitive remainders.) |
| $3\rho$ remainders (two- and three-Wilson-line pieces) | These are **nonzero as written** — they are *supposed* to cancel against the genuine $3\rho$ log terms of Groups I–III. That cancellation is **not yet demonstrated** in the uploaded files: the genuine-$3\rho$ collections and the group-wise $3\rho$ ledgers exist only in commented-out blocks of `4rho_corrected.tex`. See §5, open item (ii). |

---

## 3. Log ledger — where the BFKL and JIMWLK logs live

Classifying every cutoff log in the files ($L_J=\log\frac{k^+}{\Lambda}$ JIMWLK class,
$L_B=\log\frac{\vee}{k^+}$ BFKL class):

* **Full-interval logs.** The blue combination $[L_J+\log(\frac{\vee}{k^+}-1)] = \log\frac{\vee-k^+}{\Lambda}$
  (Groups II and III real-emission terms; Group III Rows III/V/VI extractions; Group I Row III)
  contains *both* classes: it splits as $L_J + \log\frac{\vee-k^+}{k^+}$. So for these terms
  the JIMWLK-log coefficient and the BFKL-log coefficient are **the same kernel** — your
  framework prediction (that the same structures must reassemble into JIMWLK at the
  $L_J$ end and BFKL at the $L_B$ end) is consistent with this, but it means the two
  checks are not independent for these terms. Similarly $\log\frac{\vee}{\Lambda}$
  ($\mathbb{N}_2$, $\Psi_{g\rho\rho\rho}$, parts of $\mathbb{A}^{(3)}$) $=L_J+L_B$.
* **Pure BFKL-class logs:** $\log\frac{\vee}{k^++\Lambda}$ (Group II 2d/2e/2f, Group III Rows I/II
  second structures, Group I Row IV) — these multiply $1/[(z-w)^2(z-w')^2]$-type kernels,
  i.e. exactly real-BFKL dipole-kernel structures with $z$ = unobserved-gluon position. Also
  $\log(1-\frac{k^+}{\vee})$, $\log^2(1-\frac{k^+}{\vee})$ (Group I Row II constants).
* **Pure JIMWLK-class logs:** $\log\frac{k^+-\Lambda}{\Lambda}$ (Group I Rows III/V, from the
  $p^+<k^+$ region), $\log\frac{\Lambda}{k^+-\Lambda}$, $-2L_J$ (Group I Row II $N_1$),
  $\log\frac{\Lambda}{k^+}$ ($\mathbb{A}^{(3)}$, $\mathbb{B}$-coefficients).
* **$\rho$-counting of log terms** (what must happen, per your plan):
  * $4\rho$: cancels — **verified** at the level collected in `4rho_corrected.tex`.
  * $3\rho$: the symmetrization remainders (4ρ file, "Three $\rho$" section, final eqs) must cancel
    against the genuine $3\rho$ log terms of Group I Row II($N_2$ 2a)/Row III/Row IV/Row V(part 1),
    Group II Row I (first part), Group III Rows I/II/III(part 1)/V/VI. **Open** — this is a finite,
    well-defined colour-algebra computation; see next steps.
  * $2\rho$: survivors must build $\frac{dN^{(i)}}{d^2wd^2w'dy}$ (JIMWLK, $\propto L_J$) and the
    dipole-BFKL expression ($\propto L_B$). The surviving $2\rho$ structures in the files are
    indeed all of the required forms: $[U(x')-U(w')]\otimes[U(y)-U(\cdot)]\rho\rho$ with
    $\mathcal{K}\cdot\mathcal{K}$ kernels (Group I Row V 2b, Group II Row IV), the
    $ff\,UU/N_cUU$ four-term brackets (Group II Rows I/IV), and the $z$-contracted
    Group III Rows V/VI structures. The full coefficient assembly is the step after the
    $3\rho$ cancellation. **Open.**

---

## 4. The negative-log / PV question (answer only, not implemented)

The issue is real and quantified: in Group I Row III the denominators
$p^+(x-z)^2-k^+(x-w)^2$ vanish *inside* the integration range for generic coordinates —
for $(x-w)^2>(x-z)^2$ the pole $p^+_0=k^+B/A$ sits in the $p^+>k^+$ piece, otherwise in the
$p^+<k^+$ piece, and piece 3 ($\int_\Lambda^\vee$) essentially always contains it. My numerical
check confirms your results hold **exactly in the PV sense with every such log read as
$\log\lvert\cdot\rvert$**.

What I propose to do about it (for discussion before implementing):

1. **Fix the $i\epsilon$ from LCPT once, per energy denominator.** These poles descend from
   old-fashioned energy denominators, which carry an unambiguous prescription
   ($E_{\rm interm}-E_{\rm in}-i\epsilon$ from the adiabatic switching). Propagating that
   $i\epsilon$ through the transverse Fourier transforms replaces the ill-defined
   $\log({\rm negative})$ by $\log\lvert\cdot\rvert \mp i\pi\,\theta(\cdot)$ with a *definite* sign,
   instead of the PV convention.
2. **Then classify each affected term by the Hermiticity of its $U/\rho$ structure** under the
   "$+$ c.c." operation that accompanies every Group I row. For structures that are Hermitian,
   the $\pm i\pi$ pieces cancel between the term and its conjugate and the PV answer
   ($\log\lvert\cdot\rvert$) is *exact* — nothing to fix. For anti-Hermitian structures the
   opposite happens: the $\log\lvert\cdot\rvert$ part cancels and a **real $\pi$-term (no log)
   survives**. These would be genuinely missed by the naive PV reading, so they must be
   collected term by term. (They cannot contribute to the $L_J/L_B$ coefficients — they carry
   no cutoff log — so they do not endanger the BFKL/JIMWLK identification, but they do
   affect the finite NLO remainder.)
3. Practically: this only touches Group I Row III (and the analogous $\mathbb{B}_3$ interference
   structures); everywhere else the $p^+$ denominators $p^+A+k^+B$ are positive definite.

So: **for the BFKL/JIMWLK extraction you can safely proceed with $\log\to\log\lvert\cdot\rvert$ (PV)**;
the $i\pi$ bookkeeping is a finite-term correction to be done once, at the end, with the
LCPT $i\epsilon$.

---

## 5. Open items / next steps (in order)

1. **(blocking question)** The genuine $3\rho$ collections ("Group I/II/III" $3\rho$ log terms)
   and the two-Wilson-line symmetrization ledgers exist only as *commented-out* blocks in
   `4rho_corrected.tex`. Are those blocks current (post-correction), or superseded? Once
   confirmed (or replaced), I will run the $3\rho$ cancellation check with the same
   canonicalization machinery (dummy-index relabelling + integration-variable permutation
   + $U$ orthogonality + $f$-identities).
2. **$2\rho$ assembly:** sum all surviving $2\rho$ log terms (including $\mathbb{N}_2$ and the
   virtual Group I Row II pieces) separately at order $L_J$ and $L_B$, and match against your
   JIMWLK $dN^{(1..5)}$ and dipole-BFKL targets. Prediction from the ledger in §3: the
   $1/\epsilon$ poles must cancel among Group I Row II $N_1$ ($\propto \frac{11}{3}-2L_J$),
   Group I Row V 2b-UV ($\propto 2L_J-\frac{11}{6}$), and Group II Row IV
   ($\mathcal{P}_{gg}$ pole, $\propto$ the $\xi$-integral over $[k^+/\vee, k^+/(k^++\Lambda)]$)
   — this is a strong self-consistency check to do first within the $2\rho$ assembly.
3. **Second-pass provenance check:** re-derive each row's starting expression from the
   $\Omega$-coefficient catalogue (§"Parametrization of the evolution operator" of the main
   file) — this pass verified the steps *within* each row, not the row's seed expression.
4. **FT appendix:** independently verify the two position-space Fourier identities used in
   Group I Row II 2a and Row V 2b-UV.

## 6. Minor typos found (no physics impact)

* `Group_III_corrected.tex`, Row II first display: garbled argument
  "$\bar{\mathbb{A}}^{\dagger(1)a}_i(k^+,-U(w')^{c'b}(z)$" — should read
  "$\bar{\mathbb{A}}^{\dagger(1)a}_i(k^+,-w')\,U^{c'b}(z)$".
* `Group_I_corrected.tex`, Row III $N_2$: the two red terms are identical with opposite
  signs; they cancel and can be deleted.
* `Group_I_corrected.tex`, Row V: several integration subscripts "$x',y'w'$" missing a comma.
* `Group_II_corrected.tex`, Row I: "We write it in two parts:" duplicated.

---

### Bottom line

Every derivation step I could test in the five files is **correct** — including the two
heaviest ones (Group II Row I part 2, verified to $10^{-22}$; Group I Row III with its
PV poles and arctan structures, verified to $10^{-13}$) — and the claimed $4\rho$
cancellations hold. The $3\rho$ cancellation and the final $2\rho$ → BFKL/JIMWLK matching
are not yet demonstrated in the files and are the natural next computations; nothing found
in this pass obstructs them.
