# NLO gluon production — verification of the 3ρ cancellation, JIMWLK, BFKL and finite terms (pass 2)

**Files checked:** `Evolution_three_rho_corrected.tex`, `Evolution_two_rho_corrected.tex`,
`Evolution_JIMWLK__terms_corrected.tex`, `BFKL_corrected.tex`, `Other_terms_corrected.tex`
(uploads of 2026-11-06), against the Group I/II/III files verified in pass 1.

**Method (independent, from scratch).** I built an operator-level verification engine:
- every term is represented exactly (kernel ⊗ colour ⊗ ordered ρ-product, phase
  $e^{-ik(w'-w)}$), with $\{\rho,\rho\}=AB+BA$ and
  $[\rho^a(x),\rho^b(y)]=if^{abc}\delta^{(2)}(x-y)\rho^c(x)$ implemented symbolically
  (contact terms generated automatically);
- a *symbolic canonicalizer* reduces each term to a canonical form under: relabelling of
  all contracted colour indices, relabelling of all internal position variables
  ($x,y,z,u,x'$ — only $w,w'$ are pinned by the phase), $f$-antisymmetry, kernel
  orientation signs, δ-absorption, same-point orthogonality $U^{ab}(s)U^{cb}(s)=\delta^{ac}$
  and same-point $f$-invariance. Sums of terms then cancel **iff** they cancel exactly as
  integrands — no numerics, valid for every $N_c$;
- a numeric cross-engine (random 2D points; adjoint Wilson lines as true SO(3)/8×8-adjoint
  orthogonal matrices from random SU(2)/SU(3) elements) confirms all symbolic verdicts.
- The sources I compare against are **my own transcriptions of the Group-file results**,
  which were verified at the $p^+$-integration level in pass 1 — not your collected files.

Scripts: `verification-scripts/` (`engine.py`, `symcolor.py`, `texparse.py`, `sources_LJ.py`,
`fourrho.py`, `g3cancel.py`, `perrow.py`, `singlet_bfkl.py`, `singlet_jimwlk.py`).

---

## 1. What is CONFIRMED

1. **The collected 3ρ set cancels exactly.** Parsing all blocks of
   "All three ρ term with c.c." + "Genuine terms" (three-ρ file, lines 182–255):
   40 distinct kernel⊗colour structures, every coefficient sums to zero, for arbitrary
   $N_c$. Your α/α′/β/β′ bookkeeping is internally exact.
2. **Eight of the eleven per-row collected 3ρ blocks are exactly equal to the verified
   Group-file sources** (my independent transcriptions): G-I Row IV, G-I Row V (part 1),
   G-II Row I, G-III Rows I, II, III (part 1), V, VI. The c.c. conventions came out
   consistently: the collection includes the h.c. for Groups I and III; Group II blocks are
   self-complete. (These matches also fix the Wilson-line orientation questions for those
   rows: after same-point orthogonality the "red orientation" versions are equivalent to
   the Group-file ones.)
3. **The `Evolution_two_rho` file is exactly the $\log(k^+/\Lambda)$ part of the verified
   Group results** (checked line by line against pass-1 output).
4. **BFKL final: the adjoint form (lines 570–605) equals the colour-singlet trace form
   (lines 635–644) EXACTLY** (relative difference $10^{-15}$, $N_c=2$ and $3$) — provided
   the traces are read as **adjoint** traces ($U$ adjoint, $T^e$ adjoint). The same holds
   for **JIMWLK: $N_1+N_2$ (lines 657–707) equals the trace form (lines 738–762) exactly.**
   So the singlet conversions are error-free.
5. **Finite terms (`Other_terms`) spot checks pass exactly** where I checked them against
   the verified Group results:
   - G-III Row I: the finite logs
     $\log\frac{(1-k^+/\vee)(x-z)^2+(k^+/\vee)(x-w)^2}{(x-w)^2}$ and $\log(1-k^+/\vee)$
     with their kernels are precisely the $\Lambda\to0$, $\log$-split remainder of the
     verified $\log\frac{(\vee-k^+)A+k^+B}{\Lambda A+k^+B}$ and blue-log terms. ✓
   - G-I Row III $N_1$: all three finite logs (incl. the $\lvert\cdot\rvert$ PV readings)
     match the exact $\Lambda\to0$ limits of the verified $n_2,n_3,n_5,n_6$ terms. ✓
   - $\mathcal{A}_0$ is the exact $\Lambda\to0$ limit of the verified arctan $\mathcal{A}$. ✓
   - The **red-added fourth term of $N_{3b}$ is correct**: it is exactly my verified $b_4$
     term, with $L_4=\log Z\to\log\frac{(z-w)^2}{(x-w)^2}$. ✓
   - BFKL file, G-III Row I: the $\log(\vee/k^+)$ row content is exactly the verified
     $L_B$-coefficient of the row. ✓
6. Normalisation bookkeeping: $(2\pi)^3\cdot16\pi^5=128\pi^8$ matches
   $\alpha_s^2/(8\pi^6)$ with $\alpha_s^2=g^4/16\pi^2$ — your prefactors are consistent
   with the target normalisation.

## 2. What DOES NOT check out (three localized faults)

**(A) Group I, Row III — the three-ρ file's rewrite is not equal to the source.**
The verified Group-I content (pass-1, numerically exact against the $p^+$ integral) gives
for the $\log(k^+/\Lambda)$, 3ρ part of Row III **exactly two structures**
(the $X$-type $(B-A)^{-2}$ kernels cancel pairwise between $N_2(m_5,m_6)$ and
$N_{3a}(a_3,a_4)$; $n_2$ cancels $n_5$):
$$
f^{abd}\Big[\tfrac{i}{16\pi^5}\,\mathcal K(x'-w')\!\cdot\!\mathcal K(x-w)\,
\mathcal K(x-z)\!\cdot\!\mathcal K(y-z)
\;-\;\tfrac{i}{8\pi^5}\,\frac{\mathcal K(x'-w')\!\cdot\!\mathcal K(y-z)}{(z-w)^2}\Big]
\big[U(x')-U(w')\big]\rho(x')\,U^{de}(x)\big[U^{bc}(z)-U^{bc}(y)\big]\rho^e(x)\rho^c(y).
$$
The first structure matches your collection; the second (with $1/(z-w)^2$, from the $a_2$
term) does **not** appear there — your intermediate step (three-ρ file, Row III,
"We write it as follows", eq. with kernel
$\mathcal K(x'-w')\cdot\mathcal K(y-w)\,\mathcal K(x-z)\cdot\mathcal K(w-z)$ and the
$U^{bc}(w)$ bracket) is not a legal relabelling of the source: no relabelling can put $x$
into that kernel or map $z\mapsto w$ (the phase pins $w$). Moreover the five-term list
above it carries $+2X_1-2X_2$ instead of the source's net zero.
**Important consequence:** the $1/(z-w)^2$ structure cannot be produced by any other 3ρ
source (the 4ρ reordering can never generate $(z-w)^{-2}$ since $z,w$ are not ρ-positions).
So *if* the Group-I file's Row III master formula (its eq. for the five $p^+$ pieces) is
correct, the true 3ρ sum does **not** cancel — the collected cancellation then relies on
this rewrite plus compensating changes elsewhere. The decision point is the provenance of
Row III's starting expression from $\Omega$ (the one thing neither pass has re-derived).

**(B) Group I, Row II — 4 structures differ.** After converting the $N_2$-2a logs with
$\log\frac{(a-b)^2}{(c-b)^2}=-\frac1\pi\int_z[\mathcal K(a-z)-\mathcal K(c-z)]\cdot\mathcal K(b-z)$
(exact), my source form differs from your collected GI-RowII blocks by four structures of
unit strength, all in the $f^{hbc}U^{ah}(w)\{\rho(x),\rho(y)\}$ colour channel
(your collection appears to drop/alter part of that channel's $z$-kernel pieces).

**(C) The 4ρ master equation (4rho file, the big three-group equation) is not
self-consistent as written.** Its 4ρ (fully symmetrised) sector does **not** cancel under
any c.c. convention (14–19 surviving structures), although your *collected* two- and
four-Wilson-line sums (which I verified in pass 1) do cancel. The reason is visible in the
z-Wilson-line orientations: group C of the equation has the corrected orientation
(cross-contracted index second, so $U^{d'b}(z)U^{cb}(z)=\delta^{d'c}$), while group B is
stale ($U^{bc}(z)$, giving $UU$ instead of $\delta$), and group A's $U^{bd}(z)$ is
ambiguous. Partially fixing group B reduces the survivors (14→8) but does not eliminate
them. Until this master equation is written with one consistent orientation convention,
the genuine from-4ρ 3ρ remainder (your "From four ρ's" blocks) cannot be verified from it —
and it is the third disputed input, besides (A) and (B), on which the collected 3ρ
cancellation depends.

**Summary of 2:** your collected 3ρ set is internally exact and 8/11 of its blocks are
provably source-faithful; but blocks GI-RowII, GI-RowIII and "From 4ρ" deviate from the
(pass-1-verified) sources in ways that must mutually compensate for the total to cancel.
Either those three were silently corrected relative to the Group files (then the Group
files and the 4ρ master equation need updating, and Row III's $p^+$ master formula must be
re-derived from $\Omega$), or they carry compensating errors. This is the single open knot.

## 3. Status of JIMWLK / BFKL identification

- Internal consistency of your JIMWLK and BFKL files downstream of the collection is
  **excellent**: both finals convert exactly to the colour-singlet (adjoint-trace) forms,
  and the per-row $L_B$ extractions match the verified sources where checked.
- Comparison against the target formulas you pasted in the message was **inconclusive for
  a trivial reason**: the pasted LaTeX was mangled in transit (the signs of several terms,
  especially in $dN^{(5)}$ and the $K_D$ definition, were eaten by formatting). The overall
  scale comes out right (the fitted ratio is of order $1/16\pi^5$, exactly what your
  prefactor conventions predict), but term-by-term signs cannot be fixed from the paste:
  e.g. your BFKL singlet form has the *same* sign for the
  $\mathcal K\mathcal K\,K_D\,\mathrm{Tr}$ term and the
  $K(w,w',z)K(x,y,z)\,\mathrm{Tr}_z$ term, while the pasted target has opposite signs
  — please check that one relative sign against your BFKL_JIMWLK_derivation.pdf, or upload
  the target as a .tex file and I will finish the comparison mechanically.

## 4. Finite terms — verification and simplification

Verified exactly where checked (§1.5). Given §2, a full re-organisation would be premature,
but two clean simplifications are safe now:

1. **Geometric closed form for $\mathcal{A}_0$** (exact):
   with $p=x-w$, $q=x-z$,
   $$\mathcal{A}_0=\frac{\theta}{|p\wedge q|},\qquad
   \theta=\angle\big(x-w,\;z-w\big)\in(0,\pi),$$
   i.e. $\mathcal{A}_0$ = (angle subtended at $w$ by the segment $xz$) divided by twice the
   area of the triangle $(x,w,z)$. Proof: $\arctan a+\arctan b$ with
   $a=\frac{q\cdot(q-p)}{|p\wedge q|}$, $b=\frac{p\cdot q}{|p\wedge q|}$ gives
   $\tan(\,\cdot\,)=\frac{|p\wedge q|}{(x-w)\cdot(z-w)}$ and
   $|p\wedge q|=|(x-w)\wedge(z-w)|$. This makes the finite $N_{3a},N_{3b}$ terms manifestly
   finite and rotation-invariant, and exposes the $w\to$ collinear limits.
2. The two red terms of $N_2$ in Group I Row III cancel identically — delete.
3. The $N_{2f}$ note (the $1/(z-w)^2$ coefficients $-(1-k^+/\vee)/k^+$ and
   $1/\vee-1/k^+$ cancelling) is algebraically exact: the two coefficients are equal and
   opposite identically, not just numerically.
4. Deeper simplification (combining the $\log|{\cdot}|$ structures of $N_{3a/b}$ with
   $\mathcal A_0$ into dilogarithm-free angle/log pairs, and merging the Row-IV/Row-I
   $\log[\,(1-k^+/\vee)A+(k^+/\vee)B\,]/B$ families across groups) is prepared in the
   scripts but should wait until §2 is resolved, since Row III's finite part may change.

## 5. Recommended order of fixes

1. Re-derive **Group I Row III's** five-piece $p^+$ master formula from the Ω-coefficient
   catalogue (this is the only place the $1/(z-w)^2$ 3ρ structure enters; everything hinges
   on it).
2. Rewrite the 4ρ master equation with one orientation convention (every $z$-line with the
   final-state index second), then re-run `fourrho.py` — its 4ρ sector must cancel, and its
   3ρ remainder must equal the "From four ρ's" blocks.
3. Re-check the GI-RowII collection against the four printed missing structures
   (`perrow.py` prints them).
4. Then the 3ρ cancellation and the JIMWLK extraction are proven end-to-end; the BFKL side
   follows with the same machinery (the scripts already parse your BFKL file).
