# Review of `Omega_at_g^2 (Only gluons) DRAFT` — located errors

Draft reviewed: 57 pp., "Soft Gluon Wave Function and Evolution Operator in the CGC at
Next-to-Leading Order", prepared for JHEP.

Every item below gives the equation number and the specific object. Items in **§A** I verified by
explicit re-derivation from the paper's own matrix elements (D.1)–(D.7) and its own conventions;
they are errors. Items in **§B** are internal contradictions between two places in the paper — one
of the two must change, and I say which I believe is right. **§C** is conceptual/argumentative.
**§D** is typography and prose.

Conventions I verified first, so that everything downstream rests on checked ground:

* **(D.1) is exactly right.** Computing `⟨g_i^a(k)|H_g|0⟩` from (2.20) and (B.2) gives
  `g kⁱρ^a(−k)/(4π^{3/2}|k⁺|^{3/2})`. ✓
* **(4.19) → (4.20) is right**, and **(4.20) → (5.2) is right**: matching (4.15)'s
  `(2π)^{3/2}∫d³p/(2π)³ A₂` against (4.20) gives exactly
  `A_j^b(p) = −√2 g p^j ρ^b(−p)/(√(p⁺) p²)`. ✓
* **(4.6) → (4.39) → (4.80) → (5.1) are mutually consistent.** `N² = 1 − ∫A†A` ⟹ `N = 1 − ½∫A†A`,
  and `½∫d³p/(2π)³ A†A = (g²/8π³)log(∨/Λ)∫d²p ρρ/p²`. ✓
* **(4.78) is right.** `∫_a^{1−a}dξ[ξ(1−ξ) + (1−ξ)/ξ + ξ/(1−ξ)] = 1/6 − 2ln a − 2 = −2log(Λ/k⁺) − 11/6`. ✓
* **(4.22) is right** — I initially thought the overall minus was spurious. It is not: matching
  (D.3)'s index assignment (`q ↔ l`, `p ↔ j`) onto (4.21)'s state `|g_j^b(k−p)g_k^c(p)⟩` flips the
  sign of all three tensor structures, which I checked term by term. `−igf^{abc}` is correct, and
  (4.36)'s `+igf^{acb}` is the same thing. ✓
* **(4.24) is right** (checked against D.4 including the `(2p⁺−k⁺)` factor). ✓
* **(4.53)–(4.54) are right**: the two orderings (4.50) and (4.52) combine with the decomposition
  (4.28), and the `(p⁺k²−k⁺p²)` denominators cancel exactly, leaving precisely (4.54). ✓
* **(6.9) is right**: `∫A†A p²/2p⁺ = ∫ g²ρρ/(p⁺)²`. ✓

---

## §A. Definite errors

### A1. Eq. (4.55) — missing factor `i` and missing `1/√(k⁺p⁺)`

This is the same calculation that produces (4.54), which is correct. Writing
`X = √(k⁺)/(√(p⁺)k²)`, `Y = √(p⁺)/(√(k⁺)p²)`, the sum of (4.50) and (4.52) is

```
C { X ρ^b(−p)ρ^a(k) − Y ρ^a(k)ρ^b(−p) },      C = (g²/4π³) kⁱpʲ /(p⁺k² − k⁺p²)
  = C { ½(X−Y){ρ,ρ} + ½(X+Y)[ρ^b(−p), ρ^a(k)] }
```

* The `{ρ,ρ}` piece gives `−(g²/8π³) kⁱpʲ/(√(k⁺p⁺)k²p²) {ρ^a(k),ρ^b(−p)}` = **(4.54) exactly**. ✓
* The `[ρ,ρ]` piece gives
  ```
  C·½(X+Y)·[ρ^b(−p),ρ^a(k)]  with  [ρ^b(−p),ρ^a(k)] = −i f^{abc} ρ^c(k−p)
  ==> + i (g² f^{abc}/8π³) · kⁱpʲ (p²k⁺ + p⁺k²) / ( √(k⁺p⁺) k²p² (p²k⁺ − p⁺k²) ) ρ^c(k−p)
  ```

Eq. (4.55) as printed is

```
−(g² f^{abc}/8π³) ∫dp⁺ ∫d²p |g_j^b(p)⟩ (kⁱpʲ/k²p²) ((p²k⁺+p⁺k²)/(p²k⁺−p⁺k²)) ρ^c(k−p)
```

so it is missing **(i)** the factor `i` — which is mandatory, since (4.29) itself writes the
commutator as `(i/2)f^{bad}ρ^d`; and **(ii)** the factor `1/√(k⁺p⁺)` — note (4.54) has
`∫dp⁺/√(k⁺p⁺)` while (4.55) has bare `∫dp⁺`, although both come from the same two terms.

### A2. Eq. (5.4) — A1 propagates, producing a dimensionally inconsistent bracket

(5.4) contains

```
− (g²/8π³) (kⁱpʲ/k²p²) [ (1/√(k⁺p⁺)) {ρ^a(k),ρ^b(−p)}  +  f^{abc}ρ^c(k−p) (p²k⁺+k²p⁺)/(p²k⁺−k²p⁺) ]
```

The two terms **inside the same bracket** differ by a factor `√(k⁺p⁺)` and by a factor `i`. This
is visible without any derivation: one term carries `1/√(k⁺p⁺)`, the other does not. The second
term should read `+ i f^{abc}ρ^c(k−p)(p²k⁺+k²p⁺)/(√(k⁺p⁺)(p²k⁺−k²p⁺))`.

### A3. Eq. (4.36) — the instantaneous contribution `|Ψ³_{ggρ}⟩` is missing

Eq. (4.35) states `|Ψ_{ggρ}⟩ = |Ψ^{1a}_{gg[ρ,ρ]}⟩ + |Ψ²_{ggρ}⟩ + |Ψ³_{ggρ}⟩` — three pieces.
Eq. (4.36) contains only two: the (4.31) term and the (4.22) term. The `|Ψ³_{ggρ}⟩` of (4.24) is
absent.

This is not a harmless omission, because **(5.3) does include it** — its last term,
`− ig²f^{abc}(p⁺−q⁺)ρ^a(−p−q)δ_jk / (2(2π)³√(p⁺q⁺)(p⁺+q⁺)²(q²/q⁺+p²/p⁺))`, is exactly (4.24)
rewritten with `k → p+q` (I checked the `(2p⁺−k⁺) → (p⁺−q⁺)` substitution and the sign; they
agree). So (4.36) is simply incomplete relative to what is actually used.

### A4. Eq. (5.3) — stray momentum `k` in a function of `(p,q)`

In the third group of (5.3):

```
[ 2pⁱ − (2p⁺/(p⁺+q⁺)) (pⁱ + kⁱ) ] δ_jk
```

`B₂^{cb}_{kj}(p,q)` depends only on `p` and `q`. Tracing back to (4.36) with `k → p+q` and
`p → q`, this term is `2pⁱ − (2p⁺/(p⁺+q⁺))(pⁱ+qⁱ)`. **`(pⁱ + kⁱ)` must be `(pⁱ + qⁱ)`.**
Confirmed by the fact that the *other* two structures in the same bracket, and the trailing factor
`gρ^a(−p−q)(pⁱ+qⁱ)/(4π^{3/2}|p⁺+q⁺|^{3/2})`, have already been correctly converted to `p,q`.

### A5. Eq. (5.3) — momentum/index assignment is swapped relative to Eq. (4.1)

Eq. (4.1) defines `B₂^{cb}_{kj}(p,q) a_j^{†b}(p) a_k^{†c}(q)`: index pair `(b,j)` belongs to
momentum `p`, and `(c,k)` to momentum `q`.

But every term of (5.3) uses the opposite assignment. The first two terms have `q^j p^k ρ^b(−q)ρ^c(−p)`
and `p^k q^j ρ^c(−p)ρ^b(−q)` — index `j` on momentum `q`, index `k` on momentum `p`, colour `b` on
`q`, colour `c` on `p`. The third group is the same (`q^j` sits with `δ_ik`, `p^k` with `δ_ij`).
Checked against the source (4.38), where the assignment *is* correct
(`p^k` with the gluon of momentum `p` and index `k`; `(k−p)^j` with the gluon `(k−p)` of index `j`).

So (5.3) is really `B₂^{cb}_{kj}(q,p)`. Either swap `p ↔ q` throughout (5.3), or write the
arguments in the other order. (This is invisible after symmetrisation, but as written the equation
does not match its own definition.)

### A6. Eq. (4.86) — two colour/momentum errors, both fixed by comparison with (4.88)

(4.86) and (4.88) are the two time orderings of the same process (incoming `r(d,l)`, `s(e,m)`;
outgoing `p(b,j)`, `q(c,k)`; `s` absorbed by the valence; `r → p+q`).

* **Colour factor.** (4.86) has `f^{dbe}`. The splitting `r(d) → p(b) + q(c)` must give `f^{dbc}`,
  and the valence absorption supplies `ρ^e`. As printed, `c` is a free index absent from the
  colour tensor and `e` appears twice. (4.88) correctly has `f^{dbc}`. **`f^{dbe} → f^{dbc}`.**
* **Argument of ρ.** (4.86) has `ρ^e(−s)`. The valence *absorbs* `s`, which by the paper's own
  convention (cf. (4.50), where absorbing `k` gives `ρ^a(k)`, and (4.44), where absorbing `k−p`
  gives `ρ^c(k−p)`) is `ρ^e(s)`. (4.88) correctly has `ρ^e(s)`. **`ρ^e(−s) → ρ^e(s)`.**

The rest of (4.86) — the `s^m/s²` from the `1/(s²/2s⁺)` denominator, the `1/√(p⁺q⁺r⁺s⁺)`, and the
three-term tensor structure — I checked and it is correct.

### A7. Eq. (4.94) — `g⁴` should be `g²`

`|(Ψ^{ba}_{ji}(k,p))³_{gg}⟩` is built from **two** `H_ggg` insertions (4.93), i.e. `(ig)² = −g²`.
The printed prefactor is `−g⁴ f^{ade}f^{bfh}/(128π³)`. Compare (4.84), which is the same
two-`H_ggg` topology and correctly carries `g² f^{d'ed}f^{d'cb}/(128π³)` — same `128π³`, and `g²`.
Also, this whole subsection is explicitly "at order `g²`". **`g⁴ → g²`.**

(The colour structure `f^{ade}f^{bfh}` with no shared index *is* correct here — the diagram on
p. 36 is two disconnected splittings, `k → r,s` and `p → u,v`.)

### A8. Eq. (6.4) — sign error; the printed bracket is anti-Hermitian and cannot cancel `H_g`

(6.4) prints

```
∫ d³p/(2π)³ [ − A_j^b(p) (p²/2p⁺) a_j^{†b}(p)  +  A_j^{†b}(p) (p²/2p⁺) a_j^b(p) ]
```

Two independent arguments show the relative sign is wrong:

1. **Hermiticity.** `(A E a†)† = A† E a`, so the two terms are Hermitian conjugates of one another
   and must enter with the *same* sign. As printed the bracket is anti-Hermitian, while `H_g` —
   which it is supposed to cancel — is Hermitian. An anti-Hermitian operator cannot cancel a
   Hermitian one, so (6.5) cannot follow from (6.4) as written.
2. **Direct computation.** With `Ω = 1 + gΩ⁽¹⁾`, `Ω⁽¹⁾ = A a† − A† a` (using (4.7)), and
   `Ω†H₀Ω = H₀ + [H₀, Ω⁽¹⁾]`:
   ```
   [H₀, A a† − A† a] = + A (p²/2p⁺) a†  +  A† (p²/2p⁺) a
   ```
   and with (5.2), `A_j^b(p)·p²/2p⁺ = − g p^j ρ^b(−p)/(√2 (p⁺)^{3/2})`, which is *minus* the `a†`
   coefficient of `H_g` in (2.20). So `H_g + [H₀,Ω⁽¹⁾]|_A = 0` **only** with `(+, +)`.

**The first term should be `+ A_j^b(p)(p²/2p⁺) a_j^{†b}(p)`.**

For the record, the `C`-part of (6.4) *is* correct: the h.c. of `C₁(E_r+E_q−E_p)a†a†a`, written in
the `a†aa` basis, carries `(E_p+E_q−E_r) = −(E_r−E_q−E_p)`, which reproduces the printed
`− C†(r²/2r⁺ − q²/2q⁺ − p²/2p⁺)` exactly. Only the `A` terms are wrong.

### A9. Eq. (2.23) — unbalanced brackets, duplicated prefactor, and a font break

As typeset, (2.23) opens `[` after the first prefactor/measure, then in the middle repeats
`+ (g²/8) f^{abc}f^{ade} ∫ (dp⁺/2π)(dk⁺/2π)(dr⁺/2π)dq⁺ ∫ … 1/√(p⁺k⁺r⁺q⁺)` and opens a **second**
`[`, but only **one** `]` ever closes. Two openings, one closing. Compare (2.24), which has the
same two-group structure and *is* balanced (`[…] + prefactor ∫ […]`). (2.23) is missing the closing
bracket before the repeated prefactor.

In addition, the last operator product of (2.23),
`a_i^{b†}(k)a_j^c(p)a_i^{d†}(r)a_j^{e†}(q)`, is typeset in **upright roman** while every other
operator in the equation is italic — a math-mode slip.

### A10. Eqs. (4.97) and (4.99) — inconsistent integration measure

Both are three-gluon outgoing wave functions from three-gluon incoming states, and (4.100) adds
them together. (4.97) integrates `∫d³p d³q d³u`; (4.99) integrates
`∫ (d³p/(2π)³)(d³q/(2π)³)(d³u/(2π)³)`. Their sources (4.96) and (4.98) both use the bare `d³`
measure. One of the two is wrong by `(2π)⁹`.

### A11. Duplicate bibliography entries

* **[22] and [34]** are the same paper: Altinoluk & Kovner, *Particle production at high energy and
  large transverse momentum: "The hybrid formalism" revisited*, PRD 83 (2011).
* **[47] and [50]** are the same paper: Munier, *Unitary perturbation theory on the light cone
  using adiabatic switching* ([50] additionally gives the eprint number 2510.05256).

---

## §B. Internal inconsistencies to resolve

### B1. Eq. (4.68) contradicts the paragraph immediately above it

The text on p. 29 states that "products mixing `g` and `g²` order components contribute only
beyond the perturbative accuracy considered here and can therefore be neglected." It then singles
out the overlap of `|(Ψ)_{gg}⟩` (order `g`, from `H_ggg`) with `|(Ψ)_{ggρ}⟩` (order `g²`) as "one
potentially contributing term", and argues it vanishes "due to the symmetry properties". That
overlap is `O(g³)` and is already excluded by the sentence above. Either drop (4.68) or explain
why the earlier rule does not apply.

### B2. Eq. (4.72)/(4.73) — one equation carrying two numbers, and broken cross-references

What is displayed is a single equation split across four lines, but it is given **two** numbers,
(4.72) on the first line and (4.73) on the last. The following text then says "the contributions
in the last line of Eq. (4.72)" and "The remaining contribution in Eq. (4.72)", but the last line
is (4.73). Also, this display writes `N₂` in a calligraphic font where everywhere else the symbol is blackboard-bold `ℕ₂`
(blackboard bold).

### B3. `B₃` is both matched and constrained, and the two are never reconciled

Eq. (4.10) determines `B₃†` from the `O(g²)` unitarity condition, while (5.4) determines `B₃`
independently by matching to the LCWF. The text below (5.12) says only "`B₁, D₃, E₃` are fixed by
the constraint equation" — but `B₃` appears in both places, so it is over-determined. That the two
determinations agree is a genuine, nontrivial check of the whole construction. Either state that
it has been performed, or say explicitly which route fixes `B₃` and what (4.10) is then used for.

### B4. Polarization index `k` collides with momentum `k`

(D.3) uses `l` for the polarization index of the outgoing gluon `q`, giving `δ_jl`, `δ_il`,
`δ_ij` and `q^l, p^l, k^l` — unambiguous. The main text renames it to `k` in (4.22), (4.36),
(4.46), (4.48), (4.57) and (5.5), producing expressions such as

```
( (k⁺+p⁺)/q⁺ q^k − p^k − k^k ) δ_ij           [Eq. (5.5)]
```

where `k^k` means "component `k` of the momentum `k`". This is not an error — I checked that the
contractions are consistent — but it is genuinely ambiguous on first reading, and it is
inconsistent with the appendix that the main text cites. Rename the polarization index to `l` (or
the momentum to something else) throughout.

### B5. Index labels change between (4.23) and (4.24)

(4.23) writes the state `|g_i^b(k−p) g_j^c(p)⟩` and picks up `δ_ij` from (D.4); (4.24) writes the
same state as `|g_j^b(k−p) g_k^c(p)⟩` with `δ_jk`. Harmless in isolation, but this is the kind of
relabelling that propagates.

### B6. Eq. (2.23) vs (2.24) — redundant and inconsistently labelled

(2.17) lists `H_gggg` and `H_gggg−inst.` as separate terms, and (C.6) defines only their **sum**.
(2.23) is labelled `H_gggg` alone, yet it contains momentum fractions such as
`(p⁺+k⁺)q⁺/((r⁺+q⁺)(k⁺−p⁺))` that can only come from the instantaneous `1/∂⁺` structures — the
genuine four-gluon contact vertex has no such factors. (2.24), labelled `H_gggg + H_gggg inst.`,
then covers overlapping ground. Given that (2.23) is also malformed (A9), my reading is that
(2.23) is a leftover that should be deleted, with (2.24) retained.

### B7. Non-normal-ordered operator products in (2.23)–(2.24)

Several products are written with annihilation operators to the left of creation operators, e.g.
`a_i^b(k) a_i^{c†}(p) a_j^{d†}(r) a_j^{e†}(q)` in (2.23), and
`a_i^b(p) a_i^{c†}(k) a_j^d(r) a_j^{e†}(q)` in (2.24). Since `[a,a†] ≠ 0`, these differ from their
normal-ordered forms by one-body terms. Given that normal ordering is the organising principle of
the entire paper, either normal-order these or state explicitly that the ordering as written is
intended.

### B8. Eq. (2.24), third term — squared longitudinal denominators

Terms 1 and 2 carry `1/((p⁺−k⁺)(r⁺−q⁺))`; term 3 carries `1/((p⁺+k⁺)²(r⁺+q⁺)²)`. Each `1/∂⁺` in
(C.6) supplies one power, so the squares look anomalous relative to the neighbouring terms. Worth
re-deriving.

---

## §C. Conceptual points a referee will raise

### C1. Eq. (6.7) → (6.8) is asserted, not shown — and it is the main result

The step from (6.6) to (6.7) drops `H_gggg`, `H_gg−inst.` and `H_gggg−inst.` with no intermediate
algebra, and (6.8) asserts that the `C†C` structure vanishes identically ("Evaluating the second
term explicitly, one finds … = 0") with no demonstration. This is the paper's central application.
At minimum it needs an appendix: the `O(g²)` transformation involves `[H₀,Ω⁽²⁾]`,
`[H_g+H_ggg, Ω⁽¹⁾]` and `½[[H₀,Ω⁽¹⁾],Ω⁽¹⁾]`, and a reader currently cannot verify that they have
all been included, let alone that they cancel.

### C2. The scaleless-integral argument on p. 30 needs justification in *this* scheme

`⟨(Ψ)_{gg}|(Ψ)_{gg}⟩` is set to zero because `∫d²p̃/p̃²` is scaleless, citing dimensional
regularisation [62]. But the calculation is not done in dimensional regularisation: the
longitudinal directions carry explicit cutoffs `Λ` and `∨`, and the same transverse integral
`∫d²p ρ^a(−p)ρ^a(p)/p²` is kept unintegrated elsewhere in the very same equation. Discarding an
overlapping UV/IR divergence by fiat, in a cutoff scheme, is a scheme choice that changes `ℕ`.
The UV half of it would ordinarily be absorbed into coupling renormalisation, not dropped. State
the scheme explicitly.

### C3. `E₁` and `E₂` are pure products of `C`'s — so why are they in the generator `G`?

(5.8) and (5.9) give `E₁` and `E₂` entirely as products `C·C` and `C·C†` over energy denominators,
i.e. as disconnected/reducible structures. The text introducing (5.11) says `G` is introduced "to
extract the connected contributions". But (5.12) then carries `E₁, E₂, E₃` terms explicitly. If
`G` is meant to hold only connected pieces, those products should cancel against the exponentiation
of the `C` terms and should not appear in `G`. Either the "connected" language or the content of
(5.12) needs adjusting.

### C4. `Λ`- and `∨`-dependence is raised and then dropped

§2.3 states correctly that "physical observables should ultimately be independent of the arbitrary
separation scale `∨`". The coefficients themselves are manifestly cutoff-dependent — (5.4) carries
`Θ(p⁺−k⁺−Λ)` and `Θ(k⁺−Λ−p⁺)`, and `ℕ` in (4.80) is `∝ log(∨/Λ)`. Nothing in the paper shows or
even sketches how the dependence cancels. This will be asked.

### C5. `(6.10)`: the surviving energy has no transverse propagator

`∫ d³p/(2π)³ g²ρ^b(p)ρ^b(−p)/(p⁺)²` — the `1/p²` of the Weizsäcker–Williams field has cancelled
completely against the light-cone energy `p²/2p⁺`. The arithmetic is right (I checked (6.9)), but
the `∫dp⁺/(p⁺)²` is *power* divergent at the lower cutoff, not logarithmic. Since the text calls
this "the coherent background field energy", one sentence on why a linear `1/Λ` sensitivity is the
expected answer would help the reader.

### C6. LCPT citations are thin

"light-cone perturbation theory (LCPT) [46, 47]" cites Kovchegov–Levin's book and a 2025 Munier
paper. The standard references (Lepage–Brodsky; Bjorken–Kogut–Soper; Brodsky–Pauli–Pinsky) are
absent. Also, "As discussed in Appendix A of Ref. [61]" points at a *Physics Letters B* paper
(Kovner–Lublinsky–Serino, PLB 792) — worth confirming that it has an Appendix A on adiabatic
switching, since a Munier reference would be the natural home for that statement.

Finally, in the concluding paragraph, "…requires the explicit coefficients of the evolution
operator derived here [63]" attaches [63] (Kovner–Lublinsky, *One gluon, two gluon*) to a clause
about coefficients derived *in this paper*. The citation is misplaced or belongs to a different
clause.

---

## §D. Typography and prose

**Hyphenated plurals** — "two-gluons", "three-gluons" used as plural nouns, and "a single-gluon"
used as a noun. These should be "two gluons", "three gluons", "a single gluon". Occurrences:
p. 17 ("production of two-gluons"), p. 18 ("emission of two-gluons"), p. 25 ("splitting into
two-gluons"), p. 26 ("splits into three-gluons"), p. 26 ("splits into two-gluons"), p. 27
("splits into two-gluons"), p. 28 ("first splits into two-gluons"), p. 31 ("containing
two-gluons"), p. 35 ("merge into a single-gluon"), p. 38 ("merge into a single-gluon").

**Other prose:**

| Location | Printed | Should be |
|---|---|---|
| p. 38 | "then an another incoming gluon splits" | "then another" |
| p. 30 | "the contributions … sums up to" | "sum up to" |
| p. 22 | "we get identical expression modulo the integration limits" | "an identical expression" |
| p. 48 (App. B) | "the following three eigen states" | "eigenstates" |
| p. 47 | "research assisstantship support" | "assistantship" |
| Fig. 1 caption | "the green region – Λ < p⁺ < ∨" | stray en-dash before `Λ` |

**Inconsistent hyphenation throughout** (pick one and apply globally):

* "Born–Oppenheimer" (abstract, p. 2, p. 5) vs "Born-Oppenheimer" (contents, §2.3 heading, p. 5,
  p. 46)
* "non-commuting" (abstract) vs "non commuting" (p. 11, p. 11 again, p. 46)
* "off diagonal" (abstract, p. 10) vs "off-diagonal" elsewhere
* "next-to-leading" vs "next to leading" (p. 32, p. 38, p. 50)
* "multi particle" (p. 40) vs "multi-gluon" (p. 1, p. 2)
* "Weizsäcker–Williams" vs "Weizsäcker-Williams" (p. 45, p. 46)

**Math typography:**

* (5.3), term 2: `p^k (q)^j` — the `q` is in plain italic and parenthesised, while every other
  transverse momentum in the paper is bold. Should be `p^k q^j`.
* (2.23), last term: entire operator product in upright roman (see A9).
* (2.23)/(2.24): `+ h.c.,` with a trailing comma inside the display, four times.
* (5.4): `ρ^c(−p+k)` on the first line vs `ρ^c(k−p)` on the third — same object, two notations.

---

## Priority order for fixing

1. **A8** (6.4 sign) — it breaks the paper's stated main result.
2. **A1 + A2** (4.55 and 5.4) — a coefficient of `Ω` is wrong by `i·√(k⁺p⁺)`.
3. **A3** (4.36 incomplete), **A4** (5.3 stray `k`), **A5** (5.3 argument order).
4. **A9 + B6** (2.23 malformed / redundant with 2.24).
5. **C1** (show the 6.6 → 6.8 algebra) — the largest referee-facing gap.
6. **A6, A7, A10, A11** — localized, quick.
7. **§B4** (index collision) and **§D** — presentation.
