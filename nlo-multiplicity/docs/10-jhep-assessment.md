# Addendum — Will arXiv:2607.18373 make it into JHEP?

**This replaces an earlier version written before the paper was available.** That version
estimated 45–55% and rested on three objections. Having now read v2 (57 pp., 28 Jul 2026), two of
those objections are wrong and one is weaker than I stated. The estimate goes **up**.

The paper is already formatted "Prepared for submission to JHEP", so the venue question is really
"will it survive refereeing there", not "is JHEP the right target". It is.

---

## The short answer

**~65–75%, most likely after one round of moderate revision.**

The upgrade over the blind estimate comes from four things I could not see before:

1. **The novelty claim against the prior literature is specific and correct.** The paper does not
   claim to compute the wave function for the first time. It claims that Ref. [19]
   (Lublinsky–Mulian) obtained *some* `O(g²)` and `O(g³)` contributions but "did not contain all
   contributions at order `g²`", and that the present coefficients "complete and extend" them.
   That is a well-defined, checkable claim, and it is the right one. My earlier objection — "the
   wave function is already known to higher order, so why do `O(g²)`?" — **was wrong**: partial
   knowledge at a higher order does not imply complete knowledge at a lower one.

2. **The `O(g)` trilinear coefficient `C` makes the diagonalization a real result, not a check.**
   I assumed §6 would reproduce the standard coherent-state statement. It does not. The paper
   shows that the coherent operator used in Ref. [61] is *insufficient* at order `g`: the
   three-gluon vertex generates an operator `C` (structure `a†a†a`, order `g`, carrying no `ρ`)
   at the same order, and only when `C` is included do the off-diagonal elements cancel
   identically, giving `Ω†H₁Ω = H₀` exactly (Eq. 6.5). At `O(g²)` the `C·C` contributions cancel
   (Eq. 6.8) and the sole survivor is the `ρ²` coherent background-field energy (Eq. 6.10). That
   is a correction to received wisdom, cleanly stated. It is still not a physical observable — see
   below — but it is more than internal validation.

3. **The volume and explicitness of the results.** Six coefficient families, given in closed form,
   with the `D` and `E` families expressed through `A` and `C` where possible, plus `Ω = exp(iG)`
   with `G` Hermitian by construction (Eqs. 5.11–5.12). JHEP rewards complete, reusable formula
   sets; this is one.

4. **Reception context.** DOE-funded, NC State, acknowledging Lublinsky and Skokov — i.e. the
   people whose line of work this extends — and a named independent cross-check of critical parts.
   None of that is a scientific argument, but it does mean the work has been seen by the community
   that supplies the referees, and that materially lowers the chance of a hostile report.

**What holds it below ~80%** is one surviving objection and three fixable presentation problems.

---

## The one objection that survives

**No physical observable is computed.** The paper's own conclusion is candid about this: the
construction "constitutes a necessary step towards the computation of single inclusive gluon
production at NLO", with double and triple gluon production and the multiplicity moments left to
future work. So the "so what" is entirely deferred.

For JHEP this is survivable but not free. The standard form of the report is: *"a careful and
useful technical construction, but the physics payoff is postponed; I would like to see at least
one worked consequence."* Two ways to answer it, in order of cost:

* **Cheap (recommended).** Sharpen the argument for *why the normal-ordered form matters*, which
  the introduction currently gestures at rather than makes. The point is precise and strong: the
  existing NLO literature extracts the **Hamiltonian**, i.e. the diagonal matrix element of the
  evolution operator. That suffices for evolution equations and inclusive cross sections. It does
  **not** suffice for anything requiring off-diagonal Fock structure — multiplicity distributions,
  factorial cumulants, the reduced density matrix, entanglement measures. A normal-ordered `Ω`
  gives those; a Hamiltonian cannot. One paragraph, and it converts "no result" into "a tool that
  unlocks a class of results".
* **Expensive.** Add a worked consequence. I earlier suggested holding the paper to add the
  `O(g⁴)` variance. **I retract that.** At 57 pages with a clearly staged programme already
  announced, splitting is correct; adding two more pages of a different calculation would dilute
  the paper without changing the referee's basic reading.

## Three fixable problems that cost probability

### (a) The delineation against Ref. [19] is asserted, not shown — **the highest-value fix**

The claim "[19] did not contain all contributions at order `g²`" appears once in the
introduction and once in the conclusion. It is load-bearing for the entire novelty case, and the
referee is asked to take it on trust. A referee who is an author of [19], or who has read it
recently, will want specifics.

**Fix: a short table.** One column per coefficient family (`N`, `A₁₂`, `B₁₂₃`, `C₁₂`, `D₁₂₃`,
`E₁₂₃`), one column saying "new / present in [19] / partially in [19]". Half a page. In my
judgement this single addition is worth more than anything else on this list — it converts the
central claim from an assertion into something a referee can verify in five minutes.

### (b) Cutoff and scheme dependence is raised and then dropped

§2.3 states correctly that "physical observables should ultimately be independent of the arbitrary
separation scale `∨`". But the coefficients themselves are explicitly cutoff-dependent — Eq. (5.4)
carries `Θ(p⁺ − k⁺ − Λ)` and `Θ(k⁺ − Λ − p⁺)`, and `N` in Eq. (4.80) is `∝ log(∨/Λ)` — and nothing
demonstrates or even sketches the `∨`-independence. This will be asked. It does not need a full
proof; a paragraph identifying where the `∨`-dependence sits and how it is expected to cancel
against the valence sector would defuse it. Leaving it silent invites a referee to treat it as an
unexamined gap.

### (c) Presentation defects in a paper whose value is algebraic accuracy

See `../review/draft-review-omega-g2.md` for the full itemized list. The ones that matter most for
refereeing:

* **Eq. (6.4) has a sign error** — the two `A`-terms are printed with opposite signs, making the
  bracket anti-Hermitian, so it cannot cancel the Hermitian `H_g` and Eq. (6.5) does not follow.
  This is in the section that carries the paper's application.
* **Eq. (4.55)** is missing a factor `i` and a factor `1/√(k⁺p⁺)`; the error propagates into
  **Eq. (5.4)**, where two terms inside one bracket then differ dimensionally.
* **Eq. (4.36)** omits the instantaneous contribution that Eq. (4.35) promises and Eq. (5.3)
  actually uses.
* **Eq. (2.23)** has unbalanced brackets and one operator product in upright roman.
* The step **(6.6) → (6.7) → (6.8)** is asserted with no intermediate algebra, and (6.8) — that
  the `C†C` structure vanishes identically — is stated without demonstration. This is the largest
  referee-facing gap in the paper.
* **Consolidate the internal checks.** The paper contains genuinely strong validation: unitarity
  imposed order by order, matching across four independent Fock sectors (vacuum, one-, two-,
  three-gluon), and the two diagonalization results. But these are scattered across 40 pages. A
  short "consistency checks" subsection listing them together would materially change how a
  referee experiences 57 pages of algebra — from "I cannot verify this" to "the author has
  verified this in four independent ways."

## What is *not* a problem

* **Length.** JHEP has no limit and is the right home for a long coefficient derivation.
* **`N_f = 0`.** Standard for a first paper in a programme, and explicitly justified (gluonic
  dynamics dominate at small `x`). Retract as a concern.
* **The "NLO" terminology.** I previously flagged a mismatch between "NLO" in the title and
  `O(g²)`. **Retracted.** The abstract says "up to `O(g²)`" in its second line, and the usage
  matches the established convention in this specific line of work — cf. the title of
  Ref. [19], *High Energy QCD at NLO: from light-cone wave function to JIMWLK evolution*. No
  referee in this community will be misled.
* **Single author.** Irrelevant.

---

## Recommendation

| Action | Effect on acceptance | Cost |
|---|---|---|
| Fix Eq. (6.4), (4.55)/(5.4), (4.36), (2.23) | **prerequisite** — these are errors, not polish | days |
| Add the coefficient-by-coefficient comparison table vs Ref. [19] | **largest single gain** | half a day |
| Add a "consistency checks" subsection consolidating the four validations | high | half a day |
| One paragraph: why the normal-ordered form gives what a Hamiltonian cannot | high | an hour |
| Show the (6.6) → (6.8) algebra in an appendix | high | 1–2 weeks |
| A paragraph on `Λ`/`∨`-dependence and expected cancellation | medium | a day |
| Hold the paper to add the `O(g⁴)` variance | **negative** — do not | months |

With the errors fixed and the first three additions done, I would put it at **~80%**. As it
stands, **~65–75%**.

**The most likely bad outcome is not rejection — it is a referee who asks for the Ref. [19]
delineation and the scheme-dependence discussion, and the paper spends four months in revision
that a day of work now would avoid.**

## One thing worth knowing for the follow-up

Reading the operator changes the plan in Part IX in one respect. `Ω|0⟩` is exactly
`N|0⟩ + ∫A₂a†|0⟩ + ∫∫B₂a†a†|0⟩` through `O(g²)` — the `C`, `D`, `E` families annihilate the soft
vacuum. So the multiplicity-moments paper needs **only Eqs. (5.1)–(5.3)**, which is very good news
for its feasibility, and slightly awkward news for its framing: a referee may observe that it uses
three of the six coefficient families.

The natural answer is Project E: after eikonal scattering `a → U a`, the annihilation operators no
longer kill the state and `C`, `D`, `E` all contribute. **The produced-gluon observables are what
the rest of this paper is for**, and saying so in the concluding remarks would strengthen both
papers at once.
