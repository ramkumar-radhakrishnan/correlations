# Part VIII — Publication potential: a critical assessment

You asked for honesty rather than encouragement. Here is my actual read.

## 8.0 The blunt summary

**There is one good paper here, not three.** The good paper is:

> *Non-Poissonian gluon statistics from the NLO CGC wave function: factorial cumulants, the
> negative-binomial parameter, and away-side pair correlations.*

It is a **PRD paper**. It is not a JHEP paper as currently scoped, and it is more than a short
note. Getting it to JHEP requires adding produced-particle observables (target Wilson lines), and
that roughly doubles the work.

The three things that make it publishable rather than routine:
1. `C₂` is **complete at `O(g⁴)`** with the operator that already exists — a finished, checkable
   statement, not a partial result.
2. The **displaced-squeezed-state identification** makes the result legible and quotable outside
   the small-`x` community.
3. The **away-side vs near-side discriminator** gives it a connection to the ridge literature,
   which is where the citations are.

The three things that could sink it:
1. The new effect may be numerically negligible against classical colour-charge fluctuations
   (test P3 in Part VI — run it first).
2. Gluon number in a light-cone wave function is not an observable (Part III §3.5).
3. arXiv:2605.03038 may already contain the squeezed-state framing.

---

## 8.1 Ranked projects

### **Project A — `O(g⁴)` factorial cumulants and the NBD parameter** ★★★★☆

*Compute `⟨N⟩` to `O(g⁴)` and `C₂` complete at `O(g⁴)`; evaluate in MV; extract `k_NBD` and its
quantum correction; demonstrate JIMWLK consistency of the `Y`-dependence.*

**Novelty.** Real but incremental. The wave function exists (yours); the observable is new; the
quantum-optics machinery is imported. Honest characterization: *first computation of a
Fock-space fluctuation observable from an NLO CGC evolution operator.* That is a legitimate
"first", and the completeness of `C₂` at `O(g⁴)` is what makes it a paper rather than a note.

**Technical difficulty.** Moderate. The Gaussian reduction (Part I §1.5) makes the Fock algebra
trivial — and it needs only Eqs. (5.1)–(5.3) of the paper. The real work is (i) the `ρ`-ordering
bookkeeping, (ii) the colour algebra of `⟨ρ⁴⟩`, (iii) the IR analysis. **3–5 months of focused
work** for someone who already owns the `Ω` derivation.

**Expected referee objections**, with the fix:

| Objection | Severity | Fix |
|---|---|---|
| "Gluon number in the LCWF is not an observable" | **High** | Show the `Y`- and scheme-dependence is JIMWLK; frame as a wave-function property analogous to the ugd. Mandatory |
| "Your effect is swamped by classical colour fluctuations" | **High** | Part V Model 7: quantify the ratio explicitly; use the angular discriminator |
| "The generating function is textbook quantum optics" | Medium | Cite it as such; the contribution is the QCD kernels, not the Gaussian algebra. Do not oversell |
| "How does this differ from Kovner–Lublinsky multi-gluon production?" | Medium | Prepare a one-paragraph answer. Theirs: production cross sections at LO/RFT. Yours: wave-function Fock statistics at NLO |
| "`C₂` is IR-divergent / cutoff-dependent" | Medium | Quote `k_⊥`-differentially; show which combinations are finite |
| "The `ρ`-commutator terms are dropped/mishandled" | Medium | Appendix demonstrating ordering-independence |
| "You use only three of the six coefficient families" | Medium | True, and unavoidable: `C, D, E` annihilate `\|0⟩`. Say so, and point to Project E as their consumer |
| "Why not `C₃`, `C₄`?" | Low | Because the `O(g³)` operator is needed. Say so in the text — pre-empting this converts a weakness into evidence of care |

**Additional calculations required.** The `ρ`-ordering appendix; the JIMWLK consistency check; the
IR analysis. All are within the existing formalism.

**Is numerical work sufficient?** No. A purely numerical paper here would be rejected or heavily
criticized — the community expects the analytic `O(g⁴)` expression. **You need the closed-form
`C₂` in terms of `A₂` and `B₂` with the colour algebra done, plus MV numerics.** Numerics alone is
a note; analytics alone is publishable but thin. Both together is a PRD.

**Target: PRD.** Realistic. **Verdict: do this one.**

---

### **Project B — Extend `Ω` to `O(g³)` and compute `C₃`, `C₄` completely** ★★★☆☆

**Novelty.** High if completed — nobody has the complete third factorial cumulant.

**Technical difficulty.** **Severe.** You must compute the three-gluon component of `Ω|0⟩` at
`O(g³)`, which is the genuine two-loop-level LCPT calculation that Kovner–Lublinsky–Mulian did for
the Hamiltonian. Reproducing it in normal-ordered Fock form is a substantial project on its own.
**12–24 months.**

**Expected objections.** "You have redone KLM in a different basis." That objection is fatal
unless you can show the Fock-space form contains information the Hamiltonian does not — which it
does (the Hamiltonian is the diagonal matrix element; the wave function has off-diagonal
structure), but you must argue it explicitly.

**Verdict: this is a *second* paper, and only after Project A lands.** Do not attempt it first.
The risk of spending a year and being told "this is KLM" is real.

---

### **Project C — Away-side two-gluon correlation from the squeezing kernel** ★★★★☆

*Compute `Γ(k,q)` at `O(g⁴)`, decompose into near-side (Bose) and away-side (pair-emission)
components, predict the `Δφ` structure.*

**Novelty.** Potentially the highest in the whole program — **if** it is not already implicit in
the glasma-graph literature. Verify against arXiv:1503.07126 and arXiv:1804.02910 first.

**Technical difficulty.** Low-to-moderate *given Project A* — it is the same kernels, differently
projected. Essentially free once `A₂` and `B₂` are in hand.

**Expected objections.** "Away-side correlations from momentum conservation are trivial/background."
This is the serious one. The defence must be quantitative: the *normalization*, the `Q_s`-scaling,
and the fact that this particular away-side term is `O(g⁴)` while the near-side Bose term is
`O(g⁶)` — i.e. at NLO the away side should *dominate*, which is not what a trivial
momentum-conservation argument would give.

**Verdict: fold this into Project A as its second half.** Do not publish separately — alone it is
thin, but as the physics payoff of Project A it is what makes the paper interesting.
**This is what turns Project A from ★★★ to ★★★★.**

---

### **Project D — `O(g³)` odderon contribution to the gluon multiplicity** ★★★☆☆

*Compute `2Re∫A₂^{(1)*}A₂^{(2)}`, requiring `⟨ρρρ⟩ ∝ d^{abc}`; show it vanishes in MV and evaluate
in a C-odd model.*

**Novelty.** Genuinely unexplored, and it connects to the active odderon programme.

**Technical difficulty.** Low. This is a contained calculation — **2–3 months**.

**Expected objections.** "The odderon contribution to an unobservable quantity is doubly
unobservable." Fair and hard to answer. Also: "the C-odd correlator is model-dependent."

**Verdict: a good *letter* or a section of Project A, not a standalone paper.** As a standalone it
is a short note in PRD or EPJC at best. As a section of Project A it adds genuine value at low
cost. **Recommend folding in.**

---

### **Project E — Produced-gluon cumulants (with target Wilson lines)** ★★★★★ (impact) / ★★☆☆☆ (feasibility)

*Scatter the NLO dressed state off a CGC target; compute cumulants of the produced gluon
multiplicity.*

**Novelty.** High, and it solves the observability problem outright. **It is also the natural
consumer of the paper's `C`, `D`, `E` coefficients**, which play no role in `Ω|0⟩` but contribute
the moment `a → U a`.

**Technical difficulty.** **High.** You need `⟨U U U U⟩` (quadrupole) and higher target correlators
for `C₂`, `⟨U⁶⟩` for `C₃`. The Wilson-line correlators are known in the Gaussian/MV approximation
but the bookkeeping is heavy, and the NLO real/virtual cancellation in the *production* amplitude
is a known hard problem (the NLO `pA` production literature is littered with negativity issues and
scheme subtleties).

**Expected objections.** All the standard NLO-production objections — rapidity subtraction scheme,
negative cross sections, kinematic constraints. This is a mature and contentious subfield.

**Verdict: this is the JHEP version, and it is a 2-year project.** High reward, high risk. **Do not
start here.** Reassess after Project A. Note the paper's own conclusion already announces single
inclusive gluon production at NLO as work in progress (Ref. [53]), so coordinate.

---

### **Project F — KNO scaling and its violation from the NLO wave function** ★★☆☆☆

**Novelty.** Moderate. Crowded field (dipole cascade, entanglement, nearly-Gaussian action).

**Difficulty.** Low given Project A — it is a replot.

**Objections.** "KNO in small-`x` QCD has been derived several ways already; what does the wave
function add?" Hard to answer convincingly, because the dipole-cascade derivations are closer to
the measured quantity.

**Verdict: a figure and a paragraph in Project A. Not a paper.**

---

### **Project G — Purely numerical study of `P(N)` in MV** ★☆☆☆☆

**Verdict: not publishable alone.** Without the analytic `O(g⁴)` result this is a Monte Carlo
exercise reproducing a known NBD. A referee will say exactly that. **Do not write this paper.**

---

## 8.2 Ranking

| Rank | Project | Target | Time | Risk | Recommendation |
|---|---|---|---|---|---|
| **1** | **A + C + D combined** | **PRD** | 6–9 mo | Medium | **Do this.** One substantial paper |
| 2 | E (produced gluons) | JHEP / PRD | 18–24 mo | High | Second paper, after A |
| 3 | B (`O(g³)` operator, `C₃`) | JHEP / PRD | 12–24 mo | High | Third paper, or never |
| 4 | F (KNO) | — | — | — | Section of A |
| 5 | G (numerics only) | — | — | — | **Do not** |

## 8.3 The honest bottom line

**Project A+C+D is publishable in PRD with roughly 70% confidence**, conditional on the go/no-go
test (Part VI, P3) coming out favourably. If the quantum contribution to `c₂` turns out to be
`< 1%` of the classical one across all realistic `Q_s` and `k_⊥`, then the paper becomes "we
computed a correction that does not matter", which is publishable but low-impact, and I would
downgrade to ★★☆☆☆ and suggest reframing entirely around the away-side angular structure (where
the classical contribution is absent, so a small absolute size does not matter).

**It is not a JHEP paper as scoped.** JHEP in this subfield expects either a substantially new
formalism or an observable-level calculation. Project A is neither — it is a clean, complete,
first-of-its-kind observable extracted from an existing formalism. That is PRD.

**It is more than a short note**, provided the analytic `O(g⁴)` expression, the colour algebra, the
IR analysis, the JIMWLK check, and MV numerics are all present. Drop any two of those and it
becomes a note.

**The largest single risk is not technical — it is scooping.** arXiv:2412.05085 (Kovner et al.,
Born–Oppenheimer RG) and arXiv:2605.03038 (squeezed-state radiation) are both recent and both
adjacent. The group that wrote the NLO JIMWLK papers has every tool needed to do Project A in a
few months if they decide to. **Read those two papers this week, and if the ground is clear, move
fast.**
