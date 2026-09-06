# Part VII — Literature survey: what exists, what does not

**Method note.** This survey was assembled from bibliographic search only (full texts were not
retrievable in the session where it was written). Titles, authors, and identifiers are reliable;
my characterization of each paper's *detailed* content is inference from title/abstract/known
context and should be verified by reading before any of it is repeated in a publication. Items
marked **(!) READ THIS FIRST** are ones where the overlap risk is high enough that a wrong guess
would waste months.

---

## 7.1 The direct lineage — NLO light-cone wave functions and evolution operators

| Work | Content | Relation to this program |
|---|---|---|
| Kovner, Lublinsky, Mulian, *JIMWLK at NLO* (2013–14) | The `O(g³)`/two-loop soft wave function; extraction of `H_JIMWLK^NLO` from the diagonal `S`-matrix element | **The direct ancestor.** They compute the soft wave function to the order *beyond* arXiv:2607.18373 but extract only the *Hamiltonian*, not Fock-space observables |
| [Lublinsky & Mulian, arXiv:1610.03453](https://arxiv.org/abs/1610.03453), *High Energy QCD at NLO: from light-cone wave function to JIMWLK evolution* — this is Ref. [19] of the paper | Soft LCWF to third order in `g`; Fock space of one-gluon, two-gluon and `qq̄` states | **(!) READ THIS FIRST.** The paper's own novelty claim is that [19] "did not contain all contributions at order `g²`". That claim is the load-bearing one; see `10-jhep-assessment.md` §3.1 |
| [Kovner et al., arXiv:2412.05085](https://arxiv.org/abs/2412.05085), *Born–Oppenheimer RG for High Energy Scattering: the Setup and the Wave Function* | Same Born–Oppenheimer framing as arXiv:2607.18373 | Nearest methodological neighbour; check how much of the `O(g²)` operator is already there |
| Kovner & Lublinsky, [hep-ph/0609227](https://arxiv.org/abs/hep-ph/0609227) *One gluon, two gluon: multigluon production via high energy evolution* (Ref. [63] of the paper); [hep-ph/0608258](https://arxiv.org/abs/hep-ph/0608258) *Treading on the cut*; [arXiv:0901.2560](https://arxiv.org/abs/0901.2560) *Inclusive Gluon Production in the QCD Reggeon Field Theory* | Multi-gluon **production** from the high-energy evolution operator, including Pomeron loops | **(!) READ THIS FIRST.** This is the closest existing work to "moments of the multiplicity operator from an evolution operator". It is at LO/RFT level rather than NLO-LCWF level, but a referee will absolutely ask how your work differs. You must be able to answer in one sentence |
| [Armesto, Domínguez, Kovner, Lublinsky, Skokov, JHEP 05 (2019) 025](https://link.springer.com/article/10.1007/JHEP05(2019)025), *The CGC density matrix: Lindblad evolution, entanglement entropy and Wigner functional* | Reduced density matrix of the soft sector; entanglement entropy; Wigner functional | **Highly relevant.** They already treat the soft gluon sector as a quantum state with a density matrix. Multiplicity moments are a natural observable of that density matrix — check whether they computed any |
| [arXiv:2002.02282](https://arxiv.org/abs/2002.02282), *JIMWLK Evolution, Lindblad Equation and Quantum–Classical Correspondence* | Decoherence framing of JIMWLK | Supports the "quantum vs classical fluctuation" decomposition of §1.8 |
| Kovner, Lublinsky, Serino, *Entanglement entropy, entropy production and time evolution in high energy QCD*, PLB 792 (2019) — Ref. [61] of the paper | Coherent-operator diagonalization at `O(g)` | The paper's §6 extends this by adding the `O(g)` three-gluon operator `C` |

**Assessment of this block.** The `O(g²)`/`O(g³)` wave function is *not* virgin territory. The gap
is specific and real: **nobody has taken these wave functions and computed the factorial cumulants
of the gluon number operator.** The formalism exists; the observable does not.

## 7.2 Multiplicity distributions in the CGC/Glasma — the classical result

| Work | Content | Status |
|---|---|---|
| Gelis, Lappi, McLerran, *Glittering Glasmas*, Nucl. Phys. A 828 (2009) | `P(N)` from Glasma flux tubes = **negative binomial**, with `k ∝ (N_c²−1) Q_s² S_⊥ / 2π` | **The established result.** Classical order. Your `c₂` at LO must reproduce `1/k` |
| Lappi and collaborators, various | NBD parameter from CGC, comparison to RHIC/LHC | Established |
| [arXiv:1209.4105](https://arxiv.org/abs/1209.4105), *KNO scaling from a nearly Gaussian action for small-x gluons* | KNO emerges if the effective action for colour charges at `~Q_s` is nearly Gaussian; needs saturation + running coupling | Directly relevant to Part IV §4.6 |
| [IP-Jazma, arXiv:1808.01276](https://arxiv.org/abs/1808.01276) | Critical assessment of saturation explanations of small-system collectivity | Useful for calibrating claims about phenomenological relevance |

**Assessment.** The classical/LO multiplicity distribution of the CGC is **thoroughly done**. Any
paper that produces an NBD from a CGC calculation and stops there is not publishable. The novelty
must be the `O(g⁴)` correction, the sign structure, or the angular discriminator.

## 7.3 Dipole models, KNO, and entanglement

| Work | Content |
|---|---|
| Mueller (dipole cascade); [Le, Mueller, Munier, PRD 104 (2021) 056025](https://link.aps.org/doi/10.1103/PhysRevD.104.056025) | Multiplicity distribution of dipoles from the Le–Mueller–Munier equation; generating-function formalism |
| [Kharzeev & Levin, PRD 95 (2017) 114008](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.95.114008); [arXiv:2110.04881](https://arxiv.org/abs/2110.04881) | DIS entanglement entropy `S = ln N(x,Q²)`; maximally entangled proton; `P_n` geometric |
| [arXiv:2302.01380](https://arxiv.org/abs/2302.01380), *Universality of KNO scaling in QCD at high energy and entanglement* | KNO from entanglement; comparison to H1 and ALEPH |
| [arXiv:2406.04856](https://arxiv.org/abs/2406.04856) | KNO scaling emergence in LHC jets |

**Assessment.** Generating functions and multiplicity distributions in small-`x` QCD are a *busy*
field, but almost entirely in the **dipole-cascade / partonic-counting** language, where the
"multiplicity" is a number of dipoles evolving by a Markov branching process. That is a
fundamentally different object from the Fock-space gluon number of a light-cone wave function.
**This is a genuine gap and a genuine opportunity**, but it is also a hazard: a referee from that
community will ask how your `P(N)` relates to the Mueller/Kharzeev–Levin `P_n`. Have an answer.
(Mine: theirs is a classical branching process with geometric statistics; yours is a Gaussian
quantum state with Poisson statistics at LO — they are *different* distributions and the
difference is precisely the coherence of the emitted field. Saying this clearly is itself a
contribution.)

## 7.4 Correlations, Bose enhancement, the ridge

| Work | Content | Overlap |
|---|---|---|
| [Altinoluk, Armesto, Beuf, Kovner, Lublinsky, arXiv:1503.07126](https://arxiv.org/abs/1503.07126), *Bose enhancement and the ridge* | Bose enhancement of gluons in the projectile WF → azimuthal collimation; correlations suppressed by `1/(S_⊥ × N_gluons)` | **High.** This is the `α*nα` term of §4.2. **(!) READ THIS FIRST** |
| [arXiv:1804.02910](https://arxiv.org/abs/1804.02910), *Correlations and the ridge in the CGC beyond the glasma graph approximation* | Beyond glasma graphs | High |
| [arXiv:1612.07790](https://arxiv.org/abs/1612.07790), *Exploring correlations in the CGC wave function: odd azimuthal anisotropy* | Odd harmonics from the CGC wave function | Relevant to the odderon/`⟨ρρρ⟩` project |
| [arXiv:1808.04982](https://arxiv.org/abs/1808.04982), *Energy evolution and Bose–Einstein enhancement for double parton densities* | Evolution of the Bose-enhanced double-parton density | Relevant to `Y`-dependence of `c₂` |

**Assessment.** The **near-side** (`k ≈ +q`) correlation in the projectile wave function is
well-established and thoroughly explored. **The away-side (`k ≈ −q`) correlation from the
pair-emission/squeezing kernel `B₂` of the NLO wave function is, as far as I can determine, not
computed anywhere.** That is the sharpest open item in the survey.

Caveat to that claim: away-side correlations from momentum conservation appear in many guises in
this literature, and it is possible that the `α*α*m` structure appears implicitly in a
glasma-graph calculation without being identified as a squeezing effect. **Verify by reading
1503.07126 and 1804.02910 carefully before claiming novelty.**

## 7.5 Coherent states, squeezed states, quantum optics

| Work | Content | Overlap |
|---|---|---|
| [arXiv:2605.03038](https://arxiv.org/abs/2605.03038), *Squeezed-state radiation in shockwave scattering: QCD–Gravity double copy* (2026) | Squeezed-state radiation from shockwave scattering | **(!) HIGHEST-RISK ITEM IN THE SURVEY.** Recent, directly on the squeezed-state-in-QCD theme. Read before doing anything else |
| [hep-ph/9604371](https://arxiv.org/abs/hep-ph/9604371), *Squeezed States and Particle Production in High Energy Collisions* | Squeezed-state phenomenology of multiplicity distributions | Old, phenomenological, not CGC. Low overlap but must be cited |
| [hep-ph/9902402](https://arxiv.org/abs/hep-ph/9902402), *Gluon Squeezed States in QCD Jet* | Squeezed gluon states in jets | Low overlap |
| Quantum-optics standard results | Displaced squeezed thermal states, photon-number distributions, Fredholm determinants | **This is where §2.2 comes from.** Cite it properly — the generating function is a known quantum-optics result and presenting it as new would be an easy referee kill |

**Assessment.** The coherent-state picture of the classical CGC field is folklore and is stated in
many places (including the paper's own Eq. (3.8)). The **squeezed**-state picture — that `Ω|0⟩` is
a *displaced squeezed* state and that this is what generates non-Poissonian statistics — is, in
the CGC context, either new or at most touched on in arXiv:2605.03038. **This framing is the most
marketable aspect of the project**, precisely because it makes the result legible to people
outside small-`x` QCD. It is also the aspect most vulnerable to being scooped.

## 7.6 Summary: known vs open

**Known (do not claim as new):**
1. LO CGC soft gluon state = coherent state; `dN/d²k` = WW field.
2. `P(N)` from classical CGC/Glasma = negative binomial; `k ∝ (N_c²−1)Q_s²S_⊥`.
3. Near-side Bose enhancement and HBT in the projectile wave function → the ridge.
4. NLO JIMWLK Hamiltonian and the `O(g³)` soft wave function (Kovner–Lublinsky–Mulian).
5. Multiplicity generating functions and KNO in the dipole-cascade language.
6. Entanglement-entropy-based `P_n` (Kharzeev–Levin) and its data comparisons.
7. Displaced-squeezed-state photon statistics (quantum optics; import, do not rederive).
8. CGC reduced density matrix, Lindblad/decoherence framing of JIMWLK.

**Open (candidate contributions):**
1. **Factorial cumulants of the Fock-space gluon number operator from an NLO CGC evolution
   operator.** Not done anywhere I can find. The formalism exists; the observable does not.
2. **The `O(g⁴)` quantum correction to the negative-binomial parameter `k`.** Not done.
3. **The identification of `Ω|0⟩` as a displaced *squeezed* state, and of the
   displacement–squeezing interference as the leading source of non-Poissonian statistics.**
   Possibly new; check arXiv:2605.03038.
4. **The away-side (`k ≈ −q`) two-gluon correlation from `B₂`**, and its contrast with near-side
   Bose enhancement. Sharpest open item; verify against 1503.07126.
5. **A possible sub-Poissonian (`C₂ < 0`) region** — non-classical gluon statistics in QCD. Would
   be genuinely striking; nobody has looked.
6. **The `O(g³)` odderon contribution `⟨ρρρ⟩` to the gluon multiplicity.** Small, clean, unexplored.
7. **KNO scaling and its `α_s`-suppressed violation from the NLO wave function** (as opposed to
   from the dipole cascade). Open, but crowded neighbourhood.

**Closed by prior art, do not pursue:**
- Reproducing the NBD from CGC.
- Near-side ridge from Bose enhancement.
- Three-gluon correlations (both incomplete here *and* partly covered by Kovner–Lublinsky
  multi-gluon production).

---

**Sources:**
[arXiv:2607.18373](https://arxiv.org/abs/2607.18373) ·
[arXiv:1610.03453](https://arxiv.org/abs/1610.03453) ·
[arXiv:2412.05085](https://arxiv.org/abs/2412.05085) ·
[hep-ph/0609227](https://arxiv.org/abs/hep-ph/0609227) ·
[hep-ph/0608258](https://arxiv.org/abs/hep-ph/0608258) ·
[arXiv:0901.2560](https://arxiv.org/abs/0901.2560) ·
[JHEP 05 (2019) 025](https://link.springer.com/article/10.1007/JHEP05(2019)025) ·
[arXiv:2002.02282](https://arxiv.org/abs/2002.02282) ·
[Glittering Glasmas](https://www.sciencedirect.com/science/article/abs/pii/S0375947409005077) ·
[arXiv:1209.4105](https://arxiv.org/abs/1209.4105) ·
[arXiv:1808.01276](https://arxiv.org/abs/1808.01276) ·
[PRD 104 (2021) 056025](https://link.aps.org/doi/10.1103/PhysRevD.104.056025) ·
[PRD 95 (2017) 114008](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.95.114008) ·
[arXiv:2110.04881](https://arxiv.org/abs/2110.04881) ·
[arXiv:2302.01380](https://arxiv.org/abs/2302.01380) ·
[arXiv:2406.04856](https://arxiv.org/abs/2406.04856) ·
[arXiv:1503.07126](https://arxiv.org/abs/1503.07126) ·
[arXiv:1804.02910](https://arxiv.org/abs/1804.02910) ·
[arXiv:1612.07790](https://arxiv.org/abs/1612.07790) ·
[arXiv:1808.04982](https://arxiv.org/abs/1808.04982) ·
[arXiv:2605.03038](https://arxiv.org/abs/2605.03038) ·
[hep-ph/9604371](https://arxiv.org/abs/hep-ph/9604371) ·
[hep-ph/9902402](https://arxiv.org/abs/hep-ph/9902402)
