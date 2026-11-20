# Changes in the fixed files (P1--P6)

Files in `nlo_review/fixed_files/`:

| file | base (unchanged parts) | pages |
|---|---|---|
| `Final_terms_fixed.tex/.pdf` | `final_checked/Final_terms_checked.tex` (= author's latest `paper/inputs/project/Final_terms.tex` + the earlier checked corrections, verified by diff: the checked copy differs from the author's file only by red/orange marks) | 38 |
| `Limiting_cases_final_fixed.tex/.pdf` | `limits_checked/Limiting_cases_final_1_checked.tex` (= author's latest `Limiting_cases_final.tex` + earlier checked corrections, verified by diff) | 28 |
| `Divergent_terms_fixed.tex/.pdf` | `limits_checked/Divergent_terms_checked.tex` | 8 |
| `build/` | `assemble.py` (builds the three files from the base files + `build/src/*.tex`; every insertion point is a line number of the base file, every textual replacement is asserted), `phases.py` (phase inventory), `chk/kij.py` (combined kernel of Group I Row III against the sum of the 20 author terms with `log|.|`) | |

Fix notes: `paper/notes/fix_soft.tex` (P1, P2), `fix_poles_dglap.tex` (P3, P4, P5), `fix_singlet.tex` (P6).
`notes/fix_tail.tex` (Group II Row IV tail) did not exist when the files were built; every place that depends on it carries a red NOTE "pending fix_tail".

**Marking.** Red = changed or added (each change starts with a red bracket `[Fix Pn (fix note, Step k) ...]`). Removed: orange in `Final_terms_fixed` and `Limiting_cases_final_fixed` (the replaced equations are reproduced unnumbered, in orange, directly before their replacement, under a heading "[Removed (...)]"); **teal** in `Divergent_terms_fixed`, because there orange is the author's own highlighting. Green (`\cancel`) = cancellations. Author's own marks (blue in Final_terms, dark red in Limiting_cases, orange in Divergent_terms) unchanged. The earlier red/orange marks of the checked copies are kept, except inside removed equations, where they are shown in orange with the rest.

Abbreviations used below: `N_f = (2pi)^-3 g^4/(16 pi^5) 1/k+`; `E(s) = int_0^1 dt/t (e^{-ist}-1)`; `Psi(s) = int_0^inf dt/(t(1+t)) (e^{-ist}-1)`; `Ehat(s) = int_0^inf dt/(1+t) e^{-ist}`; `J(s) = int_0^1 dt e^{-ist}`; `K(a,b,z) = K(a-z).K(b-z)`; `s = k.(z-y')` (Group III), `s = k.(x-z)` (Group I Row V); `O(y,y')` = LO structure; `N_0`, `C = N_0 alpha_s/2pi`; `B(r) = 2 gamma_E + log(r^2 mu^2/4) = log(r^2 mu^2/c_0^2)`, `c_0 = 2 e^{-gamma_E}`; `L = log(k+/Lambda)`, `ell = log(1-k+/vee)`.

---

## 1. Final_terms_fixed.tex

Base line numbers refer to `final_checked/Final_terms_checked.tex`; equation numbers to the compiled `Final_terms_fixed.pdf`.

| # | location | old | new | reason / fix-note step |
|---|---|---|---|---|
| F0 | preamble | -- | `\usepackage[normalem]{ulem}` | needed for the struck-out removed note (F5) |
| F1 | after the title note (l.35) | -- | red header: what the file is, colour code, list of the six changes, pointer to the NOTE list | task |
| F2 | Sec. 3 DGLAP, after l.229 | author's DGLAP term (kept unchanged) with `Ptilde_gg`, `[(1-xi)/xi]_+`, `log((y'-y)^2 mu^2/(4 xi^2))`, no explicit subtraction | added: **Eq. (3.2)** explicit form `N_DGLAP = -C int_0^1 dxi { P_gg(xi) F(xi) - 2N_c (1-xi)/xi F(0) }`, `F(xi) = xi^-2 int e^{-ik.(y-y')/xi} log((y'-y)^2 mu^2/(c_0^2 xi^2)) O(y,y')`, `P_gg` with `xi/(1-xi)_+` and `beta_0/2 delta(1-xi)`, `F(0) = -(4 pi/k^2) X`, `X = int_{x,x',z} O(z,z)`; calculation (`log(r^2mu^2/(c_0^2xi^2)) = 2gamma_E + log(r^2mu^2/(4xi^2))`; `F(0)` from `y'-y = xi s`, `int d^2s e^{-ik.s} log s^2 = -4pi/k^2`); NOTE on the subtraction term (see NOTE 1); scheme statement (a)-(d): bare term in full CDR has no `log xi^2` (reviewer right), counterterm with `xi^{-d}` (number-density MS-bar FF) restores it = author's form, CXY form with `xi^{-2}`; **Eq. (3.3)** `Delta = N_DGLAP - N_DGLAP^CXY = +C int dxi P_gg log xi^2 Pi(xi)`, properties, momentum-sum-rule argument (`int xi P_gg log xi^2 = -14.3225`); `mu` dependence; pole statement (collinear pole times real `P_gg+`); **Eq. (3.4)** local form `F(xi) = -log(k^2/mu^2) Pi(xi) + xi^-2 R(k/xi)` and the resulting local DGLAP form; NOTE pending fix_tail (coordinate-space IR tail at large parent position) | P4: fix_poles_dglap Steps 8, 9, 10, 11, 12; P3: Step 4; tail: Step 13 |
| F3 | Sec. 4 RC, after l.234 | author's RC term (kept) | added: `N_RC = +C beta_0/2 F(1)` exactly; `N_DGLAP + N_RC` with the two `beta_0/2 F(1)` terms cancelled in green; what the RC term is (pole collinear; coupling renormalised by UV half of `Z_3` of the produced gluon, `Z_g = Z_3^{-1/2}`; IR half = virtual `beta_0/2 delta(1-xi)` of `P_gg`; `d N_RC / d log mu^2 = +alpha_s beta_0/4pi LO`; KLSZ classification; the 11/3 virtual terms cancel between Group I Rows II and V) | P3/P4: fix_poles_dglap Steps 4, 7, 11, 12 |
| F4 | Sec. 4, after l.239 (LO+RC equation, kept) | -- | local form `LO + RC = N_0 L(k) [1 - alpha_s beta_0/4pi log(k^2/mu^2)] + C beta_0/2 R(k)` (= LO x alpha_s(k^2)/alpha_s(mu^2)); NOTE: choice coordinate/local form | P4: Step 12 |
| F5 | Sec. 5 intro, l.241 | red note of the earlier check: "Group III Rows III, V, VI keep `vee-k+` as the upper limit of the `p+` integral ..." | that note removed (orange, struck out); new red note: upper limit now `infinity`, instantaneous terms as distributions; no `1/eps` written, all poles accounted for in `Divergent_terms_fixed`; remaining `-(ell/2)(2/eps)` requires `vee -> infinity` before `eps -> 0` | P2: fix_soft (C1), Step 9; P3: fix_poles_dglap Step 6 |
| F6 | before Sec. 6 Group I (after l.242) | -- | new red subsection "Notation and conventions of the fixed finite terms": phase, c.c., operator order, PV/`log|.|`, limits; `N_f`, `K`, `K(a,b,z)`, `P`, `L^a`, `Delta`; recoil functions `E, Psi, Ehat, J, E_inf` with asymptotics, **Eq. (5.1)**; the exact split of a plus-prescription integral into JIMWLK-list term + soft remainder + regular term, **Eq. (5.2)** | fix_soft "Notation", Step 1 |
| F7 | Group I Row II, after l.250 | -- (row unchanged) | red note: the virtual 11/3 term `(alpha_s beta_0/4pi) N_0 int [4/eps + B(x-w) + B(x'-w')] O(w,w')` is present in this row (shown, green-cancelled), cancelled exactly (pole, `2gamma_E`, log, same `mu`) by the `-11/6` part of `N_2b UV` of Row V after `x->w, y->x`; UV-UV; hence no finite `beta_0 [B+B'-B(w-w')]` term | P3: fix_poles_dglap Steps 1, 2, 7 |
| F8 | Group I Row III, l.264-313 | the four equations `N_1, N_2, N_3a, N_3b` (6.4)-(6.7), incl. the recoil form of the last term of `N_3b` from the earlier check, and its notes | **removed** (orange); replaced by **Eq. (6.4)**: `(2pi)^-3 i g^4 f^{abd}/(8pi^5) 1/k+ int e^{-ik.(w-w')} { K^i(x'-w') K^j(y-z) Kcal^{ij} + [-1/2 log((z-w)^2/(x-w)^2) - E(k.(z-w))] K(y,w,z) K(x'-w').K(x-w) } L^a [U^de(x)U^bc(z) - U^de(x)U^bc(y)] rho^e(x) rho^c(y)`, `w` integral first; **Eq. (6.5)** combined kernel `Kcal^{ij}` (only positive log arguments, pole at `(x-z)^2=(x-w)^2` cancelled); calculation note: (i) PV and the explicit `log|.|` forms `log|(x-z)^2/((x-z)^2-(x-w)^2)|` (N_1, N_2), `log|((x-w)^2-(x-z)^2)/(x-w)^2|` (N_1, N_3a); (ii) sum of the 20 terms = `Kcal^{ij}` term (re-checked numerically, rel. diff <= 5e-14, `build/chk/kij.py`); (iii) third term of `N_3b` = `2i N_f l(z,w,x) K(y,w,z) P.K(x-w) C_D`, split `l = E + (l - E)`; (iv) last term of `N_3b` = `-2i N_f E K(x,w,z) P.K(y-w) C_A` moved to block S (pair A); integrability note incl. Fourier identity with the two cancelling `chi(theta)` terms (green) and condition (C3); answer to "the problematic term should cancel with something": with Group III Row VI (pair D, log part) and with the `A_0` tail of the other terms of the row after the `w` integration | P1, P5: fix_soft Steps 5, 6, 8 |
| F9 | Group I Row V, l.318-337 | first and second part of Row V with `[k+/p+]_+`, `[1/q+]_+`, `[k+/q+]_+`, `[k+/(k+-p+)]_+` | **removed** (orange); replaced by **Eq. (6.7)** `N^{I,V}_1 = (2pi)^-3 i g^4/(32pi^5) 1/k+ int e^{-ik.(z-w')} J(s) K(x'-w').K(z-x) K(y-z).K(y'-x) L^a[...]` and **Eq. (6.8)** `N^{I,V}_2` with `F_reg^{im} = int_0^1 dxi e^{-i xi s} N^{im}(xi)/(xi(y-x)^2+(1-xi)(y-z)^2)` and the `N_c 11/6 e^{-is}` term; calculation note (plus prescriptions on `f(xi)=e^{-i xi s}/D(xi)`, `int_0^1 dxi [xi(1-xi)-2] = -11/6`, prefactors); remainders `E(s)`, `e^{-is}E(-s)` moved to S (pairs B, C) | P1, P2: fix_soft Steps 1, 2, 6 |
| F10 | Group II Row I, after l.374 | -- (row unchanged) | NOTE pending fix_tail: compensation of the Group II Row IV tail by this row (second part) and Group III Row III not checked term by term | fix_poles_dglap Step 14 |
| F11 | Group II Row IV, after l.383 | -- (row unchanged) | note: (i) the second term (remnant `-(2pi)^-3 g^4 N_c/(2pi^3) 1/k^2 1/k+ X = -2N_c C (4pi/k^2) X`) goes with the subtraction `2N_c (1-xi)/xi F(0)` and must be dropped if DGLAP is written with `2N_c/xi F(0)`; (ii) NOTE pending fix_tail with the tail coefficients `+N_c/3 cos 2phi, -2N_c, +N_c/6, -11N_c/6` | P4: Step 11; tail: Step 13 |
| F12 | Group III Row III, l.391-404 | first part (`int_0^{vee-k+} dp+/2pi`, `[delta delta/(p++k+) - ([k+/(p+(p++k+))]_+ + 1/(p++k+)) delta delta]`, instantaneous `-1/k+` term) and second part | **removed** (orange); replaced by **Eq. (8.2)** `N^{III,III}_1 = (2pi)^-3 i g^4 f^{c'd'a}/(16pi^5) 1/k+ int e^{-ik.(w-y')} { Ehat(s)[K(y'-z).K(x'-y') K(y-w).K(x-z) - K(y',x,z) K(x'-y').K(y-w)] L^{c'd'b}(y') - [pi delta(s) - i P 1/s] K(y'-z).K(y-w) K^j(x-z)[K^j(x'-y') L(y') - K^j(x'-z) L(z)] } R_1^{ab}` and **Eq. (8.3)** `N^{III,III}_2` (regular brace with coefficient `(B-A)/((p++k+)A)` of `Jhat`, instantaneous `p+/D_p` term with `y'` before `p+`, closed-form last term `A/(B(A-B)) log(A/B)`); **Eq. (8.4)** `int_0^inf dp+ e^{-i p+ s/k+} = k+[pi delta(s) - i P 1/s]` (c.c. `+i`); calculation (prefactors `ig^4f/(8pi^4k+) . 1/2pi = ig^4f/(16pi^5k+)`, `x(-1/k+) . k+` for the distributional term, plus prescription with phase in the subtraction, partial fractions), finiteness of the distributional terms, condition (C7), integrability; remainders `-i N_f Psi K(y',x,z) ...` and `-2 N_f Psi P.K(x-w)[1/2 K(x,y',z) - K(w,y',z)] ...` moved to S (pairs B, C); NOTE pending fix_tail | P2: fix_soft Steps 1, 2, 6, 9; review Q4 |
| F13 | Group III Row V, l.408-410 | Row V with `int_0^{vee-k+}` and instantaneous term | **removed**; replaced by **Eq. (8.5)** (prefactor `ig^4 f^{c'ca}/(8pi^5) 1/k+`, `Ehat(s)` term, distributional term); remainder `-2i N_f Psi K(y',x,z) K(x'-y').K(y-w) f L_V R_V` moved to S (pair A) | P2: fix_soft Steps 2, 6, 9 |
| F14 | Group III Row VI, l.412-414 | Row VI, same form | **removed**; replaced by **Eq. (8.6)** (Row V with `f^{c'ca} L_V R_V -> f^{c'd'a} L^{c'd'b} R_VI`, written out); remainder (pair D) moved to S | P2: fix_soft Steps 2, 6, 9 |
| F15 | new Sec. 9 | -- | **soft-recoil block S**: colour structures `C_A, C_B, C_C, C_D` (9.1); the eight remainders, as they come from the rows and in common variables (9.2); colour identities `f L_V R_V = -C_A`, `f L R_1 = -C_B`, `f L R_2 = -C_C`, `f L R_VI|_{x<->y} = +C_D` (9.3); pair form `alpha_X V_X C_X [kappa_X(w) E(k.(z-w)) - kappa_X(w') Psi(k.(z-w'))]` (9.4) with table `alpha, V, kappa, T`; tail cancellation (green), real-minus-virtual decomposition (9.5); the block **Eq. (9.6)** `N^soft = N_f int e^{-ik.(w-w')} { -i P.K(y-w)[2C_A + C_B] F[K(x,.,z)] - P.K(x-w) C_C F[K(x,.,z) - 2K(w,.,z)] + 2i P.K(x-w) C_D F[K(y,.,z)] }`; equivalent kinematically-constrained reading (not used); conditions (C1)-(C9); not done / NOTE (`+-i eps` prescription) | P1, P2: fix_soft Steps 2, 3, 4, 6, 7, conditions |
| F16 | new Sec. 10 | -- | phase convention: inventory of all 82 phases of the file (script), all `e^{-ik.(unprimed-primed)}`; six written `e^{-ik(w-w')}` without the dot (same convention, unchanged); effect of a mixed convention (odd parts flip sign) | P5: fix_poles_dglap "Phase convention" |
| F17 | new Sec. 11 | -- | colour-singlet averaged form of all finite terms, DGLAP, RC and LO (Eqs. 11.1-11.23): assumptions (A1) colour invariance, (A2) no symmetrised three-point function; quadratic average; cubic terms = symmetrised part (zero under A2, odd in k) **plus contact terms** `i Phi delta(p_i-p_j) n /(N_c^2-1)` which do not vanish; singlet factors `D, F, Q` for every structure; every row of Groups I, II, III; DGLAP/RC/remnant/LO singlet forms; which rows are cubic; numerical checks; NOTE list of fix_singlet. For the soft-remainder rows the singlet form is written with the full row kernel (regular term + the row's own S member) | P6: fix_singlet (all sections) |
| F18 | new Sec. 12 | -- | NOTE list (11 items, see below) | task |

## 2. Limiting_cases_final_fixed.tex

Base line numbers refer to `limits_checked/Limiting_cases_final_1_checked.tex`. The finite terms are the same as in Final_terms_fixed; DGLAP/RC are not part of this file.

| # | location | old | new | reason / fix-note step |
|---|---|---|---|---|
| L0 | preamble | -- | `ulem` | as F0 |
| L1 | after l.36 (note of the earlier check) | -- | red header; states that the earlier note's statement on the upper limit `vee-k+` of Group III is superseded | task |
| L2 | after l.37 | -- | same red note as F5 (upper limit infinity, distributions, poles, order of limits) | fix_soft (C1), Step 9; fix_poles_dglap Step 6 |
| L3 | after L2 | -- | notation subsection as F6 (equations tagged N.1, N.2) | fix_soft Notation, Step 1 |
| L4 | after L3 | -- | red "DGLAP, RC and LO notation used in the notes": `O, L(q), N_0, C, c_0`, DGLAP explicit form (tag N.3), RC, `[(1-xi)/xi]_+` convention | so that the notes and the singlet section are self-contained; fix_poles_dglap Step 11 |
| L5 | Group I Row II, after l.45 | -- | as F7 | P3: Steps 1, 2, 7 |
| L6 | Group I Row III, l.59-108 | `N_1, N_2, N_3a, N_3b` and notes | removed (orange); replaced as F8 (Eqs. 1.4, 1.5) | P1, P5: fix_soft Steps 5, 6, 8 |
| L7 | Group I Row IV, after l.118 | -- | red note: Group I Row IV and Group III Row I (dark red, author's marking) cancel identically, also in singlet form; no soft remainder | fix_singlet "Rows not listed" |
| L8 | Group I Row V, l.120-139 | first and second part | removed; replaced as F9 (Eqs. 1.8, 1.9) | P1, P2: fix_soft Steps 1, 2, 6 |
| L9 | Group II Row I, after l.176 | -- | NOTE pending fix_tail (as F10) | fix_poles_dglap Step 14 |
| L10 | Group II Row IV, after l.185 | -- | remnant / subtraction note + NOTE pending fix_tail (as F11) | Steps 11, 13 |
| L11 | Group III Row III, l.198-211 | first and second part | removed; replaced as F12 (Eqs. 3.3-3.5) | P2: fix_soft Steps 1, 2, 6, 9 |
| L12 | Group III Row V, l.215-217 | Row V | removed; replaced as F13 (Eq. 3.6) | P2 |
| L13 | Group III Row VI, l.219-221 | Row VI | removed; replaced as F14 (Eq. 3.7) | P2 |
| L14 | new Sec. 4 | -- | soft-recoil block S (as F15) | P1, P2 |
| L15 | new Sec. 5 | -- | phase convention (66 phases of this file, all uniform) | P5 |
| L16 | new Sec. 6 | -- | colour-singlet section (as F17) | P6 |
| L17 | new Sec. 7 | -- | NOTE list (as F18, references to the DGLAP term point to Final_terms_fixed) | task |

## 3. Divergent_terms_fixed.tex

Base line numbers refer to `limits_checked/Divergent_terms_checked.tex`. Removed/superseded statements in **teal**.

| # | location | old | new | reason / fix-note step |
|---|---|---|---|---|
| D0 | preamble | -- | `ulem` | |
| D1 | after the title note (l.35) | -- | red header: colour code, list of changes, notation (`K, c_0, B(r), beta_0, a = alpha_s N_c/pi, O, L(q), N_0, C, Pi(xi)`), phase statement (all phases uniform) | P3, P5 |
| D2 | Group I Row II, after l.42 | -- | origin of the pole: one-rho part of `A^(3)` (I.89), `(11/3 - 2L + ell)[-2/eps + log(k^2/mu^2)]`, Fourier `int dk J_1(kr) log(k^2/mu^2) = -B(r)/r`; **Eq. (1.2)** `N^II_{11/3} = +(alpha_s beta_0/4pi) N_0 int e^{-ik.(w-w')} [4/eps + B(x-w) + B(x'-w')] O(w,w')` (c.c. incl.), = the reviewer's "missing" virtual term; region `p~ -> inf`: **UV** | P3: fix_poles_dglap Step 1 |
| D3 | Group I Row V, after l.53 | -- | after `x->w, y->x`: **Eq. (1.4)** `N_2b UV|_{-11/6} = -N^II_{11/3}`; region `z -> x`, no large-z divergence: **UV**; sign convention of the UV pole (scaleless minus IR-divergent integral) | P3: Step 2 |
| D4 | "Adding Group I", after l.60 | -- | the two equations written as one term with `{ -11/3 + 11/3 - ell }` (the 11/3's green-cancelled); UV-UV cancellation, physical reason (`Omega^dag S Omega = S` for a locally constant shock wave), consequence: no finite `beta_0[B(x-w)+B(x'-w')-B(w-w')]` LO term, emission coupling does not run with `x-w`; same `mu` required | P3: Steps 1, 2, 7 |
| D5 | after l.68 (the "= 0" equation and its bookkeeping note, kept) | -- | red note: the arithmetic of the bookkeeping is right, but (a) the DGLAP collinear pole is absorbed by the MS-bar FF counterterm with the full `P_gg`, (b) the RC pole is collinear; the coupling is renormalised by the UV half of the dressed-gluon normalisation; total with counterterms `-ell/2`, needs `vee -> inf` before `eps -> 0` | P3: Steps 3-6 |
| D6 | DGLAP note (iv), l.86 | "the pole of the DGLAP term ... is not cancelled by any other term." | that clause in teal + red "Superseded": collinear (IR), multiplies real `P_gg+`, absorbed together with the virtual `beta_0/2 delta(1-xi)` of the normalisation by the FF counterterm with full `P_gg` | P3, P4: Steps 4-6, 9 |
| D7 | after l.86 | -- | the DGLAP block of F2 (explicit form with xi->0 subtraction, NOTE on the subtraction term, scheme statement, `Delta`, poles, local form, NOTE pending fix_tail), Eqs. (3.2)-(3.4) | P4: Steps 8-12 |
| D8 | after l.90 (RC and the `log(1-k+/vee)` term) | -- | the RC block of F3 | P3, P4: Steps 4, 7, 11, 12 |
| D9 | RC section note, l.112 | "g^2 -> g^2(mu) presupposes that the 2/eps of the checked N_RC is removed by the coupling renormalisation; ..." (kept) | appended red: fulfilled, but not by the RC pole: coupling CT matched by the UV half of the dressed-gluon normalisation (`Z_3`); RC pole collinear, cancels the IR half; only in single-eps bookkeeping does the coupling CT cancel the RC pole formally | P3: Steps 3, 5, 7 |
| D10 | "Decisions recorded", l.112 | "no 1/eps poles are written (rule int d^2p/p^2 = 0), so the order of the limits vee->inf and eps->0 does not arise;" | clause "so the order of the limits ... does not arise" in teal + red "Superseded": the order does arise; total `-(ell/2)(2/eps)` zero only for `vee -> inf` first; must be stated | P3: Step 6 |
| D11 | after l.112 | -- | local form of LO+RC (as F4) | P4: Step 12 |
| D12 | new Sec. 6 | -- | **pole table and account**: regularisation (`int d^dp/p^2 = pi(2/eps_UV - 2/eps_IR)`); Step 3 normalisation of the dressed gluon, **Eq. (6.1)** `N_norm = [alpha_s beta_0/4pi - a L](2/eps_UV - 2/eps_IR) LO` with calculation; Step 4 Group II Row IV poles, **Eq. (6.2)**; Step 5 coupling and FF counterterms; Step 6 table (origin, UV/IR, coefficient, cancelled by) for Row II, Row V, normalisation UV/IR halves, coupling CT, Row IV `xi->1`, real `P_gg+`, `xi->0`, FF CT; sum UV `ell/2`, IR `-ell`; single-eps bookkeeping with green cancellations `= -ell/2`; Step 7 answers to Q2 of the review (items 1-4) | P3: Steps 3-7 |
| D13 | new Sec. 7 | -- | NOTE list (5 items) | task |

---

## 4. Calculations shown explicitly in the files (where to find them)

- Group I Row III combined kernel, PV and `log|.|`, the third and last term of `N_3b`, the Fourier identity of the tails: Final Sec. 6.3 (calculation, integrability, answer notes).
- Group I Row V plus prescriptions and the `-11/6`: Final Sec. 6 Row V calculation note.
- Group III Row III `p+` integrals, plus prescription with the phase, partial fractions, distributional identity and its prefactor: Final Sec. 8, Row III notes.
- Soft block: colour identities, pair table, tail cancellation: Final Sec. 9.
- DGLAP explicit form, `F(0)`, scheme, `Delta`, local form; RC = `-delta` part of DGLAP: Final Secs. 3, 4.
- Pole algebra (Row II/Row V relabelling, prefactor identities `g^4N_c/(32pi^4) 11/3 = g^2/(2pi^2) alpha_s beta_0/4pi`, normalisation, table, single-eps sum): Divergent_terms_fixed Secs. 1, 2, 6.
- Singlet factors and contact-term coefficients: Final Sec. 11.

## 5. NOTE questions for the author (all files)

1. **DGLAP subtraction term.** fix_poles_dglap Step 11 writes the xi->0 subtraction as `-(2N_c/xi) F(0)`; this equals the author's `-2N_c (1-xi)/xi F(0)` **plus** exactly the remnant `-(2pi)^-3 g^4 N_c/(2pi^3) (1/k^2)(1/k+) X` of Group II Row IV (checked: `int_0^1 dxi [2N_c/xi - 2N_c(1-xi)/xi] F(0) = 2N_c F(0)`, `2N_c C 4pi/k^2 = g^4 N_c/(2pi^3)/(2pi)^3/k+`), while Step 11 also keeps the remnant in Row IV (would count it twice). Author's convention kept (remnant in Row IV). Which writing for the paper?
2. **FF scheme.** The author's `log((y'-y)^2 mu^2/(4 xi^2))` is kept; it corresponds to MS-bar number-density FFs renormalised with the d-dimensional Jacobian `xi^{-d}`; CXY form differs by `Delta`. Confirm the FF definition; no independent comparison with a standard NLO collinear `pp -> gX` calculation was made.
3. **Coordinate or local form of DGLAP and RC.** The coordinate forms are each logarithmically IR divergent at large parent position; the local form is not. Choice for the paper (pending fix_tail).
4. **Group II Row IV large-distance tail** (pending fix_tail): compensation by Group II Row I (second part) and Group III Row III established only at the amplitude level (fix_poles_dglap Step 14); term-by-term check on the finite terms as written not done.
5. **Prescription of the energy denominators of B_3** (Group I Row III): principal value assumed; a `+-i eps` prescription would add `-+i pi theta((x-z)^2-(x-w)^2)`-type on-shell terms; not evaluated. Confirm PV.
6. **Order of limits.** The remaining pole `-(ell/2)(2/eps)` vanishes only if `vee -> inf` is taken before `eps -> 0`; should be stated in the paper.
7. **Cubic contact terms** (singlet section): the cubic rows leave quadratic contact terms after the projectile average (they do not vanish, contrary to the brief/reviews); whether they cancel among the rows or against the contact terms of the three-rho evolution terms was not checked.
8. **Integrability at the contact configurations** (singlet section): checked only for two kernels.
9. **Assumption (A2)**: no C-odd (`d^{abc}`) or chiral (`f^{abc} g_f`) three-point function of the projectile; otherwise the symmetrised cubic parts (odd in k) are added.
10. **Kinematically constrained reading** of the soft block (fix_soft Step 7): equivalent, not used; the author may prefer it.
11. **Not done in fix_soft**: complete numerical integration with Wilson lines; re-check of the regular terms (taken from Finite_terms_final); the non-singlet `L int d^2z/z^2 C_A` of the JIMWLK formula.
12. (Divergent_terms only) **Two-eps vs single-eps bookkeeping**: both give the same finite result; the paper should say which is used, and with single-eps that the cancellation of the RC pole by the coupling counterterm rests on the scaleless normalisation integral.

## Update after fix_tail (notes/fix_tail.tex)

The Group II Row IV large-distance tail was checked term by term: in the joint limit (w, w', z → ∞ at fixed w − z) Group II Row IV and Group II Row I (second part) each carry −11/6 + (1/6)cos2φ, Group III Row III (second part, with c.c.) carries +11/3 − (1/3)cos2φ (units 2N_c C (4π/k²)∫d²Y/Y² G_∞); the sum is zero, also in the BFKL and JIMWLK sectors. No written finite term needs correcting. All 15 "pending fix_tail" marks in the three files (and in build/src) were replaced by "settled in fix_tail" with this result; the three rows must be integrated together over the pair position. For a coloured projectile and a target not uniform at large distance a remainder ∝ total charge survives (the known non-singlet JIMWLK term).

## Update: author's answers of 2026-11-20 (NOTE questions 2, 3, 5, 10)

Verbatim answers: (FF scheme) "Yes"; (DGLAP/RC coordinate vs local form) "write it in both forms"; (principal values for the energy denominators of B_3) "yes"; (soft block, paired vs kinematically constrained form) "write it in both form". Every change is in red and carries "author's answer (date 2026-11-20)". Edited fragments in `build/src` (synced to the scratch build, re-assembled with `assemble.py`, compiled twice, regenerated `.tex` diffed against the previous ones: only the changes below moved). The order of limits (NOTE 6) and the three-rho / cubic / (A2) material (NOTEs 7-9, `singlet.tex`) are not touched.

| fragment | files | change |
|---|---|---|
| `dglap_add.tex` | Final, Divergent | Scheme note: "Scheme statement" -> **Definition** (MS-bar number-density FFs, pure pole, Jacobian xi^-d); Delta to CXY marked as a remark; the NOTE on the choice of definition replaced by the answer (comparison with a standard NLO `pp -> gX` calculation listed as a cross-check only). "Equivalent local form" -> **Both forms**: local DGLAP and RC terms `N^loc_DGLAP = C log(k^2/mu^2) int P_gg Pi`, `N^loc_RC = -(alpha_s beta_0/4pi) log(k^2/mu^2) LO`, split of the coordinate terms (new Eq. `eq:fx-locsplit`), equivalence `N_DGLAP + N_RC = N^loc_DGLAP + N^loc_RC + N_R` with `N_R = -C int {P_gg+ R(k/xi)/xi^2 - 2N_c (1-xi)/xi F(0)}` (new Eq. `eq:fx-equiv`; in the author's subtraction convention; fix_poles_dglap's N_R = this + remnant). The stale last note ("is being checked ... not settled here") replaced by **Large-distance behaviour of the two forms**: coordinate forms separately log-divergent at large parent position, sum with Group II Row I (2nd part) and Group III Row III not (fix_tail numbers), local form: N_R carries the DGLAP + RC tails (+N_c/6 - 11N_c/6). `\sloppy` added for that paragraph. |
| `lorc_add.tex` | Final, Divergent | NOTE (choice of form) replaced by the answer; local RC term and reference to the equivalence. |
| `soft_block.tex` | Final, Limiting | "Equivalent reading (not used here)" replaced by a subsection **Kinematically constrained form of the block**: `F_Lambda[kappa]` (Eq. `eq:fx-FLam`, with the short calculation), the constrained block `N^soft_kc` written out in full (Eq. `eq:fx-Skc`), equivalence `N^soft_kc = N^soft + log(k+/Lambda) N^J_S` (Eq. `eq:fx-kcequiv`), `N^JIMWLK + N^soft = [N^JIMWLK - L N^J_S] + N^soft_kc`. New **NOTE (pair D)**: the Group I member of pair D (the E part of the third term of N_3b) has no log(k+/Lambda) term of its own; that `L alpha_D V_D C_D [kappa_D(w) - kappa_D(w')]` is part of the JIMWLK list is taken from fix_soft Step 7, not re-derived (the identity itself does not depend on it). (C6) and "Not done (4)": PV stated as the prescription, NOTE removed. |
| `g1r3_new.tex` | Final, Limiting | PV calculation (i): author's answer added. |
| `notes_final.tex`, `notes_lim.tex` | Final, Limiting | Items 2, 3, 5, 10 marked **Answered** with the answer; item 4 (Row IV tail) updated to **Settled** (consistency: it still read "pending" after the fix_tail update). |
| `notes_div.tex` | Divergent | Items 2, 3 marked **Answered**. |
| `g2_notes.tex` | Final, Limiting | Consistency with fix_tail: Group II Row I note ("has not been checked term by term") and Group II Row IV note (ii) ("Pending fix_tail ... not settled here") updated to the settled result; Row IV note (ii) also states where the tail sits in the two forms. |
| `header_final.tex`, `header_lim.tex`, `div_header.tex` | all | one sentence listing the implemented answers. |

**Order of limits (coordinator's additional task, not implemented).** The proposed prescription (epsilon -> 0 first; FF counterterm integrated over the real-emission range k+/vee < xi < 1) was checked with exact algebra (`scratchpad/fix/poles/ct/s7_range_ct.py`, `s7.out`): with the counterterm restricted to the same range (plus distribution of P_gg on [0,1] applied to theta(xi > k+/vee)), the IR poles cancel exactly at finite vee (the counterterm supplies +2N_c Pi(1) ell, cancelling the -ell of the Row IV xi -> 1 end), but the UV pole (ell/2)(2/eps) of Group I Row II (A^(3) coefficient (11/3 - 2L + ell)) is untouched by an IR counterterm: total (ell/2)(2/eps) (two-eps and single-eps alike); with the plus subtraction restricted to [k+/vee, 1] instead, the total stays -(ell/2)(2/eps). The pole therefore does not cancel at finite vee; nothing was changed (NOTE 6 and all order-of-limits statements as before).

## Update: pair D in the kinematically constrained form (re-derived 2026-11-20)

`soft_block.tex` (Final, Limiting), `notes_final.tex`, `notes_lim.tex`; paper `sec_finite.tex`. The NOTE / \pending on pair D is replaced (in red in the fixed files) by the result of a check against `Group_I_checked` and `Group_III_checked`:
- only the real half of pair D, −log(k+/Λ) α_D V_D C_D κ_D(w′), is in the JIMWLK list (from Group III Row VI);
- the virtual half, +log(k+/Λ) α_D V_D C_D κ_D(w), is not. The third term of N_3b of Group I Row III carries −½log((z−w)²/(x−w)²) in place of log(k+/Λ). All log(k+/Λ) terms of Row I.III have the structure (x′−w′)·(y−z) and combine into the pair-A term, and the pair-D kernel appears there only with log(∨/k+).
- For pair D this term is therefore added and subtracted; Eq. (fx-kcequiv) / (kcequiv) holds as written.

The wording "pairs A–D" in the definition of N^J_S is adjusted accordingly. This is a structural argument from the checked formulas, not a numerical evaluation. Order-of-limits and cubic/(A2) text are not touched.
