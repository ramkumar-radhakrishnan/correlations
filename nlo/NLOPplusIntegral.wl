(* ::Package:: *)

(* ::Title:: *)
(*NLOPplusIntegral`*)

(* ::Text:: *)
(* Longitudinal (p^+) integration of the NLO gluon-production integrand    *)
(*                                                                          *)
(*   d^3 N_NLO / d^3 k = -(1/(2 Pi)^3) (I g^4 f^{c'd'a}/(4 Pi^4)) (1/k^+)  *)
(*      Integrate[dp^+/(2 Pi), {p^+, Lam, Vee - k^+}]                       *)
(*      Integrate[..., {w, y', z, y, x, x'}] E^(-I k.(y'-w))                *)
(*      E^(-I k.(p^+/k^+)(y'-z))                                            *)
(*      (y'-z)^m (x'-y')^k' (y-w)^i (x-z)^j                                 *)
(*        / ((y'-z)^2 (x'-y')^2 (y-w)^2 (x-z)^2)                            *)
(*      [ d_{k'm} d_{ij}/(p^+ + k^+) - d_{jm} d_{ik'}/p^+                   *)
(*        - d_{im} d_{jk'}/k^+ ]                                            *)
(*      [color structure 1] [color structure 2]                             *)
(*                                                                          *)
(* Kinematic window: Vee >> k^+ >> Lam > 0.                                 *)
(*                                                                          *)
(* Only the phase E^(-I p^+ omega / k^+), with                              *)
(*                                                                          *)
(*   omega = k.(y'-z)   (transverse dot product),                           *)
(*                                                                          *)
(* and the three p^+ poles depend on p^+, so the p^+ integral factorizes    *)
(* completely out of the transverse integrals and the color operators.      *)
(* Everything below is that factorized p^+ integral, done in closed form,   *)
(* plus the reassembled integrand.                                          *)

BeginPackage["NLOPplusIntegral`"];

NLOAssumptions::usage =
  "NLOAssumptions[kp, lam, vee] gives the ordering assumptions Vee > k^+ > Lam > 0 \
used when deriving the p^+ integral symbolically.";

PPlusKernel::usage =
  "PPlusKernel[omega, kp, lam, vee] gives {I1, I2, I3}, the exact values of\n\
  I1 = Integrate[dp/(2 Pi) E^(-I p omega/kp) / (p + kp), {p, lam, vee - kp}]\n\
  I2 = Integrate[dp/(2 Pi) E^(-I p omega/kp) / p,        {p, lam, vee - kp}]\n\
  I3 = Integrate[dp/(2 Pi) E^(-I p omega/kp) / kp,       {p, lam, vee - kp}]\n\
with omega = k.(y'-z). These multiply d_{k'm} d_{ij}, -d_{jm} d_{ik'} and \
-d_{im} d_{jk'} respectively. The omega -> 0 limit is handled separately.";

PPlusKernelDerive::usage =
  "PPlusKernelDerive[omega, kp, lam, vee] recomputes {I1, I2, I3} with Integrate \
under NLOAssumptions instead of quoting the closed form. Use it to re-derive / \
audit PPlusKernel.";

PPlusKernelAsymptotic::usage =
  "PPlusKernelAsymptotic[omega, kp, lam, vee] gives {I1, I2, I3} in the strongly \
ordered regime Vee >> k^+ >> Lam with |omega| Lam / k^+ << 1 and \
|omega| Vee / k^+ >> 1: Lam-dependence survives only as the rapidity logarithm in \
I2, and the Vee-dependence survives only in the non-decaying boundary oscillation \
of I3.";

PPlusKernelCollinear::usage =
  "PPlusKernelCollinear[kp, lam, vee] gives {I1, I2, I3} at omega -> 0 (the phase \
switched off). Here the Vee logarithm reappears in I1 and I2 and I3 grows \
linearly with Vee -- i.e. the phase, not the cutoff, is what tames the last term.";

RapidityLog::usage =
  "RapidityLog[omega, kp, lam, vee] gives the coefficient of the single logarithm \
carried by the 1/p^+ term, i.e. the log that multiplies -d_{jm} d_{ik'}.";

TransverseData::usage =
  "TransverseData[k, w, yp, z, y, x, xp] returns an association with the transverse \
building blocks: the six index contractions produced by the delta structures, the \
four propagator denominators, the phase E^(-I k.(y'-w)) and omega = k.(y'-z). \
Arguments may be explicit 2-component vectors (numeric or symbolic) or opaque \
symbols, in which case dot products stay inert as TDot[u, v].";

TDot::usage = "TDot[u, v] is the inert transverse (2d) dot product u.v.";

TransverseFactor::usage =
  "TransverseFactor[data, {I1, I2, I3}] contracts the p^+ kernel {I1, I2, I3} with \
the transverse numerator and divides by the four denominators. data comes from \
TransverseData.";

ColorStructure::usage =
  "ColorStructure[idx] gives f^{c'd'a} times the product of the two color brackets \
as an inert non-commutative expression built from U[a, b][pt] (adjoint Wilson line) \
and Rho[a][pt] (color charge). Repeated adjoint indices are summed.";

U::usage = "U[a, b][pt] is the adjoint Wilson line U^{ab}(pt) (inert).";
Rho::usage = "Rho[a][pt] is the color charge density rho^a(pt) (inert).";
FStruct::usage = "FStruct[a, b, c] is the structure constant f^{abc} (inert).";

NLOIntegrand::usage =
  "NLOIntegrand[k, w, yp, z, y, x, xp, kp, lam, vee, g, idx] gives the full \
d^3N_NLO/d^3k integrand after the p^+ integration has been performed, with the \
remaining transverse integrals left inert. Option \"Kernel\" -> \"Exact\" \
(default), \"Asymptotic\" or \"Collinear\" selects the p^+ kernel.";

NLOResult::usage =
  "NLOResult[...] is NLOIntegrand wrapped in the inert transverse integrals \
Integrate[..., {w, y', z, y, x, x'}].";

CheckPPlusKernel::usage =
  "CheckPPlusKernel[] numerically compares PPlusKernel, PPlusKernelDerive and \
PPlusKernelAsymptotic against direct NIntegrate over p^+ for a spread of \
parameters. Returns an association of maximum absolute deviations.";

Begin["`Private`"];

(* ------------------------------------------------------------------ *)
(* 0. Assumptions                                                       *)
(* ------------------------------------------------------------------ *)

NLOAssumptions[kp_, lam_, vee_] := {
  kp > 0, lam > 0, vee > kp + lam,
  Element[{kp, lam, vee}, Reals]
};

(* ------------------------------------------------------------------ *)
(* 1. The p^+ integral                                                  *)
(* ------------------------------------------------------------------ *)

(* Exact:                                                               *)
(*   Integrate[E^(-I b p)/(p + kp), {p, lam, vee - kp}]                  *)
(*     = E^(I b kp) (Ei[-I b vee] - Ei[-I b (kp + lam)])                 *)
(*   Integrate[E^(-I b p)/p, {p, lam, vee - kp}]                         *)
(*     = Ei[-I b (vee - kp)] - Ei[-I b lam]                              *)
(*   Integrate[E^(-I b p), {p, lam, vee - kp}]                           *)
(*     = (E^(-I b lam) - E^(-I b (vee - kp)))/(I b)                      *)
(* with b = omega/kp. Verified against Integrate and NIntegrate.         *)

PPlusKernel[omega_, kp_, lam_, vee_] :=
  If[TrueQ[PossibleZeroQ[omega]],
    PPlusKernelCollinear[kp, lam, vee],
    With[{b = omega/kp},
      (1/(2 Pi)) {
        Exp[I b kp] (ExpIntegralEi[-I b vee] - ExpIntegralEi[-I b (kp + lam)]),
        ExpIntegralEi[-I b (vee - kp)] - ExpIntegralEi[-I b lam],
        (Exp[-I b lam] - Exp[-I b (vee - kp)])/(I b kp)
      }
    ]
  ];

PPlusKernelDerive[omega_, kp_, lam_, vee_] :=
  Module[{b = omega/kp, p},
    (1/(2 Pi)) Integrate[
      Exp[-I b p] {1/(p + kp), 1/p, 1/kp},
      {p, lam, vee - kp},
      Assumptions -> Append[NLOAssumptions[kp, lam, vee], b > 0]
    ]
  ];

(* Strongly ordered regime.  Two facts drive the limits:                *)
(*   Ei[-I x] -> Ei[-I x] (small x)  = EulerGamma + Log[x] - I Pi/2,     *)
(*   Ei[-I X] -> -I Pi Sign[X]       (|X| -> Infinity).                  *)
(* Hence the Vee-dependence cancels out of I1 and I2 entirely: the       *)
(* phase, not the upper cutoff, provides the UV(-in-p^+) regulator, and  *)
(* the only surviving cutoff dependence is the Lam rapidity log in I2.   *)
(* I3 keeps a non-decaying boundary oscillation E^(-I omega vee/kp) that *)
(* averages to zero against any smooth transverse profile.               *)

PPlusKernelAsymptotic[omega_, kp_, lam_, vee_] :=
  With[{s = Sign[omega]},
    (1/(2 Pi)) {
      -Exp[I omega] (ExpIntegralEi[-I omega] + I Pi s),
      -EulerGamma - Log[Abs[omega] lam/kp] - I Pi s/2,
      (1 - Exp[-I omega (vee - kp)/kp])/(I omega)
    }
  ];

(* Phase switched off (omega -> 0). Note I3 ~ vee/kp: linearly divergent *)
(* without the phase, so the last delta structure must NOT be evaluated  *)
(* in the eikonal/collinear approximation.                               *)

PPlusKernelCollinear[kp_, lam_, vee_] :=
  (1/(2 Pi)) {Log[vee/(kp + lam)], Log[(vee - kp)/lam], (vee - kp - lam)/kp};

RapidityLog[omega_, kp_, lam_, vee_] :=
  If[TrueQ[PossibleZeroQ[omega]],
    Log[(vee - kp)/lam],
    Log[kp/(Abs[omega] lam)] - EulerGamma - I Pi Sign[omega]/2
  ];

(* ------------------------------------------------------------------ *)
(* 2. Transverse structure                                              *)
(* ------------------------------------------------------------------ *)

(* Contracting the three delta structures with                          *)
(*   (y'-z)^m (x'-y')^k' (y-w)^i (x-z)^j                                 *)
(* gives                                                                 *)
(*   d_{k'm} d_{ij} -> [(y'-z).(x'-y')] [(y-w).(x-z)]                    *)
(*   d_{jm} d_{ik'} -> [(y'-z).(x-z)]   [(y-w).(x'-y')]                  *)
(*   d_{im} d_{jk'} -> [(y'-z).(y-w)]   [(x-z).(x'-y')]                  *)

SetAttributes[TDot, Orderless];
TDot[u_List, v_List] := u . v;

dot[u_, v_] := If[ListQ[u] && ListQ[v], u . v, TDot[u, v]];
sq[u_] := dot[u, u];
sub[u_List, v_List] := u - v;
sub[u_, v_] := TSub[u, v];

TransverseData[k_, w_, yp_, z_, y_, x_, xp_] :=
  Module[{ypz, xpyp, yw, xz},
    ypz  = sub[yp, z];
    xpyp = sub[xp, yp];
    yw   = sub[y, w];
    xz   = sub[x, z];
    <|
      "ypz" -> ypz, "xpyp" -> xpyp, "yw" -> yw, "xz" -> xz,
      (* numerators paired exactly as the three bracket terms demand *)
      "N1" -> dot[ypz, xpyp] dot[yw, xz],   (* d_{k'm} d_{ij}  *)
      "N2" -> dot[ypz, xz] dot[yw, xpyp],   (* d_{jm}  d_{ik'} *)
      "N3" -> dot[ypz, yw] dot[xz, xpyp],   (* d_{im}  d_{jk'} *)
      "Denominator" -> sq[ypz] sq[xpyp] sq[yw] sq[xz],
      "Phase" -> Exp[-I dot[k, sub[yp, w]]],
      "omega" -> dot[k, ypz]
    |>
  ];

TransverseFactor[data_Association, {i1_, i2_, i3_}] :=
  data["Phase"] (data["N1"] i1 - data["N2"] i2 - data["N3"] i3) /
    data["Denominator"];

(* ------------------------------------------------------------------ *)
(* 3. Color structure (inert, untouched by the p^+ integration)         *)
(* ------------------------------------------------------------------ *)

ColorStructure[idx_Association] :=
  Module[{a, b, c, d, e, cp, dp, ep, w, yp, z, y, x, xp},
    {a, b, c, d, e} = Lookup[idx, {"a", "b", "c", "d", "e"}];
    {cp, dp, ep} = Lookup[idx, {"cp", "dp", "ep"}];
    {w, yp, z, y, x, xp} = Lookup[idx, {"w", "yp", "z", "y", "x", "xp"}];
    FStruct[cp, dp, a] *
      NonCommutativeMultiply[
        (U[ep, cp][yp] ** U[dp, b][z] ** Rho[ep][xp]
           - U[b, dp][z] ** U[cp, ep][xp] ** Rho[ep][xp]),
        (U[b, c][z] ** U[a, d][y] ** U[c, e][x] ** Rho[d][y] ** Rho[e][x]
           - U[a, d][y] ** Rho[d][y] ** Rho[b][x])
      ]
  ];

(* ------------------------------------------------------------------ *)
(* 4. Assembly                                                          *)
(* ------------------------------------------------------------------ *)

Options[NLOIntegrand] = {"Kernel" -> "Exact"};

NLOIntegrand[k_, w_, yp_, z_, y_, x_, xp_, kp_, lam_, vee_, g_, idx_Association,
    opts : OptionsPattern[]] :=
  Module[{data, omega, kern, pref},
    data = TransverseData[k, w, yp, z, y, x, xp];
    omega = data["omega"];
    kern = Switch[OptionValue["Kernel"],
      "Exact",      PPlusKernel[omega, kp, lam, vee],
      "Asymptotic", PPlusKernelAsymptotic[omega, kp, lam, vee],
      "Collinear",  PPlusKernelCollinear[kp, lam, vee],
      _, Message[NLOIntegrand::kern, OptionValue["Kernel"]]; $Failed
    ];
    If[kern === $Failed, Return[$Failed]];
    pref = -(1/(2 Pi)^3) (I g^4/(4 Pi^4)) (1/kp);
    pref TransverseFactor[data, kern] ColorStructure[idx]
  ];

NLOIntegrand::kern = "Unknown kernel `1`; use \"Exact\", \"Asymptotic\" or \"Collinear\".";

NLOResult[k_, w_, yp_, z_, y_, x_, xp_, kp_, lam_, vee_, g_, idx_Association,
    opts : OptionsPattern[NLOIntegrand]] :=
  Inactive[Integrate][
    NLOIntegrand[k, w, yp, z, y, x, xp, kp, lam, vee, g, idx, opts],
    {w, yp, z, y, x, xp}
  ];

(* ------------------------------------------------------------------ *)
(* 5. Self-check                                                        *)
(* ------------------------------------------------------------------ *)

CheckPPlusKernel[] :=
  Module[{cases, num, exact, derived, asy, dExact, dDerived, dAsy, p},
    cases = {
      {0.37, 1.0, 0.013, 260.0}, {-0.37, 1.0, 0.013, 260.0},
      {2.5, 1.0, 0.013, 260.0},  {-1.1, 3.0, 0.05, 900.0},
      {0.9, 2.0, 0.002, 5000.0}
    };
    num[{om_, kp_, lam_, vee_}] :=
      (1/(2 Pi)) NIntegrate[
        Exp[-I p om/kp] {1/(p + kp), 1/p, 1/kp},
        {p, lam, vee - kp}, Method -> "LevinRule", MaxRecursion -> 60,
        AccuracyGoal -> 12, PrecisionGoal -> 10];
    exact   = PPlusKernel @@@ cases;
    derived = PPlusKernelDerive @@@ cases;
    asy     = PPlusKernelAsymptotic @@@ cases;
    dExact   = Max[Abs[Flatten[exact - num /@ cases]]];
    dDerived = Max[Abs[Flatten[N[derived] - exact]]];
    (* asymptotics: compare I1, I2 only; I3 keeps a boundary oscillation *)
    dAsy = Max[Abs[Flatten[
      (PPlusKernelAsymptotic[#[[1]], #[[2]], 10^-8, 10^9] -
        PPlusKernel[#[[1]], #[[2]], 10^-8, 10^9])[[{1, 2}]] & /@ cases]]];
    <|"ExactVsNIntegrate" -> dExact,
      "ExactVsIntegrate" -> dDerived,
      "AsymptoticVsExact" -> dAsy|>
  ];

End[];
EndPackage[];
