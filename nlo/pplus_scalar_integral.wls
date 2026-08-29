#!/usr/bin/env wolframscript
(* Scalar form of the NLO p^+ integral:

     Integrate[Exp[-I a p/k] (A/(p + k) - B/p - 1/k), {p, Lam, V - k}]

   in the window V >> k >> Lam > 0, with a = k.(y'-z) the transverse phase,
   k = k^+, Lam = Lambda, V = \[Vee].  Run: wolframscript -file pplus_scalar_integral.wls *)

as = {a > 0, k > 0, Lam > 0, V > k + Lam};

Print["--- Integrate ---"];
res = Integrate[Exp[-I a p/k] (A/(p + k) - B/p - 1/k), {p, Lam, V - k},
   Assumptions -> as];
Print[Simplify[res, Assumptions -> as]];

(* Closed form, rearranged.  Valid for either sign of a: the arguments of
   ExpIntegralEi sit just off the imaginary axis, so no branch cut is crossed. *)
(* the outer parentheses matter: without them the newlines before - and +
   terminate the assignment and this parses as three separate expressions *)
PPlusExact[a_, k_, Lam_, V_, A_, B_] := (
  A Exp[I a] (ExpIntegralEi[-I a V/k] - ExpIntegralEi[-I a (k + Lam)/k])
  - B (ExpIntegralEi[-I a (V - k)/k] - ExpIntegralEi[-I a Lam/k])
  + I (Exp[-I a Lam/k] - Exp[-I a (V - k)/k])/a);

(* Strongly ordered limit, |a| Lam/k << 1 << |a| V/k, from
     Ei[-I x] = EulerGamma + Log[x] - I Pi/2 + O[x],
     Ei[-I X] -> -I Pi Sign[X]  as |X| -> Infinity.                      *)
PPlusAsymptotic[a_, k_, Lam_, V_, A_, B_] :=
  With[{s = Sign[a]},
    -A Exp[I a] (ExpIntegralEi[-I a] + I Pi s)
    - B (Log[k/(Abs[a] Lam)] - EulerGamma - I Pi s/2)
    + I (1 - Exp[-I a (V - k)/k])/a
  ];

Print["\n--- exact vs NIntegrate, and asymptotic vs exact ---"];
tests = {{0.7, 1.3, 0.011, 310., 1.9, -0.6}, {-0.7, 1.3, 0.011, 310., 1.9, -0.6},
         {2.2, 0.9, 0.004, 2000., 1., 1.}, {0.35, 2.0, 0.02, 800., -1.4, 0.8}};
numeric[a_, k_, Lam_, V_, A_, B_] :=
  NIntegrate[Exp[-I a p/k] (A/(p + k) - B/p - 1/k), {p, Lam, V - k},
    Method -> "LevinRule", MaxRecursion -> 60, AccuracyGoal -> 12,
    PrecisionGoal -> 10];

Print["exact - NIntegrate:      ", Chop[(PPlusExact @@ #) - (numeric @@ #) & /@ tests]];
Print["asymptotic - exact:      ",
  (PPlusAsymptotic[#[[1]], #[[2]], 10^-8, 10^9, #[[5]], #[[6]]] -
    PPlusExact[#[[1]], #[[2]], 10^-8, 10^9, #[[5]], #[[6]]]) & /@ tests];
