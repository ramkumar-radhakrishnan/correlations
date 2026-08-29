#!/usr/bin/env python3
"""Generate nlo/pplus_scalar_integral.nb.  Edit CELLS and re-run."""
import os
from make_notebook import build

CELLS = [
("Title", "p^+ integration of  E^(-I a p/k) (A/(p+k) - B/p - 1/k)"),
("Text",
 "Evaluate top to bottom (Evaluation > Evaluate Notebook).\n\n"
 "  Integrate[Exp[-I a p/k] (A/(p + k) - B/p - 1/k), {p, Lam, V - k}]\n\n"
 "in the window V >> k >> Lam > 0.  ASCII stand-ins:\n\n"
 "  p   = p^+          k   = k^+        Lam = Lambda      V = \\[Vee]\n"
 "  a   = k.(y'-z), the transverse phase, so the phase is E^(-I a p^+/k^+)\n"
 "  A, B = the coefficients of the 1/(p^+ + k^+) and 1/p^+ poles\n\n"
 "a is real and may have either sign.  Nothing here assumes |a| is small."),

("Section", "1.  Integrate"),
("Input",
 "as = {a > 0, k > 0, Lam > 0, V > k + Lam};\n"
 "\n"
 "res = Integrate[Exp[-I a p/k] (A/(p + k) - B/p - 1/k), {p, Lam, V - k},\n"
 "   Assumptions -> as];\n"
 "Simplify[res, Assumptions -> as]"),
("Text",
 "The kernel returns\n\n"
 "  A E^(I a) (-Ei[-I a (k + Lam)/k] + Ei[-I a V/k])\n"
 "   + B (Ei[-I a Lam/k] - Ei[I a (k - V)/k])\n"
 "   - I (-1 + E^(I a (k + Lam - V)/k))/(a E^(I a Lam/k))\n\n"
 "with Ei = ExpIntegralEi.  Note Ei[I a (k - V)/k] = Ei[-I a (V - k)/k]."),

("Section", "2.  Closed form"),
("Text",
 "The same thing, one term per pole:\n\n"
 "  Int dp E^(-I a p/k) A/(p+k) =  A E^(I a) [ Ei(-I a V/k) - Ei(-I a (k+Lam)/k) ]\n"
 "  Int dp E^(-I a p/k) (-B/p)  = -B [ Ei(-I a (V-k)/k) - Ei(-I a Lam/k) ]\n"
 "  Int dp E^(-I a p/k) (-1/k)  =  I [ E^(-I a Lam/k) - E^(-I a (V-k)/k) ] / a\n\n"
 "Valid for either sign of a: the ExpIntegralEi arguments sit just off the "
 "imaginary axis, so no branch cut is crossed.  The -I0 that puts them there is "
 "the usual p^+ -> p^+ - I0 prescription."),
("Text",
 "The outer parentheses in the next cell matter: without them the newlines "
 "before - and + terminate the assignment, and Mathematica reads three "
 "separate expressions instead of one sum."),
("Input",
 "PPlusExact[a_, k_, Lam_, V_, A_, B_] := (\n"
 "  A Exp[I a] (ExpIntegralEi[-I a V/k] - ExpIntegralEi[-I a (k + Lam)/k])\n"
 "  - B (ExpIntegralEi[-I a (V - k)/k] - ExpIntegralEi[-I a Lam/k])\n"
 "  + I (Exp[-I a Lam/k] - Exp[-I a (V - k)/k])/a);"),
("Text", "It agrees with what Integrate produced, identically:"),
("Input", "Simplify[PPlusExact[a, k, Lam, V, A, B] - res, Assumptions -> as]"),

("Section", "3.  The regime V >> k >> Lam"),
("Text",
 "Two expansions do the work:\n\n"
 "  Ei(-I x) = EulerGamma + Log[x] - I Pi/2 + O(x)     (small x)\n"
 "  Ei(-I X) -> -I Pi Sign[X]                          (|X| -> Infinity)\n\n"
 "so with |a| Lam/k << 1 << |a| V/k,"),
("Input",
 "PPlusAsymptotic[a_, k_, Lam_, V_, A_, B_] :=\n"
 "  With[{s = Sign[a]},\n"
 "    -A Exp[I a] (ExpIntegralEi[-I a] + I Pi s)\n"
 "    - B (Log[k/(Abs[a] Lam)] - EulerGamma - I Pi s/2)\n"
 "    + I (1 - Exp[-I a (V - k)/k])/a\n"
 "  ];\n"
 "\n"
 "PPlusAsymptotic[a, k, Lam, V, A, B]"),
("Text",
 "Three things fall out.\n\n"
 "(1) V has cancelled out of the A and B terms.  The phase, not the upper cutoff, "
 "is what regulates them once |a| V/k >> 1.\n\n"
 "(2) The only cutoff dependence left is the single logarithm Log[k/(|a| Lam)] on "
 "the B term -- the rapidity log rides on the 1/p^+ pole alone.\n\n"
 "(3) The -1/k term keeps a boundary oscillation E^(-I a (V-k)/k) of unit modulus: "
 "no pointwise V -> Infinity limit, it only averages away against a smooth "
 "transverse profile.  And at a -> 0 that term is -(V - k - Lam)/k, linearly "
 "divergent -- it is finite only by virtue of the phase, so do not drop the phase "
 "there.  For reference, the whole integral at a = 0 is:"),
("Input",
 "Integrate[A/(p + k) - B/p - 1/k, {p, Lam, V - k}, Assumptions -> as]"),

("Section", "4.  Checks"),
("Text",
 "Closed form against direct numerical integration, for both signs of a, and the "
 "asymptotic form against the exact one at Lam = 10^-8, V = 10^9.  Expect exact "
 "zeros in the first row and O(Lam, k/(a V)) ~ 10^-8 in the second."),
("Input",
 "tests = {{0.7, 1.3, 0.011, 310., 1.9, -0.6},\n"
 "         {-0.7, 1.3, 0.011, 310., 1.9, -0.6},\n"
 "         {2.2, 0.9, 0.004, 2000., 1., 1.},\n"
 "         {0.35, 2.0, 0.02, 800., -1.4, 0.8}};\n"
 "\n"
 "numeric[a_, k_, Lam_, V_, A_, B_] :=\n"
 "  NIntegrate[Exp[-I a p/k] (A/(p + k) - B/p - 1/k), {p, Lam, V - k},\n"
 "    Method -> \"LevinRule\", MaxRecursion -> 60, AccuracyGoal -> 12,\n"
 "    PrecisionGoal -> 10];\n"
 "\n"
 "<|\n"
 "  \"exact - NIntegrate\" -> Chop[(PPlusExact @@ #) - (numeric @@ #) & /@ tests],\n"
 "  \"asymptotic - exact\" ->\n"
 "    ((PPlusAsymptotic[#[[1]], #[[2]], 10^-8, 10^9, #[[5]], #[[6]]] -\n"
 "       PPlusExact[#[[1]], #[[2]], 10^-8, 10^9, #[[5]], #[[6]]]) & /@ tests)\n"
 "|>"),
("Input", "PPlusExact @@ tests[[1]]"),
("Text",
 "Exact vs asymptotic as a function of the phase a, at k = 1.3, Lam = 0.011, "
 "V = 310, A = 1.9, B = -0.6.  The two curves separate at small a, where "
 "|a| V/k >> 1 fails and the V logarithm has not yet been traded for the a "
 "logarithm; the wiggle in the exact curve is the boundary oscillation of (3)."),
("Input",
 "With[{k0 = 1.3, Lam0 = 0.011, V0 = 310., A0 = 1.9, B0 = -0.6},\n"
 "  Plot[\n"
 "    Evaluate[{Re[PPlusExact[a, k0, Lam0, V0, A0, B0]],\n"
 "              Re[PPlusAsymptotic[a, k0, Lam0, V0, A0, B0]],\n"
 "              Im[PPlusExact[a, k0, Lam0, V0, A0, B0]],\n"
 "              Im[PPlusAsymptotic[a, k0, Lam0, V0, A0, B0]]}],\n"
 "    {a, 0.02, 3},\n"
 "    PlotStyle -> {Automatic, Directive[Dashed, Thick], Automatic,\n"
 "      Directive[Dashed, Thick]},\n"
 "    PlotLegends -> {\"Re exact\", \"Re asymptotic\", \"Im exact\",\n"
 "      \"Im asymptotic\"},\n"
 "    AxesLabel -> {\"a\", None}, PlotRange -> All, ImageSize -> 560]]"),
]

if __name__ == "__main__":
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "pplus_scalar_integral.nb")
    with open(path, "w") as fh:
        fh.write(build(CELLS))
    print("wrote", path)
