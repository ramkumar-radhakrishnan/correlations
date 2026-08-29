#!/usr/bin/env python3
"""Generate nlo/NLO_pplus_integral.nb from the cell list below.

The notebook is self-contained: it re-derives and defines everything, so it does
not need NLOPplusIntegral.wl.  Edit the CELLS list and re-run to regenerate.
"""
import os

def esc(s):
    return s.replace("\\", "\\\\").replace('"', '\\"')

# (style, content).  "Section" cells open a group that runs to the next Section.
CELLS = [
("Title", "Longitudinal p^+ integration of the NLO gluon-production term"),
("Text",
 "Evaluate the notebook top to bottom (Evaluation > Evaluate Notebook).\n\n"
 "We integrate\n\n"
 "  d^3 N_NLO / d^3 k = -(1/(2 Pi)^3) (I g^4 f^{c'd'a} / (4 Pi^4)) (1/k^+)\n"
 "      Int_Lam^{Vee-k^+} dp^+/(2 Pi)  Int_{w,y',z,y,x,x'}\n"
 "      E^(-I k.(y'-w)) E^(-I k.(p^+/k^+)(y'-z))\n"
 "      (y'-z)^m (x'-y')^k' (y-w)^i (x-z)^j\n"
 "        / [ (y'-z)^2 (x'-y')^2 (y-w)^2 (x-z)^2 ]\n"
 "      [ d_{k'm} d_{ij}/(p^+ + k^+) - d_{jm} d_{ik'}/p^+ - d_{im} d_{jk'}/k^+ ]\n"
 "      [ color 1 ] [ color 2 ]\n\n"
 "over p^+, in the window Vee >> k^+ >> Lam > 0.\n\n"
 "The only p^+ dependence is the phase and the three poles.  So the p^+ integral "
 "factorizes completely out of the transverse integrals AND out of the color "
 "operators -- the Wilson lines and charges pass through untouched.  Sections 1-3 "
 "do that integral; sections 4-6 put the result back together."),

("Section", "0.  Notation"),
("Text",
 "ASCII stand-ins for the symbols in the expression:\n\n"
 "  kp    = k^+          omega = k.(y'-z)   (transverse dot product)\n"
 "  lam   = Lambda       b     = omega/kp   (the phase is E^(-I b p^+))\n"
 "  vee   = \\[Vee]        g     = coupling\n\n"
 "Ordering assumed throughout: vee > kp > lam > 0."),
("Input",
 "NLOAssumptions[kp_, lam_, vee_] :=\n"
 "  {kp > 0, lam > 0, vee > kp + lam, Element[{kp, lam, vee}, Reals]};"),

("Section", "1.  The p^+ integral, by brute force"),
("Text",
 "Pull the p^+-dependent factors out.  Each delta structure comes with one pole, "
 "so there are exactly three integrals to do."),
("Input",
 "rawIntegrals =\n"
 "  Integrate[\n"
 "    Exp[-I b p] {1/(p + kp), 1/p, 1/kp},\n"
 "    {p, lam, vee - kp},\n"
 "    Assumptions -> Append[NLOAssumptions[kp, lam, vee], b > 0]\n"
 "  ]"),
("Text",
 "i.e., with b = omega/kp,\n\n"
 "  Int dp^+ E^(-I b p^+)/(p^+ + k^+) = E^(I b k^+) [ Ei(-I b Vee) - Ei(-I b (k^+ + Lam)) ]\n"
 "  Int dp^+ E^(-I b p^+)/p^+         = Ei(-I b (Vee - k^+)) - Ei(-I b Lam)\n"
 "  Int dp^+ E^(-I b p^+)/k^+         = [ E^(-I b Lam) - E^(-I b (Vee-k^+)) ] / (I b k^+)\n\n"
 "ExpIntegralEi is evaluated just off the imaginary axis, so no branch cut is "
 "crossed for either sign of omega; the -I0 that puts it there is the usual "
 "p^+ -> p^+ - I0 prescription."),

("Section", "2.  Closed form"),
("Text",
 "PPlusKernel returns {I1, I2, I3}, each already carrying the 1/(2 Pi) of the "
 "measure.  They multiply d_{k'm} d_{ij}, -d_{jm} d_{ik'} and -d_{im} d_{jk'} "
 "respectively.  omega = 0 is handled separately (section 3)."),
("Input",
 "PPlusKernelCollinear[kp_, lam_, vee_] :=\n"
 "  (1/(2 Pi)) {Log[vee/(kp + lam)], Log[(vee - kp)/lam], (vee - kp - lam)/kp};\n"
 "\n"
 "PPlusKernel[omega_, kp_, lam_, vee_] :=\n"
 "  If[TrueQ[PossibleZeroQ[omega]],\n"
 "    PPlusKernelCollinear[kp, lam, vee],\n"
 "    With[{b = omega/kp},\n"
 "      (1/(2 Pi)) {\n"
 "        Exp[I b kp] (ExpIntegralEi[-I b vee] - ExpIntegralEi[-I b (kp + lam)]),\n"
 "        ExpIntegralEi[-I b (vee - kp)] - ExpIntegralEi[-I b lam],\n"
 "        (Exp[-I b lam] - Exp[-I b (vee - kp)])/(I b kp)\n"
 "      }\n"
 "    ]\n"
 "  ];"),
("Text", "Check it against what Integrate just produced:"),
("Input",
 "Simplify[2 Pi PPlusKernel[b kp, kp, lam, vee] - rawIntegrals,\n"
 "  Assumptions -> Append[NLOAssumptions[kp, lam, vee], b > 0]]"),

("Section", "3.  The regime Vee >> k^+ >> Lam"),
("Text",
 "Two facts drive the limits:\n\n"
 "  Ei(-I x) = EulerGamma + Log[x] - I Pi/2 + O(x)      (small x)\n"
 "  Ei(-I X) -> -I Pi Sign[X]                           (|X| -> Infinity)\n\n"
 "Valid for |omega| lam/kp << 1 << |omega| vee/kp."),
("Input",
 "PPlusKernelAsymptotic[omega_, kp_, lam_, vee_] :=\n"
 "  With[{s = Sign[omega]},\n"
 "    (1/(2 Pi)) {\n"
 "      -Exp[I omega] (ExpIntegralEi[-I omega] + I Pi s),\n"
 "      -EulerGamma - Log[Abs[omega] lam/kp] - I Pi s/2,\n"
 "      (1 - Exp[-I omega (vee - kp)/kp])/(I omega)\n"
 "    }\n"
 "  ];\n"
 "\n"
 "RapidityLog[omega_, kp_, lam_, vee_] :=\n"
 "  If[TrueQ[PossibleZeroQ[omega]],\n"
 "    Log[(vee - kp)/lam],\n"
 "    Log[kp/(Abs[omega] lam)] - EulerGamma - I Pi Sign[omega]/2\n"
 "  ];"),
("Input", "2 Pi PPlusKernelAsymptotic[omega, kp, lam, vee]"),
("Text",
 "Three things to take away, all checked numerically in section 7.\n\n"
 "(1) Vee has dropped out of I1 and I2.  The phase, not the upper cutoff, is what "
 "regulates them once |omega| Vee/k^+ >> 1.  The one surviving cutoff dependence "
 "is the single rapidity logarithm Log[k^+/(|omega| Lam)] in I2 -- so it rides "
 "entirely on the -d_{jm} d_{ik'} structure.  That is the piece a JIMWLK-type "
 "rapidity evolution has to absorb.\n\n"
 "(2) I3 has no pointwise Vee -> Infinity limit: its boundary term keeps "
 "oscillating with unit modulus.  It averages to zero only against a smooth "
 "transverse profile, so keep I3 exact rather than taking a limit under the "
 "y', z integrals.\n\n"
 "(3) Do not switch the phase off in I3.  At omega -> 0 the Vee logarithm returns "
 "in I1 and I2 and I3 -> (Vee - k^+ - Lam)/k^+ diverges LINEARLY -- that term is "
 "finite only by virtue of the phase."),
("Input", "2 Pi PPlusKernelCollinear[kp, lam, vee]"),
("Input", "RapidityLog[omega, kp, lam, vee]"),

("Section", "4.  Transverse structure"),
("Text",
 "Contracting the three delta structures with the numerator\n"
 "  (y'-z)^m (x'-y')^k' (y-w)^i (x-z)^j\n"
 "collapses it to three scalar pairings:\n\n"
 "  d_{k'm} d_{ij} -> [(y'-z).(x'-y')] [(y-w).(x-z)]      (N1)\n"
 "  d_{jm} d_{ik'} -> [(y'-z).(x-z)]   [(y-w).(x'-y')]    (N2)\n"
 "  d_{im} d_{jk'} -> [(y'-z).(y-w)]   [(x-z).(x'-y')]    (N3)\n\n"
 "Points may be given as explicit 2-component vectors (then dot products are "
 "evaluated) or as opaque symbols (then they stay inert as TDot / TSub)."),
("Input",
 "SetAttributes[TDot, Orderless];\n"
 "TDot[u_List, v_List] := u . v;\n"
 "dot[u_, v_] := If[ListQ[u] && ListQ[v], u . v, TDot[u, v]];\n"
 "sq[u_] := dot[u, u];\n"
 "sub[u_List, v_List] := u - v;\n"
 "sub[u_, v_] := TSub[u, v];\n"
 "\n"
 "TransverseData[k_, w_, yp_, z_, y_, x_, xp_] :=\n"
 "  Module[{ypz, xpyp, yw, xz},\n"
 "    ypz = sub[yp, z]; xpyp = sub[xp, yp]; yw = sub[y, w]; xz = sub[x, z];\n"
 "    <|\n"
 "      \"ypz\" -> ypz, \"xpyp\" -> xpyp, \"yw\" -> yw, \"xz\" -> xz,\n"
 "      \"N1\" -> dot[ypz, xpyp] dot[yw, xz],\n"
 "      \"N2\" -> dot[ypz, xz] dot[yw, xpyp],\n"
 "      \"N3\" -> dot[ypz, yw] dot[xz, xpyp],\n"
 "      \"Denominator\" -> sq[ypz] sq[xpyp] sq[yw] sq[xz],\n"
 "      \"Phase\" -> Exp[-I dot[k, sub[yp, w]]],\n"
 "      \"omega\" -> dot[k, ypz]\n"
 "    |>\n"
 "  ];\n"
 "\n"
 "TransverseFactor[d_Association, {i1_, i2_, i3_}] :=\n"
 "  d[\"Phase\"] (d[\"N1\"] i1 - d[\"N2\"] i2 - d[\"N3\"] i3) / d[\"Denominator\"];"),
("Input",
 "TableForm[\n"
 "  Transpose[{{\"d_{k'm} d_{ij}\", \"d_{jm} d_{ik'}\", \"d_{im} d_{jk'}\"},\n"
 "    Lookup[TransverseData[k, w, yp, z, y, x, xp], {\"N1\", \"N2\", \"N3\"}]}],\n"
 "  TableHeadings -> {None, {\"delta structure\", \"transverse contraction\"}}]"),

("Section", "5.  Color structure"),
("Text",
 "Untouched by the p^+ integration -- kept inert.  U[a,b][pt] is the adjoint "
 "Wilson line U^{ab}(pt), Rho[a][pt] the color charge rho^a(pt), FStruct the "
 "structure constant.  Repeated adjoint indices are summed."),
("Input",
 "ColorStructure[idx_Association] :=\n"
 "  Module[{a, b, c, d, e, cp, dp, ep, w, yp, z, y, x, xp},\n"
 "    {a, b, c, d, e} = Lookup[idx, {\"a\", \"b\", \"c\", \"d\", \"e\"}];\n"
 "    {cp, dp, ep} = Lookup[idx, {\"cp\", \"dp\", \"ep\"}];\n"
 "    {w, yp, z, y, x, xp} = Lookup[idx, {\"w\", \"yp\", \"z\", \"y\", \"x\", \"xp\"}];\n"
 "    FStruct[cp, dp, a] *\n"
 "      NonCommutativeMultiply[\n"
 "        (U[ep, cp][yp] ** U[dp, b][z] ** Rho[ep][xp]\n"
 "           - U[b, dp][z] ** U[cp, ep][xp] ** Rho[ep][xp]),\n"
 "        (U[b, c][z] ** U[a, d][y] ** U[c, e][x] ** Rho[d][y] ** Rho[e][x]\n"
 "           - U[a, d][y] ** Rho[d][y] ** Rho[b][x])\n"
 "      ]\n"
 "  ];"),

("Section", "6.  Assembly"),
("Text",
 "The full integrand with the p^+ integration done, transverse integrals left "
 "inert.  Option \"Kernel\" selects \"Exact\" (default), \"Asymptotic\" or "
 "\"Collinear\".\n\n"
 "Note on the measure: as written it is Int_{w,y',z,y,x}, but the integrand "
 "contains x' (in (x'-y') and rho^{e'}(x')).  NLOResult integrates over "
 "{w, y', z, y, x, x'}; if x' is meant to be tied to another point, substitute "
 "for it first."),
("Input",
 "Options[NLOIntegrand] = {\"Kernel\" -> \"Exact\"};\n"
 "\n"
 "NLOIntegrand[k_, w_, yp_, z_, y_, x_, xp_, kp_, lam_, vee_, g_,\n"
 "    idx_Association, OptionsPattern[]] :=\n"
 "  Module[{data, omega, kern, pref},\n"
 "    data = TransverseData[k, w, yp, z, y, x, xp];\n"
 "    omega = data[\"omega\"];\n"
 "    kern = Switch[OptionValue[\"Kernel\"],\n"
 "      \"Exact\", PPlusKernel[omega, kp, lam, vee],\n"
 "      \"Asymptotic\", PPlusKernelAsymptotic[omega, kp, lam, vee],\n"
 "      \"Collinear\", PPlusKernelCollinear[kp, lam, vee]];\n"
 "    pref = -(1/(2 Pi)^3) (I g^4/(4 Pi^4)) (1/kp);\n"
 "    pref TransverseFactor[data, kern] ColorStructure[idx]\n"
 "  ];\n"
 "\n"
 "NLOResult[k_, w_, yp_, z_, y_, x_, xp_, kp_, lam_, vee_, g_,\n"
 "    idx_Association, opts : OptionsPattern[NLOIntegrand]] :=\n"
 "  Inactive[Integrate][\n"
 "    NLOIntegrand[k, w, yp, z, y, x, xp, kp, lam, vee, g, idx, opts],\n"
 "    {w, yp, z, y, x, xp}\n"
 "  ];"),
("Input",
 "idx = <|\"a\" -> a, \"b\" -> b, \"c\" -> c, \"d\" -> d, \"e\" -> e,\n"
 "        \"cp\" -> cp, \"dp\" -> dp, \"ep\" -> ep,\n"
 "        \"w\" -> w, \"yp\" -> yp, \"z\" -> z, \"y\" -> y, \"x\" -> x, \"xp\" -> xp|>;\n"
 "NLOResult[k, w, yp, z, y, x, xp, kp, lam, vee, g, idx]"),

("Section", "7.  Numerical audit"),
("Text",
 "Exact kernel vs direct NIntegrate over p^+, and the asymptotic kernel vs the "
 "exact one (I1, I2 only -- I3 keeps its boundary oscillation).  Expect ~1*^-16 "
 "for the first and ~1*^-9, i.e. O(Lam, k^+/(omega Vee)), for the second."),
("Input",
 "cases = {{0.37, 1.0, 0.013, 260.0}, {-0.37, 1.0, 0.013, 260.0},\n"
 "         {2.5, 1.0, 0.013, 260.0}, {-1.1, 3.0, 0.05, 900.0},\n"
 "         {0.9, 2.0, 0.002, 5000.0}};\n"
 "\n"
 "numKernel[{om_, kp_, lam_, vee_}] :=\n"
 "  (1/(2 Pi)) NIntegrate[Exp[-I p om/kp] {1/(p + kp), 1/p, 1/kp},\n"
 "    {p, lam, vee - kp}, Method -> \"LevinRule\", MaxRecursion -> 60,\n"
 "    AccuracyGoal -> 12, PrecisionGoal -> 10];\n"
 "\n"
 "<|\n"
 "  \"ExactVsNIntegrate\" ->\n"
 "    Max[Abs[Flatten[(PPlusKernel @@@ cases) - (numKernel /@ cases)]]],\n"
 "  \"AsymptoticVsExact\" ->\n"
 "    Max[Abs[Flatten[(PPlusKernelAsymptotic[#[[1]], #[[2]], 10^-8, 10^9] -\n"
 "        PPlusKernel[#[[1]], #[[2]], 10^-8, 10^9])[[{1, 2}]] & /@ cases]]]\n"
 "|>"),
("Text",
 "And the reassembled transverse factor against a direct one-dimensional "
 "integration of the whole bracket, at an arbitrary configuration of points:"),
("Input",
 "{kv, wv, ypv, zv, yv, xv, xpv} =\n"
 "  {{0.8, -0.5}, {0.1, 0.2}, {1.3, -0.4}, {-0.6, 0.9},\n"
 "   {0.35, 1.1}, {-1.2, -0.7}, {0.9, 0.15}};\n"
 "{kp0, lam0, vee0} = {1.0, 0.01, 400.0};\n"
 "dat = TransverseData[kv, wv, ypv, zv, yv, xv, xpv];\n"
 "\n"
 "assembled = TransverseFactor[dat, PPlusKernel[dat[\"omega\"], kp0, lam0, vee0]];\n"
 "\n"
 "direct = Exp[-I kv . (ypv - wv)]/\n"
 "    (((ypv - zv) . (ypv - zv)) ((xpv - ypv) . (xpv - ypv))\n"
 "      ((yv - wv) . (yv - wv)) ((xv - zv) . (xv - zv))) *\n"
 "  NIntegrate[\n"
 "    Exp[-I p (kv . (ypv - zv))/kp0]/(2 Pi) (\n"
 "      ((ypv - zv) . (xpv - ypv)) ((yv - wv) . (xv - zv))/(p + kp0)\n"
 "      - ((ypv - zv) . (xv - zv)) ((yv - wv) . (xpv - ypv))/p\n"
 "      - ((ypv - zv) . (yv - wv)) ((xv - zv) . (xpv - ypv))/kp0),\n"
 "    {p, lam0, vee0 - kp0}, Method -> \"LevinRule\", MaxRecursion -> 60,\n"
 "    AccuracyGoal -> 12];\n"
 "\n"
 "<|\"omega\" -> dat[\"omega\"], \"assembled\" -> assembled, \"direct\" -> direct,\n"
 "  \"difference\" -> Abs[assembled - direct]|>"),
]


def build():
    out, i, n = [], 0, len(CELLS)
    while i < n:
        style, content = CELLS[i]
        if style in ("Title", "Section"):
            group = ['Cell["%s", "%s"]' % (esc(content), style)]
            i += 1
            while i < n and CELLS[i][0] != "Section" and not (
                    style == "Section" and CELLS[i][0] == "Title"):
                s2, c2 = CELLS[i]
                group.append('Cell["%s", "%s"]' % (esc(c2), s2))
                i += 1
            out.append("Cell[CellGroupData[{\n" + ",\n".join(group) + "\n}, Open]]")
        else:
            out.append('Cell["%s", "%s"]' % (esc(content), style))
            i += 1
    body = ",\n".join(out)
    return ("(* Content-type: application/vnd.wolfram.mathematica *)\n"
            "(* Wolfram Notebook File *)\n"
            "(* http://www.wolfram.com/nb *)\n\n"
            "Notebook[{\n" + body + "\n},\n"
            'WindowSize->{1100, 850},\n'
            "WindowMargins->Automatic,\n"
            'StyleDefinitions->"Default.nb"\n]\n')


if __name__ == "__main__":
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "NLO_pplus_integral.nb")
    with open(path, "w") as fh:
        fh.write(build())
    print("wrote", path)
