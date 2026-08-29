#!/usr/bin/env wolframscript
(* ::Package:: *)

(* Demo / verification driver for NLOPplusIntegral`.
   Run with:  wolframscript -file nlo/nlo_pplus_demo.wls
   (or, in a notebook:  Get["nlo/NLOPplusIntegral.wl"] and evaluate the cells.) *)

SetDirectory[DirectoryName[$InputFileName]];
Get["NLOPplusIntegral.wl"];

Print["=== 1. The p^+ integral, done symbolically ==========================="];
Print["Integrate under Vee > k^+ > Lam > 0, with b = omega/k^+ = k.(y'-z)/k^+:"];
Print[
  TableForm[
    Transpose[{
      {"1/(p^+ + k^+)", "1/p^+", "1/k^+"},
      PPlusKernelDerive[omega, kp, lam, vee] 2 Pi
    }],
    TableHeadings -> {None, {"p^+ pole", "2 Pi x Integrate[dp^+ E^(-I p^+ omega/k^+) x pole]"}}
  ]
];

Print["\nSame thing quoted in closed form by PPlusKernel (times 2 Pi):"];
Print[Simplify[2 Pi PPlusKernel[omega, kp, lam, vee]]];

Print["\n=== 2. Strongly ordered limit Vee >> k^+ >> Lam ====================="];
Print["Asymptotic kernel (times 2 Pi):"];
Print[Simplify[2 Pi PPlusKernelAsymptotic[omega, kp, lam, vee]]];
Print["\nThe rapidity logarithm sits entirely in the 1/p^+ term:"];
Print["  ", RapidityLog[omega, kp, lam, vee]];
Print["Vee has dropped out of I1 and I2: the phase E^(-I p^+ omega/k^+), not the"];
Print["upper cutoff, regulates them.  Only Log[1/Lam] survives, as it must for a"];
Print["JIMWLK-type rapidity evolution to absorb it."];

Print["\nWith the phase switched off (omega -> 0) the Vee log comes back and the"];
Print["last structure diverges linearly -- do not drop the phase there:"];
Print["  ", 2 Pi PPlusKernelCollinear[kp, lam, vee]];

Print["\n=== 3. Numerical audit =============================================="];
Print[CheckPPlusKernel[]];

Print["\n=== 4. Full integrand after the p^+ integration ====================="];
idx = <|"a" -> a, "b" -> b, "c" -> c, "d" -> d, "e" -> e,
        "cp" -> cp, "dp" -> dp, "ep" -> ep,
        "w" -> w, "yp" -> yp, "z" -> z, "y" -> y, "x" -> x, "xp" -> xp|>;
res = NLOResult[k, w, yp, z, y, x, xp, kp, lam, vee, g, idx];
Print[res];

Print["\nTransverse contractions the three delta structures collapse to:"];
Print[
  TableForm[
    {{"d_{k'm} d_{ij}", (TransverseData[k, w, yp, z, y, x, xp])["N1"]},
     {"d_{jm} d_{ik'}", (TransverseData[k, w, yp, z, y, x, xp])["N2"]},
     {"d_{im} d_{jk'}", (TransverseData[k, w, yp, z, y, x, xp])["N3"]}},
    TableHeadings -> {None, {"delta structure", "transverse contraction"}}
  ]
];

Print["\n=== 5. A fully numerical point ======================================"];
pt = {{0.8, -0.5}, {0.1, 0.2}, {1.3, -0.4}, {-0.6, 0.9},
      {0.35, 1.1}, {-1.2, -0.7}, {0.9, 0.15}};
dat = TransverseData @@ pt;
Print["omega = k.(y'-z) = ", dat["omega"]];
Print["transverse factor x p^+ integral = ",
  TransverseFactor[dat, PPlusKernel[dat["omega"], 1.0, 0.01, 400.0]]];
