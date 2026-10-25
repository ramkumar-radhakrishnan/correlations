"""Origin of log^2(1-a) in Group I Row II part 1 (Calculation_main.tex, |Psi^A_grho> + |Psi^B_grho>).

|Psi^B> after q+ = xi k+ (measured gluon q+, parent q+/xi < vee  =>  xi > a = q+/vee;
upper end (1-xi)/xi > Lambda/q+  =>  xi < 1/(1+delta), delta = Lambda/q+):
   int_a^{1/(1+delta)} dxi/(xi(1-xi)) [ (xi-2) ln xi + (xi+1) ln(1-xi) - (xi-2) P_eps ].
|Psi^A>: P_eps [11/3 + log(Lambda/q+) - 2 log(vee/Lambda)] + pi^2/3 - 67/9 + log^2(Lambda/q+) - log^2(vee/q+).
Target (Row II part 1): P_eps [11/3 + 2 Lv + log(1-a) - 2 log(vee/Lambda)] + pi^2/3 - 67/9 + log^2(1-a) + Lv log(1-a).
Then the same with the xi > a limit kept only for the 2/xi-singular terms (AABKL-like): log(1-a), log^2(1-a),
Lv log(1-a) all disappear."""
import json
import mpmath as mp

mp.mp.dps = 30
out = {}
rows = []
for a, delta in ((0.3, 1e-7), (0.05, 1e-8), (0.7, 1e-7), (0.9, 1e-9)):
    a, delta = mp.mpf(a), mp.mpf(delta)
    up = 1 / (1 + delta)
    Lv = -mp.log(a)
    lnd = mp.log(delta)                      # log(Lambda/q+)
    fin_B = mp.quad(lambda x: ((x - 2) * mp.log(x) + (x + 1) * mp.log(1 - x)) / (x * (1 - x)), [a, 0.5, up])
    pol_B = mp.quad(lambda x: -(x - 2) / (x * (1 - x)), [a, 0.5, up])
    vL = Lv - lnd                            # log(vee/Lambda)
    fin_A = mp.pi**2 / 3 - mp.mpf(67) / 9 + lnd**2 - Lv**2
    pol_A = mp.mpf(11) / 3 + lnd - 2 * vL
    fin_t = mp.pi**2 / 3 - mp.mpf(67) / 9 + mp.log(1 - a)**2 + Lv * mp.log(1 - a)
    pol_t = mp.mpf(11) / 3 + 2 * Lv + mp.log(1 - a) - 2 * vL
    # restriction xi > a only on the 2/xi-singular terms (-2 ln xi/xi, 2/xi); the rest from 0
    reg = lambda x: (-mp.log(x) / (1 - x) + mp.log(1 - x) / x + 2 * mp.log(1 - x) / (1 - x))
    fin_B0 = mp.quad(lambda x: -2 * mp.log(x) / x, [a, 1]) + mp.quad(reg, [0, 0.5, up])
    pol_B0 = mp.quad(lambda x: 2 / x, [a, 1]) + mp.quad(lambda x: 1 / (1 - x), [0, 0.5, up])
    fin_0 = mp.pi**2 / 3 - mp.mpf(67) / 9
    pol_0 = mp.mpf(11) / 3 + 2 * Lv - 2 * vL
    rows.append(dict(a=float(a), finite=float(abs(fin_A + fin_B - fin_t)), pole=float(abs(pol_A + pol_B - pol_t)),
                     finite_unrestricted=float(abs(fin_A + fin_B0 - fin_0)), pole_unrestricted=float(abs(pol_A + pol_B0 - pol_0))))
out['psiA_plus_psiB_vs_rowII'] = rows
# the four pieces of the finite integral, closed form, at a = 0.3
a = mp.mpf('0.3')
pieces = {
    '-2 ln xi/xi  -> Lv^2': (mp.quad(lambda x: -2 * mp.log(x) / x, [a, 1]), mp.log(a)**2),
    '-ln xi/(1-xi) -> Li2(1-a)': (mp.quad(lambda x: -mp.log(x) / (1 - x), [a, 1]), mp.polylog(2, 1 - a)),
    'ln(1-xi)/xi -> Li2(a)-pi^2/6': (mp.quad(lambda x: mp.log(1 - x) / x, [a, 1]), mp.polylog(2, a) - mp.pi**2 / 6),
}
out['pieces'] = {k: float(abs(v[0] - v[1])) for k, v in pieces.items()}
# 2 ln(1-xi)/(1-xi) from a to 1-d  = ln^2(1-a) - ln^2 d
d = mp.mpf('1e-6')
out['pieces']['2 ln(1-xi)/(1-xi) -> ln^2(1-a) - ln^2 delta'] = float(abs(
    mp.quad(lambda x: 2 * mp.log(1 - x) / (1 - x), [a, 0.9, 1 - d]) - (mp.log(1 - a)**2 - mp.log(d)**2)))
# compact form and small-a behaviour
for aa in (1e-2, 1e-4, 1e-6):
    aa = mp.mpf(aa)
    lhs = -mp.log(1 - aa)**2 + mp.log(aa) * mp.log(1 - aa)
    rhs = -mp.log(1 - aa) * mp.log((1 - aa) / aa)
    out.setdefault('compact', []).append(dict(a=float(aa), diff=float(abs(lhs - rhs)), value=float(lhs),
                                              value_over_a_log=float(lhs / (aa * mp.log(1 / aa)))))
json.dump(out, open('log2_origin.json', 'w'), indent=1)
print(json.dumps(out, indent=1))
