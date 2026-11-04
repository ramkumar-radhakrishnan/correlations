r"""The coefficient of log mu^2 of every source in your rows, in units of LO (alpha_s = g^2/4pi), and their sum.

Sources (your final equations):
  A3            Group I Row II part 1      (g^4 N_c/32pi^4)[-E_X][11/3 - 2L + l1] ...        and its c.c.
  shockwave     Group I Row V 2b UV        (g^4 N_c/16pi^4)[2L + 2lam - 11/6][+E_(y-z)] ...  and its c.c.
  real          Group II Row IV lines 1, 2 (g^2/4pi^2) N_c [-E_Delta] (xi integral; and g(1)[L + log(vee/k+ - 1)])
  coupling CT   g^2 of LO in MS-bar
E_r = 2/eps + log(r^2 mu^2/c0^2): every pole comes with log mu^2 with the same coefficient (times eps/2), so the
log mu^2 bookkeeping is the pole bookkeeping.  L = log(k+/Lambda), l1 = log(1 - k+/vee), lam = log(1 - Lambda/k+),
a = k+/vee.  P_user: your xi/(1-xi)_+ (subtraction g(1) on [a,1]); P_std: the MS-bar plus distribution of [0,1]
restricted to [a,1].  The relation P_user = P_std - 2N_c [l1 + log(vee/k+)] delta(1-xi) is checked numerically.
Output: mu_check.json"""
import json
import numpy as np
import sympy as sp
from scipy import integrate

Nc, L, l1, lam, a = sp.symbols('N_c L l1 lambda a', positive=True)
Lv = sp.Symbol('Lv', positive=True)    # log(vee/k+)
PG = sp.Symbol('Pstar')              # int_a^1 dxi/xi^2 P^(r)_std(xi) LO_xi  /  LO   (a symbol: the convolution)
b0 = sp.Rational(11, 3) * Nc
out = {}

# prefactor conversions to LO = 1/(2pi)^3 g^2/(2pi^2) 1/k+ int ...  (alpha_s = g^2/4pi)
g = sp.Symbol('g', positive=True)
LOpref = 1 / (2 * sp.pi)**3 * g**2 / (2 * sp.pi**2)
conv = {'A3': sp.simplify(1 / (2 * sp.pi)**3 * g**4 * Nc / (32 * sp.pi**4) / LOpref),
        'shockwave': sp.simplify(1 / (2 * sp.pi)**3 * g**4 * Nc / (16 * sp.pi**4) / LOpref),
        'real': sp.simplify(g**2 / (4 * sp.pi**2) * Nc)}
als = g**2 / (4 * sp.pi)
out['prefactor / LO prefactor'] = {k: '%s = %s' % (v, sp.simplify(v / als)) + ' alpha_s' for k, v in conv.items()}
# A3: colour [U(x')-U(w')][U(w)-U(x)] = -(LO colour after y'->w', y->w); bracket -E_X: sign (+); with c.c. -> 2 log mu^2
# shockwave: colour [U(x')-U(w')][U(x)-U(w)] after y->x, z->w = +(LO colour); bracket +E: sign (+); with c.c. -> 2
# real: bracket [-d/eps - d gamma - (d/2) log(...mu^2/4) - 1] = -E_Delta + O(eps): sign (-); no c.c.
mu = {}
mu['A3 (Row II part 1) + c.c.'] = 2 * sp.Rational(1, 4) * Nc * (sp.Rational(11, 3) - 2 * L + l1)          # units alpha_s/pi LO
mu['shockwave through the loop (Row V 2b UV) + c.c.'] = 2 * sp.Rational(1, 2) * Nc * (-sp.Rational(11, 6) + 2 * L + 2 * lam)
# real, line 1: -(1/2) [P_user-conv]; P_user-conv = Pstar - 2N_c (l1 + Lv); line 2: -(1/2) 2N_c (L + log(vee/k+ - 1))
mu['real, Group II Row IV line 1'] = -sp.Rational(1, 2) * (PG - 2 * Nc * (l1 + Lv))
mu['real, Group II Row IV line 2'] = -sp.Rational(1, 2) * 2 * Nc * (L + Lv + l1)
mu['coupling counterterm of LO (mu_R)'] = -sp.Rational(1, 2) * b0 / 2 * 2 / 2          # -(alpha_s beta0/4pi) -> units alpha_s/pi

def parts(e):
    e = sp.expand(e)
    return {'beta0 (11/3)': str(sp.simplify(e.subs({L: 0, l1: 0, lam: 0, PG: 0, Lv: 0}))),
            'x log(k+/Lambda)': str(sp.expand(e).coeff(L)), 'x log(1-k+/vee)': str(sp.expand(e).coeff(l1)),
            'x log(vee/k+)': str(sp.expand(e).coeff(Lv)),
            'x P conv': str(sp.expand(e).coeff(PG))}

out['coefficient of log mu^2, units alpha_s/pi x LO'] = {k: parts(v) for k, v in mu.items()}
groupI = mu['A3 (Row II part 1) + c.c.'] + mu['shockwave through the loop (Row V 2b UV) + c.c.']
real = mu['real, Group II Row IV line 1'] + mu['real, Group II Row IV line 2']
out['Group I (4 sources)'] = str(sp.expand(groupI))
out['Group I, beta0 part (L = l1 = lam = 0)'] = str(sp.simplify(groupI.subs({L: 0, l1: 0, lam: 0})))
out['real (1 source)'] = str(sp.expand(real))
out['evolution part, x L: A3 : shockwave : real'] = [str(sp.expand(v).coeff(L)) for v in
                                                    (mu['A3 (Row II part 1) + c.c.'], mu['shockwave through the loop (Row V 2b UV) + c.c.'], real)]
tot5 = sp.expand(groupI + real)
out['sum of the 5 sources'] = str(tot5)
out['sum + coupling counterterm'] = str(sp.expand(tot5 + mu['coupling counterterm of LO (mu_R)']))
out['expected: -(1/2)[Pstar + beta0/2] + (N_c/2) l1 + 2 N_c lam'] = str(sp.simplify(
    tot5 + mu['coupling counterterm of LO (mu_R)'] - (-sp.Rational(1, 2) * (PG + b0 / 2) + Nc / 2 * l1 + 2 * Nc * lam)))


# numerical check of P_user = P_std - [l1 + log(vee/k+)] g(1) for the xi/(1-xi)_+ term (per unit 2N_c)
def check_plus(av, gfun):
    user = integrate.quad(lambda x: (gfun(x) - gfun(1.0)) / (x * (1 - x)), av, 1, limit=400)[0]
    std = integrate.quad(lambda x: (gfun(x) / x - gfun(1.0)) / (1 - x), av, 1, limit=400)[0] + gfun(1.0) * np.log(1 - av)
    return user - (std - gfun(1.0) * (np.log(1 - av) + np.log(1 / av)))


out['P_user - (P_std - g(1)[l1 + log(vee/k+)]), numerical'] = [
    float(check_plus(av, gf)) for av in (0.05, 0.3) for gf in (lambda x: np.cos(3 / x) / x, lambda x: np.exp(-x) * x**2)]
json.dump(out, open('mu_check.json', 'w'), indent=1)
print(json.dumps(out, indent=1))
