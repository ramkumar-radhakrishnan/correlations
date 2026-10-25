"""Group II Row IV: the lower end xi -> k+/vee of int dxi/xi^2 P(xi) LO_xi[l_Delta] gives a log(vee/k+).
At small xi, LO_xi probes transverse momentum q = k/xi, and the only power-law tail comes from log(Delta^2):
    int d^2Delta e^{-i q.Delta} log(Delta^2) h(Delta) -> -4 pi h(0)/q^2,   so  LO_xi[l_Delta] ~ xi^2 and
    int_a dxi/xi^2 (2N_c/xi) xi^2 -> 2 N_c log(1/a) = 2 N_c log(vee/k+).
Check of the asymptotic for h = exp(-Delta^2/2) and h = 1/(1+Delta^2)^2 (radial):  q^2 x (2 pi int r J0(q r) log(r^2) h(r) dr)."""
import json
import mpmath as mp

mp.mp.dps = 20
out = {}
for name, h in (('gauss', lambda r: mp.e**(-r**2 / 2)), ('lorentz', lambda r: 1 / (1 + r**2)**2)):
    vals = []
    for q in (10, 30, 100):
        f = lambda r: r * mp.besselj(0, q * r) * mp.log(r**2) * h(r)
        I = 2 * mp.pi * mp.quadosc(f, [0, mp.inf], omega=q)
        vals.append(dict(q=q, q2I_over_minus4pi=float(q**2 * I / (-4 * mp.pi))))
    out[name] = vals
    print(name, vals)
json.dump(out, open('ft_asym.json', 'w'), indent=1)
