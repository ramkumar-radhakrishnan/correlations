"""Group II Row IV with general d_perp: checks.
1. polarization sum sum_{ij} V_{ij;km} V_{ij;k'm'} in d = 2..5 against the closed form
2. the tensor z integral J^{mm'} = int d^dz K_d^m(y-z) K_d^{m'}(y'-z), K_d the d-dim Fourier transform of q^m/q^2:
   (a) d = 2: the traceless part (finite) by quadrature
   (b) d = 3: the full tensor by quadrature
3. the contraction with the bracket (sympy, symbolic d) and the reduction to the earlier result at d = 2 - eps
4. the exact endpoint (plus) decomposition at xi -> 1 (quadrature)"""
import json
import numpy as np
import sympy as sp
from scipy import integrate, special

out = {}
rng = np.random.default_rng(11)

# 1. polarization sum
res = []
for d in (2, 3, 4, 5):
    I = np.eye(d)
    for t in range(2):
        p, k = rng.uniform(0.2, 2, size=2)
        P = p + k
        V = (np.einsum('km,ij->ijkm', I, I) - P / p * np.einsum('jm,ik->ijkm', I, I) - P / k * np.einsum('im,kj->ijkm', I, I))
        S = np.einsum('ijkm,ijlq->kmlq', V, V)                    # [k, m, k', m']
        F = ((d - 2 * P / p - 2 * P / k) * np.einsum('km,lq->kmlq', I, I) + ((P / p)**2 + (P / k)**2) * np.einsum('mq,kl->kmlq', I, I)
             + 2 * P**2 / (p * k) * np.einsum('ml,kq->kmlq', I, I))
        res.append(float(np.abs(S - F).max()))
out['polarization_sum_max_diff'] = max(res)


# 2. tensor integral
def Kd(r, d):
    return np.pi**(1 - d / 2) * special.gamma(d / 2) * r / np.linalg.norm(r)**d


def J_closed(D, d):
    """(pi^{2-d/2}/2) (D^2)^{1-d/2} [Gamma(d/2-1) delta - 2 Gamma(d/2) Dhat Dhat]  (no MS-bar factor)"""
    n = np.linalg.norm(D)
    Dh = D / n
    return np.pi**(2 - d / 2) / 2 * n**(2 - d) * (special.gamma(d / 2 - 1) * np.eye(d) - 2 * special.gamma(d / 2) * np.outer(Dh, Dh))


# (a) d = 2: J^{11} - J^{22} (traceless, finite) with Delta = (1, 0), y = 0.  Polar coordinates around 0,
#     r outer (split at the singular radius |Delta| = 1), theta inner; r > 4 mapped by r = 4/t.
D2 = np.array([1.0, 0.0])


def inner2(r):
    def g(th):
        z = r * np.array([np.cos(th), np.sin(th)])
        a, b = Kd(-z, 2), Kd(D2 - z, 2)
        return (a[0] * b[0] - a[1] * b[1]) * r
    return integrate.quad(g, -np.pi, np.pi, points=[0.0], limit=400, epsabs=1e-13, epsrel=1e-12)[0]


val = sum(integrate.quad(inner2, a_, b_, limit=400, epsabs=1e-12, epsrel=1e-11)[0]
          for a_, b_ in ((0, 0.5), (0.5, 1.0), (1.0, 1.5), (1.5, 4.0)))
val += integrate.quad(lambda t: inner2(4 / t) * 4 / t**2, 0, 1, limit=400, epsabs=1e-12, epsrel=1e-11)[0]
out['d2_traceless_numeric'] = val
out['d2_traceless_closed'] = -np.pi


# (b) d = 3: Delta = (0, 0, 1.3), y = 0.  J = a delta + b Dhat Dhat; trace = 3a + b, J33 = a + b.
#     Axisymmetric: spherical coordinates around 0, r outer (split at 1.3), cos(theta) inner.
D3 = np.array([0.0, 0.0, 1.3])


def inner3(r, comp):
    def g(ct):
        st = np.sqrt(max(1 - ct * ct, 0.0))
        z = r * np.array([st, 0.0, ct])
        a, b = Kd(-z, 3), Kd(D3 - z, 3)
        return 2 * np.pi * r * r * (a @ b if comp == 'trace' else a[2] * b[2])
    return integrate.quad(g, -1, 1, points=[1.0], limit=400, epsabs=1e-13, epsrel=1e-12)[0]


def J3(comp):
    v = sum(integrate.quad(lambda r: inner3(r, comp), a_, b_, limit=400, epsabs=1e-12, epsrel=1e-11)[0]
            for a_, b_ in ((0, 0.6), (0.6, 1.3), (1.3, 2.0), (2.0, 5.0)))
    return v + integrate.quad(lambda t: inner3(5 / t, comp) * 5 / t**2, 0, 1, limit=400, epsabs=1e-12, epsrel=1e-11)[0]


tr, j33 = J3('trace'), J3('33')
Jc = J_closed(D3, 3)
out['d3_trace_numeric'], out['d3_trace_closed'] = tr, float(np.trace(Jc))
out['d3_J33_numeric'], out['d3_J33_closed'] = j33, float(Jc[2, 2])
out['d3_max_rel'] = max(abs(tr - np.trace(Jc)), abs(j33 - Jc[2, 2])) / abs(np.trace(Jc))

# 3. contraction (sympy)
d, p, k, eps, E, A, B = sp.symbols('d_perp p k epsilon E A B')
P = p + k
c1 = d * p / P**2 - 2 / k          # delta_km delta_k'm'
c2 = 1 / p + p / k**2              # delta_mm' delta_kk'
c3 = 2 / k                         # delta_mk' delta_km'
XX, XD = sp.symbols('XX XD')       # X.X'  and (X.Dhat)(X'.Dhat)
# J = A delta + B Dhat Dhat ; J^{mm} = d A + B
Sigma = (c1 + c3) * (A * XX + B * XD) + c2 * (d * A + B) * XX
Sig = sp.expand(Sigma.subs({A: -sp.pi / 2 * E, B: -sp.pi}))
target = -sp.pi / 2 * (d * (p / P**2 + 1 / p + p / k**2) * E * XX + 2 * (1 / p + p / k**2) * XX + 2 * d * p / P**2 * XD)
out['contraction_minus_result'] = str(sp.simplify(Sig - sp.expand(target)))
# at d = 2 - eps with E = 2/eps + E0 (E0 finite): compare with the earlier (d_perp = 2 - eps) result
E0 = sp.symbols('E0')
red = sp.series(target.subs({d: 2 - eps, E: 2 / eps + E0}), eps, 0, 1).removeO()
earlier = -sp.pi * (XX * ((p / P**2 + p / k**2 + 1 / p) * (2 / eps + E0) - p / P**2) + 2 * p / P**2 * XD)
out['reduction_to_earlier_dperp_2_minus_eps'] = str(sp.simplify(sp.expand(red - earlier)))
# xi form: k * (p/P^2 + 1/p + p/k^2) = Phi(xi), k*(1/p + p/k^2) = xi/(1-xi) + (1-xi)/xi, xi = k/P
xi = sp.symbols('xi', positive=True)
sub = {p: k * (1 - xi) / xi}
out['xi_Phi'] = str(sp.simplify(sp.expand((k * (p / P**2 + 1 / p + p / k**2)).subs(sub)) - (xi * (1 - xi) + xi / (1 - xi) + (1 - xi) / xi)))
out['xi_c2'] = str(sp.simplify((k * (1 / p + p / k**2)).subs(sub) - (xi / (1 - xi) + (1 - xi) / xi)))
out['xi_c1'] = str(sp.simplify((k * p / P**2).subs(sub) - xi * (1 - xi)))

# 4. endpoint decomposition at xi -> 1 with the exact limits (numbers)
kp, Lam, vee, kD = 1.3, 0.07, 9.0, 0.8
x0, x1 = kp / vee, kp / (kp + Lam)
dp_, E_ = 2 - 0.3, 1.7       # any numbers for d_perp and E: the identity is algebraic
g = lambda x: np.exp(-1j * kD / x) / x**2
Phi = lambda x: x * (1 - x) + x / (1 - x) + (1 - x) / x
full = lambda x: g(x) * (dp_ * Phi(x) * E_ + 2 * (x / (1 - x) + (1 - x) / x))
q = lambda f, a, b: integrate.quad(lambda x: f(x).real, a, b, limit=400, epsrel=1e-13)[0] + 1j * integrate.quad(lambda x: f(x).imag, a, b, limit=400, epsrel=1e-13)[0]
lhs = q(full, x0, x1)
reg = lambda x: (dp_ * E_ + 2) * (g(x) - g(1)) / (1 - x) + g(x) * (dp_ * E_ * (x * (1 - x) - 1 + (1 - x) / x) + 2 * (-1 + (1 - x) / x))
rhs = q(reg, x0, x1) + (dp_ * E_ + 2) * g(1) * np.log((vee - kp) * (kp + Lam) / (vee * Lam))
out['endpoint_rel'] = float(abs(lhs - rhs) / abs(lhs))
json.dump(out, open('check_d.json', 'w'), indent=1)
print(json.dumps(out, indent=1))
