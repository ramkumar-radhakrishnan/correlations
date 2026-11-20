"""Checks of the Group II Row IV derivation (DGLAP / RC origin)."""
import sympy as sp, numpy as np
from scipy import integrate

# (a) tensor product, line 59:  xi(1-xi) * sum_ij T^{km}_{ij} T^{k'm'}_{ij}
xi, d = sp.symbols('xi d', positive=True)
a, b = 1/(1-xi), 1/xi
def T(k, m, i, j): return sp.KroneckerDelta(k,m)*sp.KroneckerDelta(i,j) - a*sp.KroneckerDelta(j,m)*sp.KroneckerDelta(i,k) - b*sp.KroneckerDelta(i,m)*sp.KroneckerDelta(k,j)
ok = True
for dim in (2, 3, 4):           # check the d-dimensional identity for several integer d
    for k, m, kp, mp in np.ndindex(dim, dim, dim, dim):
        lhs = xi*(1-xi)*sum(T(k,m,i,j)*T(kp,mp,i,j) for i in range(dim) for j in range(dim))
        D = sp.KroneckerDelta
        rhs = dim*xi*(1-xi)*D(k,m)*D(kp,mp) + 2*(D(m,kp)*D(k,mp)-D(k,m)*D(kp,mp)) + (xi/(1-xi)+(1-xi)/xi)*D(m,mp)*D(k,kp)
        if sp.simplify(lhs - rhs) != 0: ok = False
print('(a) tensor contraction, line 59, d=2,3,4:', 'OK' if ok else 'MISMATCH')

# (b) traceless (finite) part of the r-integral, line 77:
#     int d^2r [r^m (r+u)^m' - (1/2) delta^{mm'} r.(r+u)] / (r^2 (r+u)^2) = -pi (u^m u^m'/u^2 - delta/2)
u = np.array([0.7, -0.4])
def integrand(rho, phi, m, mp):
    r = rho*np.array([np.cos(phi), np.sin(phi)]); s = r + u
    val = r[m]*s[mp] - 0.5*(m == mp)*(r @ s)
    return rho*val/((r @ r)*(s @ s))
res = {}
for m, mp in [(0,0), (0,1), (1,1)]:
    v = 0
    for lo, hi in [(0, 0.3), (0.3, 3), (3, 30), (30, 3000)]:
        v += integrate.dblquad(lambda phi, rho: integrand(rho, phi, m, mp), lo, hi, 0, 2*np.pi, epsabs=1e-10, epsrel=1e-9)[0]
    res[(m, mp)] = v
uu = u/np.linalg.norm(u)
for (m, mp), v in res.items():
    print(f'(b) finite part m={m},m\'={mp}: numeric {v:+.5f}   -pi(u^m u^m\'/u^2 - delta/2) = {-np.pi*(uu[m]*uu[mp]-0.5*(m==mp)):+.5f}')

# (c) trace of the dim-reg r-integral in d=2-eps: -(pi/2) d L - pi  ->  -pi L  + O(eps)  (eps*L -> 2)
eps, L = sp.symbols('epsilon L')
tr = -(sp.pi/2)*(2-eps)*(2/eps + L) - sp.pi
print('(c) trace of r-integral:', sp.simplify(sp.series(tr, eps, 0, 1).removeO()), ' (expected -pi*(2/eps + L) with L the finite log)')

# (d) momentum sum rule of P_gg (pure glue)
x = sp.symbols('x', positive=True)
Nc = sp.symbols('N_c', positive=True)
reg = 2*Nc*sp.integrate(x*(x*(1-x) + (1-x)/x), (x, 0, 1)) + 2*Nc*sp.integrate((x*x - 1)/(1-x), (x, 0, 1))
print('(d) int xi P_gg = ', sp.simplify(reg + sp.Rational(11,6)*Nc), ' (beta0 = 11 Nc/3)')
