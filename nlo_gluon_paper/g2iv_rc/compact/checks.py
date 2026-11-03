"""Checks for G2IV_compact_RC.tex.   python3 checks.py"""
import sympy as sp

g, Nc, kp, V, k, xi, eps, aS, b0, L, r, mu, muR, Q = sp.symbols(
    'g N_c kplus vee k xi epsilon alpha_s beta_0 L r mu mu_R Q', positive=True)
pi = sp.pi

# (1) g^2/8pi^2 = alpha_s/2pi ;  g^2/4pi^2 N_c = (alpha_s/2pi) 2N_c
assert sp.simplify(g**2 / (8 * pi**2) - (g**2 / (4 * pi)) / (2 * pi)) == 0

# (2) exact bracket with d_perp = 2 - eps
dp = 2 - eps
B = -dp / eps - dp * sp.EulerGamma - dp / 2 * L - 1          # L = log(Delta^2 mu^2/4)
assert sp.simplify(B - (-2 / eps - dp / 2 * (2 * sp.EulerGamma + L))) == 0
print('bracket = -2/eps - (d/2) l_Delta : ok')

# (3) xi integral of line 2 after the y' Fourier transform (xi^2 from the FT)
I = sp.integrate((1 - xi) / xi**3 * xi**2, (xi, kp / V, 1))
assert sp.simplify(I - (sp.log(V / kp) - (1 - kp / V))) == 0
print('int (1-xi)/xi = log(vee/k+) - (1-k+/vee) : ok')

# (4) normalisation: line 2 coefficient vs H_4 coefficient (common 1/(2pi)^3 1/k+ int_{x,x',z} dropped)
line2 = (g**2 / (4 * pi**2)) * Nc * g**2 / (2 * pi**2) * (-4 * pi / k**2)     # FT of log r^2 = -4pi/k^2 (xi^2 taken out)
H4 = g**4 * Nc / (8 * pi**5) * (4 * pi**2 / k**2)                          # |FT of r^i/r^2|^2 = 4pi^2/k^2
assert sp.simplify(line2 + H4) == 0
print('line 2 = -[log(vee/k+) - (1-k+/vee)] H_4 : ok')

# (5) running coupling: (alpha_s b0/4pi)[2/eps + l(mu)] - (alpha_s b0/4pi)(2/eps) = (alpha_s b0/4pi) log(r^2 muR^2/c0^2) at mu=muR
c0 = 2 * sp.exp(-sp.EulerGamma)
lD = 2 * sp.EulerGamma + sp.log(r**2 * muR**2 / 4)
assert sp.simplify(sp.expand_log(lD - sp.log(r**2 * muR**2 / c0**2), force=True)) == 0
#   alpha_s(Q) = alpha_s(muR)[1 - alpha_s b0/(4pi) log(Q^2/muR^2)] ; Q = c0/r
aQ = aS * (1 - aS * b0 / (4 * pi) * sp.log(Q**2 / muR**2))
rc = aS * (1 + aS * b0 / (4 * pi) * sp.log(r**2 * muR**2 / c0**2))
assert sp.simplify(sp.expand_log(aQ.subs(Q, c0 / r) - rc, force=True)) == 0
#   mu_R independence at O(alpha_s^2): d/dlog muR^2 [alpha_s(muR)] = -b0 alpha_s^2/4pi
t = sp.symbols('t')   # t = log muR^2
dLO = -b0 * aS**2 / (4 * pi)
drc = aS**2 * b0 / (4 * pi) * sp.diff(sp.log(r**2) + t - sp.log(c0**2), t)
assert sp.simplify(dLO + drc) == 0
print('running coupling alpha_s(c0/|y\'-y|), mu_R independent : ok')
