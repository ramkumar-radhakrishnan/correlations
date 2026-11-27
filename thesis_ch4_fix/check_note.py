"""Independent check of the formulas exactly as printed in thesis_corrections_note.tex:
normalisation, Schroedinger residual at O(g), Fourier transform, and the thesis-style n=2,3 forms."""
import numpy as np
from math import sqrt, pi
g = 0.02
x = np.linspace(-12, 12, 6001); dx = x[1] - x[0]
P = [lambda x: 1 - g/16*(4*x**4 + 12*x**2 - 9),
     lambda x: sqrt(2)*x*(1 - g/16*(4*x**4 + 20*x**2 - 45)),
     lambda x: (1/sqrt(2))*((2*x**2 - 1) - g/16*(8*x**6 + 52*x**4 - 246*x**2 + 45)),
     lambda x: (1/sqrt(3))*((2*x**3 - 3*x) - g/16*(8*x**7 + 60*x**5 - 510*x**3 + 315*x))]
Q = [lambda k: 1 - g/16*(4*k**4 - 36*k**2 + 15),
     lambda k: sqrt(2)*k*(1 - g/16*(4*k**4 - 60*k**2 + 75)),
     lambda k: (1/sqrt(2))*((2*k**2 - 1) - g/16*(8*k**6 - 172*k**4 + 426*k**2 - 75)),
     lambda k: (1/sqrt(3))*((2*k**3 - 3*k) - g/16*(8*k**7 - 228*k**5 + 930*k**3 - 525*k))]
Ninv2 = [1 + 39/32*g**2, 1 + 315/32*g**2, 1 + 1257/32*g**2, 1 + 3645/32*g**2]
E = [0.5 + 0.75*g, 1.5 + 15/4*g, 2.5 + 39/4*g, 3.5 + 75/4*g]
k = np.linspace(-8, 8, 801)
for n in range(4):
    N = Ninv2[n]**-0.5
    psi = N*pi**-0.25*np.exp(-x**2/2)*P[n](x)
    norm = np.sum(psi**2)*dx
    # residual (H - E) psi should be O(g^2)
    d2 = np.gradient(np.gradient(psi, dx), dx)
    res = -0.5*d2 + (x**2/2 + g*x**4)*psi - E[n]*psi
    m = np.abs(x) < 5
    ft = (psi[None, :]*np.exp(-1j*np.outer(k, x))).sum(1)*dx/np.sqrt(2*pi)
    an = (-1j)**n*N*pi**-0.25*np.exp(-k**2/2)*Q[n](k)
    print('n=%d  norm-1=%.1e  max|(H-E)psi|/g=%.4f (O(g))  |FT-printed|=%.1e' %
          (n, norm-1, np.abs(res[m]).max()/g, np.abs(ft-an).max()))
# residual scaling: should drop ~4x when g halves (O(g^2))
for gg in (0.02, 0.01):
    g = gg
    psi = (1+39/32*g**2)**-0.5*pi**-0.25*np.exp(-x**2/2)*(1 - g/16*(4*x**4 + 12*x**2 - 9))
    d2 = np.gradient(np.gradient(psi, dx), dx)
    res = -0.5*d2 + (x**2/2 + g*x**4)*psi - (0.5+0.75*g)*psi
    print('g=%.3f  max|(H-E)psi| on |x|<5 = %.3e' % (g, np.abs(res[np.abs(x) < 5]).max()))
g = 0.02
# thesis-style n=2 and n=3 (Eqs. n2thesis, n3thesis) equal P_2, P_3
s2 = (1/sqrt(2))*(2*x**2-1) + g/16*(12*sqrt(2) - 14*sqrt(2)*(4*x**4-12*x**2+3) - (1/sqrt(2))*(8*x**6-60*x**4+90*x**2-15))
s3 = (1/sqrt(3))*(2*x**3-3*x) + g/16*(40*sqrt(3)*x - 12*sqrt(3)*(4*x**5-20*x**3+15*x) - (1/sqrt(3))*(8*x**7-84*x**5+210*x**3-105*x))
print('thesis-style n=2 vs P2: %.1e   n=3 vs P3: %.1e' % (np.abs(s2-P[2](x)).max(), np.abs(s3-P[3](x)).max()))
# normalisation closed form for general n
def cn(n, g):
    c = {n: 1.0, n+4: -g/16*np.sqrt((n+1)*(n+2)*(n+3)*(n+4)), n+2: -g/4*(2*n+3)*np.sqrt((n+1)*(n+2))}
    if n >= 2: c[n-2] = g/4*(2*n-1)*np.sqrt(n*(n-1))
    if n >= 4: c[n-4] = g/16*np.sqrt(n*(n-1)*(n-2)*(n-3))
    return c
print('closed-form N^-2 check n=0..9:',
      max(abs(sum(v*v for v in cn(n, g).values()) - (1 + g*g/128*(65*n**4+130*n**3+487*n**2+422*n+156))) for n in range(10)))
print('n=5 coefficients:', {m: round(v/g, 4) for m, v in sorted(cn(5, g).items())},
      'expected c1=%.4f c3=%.4f c7=%.4f c9=%.4f' % (sqrt(120)/16, 9*sqrt(20)/4, -13*sqrt(42)/4, -3*sqrt(21)/4))
# sinc ground state closed forms
from scipy.special import sici
Si2pi = sici(2*pi)[0]
print('2pi/Si(2pi) = %.5f' % (2*pi/Si2pi))
