"""The MS-bar answer written entirely in transverse position space.

  T   = int d^2r e^{iq.r}/r^2        needs a transverse UV regulator
  G   = int d^2r e^{iq.r}/(r^2 D)    D = (r-s)^2 + M^2 ,  D0 = D|_{r=0} = s^2/xibar

three equivalent regulators, and the dictionary between them.
"""
import numpy as np, sys
sys.path.insert(0, '/home/user/correlations/zint/src')
from scipy.integrate import quad
from scipy.special import k0, j0
from uv_master import Ihat
g = 0.5772156649015329
s = np.array([0.7, -0.4]); k = np.array([0.9, 0.6]); s2 = float(s@s)

def T_disc(qn, r0):                      # int_{|r|>r0} d^2r e^{iq.r}/r^2  = 2pi int_r0^inf dr J0(qr)/r
    f = lambda r: j0(qn*r)/r
    v, pts = 0.0, [r0] + list(np.arange(1, 400)*np.pi/qn + r0)
    part = []
    for a, b in zip(pts[:-1], pts[1:]):
        v += quad(f, a, b, limit=200)[0]; part.append(v)
    part = np.array(part)                 # average the oscillating partial sums
    for _ in range(6): part = 0.5*(part[:-1] + part[1:])
    return 2*np.pi*part[-1]

print("1. THREE REGULATORS FOR  T = int d^2r e^{iq.r}/r^2   AGREE, WITH r0 = m")
print("   disc  |r|>r0 :  pi log(4 e^-2gamma /(q^2 r0^2))")
print("   mass  1/(r^2+m^2) :  2 pi K0(m|q|)  ->  pi log(4 e^-2gamma /(q^2 m^2))")
print(f"   {'q':>6} {'r0 = m':>9} {'disc (quad)':>16} {'2pi K0(m q)':>16} {'pi log(4e^-2g/q^2r0^2)':>24}")
for qn in (0.6, 1.4):
    for r0 in (1e-2, 1e-3, 1e-4):
        print(f"   {qn:6.2f} {r0:9.0e} {T_disc(qn, r0):16.9f} {2*np.pi*k0(r0*qn):16.9f}"
              f" {np.pi*np.log(4*np.exp(-2*g)/(qn**2*r0**2)):24.9f}")
print("   -> the sharp disc cut and the mass regulator are the SAME to O(r0^2), with r0 = m.")
print()

print("2. MS-bar DICTIONARY.   pi log(mu^2/q^2)  ==  pi log(4 e^-2gamma/(q^2 r0^2))")
print("   =>  r0 = 2 e^-gamma / mu      (mu here is the MS-bar scale)")
for mu2 in (1.3, 4.0):
    mu = np.sqrt(mu2); r0 = 2*np.exp(-g)/mu
    for qn in (0.6, 1.4):
        print(f"   mu^2={mu2:<5} q={qn:<5} r0={r0:.6f}:  pi log(mu^2/q^2) = {np.pi*np.log(mu2/qn**2):+.9f}"
              f"   pi log(4e^-2g/(q^2 r0^2)) = {np.pi*np.log(4*np.exp(-2*g)/(qn**2*r0**2)):+.9f}")
print()

print("3. THE FULL MASTER IN POSITION SPACE:   G = (pi/D0) [ log(D0/r0^2) + Ihat(xi) ]")
print("   (Ihat carries no regulator and no mu; it is the same object as before.)")
print("   direct 2D quadrature of int_{|r|>r0} d^2r e^{iq.r}/(r^2 D)  vs  the formula:")
def G_disc(s, k, xi, r0, nth=4096):
    xb = 1-xi; M2 = xi/xb*s2; q = xb*k; 
    th = np.linspace(0, 2*np.pi, nth, endpoint=False)
    u = np.stack([np.cos(th), np.sin(th)], -1)
    def rad(rho):
        r = u*rho
        val = np.exp(1j*(r@q))/(rho**2*(np.sum((r-s)**2, -1)+M2))
        return rho*np.mean(val)*2*np.pi
    pts = [r0] + list(np.geomspace(max(r0*2, 1e-6), 2e4, 220))
    tot = 0
    for a, b in zip(pts[:-1], pts[1:]):
        if b <= a: continue
        tot += (quad(lambda x: rad(x).real, a, b, limit=200)[0]
                + 1j*quad(lambda x: rad(x).imag, a, b, limit=200)[0])
    return tot
for xi in (0.3, 0.65):
    xb = 1-xi; D0 = s2/xb
    for r0 in (1e-2, 1e-3):
        num = G_disc(s, k, xi, r0)
        ana = np.pi/D0*(np.log(D0/r0**2) + Ihat(s, k, xi))
        print(f"   xi={xi:<5} r0={r0:.0e}:  quad {num:+.9f}   formula {ana:+.9f}   rel {abs(num-ana)/abs(ana):.1e}")
print()

print("4. TRADING r0 FOR THE MS-bar mu REPRODUCES THE PUBLISHED FORM")
print("   log(D0/r0^2)  with  r0^2 = 4 e^-2gamma/mu^2  and  D0 = s^2/xibar :")
print("      = log( mu^2 s^2 / xibar ) + 2(gamma - log 2)")
for mu2 in (1.3, 4.0):
    for xi in (0.3, 0.65):
        xb = 1-xi; D0 = s2/xb; r02 = 4*np.exp(-2*g)/mu2
        print(f"   mu^2={mu2:<5} xi={xi:<5}:  log(D0/r0^2) = {np.log(D0/r02):+.12f}"
              f"   log(mu^2 s^2/xibar)+2(g-log2) = {np.log(mu2*s2/xb)+2*(g-np.log(2)):+.12f}")
