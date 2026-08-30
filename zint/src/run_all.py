"""Every numerical check quoted in the note.  Writes RESULTS.txt."""
import numpy as np, sys, io, contextlib
sys.path.insert(0, '/home/user/correlations/zint/src')
from quad2d import integrate
from masters import Kin, G_scalar, Psi, Nten, Tten, Phi, A7, A8, A9
from assemble import W_parts, T_fin, coeffs
from scipy.special import k0
import mpmath as mp

s = np.array([0.7, -0.4]); q = np.array([0.9, 0.6]); P = np.array([1.3, 0.5])
xi = 0.3; xb = 1-xi; s2 = float(s@s); M2 = xi/xb*s2; D0 = s2+M2; sp_ = np.hypot(*s)
K = Kin(s, q, M2)
Dq = lambda r: np.sum((r-s)**2, -1)+M2
r2 = lambda r: np.sum(r**2, -1); u2 = lambda r: np.sum((r-s)**2, -1)
out = []
def p(*a): out.append(" ".join(str(x) for x in a)); print(*a)

p("kinematics:  s =", s, " q =", q, " P =", P, " xi =", xi)
p("             s^2 = %.4f  M^2 = %.6f  D0 = %.6f  s.P = %.4f  q.s = %.4f"
  % (s2, M2, D0, s@P, q@s))
p("")
p("A. master integrals:  alpha-representation  vs  direct 2D quadrature")
p("%-10s %-42s %-42s %s" % ("", "alpha-representation", "2D quadrature", "rel.diff"))
def row(name, a, b):
    p("%-10s %-42s %-42s %.1e" % (name, f"{a:+.10f}", f"{b:+.10f}", abs(a-b)/abs(b)))
ps, N, T, ph = Psi(K), Nten(K), Tten(K), Phi(K)
for i in range(2):
    row(f"Psi^{i}", ps[i], integrate(lambda r, i=i: r[..., i]/(r2(r)*Dq(r)), q, split=(sp_,)))
for (i, j) in [(0,0),(0,1),(1,1)]:
    row(f"N^{i}{j}", N[i,j], integrate(lambda r,i=i,j=j: r[...,i]*r[...,j]/(r2(r)*Dq(r)), q, split=(sp_,)))
for (i, j) in [(0,0),(0,1)]:
    row(f"Tt^{i}{j}", T[i,j], integrate(lambda r,i=i,j=j:
        (r[...,i]*r[...,j]/r2(r)**2-(1.0 if i==j else 0.0)/(2*r2(r)))/Dq(r), q, split=(sp_,)))
row("A7", A7(K), integrate(lambda r: np.sum(r*(r-s),-1)/(r2(r)*u2(r)), q, split=(sp_,), ntheta=2048))
row("A8", A8(K,P), integrate(lambda r: (r@P)*np.sum(s*(r-s),-1)/(r2(r)*u2(r)), q, split=(sp_,), ntheta=2048))
row("A9", A9(K,P), integrate(lambda r: (r@s)*np.sum(P*(r-s),-1)/(r2(r)*u2(r)), q,
                             split=(sp_,), ntheta=4096, nlev=20, R0=12.0))
p("")
p("B. exact identities")
p("   trace N^ij - Phi                      : %.2e" % abs(N[0,0]+N[1,1]-ph))
p("   trace Tt^ij                           : %.2e" % abs(T[0,0]+T[1,1]))
p("   G(m) alpha-rep vs (2pi/D0)K0(mq)+(2s.Psi-Phi)/D0, m=1e-3 : %.2e"
  % (abs(G_scalar(Kin(s,q,M2,1e-6))-((2*np.pi/D0)*k0(1e-3*np.hypot(*q))+(2*(s@ps)-ph)/D0))
     / abs(G_scalar(Kin(s,q,M2,1e-6)))))
p("   A8(P=s) - A9(P=s)                     : %.2e" % abs(A8(K,s)-A9(K,s)))
p("")
p("C. the regular block  W_reg  (assembled from masters)  vs direct quadrature of T_fin/D")
for x in (0.3, 0.65):
    R = W_parts(s, q, P, x, 1e-3)
    M2x = x/(1-x)*s2
    num = integrate(lambda r: T_fin(r,s,P,x)/(np.sum((r-s)**2,-1)+M2x), q,
                    split=(sp_,), ntheta=4096, nlev=20, R0=12.0)
    p("   xi=%.2f  assembled %-40s quad %-40s rel %.1e"
      % (x, f"{R['reg']:+.8f}", f"{num:+.8f}", abs(R['reg']-num)/abs(num)))
p("")
p("D. the UV coefficient: d W / d log(1/m^2) of the *regulated original* integrand")
for x in (0.3, 0.65):
    c = coeffs(x); M2x = x/(1-x)*s2; D0x = s2+M2x
    def Tfull(r, m2, x=x, c=c):
        rr = np.sum(r**2,-1)+m2; u = r-s; uu = np.sum(u**2,-1)
        Pr, sr, sP = r@P, r@s, s@P
        return (2*x*(1-x)*Pr*sr/rr**2 + (x/(1-x)+(1-x)/x)*sP/rr + c['c3']*Pr/rr
                + c['c4']*Pr*sr/(s2*rr) + c['c5']*sP*sr/(s2*rr) + c['c6']*sP/s2
                + c['c7']*sP*np.sum(r*u,-1)/(rr*uu) + c['c8']*Pr*np.sum(s*u,-1)/(rr*uu)
                + c['c9']*sr*np.sum(P*u,-1)/(rr*uu))
    v = [integrate(lambda r, m=m: Tfull(r,m*m)/(np.sum((r-s)**2,-1)+M2x), q,
                   split=(sp_,), ntheta=4096, nlev=20, R0=12.0) for m in (1e-2, 1e-3)]
    slope = (v[1]-v[0])/np.log(1e2)   # d log(1/m^2) between m=1e-2 and m=1e-3
    p("   xi=%.2f  measured %+.6f   predicted C_UV (s.P) pi/D0 = %+.6f   rel %.1e"
      % (x, slope.real, c['CUV']*(s@P)*np.pi/D0x, abs(slope.real-c['CUV']*(s@P)*np.pi/D0x)/abs(c['CUV']*(s@P)*np.pi/D0x)))
p("")
p("E. xi -> 0 : the ln(1/M^2) of the companion cancels that of the regular block")
p("   %-8s %-26s %-26s %s" % ("xi", "xi*W_companion", "xi*W_reg", "xi*(W_comp+W_reg)"))
for x in (1e-2, 1e-3, 1e-4, 1e-5, 1e-6):
    R = W_parts(s, q, P, x, 1e-3)
    p("   %-8.0e %-26s %-26s %s" % (x, f"{x*R['companion']:+.4f}", f"{x*R['reg']:+.4f}",
                                    f"{x*(R['companion']+R['reg']):+.6f}"))
pred = -np.pi*np.exp(1j*(q@s))*(s@P)/s2
p("   coefficient of ln(1/M^2)/xi in W_reg: predicted %s" % f"{pred:+.6f}")
xs = [1e-4, 3e-5, 1e-5]; vv = [W_parts(s,q,P,x,1e-3)['reg'] for x in xs]
a = (xs[2]*vv[2]-xs[1]*vv[1])/(np.log(1/xs[2])-np.log(1/xs[1]))
p("                                        measured  %s   ratio %.5f" % (f"{a:+.6f}", (a/pred).real))
p("")
p("F. xi -> 1 at fixed q :  W_reg -> -(s.P) A7 /(2 s^2)")
p("   prediction %s" % f"{-(s@P)*A7(Kin(s,q,1.0))/(2*s2):+.7f}")
for xbv in (1e-3, 1e-4, 1e-5):
    p("   xibar=%-8.0e W_reg = %s" % (xbv, f"{W_parts(s,q,P,1-xbv,1e-3)['reg']:+.7f}"))
p("")
p("G. the p+ integrals of P_gg/(2Nc)   (mpmath, 30 digits, lambda=1e-8)")
mp.mp.dps = 30; lam = mp.mpf('1e-8'); lk = mp.log(1/lam)
I = [mp.quad(lambda t:(1-t)/t,[lam,1-lam]), mp.quad(lambda t:t/(1-t),[lam,1-lam]),
     mp.quad(lambda t:t*(1-t),[lam,1-lam])]
J = [mp.quad(lambda t:(1-t)/t*mp.log(1/(1-t)),[lam,1-lam]),
     mp.quad(lambda t:t/(1-t)*mp.log(1/(1-t)),[lam,1-lam]),
     mp.quad(lambda t:t*(1-t)*mp.log(1/(1-t)),[lam,1-lam])]
p("   sum of plain integrals   %s   2 l_k - 11/6      = %s" % (mp.nstr(sum(I),13), mp.nstr(2*lk-mp.mpf(11)/6,13)))
p("   sum of log-weighted ones %s   l_k^2/2+pi^2/6-67/36 = %s" % (mp.nstr(sum(J),13), mp.nstr(lk**2/2+mp.pi**2/6-mp.mpf(67)/36,13)))
open('/home/user/correlations/zint/RESULTS.txt','w').write("\n".join(out)+"\n")
