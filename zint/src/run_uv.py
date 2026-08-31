"""Every number in the UV-only note. -> ../RESULTS_uv.txt"""
import numpy as np, sympy as sp, mpmath as mp
out=[]
def p(*a): out.append(" ".join(str(x) for x in a)); print(*a)

s=np.array([0.7,-0.4]); P=np.array([1.3,0.5]); s2=s@s; P2=P@P; sP=s@P
p("kinematics: s = %s  P = %s   s^2 = %.4f  s.P = %.4f"%(s,P,s2,sP))
p("")
p("A. the angular average, exactly")
th=2*np.pi*np.arange(200000)/200000
R=1e-3; rv=np.stack([R*np.cos(th),R*np.sin(th)],-1)
p("   < (P.r)(s.r)/r^4 > r^2/(s.P) = %.12f     1/d = 1/2 = 0.5"
  %(((rv@P)*(rv@s)/np.sum(rv*rv,-1)**2).mean()*R**2/sP))
p("   < (P.r)/r^2 >     r  /|P|    = %.3e      (odd, vanishes)"
  %(((rv@P)/np.sum(rv*rv,-1)).mean()*R/np.hypot(*P)))
p("")
p("B. UV coefficient of each Theta_i, from the raw tensors of Eq. (1.1)")
def parts(r, xi, dperp=2.0):
    xb=1-xi; r2=np.sum(r*r,-1); u=s-r; u2=np.sum(u*u,-1)
    Pr=r@P; sr=r@s; Pu=u@P; su=u@s; ur=np.sum(u*r,-1); dd=Pr/(P2*r2)
    return [dd*dperp/(2*r2)*(s2-u2)*xi*xb,
            -(xi/xb)*( sP*(-r2)/r2 - sP*ur/(2*u2) )/(P2*r2),
            xb*dd*( np.sum(u*(-r),-1)/r2 + su/(2*s2) - (s2-u2)/(2*r2) ),
            -(xb/xi)*( Pu*(-r2)/r2 + sr*Pu/(2*s2) )/(P2*r2),
            xi*dd*( -sr/r2 - su/(2*u2) - (s2-u2)/(2*r2) ),
            -( sr*(-Pr)/r2 - Pu*sr/(2*u2) + ur*(-Pr)/r2 + sP*ur/(2*s2) )/(P2*r2)]
th=2*np.pi*np.arange(40000)/40000
for xi in (0.3,0.65):
    xb=1-xi
    for R in (1e-4,):
        rv=np.stack([R*np.cos(th),R*np.sin(th)],-1)
        c=[v.mean()*R**2/(sP/P2) for v in parts(rv,xi)]
        w=[xi*xb, xi/xb, xb, xb/xi, xi, 1.0]
        p("   xi=%.2f : Theta_i UV coefficients %s"%(xi,np.round(c,7)))
        p("            = %s  (weights)  ->  sum %.9f   C_UV = %.9f"
          %(np.round(w,5), sum(c), xi*xb+xi/xb+xb/xi))
p("")
p("C. the cancellation among Theta_3, Theta_5, Theta_6 holds in any d")
d,dp,x=sp.symbols('d d_perp xi',positive=True); xb=1-x
tot = dp*x*xb/d + x/xb + xb/x + (-2/d)*xb + (-2/d)*x + (2/d)*1
p("   sum = %s"%sp.simplify(tot))
p("   at d = d_perp : %s"%sp.simplify(tot.subs(d,dp)))
p("   Pgg/(2Nc)     : %s"%sp.simplify(x*xb+x/xb+xb/x))
p("   identical?    : %s"%(sp.simplify(tot.subs(d,dp)-(x*xb+x/xb+xb/x))==0))
p("")
p("D. the log coefficient measured on the FULL 2D integral (independent of all the algebra)")
p("   d W / d log(1/m^2) of  int d^2r e^{iq.r} bracket/D , r^2 -> r^2+m^2 :")
p("      xi=0.30  measured +7.137881   predicted C_UV (s.P) pi/D_0 = +7.138842   rel 1.3e-04")
p("      xi=0.65  measured +3.150126   predicted C_UV (s.P) pi/D_0 = +3.150493   rel 1.2e-04")
p("   (from RESULTS.txt, block D; regulator r^2 -> r^2+m^2 applied to the original integrand)")
p("")
p("E. the p+ integrals of C_UV = Pgg/(2Nc)   (mpmath, 30 digits, lambda = 1e-9)")
mp.mp.dps=30; lam=mp.mpf('1e-9'); lk=mp.log(1/lam)
C=lambda t: t*(1-t)+t/(1-t)+(1-t)/t
I0=mp.quad(lambda t:C(t),[lam,1-lam])
I2=mp.quad(lambda t:C(t)*mp.log(1/(1-t)),[lam,1-lam])
I1=mp.quad(lambda t:C(t)*mp.log(t*(1-t)),[lam,1-lam])
p("   int C_UV                 = %s    2 l_k - 11/6            = %s"%(mp.nstr(I0,14),mp.nstr(2*lk-mp.mpf(11)/6,14)))
p("   int C_UV log(1/(1-xi))   = %s    l_k^2/2 + pi^2/6 - 67/36 = %s"%(mp.nstr(I2,14),mp.nstr(lk**2/2+mp.pi**2/6-mp.mpf(67)/36,14)))
p("   int C_UV log(xi(1-xi))   = %s   -l_k^2 + 67/18 - pi^2/3   = %s   [= paper (C.10)]"%(mp.nstr(I1,14),mp.nstr(-lk**2+mp.mpf(67)/18-mp.pi**2/3,14)))
p("   the three pieces separately:")
for nm,f,pred in [("xi(1-xi)  [Theta_1]",lambda t:t*(1-t),mp.mpf(1)/6),
                  ("xi/(1-xi) [Theta_2]",lambda t:t/(1-t),lk-1),
                  ("(1-xi)/xi [Theta_4]",lambda t:(1-t)/t,lk-1)]:
    p("      int %s = %s   predicted %s"%(nm,mp.nstr(mp.quad(f,[lam,1-lam]),12),mp.nstr(pred,12)))
for nm,f,pred in [("xi(1-xi)  [Theta_1]",lambda t:t*(1-t)*mp.log(1/(1-t)),mp.mpf(5)/36),
                  ("xi/(1-xi) [Theta_2]",lambda t:t/(1-t)*mp.log(1/(1-t)),lk**2/2-1),
                  ("(1-xi)/xi [Theta_4]",lambda t:(1-t)/t*mp.log(1/(1-t)),mp.pi**2/6-1)]:
    p("      int %s log(1/(1-xi)) = %s   predicted %s"%(nm,mp.nstr(mp.quad(f,[lam,1-lam]),12),mp.nstr(pred,12)))
open('/home/user/correlations/zint/RESULTS_uv.txt','w').write("\n".join(out)+"\n")
