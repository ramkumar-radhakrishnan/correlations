"""Group II Row IV: check the plus prescription done directly in p+ (five-term form)."""
import mpmath as mp
mp.mp.dps=25
def check(K,Lam,Vee,kap,dperp,XX,RR,epsinv,Lr):
    M=epsinv+1-Lr; T=Vee-K
    S=lambda p:K+p
    ph=lambda p: mp.e**(-1j*kap*S(p)/K)          # e^{-i kappa} e^{-i kappa p+/k+}
    ph0=mp.e**(-1j*kap)
    # --- exact starting point (verified three-block form, limits Lambda .. T) ---
    I=lambda p: ph(p)*( (p/S(p)**2+p/K**2+1/p)*M*XX - dperp*p/S(p)**2*RR - (p/K**2+1/p)*XX )
    exact=mp.quad(I,[Lam,K,T])
    # --- the five terms as written ---
    t1 = mp.quad(lambda p: ph(p)*(p/S(p)**2+p/K**2)*M*XX,[0,K,T]) \
       + mp.quad(lambda p: (ph(p)-ph0)/p*M*XX,[0,K,T])                       # [1/p+]_+
    t2 = ph0*mp.log(T/Lam)*M*XX
    t3 = -mp.quad(lambda p: ph(p)*dperp*p/S(p)**2*RR,[Lam,K,T])
    t4 = -mp.quad(lambda p: ph(p)*p/K**2*XX,[0,K,T]) \
         -mp.quad(lambda p: (ph(p)-ph0)/p*XX,[0,K,T])
    t5 = -ph0*mp.log(T/Lam)*XX
    return exact, t1+t2+t3+t4+t5, t2+t5
print("="*78); print("THE FIVE-TERM p+ FORM vs THE EXACT STARTING POINT"); print("="*78)
for a in [(1.0,1e-5,60.0, 0.9,2.0, 1.7,-0.6, 12.0, 2.3),
          (2.5,2e-6,150.0,-1.4,2.0,-0.8, 1.1, -7.0,-1.1),
          (0.7,1e-5,40.0,  2.2,2.0, 2.4, 0.9,  5.0, 0.4)]:
    ex,su,ev = check(*a)
    print("   k+=%.1f kappa=%+.1f Lam=%.0e"%(a[0],a[3],a[1]))
    print("      exact      %s"%mp.nstr(ex,12))
    print("      five terms %s   diff %s  (O(Lambda), from extending 0..T)"%(mp.nstr(su,12),mp.nstr(abs(ex-su),3)))
    print("      terms 2+5  %s   = evolution sector"%mp.nstr(ev,10))
print()
print("="*78); print("WHICH LOG?  p+-prescription vs xi-prescription"); print("="*78)
K,Lam,Vee=1.0,1e-5,60.0
print("   p+ version :  delta(p+) log[(V-k+)/Lambda] = log(%.4e)"%((Vee-K)/Lam))
print("   xi version :  delta(1-xi) log[(1-xi0)/(1-xi1)] = log(%.4e)"%((1-K/Vee)*(K+Lam)/Lam))
print("   they differ by log[(V-k+)/k+] = %.4f -- a FINITE reshuffling between the"%mp.log((Vee-K)/K))
print("   delta term and the [ ]_+ remainder, because  dxibar/xibar = xi dp+/p+ .")
print("   The SUM is identical; only the split is convention.")
