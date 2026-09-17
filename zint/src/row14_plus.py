"""Row 14: what the 1/p+ master really contains, and what that means for [1/p+]_+."""
import numpy as np, mpmath as mp
mp.mp.dps=25
K=1.0
print("="*78); print("THE 1/p+ MASTER"); print("="*78)
A,B=0.62,1.45; be=K*B/A
print("   M[1/p] = (1/(A beta)) [ Ei(-i kappa p) - e^{i kappa beta} Ei(-i kappa(p+beta)) ]")
print("   as Lambda -> 0 :  Ei(-i kappa Lambda) -> gamma_E + log(|kappa| Lambda) - i pi sgn(kappa)/2")
print()
print("   %-10s %-10s %22s %22s"%("kappa","Lambda","Ei(-i kappa Lambda)","gamma+log(|k|Lam)-i pi/2"))
for kap in (0.83, 2.4):
    for Lm in (1e-6,1e-8):
        e=complex(mp.ei(-1j*kap*Lm))
        a=np.euler_gamma+np.log(abs(kap)*Lm)-1j*np.pi/2*np.sign(kap)
        print("   %-10.2f %-10.0e %22s %22s"%(kap,Lm,"%.6f%+.6fj"%(e.real,e.imag),
                                              "%.6f%+.6fj"%(a.real,a.imag)))
print()
print("   => the rapidity logarithm comes out as  log(|kappa| Lambda),  NOT  log(Lambda) alone.")
print("      kappa = k.(y'-z)/k+ , so the log carries a TRANSVERSE dependence:")
print("         log(|kappa| Lambda) = log(Lambda) + log( |k.(y'-z)| / k+ )")
print()
print("="*78); print("WHY THAT MATTERS FOR THE + PRESCRIPTION"); print("="*78)
print("   At FIXED transverse points the identity is exact and harmless:")
print("      g(p+) = e^{-i kappa p+}/(p+ A + k+ B) ,  g(0) = 1/(k+ B)")
print("      int_Lam^P dp/p g = g(0) log(P/Lam) + int_0^P dp/p [g - g(0)]   (exact)")
print()
print("   The trouble is the z integral that follows.  At large |z|:")
print("      A = (x-z)^2 ~ |z|^2  so  g(p+) ~ 1/(p+ |z|^2) -> 0 ,  but  g(0) = 1/(k+B) is")
print("      z-INDEPENDENT.  And kappa ~ |z|, so the phase damps the full integrand only for")
print("      p+ |kappa| |z| >> 1, i.e. only out to  |z|_max ~ k+/(|k| p+).")
print("   The damping scale DEPENDS ON p+ : the z and p+ integrations do not commute.")
print()
# demonstrate the damping scale
kt=np.array([0.7,-0.4]); kn=np.sqrt(kt@kt)
print("   |z|_max ~ k+/(|k| p+):")
for P in (1e-1,1e-2,1e-3,1e-4):
    print("      p+ = %-8.0e  ->  |z|_max ~ %10.1f"%(P,K/(kn*P)))
print()
print("   consequence: the log(Lambda) multiplies a z-integral that is itself log divergent,")
print("   and log|kappa| ~ log|z| against d^2z/|z|^2 gives a SECOND logarithm.  So after the")
print("   z integration the row carries  log(1/Lambda) x log(R)  and  log^2(R)  --  equivalently")
print("   a MIXED  (1/eps) log(1/Lambda)  once z is regulated in d = 2-2eps.")
