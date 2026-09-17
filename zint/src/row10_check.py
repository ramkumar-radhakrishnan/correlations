"""Row 10: A^(1)dagger x A^(3).  Check the N1/N2 split."""
import sympy as sp, numpy as np
g,pi = sp.symbols('g pi',positive=True); K,V,Lam,Nc = sp.symbols('kplus V Lambda N_c',positive=True)
I = sp.I

print("="*76); print("1.  PREFACTORS"); print("="*76)
A1  = I*g/(sp.sqrt(2)*pi*sp.sqrt(K))                      # A^(1)  from  -A = -(1/sqrt2 pi)(ig/sqrt k+)
A1d = -I*g/(sp.sqrt(2)*pi*sp.sqrt(K))                     # its conjugate
A3a = -I*g**3*Nc/(16*sp.sqrt(2)*pi**3*sp.sqrt(K))         # A^(3) one-rho piece
A3b =  I*g**3/(32*sp.sqrt(2)*sp.sqrt(K)*pi**5)            # A^(3) two-rho piece
over = 1/(2*pi)**3
print("   bra bracket = -Abar^dag + A^dag U  ->  +(ig/(sqrt2 pi sqrt k+)) [U(x') - U(w')]")
n1 = sp.simplify(over*(-A1d)*A3a)
print("   N1 = (1/(2pi)^3) x (-A^dag) x A3a =", sp.simplify(n1))
print("        quoted  (1/(2pi)^3) g^4 Nc/(32 pi^4 k+)  ->  ratio",
      sp.simplify(n1/(over*g**4*Nc/(32*pi**4*K))))
n2 = sp.simplify(over*(-A1d)*A3b)
print("   N2 = (1/(2pi)^3) x (-A^dag) x A3b =", sp.simplify(n2))
print("        quoted -(1/(2pi)^3) g^4/(64 pi^6 k+)     ->  ratio",
      sp.simplify(n2/(-over*g**4/(64*pi**6*K))))

print(); print("="*76); print("2.  THE LOG BRACKET IN N1 COLLAPSES"); print("="*76)
expr = sp.log(Lam/K) + sp.log(V/K - 1) - sp.log(V/Lam)
simp = 2*sp.log(Lam/K) + sp.log(1 - K/V)
print("   log(Lam/k+) + log(V/k+ - 1) - log(V/Lam)  ==  2 log(Lam/k+) + log(1 - k+/V) ?")
print("      symbolic difference :", sp.simplify(sp.expand_log(expr-simp, force=True)))
for kk,vv,ll in [(1.0,50.0,1e-4),(2.3,17.0,1e-6),(0.4,900.0,1e-3)]:
    a=np.log(ll/kk)+np.log(vv/kk-1)-np.log(vv/ll); b=2*np.log(ll/kk)+np.log(1-kk/vv)
    print("      k+=%-5.1f V=%-6.1f Lam=%.0e :  %14.9f  vs %14.9f   diff %.1e"%(kk,vv,ll,a,b,abs(a-b)))
print()
print("   so the second factor is   11/3 - 2 log(k+/Lambda) + log(1 - k+/V)")

print(); print("="*76); print("3.  THE V -> infinity LIMIT OF N1"); print("="*76)
tail = sp.log(1-K/V)**2 + sp.log(1-K/V)*sp.log(V/K)
print("   log^2(1-k+/V) + log(1-k+/V) log(V/k+)  -->", sp.limit(tail, V, sp.oo))
print("   log(1-k+/V)                            -->", sp.limit(sp.log(1-K/V), V, sp.oo))
print("   => ALL the V dependence of N1 vanishes as V -> infinity, leaving")
print()
print("      [ -2/eps - 2 gamma - log((x-w)^2 mu^2/4) ] x [ 11/3 - 2 log(k+/Lambda) ] - 67/9 + pi^2/3")
print()
print("   -67/9 + pi^2/3 =", float(-sp.Rational(67,9)+sp.pi**2/3),
      "  <- the standard NLO BFKL/BK constant")
print("   11/3 = b0/Nc for pure glue, so  Nc x (-2/eps)(11/3) = -2 b0/eps :")
print("      THIS 1/eps IS ULTRAVIOLET -- it renormalises g^2.  (Contrast rows 8-9.)")
print("   and  -2/eps - 2 gamma - log(r^2 mu^2/4) = -2[1/eps + log(r/r_0)],  r_0 = 2 e^-gamma/mu")
print("      exactly the position-space MS-bar dictionary used earlier in this work.")

print(); print("="*76); print("4.  THE log V IN A^(3)'s TWO-rho BRACKET CANCELS"); print("="*76)
print("   coefficient of log V, summing the four V-dependent terms:")
kv,pv,ki,pii = sp.symbols('k2 p2 ki pi_', positive=True); kp = sp.Symbol('kdotp')
kmp = kv - 2*kp + pv
c = ( ((2*pv - 2*kp + kv)*ki + kv*pii)/kmp     # from log[1 + (V/k+-1) k^2/(k-p)^2]
      - (pv*ki + kv*pii)/kmp                   # from -log((V-k+)/Lambda)
      - ki                                     # from -k^i log(V/k+)
      - 2*ki + 2*ki )                          # -2k^i log(V/Lam)  and  -2k^i log(k+p^2/(V k^2+..))
print("      ", sp.simplify(sp.expand(c)), "   <- EXACTLY ZERO")
print("   check: (2p^2-2k.p+k^2-p^2) = (k-p)^2, so that term gives +k^i, cancelling -k^i;")
print("          and the two -2k^i logs cancel because log(k+p^2/(V k^2+k+p^2)) -> -log(V k^2/(k+p^2)).")
