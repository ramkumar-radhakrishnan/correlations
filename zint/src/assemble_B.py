"""Row with B_2 in the ket instead of two A's: prefactor, colour, and the large-|z| limit."""
import sympy as sp, numpy as np
g, pi = sp.symbols('g pi', positive=True)
kp, pp = sp.symbols('kplus pplus', positive=True); Kp = pp + kp
I = sp.I

print("A. PREFACTOR")
Adag = -I*g/(sp.sqrt(2)*pi*sp.sqrt(Kp))                      # A^dag(K+)
Cpre = -g/(2*pi*sp.sqrt(2*Kp))*sp.sqrt(pp*(Kp-pp))/Kp        # C, first argument K+
Bpre = -g**2/(8*pi**2*sp.sqrt(pp*kp))                        # B_2
overall = sp.Integer(4)/(2*pi)**3                            # the 4/(2pi)^3 out front
prod = sp.simplify(overall*Adag*Cpre*Bpre)
print("   signs/i :  (-i) from A^dag, (-1) from C, (-1) from B  ->  (-i)(-1)(-1) = -i")
print("   numeric :  (4/(2pi)^3) x 1/(sqrt2 pi) x 1/(2 sqrt2 pi) x 1/(8 pi^2)")
print(f"            = {sp.simplify(overall/(sp.sqrt(2)*pi)/(2*sp.sqrt(2)*pi)/(8*pi**2))}"
      f"   vs quoted (1/(2pi)^3)(1/(8 pi^4)) = {sp.simplify(1/(2*pi)**3/(8*pi**4))}")
print("   k+ powers: B's 1/sqrt(p+ k+) cancels C's sqrt(p+ k+) -- exactly as the two A's did")
print(f"            product = {sp.simplify(prod)}")
J = (Kp/kp)**2
total = sp.simplify(prod*J)
target = -sp.Integer(1)/(2*pi)**3*I*g**4/(8*pi**4)/kp**2
print(f"   x delta Jacobian (K+/k+)^2 : {total}")
print(f"   quoted                     : {sp.simplify(target)}")
print(f"   difference                 : {sp.simplify(total-target)}    <-- 0 means correct")
print()

print("B. COLOUR: the index order inside B_2")
print("   B^{bc}_{2ji}(k+,-w; p+,-z) carries { rho^c(y), rho^b(x) } : FIRST index -> rho(x),")
print("   SECOND index -> rho(y).  So Bbar^{ca} -> U^{ce}(x) U^{ad}(y) { rho^d(y), rho^e(x) },")
print("   and  -U^{bc}(z) Bbar^{ca} = -U^{bc}(z) U^{ce}(x) U^{ad}(y) {rho^d(y),rho^e(x)}   MATCHES")
print("   term 1: U^{ac}(w) B^{bc} = U^{ac}(w) {rho^c(y), rho^b(x)}                       MATCHES")
print("   note term 1 carries NO U(z) at all.")
print()

print("C. LARGE |z| : does the colour factor vanish?")
rng = np.random.default_rng(17); N2 = 8
def rand_adj():
    A = rng.normal(size=(N2,N2)); Q,_ = np.linalg.qr(A); return Q
Uy, Uw, Ux, Uyp, Uxp, Uz = [rand_adj() for _ in range(6)]
rx, ry, rxp = [rng.normal(size=N2) for _ in range(3)]
E = np.eye(N2)
F = rng.normal(size=(N2,N2,N2))
F = (F - F.transpose(1,0,2) + F.transpose(1,2,0) - F.transpose(2,1,0)
     + F.transpose(2,0,1) - F.transpose(0,2,1))/6
anti = np.einsum('d,e->de', ry, rx) + np.einsum('e,d->de', rx, ry)   # {rho^d(y), rho^e(x)} (c-number model)
def total(Uz):
    bra = (np.einsum('ec,db,e->cdb', Uyp, Uz, rxp) - np.einsum('bd,ce,e->cdb', Uz, Uxp, rxp))
    ket = (np.einsum('ac,cb->ab', Uw, anti)                          # U^{ac}(w){rho^c(y),rho^b(x)}
           - np.einsum('bc,ce,ad,de->ab', Uz, Ux, Uy, anti))
    return np.einsum('CDA,CDB,AB->', F, bra, ket)
print(f"   full contraction at U(z) = 1      : {total(E):+.10f}")
print(f"   full contraction at a generic U(z): {total(Uz):+.10f}")
print("   the U(z) -> 1 limit is z-independent and NONZERO, so the coefficient of the")
print("   large-|z| logarithm survives -- same as the previous two rows.")
print()
print("D. which U(z) products collapse by unitarity?")
for nm, M in [("bra1 x ket1 : U^{d'b}(z) x 1   ", Uz),
              ("bra1 x ket2 : U^{d'b}(z)U^{bc}(z)", np.einsum('db,bc->dc', Uz, Uz)),
              ("bra2 x ket1 : U^{bd'}(z) x 1   ", Uz.T),
              ("bra2 x ket2 : U^{bd'}(z)U^{bc}(z)", np.einsum('bd,bc->dc', Uz, Uz))]:
    ok = np.abs(M-E).max() < 1e-10
    print(f"   {nm}  {'-> delta' if ok else '-> still z-dependent'}")
