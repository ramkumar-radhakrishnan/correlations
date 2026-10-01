import sympy as sp
from sympy import I, Rational, sqrt, simplify, expand, factor
exec(open('lcspinor.py').read().split("print(\"\\n\" + \"=\"*78)")[0])   # reuse spinor machinery

print("="*78)
print("VERTEX 2:  q(p,lam) -> q(p',lam') + g(k,m)   in LC gauge A^+=0   [massive quark]")
print("="*78)
qp = sp.Symbol('qplus', positive=True)
Pf, xi, mm = sp.symbols('P xi m', positive=True)     # parent fraction P, gluon fraction xi
k1s,k2s = sp.symbols('kappa_1 kappa_2', real=True)   # gluon transverse momentum (parent at rest transversely)

def slash_eps(mIdx, kperp, kplus):
    """LC-gauge polarisation: eps^+ = 0, eps^m = delta^{m,mIdx}, eps^- = 2 (eps.k)/k^+"""
    epsm = [sp.Integer(1) if j==mIdx else sp.Integer(0) for j in (0,1)]
    epsminus = 2*(epsm[0]*kperp[0]+epsm[1]*kperp[1])/kplus
    return Rational(1,2)*(gp*epsminus) - g[0]*epsm[0] - g[1]*epsm[1]

kperp = [k1s,k2s]
out = {}
for lam in (+1,-1):
    for lamp in (+1,-1):
        Uin  = u(Pf*qp, 0, 0, mm, lam)                       # parent: p_perp = 0
        Uout = u((Pf-xi)*qp, -k1s, -k2s, mm, lamp)           # daughter
        row = []
        for mIdx in (0,1):
            SE = slash_eps(mIdx, kperp, xi*qp)
            row.append(simplify(expand((bar(Uout)*SE*Uin)[0])))
        out[(lam,lamp)] = row

print("\n  V^m = ubar(p',lam') slash(eps*_m) u(p,lam)   [common factor sqrt(P(P-xi)) qplus stripped]")
pref = sqrt(Pf*(Pf-xi))*qp
for (lam,lamp),row in out.items():
    tag = "lam'=lam  (helicity conserving)" if lamp==lam else "lam'=-lam (MASS / helicity flip)"
    print(f"   lam={lam:+d} lam'={lamp:+d}  [{tag}]")
    for mIdx,nm in ((0,'m=1'),(1,'m=2')):
        print(f"      {nm}: {sp.factor(simplify(row[mIdx]*pref/qp))}")

print("\n  --- structure test (helicity conserving) ---")
print("  claim:  V^m  propto  [ (2(P-xi)+xi) delta^{rm} + i lam xi eps^{rm} ] kappa^r / (xi * P ... )")
epsT = sp.Matrix([[0,1],[-1,0]])
for lam in (+1,-1):
    row = out[(lam,lam)]
    for mIdx,nm in ((0,'m=1'),(1,'m=2')):
        claim = ( (2*(Pf-xi)+xi)*kperp[mIdx]
                  + I*lam*xi*sum(epsT[r,mIdx]*kperp[r] for r in (0,1)) )
        got = simplify(expand(row[mIdx]*pref/qp))
        print(f"    lam={lam:+d} {nm}:  got/claim = {sp.simplify(sp.cancel(got/claim))}")

print("\n  --- structure test (mass / helicity-flip term) ---")
for lam in (+1,-1):
    row = out[(lam,-lam)]
    got = [simplify(expand(row[i]*pref/qp)) for i in (0,1)]
    print(f"    lam={lam:+d}:  m=1: {sp.factor(got[0])}    m=2: {sp.factor(got[1])}")
