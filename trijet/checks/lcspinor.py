import sympy as sp
from sympy import I, Rational, sqrt, simplify, expand, factor

# ---------- Dirac matrices (Dirac representation) ----------
Z2 = sp.zeros(2,2); I2 = sp.eye(2)
s1 = sp.Matrix([[0,1],[1,0]]); s2 = sp.Matrix([[0,-I],[I,0]]); s3 = sp.Matrix([[1,0],[0,-1]])
def blk(A,B,C,D): return sp.Matrix(sp.BlockMatrix([[A,B],[C,D]]))
g0 = blk(I2,Z2,Z2,-I2)
g  = [blk(Z2,s,-s,Z2) for s in (s1,s2,s3)]        # gamma^1, gamma^2, gamma^3
gp = g0 + g[2]                                     # gamma^+ = gamma^0 + gamma^3
gm = g0 - g[2]                                     # gamma^-
beta = g0
alpha = [g0*g[i] for i in range(3)]                # alpha^i = gamma^0 gamma^i
def bar(sp_):  return (sp_.H)*g0

# ---------- light-cone spinors (Lepage-Brodsky) ----------
chi = {+1: sp.Matrix([1,0,1,0])/sqrt(2), -1: sp.Matrix([0,1,0,-1])/sqrt(2)}
def u(kp,k1,k2,m,lam):
    return (kp*sp.eye(4) + beta*m + alpha[0]*k1 + alpha[1]*k2)*chi[lam]/sqrt(kp)
def v(kp,k1,k2,m,lam):
    return (kp*sp.eye(4) - beta*m + alpha[0]*k1 + alpha[1]*k2)*chi[-lam]/sqrt(kp)

kp,k1,k2,m = sp.symbols('kplus k_1 k_2 m', positive=True)
print("=== spinor normalisation checks ===")
for lam in (+1,-1):
    U=u(kp,k1,k2,m,lam); V=v(kp,k1,k2,m,lam)
    print(f"  lam={lam:+d}:  ubar u = {str(simplify((bar(U)*U)[0]))}   (expect 2m)",
          f"  ubar g+ u = {str(simplify((bar(U)*gp*U)[0]))} (expect 2k+)",
          f"  vbar v = {str(simplify((bar(V)*V)[0]))} (expect -2m)")
# completeness:  sum_lam u ubar = slash(k) + m   with k^- = (k_perp^2+m^2)/k^+
km = (k1**2+k2**2+m**2)/kp
slashk = Rational(1,2)*(gp*km + gm*kp) - g[0]*k1 - g[1]*k2
S = sum((u(kp,k1,k2,m,l)*bar(u(kp,k1,k2,m,l)) for l in (+1,-1)), sp.zeros(4,4))
print("  sum_lam u ubar - (slash k + m) =", simplify(S - (slashk + m*sp.eye(4))))

print("\n" + "="*78)
print("VERTEX 1:  gamma*_T(l) -> q(k1,lam1) + qbar(k2,lam2)   [massive]")
print("="*78)
qp, z, P1, P2, mm = sp.symbols('qplus zeta P_1 P_2 m', positive=True)
eps = sp.Matrix([[0,1],[-1,0]])           # eps^{12}=+1
lam = sp.Symbol('lambda')
res = {}
for l1 in (+1,-1):
    for l2 in (+1,-1):
        U = u(z*qp, P1, P2, mm, l1)
        V = v((1-z)*qp, -P1, -P2, mm, l2)
        row = [simplify(expand((bar(U)*g[l]*V)[0])) for l in (0,1)]
        res[(l1,l2)] = row

print("  amplitude  V^l = ubar(k1,lam1) gamma^l v(k2,lam2),  normalised by sqrt(zeta(1-zeta)) qplus :")
N = sqrt(z*(1-z))*qp
for (l1,l2),row in res.items():
    r = [simplify(expand(c*N/qp)) for c in row]           # strip 1/sqrt(z(1-z)) and qplus
    tag = "lam2 = -lam1 (massless-type)" if l2==-l1 else "lam2 = +lam1 (MASS term)"
    print(f"   lam1={l1:+d} lam2={l2:+d}  [{tag}]")
    print(f"      l=1: {sp.factor(simplify(r[0]))}")
    print(f"      l=2: {sp.factor(simplify(r[1]))}")

print("\n  --- decomposition test ---")
print("  claim (helicity non-flip, lam2=-lam1):  V^l = (2 zeta - 1) delta^{li} P^i + i lam1 eps^{li} P^i   (times 1/sqrt(z zbar))")
for l1 in (+1,-1):
    row = res[(l1,-l1)]
    for l,name in ((0,'l=1'),(1,'l=2')):
        P = [P1,P2]
        claim = (2*z-1)*P[l] + I*l1*sum(eps[l,i]*P[i] for i in (0,1))
        got   = simplify(expand(row[l]*sqrt(z*(1-z))))
        # allow overall constant factor
        ratio = simplify(got/claim)
        print(f"    lam1={l1:+d} {name}: got/claim = {sp.simplify(ratio)}")
print("\n  claim (mass term, lam2=+lam1): V^l propto m, no P dependence")
for l1 in (+1,-1):
    row = res[(l1,l1)]
    got = [simplify(expand(c*sqrt(z*(1-z)))) for c in row]
    print(f"    lam1={l1:+d}:  l=1: {sp.factor(got[0])}   l=2: {sp.factor(got[1])}")
