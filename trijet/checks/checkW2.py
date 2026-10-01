import sympy as sp

# Represent S(a,b) and Q(a,b,c,d) as commuting symbols keyed by their arguments.
S = lambda a,b: sp.Symbol(f"S({a},{b})")
Q = lambda a,b,c,d: sp.Symbol(f"Q({a},{b},{c},{d})")
def Sd(a,b):                       # S(a,a) = 1  (unitarity)
    return sp.Integer(1) if a==b else S(a,b)
def Qd(a,b,c,d):                   # collapse coincident arguments
    if a==b and c==d: return sp.Integer(1)
    if a==b: return Sd(c,d)
    if c==d: return Sd(a,b)
    return Q(a,b,c,d)

def W(a,b,c, ab,bb,cb):
    """W = S_qqbarg,qqbarg - S_qqbarg - S^dag_qqbarg + 1, large Nc, as in eqs (1.5)-(1.8)."""
    SS  = Qd(a,c,cb,ab)*Qd(c,b,bb,cb)      # S_{qqbarg qqbarg}
    Sam = Sd(a,c)*Sd(c,b)                  # S_{qqbarg}(a,b,c)
    Scc = Sd(cb,ab)*Sd(bb,cb)              # S^dag_{qqbarg}(ab,bb,cb)
    return SS - Sam - Scc + 1

def comb(f, args_full, args_amp, args_conj, args_both):
    return sp.expand(f(*args_full) - f(*args_amp) - f(*args_conj) + f(*args_both))

x,y,z,xb,yb,zb = 'x','y','z','xb','yb','zb'

print("=== 1.1  antiquark x antiquark   (y,z -> y' ; ybar,zbar -> ybar')   xi -> 0 : y'->y ===")
r1 = comb(W, (x,y,z, xb,yb,zb), (x,y,y, xb,yb,zb), (x,y,z, xb,yb,yb), (x,y,y, xb,yb,yb))
print("  ", sp.simplify(r1))

print("\n=== 1.2  quark x quark          (x,z -> x' ; xbar,zbar -> xbar')  xi -> 0 : x'->x ===")
r2 = comb(W, (y,x,z, yb,xb,zb), (y,x,x, yb,xb,zb), (y,x,z, yb,xb,xb), (y,x,x, yb,xb,xb))
print("  ", sp.simplify(r2))

print("\n=== 1.3  interference           (x,z -> x' ; ybar,zbar -> ybar')  xi -> 0 ===")
r3 = comb(W, (x,y,z, xb,yb,zb), (x,y,x, xb,yb,zb), (x,y,z, xb,yb,yb), (x,y,x, xb,yb,yb))
print("  ", sp.simplify(r3))

print("\n=== UV check: does each combination vanish when z -> y (resp. x) BEFORE the xi->0 limit? ===")
# with y' kept general we cannot do this symbolically; but at z=y one has y'=y identically:
c1 = comb(W, (x,y,y, xb,yb,zb), (x,y,y, xb,yb,zb), (x,y,y, xb,yb,yb), (x,y,y, xb,yb,yb))
print("  1.1 at z=y (so y'=y):", sp.simplify(c1))
c1b= comb(W, (x,y,z, xb,yb,yb), (x,y,y, xb,yb,yb), (x,y,z, xb,yb,yb), (x,y,y, xb,yb,yb))
print("  1.1 at zbar=ybar    :", sp.simplify(c1b))
c3 = comb(W, (x,y,x, xb,yb,zb), (x,y,x, xb,yb,zb), (x,y,x, xb,yb,yb), (x,y,x, xb,yb,yb))
print("  1.3 at z=x (so x'=x):", sp.simplify(c3))

print("\n=== coincidence limit  xbar=x, ybar=y, zbar=z  of the three kernels ===")
a1,a2,b1,b2,c1_,c2 = sp.symbols('ax ay bx by cx cy', real=True)  # x, y, z as 2-vectors
X=sp.Matrix([a1,a2]); Y=sp.Matrix([b1,b2]); Z=sp.Matrix([c1_,c2])
A=(Z-X); B=(Z-Y)
Kqq   = (A.dot(A))/((A.dot(A))*(A.dot(A)))
Kaqaq = (B.dot(B))/((B.dot(B))*(B.dot(B)))
Kqaq  = (A.dot(B))/((A.dot(A))*(B.dot(B)))
tot = sp.simplify(Kaqaq + Kqq - 2*Kqaq)
bk  = sp.simplify(((X-Y).dot(X-Y))/((A.dot(A))*(B.dot(B))))
print("  K_aqaq + K_qq - 2 K_qaq - (x-y)^2/((z-x)^2 (z-y)^2) =", sp.simplify(tot-bk))
