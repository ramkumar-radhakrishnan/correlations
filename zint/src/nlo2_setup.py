"""The new row:  four Weizsacker-Williams kernels contracted through
      [ d_{k'm} d_{ij}/(p+ + k+)  -  d_{jm} d_{ik'}/p+  -  d_{im} d_{jk'}/k+ ].

   K1^m = (y'-z)^m/(y'-z)^2      K2^k' = (x'-y')^k'/(x'-y')^2
   K3^i = (y-w)^i/(y-w)^2        K4^j  = (x-z)^j/(x-z)^2

   With xi = p+/k+ the bracket is (1/k+)[ d_{k'm}d_{ij}/(1+xi) - d_{jm}d_{ik'}/xi
                                          - d_{im}d_{jk'} ], so the transverse
   structure is a sum of three products of two scalar dipole kernels:

     T1 = (K1.K2)(K3.K4)/(1+xi)      T2 = -(K1.K4)(K2.K3)/xi      T3 = -(K1.K3)(K2.K4)
"""
import sympy as sp

# --- symbolic check that the bracket really collapses to those three products ---
yp, z, xp, x, y, w = [sp.Matrix(sp.symbols(f'{n}1 {n}2', real=True))
                      for n in ('yp','z','xp','xx','yy','ww')]
xi = sp.symbols('xi', positive=True)
def K(a, b):                      # (a-b)^i/(a-b)^2
    d = a - b; return d/ (d.T*d)[0,0]
K1, K2, K3, K4 = K(yp,z), K(xp,yp), K(y,w), K(x,z)
dot = lambda u,v: (u.T*v)[0,0]

brack = 0
for m in range(2):
    for kp in range(2):
        for i in range(2):
            for j in range(2):
                d = lambda a,b: 1 if a==b else 0
                coef = (sp.Rational(1,1)/(1+xi))*d(kp,m)*d(i,j) \
                     - (1/xi)*d(j,m)*d(i,kp) - d(i,m)*d(j,kp)
                if coef != 0:
                    brack += K1[m]*K2[kp]*K3[i]*K4[j]*coef
T1 = dot(K1,K2)*dot(K3,K4)/(1+xi)
T2 = -dot(K1,K4)*dot(K2,K3)/xi
T3 = -dot(K1,K3)*dot(K2,K4)
print("bracket == T1 + T2 + T3 :", sp.simplify(brack - (T1+T2+T3)) == 0)

# --- power counting: every coincidence limit is 1/|sep|, hence integrable in 2D ---
print()
print("power counting of the transverse integrand (four WW kernels, no squared")
print("denominators anywhere):")
for nm, pair in [("y' -> z", "K1"), ("x' -> y'", "K2"), ("y -> w", "K3"), ("x -> z", "K4")]:
    print(f"   {nm:10s}: only {pair} is singular, ~ 1/|sep| ; int d^2(sep) ~ int rho drho/rho  -> finite")
print("   no pair of points appears in two kernels, so no 1/sep^2 anywhere:")
print("   the z and y' integrations are ULTRAVIOLET FINITE.")
