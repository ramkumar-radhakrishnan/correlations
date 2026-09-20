import sympy as sp

# logs of transverse quantities (all natural logs of squares)
Txy,Txyp,Txxp,Tyyp,K1,K2,lR,lmu = sp.symbols('Txy Txyp Txxp Tyyp K1 K2 lR lmu')
l1,l2,l0 = sp.symbols('l1 l2 l0')            # ln z1, ln z2, ln z0
eps,gE,pi2,ln2,lnpi = sp.symbols('eps gE pi2 ln2 lnpi')   # pi2 = pi^2
L1, L2 = l1-l0, l2-l0

lc = 2*(ln2-gE)                # ln c0^2 , c0 = 2 e^{-gE}
G  = gE + lnpi + lmu           # ln( e^{gE} pi mu^2 )
Mt = 2*ln2 + lnpi - gE + lmu   # ln( mu~^2 ),  mu~^2 = 4 pi e^{-gE} mu^2
assert sp.simplify(G - Mt + lc) == 0, "MSbar scale relation broken"

# ---- virtual  (B.23) + c.c. ,  overall factor alpha_s C_F/pi  --------------
V = sp.Rational(1,2)*((L1+L2-sp.Rational(3,2))*(2/eps + G + Txy)
                      + sp.Rational(1,2)*(l2-l1)**2 - pi2/6 + 2) \
  + sp.Rational(1,2)*((L1+L2-sp.Rational(3,2))*(2/eps + G + Txyp)
                      + sp.Rational(1,2)*(l2-l1)**2 - pi2/6 + 2)

# ---- real, in-cone  (B.24) + (1<->2) --------------------------------------
def incone(Li, li, Ki):
    return ((sp.Rational(3,4)-Li)*(2/eps) + li**2 - l0**2 - pi2/6
            + (Li-sp.Rational(3,4))*(lR + Ki - Mt - 2*li)
            + sp.Rational(1,4) + sp.Rational(3,2)*(1 - (li-ln2)))
C = incone(L1,l1,K1) + incone(L2,l2,K2)

# ---- real, out-of-cone soft-div  (B.7) + (1<->2) --------------------------
O = (L1**2 - L1*(K1 + Txxp + lR - lc)) + (L2**2 - L2*(K2 + Tyyp + lR - lc))

# ---- target: CSSV (B.25) ---------------------------------------------------
T = (-L1*(Txxp - (Txy+Txyp)/2) - L2*(Tyyp - (Txy+Txyp)/2)
     - sp.Rational(3,4)*(K1+K2+Txy+Txyp-2*lc)
     - 3*sp.Rational(1,2)*lR                    # -3 ln R = -(3/2) ln R^2
     + sp.Rational(1,2)*(l1-l2)**2
     + sp.Rational(11,2) + 3*ln2 - pi2/2)

diff = sp.expand(sp.simplify(V + C + O - T))
print("V + C + O - (B.25)  =", diff)
print()
print("pole (1/eps) coefficient of V+C+O :",
      sp.simplify(sp.expand(V+C+O).coeff(eps,-1)))
print("coefficient of l0**2 :", sp.expand(V+C+O).coeff(l0,2))
print("coefficient of l0*l1 :", sp.expand(V+C+O).coeff(l0,1).coeff(l1,1))
print("coefficient of l0*l2 :", sp.expand(V+C+O).coeff(l0,1).coeff(l2,1))
print("coefficient of lmu   :", sp.simplify(sp.expand(V+C+O).coeff(lmu,1)))
print("coefficient of lR    :", sp.simplify(sp.expand(V+C+O).coeff(lR,1)))
print("full l0-dependence   :", sp.simplify(sp.expand(V+C+O) - sp.expand(V+C+O).subs(l0,0)))
