import sympy as sp
# exact identity used in section 4.2 of trijet_divergence_separation.pdf
c = sp.symbols('x1 x2 y1 y2 z1 z2 xb1 xb2 yb1 yb2 zb1 zb2', real=True)
X,Y,Z   = sp.Matrix(c[0:2]), sp.Matrix(c[2:4]), sp.Matrix(c[4:6])
Xb,Yb,Zb= sp.Matrix(c[6:8]), sp.Matrix(c[8:10]), sp.Matrix(c[10:12])
a,b,ab,bb = Z-X, Z-Y, Zb-Xb, Zb-Yb
K = lambda u,v: u.dot(v)/(u.dot(u)*v.dot(v))
A  = a/a.dot(a)  - b/b.dot(b)          # eikonal emission amplitude
Ab = ab/ab.dot(ab) - bb/bb.dot(bb)
lhs = K(b,bb) + K(a,ab) - K(a,bb) - K(b,ab)
print("K_qbqb + K_qq - K_qqb - K_qbq  -  A.Abar =", sp.simplify(lhs - A.dot(Ab)))
sub = dict(zip(c[6:12], c[0:6]))
print("coincidence limit - (x-y)^2/((z-x)^2 (z-y)^2) =",
      sp.simplify(lhs.subs(sub) - (X-Y).dot(X-Y)/(a.dot(a)*b.dot(b))))
