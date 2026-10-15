import sympy as sp
z,x = sp.symbols('zeta xi')
e   = 1 - z - x                      # eta = 1 - zeta - xi  (momentum conservation imposed)
zb  = 1 - z                          # = eta + xi
am,ac,bm,bc,cm,cc = sp.symbols('a_m a_c b_m b_c c_m c_c')
dd = 2*(am*ac+1)*(bm*bc+cm*cc)
ee = -2*(am+ac)*(cm*bc+bm*cc)
def test(name, subs, dd_t, ee_t):
    d = sp.simplify(sp.expand(dd.subs(subs)-dd_t)); p = sp.simplify(sp.expand(ee.subs(subs)-ee_t))
    print(f"  {name:32s} dd: {'MATCH' if d==0 else sp.factor(d)}    ee: {'MATCH' if p==0 else sp.factor(p)}")

print("=== on-shell (zeta+eta+xi=1) ===")
test("(1.2)  qbar x qbar",
     {am:2*z-1, ac:2*z-1, bm:2*e+x, bc:2*e+x, cm:x, cc:x},
     8*(e**2+zb**2)*(z**2+zb**2), -8*x*(2*z-1)*(2*zb-x))
test("(1.10) q x q",
     {am:2*e-1, ac:2*e-1, bm:2*z+x, bc:2*z+x, cm:x, cc:x},
     8*(z**2+(1-e)**2)*(e**2+(1-e)**2), -8*x*(2*e-1)*(2*(1-e)-x))
for sm,sc,tag in [(1,1,'+,+'),(1,-1,'+,-'),(-1,1,'-,+'),(-1,-1,'-,-')]:
    test(f"(1.18) interference [c: {tag}]",
         {am:2*(1-e)-1, ac:2*z-1, bm:2*z+x, bc:2*e+x, cm:sm*x, cc:sc*x},
         8*(1-x-2*z*e)*(x*(1-x)+2*z*e), -8*x*(1-2*z-x)**2)
