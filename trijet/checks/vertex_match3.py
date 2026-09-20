import sympy as sp
z,x = sp.symbols('zeta xi'); e = 1-z-x; zb = 1-z
am,ac,bm,bc,cm,cc = sp.symbols('a_m a_c b_m b_c c_m c_c')
dd = 2*(am*ac+1)*(bm*bc+cm*cc); ee = -2*(am+ac)*(cm*bc+bm*cc)
def test(name, s, dt, et):
    d=sp.simplify(sp.expand(dd.subs(s)-dt)); p=sp.simplify(sp.expand(ee.subs(s)-et))
    ok = (d==0 and p==0); print(f"  {name:38s} {'MATCH' if ok else 'dd='+str(sp.factor(d))+'  ee='+str(sp.factor(p))}")

# PHYSICAL assignment:
#   photon vertex   a = 2*z_q - 1          (z_q = quark's fraction of the photon)
#   emission vertex b = 2*(P-xi) + xi ,    c = xi * lambda_emitter/lambda_1
#   massless helicity conservation at the photon vertex => lambda_qbar = -lambda_q
#     -> quark  emits :  a = 2(1-eta)-1 = 1-2eta ,  c = +xi  (emitter helicity = +lambda_1)
#     -> antiq. emits :  a = 2 zeta -1            ,  c = -xi  (emitter helicity = -lambda_1)
Aq   = {'a': 1-2*e,   'b': 2*z+x, 'c': +x}     # quark emits the gluon
Aqb  = {'a': 2*z-1,   'b': 2*e+x, 'c': -x}     # antiquark emits the gluon

print("=== physical assignment: a = 2 z_q - 1, c = xi * (emitter helicity)/lambda_1 ===")
test("(1.2)  qbar x qbar", {am:Aqb['a'],ac:Aqb['a'],bm:Aqb['b'],bc:Aqb['b'],cm:Aqb['c'],cc:Aqb['c']},
     8*(e**2+zb**2)*(z**2+zb**2), -8*x*(2*z-1)*(2*zb-x))
test("(1.10) q x q",       {am:Aq['a'], ac:Aq['a'], bm:Aq['b'], bc:Aq['b'], cm:Aq['c'], cc:Aq['c']},
     8*(z**2+(1-e)**2)*(e**2+(1-e)**2), -8*x*(2*e-1)*(2*(1-e)-x))
test("(1.18) interference (q amp x qbar conj)",
     {am:Aq['a'], ac:Aqb['a'], bm:Aq['b'], bc:Aqb['b'], cm:Aq['c'], cc:Aqb['c']},
     8*(1-x-2*z*e)*(x*(1-x)+2*z*e), -8*x*(1-2*z-x)**2)

print("\n=== the emission vertex in standard splitting variables ===")
P,w = sp.symbols('P w', positive=True)   # parent fraction P, gluon takes w = xi/P of it
print("  b = 2(P-xi)+xi = P(2-w) ,  c = xi = P w   ->  V_g^{rm} = P[(2-w) d^{rm} + i lam w eps^{rm}]")
print("  check qbar channel: P = 1-zeta, w = xi/(1-zeta):",
      sp.simplify(( (1-z)*(2-x/(1-z)) ) - (2*e+x)), "(should be 0)")
print("  check q    channel: P = 1-eta , w = xi/(1-eta):",
      sp.simplify(( (1-e)*(2-x/(1-e)) ) - (2*z+x)), "(should be 0)")
