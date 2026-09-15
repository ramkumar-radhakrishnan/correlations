import numpy as np, sympy as sp, sys, io, contextlib
sys.path.insert(0,'/home/user/correlations/zint/src')
with contextlib.redirect_stdout(io.StringIO()):
    from row4 import F, T
rng = np.random.default_rng(4)
x0=np.array([0.3,-0.2]); y0=np.array([-0.6,0.9]); z0=np.array([0.8,0.5])
w0=np.array([-0.4,-0.7]); V0=np.array([0.4,0.7]); xi=0.4
th=np.linspace(0,2*np.pi,8192,endpoint=False); u=np.stack([np.cos(th),np.sin(th)],-1)

print("3. TRIPLE COLLAPSE  x, z, w all together (the Feynman denominator's zero)")
print("   integrand ~ 1/rho^2 there, but 3 points -> 2 relative coords -> measure rho^4 drho/rho,")
print("   net rho^2 drho/rho : convergent.  Check the 1/rho^2:")
c = np.array([0.1,0.1])
for rho in (1e-2,1e-3,1e-4):
    vals=[F(c+uu*rho, y0, c+np.roll(uu,1)*rho, c-uu*rho, V0, xi) for uu in u]
    print(f"      rho={rho:.0e}:  rho^2<F> = {np.mean(vals)*rho**2:+.6e}")
print()
print("4. LARGE DISTANCE")
for nm,mv in (("|z| -> inf",'z'), ("|w| -> inf",'w')):
    print(f"   {nm}")
    for R in (1e2,1e3,1e4):
        vals=[]
        for uu in u:
            a=dict(x=x0,y=y0,z=z0,w=w0)
            a[mv]=uu*R
            vals.append(F(a['x'],a['y'],a['z'],a['w'],V0,xi))
        m=np.mean(vals)
        print(f"      R={R:.0e}:  R^2<F> = {m*R**2:+.4e}   R^3<F> = {m*R**3:+.4e}")
print()

print("5. THE p+ INTEGRALS  (there is NO phase here, so they are elementary)")
p,L,P,A,B,k = sp.symbols('pplus Lambda P A B kplus', positive=True)
I1 = sp.integrate(1/((p+k)*(p*A+k*B)), (p, L, P))
I2 = sp.integrate(1/(p*(p*A+k*B)), (p, L, P))
I3 = sp.integrate(1/(p*A+k*B), (p, L, P))
for nm, e, closed in [
  ("I1 = int dp+ / [(p+ + k+)(p+ A + k+ B)]", I1,
   sp.log((P+k)/(L+k))/(k*(B-A)) - sp.log((P*A+k*B)/(L*A+k*B))/(k*(B-A))),
  ("I2 = int dp+ / [p+ (p+ A + k+ B)]", I2,
   sp.log(P*(L*A+k*B)/(L*(P*A+k*B)))/(k*B)),
  ("I3 = int dp+ / [p+ A + k+ B]", I3,
   sp.log((P*A+k*B)/(L*A+k*B))/A)]:
    print(f"   {nm}")
    print(f"      closed form: {closed}")
    subs = {A:1.7, B:0.9, k:1.0, L:1e-6, P:37.0}
    print(f"      sympy {float(e.subs(subs)):.10f}   closed {float(closed.subs(subs)):.10f}"
          f"   diff {abs(float(e.subs(subs))-float(closed.subs(subs))):.2e}")
print()
print("   Lambda -> 0 limits:")
print("      I1 -> log( V B / [ (V-k+) A + k+ B ] ) / [ k+ (B - A) ]        FINITE")
print("      I2 -> [ log((V-k+)/Lambda) + log( k+ B / [(V-k+)A + k+ B] ) ] / (k+ B)   <- the rapidity log")
print("      I3 -> log( [(V-k+) A + k+ B] / (k+ B) ) / A                     FINITE")
for nm, closed, lim in [
  ("I1", sp.log((P+k)/(L+k))/(k*(B-A)) - sp.log((P*A+k*B)/(L*A+k*B))/(k*(B-A)),
        sp.log((P+k)*k*B/(k*((P*A+k*B))))/(k*(B-A))),
  ("I3", sp.log((P*A+k*B)/(L*A+k*B))/A, sp.log((P*A+k*B)/(k*B))/A)]:
    for lam in (1e-4,1e-6,1e-8):
        s={A:1.7,B:0.9,k:1.0,L:lam,P:37.0}
        print(f"      {nm} at Lambda={lam:.0e}: {float(closed.subs(s)):.10f}   limit {float(lim.subs(s)):.10f}")
