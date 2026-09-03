import numpy as np, sys
sys.path.insert(0,'/home/user/correlations/zint/src')
from twokernel import Kvec, closed_form_vec
rng = np.random.default_rng(7)
for _ in range(4): rng.normal(size=(4,2))
print("3. SYMMETRIES the closed form must have (the integrand is symmetric under")
print("   w<->w' and under x<->y separately, and under (w,w')<->(x,y)):")
P = rng.normal(size=(4,2)); wp,w,x,y = P
print(f"   base             {closed_form_vec(wp,w,x,y):+.12f}")
print(f"   w <-> w'         {closed_form_vec(w,wp,x,y):+.12f}")
print(f"   x <-> y          {closed_form_vec(wp,w,y,x):+.12f}")
print(f"   both             {closed_form_vec(w,wp,y,x):+.12f}")
print(f"   (w,w') <-> (x,y) {closed_form_vec(x,y,wp,w):+.12f}")
print()
print("4. REMOVABLE singularity of the formula at w' = x (0/0, the log vanishes there)")
for d in (1e-2,1e-4,1e-6,1e-8):
    xx = wp + np.array([d,0.0])
    print(f"   |x-w'|={d:7.0e}:  I = {closed_form_vec(wp,w,xx,y):+16.8f}")
print("   -> diverges as log(1/|x-w'|^2): two kernels singular at the SAME point,")
print("      integrand ~ 1/rho^2 there.  Coefficient:")
for d in (1e-4,1e-6,1e-8):
    xx = wp + np.array([d,0.0])
    K2 = ((wp-y)@(w-xx))/(((wp-y)@(wp-y))*((w-xx)@(w-xx)))
    print(f"      |x-w'|={d:7.0e}:  I/log(1/|x-w'|^2) = {closed_form_vec(wp,w,xx,y)/np.log(1/d**2):+.9f}"
          f"    -pi/2 * K(w'-y, w-x) = {-np.pi/2*K2:+.9f}")
print()
print("5. LARGE-|z| POWER COUNTING (why there is no infrared divergence)")
th = np.linspace(0,2*np.pi,8192,endpoint=False)
u = np.stack([np.cos(th),np.sin(th)],-1)
for rho in (1e1,1e2,1e3,1e4):
    z = u*rho
    val = np.mean(Kvec(wp,w,z)*Kvec(x,y,z))
    print(f"   rho={rho:8.0e}   rho^4 x <integrand> = {val*rho**4:+.6f}   rho^2 x <integrand> = {val*rho**2:+.3e}")
print("   -> each kernel falls as 1/rho^2, the product as 1/rho^4; int rho drho/rho^4 converges.")
print()
print("6. SHORT-DISTANCE POWER COUNTING at each of the four points")
for nm,p in (("z->w'",wp),("z->w",w),("z->x",x),("z->y",y)):
    for rr in (1e-2,1e-4):
        z = p + u*rr
        val = np.mean(np.abs(Kvec(wp,w,z)*Kvec(x,y,z)))
        print(f"   {nm:6s} rho={rr:.0e}:  rho x <|integrand|> = {val*rr:.6f}")
print("   -> integrand ~ 1/rho (only ONE kernel is singular at each point);")
print("      int rho drho / rho = int drho converges.  No ultraviolet divergence.")
