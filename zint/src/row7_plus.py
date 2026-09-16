"""Can the + prescription be applied to row 7?  Check BOTH subtraction terms."""
import numpy as np
kk = np.array([0.9, 0.6])

print("The three index structures contract the four kernels as:")
print("   d_im d_jk  ->  (x'-w').(z-x)  x  (y-z).(y'-x)          [the '1' term]")
print("   d_jm d_ik  ->  (x'-w').(y-z)  x  (z-x).(y'-x)          [the 1/xi term]")
print("   d_ij d_km  ->  (x'-w').(y'-x) x  (y-z).(z-x)           [the 1/(1-xi) term]")
print()
print("Both singular structures leave a BK REAL-EMISSION KERNEL:")
print("   1/xi     term:  (z-x).(y'-x)/[(z-x)^2 (y'-x)^2] = K(z,y'; x)   emission point x")
print("   1/(1-xi) term:  (y-z).(z-x)/[(y-z)^2 (z-x)^2]  = K(y,x ; z)   emission point z")
print()

th = np.linspace(0,2*np.pi,8192,endpoint=False); u = np.stack([np.cos(th),np.sin(th)],-1)
z0 = np.array([0.8,0.5]); yp0 = np.array([-0.3,0.4])
y0 = np.array([-0.6,0.9]); x0 = np.array([0.3,-0.2])

def K(a,b,e):                       # (a-e).(b-e)/[(a-e)^2 (b-e)^2]
    A=a-e; B=b-e
    return np.sum(A*B,-1)/(np.sum(A*A,-1)*np.sum(B*B,-1))

print("A. the xi -> 0 subtraction term.  Its x-integrand is K(z,y';x), phase e^{+i xi k.x} -> 1")
print(f"   {'R':>8} {'R^2<K> with phase, xi=0.3':>28} {'R^2<K> at xi = 0 (no phase)':>30}")
for R in (1e2,1e3,1e4,1e5):
    xs = u*R
    with_ph = np.mean(np.exp(1j*0.3*(xs@kk))*K(z0,yp0,xs))*R**2
    no_ph   = np.mean(K(z0,yp0,xs))*R**2
    print(f"   {R:8.0e} {abs(with_ph):28.6e} {no_ph:30.8f}")
print("   -> with the phase it decays; at xi = 0 it sits on +1/R^2 : LOG DIVERGENT")
print()
print("B. the xi -> 1 subtraction term.  Its z-integrand is K(y,x;z), phase e^{+i(1-xi)k.z} -> 1")
print(f"   {'R':>8} {'R^2<K> with phase, xi=0.7':>28} {'R^2<K> at xi = 1 (no phase)':>30}")
for R in (1e2,1e3,1e4,1e5):
    zs = u*R
    with_ph = np.mean(np.exp(1j*0.3*(zs@kk))*K(y0,x0,zs))*R**2
    no_ph   = np.mean(K(y0,x0,zs))*R**2
    print(f"   {R:8.0e} {abs(with_ph):28.6e} {no_ph:30.8f}")
print("   -> same: at xi = 1 it sits on +1/R^2 : LOG DIVERGENT")
print()
print("C. and the subtracted combination keeps the phase-less -1:")
for nm, a, b, xiv in (("xi->0, (e^{i xi k.x} - 1) K(z,y';x)", z0, yp0, 0.3),
                      ("xi->1, (e^{i(1-xi)k.z} - 1) K(y,x;z)", y0, x0, 0.3)):
    row=[]
    for R in (1e3,1e4,1e5):
        pts=u*R
        row.append(np.mean((np.exp(1j*xiv*(pts@kk))-1)*K(a,b,pts))*R**2)
    print(f"   {nm}")
    print(f"      R^2<.> = " + "  ".join(f"{v:+.6f}" for v in row) + "   -> -1 : still LOG DIVERGENT")
