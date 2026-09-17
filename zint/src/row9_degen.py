"""Row 9: the large-|z| tail of B_2^(1) IS B_2^(2)'s kernel -- pointwise, and in the prefactor."""
import numpy as np, sympy as sp
P,K = 0.37,1.0; S=P+K
def Wji(x,z,w):
    xz,xw,zw = x-z,x-w,z-w; s=lambda v: v@v
    return (np.eye(2)*(s(xz)-s(xw))/(2*s(zw))
            + (S/K)*(np.outer(xz,zw)/s(zw) - np.outer(xz,xw)/(2*s(xw)))
            + (S/P)*(np.outer(zw,xw)/s(zw) + np.outer(xz,xw)/(2*s(xz))))
Dd = lambda x,z,w: P*((x-z)@(x-z))+K*((x-w)@(x-w))
Kv = lambda a,b: (a-b)/((a-b)@(a-b))
x = np.array([0.31,-0.77]); w = np.array([1.13,0.42])
print("="*78); print("POINTWISE (fixed direction, not angular averaged)"); print("="*78)
print("   W_{ji}(x,z,w)/[S D]   vs   -(1/(2 p+ k+)) (x-z)^j (x-w)^i/[(x-z)^2 (x-w)^2]")
print("   %-10s %15s %15s %10s   %15s %10s" % ("|z|","exact[0,0]","limit[0,0]","ratio","exact[1,0]","ratio"))
dirn = np.array([np.cos(0.7), np.sin(0.7)])
for R in (1e1,1e2,1e3,1e4,1e5,1e6):
    z = R*dirn
    e = Wji(x,z,w)/(S*Dd(x,z,w))
    l = -np.outer(Kv(x,z), Kv(x,w))/(2*P*K)
    print("   %-10.0e %15.6e %15.6e %10.6f   %15.6e %10.6f"
          % (R, e[0,0], l[0,0], e[0,0]/l[0,0], e[1,0], e[1,0]/l[1,0]))
print("   -> ratio -> 1 in every component.")

print(); print("="*78); print("THE PREFACTORS COINCIDE TOO"); print("="*78)
g,pi = sp.symbols('g pi',positive=True); Ps,Ks = sp.symbols('pplus kplus',positive=True); Ss=Ps+Ks
B1pref = sp.I*g**2*sp.sqrt(Ps*Ks)/(4*pi**2*Ss)
tail   = sp.simplify(B1pref*(-Ss/(2*Ps*Ks)))      # W/D -> -(S/(2 p+ k+)) x kernel
B2pref = -g**2/(8*pi**2*sp.sqrt(Ps*Ks))
print("   B^(1) prefactor x the large-|z| weight  -S/(2 p+ k+) :", sp.simplify(tail))
print("   B^(2) prefactor                                      :", sp.simplify(B2pref))
print("   ratio                                                :", sp.simplify(tail/B2pref))
print()
print("   => at large |z| the one-rho vertex becomes")
print("        -i g^2/(8 pi^2 sqrt(p+k+))  int_x f^{cba} rho^a(x) (x-w)^i (x-z)^j/[(x-w)^2 (x-z)^2]")
print("      and the two-rho vertex is")
print("        -  g^2/(8 pi^2 sqrt(p+k+))  int_{x,y} {rho^c(y),rho^b(x)} (y-w)^i (x-z)^j/[..]")
print("      SAME kernel, SAME 1/(8 pi^2 sqrt(p+ k+)).  They differ only by")
print("        -i f^{cba} rho^a(x)   [one rho, both kernels at the SAME point x]")
print("      versus")
print("        -  {rho^c(y), rho^b(x)}   [two rho, kernels at INDEPENDENT points]")
