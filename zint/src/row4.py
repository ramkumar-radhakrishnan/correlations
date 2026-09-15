"""Row 4: two A-daggers and the gluon B_2, no C.  Prefactor, UV scan, p+ integrals."""
import sympy as sp, numpy as np

# ---------------- 1. prefactor ------------------------------------------------
g, pi = sp.symbols('g pi', positive=True); kp, pp = sp.symbols('kplus pplus', positive=True)
I = sp.I
Adk = -I*g/(sp.sqrt(2)*pi*sp.sqrt(kp))          # A^dag(k+, -w')
Adp = -I*g/(sp.sqrt(2)*pi*sp.sqrt(pp))          # A^dag(p+, -z)
Bpre = I*g**2*sp.sqrt(pp*kp)/(4*pi**2*(pp+kp))  # B_2(k+,-w; p+,-z), bracket factored out
overall = sp.Integer(2)/(2*pi)**3
print("1. PREFACTOR")
print("   signs/i : (-i)(-i)(+i) = -i")
print("   numeric : (2/(2pi)^3) x 1/(sqrt2 pi)^2 x 1/(4 pi^2) =",
      sp.simplify(overall/(2*pi**2)/(4*pi**2)), " vs quoted (1/(2pi)^3)(1/(4 pi^4)) =",
      sp.simplify(1/(2*pi)**3/(4*pi**4)))
print("   k+ powers: (1/sqrt(k+))(1/sqrt(p+)) x sqrt(p+ k+)/(p+ + k+) =",
      sp.simplify(1/sp.sqrt(kp)/sp.sqrt(pp)*sp.sqrt(pp*kp)/(pp+kp)), " = 1/(p+ + k+)")
print("   -> that 1/(p+ + k+) divides B's bracket, turning")
print("        [ d^ij(A-B)/(2(z-w)^2) + ((p+ + k+)/k+)(...) + ((p+ + k+)/p+)(...) ]")
print("      into")
print("        [ d^ij(A-B)/(2(p+ + k+)(z-w)^2) + (1/k+)(...) + (1/p+)(...) ]   <-- the quoted bracket")
tot = sp.simplify(overall*Adk*Adp*Bpre*(pp+kp))   # x (p++k+) since the bracket absorbed 1/(p++k+)
print("   full (bracket normalised) :", tot,
      "   quoted -i g^4/((2pi)^3 4 pi^4) =", sp.simplify(-I*g**4/(2*pi)**3/(4*pi**4)))
print("   difference:", sp.simplify(tot + I*g**4/(2*pi)**3/(4*pi**4)))
print()

# ---------------- 2. the transverse integrand --------------------------------
def T(i, j, x, y, z, w, xi):          # k+ = 1, p+ = xi
    A = (x-z)@(x-z); B = (x-w)@(x-w); zw = z-w; zw2 = zw@zw
    xz, xw = x-z, x-w
    d = 1.0 if i == j else 0.0
    t1 = d*(A-B)/(2*(xi+1)*zw2)
    t2 = (xw[i]*zw[j]/zw2 + xz[j]*xw[i]/(2*A))/xi
    t3 = (xz[j]*zw[i]/zw2 - xz[j]*xw[i]/(2*B))
    return t1 + t2 + t3

def F(x, y, z, w, V, xi):             # V = (x'-w')^i/(x'-w')^2 spectator
    A = (x-z)@(x-z); B = (x-w)@(x-w)
    yz = y-z; yz2 = yz@yz
    den = xi*A + B
    s = 0.0
    for i in range(2):
        for j in range(2):
            s += V[i]*yz[j]/yz2*T(i, j, x, y, z, w, xi)
    return s/den

rng = np.random.default_rng(4)
base = dict(x=np.array([0.3,-0.2]), y=np.array([-0.6,0.9]), z=np.array([0.8,0.5]),
            w=np.array([-0.4,-0.7]), V=np.array([0.4,0.7]))
xi = 0.4
th = np.linspace(0, 2*np.pi, 8192, endpoint=False)
u = np.stack([np.cos(th), np.sin(th)], -1)

print("2. UV SCAN.  rho^2 x <integrand>_theta at each coincidence.")
print("   -> 0 means integrable (no log divergence); nonzero constant means log divergent.")
def scan(name, move, centre):
    print(f"   {name}")
    for rho in (1e-2, 1e-3, 1e-4, 1e-5):
        vals = []
        for uu in u:
            kw = dict(base); kw[move] = centre + uu*rho
            vals.append(F(kw['x'], kw['y'], kw['z'], kw['w'], kw['V'], xi))
        print(f"      rho={rho:.0e}:  rho^2<F> = {np.mean(vals)*rho**2:+.6e}   <F> = {np.mean(vals):+.4e}")
scan("x -> z  (move x about z)", 'x', base['z'])
scan("x -> w  (move x about w)", 'x', base['w'])
scan("z -> w  (move z about w)", 'z', base['w'])
scan("y -> z  (move y about z)", 'y', base['z'])
