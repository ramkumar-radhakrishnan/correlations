"""Misc checks:
(1) 4rho_corrected: two-Wilson-line bracket (line 72) vanishes; four-Wilson-line bracket
    (line 82) vanishes after x<->y symmetrization of the kernel.
(2) Group II Row IV z-integral identity (line 110): verify the FINITE traceless part
    numerically in d=2:  int_z (y-z)^m (y'-z)^m' / ((y-z)^2 (y'-z)^2)
    = -(pi/2) [UV/IR log] delta^mm' - pi r^m r^m' / r^2 ,  r = y'-y.
    The traceless projection is finite: J^mm' - (tr J/2) delta^mm' -> -pi (r^m r^m'/r^2 - delta^mm'/2).
(3) Group I Row II N2: check that eq(58)-bracket minus R^i (eq 74-79) is proportional to
    log(k+/Lambda) (i.e. no residual V-dependence), so the 2a/2b split is consistent.
"""
import mpmath as mp
import sympy as sp
import random

PI = mp.pi

# ---------- (1) 4rho cancellations ----------
def check_4rho_2W():
    # U entries are commuting numbers; rho factored out front (fixed order)
    rng = random.Random(7)
    Axp, Bu, Awp, Bw = [mp.mpf(rng.uniform(0.3, 2.0)) for _ in range(4)]
    # A..=U^{ab'}(x'), B..=U^{ab}(u), etc. with a summed -> model each product as independent?
    # Since every term is U^{ab'}(arg1) U^{ab}(arg2) summed over a, model the sum over a by
    # a dot product of random vectors indexed by argument:
    N = 5
    Uxp = [mp.mpf(rng.uniform(-1, 1)) for _ in range(N)]   # U^{ab'}(x') as vector in a
    Uwp = [mp.mpf(rng.uniform(-1, 1)) for _ in range(N)]
    Uu = [mp.mpf(rng.uniform(-1, 1)) for _ in range(N)]
    Uw = [mp.mpf(rng.uniform(-1, 1)) for _ in range(N)]
    d = lambda P, Q: sum(p*q for p, q in zip(P, Q))
    total = (-d(Uxp, Uu) + d(Uwp, Uu) + d(Uw, Uxp) - d(Uu, Uxp) + d(Uxp, Uw) - d(Uwp, Uw)
             + d(Uu, Uwp) - d(Uw, Uwp) + 2*d(Uw, Uwp) + 2*d(Uxp, Uu)
             - 2*d(Uxp, Uw) - 2*d(Uu, Uwp))
    print("4rho 2-Wilson-line bracket sum =", mp.nstr(total, 5), "-> PASS" if abs(total) < 1e-20 else "-> FAIL")

def check_4rho_4W():
    # line 82 bracket, with kernel symmetric under x<->y (and dummy d<->e, rho relabel).
    # Model U^{ab'}(s) U^{ac}(t) U^{bd}(u) U^{be}(v) as X(s,t)*Y(u,v) with
    # X(s,t)=sum_a U^{ab'}(s)U^{ac}(t), Y(u,v)=sum_b U^{bd}(u)U^{be}(v) -- but d,e indices attach
    # to rho^d(x) rho^e(y). Under x<->y symmetrization we relabel d<->e, i.e. Y(u,v)->Y(v,u) AND
    # swap coordinate symbols x<->y inside arguments. Represent arguments symbolically.
    xs = sp.symbols("x y z w wp xp u", commutative=True)
    x, y, z, w, wp_, xp_, u = xs
    X = sp.Function("X")   # X(arg_of_b', arg_of_c)
    Y = sp.Function("Y")   # Y(arg_of_d, arg_of_e)
    expr = (4*X(wp_, w)*Y(z, y) + 2*X(xp_, u)*Y(x, z) - 4*Y(x, z)*X(wp_, w)
            - 2*X(xp_, u)*Y(z, y) - 2*X(wp_, u)*Y(z, y) + 2*X(wp_, u)*Y(x, z))
    # x<->y swap with d<->e relabel: Y(a,b) -> Y(b,a) after coordinate swap
    def swap(e):
        e2 = e.subs({x: y, y: x}, simultaneous=True)
        # relabel d<->e: Y(a,b)->Y(b,a)
        return e2.replace(Y, lambda a, b: sp.Function("Ytmp")(b, a)).replace(
            sp.Function("Ytmp"), lambda a, b: Y(a, b))
    sym = sp.expand(expr + swap(expr))
    print("4rho 4-Wilson-line bracket + (x<->y image) =", sym, "-> PASS" if sym == 0 else "-> FAIL")

# ---------- (2) z-integral identity ----------
def check_zintegral():
    # numeric in d=2 with IR cutoff R; traceless part must converge to -pi (r^m r^m'/r^2 - delta/2)
    y = (mp.mpf('0.17'), mp.mpf('-0.4'))
    yp = (mp.mpf('0.9'), mp.mpf('0.33'))
    r = (yp[0]-y[0], yp[1]-y[1])
    r2 = r[0]**2 + r[1]**2
    R = mp.mpf('60.0')
    def integrand(m, mprime):
        def f(t, s):  # polar around midpoint
            zx = t*mp.cos(s); zy = t*mp.sin(s)
            a = (y[0]-zx, y[1]-zy); b = (yp[0]-zx, yp[1]-zy)
            a2 = a[0]**2+a[1]**2; b2 = b[0]**2+b[1]**2
            return t*a[m]*b[mprime]/(a2*b2)
        return mp.quad(f, [0, R], [0, 2*PI], maxdegree=6)
    J = [[integrand(m, n) for n in range(2)] for m in range(2)]
    tr = J[0][0]+J[1][1]
    P = [[J[m][n] - (tr/2 if m == n else 0) for n in range(2)] for m in range(2)]
    tgt = [[-PI*(r[m]*r[n]/r2 - (mp.mpf(1)/2 if m == n else 0)) for n in range(2)] for m in range(2)]
    print("z-integral traceless part numeric:", [[mp.nstr(P[m][n], 6) for n in range(2)] for m in range(2)])
    print("                        claimed :", [[mp.nstr(tgt[m][n], 6) for n in range(2)] for m in range(2)])
    # trace: claimed = -(pi)[2/eps+...] ~ divergent log: compare numeric trace with -2pi log(R...) form:
    # J_tr ~ -2 pi log( r / (2R) ) ??? just print trace vs pi*log(4R^2/r^2) etc. for info
    print("trace numeric:", mp.nstr(tr, 8), " pi*log(4R^2/(r^2*e^{2gamma}))? ",
          mp.nstr(PI*mp.log(4*R**2/(r2*mp.exp(2*mp.euler))), 8), " pi*log(R^2/r2):",
          mp.nstr(PI*mp.log(R**2/r2), 8))

# ---------- (3) Group I Row II N2 split ----------
def check_g1row2_split():
    k1, k2, p1, p2, K, V, Lam = sp.symbols("k1 k2 p1 p2 Kp V Lam", positive=True)
    kk = sp.Matrix([k1, k2]); pp = sp.Matrix([p1, p2])
    k2s = k1**2+k2**2; p2s = p1**2+p2**2
    kdp = k1*p1+k2*p2
    kmp2 = (k1-p1)**2+(k2-p2)**2
    i = 0  # check component i=0 (and 1 by symmetry)
    def bracket_eq58(i):
        ki = [k1, k2][i]; pi = [p1, p2][i]
        return ( (-sp.log(V/Lam)+sp.log(V/K))*(k2s*pi - kdp*ki)/kmp2
                 - 2*ki*sp.log(V/Lam)
                 - 2*ki*sp.log(K*p2s/(V*k2s+K*p2s))
                 + ((2*p2s-2*kdp+k2s)*ki + k2s*pi)/kmp2*sp.log(1+(V/K-1)*k2s/kmp2)
                 - (p2s*ki+k2s*pi)/kmp2*(sp.log(V/Lam)+sp.log(1-K/V))
                 - ki*sp.log(V/K) )
    def R_eq79(i):
        ki = [k1, k2][i]; pi = [p1, p2][i]
        return ( -2*ki*sp.log(p2s/(k2s+K*p2s/V))
                 + ((2*p2s-2*kdp+k2s)*ki + k2s*pi)/kmp2*sp.log(k2s/kmp2 + (K/V)*(1-k2s/kmp2))
                 - (p2s*ki+k2s*pi)/kmp2*sp.log(1-K/V) )
    for i in (0, 1):
        delta = sp.simplify(sp.expand_log(bracket_eq58(i) - R_eq79(i), force=True))
        # delta should be  LJ * X_i  with LJ = log(K/Lam): check no V dependence:
        dV = sp.simplify(sp.diff(delta, V))
        coef = sp.simplify(delta / sp.log(K/Lam))
        print(f"component {i}: d(delta)/dV = {dV}")
        print(f"component {i}: delta/log(K/Lam) = {sp.simplify(coef)}")

if __name__ == "__main__":
    check_4rho_2W()
    check_4rho_4W()
    check_g1row2_split()
    check_zintegral()
