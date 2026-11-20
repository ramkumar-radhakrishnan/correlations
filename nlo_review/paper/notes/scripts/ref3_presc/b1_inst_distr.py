# Ref.3 / Major 5(b).  Instantaneous (p+-independent) terms of Group III, Rows III.III, III.V, III.VI.
# After the w' integral the p+ dependence is the phase only.  With t = p+/k+, r = y'-z, q = (1+t) k, the structure is
#     I = int_0^{T} dt  h(t),     h(t) = int d^2r e^{i q.r} K^m(r) [b(r) - b(0)],     T = vee/k+,
# with b(r) = K^j(x'-y') L(y') at y' = z + r (smooth model functions below).  This script
#  (1) gives the p+ distribution h(t) and the cumulative fraction of the finite part from t < t0  (where in p+ it comes from);
#  (2) shows that the large-t tail is ~ c/t^2 for smooth b, ~ t^(-1-alpha) for a Hoelder-alpha b (rough target);
#  (3) replaces the sharp cut t<T by a smooth f(t/T) for several f: the limit T->oo is f-independent, the approach is O(1/T);
#  (4) the subtraction term (the BFKL piece), int dt f(t/T)/(1+t) = ln T + c_f + O(1/T), c_f = int du/u [f(u)-theta(1-u)],
#      is the only f-dependence: a universal shift ln(vee) -> ln(vee) + c_f.
import numpy as np
from scipy.integrate import quad
from scipy.special import jv, gamma

# --- model b(r)-b(0) = (u.r) |r|^(alpha-1) exp(-r^2/(2 s^2))   (alpha=1: smooth, linear at r=0; alpha<1: Hoelder) ---
# h^m(q) = u^n int d^2r e^{iq.r} (r^m r^n / r^2) |r|^(alpha-1)... we use  K^m(r)(u.r) = u^n r^m r^n / r^2, so
# h^m = u^n [ pi delta^{mn} J0-moment - pi (2 qh^m qh^n - delta^{mn}) J2-moment ],  moments = int_0^oo dr r * g(r) * J_{0,2}(q r),
# with g(r) = |r|^(alpha-1) exp(-r^2/2s^2).
import mpmath as mp
def moments(q, alpha, s):
    """int_0^oo dr r^alpha e^{-a r^2} J_nu(q r), nu=0,2, a=1/(2 s^2): closed form (Gradshteyn 6.631.1)."""
    a = 1.0 / (2 * s * s)
    if alpha == 1.0:
        e = np.exp(-q * q / (4 * a))
        return e / (2 * a), 2.0 / (q * q) * (1 - e) - e / (2 * a)
    out = []
    mu = alpha + 1.0
    for nu in (0, 2):
        h = (nu + mu) / 2
        val = mp.gamma(h) * (q / 2) ** nu / (2 * a ** h * mp.gamma(nu + 1)) * mp.hyp1f1(h, nu + 1, -q * q / (4 * a))
        out.append(float(val))
    return tuple(out)

def h_par_perp(t, alpha=1.0, s=1.0, kabs=1.0, upar=1.0, uperp=0.5):
    """components of h along k (par) and perpendicular to k, for u = (upar, uperp) in the (khat, khat_perp) basis."""
    q = (1 + t) * kabs
    m0, m2 = moments(q, alpha, s)
    # 2 qh qh - 1 = diag(1,-1) in this basis
    hpar = upar * (np.pi * m0 - np.pi * m2)
    hperp = uperp * (np.pi * m0 + np.pi * m2)
    return hpar, hperp

def run(alpha, s, label):
    ts = np.concatenate([np.linspace(0, 20, 401)[:-1], np.geomspace(20, 4000, 200)])
    H = np.array([h_par_perp(t, alpha, s) for t in ts])
    # total by trapezoid on the grid + analytic tail from a fit c/t^(1+alpha) beyond the grid
    tot = np.trapezoid(H, ts, axis=0)
    tail_c = H[-1] * ts[-1] ** (1 + alpha)
    tot_inf = tot + tail_c * ts[-1] ** (-alpha) / alpha
    print(f"\n[{label}] alpha={alpha}, s={s}  (|k|=1, k+=1): finite part int_0^oo h dt = par {tot_inf[0]: .6f}, perp {tot_inf[1]: .6f}")
    for t0 in (1, 3, 10, 30, 100, 1000):
        i = np.searchsorted(ts, t0)
        frac = np.trapezoid(H[:i + 1], ts[:i + 1], axis=0) / tot_inf
        print(f"    fraction from p+ < {t0:5d} k+ :  par {frac[0]: .4f}   perp {frac[1]: .4f}")
    # large-t behaviour
    for t in (100.0, 400.0, 1600.0):
        hp = np.array(h_par_perp(t, alpha, s))
        print(f"    t={t:7.0f}:  t^(1+alpha) h = {hp[0]*t**(1+alpha): .5f} {hp[1]*t**(1+alpha): .5f}")
    return ts, H, tot_inf

ts, H, lim = run(1.0, 1.0, "smooth b, scale s=1/|k|")
run(1.0, 0.1, "smooth b, small scale s=0.1/|k| (p+ ~ k+/(|k| s))")
run(0.5, 1.0, "Hoelder-1/2 b (rough target)")

# (3) smooth cut-offs f(t/T)
fs = {"sharp": lambda u: (u < 1).astype(float), "exp(-u)": lambda u: np.exp(-u),
      "exp(-u^2)": lambda u: np.exp(-u * u), "1/(1+u^2)": lambda u: 1 / (1 + u * u)}
print("\n(3) smooth cut-off f(p+/vee) on the instantaneous term (alpha=1, s=1): int dt f(t/T) h(t)  [par component]")
tfine = np.concatenate([np.linspace(0, 20, 2001)[:-1], np.geomspace(20, 2e5, 3000)])
Hf = np.array([h_par_perp(t, 1.0, 1.0)[0] for t in tfine])
for T in (10.0, 100.0, 1000.0):
    row = []
    for name, f in fs.items():
        row.append(np.trapezoid(Hf * f(tfine / T), tfine) - lim[0])
    print(f"    T=vee/k+={T:6.0f}: deviation from the limit: " + "  ".join(f"{n}:{v: .2e}" for n, v in zip(fs, row)))

# (4) the subtraction (BFKL) piece: int_0^oo dt f(t/T)/(1+t) - ln T  -> c_f
print("\n(4) BFKL piece: int dt f(t/T)/(1+t) - ln T  vs  c_f = int du/u [f(u)-theta(1-u)]")
for name, f in fs.items():
    cf = quad(lambda u: (f(np.array(u)) - (u < 1)) / u, 0, 1, limit=200)[0] + quad(lambda u: f(np.array(u)) / u, 1, np.inf, limit=200)[0]
    vals = []
    for T in (1e2, 1e3, 1e4):
        v = quad(lambda t: f(np.array(t / T)) / (1 + t), 0, 50 * T, limit=2000, points=[T])[0] - np.log(T)
        vals.append(v)
    print(f"    {name:10s}: c_f={cf: .6f}   T=1e2,1e3,1e4: {vals[0]: .6f} {vals[1]: .6f} {vals[2]: .6f}")
print("    (c_f: sharp 0, exp(-u) -gamma_E, exp(-u^2) -gamma_E/2, 1/(1+u^2) 0)")
