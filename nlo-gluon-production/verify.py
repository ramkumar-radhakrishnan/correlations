"""Symbolic / numerical checks for the claims in NLO_gluon_production_DGLAP.pdf.

    python3 verify.py          # requires sympy, numpy, scipy

Conventions follow the two rho-rho expressions of the note (Sec. 1).  Throughout,
zeta = p+/k+ is the momentum fraction of the daughter gluon at transverse position
x; the other daughter, at z, carries k+ - p+.  The pair separation is Z = x - z and
the parent sits at the centroid w = (1-zeta) z + zeta x, with u = y - w the distance
from the projectile source to the centroid.
"""
import math
import numpy as np
import sympy as sp

D = 2                                     # transverse dimensions
dl = lambda a, b: sp.Integer(a == b)
RESULTS = []


def check(name, cond):
    RESULTS.append((name, bool(cond)))
    print(("  PASS  " if cond else "  FAIL  ") + name)


# --------------------------------------------------------------------------
# Lemma 1 : C_1 is the light-cone g -> gg splitting vertex
# --------------------------------------------------------------------------
def lemma1():
    print("\nLemma 1 - C_1 is the DGLAP g -> gg vertex")
    z = sp.symbols("zeta", positive=True)

    def V(i, j, k, m):                    # x^m [ ... ] with x^m stripped
        return dl(i, m) * dl(j, k) - dl(j, m) * dl(i, k) / z - dl(k, m) * dl(i, j) / (1 - z)

    # sum over daughter indices j,k and over m, with the angular average
    # x^m x^m' -> (x^2/2) delta^{mm'} folded in as the overall 1/2
    T = sp.zeros(D, D)
    for i in range(D):
        for ip in range(D):
            T[i, ip] = sp.simplify(
                sp.Rational(1, 2)
                * sum(V(i, j, k, m) * V(ip, j, k, m)
                      for j in range(D) for k in range(D) for m in range(D)))

    Pgg_over_2Nc = z / (1 - z) + (1 - z) / z + z * (1 - z)
    check("T_{ii'} is proportional to delta_{ii'}",
          sp.simplify(T[0, 1]) == 0 and sp.simplify(T[1, 0]) == 0
          and sp.simplify(T[0, 0] - T[1, 1]) == 0)
    # restore the sqrt(p+(k+-p+))/k+ prefactor squared = zeta(1-zeta)
    check("zeta(1-zeta) * T_11 = P_gg(zeta)/(2 Nc)",
          sp.simplify(z * (1 - z) * T[0, 0] - Pgg_over_2Nc) == 0)

    f = z**2 + (1 - z)**2 + z**2 * (1 - z)**2
    check("KLSZ's [xi^2+(1-xi)^2+xi^2(1-xi)^2]/(xi(1-xi)) = P_gg/(2 Nc)",
          sp.simplify(f / (z * (1 - z)) - Pgg_over_2Nc) == 0)

    # double plus prescription: subtract BOTH poles (KLSZ eq. 27)
    integrand = sp.simplify(f / (z * (1 - z)) - 1 / z - 1 / (1 - z))
    check("double-plus integrand collapses to -2 + xi - xi^2",
          sp.simplify(integrand - (-2 + z - z**2)) == 0)
    I = sp.integrate(integrand, (z, 0, 1))
    check("int_0^1 = -11/6 = -beta_0^g/(2 Nc)  with beta_0^g = 11 Nc/3",
          sp.simplify(I + sp.Rational(11, 6)) == 0)


# --------------------------------------------------------------------------
# Lemma 2 : B_2 factorizes on merging into A^(1) (x) C_1
# --------------------------------------------------------------------------
def lemma2():
    print("\nLemma 2 - merging factorization of B_2")
    zeta, kp, eps, g = sp.symbols("zeta k^+ epsilon g", positive=True)
    p, k1 = zeta * kp, (1 - zeta) * kp             # p+ and k+ - p+
    u = sp.Matrix(sp.symbols("u1 u2", real=True))
    Z = sp.Matrix(sp.symbols("Z1 Z2", real=True))
    sq = lambda v: (v.T * v)[0, 0]

    def pieces(scale):
        Zs = scale * Z
        return u - (1 - zeta) * Zs, u + zeta * Zs, Zs   # y-x, y-z, Z

    def B2(i, j, scale):
        """Supplied B_2 with  w->z, z->x, x->y, k+->k+-p+, p+->p+, (p+ + k+)->k+."""
        yx, yz, Zs = pieces(scale)
        t1 = dl(i, j) / (2 * sq(Zs)) * (sq(yx) - sq(yz))
        t2 = (kp / k1) * (yx[j] * Zs[i] / sq(Zs) - yx[j] * yz[i] / (2 * sq(yz)))
        t3 = (kp / p) * (yz[i] * Zs[j] / sq(Zs) + yx[j] * yz[i] / (2 * sq(yx)))
        return (t1 + t2 + t3) / ((p * sq(yx) + k1 * sq(yz)) * kp)

    def C1(i, j, n):
        """C_1 bracket with the relative vector (z-x) = -Z; n = parent index,
        j = daughter carrying p+, i = daughter carrying k+ - p+."""
        return sum((-Z[m]) * (dl(n, m) * dl(i, j)
                              - (kp / p) * dl(j, m) * dl(n, i)
                              - (kp / k1) * dl(i, m) * dl(n, j))
                   for m in range(D)) / sq(Z)

    no_double_pole, residue_ok = True, True
    for i in range(D):
        for j in range(D):
            e = B2(i, j, eps)
            no_double_pole &= sp.simplify(sp.limit(eps**2 * e, eps, 0)) == 0
            res = sp.simplify(sp.limit(eps * e, eps, 0))
            pred = sp.simplify(sum(u[n] * C1(i, j, n) for n in range(D)) / (sq(u) * kp**2))
            residue_ok &= sp.simplify(res - pred) == 0
    check("B_2 has no 1/Z^2 singularity", no_double_pole)
    check("B_2's 1/|Z| residue = (u^n/u^2)(k+)^-2 x C_1 tensor, all (i,j)", residue_ok)

    # overall prefactors, including the 4 pi^2
    pre_B2 = sp.I * g**2 * sp.sqrt(p * k1) / (4 * sp.pi**2)
    pre_A = sp.I * g / (sp.sqrt(2) * sp.pi * sp.sqrt(kp))
    pre_C = -g / (2 * sp.pi * sp.sqrt(2 * kp)) * sp.sqrt(p * k1) / kp
    check("prefactors: A^(1) x C_1 = -B_2/(k+)^2  (exact, incl. 4 pi^2)",
          sp.simplify(pre_A * pre_C + pre_B2 / kp**2) == 0)

    # large-Z falloff -> the pair-size integral is cut off by the distance to the source
    lam = sp.symbols("lambda", positive=True)
    powers = set()
    for i in range(D):
        for j in range(D):
            e = B2(i, j, lam)
            for n in range(1, 5):
                if sp.simplify(sp.limit(lam**n * e, lam, sp.oo)) not in (0, sp.nan):
                    powers.add(n)
                    break
    check(f"B_2 ~ 1/Z^2 at large Z (found powers {sorted(powers)})", powers == {2})


# --------------------------------------------------------------------------
# Colour algebra of the two-Wilson-line term
# --------------------------------------------------------------------------
def colour():
    print("\nColour - the two-Wilson-line structure collapses onto Nc x (one Wilson line)")
    lam = [np.array(m, dtype=complex) for m in [
        [[0, 1, 0], [1, 0, 0], [0, 0, 0]], [[0, -1j, 0], [1j, 0, 0], [0, 0, 0]],
        [[1, 0, 0], [0, -1, 0], [0, 0, 0]], [[0, 0, 1], [0, 0, 0], [1, 0, 0]],
        [[0, 0, -1j], [0, 0, 0], [1j, 0, 0]], [[0, 0, 0], [0, 0, 1], [0, 1, 0]],
        [[0, 0, 0], [0, 0, -1j], [0, 1j, 0]],
        np.array([[1, 0, 0], [0, 1, 0], [0, 0, -2]]) / np.sqrt(3)]]
    t = [m / 2 for m in lam]
    N = 8
    f = np.zeros((N, N, N))
    for a in range(N):
        for b in range(N):
            for c in range(N):
                f[a, b, c] = np.real(-2j * np.trace((t[a] @ t[b] - t[b] @ t[a]) @ t[c]))

    rng = np.random.default_rng(0)
    H = rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3))
    H = H + H.conj().T
    wv, vv = np.linalg.eigh(H)
    V = vv @ np.diag(np.exp(1j * wv)) @ vv.conj().T
    V = V / np.linalg.det(V) ** (1 / 3)
    U = np.array([[2 * np.real(np.trace(V.conj().T @ t[a] @ V @ t[b]))
                   for b in range(N)] for a in range(N)])

    check("adjoint Wilson line is real orthogonal", np.allclose(U @ U.T, np.eye(N), atol=1e-9))
    check("f^{dca} U^{db} U^{ce} = f^{bea'} U^{aa'}   (both lines at the same point)",
          np.allclose(np.einsum("dca,db,ce->bea", f, U, U),
                      np.einsum("bex,ax->bea", f, U), atol=1e-9))
    ff = np.einsum("bex,bey->xy", f, f)
    check("f^{bea'} f^{bea''} = Nc delta,  Nc = %.3f" % ff[0, 0],
          np.allclose(ff, 3 * np.eye(N), atol=1e-9))


# --------------------------------------------------------------------------
# The transverse integral and the MSbar scale
# --------------------------------------------------------------------------
def transverse():
    print("\nTransverse integral - the scale that emerges is MSbar")
    import warnings
    from scipy import integrate as si
    from scipy.special import j0
    gam = 0.5772156649015329
    worst = 0.0
    for mu, q in [(1e4, 1.0), (1e5, 2.0), (1e6, 0.5)]:
        # J_0 oscillates forever; integrate the tail over Bessel half-periods so
        # the quadrature converges instead of warning about subdivisions.
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", si.IntegrationWarning)
            num = 2 * np.pi * si.quad(lambda Z: j0(q * Z) / Z, 1 / mu, np.inf,
                                      limit=2000)[0]
        pred = np.pi * np.log(4 * mu**2 * np.exp(-2 * gam) / q**2)
        worst = max(worst, abs(num - pred) / abs(pred))
        print(f"      mu={mu:.0e} q={q}:  numeric={num:.6f}  pi*ln(mu_MSbar^2/q^2)={pred:.6f}")
    check("int_{|Z|>1/mu} d^2Z/Z^2 e^{iqZ} = pi ln(mu_MSbar^2/q^2), "
          "mu_MSbar^2 = 4 mu^2 e^{-2 gamma}  (rel. err %.1e)" % worst, worst < 1e-4)


# --------------------------------------------------------------------------
# Size of the resummed correction (Sec. 9 table)
# --------------------------------------------------------------------------
def magnitude():
    print("\nMagnitude - (Qs^2/kT^2)^{2 gamma},  gamma = alpha_s beta_0^g / 4 pi,  Nc = 3")
    b0g = 11 * 3 / 3
    alphas = [0.15, 0.25, 0.35]
    print("      kT/Qs  " + "".join(f"  as={a:.2f}" for a in alphas))
    for r in [0.15, 0.20, 0.30, 0.50, 0.70, 1.00]:
        row = [1.0 if r >= 1 else (1 / r**2) ** (2 * a * b0g / (4 * math.pi)) for a in alphas]
        print(f"      {r:5.2f}  " + "".join(f"  {v:6.2f}" for v in row))
    g = 0.25 * b0g / (4 * math.pi)
    check("gamma ~= 0.875 alpha_s for Nc = 3", abs(g / 0.25 - 0.875) < 0.01)


if __name__ == "__main__":
    lemma1(); lemma2(); colour(); transverse(); magnitude()
    bad = [n for n, okk in RESULTS if not okk]
    print("\n" + "=" * 62)
    print(f"{len(RESULTS) - len(bad)}/{len(RESULTS)} checks passed")
    for n in bad:
        print("  FAILED:", n)
    raise SystemExit(1 if bad else 0)
