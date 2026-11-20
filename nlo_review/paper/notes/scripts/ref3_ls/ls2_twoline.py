# Ref.3 / Li-Skokov comparison, two-Wilson-line sector (no on-shell denominators involved).
# LS (arXiv:2104.01879): n4 ⊃ (1/pi) int E_1 f^{ade} rho^d(P-P1) rho^e(P1) rho^{e2}(p2) U^{ac}(k-P) U^{e2c}(-k-p2) + c.c.,
#   E_1 = Gamma_1^j L^j(p2,-k),  Gamma_1^j = (1/2) (P1 x P)/(|P-P1|^2 P1^2) eps^{ij} P^i / P^2   (real).
# Paper (app_finite Eq. app_IIb, Row I.II N_2b, structure f^{hbc} U^{ah}(w) {rho^b(x),rho^c(y)}), transformed to momentum
#   space in the conventions of ls1_kernels.py (derivation in the note):
#   n4 ⊃ -(1/pi) int_{q,p} Gamma_2b^i M^{a*}_i(k) f^{hbc} U^{ha}(k-q) rho^b(p) rho^c(q-p) + c.c.,
#   Gamma_2b^i = (i/(2 pi)) R^i(q,p)/q^2,  R^i real (app_IIb).    [LS: -(1/pi) int Gamma_1^i(P=q, P1=q-p) M^*_i ... ]
# A real kernel contributes to the k-even part, i x real to the k-odd part (both reflection-odd as amplitudes).
# Tests: (a) parities; (b) is R^i/q^2 proportional to Gamma_1^i (i.e. would a phase change alone reconcile them)?
# usage: python3 -I ls2_twoline.py
import numpy as np
cr = lambda a, b: a[0] * b[1] - a[1] * b[0]
n2 = lambda a: a @ a
eps = lambda v: np.array([-v[1], v[0]])        # (eps v)^i = eps^{ji} v^j ... here: eps^{ij} P^i -> components (-P2, P1)?
def epsP(P):                                    # vector with components eps^{ij} P^i, j = 1,2 (eps^{12} = +1)
    return np.array([-P[1], P[0]])

def R0(q, p):
    return (-2 * q * np.log(n2(p) / n2(q))
            + ((2 * n2(p) - 2 * (q @ p) + n2(q)) * q + n2(q) * p) / n2(q - p) * np.log(n2(q) / n2(q - p)))
def Rvec(q, p):
    return 0.5 * (R0(q, p) / n2(p) - R0(q, q - p) / n2(q - p))
def G2b(q, p): return Rvec(q, p) / n2(q)        # times i/(2 pi)
def G1(P, P1): return 0.5 * cr(P1, P) / (n2(P - P1) * n2(P1)) * epsP(P) / n2(P)

def main():
    rng = np.random.default_rng(7)
    R = lambda: rng.normal(size=2)
    a1 = a2 = r1 = r2 = 0.0
    lon_f = []; ratio = []
    for _ in range(20000):
        q, p = R(), R()
        a1 = max(a1, np.abs(G2b(q, q - p) + G2b(q, p)).max() / np.abs(G2b(q, p)).max())       # antisymmetry p <-> q-p
        a2 = max(a2, np.abs(G1(q, p) + G1(q, q - p)).max() / np.abs(G1(q, p)).max())
        r1 = max(r1, np.abs(G2b(-q, -p) + G2b(q, p)).max() / np.abs(G2b(q, p)).max())          # reflection-odd
        r2 = max(r2, np.abs(G1(-q, -(q - p)) + G1(q, q - p)).max() / np.abs(G1(q, q - p)).max())
        g = G2b(q, p); h = G1(q, q - p)
        qh = q / np.sqrt(n2(q)); th = epsP(q) / np.sqrt(n2(q))
        lon_f.append(abs(g @ qh) / np.sqrt(n2(g)))
        ratio.append((g @ th) / (h @ th) if abs(h @ th) > 1e-12 else np.nan)
    ratio = np.array(ratio); ratio = ratio[np.isfinite(ratio)]
    print('(a) antisymmetry under p <-> q-p: paper %.1e, LS %.1e;  reflection-odd: paper %.1e, LS %.1e' % (a1, a2, r1, r2))
    print('    => paper term = i x (real, reflection-odd amplitude kernel): k-ODD;  LS E_1 = real: k-EVEN, configuration by configuration')
    print('(b) Gamma_1 is purely transverse (along eps q); R/q^2 has a longitudinal (along q) part:')
    print('    |longitudinal|/|R/q^2|: median %.3f, max %.3f' % (np.median(lon_f), np.max(lon_f)))
    print('    transverse components, ratio (R/q^2)_T / Gamma_1,T : median %.3f, quartiles [%.3f, %.3f], 5-95%% [%.2f, %.2f]' % (
        np.median(ratio), np.quantile(ratio, .25), np.quantile(ratio, .75), np.quantile(ratio, .05), np.quantile(ratio, .95)))
    print('    => not proportional (the paper kernel contains log(p^2/q^2), log(q^2/(q-p)^2); Gamma_1 is rational): no constant phase')
    print('       factor maps one onto the other.')
    # the coefficient of log(Lambda/k+) in the paper's two-rho one-gluon component (app_wavefunctions Eq. Psigrhorho)
    # versus the classical O(g^3) WW field of LS, alpha_P^(3) ~ (delta - kk/k^2) p/((k-p)^2 p^2) [rho,rho]:
    dv = 0.0
    for _ in range(2000):
        k, p = R(), R()
        lcpt = (n2(k) * p - (k @ p) * k) / (n2(k - p) * n2(p) * n2(k))     # bracket of log(Lambda/k+), /(p^2 k^2)
        ww = (p - (k @ p) / n2(k) * k) / (n2(k - p) * n2(p))                 # (delta - kk/k^2) p /(|k-p|^2 p^2)
        dv = max(dv, np.abs(lcpt - ww).max() / np.abs(ww).max())
    print('(c) the transverse classical WW^(3) kernel of LS equals the coefficient of log(Lambda/k+) in the paper\'s two-rho one-gluon')
    print('    component (Eq. Psigrhorho), up to normalisation: max rel. dev. %.1e' % dv)

if __name__ == '__main__':
    main()
