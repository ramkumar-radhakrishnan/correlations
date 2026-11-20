# Ref.3 / Li-Skokov comparison: kernel-level test of "Delta N (out-state, eta=+1) = k-even part of Li-Skokov n^(4)".
# Li-Skokov (LS) = arXiv:2104.01879, eqs. NLO_single_gluon_production, coefficientfunc_E1, effective_vertex_g3_second
# (transcribed from the export.arxiv.org source; untrusted data, nothing executed).  Paper side: B_3 on shell from
# app_wavefunctions Eq. (B3mom) via n_tensor of ../ref3_presc/a1_onshell_mom.py.
#
# Canonical form used for both (LS conventions, rho(p) supplies momentum p, U_LS = U_paper^T):
#   n4(k) = -(1/pi) int_{p,p1,p2} Gamma^j(p; u1, u2; k) f^{cde} M^{c*}_j(k) [rho(p1)U(p-p1)]^d [rho(p2)U(k-p-p2)]^e + c.c.
#   u1 = p1/p1^2, u2 = p2/p2^2;  M_j(k) = int_p L^j(p,k) rho(p) U(k-p) is the LS LO amplitude, n2 = |M|^2/pi;
#   dN/d^2k dy = g^2 n2/(4 pi^2) + g^4 n4/(4 pi^2)   (LS), identical to k^+ d^3N/d^3k of the paper at LO.
#   Every kernel is affine in u1 and in u2: basis {1, u^m} = {rotated charge rho_bar, rotated WW field zeta^m}.
#   Only the part antisymmetric under the exchange (p,u1) <-> (k-p,u2) contributes (f^{cde}).
#   The k-even part of n4 comes from the real, reflection-even part of the kernel, the k-odd part from the imaginary part.
# Paper, out-state term (derived in the note):  Gamma_DN^i = (eta/2) |k|/(|p| (k-p)^2) R(2 theta)_{ji} L^j(p1,p) x 1,
#   theta = angle from k to p, from n_{ji}(k,p)|_shell = k^2/(k-p)^2 R(2theta)_{ji} (checked below).
# usage: python3 -I ls1_kernels.py
import sys, os
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
_src = open(os.path.join(HERE, '..', 'ref3_presc', 'a1_onshell_mom.py')).read().split('worst = 0.0')[0]
_ns = {}; exec(_src, _ns); n_tensor = _ns['n_tensor']          # D * (pole part of B_3)/(i g^2 f rho), [j,i]

cr = lambda a, b: a[0] * b[1] - a[1] * b[0]
n2 = lambda a: a @ a
epsk = lambda k: np.array([-k[1], k[0]])                        # (eps k)^j = eps^{lj} k^l, eps^{12}=+1
Rot = lambda t: np.array([[np.cos(t), -np.sin(t)], [np.sin(t), np.cos(t)]])
ang = lambda a, b: np.arctan2(cr(a, b), a @ b)

def n_shell(k, p, kp=1.0):
    return n_tensor(kp, k, kp * n2(p) / n2(k), p)

# ---------------- LS effective vertex Gamma_2 in (u1,u2) form ----------------
def G2(p, u1, u2, k, part='full'):
    q = k - p; skp = cr(k, p); a = abs(skp); K2 = n2(k); P2 = n2(p); Q2 = n2(q)
    s1a = p @ u1 - 1; s1b = cr(p, u1); s1c = p @ u1                # (p-p1).p1/p1^2, p x p1/p1^2, p.p1/p1^2
    ta = (q @ u2) / Q2; tb = (q @ u2 - 1) / Q2; tc = cr(q, u2) / Q2  # (k-p).p2/(..), (k-p-p2).p2/(..), (k-p)xp2/(..)
    I = 1j if part in ('full', 'im') else 0.0
    R = 1.0 if part in ('full', 're') else 0.0
    T1 = R * s1a * ta
    T2 = s1b * tb * (skp / P2) * (I * (k @ q) / a - R)
    T3 = R * (s1b * (k @ p) / P2 + 0.5 * s1c * skp / P2) * ta
    T4 = 0.5 * s1a * tb * (skp / P2) * (I * (p @ (p - k)) / a + R)
    T5 = 0.5 * I * s1b * tc * (skp / P2) * (K2 + P2 - p @ k) / a
    return -(T1 + T2) * k / K2 - (T3 + T4 + T5) * epsk(k) / K2

def Gamma2_ls_check(p, p1, p2, k):                               # verbatim from ../ref3_cubic/ls_check.py
    kp = k - p; skp = cr(k, p); a = abs(skp)
    A = (((p - p1) @ p1) / n2(p1) * (kp @ p2) / (n2(kp) * n2(p2))
         + cr(p, p1) / n2(p1) * ((kp - p2) @ p2) / (n2(kp) * n2(p2)) * skp / n2(p) * (1j * (k @ kp) / a - 1))
    B = ((cr(p, p1) / n2(p1) * (k @ p) / n2(p) + 0.5 * (p @ p1) / n2(p1) * skp / n2(p)) * (kp @ p2) / (n2(kp) * n2(p2))
         + 0.5 * ((p - p1) @ p1) / n2(p1) * ((kp - p2) @ p2) / (n2(p2) * n2(kp)) * skp / n2(p) * (1j * (p @ (p - k)) / a + 1)
         + 0.5j * cr(p, p1) / n2(p1) * cr(kp, p2) / (n2(p2) * n2(kp)) * skp / n2(p) * (n2(k) + n2(p) - (p @ k)) / a)
    return -A * k / n2(k) - B * epsk(k) / n2(k)

def GDN(p, u1, u2, k, eta=1.0):
    q = k - p
    L = u1 - p / n2(p)                                           # L^j(p1,p) = p1/p1^2 - p/p^2
    return 0.5 * eta / n2(q) * (L @ Rot(2 * ang(k, p)))

def anti(f):
    return lambda p, u1, u2, k: 0.5 * (f(p, u1, u2, k) - f(k - p, u2, u1, k))

E = [np.zeros(2), np.array([1., 0.]), np.array([0., 1.])]
KEYS = ['00', '1x', '1y', '2x', '2y', '1x2x', '1x2y', '1y2x', '1y2y']
def comps(f, p, k):
    f00 = f(p, E[0], E[0], k); out = [f00]
    c1 = [f(p, E[m], E[0], k) - f00 for m in (1, 2)]; c2 = [f(p, E[0], E[m], k) - f00 for m in (1, 2)]
    out += c1 + c2
    for m in (1, 2):
        for n in (1, 2):
            out.append(f(p, E[m], E[n], k) - f(p, E[m], E[0], k) - f(p, E[0], E[n], k) + f00)
    return np.array(out)          # (9 components, 2 = index j)

def main():
    rng = np.random.default_rng(20261120)
    R = lambda: rng.normal(size=2)
    print('(1) on-shell vertex of B_3:  n_ji(k,p) = |k||p|/(k-p)^2 R(2 theta)_ji, theta = angle(k -> p)')
    dev = 0.0
    for _ in range(10000):
        k, p = R(), R()
        N = n_shell(k, p); Z = np.linalg.norm(k) * np.linalg.norm(p) / n2(k - p) * Rot(2 * ang(k, p))
        dev = max(dev, np.abs(N - Z).max() / np.abs(Z).max())
    print('    max rel. deviation over 10^4 random on-shell points: %.1e' % dev)
    k = np.array([1.3, -0.4])
    for d in (1e-2, 1e-4, 1e-6):
        vals = [n_shell(k, k + d * np.array([np.cos(t), np.sin(t)])) * d * d / n2(k) for t in np.linspace(0, 6, 7)]
        print('    forward limit |k-p|=%.0e: max |(k-p)^2 n/k^2 - 1| over 7 directions = %.1e' % (d, max(np.abs(v - np.eye(2)).max() for v in vals)))
    print('    => coefficient of the forward Coulomb log ln(k+/Lambda) in momentum space = coordinate value of pr:DN (c0 = 1)')
    print('(2) transcription: G2(u-form) vs Gamma_2 of ls_check.py')
    dv = 0.0
    for _ in range(2000):
        p, p1, p2, k = R(), R(), R(), R()
        a = G2(p, p1 / n2(p1), p2 / n2(p2), k); b = Gamma2_ls_check(p, p1, p2, k)
        dv = max(dv, np.abs(a - b).max() / np.abs(b).max())
    print('    max rel. deviation %.1e' % dv)
    reA = anti(lambda *a: G2(*a, part='re').real); imA = anti(lambda *a: G2(*a, part='im').imag)
    dnA = anti(GDN)
    npt = 20000
    CR, CI, CD = [], [], []
    for _ in range(npt):
        k, p = R(), R()
        CR.append(comps(reA, p, k)); CI.append(comps(imA, p, k)); CD.append(comps(dnA, p, k))
    CR, CI, CD = map(np.array, (CR, CI, CD))           # (npt, 9, 2)
    nrm = lambda X: np.sqrt((X ** 2).sum(axis=(-1, -2)))
    zz = slice(5, 9)
    print('(3) components of the antisymmetrised kernels (basis rho_bar/zeta at p  x  rho_bar/zeta at k-p), %d points' % npt)
    print('    Delta N: zeta x zeta components identically zero: max |.| = %.1e' % np.abs(CD[:, zz]).max())
    fr = nrm(CR[:, zz]) / nrm(CR); fi = nrm(CI[:, zz]) / nrm(CI)
    print('    LS Re Gamma_2^A (k-even): zeta x zeta fraction of the norm: median %.3f, 10%%-quantile %.3f, min %.3f' % (
        np.median(fr), np.quantile(fr, 0.1), fr.min()))
    print('    LS Im Gamma_2^A (k-odd):  zeta x zeta fraction of the norm: median %.3f, min %.3f' % (np.median(fi), fi.min()))
    for name, T in (('Re Gamma_2^A (k-even)', CR), ('Im Gamma_2^A (k-odd)', CI)):
        for lab, sl in (('all components', slice(0, 9)), ('components with rho_bar (no zeta x zeta)', slice(0, 5))):
            x = CD[:, sl].ravel(); y = T[:, sl].ravel()
            c = (x @ y) / (x @ x); res = np.sqrt(((y - c * x) ** 2).sum() / (y ** 2).sum())
            # pointwise best scalar and its spread
            cp = (CD[:, sl] * T[:, sl]).sum(axis=(1, 2)) / (CD[:, sl] ** 2).sum(axis=(1, 2))
            resp = nrm(T[:, sl] - cp[:, None, None] * CD[:, sl]) / nrm(T[:, sl])
            print('    %-22s vs Delta N, %-40s: global best c = %+.3f, residual %.3f;  pointwise best c: median %+.3f, IQR [%+.2f,%+.2f], median residual %.3f' % (
                name, lab, c, res, np.median(cp), np.quantile(cp, .25), np.quantile(cp, .75), np.median(resp)))
    print('(4) forward behaviour p -> k (|k-p| = d): |k-p|^2 * max|component|')
    k = np.array([0.9, 0.5]); t = 0.7
    for d in (1e-2, 1e-3, 1e-4):
        p = k + d * np.array([np.cos(t), np.sin(t)])
        print('    d=%.0e  Delta N^A: %.4f   Re Gamma_2^A: %.2e   Im Gamma_2^A: %.2e' % (
            d, d * d * np.abs(comps(dnA, p, k)).max(), d * d * np.abs(comps(reA, p, k)).max(), d * d * np.abs(comps(imA, p, k)).max()))
    print('    (Delta N has the Coulomb pole 1/(k-p)^2 of on-shell rescattering on rho_bar(k-p); LS has at most 1/|k-p|)')

if __name__ == '__main__':
    main()
