#!/usr/bin/env python3
"""
ref3_dilute / d2:
 (a) polarisation dimension: the dilute two-gluon amplitude written as  M_ji = A delta_ji + N_ji
     (N built from 2D vectors), summed over d_perp transverse polarisations:
        sum |M|^2 = sum_{2D} |A delta + N|^2 + (d_perp - 2) |A|^2,
     compared with BSV A_gluons at D = d_perp + 2  (BSV: D = 4 + 2 eps, the (D-2)/4 terms).
 (b) convolution with model distributions phi_P(q) = q^2 e^{-q^2/Q_P^2}, phi_T(q) = q^2 e^{-q^2/Q_T^2}
     (colour-singlet-like: vanish as q^2 at q->0) at a few k, real two-gluon production restricted
     to 0.1 < x < 0.9 and |p| > 0.25 (away from the soft and collinear singularities),
     computed once with the paper's dilute amplitude and once with A_gluons.
"""
import sys, numpy as np
sys.path.insert(0, '/home/user/correlations/nlo_review/paper/notes/scripts/ref3_dilute')
import d1_qmrk_integrand as d1

def pieces_split(kp, k, pp, p, l):
    """same as d1.amplitude_pieces, but every piece returned as (A, N): M = A*delta + N."""
    out = {}
    def beta_split(kp, k, pp, p):
        k2, p2 = k @ k, p @ p
        s = p + k; s2 = s @ s
        N1 = 1j * np.outer(p, k) * np.sqrt(pp * kp) / (2 * (pp * k2 + kp * p2)) * (1 / (kp * p2) - 1 / (k2 * pp))
        w = p * kp - pp * k
        D2 = s2 * np.sqrt(pp * kp) * (pp + kp) * (p2 / pp + k2 / kp)
        A2 = -1j * (s @ w) / (pp + kp) / D2
        N2 = -1j * (-np.outer(s, w) / kp - np.outer(w, s) / pp) / D2
        A3 = -1j * (pp - kp) / (2 * np.sqrt(pp * kp) * (pp + kp) ** 2 * (k2 / kp + p2 / pp))
        return A2 + A3, N1 + N2
    out['1a'] = beta_split(kp, k + l, pp, p)
    out['1b'] = beta_split(kp, k, pp, p + l)
    out['1c'] = beta_split(kp, k, pp, p)
    F = d1.amplitude_pieces(kp, k, pp, p, l)
    out['2a'] = (0.0, F['2a']); out['2b'] = (0.0, F['2b'])
    qa = k + p + l; P = kp + pp; K = k + p
    V = K / (K @ K) - qa / (qa @ qa)
    kap = (pp * k - kp * p) / P; kap2 = kap @ kap
    c = -1j * np.sqrt(pp * kp) / P ** 2 / kap2
    out['3'] = (c * (kap @ V), c * (-(P / pp) * np.outer(kap, V) - (P / kp) * np.outer(V, kap)))
    return out

def cgc_sq_d(kp, k, pp, p, l, T, dperp):
    Fs = pieces_split(kp, k, pp, p, l)
    A = sum(T[t][..., None, None] * (Fs[t][0] * np.eye(2) + Fs[t][1]) for t in Fs)   # [a,b,e,d,j,i]
    Aonly = sum(T[t] * Fs[t][0] for t in Fs)
    return np.sum(np.abs(A) ** 2) + (dperp - 2) * np.sum(np.abs(Aonly) ** 2)

def main():
    rng = np.random.default_rng(7)
    out = []
    T = d1.colour_tensors(3); CA = 72
    out.append("(a) d_perp dependence: R(d) = x(1-x) q_a^4 sum|M|^2_d / (CA * A_gluons(D=d+2))")
    for trial in range(5):
        k, p, l = rng.normal(size=2), rng.normal(size=2), rng.normal(size=2)
        x = rng.uniform(0.1, 0.9); qa = k + p + l
        row = []
        for dperp in (2.0, 1.7, 3.0, 6.0):
            S = cgc_sq_d(x, k, 1 - x, p, l, T, dperp)
            row.append(x * (1 - x) * (qa @ qa) ** 2 * S / (CA * d1.A_gluons(qa, -l, k, p, x, D=dperp + 2)))
        out.append(f"  x={x:.3f}: " + "  ".join(f"d={d}: {r:.14f}" for d, r in zip((2.0, 1.7, 3.0, 6.0), row)))
    # (b) convolution
    out.append("(b) restricted real two-gluon production, model phi's, Q_P=1, Q_T=1.5 (Monte Carlo, same points both ways)")
    QP, QT = 1.0, 1.5
    phiP = lambda q: (q @ q) * np.exp(-(q @ q) / QP ** 2)
    phiT = lambda q: (q @ q) * np.exp(-(q @ q) / QT ** 2)
    nmc = 40000
    for kk in (0.7, 1.5, 3.0):
        k = np.array([kk, 0.0])
        rng2 = np.random.default_rng(int(kk * 100))
        sA = sB = 0.0
        # importance sampling: q_a ~ gaussian width 1.2, p ~ gaussian width 2 (cut |p|>0.25), x uniform in (0.1,0.9)
        for _ in range(nmc):
            qa = rng2.normal(size=2) * 1.2; wq = 1.0 / (np.exp(-(qa @ qa) / (2 * 1.44)) / (2 * np.pi * 1.44))
            p = rng2.normal(size=2) * 2.0; wp = 1.0 / (np.exp(-(p @ p) / 8.0) / (8 * np.pi))
            if p @ p < 0.0625: continue
            x = rng2.uniform(0.1, 0.9); wx = 0.8
            qb = k + p - qa; l = -qb
            meas = wq * wp * wx / (x * (1 - x))           # d^2q_a d^2p dx/(x(1-x))  [= dp+/p+ at fixed k+]
            ff = phiP(qa) * phiT(qb) / ((qa @ qa) * (qb @ qb))
            S = d1.cgc_sq(x, k, 1 - x, p, l, T)
            sA += meas * ff * x * (1 - x) * (qa @ qa) ** 2 * S / CA
            sB += meas * ff * d1.A_gluons(qa, qb, k, p, x)
        out.append(f"  k={kk}: paper (dilute amplitude) {sA/nmc:.10e}   QMRK (BSV A_gluons) {sB/nmc:.10e}   rel.diff {abs(sA/sB-1):.1e}")
    out.append("  [common factor N (abar_s/pi)/(2 pi)^2 omitted on both sides; MC error ~1% is common to both]")
    print("\n".join(out))

if __name__ == '__main__':
    main()
