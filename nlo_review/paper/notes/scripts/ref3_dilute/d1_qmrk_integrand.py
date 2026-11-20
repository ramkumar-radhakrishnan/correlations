#!/usr/bin/env python3
"""
ref3_dilute / d1: dilute-dilute limit of the paper's two-gluon amplitude, Eq. (A2gamp),
against the QMRK gg production matrix element A_gluons of Bartels--Sabio Vera--Schwennsen
(hep-ph/0608154, Eq. eq:Agluons, D=4).

Dilute limit (derived in ref3_dilute.tex, Sec. 2):  U = 1 + u, u^{ab}(l) = f^{aeb} alpha^e(l),
keep O(u) and the part of the amplitude that is linear in rho after the anticommutator parts
have been reduced:  A_B^{(2)} + A_AA^{(1)} + A_AA^{(2)} = 1/4 [L_k, A_p] + 1/4 [L_p, A_k].
All pieces are written as   sum_t  C_t^{a b e d} F_t,ij(k,p,l)  alpha^e(l) rho^d(-q_a),
q_a = k + p + l  (projectile momentum), q_b = -l (target), observed gluon (k+, k), second (p+, p).

Momentum-space coefficients: Eq. (B2mom) of app_wavefunctions (1-rho part), Eq. (C1) of
sec_setup (coordinate form, Fourier transformed by hand in the note), Eq. (A1).

Output: ratio  R = q_a^4 * sum_{colour,pol} |M|^2 / (N_c^2 (N_c^2-1) A_gluons)  at random
kinematics; agreement means R is a kinematics-independent constant.
"""
import numpy as np, itertools, sys

rng = np.random.default_rng(12345)

# ---------------- colour ----------------
def fabc(N):
    if N == 2:
        f = np.zeros((3, 3, 3))
        for a, b, c in itertools.permutations(range(3)):
            f[a, b, c] = np.linalg.det(np.eye(3)[[a, b, c]])
        return f
    # SU(3) Gell-Mann
    l = np.zeros((8, 3, 3), complex)
    l[0][0, 1] = l[0][1, 0] = 1
    l[1][0, 1] = -1j; l[1][1, 0] = 1j
    l[2][0, 0] = 1; l[2][1, 1] = -1
    l[3][0, 2] = l[3][2, 0] = 1
    l[4][0, 2] = -1j; l[4][2, 0] = 1j
    l[5][1, 2] = l[5][2, 1] = 1
    l[6][1, 2] = -1j; l[6][2, 1] = 1j
    l[7] = np.diag([1, 1, -2]) / np.sqrt(3)
    t = l / 2
    f = np.zeros((8, 8, 8))
    for a in range(8):
        for b in range(8):
            comm = t[a] @ t[b] - t[b] @ t[a]
            for c in range(8):
                f[a, b, c] = (-2j * np.trace(comm @ t[c])).real
    return f

def colour_tensors(N):
    f = fabc(N)
    T = {}
    # index order of every tensor: [a, b, e, d]  (a observed, b second, e target, d projectile)
    T['1a'] = np.einsum('aec,cbd->abed', f, f)      # f^{aec} f^{cbd}
    T['1b'] = np.einsum('bec,acd->abed', f, f)      # f^{bec} f^{acd}
    T['1c'] = -np.einsum('abc,ced->abed', f, f)     # -f^{abc} f^{ced}
    T['2a'] = np.einsum('aec,cbd->abed', f, f)      # f^{aec} f^{cbd}
    T['2b'] = np.einsum('bec,cad->abed', f, f)      # f^{bec} f^{cad}
    T['3'] = np.einsum('cba,ced->abed', f, f)       # f^{cba} f^{ced}
    return T

# ---------------- kinematic coefficient functions (g = 1) ----------------
def beta(kp, k, pp, p):
    """1-rho part of B2^{ba}_{ji}(k+,k;p+,p) = f^{abc} rho^c(-p-k) beta_{ji}; returns beta[j,i]."""
    k2, p2 = k @ k, p @ p
    s = p + k; s2 = s @ s
    t1 = 1j * np.outer(p, k) * np.sqrt(pp * kp) / (2 * (pp * k2 + kp * p2)) * (1 / (kp * p2) - 1 / (k2 * pp))
    w = p * kp - pp * k                  # w^m = p^m k^+ - p^+ k^m
    D2 = s2 * np.sqrt(pp * kp) * (pp + kp) * (p2 / pp + k2 / kp)
    # sum_{l,m} s^l w^m [d_lm d_ji/(p+ + k+) - d_im d_jl/k+ - d_jm d_il/p+]   -> matrix [j,i]
    br = (s @ w) * np.eye(2) / (pp + kp) - np.outer(s, w) / kp - np.outer(w, s) / pp
    #        [j,i]:  -d_im d_jl s^l w^m = - s_j w_i  -> outer(s,w)[j,i] = s_j w_i  ok
    #                -d_jm d_il s^l w^m = - w_j s_i  -> outer(w,s)[j,i] = w_j s_i  ok
    t2 = -1j * br / D2
    t3 = -1j * (pp - kp) * np.eye(2) / (2 * np.sqrt(pp * kp) * (pp + kp) ** 2 * (k2 / kp + p2 / pp))
    return t1 + t2 + t3

def amplitude_pieces(kp, k, pp, p, l):
    qa = k + p + l
    P = kp + pp
    F = {}
    F['1a'] = beta(kp, k + l, pp, p)
    F['1b'] = beta(kp, k, pp, p + l)
    F['1c'] = beta(kp, k, pp, p)
    pre = 1j / (2 * np.sqrt(kp * pp))
    kl = k + l; pl = p + l
    # matrices indexed [j,i]: (vector_i)(vector_j) -> outer(vec_j, vec_i)
    F['2a'] = pre * (np.outer(p, kl) / (kl @ kl * (p @ p)) - np.outer(p, k) / (k @ k * (p @ p)))
    F['2b'] = pre * (np.outer(pl, k) / (k @ k * (pl @ pl)) - np.outer(p, k) / (k @ k * (p @ p)))
    K = k + p
    V = K / (K @ K) - qa / (qa @ qa)
    kap = (pp * k - kp * p) / P
    kap2 = kap @ kap
    # sum_{m,k} kap^m Br_{mjik} V_k = (kap.V) d_ji - (P/p+) kap_j V_i - (P/k+) kap_i V_j
    M3 = (kap @ V) * np.eye(2) - (P / pp) * np.outer(kap, V) - (P / kp) * np.outer(V, kap)
    F['3'] = -1j * np.sqrt(pp * kp) / P ** 2 * M3 / kap2
    return F

def cgc_sq(kp, k, pp, p, l, T):
    F = amplitude_pieces(kp, k, pp, p, l)
    M = sum(np.einsum('abed,ji->abedji', T[t], F[t]) for t in F)
    return np.sum(np.abs(M) ** 2)

def cgc_sq_split(kp, k, pp, p, l, T, keys):
    F = amplitude_pieces(kp, k, pp, p, l)
    M = sum(np.einsum('abed,ji->abedji', T[t], F[t]) for t in keys)
    return np.sum(np.abs(M) ** 2)

# ---------------- BSV A_gluons, D = 4 ----------------
def A_gluons(qa, qb, k1, k2, x, D=4.0):
    De = qa + qb
    assert np.allclose(De, k1 + k2)
    La = (1 - x) * k1 - x * k2
    sh = La @ La / (x * (1 - x))
    Sig = sh + De @ De
    th = 2 * qa @ k1 - qa @ qa - k1 @ k1 / x
    uh = 2 * qa @ k2 - qa @ qa - k2 @ k2 / (1 - x)
    E = (2 * (2 * x - 1) * qa @ qa + 4 * La @ qa + (1 - 2 * x) / (x * (1 - x)) * La @ La
         - 2 * x * (1 - x) * ((2 * x - 1) * De @ De + 2 * La @ De) * qa @ qa / (x * (1 - x) * De @ De + La @ La))
    qa2, qb2, k12, k22 = qa @ qa, qb @ qb, k1 @ k1, k2 @ k2
    r = (1 - x) / x; rr = x / (1 - x)
    A = qa2 * qb2 * (-1 / (th * uh) + qa2 * qb2 / (4 * th * uh * k12 * k22)
                     - 0.25 * (r / (k22 * th) + rr / (k12 * uh)) + 1 / (4 * k12 * k22)
                     + 1 / Sig * (-1 / sh * (2 + (1 / th - 1 / uh) * (r * k12 - rr * k22))
                                  + 0.25 * (Sig / sh + 1) * (r / k22 + rr / k12)
                                  - qb2 / (4 * sh) * (1 / ((1 - x) * th) + 1 / (x * uh))
                                  - qa2 / (4 * sh) * ((1 + rr * k22 / k12) / th + (1 + r * k12 / k22) / uh)))
    A += (D - 2) / 4 * ((((k1 - qa) @ (k1 - qa)) * ((k2 - qa) @ (k2 - qa)) - k12 * k22) / (th * uh)) ** 2
    A -= (D - 2) / 4 * 0.25 * ((((k2 - qa) @ (k2 - qa)) - rr * k22) / uh + E / sh) * \
         ((((k1 - qa) @ (k1 - qa)) - r * k12) / th - E / sh)
    return A

def rvec(s=1.0):
    return rng.normal(size=2) * s

def main():
    out = []
    for N in (2, 3):
        T = colour_tensors(N)
        CA = N * N * (N * N - 1)
        out.append(f"=== SU({N}) ===")
        rat = []
        for trial in range(12):
            k, p, l = rvec(), rvec(), rvec()
            x = rng.uniform(0.05, 0.95)
            kp, pp = x, 1 - x           # k+ = x P+, p+ = (1-x) P+, P+ = 1
            qa = k + p + l; qb = -l
            S = cgc_sq(kp, k, pp, p, l, T)
            Ag = A_gluons(qa, qb, k, p, x)
            R = x * (1 - x) * (qa @ qa) ** 2 * S / (CA * Ag)   # x(1-x): dp+/p+ measure and 1/k+ of dN/d^3k
            rat.append(R)
            out.append(f"x={x:.3f}  |k|={np.linalg.norm(k):.3f} |p|={np.linalg.norm(p):.3f} |l|={np.linalg.norm(l):.3f}"
                       f"  CGC={S:.6e}  A_gl={Ag:.6e}  R={R:.12f}")
        rat = np.array(rat)
        out.append(f"R: mean {rat.mean():.12f}  max rel. spread {np.max(np.abs(rat/rat.mean()-1)):.2e}")
        # MRK limits (gluon k much slower than p: x -> 0) : A_gluons -> qa^2 qb^2/(k^2 p^2)
        for x in (1e-3, 1e-5):
            k, p, l = rvec(), rvec(), rvec()
            qa = k + p + l; qb = -l
            Ag = A_gluons(qa, qb, k, p, x)
            out.append(f"MRK x={x:g}: A_gluons/(qa^2 qb^2/(k^2 p^2)) = {Ag/((qa@qa)*(qb@qb)/((k@k)*(p@p))):.8f}")
            S = cgc_sq(x, k, 1 - x, p, l, T)
            out.append(f"           CGC x(1-x) q_a^4 S/(CA qa^2 qb^2/(k^2 p^2)) = {x*(1-x)*(qa@qa)**2*S/(CA*(qa@qa)*(qb@qb)/((k@k)*(p@p))):.8f}")
    # ---- sensitivity: drop pieces / change weights; agreement must be lost ----
    T = colour_tensors(3); CA = 9 * 8
    pts = [(rvec(), rvec(), rvec(), rng.uniform(0.1, 0.9)) for _ in range(6)]
    def Rset(keys, scale=None):
        rr = []
        for k, p, l, x in pts:
            F = amplitude_pieces(x, k, 1 - x, p, l)
            if scale:
                for kk, s_ in scale.items(): F[kk] = F[kk] * s_
            M = sum(np.einsum('abed,ji->abedji', T[t], F[t]) for t in keys)
            S = np.sum(np.abs(M) ** 2); qa = k + p + l
            rr.append(x * (1 - x) * (qa @ qa) ** 2 * S / (CA * A_gluons(qa, -l, k, p, x)))
        return np.array(rr)
    allk = ['1a', '1b', '1c', '2a', '2b', '3']
    out.append("=== sensitivity (SU(3), 6 random points; R must be 1 at every point) ===")
    out.append("all pieces          : " + " ".join(f"{r:.6f}" for r in Rset(allk)))
    out.append("without commutators : " + " ".join(f"{r:.6f}" for r in Rset(['1a', '1b', '1c', '3'])))
    out.append("commutators x 2     : " + " ".join(f"{r:.6f}" for r in Rset(allk, {'2a': 2, '2b': 2})))
    out.append("without A_CA        : " + " ".join(f"{r:.6f}" for r in Rset(['1a', '1b', '1c', '2a', '2b'])))
    out.append("without B2[u rho]   : " + " ".join(f"{r:.6f}" for r in Rset(['1a', '1b', '2a', '2b', '3'])))
    # x -> 1 : observed gluon much faster than the second one (target-side MRK)
    for x in (1 - 1e-3, 1 - 1e-5):
        k, p, l = rvec(), rvec(), rvec(); qa = k + p + l; qb = -l
        Ag = A_gluons(qa, qb, k, p, x)
        S = cgc_sq(x, k, 1 - x, p, l, T)
        out.append(f"MRK x=1-{1-x:.0e}: A_gl/(qa^2qb^2/(k^2p^2)) = {Ag/((qa@qa)*(qb@qb)/((k@k)*(p@p))):.8f};"
                   f"  CGC same ratio = {x*(1-x)*(qa@qa)**2*S/(CA*(qa@qa)*(qb@qb)/((k@k)*(p@p))):.8f}")
    print("\n".join(out))

if __name__ == '__main__':
    main()
