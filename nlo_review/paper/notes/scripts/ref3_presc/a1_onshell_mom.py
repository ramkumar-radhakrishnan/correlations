# Ref.3 / Major 5(a).  On-shell (delta-function) part of B_3 in momentum space, Eq. (B3mom) of app_wavefunctions.
# All vanishing energy denominators of B_3 are D = p+ k^2 - k+ p^2  (E_g(k)-E_g(p) times 2k+p+).
# Write the pole part of B_3^{ba}_{ji}(k -> p) as  i g^2 f^{abc} rho^c(k-p) n_{ji}(k,p) / D.
# With 1/(D +- i0) = P(1/D) -/+ i pi delta(D) the extra term is  -/+ i pi i g^2 f rho n delta(D).
# The on-shell term is anti-Hermitian (a unitary rotation inside the energy shell, allowed by Omega^dagger Omega = 1)
# iff n_{ji}(k,p) = n_{ij}(p,k) on the shell (see the note, Sec. 2).  This script checks that, branch by branch.
import numpy as np
rng = np.random.default_rng(1)

def n_tensor(kp, k, pp, p):
    """D * (pole part of B_3 / (i g^2 f rho)), as a 2x2 array [j,i] (j: outgoing p, i: incoming k), evaluated on the shell."""
    k2, p2 = k @ k, p @ p
    D = pp * k2 - kp * p2            # = 0 on shell; we multiply it out analytically below
    d = np.eye(2)
    out = np.zeros((2, 2))
    # (I) instantaneous, both branches
    out += np.sqrt(pp * kp) * (kp + pp) / (pp - kp) ** 2 * d
    # (II) commutator part of the two-valence-vertex diagrams:  -k^i p^j (p^2 k+ + k^2 p+)/(k^2 p^2 sqrt(k+p+)) / D
    out += -np.outer(p, k) * (p2 * kp + k2 * pp) / (k2 * p2 * np.sqrt(kp * pp))
    # (III) triple-gluon vertex: 2 sqrt(k+p+)(k-p)^l/(k+-p+)^2 (k^m p+ - k+ p^m)/D [tensor] [branch]   (D = k^2 p+ - p^2 k+)
    q = k - p
    vm = k * pp - kp * p
    T = np.zeros((2, 2))  # [j,i]
    for j in range(2):
        for i in range(2):
            s = 0.0
            for l in range(2):
                for m in range(2):
                    t = (d[j, m] * d[i, l] / pp - d[i, m] * d[j, l] / kp - d[l, m] * d[i, j] / (pp - kp))
                    s += q[l] * vm[m] * t
            T[j, i] = s
    if pp > kp:
        br = (pp - kp) / (q @ q)
    else:
        E3 = k2 / kp - q @ q / (kp - pp) - p2 / pp
        br = 1.0 / E3
    out += 2 * np.sqrt(kp * pp) / (kp - pp) ** 2 * T * br
    return out

worst = 0.0
for trial in range(2000):
    kp = rng.uniform(0.2, 3.0)
    k = rng.normal(size=2)
    pp = kp * rng.choice([rng.uniform(0.01, 0.95), rng.uniform(1.05, 30.0)])
    ang = rng.uniform(0, 2 * np.pi)
    p = np.sqrt(pp / kp) * np.linalg.norm(k) * np.array([np.cos(ang), np.sin(ang)])   # on shell: p^2/p+ = k^2/k+
    nkp = n_tensor(kp, k, pp, p)          # [j,i]
    npk = n_tensor(pp, p, kp, k)          # [i,j] : incoming p (index j), outgoing k (index i) -> array [i,j]
    diff = np.max(np.abs(nkp - npk.T)) / np.max(np.abs(nkp))
    worst = max(worst, diff)
print("max relative violation of n_ji(k,p) = n_ij(p,k) on shell over 2000 random points:", worst)
# the three pieces are separately symmetric, so the check above does not fix their relative sign; that is fixed by the
# cancellation of the double pole 1/(p+-k+)^2 between (I) and (III) (app_wavefunctions, after Eq. (B3mom)), checked here on shell:
kp, k = 1.0, np.array([0.7, -0.3])
for d in (1e-1, 1e-2, 1e-3, -1e-2, -1e-3):
    pp = kp * (1 + d); p = np.sqrt(pp / kp) * np.linalg.norm(k) * np.array([np.cos(1.1), np.sin(1.1)])
    print("p+/k+-1 = %+.0e :  (p+/k+-1)^2 * n =" % d, (n_tensor(kp, k, pp, p) * d ** 2).round(7).tolist(),
          "  |p+/k+-1| * n =", (n_tensor(kp, k, pp, p) * abs(d)).round(5).tolist())
# separate branches
kp, k = 1.0, np.array([0.7, -0.3])
for pp in (0.3, 4.0):
    p = np.sqrt(pp / kp) * np.linalg.norm(k) * np.array([np.cos(1.1), np.sin(1.1)])
    print("p+/k+ =", pp, " n(k,p) =\n", n_tensor(kp, k, pp, p), "\n n(p,k)^T =\n", n_tensor(pp, p, kp, k).T)
