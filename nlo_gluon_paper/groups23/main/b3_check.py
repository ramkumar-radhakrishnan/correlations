"""Non-instantaneous entries of B_3 in momentum space (Calculation_main.tex, B3 momentum form, third term) compared
with the non-instantaneous one-gluon wave functions derived from m2, m3 (lcpt1.psi1 for p+ > k+, lcpt1.psi2 for p+ < k+).
If B_3 is right, both entries have the same ratio to the wave function."""
import numpy as np
import lcpt1 as L

F, PI, g = L.F, L.PI, L.g
I2 = np.eye(2)


def b3_noninst(kp, k, pp, p, v, which):
    """third term of B3^{ba}_{3ji}(k+,k;p+,p) as written; which = 'gt' (Theta(p+-k+-Lambda)) or 'lt' (Theta(k+-Lambda-p+))"""
    r = k - p
    pre = 2j * g**2 * np.sqrt(kp * pp) / (kp - pp)**2 * np.einsum('abc,c->ab', F, v)
    vec = (k * pp - kp * p) / ((k @ k) * pp - (p @ p) * kp)
    br = (np.einsum('jm,il->ijlm', I2, I2) / pp - np.einsum('im,jl->ijlm', I2, I2) / kp - np.einsum('lm,ij->ijlm', I2, I2) / (pp - kp))
    t = np.einsum('l,m,ijlm->ij', r, vec, br)
    if which == 'gt':
        f = (kp - pp) / (r @ r)
    else:
        f = 1 / ((k @ k) / kp - (r @ r) / (kp - pp) - (p @ p) / pp)
    return pre[:, :, None, None] * t[None, None] * f      # [a, b, i, j]


if __name__ == '__main__':
    rng = np.random.default_rng(4)
    for t in range(3):
        k, p, v = rng.normal(size=2), rng.normal(size=2), rng.normal(size=8)
        kp = 1.3
        pp = kp * rng.uniform(1.2, 3)
        a = L.psi1(kp, k, pp, p, v)          # rho^c(-(p-k)) = rho^c(k-p): v
        b = b3_noninst(kp, k, pp, p, v, 'gt')
        m = np.abs(a) > 1e-8 * np.abs(a).max()
        r1 = (b[m] / a[m])
        pp2 = kp * rng.uniform(0.1, 0.9)
        a2 = L.psi2(kp, k, pp2, p, v)        # rho^c(+(k-p)) = v
        b2 = b3_noninst(kp, k, pp2, p, v, 'lt')
        m2 = np.abs(a2) > 1e-8 * np.abs(a2).max()
        r2 = (b2[m2] / a2[m2])
        print('p+>k+: B3/psi = %.8f%+.8fi (spread %.1e)   p+<k+: B3/psi = %.8f%+.8fi (spread %.1e)   k+ %.2f p+ %.2f/%.2f' % (
            r1.mean().real, r1.mean().imag, np.abs(r1 - r1.mean()).max(), r2.mean().real, r2.mean().imag, np.abs(r2 - r2.mean()).max(), kp, pp, pp2))
