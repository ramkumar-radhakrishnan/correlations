"""Group II, Row I, second part (the B2^dagger B2 'N2' piece):
verify the p+ integration step, eq.(56/64) -> eqs (2a)-(2f) [lines 68-83 of Group_II_corrected.tex].

Method: pick random 2D points x,x',z,w,w'; numerically integrate the p+ integrand of eq (56)
and compare with the sum of the claimed results 2a+...+2f at the same points.
Color bracket (common factor) dropped; g^4 and 1/(2pi)^3 dropped from both sides.
"""
import random
import mpmath as mp
from common import rvec, sub, dot, sq, report

mp.mp.dps = 25
PI = mp.pi

def run(seed, K, V, Lam):
    rng = random.Random(seed)
    x, xp, z, w, wp = (rvec(rng) for _ in range(5))
    xz, xw, zw = sub(x, z), sub(x, w), sub(z, w)
    xpz, xpwp, zwp = sub(xp, z), sub(xp, wp), sub(z, wp)
    A, B, C = sq(xz), sq(xw), sq(zw)
    Ap, Bp, Cp = sq(xpz), sq(xpwp), sq(zwp)

    # ---------- LHS: eq (56) ----------
    def M(pp):
        # returns 2x2 matrix M_ij for unprimed bracket
        out = [[mp.mpf(0)]*2 for _ in range(2)]
        f1 = (pp+K)/K
        f2 = (pp+K)/pp
        for i in range(2):
            for j in range(2):
                v = mp.mpf(0)
                if i == j:
                    v += (A-B)/(2*C)
                v += f1*(xz[j]*zw[i]/C - xz[j]*xw[i]/(2*B))
                v += f2*(xw[i]*zw[j]/C + xz[j]*xw[i]/(2*A))
                out[i][j] = v
        return out

    def Mp(pp):
        out = [[mp.mpf(0)]*2 for _ in range(2)]
        f1 = (pp+K)/K
        f2 = (pp+K)/pp
        for i in range(2):
            for j in range(2):
                v = mp.mpf(0)
                if i == j:
                    v += (Ap-Bp)/(2*Cp)
                v += f1*(xpz[j]*zwp[i]/Cp - xpz[j]*xpwp[i]/(2*Bp))
                v += f2*(xpwp[i]*zwp[j]/Cp + xpz[j]*xpwp[i]/(2*Ap))
                out[i][j] = v
        return out

    def integrand(pp):
        F = 1/(pp*A + K*B)
        Fp = 1/(pp*Ap + K*Bp)
        m, mpr = M(pp), Mp(pp)
        s = sum(m[i][j]*mpr[i][j] for i in range(2) for j in range(2))
        return (pp*K/(pp+K)**2) * F * Fp * s

    lhs = (1/(4*PI**4)) * (1/(2*PI)) * mp.quad(integrand, [Lam, V-K])

    # ---------- RHS: 2a..2f ----------
    LgA = mp.log(((V-K)*A + K*B)/(Lam*A + K*B))
    LgAp = mp.log(((V-K)*Ap + K*Bp)/(Lam*Ap + K*Bp))
    Lblue = mp.log(K/Lam) + mp.log(V/K - 1)
    LV = mp.log(V/(K+Lam))
    c8 = 1/(8*PI**5*K)
    d_perp = 2
    c32 = d_perp/(32*PI**5)

    # 2a
    Bil_a = ( dot(xw, zwp)*dot(xpz, zw)/(C*Cp)
            - dot(xw, xpwp)*dot(zw, xpz)/(2*C*Bp)
            + dot(xz, xpz)*dot(xw, zwp)/(2*A*Cp)
            - dot(xpz, xz)*dot(xpwp, xw)/(4*A*Bp)
            + dot(xpwp, zw)*dot(xz, zwp)/(C*Cp)
            + dot(xpwp, zw)*dot(xz, xpz)/(2*C*Ap)
            - dot(xz, zwp)*dot(xpwp, xw)/(2*B*Cp)
            - dot(xpwp, xw)*dot(xpz, xz)/(4*B*Ap) )
    r2a = c8/(Bp*A - B*Ap) * Bil_a * (
        mp.log(((V-K)*A + K*B)/(K*Bp + (V-K)*Ap)) - mp.log((Lam*A + K*B)/(K*Bp + Lam*Ap)))

    # 2b
    Bil_b = ( dot(xpwp, xw)*dot(zw, zwp)/(C*Cp)
            + dot(xpwp, xw)*dot(xpz, zw)/(2*C*Ap)
            + dot(xw, xpwp)*dot(xz, zwp)/(2*A*Cp)
            + dot(xz, xpz)*dot(xpwp, xw)/(4*A*Ap) )
    den_b = Bp*A*B - B**2*Ap
    r2b = ( c8*Lblue/den_b * Bil_b * (A - B*Ap/Bp)
          - c8*A/den_b * Bil_b * LgA
          + c8*Ap/(Bp**2*A - B*Ap*Bp) * Bil_b * LgAp )

    # 2c
    Bil_c = ( dot(xpz, xz)*dot(zw, zwp)/(C*Cp)
            - dot(xpz, xz)*dot(xpwp, zw)/(2*C*Bp)
            - dot(xw, zwp)*dot(xz, xpz)/(2*B*Cp)
            + dot(xz, xpz)*dot(xpwp, xw)/(4*B*Bp) )
    r2c = ( c8*(Bp/Ap)/(Bp*A - B*Ap) * Bil_c * LgAp
          - c8*(B/A)/(Bp*A - B*Ap) * Bil_c * LgA )

    # 2d
    Bil_d = ( dot(xw, zw)*(Ap-Bp)/(2*C*Cp)
            + dot(xz, xw)*(Ap-Bp)/(4*A*Cp)
            + dot(xpwp, zwp)*(A-B)/(2*C*Cp)
            + dot(xpz, xpwp)*(A-B)/(4*C*Ap) )
    r2d = ( -c8*LV/((A-B)*(Bp-Ap)) * Bil_d
          + c8*A/((A-B)*(Bp*A - B*Ap)) * Bil_d * LgA
          + c8*Ap/((Bp-Ap)*(Bp*A - B*Ap)) * Bil_d * LgAp )

    # 2e
    Bil_e = ( dot(xz, zw)*(Ap-Bp)/(2*C*Cp)
            - dot(xz, xw)*(Ap-Bp)/(4*B*Cp)
            + dot(xpz, zwp)*(A-B)/(2*C*Cp)
            - dot(xpz, xpwp)*(A-B)/(4*C*Bp) )
    r2e = ( c8*LV/((A-B)*(Bp-Ap)) * Bil_e
          - c8*B/((A-B)*(Bp*A - B*Ap)) * Bil_e * LgA
          - c8*Bp/((Bp-Ap)*(Bp*A - B*Ap)) * Bil_e * mp.log((K*Bp + (V-K)*Ap)/(K*Bp + Lam*Ap)) )

    # 2f
    r2f = ( c32*LV/K * (Bp*B - A*Ap)/(C*Cp*(A-B)*(Ap-Bp))
          - c32/K * (A*B)/(C*Cp*(A-B)) * (Ap-Bp)/(Bp*A - B*Ap) * LgA
          + c32*(1/V - 1/(K+Lam)) / (C*Cp)
          + c32/K * (Bp*Ap)/(C*Cp) * (A-B)/((Ap-Bp)*(Bp*A - B*Ap)) * LgAp )

    rhs = r2a + r2b + r2c + r2d + r2e + r2f
    return lhs, rhs

if __name__ == "__main__":
    ok = True
    for seed in (1, 2, 3, 4):
        K, V, Lam = mp.mpf('1.3'), mp.mpf('47.0'), mp.mpf('1e-7')
        lhs, rhs = run(seed, K, V, Lam)
        ok &= report(f"G2 Row1 part2 seed={seed}", lhs, rhs, tol=mp.mpf('1e-8'))
    print("ALL PASS" if ok else "SOME FAIL")
