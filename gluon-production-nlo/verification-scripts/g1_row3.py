"""Group I, Row III: verify p+ integration, eq (85) [5 pieces] -> N1 (eq 89) + N2 (eq 92)
+ N3a (eq 94-113) + N3b (eq 115-121) with the arctan A (eq 123-126).
Color bracket common -> dropped. Factor i*g^4*f^{abd} and 1/(2pi)^3 dropped from both sides
(all terms carry it); i.e. we compare the real kernels multiplying i*f/(pi^4...).
PV poles handled by principal value; claimed logs evaluated as log| |.
"""
import random
import mpmath as mp
from common import rvec, sub, dot, sq, report, quad_pv

mp.mp.dps = 25
PI = mp.pi

def run(seed, K, V, Lam, verbose=False):
    rng = random.Random(seed)
    x, xp, y, z, w, wp = (rvec(rng) for _ in range(6))
    xz, xw, zw, yz = sub(x, z), sub(x, w), sub(z, w), sub(y, z)
    xpwp = sub(xp, wp)
    A, B, C, G = sq(xz), sq(xw), sq(zw), dot(xw, xz)
    Ep, F = sq(xpwp), sq(yz)
    K1 = dot(xpwp, yz)/(Ep*F)
    K2 = dot(xpwp, xw)*dot(xz, yz)/(Ep*B*A*F)

    P0 = K*B/A   # pole of (p A - K B)

    # ---------------- LHS ----------------
    c8 = 1/(8*PI**4*2*PI)   # if/8pi^4 * 1/2pi -> if/16pi^5 ; we keep 'i f' dropped
    c4 = 1/(4*PI**4*2*PI)

    def ig12(pp):
        return K1*(pp+K)/((K-pp)**2*(pp*A-K*B))

    p1 = c8*quad_pv(ig12, K+Lam, V, [P0])
    p2 = c8*quad_pv(ig12, Lam, K-Lam, [P0])

    def ig3(pp):
        return K2*(K*B+pp*A)/((pp*A-K*B)*pp*K)

    p3 = -c8*quad_pv(ig3, Lam, V, [P0])

    def N(pp):
        Vv = (K*xw[0]-pp*xz[0], K*xw[1]-pp*xz[1])
        return (dot(xpwp, zw)*dot(yz, Vv)
                - (K/pp)*dot(yz, zw)*dot(xpwp, Vv)
                - (K/(K-pp))*dot(xpwp, yz)*dot(zw, Vv))

    def ig4(pp):
        return N(pp)/(Ep*F*C*(K*B-pp*A)*(K-pp))

    p4 = -c4/K*quad_pv(ig4, Lam+K, V, [P0])

    def ig5(pp):
        Dq = (K*xw[0]-pp*xz[0])**2 + (K*xw[1]-pp*xz[1])**2
        return N(pp)*pp/(Ep*F*Dq*(K*B-pp*A)*(K-pp))

    p5 = -c4*quad_pv(ig5, Lam, K-Lam, [P0])

    lhs = p1+p2+p3+p4+p5

    # ---------------- RHS ----------------
    L = lambda u: mp.log(abs(u))
    i8 = 1/(8*PI**5)
    i16 = 1/(16*PI**5)

    # N1 (eq 89)
    n1 = i8*K1/(A-B)*(1/Lam + 1/(K-V))
    n2 = -i16/K*(mp.log(V/K-1)+mp.log(K/Lam))*K1*(A+B)/(A-B)**2
    n3 = i16/K*K1*(A+B)/(A-B)**2*L((V*A-K*B)/((Lam+K)*A-K*B))
    n4 = i8*K1/(A-B)*(1/Lam - 1/(K-Lam))
    n5 = -i16/K*mp.log(Lam/(K-Lam))*K1*(A+B)/(A-B)**2
    n6 = i16/K*K1*(A+B)/(A-B)**2*L(((K-Lam)*A-K*B)/(Lam*A-K*B))
    n7 = i16/K*(mp.log(V/K)+mp.log(K/Lam))*K2
    n8 = -i8/K*K2*L((K*B-V*A)/(K*B-Lam*A))
    N1 = n1+n2+n3+n4+n5+n6+n7+n8

    # N2 (eq 92): overall minus in front of the brace
    logX = L((B*K-A*V)/(B*K-A*(K+Lam)))
    m1 = i8*K1/(B-A)*(1/(V-K) - 1/Lam)
    m2 = -i8/K*mp.log(V/(Lam+K))*dot(xpwp, xw)*dot(yz, zw)/(Ep*F*C*B)
    # m3, m4 (red) cancel identically -> skip
    m5 = i8/K*(mp.log(K/Lam)+mp.log(V/K-1))*dot(xpwp, yz)*dot(zw, xz)*B/(Ep*F*C*(B-A)**2)
    m6 = -i8/K*(mp.log(K/Lam)+mp.log(V/K-1))*dot(xpwp, yz)*dot(zw, xw)*A/(Ep*F*C*(B-A)**2)
    m7 = i8/K*logX*dot(xpwp, zw)*dot(yz, xw)/(Ep*F*C*(B-A))
    m8 = -i8/K*logX*dot(xpwp, xw)*dot(yz, zw)*A/(Ep*F*C*(B-A)*B)
    m9 = i8/K*logX*dot(xpwp, yz)*dot(zw, xw)*A/(Ep*F*C*(B-A)**2)
    m10 = -i8/K*logX*dot(xpwp, zw)*dot(yz, xz)*B/(Ep*F*C*(B-A)*A)
    m11 = i8/K*logX*dot(xpwp, xz)*dot(yz, zw)/(Ep*F*C*(B-A))
    m12 = -i8/K*logX*dot(xpwp, yz)*dot(zw, xz)*B/(Ep*F*C*(B-A)**2)
    N2 = -(m1+m2+m5+m6+m7+m8+m9+m10+m11+m12)

    # N3a (eq 94-113)
    logY = L((B*K-A*(K-Lam))/(B*K-A*Lam))
    logKL = mp.log((K-Lam)/Lam)
    sqr = mp.sqrt(A*B-G**2)
    AA = (mp.atan(((K-Lam)*A-K*G)/(K*sqr)) - mp.atan((Lam*A-K*G)/(K*sqr)))/sqr

    a1 = i8*(1/Lam - 1/(K-Lam))*dot(xpwp, yz)/(Ep*F*(B-A))
    a2 = -i8/K*logKL*dot(xpwp, yz)/(Ep*F*C)
    a3 = -i8/K*logKL*dot(xpwp, yz)*dot(zw, xw)*A/(Ep*F*C*(B-A)**2)
    a4 = i8/K*logKL*dot(xpwp, yz)*dot(zw, xz)*B/(Ep*F*C*(B-A)**2)
    a5 = -i8/K*logY*dot(xpwp, yz)*dot(zw, xw)*A/(Ep*F*C*(B-A)**2)
    a6 = i8/K*logY*dot(xpwp, yz)*dot(zw, xz)*B/(Ep*F*C*(B-A)**2)
    a7 = -i8/K*logY*dot(xpwp, zw)*dot(yz, xw)/(Ep*F*C*(B-A))
    a8 = i8/K*logY*dot(xpwp, zw)*dot(yz, xz)*B/(Ep*F*C*(B-A)*A)
    a9 = i8/K*logY*dot(yz, zw)*dot(xpwp, xw)*A/(Ep*F*C*(B-A)*B)
    a10 = -i8/K*logY*dot(yz, zw)*dot(xpwp, xz)/(Ep*F*C*(B-A))
    a11 = i8/K*AA*dot(xpwp, zw)*dot(yz, xw)/(Ep*F*C)
    a12 = i8/K*AA*dot(yz, zw)*dot(xpwp, xz)/(Ep*F*C)
    a13 = -i8/K*AA*G*dot(xpwp, zw)*dot(yz, xz)/(Ep*F*C*A)
    N3a = a1+a2+a3+a4+a5+a6+a7+a8+a9+a10+a11+a12+a13

    # N3b (eq 115-121)
    DKL = (K*xw[0]-(K-Lam)*xz[0])**2 + (K*xw[1]-(K-Lam)*xz[1])**2
    DL = (K*xw[0]-Lam*xz[0])**2 + (K*xw[1]-Lam*xz[1])**2
    logZ = mp.log(DKL/DL)
    b1 = -i8/K*AA*G*dot(yz, zw)*dot(xpwp, xw)/(Ep*F*C*B)
    b2 = -i16/K*logZ*dot(xpwp, zw)*dot(yz, xz)/(Ep*F*C*A)
    b3 = i16/K*logZ*dot(yz, zw)*dot(xpwp, xw)/(Ep*F*C*B)
    b4 = -i16/K*logZ*dot(xpwp, yz)/(Ep*F*C)
    N3b = b1+b2+b3+b4

    rhs = N1+N2+N3a+N3b
    if verbose:
        print("  pieces: p12=", mp.nstr(p1+p2, 8), " vs n1..6=", mp.nstr(n1+n2+n3+n4+n5+n6, 8))
        print("  piece3 =", mp.nstr(p3, 8), " vs n7+n8 =", mp.nstr(n7+n8, 8),
              " vs 2*(n7+n8) =", mp.nstr(2*(n7+n8), 8))
        print("  piece4 =", mp.nstr(p4, 8), " vs N2 =", mp.nstr(N2, 8))
        print("  piece5 =", mp.nstr(p5, 8), " vs N3a+N3b =", mp.nstr(N3a+N3b, 8))
    return lhs, rhs

if __name__ == "__main__":
    for seed in (11, 12, 13):
        K, V, Lam = mp.mpf('1.3'), mp.mpf('47.0'), mp.mpf('1e-6')
        lhs, rhs = run(seed, K, V, Lam, verbose=True)
        report(f"G1 Row3 total seed={seed}", lhs, rhs, tol=mp.mpf('1e-4'))
