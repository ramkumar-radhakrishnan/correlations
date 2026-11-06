"""Check the NEW Other_terms file, Group I Row III:
    [p+ masters (eq 85, PV)]  ==  L_J*jhat + L_B*bhat + FIN(file)   as Lambda->0,
with jhat/bhat extracted numerically (two Lambda's / consistency across V handled by
transcribing bhat analytically is avoided: we instead verify
    D(Lam) := LHS(Lam,V) - FIN - jhat*log(K/Lam) - bhat*log(V/K) -> 0,
where jhat = [LHS(Lam1)-LHS(Lam2)]/log(Lam2/Lam1)  (pure -logLam coefficient)
      bhat = transcribed analytically from the printed N-term logs (verified
             independently against the LHS V-dependence by running two V's).
FIN transcribed from the uploaded 2aa8deab file, lines 63-111, with BOTH readings:
  plain log[...] (as uploaded, breaks where arguments < 0) and log|...| (PV).
"""
import random
import mpmath as mp
from common import rvec, sub, dot, sq, report
import g1_row3

mp.mp.dps = 25
PI = mp.pi

def kernels(x, xp, y, z, w, wp):
    xz, xw, zw, yz = sub(x, z), sub(x, w), sub(z, w), sub(y, z)
    xpwp = sub(xp, wp)
    A, B, C, G = sq(xz), sq(xw), sq(zw), dot(xw, xz)
    Ep, F = sq(xpwp), sq(yz)
    k = dict(A=A, B=B, C=C, G=G, Ep=Ep, F=F,
             Y=dot(xpwp, yz)/(Ep*F)*(A+B)/(A-B)**2,
             K1=dot(xpwp, yz)/(Ep*F),
             K2=dot(xpwp, xw)*dot(xz, yz)/(Ep*B*A*F),
             X1=dot(xpwp, yz)*dot(zw, xz)*B/(Ep*F*C*(B-A)**2),
             X2=dot(xpwp, yz)*dot(zw, xw)*A/(Ep*F*C*(B-A)**2),
             Xt=dot(xpwp, xw)*dot(yz, zw)/(Ep*F*C*B),
             k7=dot(xpwp, zw)*dot(yz, xw)/(Ep*F*C*(B-A)),
             k8=dot(xpwp, xw)*dot(yz, zw)*A/(Ep*F*C*(B-A)*B),
             k10=dot(xpwp, zw)*dot(yz, xz)*B/(Ep*F*C*(B-A)*A),
             k11=dot(xpwp, xz)*dot(yz, zw)/(Ep*F*C*(B-A)),
             orph=dot(xpwp, yz)/(Ep*F*C),
             q7=dot(xpwp, zw)*dot(yz, xw)/(Ep*F*C),       # A0 / N3a kernels
             q8=dot(xpwp, zw)*dot(yz, xz)*B/(Ep*F*C*A),
             q9=dot(yz, zw)*dot(xpwp, xw)*A/(Ep*F*C*B),
             q10=dot(yz, zw)*dot(xpwp, xz)/(Ep*F*C),
             r11=dot(xpwp, zw)*dot(yz, xw)/(Ep*F*C),
             r12=dot(yz, zw)*dot(xpwp, xz)/(Ep*F*C),
             r13=dot(xpwp, zw)*dot(yz, xz)/(Ep*F*C*A),    # x G
             s1=dot(yz, zw)*dot(xpwp, xw)/(Ep*F*C*B),     # x G (N3b A0)
             b2k=dot(xpwp, zw)*dot(yz, xz)/(Ep*F*C*A),
             b3k=dot(yz, zw)*dot(xpwp, xw)/(Ep*F*C*B),
             b4k=dot(xpwp, yz)/(Ep*F*C))
    return k

def FIN_file(kk, K, V, absval):
    """Finite terms of the uploaded file, GI Row III (N1+N2+N3a+N3b)."""
    A, B, C, G = kk['A'], kk['B'], kk['C'], kk['G']
    L = (lambda u: mp.log(abs(u))) if absval else (lambda u: mp.log(u))
    i8, i16 = 1/(8*PI**5), 1/(16*PI**5)
    ell = mp.log(1 - K/V)
    # ---- N1 (4 terms, lines 64-67)
    N1 = (-i16/K*ell*kk['Y']
          + i16/K*kk['Y']*L((A - (K/V)*B)/(A - B))
          + i16/K*kk['Y']*L((B - A)/B)
          - i8/K*kk['K2']*L((A - (K/V)*B)/B))
    # ---- N2 (8 terms, lines 70-76), overall minus in front of brace
    lX0 = L(((K/V)*B - A)/(B - A))
    N2 = -(i8/K*ell*kk['X1'] - i8/K*ell*kk['X2']
           + i8/K*lX0*kk['k7'] - i8/K*lX0*kk['k8'] + i8/K*lX0*kk['X2']
           - i8/K*lX0*kk['k10'] + i8/K*lX0*kk['k11'] - i8/K*lX0*kk['X1'])
    # ---- N3a (6 logW + 3 A0 terms, lines 80-100)
    lW = L((B - A)/B)
    sqr = mp.sqrt(A*B - G**2)
    A0 = (mp.atan((A - G)/sqr) + mp.atan(G/sqr))/sqr
    N3a = (-i8/K*lW*kk['X2'] + i8/K*lW*kk['X1']
           - i8/K*lW*kk['k7m'] + i8/K*lW*kk['k8m']
           + i8/K*lW*kk['q9m'] - i8/K*lW*kk['q10m']
           + i8/K*A0*kk['r11'] + i8/K*A0*kk['r12'] - i8/K*A0*G*kk['r13'])
    # ---- N3b (1 A0 + 3 logZ0 terms, lines 104-111)
    lZ0 = mp.log(C/B)
    N3b = (-i8/K*A0*G*kk['s1']
           - i16/K*lZ0*kk['b2k'] + i16/K*lZ0*kk['b3k'] - i16/K*lZ0*kk['b4k'])
    return N1 + N2 + N3a + N3b

def run(seed, K, V, absval=True, verbose=False):
    rng = random.Random(seed)
    x, xp, y, z, w, wp = (rvec(rng) for _ in range(6))
    kk = build_kk(x, xp, y, z, w, wp)
    return run_kk(seed, K, V, kk, absval, verbose)

def build_kk(x, xp, y, z, w, wp):
    kk = kernels(x, xp, y, z, w, wp)
    # extra kernels for N3a logW terms (a7..a10 descendants with 1/(B-A) prefactors):
    xz, xw, zw, yz = sub(x, z), sub(x, w), sub(z, w), sub(y, z)
    xpwp = sub(xp, wp)
    A, B, C = kk['A'], kk['B'], kk['C']
    Ep, F = kk['Ep'], kk['F']
    kk['q9m'] = dot(yz, zw)*dot(xpwp, xw)*A/(Ep*F*C*(B-A)*B)
    kk['q10m'] = dot(yz, zw)*dot(xpwp, xz)/(Ep*F*C*(B-A))
    # file N3a a7,a8 descendants:
    kk['k7m'] = dot(xpwp, zw)*dot(yz, xw)/(Ep*F*C*(B-A))
    kk['k8m'] = dot(xpwp, zw)*dot(yz, xz)*B/(Ep*F*C*(B-A)*A)
    return kk

def run_kk(seed, K, V, kk, absval=True, verbose=False):
    Lam1, Lam2 = mp.mpf('1e-6')*K, mp.mpf('1e-7')*K
    lhs1 = g1_row3.run(seed, K, V, Lam1)[0]
    lhs2 = g1_row3.run(seed, K, V, Lam2)[0]
    jhat = (lhs1 - lhs2)/mp.log(Lam2/Lam1)     # coefficient of log(K/Lam)
    # analytic jhat for cross-check: i16/K*K2 - i8/K*orph  (source attribution)
    jana = (1/(16*PI**5))/K*kk['K2'] - (1/(8*PI**5))/K*kk['orph']
    # bhat analytic: from printed N logs (see derivation in chat):
    i8, i16 = 1/(8*PI**5), 1/(16*PI**5)
    bana = (-i16*kk['Y'] + i16*kk['Y'] + i16*kk['K2'] - i8*kk['K2']   # n2,n3,n7,n8
            - (  -i8*kk['Xt']                                        # m2 (inside brace)
               + i8*kk['X1'] - i8*kk['X2']                           # m5,m6
               + i8*kk['k7'] - i8*kk['k8'] + i8*kk['X2']
               - i8*kk['k10'] + i8*kk['k11'] - i8*kk['X1']))/K       # m7..m12, overall -
    fin = FIN_file(kk, K, V, absval)
    resid = lhs1 - jhat*mp.log(K/Lam1) - bana*mp.log(V/K) - fin
    if verbose:
        print(f"  jhat num={mp.nstr(jhat,8)}  analytic={mp.nstr(jana,8)}  "
              f"rel={mp.nstr(abs(jhat-jana)/max(abs(jhat),mp.mpf('1e-30')),3)}")
    scale = max(abs(lhs1), abs(fin), mp.mpf('1e-30'))
    return resid, scale, jhat, jana

if __name__ == "__main__":
    K = mp.mpf('1.3')
    for seed in (11, 12, 13, 14, 15):
        for V in (mp.mpf('47.0'), mp.mpf('11.0')):
            r, s, jh, ja = run(seed, K, V, absval=True, verbose=(V == 47))
            print(f"seed {seed} V={V}: |resid|/scale = {mp.nstr(abs(r)/s, 4)}  (PV |log| reading)")
    # now demonstrate the plain-log reading: find where arguments go negative
    print("\n--- plain log[.] (as uploaded, no absolute values) ---")
    for seed in (11, 12, 13, 14, 15):
        try:
            r, s, jh, ja = run(seed, K, mp.mpf('47.0'), absval=False)
            tag = "complex!" if isinstance(r, mp.mpc) and abs(mp.im(r)) > 1e-20 else ""
            print(f"seed {seed}: |resid|/scale = {mp.nstr(abs(r)/s, 4)} {tag}")
        except Exception as e:
            print(f"seed {seed}: ERROR {e}")
