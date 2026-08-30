"""Full assembly  W = W_UV + W_fin(companion) + W_reg  and its numerical check."""
import numpy as np, sys
sys.path.insert(0, '/home/user/correlations/zint/src')
from masters import Kin, G_scalar, Psi, Nten, Tten, Phi, A7, A8, A9
from scipy.special import k0

def coeffs(xi):
    xb = 1-xi
    return dict(CUV=xi*xb + xi/xb + xb/xi, c1=2*xi*xb, c3=xi**2-2.5*xi+2-1/xi,
                c4=xb**2/(2*xi), c5=-1/(2*xi), c6=0.5, c7=-xi/(2*xb), c8=xi/2, c9=-0.5)

def W_parts(s, q, P, xi, m):
    xb = 1-xi; s2 = float(s@s); M2 = xi/xb*s2; D0 = s2 + M2; K = Kin(s, q, M2)
    c = coeffs(xi); sP = float(s@P); qn = float(np.hypot(*q))
    ps, N, T, ph = Psi(K), Nten(K), Tten(K), Phi(K)
    PsP, PsN, sPN, ssP = P@ps, P@N@s, s@N@P, s@ps
    S7 = sP*(A7(K)      - (ph - ssP))              / M2
    S8 = (A8(K, P) - (PsN - s2*PsP))               / M2
    S9 = (A9(K, P) - (sPN - sP*ssP))               / M2
    W_uv   = c['CUV']*sP*(2*np.pi/D0)*k0(m*qn)
    W_comp = c['CUV']*sP*(2*ssP - ph)/D0
    W_reg  = (c['c1']*(P@T@s) + c['c3']*PsP + c['c4']/s2*PsN + c['c5']*sP*ssP/s2
              + c['c6']*sP*ph/s2 + c['c7']*S7 + c['c8']*S8 + c['c9']*S9)
    return dict(UV=W_uv, companion=W_comp, reg=W_reg, total=W_uv+W_comp+W_reg,
                G=G_scalar(Kin(s, q, M2, m**2)), Gid=(2*np.pi/D0)*k0(m*qn)+(2*ssP-ph)/D0,
                pieces=dict(T=c['c1']*(P@T@s), Psi=c['c3']*PsP, N=c['c4']/s2*PsN,
                            sPsi=c['c5']*sP*ssP/s2, Phi=c['c6']*sP*ph/s2,
                            S7=c['c7']*S7, S8=c['c8']*S8, S9=c['c9']*S9))

def T_fin(r, s, P, xi):
    c = coeffs(xi); s2 = s@s
    r2 = np.sum(r**2, -1); u = r - s; u2 = np.sum(u**2, -1)
    Pr, sr, sP = r@P, r@s, s@P
    return (c['c1']*(Pr*sr/r2**2 - sP/(2*r2)) + c['c3']*Pr/r2 + c['c4']*Pr*sr/(s2*r2)
            + c['c5']*sP*sr/(s2*r2) + c['c6']*sP/s2
            + c['c7']*sP*np.sum(r*u, -1)/(r2*u2) + c['c8']*Pr*np.sum(s*u, -1)/(r2*u2)
            + c['c9']*sr*np.sum(P*u, -1)/(r2*u2))

if __name__ == "__main__":
    from quad2d import integrate
    s = np.array([0.7, -0.4]); q = np.array([0.9, 0.6]); P = np.array([1.3, 0.5])
    xi = 0.3; m = 1e-3
    R = W_parts(s, q, P, xi, m)
    print("G alpha-rep (m=1e-3) :", R['G'])
    print("G via identity       :", R['Gid'], " rel", abs(R['G']-R['Gid'])/abs(R['Gid']))
    xb = 1-xi; M2 = xi/xb*(s@s); sp_ = np.hypot(*s)
    Dq = lambda r: np.sum((r-s)**2, -1) + M2
    num = integrate(lambda r: T_fin(r, s, P, xi)/Dq(r), q, split=(sp_,), ntheta=4096, nlev=18)
    print("W_fin+W_reg assembled:", R['companion']+R['reg'])
    print("W_fin+W_reg quadrature:", num, " rel", abs(R['companion']+R['reg']-num)/abs(num))
    for k, v in R['pieces'].items():
        print(f"   {k:5s} {v:+.6f}")
