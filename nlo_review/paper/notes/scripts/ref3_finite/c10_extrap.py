# Richardson extrapolation in the window radius L (the window error is O(1/(kL)): the degree -2 tail of the
# amplitudes at |w| >> Z does not oscillate for w perpendicular to k).  Residual Z At - prediction at L -> infinity.
import sys, numpy as np
t = float(sys.argv[1]) if len(sys.argv) > 1 else 1.0
sys.argv = ['x', str(t)]
exec(open(__file__.replace('c10_extrap.py', 'c10_transition.py')).read().split("if __name__")[0])
phi = 0.6; nv = np.array([np.cos(phi), np.sin(phi)])
print('t = %.2f, phi = %.2f; residual max_ij |Z At - prediction|, extrapolated a + b/L + c/L^2 from L = (3, 4.5, 6) Z' % (t, phi))
for which in ('B', 'C', 'coh'):
    for Z in (10.0, 20.0, 40.0, 80.0):
        z = Z * nv; pred, cj = predicted(z, which)
        Ls = np.array([3.0, 4.5, 6.0]); vals = []
        for Lf in Ls:
            vals.append(Z * At_numeric(z, which, h=0.15, Lfac=Lf, sig=1.0) - pred)
        vals = np.array(vals)
        M = np.vstack([np.ones(3), 1 / (Ls * Z), 1 / (Ls * Z) ** 2]).T
        ext = np.linalg.solve(M, vals.reshape(3, -1)).reshape(3, 2, 2)[0]
        print('%-4s Z=%5.0f  |c_joint|=%.3f  residual at L=3Z: %.2e   extrapolated: %.2e   Z*extrapolated: %.3f' % (
            which, Z, np.abs(cj).max(), np.abs(vals[0]).max(), np.abs(ext).max(), Z * np.abs(ext).max()), flush=True)
