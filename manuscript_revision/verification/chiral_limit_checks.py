# Chiral-limit checks of the physical interpretation added for referee point 3 (Figs. 3 and 4).
# Requires numpy, scipy. Run: python3 chiral_limit_checks.py
# Chiral-limit ground state: minimize E = int [ 1/2 (grad phi)^2 - H.grad phi ] on a square,
# discretized on links (natural Neumann BC built in). Quadratic -> sparse linear solve.
import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spl
def solve(Hx, Hz, L=10.0, n=201):
    x = np.linspace(-L/2, L/2, n); h = x[1]-x[0]
    X, Z = np.meshgrid(x, x, indexing='ij'); N = n*n; idx = np.arange(N).reshape(n, n)
    rows=[]; cols=[]; vals=[]; b = np.zeros(N)
    def add(a, c, w, Hl):   # link a->c, weight w, field component along link
        nonlocal b
        rows.extend([a,a,c,c]); cols.extend([a,c,a,c]); vals.extend([w,-w,-w,w])
        b[c] += w*h*Hl; b[a] -= w*h*Hl
    for i in range(n-1):          # x-links
        for j in range(n):
            w = 1.0 if 0<j<n-1 else 0.5
            add(idx[i,j], idx[i+1,j], w, Hx(x[i]+h/2, x[j]))
    for i in range(n):            # z-links
        for j in range(n-1):
            w = 1.0 if 0<i<n-1 else 0.5
            add(idx[i,j], idx[i,j+1], w, Hz(x[i], x[j]+h/2))
    A = sp.csr_matrix((vals,(rows,cols)), shape=(N,N)) + 1e-10*sp.eye(N)
    phi = spl.spsolve(A.tocsc(), b).reshape(n, n); phi -= phi[n//2, n//2]
    gx, gz = np.gradient(phi, h, h)
    return X, Z, phi, gx, gz, h
def trap(f, h):
    w = np.ones(f.shape[0]); w[0]=w[-1]=0.5
    return h*h*np.einsum('i,j,ij->', w, w, f)

# --- Gaussian "vortex" field of Fig. 3: H = (-d_z chi, d_x chi), chi = 4 exp(-(x^2+z^2)/4)
chi = lambda x,z: 4*np.exp(-(x*x+z*z)/4)
Hx = lambda x,z: -(-2*z/4)*chi(x,z); Hz = lambda x,z: (-2*x/4)*chi(x,z)
X,Z,phi,gx,gz,h = solve(Hx,Hz)
HX, HZ = Hx(X,Z), Hz(X,Z); nB = HX*gx+HZ*gz
print("Gaussian: max|H| =", np.abs(np.hypot(HX,HZ)).max(), " max|grad phi| =", np.hypot(gx,gz).max())
print("  max|n.H| on boundary =", np.abs(Hx(5.0, np.linspace(-5,5,1001))).max())
print("  nB on +x axis (x=2.8):", nB[np.argmin(abs(X[:,0]-2.8)), 100], " on diagonal (2,2):", nB[np.argmin(abs(X[:,0]-2)), np.argmin(abs(X[:,0]-2))])
print("  max nB, min nB:", nB.max(), nB.min(), " uniform-field scale |H|^2 max:", (HX**2+HZ**2).max())
print("  int H.grad phi =", trap(nB,h), " int (grad phi)^2 =", trap(gx**2+gz**2,h), " int H^2 =", trap(HX**2+HZ**2,h))
r = np.hypot(X,Z); m = np.abs(nB)>0.8*np.abs(nB).max(); print("  radius of |nB| peaks ~", r[m].mean())
# --- Sharp domain wall of Fig. 4: H = (0, H0 sgn x), H0=5
H0=5.0
X,Z,phi,gx,gz,h = solve(lambda x,z: 0*x, lambda x,z: H0*np.sign(x))
nB = H0*np.sign(X)*gz
print("DW: min/max of nB/H0^2:", nB.min()/H0**2, nB.max()/H0**2, " <E>/E0 =", trap(gx**2+gz**2,h)/(H0**2*100))
print("  int H.grad phi / int (grad phi)^2 =", trap(nB,h)/trap(gx**2+gz**2,h))
j0 = 100  # z=0 row
for xv in [0.5,1,2,3,4,5]:
    i = np.argmin(abs(X[:,0]-xv)); print(f"  x={xv}: d_z phi(x,0)/H0 = {gz[i,j0]/H0:.4f}")
