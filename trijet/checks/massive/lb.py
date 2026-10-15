# Lepage-Brodsky light-cone spinors, Dirac representation. a^pm = a^0 pm a^3.
import sympy as sp
I=sp.I; R=sp.Rational
Z2=sp.zeros(2,2); I2=sp.eye(2)
s1=sp.Matrix([[0,1],[1,0]]); s2=sp.Matrix([[0,-I],[I,0]]); s3=sp.Matrix([[1,0],[0,-1]])
def blk(A,B,C,D): return sp.Matrix(sp.BlockMatrix([[A,B],[C,D]]))
g0=blk(I2,Z2,Z2,-I2); g=[blk(Z2,s,-s,Z2) for s in (s1,s2,s3)]
gp=g0+g[2]; gm=g0-g[2]; beta=g0; alpha=[g0*gi for gi in g]
chi={+1:sp.Matrix([1,0,1,0])/sp.sqrt(2), -1:sp.Matrix([0,1,0,-1])/sp.sqrt(2)}
def u(kp,kx,ky,m,lam): return (kp*sp.eye(4)+beta*m+alpha[0]*kx+alpha[1]*ky)*chi[lam]/sp.sqrt(kp)
def v(kp,kx,ky,m,lam): return (kp*sp.eye(4)-beta*m+alpha[0]*kx+alpha[1]*ky)*chi[-lam]/sp.sqrt(kp)
def bar(s): return s.H*g0
def slash_eps(r,kx,ky,kp):
    """LC-gauge (A^+=0) real Cartesian polarisation r: eps^+=0, eps_perp=e_r, eps^- = 2 e_r.k/k^+.
       slash(eps) = (1/2)(gamma^+ eps^- + gamma^- eps^+) - gamma_perp.eps_perp"""
    e=[1 if r==0 else 0, 1 if r==1 else 0]
    return R(1,2)*gp*(2*(e[0]*kx+e[1]*ky)/kp) - g[0]*e[0] - g[1]*e[1]
