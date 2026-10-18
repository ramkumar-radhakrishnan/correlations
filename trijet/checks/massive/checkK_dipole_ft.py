# Check K: Fourier transform of the dipole-only terms.
#  (1) coordinates (x,y,z) <-> (c,R,Y) and (c,Rt,X): Jacobian, phases, dipole arguments
#  (2) momentum shifts produced by each dipole S(a,b)=int d^2l/(2pi)^2 e^{il.(a-b)} D(l)
#  (3) Phi_1 - Phi_2 identity
#  (4) numerical test of the shift formula with Gaussian test functions (independent of the algebra above)
import sympy as sp, numpy as np, itertools
z,x=sp.symbols('zeta xi',positive=True); e=1-z-x; zb=1-z
cx,cy,Rx,Ry,Yx,Yy=sp.symbols('c_x c_y R_x R_y Y_x Y_y',real=True)
vec=lambda a,b: sp.Matrix([a,b])
c,Rv,Yv=vec(cx,cy),vec(Rx,Ry),vec(Yx,Yy)
# ---------- antiquark channel: c = z x + e y + xi z, R = x - y', Y = y - z,  y' = (e y + xi z)/(1-z)
X_=c+zb*Rv; Yp=c-z*Rv; Y_=Yp+x/zb*Yv; Z_=Yp-e/zb*Yv
assert sp.simplify(z*X_+e*Y_+x*Z_-c)==sp.zeros(2,1) and sp.simplify(X_-(e*Y_+x*Z_)/zb-Rv)==sp.zeros(2,1) and sp.simplify(Y_-Z_-Yv)==sp.zeros(2,1)
J=sp.Matrix([[z,e,x],[1,-e/zb,-x/zb],[0,1,-1]]).det()
print("(1) antiquark channel: Jacobian (per transverse component) =",sp.simplify(J))
k1=vec(*sp.symbols('k1x k1y',real=True)); k2=vec(*sp.symbols('k2x k2y',real=True)); k3=vec(*sp.symbols('k3x k3y',real=True))
phase=sp.expand((k1.T*X_+k2.T*Y_+k3.T*Z_)[0])
K=k1+k2+k3; kt=zb*k1-z*(k2+k3); pt=(e*k3-x*k2)/zb
target=sp.expand((K.T*c+kt.T*Rv-pt.T*Yv)[0])
print("    k1.x+k2.y+k3.z = K.c + kt.R - pt.Y  with kt=(1-zeta)k1-zeta(k2+k3), pt=(eta k3-xi k2)/(1-zeta):", sp.simplify(phase-target)==0)
args_qb={'x-z':sp.simplify(X_-Z_),'z-y':sp.simplify(Z_-Y_),"x-y'":sp.simplify(X_-Yp)}
for k_,v_ in args_qb.items(): print(f"    {k_:5s} =", list(v_))
# ---------- quark channel: x' = (z x + xi z)/(z+xi), Rt = x' - y, X = x - z
Rtx,Rty,Xx,Xy=sp.symbols('Rt_x Rt_y X_x X_y',real=True); Rt,Xv=vec(Rtx,Rty),vec(Xx,Xy)
Yq=c-(1-e)*Rt; Xp=c+e*Rt; Xq=Xp+x/(1-e)*Xv; Zq=Xp-z/(1-e)*Xv
assert sp.simplify(z*Xq+e*Yq+x*Zq-c)==sp.zeros(2,1) and sp.simplify((z*Xq+x*Zq)/(z+x)-Xp)==sp.zeros(2,1)
assert sp.simplify(Xp-Yq-Rt)==sp.zeros(2,1) and sp.simplify(Xq-Zq-Xv)==sp.zeros(2,1)
Jq=sp.Matrix([[z,e,x],[sp.Rational(1)*z/(z+x),-1,x/(z+x)],[1,0,-1]]).det()
print("    quark channel: Jacobian =",sp.simplify(Jq))
phase_q=sp.expand((k1.T*Xq+k2.T*Yq+k3.T*Zq)[0])
ktq=e*(k1+k3)-(1-e)*k2; ptq=(z*k3-x*k1)/(1-e)
print("    k1.x+k2.y+k3.z = K.c + ktq.Rt - ptq.X  with ktq=eta(k1+k3)-(1-eta)k2, ptq=(zeta k3-xi k1)/(1-eta):",
      sp.simplify(phase_q-sp.expand((K.T*c+ktq.T*Rt-ptq.T*Xv)[0]))==0)
args_q={'x-z':sp.simplify(Xq-Zq),'z-y':sp.simplify(Zq-Yq),'y-z':sp.simplify(Yq-Zq),'z-x':sp.simplify(Zq-Xq),
        "x'-y":sp.simplify(Xp-Yq),"y-x'":sp.simplify(Yq-Xp)}
for k_,v_ in args_q.items(): print(f"    {k_:5s} =", list(v_))
# ---------- (2) momentum shifts:  int dR dY e^{-i kt.R + i pt.Y} M(R,Y) e^{i l.(alpha R + beta Y)} = (2pi)^2 Phi(kt - alpha l, pt + beta l)
print("(2) momentum shifts  [Phi(kt,pt) -> Phi(kt - alpha l, pt + beta l) for S(a,b) with a-b = alpha R + beta Y]")
def coef(expr,A,B):   # a-b = alpha*A + beta*B  (A,B vectors)
    al=sp.simplify(expr[0].diff(A[0])); be=sp.simplify(expr[0].diff(B[0]))
    assert sp.simplify(expr-al*A-be*B)==sp.zeros(2,1); return al,be
for k_,v_ in args_qb.items():
    al,be=coef(v_,Rv,Yv); print(f"    antiquark: S({k_}): kt -> kt - ({al}) l ,  pt -> pt + ({be}) l")
for k_,v_ in args_q.items():
    al,be=coef(v_,Rt,Xv); print(f"    quark    : S({k_}): ktq -> ktq - ({al}) l ,  ptq -> ptq + ({be}) l")
# ---------- (3) Phi_1 - Phi_2 = -xi eta [N/(zeta zb (p^2+om^2)) + I]/D
k2s,p2s,eps2,om2,Nn,In=sp.symbols('k2 p2 eps2 om2 N I')
D=z*zb**2*(p2s+om2)+x*e*(k2s+eps2)
Phi1=(zb*Nn/(k2s+eps2)-x*e*In)/D; Phi2=Nn/(z*zb*(k2s+eps2)*(p2s+om2))
print("(3) Phi_1-Phi_2 = -xi eta [N/(zeta(1-zeta)(p^2+om^2)) + I]/D :", sp.simplify(Phi1-Phi2+x*e*(Nn/(z*zb*(p2s+om2))+In)/D)==0)
# ---------- (4) Gaussian test of the shift formula (antiquark, S(x,z)S(z,y)), independent of (2)
rng=np.random.default_rng(1)
zv,xv=0.3,0.25; ev=1-zv-xv; cc=ev/(1-zv)         # x-z = R + c Y,  z-y = -Y
Mq=np.array([[1.1,0.3],[0.3,0.8]])                 # test amplitude M(R,Y)=exp(-(a R^2 + 2 g R.Y + b Y^2)) (rotation invariant)
a_,g_,b_=Mq[0,0],Mq[0,1],Mq[1,1]; Q1,Q2=1.3,0.9
kt_=np.array([0.4,-0.7]); pt_=np.array([0.9,0.2])
def gauss4(A,Jv):  # int d^4u exp(-u.A.u + J.u)
    return np.pi**2/np.sqrt(np.linalg.det(A))*np.exp(Jv@np.linalg.solve(A,Jv)/4)
I2=np.eye(2)
# LHS: coordinate integral, u=(R,Y)
A=np.block([[(a_+Q1**2/4)*I2,(g_+cc*Q1**2/4)*I2],[(g_+cc*Q1**2/4)*I2,(b_+cc*cc*Q1**2/4+Q2**2/4)*I2]])
lhs=gauss4(A,np.concatenate([-1j*kt_,1j*pt_]))
# RHS: (2pi)^2 int d2l1 d2l2/(2pi)^4 D(l1)D(l2) Phi(kt-l1, pt+c l1-l2), Phi = (2pi)^-2 int dRdY e^{-ik.R+ip.Y}M
A0=np.block([[a_*I2,g_*I2],[g_*I2,b_*I2]])
Phi=lambda k,p: gauss4(A0,np.concatenate([-1j*k,1j*p]))/(2*np.pi)**2
Dl=lambda l,Q: 4*np.pi/Q**2*np.exp(-l@l/Q**2)
xs,ws=np.polynomial.hermite.hermgauss(14)
rhs=0
for (i1,w1),(i2,w2),(j1,u1),(j2,u2) in itertools.product(zip(xs,ws),repeat=4):
    l1=Q1*np.array([i1,i2]); l2=Q2*np.array([j1,j2])   # D(l)=4pi/Q^2 e^{-l^2/Q^2}; substitute l=Q h
    rhs+=w1*w2*u1*u2*(4*np.pi)*(4*np.pi)*Phi(kt_-l1,pt_+cc*l1-l2)
rhs*= (2*np.pi)**2/(2*np.pi)**4
print(f"(4) Gaussian test, antiquark S(x,z)S(z,y):  coordinate integral = {lhs:.12f}   momentum convolution = {rhs:.12f}   rel.diff = {abs(lhs-rhs)/abs(lhs):.1e}")
