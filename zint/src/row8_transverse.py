"""Row 8: transverse singular-region scan of the Wilson-line-reduced form."""
import numpy as np
rng = np.random.default_rng(11)

def K(a,b):        # (a-b)^m/(a-b)^2  vector
    d = a-b; return d/(d@d)

def integrand(x,xp,y,yp,z,kt,s):
    """the four kernels x the full delta-bracket, contracted.  s = p+/k+ ; k+=1."""
    Ky  = K(y,z); Kx  = K(x,y); Kyp = K(yp,z); Kxp = K(xp,yp)
    m,mp,kk,kp = Ky,Kyp,Kx,Kxp
    dperp = 2.0; P = s; Kp = 1.0; S = Kp+P
    t1 = dperp*P/S**2 * (kk@m)*(kp@mp)
    t2 = 2/Kp * ((m@kp)*(mp@kk) - (mp@kp)*(kk@m))
    t3 = (P/Kp**2 + 1/P) * (m@mp)*(kk@kp)
    ph = np.exp(-1j*(kt@(yp-y))*(1+s))
    return ph*(t1+t2+t3)

x  = np.array([0.31,-0.77]); xp = np.array([-0.52,0.19])
y0 = np.array([1.13, 0.42]); yp0= np.array([-0.31,0.88])
z0 = np.array([0.05,-0.23]); kt = np.array([0.7,-0.4]); s = 0.37

print("="*74)
print("A.  COLLAPSE SCAN:  integrand x measure,  n collapsing points -> rho^{2(n-1)} drho/rho")
print("="*74)
print("%-34s %12s %12s %12s" % ("region","rho=1e-2","rho=1e-4","rho=1e-6"))
def scan(name, mk, npts):
    vals=[]
    for rho in (1e-2,1e-4,1e-6):
        acc=0
        for th in np.linspace(0,2*np.pi,24,endpoint=False):
            acc += abs(integrand(*mk(rho,th),kt,s))
        acc/=24
        vals.append(acc*rho**(2*(npts-1)))          # x measure rho^{2(n-1)}
    print("%-34s %12.4e %12.4e %12.4e" % (name,*vals))

u=lambda th: np.array([np.cos(th),np.sin(th)])
scan("y' -> y        (2 pts)", lambda r,t:(x,xp,y0,y0+r*u(t),z0), 2)
scan("z  -> y        (2 pts)", lambda r,t:(x,xp,y0,yp0,y0+r*u(t)), 2)
scan("z  -> y'       (2 pts)", lambda r,t:(x,xp,y0,yp0,yp0+r*u(t)), 2)
scan("x  -> y        (2 pts)", lambda r,t:(y0+r*u(t),xp,y0,yp0,z0), 2)
scan("x' -> y'       (2 pts)", lambda r,t:(x,yp0+r*u(t),y0,yp0,z0), 2)
scan("y',z -> y      (3 pts)", lambda r,t:(x,xp,y0,y0+r*u(t),y0+r*u(t+2.1)), 3)
scan("x,y',z -> y    (4 pts)", lambda r,t:(y0+r*u(t+1.0),xp,y0,y0+r*u(t),y0+r*u(t+2.1)), 4)
scan("x,x',y',z -> y (5 pts)", lambda r,t:(y0+r*u(t+1.0),y0+r*u(t+4.0),y0,y0+r*u(t),y0+r*u(t+2.1)), 5)
print("\n(a constant column = log divergence;  growing to the right = power divergence;")
print(" shrinking to the right = convergent.)")

print()
print("="*74)
print("B.  LARGE-|z| (IR):  angular average of integrand x measure  rho^2 drho/rho")
print("="*74)
print("%-20s %14s" % ("|z|","avg * rho^2"))
for R in (1e1,1e2,1e3,1e4,1e5):
    acc=0
    for th in np.linspace(0,2*np.pi,400,endpoint=False):
        acc += integrand(x,xp,y0,yp0,R*u(th),kt,s)
    acc = acc/400*R**2
    print("%-20.0e %14.6f%+.6fj" % (R,acc.real,acc.imag))
print("(a constant = log divergence in the z integral at large |z|.)")

print()
print("="*74)
print("C.  behaviour of the z-INTEGRAL as y'->y  (the user's suspected region)")
print("="*74)
def Tnum(r, R=2000.0, nr=4000, nth=720):
    """T^{mm'} = int_{|z|<R} d^2z (y-z)^m (y'-z)^{m'}/((y-z)^2 (y'-z)^2), y=r/2, y'=-r/2"""
    y  =  r/2; yp = -r/2
    # log-radial grid centred on origin, with the two singular points integrable
    rr = np.geomspace(1e-7,R,nr); th = np.linspace(0,2*np.pi,nth,endpoint=False)
    RR,TH = np.meshgrid(rr,th,indexing='ij')
    Z = np.stack([RR*np.cos(TH),RR*np.sin(TH)],axis=-1)
    U = y-Z; V = yp-Z
    u2 = (U**2).sum(-1); v2=(V**2).sum(-1)
    F = U[...,:,None]*V[...,None,:]/(u2*v2)[...,None,None]
    w = (RR*np.gradient(rr)[:,None]*(2*np.pi/nth))[...,None,None]
    return (F*w).sum((0,1))
for rr in (1.0,0.3,0.1,0.03,0.01):
    r = np.array([rr,0.0])
    Tn = Tnum(r)
    pred = np.pi/2*(np.log(2000.0**2/rr**2)+1)*np.eye(2) - np.pi*np.outer(r,r)/rr**2
    print("|r|=%-6.2f  numeric T = [[%9.5f,%9.5f],[%9.5f,%9.5f]]" % (rr,*Tn.ravel()))
    print("             closed  T = [[%9.5f,%9.5f],[%9.5f,%9.5f]]" % (*pred.ravel(),))
