"""UV scan of the p+-INTEGRATED expression (the three-log form), region by region.
k+ = 1.  A = (x-z)^2, B = (x-w)^2, D = V A + B, DL = Lam A + B  (V = vee - k+)."""
import numpy as np
V, Lam, kp = 37.0, 1e-6, 1.0
vee = V + kp

def pieces(i, j, x, y, z, w, Vv):
    A=(x-z)@(x-z); B=(x-w)@(x-w); zw=z-w; zw2=zw@zw; xz,xw=x-z,x-w
    yz=y-z; K = Vv[i]*yz[j]/(yz@yz)
    d = 1.0 if i==j else 0.0
    D  = V*A + kp*B
    DL = Lam*A + kp*B
    br1 = ( xz[j]*zw[i]/(zw2*A) - xz[j]*xw[i]/(2*B*A)
          - xw[i]*zw[j]/(zw2*B) - xz[j]*xw[i]/(2*A*B) + d/(2*zw2) )
    t1 = K*br1*np.log(D/DL)
    t2 = -K*(d/(2*zw2))*np.log(vee/(Lam+kp))          # = the quoted 2nd term, in tensor form
    t3 = K*( xw[i]*zw[j]/(zw2*B) + xz[j]*xw[i]/(2*A*B) )*np.log(V/Lam)
    return t1, t2, t3

def F(x,y,z,w,Vv, which='all'):
    s = 0.0
    for i in range(2):
        for j in range(2):
            a,b,c = pieces(i,j,x,y,z,w,Vv)
            s += {'all':a+b+c,'t1':a,'t2':b,'t3':c,'t1+t2':a+b}[which]
    return s

x0=np.array([0.3,-0.2]); y0=np.array([-0.6,0.9]); z0=np.array([0.8,0.5])
w0=np.array([-0.4,-0.7]); Vv=np.array([0.4,0.7])
th=np.linspace(0,2*np.pi,8192,endpoint=False); u=np.stack([np.cos(th),np.sin(th)],-1)

def scan(name, mv, centre, whichs=('all',)):
    print(f"   {name}")
    for which in whichs:
        row=[]
        for rho in (1e-2,1e-3,1e-4,1e-5):
            vals=[]
            for uu in u:
                a=dict(x=x0,y=y0,z=z0,w=w0); a[mv]=centre+uu*rho
                vals.append(F(a['x'],a['y'],a['z'],a['w'],Vv,which))
            row.append(np.mean(vals)*rho**2)
        print(f"      {which:6s} rho^2<F> = " + "  ".join(f"{v:+.6e}" for v in row))

print("rho^2 x <integrand> : -> 0 means integrable, nonzero constant means LOG DIVERGENT")
print()
scan("x -> z", 'x', z0)
scan("x -> w", 'x', w0)
scan("y -> z", 'y', z0)
print()
print("   z -> w  -- SPLIT INTO PIECES, this is the interesting one:")
scan("z -> w", 'z', w0, whichs=('t1','t2','t1+t2','t3','all'))
