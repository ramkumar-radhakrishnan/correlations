"""Row 9: the large-|z| tail of B_2's one-rho piece, and what it turns into."""
import numpy as np
P,K = 0.37,1.0; S=P+K
def Wji(x,z,w):
    xz,xw,zw = x-z,x-w,z-w; s=lambda v: v@v
    return (np.eye(2)*(s(xz)-s(xw))/(2*s(zw))
            + (S/K)*(np.outer(xz,zw)/s(zw) - np.outer(xz,xw)/(2*s(xw)))
            + (S/P)*(np.outer(zw,xw)/s(zw) + np.outer(xz,xw)/(2*s(xz))))
Dd = lambda x,z,w: P*((x-z)@(x-z))+K*((x-w)@(x-w))
Kv = lambda a,b: (a-b)/((a-b)@(a-b))
x = np.array([0.31,-0.77]); w = np.array([1.13,0.42])
print("="*78); print("1.  LARGE-|z| LIMIT OF THE ONE-rho VERTEX"); print("="*78)
print("   claim:  W_{ji}(x,z,w)/[S D(x,z,w)]  -->  -(1/(2 p+ k+)) (x-w)^i (x-z)^j/[(x-w)^2 (x-z)^2]")
print()
print("   %-10s %16s %16s %10s" % ("|z|","exact [0,0]","limit [0,0]","ratio"))
for R in (1e1,1e2,1e3,1e4,1e5,1e6):
    ae = np.zeros((2,2)); al = np.zeros((2,2))
    for th in np.linspace(0,2*np.pi,64,endpoint=False):
        z = R*np.array([np.cos(th),np.sin(th)])
        ae += Wji(x,z,w)/(S*Dd(x,z,w))
        al += -np.outer(Kv(x,z),Kv(x,w))/(2*P*K)
    e,l = ae/64, al/64
    print("   %-10.0e %16.8e %16.8e %10.6f" % (R,e[0,0],l[0,0],e[0,0]/l[0,0]))
print()
print("   -> the gluon-exchange (one-rho) piece degenerates at large |z| into exactly")
print("      the TWO-rho piece's transverse structure, with weight -1/(2 p+ k+).")
print()
print("="*78); print("2.  HENCE A LOGARITHMIC *INFRARED* DIVERGENCE IN z"); print("="*78)
kt=np.array([0.7,-0.4]); xp=np.array([-0.52,0.19]); wp=np.array([-0.31,0.88]); y=np.array([0.62,0.41])
ph = np.exp(-1j*(kt@(wp-w)))
N1 = lambda z: ph*np.einsum('ji,j,i->',Wji(x,z,w)/S,Kv(y,z),Kv(xp,wp))/Dd(x,z,w)
N2 = lambda z: ph*np.einsum('ji,ji->',Wji(x,z,w),Wji(xp,z,wp))*(P*K/S**2)/(Dd(x,z,w)*Dd(xp,z,wp))
print("   %-10s %18s %18s" % ("|z|","N1: <.> x |z|^2","N2: <.> x |z|^2"))
for R in (1e2,1e3,1e4,1e5,1e6):
    a1=a2=0j
    for th in np.linspace(0,2*np.pi,2048,endpoint=False):
        z=R*np.array([np.cos(th),np.sin(th)]); a1+=N1(z); a2+=N2(z)
    print("   %-10.0e %18.9f %18.9f" % (R,(a1/2048*R**2).real,(a2/2048*R**2).real))
kk = Kv(x,w)@Kv(xp,wp)
print()
print("   predicted from part 1  (phase factored, real part):")
print("     N1 : -(1/(2 p+ k+)) (x-w).(x'-w')/[(x-w)^2 (x'-w')^2]  = %18.9f"
      % (-(1/(2*P*K))*kk*ph).real)
print("     N2 : +(1/(4 p+ k+)) (x-w).(x'-w')/[(x-w)^2 (x'-w')^2]  = %18.9f"
      % ((1/(4*P*K))*kk*ph).real)
