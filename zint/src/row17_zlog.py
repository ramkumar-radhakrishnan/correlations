"""Row 17: the transverse log after the p+ integration -- it is the rapidity span."""
import numpy as np
s=lambda v:v@v
def taus(x,z,w):
    A,B,C=s(x-z),s(x-w),s(z-w)
    W=(z-w)/C-(x-w)/(2*B); V=(z-w)/C+(x-z)/(2*A)
    return (A-B)/(2*C), W, V, A, B
def Mall(A,B,Ap,Bp,k,Lam,Vv):
    al,be=k*B/A,k*Bp/Ap; P=Vv-k; L=np.log
    F=lambda p:(  # returns dict of antiderivatives
        {1:L((p+al)/(p+be))/(be-al),
         2:(-al*L(p+al)+be*L(p+be))/(k*(be-al)),
         3:k*(L(p)/(al*be)+L(p+al)/(al*(al-be))+L(p+be)/(be*(be-al))),
         4:k*(L(p+k)/((al-k)*(be-k))+L(p+al)/((k-al)*(be-al))+L(p+be)/((k-be)*(al-be))),
         5:(-k*L(p+k)/((al-k)*(be-k))-al*L(p+al)/((k-al)*(be-al))-be*L(p+be)/((k-be)*(al-be))),
         6:k*(k/((al-k)*(be-k))/(p+k)+(al*be-k*k)/((al-k)**2*(be-k)**2)*L(p+k)
             -al/((k-al)**2*(be-al))*L(p+al)-be/((k-be)**2*(al-be))*L(p+be))})
    hi,lo=F(P),F(Lam)
    return {n:(hi[n]-lo[n])/(A*Ap) for n in hi}
def G(x,xp,z,w,wp,k,Lam,Vv):
    t1,W,V,A,B=taus(x,z,w); t1p,Wp,Vp,Ap,Bp=taus(xp,z,wp)
    M=Mall(A,B,Ap,Bp,k,Lam,Vv)
    #  (a,b) -> (contraction, master index)
    return ( t1*t1p*2                 *M[6]      # T1T1'   h=pk/(p+k)^2
           + t1*((xp-z)@Wp)           *M[5]      # T1T2'   h=p/(p+k)
           + ((x-z)@W)*0              *0
           + t1*(Vp@(xp-wp))          *M[4]      # T1T3'   h=k/(p+k)
           + t1p*((x-z)@W)            *M[5]      # T2T1'
           + t1p*(V@(x-w))            *M[4]      # T3T1'
           + ((x-z)@(xp-z))*(W@Wp)    *M[2]      # T2T2'   h=p/k
           + ((x-z)@Vp)*(W@(xp-wp))   *M[1]      # T2T3'   h=1
           + (V@(xp-z))*((x-w)@Wp)    *M[1]      # T3T2'
           + (V@Vp)*((x-w)@(xp-wp))   *M[3] )    # T3T3'   h=k/p
rng=np.random.default_rng(3)
x,xp,w,wp=[rng.normal(size=2) for _ in range(4)]
k=1.0; B,Bp=s(x-w),s(xp-wp); Kww=((x-w)@(xp-wp))/(B*Bp)
print("="*78); print("LARGE-|z| COEFFICIENT OF THE p+-INTEGRATED INTEGRAND"); print("="*78)
print("   x=%s  x'=%s"%(np.round(x,4),np.round(xp,4)))
print("   w=%s  w'=%s"%(np.round(w,4),np.round(wp,4)))
print("   B=(x-w)^2=%.6f   B'=(x'-w')^2=%.6f   K_ww' = %.9f"%(B,Bp,Kww))
print()
for Lam,Vv in ((1e-6,1e5),(1e-9,1e5),(1e-6,1e8),(1e-3,1e3)):
    pred=Kww/(4*k)*np.log((Vv-k)/Lam)
    row=[]
    for R in (1e3,1e4,1e5,1e6):
        acc=0.0
        for t in np.linspace(0,2*np.pi,2048,endpoint=False):
            acc+=G(x,xp,R*np.array([np.cos(t),np.sin(t)]),w,wp,k,Lam,Vv)
        row.append(acc/2048*R*R)
    print("   Lam=%-8.0e V=%-8.0e  <G>|z|^2 : "%(Lam,Vv)
          +"  ".join("%12.7f"%v for v in row)+"   predicted %12.7f"%pred)
print()
print("   prediction  =  [ (x-w).(x'-w') / ((x-w)^2 (x'-w')^2) ] * log((q+-k+)/Lambda) / (4k+)")
print("   -> the ONLY surviving transverse log is int d^2z/|z|^2 with coefficient")
print("      = (LO transverse structure) x (rapidity span).  The log|z| pieces of")
print("      T2T2' and T3T3' cancel, and the (B+B')/(B'-B) log(B'/B) leftovers of")
print("      T2T2'+T3T3' cancel against T2T3'+T3T2'.")
