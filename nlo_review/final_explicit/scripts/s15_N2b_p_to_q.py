import mpmath as mp
mp.mp.dps=40
def R(q,p):
    q2=q[0]**2+q[1]**2; p2=p[0]**2+p[1]**2; qp=q[0]*p[0]+q[1]*p[1]
    d=(q[0]-p[0])**2+(q[1]-p[1])**2
    c=(2*p2-2*qp+q2)
    return [ -2*q[i]*mp.log(p2/q2) + (c*q[i]+q2*p[i])/d*mp.log(q2/d) for i in range(2)]
def n2(v): return v[0]**2+v[1]**2
q=[mp.mpf('0.7'),mp.mpf('-0.4')]
def terms(p):
    qmp=[q[0]-p[0],q[1]-p[1]]
    A=[x/n2(p) for x in R(q,p)]; B=[x/n2(qmp) for x in R(q,qmp)]
    F=[(A[i]-B[i])/2 for i in range(2)]
    return A,B,F
print("near p=q: |r|, max|A|*r^2, max|F|*r, mean_angle F*r (x comp)")
for r in ['1e-2','1e-4','1e-6','1e-8']:
    r=mp.mpf(r); mA=0; mF=0; avg=[0,0]; N=64
    for k in range(N):
        th=2*mp.pi*k/N; p=[q[0]+r*mp.cos(th), q[1]+r*mp.sin(th)]
        A,B,F=terms(p); mA=max(mA,abs(A[0])*r**2); mF=max(mF,abs(F[0])*r)
        avg=[avg[i]+F[i]/N for i in range(2)]
    print(mp.nstr(r,3), mp.nstr(mA,6), mp.nstr(mF,6), mp.nstr(avg[0],6), mp.nstr(avg[1],6))
print("near p=0:")
for r in ['1e-2','1e-4','1e-6','1e-8']:
    r=mp.mpf(r); mA=0; mF=0; N=64; avg=[0,0]
    for k in range(N):
        th=2*mp.pi*k/N; p=[r*mp.cos(th), r*mp.sin(th)]
        A,B,F=terms(p); mA=max(mA,abs(A[0])*r**2); mF=max(mF,abs(F[0])*r)
        avg=[avg[i]+F[i]/N for i in range(2)]
    print(mp.nstr(r,3), mp.nstr(mA,6), mp.nstr(mF,6), mp.nstr(avg[0],6), mp.nstr(avg[1],6))
