"""Check the non-instantaneous part of the user's momentum-space B3 against (2pi)^3/k+ x the xi-forms
(onegrho12xi) and (onegrho22xi). Common factor i g^2 f^{abc} rho^c(k-p) stripped; tensor index [i,j]."""
import numpy as np
rng=np.random.default_rng(5); d=np.eye(2)
def T(xi,i,j,l,m): return d[i,l]*d[j,m]-d[j,l]*d[i,m]/xi-d[i,j]*d[l,m]/(1-xi)
def xiform(kp,k,pp,p):
    xi=pp/kp; out=np.zeros((2,2))
    for i in range(2):
        for j in range(2):
            s=0
            for l in range(2):
                for m in range(2):
                    if xi>1:
                        s+=-1/(4*np.pi**3)*np.sqrt(xi)/(1-xi)*(p-k)[m]*(p-xi*k)[l]/(((p-k)@(p-k))*(p@p-xi*k@k))*T(xi,i,j,l,m)
                    else:
                        s+= 1/(4*np.pi**3)*xi**1.5/(1-xi)*(k-p)[m]*(p-xi*k)[l]/((p@p-xi*k@k)*((p-xi*k)@(p-xi*k)))*T(xi,i,j,l,m)
            out[i,j]=s*8*np.pi**3/kp
    return out
def user(kp,k,pp,p):
    out=np.zeros((2,2)); q=k-p
    for i in range(2):
        for j in range(2):
            s=0
            for l in range(2):
                for m in range(2):
                    br=d[j,m]*d[i,l]/pp-d[i,m]*d[j,l]/kp-d[l,m]*d[i,j]/(pp-kp)
                    s+=2*q[l]*np.sqrt(kp*pp)/(kp-pp)**2*(k[m]*pp-kp*p[m])/(k@k*pp-p@p*kp)*br
            reg=(kp-pp)/(q@q) if pp>kp else 1/(k@k/kp-q@q/(kp-pp)-p@p/pp)
            out[i,j]=s*reg
    return out
for pp in [2.3,0.4]:
    kp=1.0; k,p=rng.normal(size=(2,2))
    print("p+/k+ =",pp," user/xi-form:",np.round(user(kp,k,pp,p)/xiform(kp,k,pp,p),4).ravel())
