import numpy as np
rng=np.random.default_rng(1)
def dot(a,b): return a@b
def sq(a): return a@a
def A0f(x,w,z):
    a=x-w; b=x-z; c=z-w
    det=np.sqrt(sq(a)*sq(b)-dot(a,b)**2)
    return np.arccos(dot(a,c)/np.sqrt(sq(a)*sq(c)))/det
def A0_author(x,w,z):
    a=x-w;b=x-z
    det=np.sqrt(sq(a)*sq(b)-dot(a,b)**2)
    return (np.arctan((sq(b)-dot(a,b))/det)+np.arctan(dot(a,b)/det))/det
def author(x,y,z,w,xp,wp):
    P=xp-wp; a=x-w; b=x-z; c=z-w; Y=y-z
    A=sq(a);B=sq(b);C=sq(c);PP=sq(P);YY=sq(Y)
    L1=np.log(abs(B/(B-A))); L2=np.log(abs((A-B)/A)); A0=A0_author(x,w,z)
    t=0
    # N1 (units i g^4 f/(8pi^5 k+)): first two have 1/16 -> factor 1/2
    t+=0.5*dot(P,Y)/(PP*YY)*(A+B)/(B-A)**2*L1
    t+=0.5*dot(P,Y)/(PP*YY)*(A+B)/(B-A)**2*L2
    t+=-dot(P,a)*dot(b,Y)/(PP*A*B*YY)*np.log(B/A)
    # N2 overall minus
    s=0
    s+= dot(P,c)*dot(Y,a)/(PP*YY*C*(A-B))*L1
    s+=-dot(P,a)*dot(Y,c)*B/(PP*YY*C*(A-B)*A)*L1
    s+= dot(P,Y)*dot(c,a)*B/(PP*YY*C*(A-B)**2)*L1
    s+=-dot(P,c)*dot(Y,b)*A/(PP*YY*C*(A-B)*B)*L1
    s+= dot(P,b)*dot(Y,c)/(PP*YY*C*(A-B))*L1
    s+=-dot(P,Y)*dot(c,b)*A/(PP*YY*C*(A-B)**2)*L1
    t-=s
    # N3a
    t+=-dot(P,Y)*dot(c,a)*B/(PP*YY*C*(A-B)**2)*L2
    t+= dot(P,Y)*dot(c,b)*A/(PP*YY*C*(A-B)**2)*L2
    t+=-dot(P,c)*dot(Y,a)/(PP*YY*C*(A-B))*L2
    t+= dot(P,c)*dot(Y,b)*A/(PP*YY*C*(A-B)*B)*L2
    t+= dot(Y,c)*dot(P,a)*B/(PP*YY*C*(A-B)*A)*L2
    t+=-dot(Y,c)*dot(P,b)/(PP*YY*C*(A-B))*L2
    t+= A0*dot(P,c)*dot(Y,a)/(PP*YY*C)
    t+= A0*dot(Y,c)*dot(P,b)/(PP*YY*C)
    t+=-A0*dot(a,b)*dot(P,c)*dot(Y,b)/(PP*YY*C*B)
    # N3b first two
    t+=-A0*dot(a,b)*dot(Y,c)*dot(P,a)/(PP*YY*C*A)
    t+=-0.5*np.log(C/A)*dot(P,c)*dot(Y,b)/(PP*YY*C*B)
    return t
def fixsoft(x,y,z,w,xp,wp):
    P=xp-wp; a=x-w; b=x-z; c=z-w; Y=y-z
    A=sq(a);B=sq(b);C=sq(c)
    KP=P/sq(P); KY=Y/sq(Y)
    A0=A0f(x,w,z)
    K=np.zeros((2,2))
    d=np.eye(2)
    K+=np.log(B/A)*(d/(2*C)-np.outer(a,b)/(A*B)+1/(C*(A-B))*(np.outer(a,c)*B/A+np.outer(c,b)*A/B-np.outer(c,a)-np.outer(b,c)))
    K+=A0/C*(np.outer(c,a)+np.outer(b,c)-dot(a,b)*(np.outer(c,b)/B+np.outer(a,c)/A))
    K+=-1/(2*C)*np.log(C/A)*np.outer(c,b)/B
    return KP@K@KY
for _ in range(6):
    x,y,z,w,xp,wp=[rng.normal(size=2)*2 for _ in range(6)]
    print(author(x,y,z,w,xp,wp), fixsoft(x,y,z,w,xp,wp), A0f(x,w,z)-A0_author(x,w,z))
