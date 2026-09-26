"""2b block by block: which blocks are UV divergent, which have p+ endpoint poles,
and what each one needs.  Residues are quoted for the eq-(1.16) form (1/k+ out front)."""
import numpy as np, sympy as sp
rng = np.random.default_rng(2026)
def sq(v): return float(v@v)

def blocks(D,x,y,z,kp,pp):
    """the six blocks of (1.16) separately, each including its own p+ weight."""
    I=np.eye(D); q=kp-pp
    T1=I*D/(2*sq(x-z))*(sq(y-x)-sq(y-z))*(pp*q/kp**2)
    T2=-(pp/q)*( np.outer(y-x,x-z)/sq(x-z) - np.outer(y-x,y-z)/(2*sq(y-z)) )
    T3=(q/kp)*I*( (y-z)@(x-z)/sq(x-z) + (y-x)@(y-z)/(2*sq(y-x)) - (sq(y-x)-sq(y-z))/(2*sq(x-z)) )
    T4=-(q/pp)*( np.outer(y-z,x-z)/sq(x-z) + np.outer(y-z,y-x)/(2*sq(y-x)) )
    T5=(pp/kp)*I*( (y-x)@(x-z)/sq(x-z) - (y-x)@(y-z)/(2*sq(y-z)) - (sq(y-x)-sq(y-z))/(2*sq(x-z)) )
    T6=-( np.outer(x-z,y-x)/sq(x-z) - np.outer(y-z,y-x)/(2*sq(y-z))
          + np.outer(x-z,y-z)/sq(x-z) + np.outer(y-x,y-z)/(2*sq(y-x)) )
    return [T1,T2,T3,T4,T5,T6]

def contract(D,x,y,z,xp,wp,kp,pp,which):
    V=(xp-wp)/sq(xp-wp); K=(z-x)/sq(z-x)
    B=blocks(D,x,y,z,kp,pp)
    M=sum(B[i] for i in which)
    return float(V@M@K)/(pp*sq(y-x)+(kp-pp)*sq(y-z))/kp

def ang(fn,r,D,m=4096):
    """Exact for the leading residue: it is a quadratic form in uhat, and averaging over
    the 2D coordinate directions +-e_i reproduces <uhat^m uhat^n> = delta^{mn}/D exactly
    while killing the odd pieces."""
    tot=0.0
    for i in range(D):
        e=np.zeros(D); e[i]=1.0
        tot+=fn(r*e)+fn(-r*e)
    return tot/(2*D)

print("="*78)
print("1.  THE UV RESIDUE OF EACH BLOCK  (coefficient of W/(z-x)^2 , W = (x'-w').(y-x)/[...])")
print("="*78)
print("   predicted:  T1 -> eta etabar ,  T2 -> eta/etabar ,  T4 -> etabar/eta ,")
print("               T3 -> -2 etabar/d_perp , T5 -> -2 eta/d_perp , T6 -> +2/d_perp .")
print()
for D in (2,3,4):
    x=np.zeros(D); y=np.array([-0.7,1.1]+[0.4]*(D-2)); xp=np.array([1.6,0.4]+[-0.3]*(D-2))
    wp=np.array([-0.5,-1.3]+[0.8]*(D-2)); kp=1.0; e=0.35; q=1-e
    V=(xp-wp)/sq(xp-wp); R=y-x; W=(V@R)/sq(R)
    pred={0:e*q, 1:e/q, 2:-2*q/D, 3:q/e, 4:-2*e/D, 5:2.0/D}
    print("   d_perp = %d"%D)
    print("     %-6s %18s %18s"%("block","numerical residue","predicted"))
    for i in range(6):
        r=1e-4
        val=r**2*ang(lambda d: contract(D,x,y,x+d,xp,wp,kp,e,[i]),r,D)*kp/W
        print("     T%-5d %18.8f %18.8f"%(i+1,val,pred[i]))
    tot=r**2*ang(lambda d: contract(D,x,y,x+d,xp,wp,kp,e,[2,4,5]),r,D)*kp/W
    print("     T3+T5+T6 %15.3e   <-- predicted 0  (because k+^2/(p+q) - k+/p+ - k+/q = 0)"%tot)
    tot2=r**2*ang(lambda d: contract(D,x,y,x+d,xp,wp,kp,e,[0,1,3]),r,D)*kp/W
    print("     T1+T2+T4 %15.8f   <-- predicted C_UV = %.8f"%(tot2, e*q+e/q+q/e))
    print()

print("="*78)
print("2.  THE SUMMARY TABLE")
print("="*78)
eta=sp.Symbol('eta',positive=True); eb=1-eta
tab=[("T1","eta*etabar","eta*etabar","regular"),
     ("T2","-eta/etabar","eta/etabar","POLE at p+ = k+"),
     ("T3","etabar","-2*etabar/d_perp","regular"),
     ("T4","-etabar/eta","etabar/eta","POLE at p+ = 0"),
     ("T5","eta","-2*eta/d_perp","regular"),
     ("T6","-1","+2/d_perp","regular")]
print("   %-6s %-14s %-20s %-18s"%("block","p+ weight","UV residue","p+ endpoint"))
for a,b,c,d in tab: print("   %-6s %-14s %-20s %-18s"%(a,b,c,d))
print()
print("   T3 + T5 + T6 :  UV residue  (2/d_perp)(1 - etabar - eta) = 0  -> UV FINITE as a group,")
print("                   and all three weights are regular -> no plus prescription either.")
print("   T1 + T2 + T4 :  UV residues are EXACTLY the three terms of C_UV , one each:")
print("                   T1 -> eta etabar ,  T2 -> eta/etabar ,  T4 -> etabar/eta .")

print()
print("="*78)
print("3.  AFTER SUBTRACTION, ARE THE p+ ENDPOINT RESIDUES UV FINITE?")
print("="*78)
print("   Only T2 and T4 have p+ poles, and they are also two of the three UV-divergent blocks.")
print("   So the two divergences sit on the SAME blocks and the order matters.  Subtract first:")
print("       counterterm(T1) = eta etabar  x W e^{i etabar k.(z-x)}/(z-x)^2 / k+")
print("       counterterm(T2) = eta/etabar  x  \"    \"")
print("       counterterm(T4) = etabar/eta  x  \"    \"")
print("   then take the endpoint residue.  Numerically, u^2 x <residue> as u = z-x -> 0 :")
x0=np.zeros(2); y0=np.array([-0.7,1.1]); xp0=np.array([1.6,0.4]); wp0=np.array([-0.5,-1.3])
kvec=np.array([0.5,-0.3]); kp=1.0
V0=(xp0-wp0)/sq(xp0-wp0); R0=y0-x0; W0=(V0@R0)/sq(R0)
def blk(z,e,which):
    return contract(2,x0,y0,z,xp0,wp0,kp,e,which)
def ct(z,e,w):
    return w*W0*np.exp(1j*(1-e)*(kvec@(z-x0)))/sq(z-x0)/kp
# eta -> 0 residue of the T4 block  =  lim eta * (T4 block) ; weight is -etabar/eta
def res_T4(z,eps=1e-7):     return eps*blk(z,eps,[3])
def res_T4_ct(z,eps=1e-7):  return eps*ct(z,eps,(1-eps)/eps)
# eta -> 1 residue of the T2 block ; weight is -eta/etabar
def res_T2(z,eps=1e-7):     return eps*blk(z,1-eps,[1])
def res_T2_ct(z,eps=1e-7):  return eps*ct(z,1-eps,(1-eps)/eps)
print()
def diff(which, side, eps, r, m=8192):
    """u^2 x < (block residue) - (counterterm residue) > , residue taken at distance eps
    from the endpoint.  w is the C_UV weight of that block."""
    e = eps if side == 0 else 1-eps
    w = (1-e)/e if side == 0 else e/(1-e)
    f = lambda d: (eps*blk(x0+d, e, [which]) - eps*ct(x0+d, e, w)).real
    ph = (np.arange(m)+0.5)*(2*np.pi/m)
    return r**2*np.mean([f(r*np.array([np.cos(t),np.sin(t)])) for t in ph])
for nm, which, side in (("T4 (eta -> 0)",3,0), ("T2 (eta -> 1)",1,1)):
    eps = 1e-6
    e = eps if side==0 else 1-eps
    w = (1-e)/e if side==0 else e/(1-e)
    r = 1e-4
    ph=(np.arange(8192)+0.5)*(2*np.pi/8192)
    A=r**2*np.mean([ (eps*blk(x0+r*np.array([np.cos(t),np.sin(t)]),e,[which])).real for t in ph])
    B=r**2*np.mean([ (eps*ct(x0+r*np.array([np.cos(t),np.sin(t)]),e,w)).real for t in ph])
    print("   %-16s  block residue = %.8f    counterterm = %.8f"%(nm,A,B))
print()
print("   and the DIFFERENCE, scanned in both the transverse distance and the extraction eps:")
print("   %-16s %-10s %12s %12s %12s"%("block","eps","r=1e-3","1e-4","1e-5"))
for nm, which, side in (("T4 (eta -> 0)",3,0), ("T2 (eta -> 1)",1,1)):
    for eps in (1e-5,1e-6,1e-7,1e-8):
        v=[diff(which,side,eps,r) for r in (1e-3,1e-4,1e-5)]
        print("   %-16s %-10.0e %12.2e %12.2e %12.2e"%(nm,eps,*v))
print()
print("   Both differences fall as |z-x|^2 while each residue separately is a CONSTANT (0.0322),")
print("   i.e. log divergent.  So the subtraction removes the endpoint residues' UV exactly, and")
print("   the plus prescription on the subtracted T2 and T4 blocks is legitimate.")
print("   (T2's numbers also shrink with the extraction eps: that is the O(eps) part of the")
print("    residue, which vanishes in the true endpoint limit.)")
print()
print("="*78)
print("4.  THE p+ INTEGRALS OF THE THREE UV TERMS, BLOCK BY BLOCK")
print("="*78)
eta,lam=sp.symbols('eta lamb',positive=True); eb=1-eta
L=sp.Symbol('L',positive=True)
import mpmath as mp
mp.mp.dps = 25
closed = {"T1": ("1/6", "-5/18"),
          "T2": ("L - 1", "-L^2 + 2"),
          "T4": ("L - 1", "2 - pi^2/3")}
wf = {"T1": lambda e: e*(1-e), "T2": lambda e: e/(1-e), "T4": lambda e: (1-e)/e}
print("   %-5s %-12s %-16s %-14s %-16s"%("block","Int w","(quadrature)","Int w log eb^2","(quadrature)"))
Lm = mp.mpf('1e-8'); Lv = -mp.log(Lm)
num = {"T1": (mp.mpf(1)/6, -mp.mpf(5)/18),
       "T2": (Lv-1, -Lv**2+2),
       "T4": (Lv-1, 2-mp.pi**2/3)}
for nm in ("T1","T2","T4"):
    q0 = mp.quad(wf[nm], [Lm, 0.5, 1-Lm])
    q1 = mp.quad(lambda e: wf[nm](e)*mp.log((1-e)**2), [Lm, 0.5, 1-Lm])
    print("   %-5s %-12s %-16s %-14s %-16s"
          %(nm, closed[nm][0], mp.nstr(q0,9)+" / "+mp.nstr(num[nm][0],9),
            closed[nm][1], mp.nstr(q1,9)+" / "+mp.nstr(num[nm][1],9)))
print("   %-5s %-12s %-16s %-14s %-16s"%("sum","2L - 11/6","","-L^2 + 67/18 - pi^2/3",""))
print("   (each cell: quadrature at Lambda/k+ = 1e-8  /  the closed form at the same Lambda)")
print()
print("   so each block's UV term is   -[1/(2pi)^3](g^4/16 pi^4) x")
print("      { (Int w) [1/eps + log(k_perp^2/mubar^2)]  +  (Int w log etabar^2) }  x  (W x colour)")
print("   and the three add up to the I0, I1 of the previous note.")

print()
print("="*78)
print("5.  THE RECIPE FOR 2b")
print("="*78)
print("""   (B)  T3 + T5 + T6 :  NOTHING SPECIAL.
          UV residue cancels within the group, and the weights etabar, eta, -1 are regular at
          both p+ endpoints.  Integrate dp+ from Lambda to k+-Lambda directly (or from 0 to k+,
          the difference is O(Lambda)), transverse integrals at d_perp = 2.  No dim reg, no
          plus prescription.  Keep the three together -- individually they ARE UV divergent.

   (A)  T1 + T2 + T4 :  subtract, block by block, the counterterm
             w(eta) x W x e^{i etabar k.(z-x)}/(z-x)^2 / k+ ,   w = eta etabar, eta/etabar, etabar/eta
          (i)  the counterterms integrate in z to -pi[1/eps + log(etabar^2 k_perp^2/mubar^2)] and
               their p+ integrals are elementary -- the table in section 4.  No plus prescription
               needed there, the closed forms are exact.
          (ii) the three remainders are UV finite.  T1's remainder has no p+ pole either, so it
               is an ordinary double integral.  T2's and T4's remainders keep one endpoint pole
               each, so use a ONE-SIDED plus prescription on each:
                    T4 :  [1/eta]_+     at p+ -> 0
                    T2 :  [1/etabar]_+  at p+ -> k+
               legitimate now, by section 3.""")
