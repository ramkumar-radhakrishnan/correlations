"""Group II Row IV: check the chain (0.21) -> (0.22) -> (0.23) -> (0.24) of the revised note."""
import sympy as sp
e,P,K,S=sp.symbols('epsilon pplus kplus S'); XX,Xr,r2,M=sp.symbols('XX Xr rsq M'); pi=sp.pi
dp=2-2*e
# contraction tables for the four bracket structures, against delta^{mm'} and against r^m r^m'
tD={'km':XX,'mm':dp*XX,'mk':XX,'mk2':XX}
tR={'km':Xr,'mm':r2*XX,'mk':Xr,'mk2':Xr}
def full(t):    # the complete (0.21) bracket, including the 2/k+ pair
    return dp*P/S**2*t['km'] + 2/K*(t['mk']-t['mk2']) + (P/K**2+1/P)*t['mm']
G21 = pi*M/dp*full(tD) - pi/r2*full(tR)              # eq (0.21) contracted
print("="*76); print("(0.21) contracted, d_perp kept:"); print("="*76)
print("   %s"%sp.simplify(sp.expand(G21)))
D = pi*( P/S**2*(M*XX - dp*Xr/r2) + (P/K**2+1/P)*(M-1)*XX )
print("   equals  pi[ p+/S^2 ( M (X.X') - d_perp (X.r)(X'.r)/r^2 ) + (p+/k+^2+1/p+)(M-1)(X.X') ] ?  %s"
      %sp.simplify(sp.expand(G21-D)))
print("   (the 2/k+ pair drops out of BOTH contractions, so (0.22) is right to omit it)")

print(); print("="*76); print("(0.22): two parts"); print("="*76)
p22a =  pi*M*(P/S**2+P/K**2+1/P)*XX                      # part 1  (delta piece, d_perp cancelled)
p22b = -pi*( dp*P/S**2*Xr/r2 + (P/K**2+1/P)*XX )          # part 2  (r^m r^m' piece, single r^2)
print("   (0.22) - (0.21) = %s"%sp.simplify(sp.expand(p22a+p22b-G21)))

print(); print("="*76); print("(0.23): three blocks"); print("="*76)
b1 =  pi*M*(P/S**2+P/K**2+1/P)*XX
b2 = -pi*dp*P/S**2*Xr/r2
b3 = -pi*(P/K**2+1/P)*XX
print("   (0.23) - (0.22) = %s"%sp.simplify(sp.expand(b1+b2+b3-(p22a+p22b))))

print(); print("="*76); print("(0.24): the merge of blocks 1 and 3"); print("="*76)
m24 = pi*(M-1)*(P/S**2+P/K**2+1/P)*XX + b2
print("   (0.24) - (0.23) = %s"%sp.simplify(sp.expand(m24-(b1+b2+b3))))
print("   -> (0.24) is SHORT by  pi (p+/S^2)(X.X') .  The identity is")
print("        M[A1+A2+A3] - [A2+A3] = (M-1)[A1+A2+A3] + A1 ,   A1 = p+/S^2")
print("      and the leftover +A1 (X.X') was dropped.")
fix = pi*(M-1)*(P/S**2+P/K**2+1/P)*XX - pi*P/S**2*(dp*Xr/r2 - XX)
print("   with the second block written as  -(p+/S^2)[ d_perp (X.rhat)(X'.rhat) - (X.X') ] :")
print("     difference = %s"%sp.simplify(sp.expand(fix-(b1+b2+b3))))

print(); print("="*76); print("(0.19) vs (0.20): the power of mu"); print("="*76)
mu,r=sp.symbols('mu r',positive=True); gE=sp.EulerGamma
N=sp.gamma(1+e)*sp.gamma(1-e)**2/sp.gamma(2-2*e)
for lab,pw in (("mu^{+2eps}, as in (0.19)",+1),("mu^{-2eps}",-1)):
    Om=sp.series(sp.expand(sp.pi**(-e)*(r**2)**(-e)*mu**(2*pw*e)*N),e,0,2).removeO()
    br=sp.simplify(sp.series(sp.expand(Om*(1/e-1)),e,0,1).removeO())
    br=br.subs(mu,sp.sqrt(mu**2*sp.exp(gE)/(4*sp.pi)))     # MSbar: mu^2 -> mu^2_MSbar e^gE/(4 pi)
    print("   %-24s -> %s"%(lab,sp.simplify(sp.logcombine(sp.expand(br),force=True))))
print("   (0.20) states  1/eps + 1 - log(4 pi^2 mu^2 r^2) : that is the mu^{-2eps} line.")
print("   So (0.18)/(0.19) should carry mu^{-2eps}, not mu^{+2eps} (also the dimensionally")
print("   consistent choice, since G ~ (r^2)^{-eps}).")
