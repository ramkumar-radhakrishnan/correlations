"""Group II Row IV: check the split of G^{mm'} into its delta piece and its r^m r^m' piece."""
import sympy as sp
eps,P,K,S=sp.symbols('epsilon pplus kplus S')
XX,Xr,r2=sp.symbols('XdotXp XrXpr rsq')      # (X.X') , (X.r)(X'.r) , r^2
c=sp.Symbol('c')                              # L = 1/eps + c ,  c = 2 - log(4 pi^2 mu^2 r^2)
L=1/eps+c; dperp=2-2*eps
Ssub={S:K+P}

# --- contraction tables -------------------------------------------------------------
# structure                     contracted with delta^{mm'}        contracted with r^m r^m'
tab={'kmk\'m\'': (XX,           Xr),
     'mm\'kk\'' : (dperp*XX,    r2*XX),
     'mk\'m\'k' : (XX,          Xr),
     'm\'k\'km' : (XX,          Xr)}
def contract(which):
    """[bracket]^{kmk'm'} X^k X'^k' contracted with 'which' part of G"""
    i = 0 if which=='delta' else 1
    return ( dperp*P/S**2*tab['kmk\'m\''][i]
            + 2/K*(tab['mk\'m\'k'][i]-tab['m\'k\'km'][i])
            + (P/K**2+1/P)*tab['mm\'kk\''][i] )

pi=sp.pi
Gdelta=sp.simplify(pi*L/2*contract('delta'))
Grr   =sp.simplify(-pi/r2*contract('rr'))
tot=sp.simplify(sp.expand(Gdelta+Grr))
# drop O(eps): eps*L -> 1 , eps*(finite) -> 0
def lim(e):
    """drop O(eps):  eps*L -> 1 ,  eps*(finite) -> 0"""
    return sp.simplify(sp.series(sp.expand(e),eps,0,1).removeO())
A=(L-1)*XX-2*Xr/r2
B=(L-2)*XX
D=pi*(P/S**2*A+(P/K**2+1/P)*B)

print("="*78); print("1.  DOES THE TWO-PIECE SPLIT REPRODUCE THE FULL CONTRACTION ?"); print("="*78)
print("   delta piece  = %s"%sp.simplify(lim(Gdelta)))
print("   r^m r^m' pc. = %s"%sp.simplify(lim(Grr)))
print("   sum - D      = %s      <-- must be 0"%sp.simplify(sp.expand(lim(tot)-lim(D))))
print()
print("   note  d_perp * L / 2 = (1-eps)(1/eps + c) = 1/eps + c - 1 = L - 1 :",
      sp.simplify(sp.expand(lim(dperp*L/2)-(L-1))))

print(); print("="*78); print("2.  THE delta PIECE, TERM BY TERM"); print("="*78)
print("   delta^{km}delta^{k'm'} x delta^{mm'} = delta^{kk'}        -> (X.X') , keeps its explicit d_perp")
print("   delta^{mm'}delta^{kk'} x delta^{mm'} = d_perp delta^{kk'} -> d_perp (X.X')")
print("   the 2/k+ pair cancels")
print("   => COMMON factor d_perp :   (pi L/2) d_perp [ p+/S^2 + p+/k+^2 + 1/p+ ] (X.X')")
print("                            =  pi (L-1) [ p+/S^2 + p+/k+^2 + 1/p+ ] (X.X')")
chk=sp.simplify(lim(Gdelta) - pi*(L-1)*(P/S**2+P/K**2+1/P)*XX)
print("   check: %s"%sp.simplify(sp.expand(chk)))

print(); print("="*78); print("3.  THE r^m r^m' PIECE, TERM BY TERM"); print("="*78)
print("   delta^{km}delta^{k'm'} x r^m r^m' = (X.r)(X'.r)   ; d_perp -> 2 (multiplies a finite term)")
print("   delta^{mm'}delta^{kk'} x r^m r^m' = r^2 (X.X')")
print("   the 2/k+ pair cancels here too")
print("   => -pi [ 2 p+/S^2 (X.r)(X'.r)/r^2 + (p+/k+^2 + 1/p+)(X.X') ]")
chk2=sp.simplify(lim(Grr) + pi*(2*P/S**2*Xr/r2+(P/K**2+1/P)*XX))
print("   check: %s"%sp.simplify(sp.expand(chk2)))

print()
print("="*78); print("4.  COMPARISON WITH WHAT YOU WROTE"); print("="*78)
print("   your piece 1 :  (g^4/8 pi^3) [p+/S^2 + p+/k+^2 + 1/p+] * L            (no d_perp)")
print("   correct      :  (g^4/4 pi^3) [p+/S^2 + p+/k+^2 + 1/p+] * (L-1)")
print("   ratio correct/yours = 2 (L-1)/L =", sp.simplify(pi*(L-1)/(pi*L/2)))
print()
print("   your piece 2 :  kernel r^m r^m'/(r^2 r^2)   -- one r^2 too many")
print("   correct      :  kernel r^m r^m'/r^2         (G's second term is -pi r^m r^m'/r^2)")
print("   (dimensions: G is dimensionless at eps=0, like the delta L term, so only ONE r^2)")
