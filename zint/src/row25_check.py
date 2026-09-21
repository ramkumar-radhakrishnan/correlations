"""Row 25: the A^(1)dag Cbar^dag x B_2 interference -- prefactor, delta, contraction, w' integration."""
import numpy as np, sympy as sp, itertools
g,pi=sp.symbols('g pi',positive=True); P,K=sp.symbols('pplus kplus',positive=True); S=P+K
I=sp.I; D=sp.Symbol('D',positive=True)          # D = p+(x-z)^2 + k+(x-w)^2
print("="*80); print("1.  PREFACTOR -- and the hidden factor in your rewritten delta"); print("="*80)
A1  = I*g/(sp.sqrt(2)*pi*sp.sqrt(S))                     # A^(1)(p+ + k+)
A1d = -I*g/(sp.sqrt(2)*pi*sp.sqrt(S))                    # its conjugate
C1d = -g/(2*pi*sp.sqrt(2*S))*sp.sqrt(P*K)/S              # Cbar^dag at k+_arg = S   (real)
B2  = I*g**2*sp.sqrt(P*K)/(4*pi**2*D*S)                  # B_2^{(1rho)} position space
prod=sp.simplify(A1d*C1d*B2)
print("   A^(1)dag x Cbar^dag x B_2  =  %s"%sp.simplify(prod))
print("   phase (-i)(-1)(+i) = %s   -> the overall MINUS sign in your row is correct"%sp.simplify((-I)*(-1)*I))
full=sp.simplify(sp.Integer(4)/(2*pi)**3*prod)
print("   (4/(2pi)^3) x that       =  %s"%full)
quoted=-1/(2*pi)**3*g**4/(4*pi**4)*(1/K)*P/(S*D)
print("   YOUR quoted prefactor    =  %s"%sp.simplify(quoted))
print("   ratio quoted/derived     =  %s"%sp.simplify(quoted/full))
print()
print("   THAT RATIO IS NOT 1 -- but it is exactly compensated by how you rewrote the delta.")
w,z,y,wp,yp=sp.symbols("w z y wp yp")
defn=yp-wp+P/S*(wp-z)                  # the delta argument straight from the C definition
yours=P/K*(yp-z)+yp-wp                 # the delta argument as you wrote it
print("   definition's argument :  %s"%sp.simplify(defn))
print("   your argument         :  %s"%sp.simplify(yours))
print("   yours / definition    =  %s"%sp.simplify(sp.factor(yours/defn)))
print("   delta^(2)(lambda v) = lambda^{-2} delta^(2)(v)  =>  your delta carries an extra k+^2/S^2,")
print("   and your prefactor carries S^2/k+^2.  Product unchanged:  net ratio = %s   CORRECT."
      %sp.simplify(quoted/full*(K/S)**2))

print(); print("="*80); print("2.  THE TENSOR CONTRACTION  W^{ji} x [C bracket]_{ji}^{mk'}"); print("="*80)
print("   Cbar_{1 j i k'} : relabel the definition's (j,k,i) -> (j,i,k') , i.e. k->i, i->k' :")
print("     [ delta_im delta_jk - (S/p+) delta_jm delta_ik - (S/k+) delta_km delta_ij ]")
print("       ->  delta_k'm delta_ij - (S/p+) delta_jm delta_ik' - (S/k+) delta_im delta_jk'")
print("   contracting over i,j with W^{ji} gives")
print("     R^{mk'} = delta_{k'm} tr W  -  (S/p+) W^{mk'}  -  (S/k+) W^{k'm}")
print()
s=lambda v:v@v; d=2; dd=np.eye(2)
def Wji(x,zz,ww,Pv,Kv):
    Sv=Pv+Kv; A,B,C=s(x-zz),s(x-ww),s(zz-ww)
    Vw=(zz-ww)/C-(x-ww)/(2*B); Vz=(zz-ww)/C+(x-zz)/(2*A)
    return np.eye(2)*(A-B)/(2*C)+(Sv/Kv)*np.outer(x-zz,Vw)+(Sv/Pv)*np.outer(Vz,x-ww)
def yours_bracket(x,zz,ww,Pv,Kv,dperp=2.0):
    Sv=Pv+Kv; A,B,C=s(x-zz),s(x-ww),s(zz-ww)
    xz,xw,zw=x-zz,x-ww,zz-ww
    R=np.zeros((2,2))                                    # R[m,k']
    R+=dd*dperp/(2*C)*(A-B)
    R+=dd*(Sv/Kv)*( (xz@zw)/C - (xz@xw)/(2*B) - (A-B)/(2*C) )
    R+=dd*(Sv/Pv)*( (xw@zw)/C + (xz@xw)/(2*A) - (A-B)/(2*C) )
    R+=-(Sv/Pv)**2*( np.outer(zw,xw)/C + np.outer(xz,zw)/(2*A) )      # <-- YOUR term 4
    R+=-(Sv/Kv)**2*( np.outer(zw,xz)/C - np.outer(xw,xz)/(2*B) )
    R+=-(Sv**2/(Pv*Kv))*( np.outer(xz,zw)/C - np.outer(xz,xw)/(2*B)
                         +np.outer(xw,zw)/C + np.outer(xw,xz)/(2*A) )
    return R
def fixed_bracket(x,zz,ww,Pv,Kv,dperp=2.0):
    Sv=Pv+Kv; A,B,C=s(x-zz),s(x-ww),s(zz-ww); xz,xw,zw=x-zz,x-ww,zz-ww
    R=yours_bracket(x,zz,ww,Pv,Kv,dperp)
    R-=-(Sv/Pv)**2*np.outer(xz,zw)/(2*A)                 # remove your piece
    R+=-(Sv/Pv)**2*np.outer(xz,xw)/(2*A)                 # put the correct one
    return R
rng=np.random.default_rng(7)
for t in range(3):
    x,zz,ww=[rng.normal(size=2)*1.2 for _ in range(3)]; Pv,Kv=(0.37,1.0) if t==0 else (rng.uniform(.3,3),rng.uniform(.3,2))
    Sv=Pv+Kv; W=Wji(x,zz,ww,Pv,Kv)
    Rex=dd*np.trace(W) - (Sv/Pv)*W - (Sv/Kv)*W.T          # R^{mk'} = d_k'm trW - (S/p)W^{mk'} - (S/k)W^{k'm}
    print("   trial %d (p+=%.3f k+=%.3f):  ||yours - exact|| = %9.5f    ||fixed - exact|| = %.2e"
          %(t+1,Pv,Kv,np.abs(yours_bracket(x,zz,ww,Pv,Kv)-Rex).max(),
            np.abs(fixed_bracket(x,zz,ww,Pv,Kv)-Rex).max()))
print()
print("   ==> ONE typo, in your term 4 (the -(S/p+)^2 group):")
print("       you wrote   (x-z)^m (z-w)^{k'} / [2 (x-z)^2]")
print("       it must be  (x-z)^m (x-w)^{k'} / [2 (x-z)^2]")
print("       because that group is -(S/p+)^2 V_z^m (x-w)^{k'} with V_z = (z-w)/(z-w)^2 + (x-z)/2(x-z)^2,")
print("       so BOTH pieces carry the SAME second factor (x-w)^{k'}.  Terms 5 and 6 are correct.")

print(); print("="*80); print("3.  THE w' INTEGRATION"); print("="*80)
print("   your delta: w'-coefficient is -1  =>  Jacobian = 1  and  w' = y' + (p+/k+)(y'-z)")
wsol=sp.solve(sp.Eq(P/K*(yp-z)+yp-wp,0),wp)[0]
print("   w' solved = %s"%sp.simplify(wsol))
print("   w' - z    = %s      =>  (w'-z)^m/(w'-z)^2 = (k+/S) (y'-z)^m/(y'-z)^2"%sp.simplify(sp.factor(wsol-z)))
print("   w' - w    = %s      =>  e^{-ik(w'-w)} = e^{-ik(y'-w)} e^{-ik (p+/k+)(y'-z)}"
      %sp.simplify(wsol-sp.Symbol('w')))
print("   net factor = 1 x (k+/S) = k+/S ;  prefactor (1/k+)(p+/S) x (k+/S) = p+/S^2 ,")
print("   so the 1/k+ out front cancels and p+/S^2 goes into the bracket.  Checking each group:")
for nm,before,after in (("d_perp",sp.Symbol('dperp')/1,sp.Symbol('dperp')*P/S**2),
                        ("S/k+",S/K,P/(K*S)),("S/p+",S/P,1/S),
                        ("-(S/p+)^2",-S**2/P**2,-1/P),("-(S/k+)^2",-S**2/K**2,-P/K**2),
                        ("-S^2/(p+k+)",-S**2/(P*K),-1/K)):
    print("     %-12s : (before) x p+/S^2 = %-22s   you wrote %-14s   diff %s"
          %(nm,sp.simplify(before*P/S**2),sp.simplify(after),sp.simplify(before*P/S**2-after)))
