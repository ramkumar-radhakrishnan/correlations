"""Group II Row IV: check the four-term evolution/DGLAP/finite form (N_c counting)."""
import sympy as sp
Nc,g,pi,K,eps,Lr=sp.symbols('N_c g pi kplus epsbarinv L_r')
XX,RR,Chat,dp,x=sp.symbols('XX RR Chat d_perp xi')
P0=1/(2*pi)**3*g**4/(8*pi**4)*Nc/K          # 1/(2pi)^3 g^4/8pi^4 Nc/k+
Phat=2*Nc*Chat                               # Phat_gg = 2 Nc Chat ,  Chat = xi xib + xib/xi + xi/[1-xi]_+
pole=Chat-x*(1-x)-(1-x)/x                    # = xi/[1-xi]_+
mine=P0*( Chat*(eps+1-Lr)*XX - dp*x*(1-x)*RR - pole*XX )   # verified ascending form
LOfac=1/(2*pi)**3*g**2/(pi**2*K)             # d^3N_LO = LOfac x N_LO   (= 2/(2pi)^3 |A^(1)|^2)
t2=g**2/(8*pi**2)*eps*(Phat/2)*LOfac         # their DGLAP block, ascending
t3=-P0*( dp*x*(1-x)*RR + pole*XX )           # their tensor block, ascending
t4w=+P0*(Phat/2)*(1-Lr)*XX                   # their last block AS WRITTEN
t4c=+P0*Chat*(1-Lr)*XX                       # corrected
print("term 2 / (P0 Chat epsbarinv XX)  =", sp.simplify(t2/(P0*Chat*eps)), "  (per unit XX) -> 1 means OK")
print("term 3 matches                   :", sp.simplify(sp.expand(t3-P0*(-dp*x*(1-x)*RR-pole*XX)))==0)
print("term 4 written / correct         =", sp.simplify(t4w/t4c), "   <-- ONE N_c TOO MANY")
print("sum with term 4 fixed - verified =", sp.simplify(sp.expand(P0*Chat*eps*XX+t3+t4c-mine)))
