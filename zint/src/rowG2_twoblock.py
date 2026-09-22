"""Group II Row IV: can blocks 1 and 3 be merged into a single (M-1) bracket?  No."""
import sympy as sp
e,P,K,S=sp.symbols('epsilon pplus kplus S'); XX,Xr,r2,M=sp.symbols('XX Xr rsq M'); pi=sp.pi
dp=2-2*e
D     = pi*( P/S**2*(M*XX - dp*Xr/r2) + (P/K**2+1/P)*(M-1)*XX )
two   = pi*( (M-1)*(P/S**2 + P/K**2 + 1/P)*XX - dp*P/S**2*Xr/r2 )
fixed = pi*( (M-1)*(P/S**2 + P/K**2 + 1/P)*XX - P/S**2*(dp*Xr/r2 - XX) )
print("naive merge  :  D - it =", sp.simplify(sp.expand(D-two)), "   <-- short by pi (p+/S^2)(X.X')")
print("fixed merge  :  D - it =", sp.simplify(sp.expand(D-fixed)))
print()
print("reason: the -1 comes from the trace of -pi rhat^m rhat^m' and only the")
print("  delta_mm' delta_kk' structure returns the SAME tensor (X.X') when contracted with it;")
print("  delta_km delta_k'm' returns (X.rhat)(X'.rhat), a different tensor, so its delta part")
print("  keeps the full M and its rhat part stays in the second block.")
print()
print("so the second block must read   -(p+/S^2) [ d_perp (X.rhat)(X'.rhat) - (X.X') ]")
print("                             =  -(p+/S^2) X^m X'^m' [ d_perp rhat^m rhat^m' - delta^mm' ]")
