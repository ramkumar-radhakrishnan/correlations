# Order of limits, follow-up to s7: does a range-restricted coupling counterterm (UV half of the dressed
# produced-gluon normalisation <gg|gg>) supply an ell = log(1-k+/vee) term that cancels the UV pole (ell/2)(2/eps)
# of Group I Row II?  Units: (2/eps) * u * LO, u = alpha_s N_c/pi (as in s1_pole_algebra.py).
# Sources: Calculation_main_1_checked.tex, Normalization: <(Psi)_gg|(Psi)_gg> = -(g^2Nc/8pi^3)(2log(Lambda/k+)+11/6) int d^2p/p^2,
# the gg state of the produced gluon k+ has daughters p+, k+-p+, range Lambda < p+ < k+-Lambda;
# A^(3) one-rho pole factor: 11/3 + log(Lambda/k+) + log(vee/k+ - 1) - log(vee/Lambda) = 11/3 - 2L + ell.
import sympy as sp
t,lam,V=sp.symbols('t lambda V',positive=True)       # t=p+/k+, lam=Lambda/k+, V=vee/k+ (>1)
L=sp.log(1/lam); ell=sp.log(1-1/V)
P=t*(1-t)+t/(1-t)+(1-t)/t                              # P_gg+/(2Nc) in t (symmetric t<->1-t)
# 0) A^(3) bracket, Row II UV pole: 11/3 - 2L + ell
br=sp.Rational(11,3)+sp.log(lam)+sp.log(V-1)-sp.log(V/lam)
print("[0] A3 pole factor - (11/3-2L+ell) =",sp.simplify(sp.expand_log(br-(sp.Rational(11,3)-2*L+ell),force=True)))
# 1) normalisation: int_lam^{1-lam} P dt = 2L - 11/6 (+O(lam))
Nfull=sp.series(sp.integrate(sp.expand(P),(t,lam,1-lam)),lam,0,1).removeO()
print("[1] int_lam^{1-lam} P dt =",sp.simplify(Nfull))
# range-restricted as the loop of A^(3): p+ < vee-k+ <=> t < V-1. For V>2 the upper limit min(1-lam,V-1)=1-lam: unchanged.
print("    restriction t < V-1 is inactive for V > 2 (daughters of the produced gluon have p+ < k+ < vee-k+)")
# 2) UV poles (units u*2/eps_UV*LO), two-eps bookkeeping
uv={'Row II':sp.Rational(11,6)-L+ell/2,'Row V':-sp.Rational(11,6)+2*L,'norm UV':(sp.Rational(11,12)-L),'coupling CT':-sp.Rational(11,12)}
print("[2] UV sum =",sp.simplify(sum(uv.values())),"; L parts:",{k:sp.simplify(sp.diff(v,lam)*(-lam)) for k,v in uv.items()})
# 3) IR side with the range-restricted FF counterterm (s7): the FF CT supplies +ell (R1) or 0 (R2) beyond the full-range table
for name,ctell in [('R1 (P_gg on [0,1] applied to theta(xi>k+/vee) Pi)',ell),('R2 (plus subtraction on [k+/vee,1])',0)]:
    IR2=(-sp.Rational(11,12)+L)+(-L-ell)+sp.Rational(11,12)+ctell        # norm IR half + Row IV xi->1 + delta part; P_gg(x)Pi cancels
    tot2=sp.simplify(sum(uv.values())+IR2)
    tot1=sp.simplify(uv['Row II']+uv['Row V']+(-L-ell+sp.Rational(11,12))+uv['coupling CT']+ctell)   # normalisation = 0
    print("[3] (a) restricted FF CT only,",name,": IR sum =",sp.simplify(IR2),"; two-eps total =",tot2,"; single-eps total =",tot1)
print("[4] (b) + range-restricted coupling CT (norm. UV half with p+<vee-k+): identical to (a) for vee>2k+, since [1] is unchanged;")
print("    leftover at finite vee, eps->0 first: (ell/2)(2/eps) [R1] or -(ell/2)(2/eps) [R2], ell=log(1-k+/vee)=-k+/vee+O((k+/vee)^2)")
