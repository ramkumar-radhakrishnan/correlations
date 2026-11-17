# Symmetric (PBW) decomposition of ordered products of rho's with [rho^a(p),rho^b(q)] = i f^{abc} rho^c(p) delta(p-q)
import sys; sys.path.insert(0,"/tmp/claude-0/-home-user-correlations/0ab97fcf-94a2-5aec-8bad-586400481685/scratchpad/agent_bfkl")
from tworho import *
from fractions import Fraction
from math import factorial
BERN=[Fraction(1),Fraction(-1,2),Fraction(1,6),Fraction(0),Fraction(-1,30)]
_cache={}
def decomp(n):
    """word of n slots (0..n-1) -> list of (coef, tuple of trees); tree = int slot | (tree,tree) meaning commutator"""
    if n in _cache: return _cache[n]
    def rec(slots):
        if len(slots)==1: return [(Fraction(1),(slots[0],))]
        x=slots[0]; out=[]
        for c,ys in rec(slots[1:]):
            k=len(ys)
            for j in range(0,k+1):
                if BERN[j]==0: continue
                for tup in itertools.permutations(range(k),j):
                    g=x
                    for i in reversed(tup): g=(ys[i],g)
                    rest=tuple(ys[i] for i in range(k) if i not in tup)
                    out.append((c*BERN[j]/factorial(j),(g,)+rest))
        return out
    _cache[n]=rec(list(range(n))); return _cache[n]
def leaves(t): return [t] if isinstance(t,int) else leaves(t[0])+leaves(t[1])
def contract(env,Tn,trees):
    """Tn: tensor with one axis per slot. returns tensor with one axis per tree (in order of trees)"""
    f=env["f"]; n=Tn.ndim; L="abcdefghijklmnopqrstuvwxyz"; nxt=[n]
    ops=[Tn]; subs=[L[:n]]
    def idx(t):
        if isinstance(t,int): return L[t]
        a=idx(t[0]); b=idx(t[1]); c=L[nxt[0]]; nxt[0]+=1
        ops.append(1j*f); subs.append(a+b+c); return c
    outs="".join(idx(t) for t in trees)
    return np.einsum(",".join(subs)+"->"+outs,*ops,optimize=True)
def comp(env,texts,pos,cc=False,labels=("x'","y","x","u"),eps=0.0):
    """degree-k symmetric component at positions pos=(p1..pk): coefficient C^{a1..ak}(p1..pk) of sigma(rho(p1)...rho(pk)),
       summed over all ways the operator's rho's (after commutators) land on (p1..pk) in this order of classes"""
    k=len(pos); X=env["X"]; save={l:X[l].copy() for l in labels}; n=env["n"]
    out=np.zeros((n,)*k,complex)
    # which labels occur?
    op0=ev(env,texts,cc); labs=sorted(set(l for key in op0 for l in key))
    for assign in itertools.product(range(k),repeat=len(labs)):
        if len(set(assign))<k: continue
        cnt={}
        for l,a in zip(labs,assign):
            m=cnt.get(a,0); cnt[a]=m+1
            X[l]=pos[a].copy()+(eps*m*np.array([0.6,-0.8]) if eps else 0)
        A=dict(zip(labs,assign))
        op=ev(env,texts,cc)
        for key,Tn in op.items():
            if len(key)<k: continue
            for c,trees in decomp(len(key)):
                if len(trees)!=k: continue
                cls=[set(A[key[s]] for s in leaves(t)) for t in trees]
                if any(len(s)!=1 for s in cls): continue
                cl=[next(iter(s)) for s in cls]
                if sorted(cl)!=list(range(k)): continue
                w=contract(env,Tn,trees)            # axes: trees order
                out+=float(c)*np.moveaxis(w,list(range(k)),cl)   # put tree t on class axis
    for l in save: X[l]=save[l]
    return out
def inv(env,texts,pos,**kw):
    """comp() already sums over all assignments of the operator's rho's to the positions, i.e. it is the
       relabelling-invariant (totally symmetric) coefficient function; normalise by the common prefactor"""
    return comp(env,texts,pos,**kw)/N0
if __name__=="__main__":
    np.seterr(all="ignore")
    B=blocks("9ea88190-BFKL.tex"); R=lambda *ls:[rhs(B[l]) for l in ls]
    env=setup(3,1); rng=np.random.default_rng(101); u=rng.normal(size=2); v=rng.normal(size=2)
    # consistency with the first implementation
    a=inv(env,R(187),[u,v],cc=True); b=Finv(env,R(187),u,v,cc=True); print("selfcheck 2-rho (GI RowIV): %.2e"%np.max(np.abs(a-b)))
    print(len(decomp(3)),len(decomp(4)))
