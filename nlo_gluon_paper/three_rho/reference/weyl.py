"""Weyl (fully symmetric) decomposition of ordered rho products.
   [rho^a(p), rho^b(q)] = i f^{abc} rho^c(p) delta(p-q)
   Term: c (complex), Us ((i,j,pos)...), fs ((a,b,c)...), rs ((idx,pos)... ordered), mp (dict pos->pos merges)
"""
from itertools import permutations
from math import factorial
import itertools
_cnt=itertools.count()
class T:
    __slots__=('c','Us','fs','rs','mp')
    def __init__(s,c,Us,fs,rs,mp=None):
        s.c=complex(c); s.Us=tuple(Us); s.fs=tuple(fs); s.rs=tuple(rs); s.mp=dict(mp or {})
    def __repr__(s):
        return f"{s.c:.4g} U{s.Us} f{s.fs} r{s.rs} mp{s.mp}"
def subst(t,q,p):
    Us=[(i,j,p if x==q else x) for i,j,x in t.Us]
    rs=[(i,p if x==q else x) for i,x in t.rs]
    mp={k:(p if v==q else v) for k,v in t.mp.items()}; mp[q]=p
    return Us,rs,mp
def commute_sort(t, target_rank):
    """bubble-sort t.rs into increasing target_rank; return (sorted_term, list_of_commutator_terms)"""
    rs=list(t.rs); comms=[]
    changed=True
    while changed:
        changed=False
        for k in range(len(rs)-1):
            (a,p),(b,q)=rs[k],rs[k+1]
            if target_rank[(a,p)]>target_rank[(b,q)]:
                # rho^a(p) rho^b(q) = rho^b(q) rho^a(p) + i f^{abc} rho^c(p) delta(p-q)
                cn='k%d'%next(_cnt)
                newrs=rs[:k]+[(cn,p)]+rs[k+2:]
                tmp=T(t.c*1j,t.Us,t.fs+((a,b,cn),),newrs,t.mp)
                Us,rr,mp=subst(tmp,q,p)
                comms.append(T(tmp.c,Us,tmp.fs,rr,mp))
                rs[k],rs[k+1]=rs[k+1],rs[k]; changed=True
    return comms
def weyl(t,nmin=2):
    """returns list of terms with symmetric rho products (rs order irrelevant)"""
    n=len(t.rs)
    if n<=1: return [t]
    rank={r:k for k,r in enumerate(t.rs)}
    out=[T(t.c,t.Us,t.fs,t.rs,t.mp)]   # Sym(P)
    lower=[]
    for perm in permutations(range(n)):
        if perm==tuple(range(n)): continue
        tp=T(t.c,t.Us,t.fs,[t.rs[k] for k in perm],t.mp)
        for cterm in commute_sort(tp,rank):
            cterm.c*= -1.0/factorial(n)
            lower.append(cterm)
    for l in lower:
        if len(l.rs)>=nmin: out+=weyl(l,nmin)
    return out
