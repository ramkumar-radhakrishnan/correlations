import sys; sys.path.insert(0,'.')
from mono import *
def canon_vec(expr):
    expr=expr.replace(" ","").replace(r"^{\prime}","'")
    terms=re.findall(r"([+-]?)([a-z]'?)",expr)
    pos=[n for s,n in terms if s!='-']; neg=[n for s,n in terms if s=='-']
    assert len(pos)==1 and len(neg)==1,expr
    a,b=pos[0],neg[0]
    # canonical orientation: sign tracked
    if (a,b)>(b,a): return (b,a),-1
    return (a,b),1
def kernel_key(s):
    s=s.replace(r"\,","").replace(r"^{\prime}","'")
    key=[]; sign=1
    for m in re.finditer(r"\\mathcal\{K\}\(([^()]*)\)\\cdot\s*\\mathcal\{K\}\(([^()]*)\)",s):
        (a,sa),(b,sb)=canon_vec(m.group(1)),canon_vec(m.group(2)); sign*=sa*sb
        key.append("K%s.K%s"%tuple(sorted(["(%s-%s)"%a,"(%s-%s)"%b])))
    for m in re.finditer(r"\\frac\{1\}\{\(([^()]*)\)\^\{2\}\}",s):
        a,_=canon_vec(m.group(1)); key.append("1/(%s-%s)^2"%a)
    return " ".join(sorted(key)),sign
def term_parts(env,s):
    """value of colour*coefficient (kernel stripped), kernel key, has U(z)"""
    key,sign=kernel_key(s)
    hasz=bool(re.search(r"U\^\{[^}]*\}\(z\)",s))
    # evaluate with kernels set to 1: strip kernels then parse
    t=re.sub(r"\\mathcal\{K\}\(([^()]*)\)\\cdot\s*\\mathcal\{K\}\(([^()]*)\)","",s.replace(r"^{\prime}","'"))
    t=re.sub(r"\\frac\{1\}\{\(([^()]*)\)\^\{2\}\}","",t)
    return key,hasz,sign*parse_term(env,t)
def collect(env,blk,acc):
    for t in block_terms(blk):
        k,h,v=term_parts(env,t); acc[(k,h)]=acc.get((k,h),0)+v
    return acc
