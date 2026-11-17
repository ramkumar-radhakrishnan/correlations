import sys; sys.path.insert(0,'.')
from mono2 import *
def grouped_terms(blk):
    """block of the form  ... \\Bigg[  KERNELS \\Big[ monomials \\Big]  KERNELS \\Big[ ... \\Big] ... \\Bigg]  -> list of monomial strings (kernel included)"""
    a=blk.index(r"\Bigg["); b=blk.rindex(r"\Bigg]")
    body=blk[a+6:b]
    body=re.sub(r"\\notag\s*\\\\\s*&?"," ",body).replace("&"," ")
    out=[]
    for m in re.finditer(r"([^\[\]]*?)\\Big\[(.*?)\\Big\]",body,re.S):
        ker=m.group(1).strip().lstrip("+").strip(); inner=m.group(2)
        assert r"\rho" not in ker,ker
        parts=re.split(r"(?=[+-]\s*(?:\d|\\,|N_|f\^|U\^|\\frac|\\rho))",inner)
        for p in parts:
            p=p.strip()
            if p: out.append(p[0]+" "+ker+" "+p[1:] if p[0] in "+-" else ker+" "+p)
    return out
def collect_grouped(env,blk,acc):
    for t in grouped_terms(blk):
        k,h,v=term_parts(env,t); acc[(k,h)]=acc.get((k,h),0)+v
    return acc
if __name__=="__main__":
    B=dict(blocks("9aa89920-Evolution_JIMWLK__terms.tex"))
    worst=0
    for seed in (5,6,7):
        env=setup(3,seed); T={}; F={}
        collect(env,B[656],T); collect(env,B[687],T)
        collect_grouped(env,B[709] if 709 in B else B[708],F)
        for k in sorted(set(T)|set(F)):
            worst=max(worst,abs(T.get(k,0)-F.get(k,0)))
            if seed==5 and abs(T.get(k,0)-F.get(k,0))>1e-9: print("MISMATCH",k,T.get(k,0),F.get(k,0))
    print("JIMWLK: listed total (2 blocks) vs final grouped formula: number of terms",len(grouped_terms(B[709] if 709 in B else B[708])),"; largest difference",worst)
