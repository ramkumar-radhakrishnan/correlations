import sys; sys.path.insert(0,'.')
from grouped import *
fn="9ea88190-BFKL.tex"
B=dict(blocks(fn))
head=lambda blk: re.sub(r"\s+"," ",blk[blk.find("=")+1:blk.find("=")+200])
cand=[l for l in sorted(B) if r"\int_{x,x',z,w,w'}" in B[l][:600] and r"\frac{g^{4}}{16\pi^{5}}" in B[l][:400]]
print("blocks with the common header:",cand)
mono=[l for l in cand if l<569]; fin=[l for l in cand if l>=569]
def is_grouped(blk): return r"\Big[" in blk
worst=0
for seed in (5,6):
    env=setup(3,seed); S={}; per={}
    for l in mono:
        try: collect(env,B[l],S)
        except AssertionError as e: print("PARSE PROBLEM in block",l,str(e)[:300]); raise
    F={}
    for l in fin:
        if is_grouped(B[l]): G={}; collect_grouped(env,B[l],G); per[l]=G
        else: G={}; collect(env,B[l],G); per[l]=G
    if seed==5:
        print("final blocks:",{l:("grouped" if is_grouped(B[l]) else "list",len(per[l])) for l in fin})
        print("%-50s %-6s %13s"%("kernel","U(z)?","sum of rows")+"".join(" %13s"%("block %d"%l) for l in fin))
        for k in sorted(set(S)|set().union(*[set(per[l]) for l in fin])):
            vals=[S.get(k,0)]+[per[l].get(k,0) for l in fin]
            if max(abs(v) for v in vals)>1e-9: print("%-50s %-6s"%(k[0],k[1])+"".join(" %+13.7f"%v for v in vals))
    for l in fin:
        for k in set(S)|set(per[l]): worst=max(worst,abs(S.get(k,0)-per[l].get(k,0)))
print("largest difference (sum of rows vs each final block):",worst)
