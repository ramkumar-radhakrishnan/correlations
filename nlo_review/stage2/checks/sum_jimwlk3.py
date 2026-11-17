import sys; sys.path.insert(0,'.')
from mono2 import *
B=dict(blocks("9aa89920-Evolution_JIMWLK__terms.tex"))
mono=[l for l in sorted(B) if l<655]
worst=0
for seed in (5,6):
    env=setup(3,seed); S={}; F={}
    for l in mono: collect(env,B[l],S)
    collect(env,B[656],F); collect(env,B[687],F)
    keys=sorted(set(S)|set(F))
    if seed==5: print("%-52s %-6s %14s %14s %10s"%("kernel","U(z)?","sum of rows","listed total","diff"))
    for k in keys:
        a=S.get(k,0); b=F.get(k,0); worst=max(worst,abs(a-b))
        if seed==5 and (abs(a)>1e-9 or abs(b)>1e-9): print("%-52s %-6s %+14.8f %+14.8f %+10.1e"%(k[0],k[1],a,b,a-b))
print("largest difference over both seeds:",worst)
