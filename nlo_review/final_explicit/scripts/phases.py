# Inventory of all k-dependent phases e^{-ik...} in a TeX file (up to \end{document}).
import re,sys
src=open(sys.argv[1]).read()
end=src.find(r'\end{document}')
body=src[:end]
def grab(s,i):
    # s[i]=='{' ; return balanced content
    d=0
    for j in range(i,len(s)):
        if s[j]=='{': d+=1
        elif s[j]=='}':
            d-=1
            if d==0: return s[i+1:j],j
    raise ValueError
out=[]
for m in re.finditer(r'e\^\{',body):
    c,j=grab(body,m.end()-1)
    if 'k' not in c: continue
    line=body[:m.start()].count('\n')+1
    cc=re.sub(r'\\,|\;|\s|\\cdot|\\left|\\right|\{|\}','',c).replace("^\\prime","'").replace("\\prime","'")
    out.append((line,cc))
for l,c in out: print(l,c)
print('total',len(out))
