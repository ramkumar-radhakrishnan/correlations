import re,sys,collections
t=open(sys.argv[1]).read()
# drop removed (orange) blocks
parts=re.split(r'(\{\\color\{orange\}\\noindent\[\\textbf\{Removed\})',t)
keep=parts[0]; rem=''
i=1
while i<len(parts):
    seg=parts[i]+parts[i+1]
    j=seg.find('}\\par\n',)  # block ends with "\n}\par\n"
    k=seg.find('\n}\\par\n')
    rem+=seg[:k]; keep+=seg[k:]
    i+=2
pat=re.compile(r'e\^\{([+-]?i[^{}]*?(?:\{[^{}]*\}[^{}]*?)*)\}')
c=collections.Counter()
for m in pat.finditer(keep):
    s=m.group(1).replace(' ','').replace('\\,','').replace('\\cdot','.')
    if 'k' in s: c[s]+=1
for k,v in sorted(c.items(),key=lambda x:-x[1]): print(v,k)
