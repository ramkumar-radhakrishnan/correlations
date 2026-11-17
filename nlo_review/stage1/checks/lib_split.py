import re
UP="/root/.claude/uploads/0ab97fcf-94a2-5aec-8bad-586400481685/"
FILES={"I":"16e62d79-Group_I.tex","II":"24c35375-Group_II.tex","III":"96b5b029-Group_III.tex"}
BLOCKS=[  # (group, title, [(start,end),...])   line numbers of the final formula(s) of each row
 ("I","Row II",[(48,50),(60,63),(64,69),(71,76)]),
 ("I","Row III",[(85,87),(88,90),(91,110),(111,119),(121,123)]),
 ("I","Row IV",[(132,135)]),
 ("I","Row V, first part",[(239,241)]),
 ("I","Row V, second part",[(243,245),(246,248),(249,263)]),
 ("II","Row I, first part",[(48,50)]),
 ("II","Row I, second part",[(64,66),(67,69),(70,72),(73,75),(76,78),(79,81)]),
 ("II","Row IV",[(250,267)]),
 ("III","Row I",[(46,48)]),
 ("III","Row II",[(57,59)]),
 ("III","Row III, first part",[(74,76)]),
 ("III","Row III, second part",[(88,92)]),
 ("III","Row V",[(111,113)]),
 ("III","Row VI",[(126,128)]),
]
def get(group,a,b):
    L=open(UP+FILES[group]).read().split("\n")
    return "\n".join(L[a-1:b])
def balanced(s,i,op='{',cl='}'):
    d=0
    for j in range(i,len(s)):
        if s[j]==op: d+=1
        elif s[j]==cl:
            d-=1
            if d==0: return j+1
    raise ValueError("unbalanced")
def norm(x): return re.sub(r'\\,|\;|\\!|\s+','',x)
def blue_groups(s):
    out=[]
    for m in re.finditer(r'\{\\color\{blue\}',s):
        e=balanced(s,m.start()); out.append((m.start(),e,s[m.end():e-1]))
    return out
SIZES=r'(?:\\left|\\Big|\\big|\\Bigg|\\bigg)?'
def log_spans(s):
    out=[]
    for m in re.finditer(r'\\log(\^\{2\})?',s):
        j=m.end()
        while s[j]==' ': j+=1
        if s[j]=='{': e=balanced(s,j)
        elif s[j:].startswith('\\frac'):
            e=balanced(s,j+5); e=balanced(s,e)
        else:
            mm=re.match(r'(\\left|\\Big|\\big|\\Bigg|\\bigg)\s*([\(\[])',s[j:])
            assert mm,(s[m.start():m.start()+60])
            size,br=mm.group(1),mm.group(2); cl={'(' :')','[':']'}[br]
            if size=='\\left':
                depth=0; e=None
                for q in re.finditer(r'\\left|\\right',s[j:]):
                    depth+= 1 if q.group()=='\\left' else -1
                    if depth==0: e=j+q.end()+1; break
            else:
                # matching same-size closing bracket
                k=s.index(size+cl,j+mm.end()); e=k+len(size)+1
        out.append((m.start(),e,s[m.start():e]))
    return out
