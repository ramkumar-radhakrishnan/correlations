import re,sys,difflib
UP="/root/.claude/uploads/0ab97fcf-94a2-5aec-8bad-586400481685/"
F={"I":"16e62d79-Group_I.tex","II":"24c35375-Group_II.tex","III":"96b5b029-Group_III.tex","O":"f10fa593-Other_terms.tex"}
def get(g,a,b):
    L=open(UP+F[g]).read().split("\n"); return "\n".join(L[a-1:b])
def norm(s):
    s=re.sub(r"\\notag|\\nonumber","",s); s=re.sub(r"\\\\(\[[^\]]*\])?","",s); s=s.replace("&","")
    s=re.sub(r"\\,|\;|\\!|\s+","",s)
    s=s.replace(r"\left","").replace(r"\right","")
    for b in (r"\Bigg",r"\bigg",r"\Big",r"\big"): s=s.replace(b,"")
    s=re.sub(r"\{\\color\{(blue|red|green)\}","{",s)
    return s
def show(g,a,b,oa,ob,ctx=28,minlen=1):
    A=norm(get(g,a,b)); B=norm(get("O",oa,ob))
    sm=difflib.SequenceMatcher(None,A,B,autojunk=False)
    print(f"##### Group {g} lines {a}-{b} ({len(A)} chars)  vs  Other_terms lines {oa}-{ob} ({len(B)} chars)   similarity {sm.ratio():.3f}")
    for tag,i1,i2,j1,j2 in sm.get_opcodes():
        if tag=="equal": continue
        print(f"  [{tag}] ...{A[max(0,i1-ctx):i1]}⟦{A[i1:i2][:200]}⟧ => ⟦{B[j1:j2][:200]}⟧{B[j2:j2+ctx]}...")
if __name__=="__main__":
    args=sys.argv[1:]
    for k in range(0,len(args),5):
        g,a,b,oa,ob=args[k:k+5]; show(g,int(a),int(b),int(oa),int(ob))
