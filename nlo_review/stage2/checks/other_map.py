import re,sys
UP="/root/.claude/uploads/0ab97fcf-94a2-5aec-8bad-586400481685/"
def align_blocks(fn):
    L=open(UP+fn).read().split("\n"); out=[]; start=None; sec=""
    for i,l in enumerate(L,1):
        if l.strip().startswith('%'): continue
        m=re.search(r"\\(sub)*section\*?\{([^}]*)",l)
        if m: sec=(sec.split("|")[0] if m.group(1) else "")+"|"+m.group(2) if m.group(1) else m.group(2)
        if '\\begin{align}' in l: start=i
        if '\\end{align}' in l and start:
            s="\n".join(L[start-1:i]); out.append((start,i,sec,s)); start=None
    return out
if __name__=="__main__":
    for st,en,sec,s in align_blocks("f10fa593-Other_terms.tex"):
        lhs=re.search(r"\\frac\{d\^\{3\}N[^=&]*",s)
        print(f"{st:4d}-{en:4d} chars={len(s):6d} sec={sec[:45]:45s} lhs={lhs.group(0)[:38] if lhs else s[13:50]!r}  plus={'_{+}' in s}  pint={'dp^{+}' in s or 'dq^{+}' in s}  nlog={s.count('log')}")
