import sys,re; sys.path.insert(0,'.')
from lib_split import *
BL=lambda x: r"{\color{blue}"+x+"}"
LL=r"\log\frac{k^{+}}{\Lambda}"; LV=r"\log\frac{\vee}{k^{+}}"
r1=r"\log\Big(1-\frac{k^{+}}{\vee}\Big)"; r2=r"\log\Big(1+\frac{\Lambda}{k^{+}}\Big)"; r3=r"\log\Big(1-\frac{\Lambda}{k^{+}}\Big)"
both=r"\Big["+BL(LL+"+"+LV)+"+"+r1+r"\Big]"
BLUE={
 norm(r"-2\,\log{\left(\frac{k^{+}}{\Lambda}\right)}"): BL(r"-2\,"+LL),
 norm(r"\left[\log{\left(\frac{\vee}{k^{+}}-1\right)}+\log{\left(\frac{k^{+}}{\Lambda}\right)}\right]"): both,
 norm(r"\left[\log{\left(\frac{k^{+}}{\Lambda}\right)}+\log{\left(\frac{\vee}{k^{+}}-1\right)}\right]"): both,
 norm(r"\log{\left[\frac{\Lambda}{k^{+}-\Lambda}\right]}"): r"\Big[-"+BL(LL)+"-"+r3+r"\Big]",
 norm(r"\left[\log{\left[\frac{\vee}{k^{+}}\right]}+\log{\left[\frac{k^{+}}{\Lambda}\right]}\right]"): BL(r"\Big["+LV+"+"+LL+r"\Big]"),
 norm(r"\log{\left[\frac{\vee}{\Lambda+k^{+}}\right]}"): r"\Big["+BL(LV)+"-"+r2+r"\Big]",
 norm(r"\log{\left[\frac{\vee}{k^{+}+\Lambda}\right]}"): r"\Big["+BL(LV)+"-"+r2+r"\Big]",
 norm(r"\left[\log{\left[\frac{\vee}{\Lambda}\right]}+\log{\left[1-\frac{k^{+}}{\vee}\right]}\right]"): both,
 norm(r"\log{\left(\frac{k^{+}-\Lambda}{\Lambda}\right)}"): r"\Big["+BL(LL)+"+"+r3+r"\Big]",
 norm(r"\log{\left[\frac{k^{+}-\Lambda}{\Lambda}\right]}"): r"\Big["+BL(LL)+"+"+r3+r"\Big]",
}
def composite(key):
    # log[((V-k)A^2 + k B^2)/(Lam A^2 + k B^2)] in its orderings -> blue log(V/k) + exact remainder
    m=re.fullmatch(r"\\log\{\\left\[\\frac\{(.*)\}\{(.*)\}\\right\]\}",key)
    if not m: return None
    num,den=m.group(1),m.group(2)
    a=re.fullmatch(r"\(\\vee-k\^\{\+\}\)\(([^()]*)\)\^\{2\}\+k\^\{\+\}\(([^()]*)\)\^\{2\}",num) 
    if a: A,B=a.group(1),a.group(2)
    else:
        a=re.fullmatch(r"k\^\{\+\}\(([^()]*)\)\^\{2\}\+\(\\vee-k\^\{\+\}\)\(([^()]*)\)\^\{2\}",num)
        if not a: return None
        B,A=a.group(1),a.group(2)
    d1=norm(rf"\Lambda({A})^{{2}}+k^{{+}}({B})^{{2}}"); d2=norm(rf"k^{{+}}({B})^{{2}}+\Lambda({A})^{{2}}")
    if den not in (d1,d2): return None
    rem=rf"\log\frac{{\big(1-\frac{{k^{{+}}}}{{\vee}}\big)({A})^{{2}}+\frac{{k^{{+}}}}{{\vee}}({B})^{{2}}}}{{({B})^{{2}}+\frac{{\Lambda}}{{k^{{+}}}}({A})^{{2}}}}"
    return r"\Big["+BL(LV)+"+"+rem+r"\Big]"
A2="(x-z)^{2}"; B2="(x-w)^{2}"; kv=r"\frac{k^{+}}{\vee}"; lk=r"\frac{\Lambda}{k^{+}}"
OTHER={
 norm(r"\log\Big(\frac{k^{+}}{\Lambda}\Big)"): BL(LL),
 norm(r"\log\Big[\frac{k^{+}}{\Lambda}\Big]"): BL(LL),
 norm(r"\log\frac{k^{+}}{\Lambda}"): BL(LL),
 norm(r"\log\frac{\vee}{k^{+}}"): BL(LV),
 norm(r"\log{\left[\frac{\vee}{\Lambda+k^{+}}\right]}"): r"\Big["+BL(LV)+"-"+r2+r"\Big]",
 norm(r"\log{\left[\frac{\left[\vee(x-z)^{2}-k^{+}(x-w)^{2}\right]}{\left[(\Lambda+k^{+})(x-z)^{2}-k^{+}(x-w)^{2}\right]}\right]}"):
     r"\Big["+BL(LV)+rf"+\log\frac{{{A2}-{kv}{B2}}}{{\big(1+{lk}\big){A2}-{B2}}}\Big]",
 norm(r"\log{\left[\frac{\left[k^{+}(x-w)^{2}-\vee(x-z)^{2}\right]}{\left[k^{+}(x-w)^{2}-\Lambda(x-z)^{2}\right]}\right]}"):
     r"\Big["+BL(LV)+rf"+\log\frac{{{kv}{B2}-{A2}}}{{{B2}-{lk}{A2}}}\Big]",
 norm(r"\log{\left[\frac{(x-w)^{2}k^{+}-(x-z)^{2}\vee}{(x-w)^{2}k^{+}-(x-z)^{2}(k^{+}+\Lambda)}\right]}"):
     r"\Big["+BL(LV)+rf"+\log\frac{{{kv}{B2}-{A2}}}{{{B2}-\big(1+{lk}\big){A2}}}\Big]",
}
UNCHANGED=[  # logs that carry no large logarithm pointwise
 r"\log{\left[\frac{(x-w)^{2}\mu^{2}_{\overline{MS}}}{4}\right]}", r"\log{\left(1-\frac{k^{+}}{\vee}\right)}", r"\log^{2}{\left[1-\frac{k^{+}}{\vee}\right]}",
 r"\log{\left[1-\frac{k^{+}}{\vee}\right]}", r"\log{\left[\frac{\vee}{k^{+}}\right]}",   # the last one only occurs multiplied by log(1-k+/V) (Group I Row II)
 r"\log\frac{(w-y)^{2}}{(x-y)^{2}}", r"\log\frac{(w-x)^{2}}{(x-y)^{2}}",
 r"\log\left(\frac{p^{2}}{k^{2}+\frac{k^{+}}{\vee}p^{2}}\right)", r"\log\left[\frac{k^{2}}{(k-p)^{2}}+\frac{k^{+}}{\vee}\left(1-\frac{k^{2}}{(k-p)^{2}}\right)\right]",
 r"\log\left(1-\frac{k^{+}}{\vee}\right)",
 r"\log{\left[\frac{\left[(k^{+}-\Lambda)(x-z)^{2}-k^{+}(x-w)^{2}\right]}{\left[\Lambda(x-z)^{2}-k^{+}(x-w)^{2}\right]}\right]}",
 r"\log{\left[\frac{(x-w)^{2}k^{+}-(x-z)^{2}(k^{+}-\Lambda)}{(x-w)^{2}k^{+}-(x-z)^{2}\Lambda}\right]}",
 r"\log{\left[\frac{\left[k^{+}(x-w)-(k^{+}-\Lambda)(x-z)\right]^{2}}{\left[k^{+}(x-w)-\Lambda(x-z)\right]^{2}}\right]}",
 r"\log\frac{(y-x)^{2}\bar\mu^{2}_{\overline{MS}}}{4}",
 r"\log{\left[\frac{(\vee-k^{+})(x-z)^{2}+k^{+}(x-w)^{2}}{k^{+}(x'-w')^{2}+(\vee-k^{+})(x'-z)^{2}}\right]}",
 r"\log{\left[\frac{\Lambda(x-z)^{2}+k^{+}(x-w)^{2}}{k^{+}(x'-w')^{2}+\Lambda(x'-z)^{2}}\right]}",
 r"\log\frac{(y'-y)^{2}\mu^{2}_{\overline{MS}}}{4\xi^{2}}", r"\log\frac{(y'-y)^{2}\mu^{2}_{\overline{MS}}}{4}", r"\log\Big(1-\frac{k^{+}}{\vee}\Big)",
]
UNCHANGED={norm(x) for x in UNCHANGED}
stats={"blue":0,"other":0,"composite":0,"unchanged":0}
def rewrite(s):
    # 1. blue groups -> placeholders
    reps=[]
    out=[]; pos=0
    for st,e,inner in blue_groups(s):
        k=norm(inner); assert k in BLUE,("unknown blue group",k)
        out.append(s[pos:st]); out.append(f"@@{len(reps)}@@"); reps.append(BLUE[k]); pos=e; stats["blue"]+=1
    out.append(s[pos:]); t="".join(out)
    # 2. remaining logs
    out=[]; pos=0
    for st,e,x in log_spans(t):
        k=norm(x)
        if k in OTHER: new=OTHER[k]; stats["other"]+=1
        elif composite(k): new=composite(k); stats["composite"]+=1
        else:
            assert k in UNCHANGED,("unclassified log",k); new=None; stats["unchanged"]+=1
        if new is not None:
            out.append(t[pos:st]); out.append(f"@@{len(reps)}@@"); reps.append(new); pos=e
    out.append(t[pos:]); t="".join(out)
    for i,r in enumerate(reps): t=t.replace(f"@@{i}@@",r)
    return t
if __name__=="__main__":
    import json
    res={}
    for g,title,spans in BLOCKS:
        res[f"{g}|{title}"]=[rewrite(get(g,a,b)) for a,b in spans]
    json.dump(res,open("split_blocks.json","w"))
    print(stats)
