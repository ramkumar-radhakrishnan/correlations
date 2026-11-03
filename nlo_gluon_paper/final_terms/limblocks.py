import re
def blocks(fn):
    t=open(fn).read()
    t=t[t.index(r'\begin{document}'):]
    t=re.sub(r'\\begin\{comment\}.*?\\end\{comment\}','',t,flags=re.S)
    t='\n'.join(l for l in t.split('\n') if not l.lstrip().startswith('%'))
    out=[];path={}
    for m in re.finditer(r'\\(section|subsection|subsubsection)\*?\{([^}]*)\}|\\begin\{align\}(.*?)\\end\{align\}',t,re.S):
        if m.group(1):
            path[m.group(1)]=m.group(2)
            if m.group(1)=='section': path.pop('subsection',None);path.pop('subsubsection',None)
            if m.group(1)=='subsection': path.pop('subsubsection',None)
        else:
            out.append((' / '.join(path.get(k,'') for k in ('section','subsection','subsubsection')), m.group(3)))
    return out
def toks(s):
    s=re.sub(r'\\notag|\\\\|&|\\,|\;|\\!|\\left|\\right|\\bigg|\\Bigg|\\big|\\Big','',s)
    return re.findall(r'\\[A-Za-z]+|[A-Za-z0-9]|\S',s)
