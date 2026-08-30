import pathlib
from playwright.sync_api import sync_playwright
CHROME = "/opt/pw-browsers/chromium-1194/chrome-linux/chrome"
out = str(pathlib.Path(__file__).resolve().parent.parent / "z_integration_UV_finite_regular.pdf")
with sync_playwright() as p:
    b = p.chromium.launch(executable_path=CHROME, args=["--no-sandbox"])
    pg = b.new_page()
    errs = []
    pg.on("pageerror", lambda e: errs.append(str(e)))
    pg.on("console", lambda m: errs.append(m.type + ": " + m.text) if m.type == "error" else None)
    pg.goto("file://" + str(pathlib.Path(__file__).resolve().parent / "document.html"))
    pg.wait_for_function("document.title === 'ready'", timeout=20000)
    pg.wait_for_timeout(800)
    # any TeX that failed to parse is left in KaTeX error colour
    bad = pg.eval_on_selector_all(".katex-error", "els => els.map(e => e.textContent)")
    raw = pg.evaluate("""() => {
        const t = document.body.innerText;
        const m = t.match(/\\$\\$?[^$]{0,60}/g);
        return m ? m.slice(0,10) : [];
    }""")
    pg.pdf(path=out, format="A4", print_background=True,
           margin={"top":"20mm","bottom":"20mm","left":"18mm","right":"18mm"},
           display_header_footer=True,
           header_template="<div></div>",
           footer_template='<div style="width:100%;font-size:8pt;color:#777;'
                           'font-family:Times New Roman,serif;text-align:center;">'
                           '<span class="pageNumber"></span> / <span class="totalPages"></span></div>')
    b.close()
print("page errors:", errs or "none")
print("katex parse errors:", bad or "none")
print("unrendered $ fragments:", raw or "none")
