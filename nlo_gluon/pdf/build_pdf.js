// Render decomposition.html (KaTeX) to ../nlo_gluon_decomposition.pdf with headless Chromium.
const path = require('path');
const { chromium } = require('playwright');
(async () => {
  const b = await chromium.launch();
  const pg = await b.newPage();
  const errs = [];
  pg.on('pageerror', e => errs.push(String(e)));
  await pg.goto('file://' + path.join(__dirname, 'decomposition.html'));
  await pg.waitForFunction("document.body.dataset.ready === '1'", null, { timeout: 30000 });
  await pg.waitForTimeout(1000);
  const bad = await pg.$$eval('.katex-error', els => els.map(e => e.textContent));
  const out = path.join(__dirname, '..', 'nlo_gluon_decomposition.pdf');
  await pg.pdf({ path: out, format: 'A4', printBackground: true,
    margin: { top: '18mm', bottom: '18mm', left: '15mm', right: '15mm' },
    displayHeaderFooter: true, headerTemplate: '<div></div>',
    footerTemplate: '<div style="width:100%;font-size:8pt;color:#777;font-family:Times New Roman,serif;text-align:center;"><span class="pageNumber"></span> / <span class="totalPages"></span></div>' });
  await b.close();
  console.log('page errors:', errs.length ? errs : 'none');
  console.log('katex errors:', bad.length ? bad : 'none');
  console.log('written:', out);
})();
