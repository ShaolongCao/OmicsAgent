/* Run with node tests/review_ui.cjs. Exercises the actual embedded review engine. */
const assert = require('node:assert/strict');
const {spawnSync} = require('node:child_process');
const path = require('node:path');
const vm = require('node:vm');
const root = path.resolve(__dirname, '..');
const fixture = {schema_version:1,run_hash:'run-a',provenance:{},selection:{eligible:{conflict:1},complete:false,missing:['agreement']},thresholds:{},sources:{},cases:[{
  id:'case-b',kind:'conflict',iteration:1,cluster:'b',n_cells:1,l1:'T',l2:'Unknown',l1_supported:true,
  original_l1:'old',proposed_l1:'T',support:['support'],counterevidence:['opposition'],unresolved:'subtype',next_check:null,decision:{},markers:[],donor_summary:[],links:{},evidence_ids:['decisions:1:b'],cells:[{cell_id:'=danger',original_l1:'old',proposed_l1:'T',proposed_l2:'Unknown'}]
}]};
fixture.provenance.test = '</script><img src=x onerror=alert(1)>& __ENGINE__';
fixture.cases[0].markers = [{level:'l1',label:'T',gene:'CD3D',mean_log_expression:2,fraction_detected:.8}];
fixture.cases[0].links = {decision:'hierarchy/decisions.csv',unsafe:'javascript:alert(1)',encoded:'%6aavascript%3aalert(1)'};
fixture.cases.push({...fixture.cases[0],id:'case-a',kind:'agreement',l2:'CD4',cells:[]});
fixture.cases.push({...fixture.cases[0],id:'case-c',kind:'insufficient',l1_supported:false,cells:[]});
const py = process.env.PYTHON || path.join(root,'.venv','Scripts','python.exe');
const result = spawnSync(py,['-c',"import json,sys;sys.path.insert(0,'src');from omicsagent.review_html import render_review;print(render_review(json.load(sys.stdin)))"],{cwd:root,input:JSON.stringify(fixture),encoding:'utf8',env:{...process.env,PYTHONIOENCODING:'utf-8'}});
assert.equal(result.status,0,result.stderr);
const html = result.stdout;
assert(!html.includes(fixture.provenance.test),'embedded JSON must not terminate script');
const embedded = JSON.parse(html.match(/<script id="report-data" type="application\/json">([\s\S]*?)<\/script>/)[1]);
assert.equal(embedded.provenance.test,fixture.provenance.test);
const context = {console};
vm.createContext(context);
vm.runInContext(html.match(/<script id="review-engine">([\s\S]*?)<\/script>/)[1],context);
const engine = context.ReviewEngine.create(fixture);
assert.equal(engine.state('case-b'),'original');
const first = engine.act('case-b','accept','markers checked',['decisions:1:b'],'2026-10-28T12:00:00.000Z');
assert.equal(first.sequence,1);
assert.equal(first.run_hash,'run-a');
assert.equal(engine.state('case-b'),'accepted');
assert.match(engine.csv(),/"'=danger","old","T","Unknown","case-b"/);
const saved = engine.exportReview();
engine.act('case-b','undo','recheck',[],'2026-10-28T12:01:00.000Z');
assert.equal(engine.state('case-b'),'original');
assert.match(engine.csv(),/"'=danger","old","old","","case-b"/);
engine.importReview(saved);
assert.equal(engine.state('case-b'),'accepted');
assert.equal(engine.exportReview().events.length,1);
const corruptions = [
  d=>{d.run_hash='other'},d=>{d.events[0].proposal_id='unknown'},
  d=>{d.events[0].action='mystery'},d=>{d.events[0].sequence=2},
  d=>{d.events[0].timestamp='tomorrow'},d=>{d.events[0].prior_state='accepted'},
  d=>{d.events[0].new_state='original'},d=>{d.events[0].evidence_refs=['invented']},
  d=>{d.events[0].run_hash='other'},d=>{d.schema_version=9},
  d=>{d.events[0].reason={}},d=>{d.events.push({...d.events[0],sequence:2})},
];
for (const corrupt of corruptions) {
  const bad = JSON.parse(JSON.stringify(saved));corrupt(bad);
  assert.throws(()=>engine.importReview(bad));
  assert.equal(engine.state('case-b'),'accepted','import rejection must be atomic');
  assert.equal(engine.exportReview().events.length,1);
}
const forged = JSON.parse(JSON.stringify(saved));forged.state={'case-b':'original'};
engine.importReview(forged);
assert.equal(engine.state('case-b'),'accepted','state must be recomputed from events');
assert.throws(()=>engine.act('case-b','accept','',[]));
assert.equal(fixture.cases[0].cells[0].original_l1,'old');
console.log('PASS: safe embedding, accept/undo, JSON roundtrip, 12 atomic rejections, forged-state rejection, original-label CSV and injection protection');

if (process.env.REVIEW_BROWSER === '1') {
  (async () => {
    const fs = require('node:fs');
    const os = require('node:os');
    const {pathToFileURL} = require('node:url');
    const {chromium} = require(process.env.PLAYWRIGHT_MODULE || 'playwright');
    const dir = fs.mkdtempSync(path.join(os.tmpdir(),'omicsagent-review-'));
    const filename = path.join(dir,'review.html');fs.writeFileSync(filename,html);
    const browser = await chromium.launch({headless:true,...(process.env.PLAYWRIGHT_CHANNEL?{channel:process.env.PLAYWRIGHT_CHANNEL}:{})});
    try {
      const page = await browser.newPage({viewport:{width:1280,height:1000},acceptDownloads:true});
      const errors=[];page.on('pageerror',error=>errors.push(error.message));
      await page.goto(pathToFileURL(filename).href);
      assert.equal(await page.locator('#select-case-b').getAttribute('aria-pressed'),'true');
      await page.locator('#select-case-c').click();
      assert.match(await page.locator('#case-panel').textContent(),/Unresolved; retaining a label does not validate it/);
      await page.locator('#select-case-b').click();
      await page.locator('#marker-evidence summary').click();
      assert.equal(await page.locator('#marker-evidence svg circle').count(),1);
      await page.locator('#sources summary').click();
      assert.equal(await page.locator('#sources a').count(),1,'unsafe source URLs are not clickable');
      assert.equal(await page.locator('#sources img').count(),0,'untrusted strings remain text');
      await page.locator('#review-reason').fill('Checked CD3D and retained uncertainty');
      await page.locator('#accept-proposal').click();
      assert.equal(await page.locator('#accept-proposal').isDisabled(),true);
      const jsonDownload = page.waitForEvent('download');await page.locator('#export-json').click();
      const jsonPath=await (await jsonDownload).path();
      const exported=JSON.parse(fs.readFileSync(jsonPath,'utf8'));
      assert.equal(exported.events[0].reason,'Checked CD3D and retained uncertainty');
      const csvDownload=page.waitForEvent('download');await page.locator('#export-csv').click();
      const csvPath=await (await csvDownload).path();
      assert.match(fs.readFileSync(csvPath,'utf8'),/"'=danger","old","T","Unknown","case-b"/);
      await page.locator('#undo-proposal').click();
      assert.equal(await page.locator('#accept-proposal').isDisabled(),false);
      await page.locator('#import-json').setInputFiles({name:'review.json',mimeType:'application/json',buffer:Buffer.from(JSON.stringify(exported))});
      await page.waitForFunction(()=>document.getElementById('status').textContent.startsWith('Review imported'));
      assert.equal(await page.locator('#accept-proposal').isDisabled(),true);
      exported.events[0].action='unknown';
      await page.locator('#import-json').setInputFiles({name:'bad.json',mimeType:'application/json',buffer:Buffer.from(JSON.stringify(exported))});
      await page.waitForFunction(()=>document.getElementById('status').textContent.startsWith('Import rejected'));
      assert.equal(await page.locator('#accept-proposal').isDisabled(),true);
      await page.setViewportSize({width:390,height:844});
      assert.equal(await page.evaluate(()=>document.documentElement.scrollWidth<=window.innerWidth),true,'mobile viewport must not overflow');
      if(process.env.REVIEW_SCREENSHOT)await page.screenshot({path:process.env.REVIEW_SCREENSHOT,fullPage:true});
      assert.deepEqual(errors,[]);
      await page.reload();
      assert.equal(await page.locator('#accept-proposal').isDisabled(),false,'reopening starts a fresh original draft');
      console.log('PASS browser: file URL, conflict default, case switching, dotplot expansion, safe sources, accept/undo, downloaded JSON/CSV, atomic import, mobile layout and fresh reload');
    } finally {await browser.close();fs.rmSync(dir,{recursive:true,force:true});}
  })().catch(error=>{console.error(error);process.exitCode=1;});
}
