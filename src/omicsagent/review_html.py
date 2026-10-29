"""Self-contained, offline scientist review with a run-scoped local draft."""

from __future__ import annotations

import json
from typing import Any


def render_review(report: dict[str, Any]) -> str:
    """Render saved evidence without modifying the analysis or its original labels."""
    data = json.dumps(report, ensure_ascii=True, allow_nan=False, separators=(",", ":"))
    for char, escaped in (("&", "\\u0026"), ("<", "\\u003c"), (">", "\\u003e")):
        data = data.replace(char, escaped)
    return _PAGE.replace("__ENGINE__", _ENGINE).replace("__REPORT_JSON__", data)


_ENGINE = r"""
'use strict';
var ReviewEngine = (() => {
  function create(report) {
    if (report.schema_version !== 1 || typeof report.run_hash !== 'string' || !report.run_hash) throw Error('Unsupported report');
    const proposals = new Map(report.cases.map(c => [c.id, c]));
    if (proposals.size !== report.cases.length) throw Error('Duplicate case IDs');
    let events = [];
    let states = new Map([...proposals.keys()].map(id => [id, 'original']));
    const copy = value => JSON.parse(JSON.stringify(value));
    const ensure = (condition, message) => {if (!condition) throw Error(message);};
    function replay(document) {
      ensure(document && document.schema_version === 1, 'Unsupported review schema');
      ensure(document.run_hash === report.run_hash, 'Review belongs to a different analysis run');
      ensure(Array.isArray(document.events), 'Review events must be an array');
      const next = new Map([...proposals.keys()].map(id => [id, 'original']));
      let last = -Infinity;
      document.events.forEach((event, i) => {
        ensure(event && typeof event === 'object', 'Invalid review event');
        ensure(event.run_hash === report.run_hash, 'Event run mismatch');
        ensure(event.sequence === i + 1, 'Event sequence is not contiguous');
        ensure(proposals.has(event.proposal_id), 'Unknown proposal ID');
        ensure(event.action === 'accept' || event.action === 'undo', 'Unknown review action');
        ensure(typeof event.timestamp === 'string' && /^\d{4}-\d\d-\d\dT\d\d:\d\d:\d\d\.\d{3}Z$/.test(event.timestamp), 'Invalid timestamp');
        const time = Date.parse(event.timestamp);
        ensure(Number.isFinite(time) && new Date(time).toISOString() === event.timestamp && time >= last, 'Invalid timestamp order');
        last = time;
        ensure(typeof event.reason === 'string', 'Reason must be text');
        const allowed = proposals.get(event.proposal_id).evidence_ids || [];
        ensure(Array.isArray(event.evidence_refs) && event.evidence_refs.every(ref => typeof ref === 'string' && allowed.includes(ref)), 'Unknown evidence reference');
        const prior = next.get(event.proposal_id);
        const target = event.action === 'accept' ? 'accepted' : 'original';
        ensure(prior !== target && event.prior_state === prior && event.new_state === target, 'Invalid review state transition');
        next.set(event.proposal_id, target);
      });
      return next;
    }
    function exportReview() {return {schema_version:1, run_hash:report.run_hash, events:copy(events)};}
    function importReview(value) {
      const document = copy(typeof value === 'string' ? JSON.parse(value) : value);
      const next = replay(document); // Validate the complete history before changing the draft.
      events = document.events;
      states = next;
    }
    function act(id, action, reason = '', evidence = [], timestamp) {
      const last = events.length ? Date.parse(events[events.length - 1].timestamp) : 0;
      const event = {sequence:events.length + 1, run_hash:report.run_hash, proposal_id:id, action,
        timestamp:timestamp || new Date(Math.max(Date.now(), last)).toISOString(), prior_state:states.get(id),
        new_state:action === 'accept' ? 'accepted' : 'original', reason, evidence_refs:evidence};
      const draft = exportReview();draft.events.push(event);importReview(draft);return copy(event);
    }
    function csv() {
      const quote = value => {
        let text = String(value == null ? '' : value);
        if (/^[\s\u0000-\u001f]*[=+@-]/.test(text) || /^[\t\r\n]/.test(text)) text = "'" + text;
        return '"' + text.replace(/"/g, '""') + '"';
      };
      const rows = [['cell_id','original_l1','reviewed_l1','reviewed_l2','proposal_id','review_state','run_hash']];
      for (const c of proposals.values()) for (const cell of c.cells || []) {
        const accepted = states.get(c.id) === 'accepted';
        rows.push([cell.cell_id,cell.original_l1,accepted ? cell.proposed_l1 : cell.original_l1,
          accepted ? cell.proposed_l2 : '',c.id,states.get(c.id),report.run_hash]);
      }
      return rows.map(row => row.map(quote).join(',')).join('\r\n') + '\r\n';
    }
    return {act, state:id => states.get(id), exportReview, importReview, csv};
  }
  return {create};
})();
"""


_PAGE = r"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<meta name="color-scheme" content="light"><title>OmicsAgent · Scientist review</title>
<style>
:root{--ink:#182d34;--muted:#526b73;--paper:#f3f5f2;--line:#d9e1dc;--teal:#176e68;--pale:#e4f0eb;--amber:#895314;--red:#8c4240}*{box-sizing:border-box}body{margin:0;background:var(--paper);color:var(--ink);font:15px/1.55 system-ui,-apple-system,Segoe UI,sans-serif}button,input,textarea{font:inherit}button,.file-button{cursor:pointer;border:1px solid var(--line);background:white;border-radius:7px;padding:9px 14px;color:var(--ink)}button:hover,.file-button:hover{border-color:var(--teal);background:#f0f6f3}button:disabled{opacity:.45;cursor:default}button:focus-visible,a:focus-visible,summary:focus-visible,textarea:focus-visible{outline:3px solid #d69d40;outline-offset:3px}h1,h2,h3,p{margin:0}h1{font-size:32px;letter-spacing:-1px;font-weight:650}h2{font-size:24px;letter-spacing:-.5px}h3{font-size:15px}a{color:var(--teal);text-underline-offset:3px}header{border-bottom:1px solid var(--line);background:#fff}.top{max-width:1260px;margin:auto;padding:22px 32px;display:flex;justify-content:space-between;align-items:center;gap:20px}.brand{font-weight:750;font-size:18px;letter-spacing:-.5px}.brand small{font-weight:450;font-size:13px;color:var(--muted);margin-left:14px}.eyebrow{text-transform:uppercase;font-size:11px;font-weight:750;letter-spacing:1.6px;color:var(--teal)}.intro{max-width:1260px;padding:32px 32px 22px;margin:auto;display:flex;justify-content:space-between;gap:30px;align-items:start}.intro p{max-width:690px;color:var(--muted);margin-top:8px}.run{font:11px ui-monospace,monospace;overflow-wrap:anywhere;max-width:230px;color:var(--muted);text-align:right}.layout{max-width:1260px;margin:0 auto;padding:0 32px 36px;display:grid;grid-template-columns:240px minmax(0,1fr);gap:24px}.cases{display:flex;flex-direction:column;gap:10px}.case-button{width:100%;text-align:left;padding:16px;display:block;background:#fbfcfa}.case-button[aria-pressed=true]{border:1px solid var(--teal);box-shadow:inset 3px 0 var(--teal);background:#eaf3ef}.case-button strong{display:block;margin:5px 0}.case-button small{display:block;color:var(--muted);font-size:12px}.kind{font-size:11px;font-weight:750;text-transform:uppercase;letter-spacing:1px}.conflict{color:var(--amber)}.agreement{color:var(--teal)}.insufficient{color:var(--red)}.coverage{font-size:12px;color:var(--muted);padding:14px 2px}.panel{border:1px solid var(--line);border-radius:12px;background:#fff;overflow:hidden}.case-head{padding:24px 26px;border-bottom:1px solid var(--line);display:flex;gap:20px;justify-content:space-between}.case-head p{color:var(--muted);font-size:13px;margin-top:5px}.badge{display:inline-block;background:var(--pale);border-radius:20px;padding:5px 10px;font-size:12px;height:fit-content;white-space:nowrap}.content{padding:24px 26px}.judgment{display:grid;grid-template-columns:1fr 1fr;gap:16px;margin:0 0 24px}.judgment>div{border:1px solid var(--line);border-radius:8px;padding:16px}.value{font-size:25px;letter-spacing:-.5px;font-weight:650;margin:3px 0}.muted{color:var(--muted);font-size:13px}.evidence{display:grid;grid-template-columns:1fr 1fr;gap:24px}.evidence h3{margin-bottom:8px}.evidence ul{padding-left:19px;margin:0;font-size:14px}.evidence li+li{margin-top:8px}.open-question{background:#f7f6f0;border-left:3px solid #c79a57;padding:16px 18px;margin-top:24px}.open-question p{margin-top:5px;font-size:14px}.next{margin-top:16px;border-top:1px solid #e5e0d1;padding-top:13px}.next strong{font-size:12px;text-transform:uppercase;letter-spacing:.8px}.next .observation{color:var(--muted);font-size:12px}.review{padding:22px 26px;background:#f6f9f6;border-top:1px solid var(--line)}.diff{display:flex;flex-wrap:wrap;gap:10px 26px;margin:12px 0;font-size:13px}.diff b{display:block;color:var(--muted);font-size:11px;text-transform:uppercase;letter-spacing:.6px}textarea{width:100%;min-height:52px;resize:vertical;border:1px solid var(--line);border-radius:6px;padding:9px;margin:8px 0;background:white}.primary{background:var(--teal);border-color:var(--teal);color:white}.primary:hover{background:#10534f;color:white}.actions{display:flex;gap:10px;align-items:center;flex-wrap:wrap}.note{font-size:12px;color:var(--muted);margin-top:10px}details{border-top:1px solid var(--line)}summary{cursor:pointer;padding:17px 26px;font-weight:600;font-size:14px}details>div{padding:0 26px 22px}.scroll{overflow:auto}table{width:100%;border-collapse:collapse;text-align:left;font-size:12px}th,td{padding:8px 10px;border-bottom:1px solid var(--line);vertical-align:top}th{font-weight:650;color:var(--muted);white-space:nowrap}td{overflow-wrap:anywhere}pre{white-space:pre-wrap;overflow-wrap:anywhere;font-size:12px;background:#f5f7f4;padding:14px;border-radius:6px}.toolbar{display:flex;gap:8px;flex-wrap:wrap}.file-button{display:inline-block;font-size:13px}.toolbar button{font-size:13px}.file-input{position:absolute;width:1px;height:1px;overflow:hidden;clip:rect(0,0,0,0)}#status{padding:12px 16px;border:1px solid var(--line);border-radius:7px;margin-top:12px;font-size:13px}#status:empty{display:none}.footer{max-width:1260px;margin:auto;padding:0 32px 30px;color:var(--muted);font-size:12px}.source-list{display:flex;flex-wrap:wrap;gap:10px 20px;padding:5px 0 15px}.plot{min-width:560px;max-width:100%;height:auto;display:block}.legend{font-size:12px;color:var(--muted);margin-bottom:12px}.history-empty{color:var(--muted);font-size:13px}.source-path{font:11px ui-monospace,monospace}.selection-warning{color:var(--red)}@media(max-width:800px){.top{align-items:start;flex-direction:column}.intro{padding-top:22px}.intro .run{display:none}.layout{grid-template-columns:1fr;gap:14px}.cases{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:7px}.case-button{padding:11px}.case-button strong{font-size:13px}.case-button small{font-size:10px}.coverage{grid-column:1/-1;padding:3px}.judgment,.evidence{gap:12px}.case-head,.content,.review{padding:20px}.case-head{flex-wrap:wrap}.brand small{display:block;margin-left:0}.intro,.layout,.top{padding-left:18px;padding-right:18px}h1{font-size:28px}}@media(max-width:480px){.evidence,.judgment{grid-template-columns:1fr}.cases{grid-template-columns:1fr}.case-button small{display:inline}.case-button strong{display:inline;margin:0 10px}.case-button .kind{display:block}h2{font-size:21px}}@media print{button,.toolbar,textarea{display:none}.layout{display:block}.cases{display:none}.panel{border:0}body{background:white}details{break-inside:avoid}}
</style></head><body>
<header><div class="top"><div class="brand">OmicsAgent <small>Evidence-aware annotation review</small></div><div class="toolbar"><button id="export-json">Export review JSON</button><label class="file-button" for="import-json">Import review JSON</label><input class="file-input" id="import-json" type="file" accept=".json,application/json"><button id="export-csv">Export reviewed cells CSV</button></div></div></header>
<section class="intro"><div><div class="eyebrow">Scientist workspace · saved evidence</div><h1>Resolve what the evidence supports.</h1><p>Compare the annotation proposal with its measured support, inspect the disagreement, and record a review choice.</p><div id="status" role="status" aria-live="polite"></div></div><div class="run" id="run-hash"></div></section>
<main class="layout"><nav class="cases" id="case-selectors" aria-label="Review cases"></nav><section class="panel" id="case-panel" aria-label="Selected case"></section></main>
<footer class="footer">Review state is a local draft in this open page. Export JSON to keep or transfer the history; reopening the report starts from the original labels. Accept and undo never modify analysis files. CSV contains only the selected case cells and starts from each cell’s original L1; subtype stays blank until acceptance. Published annotations are reference metadata, not independent ground truth.</footer>
<script id="report-data" type="application/json">__REPORT_JSON__</script>
<script id="review-engine">__ENGINE__</script>
<script id="review-ui">
'use strict';
(() => {
  const report = JSON.parse(document.getElementById('report-data').textContent);
  const engine = ReviewEngine.create(report);
  let selected = (report.cases.find(c => c.kind === 'conflict') || report.cases[0] || {}).id;
  const byId = id => document.getElementById(id);
  const el = (tag, text, cls) => {const n = document.createElement(tag);if(text !== undefined)n.textContent=String(text);if(cls)n.className=cls;return n;};
  const kinds = {agreement:'Agreement',conflict:'Evidence conflict',insufficient:'Insufficient evidence'};
  function status(text) {byId('status').textContent=text;}
  function table(rows) {
    if (!rows || !rows.length) return el('p','No measurements available.','muted');
    const columns=[...new Set(rows.flatMap(row=>Object.keys(row)))];
    const wrap=el('div',undefined,'scroll'),t=el('table'),head=el('thead'),tr=el('tr');
    columns.forEach(k=>tr.append(el('th',k.replaceAll('_',' '))));head.append(tr);t.append(head);
    const body=el('tbody');rows.forEach(row=>{const r=el('tr');columns.forEach(k=>r.append(el('td',row[k] == null ? '—' : typeof row[k] === 'object' ? JSON.stringify(row[k]) : row[k])));body.append(r);});t.append(body);wrap.append(t);return wrap;
  }
  function details(title, child, id) {const d=el('details');if(id)d.id=id;d.append(el('summary',title));const body=el('div');body.append(child);d.append(body);return d;}
  function plot(markers) {
    const wrap=el('div');wrap.append(el('p','Dot area = fraction detected (0–100%); color = mean log expression, scaled within this case. Missing measurements are omitted.','legend'));
    const rows=markers.filter(m=>Number.isFinite(m.mean_log_expression)&&Number.isFinite(m.fraction_detected));
    if (!rows.length) {wrap.append(el('p','Marker measurements unavailable.','muted'));return wrap;}
    const ns='http://www.w3.org/2000/svg',make=(tag,attrs,text)=>{const n=document.createElementNS(ns,tag);Object.entries(attrs).forEach(([k,v])=>n.setAttribute(k,String(v)));if(text!==undefined)n.textContent=text;return n;};
    const genes=[...new Set(rows.map(m=>m.gene))],labels=[...new Set(rows.map(m=>m.level+' · '+m.label))];
    const width=Math.max(560,190+genes.length*54),height=95+labels.length*42;
    const svg=make('svg',{viewBox:`0 0 ${width} ${height}`,role:'img','aria-label':'Marker expression and detection dotplot',class:'plot'});
    svg.append(make('title',{},'Marker expression and detection by label'));
    const max=Math.max(...rows.map(m=>m.mean_log_expression),.01);
    genes.forEach((gene,i)=>svg.append(make('text',{x:190+i*54,y:26,'text-anchor':'middle','font-size':11,fill:'#526b73'},gene)));
    labels.forEach((label,i)=>{svg.append(make('text',{x:170,y:57+i*42,'text-anchor':'end','font-size':11,fill:'#182d34'},label));svg.append(make('line',{x1:180,x2:width-15,y1:53+i*42,y2:53+i*42,stroke:'#edf1ed'}));});
    rows.forEach(m=>{const ratio=Math.max(0,Math.min(1,m.mean_log_expression/max)),circle=make('circle',{cx:190+genes.indexOf(m.gene)*54,cy:53+labels.indexOf(m.level+' · '+m.label)*42,r:Math.sqrt(Math.max(0,Math.min(1,m.fraction_detected)))*15,fill:`hsl(173 48% ${86-ratio*58}%)`,stroke:'#236b63','stroke-width':.5});circle.append(make('title',{},`${m.gene}; ${m.label}; mean log expression ${m.mean_log_expression}; fraction detected ${m.fraction_detected}`));svg.append(circle);});
    const scroll=el('div',undefined,'scroll');scroll.append(svg);wrap.append(scroll);return wrap;
  }
  function safeLink(path) {
    // Only evidence files within the exported report folder are navigable.
    if(typeof path!=='string'||!path||/[\u0000-\u0020\\:?#]/.test(path)||path.startsWith('/'))return null;
    try {if(decodeURIComponent(path).split('/').some(p=>p==='..'||p==='.')||/[\\:?#\u0000-\u0020]/.test(decodeURIComponent(path)))return null;}catch{return null;}
    return path;
  }
  function render() {
    const nav=byId('case-selectors');nav.replaceChildren();
    report.cases.forEach(c=>{const b=el('button',undefined,'case-button');b.id='select-'+c.id;b.setAttribute('aria-pressed',String(c.id===selected));b.append(el('span',kinds[c.kind]||c.kind,'kind '+c.kind),el('strong',c.l1),el('small',`Iteration ${c.iteration} · Cluster ${c.cluster} · ${c.n_cells} cells`));b.addEventListener('click',()=>{selected=c.id;render();});nav.append(b);});
    const coverage=el('div',undefined,'coverage');coverage.append(el('div',`Selected ${report.cases.length} cases · eligible: `+Object.entries(report.selection.eligible||{}).map(([k,n])=>`${k} ${n}`).join(', ')));
    if(!report.selection.complete)coverage.append(el('p','Missing measured case types: '+(report.selection.missing||[]).join(', '),'selection-warning'));nav.append(coverage);
    const panel=byId('case-panel');panel.replaceChildren();const c=report.cases.find(c=>c.id===selected);
    if(!c){const empty=el('div',undefined,'content');empty.append(el('h2','Final round remains unresolved'),el('p','No eligible cases were found. Earlier judgments are not reused.'));if((report.selection.skipped_pools||[]).length)empty.append(table(report.selection.skipped_pools));panel.append(empty);return;}
    const accepted=engine.state(c.id)==='accepted',head=el('div',undefined,'case-head'),heading=el('div');heading.append(el('div',kinds[c.kind]||c.kind,'eyebrow '+c.kind),el('h2',`Cluster ${c.cluster} · iteration ${c.iteration}`),el('p',`${c.n_cells} affected cells · ${c.id}`));head.append(heading,el('span',accepted?'Proposal accepted in draft':'Original labels retained','badge'));panel.append(head);
    const content=el('div',undefined,'content'),judgment=el('div',undefined,'judgment');
    [[ 'Coarse identity · L1',c.l1,c.l1_supported?'Supported by measured evidence':'Unresolved; retaining a label does not validate it'],['Subtype · L2',c.l2,c.l2==='Unknown'?'Deferred; inspect the unresolved question':'Provisional subtype judgment']].forEach(([label,value,note])=>{const box=el('div');box.append(el('div',label,'eyebrow'),el('p',value,'value'),el('p',note,'muted'));judgment.append(box);});content.append(judgment);
    const evidence=el('div',undefined,'evidence');[['Supporting evidence',c.support],['Counterevidence / limitations',c.counterevidence]].forEach(([title,items])=>{const box=el('div');box.append(el('h3',title));const ul=el('ul');(items&&items.length?items:['No measured evidence available.']).forEach(item=>ul.append(el('li',item)));box.append(ul);evidence.append(box);});content.append(evidence);
    const question=el('div',undefined,'open-question');question.append(el('h3','What remains unresolved'),el('p',c.unresolved||'No unresolved question was recorded.'));
    if(c.next_check){const names={donor:'Donor-stratified markers',cross_donor_conflict:'Query/reference subtype comparison',competing_markers:'Competing subtype markers',coverage:'Marker coverage and detection',subtype_separation:'Subtype marker co-expression'};const next=el('div',undefined,'next');next.append(el('strong','Prioritized next check · '+(names[c.next_check.kind]||'Evidence inspection')),el('p',c.next_check.reason),el('p',c.next_check.observation,'observation'));question.append(next);}content.append(question);panel.append(content);
    const review=el('div',undefined,'review');review.append(el('h3','Review the proposal'));
    const diff=el('div',undefined,'diff');[['Original L1',c.original_l1],['Proposed L1 / L2',c.proposed_l1+' / '+c.l2],['Reviewed L1',accepted?c.proposed_l1:c.original_l1]].forEach(([k,v])=>{const x=el('div');x.append(el('b',k),el('span',v));diff.append(x);});review.append(diff);
    const label=el('label','Review reason (optional)');label.htmlFor='review-reason';const reason=el('textarea');reason.id='review-reason';reason.placeholder='Record the evidence behind your choice';review.append(label,reason);
    const actions=el('div',undefined,'actions'),accept=el('button','Accept proposal','primary'),undo=el('button','Undo acceptance');accept.id='accept-proposal';undo.id='undo-proposal';accept.disabled=accepted;undo.disabled=!accepted;
    function act(action){try{engine.act(c.id,action,reason.value,c.evidence_ids||[]);render();status(action==='accept'?'Proposal accepted in the local draft. Export JSON to save the review.':'Acceptance undone in the local draft. Original labels restored.');}catch(error){status(error.message);}}
    accept.addEventListener('click',()=>act('accept'));undo.addEventListener('click',()=>act('undo'));actions.append(accept,undo);review.append(actions,el('p','Applies to the '+(c.cells||[]).length+' saved member cells. Accepting preserves any Unknown subtype and does not resolve conflicting evidence.','note'));panel.append(review);
    panel.append(details('Marker evidence · dotplot and measurements',(()=>{const box=plot(c.markers||[]);box.append(table(c.markers||[]));return box;})(),'marker-evidence'));
    panel.append(details('Donor evidence · measured denominators',table(c.donor_summary||[]),'donor-evidence'));
    panel.append(details('Affected cells · original → proposed → reviewed',table((c.cells||[]).map(cell=>({cell_id:cell.cell_id,original_l1:cell.original_l1,proposed_l1:cell.proposed_l1,proposed_l2:cell.proposed_l2,reviewed_l1:accepted?cell.proposed_l1:cell.original_l1,reviewed_l2:accepted?cell.proposed_l2:''}))),'affected-cells'));
    const params=el('div');params.append(el('h3','Recorded decision'),el('pre',JSON.stringify(c.decision,null,2)),el('h3','Thresholds'),el('pre',JSON.stringify(report.thresholds,null,2)));panel.append(details('Decision details and analysis parameters',params,'parameters'));
    const sources=el('div'),links=el('div',undefined,'source-list');Object.entries(c.links||{}).forEach(([label,path])=>{const safe=safeLink(path);if(safe){const a=el('a',label);a.setAttribute('href',safe);links.append(a);}else links.append(el('span',label+' (path unavailable)'));});sources.append(links,el('h3','Run provenance'),el('pre',JSON.stringify(report.provenance,null,2)),el('h3','Source file checksums'),table(Object.entries(report.sources||{}).map(([path,sha])=>({path,sha256:sha}))),el('p','Evidence IDs: '+(c.evidence_ids||[]).join(', '),'source-path'));panel.append(details('Source artifacts and provenance',sources,'sources'));
    const history=engine.exportReview().events.filter(e=>e.proposal_id===c.id);panel.append(details('Review history · '+history.length+' events',history.length?table(history):el('p','No review actions recorded for this case.','history-empty'),'review-history'));
  }
  function download(name,data,type){const blob=new Blob([data],{type}),url=URL.createObjectURL(blob),a=el('a');a.href=url;a.download=name;document.body.append(a);a.click();a.remove();setTimeout(()=>URL.revokeObjectURL(url),1000);}
  byId('run-hash').textContent='ANALYSIS RUN\n'+report.run_hash;
  byId('export-json').addEventListener('click',()=>{download('omicsagent-review.json',JSON.stringify(engine.exportReview(),null,2),'application/json');status('Review JSON exported. It includes the analysis run and complete action history.');});
  byId('export-csv').addEventListener('click',()=>{download('omicsagent-reviewed-cells.csv',engine.csv(),'text/csv;charset=utf-8');status('Reviewed cell labels exported. Formula-like identifiers and labels are prefixed with an apostrophe for spreadsheet safety.');});
  byId('import-json').addEventListener('change',async event=>{try{const file=event.target.files[0];if(!file)return;engine.importReview(await file.text());render();status('Review imported. Every event was validated against this run; labels were recomputed from its history.');}catch(error){status('Import rejected; draft unchanged. '+error.message);}finally{event.target.value='';}});
  render();
})();
</script></body></html>
"""
