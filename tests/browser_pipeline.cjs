// Inspect one explicitly selected existing real pipeline run. Never starts tools or calls a model.
const {chromium} = require('../.local/web-qa/node_modules/playwright');
const assert = require('node:assert/strict');
const fs = require('node:fs/promises');
const path = require('node:path');
const crypto = require('node:crypto');

(async()=>{
  const origin=process.env.CUTTAG_TEST_URL || 'http://127.0.0.1:8766';
  const id=process.env.CUTTAG_PIPELINE_REUSE_RUN;
  assert.match(id || '',/^run-[0-9a-f]{12}$/,'Set CUTTAG_PIPELINE_REUSE_RUN to an existing real pipeline case');
  const detailUrl=origin+'/#history/'+id;
  const browser=await chromium.launch({executablePath:process.env.CUTTAG_BROWSER_PATH,headless:true});
  const folder=path.resolve(__dirname,'../outputs/preview');
  await fs.mkdir(folder,{recursive:true});
  const errors=[],external=[];
  try {
    const page=await browser.newPage({viewport:{width:1440,height:1100},reducedMotion:'reduce'});
    page.on('pageerror',e=>errors.push(e.message));
    page.on('request',r=>{if(new URL(r.url()).origin!==origin)external.push(r.url());});
    await page.goto(detailUrl);
    const detail=page.locator('#history-run-content');
    await page.waitForFunction(()=>['基础分析完成','任务失败','已超时','已停止','已中断'].includes(document.querySelector('#history-run-label')?.textContent),null,{timeout:7200000});
    const record=await(await page.request.get(`${origin}/api/runs/${id}`)).json();
    assert.equal(record.status,'completed_pipeline',JSON.stringify({status:record.status,error:record.error}));
    assert.equal(record.result.coverage.scope,'whole_genome');
    assert.equal(record.actual_model_requests,0);
    assert.equal(record.execution.inputs_unchanged,true);
    assert.ok(await detail.locator('.pipeline-metrics tbody tr').count() >= 4);
    assert.match(await detail.innerText(),/未使用 IgG 对照/);
    for (const name of ['coverage.cpm.bw','narrow.peaks.bed','broad.peaks.bed','tss.matrix.gz','tss.profile.png','tss.heatmap.png']) {
      const artifact=record.artifacts.find(a=>a.name===name);
      assert.ok(artifact,name);
      const response=await page.request.get(`${origin}/api/runs/${id}/artifacts/${name}`);
      assert.equal(response.status(),200);
      const body=await response.body();
      assert.equal(body.length,artifact.bytes);
      assert.equal(crypto.createHash('sha256').update(body).digest('hex'),artifact.sha256);
    }
    const bw=await page.request.get(`${origin}/api/runs/${id}/artifacts/coverage.cpm.bw`,{headers:{Range:'bytes=0-3'}});
    assert.equal(bw.status(),206);
    assert.equal((await bw.body()).toString('hex'),'26fc8f88');
    assert.equal((await page.request.get(`${origin}/api/runs/${id}/artifacts/analysis.bam`)).status(),404);
    assert.equal((await page.request.get(`${origin}/api/runs/${id}/artifacts/.env`)).status(),404);
    for(const width of [1440,1280,1024,768,390,320]) {
      await page.setViewportSize({width,height:1000});
      await detail.locator('.result-plot').evaluateAll(images=>Promise.all(images.map(img=>img.decode())));
      assert.equal(await page.evaluate(()=>document.documentElement.scrollWidth>innerWidth),false,`Overflow ${width}`);
      await detail.locator('details[data-detail="tss"] summary').click();
      assert.equal(await page.evaluate(()=>document.documentElement.scrollWidth>innerWidth),false,`Expanded overflow ${width}`);
      await detail.locator('details[data-detail="tss"] summary').click();
    }
    await page.setViewportSize({width:1440,height:1100});
    await page.evaluate(()=>scrollTo(0,0));
    await page.screenshot({path:path.join(folder,'pipeline-desktop.png'),fullPage:true});
    const mobile=await browser.newPage({viewport:{width:390,height:844},reducedMotion:'reduce'});
    mobile.on('pageerror',e=>errors.push(e.message));
    await mobile.goto(detailUrl);
    await mobile.locator('#history-run-content .pipeline-metrics').waitFor();
    await mobile.locator('#history-run-content .result-plot').evaluateAll(images=>Promise.all(images.map(img=>img.decode())));
    await mobile.evaluate(()=>scrollTo(0,0));
    await mobile.screenshot({path:path.join(folder,'pipeline-mobile.png'),fullPage:true});
    await page.reload();
    await detail.locator('.pipeline-metrics').waitFor();
    await detail.locator('.result-plot').evaluateAll(images=>Promise.all(images.map(img=>img.decode())));
    assert.deepEqual(errors,[]);assert.deepEqual(external,[]);
    const evidence={status:'passed',run:id,source:record.source,real_tool_executions:record.actual_tool_executions,
      actual_model_requests:0,scope:'whole_genome',input_metadata_unchanged:record.execution.inputs_unchanged,
      downloads:'BW / both BED / matrix / plots: SHA256 matched',range_download:'passed',
      widths:[1440,1280,1024,768,390,320],browser_errors:errors,external_browser_resources:external,
      scientific_quality:'not_assessed'};
    await fs.writeFile(path.join(folder,'pipeline-validation.json'),JSON.stringify(evidence,null,2)+'\n');
    console.log(JSON.stringify(evidence));
  } finally {await browser.close();}
})().catch(e=>{console.error(e);process.exitCode=1});
