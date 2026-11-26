// Opt-in: two real model requests with explicitly authorized minimal summaries, two local tool rounds.
const { chromium } = require("../.local/web-qa/node_modules/playwright");
const assert = require("node:assert/strict");
const fs = require("node:fs/promises");
const path = require("node:path");
const crypto = require("node:crypto");

async function main() {
  const sampleId = process.env.CUTTAG_TEST_SAMPLE;
  assert.match(sampleId || "", /^[A-Za-z0-9_-][A-Za-z0-9_.-]{0,119}$/, "Set CUTTAG_TEST_SAMPLE to an authorized local sample");
  const origin = process.env.CUTTAG_TEST_URL || "http://127.0.0.1:8766";
  const root = path.resolve(__dirname, "..");
  const privateRoot = path.join(root, "data/private", sampleId);
  async function snapshot(folder, result = {}) {
    for (const entry of await fs.readdir(folder, {withFileTypes:true})) {
      const file = path.join(folder,entry.name);
      if(entry.isDirectory()) await snapshot(file,result);
      else {const s=await fs.stat(file);result[path.relative(privateRoot,file)]=[s.size,s.mtimeMs];}
    }
    return result;
  }
  const before = await snapshot(privateRoot);
  const browser = await chromium.launch({executablePath:process.env.CUTTAG_BROWSER_PATH,headless:true});
  try {
    const page = await browser.newPage({viewport:{width:1440,height:1000},reducedMotion:"reduce"});
    const errors=[],external=[];
    page.on("pageerror",e=>errors.push(e.message));
    page.on("console",m=>{if(m.type()==="error")errors.push(m.text());});
    page.on("request",r=>{if(new URL(r.url()).origin!==origin)external.push(r.url());});
    let id = process.env.CUTTAG_AGENT_REUSE_RUN;
    if (!id) {
      const creation = await page.request.post(origin + "/api/runs", {data: {
        sample_id:sampleId, mode:"agent", name:"Agent 验收 · 真实基线与一次重跑"
      }});
      assert.equal(creation.status(),201);
      id = (await creation.json()).id;
    }
    assert.match(id, /^run-[0-9a-f]{12}$/);
    const detailUrl = origin + "/#history/" + id;
    await page.goto(detailUrl);
    await page.waitForFunction(()=>["闭环完成","任务失败","已超时"].includes(document.querySelector("#history-run-label").textContent),null,{timeout:180000});
    const detail = page.locator("#history-run-content");
    const record=await(await page.request.get(`${origin}/api/runs/${id}`)).json();
    assert.equal(record.status,"completed_agent",JSON.stringify(record));
    assert.deepEqual(record.rounds.map(r=>r.result.bin_size),[25,50]);
    assert.equal(record.reruns_completed,1);
    assert.equal(record.actual_model_requests,2);
    assert.equal(record.actual_tool_executions,6);
    assert.equal(record.decisions[0].action,"regenerate_coverage");
    assert.equal(record.decisions[1].action,"finish");
    assert.equal(await detail.locator(".agent-comparison tbody tr").count(),2);
    assert.equal(await detail.locator(".agent-decisions .agent-decision").count(),2);
    assert.equal(await detail.locator(".agent-rounds svg").count(),2);
    for(const request of record.model_requests){
      const serialized=JSON.stringify(request.summary);
      assert.equal(serialized.includes(record.sample_id),false);
      assert.equal(serialized.includes("/Users/"),false);
      assert.equal(serialized.includes("/public/"),false);
      assert.equal(serialized.includes("raw/"),false);
      assert.equal(request.http_attempted,true);
    }
    for(const r of record.rounds){
      assert.equal(r.execution.exit_code,0);
      assert.equal(r.execution.inputs_unchanged,true);
      assert.equal(r.bam_stats.quickcheck,"passed");
      assert.ok(Number.isInteger(r.bam_stats.mapped_segments));
      const response=await page.request.get(`${origin}/api/runs/${record.id}/rounds/${r.index}/artifacts/coverage.bw`);
      const content=await response.body();
      assert.equal(response.status(),200);
      assert.equal(content.readUInt32LE(0),0x888ffc26);
      assert.equal(crypto.createHash("sha256").update(content).digest("hex"),r.result.sha256);
      assert.equal((await page.request.get(`${origin}/api/runs/${record.id}/rounds/${r.index}/artifacts/source.bam`)).status(),404);
    }
    const [download]=await Promise.all([page.waitForEvent("download"),detail.getByRole("link",{name:"下载运行记录",exact:true}).click()]);
    assert.equal(JSON.parse(await fs.readFile(await download.path(),"utf8")).id,record.id);
    await page.reload();
    await detail.locator(".agent-rounds svg").last().waitFor();
    const folder=path.join(root,"outputs/preview");
    for(const width of [1440,1280,1024,768,390,320]){
      const preview=await browser.newPage({viewport:{width,height:width<700?844:1000},reducedMotion:"reduce"});
      await preview.goto(detailUrl);
      await preview.locator("#history-run-content .agent-rounds svg").last().waitFor();
      await preview.evaluate(async()=>{await document.fonts.ready;await new Promise(r=>requestAnimationFrame(()=>requestAnimationFrame(r)));window.scrollTo(0,0);});
      assert.equal(await preview.evaluate(()=>document.documentElement.scrollWidth>innerWidth),false,`Agent overflow at ${width}`);
      if(width===1440||width===390)await preview.screenshot({path:path.join(folder,`agent-${width===1440?"desktop":"mobile"}.png`),fullPage:true});
      await preview.close();
    }
    await page.locator('[data-view="settings"]').click();
    assert.match(await page.locator("#settings-content").innerText(),/最近一次合成测试通过/);
    assert.match(await page.locator("#settings-content").innerText(),/已接入基础检查与统计/);
    assert.deepEqual(await snapshot(privateRoot),before);
    assert.deepEqual(errors,[]);
    assert.deepEqual(external,[]);
    const summary={status:"passed",run:record.id,reusedComputationForUIRecheck:!!process.env.CUTTAG_AGENT_REUSE_RUN,rounds:record.rounds.map(r=>({binSize:r.result.bin_size,bwBytes:r.result.bytes,sha256:r.result.sha256,exitCode:r.execution.exit_code})),
      modelRequests:2,samtoolsCommands:4,bamCoverageExecutions:2,reruns:1,inputMetadataUnchanged:true,modelExport:"explicitly authorized minimal summaries only",browserErrors:0,externalBrowserRequests:0,responsiveWidths:[1440,1280,1024,768,390,320],endReason:record.end_reason};
    await fs.writeFile(path.join(folder,"agent-validation.json"),JSON.stringify(summary,null,2)+"\n");
    console.log(JSON.stringify(summary));
  }finally{await browser.close();}
}
main().catch(e=>{console.error(e);process.exitCode=1});
