// Opt-in real local computation: writes regional test outputs, never changes sample inputs.
const { chromium } = require("../.local/web-qa/node_modules/playwright");
const assert = require("node:assert/strict");
const fs = require("node:fs/promises");
const path = require("node:path");
const crypto = require("node:crypto");

async function main() {
  const sampleId = process.env.CUTTAG_TEST_SAMPLE;
  assert.match(sampleId || "", /^[A-Za-z0-9_-][A-Za-z0-9_.-]{0,119}$/, "Set CUTTAG_TEST_SAMPLE to an authorized local sample");
  const origin = process.env.CUTTAG_TEST_URL || "http://127.0.0.1:8766";
  const sampleRoot = path.resolve(__dirname, "../data/private", sampleId);
  async function snapshot(folder, output = {}) {
    for (const item of await fs.readdir(folder, { withFileTypes: true })) {
      const file = path.join(folder, item.name);
      if (item.isDirectory()) await snapshot(file, output);
      else { const stat = await fs.stat(file); output[path.relative(sampleRoot, file)] = [stat.size, stat.mtimeMs]; }
    }
    return output;
  }
  const before = await snapshot(sampleRoot);
  const browser = await chromium.launch({ executablePath: process.env.CUTTAG_BROWSER_PATH, headless: true });
  try {
    const page = await browser.newPage({ viewport: { width: 1440, height: 1000 }, reducedMotion: "reduce" });
    const errors = [], external = [];
    page.on("pageerror", error => errors.push(error.message));
    page.on("console", message => { if (message.type() === "error") errors.push(message.text()); });
    page.on("request", request => { if (new URL(request.url()).origin !== origin) external.push(request.url()); });
    await page.goto(origin);
    await page.locator("#sample-select").selectOption(sampleId);
    await page.locator("#run-mode").selectOption("coverage");
    async function create(name) {
      await page.locator("#task-name").fill(name);
      const [response] = await Promise.all([
        page.waitForResponse(r => r.url() === origin + "/api/runs" && r.request().method() === "POST"),
        page.getByRole("button", { name: "开始真实计算", exact: true }).click(),
      ]);
      assert.equal(response.status(), 201);
      const record = await response.json();
      await page.waitForFunction(id => document.querySelector("#run-content").dataset.runId === id, record.id);
      return record;
    }
    const stopped = await create("工具联调 · 停止区域试跑");
    await page.getByRole("button", { name: "停止计算", exact: true }).click();
    await page.waitForFunction(() => document.querySelector("#current-run-label").textContent === "已停止", null, { timeout: 30000 });
    const stoppedRecord = await (await page.request.get(`${origin}/api/runs/${stopped.id}`)).json();
    assert.equal(stoppedRecord.status, "cancelled");
    assert.equal((await page.request.get(`${origin}/api/runs/${stopped.id}/artifacts/coverage.bw`)).status(), 404);

    const created = await create("工具联调 · 真实区域覆盖轨道");
    await page.waitForFunction(() => !["等待计算", "计算中"].includes(document.querySelector("#current-run-label").textContent), null, { timeout: 180000 });
    const record = await (await page.request.get(`${origin}/api/runs/${created.id}`)).json();
    assert.equal(record.status, "completed", JSON.stringify(record));
    assert.equal(record.actual_tool_executions, 3);
    assert.equal(record.bam_stats.quickcheck, "passed");
    assert.ok(Number.isInteger(record.bam_stats.mapped_segments));
    assert.equal(record.actual_model_requests, 0);
    assert.equal(record.execution.exit_code, 0);
    assert.equal(record.execution.tool_exit.exit_code, 0);
    assert.equal(record.execution.inputs_unchanged, true);
    assert.equal(record.execution.limits.network, "none");
    assert.equal(record.result.validation, "passed");
    assert.equal(record.result.preview.points.length, 200);
    assert.equal(record.result.start, 1000000);
    assert.equal(record.result.end, 2000000);
    await page.locator(".coverage-chart svg").waitFor();
    const [download] = await Promise.all([page.waitForEvent("download"), page.getByRole("link", { name: "下载 BW", exact: true }).click()]);
    const data = await fs.readFile(await download.path());
    assert.equal(data.readUInt32LE(0), 0x888ffc26);
    assert.equal(data.length, record.result.bytes);
    assert.equal(crypto.createHash("sha256").update(data).digest("hex"), record.result.sha256);
    await page.reload();
    await page.waitForFunction(id => document.querySelector("#run-content").dataset.runId === id && !!document.querySelector(".coverage-chart svg"), record.id);
    await page.waitForFunction(() => !document.querySelector("#toast").classList.contains("show"));
    const folder = path.resolve(__dirname, "../outputs/preview");
    await fs.mkdir(folder, { recursive: true });
    for (const [width, name] of [[1440, "desktop"], [390, "mobile"]]) {
      // A fresh mobile surface avoids a stale desktop compositor frame in full-page captures.
      const preview = width === 390 ? await browser.newPage({ viewport: { width: 390, height: 844 }, reducedMotion: "reduce" }) : page;
      if (preview !== page) { await preview.goto(origin); await preview.locator(".coverage-chart svg").waitFor(); }
      await preview.evaluate(async () => { await document.fonts.ready; await new Promise(resolve => requestAnimationFrame(() => requestAnimationFrame(resolve))); window.scrollTo(0, 0); });
      await preview.screenshot({ path: path.join(folder, `coverage-${name}.png`), fullPage: true });
      assert.equal(await preview.evaluate(() => document.documentElement.scrollWidth > innerWidth), false);
      if (preview !== page) await preview.close();
    }
    assert.deepEqual(await snapshot(sampleRoot), before);
    assert.deepEqual(errors, []);
    assert.deepEqual(external, []);
    const summary = { status: "passed", completed: record.id, cancelled: stopped.id, toolVersion: record.result.tool_version,
      region: record.parameters.region, bwBytes: data.length, inputMetadataUnchanged: true, toolExitCode: 0,
      actualToolExecutions: 3, actualModelRequests: 0, browserErrors: 0, externalRequests: 0 };
    await fs.writeFile(path.join(folder, "coverage-validation.json"), JSON.stringify(summary, null, 2) + "\n");
    console.log(JSON.stringify(summary));
  } finally { await browser.close(); }
}
main().catch(error => { console.error(error); process.exitCode = 1; });
