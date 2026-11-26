// Check the data-free startup state on an empty local checkout.
const { chromium } = require("../.local/web-qa/node_modules/playwright");
const assert = require("node:assert/strict");

async function main() {
  const origin = new URL(process.env.CUTTAG_TEST_URL || "http://127.0.0.1:8766").origin;
  const browser = await chromium.launch({ executablePath: process.env.CUTTAG_BROWSER_PATH });
  try {
    const page = await browser.newPage({ reducedMotion: "reduce" });
    const errors = [], external = [];
    page.on("pageerror", error => errors.push(error.message));
    page.on("request", request => {
      if (new URL(request.url()).origin !== origin) external.push(request.url());
    });
    await page.goto(origin);
    await page.getByText("请先准备本地样本和质控报告。", { exact: true }).waitFor();
    assert.ok(await page.locator("#start-demo").isDisabled());
    const views = [
      ["workspace", "#selected-sample-summary", /请先准备/],
      ["samples", "#samples-content", /暂无本地样本/],
      ["history", "#history-content", /暂无运行记录/],
      ["settings", "#settings-content", /模型/],
    ];
    for (const width of [1440, 1024, 768, 390, 320]) {
      await page.setViewportSize({ width, height: 900 });
      for (const [view, selector, expected] of views) {
        await page.locator(`.nav-item[data-view="${view}"]`).click();
        await page.locator(`#view-${view}`).waitFor({ state: "visible" });
        assert.match(await page.locator(selector).innerText(), expected);
        assert.equal(await page.evaluate(() => document.documentElement.scrollWidth > innerWidth),
                     false, `${view} overflows at ${width}px`);
      }
    }
    assert.deepEqual(errors, []);
    assert.deepEqual(external, []);
    console.log("Browser smoke check passed");
  } finally {
    await browser.close();
  }
}
main().catch(error => { console.error(error); process.exitCode = 1; });
