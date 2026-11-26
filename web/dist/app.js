"use strict";

const paths = {
  workspace: '<rect x="3" y="4" width="18" height="16" rx="1"/><path d="M3 9h18M9 9v11"/>',
  samples: '<path d="M4 5h16v15H4zM8 2v6M16 2v6M8 12h8M8 16h5"/>',
  history: '<circle cx="12" cy="12" r="9"/><path d="M12 7v5l3 2"/>',
  settings: '<path d="M3 7h18M3 17h18"/><circle cx="8" cy="7" r="2.5" fill="var(--soft)"/><circle cx="16" cy="17" r="2.5" fill="var(--soft)"/>',
};
const icon = (name) => `<svg viewBox="0 0 24 24" aria-hidden="true">${paths[name] || ""}</svg>`;
const escape = (value) => String(value ?? "").replace(/[&<>"']/g, (char) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" }[char]));
const $ = (selector) => document.querySelector(selector);
const state = { samples: [], runs: [], archivedRunIds: new Set(), setup: {}, selected: null, currentRun: null,
  historyRunId: null, historyRun: null, historyError: "", poll: null, busy: false,
  sampleTab: "results", libraryBaselines: new Map(), stopping: new Set(), deciding: new Set(), archiving: new Set(),
  showArchived: false, runError: "", pollError: "" };
const views = { workspace: "分析工作台", samples: "样本库", history: "运行记录", settings: "接入设置" };
const taskNames = {review:"单样本分析 · Agent", pipeline:"FASTQ 分析与模型复核",
  normalization:"预设 · CPM 与 Spike-in 标准化", igg:"预设 · IgG 对照检查", paired_agent:"配对样本分析 · Agent"};
const sampleTabs = ["results", "report", "files"];
const count = (value) => value == null ? "—" : Number(value).toLocaleString("zh-CN");
const percent = (value) => value == null ? "—" : (Number(value) * 100).toFixed(2) + "%";
const bytes = (value) => value >= 1e9 ? (value / 1e9).toFixed(2) + " GB" : value >= 1e6 ? (value / 1e6).toFixed(1) + " MB" : value >= 1000 ? Math.round(value / 1000) + " KB" : value + " B";
const stamp = (value, onlyTime = false) => new Intl.DateTimeFormat("zh-CN", { timeZone: "Asia/Shanghai", ...(onlyTime ? {} : { year: "numeric", month: "2-digit", day: "2-digit" }), hour: "2-digit", minute: "2-digit", ...(onlyTime ? { second: "2-digit" } : {}), hour12: false }).format(new Date(value));
const currentSample = () => state.samples.find((item) => item.id === state.selected) || state.samples[0];
const activeRun = (run) => ["running", "queued"].includes(run?.status);
const listedRun = (run) => activeRun(run) || !state.archivedRunIds.has(run.id);
const runLabel = (run) => ({ queued: run.mode === "demo" ? "等待演示" : "等待执行", running: ["review", "agent"].includes(run.mode) ? "Agent 运行中" : run.baseline_ready ? "质控复核中" : run.mode === "demo" ? "演示进行中" : "计算中", awaiting_rerun: "待你确认", completed_agent: "闭环完成", completed_review: run.paired_rerun_accepted === false ? "已保留原结果" : "复核完成", evidence_ready: "证据已整理", completed_pipeline: "基础分析完成", completed_spikein: "轨道已生成", completed_comparison: "两份轨道已生成", blocked_no_matched_igg: "无匹配 IgG，已停止找峰", completed: "计算完成", completed_demo: "演示完成", failed: "任务失败", timed_out: "已超时", cancelled: "已停止", interrupted: "已中断" }[run.status] || "未知状态");
const runClass = (run) => activeRun(run) ? "running" : ["failed", "timed_out"].includes(run.status) ? "failed" : run.status === "cancelled" ? "stopped" : run.status === "interrupted" ? "interrupted" : "";
const scrollToRun = () => $("#current-run-panel").scrollIntoView({ behavior: "instant", block: "nearest" });
const modelText = (value) => String(value || "").replace(/requested_bin_size/g, "目标窗口").replace(/bin_size/g, "窗口").replace(/reruns_left/g, "剩余重跑次数");
const latestResult = (mode, status) => state.runs.find(r => r.sample_id === currentSample()?.id && r.mode === mode && r.status === status);
const cpmBaseline = () => completedFor(currentSample()?.id, "pipeline", "completed_pipeline")[0];
const awaitingPairAdvice = sampleId => state.runs.find(run => run.sample_id === sampleId && run.paired_adaptive && run.status === "awaiting_rerun" && listedRun(run));
const decisionLabel = d => ({keep:"保留现有结果", request_evidence:"需要补充证据", finish:"结束任务",
  tune_igg_qvalue:`重跑带 IgG 的找峰，q = ${d.qvalue}`, prefer_rerun:"新版对本次探索目标更合适",
  compare_peak_sets:"比较窄峰和宽峰", compare_igg:"比较加入 IgG 前后的候选峰", compare_resolution:`用 ${d.bin_size} bp 聚合现有信号做比较`,
  regenerate_coverage:`将窗口调整为 ${d.bin_size} bp`})[d.action] || "未识别的建议";

function toast(message) {
  const node = $("#toast");
  node.textContent = message;
  node.classList.add("show");
  clearTimeout(toast.timer);
  toast.timer = setTimeout(() => node.classList.remove("show"), 3200);
}

async function api(path, options = {}) {
  const response = await fetch(path, { ...options, headers: { "Content-Type": "application/json", ...options.headers } });
  const data = await response.json();
  if (!response.ok) throw new Error(data.error || "本地服务暂时没有响应");
  return data;
}

function showView(route) {
  const runId = /^history\/(run-[0-9a-f]{12})$/.exec(route)?.[1] || null;
  const previousRunId = state.historyRunId;
  let view = runId ? "history" : route;
  if (!views[view]) view = "workspace";
  document.querySelectorAll(".view").forEach((node) => { node.hidden = node.id !== "view-" + view; });
  document.querySelectorAll(".nav-item").forEach((node) => {
    node.classList.toggle("active", node.dataset.view === view);
    if (node.dataset.view === view) node.setAttribute("aria-current", "page");
    else node.removeAttribute("aria-current");
  });
  document.title = views[view] + " · CUT&Tag";
  $("#history-content").hidden = !!runId;
  $("#history-create").hidden = !!runId;
  $("#history-detail").hidden = !runId;
  if (runId !== previousRunId) {
    state.historyRunId = runId;
    state.historyRun = state.runs.find(run => run.id === runId) || null;
    state.historyError = "";
    if (runId) loadHistoryRun(runId);
  }
  if (view === "history" && !runId) {
    renderHistory();
    if (previousRunId) $("#history-content").querySelector(`[data-run="${previousRunId}"]`)?.focus();
  }
  if (runId || view === "workspace") renderRun();
  if (runId && runId !== previousRunId) $("#history-back").focus();
  if (view === "samples") renderSamples();
  if (view === "settings") renderSettings();
}

async function loadHistoryRun(id) {
  try {
    const detail = await api(`/api/runs/${encodeURIComponent(id)}`);
    if (state.historyRunId === id) state.historyRun = detail;
  } catch (error) {
    if (state.historyRunId === id) state.historyError = error.message;
  }
  if (state.historyRunId === id) renderRun();
}

function go(view) {
  if (location.hash === "#" + view) showView(view);
  else location.hash = view;
}

function renderStartButton() {
  const running = state.runs.some(activeRun), mode = $("#run-mode").value;
  if (currentSample()?.pair) {
    $("#start-demo").disabled = running || state.busy || !modelReady();
    $("#start-demo").textContent = running ? "任务进行中…" : "启动 Agent";
    return;
  }
  const review = mode === "review", preset = mode === "normalization" || mode === "igg";
  const ready = preset ? modelReady() && !!presetReview(mode) && (mode === "igg" || state.setup.pipeline_runtime?.ready) :
    state.setup.pipeline_runtime?.ready && (review ? !!cpmBaseline() : currentSample()?.report_available);
  $("#start-demo").disabled = running || state.busy || !ready;
  $("#start-demo").textContent = running ? "任务进行中…" : preset ? "运行预设" : review ?
    (modelReady() ? "启动 Agent" : "查看本地质控证据") : (modelReady() ? "开始分析并复核" : "开始分析");
}

function modelAuthorized() {
  const sample = currentSample();
  return sample?.pair ? sample.pair.external_model_allowed : state.setup.qc_review_sample_ids?.includes(sample?.id);
}

function modelReady() { return state.setup.model_connected && state.setup.minimal_summary_allowed && modelAuthorized(); }

function presetReview(mode) {
  return state.runs.find(run => {
    if (run.mode !== "review" || run.status !== "completed_review" || run.sample_id !== currentSample()?.id || !listedRun(run)) return false;
    const qc = run.track_review?.sent_summaries?.[0]?.qc;
    return mode === "normalization" ? qc?.spikein?.spikein_fragments > 0 :
      qc?.experiment_context?.matched_igg === "unavailable";
  });
}

function renderMode() {
  const select = $("#run-mode"), baseline = cpmBaseline();
  const pair = currentSample()?.pair, name = $("#task-name");
  const options = $("#analysis-options");
  options.hidden = !!pair;
  if (!pair && !baseline) options.open = true;
  name.disabled = !!pair || !options.open;
  const pairedOption = select.querySelector('[value="paired_agent"]');
  pairedOption.hidden = pairedOption.disabled = !pair;
  select.querySelector('[value="pipeline"]').disabled = !!pair;
  select.querySelector('[value="review"]').disabled = !pair && !baseline;
  for (const mode of ["normalization", "igg"]) select.querySelector(`[value="${mode}"]`).disabled = !!pair || !presetReview(mode);
  if (pair) {
    select.value = "paired_agent";
    const pending = !modelAuthorized() ? "这组摘要的外部模型使用待授权。" : !modelReady() ? "请先完成模型连接测试。" : "";
    const scope = "Agent 自动整理证据，再决定复核、比较还是建议调参。你接受建议后才会重跑。";
    $("#input-note").innerHTML = `${scope}${pending} <a href="${pairUrl(currentSample())}/review-summary" download>下载拟发送摘要</a>`;
  } else {
    if (!options.open || select.selectedOptions[0]?.disabled) select.value = baseline ? "review" : "pipeline";
    const mode = select.value;
    const action = modelReady() ? "模型会根据证据自行选择下一步。" : "先在本地整理证据，完整摘要获准后可调用模型。";
    const notes = {normalization: "复用配套 BAM 生成两份全基因组轨道。上次用时约 6 分钟 20 秒。",
      igg: "核对已确认无匹配 IgG 的情况。条件缺失时保留原 BED。"};
    $("#input-note").textContent = notes[mode] ? notes[mode] + (modelReady() ? "" : " 请先确认模型接入和摘要许可。") :
      mode === "review" && baseline ? `复用已有 CPM 轨道、候选峰和质控。 ${action}` :
        `上次比对用了约 1 小时。 ${action}`;
  }
  if (Object.values(taskNames).includes(name.value)) name.value = taskNames[select.value];
  renderStartButton();
}

function renderExistingResults() {
  const current = currentSample();
  const pending = current?.pair && awaitingPairAdvice(current.id);
  $("#workspace-advice").hidden = !pending;
  $("#workspace-advice").innerHTML = pending ? pairAdviceLink(pending) : "";
  if (current?.pair) {
    $("#workspace-results").hidden = false;
    $("#workspace-results-list").innerHTML = '<li><span>全基因组轨道与两套候选峰</span><button class="table-action" data-pair-tab="results">查看结果</button></li><li><span>加样记录与 IGV 复核窗口</span><button class="table-action" data-pair-tab="report">查看本地证据</button></li>';
    return;
  }
  const items = [[cpmBaseline(), "CPM · 全基因组分析"], [latestResult("spikein", "completed_spikein"), "Spike-in · 全基因组轨道"]].filter(([run]) => run);
  $("#workspace-results").hidden = !items.length;
  $("#workspace-results-list").innerHTML = items.map(([run, title]) => `<li><span>${title}</span><button class="table-action" data-run="${escape(run.id)}" aria-label="查看${title}">查看结果</button></li>`).join("");
}

function logPreview(url, key, open, title = "工具日志") {
  return `<details class="log-preview" data-detail="${escape(key)}" data-log-url="${escape(url)}" ${open.includes(key) ? "open" : ""}>
    <summary>${escape(title)}</summary><pre tabindex="0" aria-label="${escape(title)}">正在读取日志…</pre>
    <a href="${escape(url)}" download>下载完整日志</a></details>`;
}

const logRequests = new Map();
document.addEventListener("toggle", async (event) => {
  const detail = event.target, url = detail.dataset?.logUrl;
  if (!url || !detail.open) return;
  const pre = detail.querySelector("pre");
  pre.textContent = "正在读取日志…";
  if (!logRequests.has(url)) logRequests.set(url, fetch(url, {headers: {Range: "bytes=-65536"}}).then(async response => {
    if (response.status === 416) return "日志文件暂时为空。";
    if (!response.ok) throw new Error("日志暂时无法读取，可以稍后重试。");
    const text = (await response.text()).trimEnd();
    return text || "日志文件暂时为空。";
  }).finally(() => logRequests.delete(url)));
  try { pre.textContent = await logRequests.get(url); }
  catch (error) { pre.textContent = error.message; }
}, true);

function coverageChart(result, title = "本次区域覆盖轨道") {
  const values = result.preview.points;
  const max = Math.max(1, ...values.filter((value) => value != null));
  const left = 68, right = 568, top = 20, bottom = 185;
  let gap = true;
  const line = values.map((value, index) => {
    if (value == null) { gap = true; return ""; }
    const segment = `${gap ? "M" : "L"}${(left + (index + .5) / values.length * (right - left)).toFixed(1)},${(bottom - value / max * (bottom - top)).toFixed(1)}`;
    gap = false;
    return segment;
  }).join(" ");
  return `<section class="chart-section coverage-chart"><h3>${escape(title)}</h3><p class="chart-unit">覆盖信号（未标准化）</p><div class="fragment-chart"><svg viewBox="0 0 600 245" role="img" aria-label="本次真实 bamCoverage 生成的 chr1 区域覆盖轨道">${[0, .5, 1].map((v) => { const y = bottom - v * (bottom - top); return `<line class="chart-grid" x1="${left}" x2="${right}" y1="${y}" y2="${y}"/><text class="chart-axis" x="${left - 8}" y="${y + 4}" text-anchor="end">${(max * v).toFixed(1)}</text>`; }).join("")}<path class="chart-line" d="${line}"/>${[0, .25, .5, .75, 1].map((v) => `<text class="chart-axis" x="${left + v * (right - left)}" y="208" text-anchor="${v === 1 ? "end" : "middle"}">${(1 + v).toFixed(2)}</text>`).join("")}<text class="chart-axis" x="${right}" y="237" text-anchor="end">chr1 位置（Mb）</text></svg></div><p class="source-note">来自本次 BW，按 200 段均值预览。仅为区域试跑，非全基因组结果。</p></section>`;
}

const evidenceLabels = {sequencing_qc:"测序与清理", pair_filter:"双端过滤", bam_counts:"BAM 比对统计", fragment_qc:"片段与重复",
  bw_statistics:"BigWig 局部信号", global_track:"BigWig 全基因组统计", narrow_peaks:"窄峰与 FRiP", broad_peaks:"宽峰与 FRiP",
  tss_profile:"TSS 平均信号", spikein:"Spike-in", experiment_context:"实验背景与对照", comparison:"本次比较",
  baseline:"无 IgG 候选峰与落峰比例", settings:"分析参数", igv_screen:"被筛掉的峰与 IGV 探索统计"};
const contextLabels = {matched_igg:"匹配 IgG 的情况", target_expectation:"实验靶点与预期富集位点",
  biological_replicates:"生物学重复与分组", spikein_design:"Spike-in 实验设计", target_antibody:"确切修饰与抗体",
  control_design:"实验对照设计", igv_review:"人工 IGV 复核", spikein_quantitation:"Spike-in 定量可比性"};

function renderReviewEvidence(review, open) {
  const payload = review.evidence, qc = payload.qc;
  const sent = review.sent_summaries?.length;
  if (!qc) return `<section class="reused-results"><h3>配对质控证据</h3>
    ${(sent ? review.sent_summaries : [payload]).map((summary, index) => `<details class="log-preview" data-detail="pair-qc-${index}" ${open.includes(`pair-qc-${index}`) ? "open" : ""}>
    <summary>${sent ? index || summary.comparison?.rerun ? "查看比较后请求的摘要" : "查看首轮请求的摘要" : "查看尚未发送的摘要"}</summary>
    <pre>${escape(JSON.stringify(summary, null, 2))}</pre></details>`).join("")}</section>`;
  const included = ["bam_counts", "bw_statistics", ...Object.keys(qc)].map(key => evidenceLabels[key]);
  const absent = ["sequencing_qc", "pair_filter", "tss_profile", "spikein"].filter(key => !qc[key]).map(key => evidenceLabels[key]);
  const context = qc.experiment_context;
  const igg = context.matched_igg === "unavailable" ? "当前样本已确认无匹配 IgG。" : "匹配 IgG 的情况待确认。";
  const unknown = ["target_expectation", "biological_replicates", "spikein_design"]
    .filter(key => context[key] === "unconfirmed").map(key => contextLabels[key]);
  const bam = payload.core.bam_counts, narrow = qc.narrow_peaks, broad = qc.broad_peaks;
  const highlights = [
    qc.sequencing_qc && bam.pre_filter_primary_segments > 0 ?
      `清理后 Q30 碱基比例 ${percent(qc.sequencing_qc.clean_q30)}，主比对率 ${percent(bam.pre_filter_primary_mapped_segments / bam.pre_filter_primary_segments)}。` : "",
    narrow && broad ? `窄峰 FRiP ${percent(narrow.frip)}，宽峰 FRiP ${percent(broad.frip)}。两套峰的覆盖范围不同。` : "",
    qc.spikein ? `有 ${count(qc.spikein.spikein_fragments)} 个 Spike-in 片段。${qc.spikein.same_alignment_processing_as_cpm === false ?
      "轨道使用的 BAM 处理和 CPM 不同。" : "两条轨道的 BAM 可比性尚未确认。"}` : ""
  ].filter(Boolean);
  return `<section class="reused-results"><h3>质控证据</h3>
    <p>${escape(included.join("、"))}</p>
    ${highlights.length ? `<ul class="decision-evidence">${highlights.map(row => `<li>${escape(row)}</li>`).join("")}</ul>` : ""}
    <p class="field-help">${igg}${unknown.length ? `待确认的信息包括${escape(unknown.join("、"))}。` : ""}</p>
    ${absent.length ? `<p class="field-help">缺少经核对的本地产物：${escape(absent.join("、"))}。</p>` : ""}
    <details class="log-preview" data-detail="qc-payload" ${open.includes("qc-payload") ? "open" : ""}><summary>${sent ? "查看首轮请求的证据摘要" : "查看本地尚未发送的证据摘要"}</summary>
      <pre>${escape(JSON.stringify(sent ? review.sent_summaries[0] : payload, null, 2))}</pre>
    </details>
    ${review.comparison ? `<details class="log-preview" data-detail="qc-comparison" ${open.includes("qc-comparison") ? "open" : ""}><summary>比较结果</summary>
      <pre>${escape(JSON.stringify(review.comparison, null, 2))}</pre></details>` : ""}
  </section>`;
}

function renderModelTrace(run, open = []) {
  const entry = (title, d, effect, summary, key, context) => {
    const size = summary?.parameters?.bin_size, stat = summary?.bw_statistics, bam = summary?.bam_counts;
    const hasUpstream = bam?.pre_filter_primary_segments > 0 && bam?.pre_filter_primary_mapped_segments != null;
    const regional = stat?.interval_count > 0 && stat?.nonzero_intervals != null;
    const legacyEvidence = [hasUpstream ? `<li>清理后的主比对率：${percent(bam.pre_filter_primary_mapped_segments / bam.pre_filter_primary_segments)}（过滤前 BAM）</li>` : "",
      bam?.mapped_segments != null ? `<li>过滤后的 BAM：${count(bam.mapped_segments)} 条已比对读段${bam.post_filter_duplicate_segments != null ? `，其中 ${count(bam.post_filter_duplicate_segments)} 条标记为重复` : ""}</li>` : "",
      size != null ? `<li>信号窗口：${escape(size)} bp${regional ? `；所查区域 ${percent(stat.nonzero_intervals / stat.interval_count)} 的窗口有信号` : ""}</li>` : ""].join("");
    const evidence = d.evidence_used?.length ? d.evidence_used.map(name => `<li>${escape(evidenceLabels[name] || name)}</li>`).join("") : legacyEvidence;
    const current = key === "current-track-review";
    const missing = d.missing_information || [];
    const iggCorrection = context?.matched_igg === "unavailable" && missing.some(item => /igg/i.test(item));
    return `<section class="agent-decision${current ? " current-decision" : ""}"><h4>${escape(title)}</h4>
      <p class="decision-result"><span>模型建议</span><strong>${escape(decisionLabel(d))}</strong></p>
      <dl class="decision-summary"><div><dt>实际执行</dt><dd>${escape(effect)}</dd></div>${iggCorrection ? "" : `<div><dt>理由</dt><dd>${escape(modelText(d.reason))}</dd></div>`}</dl>
      ${iggCorrection ? '<p class="field-help">当前样本已确认无匹配 IgG。此轮模型将它错列为待补证据。</p>' : ""}
      <details data-detail="${escape(key)}" ${open.includes(key) ? "open" : ""}><summary>${iggCorrection ? "引用的证据与模型原回复" : "引用的证据和缺失信息"}</summary>
        ${evidence ? `<ul class="decision-evidence">${evidence}</ul>` : ""}<p>${escape(modelText(d.assessment))}</p>
        ${iggCorrection ? `<p>${escape(d.reason)}</p>` : ""}
        ${missing.length ? `<p class="field-help">${iggCorrection ? "原回复所列信息" : run.paired_adaptive && ["exploratory_reviewed", "automated_screened"].includes(context?.igv_review) ? "生信验收时确认" : "模型认为还需确认"}</p><ul class="decision-evidence">${missing.map(item => `<li>${escape(context && Object.hasOwn(contextLabels, item) ? contextLabels[item] : item)}</li>`).join("")}</ul>` : ""}</details></section>`;
  };
  const decisions = run.decisions || [], requests = run.model_requests || [], later = run.linked_reviews || [], review = run.track_review;
  let rows = decisions.map(d => {
    const previous = run.rounds?.find(r => r.index === d.round_index), next = run.rounds?.find(r => r.index === d.round_index + 1);
    const req = requests.find(r => r.round_index === d.round_index);
    const effect = d.action === "finish" ? "未继续重跑。" : next?.status === "completed" ?
      `已重跑：${previous?.parameters?.bin_size ?? "—"} → ${next.parameters.bin_size} bp，产物检查通过。` :
      next ? "已尝试重跑，未产出通过检查的结果。" : activeRun(run) ? "尚未完成重跑。" : "未执行重跑。";
    return entry(d.round_index ? "重跑后 · 模型复评" : "基线后 · 模型建议", d, effect, req?.summary, `assessment-${d.round_index}`);
  }).join("");
  rows += requests.filter(r => r.error && !decisions.some(d => d.round_index === r.round_index)).map(r =>
    `<p class="field-help">${r.round_index ? "重跑后" : "基线后"}未取得可执行的模型建议：${escape(r.error)}。没有据此调整。</p>`).join("");
  const trackEntry = (r, title, key) => {
    const d = r.answer?.decision;
    if (!d) return `<section class="agent-decision"><h4>${escape(title)}</h4><p>${escape(r.error || (r.status === "reviewing" ?
      "尚未收到模型建议。" : r.status === "evidence_ready" ? "完整质控摘要待确认使用范围，暂无本次模型决定。" : "未取得可执行的模型建议，未调整。"))}</p></section>`;
    const effect = run.paired_adaptive ? d.action === "tune_igg_qvalue" ? (r.rerun ? "已生成新的 BED 和落峰统计。" :
      r.error ? "重跑未生成经过校验的结果。" : r.rerun_started ? "正在执行找峰和落峰统计。" :
      run.paired_rerun_accepted === false ? "已选择保留原结果。" : run.paired_rerun_accepted ? "已接受建议，等待工具开始。" : "等待你确认是否重跑。") :
      "采用已有结果，本轮无需重跑。" : r.comparison ? comparisonEffect(r.comparison) : d.action.startsWith("compare_") ?
      "比较未完成。" : "本次无需重跑，保留现有产物。";
    return entry(title, d, effect, r.summary, key, r.evidence?.qc?.experiment_context || r.evidence?.experiment_context);
  };
  if (review) {
    rows += trackEntry(review, review.status === "failed" && !review.answer ? "复核失败原因" : "下一步怎么做", "current-track-review");
    if (review.reassessment?.decision) rows += entry(run.paired_adaptive ? "重跑后 · 模型复评" : "比较后 · 模型复评",
      review.reassessment.decision, "复评完成，新旧产物均已保存。", null, "qc-reassessment", review.evidence?.qc?.experiment_context || review.evidence?.experiment_context);
    else if (review.answer && review.error) rows += `<p class="inline-error">${escape(review.error)}</p>`;
  }
  if (!rows) rows = `<p class="field-help">${(run.mode === "agent" || run.mode === "review" || run.model_review_requested) && activeRun(run) ? "尚未收到模型建议。" : run.actual_model_requests ? "有模型请求，但没有可执行的建议。" : "本次运行未调用模型，没有模型调整。"}</p>`;
  rows += later.map((r, i) => trackEntry(r, `事后复核 · ${stamp(r.at)}`, `review-${i}`)).join("");
  return `<section class="model-trace agent-decisions" aria-label="模型决策记录"><h3>模型决策记录</h3>${rows}</section>
    ${review?.evidence ? renderReviewEvidence(review, open) : ""}`;
}

function comparisonEffect(comparison) {
  if (comparison.action === "compare_peak_sets") return "已比较同一 BAM 的候选窄峰与宽峰，保留两套结果。";
  if (comparison.action === "compare_igg") return `已读取带 IgG 的 ${count(comparison.peak_count)} 个候选峰。目标 FRiP ${percent(comparison.target.frip)}，IgG FRiP ${percent(comparison.igg.frip)}。`;
  return `已用 ${comparison.bin_bp ?? comparison.window_bp} bp 对比现有区域信号，保留原轨道。`;
}

function renderAgent(run, content) {
  const open = content.dataset.runId === run.id ? Array.from(content.querySelectorAll("details[open]")).map(n => n.dataset.detail) : [];
  const rounds = run.rounds || [], decisions = run.decisions || [];
  const link = (round, name) => `/api/runs/${encodeURIComponent(run.id)}/rounds/${round}/artifacts/${name}`;
  content.dataset.runId = run.id;
  content.innerHTML = `<div class="run-topline"><h3>${escape(run.name)}</h3><p>${escape(run.sample_id)} · ${stamp(run.created_at)}</p></div>
    ${activeRun(run) ? `<p class="run-disclaimer" role="status">${escape(run.message)}</p>` : ""}
    <p class="source-note">工程目标：生成 50 bp 区域轨道，基线为 25 bp。窗口调整不代表实验质量改善。</p>
    ${rounds.length ? `<div class="table-wrap agent-comparison"><table><thead><tr><th>轮次</th><th>窗口</th><th>执行 / 文件检查</th><th>结果</th></tr></thead><tbody>${rounds.map(r => `<tr><td>${r.index ? "重跑" : "基线"}</td><td>${r.parameters.bin_size} bp</td><td>${r.status === "completed" ? "退出正常 / 通过" : escape(runLabel(r))}</td><td>${r.status === "completed" ? `<a href="${link(r.index, "coverage.bw")}" download>BW</a> <a href="${link(r.index, "tool.log")}" download>日志</a>` : "—"}</td></tr>`).join("")}</tbody></table></div>` : ""}
    ${renderModelTrace(run, open)}
    ${!activeRun(run) && !(run.status === "completed_agent" && decisions.at(-1)?.reason === run.end_reason) ? `<p class="run-disclaimer ${run.status === "failed" ? "inline-error" : ""}" role="${run.status === "failed" ? "alert" : "status"}">结束原因：${escape(modelText(run.end_reason || run.message))}</p>` : ""}
    ${state.pollError ? `<p class="inline-error" role="alert">${escape(state.pollError)}</p>` : ""}${state.runError ? `<p class="inline-error" role="alert">${escape(state.runError)}</p>` : ""}
    <div class="agent-rounds">${rounds.filter(r => r.status === "completed").map(r => coverageChart(r.result, `${r.index ? "重跑" : "基线"} · ${r.parameters.bin_size} bp`)).join("")}</div>
    ${rounds.filter(r => r.artifacts?.some(f => f.name === "tool.log")).map(r => logPreview(link(r.index, "tool.log"), `round-log-${r.index}`, open, r.index ? "重跑工具日志" : "基线工具日志")).join("")}
    <div class="run-actions">${activeRun(run) ? `<button class="button danger" data-cancel="${escape(run.id)}" ${state.stopping.has(run.id) ? "disabled" : ""}>${state.stopping.has(run.id) ? "正在停止…" : "停止任务"}</button>` : ""}<a class="button secondary" href="/api/runs/${encodeURIComponent(run.id)}/report" download>下载运行记录</a></div>
    <details data-detail="steps" ${open.includes("steps") ? "open" : ""}><summary>运行步骤</summary><ol class="timeline">${run.steps.map((step, i) => `<li class="timeline-step"><span class="step-number">${i + 1}</span><div><strong>${escape(step.label)}</strong><span class="step-state">${escape({pending:"未开始",running:"进行中",done:"完成",failed:"失败",cancelled:"已停止",skipped:"未执行"}[step.status] || step.status)}</span><p>${escape(step.detail)}</p></div></li>`).join("")}</ol></details>
    <details data-detail="events" ${open.includes("events") ? "open" : ""}><summary>事件记录</summary><ul class="event-list">${run.events.map(event => `<li><time>${stamp(event.at, true)}</time><span>${escape(event.message)}</span></li>`).join("")}</ul></details>`;
}

function renderSampleSelection() {
  const sample = currentSample();
  $("#sample-select").innerHTML = state.samples.length ? state.samples.map((item) => `<option value="${escape(item.id)}" ${item.id === state.selected ? "selected" : ""}>${escape(item.label)}</option>`).join("") : '<option value="">尚未导入样本</option>';
  $("#selected-sample-summary").textContent = !sample ? "请先准备本地样本和质控报告。" : sample.pair ?
    `匹配 IgG 为 ${sample.pair.control_id}` : cpmBaseline() ?
      (latestResult("spikein", "completed_spikein")?.matched_igg === false ? "已有全基因组分析 · 无匹配 IgG" : "已有全基因组分析") :
      !sample.report_available ? "历史报告缺失。" : "输入：" + sample.input_kind;
  renderStartButton();
}

function chartMarkup(sample) {
  const counts = sample?.fragment_distribution?.counts || [];
  const max = Math.max(0, ...counts);
  if (!max) return '<div class="empty-state"><p>暂无片段长度统计。</p></div>';
  const left = 68, right = 568, top = 18, bottom = 205;
  const line = counts.map((value, index) => `${index ? "L" : "M"}${(left + (index + .5) / counts.length * (right - left)).toFixed(1)},${(bottom - value / max * (bottom - top)).toFixed(1)}`).join(" ");
  return `<p class="chart-unit">片段数</p><svg viewBox="0 0 600 258" role="img" aria-label="历史片段长度分布，横轴为片段长度（bp），纵轴为每 10 bp 分箱的片段数">${[0, .5, 1].map((v) => { const y = bottom - v * (bottom - top); return `<line class="chart-grid" x1="${left}" x2="${right}" y1="${y}" y2="${y}"/><text class="chart-axis" x="${left - 8}" y="${y + 4}" text-anchor="end">${v ? (max * v / 10000).toFixed(1) + "万" : "0"}</text>`; }).join("")}<path class="chart-line" d="${line}"/>${[0, 200, 400, 600, 800, 1000].map((v) => `<text class="chart-axis" x="${left + v / 1000 * (right - left)}" y="227" text-anchor="${v === 1000 ? "end" : "middle"}">${v}</text>`).join("")}<text class="chart-axis" x="${right}" y="250" text-anchor="end">片段长度（bp）</text></svg>`;
}

function setSampleTab(tab, focus = false) {
  state.sampleTab = sampleTabs.includes(tab) ? tab : "results";
  document.querySelectorAll("[data-sample-tab]").forEach((node) => {
    const selected = node.dataset.sampleTab === state.sampleTab;
    node.setAttribute("aria-selected", String(selected));
    node.tabIndex = selected ? 0 : -1;
    if (selected && focus) node.focus();
  });
  sampleTabs.forEach((name) => { const panel = $("#sample-panel-" + name); if (panel) panel.hidden = name !== state.sampleTab; });
}

function completedFor(id, mode, status) {
  return state.runs.filter(r => r.sample_id === id && r.mode === mode && r.status === status && r.result?.validation === "passed");
}

function sampleResults(sample) {
  const baselines = completedFor(sample.id, "pipeline", "completed_pipeline");
  const spikein = completedFor(sample.id, "spikein", "completed_spikein");
  const chosen = baselines.find(r => r.id === state.libraryBaselines.get(sample.id)) || baselines[0];
  const file = (run, name) => `/api/runs/${encodeURIComponent(run.id)}/artifacts/${encodeURIComponent(name)}`;
  const exists = (run, name) => run.artifacts?.some(a => a.name === name);
  const download = (run, name, label) => exists(run, name) ? `<a href="${file(run, name)}" download>${escape(label)}</a>` : "";
  const picker = baselines.length > 1 ? `<div class="field library-baseline-picker"><label for="library-baseline">CPM 分析版本</label>
    <select id="library-baseline">${baselines.map(r => `<option value="${escape(r.id)}" ${r.id === chosen.id ? "selected" : ""}>${stamp(r.updated_at)} · ${escape(r.name)}</option>`).join("")}</select></div>` : "";
  let cpm = "";
  if (chosen) {
    const qc = chosen.result.qc, meta = chosen.result.coverage;
    const metrics = [["主比对率（清理后）", chosen.report_metrics?.primary_mapping_rate, percent],
      ["用于分析的配对片段", qc?.fragments, count], ["标记为重复的片段比例", qc?.duplicate_fraction, percent]].filter(([,v]) => v != null);
    const reviews = state.runs.filter(r => listedRun(r) && r.sample_id === sample.id && (r.id === chosen.id || r.baseline_id === chosen.id) &&
      (r.track_review?.status === "evidence_ready" || (r.track_review?.status === "completed_review" && r.track_review.answer?.decision)))
      .sort((a,b) => Date.parse(b.updated_at) - Date.parse(a.updated_at));
    const latest = reviews.find(r => r.track_review.status === "completed_review");
    const prepared = reviews[0]?.track_review?.status === "evidence_ready" ? reviews[0] : null;
    const decision = latest?.track_review?.reassessment?.decision || latest?.track_review?.answer?.decision;
    const compared = latest?.track_review?.comparison;
    const iggKnownAbsent = (latest?.track_review?.evidence || prepared?.track_review?.evidence)?.qc?.experiment_context?.matched_igg === "unavailable";
    const effect = compared ? comparisonEffect(compared) : latest?.track_review?.evidence ?
      "使用完整质控摘要" : "使用早期最小运行摘要";
    cpm = `<section class="library-section" aria-labelledby="library-cpm-heading">
      <div class="library-section-heading"><div><h3 id="library-cpm-heading">CPM · 全基因组分析</h3><p class="field-help">hg38 · ${escape(meta.bin_size)} bp · ${stamp(chosen.updated_at)}</p></div>
        <button class="table-action" data-run="${escape(chosen.id)}">查看完整分析</button></div>
      ${metrics.length ? `<dl class="library-metrics">${metrics.map(([label,value,format]) => `<div><dt>${label}</dt><dd>${format(value)}</dd></div>`).join("")}</dl>` : ""}
      <div class="library-files">${download(chosen, "coverage.cpm.bw", "全基因组 BigWig（IGV）")}${download(chosen, "report.md", "分析说明")}${download(chosen, "fragment-qc.json", "质控数据")}</div>
      <div class="library-detail-grid">
        <section aria-labelledby="library-peaks-heading"><h3 id="library-peaks-heading">候选峰</h3>
          <div class="table-wrap"><table class="library-peaks"><thead><tr><th scope="col">类型</th><th scope="col" class="numeric">数量</th><th scope="col" class="numeric">FRiP</th><th scope="col">文件</th></tr></thead>
          <tbody>${[["narrow","窄峰"],["broad","宽峰"]].map(([key,label]) => `<tr><th scope="row">${label}</th><td class="numeric">${count(qc?.peaks?.[key]?.count)}</td><td class="numeric">${percent(qc?.frip?.[key]?.fraction)}</td><td>${download(chosen, `${key}.peaks.bed`, `下载${label} BED`)}</td></tr>`).join("")}</tbody></table></div>
          <p class="source-note">${iggKnownAbsent ? "当前样本已确认无匹配 IgG。" : "本流程未使用匹配的 IgG。"}FRiP 是与对应候选峰重叠的片段比例，不可用来直接比较窄峰与宽峰优劣。</p>
        </section>
        ${exists(chosen,"tss.profile.png") ? `<section aria-labelledby="library-tss-heading"><h3 id="library-tss-heading">转录起始位点（TSS）周围信号</h3><img class="library-plot" src="${file(chosen,"tss.profile.png")}" alt="CPM 平均信号曲线，横轴为转录起始位置前后 3 kb" loading="lazy"><p class="source-note">用来观察平均信号，不作为所有靶点的通用合格线。</p></section>` : ""}
      </div>
      ${prepared ? `<div class="library-review"><div><h3>本地质控证据</h3><p>已整理，待确认外部使用范围。</p></div>
        <button class="table-action" data-run="${escape(prepared.id)}">查看证据</button></div>` : ""}
      ${decision ? `<div class="library-review"><div><h3>这套结果的最近一次模型复核</h3>
        <p><strong>${escape(decisionLabel(decision))}</strong> · ${escape(effect)}</p></div>
        <button class="table-action" data-run="${escape(latest.id)}">查看模型依据</button></div>` : ""}
    </section>`;
  }
  const spike = spikein.length ? `<section class="library-section library-spikein" aria-labelledby="library-spikein-heading"><h3 id="library-spikein-heading">Spike-in · 独立轨道</h3>
    <ul>${spikein.map(r => `<li><span>${stamp(r.updated_at)} · ${escape(r.result.bin_size)} bp · ${count(r.spikein_fragments)} 个 spike-in 片段</span><div class="library-files">
      ${download(r, "coverage.spikein.bw", "Spike-in BigWig")}<button class="table-action" data-run="${escape(r.id)}">查看来源</button></div></li>`).join("")}</ul>
    <p class="source-note">使用配套历史输入生成。若要和 CPM 比标准化，先确认两边的 BAM 来源和窗口一致。</p></section>` : "";
  const pending = state.runs.find(r => r.sample_id === sample.id && r.mode === "pipeline" && activeRun(r));
  const empty = !chosen && !spikein.length ? (pending ? `<div class="empty-state"><h3>分析进行中</h3><button class="table-action" data-run="${escape(pending.id)}">查看进度</button></div>` :
    `<div class="empty-state"><h3>这个样本还没有可复用的分析结果</h3><p>可以先看历史报告，或在工作台启动分析。</p><button class="button secondary" data-go="workspace">前往工作台</button></div>`) : "";
  return picker + cpm + spike + empty;
}

const pairUrl = (sample) => `/api/samples/${encodeURIComponent(sample.id)}`;

function pairAdviceLink(run) {
  const q = run.track_review?.answer?.decision?.qvalue;
  return `<section class="decision-callout" aria-label="待确认的调参建议"><h3>有一条调参建议待你确认</h3>
    <p class="decision-parameter">MACS3 找峰 · q 从 0.01 调到 ${escape(q)}</p>
    <div class="run-actions"><a class="button primary" href="#history/${escape(run.id)}">查看并确认建议</a></div></section>`;
}

function pairResults(sample) {
  const pair = sample.pair;
  const pending = awaitingPairAdvice(sample.id);
  const review = state.runs.find(run => run.paired_review && run.sample_id === sample.id && listedRun(run));
  const comparison = state.runs.filter(run => run.paired_review && run.sample_id === sample.id && run.track_review?.rerun && listedRun(run))
    .sort((a, b) => b.updated_at.localeCompare(a.updated_at))[0];
  const file = (name, label) => `<a href="${pairUrl(sample)}/artifacts/${name}" download>${label}</a>`;
  const peaks = [["baseline", "无 IgG 基线"], ["with-igg", "使用 IgG"]].map(([key, label]) =>
    `<tr><th scope="row">${label}</th><td class="numeric">${count(pair.peak_counts[key])}</td><td class="numeric">${percent(pair.fragment_overlap.target[key].fraction)}</td><td class="numeric">${percent(pair.fragment_overlap.igg[key].fraction)}</td><td>${file(`${key}.peaks.bed`, "下载 BED")}</td></tr>`).join("");
  const tracks = [["target", "目标样本"], ["igg", "IgG 对照"]].map(([role, label]) =>
    `<li><strong>${label}</strong><div class="library-files">${file(`${role}-cpm.bw`, "CPM BigWig（IGV）")}${file(`${role}-tss-profile.png`, "TSS 曲线")}${file(`${role}-tss-heatmap.png`, "TSS 热图")}${file(`${role}-spikein.bw`, "历史 Spike-in BigWig")}</div></li>`).join("");
  return `${pending ? pairAdviceLink(pending) : ""}<section class="library-section" aria-labelledby="pair-peaks-heading">
    <div class="library-section-heading"><h3 id="pair-peaks-heading">候选峰对比</h3>
      ${comparison ? `<a class="button secondary" href="#history/${escape(comparison.id)}">查看最近一次重跑对比</a>` : ""}</div>
    <p class="field-help">基于已过滤 BAM 的本地全基因组计算 · ${escape(pair.peak_settings.execution)} · q = ${escape(pair.peak_settings.qvalue)}</p>
    <div class="table-wrap"><table class="library-peaks"><thead><tr><th scope="col">方案</th><th scope="col" class="numeric">候选峰数</th><th scope="col" class="numeric">目标 FRiP</th><th scope="col" class="numeric">IgG FRiP</th><th scope="col">文件</th></tr></thead><tbody>${peaks}</tbody></table></div>
    <p class="source-note">FRiP 是落在该套候选峰上的片段比例。两套峰的范围不同，峰数和 FRiP 要结合 IGV 峰形判断。这轮试算尚未过滤 blacklist。</p>
  </section><section class="library-section library-spikein" aria-labelledby="pair-tracks-heading">
    <h3 id="pair-tracks-heading">轨道与 TSS 质控</h3><p class="field-help">本地 CPM · hg38 · ${count(pair.bin_size)} bp</p><ul>${tracks}</ul>
  </section><section class="library-section library-section-heading" aria-labelledby="pair-review-heading"><div><h3 id="pair-review-heading">模型复核</h3>
    <p class="field-help">${review ? runLabel(review) : pair.external_model_allowed ? "复用已有分析结果。" : "这组摘要的外部模型使用待授权。"}</p></div>
    ${review ? `<button class="table-action" data-run="${escape(review.id)}">查看复核记录</button>` : '<button class="table-action" data-go="workspace">前往工作台</button>'}</section>`;
}

function pairEvidence(sample) {
  const { spikein, manual_windows: windows, unconfirmed, igv_screen: screen } = sample.pair;
  const removed = screen?.removed;
  return `<section class="library-section" aria-labelledby="pair-spikein-heading">
    <div class="library-section-heading"><h3 id="pair-spikein-heading">Spike-in 加样记录</h3><a href="${pairUrl(sample)}/pair-evidence" download>下载本地证据 JSON</a></div>
    <p>两份样本均在片段化后加入，标称用量为每 ${count(spikein.cells_per_dose)} 个细胞 ${count(spikein.dose_pg)} pg。</p>
    <p class="field-help">实验端补充于 ${escape(spikein.report_date)}。${spikein.actual_cell_counts == null ? "每份样本的实际投入细胞数待确认。" : ""}</p>
    <dl class="library-metrics"><div><dt>目标历史过滤后 Spike-in 片段</dt><dd>${count(spikein.historical_fragments.target)}</dd></div>
      <div><dt>IgG 历史过滤后 Spike-in 片段</dt><dd>${count(spikein.historical_fragments.igg)}</dd></div></dl>
    ${spikein.quantitative_review_pending ? '<p class="source-note">需核对计数口径和实际投入量，再判断能否用于定量标准化。</p>' : ""}
  </section>${windows.length ? `<section class="library-section" aria-labelledby="pair-igv-heading"><h3 id="pair-igv-heading">IGV 复核窗口</h3>
    <p class="field-help">hg38 坐标。表里是与窗口重叠的候选峰数量，供人工检查峰形和背景。</p>
    <div class="table-wrap"><table><thead><tr><th scope="col">区域</th><th scope="col">IGV 坐标</th><th scope="col" class="numeric">无 IgG</th><th scope="col" class="numeric">使用 IgG</th></tr></thead><tbody>${windows.map(window =>
      `<tr><td>${escape(window.label)}</td><td>${escape(window.hg38_igv_locus)}</td><td class="numeric">${count(window.peak_overlap_counts.baseline)}</td><td class="numeric">${count(window.peak_overlap_counts["with-igg"])}</td></tr>`).join("")}</tbody></table></div>
  </section>` : ""}${removed ? `<section class="library-section" aria-labelledby="pair-igv-screen-heading"><h3 id="pair-igv-screen-heading">被筛掉的峰 · IGV 探索复核</h3>
    <p>${screen.igv_examples_reviewed ? `另有 ${count(screen.igv_examples_reviewed)} 个已查看的 IGV 窗口。` : "从两条现有 CPM 轨道自动统计。"}全基因组有 ${count(removed.total)} 个原峰与带 IgG 峰完全不重叠。</p>
    <p class="field-help">峰内平均 CPM：目标达到 IgG 两倍的有 ${count(removed.target_higher)} 个，IgG 达到目标两倍的有 ${count(removed.igg_higher)} 个，两者差距不到两倍的有 ${count(removed.similar)} 个。</p>
  </section>` : ""}<section class="library-section" aria-labelledby="pair-pending-heading"><h3 id="pair-pending-heading">待确认</h3>
    <ul>${unconfirmed.map(label => `<li>${escape(screen && label === "人工 IGV 复核结论" ? "生信同学对 IGV 示例及靶点的确认" : label)}</li>`).join("")}</ul>
    <a href="${pairUrl(sample)}/review-summary" download>下载拟发给模型的摘要</a>
  </section>`;
}

function renderSamples() {
  const sample = currentSample();
  if (!sample) {
    $("#samples-content").innerHTML = '<div class="empty-state"><h2>暂无本地样本</h2><p>准备样本和历史报告后可在此查看。</p></div>';
    return;
  }
  const m = sample.metrics || {};
  const metrics = [
    ["原始测序读段", count(m.total_reads), "fastp · 清理前"],
    ["保留测序读段", count(m.retained_reads), "fastp · 清理后"],
    ["清理前 Q30 占比", percent(m.q30_rate), "fastp · 清理前"],
    ["去重后人源片段", count(m.human_fragments), "summary.tsv"],
    ["spike-in 片段", count(m.spikein_fragments), "summary.tsv"],
    ["过滤片段去重比例", percent(m.duplicate_fraction_filtered), "summary.tsv"],
  ];
  const list = `<nav class="library-index" aria-label="本地样本列表"><ul>${state.samples.map(item => {
    const cpm = completedFor(item.id,"pipeline","completed_pipeline").length;
    const spike = completedFor(item.id,"spikein","completed_spikein").length;
    const status = item.pair ? "配对分析结果可查看" : cpm && spike ? "CPM 与 Spike-in 可查看" : cpm || spike ? "分析结果可查看" : "暂无已完成结果";
    return `<li><button data-library-sample="${escape(item.id)}" ${item.id === sample.id ? 'aria-current="true"' : ""}><strong>${escape(item.label)}</strong><small>${status}</small></button></li>`;
  }).join("")}</ul></nav>`;
  $("#samples-content").innerHTML = `<div class="library-layout">${list}<div class="library-sample"><div class="sample-heading"><h2>${escape(sample.label)}</h2>${sample.pair ? `<p class="field-help">匹配 IgG 为 ${escape(sample.pair.control_id)} · 依据为${escape(sample.pair.control_source)}</p>` : ""}</div>
    <div class="report-tabs" role="tablist" aria-label="样本详情">
      <button class="report-tab" id="sample-tab-results" role="tab" aria-controls="sample-panel-results" data-sample-tab="results">分析结果</button>
      <button class="report-tab" id="sample-tab-report" role="tab" aria-controls="sample-panel-report" data-sample-tab="report">${sample.pair ? "本地证据" : "历史报告"}</button>
      <button class="report-tab" id="sample-tab-files" role="tab" aria-controls="sample-panel-files" data-sample-tab="files">输入文件</button>
    </div>
    <section id="sample-panel-results" role="tabpanel" aria-labelledby="sample-tab-results" tabindex="0" hidden>${sample.pair ? pairResults(sample) : sampleResults(sample)}</section>
    <section id="sample-panel-report" role="tabpanel" aria-labelledby="sample-tab-report" tabindex="0" hidden>
      ${sample.pair ? pairEvidence(sample) : `<p class="source-note">来自本地历史报告，不是本次重新分析的结果；不据此判定实验是否合格。</p>
      <div class="report-layout"><div class="table-wrap"><table class="metrics-table"><caption>历史质控指标</caption><thead><tr><th scope="col">指标</th><th scope="col" class="numeric">数值</th><th scope="col">来源</th></tr></thead><tbody>${metrics.map(([label, value, source]) => `<tr><td>${label}</td><td class="numeric">${value}</td><td class="metric-source">${source}</td></tr>`).join("")}</tbody></table></div>
      <section class="chart-section" aria-labelledby="fragment-heading"><h2 id="fragment-heading">历史片段长度分布</h2><div id="fragment-chart" class="fragment-chart">${chartMarkup(sample)}</div><p class="source-note">来源：hg38.dedup.fragment_length_counts.txt<br>10 bp 分箱，仅显示 0–999 bp。</p></section></div>`}
    </section>
    <section id="sample-panel-files" role="tabpanel" aria-labelledby="sample-tab-files" tabindex="0" hidden>
      <div class="files-heading"><h2>输入文件</h2><p>${sample.files.length} 个源文件 · ${bytes(sample.total_bytes)}</p></div>
      ${sample.pair ? "" : `<p class="verification-note">${sample.verification_status === "passed" ? `基本校验：${sample.verification_checks} 项通过（已保存记录）。` : "基本校验：待检查。"}<span class="muted">仅表示传输与格式检查。</span></p>`}
      <div class="table-wrap"><table class="files-table"><thead><tr><th scope="col">文件</th><th scope="col">类型</th><th scope="col" class="numeric">大小</th><th scope="col">大小核对</th></tr></thead><tbody>${sample.files.map((file) => `<tr><td class="file-name">${escape(file.name)}</td><td>${escape({ raw: "原始数据", processed: "已处理数据", qc: "质控报告" }[file.group] || file.group)}</td><td class="numeric">${bytes(file.bytes)}</td><td class="file-status ${file.present && file.size_matches ? "" : "error"}">${!file.present ? "文件缺失" : file.size_matches ? "匹配" : "大小不匹配"}</td></tr>`).join("")}</tbody></table></div>
      <p class="table-note">对照来源清单核对大小；原始数据不提供网页下载。</p>
    </section></div></div>`;
  setSampleTab(state.sampleTab);
}

function renderPhaseProgress(title, done, total, phase) {
  const segments = Array.from({length: total}, (_, index) => {
    const status = index < done ? "phase-done" : index === done ? "phase-current" : "";
    return `<span class="${status}"></span>`;
  }).join("");
  return `<div class="phase-progress">
    <div class="phase-progress-heading"><strong>${title}</strong><span>已完成 ${done} / ${total} 个阶段</span></div>
    <div class="phase-progress-track" role="progressbar" aria-label="${title}" aria-valuemin="0"
      aria-valuemax="${total}" aria-valuenow="${done}" aria-valuetext="已完成 ${done} / ${total} 个阶段，${escape(phase)}">${segments}</div>
    <p>${escape(phase)}</p>
  </div>`;
}

function renderPipelineProgress(run) {
  const finished = run.status === "completed_pipeline" &&
    (!run.model_review_requested || run.track_review?.status === "completed_review");
  if (!activeRun(run) && !finished) return "";
  const steps = run.steps.filter(step => step.status !== "skipped");
  const total = steps.length + Number(!!run.model_review_requested);
  if (!total) return "";
  const done = finished ? total : steps.filter(step => step.status === "done").length;
  const current = steps.find(step => step.status === "running");
  const phase = finished ? "本次分析已完成" : run.status === "queued" ? "等待执行" : run.baseline_ready ? "Agent 正在复核质控结果" :
    current?.label || "正在准备输入";
  return renderPhaseProgress(run.model_review_requested ? "完整分析进度" : "基础分析进度", done, total, phase);
}

function renderPairedProgress(run) {
  if (!run.paired_adaptive) return "";
  const review = run.track_review;
  if (run.paired_rerun_accepted || review?.rerun_started) {
    if (!activeRun(run)) {
      return run.status === "completed_review" && review?.reassessment && run.rerun_figure ?
        renderPhaseProgress("调参重跑进度", 3, 3, "重跑、复评和对照图均已完成") : "";
    }
    const done = review?.reassessment ? 2 : review?.rerun ? 1 : 0;
    const phase = run.status === "queued" ? "等待重跑" : done === 2 ? "正在生成前后对照图" :
      done === 1 ? "Agent 正在复评重跑结果" : review?.rerun_started ? "正在重跑 MACS3 找峰并统计结果" : "正在核对重跑输入";
    return renderPhaseProgress("调参重跑进度", done, 3, phase);
  }
  const waiting = run.status === "awaiting_rerun";
  if (!waiting && !activeRun(run)) return "";
  const ready = run.evidence_preparing === false || !!review?.evidence;
  const phase = waiting ? "等待你确认是否重跑" : run.status === "queued" ? "等待执行" :
    ready ? "Agent 正在根据证据选择下一步" : "正在整理轨道与候选峰证据";
  return renderPhaseProgress("首轮分析进度", waiting ? 2 : Number(ready), 2, phase);
}

function renderSingleReviewProgress(run) {
  const review = run.track_review, finished = run.status === "completed_review";
  if (!activeRun(run) && !finished) return "";
  const action = review?.answer?.decision?.action;
  if (action?.startsWith("compare_") && (!finished || review.reassessment)) {
    const done = finished ? 2 : Number(!!review.comparison);
    const phase = finished ? "比较与复评已完成" : done ? "Agent 正在复评比较结果" : "正在比较已有分析结果";
    return renderPhaseProgress("比较与复评进度", done, 2, phase);
  }
  const ready = !!review?.evidence || !!review?.summary;
  const phase = finished ? "首轮分析已完成" : run.status === "queued" ? "等待执行" :
    ready ? "Agent 正在根据证据选择下一步" : "正在整理轨道与候选峰证据";
  return renderPhaseProgress("首轮分析进度", finished ? 2 : Number(ready), 2, phase);
}

function renderPipeline(run, content) {
  const summary = run.result;
  const open = content.dataset.runId === run.id ? Array.from(content.querySelectorAll("details[open]")).map(n => n.dataset.detail) : [];
  const link = name => `/api/runs/${encodeURIComponent(run.id)}/artifacts/${name}`;
  const files = new Map((run.artifacts || []).map(f => [f.name, f]));
  const fileRow = (name, label) => files.has(name) ? `<tr><td>${escape(label)}</td><td>${bytes(files.get(name).bytes)}</td><td><a href="${link(name)}" download>下载</a></td></tr>` : "";
  const stepLabels = {pending:"未开始",running:"进行中",done:"已完成",skipped:"已跳过",cancelled:"已停止",failed:"失败"};
  const qc = summary?.qc, reportMetrics = run.report_metrics || {};
  content.dataset.runId = run.id;
  content.innerHTML = `<div class="run-topline"><h3>${escape(run.name)}</h3><p>${escape(run.sample_id)} · ${run.source === "fastq" ? "原始 FASTQ" : run.source === "resume" ? "本轮生成的 BAM" : "已有去重 BAM"} · ${stamp(run.created_at)}</p></div>
    ${renderPipelineProgress(run)}
    <p class="run-disclaimer" role="status">${escape(run.message)}</p>
    ${run.error ? `<p class="inline-error" role="alert">${escape(run.error)}</p>` : ""}
    ${state.pollError ? `<p class="inline-error" role="alert">${escape(state.pollError)}</p>` : ""}
    ${renderModelTrace(run, open)}
    ${summary ? `<p class="source-note">本次真实计算。候选峰计算未使用 IgG 对照；靶点和实验设计待确认，文件检查通过不代表实验合格。</p>
      <div class="table-wrap"><table class="pipeline-metrics"><thead><tr><th>本次统计</th><th>结果</th></tr></thead><tbody>
      ${reportMetrics.raw_reads != null ? `<tr><td>清理前 / 保留读段</td><td>${count(reportMetrics.raw_reads)} / ${count(reportMetrics.retained_reads)}</td></tr><tr><td>清理前 Q30</td><td>${percent(reportMetrics.raw_q30_rate)}</td></tr><tr><td>主比对率（清理后${run.source === "resume" ? "，原 FASTQ 任务" : ""}）</td><td>${percent(reportMetrics.primary_mapping_rate)}</td></tr>` : ""}
      <tr><td>配对片段</td><td>${count(qc.fragments)}</td></tr>
      <tr><td>重复标记比例</td><td>${percent(qc.duplicate_fraction)}${run.source === "bam" ? "（输入已去重）" : "（本次保留重复）"}</td></tr>
      <tr><td>窄峰候选 / 片段 FRiP</td><td>${count(qc.peaks.narrow.count)} / ${percent(qc.frip.narrow.fraction)}</td></tr>
      <tr><td>宽峰候选 / 片段 FRiP</td><td>${count(qc.peaks.broad.count)} / ${percent(qc.frip.broad.fraction)}</td></tr></tbody></table></div>
      <p class="field-help">FRiP：与相应候选峰重叠的片段占比。当前仅描述本次结果，不按峰数或 FRiP 自动选优。</p>
      ${run.status === "completed_pipeline" && summary.validation === "passed" ?
        `<section class="rerun-visuals">${renderIgvTracks(run.id, "可在同一坐标下查看 CPM 信号、候选窄峰和宽峰。")}</section>` : ""}
      <h3 class="result-heading">分析文件</h3><div class="table-wrap"><table><thead><tr><th>文件用途</th><th>大小</th><th>操作</th></tr></thead><tbody>
      ${fileRow("coverage.cpm.bw", "全基因组 CPM 轨道 · BW")}${fileRow("narrow.peaks.bed", "窄峰候选 · BED")}${fileRow("broad.peaks.bed", "宽峰候选 · BED")}${fileRow("report.md", "结果说明")}${fileRow("fastp.html", "本次 fastp 报告（HTML 文本）")}
      </tbody></table></div>
      <section class="chart-section"><h3>片段长度分布</h3><img class="result-plot" src="${link("fragment-length.png")}" alt="本次分析的配对片段长度分布，范围 0 到 1000 bp"></section>
      <section class="chart-section"><h3>基因 TSS 周围信号</h3><p class="field-help">GENCODE v44 基因级 TSS ±3 kb，CPM。仅供描述，不适用于所有靶点的统一合格判断。</p><img class="result-plot" src="${link("tss.profile.png")}" alt="本次计算的基因 TSS 周围平均 CPM 信号曲线"></section>
      <details data-detail="tss" ${open.includes("tss") ? "open" : ""}><summary>TSS 热图、矩阵和完整文件</summary><img class="result-plot heatmap" src="${link("tss.heatmap.png")}" alt="基因 TSS 周围信号热图"><div class="table-wrap"><table><thead><tr><th>文件</th><th>大小</th><th>操作</th></tr></thead><tbody>${[...files.keys()].filter(name => !["coverage.cpm.bw","narrow.peaks.bed","broad.peaks.bed","report.md","fastp.html"].includes(name)).map(name => fileRow(name,name)).join("")}</tbody></table></div></details>
      <details data-detail="limits" ${open.includes("limits") ? "open" : ""}><summary>参数与分析限制</summary><p class="field-help">BW：25 bp、CPM、配对片段覆盖。MACS3：BAMPE、q=0.01，宽峰 cutoff=0.1。TSS：50 bp 窗口，保留缺失值。</p><ul>${summary.warnings.map(w => `<li>${escape(w)}</li>`).join("")}</ul></details>` : ""}
    <div class="run-actions">${activeRun(run) ? `<button class="button danger" data-cancel="${escape(run.id)}" ${state.stopping.has(run.id) ? "disabled" : ""}>${state.stopping.has(run.id) ? "正在停止…" : "停止计算"}</button>` : ""}<a class="button secondary" href="/api/runs/${encodeURIComponent(run.id)}/report" download>下载运行记录</a>${files.has("commands.json") ? `<a class="button secondary" href="${link("commands.json")}" download>下载工具命令记录</a>` : ""}</div>
    ${files.has("tool.log") ? logPreview(link("tool.log"), "pipeline-log", open) : ""}
    ${files.has("bowtie2.log") ? logPreview(link("bowtie2.log"), "align-log", open, "比对日志") : ""}
    ${run.steps.length ? `<details data-detail="steps" ${open.includes("steps") || activeRun(run) ? "open" : ""}><summary>运行步骤</summary><ol class="timeline">${run.steps.map((step,index) => `<li class="timeline-step ${escape(step.status)}"><span class="step-number">${index+1}</span><div><strong>${escape(step.label)}</strong><span class="step-state">${stepLabels[step.status] || escape(step.status)}</span></div></li>`).join("")}</ol></details>` : ""}
    ${run.events.length ? `<details data-detail="events" ${open.includes("events") ? "open" : ""}><summary>事件记录</summary><ul class="event-list">${run.events.map(event => `<li><time>${stamp(event.at,true)}</time><span>${escape(event.message)}</span></li>`).join("")}</ul></details>` : ""}`;
}

function renderSpikein(run, content) {
  const base = `/api/runs/${encodeURIComponent(run.id)}`;
  const open = content.dataset.runId === run.id ? Array.from(content.querySelectorAll("details[open]")).map(n => n.dataset.detail) : [];
  content.dataset.runId = run.id;
  content.innerHTML = `<div class="run-topline"><h3>${escape(run.name)}</h3><p>${escape(run.sample_id)} · ${stamp(run.created_at)}</p></div>
    ${run.error ? `<p class="inline-error" role="alert">${escape(run.error)}</p>` : `<dl class="run-summary">
      <div><dt>Spike-in 片段</dt><dd>${count(run.spikein_fragments)}</dd></div>
      <div><dt>缩放系数</dt><dd>${count(run.scale_numerator)} ÷ ${count(run.spikein_fragments)} = ${Number(run.scale_factor).toFixed(6)}</dd></div>
      <div><dt>轨道</dt><dd>hg38 全基因组 · ${escape(run.result?.bin_size)} bp · 仅 Spike-in 缩放</dd></div></dl>
      <p class="source-note">使用与 spike-in 配套的历史人源 BAM 重新生成轨道。本次轨道计算未执行有 IgG 对照的找峰。</p>`}
    ${renderModelTrace(run)}
    <div class="run-actions">${run.status === "completed_spikein" ? `<a class="button secondary" href="${base}/artifacts/coverage.spikein.bw" download>下载 BigWig</a>
      <a class="button secondary" href="${base}/artifacts/report.md" download>下载结果说明</a>` : ""}
      <a class="button secondary" href="${base}/report" download>下载运行记录</a></div>
    ${run.artifacts?.some(f => f.name === "tool.log") ? logPreview(`${base}/artifacts/tool.log`, "spikein-log", open) : ""}`;
}

function reviewActionCard(run) {
  const id = escape(run.id), review = run.track_review;
  const q = review?.answer?.decision?.qvalue;
  const busy = state.deciding.has(run.id);
  const disabled = busy ? "disabled" : "";
  if (run.status === "awaiting_rerun") {
    return `<section class="decision-callout" aria-labelledby="rerun-choice-${id}">
      <h3 id="rerun-choice-${id}">待你确认 · 重跑带 IgG 的找峰</h3>
      <p class="decision-parameter">MACS3 找峰 · q 从 0.01 调到 ${escape(q)}</p>
      <p>接受后，在同一条任务里查看新结果和模型复评。</p>
      <div class="run-actions" aria-busy="${busy}">
        <button class="button primary" data-decide-rerun="${id}" data-accept="true" ${disabled}>接受建议并重跑</button>
        <button class="button secondary" data-decide-rerun="${id}" data-accept="false" ${disabled}>保留原结果</button>
      </div></section>`;
  }
  const lastDecision = review?.reassessment?.decision || review?.answer?.decision;
  if (run.status === "completed_review" && lastDecision?.action === "request_evidence") {
    const pending = run.paired_review && awaitingPairAdvice(run.sample_id);
    return `<section class="decision-callout context-only" aria-label="这条建议需要补证据"><h3>这条建议需要补证据</h3>
      <p>模型给出了要核对的信息。重跑入口会在 Agent 提出具体的工具参数后出现。</p>
      ${run.paired_review ? `<div class="run-actions">${pending ? `<a class="button primary" href="#history/${escape(pending.id)}">查看待确认的调参建议</a>` : ""}
        <button class="button secondary" data-open-pair-agent="${escape(run.sample_id)}">前往启动 Agent</button></div>` : ""}</section>`;
  }
  return "";
}

function signedMetric(value, digits) {
  const rounded = Number(value.toFixed(digits));
  const sign = rounded > 0 ? "+" : rounded < 0 ? "−" : "";
  return sign + Math.abs(rounded).toLocaleString("zh-CN", {minimumFractionDigits: digits, maximumFractionDigits: digits});
}

function renderRerunComparison(run) {
  const review = run.track_review, before = review.comparison, after = review.rerun;
  const shifts = after.igv_changes, restored = shifts?.removed_now_overlapping, lost = shifts?.overlapping_now_missing;
  const id = escape(run.id);
  const originalBed = `/api/samples/${encodeURIComponent(run.sample_id)}/artifacts/with-igg.peaks.bed`;
  const newArtifact = name => `/api/runs/${encodeURIComponent(run.id)}/artifacts/${name}`;
  const beforeQ = review.evidence?.settings?.qvalue;
  const decision = review.reassessment?.decision;
  const canReplay = run.paired_adaptive && run.status === "completed_review" && review.peak_executions === 1;
  const metrics = [
    {label:"候选峰数", before:before.peak_count, after:after.peak_count, format:count, digits:0, unit:"个"},
    {label:"目标落峰比例", before:before.target.frip, after:after.target.frip, format:percent, digits:2, unit:"百分点", scale:100},
    {label:"IgG 落峰比例", before:before.igg.frip, after:after.igg.frip, format:percent, digits:2, unit:"百分点", scale:100},
  ];
  const beforeIgg = before.igv_screen?.overlaps_controlled?.median_igg_cpm;
  if (Number.isFinite(beforeIgg) && Number.isFinite(shifts?.median_igg_cpm_in_rerun_overlap)) {
    metrics.push({label:"重叠区间的 IgG 中位信号", before:beforeIgg, after:shifts.median_igg_cpm_in_rerun_overlap,
      format:value=>value.toFixed(4), digits:4, unit:"CPM"});
  }

  const rows = metrics.map(metric => {
    const max = Math.max(metric.before, metric.after) || 1;
    const valueCell = value => `<strong>${metric.format(value)}</strong><meter aria-hidden="true" min="0" max="${max}" value="${value}"></meter>`;
    const delta = signedMetric((metric.after - metric.before) * (metric.scale ?? 1), metric.digits);
    return `<tr><th scope="row">${metric.label}</th><td class="comparison-before">${valueCell(metric.before)}</td>
      <td class="comparison-after">${valueCell(metric.after)}</td><td class="comparison-delta"><strong>${delta}</strong><small>${metric.unit}</small></td></tr>`;
  }).join("");

  return `<section class="rerun-comparison" aria-labelledby="rerun-comparison-${id}">
    <div class="rerun-comparison-heading"><h3 id="rerun-comparison-${id}">重跑前后对比</h3>
      ${canReplay ? `<button class="button secondary" data-replay-rerun="${id}" ${state.deciding.has(run.id) ? "disabled" : ""}>用上次参数重跑</button>` : ""}</div>
    <div class="table-wrap"><table class="rerun-metrics"><thead><tr><th scope="col">指标</th>
      <th scope="col">重跑前 <small>q = ${escape(beforeQ ?? "—")}</small></th>
      <th scope="col">重跑后 <small>q = ${escape(after.qvalue)}</small></th><th scope="col">变化</th></tr></thead><tbody>${rows}</tbody></table></div>
    ${restored ? `<div class="rerun-intervals"><h4>原峰区间的保留变化</h4><p>找回 <strong>${count(restored.total)}</strong> 个区间，其中 <strong>${count(restored.target_higher)}</strong> 个目标信号偏强、<strong>${count(restored.igg_higher)}</strong> 个 IgG 信号偏强。${lost ? `此前保留的区间有 <strong>${count(lost.total)}</strong> 个在重跑后失去。` : ""}</p></div>` : ""}
    ${decision ? `<p class="rerun-verdict"><span>模型复评</span><strong>${escape(decisionLabel(decision))}</strong></p>` : ""}
    <div class="library-files rerun-downloads"><a href="${originalBed}" download>重跑前 BED</a><a href="${newArtifact("tuned.peaks.bed")}" download>重跑后 BED</a>
      <a href="${newArtifact("commands.json")}" download>工具命令</a><a href="${newArtifact("tool.log")}" download>计算日志</a></div>
  </section>`;
}

function renderIgvTracks(runId, description, showViewer = true) {
  const id = encodeURIComponent(runId);
  return `<div class="rerun-visuals-heading"><h3>IGV 轨道对照</h3>
    <a class="button secondary" href="/api/runs/${id}/igv-session" download>下载 IGV Desktop 会话</a></div>
    ${showViewer ? `<p class="field-help">${escape(description)}</p>
      <iframe class="igv-comparison" src="/igv-view.html?run=${id}" title="IGV 基因组轨道对照" loading="lazy"></iframe>` : ""}`;
}

function renderRerunImages(run, open) {
  const figure = run.rerun_figure;
  const figureUrl = figure ? `/api/runs/${encodeURIComponent(run.id)}/artifacts/${encodeURIComponent(figure.name)}` : "";
  const original = `/api/samples/${encodeURIComponent(run.sample_id)}/artifacts`;
  const originalImages = [["target-tss-heatmap.png", "目标 TSS 热图"], ["igg-tss-heatmap.png", "IgG TSS 热图"]];
  return `<section class="rerun-visuals" aria-label="找峰图像对比">
    ${renderIgvTracks(run.id, "目标和 IgG 共用纵轴，两轮候选峰各占一条轨道。可以切换代表区间，也可以缩放或输入坐标。", !activeRun(run))}
    ${figure ? `<details class="original-qc-images" data-detail="pair-static-figure" ${open.includes("pair-static-figure") ? "open" : ""}>
      <summary>静态对照图</summary>
      <p class="field-help">左边是重跑前，右边是重跑后。每行使用相同坐标与纵轴。</p>
      <a href="${figureUrl}" download>下载 PNG</a>
      <figure class="track-comparison"><div class="track-comparison-image"><img src="${figureUrl}" alt="代表区间的重跑前后对照。每行显示同一对目标和 IgG 覆盖轨道，下方标出各轮候选峰。" loading="lazy"></div></figure></details>` :
      run.track_review?.figure_error ? `<p class="field-help">${escape(run.track_review.figure_error)}</p>` : activeRun(run) ? '<p class="field-help" role="status">正在准备轨道对照图。</p>' : ""}
    <details class="original-qc-images" data-detail="pair-original-images" ${open.includes("pair-original-images") ? "open" : ""}>
      <summary>原始质控图</summary><p class="field-help">两轮使用同一对轨道。这些是原有的 TSS 热图。</p>
      <div class="original-qc-grid">${originalImages.map(([name, label]) => `<figure><figcaption>${label}</figcaption>
        <a href="${original}/${name}" download><img src="${original}/${name}" alt="${label}" loading="lazy"></a></figure>`).join("")}</div></details>
  </section>`;
}

const reviewStep = (title, body) => `<section class="review-step"><h3>${title}</h3>${body}</section>`;
const reviewMetrics = rows => `<dl class="library-metrics">${rows.map(([label, value]) =>
  `<div><dt>${label}</dt><dd>${value}</dd></div>`).join("")}</dl>`;

function renderReviewLayout(run, content, {open, progress, sections, baselineButton, extra = ""}) {
  const id = escape(run.id), review = run.track_review;
  const events = run.events?.length ? run.events : [{at: run.created_at, message: "从已有全基因组结果发起模型复核。"},
    ...(activeRun(run) ? [] : [{at: run.updated_at, message: run.message}])];
  content.dataset.runId = run.id;
  content.innerHTML = `<div class="run-topline"><h3>${escape(run.name)}</h3><p>${escape(run.sample_id)} · ${stamp(run.created_at)}</p></div>
    ${progress}
    <p class="run-disclaimer" role="status">${escape(run.message)}</p>
    <div class="review-steps">${sections.join("")}</div>
    ${review?.evidence ? renderReviewEvidence(review, open) : ""}
    ${review?.error && review.answer ? `<p class="inline-error" role="alert">${escape(review.error)}</p>` : ""}
    ${state.runError || state.pollError ? `<p class="inline-error" role="alert">${escape(state.runError || state.pollError)}</p>` : ""}
    <div class="run-actions">${activeRun(run) ? `<button class="button danger" data-cancel="${id}" ${state.stopping.has(run.id) ? "disabled" : ""}>停止任务</button>` : ""}
      ${baselineButton}
      <a class="button secondary" href="/api/runs/${encodeURIComponent(run.id)}/report" download>下载运行记录</a></div>
    <details data-detail="review-events" ${open.includes("review-events") ? "open" : ""}><summary>本次执行记录</summary>
      <ul class="event-list">${events.map(event=>`<li><time>${stamp(event.at,true)}</time><span>${escape(event.message)}</span></li>`).join("")}</ul>
      ${review?.elapsed_seconds != null ? `<p class="field-help">耗时 ${review.elapsed_seconds.toFixed(1)} 秒，模型请求 ${count(run.actual_model_requests)} 次。</p>` : ""}
    </details>${extra}`;
}

function renderPairedTrackReview(run, content, open) {
  const review = run.track_review;
  const sample = state.samples.find(item => item.id === run.sample_id);
  const pair = sample?.pair, decision = review?.answer?.decision;
  const afterReview = review?.reassessment?.decision;
  const before = review?.comparison || review?.evidence?.comparison;
  const beforeQ = review?.evidence?.settings?.qvalue ?? pair?.peak_settings?.qvalue;
  const screen = before?.igv_screen;
  const original = review?.evidence?.baseline;
  const beforeBed = `/api/samples/${encodeURIComponent(run.sample_id)}/artifacts/with-igg.peaks.bed`;
  const evidence = decision?.evidence_used || [];
  let confirmation = reviewActionCard(run);
  if (!confirmation && decision?.action === "tune_igg_qvalue") {
    confirmation = `<p>${run.paired_rerun_accepted === false ? "你选择了保留原结果。" : run.reused_plan_id ? "你接受了沿用原参数重跑。" : run.paired_rerun_accepted || review?.rerun_started ? "已接受建议，按建议参数执行。" : "任务结束，尚未执行重跑。"}</p>`;
  }
  let result = "";
  if (review?.rerun) {
    result = renderRerunComparison(run) + renderRerunImages(run, open);
  } else if (original && before && decision?.action === "compare_igg") {
    const rows = [["无 IgG", original], ["使用 IgG", before]];
    result = `<div class="table-wrap"><table class="library-peaks"><thead><tr><th>方案</th><th>候选峰数</th><th>目标落峰比例</th><th>IgG 落峰比例</th></tr></thead>
      <tbody>${rows.map(([label,row])=>`<tr><th scope="row">${label}</th><td>${count(row.peak_count)}</td><td>${percent(row.target.frip)}</td><td>${percent(row.igg.frip)}</td></tr>`).join("")}</tbody></table></div>`;
  } else if (decision && !activeRun(run) && run.status !== "awaiting_rerun" && run.status !== "failed") {
    result = `<p>${decision.action === "request_evidence" ? "本次复核的证据已保存。" : "本次使用已有找峰结果。"}</p>`;
  }
  if (afterReview) {
    result += `<details data-detail="pair-final-judgment" ${open.includes("pair-final-judgment") ? "open" : ""}><summary>查看模型复评理由</summary>
      <p>${escape(modelText(afterReview.reason))}</p></details>`;
  }
  renderReviewLayout(run, content, {
    open,
    progress: renderPairedProgress(run),
    sections: [
      reviewStep("本次分析", `<p>以已有过滤 BAM 和全基因组候选峰为起点，复核带 IgG 的 CUT&Tag 找峰。</p>
        ${reviewMetrics([["目标样本",escape(run.sample_id)], ["IgG 对照",escape(pair?.control_id || "—")]])}`),
      reviewStep("起始结果", reviewMetrics([["找峰阈值 q",escape(beforeQ ?? "—")], ["带 IgG 候选峰",count(before?.peak_count ?? pair?.peak_counts?.["with-igg"])],
        ["目标落峰比例",percent(before?.target?.frip ?? pair?.fragment_overlap?.target?.["with-igg"]?.fraction)],
        ["IgG 落峰比例",percent(before?.igg?.frip ?? pair?.fragment_overlap?.igg?.["with-igg"]?.fraction)]]) + `<a href="${beforeBed}" download>下载起始 BED</a>`),
      reviewStep("模型依据与建议", decision ? `<p class="decision-result"><strong>${escape(decisionLabel(decision))}</strong></p>
        <p>${escape(modelText(decision.reason))}</p>
        ${evidence.length ? `<ul class="decision-evidence">${evidence.map(key=>`<li>${escape(evidenceLabels[key] || key)}</li>`).join("")}</ul>` : ""}
        ${screen ? `<p class="field-help">自动统计覆盖 ${count(screen.removed.total)} 个原先筛掉的区间，其中 ${count(screen.removed.target_higher)} 个目标信号偏强，${count(screen.removed.igg_higher)} 个 IgG 信号偏强。</p>` : ""}
        ${run.reused_plan_id ? `<p><a href="#history/${escape(run.reused_plan_id)}">查看沿用的建议</a></p>` : ""}` :
        `<p>${escape(review?.error || "正在准备证据与模型判断。")}</p>`),
      confirmation ? reviewStep("人工决定", confirmation) : "",
      result ? reviewStep(review?.rerun ? "接受后的结果" : "分析结果", result) : "",
    ],
    baselineButton: `<button class="button secondary" data-pair-tab="results" data-sample-id="${escape(run.sample_id)}">查看原分析结果</button>`,
  });
}

function renderSingleComparison(review) {
  const comparison = review.comparison;
  const table = (labels, rows) => `<div class="table-wrap"><table class="library-peaks">
    <thead><tr><th scope="col">指标</th>${labels.map(label=>`<th scope="col" class="numeric">${escape(label)}</th>`).join("")}</tr></thead>
    <tbody>${rows.map(([label, first, second])=>`<tr><th scope="row">${label}</th><td class="numeric">${first}</td><td class="numeric">${second}</td></tr>`).join("")}</tbody>
  </table></div>`;
  const narrow = review.evidence?.qc?.narrow_peaks, broad = review.evidence?.qc?.broad_peaks;
  if (comparison.action === "compare_peak_sets" && narrow && broad) {
    return table(["窄峰", "宽峰"], [
      ["候选峰数", count(narrow.peaks), count(broad.peaks)],
      ["峰宽中位数（bp）", count(narrow.median_width_bp), count(broad.median_width_bp)],
      ["片段落峰比例", percent(narrow.frip), percent(broad.frip)],
    ]);
  }
  const core = review.evidence?.core || review.summary, before = core?.bw_statistics, after = comparison.statistics;
  if (before && after) {
    const mean = value => Number.isFinite(value) ? value.toFixed(4) : "—";
    return table([`${core.parameters.bin_size} bp 原窗口`, `${comparison.bin_bp ?? comparison.window_bp} bp 聚合窗口`], [
      ["有值窗口", count(before.interval_count), count(after.interval_count)],
      ["非零窗口", count(before.nonzero_intervals), count(after.nonzero_intervals)],
      ["平均信号", mean(before.preview_mean), mean(after.preview_mean)],
    ]);
  }
  return "";
}

function renderSingleTrackReview(run, content, open) {
  const baseline = state.runs.find(item => item.id === run.baseline_id && item.sample_id === run.sample_id);
  const review = run.track_review, decision = review?.answer?.decision;
  const qc = baseline?.result?.qc, evidence = review?.evidence?.qc;
  const bin = baseline?.result?.coverage?.bin_size ?? review?.evidence?.core?.parameters?.bin_size;
  const spikein = completedFor(run.sample_id, "spikein", "completed_spikein")[0];
  const noIgg = evidence?.experiment_context?.matched_igg === "unavailable" || spikein?.matched_igg === false;
  const files = new Set((baseline?.artifacts || []).map(f => f.name));
  const fileLink = name => `/api/runs/${encodeURIComponent(run.baseline_id)}/artifacts/${name}`;
  const downloadLinks = [["coverage.cpm.bw", "CPM BigWig"], ["narrow.peaks.bed", "窄峰 BED"], ["broad.peaks.bed", "宽峰 BED"],
    ["report.md", "分析报告"], ["tss.profile.png", "TSS 曲线"]].filter(([name]) => files.has(name));
  const missing = decision?.missing_information || [];
  const corrected = noIgg && missing.some(item => /igg/i.test(item));
  const proposal = decision ? `<p class="decision-result"><strong>${escape(decisionLabel(decision))}</strong></p>
    ${corrected ? '<p class="field-help">此历史回复曾把 IgG 列作待补证据。该样本已确认无匹配 IgG。</p>' : ""}
    ${decision.assessment && decision.assessment !== decision.reason ? `<p>${escape(modelText(decision.assessment))}</p>` : ""}
    <p>${escape(modelText(decision.reason))}</p>
    ${decision.evidence_used?.length ? `<p class="field-help">引用证据 · ${escape(decision.evidence_used.map(key => evidenceLabels[key] || key).join("、"))}</p>` : ""}
    ${missing.length ? `<p class="field-help">${corrected ? "历史回复所列信息" : "待确认"} · ${escape(missing.map(key => Object.hasOwn(contextLabels, key) ? contextLabels[key] : key).join("、"))}</p>` : ""}` :
    `<p>${escape(review?.error || (activeRun(run) ? "正在准备证据与模型判断。" : run.message))}</p>`;
  const after = review?.reassessment?.decision;
  let result = review?.comparison ? `<p>${escape(comparisonEffect(review.comparison))}</p>${renderSingleComparison(review)}` :
    run.status === "completed_review" ? "<p>已保留原 CPM 轨道和候选峰。</p>" : "";
  if (after) result += `<p class="decision-result"><strong>模型复评 · ${escape(decisionLabel(after))}</strong></p><p>${escape(modelText(after.reason))}</p>`;
  if (baseline && !activeRun(run) && ["coverage.cpm.bw", "narrow.peaks.bed", "broad.peaks.bed"].every(name => files.has(name))) {
    result += renderIgvTracks(run.id, "在同一坐标下查看这份样本的 CPM 信号、候选窄峰和宽峰。");
  }
  renderReviewLayout(run, content, {
    open,
    progress: renderSingleReviewProgress(run),
    sections: [
      reviewStep("本次分析", `<p>复核已有的全基因组 CPM、候选峰与质控证据。</p>
        ${reviewMetrics([["目标样本", escape(run.sample_id)], ["IgG 对照", noIgg ? "无匹配 IgG" : "匹配情况待核对"]])}`),
      reviewStep("起始结果", reviewMetrics([
        ["CPM 窗口", bin ? `${escape(bin)} bp` : "—"],
        ["配对片段数", count(qc?.fragments ?? evidence?.fragment_qc?.fragments)],
        ["候选窄峰", count(qc?.peaks?.narrow?.count ?? evidence?.narrow_peaks?.peaks)],
        ["候选宽峰", count(qc?.peaks?.broad?.count ?? evidence?.broad_peaks?.peaks)],
      ]) + (downloadLinks.length ? `<div class="library-files">${downloadLinks.map(([name, label]) =>
        `<a href="${fileLink(name)}" download>${label}</a>`).join("")}</div>` : "")),
      reviewStep("模型依据与建议", proposal),
      result ? reviewStep("分析结果", result) : "",
    ],
    baselineButton: baseline ? `<button class="button secondary" data-run="${escape(baseline.id)}">查看原分析结果</button>` : "",
    extra: files.has("tool.log") ? logPreview(fileLink("tool.log"), "baseline-tool-log", open, "原分析工具日志") : "",
  });
}

function renderTrackReview(run, content) {
  const sameRun = content.dataset.runId === run.id;
  const mayShowIgv = run.paired_review ? !!run.track_review?.rerun : !!run.baseline_id;
  const viewKey = !activeRun(run) && mayShowIgv ?
    JSON.stringify([run.updated_at, run.status, run.rerun_figure?.bytes, state.runError, state.pollError, state.deciding.has(run.id)]) : "";
  if (sameRun && viewKey && content.dataset.igvViewKey === viewKey) return;
  content.dataset.igvViewKey = viewKey;
  const open = sameRun ? Array.from(content.querySelectorAll("details[open]")).map(node => node.dataset.detail) : [];
  if (run.paired_review) renderPairedTrackReview(run, content, open);
  else renderSingleTrackReview(run, content, open);
}

function renderPreset(run, content) {
  const open = content.dataset.runId === run.id && content.querySelector('details[open]');
  const apiRun = `/api/runs/${encodeURIComponent(run.id)}`;
  const decision = run.answer?.decision;
  const pair = run.status === "completed_comparison", blocked = run.status === "blocked_no_matched_igg";
  const downloads = {"coverage.cpm.bw": "CPM BigWig", "coverage.spikein.bw": "Spike-in BigWig", "tool.log": "计算日志"};
  content.dataset.runId = run.id;
  content.innerHTML = `<div class="run-topline"><h3>${escape(run.name)}</h3><p>${escape(run.sample_id)} · ${stamp(run.created_at)} · ${run.source === "web" ? "网页发起" : "本地终端运行"}</p></div>
    ${activeRun(run) ? `<p class="run-disclaimer" role="status">${escape(run.message)}</p>` : ""}
    ${decision ? `<section class="model-trace"><h3>模型建议</h3><p><strong>${decision.action === "generate_pair" ? "用相同输入生成两份轨道" : "保留原候选峰"}</strong></p><p>${escape(decision.reason)}</p></section>` : ""}
    ${pair || blocked ? `<section class="reused-results"><h3>实际执行</h3>
      <p>${pair ? `复用同一份历史去重 BAM，生成两份 hg38 全基因组 BigWig，窗口均为 ${escape(run.bin_size)} bp。用时约 ${Math.round(run.execution_seconds)} 秒。` : "已确认无匹配 IgG，程序停止找峰并保留原窄峰和宽峰 BED。"}</p>
      <p class="field-help">模型请求 ${count(run.model_requests)} 次 · 轨道计算 ${count(run.coverage_executions)} 次 · 找峰 ${count(run.peak_executions)} 次 · FASTQ 比对 ${count(run.fastq_executions)} 次</p>
      <ul>${pair ? run.artifacts.map(a => `<li><a href="${apiRun}/artifacts/${escape(a.name)}" download>${downloads[a.name]}</a></li>`).join("") :
        run.preserved_peaks.map((a) => `<li><a href="/api/runs/${encodeURIComponent(run.baseline_id)}/artifacts/${escape(a.name)}" download>${a.name === "narrow.peaks.bed" ? "原窄峰 BED" : "原宽峰 BED"}</a></li>`).join("")}</ul></section>` : ""}
    ${run.error ? `<p class="inline-error" role="alert">${escape(run.error)}</p>` : ""}
    ${state.runError ? `<p class="inline-error" role="alert">${escape(state.runError)}</p>` : ""}
    ${run.sent_summary ? `<details class="log-preview" ${open ? "open" : ""}><summary>查看发给模型的摘要</summary><pre>${escape(JSON.stringify(run.sent_summary, null, 2))}</pre></details>` : ""}
    <div class="run-actions">${activeRun(run) ? `<button class="button danger" data-cancel="${escape(run.id)}" ${state.stopping.has(run.id) ? "disabled" : ""}>停止预设</button>` : ""}
      <a class="button secondary" href="${apiRun}/report" download>下载运行记录</a></div>`;
}

function renderRun() {
  const inHistory = !!state.historyRunId;
  const sampleId = currentSample()?.id;
  const recent = state.runs.filter(item => item.sample_id === sampleId && listedRun(item));
  const run = inHistory ? state.historyRun : state.currentRun?.sample_id === sampleId ? state.currentRun :
    recent.find(activeRun) || recent.find(item => item.mode === "review") || recent[0];
  const content = $(inHistory ? "#history-run-content" : "#run-content");
  const label = $(inHistory ? "#history-run-label" : "#current-run-label");
  renderStartButton();
  if (inHistory) $("#history-run-status").textContent = state.historyError || (run ? "" : "正在读取运行详情…");
  else $("#run-empty").hidden = !!run;
  content.hidden = !run;
  label.textContent = run ? runLabel(run) : "";
  label.className = "status " + (run ? runClass(run) : "");
  if (!run) return;
  if (run.mode === "review") return renderTrackReview(run, content);
  if (run.mode === "pipeline") return renderPipeline(run, content);
  if (run.mode === "agent") return renderAgent(run, content);
  if (run.mode === "spikein") return renderSpikein(run, content);
  if (run.mode === "preset") return renderPreset(run, content);
  const sameRun = content.dataset.runId === run.id;
  const open = sameRun ? Array.from(content.querySelectorAll("details[open]")).map((node) => node.dataset.detail) : [];
  const focused = content.contains(document.activeElement) ? document.activeElement.closest("details")?.dataset.detail : null;
  const real = run.mode === "coverage";
  const stepLabels = { pending: "未开始", running: "进行中", done: "已完成", simulated: "模拟", cancelled: "已停止", failed: "失败" };
  const steps = run.steps.map((step, index) => `<li class="timeline-step ${escape(step.status)}"><span class="step-number">${index + 1}</span><div><strong>${escape(step.label)}</strong><span class="step-state">${escape(stepLabels[step.status] || step.status)}</span><p>${escape(step.detail)}</p></div></li>`).join("");
  content.dataset.runId = run.id;
  content.innerHTML = `<div class="run-topline"><h3>${escape(run.name)}</h3><p>${escape(run.sample_id)} · ${stamp(run.created_at)}</p></div>
    ${activeRun(run) ? `<progress class="run-progress" max="100" value="${run.progress}" aria-label="任务步骤进度">${run.progress}%</progress>` : ""}
    ${real ? `<dl class="run-summary"><div><dt>真实区域试跑</dt><dd>chr1:1,000,000–2,000,000</dd></div><div><dt>计算参数</dt><dd>窗口 25 bp · 未标准化 · 仅计 R1 读段</dd></div></dl>` : `<dl class="run-summary"><div><dt>历史报告</dt><dd>人源片段 <strong>${count(run.baseline.metrics.human_fragments)}</strong></dd></div><div><dt>演示计划</dt><dd>轨道窗口 <strong>${run.plan.bin_size_before} → ${run.plan.bin_size_after} bp</strong></dd></div></dl>`}
    <p class="run-disclaimer">${escape(run.message)}</p>
    ${run.error ? `<p class="inline-error" role="alert">${escape(run.error)}</p>` : ""}
    ${state.runError ? `<p class="inline-error" role="alert">${escape(state.runError)}</p>` : ""}
    ${state.pollError ? `<p class="inline-error" role="alert">${escape(state.pollError)}</p>` : ""}
    ${renderModelTrace(run, open)}
    ${real && run.status === "completed" && run.result ? coverageChart(run.result) : ""}
    ${run.artifacts?.some(f => f.name === "tool.log") ? logPreview(`/api/runs/${encodeURIComponent(run.id)}/artifacts/tool.log`, "coverage-log", open) : ""}
    <div class="run-actions">${activeRun(run) ? `<button class="button danger" data-cancel="${escape(run.id)}" ${state.stopping.has(run.id) ? "disabled" : ""}>${state.stopping.has(run.id) ? "正在停止…" : real ? "停止计算" : "停止演示"}</button>` : ""}${real && run.status === "completed" ? `<a class="button secondary" href="/api/runs/${encodeURIComponent(run.id)}/artifacts/coverage.bw" download>下载 BW</a>` : ""}<a class="button secondary" href="/api/runs/${encodeURIComponent(run.id)}/report" download>下载运行记录</a>${run.artifacts?.some((file) => file.name === "tool.log") ? `<a class="button secondary" href="/api/runs/${encodeURIComponent(run.id)}/artifacts/tool.log" download>下载日志</a>` : ""}</div>
    <details data-detail="steps" ${open.includes("steps") ? "open" : ""}><summary>运行步骤</summary><ol class="timeline">${steps}</ol></details>
    <details data-detail="events" ${open.includes("events") ? "open" : ""}><summary>事件记录</summary><ul class="event-list">${run.events.map((event) => `<li><time>${stamp(event.at, true)}</time><span>${escape(event.message)}</span></li>`).join("")}</ul></details>`;
  if (focused) content.querySelector(`details[data-detail="${focused}"] summary`)?.focus({ preventScroll: true });
}

function renderHistory() {
  const runs = state.runs.filter(run => state.archivedRunIds.has(run.id) === state.showArchived);
  const filters = `<div class="history-filters" aria-label="切换运行记录视图">
    <button class="button secondary" data-history-view="current" aria-pressed="${!state.showArchived}">纵览</button>
    <button class="button secondary" data-history-view="archived" aria-pressed="${state.showArchived}">已删除</button>
    ${state.showArchived ? '<p class="field-help">已删除的记录可移回纵览，产物仍可查看。</p>' : ""}
  </div>`;
  if (!runs.length) {
    $("#history-content").innerHTML = filters + `<div class="empty-state"><h2>${state.showArchived ? "暂无已删除记录" : "暂无运行记录"}</h2></div>`;
    return;
  }
  const type = run => run.paired_adaptive ? "配对样本 Agent" : run.mode === "pipeline" ? (run.model_review_requested ? "FASTQ 分析与模型复核" : "基础分析 · " + (run.source === "fastq" ? "FASTQ" : "BAM")) :
    ({ review: "模型复核", spikein: "Spike-in 轨道", preset: "模型预设", agent: "旧区域联通测试", coverage: "区域工具试跑", demo: "流程演示" }[run.mode] || "分析记录");
  $("#history-content").innerHTML = filters + `<div class="table-wrap"><table class="history-table"><thead><tr><th scope="col">任务 / 样本</th><th scope="col">类型</th><th scope="col">创建时间</th><th scope="col">状态</th><th scope="col">操作</th></tr></thead><tbody>${runs.map((run) => `<tr><td><span class="history-name">${escape(run.name)}</span><span class="table-secondary">${escape(run.sample_id)}</span></td><td>${type(run)}</td><td>${stamp(run.created_at)}</td><td><span class="status ${runClass(run)}">${runLabel(run)}</span></td><td><div class="history-row-actions"><button class="table-action" data-run="${escape(run.id)}" aria-label="查看记录：${escape(run.name)}">${run.status === "awaiting_rerun" ? "查看建议" : "查看记录"}</button><button class="table-action${state.showArchived ? "" : " history-delete"}" data-archive-run="${escape(run.id)}" data-archive="${!state.showArchived}" aria-label="${state.showArchived ? "恢复" : "删除"}${escape(run.name)}" ${state.archiving.has(run.id) ? "disabled" : ""}>${state.showArchived ? "恢复" : "删除"}</button></div></td></tr>`).join("")}</tbody></table></div>`;
}

function renderSettings() {
  const setup = state.setup;
  $("#settings-content").innerHTML = `<section class="settings-section"><h2>基础分析</h2><div><dl class="settings-rows"><div><dt>工具流程</dt><dd>${setup.pipeline_runtime?.ready ? "FASTQ 清理、比对、BAM、BW、Peak、QC 图表已接入" : "尚未就绪"}</dd></div><div><dt>参考</dt><dd>当前本地匹配的 hg38 / GENCODE v44</dd></div><div><dt>标准化</dt><dd>CPM；IgG 和 spike-in 策略待补实验条件</dd></div></dl></div></section><section class="settings-section"><h2>分析工具</h2><div><dl class="settings-rows"><div><dt>bamCoverage</dt><dd>${setup.tools_connected ? "可用于区域计算" : "尚未就绪"}</dd></div><div><dt>SAMtools</dt><dd>${setup.coverage_runtime?.samtools_ready ? "已接入基础检查与统计" : "尚未就绪"}</dd></div><div><dt>输入</dt><dd>本地 BAM + BAI 索引</dd></div></dl></div></section>
    <section class="settings-section"><h2>模型服务</h2><div><dl class="settings-rows"><div><dt>连接测试</dt><dd>${setup.model_connected ? "最近一次合成测试通过" : "尚未通过"}</dd></div><div><dt>模型</dt><dd>${escape(setup.model_name || "尚未配置")}</dd></div><div><dt>API Key</dt><dd>${setup.model_key_configured ? "已填写" : "尚未填写"}</dd></div><div><dt>当前样本完整质控</dt><dd>${modelAuthorized() ? "外部使用范围已确认" : "可在本地查看，外部使用范围待确认"}</dd></div></dl><button class="button secondary" data-test-model ${!setup.model_key_configured || state.testingModel ? "disabled" : ""}>${state.testingModel ? "正在测试…" : "测试模型（合成内容）"}</button><p id="model-test-result" class="field-help" role="status">${escape(state.modelTestMessage || "")}</p><details><summary>配置说明</summary><p>现在不用填写 Key。本机私有 <code>.env</code> 已配置；网页不接收或显示 Key。</p><p>完整质控摘要只包含结构化统计和实验背景状态。可在复核记录查看拟发送的摘要。</p></details></div></section>`;
}

async function refreshRuns() {
  const data = await api("/api/runs");
  state.runs = data.runs;
  state.archivedRunIds = new Set(data.archived_run_ids);
  state.pollError = "";
  state.runs.filter((run) => !activeRun(run)).forEach((run) => state.stopping.delete(run.id));
  if (state.currentRun) {
    const item = state.runs.find((run) => run.id === state.currentRun.id);
    if (item) state.currentRun = { ...state.currentRun, ...item };
  }
  const historyRun = state.runs.find(run => run.id === state.historyRunId);
  if (historyRun) state.historyRun = { ...state.historyRun, ...historyRun };
  renderMode();
  renderExistingResults();
  renderRun();
  renderHistory();
  if (location.hash === "#samples") renderSamples();
  schedulePoll();
}

function schedulePoll() {
  clearTimeout(state.poll);
  if (state.runs.some(activeRun)) state.poll = setTimeout(() => refreshRuns().catch((error) => { state.pollError = "运行状态暂未更新：" + error.message; renderRun(); schedulePoll(); }), 900);
}

async function load() {
  try {
    const [samples, setup, runs] = await Promise.all([api("/api/samples"), api("/api/health"), api("/api/runs")]);
    $(".app-error")?.remove();
    state.samples = samples.samples;
    state.setup = setup;
    state.runs = runs.runs;
    state.archivedRunIds = new Set(runs.archived_run_ids);
    if (!state.samples.some((sample) => sample.id === state.selected)) state.selected = state.samples[0]?.id || null;
    renderMode();
    state.currentRun = state.runs.find((run) => run.id === state.currentRun?.id) || state.runs.find(listedRun) || null;
    if (state.currentRun?.mode === "pipeline") state.currentRun = await api(`/api/runs/${encodeURIComponent(state.currentRun.id)}`);
    renderSampleSelection(); renderExistingResults(); renderSamples(); renderHistory(); renderSettings(); renderRun();
    showView(location.hash.slice(1));
    schedulePoll();
  } catch (error) {
    if (!$(".app-error")) $("#main-content").insertAdjacentHTML("afterbegin", `<div class="app-error" role="alert"><strong>本地服务暂时没有响应</strong><p>${escape(error.message)}</p><button class="button secondary" data-retry>重新读取</button></div>`);
    $("#start-demo").disabled = true;
  }
}

function selectSample(id) {
  const previous = currentSample()?.id;
  state.selected = id;
  if (id !== previous) $("#analysis-options").open = false;
  $("#form-error").textContent = "";
  renderSampleSelection(); renderMode(); renderExistingResults(); renderSamples(); renderRun();
}

function openPairTab(tab) {
  state.sampleTab = tab;
  if (location.hash !== "#samples") history.pushState(null, "", "#samples");
  showView("samples");
  setSampleTab(tab, true);
}

$("#sample-select").addEventListener("change", (event) => selectSample(event.target.value));
$("#run-mode").addEventListener("change", renderMode);
$("#analysis-options").addEventListener("toggle", renderMode);
document.addEventListener("change", (event) => {
  if (event.target.id !== "library-baseline") return;
  state.libraryBaselines.set(currentSample().id, event.target.value);
  renderSamples();
  $("#library-baseline")?.focus();
});

$("#task-form").addEventListener("submit", async (event) => {
  event.preventDefault();
  if ($("#start-demo").disabled || state.busy || state.runs.some(activeRun)) return;
  state.busy = true;
  $("#form-error").textContent = "";
  renderStartButton();
  try {
    const mode = $("#run-mode").value;
    const name = $("#task-name");
    const payload = { sample_id: $("#sample-select").value, name: name.disabled ? taskNames[mode] : name.value, mode };
    if (mode === "pipeline") Object.assign(payload, { source: "fastq", model_review: true });
    else if (mode === "review") {
      if (!currentSample()?.pair) payload.baseline_id = cpmBaseline()?.id;
    }
    else if (mode !== "paired_agent") Object.assign(payload, { mode: "preset", preset: mode, review_id: presetReview(mode)?.id });
    const run = await api("/api/runs", { method: "POST", body: JSON.stringify(payload) });
    state.currentRun = run;
    state.runs = [run, ...state.runs.filter((item) => item.id !== run.id)];
    state.runError = "";
    state.pollError = "";
    renderRun();
    renderHistory();
    schedulePoll();
    toast(mode === "pipeline" ? "已创建分析任务" : mode === "paired_agent" || (mode === "review" && modelReady()) ? "已启动 Agent" :
      mode === "review" ? "已创建质控复核任务" : "已创建预设任务");
    scrollToRun();
  } catch (error) {
    $("#form-error").textContent = error.message;
    schedulePoll();
  } finally {
    state.busy = false;
    renderStartButton();
  }
});

document.addEventListener("click", async (event) => {
  const navigation = event.target.closest("[data-view], [data-go]");
  if (navigation) return go(navigation.dataset.view || navigation.dataset.go);
  if (event.target.closest("#history-back")) {
    if (history.state?.cuttagRunParent) history.back();
    else go("history");
    return;
  }
  const tab = event.target.closest("[data-sample-tab]");
  if (tab) return setSampleTab(tab.dataset.sampleTab);
  const historyFilter = event.target.closest("[data-history-view]");
  if (historyFilter) {
    state.showArchived = historyFilter.dataset.historyView === "archived";
    renderHistory();
    $("#history-content [aria-pressed='true']")?.focus();
    return;
  }
  const archiveButton = event.target.closest("[data-archive-run]");
  if (archiveButton) {
    const id = archiveButton.dataset.archiveRun;
    const archived = archiveButton.dataset.archive === "true";
    if (state.archiving.has(id)) return;
    state.archiving.add(id);
    archiveButton.disabled = true;
    try {
      const result = await api(`/api/runs/${encodeURIComponent(id)}/archive`, {
        method: "POST", body: JSON.stringify({archived})
      });
      if (result.archived_run_ids.includes(id)) state.archivedRunIds.add(id);
      else state.archivedRunIds.delete(id);
      if (archived && state.currentRun?.id === id && !activeRun(state.currentRun)) {
        state.currentRun = state.runs.find(listedRun) || null;
      }
      renderExistingResults();
      if (location.hash === "#samples") renderSamples();
      toast(archived ? "已从纵览删除，可在「已删除」恢复。" : "已恢复到纵览。");
    } catch (error) {
      toast(error.message);
    } finally {
      state.archiving.delete(id);
      renderHistory();
      if (location.hash === "#history") $("#history-content [aria-pressed='true']")?.focus({preventScroll:true});
    }
    return;
  }
  const pairTab = event.target.closest("[data-pair-tab]");
  if (pairTab) {
    if (pairTab.dataset.sampleId) selectSample(pairTab.dataset.sampleId);
    return openPairTab(pairTab.dataset.pairTab);
  }
  const tuning = event.target.closest("[data-open-pair-agent]");
  if (tuning) {
    selectSample(tuning.dataset.openPairAgent);
    $("#run-mode").value = "paired_agent";
    renderMode();
    return go("workspace");
  }
  const choice = event.target.closest("[data-decide-rerun], [data-replay-rerun]");
  if (choice) {
    const replay = !!choice.dataset.replayRerun;
    const id = choice.dataset.replayRerun || choice.dataset.decideRerun;
    if (state.deciding.has(id)) return;
    state.deciding.add(id);
    state.runError = "";
    let showProgress = false;
    renderRun();
    try {
      const action = replay ? "replay-paired-rerun" : "decide-paired-rerun";
      const updated = await api(`/api/runs/${encodeURIComponent(id)}/${action}`,
        {method:"POST", body:JSON.stringify(replay ? {} : {accept:choice.dataset.accept === "true"})});
      if (replay) {
        state.runs = [updated, ...state.runs.filter(run => run.id !== updated.id)];
        state.currentRun = updated;
        go(`history/${updated.id}`);
      } else {
        showProgress = updated.paired_rerun_accepted && activeRun(updated);
        state.runs = state.runs.map(run => run.id === id ? updated : run);
        if (state.currentRun?.id === id) state.currentRun = updated;
        if (state.historyRunId === id) state.historyRun = updated;
      }
      renderHistory();
      if (location.hash === "#samples") renderSamples();
      schedulePoll();
    } catch (error) {
      state.runError = error.message;
      toast(error.message);
    } finally {
      state.deciding.delete(id);
      renderExistingResults();
      renderRun();
      if (showProgress) {
        const panel = state.historyRunId === id ? $("#history-run-content") :
          !$("#view-workspace").hidden && $("#run-content").dataset.runId === id ? $("#run-content") : null;
        panel?.querySelector(".phase-progress")?.scrollIntoView({block:"nearest"});
      }
    }
    return;
  }
  const librarySample = event.target.closest("[data-library-sample]");
  if (librarySample) {
    state.sampleTab = "results";
    selectSample(librarySample.dataset.librarySample);
    $('.library-index [aria-current="true"]')?.focus({preventScroll:true});
    return;
  }
  const runButton = event.target.closest("[data-run]");
  if (runButton) {
    const id = runButton.dataset.run;
    if (runButton.closest("#view-history")) {
      history.pushState({cuttagRunParent: true}, "", "#history/" + id);
      showView(location.hash.slice(1));
      return;
    }
    state.currentRun = state.runs.find((run) => run.id === id);
    if (state.currentRun && currentSample()?.id !== state.currentRun.sample_id) selectSample(state.currentRun.sample_id);
    state.runError = "";
    renderRun(); go("workspace");
    // Hash navigation reveals the target before scrolling on narrow screens.
    requestAnimationFrame(scrollToRun);
    if (state.currentRun?.mode === "pipeline") {
      try {
        const detail = await api(`/api/runs/${encodeURIComponent(id)}`);
        if (state.currentRun?.id === id) { state.currentRun = detail; renderRun(); }
      } catch (error) {
        if (state.currentRun?.id === id) { state.pollError = error.message; renderRun(); }
      }
    }
    return;
  }
  const cancel = event.target.closest("[data-cancel]");
  if (cancel) {
    const id = cancel.dataset.cancel;
    state.stopping.add(id);
    state.runError = "";
    renderRun();
    try {
      const updated = await api(`/api/runs/${encodeURIComponent(id)}/cancel`, { method: "POST", body: "{}" });
      state.runs = state.runs.map((run) => run.id === id ? updated : run);
      if (state.currentRun?.id === id) state.currentRun = updated;
      if (state.historyRunId === id) state.historyRun = updated;
      if (!activeRun(updated)) state.stopping.delete(id);
      toast("已请求停止任务"); schedulePoll();
    }
    catch (error) { state.runError = error.message; state.stopping.delete(id); }
    finally { renderRun(); renderHistory(); }
  }
  if (event.target.closest("[data-retry]")) load();
  if (event.target.closest("[data-test-model]")) {
    state.testingModel = true; state.modelTestMessage = "正在请求模型，使用固定合成内容…"; renderSettings();
    try {
      const result = await api("/api/model/test", { method: "POST", body: "{}" });
      state.modelTestMessage = result.ok ? "合成测试通过，建议格式有效。" : result.error || "测试未通过";
      state.setup = await api("/api/health");
      renderMode();
    } catch (error) { state.modelTestMessage = error.message; }
    finally { state.testingModel = false; renderSettings(); }
  }
});

document.addEventListener("keydown", (event) => {
  if (!event.target.matches("[data-sample-tab]") || !["ArrowLeft", "ArrowRight", "Home", "End"].includes(event.key)) return;
  event.preventDefault();
  const index = sampleTabs.indexOf(event.target.dataset.sampleTab);
  const tab = event.key === "Home" ? sampleTabs[0] : event.key === "End" ? sampleTabs.at(-1) :
    sampleTabs[(index + (event.key === "ArrowLeft" ? -1 : 1) + sampleTabs.length) % sampleTabs.length];
  setSampleTab(tab, true);
});

window.addEventListener("hashchange", () => showView(location.hash.slice(1)));
document.querySelectorAll("[data-icon]").forEach((node) => { node.innerHTML = icon(node.dataset.icon); });
load();
