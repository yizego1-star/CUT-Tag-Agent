# 本地使用

## 基础分析

在项目根目录构建工具镜像。

```bash
docker build -t cuttag-deeptools:3.5.6 tools/deeptools
docker build -t cuttag-pipeline:2026-11-17 tools/pipeline
docker build -f tools/pipeline/Dockerfile.resume -t cuttag-pipeline-resume:2026-11-18 tools/pipeline
```

将获准使用的样本放入 `data/private/<样本ID>/`。文件名中的 `<样本ID>` 与目录名一致。

| 输入 | 本地相对路径 |
| --- | --- |
| 双端 FASTQ | `raw/<样本ID>_1.fq.gz`、`raw/<样本ID>_2.fq.gz` |
| 已过滤和去重的 BAM | `processed/<样本ID>.hg38.dedup.bam` 及同名 `.bai` |
| 可选质控报告 | `qc/fastp.json`、`qc/summary.tsv` |
| 可选文件清单和核对记录 | `provenance/source_manifest.json`、`provenance/verification.json` |

样本配套参考放入 `data/support/<样本ID>/`。流程需要匹配的 `reference/hg38/hg38.bam-compatible.chrom.sizes`、`annotation/gencode.v44.genes.bam-compatible.bed` 和完成核对后的 `provenance/reference_compatibility.json`。FASTQ 入口另需完整的 `reference/hg38/bowtie2/hg38_noalt.*.bt2`。以同一参考组装核对 BAM、索引和注释。网页支持从 FASTQ 开始。已有 BAM 可通过本地 API 发起基础分析。

用 `python3 app.py` 启动网页。独立本地任务记录保存在 `outputs/web_runs/`，分析结果保存在 `outputs/`。

例如从本地样本 `my-sample` 的已有 BAM 开始，不请求模型。

```bash
curl -fsS http://127.0.0.1:8766/api/runs \
  -H 'Content-Type: application/json' \
  -d '{"sample_id":"my-sample","mode":"pipeline","source":"bam","model_review":false}'
```

`tools.resume_pipeline` 用于从已停止的 FASTQ 任务的可验证 BAM 检查点续跑。独立 Spike-in 轨道入口见 `python3 -m tools.spikein_track --help`。

## 模型与配对复核

可将 `.env.example` 复制为本机私有 `.env`，填写自己的 HTTPS 服务地址、模型名称和 API Key。接入设置里的连通测试使用合成摘要。使用真实摘要前需确认服务和数据的适用授权，再打开 `CUTTAG_ALLOW_MINIMAL_SUMMARY`。完整质控复核另外核对 `.local/` 中的样本、服务地址、模型和字段组。字段定义见 `cuttag_web/qc_evidence.py` 与 `cuttag_web/hospital_pair_review.py`。

配对入口使用两个本地固定角色 ID，`paired-target` 和 `paired-igg`。将通过核对的输入和清单放在 `data/support/paired-case/`，共享参考放在 `data/support/shared-reference/`。运行 `python3 -m tools.prepare_hospital_pair inventory` 先检查登记，使用 `pilot` 或 `full` 生成候选峰与轨道。当前配对页面还需要本地登记的 `spikein-context.json` 和两份历史 Spike-in 轨道。输入清单、角色校验与产物约定见下列实现。

- [清单和计算入口](../tools/prepare_hospital_pair.py)
- [配对网页证据](../cuttag_web/hospital_pair.py)

## 公开 H9 开发流程

`tools/run_h9.py` 保留了从 FASTQ 到候选峰和 IgG 对照比较的专用流程。参考与 FASTQ 需另行准备。公开数据来源为 [GSE200768](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE200768)。

准备该流程所需的下载与校验清单后，可运行 `python3 -m tools.run_h9 --pairs-per-sample 100000 --compare-tools` 做限量工程检查。该入口需要上述 Docker 镜像、共享 hg38 参考及 `data/public/GSE200768_H9_H3K4me1/` 下的本地 FASTQ 和清单。输入检查及文件布局见 [H9 入口](../tools/run_h9.py)。

## 浏览器检查

`tests/browser_smoke.cjs` 检查空样本库时的导航、任务状态与响应式布局。先启动空白项目，再准备 Playwright 和 Chromium。

```bash
npm install --prefix .local/web-qa playwright
.local/web-qa/node_modules/.bin/playwright install chromium
node tests/browser_smoke.cjs
```

可用 `CUTTAG_TEST_URL` 指定测试服务，用 `CUTTAG_BROWSER_PATH` 选择已有 Chromium。其他浏览器脚本用于特定工具或既有运行记录的复查，其中早期区域试跑脚本依赖历史界面入口。
