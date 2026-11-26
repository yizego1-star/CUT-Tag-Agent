# CUT&Tag Agent

用本地网页运行 CUT&Tag 基础分析、查看质控和 IGV 轨道，并用可选的模型复核结果。

- 从双端 FASTQ 或已有 BAM 开始，用 fastp、Bowtie 2、SAMtools、deepTools 和 MACS3 生成质控、CPM 轨道和候选峰。
- 比较目标与匹配 IgG。模型可以根据结构化质控摘要建议调整 MACS3 的 q 值。人工确认后最多重跑一次，保留两轮记录。
- 在浏览器里查看 BigWig、BED 与重跑前后的对照。科学判读需要结合实验设计和人工复核。

## 本地启动

需要 Python 3.9 或更高版本。网页无需构建，实际分析另需 Docker 和自己的本地输入。

```bash
python3 app.py
```

打开 [http://127.0.0.1:8766](http://127.0.0.1:8766)。首次启动时样本库为空。添加样本、构建工具镜像和配置模型的步骤见[本地使用说明](docs/local-setup.md)。

## 代码与验证

| 路径 | 内容 |
| --- | --- |
| `app.py`、`cuttag_web/` | 本地 API、任务编排、质控与模型复核 |
| `web/dist/` | 原生网页和随附的 IGV.js |
| `tools/` | 容器配方、分析脚本和本地工具入口 |
| `tests/` | 合成输入测试与浏览器检查 |

```bash
python3 -m unittest discover -s tests -v
```

仓库提供项目源码。测序数据、参考文件、运行结果、授权记录和凭据均由使用者保存在本地。随附的 IGV.js 使用 [MIT 许可证](web/dist/vendor/igv-3.8.9-LICENSE.txt)。
