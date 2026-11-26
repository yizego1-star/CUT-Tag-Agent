import argparse
import copy
import csv
import json
import re
import shutil
import subprocess
import threading
import time
import uuid
from datetime import datetime, timedelta, timezone
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import unquote, urlsplit
from xml.etree.ElementTree import Element, SubElement, tostring
from tools.pair_comparison_figure import FIGURE_NAME, LOCUS_NAME, OUTPUT_DIR as FIGURE_DIR, build_comparison_figure, comparison_figure_info
from tools.prepare_hospital_pair import REFERENCE, REF_FILES, rerun_dir, rerun_peaks, stop_container
from tools.run_cuttag_preset import review_source, run_preset
from tools.spikein_track import cancel_tracks

from .coverage import ARTIFACTS, DockerCoverageRunner
from . import hospital_pair_review as pair_review
from .hospital_pair import CONTROL as PAIR_CONTROL, TARGET as PAIR_TARGET, HospitalPair
from .model import ModelClient, validate_qc_payload
from .agent import drive_agent
from .pipeline import PipelineRunner, MIMES as PIPELINE_ARTIFACTS, STAGES as PIPELINE_STAGES, sha256
from .track_review import compare_evidence, make_qc_evidence, review_adaptive_pair, review_baseline, review_evidence

CONTAINER_MODES = ("coverage", "agent", "pipeline")
REAL_MODES = (*CONTAINER_MODES, "review", "preset")
SPIKEIN_ARTIFACTS = {"coverage.spikein.bw": "application/octet-stream", "report.md": "text/plain; charset=utf-8", "tool.log": "text/plain; charset=utf-8"}
PRESET_NAMES = {"normalization": "预设 · CPM 与 Spike-in 标准化", "igg": "预设 · 缺匹配 IgG 时停止找峰"}
PRESET_ARTIFACTS = {key: PIPELINE_ARTIFACTS[key] for key in ("coverage.cpm.bw", "tool.log")}
PRESET_ARTIFACTS["coverage.spikein.bw"] = SPIKEIN_ARTIFACTS["coverage.spikein.bw"]
PAIR_RERUN_ARTIFACTS = {"tuned.peaks.bed": PIPELINE_ARTIFACTS["narrow.peaks.bed"],
                        "tuned.json": "application/json", "commands.json": "application/json",
                        "tool.log": PIPELINE_ARTIFACTS["tool.log"]}
PIPELINE_LOGS = ("tool.log", "bowtie2.log", "commands.json")


PROJECT_ROOT = Path(__file__).resolve().parent.parent
STEPS = [
    ("inputs", "检查输入", "确认本地样本和报告可读"),
    ("report", "读取历史报告", "整理已有质控指标"),
    ("plan", "生成演示计划", "使用固定示例，不调用模型"),
    ("rerun", "演示重跑步骤", "真实工具待接入，本次不重新计算"),
    ("finish", "保存运行记录", "保留输入、参数和步骤状态"),
]


def now():
    return datetime.now(timezone.utc).isoformat()


def read_json(path, default=None):
    try:
        if path.stat().st_size > 2 * 1024 * 1024:
            return default
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return default


def number(value):
    try:
        return float(value) if "." in str(value) else int(value)
    except (ValueError, TypeError):
        return None


def first_peak_locus(path):
    with path.open(encoding="utf-8") as stream:
        line = stream.readline()
    if not line:
        return None
    chrom, start, end = line.split()
    flank = 800
    return f"{chrom}:{max(1, int(start) + 1 - flank)}-{int(end) + flank}"


class Workbench:
    def __init__(self, project_root=PROJECT_ROOT, step_delay=0.65, coverage_runner=None, pipeline_runner=None):
        self.root = Path(project_root).resolve()
        self.data_root = self.root / "data" / "private"
        self.output_root = self.root / "outputs" / "web_runs"
        self.output_root.mkdir(parents=True, exist_ok=True)
        self.step_delay = step_delay
        self.lock = threading.RLock()
        self.runs = {}
        self.local_runs = {}
        self.cancelled = set()
        self.coverage = coverage_runner or DockerCoverageRunner(self.root)
        self.model = ModelClient(self.root)
        self.pipeline = pipeline_runner or PipelineRunner(self.root)
        self.hospital_pair = HospitalPair(self.root)
        for path in self.output_root.glob("run-*.json"):
            item = read_json(path)
            if isinstance(item, dict) and isinstance(item.get("id"), str):
                if item.get("status") in ("queued", "running"):
                    if item.get("mode") in CONTAINER_MODES and not item.get("baseline_ready"):
                        (self.pipeline if item["mode"] == "pipeline" else self.coverage).cancel(item["id"])
                    if item.get("mode") == "preset":
                        self.cancel_preset_container(item)
                    if item.get("paired_adaptive") and (item.get("evidence_preparing") or item.get("track_review", {}).get("rerun_started")):
                        stop_container('cuttag-pair-' + item['id'][4:])
                    if item.get("baseline_ready"):
                        item["status"] = "completed_pipeline"
                        item["message"] = "基础分析产物已保留；模型复核因本地服务重启而中断。"
                    else:
                        item["status"] = "interrupted"
                        item["message"] = "本地服务已重启，之前的任务未继续。"
                    if item.get("track_review", {}).get("status") == "reviewing":
                        item["track_review"].update(status="interrupted", error="本地服务已重启，模型复核未继续。")
                    self.save(item)
                if item.get("mode") == "pipeline" and item["status"] == "completed_pipeline" and "report_metrics" not in item:
                    item["report_metrics"] = self.pipeline.report_metrics(self.pipeline.output_dir(item["id"]), item.get("source"))
                    self.save(item)
                self.runs[item["id"]] = item

    def sample_path(self, sample_id):
        if not re.fullmatch(r"[A-Za-z0-9_.-]{1,120}", sample_id):
            raise ValueError("样本名称无效")
        path = (self.data_root / sample_id).resolve()
        if path.parent != self.data_root.resolve() or not path.is_dir():
            raise ValueError("没有找到这个本地样本")
        return path

    def sample(self, sample_id):
        path = self.sample_path(sample_id)
        manifest = read_json(path / "provenance" / "source_manifest.json", {})
        verification = read_json(path / "provenance" / "verification.json", {})
        fastp = read_json(path / "qc" / "fastp.json", {})
        summary = fastp.get("summary", {})
        before = summary.get("before_filtering", {})
        after = summary.get("after_filtering", {})
        metrics = {
            "total_reads": before.get("total_reads"),
            "retained_reads": after.get("total_reads"),
            "q30_rate": before.get("q30_rate"),
            "gc_content": before.get("gc_content"),
            "read_length": before.get("read1_mean_length"),
        }
        table = path / "qc" / "summary.tsv"
        if table.is_file():
            with table.open(encoding="utf-8", newline="") as stream:
                for row in csv.DictReader(stream, delimiter="\t"):
                    if row.get("sample") == sample_id:
                        for field in ("human_fragments", "spikein_fragments", "duplicate_fraction_filtered", "scale_factor"):
                            metrics[field] = number(row.get(field))
                        metrics["historical_bw_status"] = row.get("bw_status")
                        break
        files = []
        for item in manifest.get("files", []):
            relative = item.get("local_relative_path", "")
            file_path = (path / relative).resolve()
            if path not in file_path.parents:
                continue
            size = file_path.stat().st_size if file_path.is_file() else 0
            files.append({
                "name": Path(relative).name, "group": relative.split("/")[0],
                "bytes": size, "present": file_path.is_file(),
                "size_matches": size == item.get("expected_bytes"),
            })
        bins = [0] * 100
        lengths = path / "qc" / "hg38.dedup.fragment_length_counts.txt"
        if lengths.is_file():
            for line in lengths.read_text(encoding="utf-8").splitlines():
                fields = line.split()
                if len(fields) < 2:
                    continue
                try:
                    count, length = int(fields[0]), int(fields[1])
                    if count >= 0 and 0 <= length < 1000:
                        bins[length // 10] += count
                except ValueError:
                    continue
        return {
            "id": sample_id, "label": sample_id,
            "input_kind": "双端 FASTQ / BAM + 质控报告", "metrics": metrics, "files": files,
            "total_bytes": sum(item["bytes"] for item in files),
            "verification_status": verification.get("status", "not_checked"),
            "verification_checks": len(verification.get("checks", [])),
            "report_available": bool(before),
            "fragment_distribution": {"bin_size": 10, "counts": bins},
            "context_status": "实验分组和对照说明待补充",
            "source_note": "本地已有报告；不是本次重新分析的结果",
        }

    def samples(self):
        result = []
        if self.data_root.is_dir():
            for path in sorted(self.data_root.iterdir()):
                if path.is_dir():
                    try:
                        result.append(self.sample(path.name))
                    except (ValueError, OSError):
                        continue
        if self.hospital_pair.folder.is_dir():
            try:
                result = [item for item in result if item["id"] != PAIR_TARGET] + [self.paired_sample()]
            except (ValueError, OSError, KeyError, TypeError):
                pass
        return result

    def paired_sample(self):
        sample = self.hospital_pair.sample()
        sample["pair"]["external_model_allowed"] = self.model.pair_review_allowed()
        return sample

    def setup(self):
        env = {}
        path = self.root / ".env"
        if path.is_file():
            for line in path.read_text(encoding="utf-8").splitlines():
                if line.strip() and not line.lstrip().startswith("#") and "=" in line:
                    key, value = line.split("=", 1)
                    env[key.strip()] = value.strip().strip("\"'")
        runtime = self.coverage.probe()
        return {
            "mode": "local_workbench", "tools_connected": runtime["ready"], "coverage_runtime": runtime,
            "pipeline_runtime": self.pipeline.probe(),
            **self.model.status(),
            "review_sample_ids": self.review_sample_ids(),
            "qc_review_sample_ids": [sample_id for sample_id in self.review_sample_ids()
                                     if self.model.qc_review_allowed(sample_id)],
            "data_access": "local_read_only", "run_mode": "demo_coverage_agent_pipeline",
        }

    def list_runs(self):
        with self.lock:
            imported = self.completed_local_runs()
            items = copy.deepcopy(sorted({**imported, **self.runs}.values(), key=lambda item: item["created_at"], reverse=True))
            for item in items:
                self.include_run_artifacts(item)
            return items

    def archived_run_ids(self):
        entries = read_json(self.root / ".local" / "run_history_archive.json", [])
        if not isinstance(entries, list):
            return []
        return sorted({entry for entry in entries
                       if isinstance(entry, str) and re.fullmatch(r"run-[0-9a-f]{12}", entry)})

    def set_run_archived(self, run_id, archived):
        if type(archived) is not bool:
            raise ValueError('记录状态无效')
        with self.lock:
            self.get_run(run_id)
            ids = set(self.archived_run_ids())
            if archived:
                ids.add(run_id)
            else:
                ids.discard(run_id)
            path = self.root / '.local/run_history_archive.json'
            path.parent.mkdir(parents=True, exist_ok=True)
            temporary = path.with_suffix('.tmp')
            temporary.write_text(json.dumps(sorted(ids), indent=2) + '\n', encoding='utf-8')
            temporary.replace(path)
            return {'archived_run_ids': sorted(ids)}

    def include_pipeline_logs(self, item):
        """Log files remain readable even when incomplete analysis outputs are withheld."""
        if item.get("mode") != "pipeline":
            return
        folder = self.pipeline.output_dir(item["id"])
        logs = []
        for name in PIPELINE_LOGS:
            path = folder / name
            try:
                if not path.is_symlink() and path.is_file():
                    logs.append({"name": name, "bytes": path.stat().st_size})
            except OSError:
                continue
        item["artifacts"] = [a for a in item.get("artifacts", []) if a["name"] not in PIPELINE_LOGS] + logs

    def include_run_artifacts(self, item):
        self.include_pipeline_logs(item)
        if item.get('paired_adaptive') and item.get('track_review', {}).get('rerun'):
            figure = comparison_figure_info(self.root, item['id'])
            if figure:
                item['rerun_figure'] = {'name': FIGURE_NAME, 'bytes': figure[1]}

    def get_run(self, run_id):
        with self.lock:
            item = self.runs.get(run_id) or self.completed_local_runs().get(run_id)
            if not item:
                raise ValueError("没有找到这条运行记录")
            result = copy.deepcopy(item)
            self.include_run_artifacts(result)
            if item.get("mode") == "pipeline":
                reviews = []
                for path in (self.root / "outputs/agent_reviews").glob("review-*.json"):
                    if path.is_symlink() or not re.fullmatch(r"review-[0-9a-f]{12}", path.stem):
                        continue
                    review = read_json(path, {})
                    if isinstance(review, dict) and review.get("baseline_run") == run_id and review.get("model_requests_attempted"):
                        at = datetime.fromtimestamp(path.stat().st_mtime, timezone.utc).isoformat()
                        reviews.append({"id": path.stem, "at": at, **{k: review[k] for k in
                            ("status", "answer", "summary", "comparison", "error") if k in review}})
                result["linked_reviews"] = sorted(reviews, key=lambda row: row["at"])
            return result

    def completed_local_runs(self):
        imported = {}
        for subdir, filename in (("pipeline_runs", "resume-result.json"), ("spikein_runs", "record.json"),
                                 ("preset_runs", "record.json")):
            base = self.root / "outputs" / subdir
            for receipt in base.glob("run-*/" + filename):
                folder = receipt.parent
                run_id = folder.name
                if not re.fullmatch(r"run-[0-9a-f]{12}", run_id) or receipt.resolve() != receipt:
                    continue
                try:
                    modified = receipt.stat().st_mtime_ns
                    cached = self.local_runs.get(receipt)
                    if cached and cached[0] == modified:
                        imported[run_id] = cached[1]
                        continue
                    record = read_json(receipt, {})
                    if subdir == "pipeline_runs":
                        execution = record.get("execution", {})
                        prior = self.runs.get(execution.get("resumed_from"), {})
                        if prior.get("mode") != "pipeline" or record.get("status") != "completed_pipeline":
                            continue
                        summary, artifacts = self.pipeline.validate(folder, "resume")
                        finished = datetime.fromtimestamp(receipt.stat().st_mtime, timezone.utc)
                        log = {"name": "tool.log", "bytes": (folder / "tool.log").stat().st_size}
                        elapsed = execution.get("elapsed_seconds", 0)
                        item = {"id": run_id, "mode": "pipeline", "source": "resume", "status": "completed_pipeline", "progress": 100,
                                "name": "单样本基础分析 · 本地续跑", "sample_id": prior["sample_id"], "result": summary,
                                "created_at": (finished - timedelta(seconds=elapsed)).isoformat(), "updated_at": finished.isoformat(),
                                "execution": execution, "artifacts": artifacts + [log], "actual_model_requests": 0,
                                "message": "基础分析完成；复用本轮 FASTQ 已生成的 BAM。", "steps": [], "events": []}
                        item["report_metrics"] = self.pipeline.report_metrics(self.pipeline.output_dir(prior["id"]), "fastq")
                    elif subdir == "spikein_runs":
                        if (record.get("id") != run_id or record.get("status") not in ("completed_spikein", "failed")
                                or record.get("mode") != "spikein"):
                            continue
                        if record["status"] == "completed_spikein":
                            artifacts = record.get("artifacts", [])
                            if record.get("result", {}).get("validation") != "passed" or {a["name"] for a in artifacts} != set(SPIKEIN_ARTIFACTS):
                                continue
                            for artifact in artifacts:
                                path = folder / artifact["name"]
                                if not path.is_file() or path.resolve() != path or path.stat().st_size != artifact["bytes"] or sha256(path) != artifact["sha256"]:
                                    raise ValueError("轨道产物校验未通过")
                        item = record
                    else:
                        item = self.preset_item(record, folder)
                    imported[run_id] = item
                    self.local_runs[receipt] = modified, item
                except (ValueError, OSError, KeyError, TypeError, OverflowError):
                    continue
        return imported

    @staticmethod
    def preset_item(record, folder):
        if not isinstance(record, dict):
            raise ValueError("无效的预设记录")
        preset, status = record.get("preset"), record.get("status")
        allowed = {"normalization": "completed_comparison", "igg": "blocked_no_matched_igg"}
        if (preset not in allowed or record.get("id") != folder.name
                or status not in (allowed[preset], "failed", "cancelled")):
            raise ValueError("无效的预设记录")
        item = {**record, "mode": "preset", "name": PRESET_NAMES[preset], "source": "terminal", "artifacts": [],
                "updated_at": record.get("updated_at", record["created_at"])}
        if status != "completed_comparison":
            return item
        track_dir = folder / "tracks"
        path = track_dir / "record.json"
        if path.resolve() != path:
            raise ValueError("轨道记录路径无效")
        tracks = read_json(path, {})
        artifacts = tracks.get("artifacts", [])
        if (tracks.get("status") != "completed_comparison" or tracks.get("result", {}).get("validation") != "passed"
                or len(artifacts) != len(PRESET_ARTIFACTS) or {row["name"] for row in artifacts} != PRESET_ARTIFACTS.keys()):
            raise ValueError("标准化对比产物未通过检查")
        for row in artifacts:
            path = track_dir / row["name"]
            if not path.is_file() or path.resolve() != path or path.stat().st_size != row["bytes"] or sha256(path) != row["sha256"]:
                raise ValueError("标准化对比产物已改变")
        item.update(artifacts=artifacts, bin_size=tracks["result"]["cpm"]["bin_size"],
                    spikein_fragments=tracks["spikein_fragments"], scale_factor=tracks["scale_factor"])
        return item

    @staticmethod
    def cancel_preset_container(item):
        if item.get("preset") == "normalization" and re.fullmatch(r"run-[0-9a-f]{12}", item["id"]):
            try:
                cancel_tracks(item["id"])
            except (OSError, subprocess.SubprocessError):
                pass

    def save(self, item):
        path = self.output_root / (item["id"] + ".json")
        temporary = path.with_suffix(".tmp")
        temporary.write_text(json.dumps(item, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
        temporary.replace(path)

    def create_run(self, payload):
        allowed = {"sample_id", "name", "mode", "target_bin_size", "source", "baseline_id", "model_review", "preset", "review_id", "defer_rerun"}
        if not isinstance(payload, dict) or set(payload) - allowed:
            raise ValueError("任务参数不受支持")
        mode = payload.get("mode", "demo")
        if mode == "paired_agent":
            if payload.get("sample_id") != PAIR_TARGET:
                raise ValueError("自动调参目前只支持已登记的配对结果")
            return self.create_track_review(payload, adaptive=True)
        if "defer_rerun" in payload:
            raise ValueError("延后重跑只用于配对自动调参")
        if payload.get("sample_id") == PAIR_TARGET and mode != "review":
            raise ValueError("这对样本支持复核已有结果，请选择质控复核")
        if mode == "preset":
            return self.create_preset_run(payload)
        if mode not in ("demo", "coverage", "agent", "pipeline", "review") or set(payload) & {"preset", "review_id"}:
            raise ValueError("不支持此任务模式")
        if mode == "review":
            return self.create_track_review(payload)
        if mode == "pipeline":
            return self.create_pipeline_run(payload)
        if mode == "agent":
            return self.create_agent_run(payload)
        if mode == "coverage":
            return self.create_coverage_run(payload)
        if set(payload) & {"source", "baseline_id", "model_review"}:
            raise ValueError("此任务模式不支持输入起点参数")
        sample = self.sample(str(payload.get("sample_id", "")))
        if not sample["report_available"]:
            raise ValueError("需要先准备该样本的 fastp 报告")
        size = payload.get("target_bin_size", 25)
        if type(size) is not int or size not in (25, 50):
            raise ValueError("演示参数仅支持 25 或 50")
        name = str(payload.get("name", "单样本流程演示")).strip()
        if not 1 <= len(name) <= 80:
            raise ValueError("任务名称请填写 1 至 80 个字符")
        with self.lock:
            if sum(item["status"] in ("running", "queued") for item in self.runs.values()) >= 3:
                raise ValueError("请等待当前演示结束后再开始")
            run_id = "run-" + uuid.uuid4().hex[:12]
            item = {
                "id": run_id, "name": name, "sample_id": sample["id"],
                "mode": "demo", "status": "queued", "created_at": now(),
                "updated_at": now(), "progress": 0, "max_reruns": 1,
                "actual_tool_executions": 0, "actual_model_requests": 0,
                "message": "演示已创建，尚未进行真实重新分析。",
                "baseline": {"source": "historical_report", "metrics": sample["metrics"]},
                "plan": {"source": "fixed_demo", "action": "regenerate_coverage", "bin_size_before": 10, "bin_size_after": size},
                "steps": [{"id": key, "label": label, "detail": detail, "status": "pending"} for key, label, detail in STEPS],
                "events": [{"at": now(), "message": "创建本地流程演示。工具与模型均未调用。"}],
            }
            self.runs[run_id] = item
            self.save(item)
        threading.Thread(target=self._drive, args=(run_id,), daemon=True).start()
        return self.get_run(run_id)

    def create_coverage_run(self, payload):
        if set(payload) - {"sample_id", "name", "mode"}:
            raise ValueError("真实试跑使用固定区域和参数，不接受额外执行参数")
        sample = self.sample(str(payload.get("sample_id", "")))
        self.coverage.input_bam(sample["id"])
        if not self.coverage.probe()["ready"]:
            raise ValueError("本地 bamCoverage 工具尚未就绪")
        name = str(payload.get("name", "区域覆盖轨道试跑")).strip()
        if not 1 <= len(name) <= 80:
            raise ValueError("任务名称请填写 1 至 80 个字符")
        with self.lock:
            if any(item.get("mode") in REAL_MODES and item["status"] in ("queued", "running") for item in self.runs.values()):
                raise ValueError("一次只运行一个真实试跑，请等待或停止当前任务")
            run_id = "run-" + uuid.uuid4().hex[:12]
            item = {"id": run_id, "name": name, "sample_id": sample["id"], "mode": "coverage",
                    "status": "queued", "created_at": now(), "updated_at": now(), "progress": 0,
                    "max_reruns": 0, "actual_tool_executions": 0, "actual_model_requests": 0,
                    "scope": "regional_tool_test", "message": "已创建真实区域试跑，尚未调用模型。",
                    "parameters": {"tool": "bamCoverage", "region": "chr1:1000000:2000000", "bin_size": 25,
                                   "normalization": "None", "scale_factor": 1, "sam_flag_include": 64},
                    "steps": [{"id": key, "label": label, "detail": detail, "status": "pending"} for key, label, detail in [
                        ("inputs", "校验输入", "确认本地 BAM 和索引，记录输入元数据"),
                        ("coverage", "生成区域覆盖轨道", "在无网络容器中真实运行 bamCoverage"),
                        ("validate", "检查并保存结果", "校验 bigWig 可读性、区域、大小与摘要")]],
                    "events": [{"at": now(), "message": "创建真实区域试跑；不调用模型。"}]}
            self.runs[run_id] = item
            self.save(item)
        threading.Thread(target=self._drive_coverage, args=(run_id,), daemon=True).start()
        return self.get_run(run_id)

    def create_pipeline_run(self, payload):
        if set(payload) - {"sample_id", "name", "mode", "source", "model_review"} or type(payload.get("model_review", False)) is not bool:
            raise ValueError("基础分析不接受自定义命令、路径或额外参数")
        sample = self.sample(str(payload.get("sample_id", "")))
        source = payload.get("source", "fastq")
        self.pipeline.inputs(sample["id"], source)
        if not self.pipeline.probe().get("ready"):
            raise ValueError("本地基础分析工具尚未准备好")
        model_review = payload.get("model_review", False)
        name = str(payload.get("name", "单样本基础分析")).strip()
        if not 1 <= len(name) <= 80:
            raise ValueError("任务名称请填写 1 至 80 个字符")
        with self.lock:
            if any(r.get("mode") in REAL_MODES and r["status"] in ("queued", "running") for r in self.runs.values()):
                raise ValueError("请等待当前真实任务结束后再创建")
            run_id = "run-" + uuid.uuid4().hex[:12]
            item = {"id": run_id, "sample_id": sample["id"], "name": name, "mode": "pipeline", "source": source,
                    "status": "queued", "created_at": now(), "updated_at": now(), "progress": 0,
                    "max_reruns": 0, "actual_tool_executions": 0, "actual_model_requests": 0,
                    "model_review_requested": model_review,
                    "message": "已创建单样本基础分析，所有计算在本地进行。", "scope": "single_sample_basic_pipeline",
                    "steps": [{"id": key, "label": label, "detail": "本地真实计算",
                               "status": "skipped" if source == "bam" and key in ("fastp", "align", "bam") else "pending"}
                              for key, label in PIPELINE_STAGES],
                    "events": [{"at": now(), "message": "基础分析通过后进入模型复核。" if model_review else "基础分析不请求模型。"}]}
            self.runs[run_id] = item
            self.save(item)
        threading.Thread(target=self._drive_pipeline, args=(run_id,), daemon=True).start()
        return self.get_run(run_id)

    def _drive_pipeline(self, run_id):
        def phase(key, message, execution):
            with self.lock:
                item = self.runs[run_id]
                for step in item["steps"]:
                    if step["status"] == "running":
                        step["status"] = "done"
                    if step["id"] == key:
                        step["status"] = "running"
                index = next(i for i, row in enumerate(PIPELINE_STAGES) if row[0] == key)
                item.update(status="running", progress=round(index / len(PIPELINE_STAGES) * 95),
                            message=message, execution=copy.deepcopy(execution), updated_at=now())
                item["events"].append({"at": now(), "message": message})
                self.save(item)
        try:
            item = self.get_run(run_id)
            result = self.pipeline.run(run_id, item["sample_id"], item["source"], lambda: run_id in self.cancelled, phase)
        except Exception as exc:
            result = {"status": "failed", "error": str(exc)}
        with self.lock:
            item = self.runs[run_id]
            item.update(result)
            if run_id in self.cancelled:
                item.update(status="cancelled", artifacts=[a for a in item.get("artifacts", []) if a["name"] in PIPELINE_LOGS])
                item.pop("result", None)
            success = item["status"] == "completed_pipeline"
            item["message"] = "基础分析完成；候选结果与本次图表可查看和下载。" if success else {
                "cancelled": "基础分析已停止，未发布未完成结果。", "timed_out": "超过本次时间上限，已停止。"
            }.get(item["status"], "基础分析失败，请查看错误和运行日志。")
            if item["status"] == "failed":
                item["error"] = result.get("error") or result.get("execution", {}).get("error", "工具或结果检查失败")
            for step in item["steps"]:
                if success and step["status"] != "skipped":
                    step["status"] = "done"
                elif step["status"] == "running":
                    step["status"] = "cancelled" if item["status"] == "cancelled" else "failed"
            item.update(updated_at=now(), progress=100 if success else item["progress"])
            will_review = success and item.get("model_review_requested")
            if will_review:
                item.update(baseline_ready=True, status="running", progress=96, message="基础分析已完成，正在整理质控证据。")
            item["events"].append({"at": now(), "message": item["message"]})
            self.save(item)
        if will_review:
            self.drive_track_review(run_id, run_id, item["source"])

    def review_sample_ids(self):
        scope = read_json(self.root / ".local/model_summary_authorization.json", {})
        return scope.get("sample_ids", []) if isinstance(scope, dict) and isinstance(scope.get("sample_ids"), list) else []

    def latest_spikein(self, sample_id):
        with self.lock:
            return max((run for run in self.completed_local_runs().values()
                        if run.get("mode") == "spikein" and run.get("status") == "completed_spikein"
                        and run.get("sample_id") == sample_id), key=lambda run: run["created_at"], default=None)

    def sample_chat_evidence(self, sample_id, reference):
        runs = [run for run in self.list_runs() if run["sample_id"] == sample_id]
        if sample_id == PAIR_TARGET:
            current = pair_review.adaptive_payload(self.hospital_pair.sample()["pair"])
            chosen = reference or next((run for run in runs if run.get("track_review", {}).get("rerun")), {})
            review = chosen.get("track_review", {})
            saved = review.get("evidence") or {}
            basis = saved if "comparison" in saved else current
            evidence = basis
            if review.get("rerun"):
                evidence = {**basis, "comparison": {**basis["comparison"], "rerun": review["rerun"]}}
            return evidence, None, basis == current

        self.sample_path(sample_id)
        baseline_id = reference.get("baseline_id") if reference else None
        if reference and reference["mode"] == "pipeline":
            baseline_id = reference["id"]
        if baseline_id:
            baseline = self.get_run(baseline_id)
        else:
            baseline = next((run for run in runs if run["mode"] == "pipeline"
                             and run["status"] == "completed_pipeline"), {})
        if (baseline.get("sample_id") != sample_id or baseline.get("mode") != "pipeline"
                or (baseline.get("status") != "completed_pipeline" and not baseline.get("baseline_ready"))
                or baseline.get("result", {}).get("validation") != "passed"):
            raise ValueError("当前这条运行还没有可引用的完整质控结果。")
        baseline_id = baseline["id"]
        candidates = ([reference] if reference else []) + runs
        for run in candidates:
            if run.get("baseline_id") == baseline_id or run["id"] == baseline_id:
                evidence = run.get("track_review", {}).get("evidence", {})
                if evidence.get("data_source") == "cuttag_qc_evidence_v1":
                    return evidence, baseline_id, True
        evidence, _ = make_qc_evidence(self.pipeline.output_dir(baseline_id), baseline["source"], sample_id,
                                       self.latest_spikein(sample_id))
        return evidence, baseline_id, True

    def ask_sample_chat(self, payload):
        fields = {"sample_id", "run_id", "question", "history"}
        if not isinstance(payload, dict) or set(payload) - fields:
            raise ValueError("对话请求的字段无效。")
        sample_id, reference_id = payload.get("sample_id"), payload.get("run_id")
        question, history = payload.get("question"), payload.get("history", [])
        if not isinstance(sample_id, str) or not sample_id:
            raise ValueError("请先选择要询问的样本。")
        if not isinstance(question, str) or not 1 <= len(question.strip()) <= 500:
            raise ValueError("请输入 1 至 500 字的问题。")
        if (not isinstance(history, list) or len(history) > 2
                or any(not isinstance(turn, dict) or set(turn) != {"question", "reply"} for turn in history)):
            raise ValueError("对话上下文格式无效。")
        reference = None
        if reference_id is not None:
            if not isinstance(reference_id, str) or not re.fullmatch(r"run-[0-9a-f]{12}", reference_id):
                raise ValueError("运行记录无效。")
            reference = self.get_run(reference_id)
            if reference.get("sample_id") != sample_id:
                raise ValueError("运行记录和已选样本不一致。")

        labels = {sample_id: "当前样本"}
        if sample_id == PAIR_TARGET:
            labels[PAIR_CONTROL] = "匹配 IgG"
        if reference:
            labels[reference["id"]] = "当前运行"
        known_samples = {path.name for path in self.data_root.glob("*") if path.is_dir()} | {PAIR_TARGET, PAIR_CONTROL}
        other_samples = known_samples - labels.keys()
        try:
            evidence, baseline_id, allow_actions = self.sample_chat_evidence(sample_id, reference)
            answer = self.model.ask_sample(sample_id, evidence, question, history, labels, other_samples,
                                           allow_actions=allow_actions)
            decision = answer["decision"]
            result = {"sample_id": sample_id, "reply": decision["reply"]}
            proposal = decision["proposal"]
            if proposal is None:
                return result
            if sample_id == PAIR_TARGET:
                result["suggested_run"] = self.save_chat_rerun(evidence, answer, reference_id)
                return result
            self.artifact(baseline_id, "coverage.cpm.bw")
            comparison = compare_evidence(self.pipeline.output_dir(baseline_id), evidence, proposal)
            evaluated = {**evidence, "comparison": comparison}
            validate_qc_payload(evaluated)
            result["comparison"] = {"evidence": evidence, "comparison": comparison}
        except (OSError, subprocess.SubprocessError) as exc:
            raise ValueError("本地质控证据读取或比较未完成。") from exc
        try:
            follow_up = self.model.ask_sample(sample_id, evaluated, question, history, labels, other_samples,
                                              allow_actions=False)
            result["reply"] = follow_up["decision"]["reply"]
        except ValueError:
            result["notice"] = "已完成本地比较。模型复评暂未返回。"
        return result

    def save_chat_rerun(self, evidence, answer, reference_id):
        basis = {**evidence, "comparison": {key: value for key, value in evidence["comparison"].items() if key != "rerun"}}
        decision = pair_review.validate_decision(answer["decision"]["proposal"], basis, adaptive=True)
        status = "awaiting_rerun"
        review = {"status": status, "model_requests_attempted": 1, "full_pipeline_executions": 0, "peak_executions": 0,
                  "evidence_schema": pair_review.SCHEMA, "evidence": basis, "comparison": basis["comparison"],
                  "sent_summaries": [evidence], "answer": {**answer, "decision": decision}}
        run_id, at = "run-" + uuid.uuid4().hex[:12], now()
        message = "对话中提出一项找峰调整，等待你确认是否重跑。"
        item = {"id": run_id, "name": "对话建议 · 带 IgG 找峰", "mode": "review", "sample_id": PAIR_TARGET,
                "chat_reference_id": reference_id, "status": status, "created_at": at, "updated_at": at,
                "paired_review": True, "paired_adaptive": True, "actual_model_requests": 1, "track_review": review,
                "message": message, "events": [{"at": at, "message": message}]}
        with self.lock:
            self.save(item)
            self.runs[run_id] = item
        return self.get_run(run_id)

    def create_track_review(self, payload, adaptive=False):
        fields = {"sample_id", "name", "mode", "baseline_id"}
        if adaptive:
            fields.add("defer_rerun")
        if set(payload) - fields:
            raise ValueError("复核不接受目标窗口、自定义参数或路径")
        if adaptive and payload.get("defer_rerun", True) is not True:
            raise ValueError("重跑需要先接受 Agent 的建议")
        baseline_id = payload.get("baseline_id")
        sample_id = payload.get("sample_id")
        paired = sample_id == PAIR_TARGET
        if paired:
            if "baseline_id" in payload:
                raise ValueError("配对复核使用已登记的全基因组结果")
            if not self.model.pair_review_allowed():
                raise ValueError("配对样本摘要的外部模型使用范围待确认")
            model = self.model.status()
            if not (model["model_connected"] and model["minimal_summary_allowed"]):
                raise ValueError("请先确认模型连接和摘要使用范围")
            source = None
        else:
            if not isinstance(baseline_id, str) or not re.fullmatch(r"run-[0-9a-f]{12}", baseline_id):
                raise ValueError("请选择已完成的全基因组结果")
            baseline = self.get_run(baseline_id)
            if (baseline.get("mode") != "pipeline" or baseline.get("status") != "completed_pipeline"
                    or baseline.get("sample_id") != sample_id):
                raise ValueError("请选择与当前样本匹配的已完成全基因组结果")
            source = baseline["source"]
        name = str(payload.get("name", "已有结果 · 模型复核")).strip()
        if not 1 <= len(name) <= 80:
            raise ValueError("任务名称请填写 1 至 80 个字符")
        with self.lock:
            if any(r.get("mode") in REAL_MODES and r["status"] in ("queued", "running") for r in self.runs.values()):
                raise ValueError("请等待当前任务结束后再复核")
            run_id = "run-" + uuid.uuid4().hex[:12]
            item = {"id": run_id, "mode": "review", "sample_id": sample_id, "name": name, "baseline_id": baseline_id,
                    "status": "queued", "actual_model_requests": 0, "created_at": now(), "updated_at": now(),
                    "message": "复用已完成的全基因组分析；不重新比对。", "events": []}
            if paired:
                item["paired_review"] = True
            if adaptive:
                item["paired_adaptive"] = True
            self.runs[run_id] = item
            self.save(item)
        threading.Thread(target=self.drive_track_review, args=(run_id, baseline_id, source), daemon=True).start()
        return self.get_run(run_id)

    def decide_paired_rerun(self, run_id, accept):
        if type(accept) is not bool:
            raise ValueError('请选择接受建议或保留原结果')
        with self.lock:
            item = self.runs.get(run_id, {})
            advice = item.get('track_review', {}).get('answer', {}).get('decision', {})
            if (not item.get('paired_adaptive') or item.get('status') != 'awaiting_rerun'
                    or advice.get('action') != 'tune_igg_qvalue'):
                raise ValueError('这条记录当前没有可继续的调参重跑')
            if accept:
                self.check_paired_rerun_ready()
                status, message = 'queued', '已接受 Agent 的建议，正在排队重跑找峰。'
            else:
                status, message = 'completed_review', '已选择保留原结果，本次未重跑。'
                item['track_review']['status'] = status
            decided_at = now()
            item.update(status=status, message=message, updated_at=decided_at,
                        paired_rerun_accepted=accept)
            item['events'].append({'at': decided_at, 'message': message})
            self.save(item)
        if accept:
            threading.Thread(target=self.drive_track_review, args=(run_id, None, None, True), daemon=True).start()
        return self.get_run(run_id)

    def resume_paired_rerun(self, run_id):
        return self.decide_paired_rerun(run_id, True)

    def check_paired_rerun_ready(self):
        model = self.model.status()
        if not self.model.pair_review_allowed() or not (model['model_connected'] and model['minimal_summary_allowed']):
            raise ValueError('请先确认模型连接和摘要使用范围')
        if any(r.get('mode') in REAL_MODES and r['status'] in ('queued', 'running') for r in self.runs.values()):
            raise ValueError('请等待当前任务结束后再重跑')

    def replay_paired_rerun(self, source_id):
        with self.lock:
            source = self.runs.get(source_id, {})
            review = source.get('track_review', {})
            if (source.get('sample_id') != PAIR_TARGET or not source.get('paired_adaptive')
                    or source.get('status') != 'completed_review' or review.get('peak_executions') != 1):
                raise ValueError('这条记录没有已执行并可复用的找峰参数')
            self.check_paired_rerun_ready()
            evidence = review['evidence']
            pair = self.hospital_pair.sample()['pair']
            if 'igv_screen' not in evidence['comparison']:
                pair['igv_screen'] = None
            if evidence != pair_review.adaptive_payload(pair):
                raise ValueError('原记录对应的输入证据已变化，请从工作台重新做调参判断')
            advice = pair_review.validate_decision(review['answer']['decision'], evidence, adaptive=True)
            if advice['action'] != 'tune_igg_qvalue' or advice['qvalue'] != review.get('rerun', {}).get('qvalue'):
                raise ValueError('找峰参数与原执行记录不一致')

            run_id, created = 'run-' + uuid.uuid4().hex[:12], now()
            message = f"已接受沿用上次参数 q={advice['qvalue']}，正在排队重跑。"
            track = {'evidence_schema': pair_review.SCHEMA, 'evidence': evidence, 'comparison': evidence['comparison'],
                     'answer': {**review['answer'], 'decision': advice}, 'sent_summaries': [], 'model_requests_attempted': 0,
                     'status': 'reviewing', 'full_pipeline_executions': 0, 'peak_executions': 0}
            item = {'id': run_id, 'sample_id': PAIR_TARGET, 'name': 'IgG 找峰 · 沿用上次参数', 'mode': 'review',
                    'status': 'queued', 'created_at': created, 'updated_at': created, 'actual_model_requests': 0,
                    'paired_review': True, 'paired_adaptive': True, 'paired_rerun_accepted': True, 'reused_plan_id': source_id,
                    'message': message, 'events': [{'at': created, 'message': message}], 'track_review': copy.deepcopy(track)}
            self.runs[run_id] = item
            self.save(item)
        threading.Thread(target=self.drive_track_review, args=(run_id, None, None, True), daemon=True).start()
        return self.get_run(run_id)

    def drive_track_review(self, run_id, baseline_id, source, resume=False):
        def update(review):
            with self.lock:
                item = self.runs[run_id]
                if review.get("model_requests_attempted", 0) > item.get("actual_model_requests", 0):
                    summaries = review.get("sent_summaries", [])
                    message = "已发起完整质控摘要请求。"
                    if item.get("paired_adaptive") and review.get("rerun"):
                        message = "已发起新旧找峰结果的模型复评。"
                    elif not item.get("paired_adaptive") and review.get("comparison"):
                        message = "已发起对照结果复评请求。"
                    if item.get("paired_review") and len(summaries) >= 2 and summaries[-1] == summaries[-2]:
                        message = "上一次模型回复未通过检查，已用同一份摘要重试。"
                    item["events"].append({"at": now(), "message": message})
                if review.get("answer") and not item.get("track_review", {}).get("answer"):
                    item["events"].append({"at": now(), "message": "模型回复已返回并通过格式、动作范围检查。"})
                if not item.get("paired_adaptive") and review.get("comparison") and not item.get("track_review", {}).get("comparison"):
                    item["events"].append({"at": now(), "message": "已复用原产物完成模型要求的比较。"})
                if review.get("rerun_started") and not item.get("track_review", {}).get("rerun_started"):
                    item["message"] = "正在使用同一对 BAM 重跑带 IgG 的 MACS3 找峰。"
                    item["events"].append({"at": now(), "message": item["message"]})
                if review.get("rerun") and not item.get("track_review", {}).get("rerun"):
                    item["message"] = "新候选峰和落峰比例已计算，正在复评。"
                    item["events"].append({"at": now(), "message": item["message"]})
                item.update(status="running", updated_at=now(), track_review=copy.deepcopy(review),
                            actual_model_requests=review.get("model_requests_attempted", 0))
                self.save(item)
        with self.lock:
            item = self.runs[run_id]
            previous = copy.deepcopy(item["track_review"]) if resume else None
            item.update(status="running", message="正在整理已校验的本地质控证据。")
            if not resume:
                item['track_review'] = {"status": "reviewing"}
            item["events"].append({"at": now(), "message": item["message"]})
            self.save(item)
        sample_id = item["sample_id"]
        try:
            def cancelled():
                return run_id in self.cancelled

            if item.get("paired_adaptive"):
                pair = self.hospital_pair.sample()["pair"]
                if resume:
                    if 'igv_screen' not in previous['evidence']['comparison']:
                        pair['igv_screen'] = None
                else:
                    with self.lock:
                        item.update(evidence_preparing=True, updated_at=now(),
                                    message="正在自动核对轨道与两套候选峰，准备模型所需证据。")
                        item['events'].append({"at": now(), "message": item['message']})
                        self.save(item)
                    pair['igv_screen'] = self.hospital_pair.prepare_screen(pair, run_id, cancelled)
                    with self.lock:
                        item.update(evidence_preparing=False, updated_at=now(), message="本地证据已准备，正在等待 Agent 决定。")
                        item['events'].append({"at": now(), "message": item['message']})
                        self.save(item)
                result = review_adaptive_pair(pair, self.model,
                                              lambda q, stop: rerun_peaks(run_id, q, stop, inspect_igv=bool(pair.get('igv_screen'))),
                                              cancelled, update,
                                              defer=not resume, previous=previous)
                if result.get('rerun') and not cancelled() and (rerun_dir(self.root, run_id) / 'tuned.peaks.bed').is_file():
                    try:
                        build_comparison_figure(self.root, run_id, cancelled)
                    except (OSError, ValueError, KeyError, TypeError, subprocess.SubprocessError):
                        result['figure_error'] = '轨道对照图生成失败，指标和 BED 已保存。'
            elif item.get("paired_review"):
                pair = self.hospital_pair.sample()["pair"]
                result = review_evidence(pair_review.SCHEMA, lambda: (pair_review.make_payload(pair), {}),
                                         self.model.pair_review_allowed, self.model.status, self.model.review_pair,
                                         lambda _payload, decision: pair_review.compare_pair(pair, decision), cancelled, update)
            else:
                result = review_baseline(self.pipeline.output_dir(baseline_id), source, self.model, cancelled,
                                         update, sample_id=sample_id, spikein=self.latest_spikein(sample_id))
        except (OSError, KeyError, TypeError, subprocess.SubprocessError):
            result = {"model_requests_attempted": 0, **(previous or {}),
                      "status": "cancelled" if cancelled() else "failed", "error": "本地轨道证据无法读取或未在时限内生成"}
        except ValueError as exc:
            result = {"model_requests_attempted": 0, **(previous or {}), "status": "failed", "error": str(exc)}
        with self.lock:
            item = self.runs[run_id]
            item.pop('evidence_preparing', None)
            result_status = result["status"]
            message = {"completed_review": "模型复核完成。", "cancelled": "模型复核已停止。",
                       "awaiting_rerun": "Agent 已给出调参建议，等你确认是否重跑。",
                       "evidence_ready": "质控证据已在本地整理。完整摘要的外部使用范围待确认，本次未发起模型请求。"
                       }.get(result_status, "模型复核未完成，请查看决策记录。")
            if result.get('existing_comparison_selected'):
                message = '已比较两套现有候选峰。'
            if item["mode"] == "pipeline":
                item.update(status="completed_pipeline", progress=100, message="基础分析产物已保留。" + message)
            else:
                item.update(status=result_status, message=message)
            item["events"].append({"at": now(), "message": message})
            item.update(track_review=result, actual_model_requests=result["model_requests_attempted"], updated_at=now())
            self.save(item)

    def create_preset_run(self, payload):
        preset = payload.get("preset")
        if set(payload) - {"mode", "preset", "sample_id", "review_id", "name"} or preset not in PRESET_NAMES:
            raise ValueError("不支持此预设或参数")
        review_id = payload.get("review_id")
        review, summary = review_source(self.root, review_id)
        sample_id = review["sample_id"]
        if sample_id != payload.get("sample_id"):
            raise ValueError("请选择属于当前样本的已完成质控复核")
        self.sample_path(sample_id)
        validate_qc_payload(summary)
        if not self.model.qc_review_allowed(sample_id) or not self.model.status()["model_connected"]:
            raise ValueError("请先确认该样本的摘要使用范围和模型连接")
        qc = summary["qc"]
        if preset == "normalization" and not qc.get("spikein", {}).get("spikein_fragments"):
            raise ValueError("该复核缺少配套 Spike-in 计数")
        if preset == "igg" and qc["experiment_context"]["matched_igg"] != "unavailable":
            raise ValueError("该预设用于确认当前样本无匹配 IgG 的停止分支")
        name = str(payload.get("name", PRESET_NAMES[preset])).strip()
        if not 1 <= len(name) <= 80:
            raise ValueError("任务名称请填写 1 至 80 个字符")
        with self.lock:
            if any(r.get("mode") in REAL_MODES and r["status"] in ("queued", "running") for r in self.runs.values()):
                raise ValueError("请等待当前真实任务结束后再启动预设")
            run_id = "run-" + uuid.uuid4().hex[:12]
            created = now()
            item = {"id": run_id, "mode": "preset", "preset": preset, "source": "web", "name": name,
                    "sample_id": sample_id, "review_id": review_id, "baseline_id": review["baseline_id"],
                    "status": "queued", "created_at": created, "updated_at": created,
                    "model_requests": 0, "fastq_executions": 0, "peak_executions": 0,
                    "message": "准备复用已校验的质控摘要。"}
            self.runs[run_id] = item
            self.save(item)
        threading.Thread(target=self.drive_preset_run, args=(run_id, preset, review_id), daemon=True).start()
        return self.get_run(run_id)

    def drive_preset_run(self, run_id, preset, review_id):
        def update(trace):
            with self.lock:
                item = self.runs[run_id]
                item.update(copy.deepcopy({key: trace[key] for key in ("sent_summary", "answer", "model_requests") if key in trace}))
                item.update(status="running", updated_at=now(),
                            message="模型建议已返回，正在执行预设检查或计算。" if trace.get("answer") else "正在等待模型建议。")
                self.save(item)

        update({})
        try:
            record = run_preset(self.root, preset, review_id, run_id=run_id, model=self.model,
                                cancelled=lambda: run_id in self.cancelled, update=update)
            result = self.preset_item(record, self.root / "outputs/preset_runs" / run_id)
            for key in ("name", "created_at", "source"):
                result.pop(key)
        except (ValueError, OSError, TypeError, KeyError, subprocess.SubprocessError) as exc:
            result = {"status": "failed", "error": str(exc), "artifacts": []}
        with self.lock:
            item = self.runs[run_id]
            item.update(result)
            item.update(updated_at=now(), message={
                "completed_comparison": "两份标准化轨道已生成并通过检查。",
                "blocked_no_matched_igg": "已确认无匹配 IgG，保留原候选峰，未启动找峰。",
                "cancelled": "预设已停止。",
            }.get(item["status"], "预设未完成，请查看错误。"))
            self.save(item)

    def create_agent_run(self, payload):
        if set(payload) - {"sample_id", "name", "mode"}:
            raise ValueError("Agent 任务不接受额外命令、路径或参数")
        sample = self.sample(str(payload.get("sample_id", "")))
        self.coverage.input_bam(sample["id"])
        runtime, model = self.coverage.probe(), self.model.status()
        if not runtime.get("ready") or not runtime.get("samtools_ready"):
            raise ValueError("请先准备 bamCoverage 与 SAMtools 环境")
        if not model["model_connected"]:
            raise ValueError("请先在接入设置完成模型的合成连通测试")
        if not model["minimal_summary_allowed"]:
            raise ValueError("最小真实摘要的 API 使用范围尚未确认")
        name = str(payload.get("name", "单样本 Agent 闭环")).strip()
        if not 1 <= len(name) <= 80:
            raise ValueError("任务名称请填写 1 至 80 个字符")
        with self.lock:
            if any(r.get("mode") in REAL_MODES and r["status"] in ("queued", "running") for r in self.runs.values()):
                raise ValueError("请等待当前真实任务结束后再创建")
            run_id = "run-" + uuid.uuid4().hex[:12]
            item = {"id": run_id, "sample_id": sample["id"], "name": name, "mode": "agent", "status": "queued",
                    "created_at": now(), "updated_at": now(), "progress": 0, "max_reruns": 1,
                    "actual_tool_executions": 0, "actual_model_requests": 0, "rounds": [], "decisions": [], "model_requests": [],
                    "scope": "regional_engineering_agent", "message": "已创建真实 Agent 工程闭环任务。",
                    "steps": [{"id": key, "label": label, "detail": detail, "status": "pending"} for key, label, detail in [
                        ("baseline", "真实基线计算", "SAMtools 检查与统计，生成 25 bp BW"),
                        ("decision", "模型评估与建议", "仅发送已获准的最小运行摘要"),
                        ("rerun", "一次受控重跑", "程序检查建议后执行允许的覆盖窗口调整"),
                        ("review", "模型复评与结束", "重读结果，最多一次重跑，明确结束原因")]],
                    "events": [{"at": now(), "message": "创建 Agent 工程联通任务；非科学能力验收。"}]}
            self.runs[run_id] = item
            self.save(item)
        threading.Thread(target=self._drive_agent, args=(run_id,), daemon=True).start()
        return self.get_run(run_id)

    def _drive_agent(self, run_id):
        def update(key, message, partial):
            with self.lock:
                item = self.runs[run_id]
                for step in item["steps"]:
                    if step["status"] == "running" and step["id"] != key:
                        step["status"] = "done"
                    if step["id"] == key:
                        step["status"] = "running"
                item.update(partial)
                item.update(status="running", message=message, updated_at=now(), progress={"baseline": 10, "decision": 40, "rerun": 65, "review": 90}[key])
                item["events"].append({"at": now(), "message": message})
                self.save(item)
        try:
            result = drive_agent(run_id, self.runs[run_id]["sample_id"], self.coverage, self.model, lambda: run_id in self.cancelled, update)
        except Exception:
            result = {"status": "failed", "end_reason": "Agent 运行出现错误，未继续执行；请检查工具和模型环境"}
        with self.lock:
            item = self.runs[run_id]
            item.update(result)
            item.update(updated_at=now(), message=result["end_reason"], progress=100 if item["status"] == "completed_agent" else item["progress"])
            for step in item["steps"]:
                if step["status"] == "running":
                    step["status"] = "done" if item["status"] == "completed_agent" else "cancelled" if item["status"] == "cancelled" else "failed"
                elif step["status"] == "pending" and item["status"] == "completed_agent":
                    step["status"] = "skipped"
            item["events"].append({"at": now(), "message": item["message"]})
            self.save(item)

    def _drive_coverage(self, run_id):
        def phase(key, execution):
            with self.lock:
                item = self.runs[run_id]
                for step in item["steps"]:
                    if step["status"] == "running":
                        step["status"] = "done"
                    if step["id"] == key:
                        step["status"] = "running"
                item.update(status="running", execution=copy.deepcopy(execution), updated_at=now(),
                            progress={"inputs": 10, "coverage": 25, "validate": 85}[key])
                item["message"] = {"inputs": "正在准备真实试跑输入。", "coverage": "正在真实生成区域覆盖轨道；不调用模型。", "validate": "工具已退出，正在校验真实产物。"}[key]
                item["events"].append({"at": now(), "message": item["message"]})
                self.save(item)
        try:
            result = self.coverage.run(run_id, self.runs[run_id]["sample_id"], lambda: run_id in self.cancelled, phase)
        except Exception as exc:
            result = {"status": "failed", "error": str(exc)}
        with self.lock:
            item = self.runs[run_id]
            item.update(result)
            if run_id in self.cancelled:
                item["status"] = "cancelled"
                item.pop("result", None)
            status = item["status"]
            item["message"] = {"completed": "区域覆盖轨道已真实生成并通过文件检查；非全基因组结果，未调用模型。",
                               "cancelled": "真实试跑已停止；未发布未完成的覆盖轨道。",
                               "timed_out": "真实试跑超过时间限制，已停止；请查看日志。",
                               "failed": "真实试跑失败，请查看错误原因与工具日志。"}.get(status, "任务已结束。")
            if status == "failed":
                item["error"] = result.get("error") or result.get("execution", {}).get("error", "工具或产物检查失败")
            for step in item["steps"]:
                if status == "completed":
                    step["status"] = "done"
                elif step["status"] == "running":
                    step["status"] = "cancelled" if status == "cancelled" else "failed"
            item.update(updated_at=now(), progress=100 if status == "completed" else item["progress"])
            item["events"].append({"at": now(), "message": item["message"]})
            self.save(item)

    def paired_igv_inputs(self, run):
        run_id = run["id"]
        review = run.get('track_review', {})
        if not review.get('rerun'):
            raise ValueError('这条记录没有重跑结果')
        tuned, _ = self.artifact(run_id, 'tuned.peaks.bed')
        figure_meta = read_json(rerun_dir(self.root, run_id) / FIGURE_DIR / LOCUS_NAME, {})
        loci = figure_meta.get('loci', [])
        if not loci:
            loci = [row['hg38_igv_locus'] for row in self.hospital_pair.sample()['pair']['manual_windows']]
        return tuned, [locus.replace(',', '') for locus in loci], review['rerun']['qvalue']

    def igv_reference_sizes(self):
        return self.root / REFERENCE.relative_to(PROJECT_ROOT) / REF_FILES[0]

    def igv_inputs(self, run_id):
        run = self.get_run(run_id)
        sample_id = run["sample_id"]
        if run.get("paired_adaptive"):
            tuned, loci, qvalue = self.paired_igv_inputs(run)
            base = f"/api/samples/{PAIR_TARGET}/artifacts/"
            tracks = [
                (name, fmt, color, base + alias, self.hospital_pair.artifact(alias)[0])
                for alias, name, fmt, color in (
                    ("target-cpm.bw", "Target CPM", "bigwig", "#244c6d"),
                    ("igg-cpm.bw", "IgG CPM", "bigwig", "#218981"),
                    ("with-igg.peaks.bed", "Before peaks (q=0.01)", "bed", "#738095"),
                )
            ]
            tracks.append((f"After peaks (q={qvalue:g})", "bed", "#c97526",
                           f"/api/runs/{run_id}/artifacts/tuned.peaks.bed", tuned))
            view = {"sample_id": sample_id, "loci": loci, "qvalue": qvalue}
        else:
            baseline = self.get_run(run["baseline_id"]) if run.get("baseline_id") else run
            if (baseline.get("mode") != "pipeline" or baseline.get("status") != "completed_pipeline"
                    or baseline.get("sample_id") != sample_id or baseline.get("result", {}).get("validation") != "passed"):
                raise ValueError("这条记录没有可复用的全基因组轨道")
            baseline_id = baseline["id"]
            base = f"/api/runs/{baseline_id}/artifacts/"
            tracks = [
                (name, fmt, color, base + filename, self.artifact(baseline_id, filename)[0])
                for filename, name, fmt, color in (
                    ("coverage.cpm.bw", "Target CPM", "bigwig", "#244c6d"),
                    ("narrow.peaks.bed", "Narrow peaks", "bed", "#738095"),
                    ("broad.peaks.bed", "Broad peaks", "bed", "#c97526"),
                )
            ]
            loci = [first_peak_locus(path) for _, fmt, _, _, path in tracks if fmt == "bed"]
            view = {"sample_id": sample_id, "loci": list(dict.fromkeys(locus for locus in loci if locus))}
        view["tracks"] = [
            {"name": name, "format": fmt, "color": color, "url": url} for name, fmt, color, url, _ in tracks
        ]
        return view, [path for _, _, _, _, path in tracks]

    def igv_session(self, run_id):
        view, paths = self.igv_inputs(run_id)
        reference = self.igv_reference_sizes()
        genome = str(reference) if reference.is_file() else 'hg38'
        loci = view["loci"]
        session = Element('Session', genome=genome, locus=loci[0] if loci else 'All', version='8')
        resources = SubElement(session, 'Resources')
        for track, path in zip(view["tracks"], paths):
            kind = "bw" if track["format"] == "bigwig" else "bed"
            SubElement(resources, "Resource", path=str(path), name=track["name"], type=kind)
        return tostring(session, encoding='utf-8', xml_declaration=True)

    def artifact(self, run_id, name, round_index=None):
        item = self.get_run(run_id)
        if item.get('paired_adaptive'):
            if name == FIGURE_NAME and round_index is None and item.get('rerun_figure'):
                return rerun_dir(self.root, run_id) / FIGURE_DIR / name, 'image/png'
            if round_index is not None or name not in PAIR_RERUN_ARTIFACTS:
                raise ValueError('文件不可下载')
            artifacts = item.get('track_review', {}).get('artifacts', [])
            meta = next((row for row in artifacts if row['name'] == name), None)
            path = rerun_dir(self.root, run_id) / name
            if not meta or path.resolve() != path or not path.is_file() or path.stat().st_size != meta['bytes']:
                raise ValueError('调参产物不存在或大小与记录不一致')
            return path, PAIR_RERUN_ARTIFACTS[name]
        if item.get("mode") == "preset":
            if round_index is not None or item["status"] != "completed_comparison" or name not in PRESET_ARTIFACTS:
                raise ValueError("文件不可下载")
            meta = next((row for row in item["artifacts"] if row["name"] == name), None)
            path = self.root / "outputs/preset_runs" / run_id / "tracks" / name
            if not meta or not path.is_file() or path.resolve() != path or path.stat().st_size != meta["bytes"]:
                raise ValueError("产物缺失或已改变")
            return path, PRESET_ARTIFACTS[name]
        if item.get("mode") == "spikein":
            if round_index is not None or item.get("status") != "completed_spikein" or name not in SPIKEIN_ARTIFACTS:
                raise ValueError("文件不可下载")
            meta = next((a for a in item["artifacts"] if a["name"] == name), None)
            path = self.root / "outputs/spikein_runs" / run_id / name
            if not meta or not path.is_file() or path.resolve() != path or path.stat().st_size != meta["bytes"]:
                raise ValueError("产物缺失或已改变")
            return path, SPIKEIN_ARTIFACTS[name]
        if item.get("mode") == "pipeline":
            if round_index is not None or name not in PIPELINE_ARTIFACTS:
                raise ValueError("文件不可下载")
            allowed = {row["name"]: row for row in item.get("artifacts", [])}
            if name not in allowed or (item["status"] != "completed_pipeline" and not item.get("baseline_ready") and name not in PIPELINE_LOGS):
                raise ValueError("产物尚未通过检查")
            folder = self.pipeline.output_dir(run_id)
            path = (folder / name).resolve()
            if path.parent != folder or not path.is_file() or (name not in PIPELINE_LOGS and path.stat().st_size != allowed[name]["bytes"]):
                raise ValueError("产物缺失或已改变")
            return path, PIPELINE_ARTIFACTS[name]
        if item.get("mode") not in ("coverage", "agent") or name not in ARTIFACTS:
            raise ValueError("文件不可下载")
        folder = self.coverage.output_dir(run_id)
        status = item["status"]
        if item["mode"] == "agent":
            if type(round_index) is not int or round_index not in (0, 1):
                raise ValueError("需要有效的工具轮次")
            round_item = next((r for r in item.get("rounds", []) if r["index"] == round_index), None)
            if not round_item:
                raise ValueError("该轮结果尚不可用")
            status = round_item["status"]
            folder = folder / ("round-" + str(round_index))
        elif round_index is not None:
            raise ValueError("此任务没有多轮产物")
        if name in ("coverage.bw", "result.json") and status != "completed":
            raise ValueError("覆盖轨道尚未通过检查")
        path = (folder / name).resolve()
        if path.parent != folder or not path.is_file() or path.stat().st_size > 8 * 1024 * 1024:
            raise ValueError("文件尚不可用")
        return path, ARTIFACTS[name]

    def _drive(self, run_id):
        for index in range(len(STEPS)):
            with self.lock:
                item = self.runs[run_id]
                if self._stop_requested(item):
                    return
                item["status"] = "running"
                item["steps"][index]["status"] = "running"
                item["updated_at"] = now()
                self.save(item)
            time.sleep(self.step_delay)
            with self.lock:
                item = self.runs[run_id]
                if self._stop_requested(item):
                    return
                item["steps"][index]["status"] = "done" if index in (0, 1, 4) else "simulated"
                item["progress"] = round((index + 1) / len(STEPS) * 100)
                item["updated_at"] = now()
                item["events"].append({"at": now(), "message": item["steps"][index]["label"] + "：" + item["steps"][index]["detail"]})
                if index == len(STEPS) - 1:
                    item["status"] = "completed_demo"
                    item["message"] = "网页流程演示完成；未调用模型、未执行分析工具，也未生成新的 BW。"
                self.save(item)

    def _stop_requested(self, item):
        if item["id"] not in self.cancelled:
            return False
        item["status"] = "cancelled"
        item["message"] = "演示已停止，原始输入未改变。"
        item["updated_at"] = now()
        for step in item["steps"]:
            if step["status"] == "running":
                step["status"] = "cancelled"
        item["events"].append({"at": now(), "message": item["message"]})
        self.save(item)
        return True

    def cancel_run(self, run_id):
        with self.lock:
            item = self.get_run(run_id)
            if item["status"] in ("queued", "running"):
                self.cancelled.add(run_id)
        return self.get_run(run_id)


def handler_for(workbench, static_root):
    class Handler(BaseHTTPRequestHandler):
        def send_file(self, path, mime, filename=None):
            size, start, end = path.stat().st_size, 0, path.stat().st_size - 1
            partial = self.headers.get("Range")
            if partial:
                match = re.fullmatch(r"bytes=(\d*)-(\d*)", partial)
                if not match or (not match[1] and not match[2]):
                    return self.send(416, {"error": "不支持该文件范围"})
                if match[1]:
                    start = int(match[1])
                    end = min(int(match[2]), size - 1) if match[2] else size - 1
                else:
                    start = max(0, size - int(match[2]))
                if start >= size or start > end or size == 0:
                    self.send_response(416)
                    self.send_header("Content-Range", "bytes */" + str(size))
                    self.send_header("Content-Length", "0")
                    self.end_headers()
                    return
            length = max(0, end - start + 1)
            self.send_response(206 if partial else 200)
            self.send_header("Content-Type", mime)
            self.send_header("Content-Length", str(length))
            self.send_header("Accept-Ranges", "bytes")
            self.send_header("Cache-Control", "no-store")
            self.send_header("X-Content-Type-Options", "nosniff")
            self.send_header("Content-Disposition", 'attachment; filename="' + (filename or path.name) + '"')
            if partial:
                self.send_header("Content-Range", "bytes %d-%d/%d" % (start, end, size))
            self.end_headers()
            if self.command == "HEAD":
                return
            try:
                with path.open("rb") as stream:
                    stream.seek(start)
                    while length:
                        data = stream.read(min(length, 1024 * 1024))
                        if not data:
                            break
                        self.wfile.write(data)
                        length -= len(data)
            except (BrokenPipeError, ConnectionResetError):
                pass

        def do_HEAD(self):
            self.do_GET()
        def log_message(self, *_args):
            pass

        def send(self, status, body, content_type="application/json; charset=utf-8", filename=None, *, igv_view=False):
            if isinstance(body, (dict, list)):
                body = json.dumps(body, ensure_ascii=False).encode("utf-8")
            elif isinstance(body, str):
                body = body.encode("utf-8")
            self.send_response(status)
            self.send_header("Content-Type", content_type)
            self.send_header("Content-Length", str(len(body)))
            self.send_header("Cache-Control", "no-store")
            self.send_header("X-Content-Type-Options", "nosniff")
            self.send_header("Referrer-Policy", "no-referrer")
            styles = "'self' 'unsafe-inline'" if igv_view else "'self'"
            images = "'self' data: blob:" if igv_view else "'self' data:"
            self.send_header("Content-Security-Policy", f"default-src 'self'; script-src 'self'; style-src {styles}; img-src {images}; connect-src 'self'; base-uri 'none'")
            if filename:
                self.send_header("Content-Disposition", 'attachment; filename="' + filename + '"')
            self.end_headers()
            if self.command != "HEAD":
                self.wfile.write(body)

        def do_GET(self):
            path = unquote(urlsplit(self.path).path)
            try:
                if path == "/api/health":
                    return self.send(200, {"status": "ok", **workbench.setup()})
                if path == "/api/samples":
                    return self.send(200, {"samples": workbench.samples()})
                if path == f"/api/samples/{PAIR_TARGET}/pair-evidence":
                    sample = workbench.paired_sample()
                    return self.send(200, {"sample_id": sample["id"], **sample["pair"]},
                                     filename=f"{PAIR_TARGET}.pair-evidence.json")
                if path == f"/api/samples/{PAIR_TARGET}/review-summary":
                    pair = workbench.hospital_pair.sample()["pair"]
                    summary = pair_review.make_payload(pair)
                    comparison = {**summary, "comparison": pair_review.compare_pair(pair, {"action": "compare_igg"})}
                    pair_review.validate_payload(comparison)
                    adaptive = pair_review.adaptive_payload(pair)
                    pair_review.validate_payload(adaptive, adaptive=True)
                    return self.send(200, {"initial": summary, "if_compare_igg": comparison, "adaptive": adaptive},
                                     filename="paired-qc-model-summary.json")
                paired = re.fullmatch(r"/api/samples/" + re.escape(PAIR_TARGET) + r"/artifacts/([a-zA-Z0-9_.-]+)", path)
                if paired:
                    file, mime = workbench.hospital_pair.artifact(paired[1])
                    return self.send_file(file, mime, filename=paired[1])
                if path == "/api/runs":
                    return self.send(200, {"runs": workbench.list_runs(), "archived_run_ids": workbench.archived_run_ids()})
                if path == "/api/igv/hg38.chrom.sizes":
                    reference = workbench.igv_reference_sizes()
                    if not reference.is_file() or reference.resolve() != reference:
                        raise ValueError("本地 hg38 染色体长度表不可用")
                    return self.send(200, reference.read_bytes(), "text/plain; charset=utf-8")
                igv_run = re.fullmatch(r"/api/runs/(run-[0-9a-f]{12})/igv-(session|view)", path)
                if igv_run:
                    run_id, action = igv_run.groups()
                    if action == "session":
                        return self.send(200, workbench.igv_session(run_id), "application/xml; charset=utf-8",
                                         filename=run_id + "-igv.xml")
                    view, _ = workbench.igv_inputs(run_id)
                    return self.send(200, view)
                artifact = re.fullmatch(r"/api/runs/(run-[0-9a-f]{12})/artifacts/([a-zA-Z0-9_.-]+)", path)
                if artifact:
                    file, mime = workbench.artifact(artifact[1], artifact[2])
                    return self.send_file(file, mime)
                artifact = re.fullmatch(r"/api/runs/(run-[0-9a-f]{12})/rounds/([01])/artifacts/([a-zA-Z0-9_.-]+)", path)
                if artifact:
                    file, mime = workbench.artifact(artifact[1], artifact[3], int(artifact[2]))
                    return self.send_file(file, mime)
                match = re.fullmatch(r"/api/runs/(run-[0-9a-f]{12})(/report)?", path)
                if match:
                    return self.send(200, workbench.get_run(match[1]), filename=match[1] + ".json" if match[2] else None)
                assets = {
                    "/": ("index.html", "text/html"), "/index.html": ("index.html", "text/html"),
                    "/app.js": ("app.js", "text/javascript"), "/styles.css": ("styles.css", "text/css"),
                    "/favicon.svg": ("favicon.svg", "image/svg+xml"),
                    "/igv-view.html": ("igv-view.html", "text/html"), "/igv-view.js": ("igv-view.js", "text/javascript"),
                    "/igv-view.css": ("igv-view.css", "text/css"),
                    "/vendor/igv-3.8.9.esm.min.js": ("vendor/igv-3.8.9.esm.min.js", "text/javascript"),
                    "/vendor/igv-3.8.9-LICENSE.txt": ("vendor/igv-3.8.9-LICENSE.txt", "text/plain"),
                }
                if path in assets:
                    filename, mime = assets[path]
                    return self.send(200, (static_root / filename).read_bytes(), mime + "; charset=utf-8",
                                     igv_view=path == "/igv-view.html")
                self.send(404, {"error": "页面不存在"})
            except ValueError as exc:
                self.send(404, {"error": str(exc)})
            except (OSError, KeyError, TypeError):
                self.send(500, {"error": "本地文件读取失败，请检查数据准备状态"})

        def do_POST(self):
            host = self.headers.get("Host", "")
            allowed_origins = {"http://127.0.0.1:" + str(self.server.server_port), "http://localhost:" + str(self.server.server_port)}
            if "http://" + host not in allowed_origins or self.headers.get("Origin", "http://" + host) not in allowed_origins:
                return self.send(403, {"error": "仅接受本地页面的请求"})
            if not self.headers.get("Content-Type", "").startswith("application/json"):
                return self.send(415, {"error": "需要 JSON 请求"})
            try:
                length = int(self.headers.get("Content-Length", "0"))
                if not 0 < length <= 65536:
                    return self.send(413, {"error": "请求大小不受支持"})
                payload = json.loads(self.rfile.read(length))
                path = urlsplit(self.path).path
                if path == "/api/runs":
                    return self.send(201, workbench.create_run(payload))
                if path == "/api/agent-chat":
                    return self.send(200, workbench.ask_sample_chat(payload))
                if path == "/api/model/test":
                    if payload != {}:
                        raise ValueError("模型测试只接受固定的合成内容")
                    return self.send(200, workbench.model.test())
                match = re.fullmatch(r"/api/runs/(run-[0-9a-f]{12})/cancel", path)
                if match:
                    return self.send(200, workbench.cancel_run(match[1]))
                match = re.fullmatch(r"/api/runs/(run-[0-9a-f]{12})/resume-paired-rerun", path)
                if match and payload == {}:
                    return self.send(200, workbench.resume_paired_rerun(match[1]))
                match = re.fullmatch(r"/api/runs/(run-[0-9a-f]{12})/decide-paired-rerun", path)
                if match and isinstance(payload, dict) and set(payload) == {"accept"}:
                    return self.send(200, workbench.decide_paired_rerun(match[1], payload["accept"]))
                match = re.fullmatch(r"/api/runs/(run-[0-9a-f]{12})/replay-paired-rerun", path)
                if match and payload == {}:
                    return self.send(201, workbench.replay_paired_rerun(match[1]))
                match = re.fullmatch(r"/api/runs/(run-[0-9a-f]{12})/archive", path)
                if match and isinstance(payload, dict) and set(payload) == {"archived"}:
                    return self.send(200, workbench.set_run_archived(match[1], payload['archived']))
                self.send(404, {"error": "操作不存在"})
            except (ValueError, TypeError) as exc:
                self.send(400, {"error": str(exc) or "任务参数无效"})

    return Handler


def main():
    parser = argparse.ArgumentParser(description="本地 CUT&Tag 网页原型")
    parser.add_argument("--port", type=int, default=8766)
    args = parser.parse_args()
    state = Workbench()
    server = ThreadingHTTPServer(("127.0.0.1", args.port), handler_for(state, PROJECT_ROOT / "web" / "dist"))
    print("CUT&Tag workbench · local preview", flush=True)
    print("Local: http://127.0.0.1:" + str(server.server_port), flush=True)
    print("本地工具、模型测试与受限 Agent 闭环。模型状态见接入设置。Ctrl+C 停止。", flush=True)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        for item in state.list_runs():
            if item.get("mode") in REAL_MODES and item["status"] in ("queued", "running"):
                state.cancelled.add(item["id"])
                if item["mode"] in CONTAINER_MODES and not item.get("baseline_ready"):
                    (state.pipeline if item["mode"] == "pipeline" else state.coverage).cancel(item["id"])
                if item["mode"] == "preset":
                    state.cancel_preset_container(item)
                if item.get("paired_adaptive") and (item.get("evidence_preparing") or item.get("track_review", {}).get("rerun_started")):
                    stop_container('cuttag-pair-' + item['id'][4:])
        server.server_close()


if __name__ == "__main__":
    main()
