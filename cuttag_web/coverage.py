"""Local, bounded bamCoverage execution. Never uses SSH or a model service."""
import hashlib
import json
import math
import os
import re
import signal
import subprocess
import time
from pathlib import Path


IMAGE = "cuttag-deeptools:3.5.6"
REGION = "chr1:1000000:2000000"
ARTIFACTS = {"coverage.bw": "application/octet-stream", "tool.log": "text/plain; charset=utf-8",
             "result.json": "application/json", "requirements.lock": "text/plain; charset=utf-8",
             "bam-stats.json": "application/json", "samtools-quickcheck.txt": "text/plain; charset=utf-8", "samtools-idxstats.tsv": "text/tab-separated-values"}


def read_record(path):
    try:
        if path.stat().st_size > 1024 * 1024:
            return {}
        return json.loads(path.read_text())
    except (OSError, ValueError):
        return {}


class DockerCoverageRunner:
    def __init__(self, root, timeout=600):
        self.root = Path(root).resolve()
        self.timeout = timeout
        self._probe = None
        self._checked = 0

    def probe(self):
        if self._probe is not None and time.monotonic() - self._checked < 15:
            return dict(self._probe)
        info = {"ready": False, "tool": "bamCoverage", "region": REGION, "bin_size": 25}
        try:
            result = subprocess.run(["docker", "image", "inspect", IMAGE], capture_output=True, text=True, timeout=5)
            if result.returncode:
                raise ValueError("本地工具容器尚未准备好")
            metadata = json.loads(result.stdout)[0]
            image_id = metadata["Id"]
            if self._probe and self._probe.get("image_id") == image_id:
                info = dict(self._probe)
            else:
                version = subprocess.run(["docker", "run", "--rm", "--network", "none", "--read-only",
                    "--memory", "512m", "--tmpfs", "/tmp:rw,size=64m", "--entrypoint", "bamCoverage",
                    image_id, "--version"], capture_output=True, text=True, timeout=20)
                if version.returncode or version.stdout.strip() != "bamCoverage 3.5.6":
                    raise ValueError("bamCoverage 版本检查未通过")
                info.update(ready=True, image_id=image_id, architecture=metadata["Architecture"], version="3.5.6")
                stats = subprocess.run(["docker", "run", "--rm", "--network", "none", "--read-only", "--entrypoint", "samtools", image_id, "--version"], capture_output=True, text=True, errors="replace", timeout=10)
                info.update(samtools_ready=stats.returncode == 0, samtools_version=stats.stdout.splitlines()[0] if stats.returncode == 0 else None)
        except (OSError, ValueError, KeyError, IndexError, subprocess.TimeoutExpired):
            info["reason"] = "请先准备本地 bamCoverage 容器并确认 Docker 正在运行"
        self._probe, self._checked = info, time.monotonic()
        return dict(info)

    @staticmethod
    def container_name(run_id):
        if not re.fullmatch(r"run-[0-9a-f]{12}", run_id):
            raise ValueError("运行标识无效")
        return "cuttag-" + run_id

    def output_dir(self, run_id):
        self.container_name(run_id)
        base = (self.root / "outputs" / "tool_runs").resolve()
        path = (base / run_id).resolve()
        if self.root not in base.parents or path.parent != base:
            raise ValueError("输出目录越界")
        return path

    def input_bam(self, sample_id):
        if not re.fullmatch(r"[A-Za-z0-9_.-]{1,120}", sample_id):
            raise ValueError("样本名称无效")
        sample = (self.root / "data" / "private" / sample_id).resolve()
        bam = (sample / "processed" / (sample_id + ".hg38.dedup.bam")).resolve()
        if sample.parent != (self.root / "data" / "private").resolve() or sample not in bam.parents:
            raise ValueError("输入目录越界")
        index = Path(str(bam) + ".bai")
        if not bam.is_file() or not index.is_file() or index.resolve().parent != bam.parent:
            raise ValueError("缺少本地去重 BAM 或 BAI 索引")
        return bam

    def command(self, run_id, bam, folder, image_id, bin_size=25):
        if type(bin_size) is not int or bin_size not in (25, 50):
            raise ValueError("覆盖窗口只允许 25 或 50")
        return ["docker", "run", "--rm", "--name", self.container_name(run_id),
                "--network", "none", "--read-only", "--cap-drop", "ALL", "--security-opt", "no-new-privileges",
                "--cpus", "2", "--memory", "2g", "--pids-limit", "128", "--user", "1000:1000",
                "--tmpfs", "/tmp:rw,nosuid,size=256m,mode=1777",
                "--mount", "type=bind,src=" + str(bam.parent) + ",dst=/input,readonly",
                "--mount", "type=bind,src=" + str(folder) + ",dst=/output",
                image_id, bam.name, str(bin_size)]

    def cancel(self, run_id):
        name = self.container_name(run_id)
        for command in (["docker", "stop", "--time", "2", name], ["docker", "rm", "-f", name]):
            try:
                subprocess.run(command, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=8)
            except (OSError, subprocess.TimeoutExpired):
                pass

    @staticmethod
    def validate_result(folder, bin_size=25):
        result = read_record(folder / "result.json")
        track = folder / "coverage.bw"
        if result.get("validation") != "passed" or (result.get("chrom"), result.get("start"), result.get("end"), result.get("bin_size")) != ("chr1", 1000000, 2000000, bin_size):
            raise ValueError("真实产物校验报告缺失或范围不符")
        with track.open("rb") as stream:
            content = stream.read(8 * 1024 * 1024 + 1)
        if len(content) > 8 * 1024 * 1024 or len(content) < 64 or content[:4] != b"\x26\xfc\x8f\x88":
            raise ValueError("bigWig 文件格式或大小检查失败")
        if hashlib.sha256(content).hexdigest() != result.get("sha256") or len(content) != result.get("bytes"):
            raise ValueError("bigWig 与校验记录不匹配")
        points = result.get("preview", {}).get("points", [])
        if len(points) != 200 or any(value is not None and (not isinstance(value, (int, float)) or not math.isfinite(value)) for value in points):
            raise ValueError("覆盖轨道预览数据无效")
        return result

    def run(self, run_id, sample_id, should_cancel, phase, bin_size=25, round_index=None):
        runtime = self.probe()
        if not runtime["ready"]:
            raise ValueError(runtime["reason"])
        bam = self.input_bam(sample_id)
        sources = [bam, Path(str(bam) + ".bai")]
        snapshot = {str(path.relative_to(self.root)): {"bytes": path.stat().st_size, "mtime_ns": path.stat().st_mtime_ns} for path in sources}
        folder = self.output_dir(run_id)
        if round_index is not None:
            if type(round_index) is not int or round_index not in (0, 1):
                raise ValueError("最多运行基线与一次重跑")
            folder = folder / ("round-" + str(round_index))
        folder.mkdir(parents=True, exist_ok=False)
        command = self.command(run_id, bam, folder, runtime["image_id"], bin_size)
        execution = {"runtime": runtime, "argv": command, "inputs": snapshot, "timeout_seconds": self.timeout,
                     "limits": {"cpu": 2, "memory_gb": 2, "network": "none", "input_mount": "readonly"},
                     "output_directory": str(folder.relative_to(self.root)), "exit_code": None}
        phase("inputs", execution)
        if should_cancel():
            return {"status": "cancelled", "execution": execution, "actual_tool_executions": 0}
        status, process = "failed", None
        started_at = time.monotonic()
        try:
            with (folder / "tool.log").open("wb") as log:
                process = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT, start_new_session=True)
                phase("coverage", execution)
                while process.poll() is None:
                    if should_cancel() or time.monotonic() - started_at > self.timeout:
                        status = "cancelled" if should_cancel() else "timed_out"
                        self.cancel(run_id)
                        try:
                            process.wait(timeout=3)
                        except subprocess.TimeoutExpired:
                            os.killpg(process.pid, signal.SIGTERM)
                            try:
                                process.wait(timeout=2)
                            except subprocess.TimeoutExpired:
                                os.killpg(process.pid, signal.SIGKILL)
                                process.wait(timeout=2)
                        break
                    time.sleep(.1)
                execution["exit_code"] = process.returncode
            if status not in ("cancelled", "timed_out"):
                if should_cancel():
                    status = "cancelled"
                elif process.returncode == 0:
                    phase("validate", execution)
                    result = self.validate_result(folder, bin_size)
                    stats = read_record(folder / "bam-stats.json")
                    if stats.get("quickcheck") != "passed" or any(type(stats.get(k)) is not int or stats[k] < 0 for k in ("mapped_segments", "unmapped_segments", "reference_count")) or stats["reference_count"] <= 0:
                        raise ValueError("SAMtools 统计记录缺失或无效")
                    status = "completed"
                else:
                    execution["error"] = "工具未成功退出，请查看运行日志"
        except (OSError, ValueError, subprocess.SubprocessError) as exc:
            execution["error"] = str(exc)
        finally:
            if status != "completed":
                self.cancel(run_id)
            if process is not None and process.poll() is None:
                self.cancel(run_id)
                process.kill()
                process.wait(timeout=5)
        execution["elapsed_seconds"] = round(time.monotonic() - started_at, 3)
        execution["tool_start"] = read_record(folder / "tool-start.json")
        execution["tool_exit"] = read_record(folder / "tool-exit.json")
        invocations = read_record(folder / "stat-invocations.json")
        execution["stat_invocations"] = invocations if isinstance(invocations, list) else []
        execution["tool_counts"] = {"samtools": len(execution["stat_invocations"]), "bamCoverage": int(bool(execution["tool_start"]))}
        execution["inputs_unchanged"] = all(path.stat().st_size == snapshot[str(path.relative_to(self.root))]["bytes"] and path.stat().st_mtime_ns == snapshot[str(path.relative_to(self.root))]["mtime_ns"] for path in sources)
        if not execution["inputs_unchanged"]:
            status = "failed"
            execution["error"] = "输入元数据发生变化，不能确认输入保持不变"
        output = {"status": status, "execution": execution, "actual_tool_executions": sum(execution["tool_counts"].values())}
        if status == "completed":
            output["result"] = result
            output["bam_stats"] = read_record(folder / "bam-stats.json")
        output["artifacts"] = [{"name": name, "bytes": (folder / name).stat().st_size} for name in ARTIFACTS
                               if (folder / name).is_file() and (status == "completed" or name not in ("coverage.bw", "result.json"))]
        return output
