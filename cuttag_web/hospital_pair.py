"""Verified local evidence for the registered target/IgG analysis."""
import hashlib
import json
from pathlib import Path

from .hospital_pair_review import validate_igv
from .pipeline import MIMES


TARGET = "paired-target"
CONTROL = "paired-igg"
CASE_PATH = Path("outputs/hospital_pair") / TARGET
SOURCE_PATH = Path("data/support/paired-case/inputs")
SAMPLES = {"target": TARGET, "igg": CONTROL}
ARTIFACTS = {
    "target-cpm.bw": ("target", "coverage.cpm.bw"),
    "igg-cpm.bw": ("igg", "coverage.cpm.bw"),
    "baseline.peaks.bed": ("peaks", "baseline.peaks.bed"),
    "with-igg.peaks.bed": ("peaks", "with-igg.peaks.bed"),
    "target-tss-profile.png": ("target", "tss.profile.png"),
    "igg-tss-profile.png": ("igg", "tss.profile.png"),
    "target-tss-heatmap.png": ("target", "tss.heatmap.png"),
    "igg-tss-heatmap.png": ("igg", "tss.heatmap.png"),
}
SCREEN_ARTIFACTS = ("baseline.peaks.bed", "with-igg.peaks.bed", "target-cpm.bw", "igg-cpm.bw")
UNKNOWN_LABELS = {
    "exact_mark_and_antibody": "目标的确切修饰与抗体",
    "biological_replicates": "生物学重复关系",
    "biological_control_signoff": "实验记录中的对照设计",
    "manual_igv_review": "人工 IGV 复核结论",
}


def local_file(path):
    if path.resolve() != path or not path.is_file():
        raise ValueError("配对样本文件缺失或路径无效")
    return path


def read(path):
    with local_file(path).open(encoding="utf-8") as stream:
        return json.load(stream)


class HospitalPair:
    def __init__(self, root):
        self.root = root
        self.folder = root / CASE_PATH

    def evidence(self):
        source = read(self.folder / "source-materials.json")
        computed = read(self.folder / "full/agent-evidence.json")
        spike = read(self.folder / "spikein-context.json")
        confirmed = spike["actual_experiment_report"]
        if ({role: sample["sample_id"] for role, sample in source["samples"].items()} != SAMPLES
                or set(spike["samples"]) != set(SAMPLES.values()) or confirmed["applies_to_both_samples"] is not True
                or computed["status"] != "tool_artifacts_ready" or computed["scope"] != "whole_genome"):
            raise ValueError("本地配对记录与已完成的计算不匹配")
        return source, computed, spike

    def artifact(self, name, source=None):
        if name in ARTIFACTS:
            stage, filename = ARTIFACTS[name]
            folder = self.folder / "full" / stage
            size = read(folder / "stage.json")["outputs"].get(filename)
            path = folder / filename
        elif name in ("target-spikein.bw", "igg-spikein.bw"):
            role = name.split("-", 1)[0]
            sample = SAMPLES[role]
            filename = f"{sample}.spikein.bw"
            relative = SOURCE_PATH / sample / "bw" / filename
            source = source if source is not None else self.evidence()[0]
            size = source["samples"][role]["materials"][str(relative)]["bytes"]
            path = self.root / relative
        else:
            raise ValueError("没有这个配对样本产物")
        if local_file(path).stat().st_size != size:
            raise ValueError("配对样本产物大小与登记结果不一致")
        mime_key = "narrow.peaks.bed" if path.suffix == ".bed" else filename
        return path, MIMES.get(mime_key, "application/octet-stream")

    def screen_recipe(self):
        mounts, identities = [], {}
        for name in SCREEN_ARTIFACTS:
            path = self.artifact(name)[0]
            stat = path.stat()
            identities[name] = {"bytes": stat.st_size, "mtime_ns": stat.st_mtime_ns}
            mounts.append((path, "/baseline/" + path.relative_to(self.folder / "full").as_posix()))
        script = local_file(self.root / "tools/inspect_paired_igv.py")
        image = read(self.folder / "full/peaks/stage.json")["inputs"]["image"]
        receipt = {"files": identities, "screen_code": hashlib.sha256(script.read_bytes()).hexdigest(), "image": image}
        key = hashlib.sha256(json.dumps(receipt, sort_keys=True).encode()).hexdigest()[:24]
        folder = self.folder / "automatic-peak-screen" / key
        if folder.resolve() != folder:
            raise ValueError("轨道统计目录越界")
        return mounts, receipt, folder

    def igv_screen(self, peak_counts, recipe):
        _, receipt, folder = recipe
        manual = None
        record_path = self.folder / "igv-review/agent-igv-screen.json"
        if record_path.is_file():
            record = read(record_path)
            if record["status"] == "exploratory_reviewed" and record["inputs"] == receipt["files"]:
                validate_igv(record["summary"], peak_counts["baseline"])
                manual = record["summary"]
        if not folder.is_dir():
            return manual
        report_name = "interval-summary.json"
        state = read(folder / "stage.json")
        if state["inputs"] != receipt or local_file(folder / report_name).stat().st_size != state["outputs"][report_name]:
            raise ValueError("自动轨道统计与输入记录不一致")
        report = read(folder / report_name)
        if (report["baseline_peak_count"], report["controlled_peak_count"]) != (peak_counts["baseline"], peak_counts["with-igg"]):
            raise ValueError("自动轨道统计与已登记的峰表数量不一致")

        def group(name):
            data = report["groups"][name]
            return {"total": data["total"], **data["counts_by_fold"]["2"], "target_local_peak": data["target_higher_focal_2x"],
                    **{f"median_{role}_cpm": round(data[f"median_{role}_cpm"], 6) if data[f"median_{role}_cpm"] is not None else None
                       for role in SAMPLES}}

        screen = {"ratio_threshold": 2, "igv_examples_reviewed": 0, "removed": group("removed"),
                  "overlaps_controlled": group("overlap_igg_controlled")}
        validate_igv(screen, peak_counts["baseline"])
        if manual and {**manual, "igv_examples_reviewed": 0} == screen:
            return manual
        return screen

    def prepare_screen(self, pair, run_id, cancelled):
        from tools.prepare_hospital_pair import stage

        mounts, receipt, folder = self.screen_recipe()
        folder.parent.mkdir(parents=True, exist_ok=True)
        stage(folder, ["screen_igv"], mounts, ("interval-summary.json",), receipt, cancelled=cancelled,
              container="cuttag-pair-" + run_id[4:], extra_files=("tools/inspect_paired_igv.py",), timeout=120)
        current = self.screen_recipe()
        if current[1] != receipt:
            raise ValueError("自动统计过程中轨道或峰表发生变化")
        result = self.igv_screen(pair["peak_counts"], current)
        if result is None:
            raise ValueError("自动轨道统计缺失，无法完成证据整理")
        return result

    def sample(self):
        source, computed, spike = self.evidence()
        comparison = computed["peak_comparison"]
        confirmed = spike["actual_experiment_report"]
        if confirmed["addition_stage"] != "after_fragmentation":
            raise ValueError("配对样本的加样时点与已登记记录不一致")
        names = (*ARTIFACTS, "target-spikein.bw", "igg-spikein.bw")
        for name in names:
            self.artifact(name, source)
        files = []
        for role, sample in source["samples"].items():
            for name, record in sample["materials"].items():
                relative = Path(name)
                if relative.parent.parent != SOURCE_PATH / SAMPLES[role]:
                    raise ValueError("历史材料超出配对样本目录")
                path = self.root / relative
                present = path.is_file() and path.resolve() == path
                size = path.stat().st_size if present else 0
                files.append({"name": relative.name, "group": "目标" if role == "target" else "IgG",
                              "bytes": size, "present": present, "size_matches": present and size == record["bytes"]})
        windows_path = self.folder / "web-integration-preview.json"
        windows = []
        if windows_path.is_file():
            preview = read(windows_path)
            if preview["pair"] == SAMPLES and preview["historical_execution"]["peak_counts"] == comparison["peak_counts"]:
                windows = preview["manual_review_windows"]
        return {
            "id": TARGET,
            "label": TARGET,
            "input_kind": "两份已过滤 BAM · 全基因组下游分析已完成",
            "report_available": False,
            "files": files,
            "total_bytes": sum(item["bytes"] for item in files),
            "pair": {
                "control_id": CONTROL,
                "control_source": "本地登记的匹配对照记录",
                "artifacts": names,
                "peak_settings": comparison["settings"],
                "peak_counts": comparison["peak_counts"],
                "fragment_overlap": {role: comparison["fragment_metrics"][role]["frip"] for role in SAMPLES},
                "fragments": {role: computed["track_summaries"][role]["fragments"] for role in SAMPLES},
                "bin_size": computed["track_summaries"]["target"]["coverage"]["bin_size"],
                "spikein": {
                    "addition_stage": confirmed["addition_stage"],
                    "dose_pg": confirmed["nominal_spikein_pg"],
                    "cells_per_dose": confirmed["nominal_cells_per_dose"],
                    "actual_cell_counts": confirmed["actual_input_cell_count_per_sample"],
                    "report_date": spike["date"],
                    "historical_fragments": {role: int(source["samples"][role]["historical_spikein"]["spikein_fragments"])
                                             for role in SAMPLES},
                    "quantitative_review_pending": spike["quantitative_comparability"] == "pending_review",
                },
                "manual_windows": windows,
                "igv_screen": self.igv_screen(comparison["peak_counts"], self.screen_recipe()),
                "unconfirmed": [UNKNOWN_LABELS[key] for key in computed["unconfirmed"] if key in UNKNOWN_LABELS],
            },
        }
