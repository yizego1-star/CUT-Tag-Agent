"""Read-only local inputs -> isolated basic CUT&Tag outputs. No SSH or model export."""
import hashlib
import json
import math
import re
import shutil
import subprocess
import time
from pathlib import Path

from .coverage import DockerCoverageRunner, read_record

IMAGE = 'cuttag-pipeline:2026-11-17'
RESUME_IMAGE = 'cuttag-pipeline-resume:2026-11-18'
DISK_RESERVE = 50 * 10 ** 9
RESUME_START_FREE = 60 * 10 ** 9
STAGES = [('inputs', '校验输入与参考'), ('fastp', 'FASTQ 清理与质控'), ('align', '比对 hg38'),
          ('bam', '配对过滤与 BAM 整理'), ('coverage', '全基因组 CPM 轨道'), ('peaks', '候选 Peak BED'),
          ('qc', '片段质控与 FRiP'), ('deeptools', 'TSS 矩阵与图表'), ('validate', '校验并保存产物')]
MIMES = {'coverage.cpm.bw': 'application/octet-stream', 'narrow.peaks.bed': 'text/tab-separated-values',
         'broad.peaks.bed': 'text/tab-separated-values', 'tss.matrix.gz': 'application/gzip',
         'tss.regions.bed': 'text/tab-separated-values', 'tss.profile.png': 'image/png', 'tss.heatmap.png': 'image/png',
         'tss.profile.tsv': 'text/tab-separated-values', 'fragment-length.png': 'image/png',
         'fragment-length.tsv': 'text/tab-separated-values', 'fragment-qc.json': 'application/json',
         'summary.json': 'application/json', 'report.md': 'text/plain; charset=utf-8',
         'analysis.flagstat.json': 'application/json', 'idxstats.tsv': 'text/tab-separated-values',
         'quickcheck.txt': 'text/plain; charset=utf-8', 'versions.json': 'application/json',
         'commands.json': 'application/json', 'fastp.json': 'application/json', 'fastp.html': 'text/plain; charset=utf-8',
         'bowtie2.log': 'text/plain; charset=utf-8', 'alignment.flagstat.json': 'application/json',
         'pair-filter.json': 'application/json', 'markdup.txt': 'text/plain; charset=utf-8',
         'tool.log': 'text/plain; charset=utf-8', 'artifacts.json': 'application/json'}


def sha256(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for data in iter(lambda: stream.read(4 * 1024 * 1024), b''):
            digest.update(data)
    return digest.hexdigest()


class PipelineRunner:
    def __init__(self, root, timeout=7200, *, image=IMAGE):
        self.root = Path(root).resolve()
        self.timeout = timeout
        self.image = image
        self._runtime = None
        self._checked = 0

    @staticmethod
    def container_name(run_id):
        DockerCoverageRunner.container_name(run_id)
        return 'cuttag-pipeline-' + run_id

    def output_dir(self, run_id):
        self.container_name(run_id)
        base = (self.root / 'outputs' / 'pipeline_runs').resolve()
        folder = (base / run_id).resolve()
        if self.root not in base.parents or folder.parent != base:
            raise ValueError('基础分析输出路径越界')
        return folder

    def inputs(self, sample_id, source, *, checkpoint=None):
        if ((source not in ('fastq', 'bam') and not (source == 'resume' and checkpoint))
                or not isinstance(sample_id, str) or not re.fullmatch(r'[A-Za-z0-9_.-]{1,120}', sample_id)):
            raise ValueError('输入起点或样本无效')
        sample = (self.root / 'data/private' / sample_id).resolve()
        support = (self.root / 'data/support' / sample_id).resolve()
        if sample.parent != (self.root / 'data/private').resolve() or support.parent != (self.root / 'data/support').resolve():
            raise ValueError('输入目录越界')
        compatible = read_record(support / 'provenance/reference_compatibility.json')
        if compatible.get('status') != 'passed':
            raise ValueError('需要通过匹配检查的本地参考包')
        ref_files = [support / 'reference/hg38/hg38.bam-compatible.chrom.sizes',
                     support / 'annotation/gencode.v44.genes.bam-compatible.bed']
        if source == 'fastq':
            files = [sample / 'raw' / (sample_id + '_' + str(i) + '.fq.gz') for i in (1, 2)]
            ref_files += [support / 'reference/hg38/bowtie2' / ('hg38_noalt.' + suffix + '.bt2')
                          for suffix in ('1', '2', '3', '4', 'rev.1', 'rev.2')]
        elif source == 'bam':
            files = [sample / 'processed' / (sample_id + '.hg38.dedup.bam' + suffix) for suffix in ('', '.bai')]
        else:
            files = []
        for path in files + ref_files:
            if not path.is_file() or path.stat().st_size == 0 or (sample not in path.resolve().parents and support not in path.resolve().parents):
                raise ValueError('基础流程缺少有效的本地输入或参考文件')
        return sample, support, files + ref_files + list(checkpoint or ())

    def resume_inputs(self, previous_id):
        folder = self.output_dir(previous_id)
        record = read_record(self.root / 'outputs' / 'web_runs' / (previous_id + '.json'))
        commands = read_record(folder / 'commands.json')
        if (folder.name != previous_id or not isinstance(record, dict) or record.get('id') != previous_id
                or record.get('mode') != 'pipeline' or record.get('source') != 'fastq' or record.get('status') != 'failed'
                or not isinstance(record.get('sample_id'), str) or not isinstance(commands, list)):
            raise ValueError('续跑只接受本地已停止且留有记录的 FASTQ 基础分析')
        completed = [entry['argv'] for entry in commands if isinstance(entry, dict) and entry.get('exit_code') == 0
                     and isinstance(entry.get('argv'), list)]
        if not {'index', 'quickcheck', 'idxstats'}.issubset(
                {argv[1] for argv in completed if len(argv) >= 3 and argv[0] == 'samtools' and '/output/analysis.bam' in argv}):
            raise ValueError('旧任务没有已完成的 BAM 检查记录')
        bam = folder / 'analysis.bam'
        pair = (bam, bam.with_suffix('.bam.bai'))
        if any(p.is_symlink() or not p.is_file() or p.stat().st_size == 0 for p in pair):
            raise ValueError('旧任务的 BAM 或索引缺失')
        return record['sample_id'], pair

    def probe(self):
        if self._runtime is not None and time.monotonic() - self._checked < 15:
            return dict(self._runtime)
        info = {'ready': False, 'source_options': ['fastq', 'bam']}
        try:
            p = subprocess.run(['docker', 'image', 'inspect', self.image], capture_output=True, text=True, timeout=5)
            if p.returncode:
                raise ValueError('本地基础分析镜像尚未准备好')
            metadata = json.loads(p.stdout)[0]
            if self._runtime and self._runtime.get('image_id') == metadata['Id']:
                info = dict(self._runtime)
            else:
                check = subprocess.run(['docker', 'run', '--rm', '--network', 'none', '--read-only', '--entrypoint', 'cat',
                                        metadata['Id'], '/opt/pipeline-versions.txt'], capture_output=True, text=True, timeout=10)
                if check.returncode or not all(tool in check.stdout for tool in ('fastp', 'bowtie2', 'macs', 'samtools')):
                    raise ValueError('基础工具版本探测失败')
                info.update(ready=True, image_id=metadata['Id'], architecture=metadata['Architecture'], packages=check.stdout.strip())
        except (OSError, ValueError, KeyError, IndexError, subprocess.TimeoutExpired):
            info['reason'] = '请准备本地基础分析镜像并启动 Docker'
        self._runtime, self._checked = info, time.monotonic()
        return dict(info)

    def command(self, run_id, sample, support, folder, image_id, source, *, checkpoint=None):
        inputs = ['--mount', 'type=bind,src=' + str(sample) + ',dst=/input,readonly']
        if checkpoint:
            inputs = [arg for p in checkpoint for arg in ('--mount', 'type=bind,src=' + str(p) + ',dst=/resume/' + p.name + ',readonly')]
        return ['docker', 'run', '--rm', '--pull', 'never', '--name', self.container_name(run_id), '--network', 'none', '--read-only',
                '--cap-drop', 'ALL', '--security-opt', 'no-new-privileges', '--cpus', '4', '--memory', '6g',
                '--pids-limit', '256', '--user', '1000:1000', '--tmpfs', '/tmp:rw,nosuid,size=256m,mode=1777',
                '--env', 'TMPDIR=/output/.scratch', *inputs,
                '--mount', 'type=bind,src=' + str(support) + ',dst=/support,readonly',
                '--mount', 'type=bind,src=' + str(folder) + ',dst=/output', image_id, source, sample.name]

    def cancel(self, run_id):
        name = self.container_name(run_id)
        for command in (['docker', 'stop', '--time', '2', name], ['docker', 'rm', '-f', name]):
            try:
                subprocess.run(command, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=8)
            except (OSError, subprocess.TimeoutExpired):
                pass

    @staticmethod
    def report_metrics(folder, source):
        if source != 'fastq':
            return {}
        fastp = read_record(folder / 'fastp.json').get('summary', {})
        before, after = fastp.get('before_filtering', {}), fastp.get('after_filtering', {})
        metrics = {}
        for key, value in {'raw_reads': before.get('total_reads'), 'retained_reads': after.get('total_reads'),
                           'raw_q30_rate': before.get('q30_rate')}.items():
            if type(value) in (int, float) and math.isfinite(value) and value >= 0:
                metrics[key] = value
        flagstat = read_record(folder / 'alignment.flagstat.json')
        groups = [flagstat.get(name, {}) for name in ('QC-passed reads', 'QC-failed reads')]
        primary = sum(g.get('primary', 0) for g in groups)
        mapped = sum(g.get('primary mapped', 0) for g in groups)
        if primary > 0 and 0 <= mapped <= primary:
            metrics.update(primary_alignments=primary, primary_mapped=mapped, primary_mapping_rate=mapped / primary)
        return metrics

    @staticmethod
    def validate(folder, source):
        summary = read_record(folder / 'summary.json')
        artifacts = read_record(folder / 'artifacts.json')
        if (summary.get('validation') != 'passed' or summary.get('source') != source
                or summary.get('scientific_quality') != 'not_assessed' or artifacts.get('validation') != 'passed'):
            raise ValueError('基础分析结果未通过检查')
        required = {'coverage.cpm.bw', 'narrow.peaks.bed', 'broad.peaks.bed', 'tss.matrix.gz',
                    'tss.profile.png', 'tss.heatmap.png', 'fragment-length.png', 'fragment-qc.json',
                    'summary.json', 'report.md', 'versions.json', 'commands.json'}
        if source == 'fastq':
            required |= {'fastp.json', 'alignment.flagstat.json', 'pair-filter.json'}
        rows = artifacts.get('artifacts', [])
        if not isinstance(rows, list) or not required.issubset({a.get('name') for a in rows}):
            raise ValueError('缺少应交付的基础分析产物')
        seen = set()
        for artifact in rows:
            name = artifact.get('name')
            if name not in MIMES or name in seen:
                raise ValueError('产物清单包含未支持或重复文件')
            seen.add(name)
            path = (folder / name).resolve()
            if path.parent != folder.resolve() or not path.is_file() or path.stat().st_size != artifact.get('bytes'):
                raise ValueError('产物文件缺失或大小不匹配')
            if sha256(path) != artifact.get('sha256'):
                raise ValueError('产物摘要不匹配')
            if name.endswith('.png') and path.read_bytes()[:8] != b'\x89PNG\r\n\x1a\n':
                raise ValueError('图表文件格式不正确')
        with (folder / 'coverage.cpm.bw').open('rb') as stream:
            if stream.read(4) != b'\x26\xfc\x8f\x88':
                raise ValueError('全基因组 bigWig 文件格式不正确')
        track = summary.get('coverage', {})
        if track.get('scope') != 'whole_genome' or track.get('reference_count') != 195 or track.get('normalization') != 'CPM':
            raise ValueError('轨道范围或标准化不符合约定')
        points = track.get('preview', {}).get('points', [])
        if len(points) != 200 or any(v is not None and (type(v) not in (int, float) or not math.isfinite(v)) for v in points):
            raise ValueError('轨道预览数据无效')
        return summary, rows

    def run(self, run_id, sample_id, source, cancelled, phase, *, resume_from=None):
        checkpoint = None
        if resume_from is not None:
            original_sample, checkpoint = self.resume_inputs(resume_from)
            if original_sample != sample_id or source != 'resume':
                raise ValueError('续跑来源与样本不匹配')
        runtime = self.probe()
        if not runtime.get('ready'):
            raise ValueError(runtime['reason'])
        sample, support, inputs = self.inputs(sample_id, source, checkpoint=checkpoint)
        folder = self.output_dir(run_id)
        folder.parent.mkdir(parents=True, exist_ok=True)
        free = shutil.disk_usage(folder.parent).free
        minimum = RESUME_START_FREE if checkpoint else 30 * 1024 ** 3
        if free < minimum:
            raise ValueError('本地续跑需要至少 60 GB 可用空间' if checkpoint else '本地可用磁盘不足 30 GiB，请先安排计算空间')
        before = {str(p.relative_to(self.root)): [p.stat().st_size, p.stat().st_mtime_ns] for p in inputs}
        folder.mkdir(exist_ok=False)
        argv = self.command(run_id, sample, support, folder, runtime['image_id'], source, checkpoint=checkpoint)
        execution = {'argv': argv, 'runtime': runtime, 'source': source, 'inputs': before,
                     'timeout_seconds': self.timeout, 'limits': {'cpus': 4, 'memory_gb': 6, 'network': 'none', 'input_mounts': 'readonly'},
                     'output_directory': str(folder.relative_to(self.root)), 'exit_code': None}
        if checkpoint:
            execution.update(resumed_from=resume_from, disk_reserve_bytes=DISK_RESERVE, min_free_bytes_observed=free)
        phase('inputs', '准备基础分析输入', execution)
        started, status, proc, last = time.monotonic(), 'failed', None, ''
        scratch = folder / '.scratch'
        try:
            if cancelled():
                status = 'cancelled'
            else:
                scratch.mkdir()
                with (folder / 'tool.log').open('wb') as stream:
                    proc = subprocess.Popen(argv, stdout=stream, stderr=subprocess.STDOUT, start_new_session=True)
                    while proc.poll() is None:
                        low_space = False
                        if checkpoint:
                            free = shutil.disk_usage(folder).free
                            execution['min_free_bytes_observed'] = min(execution['min_free_bytes_observed'], free)
                            low_space = free <= DISK_RESERVE
                        stop_requested = cancelled()
                        expired = time.monotonic() - started > self.timeout
                        if low_space or stop_requested or expired:
                            status = 'failed' if low_space else ('cancelled' if stop_requested else 'timed_out')
                            if low_space:
                                execution['error'] = '可用磁盘空间已接近 50 GB，已触发保护停止。'
                            self.cancel(run_id)
                            try:
                                proc.wait(timeout=5)
                            except subprocess.TimeoutExpired:
                                proc.kill()
                                proc.wait(timeout=5)
                            break
                        current = read_record(folder / 'progress.json')
                        key = current.get('stage')
                        if key != last and key in dict(STAGES):
                            phase(key, current.get('message', dict(STAGES)[key]), execution)
                            last = key
                        time.sleep(.5)
                    execution['exit_code'] = proc.returncode
                if status not in ('cancelled', 'timed_out') and 'error' not in execution:
                    if cancelled():
                        status = 'cancelled'
                    elif proc.returncode == 0:
                        phase('validate', '正在核对完整产物及摘要', execution)
                        summary, artifacts = self.validate(folder, source)
                        status = 'completed_pipeline'
                    else:
                        execution['error'] = read_record(folder / 'failure.json').get('error', '工具失败，请查看运行日志')
        except (OSError, ValueError, subprocess.SubprocessError) as exc:
            execution['error'] = str(exc)
        finally:
            if status != 'completed_pipeline':
                self.cancel(run_id)
            if proc is not None and proc.poll() is None:
                proc.kill()
                proc.wait(timeout=5)
            try:
                if scratch.is_dir():
                    shutil.rmtree(scratch)
            except OSError as exc:
                execution['scratch_cleanup_error'] = str(exc)
        execution['elapsed_seconds'] = round(time.monotonic() - started, 3)
        execution['inputs_unchanged'] = all(p.is_file() and [p.stat().st_size, p.stat().st_mtime_ns] == before[str(p.relative_to(self.root))] for p in inputs)
        if not execution['inputs_unchanged']:
            status = 'failed'
            execution['error'] = '输入元数据发生变化，结果未发布'
        commands = read_record(folder / 'commands.json')
        result = {'status': status, 'execution': execution, 'actual_tool_executions': len(commands) if isinstance(commands, list) else 0,
                  'actual_model_requests': 0, 'artifacts': []}
        if status == 'completed_pipeline':
            result.update(report_metrics=self.report_metrics(folder, source), result=summary, artifacts=artifacts + [{'name': 'tool.log', 'bytes': (folder / 'tool.log').stat().st_size}])
        else:
            result['artifacts'] = [{'name': name, 'bytes': (folder / name).stat().st_size} for name in ('tool.log', 'commands.json', 'bowtie2.log') if (folder / name).is_file()]
        return result
