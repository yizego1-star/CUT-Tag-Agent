"""Build a local paired CUT&Tag case from verified target and IgG BAMs."""
import argparse
from contextlib import contextmanager
import fcntl
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess
import time
from uuid import uuid4

from cuttag_web.pipeline import RESUME_IMAGE


ROOT = Path(__file__).resolve().parent.parent
INTAKE = ROOT / 'data/support/paired-case'
REFERENCE = ROOT / 'data/support/shared-reference'
OUTPUT = ROOT / 'outputs/hospital_pair/paired-target'
SAMPLES = {'target': 'paired-target', 'igg': 'paired-igg'}
PILOT_REGION = 'chr1:1000001-2000000'  # samtools uses 1-based inclusive coordinates
TRACK_REGION = 'chr1:1000000:2000000'
REF_FILES = ('reference/hg38/hg38.bam-compatible.chrom.sizes',
             'annotation/gencode.v44.genes.bam-compatible.bed')
TSS_FILES = ('tss.matrix.gz', 'tss.profile.tsv', 'tss.profile.png', 'tss.heatmap.png', 'tss.regions.bed')
WORKER_FILES = ('tools/prepare_hospital_pair.py', 'tools/pipeline/worker.py', 'cuttag_web/pipeline.py',
                'cuttag_web/coverage.py')


def load(path):
    return json.loads(path.read_text())


def write(path, value):
    temp = path.with_name(path.name + '.tmp')
    temp.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + '\n')
    temp.replace(path)


def local_path(path):
    if not path.is_file() or path.resolve() != path or ROOT not in path.parents:
        raise ValueError(f'项目输入缺失或路径越界 {path}')
    return path


def file_identity(path):
    stat = local_path(path).stat()
    return {'bytes': stat.st_size, 'mtime_ns': stat.st_mtime_ns}


def inventory():
    provenance = INTAKE / 'provenance'
    pair = load(local_path(provenance / 'pair_matching_and_qc_evidence.json'))
    checks = load(local_path(provenance / 'paired_material_validation.json'))
    records = load(local_path(provenance / 'paired_download_manifest.json'))['files']
    context = load(local_path(provenance / 'experiment_context_evidence.json'))
    if (pair['target'], pair['control']) != tuple(SAMPLES.values()) or not pair['matched_control_confirmed']:
        raise ValueError('本地登记的对照关系与本案例不一致')
    if checks['status'] != 'passed':
        raise ValueError('配对数据包没有通过已有基础校验')
    by_sample = {row['sample']: row for row in pair['samples']}
    samples = {}
    for role, sample in SAMPLES.items():
        prefix = f'inputs/{sample}/'
        items = {}
        for row in records:
            name = row['local_relative_path']
            if not name.startswith(prefix):
                continue
            path = INTAKE / name
            identity = file_identity(path)
            if identity['bytes'] != row['actual_bytes'] or row['status'] != 'downloaded':
                raise ValueError(f'配对材料和已校验清单不一致 {path}')
            items[str(path.relative_to(ROOT))] = {**identity, 'recorded_sha256': row['sha256']}
        historical = by_sample[sample]
        samples[role] = {'sample_id': sample, 'input_stage': 'historical_filtered_bam', 'materials': items,
                         'fastp_before': historical['fastp_before'], 'fastp_after': historical['fastp_after'],
                         'historical_spikein': historical['summary']}
    return {
        'status': 'source_materials_verified',
        'scope': 'local_paired_case',
        'reference': {name: file_identity(REFERENCE / name) for name in REF_FILES},
        'samples': samples,
        'control_assignment': 'local_pair_record',
        'peak_mode': 'exploratory_narrow',
        'unconfirmed': ['exact_mark_and_antibody', 'biological_replicates', 'spikein_dose_and_addition_stage',
                        'biological_control_signoff', 'manual_igv_review'],
        'review_windows': context['review_windows'],
        'not_locally_available': ['raw_fastq_for_this_pair'],
        'spikein_quantitative_verdict': 'pending_experiment_record',
        'external_model_allowed': False,
    }


def mount(source, destination, writable=False):
    return ['--mount', f'type=bind,src={source},dst={destination}' + ('' if writable else ',readonly')]


def stage(folder, args, inputs, expected, receipt, *, cancelled=lambda: False, container=None, extra_files=(), timeout=4*3600):
    if folder.exists():
        saved = load(folder / 'stage.json')
        if saved['inputs'] != receipt or set(saved['outputs']) != set(expected):
            raise ValueError(f'{folder.name} 运行参数或输入已变，请保留已有结果并使用独立目录')
        for name, value in saved['outputs'].items():
            if (folder / name).stat().st_size != value:
                raise ValueError(f'阶段产物尺寸已变化 {folder / name}')
        print(f'复用 {folder.name}', flush=True)
        return
    attempt = folder.with_name(folder.name + '.attempt-' + uuid4().hex[:8])
    attempt.mkdir()
    scratch = attempt / '.scratch'
    scratch.mkdir()
    container = container or 'cuttag-pair-' + uuid4().hex[:12]
    command = ['docker', 'run', '--rm', '--name', container, '--pull', 'never', '--network', 'none',
               '--read-only', '--cpus', '4', '--memory', '6g', '--cap-drop', 'ALL', '--user', '1000:1000',
               '--env', 'TMPDIR=/output/.scratch', '--env', 'MPLCONFIGDIR=/output/.scratch/matplotlib',
               '--security-opt', 'no-new-privileges', '--env', 'PYTHONPATH=/code', '--entrypoint', 'python']
    for name in (*WORKER_FILES, *extra_files):
        command += mount(ROOT / name, '/code/' + name)
    for name in REF_FILES:
        command += mount(REFERENCE / name, '/support/' + name)
    for path, destination in inputs:
        command += mount(path, destination)
    command += mount(attempt, '/output', writable=True)
    command += [receipt['image'], '-B', '-m', 'tools.prepare_hospital_pair', '--worker', *args]
    print(f'开始 {folder.name}，日志 {attempt / "tool.log"}', flush=True)
    start = time.monotonic()
    proc = None
    try:
        with (attempt / 'tool.log').open('wb') as log:
            if cancelled():
                raise InterruptedError('已停止任务')
            proc = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT)
            while proc.poll() is None:
                if cancelled():
                    raise InterruptedError('已停止任务')
                if time.monotonic() - start > timeout:
                    raise subprocess.TimeoutExpired(command, timeout)
                try:
                    proc.wait(timeout=1)
                except subprocess.TimeoutExpired:
                    continue
    except (KeyboardInterrupt, InterruptedError, subprocess.TimeoutExpired):
        if proc is not None:
            stop_container(container)
            if proc.poll() is None:
                proc.kill()
            proc.wait(timeout=10)
        raise
    if proc.returncode:
        raise ValueError(f'{folder.name} 退出码 {proc.returncode}，请查 {attempt / "tool.log"}')
    outputs = {name: (attempt / name).stat().st_size for name in expected}
    if any(size == 0 and not name.endswith('.peaks.bed') for name, size in outputs.items()):
        raise ValueError(f'{folder.name} 必要产物为空')
    shutil.rmtree(scratch)
    elapsed = round(time.monotonic() - start, 2)
    write(attempt / 'stage.json', {'inputs': receipt, 'outputs': outputs, 'seconds': elapsed})
    attempt.rename(folder)
    print(f'完成 {folder.name}，{elapsed} 秒', flush=True)


def stop_container(name):
    try:
        subprocess.run(['docker', 'rm', '-f', name], capture_output=True, timeout=10, check=False)
    except (OSError, subprocess.TimeoutExpired):
        pass


def source_mounts():
    for role, sample in SAMPLES.items():
        bam = INTAKE / f'inputs/{sample}/bam/{sample}.hg38.filtered.bam'
        yield local_path(bam), f'/alignments/{role}.bam'
        yield local_path(bam.with_suffix('.bam.bai')), f'/alignments/{role}.bam.bai'


def rerun_dir(root, run_id):
    if not re.fullmatch(r'run-[0-9a-f]{12}', run_id):
        raise ValueError('调参任务编号无效')
    folder = root / 'outputs/paired_reruns' / run_id
    if folder.resolve() != folder:
        raise ValueError('调参产物目录越界')
    return folder


def rerun_peaks(run_id, qvalue, cancelled, *, inspect_igv=False):
    from cuttag_web.hospital_pair_review import QVALUES, validate_metrics

    if qvalue not in QVALUES:
        raise ValueError('不支持的调参阈值')
    baseline = load(local_path(OUTPUT / 'full/peaks/stage.json'))['inputs']
    inputs = list(source_mounts())
    bam_inputs = {dest: file_identity(path) for path, dest in inputs}
    reference = {name: file_identity(REFERENCE / name) for name in REF_FILES}
    if bam_inputs != baseline['bam_inputs'] or reference != baseline['reference']:
        raise ValueError('调参与基线的 BAM、索引或参考文件不一致')
    image = subprocess.check_output(['docker', 'image', 'inspect', '--format', '{{.Id}}', RESUME_IMAGE],
                                    text=True, timeout=10).strip()
    if image != baseline['image']:
        raise ValueError('调参与基线的 MACS3 容器镜像不一致')
    folder = rerun_dir(ROOT, run_id)
    folder.parent.mkdir(parents=True, exist_ok=True)
    receipt = {'image': image, 'bam_inputs': bam_inputs, 'reference': reference, 'qvalue': qvalue,
               'implementation': {name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest()
                                  for name in (*WORKER_FILES, 'cuttag_web/hospital_pair_review.py')}}
    names = ('tuned.peaks.bed', 'tuned.json', 'commands.json', 'tool.log')
    args = ['tune', str(qvalue)]
    extra_files = ['cuttag_web/hospital_pair_review.py']
    if inspect_igv:
        args.append('inspect_igv')
        script = 'tools/inspect_paired_igv.py'
        extra_files.append(script)
        receipt['implementation'][script] = hashlib.sha256((ROOT / script).read_bytes()).hexdigest()
        full = OUTPUT / 'full'
        evidence_names = ('target/coverage.cpm.bw', 'igg/coverage.cpm.bw', 'peaks/baseline.peaks.bed', 'peaks/with-igg.peaks.bed')
        receipt['igv_inputs'] = {name: file_identity(full / name) for name in evidence_names}
        inputs += [(full / name, '/baseline/' + name) for name in evidence_names]
    stage(folder, args, inputs, names, receipt, cancelled=cancelled,
          container='cuttag-pair-' + run_id[4:], extra_files=extra_files)
    measured = {dest: file_identity(path) for path, dest in inputs}
    expected_inputs = {**bam_inputs, **{'/baseline/' + name: value for name, value in receipt.get('igv_inputs', {}).items()}}
    if (measured != expected_inputs
            or {name: file_identity(REFERENCE / name) for name in REF_FILES} != reference):
        raise ValueError('调参过程中 BAM、索引或参考发生变化')
    result = load(local_path(folder / 'tuned.json'))
    validate_metrics({key: result[key] for key in ('peak_count', 'target', 'igg')})
    with (folder / 'tuned.peaks.bed').open() as bed:
        count = sum(1 for _ in bed)
    commands = load(folder / 'commands.json')
    peak = [row for row in commands if row['argv'][1:2] == ['callpeak']]
    result_fields = {'qvalue', 'peak_count', 'target', 'igg'} | ({'igv_changes'} if inspect_igv else set())
    if (set(result) != result_fields or result['qvalue'] != qvalue
            or count != result['peak_count'] or len(peak) != 1 or any(row['exit_code'] != 0 for row in commands)
            or peak[0]['argv'][peak[0]['argv'].index('-q') + 1] != str(qvalue)
            or peak[0]['argv'][peak[0]['argv'].index('-c') + 1] != '/alignments/igg.bam'):
        raise ValueError('新候选峰和成功执行的工具记录不一致')
    return result, [{'name': name, 'bytes': local_path(folder / name).stat().st_size} for name in names]


def prepare(scope):
    metadata = inventory()
    OUTPUT.mkdir(mode=0o700, parents=True, exist_ok=True)
    with (OUTPUT / '.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        write(OUTPUT / 'source-materials.json', metadata)
        if scope == 'inventory':
            return OUTPUT
        if scope == 'full' and shutil.disk_usage(OUTPUT).free < 30 * 1024 ** 3:
            raise ValueError('全基因组分析需要至少 30 GiB 本地可用空间')
        dest = OUTPUT / scope
        dest.mkdir(mode=0o700, exist_ok=True)
        image = subprocess.check_output(['docker', 'image', 'inspect', '--format', '{{.Id}}', RESUME_IMAGE],
                                        text=True, timeout=10).strip()
        inputs = list(source_mounts())
        receipt = {'image': image, 'scope': scope, 'reference': metadata['reference'],
                   'bam_inputs': {target: file_identity(path) for path, target in inputs},
                   'implementation': {name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest()
                                      for name in WORKER_FILES}}
        if scope == 'pilot':
            alignment = dest / 'alignments'
            stage(alignment, ['subset'], inputs, tuple(f'{role}.bam{ext}' for role in SAMPLES for ext in ('', '.bai')),
                  receipt)
            inputs = [(alignment, '/alignments')]
            receipt['pilot_alignments'] = load(alignment / 'stage.json')['outputs']
        for role in SAMPLES:
            expected = ('coverage.cpm.bw', 'track.json', 'flagstat.json', 'idxstats.tsv')
            if scope == 'full':
                expected += TSS_FILES
            stage(dest / role, ['track', role, scope], inputs, expected, receipt)
        stage(dest / 'peaks', ['peaks', scope], inputs,
              ('comparison.json', 'baseline.peaks.bed', 'with-igg.peaks.bed'), receipt)
        tracks = {role: load(dest / role / 'track.json') for role in SAMPLES}
        comparison = load(dest / 'peaks/comparison.json')
        if any(tracks[role]['fragments'] != comparison['fragment_metrics'][role]['fragments'] for role in SAMPLES):
            raise ValueError('轨道与候选峰的配对片段统计不一致')
        write(dest / 'agent-evidence.json', {
            'status': 'pilot_ready' if scope == 'pilot' else 'tool_artifacts_ready',
            'scope': PILOT_REGION if scope == 'pilot' else 'whole_genome',
            'source_materials': '../source-materials.json',
            'track_summaries': tracks,
            'peak_comparison': comparison,
            'can_compare_peak_workflow': True,
            'can_approve_spikein_standardization': False,
            'scientific_verdict': 'pending',
            'unconfirmed': metadata['unconfirmed'],
            'external_model_allowed': False,
        })
        return dest


def worker_main(args):
    import pyBigWig
    import pysam
    from tools.pipeline import worker

    if args == ['screen_igv']:
        from tools.inspect_paired_igv import inspect
        worker.write('interval-summary.json', inspect(Path('/baseline')))
        return

    if len(args) == 3 and args[0] == 'plot_pair':
        from tools.pair_comparison_figure import FIGURE_NAME, LOCUS_NAME, display_comparison
        loci = display_comparison(Path('/baseline'), Path('/tuned.peaks.bed'), tuple(map(float, args[1:])),
                                  worker.OUT / FIGURE_NAME)
        worker.write(LOCUS_NAME, {'loci': loci})
        return

    chroms = worker.sizes()
    bams = {role: f'/alignments/{role}.bam' for role in SAMPLES}

    @contextmanager
    def checked_bam(path):
        worker.run(['samtools', 'quickcheck', '-v', path])
        with pysam.AlignmentFile(path, 'rb') as bam:
            if dict(zip(bam.references, bam.lengths)) != chroms or not bam.check_index():
                raise ValueError('BAM 参考或索引与本案例不匹配')
            yield bam

    if args == ['subset']:
        for role, path in bams.items():
            with checked_bam(path):
                pass
            output = str(worker.OUT / f'{role}.bam')
            worker.run(['samtools', 'view', '--fetch-pairs', '-b', '-o', output, path, PILOT_REGION])
            worker.run(['samtools', 'index', output])
        return

    if args[0] == 'track':
        _, role, scope = args
        bam = bams[role]
        with checked_bam(bam) as expected:
            mapped = expected.mapped
        worker.run(['samtools', 'flagstat', '-O', 'json', bam], 'flagstat.json')
        worker.run(['samtools', 'idxstats', bam], 'idxstats.tsv')
        counts = load(worker.OUT / 'flagstat.json')
        good, bad = counts['QC-passed reads'], counts['QC-failed reads']
        if (not mapped or mapped % 2 or bad['total'] or good['read1'] != mapped // 2
                or good['read2'] != mapped // 2 or any(good[key] != mapped for key in
                    ('total', 'primary', 'mapped', 'properly paired', 'with itself and mate mapped'))
                or good['with mate mapped to a different chr']):
            raise ValueError('BAM 实测配对计数不符合按片段计算 CPM 的约定')
        fragments = mapped // 2
        if scope == 'full':
            coverage = worker.coverage(bam, chroms, fragments=fragments)
            tss = worker.tss()
            coverage.pop('preview')
        else:
            track = str(worker.OUT / 'coverage.cpm.bw')
            worker.run(['bamCoverage', '-b', bam, '-o', track, '--region', TRACK_REGION, '--binSize', '25',
                        '--scaleFactor', str(1_000_000 / fragments), '--extendReads', '--samFlagInclude', '64',
                        '--samFlagExclude', '2828', '--numberOfProcessors', '3'])
            with pyBigWig.open(track) as bw:
                if bw.chroms().get('chr1') != chroms['chr1'] or not bw.header()['nBasesCovered']:
                    raise ValueError('限量轨道校验失败')
            coverage, tss = {'scope': PILOT_REGION, 'normalization': 'CPM_within_pilot', 'bin_size': 25}, None
        worker.write('track.json', {'fragments': fragments, 'coverage': coverage, 'tss': tss})
        return

    if args[0] == 'tune':
        from cuttag_web.hospital_pair_review import QVALUES
        if args[2:] not in ([], ['inspect_igv']) or len(args) < 2 or args[1] not in tuple(map(str, QVALUES)):
            raise ValueError('不支持的调参阈值')
        qvalue = float(args[1])
        for bam in bams.values():
            with checked_bam(bam):
                pass
        (worker.OUT / 'peaks').mkdir()
        result = worker.call_peaks(bams['target'], chroms, 'tuned', controls=[bams['igg']], qvalue=qvalue)
        worker.stage('qc', '正在扫描目标和 IgG 的配对片段，统计新峰内片段比例')
        metrics = {role: worker.fragment_counts(bam, {'tuned': result}) for role, bam in bams.items()}
        summary = {'qvalue': qvalue, 'peak_count': result['count'],
                   **{role: {'fragments': row['fragments'], 'frip': row['frip']['tuned']['fraction']}
                      for role, row in metrics.items()}}
        if args[2:] == ['inspect_igv']:
            from tools.inspect_paired_igv import changes_after_rerun
            worker.stage('igv', '正在核对重跑找回了哪些原目标和 IgG 信号区间')
            summary['igv_changes'] = changes_after_rerun(Path('/baseline'), worker.OUT / 'tuned.peaks.bed')
        worker.write('tuned.json', summary)
        return

    if args[0] != 'peaks':
        raise ValueError('未知的生信计算阶段')
    (worker.OUT / 'peaks').mkdir()
    results = {
        'baseline': worker.call_peaks(bams['target'], chroms, 'baseline'),
        'with-igg': worker.call_peaks(bams['target'], chroms, 'with-igg', controls=[bams['igg']]),
    }
    worker.write('comparison.json', {
        'scope': 'whole_genome' if args[1] == 'full' else PILOT_REGION,
        'settings': {'mode': 'exploratory_narrow', 'qvalue': 0.01, 'execution': 'MACS3 BAMPE',
                     'keep_dup': 'all', 'blacklist_filter': 'not_applied', 'igg_used_only_in': 'with-igg'},
        'peak_counts': {mode: result['count'] for mode, result in results.items()},
        'fragment_metrics': {role: worker.fragment_counts(path, results) for role, path in bams.items()},
    })


def main():
    import sys
    if sys.argv[1:2] == ['--worker']:
        worker_main(sys.argv[2:])
        return
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('scope', choices=('inventory', 'pilot', 'full'))
    args = parser.parse_args()
    print(prepare(args.scope), flush=True)


if __name__ == '__main__':
    main()
