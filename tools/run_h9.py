"""Run the public H9 workflow locally and keep each validated BAM for later review."""
import argparse
import fcntl
import json
import os
from pathlib import Path
import shutil
import signal
import subprocess
import sys
import time
from uuid import uuid4

from cuttag_web.h9_review import CONTEXT
from cuttag_web.model import ModelClient
from cuttag_web.pipeline import RESUME_IMAGE
from tools.h9_samples import CONTROLS, GEO_SAMPLES, SAMPLES, TARGETS


ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / 'data/public/GSE200768_H9_H3K4me1'
REFERENCE = ROOT / 'data/support/shared-reference/reference/hg38'
SIZE_TABLE = REFERENCE / 'hg38.bam-compatible.chrom.sizes'
INDEX = REFERENCE / 'bowtie2'
OUTPUT = ROOT / 'outputs/h9_development'
MIB = 1024 ** 2


def read(path):
    return json.loads(path.read_text())


def write(path, value):
    temp = path.with_suffix('.tmp')
    temp.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + '\n')
    temp.replace(path)


def local_inputs():
    checks = read(DATA / 'metadata/fastq_validation.json')
    manifest = read(DATA / 'metadata/raw_download_manifest.json')
    metadata = read(DATA / 'metadata/selected_sample_records.json')
    records = {row['accession']: row['fields'] for row in metadata['samples']}
    verified = {row['file']: row for row in checks['files']}
    listed = {row['local_relative_path']: row for row in manifest['files']}
    if (checks['status'] != 'passed' or checks['checked_files'] != 8 or manifest['source_is_public'] is not True
            or set(verified) != set(listed)):
        raise ValueError('H9 输入需要已校验的公开下载清单')
    result = {}
    for sample in SAMPLES:
        accession = GEO_SAMPLES[sample]
        role = 'target' if sample in TARGETS else 'igg'
        label = 'H3K4me1' if role == 'target' else 'IgG'
        replicate = (TARGETS if role == 'target' else CONTROLS).index(sample) + 1
        fields = records[accession]
        antibody = ('H3K4me1(Abcam, ab8895)' if role == 'target' else 'IgG(Santa Cruz Biotechnology, sc-2027)')
        if (fields['!Sample_title'] != [f'H9, {label} Rep{replicate}']
                or fields['!Sample_source_name_ch1'] != ['H9']
                or 'antibody: ' + antibody not in fields['!Sample_characteristics_ch1']
                or 'Assembly: hg38' not in fields['!Sample_data_processing']):
            raise ValueError('公开 H9 样本角色、重复、抗体或参考版本不匹配')
        paths = [DATA / 'raw' / sample / f'{sample}_{mate}.fastq.gz' for mate in (1, 2)]
        counts = []
        for path in paths:
            relative = str(path.relative_to(DATA))
            row, check = listed[relative], verified[relative]
            if (path.resolve() != path or row['geo_accession'] != accession or row['run'] != sample
                    or check['ena_md5_match'] is not True or check['gzip_crc_checked_to_eof'] is not True
                    or path.stat().st_size != row['bytes']):
                raise ValueError('公开 H9 FASTQ 与已校验清单不一致')
            counts.append(check['records'])
        if counts[0] != counts[1]:
            raise ValueError('公开 H9 双端读段数不一致')
        result[sample] = {'role': role, 'replicate': replicate, 'geo': accession, 'read_pairs': counts[0],
                          'antibody': antibody,
                          'files': {str(path.relative_to(ROOT)): [path.stat().st_size, path.stat().st_mtime_ns]
                                    for path in paths}}
    return result


def mount(source, target, writable=False):
    return ['--mount', f'type=bind,src={source},dst={target}' + ('' if writable else ',readonly')]


def stage(folder, mounts, args, receipt, *, native_input=None):
    if folder.exists():
        saved = read(folder / 'stage.json')
        if saved['inputs'] != receipt:
            raise ValueError(f'{folder.name} 输入已变化，需要独立输出目录')
        for name, size in saved['outputs'].items():
            if (folder / name).stat().st_size != size:
                raise ValueError(f'{folder.name} 已登记的输出损坏')
        print(f'复用 {folder.name}', flush=True)
        return
    if shutil.disk_usage(folder.parent).free < 30 * 1024 * MIB:
        raise ValueError('H9 本地分析需要至少 30 GiB 可用空间')
    attempt = folder.with_name(folder.name + '.attempt-' + uuid4().hex[:8])
    attempt.mkdir()
    scratch = attempt / '.scratch'
    scratch.mkdir()
    container = 'cuttag-h9-' + uuid4().hex[:12]
    command = ['docker', 'run', '--rm', '--name', container, '--pull', 'never', '--network', 'none', '--read-only', '--cpus', '4',
               '--memory', '6g', '--cap-drop', 'ALL', '--security-opt', 'no-new-privileges', '--user', '1000:1000',
               '--env', 'TMPDIR=/output/.scratch', '--env', 'MPLCONFIGDIR=/output/.scratch/matplotlib',
               '--env', 'PYTHONPATH=/code', '--entrypoint', 'python']
    for file in ('tools/h9_worker.py', 'tools/h9_samples.py', 'tools/pipeline/worker.py'):
        command += mount(ROOT / file, '/code/' + file)
    command += mount(SIZE_TABLE, '/support/reference/hg38/' + SIZE_TABLE.name)
    command += mount(attempt, '/output', writable=True) + mounts
    command += [receipt['image_id'], '-B', '-m', 'tools.h9_worker', *args]
    env = None
    if native_input is not None:
        from tools.h9_native import isolated_env
        env = isolated_env(scratch)
        command = [sys.executable, '-B', '-m', 'tools.h9_worker', '--native-root', str(ROOT),
                   '--output', str(attempt), '--input', str(native_input), *args]
    print(f'开始 {folder.name}，日志 {attempt / "tool.log"}', flush=True)
    start = time.monotonic()
    try:
        with (attempt / 'tool.log').open('wb') as log:
            if native_input is None:
                code = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, timeout=4 * 3600).returncode
            else:
                with subprocess.Popen(command, env=env, stderr=subprocess.STDOUT, stdout=log,
                                      start_new_session=True) as process:
                    try:
                        code = process.wait(timeout=24 * 3600)
                    except BaseException:
                        os.killpg(process.pid, signal.SIGTERM)
                        try:
                            process.wait(timeout=8)
                        except subprocess.TimeoutExpired:
                            os.killpg(process.pid, signal.SIGKILL)
                        raise
    except (KeyboardInterrupt, subprocess.TimeoutExpired):
        if native_input is None:
            subprocess.run(['docker', 'rm', '-f', container], capture_output=True, timeout=10, check=False)
        raise
    if code:
        raise RuntimeError(f'{folder.name} 退出码 {code}，查看 {attempt / "tool.log"}')
    if args[0] == 'prepare':
        wanted = ('prepared.json', 'analysis.bam', 'analysis.bam.bai')
        info = read(attempt / 'prepared.json')
        expected_pairs = receipt['pair_limit'] or receipt['sample']['read_pairs']
        if info['sample'] != args[1] or info['raw_read_pairs'] != expected_pairs:
            raise ValueError('H9 处理结果与选定样本或读段数不符')
        shutil.rmtree(attempt / 'work')
    elif args[0] == 'track':
        wanted = ('track.json', 'coverage.cpm.bw')
    else:
        wanted = ('metrics.json', *(sample + '.peaks.bed' for sample in TARGETS))
        metrics = read(attempt / 'metrics.json')
        if metrics['controls'] != (list(CONTROLS) if args[0] == 'compare-igg' else []):
            raise ValueError('H9 峰结果的对照配置不符')
    outputs = {name: (attempt / name).stat().st_size for name in wanted}
    if any(size == 0 for name, size in outputs.items() if not name.endswith('.bed')):
        raise ValueError(f'{folder.name} 的必要文件为空')
    shutil.rmtree(scratch)
    elapsed = round(time.monotonic() - start, 2)
    write(attempt / 'stage.json', {'inputs': receipt, 'outputs': outputs, 'seconds': elapsed})
    attempt.rename(folder)
    print(f'完成 {folder.name}，{elapsed} 秒', flush=True)


def model_metrics(path):
    metrics = read(path)
    return {'peak_mode': metrics['peak_mode'], 'qvalue': metrics['qvalue'], 'control_count': len(metrics['controls']),
            'peak_counts': [metrics['peak_counts'][sample] for sample in TARGETS],
            'target_frip': [metrics['target_frip'][sample]['fraction'] for sample in TARGETS],
            'igg_overlap': [metrics['igg_overlap'][control][target]['fraction']
                            for control in CONTROLS for target in TARGETS],
            'replicate_peak_jaccard': metrics['replicate_peak_jaccard']}


def review_model(output, scope, compare):
    model = ModelClient(ROOT)
    attempts = output / 'reviews'
    attempts.mkdir(exist_ok=True)
    path = attempts / (uuid4().hex + '.json')
    record = {'status': 'running', 'input_scope': scope, 'model_requests': 0, 'tool_comparisons': 0, 'rounds': []}
    write(path, record)
    print(f'模型决策记录 {path}', flush=True)
    baseline = model_metrics(output / 'baseline/metrics.json')
    pair_counts = [read(output / 'samples' / sample / 'prepared.json')['raw_read_pairs'] for sample in SAMPLES]

    def request(round_index, with_igg=None):
        payload = {'data_source': 'public_h9_summary', 'input_scope': scope, 'read_pairs_per_sample': pair_counts,
                   'round_index': round_index, 'experiment_context': CONTEXT, 'baseline': baseline, 'with_igg': with_igg}
        entry = {'summary': payload}
        record['rounds'].append(entry)

        def sent():
            record['model_requests'] += 1
            write(path, record)

        entry['answer'] = model.review_public_h9(payload, on_request=sent)
        write(path, record)
        return entry['answer']['decision']['action']

    try:
        action = request(0)
        if action == 'compare_igg':
            reused = (output / 'with-igg').exists()
            compare()
            record.update(tool_comparisons=int(not reused), reused_comparison=reused)
            action = request(1, model_metrics(output / 'with-igg/metrics.json'))
        record.update(status='completed', final_action=action)
    except (ValueError, OSError, RuntimeError) as exc:
        record.update(status='failed', error=str(exc))
        raise
    finally:
        write(path, record)


def run(pair_limit, prepare_only=False, *, model=False, compare_tools=False, native_recipe=None):
    sources = local_inputs()
    output = OUTPUT / ('full' if pair_limit is None else f'prefix-{pair_limit}')
    samples = output / 'samples'
    samples.mkdir(parents=True, exist_ok=True)
    with (output / '.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        docker_image = None if native_recipe else subprocess.check_output(
            ['docker', 'image', 'inspect', '--format', '{{.Id}}', RESUME_IMAGE], text=True).strip()
        recipe = {'image_id': docker_image, 'pair_limit': pair_limit}
        if native_recipe:
            recipe['native'] = native_recipe

        def process(folder, mounts, args, receipt, input_path):
            stage(folder, mounts, args, receipt, native_input=input_path if native_recipe else None)

        scope = 'full' if pair_limit is None else 'prefix_smoke'
        write(output / 'experiment.json', {'context': CONTEXT, 'samples': sources, 'input_scope': scope})
        for sample in SAMPLES:
            args = ['prepare', sample]
            if pair_limit is not None:
                args += ['--pair-limit', str(pair_limit)]
            inputs_path = DATA / 'raw' / sample
            process(samples / sample, mount(inputs_path, '/reads') + mount(INDEX, '/support/reference/hg38/bowtie2'),
                    args, {**recipe, 'sample': sources[sample]}, inputs_path)
        if not prepare_only:
            bams = {sample: read(samples / sample / 'stage.json') for sample in SAMPLES}
            tracks = output / 'tracks'
            tracks.mkdir(exist_ok=True)
            for sample in SAMPLES:
                process(tracks / sample, mount(samples / sample, '/sample'), ['track'],
                        {**recipe, 'bam': bams[sample]}, samples / sample)
            inputs = {**recipe, 'samples': bams}
            mounts = mount(samples, '/samples')
            process(output / 'baseline', mounts, ['baseline'], inputs, samples)

            def compare():
                process(output / 'with-igg', mounts, ['compare-igg'], inputs, samples)

            if model:
                review_model(output, scope, compare)
            elif compare_tools:
                compare()
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--pairs-per-sample', type=int, help='限量工程联调，按原序读取每份 FASTQ 的前 N 对')
    parser.add_argument('--prepare-only', action='store_true', help='只生成并保留 BAM 和索引')
    review = parser.add_mutually_exclusive_group()
    review.add_argument('--model', action='store_true', help='使用已配置模型复核，仅发送公开 H9 数值摘要')
    review.add_argument('--compare-tools', action='store_true', help='直接实测带 IgG 的工具动作，跳过模型')
    args = parser.parse_args()
    if args.pairs_per_sample is not None and not 1 <= args.pairs_per_sample <= 1_000_000:
        parser.error('--pairs-per-sample 应在 1 到 1,000,000 范围内')
    if args.prepare_only and (args.model or args.compare_tools):
        parser.error('--prepare-only 只执行样本准备')
    print(run(args.pairs_per_sample, args.prepare_only, model=args.model, compare_tools=args.compare_tools), flush=True)


if __name__ == '__main__':
    main()
