"""Produce a separate, historically matched, spike-in scaled hg38 track locally."""
import argparse
import csv
from datetime import datetime, timezone
import json
import math
from pathlib import Path, PurePosixPath
import re
import shutil
import subprocess
import time
import uuid

from cuttag_web.coverage import read_record
from cuttag_web.pipeline import DISK_RESERVE, RESUME_IMAGE, RESUME_START_FREE, sha256


RUN_TRACK = r"""
import json, math, sys, subprocess
import pyBigWig, pysam

bam = '/human/analysis.bam'
compare = sys.argv[3] == 'compare'
bin_size = 25 if compare else 10
with pysam.AlignmentFile(bam, 'rb') as f:
    f.check_index()
    chroms = dict(zip(f.references, f.lengths))
    if sum(s.mapped for s in f.get_index_statistics()) != 2 * int(sys.argv[2]):
        raise ValueError('Historical BAM count differs from the recorded fragment count')
with open('/chrom.sizes') as s:
    assert chroms == {str(c): int(n) for c, n in map(str.split, s)}
subprocess.run(['samtools', 'quickcheck', '-v', bam], check=True)
base = ['bamCoverage', '-b', bam, '--outFileFormat', 'bigwig', '--binSize', str(bin_size),
        '--extendReads', '--samFlagInclude', '64', '--numberOfProcessors', '3']
if compare:
    base.append('--exactScaling')

def track(kind, normalization):
    path = '/output/coverage.' + kind + '.bw'
    cmd = base + ['-o', path, *normalization]
    print('Generating ' + kind, flush=True)
    subprocess.run(cmd, check=True)
    with pyBigWig.open(path) as bw:
        header = bw.header()
        points = bw.stats('chr1', 1000000, 2000000, nBins=200, type='mean', exact=True)
        if (not bw.isBigWig() or bw.chroms() != chroms or not header['nBasesCovered']
                or not any(v is not None and math.isfinite(v) and v > 0 for v in points)):
            raise ValueError(kind + ' BigWig verification failed')
    return {'reference_count': len(chroms), 'bin_size': bin_size, 'validation': 'passed',
            'preview': points, 'command': cmd, 'bases_covered': header['nBasesCovered']}

spikein = track('spikein', ['--normalizeUsing', 'None', '--scaleFactor', sys.argv[1]])
result = spikein
if compare:
    cpm = track('cpm', ['--normalizeUsing', 'CPM'])
    if cpm['bases_covered'] != spikein['bases_covered']:
        raise ValueError('The two BigWigs cover different base counts')
    ratios = [s / c for c, s in zip(cpm['preview'], spikein['preview']) if c and s]
    if not ratios or not all(math.isclose(r, ratios[0], rel_tol=1e-5) for r in ratios):
        raise ValueError('The local signal differs by more than a scaling factor')
    result = {'cpm': cpm, 'spikein': spikein, 'validation': 'passed',
              'local_spikein_to_cpm_ratio': sum(ratios) / len(ratios), 'compared_preview_bins': len(ratios)}
with open('/output/track.json', 'w') as file:
    json.dump(result, file, allow_nan=False)
"""


def matched_inputs(root, sample):
    if not re.fullmatch(r'[A-Za-z0-9_-][A-Za-z0-9_.-]{0,119}', sample):
        raise ValueError('样本名无效')
    private = root / 'data/private' / sample
    support = root / 'data/support' / sample
    original = read_record(private / 'provenance/source_manifest.json')
    supplement = read_record(support / 'provenance/spikein_sample_manifest.json')
    matches = read_record(support / 'provenance/spikein_sample_match.json')
    if any(item.get('sample_id') != sample for item in (original, supplement, matches)) or matches.get('status') != 'matched_historical_sample_files':
        raise ValueError('缺少当前样本的 spike-in 匹配记录')
    bam = private / 'processed' / (sample + '.hg38.dedup.bam')
    bai = bam.with_suffix('.bam.bai')
    bed = support / 'sample_spikein/bed' / (sample + '.spikein.fragments.bed')
    chroms = support / 'reference/hg38/hg38.bam-compatible.chrom.sizes'
    source = PurePosixPath(supplement['source_result_directory'])
    for path, manifest, base in ((bam, original, private), (bai, original, private), (bed, supplement, support)):
        item = next((x for x in manifest['files'] if x['local_relative_path'] == str(path.relative_to(base))), {})
        if (path.resolve() != path or not path.is_file() or item.get('expected_bytes') != path.stat().st_size
                or PurePosixPath(item.get('remote_path', '/')).parent.parent != source):
            raise ValueError('人源 BAM 与 spike-in 不属于同一次历史分析，或文件已变化')
        if path == bed and sha256(path) != item.get('sha256'):
            raise ValueError('spike-in 片段表校验失败')
    with (private / 'qc/summary.tsv').open(newline='') as f:
        row = next((r for r in csv.DictReader(f, delimiter='\t') if r['sample'] == sample), None)
    with bed.open() as f:
        count = sum(1 for line in f if line.strip())
    if not row or not count or int(row['spikein_fragments']) != count:
        raise ValueError('spike-in 片段计数缺失或与历史记录不符')
    factor = 10000 / count
    if not math.isclose(factor, float(row['scale_factor']), rel_tol=1e-9):
        raise ValueError('spike-in 计算系数与历史记录不符')
    return bam, bai, bed, chroms, count, factor, int(row['human_fragments'])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('sample')
    parser.add_argument('--no-matched-igg', action='store_true', help='仅在已核对没有匹配 IgG 时启用')
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    base = root / 'outputs/spikein_runs'
    base.mkdir(exist_ok=True)
    folder = base / ('run-' + uuid.uuid4().hex[:12])
    record = run_tracks(args.sample, folder, matched_inputs(root, args.sample), no_matched_igg=args.no_matched_igg)
    return int(record['status'] != 'completed_spikein')


def cancel_tracks(run_id):
    subprocess.run(['docker', 'rm', '-f', 'cuttag-spikein-' + run_id],
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=10)


def run_tracks(sample, folder, inputs, *, compare=False, no_matched_igg=False, cancelled=lambda: False):
    bam, bai, bed, chroms, count, factor, human = inputs
    if shutil.disk_usage(folder.parent).free < RESUME_START_FREE:
        raise ValueError('生成轨道前需要至少 60 GB 空闲空间')
    run_id = folder.parent.name if compare else folder.name
    folder.mkdir()
    scratch = folder / '.scratch'
    scratch.mkdir()
    name = 'cuttag-spikein-' + run_id
    before = [(p.stat().st_size, p.stat().st_mtime_ns) for p in (bam, bai, bed, chroms)]
    stamp = lambda: datetime.now(timezone.utc).isoformat()
    record = {'id': run_id, 'sample_id': sample, 'mode': 'normalization_comparison' if compare else 'spikein',
              'source': 'matched_historical_bam',
              'name': 'Spike-in 轨道 · 历史配套 BAM', 'created_at': stamp(), 'status': 'failed', 'artifacts': [],
              'actual_model_requests': 0, 'human_fragments': human,
              'spikein_fragments': count, 'scale_factor': factor, 'scale_numerator': 10000}
    if not compare and no_matched_igg:
        record['matched_igg'] = False
    bind = lambda path, target, access='': ['--mount', 'type=bind,src=' + str(path) + ',dst=' + target + access]
    argv = ['docker', 'run', '--rm', '--pull', 'never', '--name', name, '--network', 'none', '--read-only',
            '--cap-drop', 'ALL', '--security-opt', 'no-new-privileges', '--cpus', '4', '--memory', '6g',
            '--pids-limit', '256', '--user', '1000:1000', '--env', 'TMPDIR=/output/.scratch',
            '--env', 'MPLCONFIGDIR=/output/.scratch/matplotlib',
            *bind(bam, '/human/analysis.bam', ',readonly'), *bind(bai, '/human/analysis.bam.bai', ',readonly'),
            *bind(chroms, '/chrom.sizes', ',readonly'), *bind(folder, '/output'), '--entrypoint', 'python',
            RESUME_IMAGE, '-I', '-c', RUN_TRACK, format(factor, '.12g'), str(human), 'compare' if compare else 'single']
    print(folder, flush=True)
    started = time.monotonic()
    try:
        with (folder / 'tool.log').open('wb') as log:
            with subprocess.Popen(argv, stdout=log, stderr=subprocess.STDOUT) as process:
                try:
                    while process.poll() is None:
                        if cancelled():
                            raise ValueError('任务已停止')
                        if shutil.disk_usage(folder).free <= DISK_RESERVE:
                            raise ValueError('磁盘剩余接近 50 GB，已停止')
                        if time.monotonic() - started > 1200:
                            raise ValueError('本地轨道生成超时（20 分钟）')
                        time.sleep(1)
                except (ValueError, OSError, KeyboardInterrupt):
                    try:
                        cancel_tracks(run_id)
                    finally:
                        if process.poll() is None:
                            process.kill()
                    raise
                if process.returncode != 0:
                    raise ValueError('生成轨道失败，请查看工具日志')
        track = read_record(folder / 'track.json')
        if track.get('validation') != 'passed' or before != [(p.stat().st_size, p.stat().st_mtime_ns) for p in (bam, bai, bed, chroms)]:
            raise ValueError('结果校验未通过，或本地输入发生变化')
        record['result'] = track
        names = ['coverage.spikein.bw', 'tool.log']
        if compare:
            record['status'] = 'completed_comparison'
            names.append('coverage.cpm.bw')
        else:
            record['status'] = 'completed_spikein'
            (folder / 'report.md').write_text(
                f'# Spike-in 单样本轨道\n\n真实本地计算。使用同轮历史去重人源 BAM 和已核对的 spike-in 片段表。\n\n'
                f'- {count:,} 个历史 spike-in 片段；缩放系数 10000 / {count} = {factor:.12g}。\n'
                '- 10 bp 全基因组轨道，只应用 spike-in 系数。原 CPM 结果保留。\n'
                '- 本次轨道计算未执行有 IgG 对照的找峰。\n'
                '- 未调用模型；这次按可追溯的历史做法重算轨道，不判断实验合格或跨样本可比性。\n', encoding='utf-8')
            names.append('report.md')
        record['artifacts'] = [{'name': n, 'bytes': (folder / n).stat().st_size, 'sha256': sha256(folder / n)}
                               for n in names]
    except (OSError, ValueError, subprocess.SubprocessError, KeyboardInterrupt) as exc:
        record.update(status='failed', error=str(exc) or '用户停止', artifacts=[])
    finally:
        shutil.rmtree(scratch, ignore_errors=True)
        record.update(updated_at=stamp(), elapsed_seconds=round(time.monotonic() - started, 3))
        temp = folder / 'record.tmp'
        temp.write_text(json.dumps(record, ensure_ascii=False, indent=2, allow_nan=False) + '\n')
        temp.replace(folder / 'record.json')
    return record


if __name__ == '__main__':
    try:
        raise SystemExit(main())
    except (ValueError, OSError, KeyError) as exc:
        raise SystemExit(str(exc))
