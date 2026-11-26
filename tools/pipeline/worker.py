"""Single-sample CUT&Tag tools. All argv and mount paths are controlled locally."""
import argparse
import bisect
import hashlib
import importlib.metadata
import itertools
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import time

import numpy as np
import pyBigWig
import pysam

OUT = Path('/output')
SUPPORT = Path('/support')
MACS = '/usr/bin/macs3'
COMMANDS = []


def write(name, value):
    path = OUT / name
    temp = path.with_suffix(path.suffix + '.tmp')
    temp.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + '\n')
    temp.replace(path)


def stage(key, message):
    write('progress.json', {'stage': key, 'message': message})
    print(message, flush=True)


def run(argv, output=None):
    entry = {'argv': argv, 'started': time.time(), 'exit_code': None}
    COMMANDS.append(entry)
    write('commands.json', COMMANDS)
    stream = (OUT / output).open('wb') if output else None
    try:
        code = subprocess.call(argv, stdout=stream)
        entry.update(exit_code=code, seconds=round(time.time() - entry['started'], 3))
        write('commands.json', COMMANDS)
        if code:
            raise RuntimeError('Tool exited with code %s: %s' % (code, argv[0]))
    finally:
        if stream:
            stream.close()


def align(r1, r2):
    argv = ['bowtie2', '--local', '--very-sensitive-local', '--no-mixed', '--no-discordant', '--dovetail',
            '--phred33', '-I', '10', '-X', '1000', '-p', '3', '-x', str(SUPPORT / 'reference/hg38/bowtie2/hg38_noalt'),
            '-1', r1, '-2', r2]
    sorter = ['samtools', 'sort', '-n', '-@', '1', '-m', '384M', '-T', str(OUT / 'work/align'), '-o', str(OUT / 'work/name.bam'), '-']
    entries = [{'argv': a, 'started': time.time(), 'exit_code': None} for a in (argv, sorter)]
    COMMANDS.extend(entries)
    write('commands.json', COMMANDS)
    with (OUT / 'bowtie2.log').open('wb') as log:
        first = subprocess.Popen(argv, stdout=subprocess.PIPE, stderr=log)
        second = subprocess.Popen(sorter, stdin=first.stdout)
        first.stdout.close()
        try:
            codes = [None, second.wait()]
            codes[0] = first.wait()
        finally:
            for proc in (first, second):
                if proc.poll() is None:
                    proc.kill()
                    proc.wait()
        for entry, code in zip(entries, codes):
            entry.update(exit_code=code, seconds=round(time.time() - entry['started'], 3))
        write('commands.json', COMMANDS)
        if any(codes):
            raise RuntimeError('Alignment/sort pipeline failed; see bowtie2.log and tool.log')


def filter_pairs():
    retained, total = 0, 0
    with pysam.AlignmentFile(str(OUT / 'work/name.bam'), 'rb') as inp:
        with pysam.AlignmentFile(str(OUT / 'work/filtered.bam'), 'wb', header=inp.header) as out:
            for _, items in itertools.groupby(inp, lambda a: a.query_name):
                records = list(items)
                total += 1
                if len(records) != 2:
                    continue
                a, b = records
                if any(r.flag & 2828 or not r.is_proper_pair or r.mapping_quality < 30 for r in records):
                    continue
                if a.is_read1 == b.is_read1 or a.reference_id != b.reference_id or a.reference_name in ('chrM', 'chrEBV'):
                    continue
                start, end = min(a.reference_start, b.reference_start), max(a.reference_end, b.reference_end)
                if start < 0 or not 0 < end - start <= 1000:
                    continue
                out.write(a)
                out.write(b)
                retained += 1
    if not retained:
        raise RuntimeError('No fragments retained by the fixed paired-read filter')
    write('pair-filter.json', {'pairs_examined': total, 'pairs_retained': retained, 'min_both_mate_mapq': 30,
                             'excluded_chromosomes': ['chrM', 'chrEBV'], 'max_fragment_length': 1000,
                             'duplicate_policy': 'mark_and_retain'})


def prepare_fastq(sample, *, input_reads=None, pair_limit=None):
    stage('fastp', '正在清理双端 FASTQ，并生成本次 fastp 报告')
    r1, r2 = (str(OUT / f'work/clean.R{mate}.fq.gz') for mate in (1, 2))
    input_r1, input_r2 = input_reads or ('/input/raw/' + sample + '_1.fq.gz', '/input/raw/' + sample + '_2.fq.gz')
    limit = ['--reads_to_process', str(pair_limit)] if pair_limit is not None else []
    run(['fastp', '-i', input_r1, '-I', input_r2,
         '-o', r1, '-O', r2, '--detect_adapter_for_pe', '--length_required', '20', '--thread', '3',
         '--html', str(OUT / 'fastp.html'), '--json', str(OUT / 'fastp.json'), *limit])
    stage('align', '正在将清理后的数据比对到匹配的 hg38 索引')
    align(r1, r2)
    run(['samtools', 'quickcheck', '-v', str(OUT / 'work/name.bam')])
    run(['samtools', 'flagstat', '-O', 'json', str(OUT / 'work/name.bam')], 'alignment.flagstat.json')
    stage('bam', '正在进行双端过滤、重复标记、坐标排序和索引')
    filter_pairs()
    run(['samtools', 'fixmate', '-m', str(OUT / 'work/filtered.bam'), str(OUT / 'work/fixmate.bam')])
    run(['samtools', 'sort', '-@', '2', '-m', '512M', '-T', str(OUT / 'work/coordinate'),
         '-o', str(OUT / 'work/sorted.bam'), str(OUT / 'work/fixmate.bam')])
    bam = str(OUT / 'analysis.bam')
    run(['samtools', 'markdup', '-c', '-s', '-f', str(OUT / 'markdup.txt'), str(OUT / 'work/sorted.bam'), bam])
    run(['samtools', 'index', bam])
    return bam


def sizes():
    with (SUPPORT / 'reference/hg38/hg38.bam-compatible.chrom.sizes').open() as stream:
        return {f[0]: int(f[1]) for line in stream if (f := line.split())}


def peaks(bam, chroms):
    stage('peaks', '正在生成无 IgG 对照的候选窄峰和宽峰；不自动选优')
    return {mode: call_peaks(bam, chroms, mode, mode=mode) for mode in ('narrow', 'broad')}


def call_peaks(bam, chroms, name, *, mode='narrow', controls=(), qvalue=0.01):
    args = [MACS, 'callpeak', '-t', bam, '-f', 'BAMPE', '-g', 'hs', '--keep-dup', 'all',
            '-q', str(qvalue), '-n', name, '--outdir', str(OUT / 'peaks')]
    if mode == 'broad':
        args += ['--broad', '--broad-cutoff', '0.1']
    if controls:
        args += ['-c', *controls]
    run(args)
    extension = 'broadPeak' if mode == 'broad' else 'narrowPeak'
    path = OUT / 'peaks' / (name + '_peaks.' + extension)
    intervals, count = {}, 0
    with path.open() as stream, (OUT / (name + '.peaks.bed')).open('w') as bed:
        for line in stream:
            f = line.rstrip().split('\t')
            chrom, start, end = f[0], int(f[1]), int(f[2])
            if chrom not in chroms or not 0 <= start < end <= chroms[chrom]:
                raise ValueError('Peak BED coordinates do not match the reference')
            bed.write('%s\t%d\t%d\n' % (chrom, start, end))
            intervals.setdefault(chrom, []).append((start, end))
            count += 1
    merged = {}
    for chrom, rows in intervals.items():
        regions = []
        for start, end in sorted(rows):
            if regions and start <= regions[-1][1]:
                regions[-1] = (regions[-1][0], max(end, regions[-1][1]))
            else:
                regions.append((start, end))
        merged[chrom] = {'starts': [a for a, _ in regions], 'ends': [b for _, b in regions]}
    return {'count': count, 'regions': merged}


def fragment_counts(bam, peak_sets):
    histogram = [0] * 1001
    fragments, duplicate_fragments = 0, 0
    hits = dict.fromkeys(peak_sets, 0)
    with pysam.AlignmentFile(bam, 'rb') as inp:
        for r in inp:
            if not r.is_read1 or not r.is_proper_pair or r.flag & 2828:
                continue
            start = min(r.reference_start, r.next_reference_start)
            end = start + abs(r.template_length)
            if start < 0 or end <= start:
                continue
            fragments += 1
            duplicate_fragments += int(r.is_duplicate)
            if end - start <= 1000:
                histogram[end - start] += 1
            for mode, peak_set in peak_sets.items():
                region = peak_set['regions'].get(r.reference_name)
                if region:
                    idx = bisect.bisect_left(region['starts'], end) - 1
                    if idx >= 0 and region['ends'][idx] > start:
                        hits[mode] += 1
    if not fragments:
        raise ValueError('No paired fragments available for metrics')
    return {'fragments': fragments, 'marked_duplicate_fragments': duplicate_fragments,
               'duplicate_fraction': duplicate_fragments / fragments, 'fragment_length_counts': histogram,
               'frip': {mode: {'overlapping_fragments': hits[mode], 'total_fragments': fragments,
                               'fraction': hits[mode] / fragments} for mode in hits},
               'frip_definition': 'paired fragments overlapping at least one candidate peak; duplicate flags retained; no quality threshold'}


def fragment_metrics(bam, peak_sets):
    metrics = fragment_counts(bam, peak_sets)
    histogram = metrics['fragment_length_counts']
    write('fragment-qc.json', metrics)
    with (OUT / 'fragment-length.tsv').open('w') as stream:
        stream.write('length_bp\tfragments\n')
        stream.writelines('%d\t%d\n' % row for row in enumerate(histogram))
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(7, 3.2), layout='constrained')
    ax.plot(range(len(histogram)), histogram, color='#2463a7', linewidth=1)
    ax.set(xlabel='Fragment length (bp)', ylabel='Fragments', xlim=(0, 1000))
    ax.spines[['top', 'right']].set_visible(False)
    fig.savefig(OUT / 'fragment-length.png', dpi=150)
    plt.close(fig)
    return metrics


def coverage(bam, chroms, *, fragments=None):
    stage('coverage', '正在生成全基因组 CPM 覆盖轨道（配对片段，25 bp）')
    scaling = (['--normalizeUsing', 'CPM', '--exactScaling'] if fragments is None
               else ['--scaleFactor', str(1_000_000 / fragments)])
    track = str(OUT / 'coverage.cpm.bw')
    run(['bamCoverage', '-b', bam, '-o', track, '--binSize', '25', *scaling,
         '--extendReads', '--samFlagInclude', '64', '--samFlagExclude', '2828', '--numberOfProcessors', '3'])
    with pyBigWig.open(track) as bw:
        if not bw.isBigWig() or bw.chroms() != chroms:
            raise ValueError('Whole-genome bigWig reference dictionary mismatch')
        header = bw.header()
        if not header.get('nBasesCovered') or not all(math.isfinite(header[k]) for k in ('minVal', 'maxVal', 'sumData', 'sumSquared')):
            raise ValueError('Whole-genome bigWig has no valid finite signal')
        preview = bw.stats('chr1', 1000000, 2000000, nBins=200, type='mean', exact=True)
        if any(v is not None and not math.isfinite(v) for v in preview):
            raise ValueError('Invalid preview signal')
    return {'scope': 'whole_genome', 'reference_count': len(chroms), 'header': header,
            'normalization': 'CPM', 'bin_size': 25, 'paired_fragment_coverage': True,
            'preview': {'chrom': 'chr1', 'start': 1000000, 'end': 2000000, 'points': preview}, 'validation': 'passed'}


def tss():
    stage('deeptools', '正在生成 GENCODE v44 基因级 TSS ±3 kb 矩阵和图；不判定合格线')
    matrix = str(OUT / 'tss.matrix.gz')
    run(['computeMatrix', 'reference-point', '--referencePoint', 'TSS', '-S', str(OUT / 'coverage.cpm.bw'),
         '-R', str(SUPPORT / 'annotation/gencode.v44.genes.bam-compatible.bed'), '-b', '3000', '-a', '3000',
         '--binSize', '50', '--numberOfProcessors', '3', '--sortRegions', 'keep',
         '-o', matrix, '--outFileSortedRegions', str(OUT / 'tss.regions.bed')])
    run(['plotProfile', '-m', matrix, '-out', str(OUT / 'tss.profile.png'),
         '--outFileNameData', str(OUT / 'tss.profile.tsv'), '--samplesLabel', 'CPM', '--plotTitle', 'Gene TSS signal (descriptive)',
         '--yAxisLabel', 'CPM', '--colors', '#2463a7'])
    run(['plotHeatmap', '-m', matrix, '-out', str(OUT / 'tss.heatmap.png'), '--colorMap', 'Blues',
         '--samplesLabel', 'CPM', '--plotTitle', 'Gene TSS signal (descriptive)', '--whatToShow', 'heatmap and colorbar'])
    import gzip
    with gzip.open(OUT / 'tss.matrix.gz', 'rt') as stream:
        header = stream.readline()
        if not header.startswith('@'):
            raise ValueError('deepTools matrix header missing')
        spec = json.loads(header[1:])
        rows = 0
        for line in stream:
            fields = line.split('\t')
            if len(fields) != 126:
                raise ValueError('deepTools matrix dimensions mismatch')
            rows += 1
    if rows != 62700 or spec.get('sample_boundaries') != [0, 120]:
        raise ValueError('Unexpected gene TSS matrix dimensions')
    return {'annotation': 'GENCODE v44 gene-level TSS', 'flank_bp': 3000, 'bin_size': 50,
            'genes': rows, 'columns': 120, 'missing_signal': 'NaN retained', 'quality_verdict': 'not_assessed'}


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for data in iter(lambda: stream.read(4 * 1024 * 1024), b''):
            h.update(data)
    return h.hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('source', choices=['fastq', 'bam', 'resume'])
    parser.add_argument('sample')
    args = parser.parse_args()
    os.chdir(OUT)
    (OUT / 'work').mkdir()
    (OUT / 'peaks').mkdir()
    stage('inputs', '正在检查本地参考与输入')
    versions = {'packages': Path('/opt/pipeline-versions.txt').read_text(),
                'deeptools': importlib.metadata.version('deeptools'), 'pysam': pysam.__version__}
    write('versions.json', versions)
    if args.source == 'fastq':
        bam = prepare_fastq(args.sample)
    elif args.source == 'resume':
        if Path(tempfile.gettempdir()) != OUT / '.scratch':
            raise ValueError('续跑临时目录不可写，请检查任务磁盘挂载')
        bam = '/resume/analysis.bam'
    else:
        bam = '/input/processed/' + args.sample + '.hg38.dedup.bam'
    run(['samtools', 'quickcheck', '-v', bam], 'quickcheck.txt')
    run(['samtools', 'flagstat', '-O', 'json', bam], 'analysis.flagstat.json')
    run(['samtools', 'idxstats', bam], 'idxstats.tsv')
    chroms = sizes()
    with pysam.AlignmentFile(bam, 'rb') as inp:
        if dict(zip(inp.references, inp.lengths)) != chroms:
            raise ValueError('BAM and local reference do not match')
        if args.source == 'resume':
            inp.check_index()
    track = coverage(bam, chroms)
    peak_sets = peaks(bam, chroms)
    stage('qc', '正在统计配对片段、候选峰 FRiP 并检查结果')
    qc = fragment_metrics(bam, peak_sets)
    tss_info = tss()
    stage('validate', '正在核对产物坐标、矩阵和文件摘要')
    qc['peaks'] = {k: {'count': v['count'], 'qvalue': 0.01, 'control': None,
                       'broad_cutoff': 0.1 if k == 'broad' else None} for k, v in peak_sets.items()}
    warnings = ['准确靶点、抗体、分组和生物学重复未确认。', '未使用匹配 IgG，对照相关策略未验证。',
                'spike-in 实验设计未确认，本次只使用 CPM。', '窄峰和宽峰均为探索性候选，不按峰数自动选优。',
                'TSS 图为基因级描述性结果，不作为所有靶点通用质控线。', '未执行 blacklist 过滤或 Motif 富集。']
    if args.source == 'bam':
        warnings += ['起点是已有去重 BAM；没有重新清理或比对 FASTQ；上游指标仍是历史来源。']
    elif args.source == 'resume':
        warnings += ['本轮从前次 FASTQ 分析生成的 BAM 续跑。复用已标记并保留重复片段的 BAM；本轮没有重新清理或比对 FASTQ。']
    else:
        warnings += ['本次重复片段仅标记并保留；与历史去重 BAM 的处理策略不同。']
    summary = {'validation': 'passed', 'source': args.source, 'scope': 'single_sample_basic_pipeline',
               'coverage': track, 'qc': qc, 'tss': tss_info, 'warnings': warnings,
               'scientific_quality': 'not_assessed', 'model_requests': 0, 'versions': versions}
    write('summary.json', summary)
    (OUT / 'report.md').write_text('# 单样本基础分析结果\n\n真实本地工具计算。文件检查通过不等于实验合格。\n\n'
        + '\n'.join('- ' + w for w in warnings) + '\n\n'
        + 'coverage.cpm.bw：全基因组 CPM 配对片段轨道。\n\n'
        + 'narrow.peaks.bed / broad.peaks.bed：无 IgG 的两种候选区域，BED3，供后续 Motif 流程选择。\n\n'
        + 'tss.matrix.gz / tss.profile.png / tss.heatmap.png：GENCODE v44 基因级 TSS ±3 kb。\n\n'
        + 'fragment-qc.json：片段数、重复标记和按片段重叠定义的 FRiP；不设置通用合格线。\n\n'
        + 'versions.json / commands.json / tool.log：工具版本、实际参数和运行日志。\n')
    names = ['coverage.cpm.bw', 'narrow.peaks.bed', 'broad.peaks.bed', 'tss.matrix.gz', 'tss.regions.bed',
             'tss.profile.png', 'tss.profile.tsv', 'tss.heatmap.png', 'fragment-length.png', 'fragment-length.tsv',
             'fragment-qc.json', 'summary.json', 'report.md', 'analysis.flagstat.json', 'idxstats.tsv', 'quickcheck.txt',
             'versions.json', 'commands.json']
    if args.source == 'fastq':
        names += ['fastp.json', 'fastp.html', 'bowtie2.log', 'alignment.flagstat.json', 'pair-filter.json', 'markdup.txt']
    artifacts = []
    for name in names:
        path = OUT / name
        artifacts.append({'name': name, 'bytes': path.stat().st_size, 'sha256': digest(path)})
    write('artifacts.json', {'validation': 'passed', 'artifacts': artifacts})
    stage('done', '基础分析完成；请结合实验背景审阅候选结果')


if __name__ == '__main__':
    try:
        main()
    except Exception as exc:
        write('failure.json', {'error': str(exc), 'tool_executions': len(COMMANDS)})
        raise
