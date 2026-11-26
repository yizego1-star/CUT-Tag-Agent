"""Container entry point for the public H9 development workflow."""
import argparse
import json
import math
import os
from pathlib import Path

import pyBigWig
import pysam

from tools.pipeline import worker
from tools.h9_samples import CONTROLS, SAMPLES, TARGETS


REVIEW_REGION = 'chr3:181649967:181792713'


def check_bam(path, chroms):
    worker.run(['samtools', 'quickcheck', '-v', path])
    with pysam.AlignmentFile(path, 'rb') as bam:
        if dict(zip(bam.references, bam.lengths)) != chroms or not bam.check_index():
            raise ValueError('H9 BAM 参考或索引不匹配')


def prepare(sample, pair_limit, reads_root=Path('/reads')):
    out = worker.OUT
    (out / 'work').mkdir()
    reads = tuple(str(reads_root / f'{sample}_{mate}.fastq.gz') for mate in (1, 2))
    bam = worker.prepare_fastq(sample, input_reads=reads, pair_limit=pair_limit)
    chroms = worker.sizes()
    check_bam(bam, chroms)
    qc = json.loads((out / 'fastp.json').read_text())['summary']
    fragments = json.loads((out / 'pair-filter.json').read_text())['pairs_retained']
    reads_examined = qc['before_filtering']['total_reads']
    if pair_limit is not None and reads_examined != 2 * pair_limit:
        raise ValueError('限量联调实际读取的 FASTQ 数量偏离配置')
    worker.write('prepared.json', {'sample': sample, 'pair_limit': pair_limit, 'raw_read_pairs': reads_examined // 2,
                                  'clean_read_pairs': qc['after_filtering']['total_reads'] // 2,
                                  'human_fragments': fragments})
    worker.stage('done', '已保存可复用的 BAM 和索引')


def track(sample_root=Path('/sample')):
    info = json.loads((sample_root / 'prepared.json').read_text())
    chroms = worker.sizes()
    bam = str(sample_root / 'analysis.bam')
    check_bam(bam, chroms)
    fragments = info['human_fragments']
    # Preparation retains one passing R1 per fragment, including marked duplicates.
    with pysam.AlignmentFile(bam, 'rb') as inp:
        if fragments <= 0 or inp.mapped != 2 * fragments or inp.unmapped:
            raise ValueError('H9 BAM 索引计数与 CPM 片段分母不符')
    worker.stage('coverage', '正在生成 H9 CPM 轨道')
    track_path = worker.OUT / 'coverage.cpm.bw'
    if info['pair_limit'] is None:
        worker.coverage(bam, chroms, fragments=fragments)
    else:
        # The prepared BAM retains exactly one passing R1 for each retained pair.
        worker.run(['bamCoverage', '-b', bam, '-o', str(track_path), '--region', REVIEW_REGION,
                    '--binSize', '25', '--scaleFactor', str(1_000_000 / fragments), '--extendReads',
                    '--samFlagInclude', '64', '--samFlagExclude', '2828', '--numberOfProcessors', '3'])
        with pyBigWig.open(str(track_path)) as bw:
            if bw.chroms().get('chr3') != chroms['chr3'] or not math.isfinite(bw.header()['sumData']):
                raise ValueError('H9 预览轨道校验失败')
    worker.write('track.json', {'track': track_path.name, 'normalization': 'CPM', 'fragments': fragments,
                               'scope': 'whole_genome' if info['pair_limit'] is None else REVIEW_REGION})
    worker.stage('done', '已保存轨道')


def replicate_jaccard(first, second):
    intersection = total = 0
    for chrom in first.keys() | second.keys():
        a, b = (list(zip(block['starts'], block['ends'])) for block in
                (value.get(chrom, {'starts': [], 'ends': []}) for value in (first, second)))
        total += sum(end - start for start, end in a + b)
        i = j = 0
        while i < len(a) and j < len(b):
            x, y = a[i]
            u, v = b[j]
            intersection += max(0, min(y, v) - max(x, u))
            if y <= v:
                i += 1
            else:
                j += 1
    return intersection / (total - intersection) if total > intersection else None


def analyze(use_igg, samples_root=Path('/samples')):
    (worker.OUT / 'peaks').mkdir()
    chroms = worker.sizes()
    bams = {sample: str(samples_root / sample / 'analysis.bam') for sample in SAMPLES}
    for bam in bams.values():
        check_bam(bam, chroms)
    worker.stage('peaks', '正在使用两份 H9 IgG 重新找窄峰' if use_igg else '正在生成 H9 候选窄峰基线')
    controls = [bams[sample] for sample in CONTROLS] if use_igg else []
    peaks = {sample: worker.call_peaks(bams[sample], chroms, sample, controls=controls) for sample in TARGETS}
    worker.stage('qc', '正在计算目标和 IgG 落峰比例及重复峰重叠')
    counts = {sample: worker.fragment_counts(bam, peaks)['frip'] for sample, bam in bams.items()}
    worker.write('metrics.json', {
        'peak_mode': 'narrow', 'qvalue': 0.01, 'controls': list(CONTROLS) if use_igg else [],
        'peak_counts': {sample: peaks[sample]['count'] for sample in TARGETS},
        'target_frip': {sample: counts[sample][sample] for sample in TARGETS},
        'igg_overlap': {sample: counts[sample] for sample in CONTROLS},
        'replicate_peak_jaccard': replicate_jaccard(*(peaks[sample]['regions'] for sample in TARGETS)),
    })
    worker.stage('done', '候选窄峰和证据已保存')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--native-root', type=Path)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--input', type=Path)
    sub = parser.add_subparsers(dest='command', required=True)
    prep = sub.add_parser('prepare')
    prep.add_argument('sample', choices=SAMPLES)
    prep.add_argument('--pair-limit', type=int)
    sub.add_parser('track')
    sub.add_parser('baseline')
    sub.add_parser('compare-igg')
    args = parser.parse_args()
    if args.native_root:
        root = args.native_root
        if (not args.output or not args.input or root.resolve() != root or
                any(path.resolve() != path or root not in path.parents for path in (args.output, args.input))):
            parser.error('原生输入输出路径必须位于项目目录内')
        worker.OUT = args.output
        worker.SUPPORT = root
        worker.MACS = str(Path(os.path.abspath(os.sys.executable)).parent / 'macs3')
    os.chdir(worker.OUT)
    if args.command == 'prepare':
        prepare(args.sample, args.pair_limit, args.input or Path('/reads'))
    elif args.command == 'track':
        track(args.input or Path('/sample'))
    else:
        analyze(args.command == 'compare-igg', args.input or Path('/samples'))


if __name__ == '__main__':
    main()
