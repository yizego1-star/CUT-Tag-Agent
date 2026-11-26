"""Collect verified, decision-relevant CUT&Tag summaries from local outputs."""
import math
from statistics import median

from .coverage import read_record


SCHEMA = 'cuttag_qc_evidence_v1'
GROUPS = {
    'sequencing_qc': {'raw_reads', 'clean_reads', 'raw_q20', 'clean_q20', 'raw_q30', 'clean_q30',
                      'raw_gc', 'clean_gc', 'clean_read1_mean_bp', 'clean_read2_mean_bp'},
    'pair_filter': {'pairs_examined', 'pairs_retained', 'min_mate_mapq', 'max_fragment_bp'},
    'fragment_qc': {'fragments', 'marked_duplicate_fragments', 'duplicate_fraction',
                    'length_p10_bp', 'length_median_bp', 'length_p90_bp', 'length_mode_bp'},
    'global_track': {'bases_covered', 'signal_min', 'signal_max', 'signal_mean_per_covered_base'},
    'narrow_peaks': {'peaks', 'reference_sequences_with_peaks', 'median_width_bp', 'summed_width_bp',
                     'qvalue', 'overlapping_fragments', 'total_fragments', 'frip'},
    'broad_peaks': {'peaks', 'reference_sequences_with_peaks', 'median_width_bp', 'summed_width_bp',
                    'qvalue', 'broad_cutoff', 'overlapping_fragments', 'total_fragments', 'frip'},
    'tss_profile': {'annotated_genes', 'flank_bp', 'bin_bp', 'center_500bp_mean_cpm', 'outer_500bp_each_end_mean_cpm',
                    'max_mean_cpm'},
    'spikein': {'spikein_fragments', 'human_fragments', 'scale_numerator', 'scale_factor',
                'bin_bp', 'same_alignment_processing_as_cpm', 'same_bin_as_cpm'},
    'experiment_context': {'matched_igg', 'target_expectation', 'biological_replicates', 'spikein_design',
                           'peak_control', 'blacklist_filter', 'reference_assembly', 'tss_annotation'},
}
AUTH_GROUPS = frozenset(GROUPS) | {'tool_parameters', 'execution_status', 'bam_counts', 'bw_regional_statistics',
                                 'local_comparison_statistics'}
CONTEXT = {
    'matched_igg': {'unavailable', 'unknown'},
    'target_expectation': {'unconfirmed'},
    'biological_replicates': {'unconfirmed'},
    'spikein_design': {'unconfirmed'},
    'peak_control': {'not_used'},
    'blacklist_filter': {'not_applied', 'unknown'},
    'reference_assembly': {'hg38'},
    'tss_annotation': {'gencode_v44_gene_level', 'unknown'},
}


def weighted_percentile(counts, fraction):
    rank = max(1, math.ceil(sum(counts) * fraction))
    for index, count in enumerate(counts):
        rank -= count
        if rank <= 0:
            return index
    raise ValueError('片段长度分布为空')


def peak_shape(path):
    widths, references = [], set()
    with path.open() as stream:
        for row in stream:
            chrom, start, end = row.split()
            start, end = int(start), int(end)
            if start < 0 or end <= start:
                raise ValueError('峰坐标范围无效')
            references.add(chrom)
            widths.append(end - start)
    return {'peaks': len(widths), 'reference_sequences_with_peaks': len(references),
            'median_width_bp': median(widths) if widths else 0,
            'summed_width_bp': sum(widths)}


def tss_profile(path, tss):
    with path.open() as stream:
        rows = list(stream)
    if len(rows) != 3:
        raise ValueError('TSS 平均信号表格式无效')
    head, group, *points = rows[2].rstrip('\n').split('\t')
    flank, step = tss['flank_bp'], tss['bin_size']
    bins = len(points)
    if ((head, group) != ('CPM', 'genes') or bins != tss['columns'] or bins * step != 2 * flank
            or step <= 0 or 250 % step or flank < 750):
        raise ValueError('TSS 平均信号表与窗口设置不一致')
    values = [float(value) for value in points]
    if any(math.isinf(value) or value < 0 for value in values):
        raise ValueError('TSS 平均信号值无效')
    center, half, edge = bins // 2, 250 // step, 500 // step

    def mean(window):
        observed = [value for value in window if math.isfinite(value)]
        return sum(observed) / len(observed) if observed else None

    return {'annotated_genes': tss['genes'], 'flank_bp': flank, 'bin_bp': step,
            'center_500bp_mean_cpm': mean(values[center - half:center + half]),
            'outer_500bp_each_end_mean_cpm': mean(values[:edge] + values[-edge:]),
            'max_mean_cpm': max((value for value in values if math.isfinite(value)), default=None)}


def collect(folder, base, artifacts, upstream, spikein, matched_igg):
    """Return narrow structured data and local-only provenance. Never export entire reports."""
    names = {item['name'] for item in artifacts}
    qc = base['qc']
    fragments, marked, lengths = qc['fragments'], qc['marked_duplicate_fragments'], qc['fragment_length_counts']
    if (type(fragments) is not int or fragments <= 0 or type(marked) is not int or not 0 <= marked <= fragments
            or not isinstance(lengths, list) or not lengths or any(type(n) is not int or n < 0 for n in lengths)
            or sum(lengths) != fragments):
        raise ValueError('片段长度分布与质控计数不一致')
    header = base['coverage']['header']
    covered = header['nBasesCovered']
    groups = {
        'fragment_qc': {'fragments': fragments, 'marked_duplicate_fragments': marked, 'duplicate_fraction': marked / fragments,
                        'length_p10_bp': weighted_percentile(lengths, .1), 'length_median_bp': weighted_percentile(lengths, .5),
                        'length_p90_bp': weighted_percentile(lengths, .9), 'length_mode_bp': lengths.index(max(lengths))},
        'global_track': {'bases_covered': covered, 'signal_min': header['minVal'], 'signal_max': header['maxVal'],
                         'signal_mean_per_covered_base': header['sumData'] / covered if covered else None},
        'experiment_context': {'matched_igg': 'unavailable' if matched_igg is False else 'unknown',
                               'target_expectation': 'unconfirmed', 'biological_replicates': 'unconfirmed',
                               'spikein_design': 'unconfirmed', 'peak_control': 'not_used',
                               'reference_assembly': 'hg38',
                               'tss_annotation': 'gencode_v44_gene_level' if base.get('tss', {}).get('annotation') == 'GENCODE v44 gene-level TSS' else 'unknown',
                               'blacklist_filter': 'not_applied' if '未执行 blacklist 过滤或 Motif 富集。' in base.get('warnings', []) else 'unknown'},
    }
    sources = {'fragment_qc': str(folder / 'summary.json'), 'global_track': str(folder / 'coverage.cpm.bw'),
               'experiment_context': '当前样本的实验背景记录与分析配置'}
    for kind in ('narrow', 'broad'):
        name = f'{kind}.peaks.bed'
        shape = peak_shape(folder / name)
        frip, peaks = qc['frip'][kind], qc['peaks'][kind]
        overlap = frip['overlapping_fragments']
        if (shape['peaks'] != peaks['count'] or frip['total_fragments'] != fragments
                or type(overlap) is not int or not 0 <= overlap <= fragments or peaks['control'] is not None):
            raise ValueError('峰表与同一分析 BAM 的 FRiP 记录不匹配')
        group = f'{kind}_peaks'
        groups[group] = {**shape, 'qvalue': peaks['qvalue'], 'overlapping_fragments': overlap,
                         'total_fragments': fragments, 'frip': overlap / fragments}
        if kind == 'broad':
            groups[group]['broad_cutoff'] = peaks['broad_cutoff']
        sources[group] = str(folder / name)
    if 'tss.profile.tsv' in names:
        path = folder / 'tss.profile.tsv'
        groups['tss_profile'] = tss_profile(path, base['tss'])
        sources['tss_profile'] = str(path)
    if 'fastp.json' in upstream:
        path = upstream['fastp.json']
        report = read_record(path)['summary']
        before, after = report['before_filtering'], report['after_filtering']
        groups['sequencing_qc'] = {'raw_reads': before['total_reads'], 'clean_reads': after['total_reads'],
                                  'raw_q20': before['q20_rate'], 'clean_q20': after['q20_rate'],
                                  'raw_q30': before['q30_rate'], 'clean_q30': after['q30_rate'],
                                  'raw_gc': before['gc_content'], 'clean_gc': after['gc_content'],
                                  'clean_read1_mean_bp': after['read1_mean_length'], 'clean_read2_mean_bp': after['read2_mean_length']}
        sources['sequencing_qc'] = str(path)
    if 'pair-filter.json' in upstream:
        path = upstream['pair-filter.json']
        pairs = read_record(path)
        if pairs['pairs_retained'] != fragments:
            raise ValueError('双端过滤记录与分析片段数不匹配')
        groups['pair_filter'] = {'pairs_examined': pairs['pairs_examined'], 'pairs_retained': pairs['pairs_retained'],
                                 'min_mate_mapq': pairs['min_both_mate_mapq'], 'max_fragment_bp': pairs['max_fragment_length']}
        sources['pair_filter'] = str(path)
    if spikein:
        result = spikein['result']
        values = {key: spikein[key] for key in ('spikein_fragments', 'human_fragments', 'scale_numerator', 'scale_factor')}
        if values['spikein_fragments'] <= 0 or not math.isclose(values['scale_numerator'] / values['spikein_fragments'], values['scale_factor']):
            raise ValueError('Spike-in 计数和缩放系数不一致')
        groups['spikein'] = {**values, 'bin_bp': result['bin_size'],
                            'same_alignment_processing_as_cpm': None if base['source'] == 'bam' else False,
                            'same_bin_as_cpm': result['bin_size'] == base['coverage']['bin_size']}
        sources['spikein'] = str(folder.parents[2] / 'outputs' / 'spikein_runs' / spikein['id'] / 'record.json')
    return groups, sources


def validate_groups(groups):
    if not isinstance(groups, dict) or set(groups) - GROUPS.keys() or 'experiment_context' not in groups:
        raise ValueError('质控摘要包含范围外证据')
    for group, values in groups.items():
        if not isinstance(values, dict) or set(values) != GROUPS[group]:
            raise ValueError('质控摘要包含范围外字段')
        for key, value in values.items():
            if group == 'experiment_context':
                valid = isinstance(value, str) and value in CONTEXT[key]
            elif key == 'same_alignment_processing_as_cpm':
                valid = value is None or type(value) is bool
            elif key == 'same_bin_as_cpm':
                valid = type(value) is bool
            elif key in ('center_500bp_mean_cpm', 'outer_500bp_each_end_mean_cpm', 'max_mean_cpm',
                         'signal_mean_per_covered_base') and value is None:
                valid = True
            else:
                valid = type(value) in (int, float) and 0 <= value <= 10 ** 16
            if not valid:
                raise ValueError('质控摘要含无效值')
