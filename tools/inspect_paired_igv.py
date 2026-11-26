"""Compare local peak intervals and CPM tracks for an IGV review of the filtered peaks."""
from bisect import bisect_right
from collections import defaultdict
import json
from math import ceil
from pathlib import Path
from statistics import median
import sys

import pyBigWig


def read_bed(path):
    by_chrom = defaultdict(list)
    with path.open() as stream:
        for line in stream:
            chrom, start, end = line.rstrip().split('\t')
            by_chrom[chrom].append((int(start), int(end)))
    return {chrom: sorted(rows) for chrom, rows in by_chrom.items()}


def merge(rows):
    merged = []
    for start, end in rows:
        if merged and start <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))
    return merged


def track_windows(track, chrom, start, end):
    segments = track.intervals(chrom, start, end) or ()
    return segments, [row[1] for row in segments]


def area(window, start, end):
    segments, ends = window
    count = signal = 0
    for i in range(bisect_right(ends, start), len(segments)):
        left, right, value = segments[i]
        if left >= end:
            break
        overlap = min(right, end) - max(left, start)
        count += overlap
        signal += overlap * value
    return count, signal


def compare_rows(target, igg, original, controlled):
    stats = defaultdict(list)
    for chrom, peaks in original.items():
        length = target.chroms(chrom)
        if length != igg.chroms(chrom) or not length:
            raise ValueError('CPM 轨道缺少匹配的染色体')
        blocks = merge(controlled.get(chrom, ()))
        ends = [end for _, end in blocks]
        for i in range(0, len(peaks), 320):
            batch = peaks[i:i+320]
            low = max(0, min(2 * start - end for start, end in batch))
            high = min(length, max(2 * end - start for start, end in batch))
            tw = track_windows(target, chrom, low, high)
            iw = track_windows(igg, chrom, low, high)
            for start, end in batch:
                if not 0 <= start < end <= length:
                    raise ValueError('候选峰坐标超出参考')
                overlap = bisect_right(ends, start)
                group = 'overlap_igg_controlled' if overlap < len(blocks) and blocks[overlap][0] < end else 'removed'
                t_bp, t_sum = area(tw, start, end)
                i_bp, i_sum = area(iw, start, end)
                row = {'chrom': chrom, 'start': start, 'end': end}
                stats[group].append(row)
                size = end - start
                if t_bp < .9 * size or i_bp < .9 * size or max(t_sum, i_sum) <= 0:
                    continue
                t, g = t_sum / t_bp, i_sum / i_bp
                row.update(target=t, igg=g)
                flank_left, flank_right = max(0, start-size), min(length, end+size)
                l_bp, l_sum = area(tw, flank_left, start)
                r_bp, r_sum = area(tw, end, flank_right)
                if l_bp + r_bp >= .9 * (flank_right - flank_left - size):
                    flank = (l_sum + r_sum) / (l_bp + r_bp) if l_bp + r_bp else 0
                    row['target_focal'] = t > 0 and t >= 2 * flank
    return stats


def label(row, fold):
    if 'target' not in row:
        return 'unscored'
    if row['target'] >= fold * row['igg']:
        return 'target_higher'
    if row['igg'] >= fold * row['target']:
        return 'igg_higher'
    return 'similar'


def summarize(rows):
    ratios = {str(fold): {category: sum(label(row, fold) == category for row in rows)
                          for category in ('target_higher', 'igg_higher', 'similar', 'unscored')}
              for fold in (1.5, 2, 3)}
    scored = [row for row in rows if 'target' in row]
    return {'total': len(rows), 'counts_by_fold': ratios,
            'median_target_cpm': median(row['target'] for row in scored) if scored else None,
            'median_igg_cpm': median(row['igg'] for row in scored) if scored else None,
            'target_higher_focal_2x': sum(label(row, 2) == 'target_higher' and row.get('target_focal', False) for row in rows)}


def examples(stats, target, igg):
    for category in ('target_higher', 'igg_higher', 'similar'):
        rows = [row for row in stats['removed'] if label(row, 2) == category and
                (category != 'target_higher' or row.get('target_focal'))]
        if not rows:
            continue
        rows.sort(key=lambda row: max(row['target'], row['igg']))
        for num, quantile in enumerate((.5, .8), 1):
            row = rows[int(quantile*(len(rows)-1))]
            chrom, start, end = (row[key] for key in ('chrom', 'start', 'end'))
            left, right = max(0, start-900), min(target.chroms(chrom), end+900)
            maximum = max(bw.stats(chrom, left, right, type='max', exact=True)[0] or 0 for bw in (target, igg))
            yield {'label': f'{category}-{num}', 'locus': f'{chrom}:{left+1}-{right}',
                   'max_cpm': max(1, ceil(1.1 * maximum)), 'description': 'original peak without overlap', **row}


def changes_after_rerun(folder, tuned_bed):
    original = read_bed(folder / 'peaks/baseline.peaks.bed')
    controlled = read_bed(folder / 'peaks/with-igg.peaks.bed')
    tuned = {chrom: merge(rows) for chrom, rows in read_bed(tuned_bed).items()}
    ends = {chrom: [end for _, end in rows] for chrom, rows in tuned.items()}
    with pyBigWig.open(str(folder / 'target/coverage.cpm.bw')) as target, pyBigWig.open(str(folder / 'igg/coverage.cpm.bw')) as igg:
        stats = compare_rows(target, igg, original, controlled)

    changes = {}
    overlapping_igg = []
    for source, key, expected_overlap in (('removed', 'removed_now_overlapping', True),
                                          ('overlap_igg_controlled', 'overlapping_now_missing', False)):
        counts = dict.fromkeys(('total', 'target_higher', 'igg_higher', 'similar', 'unscored', 'target_local_peak'), 0)
        for row in stats.get(source, ()):
            chrom = row['chrom']
            blocks = tuned.get(chrom, ())
            i = bisect_right(ends.get(chrom, ()), row['start'])
            overlaps = i < len(blocks) and blocks[i][0] < row['end']
            if overlaps and 'igg' in row:
                overlapping_igg.append(row['igg'])
            if overlaps == expected_overlap:
                category = label(row, 2)
                counts['total'] += 1
                counts[category] += 1
                counts['target_local_peak'] += int(category == 'target_higher' and row.get('target_focal', False))
        changes[key] = counts
    changes['median_igg_cpm_in_rerun_overlap'] = round(median(overlapping_igg), 6) if overlapping_igg else None
    return changes


def inspect(folder):
    original = read_bed(folder / 'peaks/baseline.peaks.bed')
    controlled = read_bed(folder / 'peaks/with-igg.peaks.bed')
    with pyBigWig.open(str(folder / 'target/coverage.cpm.bw')) as target, pyBigWig.open(str(folder / 'igg/coverage.cpm.bw')) as igg:
        stats = compare_rows(target, igg, original, controlled)
        return {'method': 'peak mean CPM with >=90% coverage on both tracks; overlap means at least 1 bp',
                'baseline_peak_count': sum(map(len, original.values())),
                'controlled_peak_count': sum(map(len, controlled.values())),
                'groups': {group: summarize(stats[group]) for group in ('removed', 'overlap_igg_controlled')},
                'igv_examples': list(examples(stats, target, igg))}


if __name__ == '__main__':
    print(json.dumps(inspect(Path(sys.argv[1])), allow_nan=False))
