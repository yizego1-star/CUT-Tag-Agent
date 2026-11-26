"""Render a matched before/after view from the local BigWig and BED files."""
from bisect import bisect_right
import hashlib


FIGURE_NAME = 'peak-comparison.png'
LOCUS_NAME = 'igv-locus.json'
OUTPUT_DIR = 'peak-comparison'
SCRIPT = 'tools/pair_comparison_figure.py'
TRACKS = ('target/coverage.cpm.bw', 'igg/coverage.cpm.bw', 'peaks/baseline.peaks.bed', 'peaks/with-igg.peaks.bed')


def build_comparison_figure(root, run_id, cancelled=lambda: False):
    from cuttag_web.hospital_pair import HospitalPair, local_file, read
    from cuttag_web.hospital_pair_review import QVALUES
    from tools.prepare_hospital_pair import file_identity, rerun_dir, stage

    pair = HospitalPair(root)
    rerun = rerun_dir(root, run_id)
    previous = pair.evidence()[1]['peak_comparison']['settings']['qvalue']
    tuned = read(rerun / 'tuned.json')['qvalue']
    bed = local_file(rerun / 'tuned.peaks.bed')
    if previous != .01 or tuned not in QVALUES or bed.stat().st_size != read(rerun / 'stage.json')['outputs'][bed.name]:
        raise ValueError('找峰参数或 BED 与执行记录不一致')
    aliases = ('target-cpm.bw', 'igg-cpm.bw', 'baseline.peaks.bed', 'with-igg.peaks.bed')
    inputs = [(pair.artifact(alias)[0], '/baseline/' + name) for name, alias in zip(TRACKS, aliases)]
    inputs.append((bed, '/tuned.peaks.bed'))
    extras = (SCRIPT, 'tools/inspect_paired_igv.py')
    receipt = {'image': read(pair.folder / 'full/peaks/stage.json')['inputs']['image'],
               'input_files': {dest: file_identity(path) for path, dest in inputs},
               'qvalues': [previous, tuned], 'scripts': {name: hashlib.sha256((root / name).read_bytes()).hexdigest() for name in extras}}
    stage(rerun / OUTPUT_DIR, ['plot_pair', str(previous), str(tuned)], inputs, (FIGURE_NAME, LOCUS_NAME), receipt,
          cancelled=cancelled, container='cuttag-pair-' + run_id[4:], extra_files=extras, timeout=120)


def comparison_figure_info(root, run_id):
    from cuttag_web.hospital_pair import local_file, read
    from tools.prepare_hospital_pair import rerun_dir

    folder = rerun_dir(root, run_id) / OUTPUT_DIR
    try:
        path = local_file(folder / FIGURE_NAME)
        size = read(folder / 'stage.json')['outputs'][FIGURE_NAME]
        with path.open('rb') as stream:
            if stream.read(8) != b'\x89PNG\r\n\x1a\n' or path.stat().st_size != size:
                return None
    except (OSError, ValueError, KeyError, TypeError):
        return None
    return path, size


def select_windows(stats, tuned, label):
    ends = {chrom: [end for _, end in spans] for chrom, spans in tuned.items()}

    def overlaps(row):
        spans = tuned.get(row['chrom'], ())
        i = bisect_right(ends.get(row['chrom'], ()), row['start'])
        return i < len(spans) and spans[i][0] < row['end']

    restored = [row for row in stats['removed'] if 'target' in row and overlaps(row)]
    earlier = [row for row in stats['overlap_igg_controlled'] if 'target' in row]
    target = [r for r in restored if label(r, 2) == 'target_higher']
    focal = [r for r in target if r.get('target_focal')]
    igg = [r for r in restored if label(r, 2) == 'igg_higher']
    shared = [r for r in restored if label(r, 2) == 'similar']
    lost = [row for row in earlier if not overlaps(row)]
    retained = [row for row in earlier if overlaps(row)]
    groups = [('Target-biased region recovered', focal or target),
              ('IgG-biased region recovered', igg) if igg else ('Similar signals recovered', shared),
              ('Region lost after rerun', lost) if lost else ('Region retained in both runs', retained)]
    for title, rows in groups:
        if rows:
            rows.sort(key=lambda row: (max(row['target'], row['igg']), row['chrom'], row['start']))
            yield title, rows[len(rows)//2]


def display_comparison(folder, tuned_bed, qvalues, output):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch
    from matplotlib.ticker import FuncFormatter, MaxNLocator
    import numpy as np
    import pyBigWig
    from tools.inspect_paired_igv import compare_rows, label, merge, read_bed

    original = read_bed(folder / 'peaks/baseline.peaks.bed')
    before = read_bed(folder / 'peaks/with-igg.peaks.bed')
    after = read_bed(tuned_bed)
    navy, teal, old_color, new_color = '#244c6d', '#218981', '#738095', '#c97526'
    with pyBigWig.open(str(folder / TRACKS[0])) as target, pyBigWig.open(str(folder / TRACKS[1])) as igg:
        stats = compare_rows(target, igg, original, before)
        windows = list(select_windows(stats, {c: merge(rows) for c, rows in after.items()}, label))
        if not windows:
            raise ValueError('两轮峰表中没有可绘制的共同或变化区间')
        figure = plt.figure(figsize=(11, 2.45 * len(windows) + 1.05), facecolor='white')
        grid = figure.add_gridspec(len(windows), 2, left=.095, right=.98, top=.88, bottom=.11,
                                  wspace=.15, hspace=.48)
        loci = []
        for row_index, (title, row) in enumerate(windows):
            chrom, start, end = (row[key] for key in ('chrom', 'start', 'end'))
            flank = max(600, min(1000, (end-start)//2))
            low, high = max(0, start-flank), min(target.chroms(chrom), end+flank)
            loci.append(f'{chrom}:{low+1}-{high}')
            track_values = [np.array(track.values(chrom, low, high), dtype=float) for track in (target, igg)]
            ymax = max(.1, max(np.nanmax(values) for values in track_values)) * 1.08
            x = np.arange(low, high) + 1
            for col, (peak_set, color) in enumerate(((before, old_color), (after, new_color))):
                sub = grid[row_index, col].subgridspec(3, 1, height_ratios=(3, 3, 1.25), hspace=.04)
                axes = [figure.add_subplot(sub[0])]
                axes.extend(figure.add_subplot(sub[i], sharex=axes[0]) for i in (1, 2))
                for ax in axes:
                    ax.axvspan(start+1, end, color='#eaf0f3', lw=0, zorder=0)
                    ax.spines[['top', 'right']].set_visible(False)
                    ax.spines[['left', 'bottom']].set_color('#d5dde4')
                    ax.set_xlim(low+1, high)
                for ax, values, stroke, name in zip(axes, track_values, (navy, teal), ('Target', 'IgG')):
                    ax.fill_between(x, 0, values, color=stroke, step='post', linewidth=0)
                    ax.set_ylim(0, ymax)
                    ax.grid(axis='y', color='#eef0f2', lw=.6)
                    ax.set_axisbelow(True)
                    ax.yaxis.set_major_locator(MaxNLocator(3))
                    ax.tick_params(axis='y', labelsize=7, length=2, colors='#48586b')
                    ax.tick_params(axis='x', bottom=False, labelbottom=False)
                    if col == 0:
                        ax.set_ylabel(name, fontsize=8, color=stroke, rotation=0, labelpad=21, va='center')
                    else:
                        ax.tick_params(axis='y', left=False, labelleft=False)
                axes[0].set_title(title if col == 0 else f'{chrom}:{low+1:,}–{high:,}', fontsize=9,
                                  fontweight='semibold', loc='left', pad=4, color='#23354b')
                peaks = [(max(low+1,a+1), min(b,high)-max(low+1,a+1)+1)
                         for a,b in peak_set.get(chrom, ()) if b > low and a < high]
                axes[-1].broken_barh(peaks, (.25,.5), facecolors=color)
                axes[-1].set_ylim(0,1)
                axes[-1].set_yticks([])
                axes[-1].xaxis.set_major_locator(MaxNLocator(4, prune='both'))
                axes[-1].xaxis.set_major_formatter(FuncFormatter(lambda value, _: f'{int(value):,}'))
                axes[-1].tick_params(axis='x', labelsize=7, length=2, colors='#48586b')
                if col == 0:
                    axes[-1].set_ylabel('Peaks', fontsize=8, color='#48586b', rotation=0, labelpad=21, va='center')
        for column, title in enumerate(('Before', 'After')):
            figure.text(.095 + .475 * column, .946, f'{title}   q = {qvalues[column]:g}', fontsize=13,
                        fontweight='semibold', color='#23354b')
        figure.legend(handles=[Patch(color=navy, label='Target CPM'), Patch(color=teal, label='IgG CPM'),
                               Patch(color='#eaf0f3', label='Reviewed interval'), Patch(color=old_color, label='Before peaks'),
                               Patch(color=new_color, label='After peaks')], loc='lower center', ncol=5, frameon=False, fontsize=8)
        figure.savefig(output, dpi=180, facecolor='white')
        plt.close(figure)
        return loci
