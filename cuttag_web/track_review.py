"""Review an existing CPM result without a prescribed decision or a new alignment."""
import json
import math
import subprocess
import time

from .coverage import REGION, read_record
from . import hospital_pair_review as pair_review
from .model import validate_qc_decision, validate_qc_payload, validate_track_review
from .pipeline import IMAGE, PipelineRunner
from .qc_evidence import SCHEMA, collect


READ_REGION = """
import json, sys
import numpy as np
import pyBigWig

size = int(sys.argv[1])
with pyBigWig.open('/track.bw') as track:
    values = track.values('chr1', 1000000, 2000000, numpy=True).reshape(-1, size)
present = ~np.isnan(values)
counts = present.sum(axis=1)
averages = np.nansum(values, axis=1, dtype=float)[counts > 0] / counts[counts > 0]
if not averages.size or not np.isfinite(averages).all():
    raise ValueError('No finite signal in the selected region')
print(json.dumps(dict(interval_count=int(averages.size), nonzero_intervals=int(np.count_nonzero(averages)),
                      preview_mean=float(np.mean(averages)), preview_max=float(np.max(averages)),
                      preview_median=float(np.median(averages)), preview_p95=float(np.percentile(averages, 95))), allow_nan=False))
"""


def read_region(track, size):
    if type(size) is not int or size not in (25, 50, 100):
        raise ValueError('不支持这个区域统计窗口')
    cmd = ['docker', 'run', '--rm', '--pull', 'never', '--network', 'none', '--read-only', '--cap-drop', 'ALL',
           '--security-opt', 'no-new-privileges', '--user', '1000:1000', '--cpus', '1', '--memory', '512m',
           '--mount', f'type=bind,src={track},dst=/track.bw,readonly',
           '--entrypoint', 'python', IMAGE, '-I', '-c', READ_REGION, str(size)]
    process = subprocess.run(cmd, capture_output=True, text=True, timeout=20)
    if process.returncode:
        raise ValueError('本地区域统计读取失败，未重跑测序分析')
    return json.loads(process.stdout)


def flagstat_counts(path, *fields):
    report = read_record(path)
    if not isinstance(report, dict):
        return {}
    groups = [report.get(k, {}) for k in ('QC-passed reads', 'QC-failed reads')]
    counts = {}
    for field in fields:
        values = [group.get(field) for group in groups if isinstance(group, dict)]
        if len(values) == 2 and all(type(value) is int and value >= 0 for value in values):
            counts[field] = sum(values)
    return counts


def original_alignment(folder, source, artifact_names, sample_id=None):
    if source == 'fastq' and 'alignment.flagstat.json' in artifact_names:
        return folder / 'alignment.flagstat.json'
    if source != 'resume':
        return None
    receipt = read_record(folder / 'resume-result.json')
    if not isinstance(receipt, dict) or receipt.get('status') != 'completed_pipeline':
        return None
    execution = receipt.get('execution', {})
    parent = execution.get('resumed_from')
    if execution.get('inputs_unchanged') is not True or not isinstance(parent, str):
        return None
    runner = PipelineRunner(folder.parents[2])
    try:
        original_sample, checkpoint = runner.resume_inputs(parent)
    except (OSError, ValueError):
        return None
    if sample_id is not None and original_sample != sample_id:
        return None
    # Reuse the caller's verified BAM provenance, never historic or other-sample counts.
    inputs = execution.get('inputs', {})
    for path in checkpoint:
        stat = path.stat()
        if inputs.get(str(path.relative_to(runner.root))) != [stat.st_size, stat.st_mtime_ns]:
            return None
    upstream = runner.output_dir(parent)
    commands = read_record(upstream / 'commands.json')
    logged = ['samtools', 'flagstat', '-O', 'json', '/output/work/name.bam']
    if isinstance(commands, list) and any(isinstance(c, dict) and c.get('argv') == logged and c.get('exit_code') == 0 for c in commands):
        path = upstream / 'alignment.flagstat.json'
        if path.resolve() == path and path.is_file():
            return path
    return None


def make_summary(folder, source='resume', *, verified=None, sample_id=None):
    base, artifacts = verified or PipelineRunner.validate(folder, source)
    artifact_names = {a['name'] for a in artifacts}
    if base['coverage']['bin_size'] != 25 or 'analysis.flagstat.json' not in artifact_names:
        raise ValueError('仅使用已验收的 25 bp CPM 覆盖基线与对应 BAM 统计')
    final = flagstat_counts(folder / 'analysis.flagstat.json', 'total', 'mapped', 'duplicates')
    total, mapped = final.get('total', 0), final.get('mapped', 0)
    if not 0 < mapped <= total:
        raise ValueError('BAM 计数记录无效')
    bam_counts = {'mapped_segments': mapped, 'unmapped_segments': total - mapped,
                  'reference_count': base['coverage']['reference_count']}
    if 'duplicates' in final and final['duplicates'] <= total:
        bam_counts['post_filter_duplicate_segments'] = final['duplicates']
    alignment = original_alignment(folder, source, artifact_names, sample_id)
    if alignment:
        original = flagstat_counts(alignment, 'primary', 'primary mapped')
        primary, primary_mapped = original.get('primary', 0), original.get('primary mapped', -1)
        if primary and 0 <= primary_mapped <= primary:
            bam_counts.update(pre_filter_primary_segments=primary, pre_filter_primary_mapped_segments=primary_mapped)
    return {'data_source': 'authorized_minimal_summary',
            'parameters': {'region': REGION, 'bin_size': 25, 'normalization': 'CPM'},
            'execution': {'exit_code': 0, 'output_valid': True},
            'bam_counts': bam_counts,
            'bw_statistics': read_region(folder / 'coverage.cpm.bw', 25)}


def make_qc_evidence(folder, source, sample_id, spikein=None):
    verified = PipelineRunner.validate(folder, source)
    base, artifacts = verified
    names = {row['name'] for row in artifacts}
    core = make_summary(folder, source, verified=verified, sample_id=sample_id)
    parent = original_alignment(folder, source, names, sample_id)
    upstream = {}
    if parent:
        commands = read_record(parent.parent / 'commands.json')
        fastp_logged = isinstance(commands, list) and any(
            isinstance(row, dict) and row.get('exit_code') == 0 and isinstance(row.get('argv'), list)
            and row['argv'] and row['argv'][0] == 'fastp' and '/output/fastp.json' in row['argv'] for row in commands)
        for name in ('fastp.json', 'pair-filter.json'):
            path = parent.parent / name
            if (path.is_file() and path.resolve() == path
                    and ((source == 'fastq' and name in names)
                         or (source == 'resume' and (name == 'pair-filter.json' or fastp_logged)))):
                upstream[name] = path
    if spikein and spikein.get('sample_id') != sample_id:
        raise ValueError('Spike-in 证据与当前分析样本不一致')
    groups, provenance = collect(folder, base, artifacts, upstream, spikein, (spikein or {}).get('matched_igg'))
    if 'sequencing_qc' in groups and groups['sequencing_qc']['clean_reads'] != core['bam_counts'].get('pre_filter_primary_segments'):
        raise ValueError('FASTQ 清理计数与同一次流程的主要比对记录不匹配')
    payload = {'data_source': SCHEMA, 'core': core, 'qc': groups}
    validate_qc_payload(payload)
    provenance.update(bam_counts=[str(folder / 'analysis.flagstat.json')] + ([str(parent)] if parent else []),
                      bw_statistics=str(folder / 'coverage.cpm.bw'))
    return payload, provenance


def compare_evidence(folder, payload, decision):
    action = decision['action']
    if action == 'compare_peak_sets':
        narrow, broad = (payload['qc'][name] for name in ('narrow_peaks', 'broad_peaks'))
        return {'action': action, 'broad_minus_narrow_peaks': broad['peaks'] - narrow['peaks'],
                'broad_minus_narrow_median_bp': broad['median_width_bp'] - narrow['median_width_bp'],
                'broad_minus_narrow_frip': broad['frip'] - narrow['frip']}
    if action == 'compare_resolution':
        return {'action': action, 'bin_bp': decision['bin_size'],
                'statistics': read_region(folder / 'coverage.cpm.bw', decision['bin_size'])}
    return None


def review_qc_baseline(folder, source, sample_id, spikein, model, cancelled, update):
    def request(payload, sent):
        answer = model.review_qc(payload, sample_id, on_request=sent)
        return {**answer, 'decision': validate_qc_decision(answer['decision'], payload)}

    return review_evidence(SCHEMA, lambda: make_qc_evidence(folder, source, sample_id, spikein),
                           lambda: model.qc_review_allowed(sample_id), model.status, request,
                           lambda payload, decision: compare_evidence(folder, payload, decision), cancelled, update)


def review_evidence(schema, load, allowed, status, request, compare, cancelled, update):
    result = {'status': 'reviewing', 'model_requests_attempted': 0, 'full_pipeline_executions': 0,
              'evidence_schema': schema, 'sent_summaries': []}
    started = time.monotonic()

    def finish(status):
        result.update(status=status, elapsed_seconds=round(time.monotonic() - started, 3))
        return result

    def decide(payload):
        def sent():
            if cancelled():
                raise ValueError('任务已停止')
            result['sent_summaries'].append(payload)
            result['model_requests_attempted'] = len(result['sent_summaries'])
            update(result)

        return request(payload, sent)

    try:
        if cancelled():
            return finish('cancelled')
        payload, provenance = load()
        result.update(evidence=payload, evidence_sources=provenance)
        if cancelled():
            return finish('cancelled')
        if not allowed():
            return finish('evidence_ready')
        if not status()['model_connected']:
            raise ValueError('模型连接尚未通过验证')
        result['answer'] = decide(payload)
        update(result)
        if cancelled():
            return finish('cancelled')
        comparison = compare(payload, result['answer']['decision'])
        if comparison:
            result['comparison'] = comparison
            update(result)
            if cancelled():
                return finish('cancelled')
            result['reassessment'] = decide({**payload, 'comparison': comparison})
            update(result)
        return finish('cancelled' if cancelled() else 'completed_review')
    except (OSError, ValueError, KeyError, TypeError, OverflowError, subprocess.SubprocessError) as exc:
        result['error'] = str(exc)
        return finish('cancelled' if cancelled() else 'failed')


def review_adaptive_pair(pair, model, rerun, cancelled, update, *, defer=False, previous=None):
    result = previous or {'evidence_schema': pair_review.SCHEMA, 'sent_summaries': [], 'model_requests_attempted': 0,
                          'status': 'reviewing', 'full_pipeline_executions': 0, 'peak_executions': 0}
    started = time.monotonic()
    elapsed = result.get('elapsed_seconds', 0)

    def finish(status):
        result.update(status=status, elapsed_seconds=round(elapsed + time.monotonic() - started, 3))
        return result

    def decide(payload):
        def sent():
            if cancelled():
                raise InterruptedError('任务已停止')
            result['sent_summaries'].append(payload)
            result['model_requests_attempted'] = len(result['sent_summaries'])
            update(result)

        return model.review_pair(payload, on_request=sent, adaptive=True)

    try:
        if cancelled():
            return finish('cancelled')
        payload = pair_review.adaptive_payload(pair)
        comparison = payload['comparison']
        pair_review.validate_payload(payload, adaptive=True)
        if previous and result.get('evidence') != payload:
            raise ValueError('配对基线与 Agent 做决定时的证据已不一致')
        result.update(status='reviewing', evidence=payload, comparison=comparison)
        if not model.pair_review_allowed():
            return finish('evidence_ready')
        if not previous:
            result['answer'] = decide(payload)
            update(result)
        decision = pair_review.validate_decision(result['answer']['decision'], payload, adaptive=True)
        if cancelled():
            return finish('cancelled')
        if decision['action'] == 'compare_igg':
            result['existing_comparison_selected'] = True
        if decision['action'] != 'tune_igg_qvalue':
            return finish('completed_review')
        if defer:
            return finish('awaiting_rerun')
        result['rerun_started'] = True
        update(result)
        tuned, files = rerun(decision['qvalue'], cancelled)
        evaluated = {**payload, 'comparison': {**comparison, 'rerun': tuned}}
        pair_review.validate_payload(evaluated, adaptive=True)
        result.update(rerun=tuned, artifacts=files, peak_executions=1)
        update(result)
        if cancelled():
            return finish('cancelled')
        result['reassessment'] = decide(evaluated)
        update(result)
        return finish('cancelled' if cancelled() else 'completed_review')
    except (OSError, ValueError, KeyError, TypeError, OverflowError, subprocess.SubprocessError) as exc:
        result['error'] = ('调参工具未在时限内完成' if isinstance(exc, subprocess.TimeoutExpired) else str(exc))
        return finish('cancelled' if cancelled() else 'failed')


def review_baseline(folder, source, model, cancelled=lambda: False, update=lambda _: None, *, sample_id=None, spikein=None):
    if sample_id is not None:
        return review_qc_baseline(folder, source, sample_id, spikein, model, cancelled, update)
    result = {'model_requests_attempted': 0, 'status': 'reviewing', 'full_pipeline_executions': 0}
    started = time.monotonic()

    def sent():
        if cancelled():
            raise ValueError('任务已停止，未发送模型请求')
        result['model_requests_attempted'] = 1
        update(result)

    try:
        if not cancelled():
            summary = make_summary(folder, source)
            result['summary'] = summary
            if not model.status()['model_connected']:
                raise ValueError('当前模型配置未通过连接验证，未发送请求')
            answer = model.review_track(summary, on_request=sent)
            choice = validate_track_review(answer['decision'], 25)
            result['answer'] = answer
            update(result)
            if choice['action'] == 'compare_resolution' and not cancelled():
                compared = read_region(folder / 'coverage.cpm.bw', choice['bin_size'])
                base = summary['bw_statistics']
                same_mean = math.isclose(base['preview_mean'], compared['preview_mean'], rel_tol=1e-6, abs_tol=1e-8)
                result['comparison'] = {'source': 'existing_bigwig_region_means', 'window_bp': choice['bin_size'],
                                        'statistics': compared, 'mean_matches_baseline': same_mean}
        result['status'] = 'cancelled' if cancelled() else 'completed_review'
    except (OSError, ValueError, KeyError, TypeError, subprocess.SubprocessError) as exc:
        result.update(status='cancelled' if cancelled() else 'failed', error=str(exc))
    result['elapsed_seconds'] = round(time.monotonic() - started, 3)
    return result
